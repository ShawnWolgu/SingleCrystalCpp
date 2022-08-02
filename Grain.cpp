#include "singleX.h"

Grain::Grain(): deform_grad(Matrix3d::Identity()), deform_grad_elas(deform_grad), deform_grad_plas(Matrix3d::Zero()), stress_tensor(Matrix3d::Zero()),
            strain_tensor(Matrix3d::Zero()), orientation(Matrix3d::Identity()), elastic_modulus(Matrix6d::Identity()), elastic_modulus_ref(Matrix6d::Identity()) {};

Grain::Grain(Matrix6d elastic_mod, Matrix3d lat_vecs, vector<Slip> s, Matrix3d orient_Mat){
    deform_grad = Matrix3d::Identity(), deform_grad_elas = deform_grad, deform_grad_plas = stress_tensor = Matrix3d::Zero(),
    strain_tensor = Matrix3d::Zero(), orientation = orient_Mat, elastic_modulus_ref = elastic_mod,  elastic_modulus = rotate_6d_stiff_modu(elastic_mod,orientation);
    lattice_vec = lat_vecs;
    for (Slip &slip_component : s) slip_sys.push_back(slip_component); 
    for (Slip &slip_component : slip_sys) slip_component.update_status(*this);
    strain_modi_tensor = Matrix6d::Zero();
    strain_modi_tensor.diagonal() << 1,1,1,2,2,2;
}

Matrix3d Grain::get_vel_grad_plas(Matrix3d stress_3d){
    Matrix3d vel_grad_plas = Matrix3d::Zero();
    Matrix3d stress_grain = orientation * stress_3d * orientation.transpose();
    for (Slip &slip_component : slip_sys) {
        slip_component.cal_strain(*this, stress_grain);
        vel_grad_plas += slip_component.dL_tensor();
    }
    return orientation.transpose() * vel_grad_plas * orientation;
}

void Grain::calc_slip_ddgamma_dtau(Matrix3d stress_3d){
    Matrix3d stress_grain = orientation * stress_3d * orientation.transpose();
    for (Slip &slip_component : slip_sys) {
        slip_component.cal_ddgamma_dtau(*this, stress_grain);
    }
}

Matrix6d Grain::get_dp_grad(){
    Matrix6d dp_grad = Matrix6d::Zero();
    for (Slip &slip_component : slip_sys) {
        dp_grad += slip_component.ddp_dsigma();
    }
    return rotate_6d_stiff_modu(dp_grad, orientation)*strain_modi_tensor;
}

Matrix<double, 3, 6> Grain::get_wp_grad(){
    Matrix6d wp_grad = Matrix6d::Zero();
    for (Slip &slip_component : slip_sys) {
        wp_grad += slip_component.dwp_dsigma();
    }
    Matrix6d return_mat = rotate_6d_stiff_modu(wp_grad, orientation)*strain_modi_tensor;
    return return_mat(Eigen::seq(3,last),all);
}

void Grain::update_status(Matrix3d vel_bc_tensor, Matrix3d vel_grad_flag, Matrix3d stress_incr, Matrix3d dstress_flag){
    Matrix3d vel_grad_elas = Matrix3d::Zero(), vel_grad_plas = Matrix3d::Zero(), spin_elas = Matrix3d::Zero();
    // update strain_rate
    if (vel_bc_tensor != Matrix3d::Zero()) strain_rate = vel_bc_tensor.cwiseAbs().maxCoeff() / timestep;
    // update elastic modulus
    elastic_modulus = rotate_6d_stiff_modu(elastic_modulus_ref,orientation);
    solve_Lsig_iteration(vel_bc_tensor, vel_grad_flag, stress_incr, dstress_flag);
    vel_grad_plas = get_vel_grad_plas(stress_tensor + stress_incr);
    vel_grad_elas = vel_bc_to_vel_grad(vel_bc_tensor) - vel_grad_plas;
    // end iteration
    stress_tensor = stress_tensor + stress_incr;
    strain_tensor = strain_tensor + 0.5 * (vel_grad_elas + vel_grad_plas + vel_grad_elas.transpose() + vel_grad_plas.transpose());
    deform_grad = (vel_bc_to_vel_grad(vel_bc_tensor)+Matrix3d::Identity()) * deform_grad;
    deform_grad_elas = (vel_grad_elas + Matrix3d::Identity()) * deform_grad_elas;
    deform_grad_plas = deform_grad_elas.inverse() * deform_grad;
    for (Slip &slip_component : slip_sys) slip_component.update_status(*this);
    spin_elas = 0.5 * (vel_grad_elas - vel_grad_elas.transpose());
    orientation = orientation * Rodrigues(spin_elas).transpose(); 
}

void Grain::solve_Lsig_iteration(Matrix3d &vel_bc_tensor, Matrix3d &vel_grad_flag, Matrix3d &stress_incr, Matrix3d &dstress_flag){
    Matrix<double, 15, 1> params, flags;
    params << tensor_trans_order_9(vel_bc_tensor), tensor_trans_order(stress_incr);
    flags << tensor_trans_order_9(vel_grad_flag), tensor_trans_order(dstress_flag);
    vector<int> unknown_idx(6,0), known_idx(9,0);
    flag_to_idx(flags, known_idx, unknown_idx);
    Matrix<double,9,1> known_params = params(known_idx);
    Vector6d unknown_params = params(unknown_idx), stress_6d = tensor_trans_order(stress_tensor), y_vec = Vector6d::Zero();
    C_ij_pri = get_C_ij_pri(stress_6d);
    Sigma_ik = get_Sigma_ik(stress_6d);
    double coeff = 1/1.2;

    do{
	coeff = 1/1.2;	y_vec = Vector6d::Zero();
   	Vector6d x_iter_save = Vector6d::Zero(), x_iter_step = Vector6d::Zero(), y_vec_last = Vector6d::Zero();
   	Matrix6d dfx_matrix;
        do{
            x_iter_save = unknown_params;
            y_vec = calc_fx(vel_bc_tensor, stress_incr);
            if (y_vec_last != Vector6d::Zero() && y_vec.norm()/y_vec_last.norm() > 1){
                coeff = 0.5 * coeff;
                unknown_params = x_iter_save - coeff * x_iter_step;
            }
            else{
                coeff = 1.2 * coeff;
    		dfx_matrix = calc_dfx(vel_bc_tensor, stress_incr, unknown_idx);
    	    	x_iter_step = -dfx_matrix.inverse()*y_vec;
                y_vec_last = y_vec;
                unknown_params = x_iter_save + coeff * x_iter_step;
            }
	    params_convert_to_matrix(params, unknown_params, unknown_idx, vel_bc_tensor, stress_incr);
        } while ((unknown_params - x_iter_save).norm() > 1e-7);
        if (abs(y_vec.norm()) > 1){
    	    if (dtime < 1e-20) {cout << "Severe disconvergence!! break!" << endl; exit(0);}
    	    else{
    	        //cout << "Severe disconvergence!! Will retry in a smaller step." << endl;
    	        dtime = dtime /2;
    	    }
        }
    } while(abs(y_vec.norm()) > 1);
    if (abs(y_vec.norm()) > 1e-1){
	cout << "End-of-step y_vec, coeff = " << endl;
        cout << y_vec.norm() << "  " << coeff << endl;
    }
}

Vector6d Grain::calc_fx(Matrix3d &vel_bc_tensor, Matrix3d &stress_incr){
    Matrix<double, 6, 9> left_M;
    left_M << C_ij_pri*strain_modi_tensor, Sigma_ik;
    Matrix3d vel_grad_elas = vel_bc_tensor - get_vel_grad_plas(stress_tensor + stress_incr);
    Vector6d return_vec = left_M * bc_modi_matrix * tensor_trans_order_9(vel_grad_elas) - tensor_trans_order(stress_incr);
    return return_vec;
}

Matrix6d Grain::calc_dfx(Matrix3d &vel_bc_tensor, Matrix3d &stress_incr,vector<int> unknown_idx){
    calc_slip_ddgamma_dtau(stress_tensor+stress_incr);
    dwp_by_dsigma = get_wp_grad();
    ddp_by_dsigma = get_dp_grad();
    Matrix6d dfx_ddsigma = -Matrix6d::Identity() - elastic_modulus*strain_modi_tensor*ddp_by_dsigma - Sigma_ik*dwp_by_dsigma;
    Matrix<double, 6, 15> dfx_dX;
    Matrix<double,15,15> trans_compo;
    dfx_dX << C_ij_pri*strain_modi_tensor, Sigma_ik, dfx_ddsigma;
    trans_compo << bc_modi_matrix, Matrix<double,9,6>::Zero(), Matrix<double,6,9>::Zero(), Matrix6d::Identity();
    dfx_dX = dfx_dX * trans_compo;
    return dfx_dX(all, unknown_idx);
}

Matrix<double, 9, 9> Grain::L_compo_to_wd(){
    Matrix<double, 9, 9> mid_M;
    mid_M << Matrix3d::Identity(), Matrix3d::Zero(), Matrix3d::Zero(),
             Matrix3d::Zero(), 0.5*Matrix3d::Identity(), 0.5*Matrix3d::Identity(),
             Matrix3d::Zero(), 0.5*Matrix3d::Identity(), -0.5*Matrix3d::Identity();
    return mid_M;
}

Matrix6d Grain::get_C_ij_pri(Vector6d &stress_6d){
    Matrix6d C_ij_pri;
    C_ij_pri << stress_6d, stress_6d, stress_6d, Vector6d::Zero(), Vector6d::Zero(), Vector6d::Zero();
    C_ij_pri = elastic_modulus + C_ij_pri;
    return C_ij_pri;
}

Matrix<double, 6, 3> Grain::get_Sigma_ik(Vector6d &stress_6d){
    Matrix<double, 6, 3> Sigma_ik;
    Sigma_ik << 0, 2*stress_6d(4), 2*stress_6d(5), 2*stress_6d(3), 0, -2*stress_6d(5),
                -2*stress_6d(3), -2*stress_6d(4), 0, stress_6d(2)-stress_6d(1), -stress_6d(5), -stress_6d(4),
                -stress_6d(5), stress_6d(2)-stress_6d(0), stress_6d(3), stress_6d(4), stress_6d(3), stress_6d(1)-stress_6d(0);
    return Sigma_ik; 
}

void Grain::print_stress_strain(ofstream &os){
    os << strain_tensor(0,0) << ',' << strain_tensor(1,1) << ','  << strain_tensor(2,2) << ',' << strain_tensor(0,1) << ',' 
       << strain_tensor(0,2) << ',' << strain_tensor(1,2) << ','  << stress_tensor(0,0) << ',' << stress_tensor(1,1) << ',' 
       << stress_tensor(2,2) << ',' << stress_tensor(0,1) << ','  << stress_tensor(0,2) << ',' << stress_tensor(1,2)  << endl;
}

void Grain::print_stress_strain_screen(){
    cout << strain_tensor(0,0) << '\t' << strain_tensor(1,1) << '\t'  << strain_tensor(2,2) << '\t' << strain_tensor(0,1) << '\t' 
       << strain_tensor(0,2) << '\t' << strain_tensor(1,2) << '\t'  << stress_tensor(0,0) << '\t' << stress_tensor(1,1) << '\t' 
       << stress_tensor(2,2) << '\t' << stress_tensor(0,1) << '\t'  << stress_tensor(0,2) << '\t' << stress_tensor(1,2) << '\t'  << endl;
}

void Grain::print_dislocation(ofstream &os){
    os << strain_tensor(0,0) << ',' << strain_tensor(1,1) << ','  << strain_tensor(2,2);
    for (Slip &slip_component : slip_sys) os << ',' << slip_component.SSD_density;
    os << endl;
}

void Grain::print_crss(ofstream &os){
    os << strain_tensor(0,0) << ',' << strain_tensor(1,1) << ','  << strain_tensor(2,2) << ',' << (slip_sys)[0].acc_strain;
    for (Slip &slip_component : slip_sys) os << ',' << slip_component.crss;
    os << endl;
}

void Grain::print_accstrain(ofstream &os){
    os << strain_tensor(0,0) << ',' << strain_tensor(1,1) << ','  << strain_tensor(2,2);
    for (Slip &slip_component : slip_sys) os << ',' << slip_component.acc_strain;
    os << endl;
}

void Grain::print_euler(ofstream &os){
    Vector3d euler_angle = Euler_trans(orientation);
    os << euler_angle(0) << ',' << euler_angle(1) << ',' << euler_angle(2);
    os << endl;
}
