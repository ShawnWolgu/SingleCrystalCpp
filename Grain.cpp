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

Matrix6d Grain::get_dp_grad(Matrix3d stress_3d){
    Matrix6d dp_grad = Matrix6d::Zero();
    Matrix3d stress_grain = orientation * stress_3d * orientation.transpose();
    for (Slip &slip_component : slip_sys) {
        slip_component.cal_ddgamma_dtau(*this, stress_grain);
        dp_grad += slip_component.ddgamma_dsigma();
    }
    return rotate_6d_stiff_modu(dp_grad, orientation);
}

void Grain::update_status(Matrix3d L_dt_tensor, Matrix3d vel_grad_flag, Matrix3d stress_incr, Matrix3d dstress_flag){
    Matrix3d vel_grad_elas = Matrix3d::Zero(), vel_grad_plas = Matrix3d::Zero(), spin_elas = Matrix3d::Zero();
    // update strain_rate
    if (L_dt_tensor != Matrix3d::Zero()) strain_rate = L_dt_tensor.cwiseAbs().maxCoeff() / timestep;
    // update elastic modulus
    elastic_modulus = rotate_6d_stiff_modu(elastic_modulus_ref,orientation);
    solve_Lsig_iteration(L_dt_tensor, vel_grad_flag, stress_incr, dstress_flag);
    L_dt_tensor(0,1) = 0, L_dt_tensor(1,0) = 0; // fix Omega12 = 0
    vel_grad_plas = get_vel_grad_plas(stress_tensor + stress_incr);
    vel_grad_elas = L_dt_tensor - vel_grad_plas;
    // end iteration

    stress_tensor = stress_tensor + stress_incr;
    strain_tensor = strain_tensor + 0.5 * (vel_grad_elas + vel_grad_plas + vel_grad_elas.transpose() + vel_grad_plas.transpose());
    deform_grad_elas = (vel_grad_elas + Matrix3d::Identity()) * deform_grad_elas;
    deform_grad_plas = deform_grad_elas.inverse() * deform_grad;
    for (Slip &slip_component : slip_sys) slip_component.update_status(*this);
    spin_elas = 0.5 * (vel_grad_elas - vel_grad_elas.transpose());
    orientation = orientation * Rodrigues(spin_elas).transpose(); 
}

void Grain::solve_Lsig_iteration(Matrix3d &L_dt_tensor, Matrix3d &vel_grad_flag, Matrix3d &stress_incr, Matrix3d &dstress_flag){
    Matrix<double, 15, 1> params, flags;
    params << tensor_trans_order_9(L_dt_tensor), tensor_trans_order(stress_incr);
    flags << tensor_trans_order_9(vel_grad_flag), tensor_trans_order(dstress_flag);
    vector<int> unknown_idx(6,0), known_idx(9,0);
    flag_to_idx(flags, known_idx, unknown_idx);
    Matrix<double,9,1> known_params = params(known_idx);
    Vector6d unknown_params = params(unknown_idx);
    
    int iter_num = 0;
    double coeff = 0.5;
    Vector6d x_iter_save = Vector6d::Zero(), x_iter_step = Vector6d::Zero();
    Vector6d y_vec = Vector6d::Zero(), y_vec_last = Vector6d::Zero();
    do{
        x_iter_save = unknown_params;
        //Solve stress_incr, L_dt_tensor
        ddp_by_dsigma = get_dp_grad(stress_tensor+stress_incr);
        y_vec = calcfx_and_update(unknown_params, known_params, unknown_idx, known_idx, right_vec(stress_incr));
        
        
        //cout << "coeff = " << coeff << endl;
        //cout << "x = " << endl << x_iter_save.transpose() << endl;
        //cout << "y = " << endl << y_vec.transpose() << endl;
        //cout << "y_norm = " << endl << y_vec.norm() << endl;
        //cout << "y_last = " << endl << y_vec_last.transpose() << endl;
        //cout << "y*y_last = " << endl << y_vec.dot(y_vec_last) << endl;

        //if (y_vec_last != Vector6d::Zero() && y_vec.norm()/y_vec_last.norm() > 1){
        if (y_vec.dot(y_vec_last) < -1e-20){
            coeff = 0.5 * coeff;
            unknown_params = x_iter_save - coeff * x_iter_step;
            continue;
        }
        else{
            coeff = 1.2 * coeff;
            x_iter_step = unknown_params - x_iter_save;
            y_vec_last = y_vec;
            unknown_params = x_iter_save + coeff * x_iter_step;
            int unk_par_idx = 0;
            for (auto &unk_idx : unknown_idx) {params(unk_idx) = unknown_params(unk_par_idx); ++unk_par_idx;}
            params_convert_to_matrix(params, L_dt_tensor, stress_incr);
        }
        ++ iter_num;
        if(iter_num > 10000) {cout << (unknown_params - x_iter_save).norm() << endl << "Not converged... but still going" << endl; break;}
    } while ((unknown_params - x_iter_save).norm() > 1e-5);
    if (abs(y_vec.norm()) > 1e-1){
        cout << "End-of-step y_vec, coeff = " << endl;
        cout << y_vec.norm() << "  " << coeff << endl;
    }
}

Vector6d Grain::calcfx_and_update(Vector6d &unknown_params, Matrix<double,9,1> known_params, vector<int> unknown_idx,vector<int> known_idx, Vector6d right_vector){
    Matrix<double, 6, 15> left_x_mid = left_matrix() * mid_matrix();
    Matrix6d solveA;
    Matrix<double,6,9> solveB;
    solveA << left_x_mid(all,unknown_idx);
    solveB << left_x_mid(all,known_idx);
    Vector6d return_vec = solveA * unknown_params + solveB * known_params - right_vector;
    unknown_params = solveA.inverse() * (-1*solveB*known_params + right_vector);
    return return_vec;
}

Matrix<double, 6, 15> Grain::left_matrix(){
    Vector6d stress_6d = tensor_trans_order(stress_tensor);
    Matrix6d C_ij_pri = get_C_ij_pri(stress_6d);
    Sigma_ik = get_Sigma_ik(stress_6d);

    Matrix6d minus_delta_il = -Matrix6d::Identity(); 
    CNPN = elastic_modulus * strain_modi_tensor * ddp_by_dsigma * strain_modi_tensor;

    Matrix<double, 6, 15> left_M;
    left_M << C_ij_pri, Sigma_ik, minus_delta_il-CNPN;
    return left_M;
}

Matrix<double, 15, 15> Grain::mid_matrix(){
    Matrix<double, 15, 15> mid_M;
    mid_M << Matrix3d::Identity(), Matrix3d::Zero(), Matrix3d::Zero(), Matrix<double,3,6>::Zero(),
             Matrix3d::Zero(), Matrix3d::Identity(), Matrix3d::Identity(), Matrix<double,3,6>::Zero(),
             Matrix3d::Zero(), 0.5*Matrix3d::Identity(), -0.5*Matrix3d::Identity(), Matrix<double,3,6>::Zero(),
             Matrix<double,6,3>::Zero(), Matrix<double,6,6>::Zero(), Matrix6d::Identity(); 
    return mid_M;
}

Vector6d Grain::right_vec(Matrix3d &stress_incr){
    vector<int> slice = {3,4,5};
    Matrix3d dplusw = get_vel_grad_plas(stress_tensor+stress_incr);
    Vector6d dpn = tensor_trans_order((Matrix3d)(0.5*dplusw + 0.5*dplusw.transpose()));
    Vector3d dwn = tensor_trans_order((Matrix3d)(0.5*dplusw - 0.5*dplusw.transpose()))(slice);
    Vector6d result = elastic_modulus * strain_modi_tensor * dpn - CNPN * tensor_trans_order(stress_incr) + Sigma_ik * dwn;
    return result;
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
