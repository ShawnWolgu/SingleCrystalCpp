#include "singleX.h"

Grain::Grain(): deform_grad(Matrix3d::Identity()), deform_grad_elas(deform_grad), deform_grad_plas(Matrix3d::Zero()), stress_tensor(Matrix3d::Zero()),
            strain_tensor(Matrix3d::Zero()), orientation(Matrix3d::Identity()), elastic_modulus(Matrix6d::Identity()) {};
    
void Grain::set_slip_sys(vector<Slip> &s){
    slip_sys = &s;
}

Matrix3d Grain::get_vel_grad_plas(Matrix3d stress_3d){
    Matrix3d vel_grad_plas = Matrix3d::Zero();
    for (Slip &slip_component : *slip_sys) {
        //slip_component.update_strain(stress_tensor + stress_incr, (vel_grad_elas + Matrix3d::Identity()) * deform_grad_elas);
        slip_component.cal_strain(*this);
        vel_grad_plas += slip_component.dL_tensor();
    }
    return vel_grad_plas;
}

void Grain::update_status(Matrix3d L_dt_tensor, Matrix3d vel_grad_flag, Matrix3d stress_incr, Matrix3d dstress_flag){
    Matrix3d vel_grad_elas = Matrix3d::Zero(), vel_grad_plas = Matrix3d::Zero(), stress_iter_save = Matrix3d::Zero(), spin_elas = Matrix3d::Zero();

    // begin iteration
    int iter_num = 0;
    do{
        stress_iter_save = stress_incr;
        //Solve stress_incr, L_dt_tensor
        cal_matrix_dpwp_by_sigma(stress_incr);
        solve_L_dsigma(L_dt_tensor, vel_grad_flag, stress_incr, dstress_flag);
        ++ iter_num;
        if(iter_num > 1000) {cout << (stress_incr - stress_iter_save).norm() << endl << "Not converged... but still going" << endl; break;}
    } while ((stress_incr - stress_iter_save).norm() > 1e-1);
    vel_grad_plas = get_vel_grad_plas(stress_tensor + stress_incr);
    vel_grad_elas = L_dt_tensor - vel_grad_plas;
    // end iteration

    stress_tensor = stress_tensor + stress_incr;
    strain_tensor = strain_tensor + 0.5 * (vel_grad_elas + vel_grad_plas + vel_grad_elas.transpose() + vel_grad_plas.transpose());
    deform_grad_elas = (vel_grad_elas + Matrix3d::Identity()) * deform_grad_elas;
    deform_grad_plas = deform_grad_elas.inverse() * deform_grad;
    for (Slip &slip_component : *slip_sys) slip_component.update_status(*this);
    spin_elas = 0.5 * (vel_grad_elas - vel_grad_elas.transpose());
    orientation = (spin_elas + Matrix3d::Identity()) * orientation * (spin_elas + Matrix3d::Identity()).transpose();
}

Matrix<double, 6, 15> Grain::left_matrix(){
    Vector6d stress_6d = tensor_trans_order(stress_tensor);
    Matrix6d C_ij_pri = get_C_ij_pri(stress_6d);
    Matrix<double, 6, 3> Sigma_ik = get_Sigma_ik(stress_6d);

    Matrix6d minus_delta_il = -Matrix6d::Identity(); 
    CP_SigQ = elastic_modulus * dp_by_sigma + Sigma_ik * wp_by_sigma;

    Matrix<double, 6, 15> left_M;
    left_M << C_ij_pri, Sigma_ik, minus_delta_il-CP_SigQ;
    //cout << left_M << endl;
    return left_M;
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

void Grain::cal_matrix_dpwp_by_sigma(Matrix3d stress_incr){
    dp_by_sigma = Matrix6d::Zero();
    wp_by_sigma = Matrix<double, 3, 6>::Zero();
    Matrix3d sigma_3d = stress_tensor + stress_incr, vel_grad_plas = Matrix3d::Zero();
    Vector6d sigma_6d = tensor_trans_order(sigma_3d), dp_6d, wp_6d; 
    for (int j = 0; j != 6; ++j){
        vel_grad_plas = get_vel_grad_plas(tensor_trans_order(get_vec_only_ith(sigma_6d, j)));
        dp_6d = tensor_trans_order((Matrix3d)(0.5 * vel_grad_plas + 0.5 * vel_grad_plas.transpose()));
        wp_6d = tensor_trans_order((Matrix3d)(0.5 * vel_grad_plas - 0.5 * vel_grad_plas.transpose()));
        for (int i = 0; i != 6; ++i){
            dp_by_sigma(i,j) = dp_6d(i) / (sigma_6d(j)+1e-2);
            if(i >= 3) wp_by_sigma(i-3,j) = wp_6d(i) / (sigma_6d(j)+1e-2);
        }
    }
}

Matrix<double, 15, 15> Grain::mid_matrix(){
    Matrix<double, 15, 15> mid_M;
    mid_M << Matrix3d::Identity(), Matrix3d::Zero(), Matrix3d::Zero(), Matrix<double,3,6>::Zero(),
             Matrix3d::Zero(), 0.5*Matrix3d::Identity(), 0.5*Matrix3d::Identity(), Matrix<double,3,6>::Zero(),
             Matrix3d::Zero(), 0.5*Matrix3d::Identity(), -0.5*Matrix3d::Identity(), Matrix<double,3,6>::Zero(),
             Matrix<double,6,3>::Zero(), Matrix<double,6,6>::Zero(), Matrix6d::Identity(); 
    return mid_M;
}

void Grain::solve_L_dsigma(Matrix3d &vel_grad, Matrix3d &vel_grad_flag, Matrix3d &stress_incr, Matrix3d &dstress_flag){
    Matrix<double, 15, 1> params, flags;
    params << tensor_trans_order_9(vel_grad), tensor_trans_order(stress_incr);
    flags << tensor_trans_order_9(vel_grad_flag), tensor_trans_order(dstress_flag);

    Matrix<double, 6, 15> left_x_mid = left_matrix() * mid_matrix();
    vector<int> unknown_idx(6,0), known_idx(9,0);
    flag_to_idx(flags, known_idx, unknown_idx);

    Matrix6d solveA;
    Matrix<double,6,9> solveB;
    Vector6d unknown_params, right_vec = CP_SigQ * tensor_trans_order(stress_tensor);
    Matrix<double,9,1> known_params;
    solveA << left_x_mid(all,unknown_idx);
    solveB << left_x_mid(all,known_idx);
    
    known_params = params(known_idx);
    unknown_params = solveA.inverse() * (-1*solveB*known_params + right_vec);

    int unk_par_idx = 0;
    for (auto &unk_idx : unknown_idx) {params(unk_idx) = unknown_params(unk_par_idx); ++unk_par_idx;}
    params_convert_to_matrix(params, vel_grad, stress_incr);
}

void Grain::print_stress_strain(ofstream &os){
    os << strain_tensor(0,0) << ',' << strain_tensor(1,1) << ','  << strain_tensor(2,2) << ',' << strain_tensor(0,1) << ',' 
       << strain_tensor(0,2) << ',' << strain_tensor(1,2) << ','  << stress_tensor(0,0) << ',' << stress_tensor(1,1) << ',' 
       << stress_tensor(2,2) << ',' << stress_tensor(0,1) << ','  << stress_tensor(0,2) << ',' << stress_tensor(1,2) << ','  << endl;
}

void Grain::print_stress_strain_screen(ofstream &os){
    os << strain_tensor(0,0) << ',' << strain_tensor(1,1) << ','  << strain_tensor(2,2) << ',' << strain_tensor(0,1) << ',' 
       << strain_tensor(0,2) << ',' << strain_tensor(1,2) << ','  << stress_tensor(0,0) << ',' << stress_tensor(1,1) << ',' 
       << stress_tensor(2,2) << ',' << stress_tensor(0,1) << ','  << stress_tensor(0,2) << ',' << stress_tensor(1,2) << ','  << endl;
    cout << strain_tensor(0,0) << '\t' << strain_tensor(1,1) << '\t'  << strain_tensor(2,2) << '\t' << strain_tensor(0,1) << '\t' 
       << strain_tensor(0,2) << '\t' << strain_tensor(1,2) << '\t'  << stress_tensor(0,0) << '\t' << stress_tensor(1,1) << '\t' 
       << stress_tensor(2,2) << '\t' << stress_tensor(0,1) << '\t'  << stress_tensor(0,2) << '\t' << stress_tensor(1,2) << '\t'  << endl;
}

void Grain::print_dislocation(ofstream &os){
    os << strain_tensor(0,0) << ',' << strain_tensor(1,1) << ','  << strain_tensor(2,2) << ',';
    for (Slip &slip_component : *slip_sys) os << slip_component.SSD_density << ',';
    os << endl;
}

void Grain::print_crss(ofstream &os){
    os << strain_tensor(0,0) << ',' << strain_tensor(1,1) << ','  << strain_tensor(2,2) << ',' << (*slip_sys)[0].acc_strain << ",";
    for (Slip &slip_component : *slip_sys) os << slip_component.crss << ',';
    os << endl;
}