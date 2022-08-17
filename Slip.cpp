#include "singleX.h"

Slip::Slip() = default;

Slip::Slip(Vector6d &slip_info, vector<double> &hardens, Matrix3d lattice_vec) {
    harden_params = hardens;
    for (int temp_idx=0; temp_idx<6; ++temp_idx){
        if (temp_idx < 3) plane_norm_disp(temp_idx) = slip_info(temp_idx);
        else {
            if(temp_idx < 6) burgers_vec(temp_idx-3) = slip_info(temp_idx);
        }
    }
    plane_norm_disp = (plane_norm_disp.transpose() * lattice_vec).transpose();
    burgers_vec = (burgers_vec.transpose() * lattice_vec).transpose();

    plane_norm = plane_norm_disp / plane_norm_disp.norm();
    schmidt = burgers_vec/burgers_vec.norm() * plane_norm.transpose();

    switch (flag_harden)
    {
    case 0:
        crss = harden_params[0];
        break;
    case 1:
        update_params = {0,0,0};
        crss = harden_params[0];
        SSD_density = harden_params[1];
        update_params[0] = harden_params[2]; // substructure dislocation density
        break;
    case 2:
        update_params = {0,0,0,0,0,0,0};
        SSD_density = harden_params[0];
        break;
    default:
        crss = harden_params[0];
        break;
    }

    acc_strain = 0;
    strain_rate_slip = 0.00;
    ddgamma_dtau = 0.00;
    shear_modulus = 0;
    disl_vel = 0.0;
};

Matrix3d Slip::dL_tensor() {return schmidt * strain_rate_slip * dtime;}

Matrix3d Slip::dstrain_tensor() {return 0.5 * (schmidt + schmidt.transpose()) * strain_rate_slip * dtime;}

Matrix3d Slip::drotate_tensor() {return 0.5 * (schmidt - schmidt.transpose()) * strain_rate_slip * dtime;}

Matrix6d Slip::ddp_dsigma() {
    Vector6d symSchmidt_6d = strain_modi_tensor * tensor_trans_order((Matrix3d)(0.5*(schmidt+schmidt.transpose())));
    Matrix6d symSchmidt_66 = symSchmidt_6d * symSchmidt_6d.transpose();
    return symSchmidt_66 * ddgamma_dtau * dtime;
}

Matrix6d Slip::dwp_dsigma() {
    Vector6d symSchmidt_6d = strain_modi_tensor * tensor_trans_order((Matrix3d)(0.5*(schmidt+schmidt.transpose())));
    Vector6d asymSchmidt_6d = strain_modi_tensor * tensor_trans_order((Matrix3d)(0.5*(schmidt-schmidt.transpose())));
    Matrix6d symSchmidt_66 = asymSchmidt_6d * symSchmidt_6d.transpose();
    return symSchmidt_66 * ddgamma_dtau * dtime;
}

double Slip::cal_rss(Matrix3d stress_tensor){
    return (stress_tensor.cwiseProduct(schmidt)).sum();
}

void Slip::cal_shear_modulus(Matrix6d elastic_modulus){
    Matrix3d slip_rotation;
    Vector3d trav_direc = burgers_vec.cross(plane_norm);
    slip_rotation << (burgers_vec/burgers_vec.norm()), plane_norm, trav_direc / trav_direc.norm();
    shear_modulus = rotate_6d_stiff_modu(elastic_modulus, slip_rotation.transpose())(3,3);
}

void Slip::cal_strain(Grain &grain, Matrix3d stress_tensor){
    /* Select different model for shear strain rate, controlled by flag_harden.
     * 0 : Voce Hardening; 1 : Dislocation Density Hardening; 2 : Dislocation Velocity Model;
     * Power model will be used in case 0 and 1, while Velocity model used in case 2;
     */
    switch (flag_harden)
    {
    case 0:
        cal_strain_pow(stress_tensor);
        break;
    case 1:
        cal_strain_ddhard(stress_tensor, grain.strain_rate);
        break;
    case 2:
        cal_strain_disvel(stress_tensor);
        break;
    default:
        cal_strain_pow(stress_tensor);
        break;
    }
}

void Slip::cal_strain_pow(Matrix3d stress_tensor){
    double rss_slip = cal_rss(stress_tensor);       
    if(abs(rss_slip) > 0.5 * crss){
        strain_rate_slip = ref_strain_rate * pow(abs(rss_slip / crss), 1/rate_sen) * sign(rss_slip); 
    }
    //cout << "rss, dgamma : " << rss_slip << "   " << strain_rate_slip*dtime << endl;
}   

void Slip::cal_strain_ddhard(Matrix3d stress_tensor, double strain_rate){
    double rss_slip = cal_rss(stress_tensor);  
    double burgers = update_params[1], disl_rem_grad = update_params[2];  
    double c_multi = harden_params[3], ref_rate = harden_params[4], H_activation = harden_params[5],
           drag_stress = harden_params[6];
    const double chi = 0.9, k_sub = 0.086, q = 4;
    if(abs(rss_slip) > 0.5 * crss){
        strain_rate_slip = ref_strain_rate * pow(abs(rss_slip / crss), 1/rate_sen) * sign(rss_slip); 
    }
    disl_rem_grad = c_multi * chi * burgers / H_activation * (1 - (k_boltzmann * temperature)/(drag_stress*1e6 * pow(burgers,3)) * log(strain_rate/ref_rate)) * SSD_density;
    update_params[2] = disl_rem_grad;
}   

void Slip::cal_strain_disvel(Matrix3d stress_tensor){
    double burgers = update_params[0];
    double rss_slip = cal_rss(stress_tensor);
    disl_vel = disl_velocity(rss_slip, harden_params, update_params);
    strain_rate_slip = SSD_density * burgers * disl_vel * sign(rss_slip);
    //if(abs(rss_slip)>10){
    //cout << "rss, back_stress , dgamma : " << rss_slip << "   " << update_params[3] << "   " << strain_rate_slip*dtime << endl;}
}

void Slip::update_status(Grain &grain){
    /* Select different model for shear strain rate, controlled by flag_harden.
     * 0 : Voce Hardening; 1 : Dislocation Density Hardening; 2 : Dislocation Velocity Model;
     */
    //Update Schmidt here.
    //Vector3d update_bv =  burgers_vec;
    Vector3d update_bv = grain.orientation * grain.deform_grad_elas * grain.orient_ref.transpose() * burgers_vec;
    Vector3d update_nv = grain.orientation * grain.deform_grad_elas.inverse().transpose() * grain.orient_ref.transpose() * plane_norm;
    schmidt = update_bv / update_bv.norm() * update_nv.transpose();
    switch (flag_harden)
    {
    case 0:
        update_voce(grain.slip_sys);
        break;
    case 1:
        update_ddhard(grain.slip_sys,update_bv.norm());
        break;
    case 2:
        update_disvel(grain.slip_sys,update_bv.norm());
        break;
    default:
        update_voce(grain.slip_sys);
        break;
    }
}

void Slip::update_voce(vector<Slip> &slip_sys){
    /*
     * Update crss and acc_strain.
     */
    double Gamma = 0;
    for(Slip &isys : slip_sys){
        Gamma += abs(isys.acc_strain);
    }
    acc_strain += abs(strain_rate_slip) * dtime;
    double tau_0 = harden_params[0], tau_1 = harden_params[1], h_0 = harden_params[2], h_1 = harden_params[3];
    double dtau_by_dGamma = h_1 + (abs(h_0/tau_1)*tau_1 - h_1) * exp(-Gamma*abs(h_0/tau_1)) + abs(h_0/tau_1)*h_1*Gamma*exp(-Gamma*abs(h_0/tau_1));
    for(Slip &isys : slip_sys){
        crss += abs(isys.strain_rate_slip) * dtime * 1 * dtau_by_dGamma;
	//acc_strain += abs(isys.strain_rate_slip)*dtime;
    }
}

void Slip::update_ddhard(vector<Slip> &slip_sys, double bv_norm){
    /*
     * Update SSD_density and dislocation density hardening model parameters.
     */
    double tau_forest, tau_subs;
    double subs_density = update_params[0], burgers = update_params[1], disl_rem_grad = update_params[2]; 
    double tau_0 = harden_params[0], c_multi = harden_params[3], ref_rate = harden_params[4], H_activation = harden_params[5],
           drag_stress = harden_params[6];
    const double chi = 0.9, k_sub = 0.086, q = 4;

    burgers = bv_norm * 1e-10;
    tau_forest = chi * burgers * shear_modulus * sqrt(SSD_density);
    tau_subs = k_sub * shear_modulus * burgers * sqrt(subs_density) * log10(1/burgers/sqrt(subs_density));

    crss = tau_0 + tau_forest + tau_subs;
    SSD_density += (c_multi * sqrt(SSD_density) - disl_rem_grad) * abs(strain_rate_slip) * dtime;
    for(Slip &isys : slip_sys){
        // harden_params[7]: fraction to recovery, update_params[1]: burgers, update_params[2]: disl_rem_grad;
        subs_density += q * isys.harden_params[7] * isys.update_params[1] * sqrt(subs_density) *  abs(isys.strain_rate_slip) * dtime * isys.update_params[2];
        acc_strain += abs(isys.strain_rate_slip) * dtime;
    }
    update_params[0] = subs_density, update_params[1] = burgers; 
}

void Slip::update_disvel(vector<Slip> &slip_sys, double bv_norm){
    /*
     * harden parameters: 0: SSD_density,
     * 1: freq_Debye, 2: c_length, 3: kink_energy_ref, 4: temperature_ref,
     * 5: Peierls_stress, 6: expo_kinkeng, 7: wave_speed, 8: c_drag, 9: c_backstress,
     * 10: c_multi, 11:c_annih, 12:HP_stress.
     * 
     * update parameters:
     * 0: burgers, 1: disl_density_for, 2: disl_density_para, 3: back_stress,
     * 4: barrier_distance
     */
    double Peierls_stress = harden_params[5], c_backstress = harden_params[9], c_multi = harden_params[10], c_annih = harden_params[11], HP_stress = harden_params[12];
    double burgers, disl_density_for, disl_density_para, back_stress, barrier_distance, cosine_n_m;

    SSD_density += (c_multi * sqrt(SSD_density) - c_annih * SSD_density) * abs(strain_rate_slip) * dtime;
    disl_density_for = disl_density_para = 0;
    for(Slip &isys : slip_sys){
	Vector3d t_vector = isys.plane_norm.cross(isys.burgers_vec);
        cosine_n_m =  plane_norm.transpose() * (t_vector / t_vector.norm());
        disl_density_for += isys.SSD_density * abs(cosine_n_m);
        disl_density_para += isys.SSD_density * sqrt(1-cosine_n_m*cosine_n_m);
    }
    
    //disl_density_para = SSD_density;
    burgers = bv_norm * 1e-10;
    //burgers = burgers_vec.norm() * 1e-10;
    back_stress = c_backstress * shear_modulus * burgers * sqrt(disl_density_para) + HP_stress;
    crss = Peierls_stress + back_stress; 
    barrier_distance = plane_norm_disp.cross(burgers_vec).norm() * 1e-10;
    acc_strain += abs(strain_rate_slip) * dtime;
    //cout << strain_rate_slip << endl;
    update_params[0] = burgers, update_params[1] = disl_density_for, update_params[2] = disl_density_para;
    update_params[3] = back_stress, update_params[4] = barrier_distance;
}

void Slip::cal_ddgamma_dtau(Matrix3d stress_tensor){
    /* Select different model for the gradient of shear strain rate by rss, controlled by flag_harden.
     * Note the slip rate will also be update in this function.
     * 0 : Voce Hardening; 1 : Dislocation Density Hardening; 2 : Dislocation Velocity Model;
     * Power model will be used in case 0 and 1, while Velocity model used in case 2;
     */
    switch (flag_harden)
    {
    case 0:
        cal_ddgamma_dtau_pow(stress_tensor);
        break;
    case 1:
        cal_ddgamma_dtau_ddhard(stress_tensor);
        break;
    case 2:
        cal_ddgamma_dtau_disvel(stress_tensor);
        break;
    default:
        cal_ddgamma_dtau_pow(stress_tensor);
        break;
    }
}

void Slip::cal_ddgamma_dtau_pow(Matrix3d stress_tensor){
    double rss_slip = cal_rss(stress_tensor);       
    if(abs(rss_slip) > 0.5 * crss){
        ddgamma_dtau = ref_strain_rate * pow(abs(rss_slip / crss), 1/rate_sen-1) * sign(rss_slip) / rate_sen / crss * sign(rss_slip); 
        strain_rate_slip = ref_strain_rate * pow(abs(rss_slip / crss), 1/rate_sen) * sign(rss_slip); 
    }
}   

void Slip::cal_ddgamma_dtau_ddhard(Matrix3d stress_tensor){
    double rss_slip = cal_rss(stress_tensor);  
    if(abs(rss_slip) > 0.5 * crss){
        ddgamma_dtau = ref_strain_rate * pow(abs(rss_slip / crss), 1/rate_sen-1) * sign(rss_slip) / rate_sen / crss * sign(rss_slip); 
        strain_rate_slip = ref_strain_rate * pow(abs(rss_slip / crss), 1/rate_sen) * sign(rss_slip); 
    }
}   

void Slip::cal_ddgamma_dtau_disvel(Matrix3d stress_tensor){
    double burgers = update_params[0];
    double rss_slip = cal_rss(stress_tensor);
    vector<double> dvel_and_vel = disl_velocity_grad(rss_slip, harden_params, update_params);
    ddgamma_dtau = SSD_density * burgers * sign(rss_slip) * dvel_and_vel[0];
    strain_rate_slip = SSD_density * burgers * dvel_and_vel[1] * sign(rss_slip);
}
