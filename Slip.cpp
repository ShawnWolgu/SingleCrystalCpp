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
    shear_modulus = 0;
};

Matrix3d Slip::dL_tensor() {return schmidt * strain_rate_slip * dtime;}

Matrix3d Slip::dstrain_tensor() {return 0.5 * (schmidt + schmidt.transpose()) * strain_rate_slip * dtime;}

Matrix3d Slip::drotate_tensor() {return 0.5 * (schmidt - schmidt.transpose()) * strain_rate_slip * dtime;}

double Slip::cal_rss(Matrix3d stress_tensor, Matrix3d deform_grad_elas){
    return (deform_grad_elas * (burgers_vec / burgers_vec.norm())).transpose() * stress_tensor * (plane_norm.transpose() * deform_grad_elas.inverse()).transpose();
}

void Slip::cal_shear_modulus(Matrix6d elastic_modulus){
    Matrix3d stress_tensor = calc_stress(schmidt, elastic_modulus);
    shear_modulus = cal_rss(stress_tensor, Matrix3d::Identity());
}

void Slip::cal_strain(Grain &grain){
    /* Select different model for shear strain rate, controlled by flag_harden.
     * 0 : Voce Hardening; 1 : Dislocation Density Hardening; 2 : Dislocation Velocity Model;
     * Power model will be used in case 0 and 1, while Velocity model used in case 2;
     */
    switch (flag_harden)
    {
    case 0:
        cal_strain_pow(grain.stress_tensor, grain.deform_grad_elas);
        break;
    case 1:
        cal_strain_ddhard(grain.stress_tensor, grain.deform_grad_elas, grain.strain_rate);
        break;
    case 2:
        cal_strain_disvel(grain.stress_tensor, grain.deform_grad_elas);
        break;
    default:
        cal_strain_pow(grain.stress_tensor, grain.deform_grad_elas);
        break;
    }
}

void Slip::cal_strain_pow(Matrix3d stress_tensor, Matrix3d deform_grad_elas){
    double rss_slip = cal_rss(stress_tensor, deform_grad_elas);       
    if(abs(rss_slip) > 0.5 * crss){
        strain_rate_slip = ref_strain_rate * pow(abs(rss_slip / crss), 1/rate_sen) * rss_slip / abs(rss_slip); 
    }
}   

void Slip::cal_strain_ddhard(Matrix3d stress_tensor, Matrix3d deform_grad_elas, double strain_rate){
    double rss_slip = cal_rss(stress_tensor, deform_grad_elas);  
    double burgers = update_params[1], disl_rem_grad = update_params[2];  
    double c_multi = harden_params[3], ref_rate = harden_params[4], H_activation = harden_params[5],
           drag_stress = harden_params[6];
    const double chi = 0.9, k_sub = 0.086, q = 4;
    if(abs(rss_slip) > 0.5 * crss){
        strain_rate_slip = ref_strain_rate * pow(abs(rss_slip / crss), 1/rate_sen) * rss_slip / abs(rss_slip); 
    }
    disl_rem_grad = c_multi * chi * burgers / H_activation * (1 - (k_boltzmann * temperature)/(drag_stress*1e6 * pow(burgers,3)) * log(strain_rate/ref_rate)) * SSD_density;
    update_params[2] = disl_rem_grad;
}   

void Slip::cal_strain_disvel(Matrix3d stress_tensor, Matrix3d deform_grad_elas){
    double rss_slip = 0.0;
    double disl_velocity, velocity_move, time_nucleation, temp_value;
    double c_ath = harden_params[1], c_cross = harden_params[2], E_activation = harden_params[3], 
           c_act = harden_params[4], c_l = harden_params[5], velocity_ref = harden_params[6], drag_coeff = harden_params[7];
    double burgers = update_params[0], tau_ath = update_params[1], tau_th = update_params[2], peierls_distance = update_params[3], 
           disl_density_perp = update_params[4], disl_density_for = update_params[5], freq_kink = update_params[6];
    
    rss_slip = cal_rss(stress_tensor, deform_grad_elas); 

    //cal shear rate
    if(abs(rss_slip)-tau_ath > 0){
        temp_value = drag_coeff * k_boltzmann * temperature / (2 * pow(burgers,3) * (abs(rss_slip)-tau_ath)) * 1e-6;
        velocity_move = velocity_ref * (sqrt(1+temp_value*temp_value)-temp_value);
        time_nucleation = 1 / (freq_kink * sinh((abs(rss_slip)-tau_ath)/tau_th));
        disl_velocity = peierls_distance / (time_nucleation + peierls_distance / velocity_move);
        strain_rate_slip = SSD_density * burgers * disl_velocity * rss_slip / abs(rss_slip);
    }
    else strain_rate_slip = 0.0;
}

void Slip::update_status(Grain &grain){
    /* Select different model for shear strain rate, controlled by flag_harden.
     * 0 : Voce Hardening; 1 : Dislocation Density Hardening; 2 : Dislocation Velocity Model;
     */
    switch (flag_harden)
    {
    case 0:
        update_voce(grain.slip_sys);
        break;
    case 1:
        update_ddhard(grain.deform_grad_elas, grain.slip_sys);
        break;
    case 2:
        update_disvel(grain.deform_grad_elas, grain.slip_sys);
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
    double tau_0 = harden_params[0], tau_1 = harden_params[1], h_0 = harden_params[2], h_1 = harden_params[3];
    double dtau_by_dGamma = h_1 + (abs(h_0/tau_1)*tau_1 - h_1) * exp(-acc_strain*abs(h_0/tau_1)) + abs(h_0/tau_1)*h_1*acc_strain*exp(-acc_strain*abs(h_0/tau_1));
    for(Slip &isys : slip_sys){
        crss += abs(isys.strain_rate_slip) * dtime * 1 * dtau_by_dGamma;
        acc_strain += abs(isys.strain_rate_slip) * dtime;
    }
}

void Slip::update_ddhard(Matrix3d deform_grad_elas, vector<Slip> &slip_sys){
    /*
     * Update SSD_density and dislocation density hardening model parameters.
     */
    double tau_forest, tau_subs;
    double subs_density = update_params[0], burgers = update_params[1], disl_rem_grad = update_params[2]; 
    double tau_0 = harden_params[0], c_multi = harden_params[3], ref_rate = harden_params[4], H_activation = harden_params[5],
           drag_stress = harden_params[6];
    const double chi = 0.9, k_sub = 0.086, q = 4;

    burgers = (deform_grad_elas * burgers_vec).norm() * 1e-10;
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

void Slip::update_disvel(Matrix3d deform_grad_elas, vector<Slip> &slip_sys){
    /*
     * Update SSD_density and dislocation velocity model parameters.
     */
    double burgers, tau_ath, tau_th, peierls_distance, disl_density_perp, disl_density_for, freq_kink;
    double cosine_n_m;
    double c_ath = harden_params[1], c_cross = harden_params[2], E_activation = harden_params[3], 
           c_act = harden_params[4], c_l = harden_params[5], c_multi = harden_params[8], c_annih = harden_params[9];

    disl_density_for = disl_density_perp = 0;
    for(Slip &isys : slip_sys){
        cosine_n_m =  plane_norm.transpose() * (isys.burgers_vec / isys.burgers_vec.norm());
        disl_density_for += isys.SSD_density * abs(cosine_n_m);
        disl_density_perp += isys.SSD_density * sqrt(1-cosine_n_m*cosine_n_m);
    }
    
    burgers = (deform_grad_elas * burgers_vec).norm() * 1e-10;
    tau_ath = c_ath * shear_modulus * burgers * sqrt(disl_density_perp); 
    tau_th = (k_boltzmann*temperature + c_cross * E_activation * burgers* sqrt(disl_density_for))/(c_act * pow(burgers,3) * 20) / 1e6;
    freq_kink = 2*debye_freq/(400*burgers) * c_l/sqrt(disl_density_for) * exp(-1* E_activation/(k_boltzmann * temperature));
    peierls_distance = plane_norm_disp.cross(burgers_vec).norm() * 1e-10;

    SSD_density += (c_multi * sqrt(SSD_density) - c_annih * SSD_density) * abs(strain_rate_slip) * dtime;
    update_params[0] = burgers, update_params[1] = tau_ath, update_params[2] = tau_th, update_params[3] = peierls_distance, 
    update_params[4] = disl_density_perp, update_params[5] = disl_density_for, update_params[6] = freq_kink;
}
