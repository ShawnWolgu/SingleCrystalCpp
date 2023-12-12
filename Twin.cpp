#include "singleX.h"

Twin::Twin() = default;

Twin::Twin(int number, Vector6d &slip_info, vector<double> &hardens, vector<double> &latents, Matrix3d lattice_vec, double f_active){
    Vector3d plane_norm_disp;
    num = number, type = twin;
    harden_params = hardens; latent_params = latents;
    flag_active = !(f_active < 1e-20);
    for (int temp_idx=0; temp_idx<6; ++temp_idx){
        if (temp_idx < 3) plane_norm_disp(temp_idx) = slip_info(temp_idx);
        else {
            if(temp_idx < 6) burgers_vec(temp_idx-3) = slip_info(temp_idx);
        }
    }
    burgers_vec = (burgers_vec.transpose() * lattice_vec).transpose();
    plane_norm = get_plane_norm(plane_norm_disp, lattice_vec);
    cout << "Plane normal of Twin system " << num << " is " << plane_norm.transpose() << endl;
    schmidt = burgers_vec/burgers_vec.norm() * plane_norm.transpose();
    crss = harden_params[0];
    rate_sen = harden_params[5];
    acc_strain = 0;
    shear_rate = 0.00;
    ddgamma_dtau = 0.00;
    shear_modulus = 0;
    twin_frac = 0.0;
}

void Twin::cal_strain(Grain &grain, Matrix3d stress_tensor){
    double rss_matrix = cal_rss(stress_tensor);
    double rate_ = ref_strain_rate * pow(abs(rss_matrix / crss), 1/rate_sen)* sign(rss_matrix) * equiv_frac; 
    switch (status) {
    inactive:
        shear_rate = (rss_matrix > 0) ? rate_ : 0.0; // rate_TN;
        break;
    growth:
        shear_rate = rate_;
        break;
    saturated:
        shear_rate = (rss_matrix < 0) ? rate_ : 0.0; // rate_MP;
        break;
    default:
        shear_rate = (rss_matrix > 0) ? rate_ : 0.0;
        break;
    }
    rss = rss_matrix;
}

void Twin::cal_ddgamma_dtau(Matrix3d stress_tensor){
    double rss_matrix = cal_rss(stress_tensor);       
    double gradient_ = ref_strain_rate * pow(abs(rss_matrix / crss), 1/rate_sen-1) * sign(rss_matrix) / rate_sen / crss * sign(rss_matrix) * equiv_frac; 
    double rate_ = ref_strain_rate * pow(abs(rss_matrix / crss), 1/rate_sen) * sign(rss_matrix) * equiv_frac; 
    switch (status){
    inactive:
        shear_rate = (rss_matrix > 0) ? rate_ : 0.0; // rate_TN;
        ddgamma_dtau = (rss_matrix > 0) ? gradient_ : 0.0; // rate_TN;
        break;
    growth:
        shear_rate = rate_;
        ddgamma_dtau = gradient_;
        break;
    saturated:
        shear_rate = (rss_matrix < 0) ? rate_ : 0.0; // rate_MP;
        ddgamma_dtau = (rss_matrix < 0) ? gradient_ : 0.0; // rate_MP;
        break;
    default:
        shear_rate = (rss_matrix > 0) ? rate_ : 0.0;
        ddgamma_dtau = (rss_matrix > 0) ? gradient_ : 0.0;
        break;
    }
}

void Twin::update_status(Grain &grain){
    double Gamma = 0;
    double tau_0 = harden_params[0], tau_1 = harden_params[1], h_0 = harden_params[2], h_1 = harden_params[3], \
           lower_bound = harden_params[6], upper_bound = harden_params[7];

    for(auto &isys : grain.mode_sys) Gamma += abs(isys->acc_strain);
    if (grain.twin_frac > upper_bound) status = saturated;

    double dtau_by_dGamma = h_1 + (abs(h_0/tau_1)*tau_1 - h_1) * exp(-Gamma*abs(h_0/tau_1)) + abs(h_0/tau_1)*h_1*Gamma*exp(-Gamma*abs(h_0/tau_1));
    for(auto &isys : grain.mode_sys)
        crss += abs(isys->shear_rate) * dtime * grain.lat_hard_mat(num,isys->num) * dtau_by_dGamma;
    
    equiv_frac = (1-grain.twin_frac) + twin_frac;
}

void Twin::update_ssd(Matrix3d dstrain, Matrix3d orientation){
    double lower_bound = harden_params[6], upper_bound = harden_params[7];
    acc_strain += shear_rate * dtime;
    twin_frac += shear_rate / harden_params[4] * dtime;
    SSD_density = twin_frac;
    if (twin_frac < lower_bound) status = inactive;
    else if (twin_frac > upper_bound) status = saturated;
    else status = growth;
}
