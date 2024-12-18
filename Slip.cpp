#include "singleX.h"
#include <vector>

Slip::Slip() = default;

Slip::Slip(int number, Vector6d &slip_info, vector<double> &hardens, vector<double> &latents, Matrix3d lattice_vec, double f_active) {
    Vector3d plane_norm_disp;
    num = number, type = slip;
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
    cout << "Plane normal of slip system " << num << " is " << plane_norm.transpose() << endl;
    schmidt = burgers_vec/burgers_vec.norm() * plane_norm.transpose();
    switch (flag_harden)
    {
        case 0:
            crss = harden_params[0];
            break;
        case 1:
            crss = harden_params[5];
            update_params = {0,0,0,0,0};
            SSD_density = harden_params[0];
            if (flag_active == true) SSD_density = SSD_density*f_active; 
            rho_init = rho_mov = rho_H = SSD_density;
            break;
        default:
            crss = harden_params[0];
            break;
    }
    acc_strain = 0;
    shear_rate = 0.00;
    ddgamma_dtau = 0.00;
    shear_modulus = 0;
    disl_vel = 0.0;
    rss = 0.0;
};

void Slip::cal_strain(Grain &grain, Matrix3d stress_tensor) {
    /* 
     * Select different model for shear strain rate, controlled by flag_harden.
     * 0 : Voce Hardening; 1 : Dislocation Velocity Model;
     * Power model will be used in case 0, while Velocity model used in case 1;
     */
    if (flag_active == false){
        shear_rate = 0.0;
        return;
    }
    switch (flag_harden)
    {
        case 0:
            cal_strain_pow(stress_tensor);
            break;
        case 1:
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
        shear_rate = ref_strain_rate * pow(abs(rss_slip / crss), 1/rate_sen)* sign(rss_slip); 
    }
    rss = rss_slip;
}   

void Slip::cal_strain_disvel(Matrix3d stress_tensor){
    double burgers = update_params[0];
    double rss_slip = cal_rss(stress_tensor);
    disl_vel = disl_velocity(rss_slip);
    shear_rate = abs(rho_mov * burgers * disl_vel) * sign(rss_slip);
    rss = rss_slip;
}

void Slip::update_status(Grain &grain){
    /* Select different model for shear strain rate, controlled by flag_harden. */
    /* 0 : Voce Hardening; 1 : Dislocation Velocity Model; */
    /* Update Schmidt here. */
    /* Vector3d update_bv =  burgers_vec; */
    /* Vector3d update_nv =  plane_norm; */
    Vector3d update_bv = grain.orientation * grain.deform_grad_elas * grain.orient_ref.transpose() * burgers_vec;
    Vector3d update_nv = grain.orientation * grain.deform_grad_elas.inverse().transpose() * grain.orient_ref.transpose() * plane_norm;
    schmidt = update_bv / update_bv.norm() * update_nv.transpose();
    switch (flag_harden)
    {
        case 0:
            update_voce(grain.mode_sys, grain.lat_hard_mat);
            break;
        case 1:
            update_disvel(grain.mode_sys, grain.lat_hard_mat, update_bv.norm());
            break;
        default:
            update_voce(grain.mode_sys, grain.lat_hard_mat);
            break;
    }
}
    
/* Update crss */
void Slip::update_voce(vector<PMode*> mode_sys, MatrixXd lat_hard_mat){
    double Gamma = 0;
    for(auto &isys : mode_sys){
        Gamma += abs(isys->acc_strain);
    }
    double tau_0 = harden_params[0], tau_1 = harden_params[1], h_0 = harden_params[2], h_1 = harden_params[3];
    double dtau_by_dGamma = h_1 + (abs(h_0/tau_1)*tau_1 - h_1) * exp(-Gamma*abs(h_0/tau_1)) + abs(h_0/tau_1)*h_1*Gamma*exp(-Gamma*abs(h_0/tau_1));
    for(auto &isys : mode_sys)
    crss += abs(isys->shear_rate) * dtime * lat_hard_mat(num,isys->num) * dtau_by_dGamma;
}

void Slip::update_disvel(vector<PMode*> mode_sys, MatrixXd lat_hard_mat, double bv_norm){
    /*
     * [velocity parameters] 
     *  1. MFP control coeffient, 2. reference frequency, 3. activation energy, 4. slip resistance, 5. energy exponent
     *  6. saturated speed, 7. drag coefficient
     * [hardening parameters] 
     *  8. forest hardening coefficient
     * [DD evolution parameters] 
     *  0. SSD_density, 9. nucleation coefficient, 10. nucleation threshold stress, 11. multiplication coefficient
     *  12. drag stress D, 13. reference strain rate, 14. c/g 
     *
     * update parameters:
     * 0: burgers, 1: mean_free_path, 2: disl_density_resist, 3: forest_stress
     */
    double burgers = bv_norm * 1e-10, joint_factor = 0.0;
    double c_mfp = harden_params[1], resistance_slip = harden_params[4], c_forest = harden_params[8], HP_stress = 0;
    double disl_density_for, disl_density_resist, joint_density, forest_stress, mean_free_path;
    disl_density_for = disl_density_resist = joint_density = 0;
    for(auto &isys : mode_sys){
        disl_density_for += isys->SSD_density;
        disl_density_resist += isys->rho_H * lat_hard_mat(num,isys->num);
        if(isys->num != num) joint_density += lat_hard_mat(num,isys->num) * sqrt(isys->rho_H-isys->rho_init) * sqrt(rho_H-rho_init);
    }
    double crss_factor = joint_factor*joint_density+disl_density_resist;
    forest_stress = c_forest * shear_modulus * burgers * sqrt(crss_factor);
    mean_free_path = c_mfp / sqrt(disl_density_for);
    crss = forest_stress;
    acc_strain += abs(shear_rate) * dtime;
    update_params[0] = burgers, update_params[1] = mean_free_path, \
    update_params[2] = disl_density_resist, update_params[3] = forest_stress;
}

void Slip::cal_ddgamma_dtau(Matrix3d stress_tensor){
    /* Select different model for shear strain rate gradient calculation, controlled by flag_harden. */
    /* 0 : Voce Hardening; 1 : Dislocation Velocity Model; */
    /* Update Schmidt here. */
    /* Vector3d update_bv =  burgers_vec; */
    /* Vector3d update_nv =  plane_norm; */
    if (flag_active){
        ddgamma_dtau = 0.0; shear_rate = 0.0;
        return;
    }
    switch (flag_harden)
    {
        case 0:
            cal_ddgamma_dtau_pow(stress_tensor);
            break;
        case 1:
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
        shear_rate = ref_strain_rate * pow(abs(rss_slip / crss), 1/rate_sen) * sign(rss_slip); 
    }
}   

void Slip::cal_ddgamma_dtau_disvel(Matrix3d stress_tensor){
    double burgers = update_params[0];
    double rss_slip = cal_rss(stress_tensor);
    vector<double> dvel_and_vel = disl_velocity_grad(rss_slip);
    ddgamma_dtau = rho_mov * burgers * sign(rss_slip) * dvel_and_vel[0];
    shear_rate = rho_mov * burgers * dvel_and_vel[1] * sign(rss_slip);
}

/* Some functions not used in the current version */

/* [Not Used]Update LH parameters */
void Slip::update_lhparams(Matrix3d dstrain){
    double lh_coeff = 0.; // Should be add as field var.
    if (flag_harden == 2){ 
        double ref_srate = 1e-3, exp_lh = -0.1;
        lh_coeff = pow(calc_equivalent_value(dstrain)/dtime/ref_srate, exp_lh);
        if (lh_coeff > 2) lh_coeff = 2;
    }
    else{}
}

/* [Not Used]Update dislocation densities under cross slip processes  */
void Slip::update_cross_slip(vector<PMode> &slip_sys, Matrix3d stress_tensor){
    vector<double> cross_params = {0,0,0,0}; // Should be add as field var.
    double cross_in = 0, cross_out = 0; // Should be add as field var.
    if (flag_harden == 2){ 
        double burgers = update_params[0], para = 0.98, rss_slip = cal_rss(stress_tensor), back_stress = update_params[3];
        double nu_cross = cross_params[0], phi = cross_params[1], cross_stress = cross_params[2], volume_cross = cross_params[3]*pow(burgers,3);
        for(auto &isys : slip_sys){
            if ((isys.num != num) && (abs(cal_cosine(isys.burgers_vec,burgers_vec))>para)){
                double rss_isys = isys.cal_rss(stress_tensor);
                double exp_term_in = (rss_slip - (crss + back_stress))/(k_boltzmann * temperature) * volume_cross;
                double exp_term_out = (rss_isys - (isys.crss + isys.update_params[3]))/(k_boltzmann * temperature) * volume_cross;
                cross_in += nu_cross * phi * isys.SSD_density * exp(exp_term_in);
                cross_out += nu_cross * phi * SSD_density * exp(exp_term_out);
            }
        }
    }
}

/* [Not Used]Update movable dislocation densities under Orowan processes  */
void Slip::update_rho_mov(vector<PMode> &mode_sys){
    if (flag_harden == 2){ 
        double burgers = update_params[0], para = 0.98, coeff = 0.;
        rho_mov = SSD_density;
        for(auto &isys : mode_sys){
            if ((isys.num != num) && (abs(cal_cosine(isys.burgers_vec,burgers_vec))>para)){
                rho_mov += isys.SSD_density * coeff;
            }
        }
    }
}

/* [Not Used]Update dislocation densities created under surface nucleation processes  */
void Slip::update_surface_nuc(Matrix3d stress_tensor){
    vector<double> surf_params = {0,0,0,0,0,0}; // Should be add as field var.
    double dSSD_surface = 0; // Should be add as field var.
    if (flag_harden == 2){ 
        double burgers = update_params[0], rss_slip = cal_rss(stress_tensor), back_stress = update_params[3], expo_alpha = harden_params[6];
        double energy_nuc = surf_params[0] * eV_to_J, c_tau = surf_params[1], freq_surfnuc = surf_params[2], distance_plane = surf_params[3] * 1e-10, grain_diameter = surf_params[4], shape_param = surf_params[5];
        double ssd_term = pow((1-SSD_density/rho_sat),3);
        double exp_term = energy_nuc * (1-pow(abs(rss_slip/(back_stress*c_tau)),expo_alpha));
        exp_term = min(exp_term,500*k_boltzmann*temperature);
        dSSD_surface = ssd_term * shape_param / (distance_plane * grain_diameter) * freq_surfnuc * exp(-exp_term/(k_boltzmann*temperature));
    }
}
