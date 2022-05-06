#include "singleX.h"

Slip::Slip() = default;

Slip::Slip(Vector3d temp_bv, Vector3d temp_pn): burgers_vec(temp_bv), plane_norm(temp_pn/temp_pn.norm()), plane_norm_disp(temp_pn),
        schmidt(temp_bv * (temp_pn/temp_pn.norm()).transpose()) {
            crss = acc_strain = 0;
            strain_rate_slip = 0.01;
};

Slip::Slip(istream &is) {
    is >> burgers_vec(0) >> burgers_vec(1) >> burgers_vec(2) >> plane_norm_disp(0) >> plane_norm_disp(1) >> plane_norm_disp(2);
    plane_norm = plane_norm_disp / plane_norm_disp.norm();
    schmidt = burgers_vec * plane_norm.transpose();
    crss = acc_strain = 0;
    strain_rate_slip = 0.01;
};

Slip::Slip(std::string &is) {
    stringstream stream(is);
    double temp[16]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    int temp_idx = 0;
    while (!stream.eof()) stream >> temp[temp_idx++];
    for (temp_idx=0; temp_idx<16; ++temp_idx){
        if (temp_idx < 3) plane_norm_disp(temp_idx) = temp[temp_idx];
        else {
            if(temp_idx < 6) burgers_vec(temp_idx-3) = temp[temp_idx];
            else harden_params.push_back(temp[temp_idx]);
            }
    }
    plane_norm = plane_norm_disp / plane_norm_disp.norm();
    schmidt = burgers_vec/burgers_vec.norm() * plane_norm.transpose();
    SSD_density = crss = harden_params[0];
    acc_strain = 0;
    strain_rate_slip = 0.01;
    shear_modulus = 0;
};

Matrix3d Slip::dL_tensor() {return schmidt * strain_rate_slip * dtime;}

Matrix3d Slip::dstrain_tensor() {return 0.5 * (schmidt + schmidt.transpose()) * strain_rate_slip * dtime;}

Matrix3d Slip::drotate_tensor() {return 0.5 * (schmidt - schmidt.transpose()) * strain_rate_slip * dtime;}

double Slip::cal_rss(Matrix3d stress_tensor, Matrix3d deform_grad_elas){
    return (deform_grad_elas * (burgers_vec / burgers_vec.norm())).transpose() * stress_tensor * (plane_norm.transpose() * deform_grad_elas.inverse()).transpose();
}

void Slip::cal_shear_modulus(Matrix6d elastic_modulus){
    Matrix3d stress_tensor = update_stress(schmidt, elastic_modulus);
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
        cal_strain_pow(grain.stress_tensor, grain.deform_grad_elas);
        break;
    case 2:
        cal_strain_disvel(grain.stress_tensor, grain.deform_grad_elas, *grain.slip_sys, grain.lattice_cons);
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

void Slip::cal_strain_disvel(Matrix3d stress_tensor, Matrix3d deform_grad_elas, vector<Slip> &slip_sys, double lattice_cons){
    double rss_slip = 0.0;
    double disl_velocity, peierls_distance;
    double tau_ath, tau_th, disl_density_perp, disl_density_for, burgers, cosine_n_m, velocity_move, time_nucleation, freq_kink, temp_value;
    double c_ath = harden_params[1], c_cross = harden_params[2], E_activation = harden_params[3], 
           c_act = harden_params[4], c_l = harden_params[5], velocity_ref = harden_params[6], drag_coeff = harden_params[7];
    burgers = (deform_grad_elas * burgers_vec).norm() * 1e-10 * lattice_cons;

    //calculate disl_density_perp, disl_density_for
    disl_density_for = disl_density_perp = 0;
    for(Slip &isys : slip_sys){
        cosine_n_m =  plane_norm.transpose() * (isys.burgers_vec / isys.burgers_vec.norm());
        disl_density_for += isys.SSD_density * abs(cosine_n_m);
        disl_density_perp += isys.SSD_density * sqrt(1-cosine_n_m*cosine_n_m);
    }
    
    //cal rss and tau_ath, tau_th
    rss_slip = cal_rss(stress_tensor, deform_grad_elas); 
    tau_ath = c_ath * shear_modulus * burgers * sqrt(disl_density_perp); 
    tau_th = (k_boltzmann*temperature + c_cross * E_activation * burgers* sqrt(disl_density_for))/(c_act * pow(burgers,3) * 20) / 1e6;

    //cal shear rate
    if(abs(rss_slip)-tau_ath > 0){
        temp_value = drag_coeff * k_boltzmann * temperature / (2*burgers * burgers * burgers *(abs(rss_slip)-tau_ath)) * 1e-6;
        velocity_move = velocity_ref * (sqrt(1+temp_value*temp_value)-temp_value);
        freq_kink = 2*debye_freq*burgers/(400*burgers) * c_l/sqrt(disl_density_for);
        time_nucleation = 1 / (freq_kink * exp(-1* E_activation/(k_boltzmann * temperature))*sinh((abs(rss_slip)-tau_ath)/tau_th));
        peierls_distance = plane_norm_disp.cross(burgers_vec).norm() * 1e-10;
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
        update_voce(*grain.slip_sys);
        break;
    case 1:
        update_voce(*grain.slip_sys);
        break;
    case 2:
        update_disvel();
        break;
    default:
        update_voce(*grain.slip_sys);
        break;
    }
}

void Slip::update_voce(vector<Slip> &slip_sys){
    double tau_0 = harden_params[0], tau_1 = harden_params[1], h_0 = harden_params[2], h_1 = harden_params[3];
    double dtau_by_dGamma = h_1 + (abs(h_0/tau_1)*tau_1 - h_1) * exp(-acc_strain*abs(h_0/tau_1)) + abs(h_0/tau_1)*h_1*acc_strain*exp(-acc_strain*abs(h_0/tau_1));
    for(Slip &isys : slip_sys){
        crss += abs(isys.strain_rate_slip) * dtime * 1 * dtau_by_dGamma;
        acc_strain += abs(isys.strain_rate_slip) * dtime;
    }
    //acc_strain += strain_rate_slip * dtime;
}

void Slip::update_disvel(){
    double c_multi = harden_params[8], c_annih = harden_params[9];
    SSD_density += (c_multi * sqrt(SSD_density) - c_annih * SSD_density) * abs(strain_rate_slip) * dtime;
}
