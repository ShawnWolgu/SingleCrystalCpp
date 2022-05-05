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
            else dd_hardening_params.push_back(temp[temp_idx]);
            }
    }
    plane_norm = plane_norm_disp / plane_norm_disp.norm();
    schmidt = burgers_vec/burgers_vec.norm() * plane_norm.transpose();
    SSD_density = dd_hardening_params[0];
    crss = acc_strain = 0;
    strain_rate_slip = 0.01;
    shear_modulus = 0;
};

Matrix3d Slip::dL_tensor() {return schmidt * strain_rate_slip * dtime;}

Matrix3d Slip::dstrain_tensor() {return 0.5 * (schmidt + schmidt.transpose()) * strain_rate_slip * dtime;}

Matrix3d Slip::drotate_tensor() {return 0.5 * (schmidt - schmidt.transpose()) * strain_rate_slip * dtime;}

void Slip::update_acc_strain() {acc_strain += strain_rate_slip * dtime;}

double Slip::cal_rss(Matrix3d stress_tensor, Matrix3d deform_grad_elas){
    return (deform_grad_elas * (burgers_vec / burgers_vec.norm())).transpose() * stress_tensor * (plane_norm.transpose() * deform_grad_elas.inverse()).transpose();
}

void Slip::update_strain(Matrix3d stress_tensor, Matrix3d deform_grad_elas){
    double rss_slip = 0.0;
    double linear_strain_rate = 0.0;

    rss_slip = cal_rss(stress_tensor, deform_grad_elas);       

    if(abs(rss_slip) > 0.5 * crss){
        strain_rate_slip = ref_strain_rate * pow(abs(rss_slip / crss), 1/m) * rss_slip / abs(rss_slip); 
        linear_strain_rate = ref_strain_rate * rss_slip / crss * 2e10;
        if (abs(linear_strain_rate) < abs(strain_rate_slip)) {
            strain_rate_slip = linear_strain_rate;
        }
    }
}

void Slip::cal_shear_modulus(Matrix6d elastic_modulus){
    Matrix3d stress_tensor = update_stress(schmidt, elastic_modulus);
    shear_modulus = cal_rss(stress_tensor, Matrix3d::Identity());
}

void Slip::cal_strain_ddh(Matrix3d stress_tensor, Matrix3d deform_grad_elas, vector<Slip> &slip_sys, double lattice_cons){
    double rss_slip = 0.0;
    double disl_velocity, peierls_distance;
    double tau_ath, tau_th, disl_density_perp, disl_density_for, burgers, cosine_n_m, velocity_move, time_nucleation, freq_kink, temp_value;
    double c_ath = dd_hardening_params[1], c_cross = dd_hardening_params[2], E_activation = dd_hardening_params[3], 
           c_act = dd_hardening_params[4], c_l = dd_hardening_params[5], velocity_ref = dd_hardening_params[6], drag_coeff = dd_hardening_params[7];
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

void Slip::update_dislocation(){
    double c_multi = dd_hardening_params[8], c_annih = dd_hardening_params[9];
    SSD_density += (c_multi * sqrt(SSD_density) - c_annih * SSD_density) * abs(strain_rate_slip) * dtime;
}
