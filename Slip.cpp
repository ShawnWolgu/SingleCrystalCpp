#include "singleX.h"

Slip::Slip() = default;

Slip::Slip(int number, Vector6d &slip_info, vector<double> &hardens, vector<double> &latents, vector<double> &surf, Matrix3d lattice_vec) {
    num = number;
    harden_params = hardens; latent_params = latents; surf_params = surf;
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
	crss = harden_params[5];
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
        strain_rate_slip = ref_strain_rate * pow(abs(rss_slip / crss), 1/rate_sen)* sign(rss_slip); 
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
    disl_vel = disl_velocity(rss_slip);
    strain_rate_slip = abs(SSD_density * burgers * disl_vel) * sign(rss_slip);
    //strain_rate_slip = abs(disl_vel) * sign(rss_slip);
}

void Slip::update_status(Grain &grain){
    /* Select different model for shear strain rate, controlled by flag_harden.
     * 0 : Voce Hardening; 1 : Dislocation Density Hardening; 2 : Dislocation Velocity Model;
     */
    //Update Schmidt here.
    //Vector3d update_bv =  burgers_vec;
    //Vector3d update_nv =  plane_norm;
    Vector3d update_bv = grain.orientation * grain.deform_grad_elas * grain.orient_ref.transpose() * burgers_vec;
    Vector3d update_nv = grain.orientation * grain.deform_grad_elas.inverse().transpose() * grain.orient_ref.transpose() * plane_norm;
    schmidt = update_bv / update_bv.norm() * update_nv.transpose();
    switch (flag_harden)
    {
    case 0:
        update_voce(grain.slip_sys, grain.lat_hard_mat);
        break;
    case 1:
        update_ddhard(grain.slip_sys, grain.lat_hard_mat, update_bv.norm());
        break;
    case 2:
        update_disvel(grain.slip_sys, grain.lat_hard_mat, update_bv.norm());
        break;
    default:
        update_voce(grain.slip_sys, grain.lat_hard_mat);
        break;
    }
}

void Slip::update_voce(vector<Slip> &slip_sys, MatrixXd lat_hard_mat){
    /*
     * Update crss and acc_strain.
     */
    double Gamma = 0;
    for(Slip &isys : slip_sys){
        Gamma += abs(isys.acc_strain);
    }
    double tau_0 = harden_params[0], tau_1 = harden_params[1], h_0 = harden_params[2], h_1 = harden_params[3];
    double dtau_by_dGamma = h_1 + (abs(h_0/tau_1)*tau_1 - h_1) * exp(-Gamma*abs(h_0/tau_1)) + abs(h_0/tau_1)*h_1*Gamma*exp(-Gamma*abs(h_0/tau_1));
    for(Slip &isys : slip_sys){
        crss += abs(isys.strain_rate_slip) * dtime * lat_hard_mat(num,isys.num) * dtau_by_dGamma;
	//acc_strain += abs(isys.strain_rate_slip)*dtime;
    }
}

void Slip::update_ddhard(vector<Slip> &slip_sys, MatrixXd lat_hard_mat, double bv_norm){
    /*
     * Update SSD_density and dislocation density hardening model parameters.
     * Warning: This function is not maintained in recent versions, maybe not accurate!
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

void Slip::update_ssd(){
    if (flag_harden == 0 || flag_harden == 1){
	acc_strain += abs(strain_rate_slip) * dtime;
    }
    if (flag_harden == 2){ 
    	double c_multi = harden_params[10], c_annih = harden_params[11];
    	SSD_density += (c_multi * sqrt(SSD_density) - c_annih * SSD_density) * abs(strain_rate_slip) * dtime + (dSSD_surface) * dtime;
    }
}

void Slip::update_cross_slip(vector<Slip> &slip_sys, Matrix3d stress_tensor){
    if (flag_harden == 2){ 
    	double burgers = update_params[0], para = 0.98, rss_slip = cal_rss(stress_tensor), back_stress = update_params[3];
    	double nu_cross = cross_params[0], phi = cross_params[1], cross_stress = cross_params[2], volume_cross = cross_params[3]*pow(burgers,3);
    	for(Slip &isys : slip_sys){
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

void Slip::update_surface_nuc(Matrix3d stress_tensor){
    if (flag_harden == 2){ 
    	double burgers = update_params[0], rss_slip = cal_rss(stress_tensor), back_stress = update_params[3], expo_alpha = harden_params[6];
    	double energy_nuc = surf_params[0] * eV_to_J, c_tau = surf_params[1], freq_surfnuc = surf_params[2], distance_plane = surf_params[3] * 1e-10, grain_diameter = surf_params[4], shape_param = surf_params[5];
	double ssd_sat = pow((harden_params[10]/harden_params[11]),2);
	double ssd_term = pow((1-SSD_density/ssd_sat),3);
	double exp_term = energy_nuc * (1-pow(abs(rss_slip/(back_stress*c_tau)),expo_alpha));
	exp_term = min(exp_term,500*k_boltzmann*temperature);
	dSSD_surface = ssd_term * shape_param / (distance_plane * grain_diameter) * freq_surfnuc * exp(-exp_term/(k_boltzmann*temperature));
    }
}

void Slip::update_disvel(vector<Slip> &slip_sys, MatrixXd lat_hard_mat, double bv_norm){
    /*
     * harden parameters: 0: SSD_density,
     * 1: freq_Debye, 2: c_length, 3: kink_energy_ref, 4: temperature_ref,
     * 5: Peierls_stress, 6: expo_kinkeng, 7: wave_speed, 8: c_drag, 9: c_backstress,
     * 10: c_multi, 11:c_annih, 12:HP_stress.
     * 
     * update parameters:
     * 0: burgers, 1: disl_density_for, 2: disl_density_para, 3: back_stress,
     * 4: barrier_distance
     *
     * cross slip parameters:
     * 0: nu_cross, 1: phi, 2: cross_stress, 3: c_volume_cross
     */
    double Peierls_stress = harden_params[5], c_backstress = harden_params[9], HP_stress = harden_params[12];
    double burgers, disl_density_for, disl_density_resist, back_stress, barrier_distance, cosine_n_m, ref_strain;
    disl_density_for = disl_density_resist = 0;
    for(Slip &isys : slip_sys){
	Vector3d t_vector = isys.plane_norm.cross(isys.burgers_vec);
        cosine_n_m =  plane_norm.transpose() * (t_vector / t_vector.norm());
        disl_density_for += isys.SSD_density;// * abs(cosine_n_m);
        disl_density_resist += isys.SSD_density * lat_hard_mat(num,isys.num);// * sqrt(1-cosine_n_m*cosine_n_m);
    }
    burgers = bv_norm * 1e-10;
    //burgers = burgers_vec.norm() * 1e-10;
    ref_strain = burgers * SSD_density / (2 * sqrt(disl_density_for));
    back_stress = c_backstress * shear_modulus * burgers * sqrt(disl_density_resist);// + HP_stress
    //CRSS iteration:
    double crss_0 = crss, f_g = 0, f_g_grad = 0, crss_norm = crss/back_stress, dg = 0;
    do{
	f_g = crss_0 - crss + abs(strain_rate_slip) * (0.01*back_stress/ref_strain) * pow(crss_norm,3) * (cosh(pow(crss_norm,-2))-1) * dtime;
	f_g_grad = -1 + 0.01*back_stress/ref_strain*abs(strain_rate_slip)*dtime*(3*pow(crss,2)/pow(back_stress,3)*(cosh(pow(crss_norm,-2))-1)-pow(crss_norm,3)*sinh(pow(crss_norm,-2))*2*back_stress/pow(crss,2));
	if (dg*f_g/f_g_grad <= 0) dg = -f_g/f_g_grad;
	else dg = -0.5 *dg;
	crss += dg;
	crss_norm = crss/back_stress;
    } while(abs(f_g)>1e-1 && abs(dg)<1e-3);
    //crss += abs(strain_rate_slip) * (0.01 * back_stress/ref_strain) * pow(crss_norm,3) * (cosh(pow(crss_norm,-2))-1) * dtime;
    barrier_distance = plane_norm_disp.cross(burgers_vec).norm() * 1e-10;
    acc_strain += abs(strain_rate_slip) * dtime;
    //cout << strain_rate_slip << endl;
    update_params[0] = burgers, update_params[1] = disl_density_for, update_params[2] = disl_density_resist;
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
    vector<double> dvel_and_vel = disl_velocity_grad(rss_slip, crss, harden_params, update_params);
    //ddgamma_dtau = dvel_and_vel[0] * sign(rss_slip);
    //strain_rate_slip = dvel_and_vel[1] * sign(rss_slip);
    ddgamma_dtau = SSD_density * burgers * sign(rss_slip) * dvel_and_vel[0];
    strain_rate_slip = SSD_density * burgers * dvel_and_vel[1] * sign(rss_slip);
}
