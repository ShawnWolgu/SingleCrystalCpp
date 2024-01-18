#include "singleX.h"
#include <vector>

void Slip::update_ssd_old(Matrix3d strain_rate, Matrix3d orientation){
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
    if (flag_harden == 0) acc_strain += abs(shear_rate) * dtime;
    if (flag_harden == 1){ 
        double c_forest = harden_params[8], c_nuc = harden_params[9], tau_nuc = harden_params[10],\
               c_multi = harden_params[11], c_annih = 0.,\
               D = harden_params[12] * 1e6, ref_srate = harden_params[13], gg = c_forest/harden_params[14],\
               burgers = update_params[0], mfp = update_params[1], forest_stress = update_params[3]; 
        /* double equi_strain_rate = calc_equivalent_value(strain_rate); */
        /* double equi_strain_rate = calc_first_principal(strain_rate); */
        /* c_multi = (custom_var > 1) ? 1/log(exp(1) * custom_var)*c_multi : c_multi; */
        double equi_strain_rate = strain_rate(2,2);
        rho_sat = c_forest * burgers / gg * (1-k_boltzmann * temperature/D/pow(burgers,3) * log(abs(equi_strain_rate)/ref_srate));
        rho_sat = max(pow(1/rho_sat,2), 0.5*SSD_density);
        double term_nuc = c_nuc * max(abs(rss)-tau_nuc,0.) / (shear_modulus * burgers * burgers);
        double term_multi = c_multi / mfp; 
        c_annih = (term_multi + term_nuc) / rho_sat;
        SSD_density += (term_multi + term_nuc - c_annih * SSD_density) * abs(shear_rate) * dtime;
        rho_mov = SSD_density;
        if(SSD_density < rho_init) rho_init = SSD_density;
    }
}

/* Update acc_strain or SSD_density */
void Slip::update_rho_hard_old(vector<PMode*> mode_sys){
    rho_H = SSD_density;
    double coplanar = 0; int active_num = 0;
    bool open = true;
    for (auto &isys : mode_sys) {
        if (isys->type != slip) {open = false; continue;}
        active_num += isys->SSD_density > isys->rho_init ? 1 : 0;
        if (isys->num == num) continue;
        int inter_mode = get_interaction_mode(burgers_vec, plane_norm, isys->burgers_vec, isys->plane_norm);
        if (inter_mode != 2) continue;
        if (coplanar == 0.) coplanar += isys->SSD_density;
        else coplanar = min(coplanar, isys->SSD_density);
    }
    open = open && (active_num > 0);
    if (open) rho_H += coplanar;
}

void Slip::update_ssd(Matrix3d strain_rate, Matrix3d orientation){
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
    if (flag_harden == 0) acc_strain += abs(shear_rate) * dtime;
    if (flag_harden == 1){ 
        double c_forest = harden_params[8], c_nuc = harden_params[9], tau_nuc = harden_params[10],\
               c_multi = harden_params[11], c_annih = 0.,\
               D = harden_params[12] * 1e6, ref_srate = harden_params[13], gg = c_forest/harden_params[14],\
               burgers = update_params[0], mfp = update_params[1], forest_stress = update_params[3]; 
        SSD_density = rho_H;
        /* double equi_strain_rate = strain_rate(2,2); */
        /* double equi_strain_rate = calc_first_principal(strain_rate); */
        double equi_strain_rate = calc_equivalent_value(strain_rate);
        rho_sat = c_forest * burgers / gg * (1-k_boltzmann * temperature/D/pow(burgers,3) * log(abs(equi_strain_rate)/ref_srate));
        rho_sat = max(pow(1/rho_sat,2), 0.5*SSD_density);
        double term_nuc = c_nuc * max(abs(rss)-tau_nuc,0.) / (shear_modulus * burgers * burgers);
        double term_multi = c_multi / mfp; 
        c_annih = (term_multi + term_nuc) / rho_sat;
        SSD_density += (term_multi + term_nuc - c_annih * SSD_density) * abs(shear_rate) * dtime;
        rho_mov = SSD_density;
        if(SSD_density < rho_init) rho_init = SSD_density;
    }
}

void Slip::update_rho_hard(vector<PMode*> mode_sys){
    double d_term_coplanar = 0, coeff_coplanar = harden_params[15];
    vector<double> coplanar_indices;
    for (auto &isys : mode_sys) {
        if (isys->type != slip) continue;
        if (isys->num == num) continue;
        int inter_mode = get_interaction_mode(burgers_vec, plane_norm, isys->burgers_vec, isys->plane_norm);
        if (inter_mode != 2) continue;
        coplanar_indices.push_back(isys->num);
    }
    for (auto &index_1 : coplanar_indices) {
        for (auto &index_2 : coplanar_indices) {
            if (index_1 >= index_2) continue;
            double plus_term = mode_sys[index_1]->SSD_density * mode_sys[index_1]->disl_vel * sqrt(mode_sys[index_2]->SSD_density) + \
                               mode_sys[index_2]->SSD_density * mode_sys[index_2]->disl_vel * sqrt(mode_sys[index_1]->SSD_density);
            d_term_coplanar += plus_term;
        }
        double minus_term = mode_sys[index_1]->SSD_density * mode_sys[index_1]->disl_vel * sqrt(SSD_density) + \
                            SSD_density * disl_vel * sqrt(mode_sys[index_1]->SSD_density);
        d_term_coplanar -= minus_term;
    }
    d_term_coplanar = d_term_coplanar * coeff_coplanar * dtime;
    if (d_term_coplanar + SSD_density < 0) d_term_coplanar = -SSD_density * 0.5;
    /* custom_var = d_term_coplanar; */
    rho_H = SSD_density + d_term_coplanar;
    rho_init = min(rho_H, rho_init);
}
