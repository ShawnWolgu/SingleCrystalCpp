#include "singleX.h"
#include <cmath>

double waiting_time(double stress_eff, double resistance_slip, double act_energy_r, double frequency_r, double energy_expo, double temperature);
double running_time(double stress_eff, double c_drag, double speed_sat, double mean_free_path, double burgers, double temperature);
vector<double> waiting_time_grad(double stress_eff, double resistance_slip, double act_energy_r, double frequency_r, double energy_expo, double temperature);
vector<double> running_time_grad(double stress_eff, double c_drag, double speed_sat, double mean_free_path, double burgers, double temperature);

double Slip::disl_velocity(double rss){
    /*
     * [velocity parameters] 
     *  1. MFP control coeffient, 2. reference frequency, 3. activation energy, 4. slip resistance, 5. energy exponent
     *  6. saturated speed, 7. drag coefficient
     * [hardening parameters] 
     *  8. forest hardening coefficient
     * [DD evolution parameters] 
     *  0. SSD_density, 9. multiplication coefficient, 10. drag stress D, 11. reference strain rate, 12. c/g 
     *
     * update parameters:
     * 0: burgers, 1: mean_free_path, 2: disl_density_resist, 3: forest_stress,
     */
    double frequency_r = harden_params[2], act_energy_r = harden_params[3], resistance_slip = harden_params[4], \
           energy_expo = harden_params[5], speed_sat = harden_params[6], c_drag = harden_params[7];
    double burgers = update_params[0], mean_free_path = update_params[1], forest_stress = update_params[3];
    double stress_eff = abs(rss) - forest_stress;
    if (stress_eff >= 0.0) {
        t_wait = waiting_time(stress_eff, resistance_slip, act_energy_r, frequency_r, energy_expo, temperature);
        t_run = running_time(stress_eff, c_drag, speed_sat, mean_free_path, burgers, temperature);
        return mean_free_path / (t_wait + t_run);
    }
    else{
        t_wait = waiting_time(stress_eff, resistance_slip, act_energy_r, frequency_r, energy_expo, temperature);
        return mean_free_path / t_wait;
    }
}

vector<double> Slip::disl_velocity_grad(double rss){
    /*
     * [velocity parameters] 
     *  1. MFP control coeffient, 2. reference frequency, 3. activation energy, 4. slip resistance, 5. energy exponent
     *  6. saturated speed, 7. drag coefficient
     * [hardening parameters] 
     *  8. forest hardening coefficient
     * [DD evolution parameters] 
     *  0. SSD_density, 9. multiplication coefficient, 10. drag stressD, 11. reference strain rate, 12. c/g 
     *
     * update parameters:
     * 0: burgers, 1: mean_free_path, 2: disl_density_resist, 3: forest_stress,
     */
    double frequency_r = harden_params[2], act_energy_r = harden_params[3], resistance_slip = harden_params[4], \
           energy_expo = harden_params[5], speed_sat = harden_params[6], c_drag = harden_params[7];
    double burgers = update_params[0], mean_free_path = update_params[1], forest_stress = update_params[3];
    double stress_eff = abs(rss) - forest_stress;
    vector<double> result;
    if (stress_eff >= 0.0) {
        vector<double> dtwait_drss = waiting_time_grad(stress_eff, resistance_slip, act_energy_r, frequency_r, energy_expo, temperature);
        vector<double> dtrun_drss = running_time_grad(stress_eff, c_drag, speed_sat, mean_free_path, burgers, temperature);
        double dvel_dtau = -sign(rss) * mean_free_path / pow((dtwait_drss[1] + dtrun_drss[1]),2) * (dtwait_drss[0] + dtrun_drss[0]);
        double velocity = mean_free_path / (dtwait_drss[1] + dtrun_drss[1]);
        result = {dvel_dtau,velocity};
    }
    else{
        vector<double> dtwait_drss = waiting_time_grad(stress_eff, resistance_slip, act_energy_r, frequency_r, energy_expo, temperature);
        double dvel_dtau = -sign(rss) * mean_free_path / pow(dtwait_drss[1],2) * dtwait_drss[0];
        double velocity = mean_free_path / dtwait_drss[1];
        result = {dvel_dtau,velocity};
    }
    return result;
}

double waiting_time(double stress_eff, double resistance_slip, double act_energy_r, double frequency_r, double energy_expo, double temperature){
    stress_eff = stress_eff * MPa_to_Pa, resistance_slip = resistance_slip * MPa_to_Pa, act_energy_r = act_energy_r * eV_to_J;
    double act_energy = 0.;
    if (stress_eff >= 0.0) act_energy = act_energy_r * (1-pow((stress_eff/resistance_slip), energy_expo)); // activation energy
    else act_energy = act_energy_r * (1+pow((abs(stress_eff)/resistance_slip), energy_expo)); // activation energy
    act_energy = min(act_energy, 500 * k_boltzmann * temperature); // avoid too large activation energy;
    act_energy = max(act_energy, -500 * k_boltzmann * temperature); // avoid too small activation energy;
    return 1 / frequency_r * exp(act_energy / (k_boltzmann * temperature));
}

vector<double> waiting_time_grad(double stress_eff, double resistance_slip, double act_energy_r, double frequency_r, double energy_expo, double temperature){
    /* Return a vector: 0. dtw/dtau, 1. tw. */
    stress_eff = stress_eff * MPa_to_Pa, resistance_slip = resistance_slip * MPa_to_Pa, act_energy_r = act_energy_r * eV_to_J;
    double act_energy = 0.;
    if (stress_eff >= 0.0) act_energy = act_energy_r * (1-pow((stress_eff/resistance_slip), energy_expo)); // activation energy
    else act_energy = act_energy_r * (1+pow((abs(stress_eff)/resistance_slip), energy_expo)); // activation energy
    act_energy = min(act_energy, 500 * k_boltzmann * temperature); // avoid too large activation energy;
    act_energy = max(act_energy, -500 * k_boltzmann * temperature); // avoid too small activation energy;
    double waiting_time = 1 / frequency_r * exp(act_energy / (k_boltzmann * temperature));
    double grad_term = 0;
    if (act_energy != 500*k_boltzmann*temperature && act_energy != -500*k_boltzmann*temperature) 
        grad_term = - act_energy_r * energy_expo * pow((abs(stress_eff)/resistance_slip), energy_expo-1) / (k_boltzmann * temperature);
    vector<double> result ={ grad_term*waiting_time*MPa_to_Pa/resistance_slip, waiting_time };
    return result;
}

double running_time(double stress_eff, double c_drag, double speed_sat, double mean_free_path, double burgers, double temperature){
    stress_eff = stress_eff * MPa_to_Pa;
    double coeff_B = (c_drag * k_boltzmann * temperature) / (speed_sat * pow(burgers,2));
    /* double coeff_B = (c_drag * k_boltzmann * temperature) / (speed_sat * burgers * mean_free_path); */
    double v_norm = 2 * burgers * stress_eff / (coeff_B * speed_sat);
    double velocity = speed_sat * (sqrt(1 + pow(v_norm, -2)) - 1/v_norm);
    return mean_free_path / max(velocity, 1e-40);
}

vector<double> running_time_grad(double stress_eff, double c_drag, double speed_sat, double mean_free_path, double burgers, double temperature){
    /* Return a vector: 0. dtr/dtau, 1. tr. */
    stress_eff = stress_eff * MPa_to_Pa;
    double coeff_B = (c_drag * k_boltzmann * temperature) / (speed_sat * burgers * burgers);
    /* double coeff_B = (c_drag * k_boltzmann * temperature) / (speed_sat * burgers * mean_free_path); */
    double v_norm = 2 * burgers * stress_eff / (coeff_B * speed_sat);
    double velocity = speed_sat * (sqrt(1 + 1/pow(v_norm,2))-1/v_norm); 
    double gradient = -1 * (2*burgers*mean_free_path)/(pow(velocity,2)*coeff_B) * pow(v_norm,-2) * (1-1/(v_norm * sqrt(1+pow(v_norm,-2))));
    vector<double> result = {gradient*MPa_to_Pa, mean_free_path/max(velocity,1e-40)};
    return result;
}
