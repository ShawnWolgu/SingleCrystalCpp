#include "singleX.h"

double waiting_time(double rss, double freq_Debye, double c_length, double burgers, double disl_density_for, double kink_energy_ref, double back_stress,\
                    double Peierls_stress, double expo_kinkeng, double temperature_ref);
double running_time(double rss, double c_drag, double wave_speed, double barrier_distance, double burgers, double back_stress);
vector<double> waiting_time_grad(double rss, double freq_Debye, double c_length, double burgers, double disl_density_for, double kink_energy_ref, double back_stress,\
                    double Peierls_stress, double expo_kinkeng, double temperature_ref);
vector<double> running_time_grad(double rss, double c_drag, double wave_speed, double barrier_distance, double burgers, double back_stress);


double Slip::disl_velocity(double rss){
    /*
     * harden parameters: 0: SSD_density,
     * 1: freq_Debye, 2: c_length, 3: kink_energy_ref, 4: temperature_ref,
     * 5: Peierls_stress, 6: expo_kinkeng, 7: wave_speed, 8: c_drag, 9: c_backstress,
     * 10: c_multi, 11:c_annih, 12:HP_stress.
     * 
     * update parameters:
     * 0: burgers, 1: disl_density_for, 2: disl_density_perp, 3: back_stress,
     * 4: barrier_distance
     */
    double freq_Debye = harden_params[1], c_length = harden_params[2], kink_energy_ref = harden_params[3],\
           temperature_ref = harden_params[4], Peierls_stress = harden_params[5], expo_kinkeng = harden_params[6],\
           wave_speed = harden_params[7], c_drag = harden_params[8];
    double burgers = update_params[0], disl_density_for = update_params[1],\
           back_stress = update_params[3], barrier_distance = update_params[4];
    if(abs(rss)-back_stress > 1e-20){
        t_wait = waiting_time(rss, freq_Debye, c_length, burgers, disl_density_for, kink_energy_ref, back_stress,\
                    Peierls_stress, expo_kinkeng, temperature_ref); 
        t_run = running_time(rss, c_drag, wave_speed, barrier_distance, burgers, back_stress);
//	cout << "t_wait, t_run = " << t_wait << "," << t_run << endl;
        return barrier_distance / (t_wait + t_run);
    }
    else{return 0.0;}
}

double waiting_time(double rss, double freq_Debye, double c_length, double burgers, double disl_density_for, double kink_energy_ref, double back_stress,\
                    double Peierls_stress, double expo_kinkeng, double temperature_ref){
    rss = rss* MPa_to_Pa, back_stress = back_stress*MPa_to_Pa, Peierls_stress = Peierls_stress*MPa_to_Pa, kink_energy_ref = kink_energy_ref*eV_to_J;
    double stress_norm = min((abs(rss) - back_stress)/Peierls_stress,1.00);
    double kink_energy = kink_energy_ref;
    double kink_volume = 400 * pow(burgers,3);
    double freq = freq_Debye * pow(burgers,2) * (c_length / sqrt(disl_density_for))/abs(kink_volume+1e-30) * exp(-kink_energy/(k_boltzmann * temperature)) * sinh(stress_norm);
//    cout << "Qs, Vs, Freq = " << kink_energy / eV_to_J << ", " << kink_volume * 1e30 << ", " <<  freq_Debye * pow(burgers,2) * (c_length / sqrt(disl_density_for))/kink_volume << endl;
    return 1 / max(freq, 1e-20);
}
double waiting_time_self(double rss, double freq_Debye, double c_length, double burgers, double disl_density_for, double kink_energy_ref, double back_stress,\
                    double Peierls_stress, double expo_kinkeng, double temperature_ref){
    rss = rss* MPa_to_Pa, back_stress = back_stress*MPa_to_Pa, Peierls_stress = Peierls_stress*MPa_to_Pa, kink_energy_ref = kink_energy_ref*eV_to_J;
    double stress_norm = min((abs(rss) - back_stress)/Peierls_stress,1.00);
    double kink_energy = kink_energy_ref * (1-temperature/temperature_ref) * pow(1-stress_norm,expo_kinkeng);
    double kink_volume = expo_kinkeng * kink_energy_ref * (1-temperature/temperature_ref) * pow(1-stress_norm,expo_kinkeng-1)/Peierls_stress;
    double freq = freq_Debye * pow(burgers,2) * (c_length / sqrt(disl_density_for))/abs(kink_volume+1e-30) * exp(-kink_energy/(k_boltzmann * temperature));
//    cout << "Qs, Vs, Freq = " << kink_energy / eV_to_J << ", " << kink_volume * 1e30 << ", " <<  freq_Debye * pow(burgers,2) * (c_length / sqrt(disl_density_for))/kink_volume << endl;
    return 1 / max(freq, 1e-20);
}

double running_time(double rss, double c_drag, double wave_speed, double barrier_distance, double burgers, double back_stress){
    rss = rss* MPa_to_Pa, back_stress = back_stress*MPa_to_Pa;
    double coeff_B = (c_drag * k_boltzmann * temperature) / (wave_speed * burgers * burgers);
    double velocity = coeff_B * wave_speed/(2*(abs(rss) - back_stress) * burgers);
    velocity = wave_speed*(sqrt(velocity*velocity+1)-velocity);
    return barrier_distance / velocity;
}

vector<double> disl_velocity_grad(double rss, vector<double> harden_params, vector<double> update_params){
    /*
     * Also calculate the velocity.
     * harden parameters: 0: SSD_density,
     * 1: freq_Debye, 2: c_length, 3: kink_energy_ref, 4: temperature_ref,
     * 5: Peierls_stress, 6: expo_kinkeng, 7: wave_speed, 8: c_drag, 9: c_backstress,
     * 10: c_multi, 11:c_annih, 12:HP_stress.
     * 
     * update parameters:
     * 0: burgers, 1: disl_density_for, 2: disl_density_perp, 3: back_stress,
     * 4: barrier_distance
     */
    double freq_Debye = harden_params[1], c_length = harden_params[2], kink_energy_ref = harden_params[3],\
           temperature_ref = harden_params[4], Peierls_stress = harden_params[5], expo_kinkeng = harden_params[6],\
           wave_speed = harden_params[7], c_drag = harden_params[8];
    double burgers = update_params[0], disl_density_for = update_params[1],\
           back_stress = update_params[3], barrier_distance = update_params[4];
    if(abs(rss)-back_stress > 1e-20){
        vector<double> dtwait_drss = waiting_time_grad(rss, freq_Debye, c_length, burgers, disl_density_for, kink_energy_ref, back_stress,\
                    Peierls_stress, expo_kinkeng, temperature_ref);
        vector<double> dtrun_drss = running_time_grad(rss, c_drag, wave_speed, barrier_distance, burgers, back_stress);
        double dvel_dtau = -barrier_distance / pow((dtwait_drss[1] + dtrun_drss[1]),2) * (dtwait_drss[0] + dtrun_drss[0]);
        double velocity = barrier_distance / (dtwait_drss[1] + dtrun_drss[1]);
        vector<double> result = {dvel_dtau,velocity};
        return result;
    }
    else{
        vector<double> result = {0.00,0.00};
        return result;
        }
}

vector<double> waiting_time_grad_self(double rss, double freq_Debye, double c_length, double burgers, double disl_density_for, double kink_energy_ref, double back_stress,\
                    double Peierls_stress, double expo_kinkeng, double temperature_ref){
    /* Return a vector: 0. dtw/dtau, 1. tw. */
    rss = rss* MPa_to_Pa, back_stress = back_stress*MPa_to_Pa, Peierls_stress = Peierls_stress*MPa_to_Pa, kink_energy_ref = kink_energy_ref*eV_to_J;
    double stress_norm = min((abs(rss) - back_stress)/Peierls_stress,1.00);
    double dvolume_dtau = -1 * sign(rss) * expo_kinkeng * (expo_kinkeng - 1) * kink_energy_ref * (1-temperature/temperature_ref) *\
                          pow(1-stress_norm,expo_kinkeng-2) / pow(Peierls_stress,2);
    double kink_energy = kink_energy_ref * (1-temperature/temperature_ref) * pow(1-stress_norm,expo_kinkeng);
    double kink_volume = expo_kinkeng * kink_energy_ref * (1-temperature/temperature_ref) * pow(1-stress_norm,expo_kinkeng-1)/Peierls_stress;
    double freq_const = freq_Debye * pow(burgers,2) * (c_length / sqrt(disl_density_for)); 
    double dtw_dtau = 1/freq_const * exp(kink_energy/k_boltzmann/temperature) * (dvolume_dtau - kink_volume * kink_volume / k_boltzmann/temperature);
    double freq = freq_const/abs(kink_volume+1e-30) * exp(-kink_energy/(k_boltzmann * temperature));
    vector<double> result ={ dtw_dtau*MPa_to_Pa, 1/max(freq, 1e-25)};
    return result;
}

vector<double> waiting_time_grad(double rss, double freq_Debye, double c_length, double burgers, double disl_density_for, double kink_energy_ref, double back_stress,\
                    double Peierls_stress, double expo_kinkeng, double temperature_ref){
    /* Return a vector: 0. dtw/dtau, 1. tw. */
    rss = rss* MPa_to_Pa, back_stress = back_stress*MPa_to_Pa, Peierls_stress = Peierls_stress*MPa_to_Pa, kink_energy_ref = kink_energy_ref*eV_to_J;
    double stress_norm = min((abs(rss) - back_stress)/Peierls_stress,1.00);
    double kink_energy = kink_energy_ref;
    double kink_volume = 400 * pow(burgers,3);
    double freq_const = freq_Debye * pow(burgers,2) * (c_length / sqrt(disl_density_for)) / abs(kink_volume + 1e-30); 
    double freq = freq_const * exp(-kink_energy/(k_boltzmann * temperature)) * sinh(stress_norm);
    double dtw_dtau = 1/freq_const * exp(kink_energy/k_boltzmann/temperature) * (-1 * cosh(stress_norm) / pow(sinh(stress_norm),2) * sign(rss) / Peierls_stress);
    vector<double> result ={ dtw_dtau*MPa_to_Pa, 1/max(freq, 1e-25)};
    return result;
}

vector<double> running_time_grad(double rss, double c_drag, double wave_speed, double barrier_distance, double burgers, double back_stress){
    /* Return a vector: 0. dtr/dtau, 1. tr. */
    rss = rss* MPa_to_Pa, back_stress = back_stress*MPa_to_Pa;
    double coeff_B = (c_drag * k_boltzmann * temperature) / (wave_speed * burgers * burgers);
    double vs_by_vel_m = (coeff_B * wave_speed) / (2*(abs(rss) - back_stress)*burgers);
    double dvm_dtau = sign(rss) * 2 * burgers / coeff_B;
    double dvd_dtau = dvm_dtau * pow(vs_by_vel_m,2) * (1 - vs_by_vel_m / sqrt(1+pow(vs_by_vel_m,2)));
    double velocity = wave_speed*(sqrt(vs_by_vel_m*vs_by_vel_m+1)-vs_by_vel_m);
    velocity = max(1e-10,velocity);
    vector<double> result = { -barrier_distance/pow(velocity,2)*dvd_dtau*MPa_to_Pa, barrier_distance/velocity};
    return result;
}
