#include "singleX.h"

double waiting_time(double rss, double freq_Debye, double c_length, double burgers, double disl_density_for, double kink_energy_ref, double back_stress,\
                    double Peierls_stress, double expo_kinkeng, double temperature_ref);
double running_time(double rss, double c_drag, double wave_speed, double barrier_distance, double burgers, double back_stress);
vector<double> waiting_time_grad(double rss, double freq_Debye, double c_length, double burgers, double disl_density_for, double kink_energy_ref, double back_stress,\
                    double Peierls_stress, double expo_kinkeng, double temperature_ref);
vector<double> running_time_grad(double rss, double c_drag, double wave_speed, double barrier_distance, double burgers, double back_stress);


double disl_velocity(double rss, vector<double> harden_params, vector<double> update_params){
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
    double freq_Debye = harden_params[1], c_length = harden_params[2], kink_energy_ref = harden_params[3] * eV_to_J,\
           temperature_ref = harden_params[4], Peierls_stress = harden_params[5], expo_kinkeng = harden_params[6],\
           wave_speed = harden_params[7], c_drag = harden_params[8];
    double burgers = update_params[0], disl_density_for = update_params[1],\
           back_stress = update_params[3], barrier_distance = update_params[4];
    if(abs(rss)-back_stress > 1e-20){
        double t_wait = waiting_time(rss, freq_Debye, c_length, burgers, disl_density_for, kink_energy_ref, back_stress,\
                    Peierls_stress, expo_kinkeng, temperature_ref); 
        double t_run = running_time(rss, c_drag, wave_speed, barrier_distance, burgers, back_stress);
        return barrier_distance / (t_wait + t_run);
    }
    else{return 0.0;}
}

double waiting_time(double rss, double freq_Debye, double c_length, double burgers, double disl_density_for, double kink_energy_ref, double back_stress,\
                    double Peierls_stress, double expo_kinkeng, double temperature_ref){
    rss = rss* MPa_to_Pa, back_stress = back_stress*MPa_to_Pa, Peierls_stress = Peierls_stress*MPa_to_Pa, kink_energy_ref = kink_energy_ref*eV_to_J;
    double stress_norm = min((abs(rss) - back_stress)/Peierls_stress,1 - 1e-3);
    double kink_energy = kink_energy_ref * (1-temperature/temperature_ref) * pow(1-stress_norm,expo_kinkeng);
    double kink_volume = expo_kinkeng * kink_energy_ref * (1-temperature/temperature_ref) * pow(1-stress_norm,expo_kinkeng-1)/Peierls_stress;
    double freq = freq_Debye * pow(burgers,2) * (c_length / sqrt(disl_density_for))/kink_volume * exp(-kink_energy/(k_boltzmann * temperature));
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

vector<double> waiting_time_grad(double rss, double freq_Debye, double c_length, double burgers, double disl_density_for, double kink_energy_ref, double back_stress,\
                    double Peierls_stress, double expo_kinkeng, double temperature_ref){
    /* Return a vector: 0. dtw/dtau, 1. tw. */
    rss = rss* MPa_to_Pa, back_stress = back_stress*MPa_to_Pa, Peierls_stress = Peierls_stress*MPa_to_Pa, kink_energy_ref = kink_energy_ref*eV_to_J;
    double stress_norm = min((abs(rss) - back_stress)/Peierls_stress,1 - 1e-3);
    double dvolume_dtau = sign(rss) * expo_kinkeng * (expo_kinkeng - 1) * kink_energy_ref * (1-temperature/temperature_ref) *\
                          pow(1-stress_norm,expo_kinkeng-2) / pow(Peierls_stress,2);
    double kink_energy = kink_energy_ref * (1-temperature/temperature_ref) * pow(1-stress_norm,expo_kinkeng);
    double kink_volume = expo_kinkeng * kink_energy_ref * (1-temperature/temperature_ref) * pow(1-stress_norm,expo_kinkeng-1)/Peierls_stress;
    double freqterm = freq_Debye * pow(burgers,2) * (c_length / sqrt(disl_density_for)) * (-dvolume_dtau/pow(kink_volume,2) + 1/(k_boltzmann*temperature)) *\
                  exp(-kink_energy/(k_boltzmann * temperature));
    double freq = freq_Debye * pow(burgers,2) * (c_length / sqrt(disl_density_for))/kink_volume * exp(-kink_energy/(k_boltzmann * temperature));
    vector<double> result ={ -1/pow(freqterm,2)*MPa_to_Pa, 1/max(freq, 1e-20)};
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
