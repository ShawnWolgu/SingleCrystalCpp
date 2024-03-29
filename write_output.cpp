#include "singleX.h"

void custom_output_initialization(Grain &grain);
void custom_output_initialization();
void print_custom(Grain &grain);

void title_output(ofstream &outf, string leads, string subs, int length);
void title_output(ofstream &outf, string leads, string subs1, string subs2, int length);
double cal_std_013467910(vector<double> vd);
double cal_mean_013467910(vector<double> vd);

void outfile_initialization(){
    cout << "Output Files Initialing " << endl;
    stress_step_file << "time,e11,e22,e33,e12,e13,e23,s11,s22,s33,s12,s13,s23" << endl;
    stress_file << "time,e11,e22,e33,e12,e13,e23,s11,s22,s33,s12,s13,s23" << endl;
    disloc_file << "time,e11,e22,e33,dd1,dd2,dd3,...," << endl;
    crss_file << "e11,e22,e33,acc_strain,crss1,crss2,crss3,..," <<endl;
    rss_file << "e11,e22,e33,acc_strain,rss1,rss2,rss3,..," <<endl;
    accstrain_file << "e11,e22,e33,acc_strain1,acc_strain2,acc_strain3,...," <<endl;
    disloc_step_file << "time,e11,e22,e33,dd1,dd2,dd3,...," << endl; 
    euler_file << "phi1,PHI,phi2" << endl;
    schmidt_file << "e11,e22,e33,sf1,sf2,sf3,...," << endl;
    disvel_file << "e11,e22,e33,vel1,vel2,vel3,...," << endl;
    custom_output_initialization();
    cout << "Finish Initialization." << endl;
}

void outfile_initialization(Grain &grain){
    cout << "Output Files Initialing (Grain) ...." << endl;
    int slip_num =  grain.mode_sys.size();

    stress_step_file << "time,e11,e22,e33,e12,e13,e23,s11,s22,s33,s12,s13,s23,pe" << endl;
    grain.print_stress_strain(stress_step_file);
    stress_file << "time,e11,e22,e33,e12,e13,e23,s11,s22,s33,s12,s13,s23,pe" << endl;
    grain.print_stress_strain(stress_file);
    cout << "This grain has " << slip_num << " slip systems" << endl;
    title_output(crss_file,"e11,e22,e33,acc_strain,","crss", slip_num);
    grain.print_crss(crss_file);
    title_output(rss_file,"e11,e22,e33,acc_strain,","crss", slip_num);
    grain.print_rss(rss_file);
    title_output(accstrain_file,"e11,e22,e33,","acc_strain", slip_num);
    grain.print_accstrain(accstrain_file);
    if (flag_harden == 1) {
        title_output(time_step_file, "time,e11,e22,e33,", "tw","tr", slip_num);
        grain.print_time(time_step_file);
        title_output(disloc_file, "time,e11,e22,e33,", "ssd", slip_num);
        grain.print_dislocation(disloc_file);
        title_output(disloc_step_file, "time,e11,e22,e33,", "ssd", slip_num);
        grain.print_dislocation(disloc_step_file);
        custom_output_initialization(grain);
    };
    euler_file << "phi1,PHI,phi2" << endl;
    grain.print_euler(euler_file);
    title_output(schmidt_file, "e11,e22,e33,", "sf", slip_num);
    grain.print_schmidt(schmidt_file);
    title_output(disvel_file, "e11,e22,e33,", "vel", slip_num);
    grain.print_disvel(disvel_file);
    cout << "Finish Initialization." << endl;
}

void outfile_close(){
    stress_step_file.close(); 
    disloc_step_file.close();
    stress_file.close();
    disloc_file.close();
    crss_file.close();
    rss_file.close();
    euler_file.close();
    custom_output_file.close();
    accstrain_file.close();
}

void substep_output(Grain &grain){
    grain.print_stress_strain(stress_step_file);
    grain.print_dislocation(disloc_step_file);
    if (flag_harden == 1) grain.print_time(time_step_file);
    grain.print_crss(crss_file);
    grain.print_disvel(disvel_file);
    grain.print_schmidt(schmidt_file);
    print_custom(grain);
}

void grain_output(Grain &grain){
    grain.print_crss(crss_file);
    grain.print_stress_strain(stress_file);
    grain.print_dislocation(disloc_file);
    grain.print_stress_strain_screen();
    grain.print_euler(euler_file);
    grain.print_accstrain(accstrain_file);
    grain.print_rss(rss_file);
    grain.print_schmidt(schmidt_file);
    grain.print_disvel(disvel_file);
    if (flag_harden == 1) grain.print_time(time_step_file);
    print_custom(grain);
}

void title_output(ofstream &outf, string leads, string subs, int length){
    outf << leads;
    for(int i = 0; i != length; ++i){
        string dtitle = subs + to_string(i);
        outf << dtitle << ",";
    }
    outf << endl;
}

void title_output(ofstream &outf, string leads, string subs1, string subs2, int length){
    outf << leads;
    for(int i = 0; i != length; ++i){
        string dtitle = subs1 + to_string(i);
        outf << dtitle << ",";
        dtitle = subs2 + to_string(i);
        outf << dtitle << ",";
    }
    outf << endl;
}

// custom output settings:
void custom_output_initialization(){
    //title_output(custom_output_file, "e11,e22,e33,", "rate", 12);
    //custom_output_file << "e11,e22,e33,U11,U12,U13,U21,U22,U23,U31,U32,U33" << endl;
    custom_output_file << "e11,e22,e33,std_ssd,std_b,std_rate" << endl;
}

void custom_output_initialization(Grain &grain){
    title_output(custom_output_file, "e11,e22,e33,", "rhosat", 12);
    //custom_output_file << "e11,e22,e33,U11,U12,U13,U21,U22,U23,U31,U32,U33" << endl;
    //custom_output_file << "e11,e22,e33,std_ssd,std_b,std_rate" << endl;
    print_custom(grain);
}

void print_custom(Grain &grain){
    custom_output_file << grain.strain_tensor(0,0) << ',' << grain.strain_tensor(1,1) << ',' << grain.strain_tensor(2,2);// << ',' << grain.slip_sys[4].update_params[0];
    //Matrix3d U = (grain.orientation.transpose() * grain.orient_ref).inverse()*grain.deform_grad_elas;
    for (auto &mode_component : grain.mode_sys) {
	custom_output_file << ',' << (mode_component->custom_var);// cross_in - slip_component.cross_out;
    }
    custom_output_file << endl;
}

double cal_std_013467910(vector<double> vd){
    double mean = cal_mean_013467910(vd);
    int idx[8] = {0,1,3,4,6,7,9,10};
    double add_number = 0.0;
    for (auto &iid : idx) add_number += pow((abs(vd[iid])-mean),2);
    return sqrt(0.125*add_number)/mean;
}

double cal_mean_013467910(vector<double> vd){
    return 0.125 * (abs(vd[0]) + abs(vd[1]) + abs(vd[3])+ abs(vd[4]) + abs(vd[6]) + abs(vd[7]) + abs(vd[9]) + abs(vd[10]));
}
