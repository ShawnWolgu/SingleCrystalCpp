#include "singleX.h"

void custom_output_initialization();
void print_custom(Grain &grain);

void title_output(ofstream &outf, string leads, string subs, int length);
void title_output(ofstream &outf, string leads, string subs1, string subs2, int length);

void outfile_initialization(){
    cout << "Output Files Initialing " << endl;
    stress_step_file << "time,e11,e22,e33,e12,e13,e23,s11,s22,s33,s12,s13,s23" << endl;
    stress_file << "time,e11,e22,e33,e12,e13,e23,s11,s22,s33,s12,s13,s23" << endl;
    disloc_file << "time,e11,e22,e33,dd1,dd2,dd3,...," << endl;
    crss_file << "e11,e22,e33,acc_strain,crss1,crss2,crss3,..," <<endl;
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
    int slip_num =  grain.slip_sys.size();

    stress_step_file << "time,e11,e22,e33,e12,e13,e23,s11,s22,s33,s12,s13,s23" << endl;
    grain.print_stress_strain(stress_step_file);
    stress_file << "time,e11,e22,e33,e12,e13,e23,s11,s22,s33,s12,s13,s23" << endl;
    grain.print_stress_strain(stress_file);
    cout << "This grain has " << slip_num << " slip systems" << endl;
    title_output(disloc_file, "time,e11,e22,e33,", "ssd", slip_num);
    grain.print_dislocation(disloc_file);
    title_output(crss_file,"e11,e22,e33,acc_strain,","crss", slip_num);
    grain.print_crss(crss_file);
    title_output(accstrain_file,"e11,e22,e33,","acc_strain", slip_num);
    grain.print_accstrain(accstrain_file);
    title_output(disloc_step_file, "time,e11,e22,e33,", "ssd", slip_num);
    grain.print_dislocation(disloc_step_file);
    if (flag_harden == 2) {
        title_output(time_step_file, "time,e11,e22,e33,", "tw","tr", slip_num);
        grain.print_time(time_step_file);
    };
    euler_file << "phi1,PHI,phi2" << endl;
    grain.print_euler(euler_file);
    title_output(schmidt_file, "e11,e22,e33,", "sf", slip_num);
    grain.print_schmidt(schmidt_file);
    title_output(disvel_file, "e11,e22,e33,", "vel", slip_num);
    grain.print_disvel(disvel_file);
    custom_output_initialization();
    cout << "Finish Initialization." << endl;
}

void outfile_close(){
    stress_step_file.close(); 
    disloc_step_file.close();
    stress_file.close();
    disloc_file.close();
    crss_file.close();
    euler_file.close();
    custom_output_file.close();
    accstrain_file.close();
}

void substep_output(Grain &grain){
    grain.print_stress_strain(stress_step_file);
    grain.print_dislocation(disloc_step_file);
    if (flag_harden == 2) grain.print_time(time_step_file);
}

void grain_output(Grain &grain){
    grain.print_crss(crss_file);
    grain.print_stress_strain(stress_file);
    grain.print_dislocation(disloc_file);
    grain.print_stress_strain_screen();
    grain.print_euler(euler_file);
    grain.print_accstrain(accstrain_file);
    grain.print_schmidt(schmidt_file);
    grain.print_disvel(disvel_file);
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
    title_output(custom_output_file, "e11,e22,e33,", "refg", 12);
    //custom_output_file << "e11,e22,e33,U11,U12,U13,U21,U22,U23,U31,U32,U33" << endl;
}

void print_custom(Grain &grain){
    custom_output_file << grain.strain_tensor(0,0) << ',' << grain.strain_tensor(1,1) << ',' << grain.strain_tensor(2,2);// << ',' << grain.slip_sys[4].update_params[0];
    //Matrix3d U = (grain.orientation.transpose() * grain.orient_ref).inverse()*grain.deform_grad_elas;
    for (Slip &slip_component : grain.slip_sys) custom_output_file << ',' << slip_component.strain_rate_slip;
    custom_output_file << endl;
}
