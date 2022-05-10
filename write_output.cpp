#include "singleX.h"

void title_output(ofstream &outf, string leads, string subs, int length);

void outfile_initialization(){
    stress_step_file << "e11,e22,e33,e12,e13,e23,s11,s22,s33,s12,s13,s23" << endl;
    stress_file << "e11,e22,e33,e12,e13,e23,s11,s22,s33,s12,s13,s23" << endl;
    disloc_file << "e11,e22,e33,dd1,dd2,dd3,...," << endl;
    crss_file << "e11,e22,e33,acc_strain,crss1,crss2,crss3,..," <<endl;
    disloc_step_file << "e11,e22,e33,dd1,dd2,dd3,...," << endl;  
}

void outfile_initialization(Grain &grain){
    int slip_num =  grain.slip_sys.size();

    stress_step_file << "e11,e22,e33,e12,e13,e23,s11,s22,s33,s12,s13,s23" << endl;
    grain.print_stress_strain(stress_step_file);
    stress_file << "e11,e22,e33,e12,e13,e23,s11,s22,s33,s12,s13,s23" << endl;
    grain.print_stress_strain(stress_file);
    title_output(disloc_file, "e11,e22,e33,", "ssd", slip_num);
    grain.print_dislocation(disloc_file);
    title_output(crss_file,"e11,e22,e33,acc_strain,","crss", slip_num);
    grain.print_crss(crss_file);
    title_output(disloc_step_file, "e11,e22,e33,", "ssd", slip_num);
    grain.print_dislocation(disloc_step_file);
}

void outfile_close(){
    stress_step_file.close(); 
    disloc_step_file.close();
    stress_file.close();
    disloc_file.close();
    crss_file.close();
}

void substep_output(Grain &grain){
    grain.print_stress_strain(stress_step_file);
    grain.print_dislocation(disloc_step_file);
}

void grain_output(Grain &grain){
    grain.print_crss(crss_file);
    grain.print_stress_strain(stress_file);
    grain.print_dislocation(disloc_file);
    grain.print_stress_strain_screen();
}

void title_output(ofstream &outf, string leads, string subs, int length){
    outf << leads;
    for(int i; i != length; ++i){
        string dtitle = subs + to_string(i);
        outf << dtitle << ",";
    }
    outf << endl;
}