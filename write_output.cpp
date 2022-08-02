#include "singleX.h"

void custom_output_initialization();
void print_custom(Grain &grain);

void title_output(ofstream &outf, string leads, string subs, int length);

void outfile_initialization(){
    cout << "Output Files Initialing " << endl;
    stress_step_file << "e11,e22,e33,e12,e13,e23,s11,s22,s33,s12,s13,s23" << endl;
    stress_file << "e11,e22,e33,e12,e13,e23,s11,s22,s33,s12,s13,s23" << endl;
    disloc_file << "e11,e22,e33,dd1,dd2,dd3,...," << endl;
    crss_file << "e11,e22,e33,acc_strain,crss1,crss2,crss3,..," <<endl;
    accstrain_file << "e11,e22,e33,acc_strain1,acc_strain2,acc_strain3,...," <<endl;
    disloc_step_file << "e11,e22,e33,dd1,dd2,dd3,...," << endl; 
    euler_file << "phi1,PHI,phi2" << endl;
    custom_output_initialization();
    cout << "Finish Initialization." << endl;
}

void outfile_initialization(Grain &grain){
    cout << "Output Files Initialing (Grain) ...." << endl;
    int slip_num =  grain.slip_sys.size();

    stress_step_file << "e11,e22,e33,e12,e13,e23,s11,s22,s33,s12,s13,s23" << endl;
    grain.print_stress_strain(stress_step_file);
    stress_file << "e11,e22,e33,e12,e13,e23,s11,s22,s33,s12,s13,s23" << endl;
    grain.print_stress_strain(stress_file);
    cout << "This grain has " << slip_num << " slip systems" << endl;
    title_output(disloc_file, "e11,e22,e33,", "ssd", slip_num);
    grain.print_dislocation(disloc_file);
    title_output(crss_file,"e11,e22,e33,acc_strain,","crss", slip_num);
    grain.print_crss(crss_file);
    title_output(accstrain_file,"e11,e22,e33,","acc_strain", slip_num);
    grain.print_accstrain(accstrain_file);
    title_output(disloc_step_file, "e11,e22,e33,", "ssd", slip_num);
    grain.print_dislocation(disloc_step_file);
    euler_file << "phi1,PHI,phi2" << endl;
    grain.print_euler(euler_file);
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
}

void grain_output(Grain &grain){
    grain.print_crss(crss_file);
    grain.print_stress_strain(stress_file);
    grain.print_dislocation(disloc_file);
    grain.print_stress_strain_screen();
    grain.print_euler(euler_file);
    grain.print_accstrain(accstrain_file);
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

// custom output settings:
void custom_output_initialization(){
    custom_output_file << "e11,e22,e33,teff1,teff2,teff3,teff4,teff5,teff6,teff7,teff8,teff9,teff10,teff11,teff12" << endl;
}

void print_custom(Grain &grain){
    custom_output_file << grain.strain_tensor(0,0) << ',' << grain.strain_tensor(1,1) << ',' << grain.strain_tensor(2,2) << ',';
    Matrix3d stress_t = grain.orientation * grain.stress_tensor * grain.orientation.transpose();
    Matrix3d elastic_f = grain.orientation * grain.deform_grad_elas * grain.orientation.transpose();
    for (Slip &slip_component : grain.slip_sys){
	custom_output_file << ',' << slip_component.cal_rss(stress_t,elastic_f) - slip_component.update_params[3];
    }
    custom_output_file << endl;
}
// code for printing dislocation velocity
//void custom_output_initialization(){
//    custom_output_file << "e11,e22,e33,v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12" << endl;
//}
//
//void print_custom(Grain &grain){
//    custom_output_file << grain.strain_tensor(0,0) << ',' << grain.strain_tensor(1,1) << ',' << grain.strain_tensor(2,2) << ',';
//    for (Slip &slip_component : grain.slip_sys) custom_output_file << ',' << slip_component.disl_vel;
//    custom_output_file << endl;
//}
