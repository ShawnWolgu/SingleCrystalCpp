#include "singleX.h"

double timestep = 0.001, substep = 0.001, dtime = substep*timestep, m = 0.05, max_strain = 0.2, temperature = 298;
int flag_harden = 0;

void flag_to_idx(Matrix<double, 15, 1> flag, vector<int> &known_idx, vector<int> &unknown_idx){
    int unknown_cur = 0, known_cur = 0;
    for (int idx = 0;idx != 15; ++idx) {
        if(flag(idx) == 0) {unknown_idx[unknown_cur] = idx; ++unknown_cur;}
        else {known_idx[known_cur] = idx; ++known_cur;}
    }
    if (unknown_cur != 6 || known_cur != 9) throw "ERROR: wrong number of unknown parameters in L or dsigma!";
}

void step_config(string in_str){
    stringstream stream(in_str);
    double temp[3]={0,0,0};
    int temp_idx = 0;
    while (!stream.eof()) stream >> temp[temp_idx++];
    timestep = 0.001, substep = 0.001, dtime = substep*timestep, max_strain = 0.2;

    if (temp[0] != 0) {
        timestep = temp[0];
    }

    if (temp[1] != 0) {
        substep = temp[1];
    }

    if (temp[2] != 0) {
        max_strain = temp[2];
    }
    dtime = substep*timestep;
}

Matrix3d load_matrix_input(ifstream &load_file){
    if(load_file.eof()){
        throw "ERROR in load file!";
    } 
    Matrix3d temp_matrix;
    string matrix_str;
    double temp[9]={0,0,0,0,0,0,0,0,0};
    double *temp_ar = temp;
    int line_count = 0;
    while(!load_file.eof() && line_count != 3){  
        getline(load_file,matrix_str);
        if(matrix_str == "") continue;
        stringstream stream(matrix_str);
        int temp_idx = 0;      
        for(; !stream.eof() && temp_idx != 3; temp_idx++) stream >> *(temp_ar++);
        if (temp_idx !=3) throw "ERROR in load file!";
        line_count++;
    }
    if (line_count != 3) throw "ERROR in load file!";
    temp_matrix << temp[0], temp[1], temp[2], temp[3], temp[4], temp[5], temp[6], temp[7], temp[8];
    cout << temp_matrix << endl;
    return temp_matrix;
}

int main(int argc, char **argv)
{
    vector<Slip> fcc_110_111;
    Matrix3d strain_tensor = Matrix3d::Zero(), strain_plastic = Matrix3d::Zero(), strain_elastic =  Matrix3d::Zero(), stress_tensor = Matrix3d::Zero();

    flag_harden = atoi(argv[1]);

    //[Slip system setting]: wait to be extracted:
    const Matrix6d elastic_modulus{
        {251700, 107900, 107900, 0, 0, 0},
        {107900, 251700, 107900, 0, 0, 0},
        {107900, 107900, 251700, 0, 0, 0},
        {0, 0, 0, 71920, 0, 0},
        {0, 0, 0, 0, 71920, 0},
        {0, 0, 0, 0, 0, 71920}};
    ifstream slip_file("FCC_111.txt",ifstream::in);
    string temp_slip;
    while (!slip_file.eof()){
        getline(slip_file,temp_slip);
        Slip tmp(temp_slip);
        tmp.cal_shear_modulus(elastic_modulus);
        fcc_110_111.push_back(tmp);
    }
    slip_file.close();
    //End of [Slip system setting]

    //[Grain initialization]
    Grain testgrain;
    ofstream grain_info("grain_stepinfo.csv",ofstream::out);
    ofstream grain_str("grain_str.csv",ofstream::out);
    ofstream disloc_file("strain_disl.csv",ofstream::out);
    ofstream crss_file("slip_crss.csv",ofstream::out);
    grain_info << "e11,e22,e33,e12,e13,e23,s11,s22,s33,s12,s13,s23" << endl;
    grain_str << "e11,e22,e33,e12,e13,e23,s11,s22,s33,s12,s13,s23" << endl;
    testgrain.elastic_modulus = elastic_modulus;
    testgrain.lattice_vec = Matrix3d::Identity() * 3.88;
    testgrain.set_slip_sys(fcc_110_111);
    testgrain.print_stress_strain(grain_info);
    testgrain.print_stress_strain(grain_str);
    //End of [Grain initialization]

    //[Loading]
    ifstream load_file("Load.txt",ifstream::in);
    string step_conf_string;
    getline(load_file, step_conf_string);
    // args: 1: timestep, 2: substep, 3: max_strain
    step_config(step_conf_string);
    Matrix3d vel_grad_tensor = load_matrix_input(load_file), vel_grad_flag = load_matrix_input(load_file),
             stress_incr = load_matrix_input(load_file), dstress_flag = load_matrix_input(load_file);
    
    //Start step
    int step_ctrl = 0;
    while(testgrain.strain_tensor.cwiseAbs().maxCoeff() < max_strain){
        step_ctrl = 0;
        while(step_ctrl < 1/substep){
            if (step_ctrl == 0) testgrain.update_status(vel_grad_tensor * timestep, vel_grad_flag, stress_incr, dstress_flag);
            else testgrain.update_status(Matrix3d::Zero(), vel_grad_flag, Matrix3d::Zero(), dstress_flag);
            testgrain.print_stress_strain(grain_info);
            ++step_ctrl;
        }
        testgrain.print_stress_strain_screen(grain_str);
        testgrain.print_dislocation(disloc_file);
        testgrain.print_crss(crss_file);
    }

    grain_info.close();
    return 0;
}