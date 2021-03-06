#include "singleX.h"

double timestep = 0.001, substep = 0.001, dtime = substep*timestep, m = 0.05, max_strain = 0.2, temperature = 298, outputstep = 0.001;
int flag_harden = 0;
string sxfile_path = "SingleX.txt", loadfile_path = "Load.txt", configure_path = "Config.txt";

ofstream stress_file("stress_grain.csv",ofstream::out), disloc_file("disloc_grain.csv",ofstream::out), crss_file("crss_grain.csv",ofstream::out), accstrain_file("accstrain.csv",ofstream::out);
ofstream stress_step_file("stress_step.csv",ofstream::out), disloc_step_file("disloc_step.csv",ofstream::out), euler_file("euler_angle_grain.csv",ofstream::out);
ofstream custom_output_file("custom_out.csv",ofstream::out);
