#include "singleX.h"

int main(int argc, char **argv)
{
    Matrix3d strain_tensor = Matrix3d::Zero(), strain_plastic = Matrix3d::Zero(), strain_elastic =  Matrix3d::Zero(), stress_tensor = Matrix3d::Zero();
    bool flag_subprint = false, flag_test = false, flag_adaptive = false;
    //[Read cmd args]
    try {  
	TCLAP::CmdLine cmd("Command description message", ' ', "0.9");
	TCLAP::ValueArg<double> stepArg("s","step","Timestep.",false,0,"double");
	TCLAP::ValueArg<double> substepArg("b","substep","Substep in a timestep.",false, 0,"double");
	TCLAP::ValueArg<double> endArg("e","endstrain","Set the final strain of this simulation.",false, 0,"double");
	TCLAP::ValueArg<double> outputArg("o","outputstep","Control the output intervals.",false, 0.001,"double");
	TCLAP::MultiArg<double> eulerArg("E","euler","Euler Angle.",false, "double");
	TCLAP::ValueArg<string> loadStr("l","load","Set the Load file path", false, "", "string");
	TCLAP::ValueArg<string> sxStr("x","sx","Set the SX file path", false, "", "string");
	TCLAP::SwitchArg subprintSwitch("p","subprint","Print the substep informations", cmd, false);
	TCLAP::SwitchArg adaptiveSwitch("a","adaptive","Adaptive the loading step during calculation", cmd, false);
	TCLAP::SwitchArg testSwitch("t","test","Only test the input configrations", cmd, false);

	cmd.add(stepArg); cmd.add(substepArg); cmd.add(outputArg); cmd.add(endArg); cmd.add(eulerArg);
	cmd.add(loadStr); cmd.add(sxStr);

	cmd.parse( argc, argv );
	timestep = stepArg.getValue(), substep = substepArg.getValue(), outputstep = outputArg.getValue(), max_strain = endArg.getValue();
	loadfile_path = loadStr.getValue(); sxfile_path = sxStr.getValue();
	flag_subprint = subprintSwitch.getValue(); flag_test = testSwitch.getValue(); flag_adaptive = adaptiveSwitch.getValue();
	if ( eulerArg.getValue().size() == 3) euler_line_input = eulerArg.getValue();
	if ( flag_subprint ) cout << "Enable printing substep informations. " << endl;
	if ( flag_adaptive ) cout << "Enable adaptive loading step." << endl;

	} catch (TCLAP::ArgException &e)  // catch exceptions
	{ std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }
    //[Read Configuration File]
    cout << "##### SXCpp Configuration #####" << endl;
    set_config();
    //[Grain initialization]
    cout << "#### Grain Initialization #####" << endl;
    Grain testgrain = read_grain();

    cout << "##### Create Output Files #####" << endl;
    outfile_initialization(testgrain);

    //[Load Setting]
    bc_modi_matrix << Matrix3d::Identity(), Matrix3d::Zero(), Matrix3d::Zero(),Matrix3d::Zero(), 0.5*Matrix3d::Identity(), 0.5*Matrix3d::Identity(),Matrix3d::Zero(), 0.5*Matrix3d::Identity(), -0.5*Matrix3d::Identity();
    vel_to_dw_matrix << bc_modi_matrix;
    Matrix3d vel_grad_tensor, vel_grad_flag, stress_incr, dstress_flag;
    cout << "##### Create Output Files #####" << endl;
    read_load(vel_grad_tensor, vel_grad_flag, stress_incr, dstress_flag);
	
    cout << "Settings: " << endl << "Timestep: " << timestep << endl << "Substep: " << substep << endl \
	 << "dtime: " << dtime << endl << "outputstep: " << outputstep << endl << "max_strain: " << max_strain << endl;
    
    if (flag_test == false){
    //[Start Loading]
    	srand(time(0));
    	if (flag_adaptive == true) adaptive_step_load_sx(testgrain, vel_grad_tensor, vel_grad_flag, stress_incr, dstress_flag, flag_subprint);
    	else singleXloading(testgrain, vel_grad_tensor, vel_grad_flag, stress_incr, dstress_flag, flag_subprint);
    	outfile_close(); cout << "finish!" << endl;
    }
    return 0;
}

