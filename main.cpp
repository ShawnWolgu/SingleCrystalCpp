#include "singleX.h"


int main(int argc, char **argv)
{
    Matrix3d strain_tensor = Matrix3d::Zero(), strain_plastic = Matrix3d::Zero(), strain_elastic =  Matrix3d::Zero(), stress_tensor = Matrix3d::Zero();
    bool flag_subprint = false;
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
    //[Read cmd args]
    try {  
	TCLAP::CmdLine cmd("Command description message", ' ', "0.9");
	TCLAP::ValueArg<double> stepArg("s","step","Timestep.",false,timestep,"double");
	TCLAP::ValueArg<double> substepArg("b","substep","Substep in a timestep.",false, substep,"double");
	TCLAP::ValueArg<double> endArg("e","endstrain","Set the final strain of this simulation.",false, max_strain,"double");
	TCLAP::ValueArg<double> outputArg("o","outputstep","Control the output intervals.",false, 0.001,"double");

	cmd.add(stepArg);
	cmd.add(substepArg);
	cmd.add(outputArg);
	cmd.add(endArg);

	TCLAP::SwitchArg subprintSwitch("p","subprint","Print the substep informations", cmd, false);

	cmd.parse( argc, argv );

	timestep = stepArg.getValue();
	substep = substepArg.getValue();
	outputstep = outputArg.getValue();
	max_strain = endArg.getValue();
	dtime = timestep*substep;
	bool flag_subprint = subprintSwitch.getValue();

	if ( flag_subprint )
	{
		cout << "Enable printing substep informations. " << endl;
	}

	cout << "Setting via cmd args: " << endl << "Timestep: " << timestep << endl \
		<< "Substep: " << substep << endl << "outputstep: " << outputstep << endl \
		<< "max_strain: " << max_strain << endl;

	} catch (TCLAP::ArgException &e)  // catch exceptions
	{ std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }
    //[Start Loading]
    singleXloading(testgrain, vel_grad_tensor, vel_grad_flag, stress_incr, dstress_flag, flag_subprint);
    
    outfile_close();
    cout << "finish!" << endl;
    return 0;
}

