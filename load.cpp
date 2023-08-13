#include "singleX.h"

void singleXloading(Grain &grain, Matrix3d vel_grad_tensor, Matrix3d vel_grad_flag, Matrix3d stress_incr, Matrix3d dstress_flag, bool flag_subprint){
    norm_time = 0.0;
    double step_time = 0.0, dtime_save = 0.0;
    double op_control = outputstep;
    while(grain.strain_tensor.cwiseAbs().maxCoeff() < max_strain){
        step_time = 0.0;
        while(step_time < timestep){
	    dtime_save = dtime = min(dtime, timestep-step_time);
            if (step_time == 0.0) grain.update_status(vel_grad_tensor * timestep, vel_grad_flag, stress_incr, dstress_flag);
            else grain.update_status(Matrix3d::Zero(), vel_grad_flag, Matrix3d::Zero(), dstress_flag);
	    norm_time += dtime/timestep;
	    if(flag_subprint) substep_output(grain);
            step_time += dtime;
	    if (dtime == dtime_save) dtime = dtime * 1.1; 
        }
        if(grain.strain_tensor.cwiseAbs().maxCoeff() > op_control){
            op_control += outputstep;
	    grain_output(grain);
        }
    }
}

void adaptive_step_load_sx(Grain &grain, Matrix3d vel_grad_tensor, Matrix3d vel_grad_flag, Matrix3d stress_incr, Matrix3d dstress_flag, bool flag_subprint){
    norm_time = 0.0;
    double acc_time = 0.0, dtime_save = 0.0, op_control = outputstep, end_time = max_strain / grain.strain_tensor.cwiseAbs().maxCoeff();
    while(grain.strain_tensor.cwiseAbs().maxCoeff() < max_strain){
        acc_time = 0.0;
        dtime_save = dtime = min(dtime, end_time-acc_time);
        grain.update_status_adaptive(vel_grad_tensor, vel_grad_flag, stress_incr, dstress_flag);
        norm_time += dtime/timestep;
        if(flag_subprint) substep_output(grain);
        acc_time += dtime;
        if (dtime == dtime_save) dtime = min(dtime * 5, timestep); 
        //cout << dtime << endl;
        if(grain.strain_tensor.cwiseAbs().maxCoeff() > op_control){
            op_control += outputstep;
            grain_output(grain);
        }
    }
}
