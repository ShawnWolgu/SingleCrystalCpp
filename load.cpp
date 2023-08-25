#include "singleX.h"

void adaptive_step_load_sx(Grain &grain, Matrix3d vel_grad_tensor, Matrix3d vel_grad_flag, Matrix3d stress_incr, Matrix3d dstress_flag, bool flag_subprint){
    norm_time = 0.0;
    double dtime_save = 0.0, op_control = outputstep, time_left = timestep;
    while(grain.strain_tensor.cwiseAbs().maxCoeff() < max_strain){
        dtime_save = dtime = min(dtime, time_left);
        grain.update_status(vel_grad_tensor, vel_grad_flag, stress_incr, dstress_flag);
        time_left -= dtime;
        norm_time += dtime/timestep;
        if(flag_subprint) substep_output(grain);
        if(grain.strain_tensor.cwiseAbs().maxCoeff() >= op_control){
            op_control += outputstep;
            grain_output(grain);
        }
        if(dtime == dtime_save) dtime = min(dtime * 5, timestep); 
        if(time_left <= 1e-20) time_left = timestep;
    }
}
