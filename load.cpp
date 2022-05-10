#include "singleX.h"

void singleXloading(Grain &grain, Matrix3d vel_grad_tensor, Matrix3d vel_grad_flag, Matrix3d stress_incr, Matrix3d dstress_flag){
    int step_ctrl = 0;
    while(grain.strain_tensor.cwiseAbs().maxCoeff() < max_strain){
        step_ctrl = 0;
        while(step_ctrl < 1/substep){
            if (step_ctrl == 0) grain.update_status(vel_grad_tensor * timestep, vel_grad_flag, stress_incr, dstress_flag);
            else grain.update_status(Matrix3d::Zero(), vel_grad_flag, Matrix3d::Zero(), dstress_flag);
            substep_output(grain);
            ++step_ctrl;
        }
        grain_output(grain);
    }
}