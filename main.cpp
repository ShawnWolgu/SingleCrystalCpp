#include "singleX.h"

int main(int argc, char **argv)
{
    Matrix3d strain_tensor = Matrix3d::Zero(), strain_plastic = Matrix3d::Zero(), strain_elastic =  Matrix3d::Zero(), stress_tensor = Matrix3d::Zero();

    //[Grain initialization]
    Grain testgrain = read_grain();
    outfile_initialization(testgrain);

    //[Load Setting]
    Matrix3d vel_grad_tensor, vel_grad_flag, stress_incr, dstress_flag;
    read_load(vel_grad_tensor, vel_grad_flag, stress_incr, dstress_flag);

    //[Start Loading]
    singleXloading(testgrain, vel_grad_tensor, vel_grad_flag, stress_incr, dstress_flag);
    
    outfile_close();
    cout << "finish!" << endl;
    return 0;
}
