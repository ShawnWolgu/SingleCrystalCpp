#include "singleX.h"

Vector6d tensor_trans_order(Matrix3d tensor){
    /** 
     * 3*3 Matrix --> 1*6 Vector: [00,11,22,12,02,01] 
     */
    return Vector6d{{tensor(0,0), tensor(1,1), tensor(2,2), tensor(1,2), tensor(0,2), tensor(0,1)}};
}

Matrix<double,9,1> tensor_trans_order_9(Matrix3d tensor){
    /** 
     * 3*3 Matrix --> 9*1 Vector: [00,11,22,12,02,01,21,20,10] 
     */
    return Matrix<double,9,1> {{tensor(0,0), tensor(1,1), tensor(2,2), tensor(1,2), tensor(0,2), tensor(0,1), tensor(2,1), tensor(2,0), tensor(1,0)}};
}

Matrix3d tensor_trans_order(Vector6d tensor){
    /** 
     * 1*6 Vector --> 3*3 Matrix: [[0,5,4],[5,1,3],[4,3,2]] 
     */
    return Matrix3d{{tensor(0), tensor(5), tensor(4)}, {tensor(5), tensor(1), tensor(3)}, {tensor(4), tensor(3), tensor(2)}};
}

Matrix3d tensor_trans_order_9(Matrix<double,9,1> tensor){
    /** 
     * 9*1 Vector --> 3*3 Matrix: [[0,5,4],[8,1,3],[7,6,2]] 
     */
    return Matrix3d{{tensor(0), tensor(5), tensor(4)}, {tensor(8), tensor(1), tensor(3)}, {tensor(7), tensor(6), tensor(2)}};
}

void params_convert_to_matrix(Matrix<double, 15, 1> &params, Matrix3d &vel_grad_elas, Matrix3d &stress_incr){
    Matrix<double,9,1> temp_vel_grad_elas = params(Eigen::seq(0,8));
    Vector6d temp_stress_incr = params(Eigen::seq(9,14));
    vel_grad_elas = tensor_trans_order_9(temp_vel_grad_elas);
    stress_incr = tensor_trans_order(temp_stress_incr);
}

Matrix3d calc_stress(Matrix3d strain_elastic, Matrix6d elastic_modulus){
    Vector6d eps_elastic, sig_vector;
    Matrix3d stress_tensor = Matrix3d::Identity();
    eps_elastic = tensor_trans_order(strain_elastic);
    sig_vector = elastic_modulus * eps_elastic;
    stress_tensor = tensor_trans_order(sig_vector);
    return stress_tensor;
}

Vector6d get_vec_only_ith(Vector6d &vector_base, int i){
    Vector6d vec_ith = Vector6d::Zero();
    vec_ith(i) = vector_base(i) + 1e-2;
    return vec_ith;
};

ostream &operator<<(ostream &os, const Slip &slip){
    os << "Slip system: [" << slip.burgers_vec.transpose() << "](" << slip.plane_norm_disp.transpose() << "), CRSS = " << slip.crss << endl
        << "Schmidt Matrix = " << endl << slip.schmidt;
    return os;
}
