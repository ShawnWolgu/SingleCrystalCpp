#include "singleX.h"

Matrix6d cal_rotation_trans_6d_for_stiff(Matrix3d M);
Matrix6d cal_rotation_trans_6d_for_compl(Matrix3d M);

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
    Matrix6d strain_modi_tensor = Matrix6d::Zero();
    strain_modi_tensor.diagonal() << 1,1,1,2,2,2;
    eps_elastic = tensor_trans_order(strain_elastic);
    sig_vector = elastic_modulus * strain_modi_tensor * eps_elastic;
    stress_tensor = tensor_trans_order(sig_vector);
    return stress_tensor;
}

Vector6d get_vec_only_ith(Vector6d &vector_base, int i){
    Vector6d vec_ith = Vector6d::Zero();
    vec_ith(i) = vector_base(i) + 1e-2;
    return vec_ith;
}

void flag_to_idx(Matrix<double, 15, 1> flag, vector<int> &known_idx, vector<int> &unknown_idx){
    int unknown_cur = 0, known_cur = 0;
    for (int idx = 0;idx != 15; ++idx) {
        if(flag(idx) == 0) {unknown_idx[unknown_cur] = idx; ++unknown_cur;}
        else {known_idx[known_cur] = idx; ++known_cur;}
    }
    if (unknown_cur != 6 || known_cur != 9) throw "ERROR: wrong number of unknown parameters in L or dsigma!";
}

Matrix6d rotate_6d_stiff_modu(Matrix6d modulus, Matrix3d rotate_matrix)
{
    Matrix6d M66 = cal_rotation_trans_6d_for_stiff(rotate_matrix);
    return M66*modulus*M66.transpose();
}

Matrix6d cal_rotation_trans_6d_for_stiff(Matrix3d M)
{
	Matrix6d M66;
    double xx,yy,zz,yz,xz,xy,yx,zx,zy;
    xx = M(0,0); 
    yy = M(1,1);
    zz = M(2,2);
    yz = M(1,2);
    xz = M(0,2);
    xy = M(0,1);
	yx = M(1,0);
	zx = M(2,0);
	zy = M(2,1);

    for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
    M66(i,j) = M(i,j)*M(i,j);

    M66(0,3) = 2*xy*xz;
    M66(0,4) = 2*xx*xz;
    M66(0,5) = 2*xx*xy;

    M66(1,3) = 2*yy*yz;
    M66(1,4) = 2*yx*yz;
    M66(1,5) = 2*yx*yy;

    M66(2,3) = 2*zy*zz;
    M66(2,4) = 2*zx*zz;
    M66(2,5) = 2*zx*zy;

    M66(3,0) = yx*zx;
    M66(3,1) = yy*zy;
    M66(3,2) = yz*zz;

    M66(4,0) = xx*zx;
    M66(4,1) = xy*zy;
    M66(4,2) = xz*zz;

    M66(5,0) = xx*yx;
    M66(5,1) = xy*yy;
    M66(5,2) = xz*yz;

    M66(3,3) = yy*zz+yz*zy;
    M66(3,4) = yx*zz+yz*zx;
    M66(3,5) = yx*zy+yy*zx;

    M66(4,3) = xy*zz+xz*zy;
    M66(4,4) = xx*zz+xz*zx;
    M66(4,5) = xx*zy+xy*zx;

    M66(5,3) = xy*yz+xz*yy;
    M66(5,4) = xx*yz+xz*yx;
    M66(5,5) = xx*yy+xy*yx;

	return M66;
}

Matrix6d rotate_6d_compl_modu(Matrix6d modulus, Matrix3d rotate_matrix)
{
    Matrix6d M66 = cal_rotation_trans_6d_for_compl(rotate_matrix);
    return M66*modulus*M66.transpose();
}

Matrix6d cal_rotation_trans_6d_for_compl(Matrix3d M)
{
	Matrix6d M66;
    double xx,yy,zz,yz,xz,xy,yx,zx,zy;
    xx = M(0,0); 
    yy = M(1,1);
    zz = M(2,2);
    yz = M(1,2);
    xz = M(0,2);
    xy = M(0,1);
	yx = M(1,0);
	zx = M(2,0);
	zy = M(2,1);

    for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
    M66(i,j) = M(i,j)*M(i,j);

    M66(0,3) = xy*xz;
    M66(0,4) = xx*xz;
    M66(0,5) = xx*xy;

    M66(1,3) = yy*yz;
    M66(1,4) = yx*yz;
    M66(1,5) = yx*yy;

    M66(2,3) = zy*zz;
    M66(2,4) = zx*zz;
    M66(2,5) = zx*zy;

    M66(3,0) = 2*yx*zx;
    M66(3,1) = 2*yy*zy;
    M66(3,2) = 2*yz*zz;

    M66(4,0) = 2*xx*zx;
    M66(4,1) = 2*xy*zy;
    M66(4,2) = 2*xz*zz;

    M66(5,0) = 2*xx*yx;
    M66(5,1) = 2*xy*yy;
    M66(5,2) = 2*xz*yz;

    M66(3,3) = yy*zz+yz*zy;
    M66(3,4) = yx*zz+yz*zx;
    M66(3,5) = yx*zy+yy*zx;

    M66(4,3) = xy*zz+xz*zy;
    M66(4,4) = xx*zz+xz*zx;
    M66(4,5) = xx*zy+xy*zx;

    M66(5,3) = xy*yz+xz*yy;
    M66(5,4) = xx*yz+xz*yx;
    M66(5,5) = xx*yy+xy*yx;

	return M66;
}

Matrix3d Euler_trans(Vector3d euler_vector){
    Matrix3d Mout;
    double PH,TH,TM;
    euler_vector = euler_vector / 180 * pi; // in radians
    PH = euler_vector(0);
    TH = euler_vector(1);
    TM = euler_vector(2);

    double SPH,STH,STM;
    double CPH,CTH,CTM;
    SPH=sin(PH); CPH=cos(PH);
    STH=sin(TH); CTH=cos(TH);
    STM=sin(TM); CTM=cos(TM);
    Mout(0,0)=CTM*CPH-SPH*STM*CTH;
    Mout(1,0)=-STM*CPH-SPH*CTM*CTH;
    Mout(2,0)=SPH*STH;
    Mout(0,1)=CTM*SPH+CPH*STM*CTH;
    Mout(1,1)=-SPH*STM+CPH*CTM*CTH;
    Mout(2,1)=-STH*CPH;
    Mout(0,2)=STH*STM;
    Mout(1,2)=CTM*STH;
    Mout(2,2)=CTH;

    return Mout;
}

Vector3d Euler_trans(Matrix3d euler_matrix){
    Vector3d Vout;
    double PH,TH,TM;
    double STH;

    TH = acos(euler_matrix(2,2));
    if(abs(euler_matrix(2,2)) >= 0.99)
    {
        TM = 0;
        PH = atan2(euler_matrix(0,1),euler_matrix(0,0));
    }
    else
    {
        STH = sin(TH);
        TM = atan2(euler_matrix(0,2)/STH,euler_matrix(1,2)/STH);
        PH = atan2(euler_matrix(2,0)/STH,-euler_matrix(2,1)/STH); 
    }
    Vout << PH,TH,TM;
    Vout = Vout * 180 / pi; // in degree
    return Vout;
}

Matrix3d Rodrigues(Matrix3d spin_elas){
    Matrix3d exp_skew;
    Vector3d omegas;
    omegas << spin_elas(2,1), spin_elas(0,2), spin_elas(1,0);
    double omega_norm = omegas.norm();
    if(omega_norm < 1e-30)
        exp_skew = Matrix3d::Identity();
    else
    {
        Matrix3d term_1 = spin_elas/omega_norm * sin(omega_norm);
        Matrix3d term_2 = (spin_elas/omega_norm) * (spin_elas/omega_norm) * (1.0 - cos(omega_norm));
        exp_skew =  Matrix3d::Identity()+ term_1 + term_2;
    }
    return exp_skew;
}

ostream &operator<<(ostream &os, const Slip &slip){
    os << "Slip system: [" << slip.burgers_vec.transpose() << "](" << slip.plane_norm_disp.transpose() << "), CRSS = " << slip.crss << endl
        << "Schmidt Matrix = " << endl << slip.schmidt;
    return os;
}

int sign(double x){
    return (0 < x) - (x < 0);
}

void stress_incr_modify(Matrix3d stress_tensor, Matrix3d &stress_incr, Matrix3d stress_iter_save){
    Vector6d incr_6d = tensor_trans_order(stress_incr), stress_6d = tensor_trans_order(stress_tensor), incr_0_6d = tensor_trans_order(stress_iter_save);
    Vector6d direction = incr_6d - incr_0_6d;
    double step = direction.norm();
    direction = direction / step;
    double ratio = incr_6d.dot(-stress_6d / stress_6d.norm())/stress_6d.norm();
    if(ratio > 0.1){
        incr_6d = step*0.1/ratio * direction + incr_0_6d;
        stress_incr = tensor_trans_order(incr_6d);
    }
}
