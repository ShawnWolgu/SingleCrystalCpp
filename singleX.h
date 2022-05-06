#ifndef SINGX
#define SINGX

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <time.h> 
#include <math.h>
using namespace std;
using Eigen::Matrix3d, Eigen::Vector3d, Eigen::Matrix, Eigen::MatrixXd, Eigen::all;
typedef Matrix<double, 6, 1> Vector6d;
typedef Matrix<double, 6, 6> Matrix6d;

extern double timestep, substep, dtime, m, max_strain, temperature;
extern int flag_harden;

Vector6d tensor_trans_order(Matrix3d tensor);
Matrix<double,9,1> tensor_trans_order_9(Matrix3d tensor);
Matrix3d tensor_trans_order(Vector6d tensor);
Matrix3d tensor_trans_order_9(Matrix<double,9,1> tensor);
Matrix3d update_stress(Matrix3d strain_elastic, Matrix6d elastic_modulus);
Vector6d get_vec_only_ith(Vector6d &vector_base, int i); 
void flag_to_idx(Matrix<double, 15, 1> flag, vector<int> &known_idx, vector<int> &unknown_idx);
void params_convert_to_matrix(Matrix<double, 15, 1> &params, Matrix3d &vel_grad_elas, Matrix3d &stress_incr);

class Grain;

class Slip {
    public:
        Vector3d burgers_vec, plane_norm, plane_norm_disp;
        Matrix3d schmidt;
        vector<double> harden_params;
        double ref_strain_rate = 0.0001, rate_sen = m, crss, acc_strain, strain_rate_slip, shear_modulus, SSD_density;
        const double k_boltzmann = 1.380649e-23, debye_freq = 9.13e13;
        Slip();
        Slip(Vector3d temp_bv, Vector3d temp_pn);
        Slip(istream &is);
        Slip(string &is);
        Matrix3d dL_tensor();
        Matrix3d dstrain_tensor();
        Matrix3d drotate_tensor();
        void cal_strain(Grain &grain);
        void cal_strain_pow(Matrix3d stress_tensor, Matrix3d deform_grad_elas);
        void cal_strain_disvel(Matrix3d stress_tensor, Matrix3d deform_grad_elas, vector<Slip> &slip_sys, double lattice_cons);
        void cal_shear_modulus(Matrix6d elastic_modulus);
        void update_status(Grain &grain);
        void update_disvel();
        void update_voce(vector<Slip> &slip_sys);
        double cal_rss(Matrix3d stress_tensor, Matrix3d deform_grad_elas);
};


class Grain{
    public:
        Matrix3d deform_grad, deform_grad_elas, deform_grad_plas, stress_tensor, strain_tensor, orientation;
        Matrix6d elastic_modulus;
        Matrix<double, 6, 6> dp_by_sigma;
        Matrix<double, 3, 6> wp_by_sigma;
        vector<Slip> *slip_sys;
        double lattice_cons;
        Grain();
        void set_slip_sys(vector<Slip> &s);
        void update_status(Matrix3d L_dt_tensor, Matrix3d vel_grad_flag, Matrix3d stress_incr, Matrix3d dstress_flag);
        void print_stress_strain(ofstream &os);
        void print_stress_strain_screen(ofstream &os);
        void print_dislocation(ofstream &os);
        void print_crss(ofstream &os);
        Matrix3d get_vel_grad_plas(Matrix3d stress_incr);
    private:
        Matrix6d CP_SigQ;
        Matrix<double, 6, 15> left_matrix();
        Matrix<double, 15, 15> mid_matrix();
        Matrix6d get_C_ij_pri(Vector6d &stress_6d);
        Matrix<double, 6, 3> get_Sigma_ik(Vector6d &stress_6d);
        void cal_matrix_dpwp_by_sigma(Matrix3d stress_incr);
        void solve_L_dsigma(Matrix3d &vel_grad_elas, Matrix3d &vel_grad_flag, Matrix3d &stress_incr, Matrix3d &dstress_flag);
};

#endif

