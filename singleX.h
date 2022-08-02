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
#include <tclap/CmdLine.h>
using namespace std;
using Eigen::Matrix3d, Eigen::Vector3d, Eigen::Matrix, Eigen::MatrixXd, Eigen::all, Eigen::last;
typedef Matrix<double, 6, 1> Vector6d;
typedef Matrix<double, 6, 6> Matrix6d;

// [constants]
#define pi 3.14159265358979323846
#define k_boltzmann 1.380649e-23
#define eV_to_J 1.60217662e-19
#define MPa_to_Pa 1e6

// [global variables]
extern double timestep, substep, dtime, m, max_strain, temperature, outputstep;
extern int flag_harden;
extern ofstream stress_file, disloc_file, crss_file, stress_step_file, disloc_step_file, disloc_step_file, euler_file, custom_output_file, accstrain_file;
extern string sxfile_path, loadfile_path, configure_path;
extern Matrix<double,9,9> bc_modi_matrix;
// [classes]
class Grain;
class Slip;

// [file read]
void set_config();
Grain read_grain();
void read_load(Matrix3d &vel_grad_tensor, Matrix3d &vel_grad_flag, Matrix3d &stress_incr, Matrix3d &dstress_flag);

// [file write]
void outfile_initialization();
void outfile_initialization(Grain &grain);
void substep_output(Grain &grain);
void grain_output(Grain &grain);
void outfile_close();

// [load apply]
void singleXloading(Grain &grain, Matrix3d vel_grad_tensor, Matrix3d vel_grad_flag, Matrix3d stress_incr, Matrix3d dstress_flag, bool flag_subprint);

// [functions]
void flag_to_idx(Matrix<double, 15, 1> flag, vector<int> &known_idx, vector<int> &unknown_idx);
void params_convert_to_matrix(Matrix<double, 15, 1> &params, Vector6d &unknown_params, vector<int> &unknown_idx, Matrix3d &vel_grad_elas, Matrix3d &stress_incr);
int sign(double x);
Vector3d Euler_trans(Matrix3d euler_matrix);
Vector6d tensor_trans_order(Matrix3d tensor);
Vector6d get_vec_only_ith(Vector6d &vector_base, int i); 
Matrix3d tensor_trans_order(Vector6d tensor);
Matrix3d tensor_trans_order_9(Matrix<double,9,1> tensor);
Matrix3d calc_stress(Matrix3d strain_elastic, Matrix6d elastic_modulus);
Matrix3d Euler_trans(Vector3d euler_vector);
Matrix3d Rodrigues(Matrix3d spin_elas);
Matrix3d vel_bc_to_vel_grad(Matrix3d vel_bc_tensor);
Matrix6d rotate_6d_stiff_modu(Matrix6d modulus, Matrix3d rotate_matrix);
Matrix6d rotate_6d_compl_modu(Matrix6d modulus, Matrix3d rotate_matrix);
Matrix<double,9,1> tensor_trans_order_9(Matrix3d tensor);

// [dislocation velocity model]
double disl_velocity(double rss, vector<double> harden_params, vector<double> update_params);
vector<double> disl_velocity_grad(double rss, vector<double> harden_params, vector<double> update_params);

// [class members]
class Slip {
    public:
        Vector3d burgers_vec, plane_norm, plane_norm_disp;
        Matrix3d schmidt;
        vector<double> harden_params, update_params;
        double ref_strain_rate = 0.001, rate_sen = m, strain_rate_slip, ddgamma_dtau, shear_modulus, SSD_density, crss, acc_strain, disl_vel;
        const double debye_freq = 9.13e13;
        Slip();
        Slip(Vector6d &slip_info, vector<double> &hardens, Matrix3d lattice_vec);
        Matrix3d dL_tensor();
        Matrix3d dstrain_tensor();
        Matrix3d drotate_tensor();
        Matrix6d ddp_dsigma();
        Matrix6d dwp_dsigma();
        void cal_strain(Grain &grain, Matrix3d stress_tensor);
        void cal_ddgamma_dtau(Grain &grain, Matrix3d stress_tensor);
        void update_status(Grain &grain);
        void cal_shear_modulus(Matrix6d elastic_modulus);
    private:
        void cal_strain_pow(Matrix3d stress_tensor, Matrix3d deform_grad_elas);
        void cal_strain_ddhard(Matrix3d stress_tensor, Matrix3d deform_grad_elas, double strain_rate);
        void cal_strain_disvel(Matrix3d stress_tensor, Matrix3d deform_grad_elas);
        void cal_ddgamma_dtau_pow(Matrix3d stress_tensor, Matrix3d deform_grad_elas);
        void cal_ddgamma_dtau_ddhard(Matrix3d stress_tensor, Matrix3d deform_grad_elas, double strain_rate);
        void cal_ddgamma_dtau_disvel(Matrix3d stress_tensor, Matrix3d deform_grad_elas);
        void update_ddhard(Matrix3d deform_grad_elas, vector<Slip> &slip_sys);
        void update_disvel(Matrix3d deform_grad_elas, vector<Slip> &slip_sys);
        void update_voce(vector<Slip> &slip_sys);
        double cal_rss(Matrix3d stress_tensor, Matrix3d deform_grad_elas);
        void cal_strain_disvel_old(Matrix3d stress_tensor, Matrix3d deform_grad_elas);
        void update_disvel_old(Matrix3d deform_grad_elas, vector<Slip> &slip_sys);
};


class Grain{
    public:
        Matrix3d lattice_vec, deform_grad, deform_grad_elas, deform_grad_plas, stress_tensor, strain_tensor, orientation;
        Matrix6d elastic_modulus, elastic_modulus_ref;
        vector<Slip> slip_sys;
        double strain_rate = 1e-3;
        Grain();
        Grain(Matrix6d elastic_mod, Matrix3d lat_vecs, vector<Slip> s, Matrix3d orient_Mat);
        void update_status(Matrix3d L_dt_tensor, Matrix3d vel_grad_flag, Matrix3d stress_incr, Matrix3d dstress_flag);
        void print_stress_strain(ofstream &os);
        void print_stress_strain_screen();
        void print_dislocation(ofstream &os);
        void print_crss(ofstream &os);
        void print_accstrain(ofstream &os);
        void print_euler(ofstream &os);
        Matrix3d get_vel_grad_plas(Matrix3d stress_incr);
    private:
        Vector6d solve_L_dsigma(Matrix3d &vel_grad_elas, Matrix3d &vel_grad_flag, Matrix3d &stress_incr, Matrix3d &dstress_flag);
        Vector6d calc_fx(Matrix3d &L_dt_tensor, Matrix3d &stress_incr);
        Matrix<double, 3, 6> dwp_by_dsigma;
        Matrix<double, 6, 3> Sigma_ik;
        Matrix<double, 9, 9> L_compo_to_wd();
        Matrix6d ddp_by_dsigma, C_ij_pri, strain_modi_tensor;
        Matrix6d get_C_ij_pri(Vector6d &stress_6d);
        Matrix6d get_dp_grad();
        Matrix6d calc_dfx(Matrix3d &L_dt_tensor, Matrix3d &stress_incr,vector<int> unknown_idx);
        Matrix<double, 3, 6> get_wp_grad();
        Matrix<double, 6, 3> get_Sigma_ik(Vector6d &stress_6d);
	void calc_slip_ddgamma_dtau(Matrix3d stress_3d);
        void solve_Lsig_iteration(Matrix3d &L_dt_tensor, Matrix3d &vel_grad_flag, Matrix3d &stress_incr, Matrix3d &dstress_flag);
};

#endif

