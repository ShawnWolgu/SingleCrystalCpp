#ifndef SINGX
#define SINGX

#include <cstdlib>
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
using Eigen::Matrix3d, Eigen::Vector3d, Eigen::Matrix, Eigen::MatrixXd, Eigen::VectorXd, Eigen::all, Eigen::last;
typedef Matrix<double, 6, 1> Vector6d;
typedef Matrix<double, 6, 6> Matrix6d;

#define pi 3.14159265358979323846
#define k_boltzmann 1.380649e-23
#define eV_to_J 1.60217662e-19
#define MPa_to_Pa 1e6

// [global variables]
extern double timestep, substep, dtime, m, max_strain, temperature, outputstep, norm_time;
extern int flag_harden;
extern ofstream stress_file, disloc_file, crss_file, rss_file, stress_step_file, disloc_step_file, disloc_step_file, euler_file, custom_output_file, accstrain_file, schmidt_file, disvel_file,time_step_file;
extern string sxfile_path, loadfile_path, configure_path;
extern Matrix6d strain_modi_tensor;
extern Matrix<double,9,9> bc_modi_matrix, vel_to_dw_matrix;
extern vector<double> euler_line_input;

// [classes]
class Grain;
class PMode;
class Slip;
class Twin;
enum mode_type {slip, twin, undefined};

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
void adaptive_step_load_sx(Grain &grain, Matrix3d vel_grad_tensor, Matrix3d vel_grad_flag, Matrix3d stress_incr, Matrix3d dstress_flag, bool flag_subprint);

// [functions]
void flag_to_idx(Matrix<double, 15, 1> flag, vector<int> &known_idx, vector<int> &unknown_idx);
void params_convert_to_matrix(Matrix<double, 15, 1> &params, Vector6d &unknown_params, vector<int> &unknown_idx, Matrix3d &vel_grad_elas, Matrix3d &stress_incr);
void cut_precision(Matrix3d &mat, int prec);
int sign(double x);
int heaviside(double x);
int get_interaction_mode(Vector3d burgers_i, Vector3d plane_i, Vector3d burgers_j, Vector3d plane_j);
double cal_cosine(Vector3d vec_i, Vector3d vec_j);
double calc_relative_error(Vector6d &v1, Vector6d &v2);
double calc_relative_error(double x, double y);
double calc_equivalent_value(Matrix3d mat);
double calc_first_principal(Matrix3d mat);
double relative_std(VectorXd &vec);
Vector3d get_plane_norm(Vector3d &plane_norm_disp, Matrix3d &lattice_vec);
Vector3d Euler_trans(Matrix3d euler_matrix);
Vector6d tensor_trans_order(Matrix3d tensor);
Matrix3d tensor_trans_order(Vector6d tensor);
Matrix3d tensor_trans_order_9(Matrix<double,9,1> tensor);
Matrix3d Euler_trans(Vector3d euler_vector);
Matrix3d Rodrigues(Matrix3d spin_elas);
Matrix3d vel_bc_to_vel_grad(Matrix3d vel_bc_tensor);
Matrix3d tensor_rot_to_CryCoord(Matrix3d tensor, Matrix3d orientation);
Matrix3d tensor_rot_to_RefCoord(Matrix3d tensor, Matrix3d orientation);
Matrix6d rotate_6d_stiff_modu(Matrix6d modulus, Matrix3d rotate_matrix);
Matrix6d rotate_6d_compl_modu(Matrix6d modulus, Matrix3d rotate_matrix);
Matrix<double,9,1> tensor_trans_order_9(Matrix3d tensor);
Matrix<double,9,1> vel_to_dw(Matrix3d tensor);

// [dislocation velocity model]

// [class members]
class PMode {
public:
    Vector3d burgers_vec, plane_norm;
    int num = -1;
    mode_type type = undefined;
    bool flag_active;
    double SSD_density, crss, acc_strain, disl_vel, custom_var, rss = 0.0, t_wait = 0.0, t_run = 0.0, rho_init=0.0, rho_H = 0.0;
    double ref_strain_rate = 0.001, rate_sen = m, shear_rate, ddgamma_dtau, shear_modulus;
    /*
     * [velocity parameters] 
     *  1. MFP control coeffient, 2. reference frequency, 3. activation energy, 4. slip resistance, 5. energy exponent
     *  6. saturated speed, 7. drag coefficient
     *
     * [hardening parameters] 
     *  8. forest hardening coefficient
     *
     * [DD evolution parameters] 
     *  0. SSD_density, 9. nucleation coefficient, 10. multiplication coefficient, 11. drag stress D, 12. reference strain rate, 13. c/g 
     */
    /* update_params: 0: burgers, 1: mean_free_path, 2: disl_density_resist, 3: forest_stress */
    /* 
     * [Twin parameters] 
     * 0. tau_0, 1. tau_1, 2. h_0, 3. h_1, 4. twin_strain, 5. SRS, 6. low_threshold, 7. high_threshold*/
    vector<double> harden_params, update_params, latent_params;
    Matrix3d schmidt;
    PMode();
    PMode(int slip_num, Vector6d &slip_info, vector<double> &hardens, vector<double> &latents, Matrix3d lattice_vec, double f_active);
    void cal_shear_modulus(Matrix6d elastic_modulus);
    double cal_rss(Matrix3d stress_tensor);
    virtual void cal_strain(Grain &grain, Matrix3d stress_tensor) {};
    virtual void cal_ddgamma_dtau(Matrix3d stress_tensor) {};
    virtual void update_status(Grain &grain) {};
    virtual void update_ssd(Matrix3d dstrain, Matrix3d orientation) {};
    virtual void update_rho_hard(vector<PMode*> mode_sys) {};
    Matrix3d dL_tensor();
    Matrix3d dstrain_tensor();
    Matrix3d drotate_tensor();
    Matrix6d ddp_dsigma();
    Matrix6d dwp_dsigma();
};

class Slip : public PMode{
public:
    double rho_sat = 0.0, rho_mov = 0.0;
    Slip();
    Slip(int slip_num, Vector6d &slip_info, vector<double> &hardens, vector<double> &latents, Matrix3d lattice_vec, double f_active);
    void cal_strain(Grain &grain, Matrix3d stress_tensor) override;
    void cal_ddgamma_dtau(Matrix3d stress_tensor) override;
    void update_status(Grain &grain) override;
    void update_ssd(Matrix3d dstrain, Matrix3d orientation) override;
    void update_rho_hard(vector<PMode*> mode_sys) override;
private:
    void cal_strain_pow(Matrix3d stress_tensor);
    void cal_strain_disvel(Matrix3d stress_tensor);
    void cal_ddgamma_dtau_pow(Matrix3d stress_tensor);
    void cal_ddgamma_dtau_disvel(Matrix3d stress_tensor);
    void update_disvel(vector<PMode*> mode_sys, MatrixXd lat_hard_mat, double bv);
    void update_voce(vector<PMode*> mode_sys, MatrixXd lat_hard_mat);
    double disl_velocity(double rss);
    vector<double> disl_velocity_grad(double rss);
    //Unused functions
    void update_ssd_old(Matrix3d dstrain, Matrix3d orientation);
    void update_rho_hard_old(vector<PMode*> mode_sys);
    void update_lhparams(Matrix3d dstrain);
    void update_cross_slip(vector<PMode> &mode_sys, Matrix3d stress_tensor);
    void update_rho_mov(vector<PMode> &mode_sys);
    void update_surface_nuc(Matrix3d stress_tensor);
};

class Twin : public PMode{
public:
    double twin_frac = 0.0, equiv_frac = 0.0;
    Twin();
    Twin(int slip_num, Vector6d &slip_info, vector<double> &hardens, vector<double> &latents, Matrix3d lattice_vec, double f_active);
    void cal_strain(Grain &grain, Matrix3d stress_tensor) override;
    void cal_ddgamma_dtau(Matrix3d stress_tensor) override;
    void update_status(Grain &grain) override;
    void update_ssd(Matrix3d dstrain, Matrix3d orientation) override;
private:
    enum twin_status {inactive, growth, saturated};
    twin_status status = inactive;
};

class Grain{
public:
    Matrix3d lattice_vec, deform_grad, deform_grad_elas, deform_grad_plas, stress_tensor, strain_tensor, orientation, orient_ref;
    Matrix6d elastic_modulus, elastic_modulus_ref;
    MatrixXd lat_hard_mat;
    vector<PMode*> mode_sys;
    double strain_rate = 1e-3, twin_frac = 0.0;
    Grain();
    Grain(Matrix6d elastic_mod, Matrix3d lat_vecs, vector<PMode*> s, MatrixXd latent_matrix, Matrix3d orient_Mat);
    void update_status(Matrix3d L_dt_tensor, Matrix3d vel_grad_flag, Matrix3d stress_incr, Matrix3d dstress_flag);
    void print_stress_strain(ofstream &os);
    void print_stress_strain_screen();
    void print_dislocation(ofstream &os);
    void print_crss(ofstream &os);
    void print_rss(ofstream &os);
    void print_accstrain(ofstream &os);
    void print_euler(ofstream &os);
    void print_schmidt(ofstream &os);
    void print_disvel(ofstream &os);
    void print_time(ofstream &os);
    Matrix3d get_vel_grad_plas(Matrix3d stress_incr);
private:
    Vector6d solve_L_dsigma(Matrix3d &vel_grad_elas, Matrix3d &vel_grad_flag, Matrix3d &stress_incr, Matrix3d &dstress_flag);
    Vector6d calc_fx(Matrix3d &L_dt_tensor, Matrix3d &stress_incr);
    Matrix<double, 3, 6> dwp_by_dsigma;
    Matrix<double, 6, 3> Sigma_ik;
    Matrix6d ddp_by_dsigma, C_ij_pri;
    Matrix6d get_C_ij_pri(Vector6d &stress_6d);
    Matrix6d get_dp_grad();
    Matrix<double, 6, 15> calc_dfx(Matrix3d &L_dt_tensor, Matrix3d &stress_incr);
    Matrix<double, 3, 6> get_wp_grad();
    Matrix<double, 6, 3> get_Sigma_ik(Vector6d &stress_6d);
    void calc_slip_ddgamma_dtau(Matrix3d stress_3d);
    void solve_iteration(Matrix3d &L_dt_tensor, Matrix3d &vel_grad_flag, Matrix3d &stress_incr, Matrix3d &dstress_flag);
};

#endif

