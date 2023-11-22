#include "singleX.h"
#include <vector>

PMode::PMode() = default;

Matrix3d PMode::dL_tensor(double twin_frac) {return schmidt * shear_rate * dtime;}

Matrix6d PMode::ddp_dsigma() {
    Vector6d symSchmidt_6d = strain_modi_tensor * tensor_trans_order((Matrix3d)(0.5*(schmidt+schmidt.transpose())));
    Matrix6d symSchmidt_66 = symSchmidt_6d * symSchmidt_6d.transpose();
    return symSchmidt_66 * ddgamma_dtau * dtime;
}

Matrix6d PMode::dwp_dsigma() {
    Vector6d symSchmidt_6d = strain_modi_tensor * tensor_trans_order((Matrix3d)(0.5*(schmidt+schmidt.transpose())));
    Vector6d asymSchmidt_6d = strain_modi_tensor * tensor_trans_order((Matrix3d)(0.5*(schmidt-schmidt.transpose())));
    Matrix6d symSchmidt_66 = asymSchmidt_6d * symSchmidt_6d.transpose();
    return symSchmidt_66 * ddgamma_dtau * dtime;
}

double PMode::cal_rss(Matrix3d stress_tensor){
    return (stress_tensor.cwiseProduct(schmidt)).sum();
}

void PMode::cal_shear_modulus(Matrix6d elastic_modulus){
    Matrix3d slip_rotation;
    Vector3d trav_direc = burgers_vec.cross(plane_norm);
    slip_rotation << (burgers_vec/burgers_vec.norm()), plane_norm, trav_direc / trav_direc.norm();
    shear_modulus = rotate_6d_stiff_modu(elastic_modulus, slip_rotation.transpose())(3,3);
    if (type == slip){
        cout << "Shear modulus of slip system " << num << " is " << shear_modulus << endl;
    }
    else if (type == twin){
        cout << "Shear modulus of twin system " << num << " is " << shear_modulus << endl;
    }
    else{
        cout << "Shear modulus of deformation system " << num << " is " << shear_modulus << endl;
    }
}

Matrix3d PMode::dstrain_tensor() {return 0.5 * (schmidt + schmidt.transpose()) * shear_rate * dtime;}

Matrix3d PMode::drotate_tensor() {return 0.5 * (schmidt - schmidt.transpose()) * shear_rate * dtime;}

