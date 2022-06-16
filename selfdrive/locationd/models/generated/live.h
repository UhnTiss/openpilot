#pragma once
#include "rednose/helpers/common_ekf.h"
extern "C" {
void live_update_4(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_9(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_10(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_12(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_35(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_32(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_13(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_14(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_33(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_H(double *in_vec, double *out_8674808526469126384);
void live_err_fun(double *nom_x, double *delta_x, double *out_826382423215298085);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_7532358635846209890);
void live_H_mod_fun(double *state, double *out_2474155941233457423);
void live_f_fun(double *state, double dt, double *out_5893392617996069488);
void live_F_fun(double *state, double dt, double *out_4522217942656145191);
void live_h_4(double *state, double *unused, double *out_4900326243380033195);
void live_H_4(double *state, double *unused, double *out_1849829052632691181);
void live_h_9(double *state, double *unused, double *out_3074700386960502334);
void live_H_9(double *state, double *unused, double *out_1608639406003100536);
void live_h_10(double *state, double *unused, double *out_294220397151957808);
void live_H_10(double *state, double *unused, double *out_7113328125087702607);
void live_h_12(double *state, double *unused, double *out_1340863276117353804);
void live_H_12(double *state, double *unused, double *out_3169627355399270614);
void live_h_35(double *state, double *unused, double *out_3634549000190945085);
void live_H_35(double *state, double *unused, double *out_5915190387724284323);
void live_h_32(double *state, double *unused, double *out_5880204314189519423);
void live_H_32(double *state, double *unused, double *out_5602568692836839295);
void live_h_13(double *state, double *unused, double *out_5340999065095199952);
void live_H_13(double *state, double *unused, double *out_6684707275772634272);
void live_h_14(double *state, double *unused, double *out_3074700386960502334);
void live_H_14(double *state, double *unused, double *out_1608639406003100536);
void live_h_33(double *state, double *unused, double *out_4252459005722954991);
void live_H_33(double *state, double *unused, double *out_9065747392363141927);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}