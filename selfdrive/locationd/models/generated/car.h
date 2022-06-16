#pragma once
#include "rednose/helpers/common_ekf.h"
extern "C" {
void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_err_fun(double *nom_x, double *delta_x, double *out_5754959651653696277);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_5191293492570975441);
void car_H_mod_fun(double *state, double *out_2119171891377142393);
void car_f_fun(double *state, double dt, double *out_6621260549472542149);
void car_F_fun(double *state, double dt, double *out_8310282395611039980);
void car_h_25(double *state, double *unused, double *out_4214223911070908245);
void car_H_25(double *state, double *unused, double *out_4592393826364716439);
void car_h_24(double *state, double *unused, double *out_7964186168728422677);
void car_H_24(double *state, double *unused, double *out_2419744227359216873);
void car_h_30(double *state, double *unused, double *out_4649887513822112941);
void car_H_30(double *state, double *unused, double *out_7110726784871965066);
void car_h_26(double *state, double *unused, double *out_6073984429677187759);
void car_H_26(double *state, double *unused, double *out_7896919796125517040);
void car_h_27(double *state, double *unused, double *out_8589790291557686536);
void car_H_27(double *state, double *unused, double *out_9112423217653643333);
void car_h_29(double *state, double *unused, double *out_8454989324001408112);
void car_H_29(double *state, double *unused, double *out_7620958129186357250);
void car_h_28(double *state, double *unused, double *out_5575389940908204572);
void car_H_28(double *state, double *unused, double *out_2538559112116826676);
void car_h_31(double *state, double *unused, double *out_8065406122667509068);
void car_H_31(double *state, double *unused, double *out_7270711693892165564);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}