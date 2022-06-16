#include "car.h"

namespace {
#define DIM 9
#define EDIM 9
#define MEDIM 9
typedef void (*Hfun)(double *, double *, double *);

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}
const static double MAHA_THRESH_25 = 3.8414588206941227;
const static double MAHA_THRESH_24 = 5.991464547107981;
const static double MAHA_THRESH_30 = 3.8414588206941227;
const static double MAHA_THRESH_26 = 3.8414588206941227;
const static double MAHA_THRESH_27 = 3.8414588206941227;
const static double MAHA_THRESH_29 = 3.8414588206941227;
const static double MAHA_THRESH_28 = 3.8414588206941227;
const static double MAHA_THRESH_31 = 3.8414588206941227;

/******************************************************************************
 *                      Code generated with SymPy 1.10.1                      *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_5754959651653696277) {
   out_5754959651653696277[0] = delta_x[0] + nom_x[0];
   out_5754959651653696277[1] = delta_x[1] + nom_x[1];
   out_5754959651653696277[2] = delta_x[2] + nom_x[2];
   out_5754959651653696277[3] = delta_x[3] + nom_x[3];
   out_5754959651653696277[4] = delta_x[4] + nom_x[4];
   out_5754959651653696277[5] = delta_x[5] + nom_x[5];
   out_5754959651653696277[6] = delta_x[6] + nom_x[6];
   out_5754959651653696277[7] = delta_x[7] + nom_x[7];
   out_5754959651653696277[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_5191293492570975441) {
   out_5191293492570975441[0] = -nom_x[0] + true_x[0];
   out_5191293492570975441[1] = -nom_x[1] + true_x[1];
   out_5191293492570975441[2] = -nom_x[2] + true_x[2];
   out_5191293492570975441[3] = -nom_x[3] + true_x[3];
   out_5191293492570975441[4] = -nom_x[4] + true_x[4];
   out_5191293492570975441[5] = -nom_x[5] + true_x[5];
   out_5191293492570975441[6] = -nom_x[6] + true_x[6];
   out_5191293492570975441[7] = -nom_x[7] + true_x[7];
   out_5191293492570975441[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_2119171891377142393) {
   out_2119171891377142393[0] = 1.0;
   out_2119171891377142393[1] = 0;
   out_2119171891377142393[2] = 0;
   out_2119171891377142393[3] = 0;
   out_2119171891377142393[4] = 0;
   out_2119171891377142393[5] = 0;
   out_2119171891377142393[6] = 0;
   out_2119171891377142393[7] = 0;
   out_2119171891377142393[8] = 0;
   out_2119171891377142393[9] = 0;
   out_2119171891377142393[10] = 1.0;
   out_2119171891377142393[11] = 0;
   out_2119171891377142393[12] = 0;
   out_2119171891377142393[13] = 0;
   out_2119171891377142393[14] = 0;
   out_2119171891377142393[15] = 0;
   out_2119171891377142393[16] = 0;
   out_2119171891377142393[17] = 0;
   out_2119171891377142393[18] = 0;
   out_2119171891377142393[19] = 0;
   out_2119171891377142393[20] = 1.0;
   out_2119171891377142393[21] = 0;
   out_2119171891377142393[22] = 0;
   out_2119171891377142393[23] = 0;
   out_2119171891377142393[24] = 0;
   out_2119171891377142393[25] = 0;
   out_2119171891377142393[26] = 0;
   out_2119171891377142393[27] = 0;
   out_2119171891377142393[28] = 0;
   out_2119171891377142393[29] = 0;
   out_2119171891377142393[30] = 1.0;
   out_2119171891377142393[31] = 0;
   out_2119171891377142393[32] = 0;
   out_2119171891377142393[33] = 0;
   out_2119171891377142393[34] = 0;
   out_2119171891377142393[35] = 0;
   out_2119171891377142393[36] = 0;
   out_2119171891377142393[37] = 0;
   out_2119171891377142393[38] = 0;
   out_2119171891377142393[39] = 0;
   out_2119171891377142393[40] = 1.0;
   out_2119171891377142393[41] = 0;
   out_2119171891377142393[42] = 0;
   out_2119171891377142393[43] = 0;
   out_2119171891377142393[44] = 0;
   out_2119171891377142393[45] = 0;
   out_2119171891377142393[46] = 0;
   out_2119171891377142393[47] = 0;
   out_2119171891377142393[48] = 0;
   out_2119171891377142393[49] = 0;
   out_2119171891377142393[50] = 1.0;
   out_2119171891377142393[51] = 0;
   out_2119171891377142393[52] = 0;
   out_2119171891377142393[53] = 0;
   out_2119171891377142393[54] = 0;
   out_2119171891377142393[55] = 0;
   out_2119171891377142393[56] = 0;
   out_2119171891377142393[57] = 0;
   out_2119171891377142393[58] = 0;
   out_2119171891377142393[59] = 0;
   out_2119171891377142393[60] = 1.0;
   out_2119171891377142393[61] = 0;
   out_2119171891377142393[62] = 0;
   out_2119171891377142393[63] = 0;
   out_2119171891377142393[64] = 0;
   out_2119171891377142393[65] = 0;
   out_2119171891377142393[66] = 0;
   out_2119171891377142393[67] = 0;
   out_2119171891377142393[68] = 0;
   out_2119171891377142393[69] = 0;
   out_2119171891377142393[70] = 1.0;
   out_2119171891377142393[71] = 0;
   out_2119171891377142393[72] = 0;
   out_2119171891377142393[73] = 0;
   out_2119171891377142393[74] = 0;
   out_2119171891377142393[75] = 0;
   out_2119171891377142393[76] = 0;
   out_2119171891377142393[77] = 0;
   out_2119171891377142393[78] = 0;
   out_2119171891377142393[79] = 0;
   out_2119171891377142393[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_6621260549472542149) {
   out_6621260549472542149[0] = state[0];
   out_6621260549472542149[1] = state[1];
   out_6621260549472542149[2] = state[2];
   out_6621260549472542149[3] = state[3];
   out_6621260549472542149[4] = state[4];
   out_6621260549472542149[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_6621260549472542149[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_6621260549472542149[7] = state[7];
   out_6621260549472542149[8] = state[8];
}
void F_fun(double *state, double dt, double *out_8310282395611039980) {
   out_8310282395611039980[0] = 1;
   out_8310282395611039980[1] = 0;
   out_8310282395611039980[2] = 0;
   out_8310282395611039980[3] = 0;
   out_8310282395611039980[4] = 0;
   out_8310282395611039980[5] = 0;
   out_8310282395611039980[6] = 0;
   out_8310282395611039980[7] = 0;
   out_8310282395611039980[8] = 0;
   out_8310282395611039980[9] = 0;
   out_8310282395611039980[10] = 1;
   out_8310282395611039980[11] = 0;
   out_8310282395611039980[12] = 0;
   out_8310282395611039980[13] = 0;
   out_8310282395611039980[14] = 0;
   out_8310282395611039980[15] = 0;
   out_8310282395611039980[16] = 0;
   out_8310282395611039980[17] = 0;
   out_8310282395611039980[18] = 0;
   out_8310282395611039980[19] = 0;
   out_8310282395611039980[20] = 1;
   out_8310282395611039980[21] = 0;
   out_8310282395611039980[22] = 0;
   out_8310282395611039980[23] = 0;
   out_8310282395611039980[24] = 0;
   out_8310282395611039980[25] = 0;
   out_8310282395611039980[26] = 0;
   out_8310282395611039980[27] = 0;
   out_8310282395611039980[28] = 0;
   out_8310282395611039980[29] = 0;
   out_8310282395611039980[30] = 1;
   out_8310282395611039980[31] = 0;
   out_8310282395611039980[32] = 0;
   out_8310282395611039980[33] = 0;
   out_8310282395611039980[34] = 0;
   out_8310282395611039980[35] = 0;
   out_8310282395611039980[36] = 0;
   out_8310282395611039980[37] = 0;
   out_8310282395611039980[38] = 0;
   out_8310282395611039980[39] = 0;
   out_8310282395611039980[40] = 1;
   out_8310282395611039980[41] = 0;
   out_8310282395611039980[42] = 0;
   out_8310282395611039980[43] = 0;
   out_8310282395611039980[44] = 0;
   out_8310282395611039980[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_8310282395611039980[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_8310282395611039980[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_8310282395611039980[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_8310282395611039980[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_8310282395611039980[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_8310282395611039980[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_8310282395611039980[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_8310282395611039980[53] = -9.8000000000000007*dt;
   out_8310282395611039980[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_8310282395611039980[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_8310282395611039980[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8310282395611039980[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8310282395611039980[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_8310282395611039980[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_8310282395611039980[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_8310282395611039980[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8310282395611039980[62] = 0;
   out_8310282395611039980[63] = 0;
   out_8310282395611039980[64] = 0;
   out_8310282395611039980[65] = 0;
   out_8310282395611039980[66] = 0;
   out_8310282395611039980[67] = 0;
   out_8310282395611039980[68] = 0;
   out_8310282395611039980[69] = 0;
   out_8310282395611039980[70] = 1;
   out_8310282395611039980[71] = 0;
   out_8310282395611039980[72] = 0;
   out_8310282395611039980[73] = 0;
   out_8310282395611039980[74] = 0;
   out_8310282395611039980[75] = 0;
   out_8310282395611039980[76] = 0;
   out_8310282395611039980[77] = 0;
   out_8310282395611039980[78] = 0;
   out_8310282395611039980[79] = 0;
   out_8310282395611039980[80] = 1;
}
void h_25(double *state, double *unused, double *out_4214223911070908245) {
   out_4214223911070908245[0] = state[6];
}
void H_25(double *state, double *unused, double *out_4592393826364716439) {
   out_4592393826364716439[0] = 0;
   out_4592393826364716439[1] = 0;
   out_4592393826364716439[2] = 0;
   out_4592393826364716439[3] = 0;
   out_4592393826364716439[4] = 0;
   out_4592393826364716439[5] = 0;
   out_4592393826364716439[6] = 1;
   out_4592393826364716439[7] = 0;
   out_4592393826364716439[8] = 0;
}
void h_24(double *state, double *unused, double *out_7964186168728422677) {
   out_7964186168728422677[0] = state[4];
   out_7964186168728422677[1] = state[5];
}
void H_24(double *state, double *unused, double *out_2419744227359216873) {
   out_2419744227359216873[0] = 0;
   out_2419744227359216873[1] = 0;
   out_2419744227359216873[2] = 0;
   out_2419744227359216873[3] = 0;
   out_2419744227359216873[4] = 1;
   out_2419744227359216873[5] = 0;
   out_2419744227359216873[6] = 0;
   out_2419744227359216873[7] = 0;
   out_2419744227359216873[8] = 0;
   out_2419744227359216873[9] = 0;
   out_2419744227359216873[10] = 0;
   out_2419744227359216873[11] = 0;
   out_2419744227359216873[12] = 0;
   out_2419744227359216873[13] = 0;
   out_2419744227359216873[14] = 1;
   out_2419744227359216873[15] = 0;
   out_2419744227359216873[16] = 0;
   out_2419744227359216873[17] = 0;
}
void h_30(double *state, double *unused, double *out_4649887513822112941) {
   out_4649887513822112941[0] = state[4];
}
void H_30(double *state, double *unused, double *out_7110726784871965066) {
   out_7110726784871965066[0] = 0;
   out_7110726784871965066[1] = 0;
   out_7110726784871965066[2] = 0;
   out_7110726784871965066[3] = 0;
   out_7110726784871965066[4] = 1;
   out_7110726784871965066[5] = 0;
   out_7110726784871965066[6] = 0;
   out_7110726784871965066[7] = 0;
   out_7110726784871965066[8] = 0;
}
void h_26(double *state, double *unused, double *out_6073984429677187759) {
   out_6073984429677187759[0] = state[7];
}
void H_26(double *state, double *unused, double *out_7896919796125517040) {
   out_7896919796125517040[0] = 0;
   out_7896919796125517040[1] = 0;
   out_7896919796125517040[2] = 0;
   out_7896919796125517040[3] = 0;
   out_7896919796125517040[4] = 0;
   out_7896919796125517040[5] = 0;
   out_7896919796125517040[6] = 0;
   out_7896919796125517040[7] = 1;
   out_7896919796125517040[8] = 0;
}
void h_27(double *state, double *unused, double *out_8589790291557686536) {
   out_8589790291557686536[0] = state[3];
}
void H_27(double *state, double *unused, double *out_9112423217653643333) {
   out_9112423217653643333[0] = 0;
   out_9112423217653643333[1] = 0;
   out_9112423217653643333[2] = 0;
   out_9112423217653643333[3] = 1;
   out_9112423217653643333[4] = 0;
   out_9112423217653643333[5] = 0;
   out_9112423217653643333[6] = 0;
   out_9112423217653643333[7] = 0;
   out_9112423217653643333[8] = 0;
}
void h_29(double *state, double *unused, double *out_8454989324001408112) {
   out_8454989324001408112[0] = state[1];
}
void H_29(double *state, double *unused, double *out_7620958129186357250) {
   out_7620958129186357250[0] = 0;
   out_7620958129186357250[1] = 1;
   out_7620958129186357250[2] = 0;
   out_7620958129186357250[3] = 0;
   out_7620958129186357250[4] = 0;
   out_7620958129186357250[5] = 0;
   out_7620958129186357250[6] = 0;
   out_7620958129186357250[7] = 0;
   out_7620958129186357250[8] = 0;
}
void h_28(double *state, double *unused, double *out_5575389940908204572) {
   out_5575389940908204572[0] = state[0];
}
void H_28(double *state, double *unused, double *out_2538559112116826676) {
   out_2538559112116826676[0] = 1;
   out_2538559112116826676[1] = 0;
   out_2538559112116826676[2] = 0;
   out_2538559112116826676[3] = 0;
   out_2538559112116826676[4] = 0;
   out_2538559112116826676[5] = 0;
   out_2538559112116826676[6] = 0;
   out_2538559112116826676[7] = 0;
   out_2538559112116826676[8] = 0;
}
void h_31(double *state, double *unused, double *out_8065406122667509068) {
   out_8065406122667509068[0] = state[8];
}
void H_31(double *state, double *unused, double *out_7270711693892165564) {
   out_7270711693892165564[0] = 0;
   out_7270711693892165564[1] = 0;
   out_7270711693892165564[2] = 0;
   out_7270711693892165564[3] = 0;
   out_7270711693892165564[4] = 0;
   out_7270711693892165564[5] = 0;
   out_7270711693892165564[6] = 0;
   out_7270711693892165564[7] = 0;
   out_7270711693892165564[8] = 1;
}
#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}




}
extern "C" {

void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
}
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<2, 3, 0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
}
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
}
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
}
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
}
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
}
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
}
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_31, H_31, NULL, in_z, in_R, in_ea, MAHA_THRESH_31);
}
void car_err_fun(double *nom_x, double *delta_x, double *out_5754959651653696277) {
  err_fun(nom_x, delta_x, out_5754959651653696277);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_5191293492570975441) {
  inv_err_fun(nom_x, true_x, out_5191293492570975441);
}
void car_H_mod_fun(double *state, double *out_2119171891377142393) {
  H_mod_fun(state, out_2119171891377142393);
}
void car_f_fun(double *state, double dt, double *out_6621260549472542149) {
  f_fun(state,  dt, out_6621260549472542149);
}
void car_F_fun(double *state, double dt, double *out_8310282395611039980) {
  F_fun(state,  dt, out_8310282395611039980);
}
void car_h_25(double *state, double *unused, double *out_4214223911070908245) {
  h_25(state, unused, out_4214223911070908245);
}
void car_H_25(double *state, double *unused, double *out_4592393826364716439) {
  H_25(state, unused, out_4592393826364716439);
}
void car_h_24(double *state, double *unused, double *out_7964186168728422677) {
  h_24(state, unused, out_7964186168728422677);
}
void car_H_24(double *state, double *unused, double *out_2419744227359216873) {
  H_24(state, unused, out_2419744227359216873);
}
void car_h_30(double *state, double *unused, double *out_4649887513822112941) {
  h_30(state, unused, out_4649887513822112941);
}
void car_H_30(double *state, double *unused, double *out_7110726784871965066) {
  H_30(state, unused, out_7110726784871965066);
}
void car_h_26(double *state, double *unused, double *out_6073984429677187759) {
  h_26(state, unused, out_6073984429677187759);
}
void car_H_26(double *state, double *unused, double *out_7896919796125517040) {
  H_26(state, unused, out_7896919796125517040);
}
void car_h_27(double *state, double *unused, double *out_8589790291557686536) {
  h_27(state, unused, out_8589790291557686536);
}
void car_H_27(double *state, double *unused, double *out_9112423217653643333) {
  H_27(state, unused, out_9112423217653643333);
}
void car_h_29(double *state, double *unused, double *out_8454989324001408112) {
  h_29(state, unused, out_8454989324001408112);
}
void car_H_29(double *state, double *unused, double *out_7620958129186357250) {
  H_29(state, unused, out_7620958129186357250);
}
void car_h_28(double *state, double *unused, double *out_5575389940908204572) {
  h_28(state, unused, out_5575389940908204572);
}
void car_H_28(double *state, double *unused, double *out_2538559112116826676) {
  H_28(state, unused, out_2538559112116826676);
}
void car_h_31(double *state, double *unused, double *out_8065406122667509068) {
  h_31(state, unused, out_8065406122667509068);
}
void car_H_31(double *state, double *unused, double *out_7270711693892165564) {
  H_31(state, unused, out_7270711693892165564);
}
void car_predict(double *in_x, double *in_P, double *in_Q, double dt) {
  predict(in_x, in_P, in_Q, dt);
}
void car_set_mass(double x) {
  set_mass(x);
}
void car_set_rotational_inertia(double x) {
  set_rotational_inertia(x);
}
void car_set_center_to_front(double x) {
  set_center_to_front(x);
}
void car_set_center_to_rear(double x) {
  set_center_to_rear(x);
}
void car_set_stiffness_front(double x) {
  set_stiffness_front(x);
}
void car_set_stiffness_rear(double x) {
  set_stiffness_rear(x);
}
}

const EKF car = {
  .name = "car",
  .kinds = { 25, 24, 30, 26, 27, 29, 28, 31 },
  .feature_kinds = {  },
  .f_fun = car_f_fun,
  .F_fun = car_F_fun,
  .err_fun = car_err_fun,
  .inv_err_fun = car_inv_err_fun,
  .H_mod_fun = car_H_mod_fun,
  .predict = car_predict,
  .hs = {
    { 25, car_h_25 },
    { 24, car_h_24 },
    { 30, car_h_30 },
    { 26, car_h_26 },
    { 27, car_h_27 },
    { 29, car_h_29 },
    { 28, car_h_28 },
    { 31, car_h_31 },
  },
  .Hs = {
    { 25, car_H_25 },
    { 24, car_H_24 },
    { 30, car_H_30 },
    { 26, car_H_26 },
    { 27, car_H_27 },
    { 29, car_H_29 },
    { 28, car_H_28 },
    { 31, car_H_31 },
  },
  .updates = {
    { 25, car_update_25 },
    { 24, car_update_24 },
    { 30, car_update_30 },
    { 26, car_update_26 },
    { 27, car_update_27 },
    { 29, car_update_29 },
    { 28, car_update_28 },
    { 31, car_update_31 },
  },
  .Hes = {
  },
  .sets = {
    { "mass", car_set_mass },
    { "rotational_inertia", car_set_rotational_inertia },
    { "center_to_front", car_set_center_to_front },
    { "center_to_rear", car_set_center_to_rear },
    { "stiffness_front", car_set_stiffness_front },
    { "stiffness_rear", car_set_stiffness_rear },
  },
  .extra_routines = {
  },
};

ekf_init(car);
