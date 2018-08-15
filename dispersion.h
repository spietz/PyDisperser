#ifndef _dispersion_h
#define _dispersion_h

#include <iostream>
#include <vector>
#include <cassert>
#include <math.h>
#include <random>
/* For parsing member-function: 
   Using bind with pointers to members 
   http://www.boost.org/doc/libs/1_34_0/libs/bind/bind.html#with_functions*/
#include <boost/bind.hpp> 	
#include "num_int.h"

class Disperser{

private:

  std::default_random_engine generator;

  /* Standard normal distribution */
  std::normal_distribution<double> N;

  /* Standard normal distribution */
  std::uniform_real_distribution<double> U;

  /* Reynold's number */
  double Re;


public:

  /* for gaussian distributions */
  double u_f, v_m;
  double  displacement;

  /* flags */
  bool bouncing_particles, direction_likelyhood, bursting_process;
  bool gaussian_u, gaussian_v;

  double fall_vel;

  /* Lower and upper limit */
  double lower_y, upper_y;

  Disperser(double Re_in, double lower_y_in = 0.0, double upper_y_in = 1.0);

  int get_sign(double y);

  double get_f();

  double get_u(double y);

  double get_v(double y);

  double get_l(double y);

  void update(const int NP, double* X, double* Y, double tF, int max_step);

};

#endif
