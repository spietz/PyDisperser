// -*- coding: utf-8 -*-


#include "dispersion.h"

#include<iostream>

Disperser::Disperser(double Re_in, double lower_y_in, double upper_y_in){

  /* default u_f, v_m for gaussian distributions */
  u_f = 0.0;
  v_m = 0.0;

  /* default "flags" */
  bouncing_particles = false;
  direction_likelyhood = false;
  bursting_process = false;
  gaussian_u = false;
  gaussian_v = false;
  displacement = 0.055;
  fall_vel = 0;

  Re = Re_in;
  lower_y = lower_y_in;
  upper_y = upper_y_in;

  /* Init generator with random seed */
  std::default_random_engine generator( std::random_device{}() );

  /* Init standard normal distribution */
  std::normal_distribution<double> N(0.0,1.0);

  /* Init standard uniform distribution */
  std::uniform_real_distribution<double> U(0.0,1.0);

};

int Disperser::get_sign(double y){

  /* 50% change of getting -1 1 independent of y*/
  if (!direction_likelyhood)
    return ((int)round(U(generator))*2-1);

  /* Instead of having 50 chance of +-1, we try to follow the physics by */
  /* saying that close to the bottom, the particle is more likely to go up. In */
  /* the midle 50 and in the top, it is more likely to go dowm. */
  if (U(generator) < (pow(2.0*(y-0.5-displacement),7)*0.5 + 0.5))
    return -1;
  return 1;

};

double Disperser::get_f(){
  /* Reynold's number functions */
  assert(Re > 0.0);

  return 2.5 * log(Re) + 5.0;
};

double Disperser::get_u(double y){
  /* Horizontal particle velocity */
  /* assert(y > lower_y && y < upper_y); */

  return 2.5 * log(y + exp(-2.0)/Re) + Disperser::get_f();
};

double Disperser::get_v(double y){
  /* Vertical particle velocity */
  double v_f = (-0.017 * y + 0.04) * Disperser::get_f();
  if (gaussian_v){
    double v = v_m + v_f*N(generator);
    if (std::abs(v) < 0.1*v_f)
      return std::copysign(0.1*v_f,v_f);
    return v;
  }
  return (double)get_sign(y)*v_f;
};

double Disperser::get_l(double y){
  /* Lenght scale of turbulence */
  if (y < 0.2)
    return (y + 0.001);
  return 0.2;
};

void Disperser::update(const int NP, double* X, double* Y, double tF, int max_step){

  double t, dt, t_burst_count, x, y, dy, v, T1 = 5.0/get_f();
  int step;

  for(int n=0; n<NP; ++n){  /* loop particles */

    x = X[n];
    y = Y[n];

    t = 0;      /* single particle time */
    dt = 0;
    t_burst_count = 0;
    step = 0;
    while((t < 0.999*tF) && (step < max_step)){ /* step time until final time */
      ++step;

      /* bursting process */
      if(y < 50./Re && bursting_process){
        t_burst_count += dt;
        if (t_burst_count > T1){
          /* std::cout << "particle " << n << " bursts!" << std::endl; */
          t_burst_count = 0;    /* reset counter */
          y = 100./Re;      /* move up */
        }
      }

      /* get vertical velocity */
      v = get_v(y);

      /* get time step size */
      dt = get_l(y)/std::abs(v);

      /* ensure incremented time doesn't exceed final*/
      if(t + dt > tF)
        dt = tF - t;

      /* increment time */
      t += dt;

      /* get y-displacement */
      dy = (v - fall_vel) * dt;

      /* no bouncing */
      if (!bouncing_particles){

        /* ensure position is within boundaries */
        if (y+dy < lower_y)
          dy = lower_y - y;
        else if (y+dy > upper_y)
          dy = upper_y - y;

        /* update x-coordinate depending on magnitude of y-displacement */
        if (std::abs(dy) < 0.01){ /* variation of u in y neglible */
          x += Disperser::get_u(y)*dt;
        }
        else{ /* variation of u in y significant -> use numerical integration */
          x += (1./dy *  integrate(boost::bind(&Disperser::get_u, this, _1),
                                   y, y+dy, 10, trapezium()) * dt);
        }


        /* update y-coordinate */
        y += dy;

      }

      /*  bouncing */
      else{

        if(y+dy < lower_y){     /* particle hit the wall */
          x += std::abs((1./dy *  integrate(boost::bind(&Disperser::get_u, this, _1),
                                            y, lower_y, 5, trapezium()) * dt));
          x += std::abs((1./dy *  integrate(boost::bind(&Disperser::get_u, this, _1),
                                            lower_y, -y-dy, 5, trapezium()) * dt));
          y = -y-dy;
        }
        else if(y+dy > upper_y){    /* particle hit free surface */
          x += std::abs((1./dy *  integrate(boost::bind(&Disperser::get_u, this, _1),
                                            y, upper_y, 5, trapezium()) * dt));
          x += std::abs((1./dy *  integrate(boost::bind(&Disperser::get_u, this, _1),
                                            upper_y, 2-y-dy, 5, trapezium()) * dt));
          y = 2-y-dy;
        }
        else{           /* particle did not hit either of them */
          x += (1./dy *  integrate(boost::bind(&Disperser::get_u, this, _1),
                                   y, y+dy, 10, trapezium()) * dt);
          y += dy;
        }

      }


    }

    if(gaussian_u)
      x += u_f*N(generator)*dt;

    X[n] = x;
    Y[n] = y;

  }


};
