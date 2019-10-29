#ifndef DOUBLE_PENDULUM_H
#define DOUBLE_PENDULUM_H

struct pendulum_params
{
	double m1;
	double m2;
	double l1;
	double l2;
	double g;
};

int
pendulum_differential
    (double* state,
     double r,
     void* params,
     double r_tol,
     double a_tol,
     double* diff);

int
pendulum_terminate
    (double* state,
     double r,
     void* params,
     double r_tol,
     double a_tol);


#endif /* DOUBLE_PENDULUM_H */
