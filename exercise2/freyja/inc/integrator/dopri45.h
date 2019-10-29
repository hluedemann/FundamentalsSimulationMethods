#ifndef DOPRI45_H
#define DOPRI45_H

/*
 * =====================================================================================
 *
 *       Filename:  dopri45.h
 *
 *    Description:  Adaptation of the 4(5)-order adaptive method by Dormand & Prince (1980)
 *
 *        Version:  1.0
 *        Created:  29/08/19 01:26:01
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Freyja Walberg (), gx231@stud.uni-heidelberg.de
 *   Organization:  Universität Heidelberg
 *
 * =====================================================================================
 */
/**
 * @file
 * @brief Adaptation of the 4(5)-order adaptive method by Dormand & Prince (1980)
 */

#include"integrator.h"

/**
 * @brief allocates the DoPri-integrator 
 *
 * @param n Dimension of the ODE to be solved
 * @param r_tol relative tolerance, not directly used by the integrator 
 * @param a_tol absolute tolerance, also used as minimum stepsize 
 * @param s_tol error tolerance for each step
 * @param target_err factor by which the error should be smaller than \p s_tol 
 * @param max_inc maximum increase factor of the stepsize, \p max_inc > 1
 * @param max_dec maximum decrease factor of the stepsize, \p max_dec < 1
 * @param max_step maximum step size for the integrator 
 * @param rel_err use relative error instead of absolute error for error estimate
 *
 * @return pointer to the integrator struct, returns NULL on error
 */
integrator_t*
dopri_alloc
    (size_t n,
     double r_tol,
     double a_tol,
     double s_tol,
     double target_err,
     double max_inc,
     double max_dec,
     double max_step,
     int rel_err);

/**
 * @brief Frees the memory allocated by dopri_alloc() 
 *
 * @param itg pointer to the integrator struct
 */
void
dopri_free
    (integrator_t* itg);


int
dopri_step
    (size_t n,
     double r_tol,
     double a_tol,
     double* Y,
     double* t,
     void* params,
     double* step_size,
     diff_function_f,
     jacobian_f jac,
     void* internal_params);

typedef struct dopri_params {
    double* k1;
    double* k2;
    double* k3;
    double* k4;
    double* k5;
    double* k6;
    double* k7;
    double t_last;
    double min_step;
    double max_step;
    double max_inc;
    double max_dec;
    double s_tol;
    double target_err;
    int rel_err;
}dopri_params_t;

#endif /* DOPRI45_H */
