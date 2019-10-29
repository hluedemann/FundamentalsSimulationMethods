#ifndef RK2_H
#define RK2_H

/*
 * =====================================================================================
 *
 *       Filename:  rk2.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  29/10/19 01:52:46
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
 * @brief Implementation of a error-corrector runge-kutta method for the integrator interface.
 */

#include"integrator.h"

/**
 * @brief allocates the runge-kutta integrator 
 *
 * @param n     Dimension of the ODE to be solved
 * @param r_tol relative tolerance, not directly used by the integrator
 * @param a_tol absolute tolerance, not directly used by the integrator
 *
 * @return pointer to the integrator struct, returns NULL on error
 */
integrator_t* 
rk2_alloc
    (size_t n,
     double r_tol,
     double a_tol);

/**
 * @brief Gives back the memory allocated by rk2_alloc() 
 *
 * @param itg pointer to the integrator struct
 */
void
rk2_free
    (integrator_t* itg);

int 
rk2_step
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

typedef struct rk2_params {
    double* k1;
    double* k2;
    double* y_tmp;
} rk2_params_t;

#endif /* RK2_H */
