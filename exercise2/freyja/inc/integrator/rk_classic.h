#ifndef RK_CLASSIC_H
#define RK_CLASSIC_H

/*
 * =====================================================================================
 *
 *       Filename:  rk_classic.h
 *
 *    Description:  Implementation of the classic fourth-order Runge-Kutta method for the integrator interface 
 *
 *        Version:  1.0
 *        Created:  20/08/19 18:33:43
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
 * @brief Implementation of the classic fourth-order Runge-Kutta method for the integrator interface.
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
rk_alloc
    (size_t n,
     double r_tol,
     double a_tol);

/**
 * @brief Gives back the memory allocated by rk_alloc() 
 *
 * @param itg pointer to the integrator struct
 */
void
rk_free
    (integrator_t* itg);

int 
rk_step
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

typedef struct rk_params {
    double* k1;
    double* k2;
    double* k3;
    double* k4;
    double* y_tmp;
} rk_params_t;

#endif /* RK_CLASSIC_H */
