/*
 * =====================================================================================
 *
 *       Filename:  rk_classic.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  20/08/19 19:59:36
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Freyja Walberg (), gx231@stud.uni-heidelberg.de
 *   Organization:  Universität Heidelberg
 *
 * =====================================================================================
 */


#include"rk_classic.h"
#include<stdlib.h>


integrator_t* 
rk_alloc
    (size_t n,
     double r_tol,
     double a_tol)
{
    integrator_t* itg = malloc(sizeof(integrator_t));
    itg->n = n;
    itg->r_tol = r_tol;
    itg->a_tol = a_tol;
    itg->Y = malloc(n*sizeof(double));
    rk_params_t* params = malloc(sizeof(rk_params_t));
    params->k1 = malloc(n*sizeof(double));
    params->k2 = malloc(n*sizeof(double));
    params->k3 = malloc(n*sizeof(double));
    params->k4 = malloc(n*sizeof(double));
    params->y_tmp = malloc(n*sizeof(double));
    itg->internal_params = params;
    itg->single_step_f = rk_step;
    return itg;
}

void
rk_free
    (integrator_t* itg)
{
    rk_params_t* params = itg->internal_params;
    free(params->k1);
    free(params->k2);
    free(params->k3);
    free(params->k4);
    free(params->y_tmp);
    free(params);
    free(itg->Y);
    free(itg);
}

int 
rk_step
    (size_t n,
     double r_tol,
     double a_tol,
     double* Y,
     double* t,
     void* params,
     double* step_size,
     diff_function_f f,
     jacobian_f jac,
     void* internal_params)
{
    double h = *step_size;
    rk_params_t* itg_params = internal_params;
    int rval = f(Y, *t, params, r_tol, a_tol, itg_params->k1);
    if(rval < 0) 
        return -1;

    for(size_t i = 0; i < n; ++i)
    {
        itg_params->y_tmp[i] = Y[i] + h/2*itg_params->k1[i];
    }
    rval = f(itg_params->y_tmp, *t + h/2, params, r_tol, a_tol, itg_params->k2);
    if(rval < 0) 
        return -1;

    for(size_t i = 0; i < n; ++i)
    {
        itg_params->y_tmp[i] = Y[i] + h/2*itg_params->k2[i];
    }
    rval = f(itg_params->y_tmp, *t + h/2, params, r_tol, a_tol, itg_params->k3);
    if(rval < 0) 
        return -1;

    for(size_t i = 0; i < n; ++i)
    {
        itg_params->y_tmp[i] = Y[i] + h*itg_params->k3[i];
    }
    rval = f(itg_params->y_tmp, *t + h, params, r_tol, a_tol, itg_params->k4);
    if(rval < 0) 
        return -1;

    for(size_t i = 0; i < n; ++i)
    {
        double k = itg_params->k1[i] + 2*itg_params->k2[i] + 2*itg_params->k3[i] + itg_params->k4[i];
        k /= 6;
        k *= h;
        Y[i] += k;
    }
    *t += h;
    return 0;
}




