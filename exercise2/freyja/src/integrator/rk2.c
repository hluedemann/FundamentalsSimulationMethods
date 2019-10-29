/*
 * =====================================================================================
 *
 *       Filename:  rk2.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  29/10/19 01:54:56
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Freyja Walberg (), gx231@stud.uni-heidelberg.de
 *   Organization:  Universität Heidelberg
 *
 * =====================================================================================
 */



#include"rk2.h"
#include<stdlib.h>


integrator_t* 
rk2_alloc
    (size_t n,
     double r_tol,
     double a_tol)
{
    integrator_t* itg = malloc(sizeof(integrator_t));
    itg->n = n;
    itg->r_tol = r_tol;
    itg->a_tol = a_tol;
    itg->Y = malloc(n*sizeof(double));
    rk2_params_t* params = malloc(sizeof(rk2_params_t));
    params->k1 = malloc(n*sizeof(double));
    params->k2 = malloc(n*sizeof(double));
    params->y_tmp = malloc(n*sizeof(double));
    itg->internal_params = params;
    itg->single_step_f = rk2_step;
    return itg;
}

void
rk2_free
    (integrator_t* itg)
{
    rk2_params_t* params = itg->internal_params;
    free(params->k1);
    free(params->k2);
    free(params->y_tmp);
    free(params);
    free(itg->Y);
    free(itg);
}

int 
rk2_step 
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
    rk2_params_t* itg_params = internal_params;
    int rval = f(Y, *t, params, r_tol, a_tol, itg_params->k1);
    if(rval < 0) 
        return -1;

    for(size_t i = 0; i < n; ++i)
    {
        itg_params->y_tmp[i] = Y[i] + h*itg_params->k1[i];
    }
    rval = f(itg_params->y_tmp, *t + h, params, r_tol, a_tol, itg_params->k2);
    if(rval < 0) 
        return -1;


    for(size_t i = 0; i < n; ++i)
    {
        double k = itg_params->k1[i] + itg_params->k2[i];
        k /= 2;
        k *= h;
        Y[i] += k;
    }
    *t += h;
    return 0;
}




