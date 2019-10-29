/*
 * =====================================================================================
 *
 *       Filename:  integrator.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  20/08/19 19:50:04
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Freyja Walberg (), gx231@stud.uni-heidelberg.de
 *   Organization:  Universität Heidelberg
 *
 * =====================================================================================
 */


#include"integrator.h"

#include<string.h>
#include<stdlib.h>
#include<stdio.h>

#define BUFF 1024

int 
itg_initialize
    (integrator_t* itg,
     double* Y,
     double t,
     void* params,
     double step_size,
     diff_function_f f,
     termination_condition_f cond,
     jacobian_f jac)
{
    memcpy(itg->Y, Y, itg->n*sizeof(double));
    itg->t = t;
    itg->params = params;
    itg->step_size = step_size;
    itg->f = f;
    itg->jac = jac;
    itg->cond = cond;
    return 0;
}

int 
itg_single_step 
    (integrator_t* itg)
{
    return itg->single_step_f(
            itg->n,
            itg->r_tol,
            itg->a_tol,
            itg->Y,
            &(itg->t),
            itg->params,
            &(itg->step_size),
            itg->f,
            itg->jac,
            itg->internal_params);
}

int 
itg_integrate
    (integrator_t* itg,
     double** y_arr,
     double** t_arr,
     size_t* iter,
     size_t maxiter)
{
    if(itg->cond(itg->Y,
                 itg->t,
                 itg->params,
                 itg->r_tol,
                 itg->a_tol) == 1)
    {
        return 0;   //already done
    }
    size_t buff_size = BUFF;
    if(maxiter < buff_size) 
        buff_size = maxiter;
    double* Y;
    if(y_arr)
    {
        Y = malloc(buff_size*itg->n*sizeof(double));
        if(!Y) 
        {
            fprintf(stderr, "Dynamic allocation failed\n");
            return -1;
        }
    }
    double* t;
    if(t_arr)
    {
        t = malloc(buff_size*sizeof(double));
        if(!t) 
        {
            fprintf(stderr, "Dynamic allocation failed\n");
            return -1;
        }
    }
    
        
    
    size_t i;
    if(y_arr)
    {
        memcpy(Y, itg->Y, itg->n*sizeof(double));
    }
    if(t_arr)
    {
        t[0] = itg->t;
    }
    for(i = 1; i <= maxiter; ++i)
    {
        int rval = itg_single_step(itg);
        if(rval < 0) return -1;     //return on error
        
        if(i == buff_size)      //check if buffer is full
        {
            buff_size *= 2;
            if(maxiter < buff_size) 
                buff_size = maxiter + 1;
            if(y_arr)
            {
                Y = realloc(Y, buff_size*itg->n*sizeof(double));
                if(!Y) 
                {
                    fprintf(stderr, "Dynamic allocation failed\n");
                    return -1;
                }
            }
            if(t_arr)
            {
                t = realloc(t, buff_size*sizeof(double));
                if(!t) 
                {
                    fprintf(stderr, "Dynamic allocation failed\n");
                    return -1;
                }
            }
        }
        if(y_arr)
        {
            memcpy(Y+itg->n*i, itg->Y, itg->n*sizeof(double));
        }
        if(t_arr)
        {
            t[i] = itg->t;
        }
        int cond = itg->cond(itg->Y, itg->t, itg->params, itg->r_tol, itg->a_tol);
        if(cond < 0) 
            return -1;      //return on error
        if(cond) 
            break;          //if termination condition returns 1, terminate
    }
    if(y_arr)
    {
        Y = realloc(Y, (i+1)*itg->n*sizeof(double));
        *y_arr = Y;
    }
    if(t_arr)
    {
        t = realloc(t, (i+1)*sizeof(double));
        *t_arr = t;
    }
    if(iter)
    {
        *iter = i;
    }
    return 0;
}

