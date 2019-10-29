/*
 * =====================================================================================
 *
 *       Filename:  dopri45.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  29/08/19 02:03:21
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Freyja Walberg (), gx231@stud.uni-heidelberg.de
 *   Organization:  Universität Heidelberg
 *
 * =====================================================================================
 */

#include"dopri45.h"

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>

#define MIN_S_TOL 1e-5


/* definition of the Butcher tableau for DoPri: */

static const double c[] = {0.0, 1.0/5.0, 3.0/10.0, 4.0/5.0, 8.0/9.0, 1.0, 1.0};

static const double a2[] = {1.0/5.0};
static const double a3[] = {3.0/40.0, 9.0/40.0};
static const double a4[] = {44.0/45.0, -56.0/15.0, 32.0/9.0};
static const double a5[] = {19372.0/6561.0, -25360.0/2187.0, 64448.0/6561.0, -212.0/729.0};
static const double a6[] = {9017.0/3168.0, -355.0/33.0, 46732.0/5247.0, 49.0/176.0, -5103.0/18656.0};
static const double a7[] = {35.0/384.0, 0, 500.0/1113.0, 125.0/192.0, -2187.0/6784.0, 11.0/84.0};

static const double b[] = {35.0/384.0, 0, 500.0/1113.0, 125.0/192.0, -2187.0/6784.0, 11.0/84.0, 0};
static const double b_[] = {5179.0/57600.0, 0, 7571.0/16695.0, 393.0/640.0, -92097.0/339200.0, 187.0/2100.0, 1.0/40.0};





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
     int rel_err)
{
    integrator_t* itg = malloc(sizeof(integrator_t));
    if(!itg) 
        return NULL;
    itg->n = n;
    itg->r_tol = r_tol;
    itg->a_tol = a_tol;
    itg->Y = malloc(n*sizeof(double));
    if(!itg->Y)
    {
        free(itg);
        return NULL;
    }
    itg->params = NULL;
    itg->f = NULL;
    itg->cond = NULL;
    itg->jac = NULL;
    itg->single_step_f = dopri_step;

    dopri_params_t* params = malloc(sizeof(dopri_params_t));
    if(!params)
    {
        free(itg->Y);
        free(itg);
        return NULL;
    }
    params->k1 = malloc(n*sizeof(double));
    params->k2 = malloc(n*sizeof(double));
    params->k3 = malloc(n*sizeof(double));
    params->k4 = malloc(n*sizeof(double));
    params->k5 = malloc(n*sizeof(double));
    params->k6 = malloc(n*sizeof(double));
    params->k7 = malloc(n*sizeof(double));
    
    if(!(params->k1 &&
         params->k2 &&
         params->k3 &&
         params->k4 &&
         params->k5 &&
         params->k6 &&
         params->k7))
    {
        free(params->k1);
        free(params->k2);
        free(params->k3);
        free(params->k4);
        free(params->k5);
        free(params->k6);
        free(params->k7);
        free(params);
        free(itg->Y);
        free(itg);
        return NULL;
    }
    params->t_last = NAN;
    params->min_step = a_tol;
    params->max_step = max_step;
    params->max_inc = max_inc;
    params->max_dec = max_dec;
    params->s_tol = s_tol;
    if(s_tol < MIN_S_TOL)
    {
        fprintf(stderr, "Warning: provided error tolerance, increased to safe minimum\n");
        params->s_tol = MIN_S_TOL;
    }
    params->target_err = target_err;
    params->rel_err = rel_err;
    itg->internal_params = params;
    return itg;
}





    
    




void
dopri_free
    (integrator_t* itg)
{
    dopri_params_t* params = itg->internal_params;
    free(params->k1);
    free(params->k2);
    free(params->k3);
    free(params->k4);
    free(params->k5);
    free(params->k6);
    free(params->k7);
    free(params);
    free(itg->Y);
    free(itg);
}


int
dopri_step
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
    dopri_params_t* dp_params = internal_params;
    if(*t != dp_params->t_last)  //calculate first k if the times do not match (first step of the integration
    {
        int ret = f(Y, *t, params, r_tol, a_tol, dp_params->k1);
        if(ret < 0)
            return ret;
        dp_params->t_last = *t;
    }

    while(1)
    {
        double y[n];
        int ret;
        //K2:
        memcpy(y, Y, n*sizeof(double));
        for(size_t i = 0; i < n; ++i)
        {
            y[i] += *step_size * (a2[0]*dp_params->k1[i]);
        }
        ret = f(y, *t + *step_size*c[1], params, r_tol, a_tol, dp_params->k2);
        if(ret < 0)
        {
            return ret;
        }

        //K3:
        memcpy(y, Y, n*sizeof(double));
        for(size_t i = 0; i < n; ++i)
        {
            y[i] += *step_size * (a3[0]*dp_params->k1[i]);
            y[i] += *step_size * (a3[1]*dp_params->k2[i]);
        }
        ret = f(y, *t + *step_size*c[2], params, r_tol, a_tol, dp_params->k3);
        if(ret < 0)
        {
            return ret;
        }

        //K4:
        memcpy(y, Y, n*sizeof(double));
        for(size_t i = 0; i < n; ++i)
        {
            y[i] += *step_size * (a4[0]*dp_params->k1[i]);
            y[i] += *step_size * (a4[1]*dp_params->k2[i]);
            y[i] += *step_size * (a4[2]*dp_params->k3[i]);
        }
        ret = f(y, *t + *step_size*c[3], params, r_tol, a_tol, dp_params->k4);
        if(ret < 0)
        {
            return ret;
        }

        //K5:
        memcpy(y, Y, n*sizeof(double));
        for(size_t i = 0; i < n; ++i)
        {
            y[i] += *step_size * (a5[0]*dp_params->k1[i]);
            y[i] += *step_size * (a5[1]*dp_params->k2[i]);
            y[i] += *step_size * (a5[2]*dp_params->k3[i]);
            y[i] += *step_size * (a5[3]*dp_params->k4[i]);
        }
        ret = f(y, *t + *step_size*c[4], params, r_tol, a_tol, dp_params->k5);
        if(ret < 0)
        {
            return ret;
        }

        //K6:
        memcpy(y, Y, n*sizeof(double));
        for(size_t i = 0; i < n; ++i)
        {
            y[i] += *step_size * (a6[0]*dp_params->k1[i]);
            y[i] += *step_size * (a6[1]*dp_params->k2[i]);
            y[i] += *step_size * (a6[2]*dp_params->k3[i]);
            y[i] += *step_size * (a6[3]*dp_params->k4[i]);
            y[i] += *step_size * (a6[4]*dp_params->k5[i]);
        }
        ret = f(y, *t + *step_size*c[5], params, r_tol, a_tol, dp_params->k6);
        if(ret < 0)
        {
            return ret;
        }

        //K7:
        memcpy(y, Y, n*sizeof(double));
        for(size_t i = 0; i < n; ++i)
        {
            y[i] += *step_size * (a7[0]*dp_params->k1[i]);
            y[i] += *step_size * (a7[1]*dp_params->k2[i]);
            y[i] += *step_size * (a7[2]*dp_params->k3[i]);
            y[i] += *step_size * (a7[3]*dp_params->k4[i]);
            y[i] += *step_size * (a7[4]*dp_params->k5[i]);
            y[i] += *step_size * (a7[5]*dp_params->k6[i]);
        }
        ret = f(y, *t + *step_size*c[6], params, r_tol, a_tol, dp_params->k7);
        if(ret < 0)
        {
            return ret;
        }

        double norm(size_t, double*);
        double err[n];

        memcpy(y, Y, n*sizeof(double));
        for(size_t i = 0; i < n; ++i)
        {
            double sum = b[0]*dp_params->k1[i];
            sum += b[1]*dp_params->k2[i];
            sum += b[2]*dp_params->k3[i];
            sum += b[3]*dp_params->k4[i];
            sum += b[4]*dp_params->k5[i];
            sum += b[5]*dp_params->k6[i];
            sum += b[6]*dp_params->k7[i];
            y[i] += *step_size *sum;
            
            
            double sum_ =  (b[0] - b_[0])*dp_params->k1[i];
            sum_ += (b[1] - b_[1])*dp_params->k2[i];
            sum_ += (b[2] - b_[2])*dp_params->k3[i];
            sum_ += (b[3] - b_[3])*dp_params->k4[i];
            sum_ += (b[4] - b_[4])*dp_params->k5[i];
            sum_ += (b[5] - b_[5])*dp_params->k6[i];
            sum_ += (b[6] - b_[6])*dp_params->k7[i];
            if(dp_params->rel_err)
            {
                err[i] = sum_ / sum;
                if(sum < a_tol)
                    err[i] = 0;
            }
            else
            {
                err[i] = *step_size * sum_;
            }

        }
        double e = norm(n, err);
        double delta = pow(dp_params->target_err*dp_params->s_tol/e, 1/5);
        if(delta < dp_params->max_dec)
        {
            *step_size *= dp_params->max_dec;
        }
        else if(delta > dp_params->max_inc)
        {
            *step_size += dp_params->max_inc;
        }
        else
        {
            *step_size *= delta;
        }
        if(*step_size > dp_params->max_step)
        {
            *step_size = dp_params->max_step;
        }
        if(*step_size < dp_params->min_step)
        {
            fprintf(stderr, "error: step size to small\n");
            return -1;
        }
        if(e < dp_params->s_tol)
        {
            memcpy(Y, y, n*sizeof(double));
            memcpy(dp_params->k1, dp_params->k7, n*sizeof(double));
            *t += *step_size;
            dp_params->t_last = *t;
            break;
        }

    }

    return 0;







}


double
norm
    (size_t n,
     double* x)
{
    double sum = 0;
    for(size_t i = 0; i < n; ++i)
    {
        sum += x[i]*x[i];
    }
    return sqrt(sum);
}





