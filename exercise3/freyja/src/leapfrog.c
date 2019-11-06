/*
 * =====================================================================================
 *
 *       Filename:  leapfrog.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  05/11/19 19:57:15
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Freyja Walberg (), gx231@stud.uni-heidelberg.de
 *   Organization:  Universität Heidelberg
 *
 * =====================================================================================
 */

#include<stdlib.h>
#include"leapfrog.h"



int
leapfrog 
    (size_t n,
     size_t steps,
     double* position,
     double* velocity,
     double dt,
     system_function_f acceleration,
     void* param)
{
    double* a0 = malloc(sizeof(double)*n);
    double* a1 = malloc(sizeof(double)*n);
    int ret = acceleration(n, position, a0, param);
    if(ret < 0)
        return ret;
    for(size_t i = 0; i < steps; ++i)
    {
        for(size_t j = 0; j < n; ++j)
        {
            position[n*(i+1) + j] = position[n*i + j] + velocity[n*i + j]*dt + 0.5*dt*dt*a0[j];
        }
        ret = acceleration(n, position + n*(i+1), a1, param);
        if(ret < 0)
            return ret;
        for(size_t j = 0; j < n; ++j)
        {
            velocity[n*(i + 1) + j] = velocity[n*i + j] + 0.5*(a0[j] + a1[j])*dt;
        }
        double* tmp = a0;
        a0 = a1;
        a1 = tmp;
        
        
    }
    free(a0);
    free(a1);
    return 0;
}
