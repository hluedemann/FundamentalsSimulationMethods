/*
 * =====================================================================================
 *
 *       Filename:  n_grav.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  05/11/19 21:26:42
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Freyja Walberg (), gx231@stud.uni-heidelberg.de
 *   Organization:  Universität Heidelberg
 *
 * =====================================================================================
 */
#include<stdio.h>
#include<math.h>
#include<string.h>

#include"n_grav.h"


int
n_grav
    (size_t n,
     double* x,
     double* a,
     void* param)
{
    struct n_grav_params* params = param;
    
    if(n%params->dim)
    {
        fprintf(stderr, "dimensionality does not match!\n");
        return -1;
    }
    memset(a, 0, sizeof(double) * n);
    for(size_t i = 0; i < n/params->dim; ++i)
    {
        
        for(size_t j = 0; j < n/params->dim; ++j)
        {
            if(i == j)
                continue;
            double diff, tmp = 0;
            for(size_t k = 0; k < params->dim; ++k)
            {
                diff = x[params->dim*j + k] - x[params->dim*i + k];
                tmp += diff*diff;
            }
            tmp = pow(tmp, 1.5);
            for(size_t k = 0; k < params->dim; ++k)
            {
                a[params->dim*i + k] += (x[params->dim*j + k] - x[params->dim*i + k])/tmp;
            }
        }
    }
    return 0;
}
