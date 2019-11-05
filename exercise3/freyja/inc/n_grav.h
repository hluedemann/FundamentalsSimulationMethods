#ifndef N_GRAV_H
#define N_GRAV_H

/*
 * =====================================================================================
 *
 *       Filename:  n_grav.h
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

struct n_grav_params
{
    size_t dim;
    double* masses;
};

int
n_grav
    (size_t n,
     double* x,
     double* a,
     void* param);

#endif /* N_GRAV_H */
