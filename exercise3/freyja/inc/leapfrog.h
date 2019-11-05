#ifndef LEAPFROG_H
#define LEAPFROG_H

/*
 * =====================================================================================
 *
 *       Filename:  leapfrog.h
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

typedef int (*system_function_f)
                (size_t n,
                 double* position,
                 double* acceleration,
                 void* param);

int
leapfrog 
    (size_t n,
     size_t steps,
     double* position,
     double* velocity,
     double dt,
     system_function_f acceleration,
     void* param);
     


#endif /* LEAPFROG_H */
