/*
 * =====================================================================================
 *
 *       Filename:  main.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  29/10/19 01:44:23
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Freyja Walberg (), gx231@stud.uni-heidelberg.de
 *   Organization:  Universität Heidelberg
 *
 * =====================================================================================
 */

#include<math.h>
#include<stdio.h>

#include<rk2.h>
#include<rk_classic.h>

#include"double_pendulum.h"

#define RTOL 1e-5
#define ATOL 1e-20

#define DELTA_T 0.05
#define STEPS 2000
int
main()
{
	struct pendulum_params params;
	params.g = 1;
	params.m1 =	0.5;
	params.m2 =	1;
	params.l1 = 2;
	params.l2 = 1;

	integrator_t *itg = rk_alloc(4, RTOL, ATOL);

	double initial_state[] = {50.0/180.0*M_PI, -120.0/180.0*M_PI, 0, 0};
	itg_initialize(itg, initial_state, 0, &params, DELTA_T, pendulum_differential, pendulum_terminate, NULL);
	double *y;
	double *t;
	itg_integrate(itg, &y, &t, NULL, STEPS);

	for(size_t i = 0; i < STEPS; ++i)
	{
		printf("%e %e %e %e %e\n", t[i], y[4*i], y[4*i+1], y[4*i+2], y[4*i+3]);
	}

}
