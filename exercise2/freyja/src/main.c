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
#include<stdlib.h>

#include<rk2.h>
#include<rk_classic.h>

#include"double_pendulum.h"
#include"pendulum_energy.h"

#include"gnuplot_i.h"

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


	double initial_state[] = {50.0/180.0*M_PI, -120.0/180.0*M_PI, 0, 0};
	double e0;
	pendulum_energy(1, initial_state, &params, &e0);


	integrator_t *itg = rk2_alloc(4, RTOL, ATOL);
	itg_initialize(itg, initial_state, 0, &params, DELTA_T, pendulum_differential, pendulum_terminate, NULL);
	double *y_rk2;
	double *t_rk2;
	itg_integrate(itg, &y_rk2, &t_rk2, NULL, STEPS);
	rk2_free(itg);

	double e_rk2[STEPS + 1];
	pendulum_energy(STEPS+1, y_rk2, &params, e_rk2);

	itg = rk_alloc(4, RTOL, ATOL);
	itg_initialize(itg, initial_state, 0, &params, DELTA_T, pendulum_differential, pendulum_terminate, NULL);
	double *y_rk4;
	double *t_rk4;
	itg_integrate(itg, &y_rk4, &t_rk4, NULL, STEPS);
	rk_free(itg);

	double e_rk4[STEPS + 1];
	pendulum_energy(STEPS+1, y_rk4, &params, e_rk4);



	for(size_t i = 0; i < STEPS + 1; ++i)
	{
		e_rk4[i] = 1 - e_rk4[i]/e0; 
		e_rk2[i] = 1 - e_rk2[i]/e0; 
	}

	gnuplot_ctrl* h;
	h = gnuplot_init();
	gnuplot_plot_xy(h, t_rk2, e_rk2, STEPS +1, "rk2");
	gnuplot_plot_xy(h, t_rk4, e_rk4, STEPS +1, "rk4");

	int tmp;
	scanf("  %i", &tmp);
	//sleep(5);
	gnuplot_close(h);

}
