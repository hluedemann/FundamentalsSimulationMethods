/*
 * =====================================================================================
 *
 *       Filename:  main.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  05/11/19 21:49:11
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Freyja Walberg (), gx231@stud.uni-heidelberg.de
 *   Organization:  Universität Heidelberg
 *
 * =====================================================================================
 */
#include<stdio.h>

#include"leapfrog.h"
#include"n_grav.h"
#include"gnuplot_i.h"


#define TIME 10
#define TIMESTEP 0.01

int
main()
{
	struct n_grav_params params;

	params.dim = 3;
	params.masses = malloc(sizeof(double)*2);
	params.masses[0] = 1;
	params.masses[1] = 1;

	size_t steps  = (size_t) TIME/TIMESTEP;
	double pos[6*(steps + 1)];
	double vel[6*(steps + 1)];
	memset(pos, 0, sizeof(double)*6);
	memset(vel, 0, sizeof(double)*6);
	pos[0] = -0.5;
	pos[3] =  0.5;
	vel[1] = -0.5;
	vel[4] =  0.5;

	leapfrog(6, steps, pos, vel, TIMESTEP, n_grav, &params);

	gnuplot_ctrl* h;
	h = gnuplot_init();

	double x1[steps + 1];
	double x2[steps + 1];
	double y1[steps + 1];
	double y2[steps + 1];
	for(size_t i = 0; i <= steps; ++i)
	{
		x1[i] = pos[6*i];
		y1[i] = pos[6*i + 1];
		x2[i] = pos[6*i + 3];
		y2[i] = pos[6*i + 4];
	}
	gnuplot_plot_xy(h, x1, y1, steps + 1, "sun 1");
	gnuplot_plot_xy(h, x2, y2, steps + 1, "sun 2");

	int tmp; 
	scanf(" %i", &tmp);
	gnuplot_close(h);

	

}
