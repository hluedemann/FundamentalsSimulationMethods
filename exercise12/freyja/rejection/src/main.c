/*
 * =====================================================================================
 *
 *       Filename:  main.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  21/01/20 11:31:43
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Freyja Walberg (), gx231@stud.uni-heidelberg.de
 *   Organization:  Universität Heidelberg
 *
 * =====================================================================================
 */

#include"monte_carlo.h"
#include"distr.h"
#include"gnuplot_i.h"


#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define N 1000000
#define A 0.0
#define B 5.0

#define N_BINS 100
#define N_POINTS 1000
int
main()
{
	double p0 = 1/12.8101;
	double *x = calloc(sizeof(double), N);
	double *x_env = calloc(sizeof(double), N);;
	
	lin_envelop_params_t params;
	params.segments = 4;
	params.y = calloc(sizeof(double), 5);
	params.x = calloc(sizeof(double), 5);

	params.x[0] = 0;
	params.y[0] = 0.01;

	params.x[1] = 1.8;
	params.y[1] = 0.15;

	params.x[2] = 2.35;
	params.y[2] = 2.5;

	params.x[3] = 3;
	params.y[3] = 0.1;

	params.x[4] = 5;
	params.y[4] = 0.02;

	double I = 0;
	for(size_t i = 0; i < params.segments; ++i)
	{
		double segment = params.y[i] + params.y[i+1];
		segment *= params.x[i+1] - params.x[i];
		I += segment;
	}

	size_t samples = monte_carlo(N, distr, &p0, 2.4, A, B, x); 
	size_t samples_env = monte_carlo_enveloped(N, distr, &p0, envelop, inv_envelop, &params, I, x_env); 
	double rejection_rate = 1 - (double) N / ( (double) samples);
	double rejection_rate_env = 1 - (double) N / ( (double) samples_env);

	double *bins = calloc(sizeof(double), N_BINS);
	double *bins_env = calloc(sizeof(double), N_BINS);
	for(size_t i = 0; i < N; ++i)
	{
		double tmp = x[i] - A;
		tmp *= N_BINS;
		tmp /= B-A;
		size_t bin = (size_t) floor(tmp);
		bins[bin] += 1;
		double tmp_env = x_env[i] - A;
		tmp_env *= N_BINS;
		tmp_env /= B-A;
		size_t bin_env = (size_t) floor(tmp_env);
		bins_env[bin_env] += 1;
	}
	double *X = calloc(sizeof(double), N_BINS);
	X[0] = A;
	double dx = ((double) (B - A)) /((double) N_BINS);
	for(size_t i = 1; i < N_BINS; ++i)
	{
		X[i] = X[i-1] + dx;
	}
	printf("rejection rate: %f\nrejection rate envelop: %f\n", rejection_rate, rejection_rate_env);

	double *points_x = calloc(sizeof(double), N_POINTS);
	double *points_y = calloc(sizeof(double), N_POINTS);;
	points_x[0] = A;
	dx = ((double) (B - A)) /((double) N_POINTS - 1);
	for(size_t i = 1; i < N_POINTS; ++i)
	{
		points_x[i] = points_x[i-1] + dx;
		points_y[i] = distr(points_x[i], &p0);
	}

	for(size_t i = 0; i < N_BINS; ++i)
	{
		bins[i] /= (double) N;
		bins_env[i] /= (double) N;
	}

	gnuplot_ctrl *h;
	h = gnuplot_init();
	gnuplot_cmd(h, "set logscale y");
	gnuplot_plot_xy(h, points_x, points_y, N_POINTS, "target distribution");
	gnuplot_plot_xy(h, X, bins, N_BINS, "rejection sampling");
	gnuplot_plot_xy(h, X, bins_env, N_BINS, "enveloped rejection sampling");

	int tmp;
	scanf(" %i", &tmp);

	gnuplot_close(h);
}
