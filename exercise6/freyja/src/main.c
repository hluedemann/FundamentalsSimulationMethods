/*
 * =====================================================================================
 *
 *       Filename:  main.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  24/11/19 17:37:47
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Freyja Walberg (), gx231@stud.uni-heidelberg.de
 *   Organization:  Universität Heidelberg
 *
 * =====================================================================================
 */

#include"tridiag.h"
#include"heat.h"
#include"gnuplot_i.h"

#include"stdio.h"
#include<stdlib.h>

#define N1 100
#define N2 1000

int
main()
{
	/* 
	double A[3*N-2] = {1,2,2,3,4,4,5};
	double x[N] = {6,7,8};
	double y[N];

	tridiag_multiply(N, A, x, y);
	for(size_t i = 0; i < N; ++i)
	{
		printf("x_%zu = %e, y_%zu = %e\n", i, x[i], i, y[i]);
	}
	printf("\n\n");
	
	tridiag_solve(N, A, y);
	for(size_t i = 0; i < N; ++i)
	{
		printf("x_%zu = %e, y_%zu = %e\n", i, x[i], i, y[i]);
	}
	printf("\n\n");
	printf("%f %f 0\n%f %f %f\n0 %f %f\n\n", A[0], A[1], A[2], A[3], A[4], A[5], A[6] );
	*/

	double* T1 = malloc(sizeof(double)*N1);
	double* T2 = malloc(sizeof(double)*N2);

	double res1 = heat_transfer(N1, 1, 1, 1, 1, 1, T1);
	double res2 = heat_transfer(N2, 1, 1, 1, 1, 1, T2);
	
	double* x1 = malloc(sizeof(double)*N1);
	double* x2 = malloc(sizeof(double)*N2);
	
	double h1 = 2.0/(N1-1);
	double h2 = 2.0/(N2-1);
	for(size_t i = 0; i < N1; ++i)
	{
		x1[i] = -1+i*h1;
	}
	for(size_t i = 0; i < N2; ++i)
	{
		x2[i] = -1+i*h2;
	}

	gnuplot_ctrl* g;
	g = gnuplot_init();
	gnuplot_cmd(g, "set terminal png");
	gnuplot_cmd(g, "set output \"heat_diffusion.png\"");

	gnuplot_setstyle(g, "lines");
	gnuplot_plot_xy(g, x1, T1, N1, "N = 100");
	gnuplot_plot_xy(g, x2, T2, N2, "N = 1000");
	int a;
	scanf(" %i", &a);

	gnuplot_close(g);
	printf("Residual for %i points: %e\nResidual for %i points: %e\n", N1, res1, N2, res2);
}
