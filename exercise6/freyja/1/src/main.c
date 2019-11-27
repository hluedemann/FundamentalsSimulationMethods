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

#define N_J1 8
#define N_J2 100
#define ITER 30

#define BUF 128

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
	gnuplot_cmd(g, "set key outside");

	gnuplot_setstyle(g, "lines");
	gnuplot_plot_xy(g, x1, T1, N1, "N = 100");
	gnuplot_plot_xy(g, x2, T2, N2, "N = 1000");


	//gnuplot_cmd(g, "set terminal png");
	//gnuplot_cmd(g, "set terminal tikz createstyle");
	gnuplot_cmd(g, "set terminal tikz");
	//gnuplot_cmd(g, "set terminal epslatex color colortext");
	gnuplot_cmd(g, "set output \"latex/heat_diffusion.tex\"");

	gnuplot_cmd(g, "replot");

	gnuplot_close(g);
	printf("Residual for %i points: %e\nResidual for %i points: %e\n", N1, res1, N2, res2);

	double* Tj1 = malloc(sizeof(double)*(ITER+1)*N_J1+1);
	double* Tj2 = malloc(sizeof(double)*(ITER+1)*N_J2);


	double resj1 = heat_transfer_jacobi(N_J1, 1, 1, 1, 1, 1, Tj1, ITER);
	double resj2 = heat_transfer_jacobi(N_J2, 1, 1, 1, 1, 1, Tj2, ITER);
	printf("Jacobi method with %i iterations:\nResidual for %i points: %e\nResidual for %i points: %e\n",ITER, N_J1, resj1, N_J2, resj2); 
	double* xj1 = malloc(sizeof(double)*N_J1);
	double* xj2 = malloc(sizeof(double)*N_J2);
	
	double hj1 = 2.0/(N_J1-1);
	double hj2 = 2.0/(N_J2-1);
	for(size_t i = 0; i < N_J1; ++i)
	{
		xj1[i] = -1+i*hj1;
	}
	for(size_t i = 0; i < N_J2; ++i)
	{
		xj2[i] = -1+i*hj2;
	}

	gnuplot_ctrl *g1, *g2;
	g1 = gnuplot_init();
	g2 = gnuplot_init();
	gnuplot_cmd(g1, "set key outside");
	gnuplot_cmd(g2, "set key outside");
	

	char* buf = malloc(BUF);
	for(size_t i = 0; i <= ITER; ++i)
	{
		sprintf(buf, "Iteration %zu", i);
		gnuplot_plot_xy(g1, xj1, Tj1 + N_J1*i, N_J1, buf);
		gnuplot_plot_xy(g2, xj2, Tj2 + N_J2*i, N_J2, buf);
		
	}
	gnuplot_cmd(g1, "set terminal tikz");
	gnuplot_cmd(g1, "set output \"latex/heat_diffusion_jacobi_1.tex\"");
	gnuplot_cmd(g1, "replot");
	gnuplot_cmd(g2, "set terminal tikz");
	gnuplot_cmd(g2, "set output \"latex/heat_diffusion_jacobi_2.tex\"");
	gnuplot_cmd(g2, "replot");

	gnuplot_close(g1);
	gnuplot_close(g2);

}
