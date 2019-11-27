/*
 * =====================================================================================
 *
 *       Filename:  main.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  27/11/19 11:36:34
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Freyja Walberg (), gx231@stud.uni-heidelberg.de
 *   Organization:  Universität Heidelberg
 *
 * =====================================================================================
 */


#include"nr.h"
#include"gradient_matrix.h"
#include"gnuplot_i.h"

#define N 100
#define M 1
#define RHO 1.0

#define ITOL 2
#define TOL 1e-5

#define MAXITER 10
int
main()
{
	alloc_gradient_matrix(N);
	double* x = calloc(sizeof(double),N*N*N + 1);
	double* b = calloc(sizeof(double),N*N*N + 1);

	for(unsigned long i = 0; i < N; ++i)
	{
		for(unsigned long j = 0; j < N; ++j)
		{
			for(unsigned long k = 0; k < N; ++k)
			{
				if( (i >= N/2 -M) && (i < N/2 + M) &&
					(j >= N/2 -M) && (j < N/2 + M) &&
					(k >= N/2 -M) && (k < N/2 + M) )
				{
					b[N*N*i + N*j + k + 1] = -RHO;
				}
			}
		}
	}
	int iter;
	double err;
	linbcg(N*N*N, b, x, ITOL, TOL, MAXITER, &iter, &err);
	printf("finished with error %e after %i iterations\n", err, iter);
	double* cut = malloc(sizeof(double)*N);


	for(unsigned long i = 0; i < N; ++i)
	{
		cut[i] = x[i*(N*N+N+1) + 1];
	}

	
	gnuplot_ctrl* g;
	g = gnuplot_init();
	gnuplot_cmd(g, "set key outside");

	gnuplot_setstyle(g, "lines");
	gnuplot_plot_x(g, cut, N, "potential along diagonal");


	gnuplot_cmd(g, "set terminal tikz");
	gnuplot_cmd(g, "set output \"poisson.tex\"");
	

	gnuplot_cmd(g, "replot");

	gnuplot_close(g);




}



