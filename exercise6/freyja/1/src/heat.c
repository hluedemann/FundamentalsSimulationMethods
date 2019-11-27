/*
 * =====================================================================================
 *
 *       Filename:  heat.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  24/11/19 18:47:56
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Freyja Walberg (), gx231@stud.uni-heidelberg.de
 *   Organization:  Universität Heidelberg
 *
 * =====================================================================================
 */

#include"heat.h"
#include"tridiag.h"
#include<stdlib.h>
#include<math.h>

double 
heat_transfer
	(size_t N,
	 double D,
	 double eps,
	 double L,
	 double T_l,
	 double T_r,
	 double* T)
{
	double h = 2*L/(N-1);
	double* A = malloc(sizeof(double)*(3*N-2));
	double* b = malloc(sizeof(double)*N);
	double* b_ = malloc(sizeof(double)*N);

	A[0] = 2;
	A[1] = -1;
	b[0] = h*h*eps/D + T_l;
	T[0] = h*h*eps/D + T_l;
	for(size_t i = 1; i < N-1; ++i)
	{
		A[3*i-1] = -1;
		A[3*i] = 2;
		A[3*i+1] = -1;
		b[i] = h*h*eps/D;
		T[i] = h*h*eps/D;
	}
	A[3*N-4] = -1;
	A[3*N -3] = 2;
	b[N-1] = h*h*eps/D + T_r;
	T[N-1] = h*h*eps/D + T_r;

	tridiag_solve(N, A, T);
	tridiag_multiply(N, A, T, b_);
	
	double res = 0;
	for(size_t i = 0; i < N; ++i)
	{
		double tmp = b[i] - b_[i];
		res += tmp*tmp;
	}
	return sqrt(res);
	
}



double
heat_transfer_jacobi
	(size_t N,
	 double D,
	 double eps,
	 double L,
	 double T_l,
	 double T_r,
	 double* T,
	 size_t iter)
{
	double h = 2*L/(N-1);
	double* A = malloc(sizeof(double)*(3*N-2));
	double* b = malloc(sizeof(double)*N);
	double* b_ = malloc(sizeof(double)*N);

	A[0] = 2;
	A[1] = -1;
	b[0] = h*h*eps/D + T_l;
	T[0] = h*h*eps/D + T_l;
	for(size_t i = 1; i < N-1; ++i)
	{
		A[3*i-1] = -1;
		A[3*i] = 2;
		A[3*i+1] = -1;
		b[i] = h*h*eps/D;
		T[i] = h*h*eps/D;
	}
	A[3*N-4] = -1;
	A[3*N -3] = 2;
	b[N-1] = h*h*eps/D + T_r;
	T[N-1] = h*h*eps/D + T_r;

	for(size_t i = 0; i < iter; ++i)
	{
		tridiag_solve_jacobi(N, A, b, T+N*i, T+N*(i+1));
	}

	
	tridiag_multiply(N, A, T, b_);
	
	double res = 0;
	for(size_t i = 0; i < N; ++i)
	{
		double tmp = b[i] - b_[i];
		res += tmp*tmp;
	}
	return sqrt(res);
	
}
