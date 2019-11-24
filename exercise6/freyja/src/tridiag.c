/*
 * =====================================================================================
 *
 *       Filename:  tridiag.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  24/11/19 17:45:45
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Freyja Walberg (), gx231@stud.uni-heidelberg.de
 *   Organization:  Universität Heidelberg
 *
 * =====================================================================================
 */

#include<stdio.h>
#include<math.h>
#include"tridiag.h"


#define TRIDIAG_TOL 1e-322	//if numbers get too small A is almost singular

void
tridiag_multiply
	(size_t N,
	 double* A,
	 double* x,
	 double* y)
{
	y[0] = A[0]*x[0] + A[1]*x[1];
	for(size_t i = 1; i < N - 1; ++i)
	{
		y[i] = A[3*i-1]*x[i-1] + A[3*i]*x[i] + A[3*i+1]*x[i+1];
	}
	y[N-1] = A[3*N-4]*x[N-2] + A[3*N-3]*x[N-1];
	
}



int
tridiag_solve
	(size_t N,
	 double* A,
	 double* b)
{
	for(size_t i = 1; i < N; ++i)
	{
		if(fabs(A[3*(i-1)]) < TRIDIAG_TOL)
		{
			fprintf(stderr, "Warning: Near-Singular Matrix");
			return -1;			
		}
		A[3*i] -= A[3*i-2]*A[3*i-1]/A[3*(i-1)];
		b[i] -= b[i-1]*A[3*i-1]/A[3*(i-1)];
	}
	b[N-1] /= A[3*(N-1)];
	for(size_t i = 2; i <= N; ++i)
	{
		b[N-i] -= b[N-i+1]*A[3*(N-i)+1];
		b[N-i] /= A[3*(N-i)];
	}


	return 0;
}
