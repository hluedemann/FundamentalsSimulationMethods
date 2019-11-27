/*
 * =====================================================================================
 *
 *       Filename:  gradient_matrix.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  26/11/19 22:01:46
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Freyja Walberg (), gx231@stud.uni-heidelberg.de
 *   Organization:  Universität Heidelberg
 *
 * =====================================================================================
 */

#include"gradient_matrix.h"
#include<stdlib.h>

unsigned long* ija;
double* sa;

void
alloc_gradient_matrix
	(unsigned long N)
{
	unsigned long n = N*N*(7*N-6) + 2;
	ija = malloc(n*sizeof(unsigned long));
	sa  = malloc(n*sizeof(double));
	for(unsigned long i = 0; i < N*N*N + 1; ++i)
	{
		sa[i] = 6;
	}
	for(unsigned long i = N*N*N + 1; i < n; ++i)
	{
		sa[i] = -1;
	}

	long curr_index = N*N*N+2;
	for(unsigned long i = 0; i < N; ++i)
	{
		for(unsigned long j = 0; j < N; ++j)
		{
			for(unsigned long k = 0; k < N; ++k)
			{
				ija[N*N*i + N*j + k + 1] = curr_index;
				
				if(i != 0)
				{
					ija[curr_index] = N*N*(i-1) + N*j + k;
					++curr_index;
				}
				if(j != 0)
				{
					ija[curr_index] = N*N*i + N*(j-1) + k;
					++curr_index;
				}
				if(k != 0)
				{
					ija[curr_index] = N*N*i + N*j + k - 1;
					++curr_index;
				}

				if(k != N-1)
				{
					ija[curr_index] = N*N*i + N*j + k + 1;
				}
				if(j != N-1)
				{
					ija[curr_index] = N*N*i + N*(j+1) + k;
				}
				if(i != N-1)
				{
					ija[curr_index] = N*N*(i+1) + N*j + k;
				}

				
			}
		}
	}
	ija[N*N*N + 1] = curr_index;


}

void
free_gradient_matrix()
{
	free(ija);
	free(sa);
}
