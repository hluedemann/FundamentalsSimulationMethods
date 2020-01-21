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

#include<stdio.h>
#include<math.h>

#define N 1000000
#define A 0.0
#define B 5.0

#define N_BINS 100

int
main()
{
	double p0 = 1/12.8101;
	double x[N];
	
	size_t samples = monte_carlo(N, distr, &p0, 2.4, A, B, x); 

	double rejection_rate = 1 - (double) N / ( (double) samples);

	size_t bins[N_BINS] = {};
	for(size_t i = 0; i < N; ++i)
	{
		double tmp = x[i] - A;
		tmp *= N_BINS;
		tmp /= B-A;
		size_t bin = (size_t) floor(tmp);
		bins[bin] += 1;
	}
	printf("rejection rate: %f\n\nBinning:\n", rejection_rate);
	int bin_width = 3;
	for(size_t i = 0; i < N_BINS; ++i)
	{
		printf("Bin %*zu: %zu\n", bin_width, i, bins[i]);
	}

}
