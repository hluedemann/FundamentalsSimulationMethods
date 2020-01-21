/*
 * =====================================================================================
 *
 *       Filename:  monte_carlo.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  21/01/20 11:28:19
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Freyja Walberg (), gx231@stud.uni-heidelberg.de
 *   Organization:  Universität Heidelberg
 *
 * =====================================================================================
 */

#include"monte_carlo.h"
#include"random_gen.h"
size_t
monte_carlo
	(size_t n,
	 distribution_f dist,
	 void* param,
	 double dist_max,
	 double a,
	 double b,
	 double* samples)
{
	size_t draws = 0;
	for(size_t i = 0; i < n; ++i)
	{
		double x, y;
		while(1)
		{
			++draws;
			x = random_gen(a, b);
			y = random_gen(0, dist_max);
			if(y < dist(x, param))
			{
				break;
			}
		}
		samples[i] = x;
	}
	return draws;
}
