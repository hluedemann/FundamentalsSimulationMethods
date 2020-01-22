/*
 * =====================================================================================
 *
 *       Filename:  random_gen.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  21/01/20 10:24:00
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Freyja Walberg (), gx231@stud.uni-heidelberg.de
 *   Organization:  Universität Heidelberg
 *
 * =====================================================================================
 */

#include"random_gen.h"
#include<limits.h>
//#ifndef __RDRND__
#include<stdlib.h>
#include<time.h>
static int init = 0;
//#endif

double
random_gen
	(double a,
	 double b)
{
	unsigned long long r;
	double rnd;
/*
#ifdef __RDRND__ 
	__builtin_ia32_rdrand64_step(&r);
	rnd = (double)r/( (double) ULLONG_MAX);
#else
*/
	if(init == 0)
	{
		srand(time(0));
		init = 1;
	}
	r = rand();
	rnd = (double) r / ( (double) RAND_MAX);
//#endif
	rnd *= b - a;
	rnd += a;
	return rnd;
	
	
	
}
