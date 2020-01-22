/*
 * =====================================================================================
 *
 *       Filename:  distr.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  21/01/20 11:45:12
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Freyja Walberg (), gx231@stud.uni-heidelberg.de
 *   Organization:  Universität Heidelberg
 *
 * =====================================================================================
 */

#include"distr.h"

#include<math.h>

double
distr
	(double x,
	 void* param)
{
	double y = *(double*) param;
	double den = pow(x - 2, 4) + pow(sin(x - 3), 8);
	return y / den;
}

double
envelop
	(double x,
	 void* param)
{
	lin_envelop_params_t *params = param;
	double y = 0;
	for(size_t i = 0; i < params->segments; ++i)
	{
		if(x < params->x[i+1])
		{
			y = params->y[i+1] - params->y[i];
			y /= params->x[i+1] - params->x[i];
			y *= x - params->x[i];
			y += params->y[i];
			break;
		}
	}
	return y/2;
}


double
inv_envelop
	(double I,
	 void* param)
{
	lin_envelop_params_t *params = param;
	double x = 0;
	for(size_t i = 0; i < params->segments; ++i)
	{
		double segment_area = 0.5*(params->y[i] + params->y[i+1]) * (params->x[i+1] - params->x[i]);
		if(I > segment_area)
		{
			I -= segment_area;
			continue;
		}
		double alpha = params->y[i+1] - params->y[i];
		alpha /= params->x[i+1] - params->x[i];
		x = params->x[i];
		x -= params->y[i]/alpha;
		double root = params->y[i] / alpha;
		root *= root;
		root += 2*I/alpha;
		root = sqrt(root);
		x += root;
	}
	return x; 
}
