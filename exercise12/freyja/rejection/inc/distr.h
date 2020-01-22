#ifndef DISTR_H
#define DISTR_H

/*
 * =====================================================================================
 *
 *       Filename:  distr.h
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
#include<stddef.h>
typedef struct lin_envelop_params
{
	size_t segments;
	double *x;
	double *y;
}lin_envelop_params_t;

double
distr
	(double x,
	 void* param);

double
envelop
	(double x,
	 void* param);
	
double
inv_envelop
	(double x,
	 void* param);

#endif /* DISTR_H */
