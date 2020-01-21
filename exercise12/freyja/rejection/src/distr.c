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
