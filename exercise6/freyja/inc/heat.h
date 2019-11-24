#ifndef HEAT_H
#define HEAT_H

/*
 * =====================================================================================
 *
 *       Filename:  heat.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  24/11/19 18:47:43
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Freyja Walberg (), gx231@stud.uni-heidelberg.de
 *   Organization:  Universität Heidelberg
 *
 * =====================================================================================
 */

#include<stddef.h>

double
heat_transfer
	(size_t N,
	 double D,
	 double eps,
	 double L,
	 double T_l,
	 double T_r,
	 double* T);
	


#endif /* HEAT_H */
