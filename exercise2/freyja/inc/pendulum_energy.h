#ifndef PENDULUM_ENERGY_H
#define PENDULUM_ENERGY_H

/*
 * =====================================================================================
 *
 *       Filename:  pendulum_energy.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  29/10/19 14:18:26
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Freyja Walberg (), gx231@stud.uni-heidelberg.de
 *   Organization:  Universität Heidelberg
 *
 * =====================================================================================
 */
#include"double_pendulum.h"

int
pendulum_energy
	(size_t n,
	 double* states,
	 struct pendulum_params* params,
	 double* energy);

#endif /* PENDULUM_ENERGY_H */
