/*
 * =====================================================================================
 *
 *       Filename:  pendulum_energy.c
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

#include<math.h>
#include<stddef.h>

#include"pendulum_energy.h"
#include"double_pendulum.h"


int
pendulum_energy
	(size_t n,
	 double* states,
	 struct pendulum_params* params,
	 double* energy)
{
	double diff[4];
	double state[4];
	for(size_t i = 0; i < n; ++i)
	{
		state[0] = states[4*i];
		state[1] = states[4*i + 1];
		state[2] = states[4*i + 2];
		state[3] = states[4*i + 3];

		pendulum_differential(state, 0, params, 0, 0, diff);

		double tmp = params->l1*diff[0];
		tmp *= tmp/2;
		tmp += params->m1 + params->m2;
		energy[i] = tmp;
		

		tmp = params->l2*diff[1];
		tmp *= tmp;
		tmp += params->m2 /2;
		energy[i] += tmp;

		tmp = params->m2*params->l1*params->l2;
		tmp *= diff[0]*diff[1];
		tmp *= cos(state[0] - state[1]);
		energy[i] += tmp;

		tmp = params->l1*(1 - cos(state[0]));
		tmp *= params->m1 + params->m2;
		tmp += params->m2*params->l2*(1-cos(state[1]));
		tmp *= params->g;

		energy[i] += tmp;
		
		
	}
	return 0;
}
