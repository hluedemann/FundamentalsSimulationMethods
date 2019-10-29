/*
 * =====================================================================================
 *
 *       Filename:  double_pendulum.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  29/10/19 01:23:34
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Freyja Walberg (), gx231@stud.uni-heidelberg.de
 *   Organization:  Universität Heidelberg
 *
 * =====================================================================================
 */

#include"double_pendulum.h"
#include<math.h>

int
pendulum_differential
    (double* state,
     double r,
     void* params,
     double r_tol,
     double a_tol,
     double* diff)
{
	struct pendulum_params *param = params;
	double cosine = cos(state[0] - state[1]);
	double cos_sq = cosine*cosine*param->m2;

	double tmp;

	cos_sq /= param->m2 + param->m1;

	diff[0] = state[2]/(param->l1*param->l1*(param->m1+param->m2));
	diff[0] *= 1 + cos_sq;
	diff[0] /= 1 - cos_sq;
	tmp = state[3]*cosine;
	tmp /= (param->l1*param->l2*(param->m1+param->m2));
	tmp /= 1 - cos_sq;
	diff[0] -= tmp;

	diff[1] = state[2]*param->m2*param->l2/param->l1*cosine;
	diff[1] /= param->m1 + param->m2;
	diff[1] = state[3] - diff[1];
	diff[1] /= param->m2*param->l2*param->l2;
	diff[1] /= 1 - cos_sq;


	tmp = sin(state[0] - state[1]);
	tmp *= param->m2*param->l1*param->l1;
	tmp *= diff[0] * diff[1];
	
	diff[2] = -sin(state[0])*param->g*param->l1;
	diff[2] *= param->m1 + param->m2;
	diff[2] -= tmp;

	diff[3] = -sin(state[1])*param->g*param->l2*param->m2;
	diff[3] += tmp;

	return 0;


}

int
pendulum_terminate
    (double* state,
     double r,
     void* params,
     double r_tol,
     double a_tol)
{
	return 0;
}
