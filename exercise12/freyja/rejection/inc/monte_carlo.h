#ifndef MONTE_CARLO_H
#define MONTE_CARLO_H

/*
 * =====================================================================================
 *
 *       Filename:  monte_carlo.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  21/01/20 10:59:02
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Freyja Walberg (), gx231@stud.uni-heidelberg.de
 *   Organization:  Universität Heidelberg
 *
 * =====================================================================================
 */
typedef double 
	  (*distribution_f) 
			(double,
			 void*);		

#include<stddef.h>
/**
 * @brief Generates a sample of points following the provided distribution in the intervall [a,b], generated with the rejection method
 *
 * @param n number of sampling points
 * @param dist the distribution of points
 * @param param parameters for the distribution 
 * @param dist_max maximum value the distribution takes on [a,b]
 * @param a lower limit of the intervall
 * @param b upper limit of the intervall
 * @param samples array of at least size n, samles will be stored in here
 *
 * @return number of total samples drawn
 */
size_t
monte_carlo
	(size_t n,
	 distribution_f dist,
	 void* param,
	 double dist_max,
	 double a,
	 double b,
	 double* samples);


/**
 * @brief Generates a sample of points following the provided distribution by using a known envelop 
 *
 * @param n number of sampling points
 * @param dist the distribution of points
 * @param param parameters for the distribution 
 * @param dist_max maximum value of the distribution on the intervall given by the envelop 
 * @param envelop envelop function
 * @param inv_envelop inverse of the integral of the envelop function
 * @param envelop_params parameters for the envelop and its inverse
 * @param I value of the integral over the envelop 
 * @param samples array of at least size n, samles will be stored in here
 *
 * @return number of total samples drawn
 */
size_t
monte_carlo_enveloped
	(size_t n,
	 distribution_f dist,
	 void* param,
	 distribution_f envelop,
	 distribution_f inv_envelop,
	 void* envelop_params,
	 double I,
	 double* samples);
	

#endif /* MONTE_CARLO_H */
