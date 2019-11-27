#ifndef GRADIENT_MATRIX_H
#define GRADIENT_MATRIX_H

/*
 * =====================================================================================
 *
 *       Filename:  gradient_matrix.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  26/11/19 22:01:46
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Freyja Walberg (), gx231@stud.uni-heidelberg.de
 *   Organization:  Universität Heidelberg
 *
 * =====================================================================================
 */

extern unsigned long* ija;
extern double* sa;

void
alloc_gradient_matrix
	(unsigned long N);

void
free_gradient_matrix();


#endif /* MATRIX_H */
