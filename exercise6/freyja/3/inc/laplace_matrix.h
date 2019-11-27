#ifndef LAPLACE_MATRIX_H
#define LAPLACE_MATRIX_H

/*
 * =====================================================================================
 *
 *       Filename:  laplace_matrix.h
 *
 *    Description:  allocator for the finite difference laplace operator as matrix
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
alloc_laplace_matrix
	(unsigned long N);

void
free_laplace_matrix();


#endif /* MATRIX_H */
