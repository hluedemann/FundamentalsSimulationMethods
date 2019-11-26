#ifndef TRIDIAG_H
#define TRIDIAG_H

/*
 * =====================================================================================
 *
 *       Filename:  tridiag.h
 *
 *    Description:  Tridiagonal vector multiply and  LSE-solver
 *
 *        Version:  1.0
 *        Created:  24/11/19 17:28:33
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Freyja Walberg (), gx231@stud.uni-heidelberg.de
 *   Organization:  Universität Heidelberg
 *
 * =====================================================================================
 */
/* Tridiagonal Matrix: the nonzero elements for all lines lie sequential, array has length 3N-2*/

#include<stddef.h>

/* calculates y <- A*x
 * x and y need at least N elements, A at least 3N-2
 */
void
tridiag_multiply
	(size_t N,
	 double* A,
	 double* x,
	 double* y);
/* solves the LSE A*x = b
 * b needs at least N elements, A needs at least 3N-2 Elements
 * the content of A gets overwritten, b stores the result x
 */
int
tridiag_solve
	(size_t N,
	 double* A,
	 double* b);

/* performes a single step of the jacobi method for the LSE A*x = b
 * A and b are conserved, x is overwritten with the new solution
 */
int
tridiag_solve_jacobi
	(size_t N,
	 double* A,
	 double* b,
	 double* x,
	 double* x1);

#endif /* TRIDIAG_H */
