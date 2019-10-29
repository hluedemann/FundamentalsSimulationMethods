#ifndef INTEGRATOR_H
#define INTEGRATOR_H



/*
 * =====================================================================================
 *
 *       Filename:  integrator.h
 *
 *    Description:  Interface for solvers of systems of first order ODE's  
 *
 *        Version:  1.0
 *        Created:  20/08/19 19:32:45
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Freyja Walberg (), gx231@stud.uni-heidelberg.de
 *   Organization:  Universität Heidelberg
 *
 * =====================================================================================
 */


/**
 * @file
 * @brief Interface for solvers of systems of first order ODE's 
 */
#include<stddef.h>


/**
 * @brief function pointer for the time derivative ODE, to be supplied by the user
 *
 * @param y state of the system
 * @param t current time of the system
 * @param params optional: additional parameters to be supplied
 * @param r_tol relative tolerance
 * @param a_tol absolute tolerance
 * @param y_diff current derivative with respect to timne to be stored here
 *
 * @return value smaller than 0 on error, 0 for success
 */
typedef int (*diff_function_f)
        (double* y,
         double t,
         void* params,
         double r_tol,
         double a_tol,
         double* y_diff);

/**
 * @brief ifunction pointer to the jacobian matrix for use in implicit solvers
 *
 * @param y state of the system
 * @param t current time of the system
 * @param params optional: additional parameters to be supplied
 * @param r_tol relative tolerance
 * @param a_tol absolute tolerance
 * @param jacobian jacobian matrix of the system function in row major
 *
 * @return value smaller than 0 for error, 0 on success
 */
typedef int (*jacobian_f)
        (double* y,
         double t,
         void* params,
         double r_tol,
         double a_tol,
         double* jacobian);

/**
 * @brief function pointer for the termination condition of the system, to be supplied by the user
 *
 * @param y state of the system
 * @param t current time of the system
 * @param params optional: additional parameters to be supplied
 * @param r_tol relative tolerance
 * @param a_tol absolute tolerance
 *
 * @return value smaller than 0 on error, value larger than 0 for normal termination, 0 otherwise
 */
typedef int (*termination_condition_f)
        (double* y,
         double t,
         void* params,
         double r_tol,
         double a_tol);





/**
 * @brief ODE integrator interface
 *
 * Interface for solvers of first order ODE solvers
 */
typedef struct integrator
{
    size_t n;       /**< size of the ODE system */
    double r_tol;   /**< relative tolerance used for the RHS of the ODE, changing might lead to undefined behaviour, reference the integrator information */
    double a_tol;   /**< absolute tolerance used for the RHS of the ODE, changing might lead to undefined behaviour, reference the integrator information */
    
    
    double* Y;      /**< current state of the integrator, at least n*sizeof(double) bytes are allocated */
    double  t;      /**< current time of the integrator */
    void* params;   
    double step_size;

    diff_function_f f;
    termination_condition_f cond;
    jacobian_f jac;     /* for use in some implicit solvers*/

    int (*single_step_f)
        (size_t n,
         double r_tol,
         double a_tol,
         double* Y,
         double* t,
         void* params,
         double* step_size,
         diff_function_f f,
         jacobian_f jac,
         void* internal_params);

    void* internal_params;  //internal parameters of the integrator, if a_tol and r_tol are used directly, writing them in here might be useful

}integrator_t;


/**
 * @brief Initialization function of the integrator, can be called multiple times without reallocation
 * 
 * @param itg integrator struct, allocated by specific solver
 * @param Y array of the initial state, values are copied into the integrator
 * @param t initial time of the integrator
 * @param params parameters of the ODE, value of the pointer is stored, content is not copied
 * @param step_size step-size parameter for the solver, exact behaviour depends on the used solver
 * @param f RHS of the ODE of type diff_function_f 
 * @param cond termination condition for the solver of type termination_condition_f, ignored if NULL
 *
 * @return -1 on error, 0 otherwise
 */
int 
itg_initialize
    (integrator_t* itg,
     double* Y,
     double t,
     void* params,
     double step_size,
     diff_function_f f,
     termination_condition_f cond,
     jacobian_f jac);

/**
 * @brief Advances the system by a single step
 *
 * @param itg integrator to bbe advanced
 *
 * @return 0 on sucsess, error code of the integrator otherwise
 */
int 
itg_single_step 
    (integrator_t* itg);

    /**
     * @brief integrates the system for \p maxiter steps, or untill the termination condition is met, whichever happens first
     *
     * @param itg integrator of the system
     * @param yarr pointer for the state of intermediate steps to be saved. Array is dynamically allocated, set to NULL to reduce overhead.
     * @param t_arr pointer for the time of intermediate steps to be saved. Array is dynamically allocated, set to NULL to reduce overhead.
     * @param iter stores the number of actual iteration steps done in here, set to NULL to ignore
     * @param maxiter maximum number of iteration steps
     *
     * @return -1 on error, 0 otherwise
     */
int 
itg_integrate
    (integrator_t* itg,
     double** yarr,
     double** t_arr,
     size_t* iter,
     size_t maxiter);





#endif /* INTEGRATOR_H */
