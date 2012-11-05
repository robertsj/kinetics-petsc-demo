#ifndef FUNCTIONS_HH_
#define FUNCTIONS_HH_

#include "constant_data.hh"
#include "data.hh"
#include "petsc.h"

/**
 *  @brief  Write a PETSc matrix to binary
 *  @param  M     PETSc matrix object
 *  @param  name  name of binary output file
 */
void write_matrix(Mat M, const char* name);

/**
 *  @brief  Write a PETSc matrix to binary
 *  @param  M     PETSc matrix object
 *  @param  name  name of binary output file
 */
void write_vector(Vec V, const char* name);

/**
 *  @brief  Update the material for a given time
 *  @param  data  reference to data object being updated
 *  @param  time  time at which to evaluate data
 */
void update_material(Data &data, double time);

/**
 *  @brief  Fill diffusion operators for steady state calculation
 *
 *  In particular, this updates the L and F operators stored in
 *  the Data object.
 *
 *  @param  data  reference to data object being updated
 *  @return       error code
 */
int fill_diffusion_operators(Data &data);

/**
 *  @brief  Fill time-dependent kinetics operator at a function of time
 *
 *  In particular, this updates the L and F operators stored in
 *  the Data object.
 *
 *  @param  A     time-dependent kinetics operator
 *  @param  data  reference to data object being updated
 *  @param  time  time at which to evaluate data
 *  @return       error code
 */
int fill_time_dependent_operator(Mat A, Data &data, double time);

/**
 *  @brief  Computes the right hand side operator for x_t = A(t)*x
 *
 *  For our linear problem, this routine updates the matrix A as
 *  for the given time.  For general nonlinear problems, the
 *  right hand side is an arbitrary function of both time and
 *  the solution vector.
 *
 *  @param  ts    PETSc time stepper context
 *  @param  time  time at which to evaluate data
 *  @param  x     solution vector
 *  @param  A     pointer to time dependent operator
 *  @param  B     pointer to preconditioner (which may be A)
 *  @param  flag  unused flag
 *  @param  ctx   user-defined context (in our case, a Data object)
 *  @return       error code
 */
int right_hand_side (TS ts, double time, Vec x, Mat *A, Mat *B,
                     MatStructure *flag, void *ctx);

/**
 *  @brief  Routine to monitor the solution at each time step
 *
 *  By defining one's own monitor function, it becomes possible
 *  to extract information at each time step that would otherwise
 *  be unavailable.  In this function, we compute the power and
 *  print it along with the time and step size.  Piped to file
 *  via standard output, this can be used for plotting in Python.
 *
 *  @param  ts    PETSc time stepper context
 *  @param  step  time step index
 *  @param  time  actual time at this step
 *  @param  x     solution vector
 *  @param  ctx   user-defined context (in our case, a Data object)
 *  @return       error code
 */
int monitor(TS ts, int step, double time, Vec x, void *ctx);

/**
 *  @brief  Compute the steady state flux and precursor concentration
 *  @param  x     solution vector
 *  @param  data  reference to problem data
 */
void compute_steady_state(Vec x, Data &data);

#endif // FUNCTIONS_HH_ 
