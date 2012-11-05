// Demo of PETSc TS for two group slab reactor kinetics problem

#include "constant_data.hh"
#include "data.hh"
#include "functions.hh"
#include <iostream>

using std::cout;
using std::endl;

/**
 * @brief Function F(u, u_t, t) = 0
 * @param t     time at step/stage being solved
 * @param u     state vector
 * @param u_t   time derivative of state vector
 * @param F     function vector
 * @param ctx   user-defined context for matrix evaluation routine
 */
int DAE_F(TS ts, double time, Vec u, Vec u_t, Vec F, void* ctx)
{
  // Our ODE has the form
  //   u_t = A(t)*u.
  // Hence, the DAE has the form
  //   F(u, u_t, t) = u_t - A(t)*u = 0
  Data data = *((Data*)ctx);
  fill_time_dependent_operator(data.A, data, time);
  // F <-- A(t)*u
  int ierr = MatMult(data.A, u, F);
  // F <-- u_t - A(t)*u
  ierr = VecAYPX(F, -1.0, u_t);
  return ierr;
}

// Jacobian F_u + a * F_u', where F is DAE_F
int DAE_J(TS ts, double time, Vec u, Vec u_t, double a,
          Mat *J, Mat *P, MatStructure *flag, void *ctx)
{
  int ierr = 0;
  // For our case, this Jacobian is simply
  //   J = a*I - A
  Data data = *((Data*)ctx);
  fill_time_dependent_operator(data.J, data, time);
  ierr = MatScale(data.J, -1.0); // J <-- -A
  ierr = MatShift(data.J, a);    // J <-- -A + a*I
  //cout << "time: " << time << endl;
  return ierr;
}


int main(int argc, char** argv)
{
  PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);

  // Command line options
  int fmm = 1; // fine mesh multiplier (default 1 per assembly)
  double dt = 0.1;
  PetscOptionsGetInt(PETSC_NULL,  "-fmm", &fmm, PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, "-dt",  &dt,  PETSC_NULL);

  TS ts;                      // PETSc time stepping context
  Vec x, temp;                // Solution and temp vectors
  Data data(fmm);             // Problem data
  int N = data.n * (2 + npc); // Problem size

  // Create the unknown vector
  VecCreateSeq(PETSC_COMM_SELF, N, &x);

  // Create time stepping solver context and time range
  TSCreate(PETSC_COMM_SELF, &ts);
  TSSetProblemType(ts, TS_LINEAR);
  TSSetInitialTimeStep(ts, 0.0, dt);
  TSSetDuration(ts, 1e6, 4.0); // default end time of 50 seconds
  TSSetExactFinalTime(ts, PETSC_TRUE); // ensure we get the final time requested

  // Set our own monitor to display stuff at each time
  TSMonitorSet(ts, monitor, &data, PETSC_NULL);

  // Create matrix data structure
  MatCreateSeqAIJ(PETSC_COMM_SELF, N, N, 5 + npc, PETSC_NULL, &data.A);
  MatCreateSeqAIJ(PETSC_COMM_SELF, N, N, 5 + npc, PETSC_NULL, &data.J);

  // ODE
  // Set the right hand side function.  It all happens in Jacobian for
  // the case of a linear problem, i.e. u_t = A(u,t)*u with A(u,t) = A(t)
//  TSSetRHSFunction(ts, PETSC_NULL, TSComputeRHSFunctionLinear, &data);
//  TSSetRHSJacobian(ts, data.A, data.A, right_hand_side, &data);

  // DAE functions
  TSSetIFunction(ts, PETSC_NULL, DAE_F, &data);
  TSSetIJacobian(ts, data.J, data.J, DAE_J, &data);

  // Final setup.
  TSSetFromOptions(ts);

  // Create diffusion operators, solve for steady state, and destroy the operators.
  MatCreateSeqAIJ(PETSC_COMM_SELF, 2*data.n, 2*data.n, 5, PETSC_NULL, &data.L);
  MatCreateSeqAIJ(PETSC_COMM_SELF, 2*data.n, 2*data.n, 5, PETSC_NULL, &data.F);
  compute_steady_state(x, data);
  MatDestroy(&data.L);
  MatDestroy(&data.F);

  // Initialize time dependent diffusion operator
  fill_time_dependent_operator(data.A, data, 0.0);

  // Run the time stepper
  TSSolve(ts, x, PETSC_NULL);

  // Final power
  double power = 0;
  double *x_a;
  VecGetArray(x, &x_a);
  for (int i = 0; i < data.n; ++i)
  {
    int m = data.matid[i];
    power += (data.dx/data.keff)*(x_a[i]*nsf[m][0]+x_a[i+data.n]*nsf[m][1]);
  }
  VecRestoreArray(x, &x_a);
  PetscPrintf(PETSC_COMM_SELF, "FINAL POWER  %-11g\n", power);

  // Free work space and finalize.
  TSDestroy(&ts);
  MatDestroy(&data.A);
  MatDestroy(&data.J);
  VecDestroy(&x);
  PetscFinalize();
}

