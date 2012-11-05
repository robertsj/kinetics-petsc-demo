#include "functions.hh"
#include <cmath>
#include <iostream>

//----------------------------------------------------------------------------//
void write_matrix(Mat M, const char* name)
{
  PetscViewer view;
  PetscViewerBinaryOpen(PETSC_COMM_SELF, name, FILE_MODE_WRITE, &view);
  MatView(M, view);
  PetscViewerDestroy(&view);
}

//----------------------------------------------------------------------------//
void write_vector(Vec V, const char* name)
{
  PetscViewer view;
  PetscViewerBinaryOpen(PETSC_COMM_SELF, name, FILE_MODE_WRITE, &view);
  VecView(V, view);
  PetscViewerDestroy(&view);
}

//----------------------------------------------------------------------------//
void update_material(Data &data, double time)
{
  for (int g = 0; g < 2; ++g)
  {
    data.tdsr[0][g] = sr[0][g];
    data.tdsr[1][g] = sr[1][g];
    double sa_in  = sr[2][g]; // "rod in"  removal cross section
    double sa_out = sr[0][g]; // "rod out" removal cross section
    if  (time <= 2.0)
      data.tdsr[2][g] = sa_in;
    else if (time <= 4.0)
      data.tdsr[2][g] = 0.5*(4.0-time)*sa_in + 0.5*(time-2.0)*sa_out;
    else if (time <= 10.0)
      data.tdsr[2][g] = sa_out;
    else if (time <= 12.0)
      data.tdsr[2][g] = 0.5*(time-10.0)*sa_in+0.5*(12.0-time)*sa_out;
    else
      data.tdsr[2][g] = sa_in;
  }
}

//----------------------------------------------------------------------------//
int fill_diffusion_operators(Data &data)
{
  // Get initial cross sections
  update_material(data, 0.0);

  // Ease notation
  double dx = data.dx;
  Mat L = data.L;
  Mat F = data.F;

  // Fill matrices
  for (int g = 0; g < 2; ++g)
  {
    for (int cell = 0; cell < data.n; ++cell)
    {
      int m = data.matid[cell];             // material id
      double leak = D[m][g]/dx;             // leakage term
      double removal = data.tdsr[m][g]*dx;  // removal term
      double fission = nsf[m][g]*dx;        // fission term
      double inscatter = -s21[m]*dx;        // scatter from 1 to 2
      int i = cell + g * data.n;
      if (cell == 0)
      {
        MatSetValue(L, i, i    , 3.0 * leak + removal,  INSERT_VALUES);
        MatSetValue(L, i, i + 1, -leak,                 INSERT_VALUES);
      }
      else if (cell == data.n - 1)
      {
        MatSetValue(L, i, i - 1, -leak,                 INSERT_VALUES);
        MatSetValue(L, i, i    , 3.0 * leak + removal,  INSERT_VALUES);
      }
      else
      {
        MatSetValue(L, i, i - 1, -leak,                 INSERT_VALUES);
        MatSetValue(L, i, i    , 2.0 * leak + removal,  INSERT_VALUES);
        MatSetValue(L, i, i + 1, -leak,                 INSERT_VALUES);
      }
      if (g == 0)
        MatSetValue(L, i + data.n, i, inscatter, INSERT_VALUES);
      MatSetValue(F, cell, i, fission, INSERT_VALUES);
    }
  }
  // Assemble matrices
  MatAssemblyBegin(L, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(L, MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(F, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(F, MAT_FINAL_ASSEMBLY);
  // write_matrix(L, "L.out");
  // write_matrix(F, "F.out");
  return 0;
}

//----------------------------------------------------------------------------//
int fill_time_dependent_operator(Mat A, Data &data, double time)
{
  // Get time-dependent cross section
  update_material(data, time);
  double dx = data.dx;

  // Fill matrix
  for (int g = 0; g < 2; ++g)
  {
    for (int cell = 0; cell < data.n; ++cell)
    {
      int m = data.matid[cell];                     // material id
      double v = velocity[g];                       // velocity
      double leak = D[m][g]/dx;                     // leakage term
      double removal = -data.tdsr[m][g]*dx;         // removal term
      double fission = nsf[m][g] * dx / data.keff;  // fission term
      double inscatter = s21[m]*dx;                 // scatter from 1 to 2
      double diag = removal;
      if (g == 0) diag += fission * (1.0 - beta_tot);
      int i = cell + data.n * g;
      if (cell == 0)
      {
        MatSetValue(A, i, i    , v * (diag - 3.0 * leak), INSERT_VALUES);
        MatSetValue(A, i, i + 1, v * leak,                INSERT_VALUES);
      }
      else if (cell == data.n - 1)
      {
        MatSetValue(A, i, i - 1, v * leak,                INSERT_VALUES);
        MatSetValue(A, i, i    , v * (diag - 3.0 * leak), INSERT_VALUES);
      }
      else
      {
        MatSetValue(A, i, i - 1, v * leak,                INSERT_VALUES);
        MatSetValue(A, i, i    , v * (diag - 2.0 * leak), INSERT_VALUES);
        MatSetValue(A, i, i + 1, v * leak,                INSERT_VALUES);
      }
      if (g == 0)
        MatSetValue(A, cell + data.n, cell, velocity[1] * inscatter, INSERT_VALUES);
      if (g == 1)
      {
        MatSetValue(A, cell, cell + data.n,
                    velocity[0] * (1.0 - beta_tot) * fission, INSERT_VALUES);
      }
      for (int j = 0; j < npc; ++j)
      {
        int prec_idx = cell + (j + 2) * data.n;
        MatSetValue(A, prec_idx, i, beta[j] * fission, INSERT_VALUES);
        if (g == 0)
        {
          MatSetValue(A, cell,     prec_idx, v * lambda[j],  INSERT_VALUES);
          MatSetValue(A, prec_idx, prec_idx, -lambda[j],     INSERT_VALUES);
        }
      }
    }
  }
  // Assemble matrices
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
  // write_matrix(A, "A.out");
  return 0;
}

//----------------------------------------------------------------------------//
int right_hand_side(TS ts, double time, Vec x, Mat *A, Mat *B,
                    MatStructure *flag, void *ctx)
{
  Data data = *((Data*) ctx); // cast void*-as-Data* and then dereference
  fill_time_dependent_operator(*A, data, time);
  MatSetOption(*A, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);
  return 0;
}

//----------------------------------------------------------------------------//
int monitor(TS ts, int step, double time, Vec x, void *ctx)
{
  // Sum the solution as a measure of power and record it.
  Data data = *((Data*) ctx);
  double power = 0;
  double *x_a;
  VecGetArray(x, &x_a);
  for (int i = 0; i < data.n; ++i)
  {
    int m = data.matid[i];
    power += (data.dx/data.keff)*(x_a[i]*nsf[m][0]+x_a[i+data.n]*nsf[m][1]);
  }
  VecRestoreArray(x, &x_a);
  data.power.push_back(power);

  // Display the time, time step, and power
  double dt;
  TSGetTimeStep(ts, &dt);
  PetscPrintf(PETSC_COMM_SELF, " %-11g  %-11g  %-11g\n", time, dt, power);

  return 0;
}

//----------------------------------------------------------------------------//
void compute_steady_state(Vec x, Data &data)
{
  // Setup operators
  fill_diffusion_operators(data);

  // Create temporary storage for flux and source
  Vec phi, temp, Q;
  VecCreateSeq(PETSC_COMM_SELF, 2*data.n, &phi);
  VecSet(phi, 1.0 / std::sqrt(2*data.n));
  VecDuplicate(phi, &temp);
  VecDuplicate(phi, &Q);

  // Create solver for inverting L
  KSP ksp;
  KSPCreate(PETSC_COMM_SELF, &ksp);
  KSPSetOperators(ksp, data.L, data.L, SAME_NONZERO_PATTERN);

  // Eigenvalue and residual norm of flux
  double keff = 1.0, resnorm = 1.0;

  // Do power iterations
  for (int i = 0; i < 5000; ++i)
  {
    VecCopy(phi, temp);
    MatMult(data.F, phi, Q);          // Q <-- F * x
    KSPSolve(ksp, Q, phi);            // phi  <-- inv(L)*phi = inv(L)*(F*phi)
    VecNorm(phi, NORM_2, &keff);      // k = |A*phi| / |phi| = |inv(L)*F*phi|
    VecScale(phi, 1 / keff);          // normalize flux
    VecAXPY(temp, -1.0, phi);         // x0 <-- inv(L)*(F*x) - x
    VecNorm(temp, NORM_2, &resnorm);  // resnorm = |x0 - x|
    if (resnorm < 1e-12) break;
  }
  data.keff = keff;

  // Normalize power to unity, and scale flux accordingly
  VecNorm(Q, NORM_1, &resnorm);
  VecScale(Q, keff/resnorm);
  VecScale(phi, keff/resnorm);

  // Get the unknown array
  double *x_a;
  VecGetArray(x, &x_a);

  // Insert unknown array into temporary vector, and copy flux into it.
  VecPlaceArray(temp, x_a);
  VecCopy(phi, temp);
  VecResetArray(temp);

  // Compute precursor concentration.
  for (int i = 0; i < npc; ++i)
  {
    VecPlaceArray(temp, x_a + (i + 2) * data.n);
    VecAXPY(temp,  beta[i] / (lambda[i] * keff), Q);
    VecResetArray(temp);
  }

  // Restore the unknown array
  VecRestoreArray(x, &x_a);

  // Clean up
  VecDestroy(&phi);
  VecDestroy(&temp);
  VecDestroy(&Q);
  KSPDestroy(&ksp);
}
