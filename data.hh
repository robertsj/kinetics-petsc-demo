#ifndef DATA_HH_
#define DATA_HH_

#include "constant_data.hh"
#include "petsc.h"

/**
 *  @struct Data
 *  @brief  Container for problem- and time-dependent data
 *
 *  A fundamental design feature of PETSc is use of "contexts", which
 *  are pointers to an object.  With respect to user data, contexts
 *  provide a way to access that data from with call back functions.
 *  For example, in the time stepping framework, a user must define
 *  a right hand side function.  When setting that function, the user
 *  can also set a context that provides access to needed to compute
 *  the right hand side (as opposed to storing all that data directly
 *  within the function itself, which could be impossible).
 */
struct Data
{
  // Constructor
  Data(int val = 1)
    : fmm(val), n(nfm * fmm), dx(10.0/fmm), keff(1.0), matid(n, 0)
  {
    // Set the material id for each fine mesh
    int cell = 0;
    for (int r = 0; r < number_regions; ++r)
      for (int i = 0; i < fmm*fm[r]; ++i, ++cell)
        matid[cell] = mt[r];
  }
  int fmm;            // fine mesh multiplier
  int n;              // total number of meshes
  double dx;          // mesh width
  double keff;        // eigenvalue
  vec_int matid;      // material id
  double tdsr[3][2];  // time dependent removal cross section
  vec_dbl power;      // power at each time step
  Mat L, F;           // diffusion loss and gain operators
  Mat A, J;           // kinetics operator and DAE jacobian
};

#endif // DATA_HH_
