/*****************************************************************************
 *
 *  psi_solver.c
 *
 *  Factory function for the solver object.
 *
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  (c) 2023 The University of Edinburgh
 *
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *
 *****************************************************************************/

#include <assert.h>

/* Available implementations ... */
#include "psi_petsc.h"
#include "psi_sor.h"
/*CHANGE INIT - 20260512 add FFT solver to psi_solver_create */
#include "psi_fft.h"
/*CHANGE END - 20260512 */

/*****************************************************************************
 *
 *  psi_solver_create
 *
 *  Returns pointer to the abstract solver type.
 *
 *  We must allow that PETSc is not available, so the return value
 *  must be checked by the caller.
 *
 *****************************************************************************/

int psi_solver_create(psi_t * psi, psi_solver_t ** solver) {

  int ifail = 0;

  assert(solver && *solver == NULL);

  switch (psi->solver.psolver) {

  case (PSI_POISSON_SOLVER_PETSC):
    {
      psi_solver_petsc_t * petsc = NULL;
      ifail = psi_solver_petsc_create(psi, &petsc);
      if (ifail == 0) *solver = (psi_solver_t *) petsc;
    }
    break;

  case (PSI_POISSON_SOLVER_SOR):
    {
      psi_solver_sor_t * sor = NULL;
      ifail = psi_solver_sor_create(psi, &sor);
      if (ifail == 0) *solver = (psi_solver_t *) sor;
    }
    break;

  /*CHANGE INIT - 20260512 FFT solver selectable from input for fe_electro */
  case (PSI_POISSON_SOLVER_FFT):
    {
      psi_solver_fft_t * fft = NULL;
      /*CHANGE INIT - 20260630 FFT Laplacian variant selectable from input via
       * electrokinetics_fft_laplacian (default discrete = matches PETSc).
       * ANALYTIC uses continuum k^2 and differs from PETSc by O(h^2). */
      psi_fft_laplacian_t lap = PSI_FFT_LAPLACIAN_DISCRETE;
      if (psi->solver.fft_laplacian == PSI_FFT_LAPLACIAN_OPT_ANALYTIC) {
        lap = PSI_FFT_LAPLACIAN_ANALYTIC;
      }
      ifail = psi_solver_fft_create_opt(psi, lap, &fft);
      /*CHANGE END - 20260630 */
      if (ifail == 0) *solver = (psi_solver_t *) fft;
    }
    break;
  /*CHANGE END - 20260512 */

  default:
    ifail = -1;
  }

  return ifail;
}

/*****************************************************************************
 *
 *  psi_solver_var_epsilon_create
 *
 *  The equivalent for the version allowing electric contrast.
 *  The same comments apply.
 *
 *****************************************************************************/

int psi_solver_var_epsilon_create(psi_t * psi, var_epsilon_t user,
				  psi_solver_t ** solver) {

  int ifail = 0;

  assert(solver && *solver == NULL);

  switch (psi->solver.psolver) {

  case (PSI_POISSON_SOLVER_PETSC):
    {
      psi_solver_petsc_t * petsc = NULL;
      ifail = psi_solver_petsc_var_epsilon_create(psi, user, &petsc);
      if (ifail == 0) *solver = (psi_solver_t *) petsc;
    }
    break;

  case (PSI_POISSON_SOLVER_SOR):
    {
      psi_solver_sor_t * sor = NULL;
      ifail = psi_solver_sor_var_epsilon_create(psi, user, &sor);
      *solver = (psi_solver_t *) sor;
    }
    break;

  default:
    ifail = -1;
  }

  return ifail;
}
