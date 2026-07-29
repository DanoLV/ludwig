/*****************************************************************************
 *
 *  psi_solver_options.h
 *
 *  Poisson solver options.
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  (c) 2023 The University of Edinburgh
 *
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *
 *****************************************************************************/

#ifndef LUDWIG_PSI_SOLVER_OPTIONS_H
#define LUDWIG_PSI_SOLVER_OPTIONS_H

#include "util_json.h"
#include "util_petsc.h"

/* Poisson solver method */

typedef enum psi_poisson_solver_enum_s {
  PSI_POISSON_SOLVER_INVALID = 0,
  PSI_POISSON_SOLVER_SOR = 1,
  PSI_POISSON_SOLVER_PETSC = 2,
  PSI_POISSON_SOLVER_NONE = 3,
  /*CHANGE INIT - 20260512 add FFT and EWALD solver ids */
  PSI_POISSON_SOLVER_FFT = 4,    /* cuFFT-based Poisson solver (fe_electro) */
  PSI_POISSON_SOLVER_EWALD = 5   /* Ewald sum GPU solver (fe_electro_ewald) */
  /*CHANGE END - 20260512 */
} psi_poisson_solver_enum_t;

/*CHANGE INIT - 20260630 FFT Laplacian variant selectable from input */
/* Reciprocal-space operator used by the FFT Poisson solver. "discrete" uses
   the eigenvalue of the same finite-difference stencil as PETSc/SOR (so the
   solutions match); "analytic" uses the continuum k^2 (differs by O(h^2)). */

typedef enum psi_fft_laplacian_enum_s {
  PSI_FFT_LAPLACIAN_OPT_DISCRETE = 0,  /* discrete stencil eigenvalue (=PETSc) */
  PSI_FFT_LAPLACIAN_OPT_ANALYTIC = 1   /* continuum k^2 */
} psi_fft_laplacian_enum_t;
/*CHANGE END - 20260630 */

/* This is intended to be general; some components might not be relevant
   in all specific cases. */

typedef struct psi_solver_options_s psi_solver_options_t;

struct psi_solver_options_s {

  psi_poisson_solver_enum_t psolver;   /* Poisson solver id */
  /*CHANGE INIT - 20260630 FFT Laplacian variant (only used when psolver==FFT) */
  psi_fft_laplacian_enum_t fft_laplacian;
  /*CHANGE END - 20260630 */
  int maxits;                          /* Maximum iterations in solver */
  int verbose;                         /* Level of verbosity */
  int nfreq;                           /* Frequency of report */
  int nstencil;                        /* Stencil option */

  double reltol;                       /* Relative tolerance */
  double abstol;                       /* Absolute tolerance */
};

const char * psi_poisson_solver_to_string(psi_poisson_solver_enum_t mytype);
psi_poisson_solver_enum_t psi_poisson_solver_from_string(const char * str);
psi_poisson_solver_enum_t psi_poisson_solver_default(void);

psi_solver_options_t psi_solver_options_default(void);
psi_solver_options_t psi_solver_options_type(psi_poisson_solver_enum_t mytype);

int psi_solver_options_to_json(const psi_solver_options_t * opts, cJSON ** js);
int psi_solver_options_from_json(const cJSON * json, psi_solver_options_t * p);

#endif
