/*****************************************************************************
 *
 *  psi_refined.h
 *
 *  Refined grid Poisson solver for improved electrostatic resolution.
 *
 *  This module implements a refined grid approach where:
 *  1. A finer grid (by integer factor) is created for Poisson solving
 *  2. Liquid charges from LB nodes are assigned to coincident refined nodes
 *  3. Particle charges are distributed using Peskin delta on the refined grid
 *  4. Poisson equation is solved on the refined grid using FFT
 *  5. Potential is transferred back to LB grid at coincident nodes
 *  6. Electric field at particles is interpolated using Peskin on refined grid
 *
 *  CHANGE INIT - 20260121 Refined grid Poisson solver
 *
 *****************************************************************************/

#ifndef LUDWIG_PSI_REFINED_H
#define LUDWIG_PSI_REFINED_H

#include "pe.h"
#include "coords.h"
#include "psi.h"
#include "colloids.h"
#include "hydro.h"

typedef struct psi_refined_s psi_refined_t;

/* Configuration options */
typedef struct psi_refined_options_s {
  int refinement_factor;         /* Integer refinement factor (e.g., 2, 4, 8) */
  double epsilon;                /* Permittivity */
  double beta;                   /* Boltzmann factor 1/(k_B T) */
  double e;                      /* Unit charge */
} psi_refined_options_t;

/* Default options */
#define PSI_REFINED_OPTIONS_DEFAULT { \
  .refinement_factor = 2,        \
  .epsilon = 1.0,                \
  .beta = 1.0,                   \
  .e = 1.0                       \
}

/* Refined grid structure */
struct psi_refined_s {
  pe_t * pe;                     /* Parallel environment */
  cs_t * cs;                     /* Original coordinate system (LB grid) */

  /* Refinement parameters */
  int rfactor;                   /* Refinement factor */

  /* Refined grid dimensions */
  int nx_ref, ny_ref, nz_ref;    /* Refined grid local size */
  int ntotal_ref[3];             /* Refined grid total size */
  double dx_ref;                 /* Refined grid spacing (1.0/rfactor) */

  /* LB grid dimensions (for reference) */
  int nx_lb, ny_lb, nz_lb;       /* LB grid local size */
  int ntotal_lb[3];              /* LB grid total size */

  /* Refined grid data (host arrays) */
  double * rho_ref;              /* Charge density on refined grid */
  double * psi_ref;              /* Potential on refined grid */
  double * efield_ref;           /* Electric field on refined grid (3 components) */

  /* FFT work arrays (device memory if using CUDA) */
  void * rho_real_d;             /* Real-space charge density */
  void * rho_complex_d;          /* Fourier-space charge density */
  void * psi_complex_d;          /* Fourier-space potential */
  void * psi_real_d;             /* Real-space potential */
  void * laplacian_eigenval_d;   /* Discrete Laplacian eigenvalues */

  /* cuFFT plans */
  void * plan_forward;           /* R2C FFT plan */
  void * plan_backward;          /* C2R FFT plan */

  /* Physical parameters */
  double epsilon;                /* Permittivity */
  double beta;                   /* Boltzmann factor */
  double e;                      /* Unit charge */
  double lx, ly, lz;             /* System dimensions */

  /* Reference to original psi */
  psi_t * psi_lb;                /* Original LB psi object */

  /* Status flag */
  int is_initialised;
};

/* Creation and destruction */
int psi_refined_create(pe_t * pe, cs_t * cs, psi_t * psi,
                        const psi_refined_options_t * opts,
                        psi_refined_t ** pref);
int psi_refined_free(psi_refined_t ** pref);

/* Main solve workflow */
int psi_refined_solve(psi_refined_t * pref, colloids_info_t * cinfo,
                       int ntimestep);

/* Individual steps (for debugging/testing) */

/* Step 1: Transfer fluid charges from LB to refined grid */
int psi_refined_transfer_charges_from_lb(psi_refined_t * pref);

/* Step 2: Distribute particle charges on refined grid using Peskin */
int psi_refined_distribute_particle_charges(psi_refined_t * pref,
                                             colloids_info_t * cinfo);

/* Step 3: Solve Poisson equation on refined grid */
int psi_refined_fft_solve(psi_refined_t * pref);

/* Step 4: Transfer potential from refined grid to LB grid */
int psi_refined_transfer_potential_to_lb(psi_refined_t * pref);

/* Step 5: Compute electric field on refined grid */
int psi_refined_compute_efield(psi_refined_t * pref);

/* Step 6: Interpolate electric field to particles using Peskin */
int psi_refined_interpolate_efield_to_particles(psi_refined_t * pref,
                                                  colloids_info_t * cinfo);

/* Utility functions */

/* Convert LB grid index to refined grid index */
int psi_refined_lb_to_ref_index(psi_refined_t * pref,
                                 int i_lb, int j_lb, int k_lb);

/* Get refined grid position from LB position */
void psi_refined_lb_to_ref_position(psi_refined_t * pref,
                                     const double r_lb[3],
                                     double r_ref[3]);

/* Get electric field at a position on refined grid (with interpolation) */
int psi_refined_efield_at_position(psi_refined_t * pref,
                                    const double r_ref[3],
                                    double efield[3]);

/* Reset refined grid arrays */
int psi_refined_reset(psi_refined_t * pref);

/* Info output */
int psi_refined_info(psi_refined_t * pref);

/* Peskin delta function (scaled for refined grid) */
double psi_refined_d_peskin(double r, double dx);

#endif /* LUDWIG_PSI_REFINED_H */
