/*****************************************************************************
 *
 *  psi_fft.h
 *
 *  Poisson solver using cuFFT for GPU acceleration.
 *
 *  Solves: nabla^2 psi = -rho_elec / epsilon
 *
 *  Using spectral method:
 *    1. FFT(rho) -> rho_hat(k)
 *    2. psi_hat(k) = rho_hat(k) / (epsilon * k^2)
 *    3. IFFT(psi_hat) -> psi(r)
 *
 *  This solver assumes periodic boundary conditions and uniform permittivity.
 *  It is O(N log N) compared to O(N^2) or worse for iterative methods.
 *
 *  CHANGE INIT - 20260119 New FFT-based Poisson solver using cuFFT
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *****************************************************************************/

#ifndef LUDWIG_PSI_FFT_H
#define LUDWIG_PSI_FFT_H

#include "psi_solver.h"
#include "colloids.h"

typedef struct psi_solver_fft_s psi_solver_fft_t;

/* Laplacian operator type */
typedef enum {
  PSI_FFT_LAPLACIAN_ANALYTIC = 0,  /* Continuum k^2 operator */
  PSI_FFT_LAPLACIAN_DISCRETE = 1   /* Discrete stencil eigenvalue */
} psi_fft_laplacian_t;

/* Influence function type for P3M */
typedef enum {
  PSI_FFT_INFLUENCE_SIMPLE = 0,    /* Simple 1/(ε·λ(k)) */
  PSI_FFT_INFLUENCE_HOCKNEY = 1,   /* Optimal Hockney Eq. 8-22 */
  PSI_FFT_INFLUENCE_EWALD = 2      /* Ewald sum with Gaussian screening */
} psi_fft_influence_t;

/*CHANGE INIT - 20260422 Ewald deconvolution kernel selection */
/* Assignment kernel used to distribute charge onto the mesh.
 * Selects which P̂(k) is divided out in the Ewald influence function. */
typedef enum {
  PSI_FFT_DECONV_NONE    = 0,   /* No deconvolution (simple Gaussian screening) */
  PSI_FFT_DECONV_PESKIN   = 1,   /* Peskin 4-point kernel */
  PSI_FFT_DECONV_BSPLINE4 = 2,   /* Cubic B-spline, support [-2,2] */
  PSI_FFT_DECONV_BSPLINE6 = 3,   /* Quintic B-spline, support [-3,3] */
  PSI_FFT_DECONV_KB4      = 4,   /* Kaiser-Bessel, support [-2,2] */
  PSI_FFT_DECONV_PESKIN6  = 5    /* Peskin 6-point C^3 kernel (Bao et al. 2016) */
} psi_fft_deconv_t;
/*CHANGE END - 20260422 Ewald deconvolution kernel selection */

struct psi_solver_fft_s {
  psi_solver_t super;           /* Superclass block */
  psi_t * psi;                  /* Retain a reference to psi_t */

  /* FFT work arrays (device memory) */
  void * rho_real_d;            /* Real-space charge density */
  void * rho_complex_d;         /* Fourier-space charge density */
  void * psi_complex_d;         /* Fourier-space potential */
  void * psi_real_d;            /* Real-space potential */

  /* Discrete Laplacian eigenvalues (device memory) */
  void * laplacian_eigenval_d;  /* -lambda(k) for discrete Laplacian */

  /* Optimal influence function Ĝ_opt(k) (Hockney Eq. 8-22) */
  void * G_opt_d;               /* Precomputed Ĝ_opt(k) for P3M */
  psi_fft_influence_t influence_type;  /* Simple, Hockney optimal, or Ewald */
  double shape_a;               /* S2 shape parameter (typically = H = 1.0) */
  int n_alias_max;              /* Max alias index for Hockney sums */

  /* Ewald sum parameters */
  double ewald_alpha;           /* Ewald splitting parameter */
  double ewald_rcut;            /* Real-space cutoff for Ewald sum */

  /* cuFFT plans */
  void * plan_forward;          /* R2C plan */
  void * plan_backward;         /* C2R plan */

  /* Grid dimensions */
  int nx, ny, nz;               /* Local grid size */
  int ntotal[3];                /* Total grid size */

  /* Physical parameters */
  double epsilon;               /* Permittivity */
  double beta;                  /* Boltzmann factor: 1/kT */
  double lx, ly, lz;            /* System dimensions */

  /* Laplacian operator type */
  psi_fft_laplacian_t laplacian_type;  /* Analytic or discrete */

  /* Initialized flag */
  int is_initialised;
};

/* Creation and destruction */
int psi_solver_fft_create(psi_t * psi, psi_solver_fft_t ** solver);
int psi_solver_fft_create_opt(psi_t * psi, psi_fft_laplacian_t laplacian_type,
                               psi_solver_fft_t ** solver);
int psi_solver_fft_create_hockney(psi_t * psi, psi_fft_laplacian_t laplacian_type,
                                   double shape_a, int n_alias_max,
                                   psi_solver_fft_t ** solver);
/*CHANGE INIT - 20260422 deconv parameter in ewald */
int psi_solver_fft_create_ewald(psi_t * psi, psi_fft_laplacian_t laplacian_type,
                                   double alpha, double rcut,
                                   psi_fft_deconv_t deconv,
                                   psi_solver_fft_t ** psolver);
/*CHANGE END - 20260422 deconv parameter in ewald */                                   
int psi_solver_fft_free(psi_solver_fft_t ** solver);

/* Enable/disable Hockney optimal influence function after creation */
int psi_solver_fft_set_influence_hockney(psi_solver_fft_t * solver,
                                          double shape_a, int n_alias_max);

/* Enable Ewald sum influence function after creation
 * Parameters:
 *   alpha  - Ewald splitting parameter (controls real/reciprocal space division)
 *   rcut   - Real-space cutoff radius for short-range correction
 *   deconv - Assignment kernel whose P̂(k) is divided out (deconvolution)
 * Typical choice: alpha = 2*pi/L, rcut = L/2 for balanced computation */
/*CHANGE INIT - 20260422 deconv parameter in set_influence_ewald */
int psi_solver_fft_set_influence_ewald(psi_solver_fft_t * solver,
                                        double alpha, double rcut,
                                        psi_fft_deconv_t deconv);
/*CHANGE END - 20260422 deconv parameter in set_influence_ewald */

/* Main solve function */
int psi_solver_fft_solve(psi_solver_fft_t * solver, int ntimestep);

/* Solve with subgrid particle charges */
int psi_solver_fft_solve_with_subgrid(psi_solver_fft_t * solver,
                                       colloids_info_t * cinfo,
                                       int ntimestep);

/* Utility functions */
int psi_solver_fft_info(psi_solver_fft_t * solver);

/*CHANGE INIT - 20260710 spectral (FFT phase) shift of a scalar lattice field */
int psi_fft_shift_field(psi_t * psi, double * field_data,
                        double ux, double uy, double uz);
int psi_fft_kill_nyquist(psi_t * psi, double * field_data,
                         int killx, int killy, int killz);
double psi_fft_nyquist_x_norm(psi_t * psi, const double * field_data);
int psi_fft_shift_field_keepnyq(psi_t * psi, double * field_data,
                                double ux, double uy, double uz);
/*CHANGE END - 20260710 */

#endif /* LUDWIG_PSI_FFT_H */
