/*****************************************************************************
 *
 *  psi_fft_pn.h
 *
 *  PN (Particle-Node) Level 3 correction for FFT-based Poisson solver.
 *
 *  This implements short-range corrections to eliminate errors introduced
 *  by the FFT mesh discretization near charged subgrid particles.
 *
 *  The correction works by computing:
 *  1. The difference between exact Coulomb potential and FFT mesh potential
 *  2. Modifying psi field near particles for correct Nernst-Planck transport
 *  3. Computing field correction on particles for accurate forces
 *  4. Applying Newton III reaction forces to fluid for momentum conservation
 *
 *  The correction kernel depends on whether analytic (k^2) or discrete
 *  Laplacian is used in the FFT solver.
 *
 *****************************************************************************/

#ifndef LUDWIG_PSI_FFT_PN_H
#define LUDWIG_PSI_FFT_PN_H

#include "psi.h"
#include "psi_fft.h"
#include "colloids.h"
#include "hydro.h"

/* PN correction parameters */
#define PSI_FFT_PN_RC       2.5    /* Cutoff radius in lattice units */
#define PSI_FFT_PN_NPTS_MAX 216    /* Max nodes within cutoff (6^3) */

/* Structure to store precomputed Green function differences */
typedef struct psi_fft_pn_green_s {
  int npoints;                    /* Number of points stored */
  int offsets[PSI_FFT_PN_NPTS_MAX][3];  /* Integer offsets (di, dj, dk) */
  double G_coulomb[PSI_FFT_PN_NPTS_MAX]; /* Coulomb Green function 1/(4*pi*eps*r) */
  double G_fft[PSI_FFT_PN_NPTS_MAX];     /* FFT Green function (inverse transform of 1/lambda) */
  double G_diff[PSI_FFT_PN_NPTS_MAX];    /* Difference: G_coulomb - G_fft */
} psi_fft_pn_green_t;

/* PN correction solver structure */
typedef struct psi_fft_pn_s {
  psi_solver_fft_t * fft_solver;  /* Reference to FFT solver */
  psi_t * psi;                    /* Reference to psi structure */
  psi_fft_pn_green_t green;       /* Precomputed Green functions */
  double epsilon;                 /* Permittivity */
  double beta;                    /* Boltzmann factor 1/kT */
  double cutoff;                  /* Cutoff radius */
  double G_solver_zero;           /* Cached G_solver(0) = (1/N³) Σ_k 1/(ε·λ(k)) */
  int is_initialized;             /* Initialization flag */
} psi_fft_pn_t;

/* Creation and destruction */
int psi_fft_pn_create(psi_solver_fft_t * fft_solver, psi_t * psi,
                       psi_fft_pn_t ** pn);
int psi_fft_pn_free(psi_fft_pn_t ** pn);

/* Precompute Green functions for given FFT kernel */
int psi_fft_pn_compute_green(psi_fft_pn_t * pn);

/* Main correction function: modifies psi and computes forces */
int psi_fft_pn_correction(psi_fft_pn_t * pn, colloids_info_t * cinfo,
                           hydro_t * hydro);

/* Utility: apply potential correction to psi near particles */
int psi_fft_pn_correct_potential(psi_fft_pn_t * pn, colloids_info_t * cinfo);

/* Utility: compute field correction on particles and reaction forces */
int psi_fft_pn_correct_field(psi_fft_pn_t * pn, colloids_info_t * cinfo,
                              hydro_t * hydro);

/* Info output */
int psi_fft_pn_info(psi_fft_pn_t * pn);

/* Subtract Coulomb self-field from Esub (call AFTER subgrid_update_Esub) */
int psi_fft_pn_subtract_self_field(psi_fft_pn_t * pn, colloids_info_t * cinfo);

/* Diagnostic: test self-field compensation with isolated particle charge */
int psi_fft_pn_test_self_field(psi_fft_pn_t * pn, colloids_info_t * cinfo);

/* Reference functions: compute Green function difference without interpolation
 * These are expensive O(N^3) but exact - use for validation */
double psi_fft_pn_compute_G_diff_direct(psi_fft_pn_t * pn,
                                         double rx, double ry, double rz);
void psi_fft_pn_compute_grad_G_diff_direct(psi_fft_pn_t * pn,
                                            double rx, double ry, double rz,
                                            double gradG[3]);

/* Mesh reference potential φᵐ(r) for S2-shaped particle (Hockney Eq. 8-75)
 * This is the potential arising from mesh discretization of a point charge.
 * Parameters: epsilon = permittivity, beta = 1/kT, r = distance, a = mesh spacing */
double psi_fft_pn_phi_mesh_reference(double epsilon, double beta, double r, double a);

/* Apply mesh reference potential correction near subgrid particles (Hockney Eq. 8-75)
 * Corrects psi field by adding: q_p * [G_Coulomb(r) - φᵐ(r)] for nodes within cutoff */
int psi_fft_pn_correct_potential_hockney(psi_fft_pn_t * pn, colloids_info_t * cinfo);

/*****************************************************************************
 *  Hockney P3M Optimal Influence Function (Eq. 8-22)
 *
 *  These functions implement the optimal influence function from Hockney &
 *  Eastwood "Computer Simulation Using Particles" Chapter 8.
 *
 *  Key equations:
 *    - Eq. 8-3:  S2 shape function S(r) = (48/πa⁴)(a/2 - r) for r < a/2
 *    - Eq. 8-4:  Reference force transform R̂(k) = -ikŜ²/k²
 *    - Eq. 8-9:  Assignment function transform Û(k) = [sin(kH/2)/(kH/2)]^(p+1)
 *    - Eq. 8-22: Optimal influence function Ĝ_opt(k)
 *
 *  The optimal Ĝ minimizes the error Q between mesh force and reference force.
 *****************************************************************************/

/* S2 shape function transform Ŝ(k) (Hockney Eq. 8-41)
 * Parameters: k = wavenumber magnitude, a = shape parameter */
double psi_fft_pn_S2_shape_transform(double k, double a);

/* Peskin/TSC assignment function transform Û(k) (Hockney Eq. 8-9 with p=2)
 * Parameters: kx,ky,kz = wavenumber components, H = mesh spacing */
double psi_fft_pn_peskin_transform(double kx, double ky, double kz, double H);

/* Reference potential transform Ŝ²(k)/k² for S2 shape (derived from Eq. 8-4)
 * Parameters: k = wavenumber magnitude, a = shape parameter */
double psi_fft_pn_reference_potential_transform(double k, double a);

/* Optimal influence function Ĝ_opt(k) (Hockney Eq. 8-22)
 * Minimizes mesh force error for given assignment (U), reference (R), and
 * differential operator (D).
 *
 * Parameters:
 *   kx, ky, kz   - wavenumber components
 *   H            - mesh spacing
 *   a            - S2 shape parameter
 *   epsilon      - permittivity
 *   D_hat_sq     - |D̂(k)|² from discrete Laplacian stencil (λ(k) eigenvalue)
 *   n_alias_max  - maximum alias index for sums (typically 2)
 *
 * Returns: Ĝ_opt(k) - optimal Green function in Fourier space
 */
double psi_fft_pn_optimal_influence_function(double kx, double ky, double kz,
                                              double H, double a,
                                              double epsilon,
                                              double D_hat_sq,
                                              int n_alias_max);

/*****************************************************************************
 *  Ewald Sum Corrections
 *
 *  The Ewald sum splits the Coulomb potential:
 *    φ(r) = φ_recip(r) + φ_real(r) - φ_self
 *
 *  The FFT computes φ_recip with Gaussian screening:
 *    Ĝ_ewald(k) = (4π/k²) · exp(-k²/(4α²)) / ε
 *
 *  Short-range correction adds the screened real-space part:
 *    φ_real(r) = erfc(α·r) / (4πε·r)
 *
 *  Self-energy correction removes self-interaction:
 *    φ_self = α / (√π · ε)
 *
 *****************************************************************************/

/* Ewald real-space potential: erfc(α·r) / (4πε·r)
 * Parameters: alpha = Ewald parameter, epsilon = permittivity, r = distance */
double psi_fft_pn_ewald_real_potential(double alpha, double epsilon, double r);

/* Ewald real-space field: -d/dr[erfc(α·r)/(4πε·r)]
 * Returns magnitude of radial field (multiply by r_hat for vector)
 * Parameters: alpha = Ewald parameter, epsilon = permittivity, r = distance */
double psi_fft_pn_ewald_real_field(double alpha, double epsilon, double r);

/* Ewald self-energy: α / (√π · ε)
 * Parameters: alpha = Ewald parameter, epsilon = permittivity */
double psi_fft_pn_ewald_self_energy(double alpha, double epsilon);

/* Apply Ewald real-space correction to potential near subgrid particles
 * Adds: q_p * erfc(α·r) / (4πε·r) for nodes within cutoff
 * Parameters: pn = PN correction structure, cinfo = colloid info */
int psi_fft_pn_correct_potential_ewald(psi_fft_pn_t * pn, colloids_info_t * cinfo);

/* Apply Ewald real-space correction to electric field near subgrid particles
 * This modifies the field at lattice nodes within the Ewald cutoff.
 * Parameters: pn = PN correction structure, cinfo = colloid info */
int psi_fft_pn_correct_field_ewald(psi_fft_pn_t * pn, colloids_info_t * cinfo,
                                    hydro_t * hydro);

/*****************************************************************************
 *  Complete Ewald Real-Space Correction
 *
 *  This function performs ALL Ewald real-space corrections in one call:
 *
 *  1. POTENTIAL CORRECTION (for Nernst-Planck ion transport):
 *     Adds erfc(αr)/(εr) contribution to ψ at lattice nodes near particles
 *
 *  2. FIELD CORRECTION ON PARTICLE (for subgrid particle force):
 *     Adds -∇[erfc(αr)/(εr)] from nearby lattice charges to pc->Esub
 *
 *  3. REACTION FORCE ON FLUID (for momentum conservation):
 *     Applies Newton III reaction force to fluid nodes via hydro_f_local_add
 *
 *  Call this AFTER psi_solver_fft_solve() and BEFORE subgrid force computation.
 *
 *  Parameters:
 *    pn     - PN correction structure (must have Ewald FFT solver)
 *    cinfo  - colloid information
 *    hydro  - hydrodynamics for reaction forces (can be NULL to skip)
 *
 *  Returns: 0 on success
 *****************************************************************************/
int psi_fft_pn_ewald_correction_full(psi_fft_pn_t * pn, colloids_info_t * cinfo,
                                      hydro_t * hydro);

/*****************************************************************************
 *  Complete Ewald Sum for Potential and Forces (Traditional Ewald, NOT FFT)
 *
 *  This function computes EVERYTHING using traditional Ewald summation
 *  (Fourier sum over k-vectors + real-space corrections):
 *
 *  1. POTENTIAL on all lattice nodes (stored in psi):
 *     - Fourier: sum over k-vectors with Gaussian screening
 *     - Real-space: particle->node + node->node within cutoff
 *
 *  2. FORCE on all lattice nodes (applied to hydro via hydro_f_local_add):
 *     - Fourier: sum over k-vectors
 *     - Real-space: particle->node + node->node within cutoff
 *
 *  3. FIELD (Esub) and FORCE (fex) on subgrid particles:
 *     - Fourier: sum over k-vectors evaluated at particle position
 *     - Real-space: node->particle + particle->particle within cutoff
 *
 *  This function REPLACES:
 *     - psi_solver_fft_solve (for potential)
 *     - subgrid_update_Esub (for particle field)
 *     - psi_force (for fluid force)
 *
 *  Call this INSTEAD of the FFT solver when you want pure Ewald.
 *
 *  Parameters:
 *    pn     - PN correction structure (uses Ewald alpha and rcut from solver)
 *    cinfo  - colloid information
 *    hydro  - hydrodynamics for fluid forces (can be NULL to skip)
 *
 *  Returns: 0 on success
 *****************************************************************************/
int psi_fft_pn_ewald_force_potential_full(psi_fft_pn_t * pn, colloids_info_t * cinfo,
                                           hydro_t * hydro);

#endif /* LUDWIG_PSI_FFT_PN_H */
