/*****************************************************************************
 *
 *  psi_fft_pn.c
 *
 *  PN (Particle-Node) Level 3 correction for FFT-based Poisson solver.
 *
 *  This implements P3M-style short-range corrections for subgrid particles
 *  to eliminate errors from FFT mesh discretization.
 *
 *  The key insight is that the FFT solver computes a potential that differs
 *  from the exact Coulomb potential near charges. The difference depends on
 *  the Green function of the FFT kernel (analytic k^2 or discrete stencil).
 *
 *  G_FFT(r) = (1/N^3) * sum_k [exp(i*k·r) / (epsilon * lambda(k))]
 *
 *  The correction is:
 *  delta_phi(r) = G_Coulomb(r) - G_FFT(r)
 *
 *  where G_Coulomb(r) = 1/(4*pi*epsilon*r) for 3D Coulomb.
 *
 *  OPTIMIZATION: Green functions are precomputed at integer lattice offsets
 *  and interpolated for fractional positions. This avoids expensive O(N^3)
 *  Fourier sums at runtime.
 *
 *****************************************************************************/

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "util.h"
#include "pe.h"
#include "coords.h"
#include "memory.h"
#include "field.h"
#include "psi_fft_pn.h"
#include "util_sum.h"
#include "subgrid.h"

 /* PI constant */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Internal helper to find Green function difference at a given offset */
static double psi_fft_pn_green_diff_lookup(psi_fft_pn_t* pn, int di, int dj, int dk);

/* Compute G_solver at r=0 (finite self-energy of the mesh) */
static double psi_fft_pn_compute_G_solver_zero(psi_fft_pn_t* pn);

/*****************************************************************************
 *
 *  psi_fft_pn_phi_mesh_reference
 *
 *  REFERENCE FUNCTION: Compute the mesh reference potential φᵐ(r) for an
 *  S2-shaped (TSC) mesh particle, as defined in Hockney & Eastwood,
 *  "Computer Simulation Using Particles", Eq. 8-75, page 300.
 *
 *  This is the potential arising from the mesh discretization of a point
 *  charge using the S2 (Triangular Shaped Cloud) assignment function.
 *
 *  φᵐ(r) = (1/(4πε₀beta)) × f(ξ)
 *
 *  where ξ = 2r/a, and a is the mesh spacing (assumed = 1 in lattice units).
 *
 *  For 0 ≤ ξ ≤ 1:
 *    f(ξ) = (1/70a)[208 - 112ξ² + 56ξ⁴ - 14ξ⁵ - 8ξ⁶ + 3ξ⁷]
 *
 *  For 1 < ξ ≤ 2:
 *    f(ξ) = (1/70a)[12/ξ + 128 + 224ξ - 448ξ² + 280ξ³ - 56ξ⁴ - 14ξ⁵ + 8ξ⁶ - ξ⁷]
 *
 *  For ξ > 2:
 *    f(ξ) = 1/r  (reduces to Coulomb potential)
 *
 *  Parameters:
 *    epsilon - permittivity (ε)
 *    beta    - Boltzmann factor (1/kT)
 *    r       - separation distance
 *    a       - mesh spacing (default 1.0 in lattice units)
 *
 *  Returns:
 *    The mesh reference potential φᵐ(r)
 *
 *****************************************************************************/

double psi_fft_pn_phi_mesh_reference(double epsilon, double beta, double r, double a) {

  if (r < 1.0e-10) {
    /* At r=0, return the finite self-potential from ξ=0 limit */
    /* f(0) = (1/70a) × 208 = 208/(70a) */
    return 208.0 * beta / (70.0 * a  * epsilon);
  }

  double xi = 2.0 * r / a;  /* ξ = 2r/a (Eq. 8-73) */
  double phi_m;

  if (xi <= 1.0) {
    /* Region 0 ≤ ξ ≤ 1 */
    double xi2 = xi * xi;
    double xi4 = xi2 * xi2;
    double xi5 = xi4 * xi;
    double xi6 = xi4 * xi2;
    double xi7 = xi6 * xi;

    phi_m = (1.0 / (70.0 * a)) * (208.0 - 112.0 * xi2 + 56.0 * xi4
                                  - 14.0 * xi5 - 8.0 * xi6 + 3.0 * xi7);
  }
  else if (xi <= 2.0) {
    /* Region 1 < ξ ≤ 2 */
    double xi2 = xi * xi;
    double xi3 = xi2 * xi;
    double xi4 = xi2 * xi2;
    double xi5 = xi4 * xi;
    double xi6 = xi4 * xi2;
    double xi7 = xi6 * xi;

    phi_m = (1.0 / (70.0 * a)) * (12.0 / xi + 128.0 + 224.0 * xi
                                  - 448.0 * xi2 + 280.0 * xi3 - 56.0 * xi4
                                  - 14.0 * xi5 + 8.0 * xi6 - xi7);
  }
  else {
    /* Region ξ > 2: reduces to Coulomb 1/r */
    phi_m = 1.0 / r;
  }

  /* Include 1/(kT*ε) factor */
  return phi_m * beta / (epsilon);
}

/*****************************************************************************
 *
 *  psi_fft_pn_S2_shape_transform
 *
 *  Compute the Fourier transform Ŝ(k) of the S2 (linearly decreasing density)
 *  particle shape function, as defined in Hockney & Eastwood Eq. 8-41.
 *
 *  S2: S(r) = (48/πa⁴)(a/2 - r) for r < a/2, 0 otherwise  (Eq. 8-3)
 *
 *  The transform is (Eq. 8-41):
 *    Ŝ(k) = (12/(ka/2)⁴) * [2 - 2cos(ka/2) - (ka/2)sin(ka/2)]
 *
 *  For small k, this approaches 1.0 (normalized).
 *  For k->0, use Taylor expansion to avoid numerical issues.
 *
 *  Parameters:
 *    k - wavenumber magnitude |k|
 *    a - particle shape parameter (mesh spacing, typically 1.0)
 *
 *  Returns:
 *    Ŝ(k) - the transform of the S2 shape function
 *
 *****************************************************************************/

double psi_fft_pn_S2_shape_transform(double k, double a) {

  double ka2 = k * a / 2.0;  /* ka/2 */

  if (ka2 < 1.0e-6) {
    /* Taylor expansion for small ka/2: Ŝ ≈ 1 - (ka/2)²/10 + ... */
    return 1.0 - ka2 * ka2 / 10.0;
  }

  double ka2_2 = ka2 * ka2;
  double ka2_4 = ka2_2 * ka2_2;

  /* Eq. 8-41: Ŝ(k) = (12/(ka/2)⁴) * [2 - 2cos(ka/2) - (ka/2)sin(ka/2)] */
  double S_hat = (12.0 / ka2_4) * (2.0 - 2.0 * cos(ka2) - ka2 * sin(ka2));

  return S_hat;
}

/*****************************************************************************
 *
 *  psi_fft_pn_peskin_transform
 *
 *  Compute the Fourier transform Û(k) of the Peskin delta function used
 *  for charge assignment. This corresponds to the TSC (Triangular Shaped
 *  Cloud) assignment function in Hockney notation (p=2 in Eq. 8-9).
 *
 *  For the TSC/Peskin scheme in 1D:
 *    Û(kᵢ) = [sin(kᵢH/2) / (kᵢH/2)]^(p+1) = [sin(kᵢH/2) / (kᵢH/2)]³
 *
 *  In 3D, the total is:
 *    Û(k) = Û(kₓ) × Û(kᵧ) × Û(k_z)
 *
 *  For small k, approaches 1.0.
 *
 *  Parameters:
 *    kx, ky, kz - wavenumber components
 *    H          - mesh spacing (typically 1.0)
 *
 *  Returns:
 *    Û(k) - the transform of the Peskin/TSC assignment function
 *
 *****************************************************************************/

/* CHANGE INIT - Peskin continuous FT
 * 1D continuous Fourier transform of the Peskin 4-point kernel:
 *   P̂_1d(k) = ∫_{-2}^{2} φ_peskin(r) cos(k·r) dr  (real+symmetric)
 * Computed analytically by integrating the two-piece formula over [-2,-1], [-1,0], [0,1], [1,2].
 * The 3D transform is separable: P̂(kx,ky,kz) = P̂_1d(kx)·P̂_1d(ky)·P̂_1d(kz).
 */
static double peskin_transform_1d(double k) {

  /* Continuous FT of Peskin kernel: P̂(k) = ∫_{-2}^{2} φ(r) cos(kr) dr
   * Midpoint rule over [0,1] and [1,2], ×2 for symmetry. 1000 pts per piece. */
  if (fabs(k) < 1.0e-10) return 1.0;

  const int N = 1000;
  double sum = 0.0;
  double dr = 1.0 / N;

  for (int i = 0; i < N; i++) {
    double r = (i + 0.5) * dr;  /* [0,1] midpoints */
    double phi = 0.125 * (3.0 - 2.0*r + sqrt(1.0 + 4.0*r - 4.0*r*r));
    sum += phi * cos(k * r) * dr;
  }
  for (int i = 0; i < N; i++) {
    double r = 1.0 + (i + 0.5) * dr;  /* [1,2] midpoints */
    double phi = 0.125 * (5.0 - 2.0*r - sqrt(-7.0 + 12.0*r - 4.0*r*r));
    sum += phi * cos(k * r) * dr;
  }
  return 2.0 * sum;
}
/* CHANGE END - Peskin continuous FT */

double psi_fft_pn_peskin_transform(double kx, double ky, double kz, double H) {

  /*CHANGE INIT - use continuous FT of Peskin kernel (not sinc^p approximation) */
  // double U_hat = 1.0;
  // double kH2;
  // double sinc;
  // /* X component */
  // kH2 = kx * H / 2.0;
  // sinc = (fabs(kH2) < 1.0e-10) ? 1.0 : sin(kH2) / kH2;
  // U_hat *= sinc * sinc * sinc;  /* p=2 -> (p+1)=3 - incorrect for Ludwig Peskin */
  // /* Y component */
  // kH2 = ky * H / 2.0;
  // sinc = (fabs(kH2) < 1.0e-10) ? 1.0 : sin(kH2) / kH2;
  // U_hat *= sinc * sinc * sinc;
  // /* Z component */
  // kH2 = kz * H / 2.0;
  // sinc = (fabs(kH2) < 1.0e-10) ? 1.0 : sin(kH2) / kH2;
  // U_hat *= sinc * sinc * sinc;
  // return U_hat;

  return peskin_transform_1d(kx / H) * peskin_transform_1d(ky / H) * peskin_transform_1d(kz / H);
  /*CHANGE END - use continuous FT of Peskin kernel */
}

/*****************************************************************************
 *
 *  psi_fft_pn_reference_force_transform
 *
 *  Compute the Fourier transform R̂(k) of the reference interparticle force
 *  for an S2-shaped charge distribution (Hockney Eq. 8-4).
 *
 *  R̂(k) = -ik Ŝ² / k²
 *
 *  For the potential (not force), we need:
 *    φ̂_ref(k) = Ŝ²(k) / k²
 *
 *  This is the potential in k-space from an S2-shaped charge.
 *
 *  Parameters:
 *    k - wavenumber magnitude |k|
 *    a - S2 shape parameter (mesh spacing)
 *
 *  Returns:
 *    Ŝ²(k)/k² - reference potential transform (excluding 1/ε factor)
 *
 *****************************************************************************/

double psi_fft_pn_reference_potential_transform(double k, double a) {

  if (k < 1.0e-10) {
    /* At k=0, this diverges but we handle it separately */
    return 0.0;
  }

  double S_hat = psi_fft_pn_S2_shape_transform(k, a);
  double S_hat_sq = S_hat * S_hat;

  return S_hat_sq / (k * k);
}

/*****************************************************************************
 *
 *  psi_fft_pn_optimal_influence_function
 *
 *  Compute the optimal influence function Ĝ_opt(k) from Hockney & Eastwood
 *  Eq. 8-22, which minimizes the error Q in the mesh-calculated force.
 *
 *  For the conventional (non-interlaced) PM calculation:
 *
 *                D̂(k) · Σₙ Û²(kₙ) R̂(kₙ)
 *    Ĝ_opt(k) = ─────────────────────────────
 *                |D̂(k)|² [Σₙ Û²(kₙ)]²
 *
 *  where:
 *    - D̂(k) is the Fourier transform of the gradient difference operator
 *    - Û(k) is the Fourier transform of the charge assignment function (Peskin/TSC)
 *    - R̂(k) is the Fourier transform of the reference force (S2 shape)
 *    - kₙ = k + n·k_g are the aliased wavenumbers (k_g = 2π/H)
 *    - The sum over n includes aliases: n = (n₁, n₂, n₃) with nᵢ ∈ {...,-1,0,1,...}
 *
 *  For the potential (not force), we use R̂ = Ŝ²/k² instead of -ikŜ²/k².
 *
 *  This function computes Ĝ_opt for a single k-vector, including alias sums
 *  up to |nᵢ| ≤ n_alias_max.
 *
 *  Parameters:
 *    kx, ky, kz   - wavenumber components (in principal zone)
 *    H            - mesh spacing (1.0 in lattice units)
 *    a            - S2 shape parameter (typically = H)
 *    epsilon      - permittivity
 *    D_hat_sq     - |D̂(k)|² from your discrete Laplacian stencil
 *    n_alias_max  - maximum alias index to sum over (typically 2-3)
 *
 *  Returns:
 *    Ĝ_opt(k) - the optimal influence function value
 *
 *****************************************************************************/

double psi_fft_pn_optimal_influence_function(double kx, double ky, double kz,
                                              double H, double a,
                                              double epsilon,
                                              double D_hat_sq,
                                              int n_alias_max) {

  /* Check for k=0 mode */
  double k_mag = sqrt(kx * kx + ky * ky + kz * kz);
  if (k_mag < 1.0e-10) {
    /* At k=0, the influence function is undefined (set potential to zero) */
    return 0.0;
  }

  double k_g = 2.0 * M_PI / H;  /* Reciprocal lattice spacing */

  /* Sums for the optimal influence function (Eq. 8-22) */
  double sum_U2 = 0.0;           /* Σₙ Û²(kₙ) */
  double sum_U2_phi = 0.0;       /* Σₙ Û²(kₙ) φ̂_ref(kₙ) with φ̂_ref = Ŝ²/k² */

  /* Sum over aliases */
  for (int nx = -n_alias_max; nx <= n_alias_max; nx++) {
    for (int ny = -n_alias_max; ny <= n_alias_max; ny++) {
      for (int nz = -n_alias_max; nz <= n_alias_max; nz++) {

        /* Aliased wavenumber: kₙ = k + n·k_g */
        double kx_n = kx + nx * k_g;
        double ky_n = ky + ny * k_g;
        double kz_n = kz + nz * k_g;

        double k_n_mag = sqrt(kx_n * kx_n + ky_n * ky_n + kz_n * kz_n);

        /* Û²(kₙ) - squared Peskin/TSC transform at aliased k */
        double U_hat_n = psi_fft_pn_peskin_transform(kx_n, ky_n, kz_n, H);
        double U2_n = U_hat_n * U_hat_n;

        sum_U2 += U2_n;

        /* For potential: φ̂_ref(k) = Ŝ²(k)/k²
         * This is the Fourier transform of the reference potential
         * for an S2-shaped charge distribution.
         */
        if (k_n_mag > 1.0e-10) {
          double S_hat_n = psi_fft_pn_S2_shape_transform(k_n_mag, a);
          double phi_ref_n = S_hat_n * S_hat_n / (k_n_mag * k_n_mag);
          sum_U2_phi += U2_n * phi_ref_n;
        }
      }
    }
  }

  if (sum_U2 < 1.0e-20) {
    return 0.0;
  }

  /* Optimal influence function for POTENTIAL (adapted from Eq. 8-22):
   *
   * For the potential (not force), the optimal influence function is:
   *
   * Ĝ_opt(k) = [Σₙ Û²(kₙ) φ̂_ref(kₙ)] / [(Σₙ Û²(kₙ))² × ε]
   *
   * where φ̂_ref(k) = Ŝ²(k)/k² is the reference potential transform.
   *
   * This does NOT include the Laplacian eigenvalue λ(k) because:
   * - The standard FFT Poisson solver uses Ĝ = 1/(ε·k²)
   * - The optimal version replaces this entirely with the alias-corrected form
   * - The aliasing sums already account for the discrete mesh effects
   *
   * Note: D_hat_sq is passed but not used for the potential formulation.
   * It would be needed for the force formulation (Eq. 8-22 directly).
   */
  (void) D_hat_sq;  /* Suppress unused warning */

  double G_opt = sum_U2_phi / (sum_U2 * sum_U2 * epsilon);

  return G_opt;
}

/*****************************************************************************
 *
 *  psi_fft_pn_compute_G_diff_direct
 *
 *  REFERENCE FUNCTION: Compute G_Coulomb - G_FFT for any separation vector
 *  (rx, ry, rz) using the full Fourier sum without any interpolation or
 *  lookup table. This is expensive O(N^3) but exact.
 *
 *  Use this for validation and comparison with optimized versions.
 *
 *****************************************************************************/

double psi_fft_pn_compute_G_diff_direct(psi_fft_pn_t* pn,
                                         double rx, double ry, double rz) {

  double dist = sqrt(rx * rx + ry * ry + rz * rz);

  if (dist < 0.001) return 0.0;  /* Singularity protection */

  /* Coulomb Green function */
  double G_coulomb = 1.0 / (4.0 * M_PI * pn->epsilon * pn->beta * dist);
  // double G_coulomb = 1.0 / (pn->epsilon * pn->beta * dist);

  /* FFT Green function: sum_k cos(k·r) / (epsilon * lambda(k)) / N^3 */
  int nx = pn->fft_solver->nx;
  int ny = pn->fft_solver->ny;
  int nz = pn->fft_solver->nz;
  double ntot_inv = 1.0 / ((double)nx * ny * nz);
  stencil_t* stencil = pn->psi->stencil;
  psi_fft_laplacian_t lap_type = pn->fft_solver->laplacian_type;

  double G_fft = 0.0;

  for (int ix = 0; ix < nx; ix++) {
    double kx_idx = (ix <= nx / 2) ? (double)ix : (double)(ix - nx);
    double kx = 2.0 * M_PI * kx_idx / (double)nx;

    for (int iy = 0; iy < ny; iy++) {
      double ky_idx = (iy <= ny / 2) ? (double)iy : (double)(iy - ny);
      double ky = 2.0 * M_PI * ky_idx / (double)ny;

      for (int iz = 0; iz < nz; iz++) {
        double kz_idx = (iz <= nz / 2) ? (double)iz : (double)(iz - nz);
        double kz = 2.0 * M_PI * kz_idx / (double)nz;

        double lambda;
        if (lap_type == PSI_FFT_LAPLACIAN_DISCRETE) {
          lambda = 0.0;
          for (int sp = 0; sp < stencil->npoints; sp++) {
            double kdotc = kx * stencil->cv[sp][X] +
              ky * stencil->cv[sp][Y] +
              kz * stencil->cv[sp][Z];
            lambda += stencil->wlaplacian[sp] * cos(kdotc);
          }
        }
        else {
          lambda = kx * kx + ky * ky + kz * kz;
        }

        if (lambda < 1.0e-15) continue;

        double kdotr = kx * rx + ky * ry + kz * rz;
        G_fft += cos(kdotr) / (pn->epsilon * lambda);
      }
    }
  }

  G_fft *= ntot_inv;

  return G_coulomb - G_fft;
}

/*****************************************************************************
 *
 *  psi_fft_pn_compute_grad_G_diff_direct
 *
 *  REFERENCE FUNCTION: Compute the gradient of G_diff using finite differences
 *  of the direct Fourier sum. No interpolation or lookup.
 *
 *  gradG[i] = -dG_diff/dr[i] (negative because E = -grad(phi))
 *
 *****************************************************************************/

void psi_fft_pn_compute_grad_G_diff_direct(psi_fft_pn_t* pn,
                                            double rx, double ry, double rz,
                                            double gradG[3]) {

  double h = 0.01;  /* Small step for finite differences */

  /* Central differences for each component */
  double Gp, Gm;

  /* dG/dx */
  Gp = psi_fft_pn_compute_G_diff_direct(pn, rx + h, ry, rz);
  Gm = psi_fft_pn_compute_G_diff_direct(pn, rx - h, ry, rz);
  gradG[X] = -(Gp - Gm) / (2.0 * h);

  /* dG/dy */
  Gp = psi_fft_pn_compute_G_diff_direct(pn, rx, ry + h, rz);
  Gm = psi_fft_pn_compute_G_diff_direct(pn, rx, ry - h, rz);
  gradG[Y] = -(Gp - Gm) / (2.0 * h);

  /* dG/dz */
  Gp = psi_fft_pn_compute_G_diff_direct(pn, rx, ry, rz + h);
  Gm = psi_fft_pn_compute_G_diff_direct(pn, rx, ry, rz - h);
  gradG[Z] = -(Gp - Gm) / (2.0 * h);
}

/*****************************************************************************
 *
 *  psi_fft_pn_compute_G_solver_zero
 *
 *  Compute G_solver at r=0. This is the finite self-energy of the FFT mesh.
 *
 *  G_solver(0) = (1/N^3) * sum_k [1 / (epsilon * lambda(k))]  for k != 0
 *
 *  Unlike G_Coulomb(0) which diverges, this is finite because the mesh
 *  provides a natural regularization. This value is needed when computing
 *  the Peskin-distributed potential at a node that coincides with a
 *  Peskin support node (same source and target).
 *
 *****************************************************************************/

static double psi_fft_pn_compute_G_solver_zero(psi_fft_pn_t* pn) {

  int nx = pn->fft_solver->nx;
  int ny = pn->fft_solver->ny;
  int nz = pn->fft_solver->nz;
  double ntot_inv = 1.0 / ((double)nx * ny * nz);
  double epsilon = pn->epsilon;
  stencil_t* stencil = pn->psi->stencil;
  psi_fft_laplacian_t lap_type = pn->fft_solver->laplacian_type;

  double G_solver_0 = 0.0;

  for (int ix = 0; ix < nx; ix++) {
    double kx_idx = (ix <= nx / 2) ? (double)ix : (double)(ix - nx);
    double kx = 2.0 * M_PI * kx_idx / (double)nx;

    for (int iy = 0; iy < ny; iy++) {
      double ky_idx = (iy <= ny / 2) ? (double)iy : (double)(iy - ny);
      double ky = 2.0 * M_PI * ky_idx / (double)ny;

      for (int iz = 0; iz < nz; iz++) {
        double kz_idx = (iz <= nz / 2) ? (double)iz : (double)(iz - nz);
        double kz = 2.0 * M_PI * kz_idx / (double)nz;

        double lambda;
        if (lap_type == PSI_FFT_LAPLACIAN_DISCRETE) {
          lambda = 0.0;
          for (int sp = 0; sp < stencil->npoints; sp++) {
            double kdotc = kx * stencil->cv[sp][X] +
              ky * stencil->cv[sp][Y] +
              kz * stencil->cv[sp][Z];
            lambda += stencil->wlaplacian[sp] * cos(kdotc);
          }
        }
        else {
          lambda = kx * kx + ky * ky + kz * kz;
        }

        /* Skip k=0 mode */
        if (lambda < 1.0e-10) continue;

        /* At r=0, exp(i*k·r) = 1, so phase factor is 1 */
        G_solver_0 += 1.0 / (epsilon * lambda);
      }
    }
  }

  G_solver_0 *= ntot_inv;

  return G_solver_0;
}

/*****************************************************************************
 *
 *  psi_fft_pn_create
 *
 *  Create the PN correction structure and precompute Green functions.
 *
 *****************************************************************************/

int psi_fft_pn_create(psi_solver_fft_t* fft_solver, psi_t* psi,
                       psi_fft_pn_t** pn) {

  psi_fft_pn_t* obj = NULL;

  assert(fft_solver);
  assert(psi);
  assert(pn);

  obj = (psi_fft_pn_t*)calloc(1, sizeof(psi_fft_pn_t));
  if (obj == NULL) {
    pe_fatal(psi->pe, "psi_fft_pn_create: calloc failed\n");
    return -1;
  }

  obj->fft_solver = fft_solver;
  obj->psi = psi;
  obj->cutoff = PSI_FFT_PN_RC;
  obj->is_initialized = 0;

  psi_epsilon(psi, &obj->epsilon);
  psi_beta(psi, &obj->beta);

  /* Initialize Green function arrays */
  obj->green.npoints = 0;
  memset(obj->green.offsets, 0, sizeof(obj->green.offsets));
  memset(obj->green.G_coulomb, 0, sizeof(obj->green.G_coulomb));
  memset(obj->green.G_fft, 0, sizeof(obj->green.G_fft));
  memset(obj->green.G_diff, 0, sizeof(obj->green.G_diff));

  /* Precompute Green functions for integer offsets */
  psi_fft_pn_compute_green(obj);

  /* Precompute and cache G_solver(0) - the finite self-energy of the mesh */
  obj->G_solver_zero = psi_fft_pn_compute_G_solver_zero(obj);

  pe_info(psi->pe, "psi_fft_pn: G_solver(0) = %14.7e\n", obj->G_solver_zero);

  obj->is_initialized = 1;
  *pn = obj;

  return 0;
}

/*****************************************************************************
 *
 *  psi_fft_pn_free
 *
 *****************************************************************************/

int psi_fft_pn_free(psi_fft_pn_t** pn) {

  if (pn && *pn) {
    free(*pn);
    *pn = NULL;
  }

  return 0;
}

/*****************************************************************************
 *
 *  psi_fft_pn_green_diff_lookup
 *
 *  Look up the precomputed G_diff for integer offset (di, dj, dk).
 *  Returns 0 if offset not found.
 *
 *****************************************************************************/

static double psi_fft_pn_green_diff_lookup(psi_fft_pn_t* pn, int di, int dj, int dk) {

  for (int p = 0; p < pn->green.npoints; p++) {
    if (pn->green.offsets[p][X] == di &&
        pn->green.offsets[p][Y] == dj &&
        pn->green.offsets[p][Z] == dk) {
      return pn->green.G_diff[p];
    }
  }

  return 0.0;  /* Not found or beyond cutoff */
}

/*****************************************************************************
 *
 *  psi_fft_pn_compute_green
 *
 *  Precompute the Green function difference for integer points within cutoff.
 *
 *  For each integer offset (di, dj, dk) with |r| <= cutoff:
 *  1. G_Coulomb(r) = 1 / (4*pi*epsilon*|r|)
 *  2. G_FFT(r) = (1/N^3) * sum_k [exp(i*k·r) / (epsilon * lambda(k))]
 *  3. G_diff = G_Coulomb - G_FFT
 *
 *****************************************************************************/

int psi_fft_pn_compute_green(psi_fft_pn_t* pn) {

  int nx, ny, nz;
  int npoints = 0;
  double epsilon;
  double rc = pn->cutoff;
  double rc_sq = rc * rc;
  int irc = (int)ceil(rc);
  stencil_t* stencil = NULL;
  psi_fft_laplacian_t lap_type;

  assert(pn);

  nx = pn->fft_solver->nx;
  ny = pn->fft_solver->ny;
  nz = pn->fft_solver->nz;
  epsilon = pn->epsilon;
  lap_type = pn->fft_solver->laplacian_type;
  stencil = pn->psi->stencil;

  /* Loop over all integer offsets within cutoff */
  for (int di = -irc; di <= irc; di++) {
    for (int dj = -irc; dj <= irc; dj++) {
      for (int dk = -irc; dk <= irc; dk++) {

        double r_sq = di * di + dj * dj + dk * dk;

        /* Skip origin and points beyond cutoff */
        if (r_sq < 1.0e-10) continue;
        if (r_sq > rc_sq) continue;

        if (npoints >= PSI_FFT_PN_NPTS_MAX) {
          pe_info(pn->psi->pe, "Warning: psi_fft_pn cutoff too large\n");
          break;
        }

        double r = sqrt(r_sq);

        /* Store offset */
        pn->green.offsets[npoints][X] = di;
        pn->green.offsets[npoints][Y] = dj;
        pn->green.offsets[npoints][Z] = dk;

        /* Coulomb Green function: G_C = 1/(4*pi*epsilon*r) */
        pn->green.G_coulomb[npoints] = 1.0 / (4.0 * M_PI * epsilon * r);

        /* FFT Green function: inverse Fourier transform of 1/(epsilon*lambda) */
        double G_fft = 0.0;
        double ntot_inv = 1.0 / ((double)nx * ny * nz);

        for (int ix = 0; ix < nx; ix++) {
          for (int iy = 0; iy < ny; iy++) {
            for (int iz = 0; iz < nz; iz++) {

              /* Wave numbers */
              double kx_idx = (ix <= nx / 2) ? (double)ix : (double)(ix - nx);
              double ky_idx = (iy <= ny / 2) ? (double)iy : (double)(iy - ny);
              double kz_idx = (iz <= nz / 2) ? (double)iz : (double)(iz - nz);

              /* Wave vector components */
              double kx = 2.0 * M_PI * kx_idx / (double)nx;
              double ky = 2.0 * M_PI * ky_idx / (double)ny;
              double kz = 2.0 * M_PI * kz_idx / (double)nz;

              /* Compute eigenvalue based on kernel type */
              double lambda;

              if (lap_type == PSI_FFT_LAPLACIAN_DISCRETE) {
                /* Discrete stencil eigenvalue */
                lambda = 0.0;
                for (int p = 0; p < stencil->npoints; p++) {
                  double kdotc = kx * stencil->cv[p][X] +
                    ky * stencil->cv[p][Y] +
                    kz * stencil->cv[p][Z];
                  lambda += stencil->wlaplacian[p] * cos(kdotc);
                }
              }
              else {
                /* Analytic k^2 eigenvalue */
                lambda = kx * kx + ky * ky + kz * kz;
              }

              /* Skip k=0 mode */
              if (lambda < 1.0e-10) continue;

              /* Phase factor: exp(i*k·r) = cos(k·r) for real G */
              double kdotr = kx * di + ky * dj + kz * dk;
              double phase = cos(kdotr);

              G_fft += phase / (epsilon * lambda);
            }
          }
        }

        G_fft *= ntot_inv;
        pn->green.G_fft[npoints] = G_fft;

        /* Correction: G_diff = G_Coulomb - G_FFT */
        pn->green.G_diff[npoints] = pn->green.G_coulomb[npoints] - G_fft;

        npoints++;
      }
    }
  }

  pn->green.npoints = npoints;

  pe_info(pn->psi->pe, "psi_fft_pn: computed Green functions for %d points\n",
          npoints);

  return 0;
}

/*****************************************************************************
 *
 *  psi_fft_pn_correction
 *
 *  Main correction function. Corrects the potential ψ near subgrid particles.
 *  This must be called BEFORE subgrid_update_Esub().
 *
 *  The FFT solver produces a potential from the Peskin-distributed charge.
 *  This correction replaces it with the exact Coulomb potential from a
 *  point charge.
 *
 *  After this correction:
 *  1. Ions see the correct Coulomb potential (for Nernst-Planck transport)
 *  2. The potential is ready for Esub calculation
 *
 *  NOTE: After calling subgrid_update_Esub(), you must call
 *  psi_fft_pn_subtract_self_field() to remove the Coulomb self-field
 *  from Esub.
 *
 *****************************************************************************/

int psi_fft_pn_correction(psi_fft_pn_t* pn, colloids_info_t* cinfo,
                           hydro_t* hydro) {

  assert(pn);
  assert(cinfo);

  (void)hydro;  /* Not used when only correcting potential */

  if (cinfo->nsubgrid == 0) return 0;

  /* Correct potential in psi structure near particles */
  psi_fft_pn_correct_potential(pn, cinfo);

  return 0;
}

/*****************************************************************************
 *
 *  psi_fft_pn_correct_potential_hockney
 *
 *  Apply mesh reference potential correction using Hockney & Eastwood Eq. 8-75.
 *
 *  The FFT solver computes the potential from the particle's charge distributed
 *  via Peskin delta function. The resulting potential at nearby nodes is the
 *  mesh reference potential φᵐ(r), not the exact Coulomb potential.
 *
 *  This function corrects the potential at nodes within the cutoff radius:
 *
 *    Δψ(r_n) = q_p * [ G_Coulomb(r_n - r_p) - φᵐ(|r_n - r_p|) ]
 *
 *  where:
 *    - G_Coulomb(r) = 1/(4πεβ|r|) is the exact Coulomb Green function
 *    - φᵐ(r) is the S2 mesh reference potential (Hockney Eq. 8-75)
 *    - q_p is the particle charge
 *    - r_n is the node position
 *    - r_p is the particle position
 *
 *  The mesh spacing a = 1.0 (lattice units).
 *
 *  After this correction:
 *    1. Ions see the correct Coulomb potential (for Nernst-Planck transport)
 *    2. The potential is ready for Esub calculation
 *
 *****************************************************************************/

int psi_fft_pn_correct_potential_hockney(psi_fft_pn_t* pn, colloids_info_t* cinfo) {

  int ncell[3];
  int nlocal[3], offset[3];
  int ic, jc, kc;
  colloid_t* pc;
  psi_t* psi;
  int irc;
  double epsilon;
  double beta;
  double a;  /* S2 shape parameter from FFT solver */

  assert(pn);
  assert(cinfo);

  if (cinfo->nsubgrid == 0) return 0;

  psi = pn->psi;
  irc = (int)ceil(pn->cutoff);
  epsilon = pn->epsilon;
  beta = pn->beta;

  /* Get shape_a from FFT solver - must match the value used in optimal influence function */
  if (pn->fft_solver && pn->fft_solver->influence_type == PSI_FFT_INFLUENCE_HOCKNEY) {
    a = pn->fft_solver->shape_a;
  } else {
    /* Fallback to default if solver doesn't have Hockney influence */
    a = 1.0;  /* Default: H = 1 lattice unit */
  }

  cs_nlocal(cinfo->cs, nlocal);
  cs_nlocal_offset(cinfo->cs, offset);
  colloids_info_ncell(cinfo, ncell);

  /* Loop over all particles */
  for (ic = 0; ic <= ncell[X] + 1; ic++) {
    for (jc = 0; jc <= ncell[Y] + 1; jc++) {
      for (kc = 0; kc <= ncell[Z] + 1; kc++) {

        colloids_info_cell_list_head(cinfo, ic, jc, kc, &pc);

        for (; pc; pc = pc->next) {

          if (pc->s.bc != COLLOID_BC_SUBGRID) continue;

          /* Net charge of particle */
          double q_p = pc->s.q0 - pc->s.q1;

          /* Particle position in local coordinates */
          double r_p[3];
          r_p[X] = pc->s.r[X] - (double)offset[X];
          r_p[Y] = pc->s.r[Y] - (double)offset[Y];
          r_p[Z] = pc->s.r[Z] - (double)offset[Z];

          /* Base lattice node for search */
          int i0 = (int)floor(r_p[X]);
          int j0 = (int)floor(r_p[Y]);
          int k0 = (int)floor(r_p[Z]);

          printf("Correcting potential for particle -------------------------------------------------\n");
          /* Loop over nearby nodes within cutoff to apply correction */
          for (int di_n = -irc; di_n <= irc + 1; di_n++) {
            for (int dj_n = -irc; dj_n <= irc + 1; dj_n++) {
              for (int dk_n = -irc; dk_n <= irc + 1; dk_n++) {

                int node_i = i0 + di_n;
                int node_j = j0 + dj_n;
                int node_k = k0 + dk_n;

                /* Check if within local domain */
                if (node_i < 1 || node_i > nlocal[X]) continue;
                if (node_j < 1 || node_j > nlocal[Y]) continue;
                if (node_k < 1 || node_k > nlocal[Z]) continue;

                /* Vector from particle to this node: r_n - r_p */
                double r_pn[3];
                r_pn[X] = (double)node_i - r_p[X];
                r_pn[Y] = (double)node_j - r_p[Y];
                r_pn[Z] = (double)node_k - r_p[Z];

                double dist_pn = sqrt(r_pn[X] * r_pn[X] + r_pn[Y] * r_pn[Y] + r_pn[Z] * r_pn[Z]);
                if (dist_pn > pn->cutoff) continue;

                /* G_Coulomb(r_n - r_p) - exact point charge potential */
                double G_Coulomb_pn = 0.0;
                if (dist_pn > 1.0e-10) {
                  G_Coulomb_pn = 1.0 * beta / (epsilon * dist_pn);
                }

                /* φᵐ(r) - mesh reference potential (Hockney Eq. 8-75) */
                double phi_mesh = psi_fft_pn_phi_mesh_reference(epsilon, beta, dist_pn, a);

                /* Potential correction: Δψ = q_p * [G_Coulomb - φᵐ] */
                double delta_psi = q_p * (G_Coulomb_pn - phi_mesh);

                /* Add correction to psi field */
                int index = cs_index(cinfo->cs, node_i, node_j, node_k);
                int addr = addr_rank0(psi->nsites, index);
                printf("Node (%d,%d,%d): r_pn=%e, G_Coulomb=%e, phi_mesh_ref=%e, phi_mesh_poisson=%e, delta_psi=%e\n",
                       node_i, node_j, node_k, dist_pn, G_Coulomb_pn, phi_mesh, psi->psi->data[addr], delta_psi);
                psi->psi->data[addr] += delta_psi;
              }
            }
          }
        }
      }
    }
  }

  return 0;
}

/*****************************************************************************
 *
 *  psi_fft_pn_correct_potential
 *
 *  Modify the potential psi near each charged particle.
 *
 *  The FFT solver computes the potential from the particle's charge
 *  distributed via Peskin delta. This differs from the Coulomb potential
 *  of a point charge. The correction replaces the distributed potential
 *  with the point-charge Coulomb potential.
 *
 *  For each node r_n near the particle at r_p with charge q_p:
 *
 *    Δψ(r_n) = q_p * [ G_Coulomb(r_n - r_p) - Σ_j w_j * G_solver(r_n - r_j) ]
 *
 *  where:
 *    - G_Coulomb(r) = 1/(4πε|r|) is the exact Coulomb Green function
 *    - G_solver(r) = G_Coulomb(r) - G_diff(r) is the FFT solver Green function
 *    - w_j = d_peskin(r_p - r_j) are Peskin weights
 *    - The sum is over nodes j in the Peskin support (typically 64 nodes)
 *
 *  This correction ensures:
 *    1. Ions see the correct Coulomb potential from the particle (for N-P)
 *    2. The particle sees the correct field when Esub is interpolated
 *    3. Self-field is automatically eliminated (point charge has no self-field)
 *
 *****************************************************************************/
int psi_fft_pn_correct_potential(psi_fft_pn_t* pn, colloids_info_t* cinfo) {

  int ncell[3];
  int nlocal[3], offset[3];
  int ic, jc, kc;
  colloid_t* pc;
  psi_t* psi;
  int irc;
  double epsilon;
  double beta;

  assert(pn);
  assert(cinfo);

  psi = pn->psi;
  irc = (int)ceil(pn->cutoff);
  epsilon = pn->epsilon;
  beta = pn->beta;

  cs_nlocal(cinfo->cs, nlocal);
  cs_nlocal_offset(cinfo->cs, offset);
  colloids_info_ncell(cinfo, ncell);

  /* Loop over all particles */
  for (ic = 0; ic <= ncell[X] + 1; ic++) {
    for (jc = 0; jc <= ncell[Y] + 1; jc++) {
      for (kc = 0; kc <= ncell[Z] + 1; kc++) {

        colloids_info_cell_list_head(cinfo, ic, jc, kc, &pc);

        for (; pc; pc = pc->next) {

          if (pc->s.bc != COLLOID_BC_SUBGRID) continue;

          /* Net charge of particle */
          double q_p = pc->s.q0 - pc->s.q1;

          /* Particle position in local coordinates */
          double r_p[3];
          r_p[X] = pc->s.r[X] - (double)offset[X];
          r_p[Y] = pc->s.r[Y] - (double)offset[Y];
          r_p[Z] = pc->s.r[Z] - (double)offset[Z];

          /* Base lattice node for Peskin support */
          int i0 = (int)floor(r_p[X]);
          int j0 = (int)floor(r_p[Y]);
          int k0 = (int)floor(r_p[Z]);

          /* Precompute Peskin weights for this particle */
          /* Peskin support: nodes from (i0-1, j0-1, k0-1) to (i0+2, j0+2, k0+2) */
          double w_peskin[4][4][4];

          for (int di = -1; di <= 2; di++) {
            for (int dj = -1; dj <= 2; dj++) {
              for (int dk = -1; dk <= 2; dk++) {
                double rx = r_p[X] - (double)(i0 + di);
                double ry = r_p[Y] - (double)(j0 + dj);
                double rz = r_p[Z] - (double)(k0 + dk);
                w_peskin[di + 1][dj + 1][dk + 1] = d_peskin(rx) * d_peskin(ry) * d_peskin(rz);
              }
            }
          }

          /* Loop over nearby nodes within cutoff to apply correction */
          for (int di_n = -irc; di_n <= irc + 1; di_n++) {
            for (int dj_n = -irc; dj_n <= irc + 1; dj_n++) {
              for (int dk_n = -irc; dk_n <= irc + 1; dk_n++) {

                int node_i = i0 + di_n;
                int node_j = j0 + dj_n;
                int node_k = k0 + dk_n;

                /* Check if within local domain */
                if (node_i < 1 || node_i > nlocal[X]) continue;
                if (node_j < 1 || node_j > nlocal[Y]) continue;
                if (node_k < 1 || node_k > nlocal[Z]) continue;

                /* Vector from particle to this node: r_n - r_p */
                double r_pn[3];
                r_pn[X] = (double)node_i - r_p[X];
                r_pn[Y] = (double)node_j - r_p[Y];
                r_pn[Z] = (double)node_k - r_p[Z];

                double dist_pn = sqrt(r_pn[X] * r_pn[X] + r_pn[Y] * r_pn[Y] + r_pn[Z] * r_pn[Z]);
                if (dist_pn > pn->cutoff) continue;

                /* Term 1: G_Coulomb(r_n - r_p) - exact point charge potential */
                double G_Coulomb_pn = 0.0;
                if (dist_pn > 0.001) {
                  G_Coulomb_pn = 1.0 / (4.0 * M_PI * epsilon * beta * dist_pn);
                }

                // /* Term 2: Σ_j w_j * G_solver(r_n - r_j)
                //  * Sum over Peskin support nodes j where charge was distributed
                //  * G_solver = G_Coulomb - G_diff
                //  */
                // double G_solver_sum = 0.0;

                // for (int di_j = -1; di_j <= 2; di_j++) {
                //   for (int dj_j = -1; dj_j <= 2; dj_j++) {
                //     for (int dk_j = -1; dk_j <= 2; dk_j++) {

                //       double w_j = w_peskin[di_j + 1][dj_j + 1][dk_j + 1];
                //       // if (w_j < 1.0e-12) continue;

                //       /* Position of Peskin node j */
                //       int src_i = i0 + di_j;
                //       int src_j = j0 + dj_j;
                //       int src_k = k0 + dk_j;

                //       /* Vector from source node j to target node n: r_n - r_j */
                //       double r_jn[3];
                //       r_jn[X] = (double)(node_i - src_i);
                //       r_jn[Y] = (double)(node_j - src_j);
                //       r_jn[Z] = (double)(node_k - src_k);

                //       double dist_jn = sqrt(r_jn[X] * r_jn[X] + r_jn[Y] * r_jn[Y] + r_jn[Z] * r_jn[Z]);

                //       /* G_solver(r_n - r_j) = G_Coulomb(r_n - r_j) - G_diff(r_n - r_j) */
                //       double G_solver_jn = 0.0;

                //       if (dist_jn > 0.001) {
                //         double G_Coulomb_jn = 1.0 / (4.0 * M_PI * epsilon * beta * dist_jn);
                //         // double G_Coulomb_jn = 1.0 / (epsilon * beta * dist_jn);
                //         double G_diff_jn = psi_fft_pn_compute_G_diff_direct(pn, r_jn[X], r_jn[Y], r_jn[Z]);
                //         G_solver_jn = G_Coulomb_jn - G_diff_jn;
                //       }
                //       else {
                //         /* Same node: j = n. G_solver(0) is finite (self-energy of mesh)
                //          * G_solver(0) = (1/N³) Σ_k 1/(ε λ(k)) which is finite.
                //          * Use the cached value computed at initialization.
                //          */
                //         G_solver_jn = pn->G_solver_zero;
                //       }

                //       G_solver_sum += w_j * G_solver_jn;
                //     }
                //   }
                // }

                // /* Potential correction: replace distributed potential with point charge */
                // double delta_psi = q_p * (G_Coulomb_pn - G_solver_sum);

                double delta_psi = q_p * G_Coulomb_pn;

                /* Add correction to psi field */
                int index = cs_index(cinfo->cs, node_i, node_j, node_k);
                int addr = addr_rank0(psi->nsites, index);
                psi->psi->data[addr] += delta_psi;
              }
            }
          }
        }
      }
    }
  }

  // /* Loop over all particles */
  // for (ic = 0; ic <= ncell[X] + 1; ic++) {
  //   for (jc = 0; jc <= ncell[Y] + 1; jc++) {
  //     for (kc = 0; kc <= ncell[Z] + 1; kc++) {

  //       colloids_info_cell_list_head(cinfo, ic, jc, kc, &pc);

  //       for (; pc; pc = pc->next) {

  //         if (pc->s.bc != COLLOID_BC_SUBGRID) continue;

  //         /* Net charge of particle */
  //         double q_p = pc->s.q0 - pc->s.q1;

  //         /* Particle position in local coordinates */
  //         double r_p[3];
  //         r_p[X] = pc->s.r[X] - (double)offset[X];
  //         r_p[Y] = pc->s.r[Y] - (double)offset[Y];
  //         r_p[Z] = pc->s.r[Z] - (double)offset[Z];

  //         /* Base lattice node for Peskin support */
  //         int i0 = (int)floor(r_p[X]);
  //         int j0 = (int)floor(r_p[Y]);
  //         int k0 = (int)floor(r_p[Z]);

  //         /* Precompute Peskin weights for this particle */
  //         /* Peskin support: nodes from (i0-1, j0-1, k0-1) to (i0+2, j0+2, k0+2) */
  //         double w_peskin[4][4][4];
  //         // double sum_w = 0.0;
  //         for (int di = -1; di <= 2; di++) {
  //           for (int dj = -1; dj <= 2; dj++) {
  //             for (int dk = -1; dk <= 2; dk++) {
  //               double rx = r_p[X] - (double)(i0 + di);
  //               double ry = r_p[Y] - (double)(j0 + dj);
  //               double rz = r_p[Z] - (double)(k0 + dk);
  //               w_peskin[di + 1][dj + 1][dk + 1] = d_peskin(rx) * d_peskin(ry) * d_peskin(rz);
  //               // sum_w += w_peskin[di + 1][dj + 1][dk + 1];
  //             }
  //           }
  //         }
  //         // assert(sum_w == 1.0);

  //         /* Loop over nearby nodes within cutoff to apply correction */
  //         for (int di_n = -irc; di_n <= irc + 1; di_n++) {
  //           for (int dj_n = -irc; dj_n <= irc + 1; dj_n++) {
  //             for (int dk_n = -irc; dk_n <= irc + 1; dk_n++) {

  //               int node_i = i0 + di_n;
  //               int node_j = j0 + dj_n;
  //               int node_k = k0 + dk_n;

  //               /* Check if within local domain */
  //               if (node_i < 1 || node_i > nlocal[X]) continue;
  //               if (node_j < 1 || node_j > nlocal[Y]) continue;
  //               if (node_k < 1 || node_k > nlocal[Z]) continue;

  //               /* Vector from particle to this node: r_n - r_p */
  //               double r_pn[3];
  //               r_pn[X] = (double)node_i - r_p[X];
  //               r_pn[Y] = (double)node_j - r_p[Y];
  //               r_pn[Z] = (double)node_k - r_p[Z];

  //               double dist_pn = sqrt(r_pn[X] * r_pn[X] + r_pn[Y] * r_pn[Y] + r_pn[Z] * r_pn[Z]);
  //               if (dist_pn > pn->cutoff) continue;

  //               /* Term 1: G_Coulomb(r_n - r_p) - exact point charge potential */
  //               double G_Coulomb_pn = 0.0;
  //               if (dist_pn > 0.001) {
  //                 G_Coulomb_pn = 1.0 / (4.0 * M_PI * epsilon * beta * dist_pn);
  //                 // G_Coulomb_pn = 1.0 / (epsilon * beta * dist_pn);
  //               }

  //               /* Term 2: Σ_j w_j * G_solver(r_n - r_j)
  //                * Sum over Peskin support nodes j where charge was distributed
  //                * G_solver = G_Coulomb - G_diff
  //                */
  //               double G_solver_sum = 0.0;

  //               for (int di_j = -1; di_j <= 2; di_j++) {
  //                 for (int dj_j = -1; dj_j <= 2; dj_j++) {
  //                   for (int dk_j = -1; dk_j <= 2; dk_j++) {

  //                     double w_j = w_peskin[di_j + 1][dj_j + 1][dk_j + 1];
  //                     // if (w_j < 1.0e-12) continue;

  //                     /* Position of Peskin node j */
  //                     int src_i = i0 + di_j;
  //                     int src_j = j0 + dj_j;
  //                     int src_k = k0 + dk_j;

  //                     /* Vector from source node j to target node n: r_n - r_j */
  //                     double r_jn[3];
  //                     r_jn[X] = (double)(node_i - src_i);
  //                     r_jn[Y] = (double)(node_j - src_j);
  //                     r_jn[Z] = (double)(node_k - src_k);

  //                     double dist_jn = sqrt(r_jn[X] * r_jn[X] + r_jn[Y] * r_jn[Y] + r_jn[Z] * r_jn[Z]);

  //                     /* G_solver(r_n - r_j) = G_Coulomb(r_n - r_j) - G_diff(r_n - r_j) */
  //                     double G_solver_jn = 0.0;

  //                     if (dist_jn > 0.001) {
  //                       double G_Coulomb_jn = 1.0 / (4.0 * M_PI * epsilon * beta * dist_jn);
  //                       // double G_Coulomb_jn = 1.0 / (epsilon * beta * dist_jn);
  //                       double G_diff_jn = psi_fft_pn_compute_G_diff_direct(pn, r_jn[X], r_jn[Y], r_jn[Z]);
  //                       G_solver_jn = G_Coulomb_jn - G_diff_jn;
  //                     }
  //                     else {
  //                       /* Same node: j = n. G_solver(0) is finite (self-energy of mesh)
  //                        * G_solver(0) = (1/N³) Σ_k 1/(ε λ(k)) which is finite.
  //                        * Use the cached value computed at initialization.
  //                        */
  //                       G_solver_jn = pn->G_solver_zero;
  //                     }

  //                     G_solver_sum += w_j * G_solver_jn;
  //                   }
  //                 }
  //               }

  //               /* Potential correction: replace distributed potential with point charge */
  //               double delta_psi = q_p * (G_Coulomb_pn - G_solver_sum);

  //               /* Add correction to psi field */
  //               int index = cs_index(cinfo->cs, node_i, node_j, node_k);
  //               int addr = addr_rank0(psi->nsites, index);
  //               psi->psi->data[addr] += delta_psi;
  //             }
  //           }
  //         }
  //       }
  //     }
  //   }
  // }

  return 0;
}

/*****************************************************************************
 *
 *  psi_fft_pn_correct_field
 *
 *  Compute the field correction on each particle and apply reaction forces.
 *
 *  The field correction at particle position from each node charge is:
 *    E_corr = q_node * (-grad G_diff) evaluated at particle position
 *
 *  The particle feels: F_particle = q_particle * E_corr
 *  Newton III: F_fluid_node = -F_particle (applied directly to node)
 *
 *****************************************************************************/

int psi_fft_pn_correct_field(psi_fft_pn_t* pn, colloids_info_t* cinfo,
                              hydro_t* hydro) {

  int ncell[3];
  int nlocal[3], offset[3];
  int ic, jc, kc;
  colloid_t* pc;
  psi_t* psi;
  double kt, eunit;
  int irc;

  assert(pn);
  assert(cinfo);

  psi = pn->psi;
  irc = (int)ceil(pn->cutoff);

  /* Get kT and eunit for force conversion */
  psi_unit_charge(psi, &eunit);
  psi_beta(psi, &kt);
  kt = 1.0 / kt;  /* Convert beta to kT */

  cs_nlocal(cinfo->cs, nlocal);
  cs_nlocal_offset(cinfo->cs, offset);
  colloids_info_ncell(cinfo, ncell);

  /* Loop over all particles */
  for (ic = 0; ic <= ncell[X] + 1; ic++) {
    for (jc = 0; jc <= ncell[Y] + 1; jc++) {
      for (kc = 0; kc <= ncell[Z] + 1; kc++) {

        colloids_info_cell_list_head(cinfo, ic, jc, kc, &pc);

        for (; pc; pc = pc->next) {

          if (pc->s.bc != COLLOID_BC_SUBGRID) continue;

          double q_p = pc->s.q0 - pc->s.q1;
          // if (fabs(q_p) < 1.0e-12) continue;

          /* Klein sums for field correction at particle */
          klein_t E_corr_k[3];
          E_corr_k[X] = klein_zero();
          E_corr_k[Y] = klein_zero();
          E_corr_k[Z] = klein_zero();

          /* Particle position in local coordinates */
          double r0[3];
          r0[X] = pc->s.r[X] - (double)offset[X];
          r0[Y] = pc->s.r[Y] - (double)offset[Y];
          r0[Z] = pc->s.r[Z] - (double)offset[Z];

          int i0 = (int)floor(r0[X]);
          int j0 = (int)floor(r0[Y]);
          int k0 = (int)floor(r0[Z]);

          /* Loop over nearby nodes */
          for (int di = -irc; di <= irc + 1; di++) {
            for (int dj = -irc; dj <= irc + 1; dj++) {
              for (int dk = -irc; dk <= irc + 1; dk++) {

                int node_i = i0 + di;
                int node_j = j0 + dj;
                int node_k = k0 + dk;

                if (node_i < 1 || node_i > nlocal[X]) continue;
                if (node_j < 1 || node_j > nlocal[Y]) continue;
                if (node_k < 1 || node_k > nlocal[Z]) continue;

                /* Separation: node -> particle (for field direction) */
                double r_np[3];
                r_np[X] = r0[X] - (double)node_i;
                r_np[Y] = r0[Y] - (double)node_j;
                r_np[Z] = r0[Z] - (double)node_k;

                double dist = sqrt(r_np[X] * r_np[X] + r_np[Y] * r_np[Y] + r_np[Z] * r_np[Z]);
                if (dist < 0.01 || dist > pn->cutoff) continue;

                /* Get charge density at node */
                int index = cs_index(cinfo->cs, node_i, node_j, node_k);
                double rho0, rho1;
                psi_rho(psi, index, 0, &rho0);
                psi_rho(psi, index, 1, &rho1);
                double q_node = rho0 - rho1;

                if (fabs(q_node) < 1.0e-14) continue;

                /* Compute gradient of G_diff using direct method (reference)
                 * E = -grad(G), so E_corr = q_node * (-grad G_diff)
                 * We use r_np = particle - node, which is -r in G_diff convention
                 * G_diff is computed for particle->node, so we negate for field direction
                 */
                double gradG[3];
                psi_fft_pn_compute_grad_G_diff_direct(pn, -r_np[X], -r_np[Y], -r_np[Z], gradG);

                /* Field at particle from node charge (note sign handling) */
                double E_from_node[3];
                E_from_node[X] = q_node * gradG[X];
                E_from_node[Y] = q_node * gradG[Y];
                E_from_node[Z] = q_node * gradG[Z];

                /* Accumulate field correction */
                klein_add_double(&E_corr_k[X], E_from_node[X]);
                klein_add_double(&E_corr_k[Y], E_from_node[Y]);
                klein_add_double(&E_corr_k[Z], E_from_node[Z]);

                /* Newton III: reaction force on fluid node */
                /* F_particle = q_p * E_corr, F_node = -F_particle */
                if (hydro != NULL) {
                  double F_on_node[3];
                  F_on_node[X] = -q_p * E_from_node[X] * kt / eunit;
                  F_on_node[Y] = -q_p * E_from_node[Y] * kt / eunit;
                  F_on_node[Z] = -q_p * E_from_node[Z] * kt / eunit;

                  hydro_f_local_add(hydro, index, F_on_node);
                }
              }
            }
          }

          /* Add field correction to particle's Esub */
          pc->Esub[X] += klein_sum(&E_corr_k[X]);
          pc->Esub[Y] += klein_sum(&E_corr_k[Y]);
          pc->Esub[Z] += klein_sum(&E_corr_k[Z]);
        }
      }
    }
  }

  return 0;
}

/*****************************************************************************
 *
 *  psi_fft_pn_info
 *
 *****************************************************************************/

int psi_fft_pn_info(psi_fft_pn_t* pn) {

  assert(pn);

  pe_info(pn->psi->pe, "\n");
  pe_info(pn->psi->pe, "FFT PN Correction\n");
  pe_info(pn->psi->pe, "-----------------\n");
  pe_info(pn->psi->pe, "Cutoff radius:    %8.4f\n", pn->cutoff);
  pe_info(pn->psi->pe, "Green fn points:  %d\n", pn->green.npoints);
  pe_info(pn->psi->pe, "Kernel type:      %s\n",
          pn->fft_solver->laplacian_type == PSI_FFT_LAPLACIAN_DISCRETE ?
          "discrete" : "analytic");
  pe_info(pn->psi->pe, "\n");

  return 0;
}

/*****************************************************************************
 *
 *  psi_fft_pn_subtract_self_field
 *
 *  Subtract the Coulomb self-field from each particle's Esub.
 *
 *  After correcting the potential to Coulomb, when the particle interpolates
 *  E = -∇ψ from the surrounding nodes, it sees a non-zero self-field because
 *  ∇G_Coulomb ≠ 0 at the interpolation nodes.
 *
 *  The self-field is:
 *    E_self = q_p * Σ_i w_i * (-∇G_Coulomb)(r_i - r_p)
 *
 *  where:
 *    - w_i = d_peskin(r_p - r_i) are Peskin interpolation weights
 *    - r_i are the interpolation nodes
 *    - r_p is the particle position
 *
 *  We subtract this from Esub so the particle doesn't feel its own field.
 *
 *  IMPORTANT: This must be called AFTER subgrid_update_Esub().
 *
 *****************************************************************************/

 /*****************************************************************************
  *
  *  psi_fft_pn_test_self_field
  *
  *  DIAGNOSTIC FUNCTION: Test if the self-field is correctly compensated.
  *
  *  This function:
  *  1. Saves current psi and rho values
  *  2. Clears psi and rho (both species)
  *  3. Distributes ONLY the particle's charge via Peskin delta to rho[0]
  *     (The FFT solver computes rho_elec = rho[0] - rho[1], so we put
  *      positive charge in species 0)
  *  4. Solves Poisson via FFT for this isolated charge
  *  5. Computes E = -∇ψ at nodes and interpolates to particle (E_from_FFT)
  *  6. Computes analytic Coulomb self-field E_self_Coulomb
  *  7. Reports comparison: if working correctly, E_from_FFT shows what
  *     the FFT solver produces from the distributed charge
  *  8. Restores original psi and rho values
  *
  *  Call this BEFORE any corrections to see the raw FFT self-field.
  *
  *****************************************************************************/

int psi_fft_pn_test_self_field(psi_fft_pn_t* pn, colloids_info_t* cinfo) {

  int ncell[3];
  int nlocal[3], offset[3];
  int icc, jcc, kcc;
  colloid_t* pc;
  double epsilon;
  double beta;
  psi_t* psi;
  stencil_t* stencil;
  cs_t* cs;

  assert(pn);
  assert(cinfo);

  if (cinfo->nsubgrid == 0) return 0;

  psi = pn->psi;
  epsilon = pn->epsilon;
  beta = pn->beta;
  stencil = psi->stencil;
  cs = psi->cs;

  cs_nlocal(cs, nlocal);
  cs_nlocal_offset(cs, offset);
  colloids_info_ncell(cinfo, ncell);

  int nsites = psi->nsites;
  int nk = psi->nk;

  /* Allocate backup arrays for psi and rho (all species) */
  double* psi_backup = (double*)malloc(nsites * sizeof(double));
  double* rho_backup = (double*)malloc(nsites * nk * sizeof(double));
  assert(psi_backup);
  assert(rho_backup);

  /* Save current psi and rho */
  for (int index = 0; index < nsites; index++) {
    psi_backup[index] = psi->psi->data[addr_rank0(nsites, index)];
    for (int n = 0; n < nk; n++) {
      int irho = addr_rank1(nsites, nk, index, n);
      rho_backup[index * nk + n] = psi->rho->data[irho];
    }
  }

  pe_info(psi->pe, "\n");
  pe_info(psi->pe, "=== PSI_FFT_PN SELF-FIELD TEST ===\n");
  pe_info(psi->pe, "Testing each particle with isolated charge distribution\n");
  pe_info(psi->pe, "\n");

  /* Loop over all particles */
  for (icc = 0; icc <= ncell[X] + 1; icc++) {
    for (jcc = 0; jcc <= ncell[Y] + 1; jcc++) {
      for (kcc = 0; kcc <= ncell[Z] + 1; kcc++) {

        colloids_info_cell_list_head(cinfo, icc, jcc, kcc, &pc);

        for (; pc; pc = pc->next) {

          if (pc->s.bc != COLLOID_BC_SUBGRID) continue;

          double q_p = pc->s.q0 - pc->s.q1;

          for (int i = 0; i < 1; i++)
          {

            /* Particle position in local coordinates */
            double r_p[3];
            r_p[X] = pc->s.r[X] - (double)offset[X] + 0.1 * (double)i;
            r_p[Y] = pc->s.r[Y] - (double)offset[Y] + 0.1 * (double)i;
            r_p[Z] = pc->s.r[Z] - (double)offset[Z] + 0.1 * (double)i;

            int i0 = (int)floor(r_p[X]);
            int j0 = (int)floor(r_p[Y]);
            int k0 = (int)floor(r_p[Z]);

            pe_info(psi->pe, "Particle %d: q = %.6e, pos = (%.4f, %.4f, %.4f)\n",
                    pc->s.index, q_p, r_p[X], r_p[Y], r_p[Z]);
            pe_info(psi->pe, "Base node: (%d, %d, %d)\n", i0, j0, k0);

            /* ============================================================
             * STEP 1: Clear psi and rho (all species)
             * ============================================================ */
            for (int index = 0; index < nsites; index++) {
              psi->psi->data[addr_rank0(nsites, index)] = 0.0;
              for (int n = 0; n < nk; n++) {
                int irho = addr_rank1(nsites, nk, index, n);
                psi->rho->data[irho] = 0.0;
              }
            }

            /* ============================================================
             * STEP 2: Distribute particle charge to rho[0] via Peskin
             *         FFT solver computes rho_elec = rho[0] - rho[1]
             *         so we put positive charge in species 0
             * ============================================================ */
            double total_charge_distributed = 0.0;
            int ni = 0;
            double sum_w = 0.0;

            // Peskin - 2 neigbours
            int i_min = imax(0, (int)floor(r_p[X] - 1.0));
            int i_max = imin(nlocal[X] + 1, (int)ceil(r_p[X] + 1.0));
            int j_min = imax(0, (int)floor(r_p[Y] - 1.0));
            int j_max = imin(nlocal[Y] + 1, (int)ceil(r_p[Y] + 1.0));
            int k_min = imax(0, (int)floor(r_p[Z] - 1.0));
            int k_max = imin(nlocal[Z] + 1, (int)ceil(r_p[Z] + 1.0));

            double w = 0.0;

            for (int i = i_min; i <= i_max; i++) {
              for (int j = j_min; j <= j_max; j++) {
                for (int k = k_min; k <= k_max; k++) {

                  int index = cs_index(cinfo->cs, i, j, k);
                  double r[3];
                  r[X] = r_p[X] - 1.0 * i;
                  r[Y] = r_p[Y] - 1.0 * j;
                  r[Z] = r_p[Z] - 1.0 * k;

                  /* Peskin weight */
                  w = d_peskin(r[X]) * d_peskin(r[Y]) * d_peskin(r[Z]);
                  // if (w > 1) pe_info(psi->pe, "Peskin > 1 en nodo: %d, %d, %d\n",
                  //                                                   i, j, k);
                  sum_w += w;

                  /* Put charge in species 0 (positive species) */
                  /* rho_elec = e * (z0*rho0 + z1*rho1) = e * (rho0 - rho1) */
                  /* For unit charge: rho0 = q_p * w */
                  int irho0 = addr_rank1(nsites, nk, index, 0);
                  psi->rho->data[irho0] += q_p * w;
                  total_charge_distributed += q_p * w;
                  ni++;
                }
              }
            }

            pe_info(psi->pe, "Total charge distributed: %.10e (should be %.10e)\n",
                    total_charge_distributed, q_p);

            /* ============================================================
             * STEP 3: Sync rho to device and solve Poisson via FFT
             *         The FFT solver reads rho from device memory
             * ============================================================ */
            field_memcpy(psi->rho, tdpMemcpyHostToDevice);
            psi_solver_fft_solve(pn->fft_solver, 0);
            /* Note: psi_solver_fft_solve already syncs psi back to host */

            /* ============================================================
             * STEP 4: Compute E = -∇ψ and interpolate to particle
             *         This is what subgrid_update_Esub would compute
             * ============================================================ */
            double E_from_FFT[3] = { 0.0, 0.0, 0.0 };

            for (int di = -1; di <= 2; di++) {
              for (int dj = -1; dj <= 2; dj++) {
                for (int dk = -1; dk <= 2; dk++) {

                  int node_i = i0 + di;
                  int node_j = j0 + dj;
                  int node_k = k0 + dk;

                  if (node_i < 1 || node_i > nlocal[X]) continue;
                  if (node_j < 1 || node_j > nlocal[Y]) continue;
                  if (node_k < 1 || node_k > nlocal[Z]) continue;

                  /* Peskin weight */
                  double w_i = d_peskin(r_p[X] - (double)node_i)
                    * d_peskin(r_p[Y] - (double)node_j)
                    * d_peskin(r_p[Z] - (double)node_k);

                  if (w_i < 1.0e-12) continue;

                  /* Compute E = -grad(psi) at this node using stencil */
                  double e_node[3] = { 0.0, 0.0, 0.0 };

                  for (int p = 1; p < stencil->npoints; p++) {
                    int8_t cx = stencil->cv[p][X];
                    int8_t cy = stencil->cv[p][Y];
                    int8_t cz = stencil->cv[p][Z];

                    int index1 = cs_index(cs, node_i + cx, node_j + cy, node_k + cz);
                    double psi_val = psi->psi->data[addr_rank0(nsites, index1)];

                    e_node[X] -= stencil->wgradients[p] * cx * psi_val;
                    e_node[Y] -= stencil->wgradients[p] * cy * psi_val;
                    e_node[Z] -= stencil->wgradients[p] * cz * psi_val;
                  }

                  /* Interpolate to particle */
                  E_from_FFT[X] += w_i * e_node[X];
                  E_from_FFT[Y] += w_i * e_node[Y];
                  E_from_FFT[Z] += w_i * e_node[Z];
                }
              }
            }

            pe_info(psi->pe, "E_from_FFT (isolated charge):   (%.10e, %.10e, %.10e)\n",
                    E_from_FFT[X], E_from_FFT[Y], E_from_FFT[Z]);

            /* ============================================================
             * STEP 5: Compute analytic Coulomb self-field
             *         E_self = q_p * Σ_i w_i * (-∇G_Coulomb)(r_i - r_p)
             * ============================================================ */
            double E_self_Coulomb[3] = { 0.0, 0.0, 0.0 };

            for (int di = -1; di <= 2; di++) {
              for (int dj = -1; dj <= 2; dj++) {
                for (int dk = -1; dk <= 2; dk++) {

                  int node_i = i0 + di;
                  int node_j = j0 + dj;
                  int node_k = k0 + dk;

                  double r_ip[3];
                  r_ip[X] = (double)node_i - r_p[X];
                  r_ip[Y] = (double)node_j - r_p[Y];
                  r_ip[Z] = (double)node_k - r_p[Z];

                  double dist = sqrt(r_ip[X] * r_ip[X] + r_ip[Y] * r_ip[Y] + r_ip[Z] * r_ip[Z]);

                  if (dist < 0.01) continue;

                  double w_i = d_peskin(r_p[X] - (double)node_i)
                    * d_peskin(r_p[Y] - (double)node_j)
                    * d_peskin(r_p[Z] - (double)node_k);

                  if (w_i < 1.0e-12) continue;

                  /* -∇G_Coulomb(r_i - r_p) = (r_i - r_p) / (4πε|r|³) */
                  double dist_cubed = dist * dist * dist;
                  double grad_G[3];
                  grad_G[X] = r_ip[X] / (4.0 * M_PI * epsilon * beta * dist_cubed);
                  grad_G[Y] = r_ip[Y] / (4.0 * M_PI * epsilon * beta * dist_cubed);
                  grad_G[Z] = r_ip[Z] / (4.0 * M_PI * epsilon * beta * dist_cubed);

                  E_self_Coulomb[X] += w_i * grad_G[X];
                  E_self_Coulomb[Y] += w_i * grad_G[Y];
                  E_self_Coulomb[Z] += w_i * grad_G[Z];
                }
              }
            }

            E_self_Coulomb[X] *= q_p;
            E_self_Coulomb[Y] *= q_p;
            E_self_Coulomb[Z] *= q_p;

            pe_info(psi->pe, "E_self_Coulomb (analytic):      (%.10e, %.10e, %.10e)\n",
                    E_self_Coulomb[X], E_self_Coulomb[Y], E_self_Coulomb[Z]);

            /* ============================================================
             * STEP 6: Compare results
             * ============================================================ */
            double E_diff[3];
            E_diff[X] = E_from_FFT[X] - E_self_Coulomb[X];
            E_diff[Y] = E_from_FFT[Y] - E_self_Coulomb[Y];
            E_diff[Z] = E_from_FFT[Z] - E_self_Coulomb[Z];

            pe_info(psi->pe, "\n");
            pe_info(psi->pe, "E_from_FFT - E_self_Coulomb:    (%.10e, %.10e, %.10e)\n",
                    E_diff[X], E_diff[Y], E_diff[Z]);
            pe_info(psi->pe, "  (This shows FFT self-field differs from Coulomb self-field)\n");

            /* Also report what we currently subtract */
            pe_info(psi->pe, "\n");
            pe_info(psi->pe, "What psi_fft_pn_subtract_self_field subtracts: E_self_Coulomb\n");
            pe_info(psi->pe, "After subtraction, residual would be: E_from_FFT - E_self_Coulomb\n");
            pe_info(psi->pe, "  = (%.10e, %.10e, %.10e)\n",
                    E_diff[X], E_diff[Y], E_diff[Z]);

            pe_info(psi->pe, "\n");
            pe_info(psi->pe, "CONCLUSION: To get Esub = 0, we should subtract E_from_FFT,\n");
            pe_info(psi->pe, "            NOT E_self_Coulomb!\n");
            pe_info(psi->pe, "\n");
          }
        }
      }
    }
  }

  /* ============================================================
   * STEP 7: Restore original psi and rho (host memory)
   * ============================================================ */
  for (int index = 0; index < nsites; index++) {
    psi->psi->data[addr_rank0(nsites, index)] = psi_backup[index];
    for (int n = 0; n < nk; n++) {
      int irho = addr_rank1(nsites, nk, index, n);
      psi->rho->data[irho] = rho_backup[index * nk + n];
    }
  }

  /* Sync restored data back to device */
  field_memcpy(psi->psi, tdpMemcpyHostToDevice);
  field_memcpy(psi->rho, tdpMemcpyHostToDevice);

  free(psi_backup);
  free(rho_backup);

  pe_info(psi->pe, "=== END SELF-FIELD TEST ===\n");
  pe_info(psi->pe, "\n");

  return 0;
}

int psi_fft_pn_subtract_self_field(psi_fft_pn_t* pn, colloids_info_t* cinfo) {

  int ncell[3];
  int nlocal[3], offset[3];
  int ic, jc, kc;
  colloid_t* pc;
  double epsilon;
  double beta;

  assert(pn);
  assert(cinfo);

  if (cinfo->nsubgrid == 0) return 0;

  epsilon = pn->epsilon;
  beta = pn->beta;

  cs_nlocal(cinfo->cs, nlocal);
  cs_nlocal_offset(cinfo->cs, offset);
  colloids_info_ncell(cinfo, ncell);

  /* Loop over all particles */
  for (ic = 0; ic <= ncell[X] + 1; ic++) {
    for (jc = 0; jc <= ncell[Y] + 1; jc++) {
      for (kc = 0; kc <= ncell[Z] + 1; kc++) {

        colloids_info_cell_list_head(cinfo, ic, jc, kc, &pc);

        for (; pc; pc = pc->next) {

          if (pc->s.bc != COLLOID_BC_SUBGRID) continue;

          /* Net charge of particle */
          double q_p = pc->s.q0 - pc->s.q1;

          /* Particle position in local coordinates */
          double r_p[3];
          r_p[X] = pc->s.r[X] - (double)offset[X];
          r_p[Y] = pc->s.r[Y] - (double)offset[Y];
          r_p[Z] = pc->s.r[Z] - (double)offset[Z];

          /* Base lattice node for Peskin support */
          int i0 = (int)floor(r_p[X]);
          int j0 = (int)floor(r_p[Y]);
          int k0 = (int)floor(r_p[Z]);

          /* Accumulate self-field */
          double E_self[3] = { 0.0, 0.0, 0.0 };

          /* Loop over Peskin interpolation nodes */
          for (int di = -1; di <= 2; di++) {
            for (int dj = -1; dj <= 2; dj++) {
              for (int dk = -1; dk <= 2; dk++) {

                int node_i = i0 + di;
                int node_j = j0 + dj;
                int node_k = k0 + dk;

                /* Vector from particle to interpolation node: r_i - r_p */
                double r_ip[3];
                r_ip[X] = (double)node_i - r_p[X];
                r_ip[Y] = (double)node_j - r_p[Y];
                r_ip[Z] = (double)node_k - r_p[Z];

                double dist = sqrt(r_ip[X] * r_ip[X] + r_ip[Y] * r_ip[Y] + r_ip[Z] * r_ip[Z]);

                if (dist < 0.01) continue;  /* Skip if too close (shouldn't happen for Peskin) */

                /* Peskin weight for this node */
                double w_i = d_peskin(r_p[X] - (double)node_i)
                  * d_peskin(r_p[Y] - (double)node_j)
                  * d_peskin(r_p[Z] - (double)node_k);

                /* -∇G_Coulomb(r_i - r_p) = (r_i - r_p) / (4πε|r_i - r_p|³)
                 * This is the electric field at r_i from a charge at r_p
                 */
                double dist_cubed = dist * dist * dist;
                double grad_G[3];
                grad_G[X] = r_ip[X] * beta / (4.0 * M_PI * epsilon * beta * dist_cubed);
                grad_G[Y] = r_ip[Y] * beta / (4.0 * M_PI * epsilon * beta * dist_cubed);
                grad_G[Z] = r_ip[Z] * beta / (4.0 * M_PI * epsilon * beta * dist_cubed);

                /* Accumulate weighted field contribution */
                E_self[X] += w_i * grad_G[X];
                E_self[Y] += w_i * grad_G[Y];
                E_self[Z] += w_i * grad_G[Z];
              }
            }
          }

          /* Multiply by particle charge to get self-field */
          E_self[X] *= q_p;
          E_self[Y] *= q_p;
          E_self[Z] *= q_p;

          /* Subtract self-field from Esub */
          pc->Esub[X] -= E_self[X];
          pc->Esub[Y] -= E_self[Y];
          pc->Esub[Z] -= E_self[Z];
        }
      }
    }
  }

  return 0;
}

/*****************************************************************************
 *
 *  psi_fft_pn_ewald_real_potential
 *
 *  Ewald real-space potential contribution:
 *    φ_real(r) = erfc(α·r) / (4πε·r)
 *
 *  This is the short-range screened potential that must be added to the
 *  FFT result (which only contains the long-range Gaussian-screened part).
 *
 *****************************************************************************/

double psi_fft_pn_ewald_real_potential(double alpha, double epsilon, double r) {

  if (r < 1.0e-10) return 0.0;  /* Avoid singularity at r=0 */

  // double phi_real = erfc(alpha * r) / (4.0 * M_PI * epsilon * r);
  double phi_real = erfc(alpha * r) / (epsilon * r);

  return phi_real;
}

/*****************************************************************************
 *
 *  psi_fft_pn_ewald_real_field
 *
 *  Ewald real-space electric field magnitude (radial component):
 *    E_real(r) = -d/dr[erfc(α·r)/(4πε·r)]
 *             = [erfc(α·r)/r² + (2α/√π)·exp(-α²r²)/r] / (4πε)
 *
 *  Returns the magnitude; multiply by r̂ = r/|r| for the vector field.
 *
 *****************************************************************************/

double psi_fft_pn_ewald_real_field(double alpha, double epsilon, double r) {

  if (r < 1.0e-10) return 0.0;  /* Avoid singularity at r=0 */

  double ar = alpha * r;
  double ar2 = ar * ar;
  double r2 = r * r;

  /* E = [erfc(αr)/r² + (2α/√π)·exp(-α²r²)/r] / (4πε) */
  double term1 = erfc(ar) / r2;
  double term2 = (2.0 * alpha / sqrt(M_PI)) * exp(-ar2) / r;

  double E_real = (term1 + term2) / (4.0 * M_PI * epsilon);

  return E_real;
}

/*****************************************************************************
 *
 *  psi_fft_pn_ewald_self_energy
 *
 *  Ewald self-energy correction per unit charge squared:
 *    φ_self = α / (√π · ε)
 *
 *  The total self-energy for a particle with charge q is:
 *    E_self = q² · φ_self / 2
 *
 *****************************************************************************/

double psi_fft_pn_ewald_self_energy(double alpha, double epsilon) {

  return alpha / (sqrt(M_PI) * epsilon);
}

/*****************************************************************************
 *
 *  psi_fft_pn_correct_potential_ewald
 *
 *  Apply Ewald real-space correction to potential near subgrid particles.
 *
 *  The FFT solver with Ewald influence function computes:
 *    φ_recip = IFFT[ Ĝ_ewald(k) · ρ̂(k) ]
 *
 *  where Ĝ_ewald(k) = (4π/k²)·exp(-k²/4α²)/ε contains only the
 *  long-range Gaussian-screened part.
 *
 *  This function adds the short-range real-space contribution:
 *    φ_real(r) = Σ_j q_j · erfc(α·|r-r_j|) / (4πε·|r-r_j|)
 *
 *  for lattice nodes within the cutoff distance of each particle.
 *
 *****************************************************************************/

int psi_fft_pn_correct_potential_ewald(psi_fft_pn_t * pn, colloids_info_t * cinfo) {

  int ncell[3];
  int nlocal[3], offset[3];
  int ic, jc, kc;
  colloid_t * pc;
  psi_t * psi;
  int irc;
  double epsilon;
  double alpha;
  double rcut;

  assert(pn);
  assert(cinfo);

  if (cinfo->nsubgrid == 0) return 0;

  psi = pn->psi;
  epsilon = pn->epsilon;

  /* Get Ewald parameters from FFT solver */
  if (pn->fft_solver == NULL ||
      pn->fft_solver->influence_type != PSI_FFT_INFLUENCE_EWALD) {
    /* Not using Ewald - nothing to do */
    return 0;
  }

  alpha = pn->fft_solver->ewald_alpha;
  rcut = pn->fft_solver->ewald_rcut;
  irc = (int)ceil(rcut);

  printf("DEBUG Ewald correction: alpha=%g, rcut=%g, irc=%d, epsilon=%g\n",
         alpha, rcut, irc, epsilon);

  cs_nlocal(cinfo->cs, nlocal);
  cs_nlocal_offset(cinfo->cs, offset);
  colloids_info_ncell(cinfo, ncell);

  /* Loop over all particles */
  for (ic = 0; ic <= ncell[X] + 1; ic++) {
    for (jc = 0; jc <= ncell[Y] + 1; jc++) {
      for (kc = 0; kc <= ncell[Z] + 1; kc++) {

        colloids_info_cell_list_head(cinfo, ic, jc, kc, &pc);

        for (; pc; pc = pc->next) {

          if (pc->s.bc != COLLOID_BC_SUBGRID) continue;

          /* Net charge of particle */
          double q_p = pc->s.q0 - pc->s.q1;

          /* Particle position in local coordinates */
          double r_p[3];
          r_p[X] = pc->s.r[X] - (double)offset[X];
          r_p[Y] = pc->s.r[Y] - (double)offset[Y];
          r_p[Z] = pc->s.r[Z] - (double)offset[Z];

          /* Base lattice node for search */
          int i0 = (int)floor(r_p[X]);
          int j0 = (int)floor(r_p[Y]);
          int k0 = (int)floor(r_p[Z]);

          /* Loop over nearby nodes within cutoff */
          for (int di = -irc; di <= irc + 1; di++) {
            for (int dj = -irc; dj <= irc + 1; dj++) {
              for (int dk = -irc; dk <= irc + 1; dk++) {

                int node_i = i0 + di;
                int node_j = j0 + dj;
                int node_k = k0 + dk;

                /* Check if within local domain */
                if (node_i < 1 || node_i > nlocal[X]) continue;
                if (node_j < 1 || node_j > nlocal[Y]) continue;
                if (node_k < 1 || node_k > nlocal[Z]) continue;

                /* Vector from particle to this node */
                double r_vec[3];
                r_vec[X] = (double)node_i - r_p[X];
                r_vec[Y] = (double)node_j - r_p[Y];
                r_vec[Z] = (double)node_k - r_p[Z];

                double r = sqrt(r_vec[X]*r_vec[X] + r_vec[Y]*r_vec[Y] + r_vec[Z]*r_vec[Z]);

                if (r > rcut) continue;
                if (r < 1.0e-10) continue;  /* Skip self-point */

                /* Ewald real-space potential correction */
                double phi_real = psi_fft_pn_ewald_real_potential(alpha, epsilon, r);

                /* Add correction to psi: φ_corr = q_p · φ_real(r) */
                int index = cs_index(psi->cs, node_i, node_j, node_k);
                double old_psi = psi->psi->data[addr_rank0(psi->psi->nsites, index)];
                psi->psi->data[addr_rank0(psi->psi->nsites, index)] += q_p * phi_real;

                /* Debug first few corrections */
                static int debug_count = 0;
                if (debug_count < 10) {
                  printf("DEBUG Ewald corr: r=%g, phi_real=%g, q_p=%g, old_psi=%g, new_psi=%g\n",
                         r, phi_real, q_p, old_psi,
                         psi->psi->data[addr_rank0(psi->psi->nsites, index)]);
                  debug_count++;
                }
              }
            }
          }
        }
      }
    }
  }

  return 0;
}

/*****************************************************************************
 *
 *  psi_fft_pn_correct_field_ewald
 *
 *  Apply Ewald real-space correction to electric field.
 *
 *  This is needed because the field E = -∇φ computed from the FFT potential
 *  only contains the long-range part. The short-range Ewald contribution
 *  must be added:
 *
 *    E_real(r) = Σ_j q_j · [erfc(α|r-r_j|)/|r-r_j|² +
 *                           (2α/√π)exp(-α²|r-r_j|²)/|r-r_j|] · (r-r_j)/|r-r_j|
 *                / (4πε)
 *
 *  This affects:
 *  1. Forces on subgrid particles (computed via field interpolation)
 *  2. Nernst-Planck ion fluxes (drift term)
 *  3. Electrokinetic stress on fluid
 *
 *  Currently this function corrects the field stored in Esub for particles.
 *  For the fluid/ions, the potential correction should be sufficient if
 *  the discrete gradient stencil is used consistently.
 *
 *****************************************************************************/

int psi_fft_pn_correct_field_ewald(psi_fft_pn_t * pn, colloids_info_t * cinfo,
                                    hydro_t * hydro) {

  int ncell[3];
  int nlocal[3], offset[3];
  int ic, jc, kc;
  colloid_t * pc;
  double alpha;
  double rcut;
  double epsilon;

  assert(pn);
  assert(cinfo);
  (void) hydro;  /* Not used yet - could add fluid force corrections */

  if (cinfo->nsubgrid == 0) return 0;

  /* Get Ewald parameters from FFT solver */
  if (pn->fft_solver == NULL ||
      pn->fft_solver->influence_type != PSI_FFT_INFLUENCE_EWALD) {
    return 0;
  }

  alpha = pn->fft_solver->ewald_alpha;
  rcut = pn->fft_solver->ewald_rcut;
  epsilon = pn->epsilon;

  cs_nlocal(cinfo->cs, nlocal);
  cs_nlocal_offset(cinfo->cs, offset);
  colloids_info_ncell(cinfo, ncell);

  /* Loop over all target particles to compute field correction at their positions */
  for (ic = 0; ic <= ncell[X] + 1; ic++) {
    for (jc = 0; jc <= ncell[Y] + 1; jc++) {
      for (kc = 0; kc <= ncell[Z] + 1; kc++) {

        colloid_t * pc_target;
        colloids_info_cell_list_head(cinfo, ic, jc, kc, &pc_target);

        for (; pc_target; pc_target = pc_target->next) {

          if (pc_target->s.bc != COLLOID_BC_SUBGRID) continue;

          double r_target[3];
          r_target[X] = pc_target->s.r[X];
          r_target[Y] = pc_target->s.r[Y];
          r_target[Z] = pc_target->s.r[Z];

          double E_corr[3] = {0.0, 0.0, 0.0};

          /* Sum contributions from all source particles within cutoff */
          for (int ic2 = 0; ic2 <= ncell[X] + 1; ic2++) {
            for (int jc2 = 0; jc2 <= ncell[Y] + 1; jc2++) {
              for (int kc2 = 0; kc2 <= ncell[Z] + 1; kc2++) {

                colloids_info_cell_list_head(cinfo, ic2, jc2, kc2, &pc);

                for (; pc; pc = pc->next) {

                  if (pc->s.bc != COLLOID_BC_SUBGRID) continue;
                  if (pc == pc_target) continue;  /* Skip self */

                  double q_source = pc->s.q0 - pc->s.q1;

                  /* Vector from source to target */
                  double r_vec[3];
                  r_vec[X] = r_target[X] - pc->s.r[X];
                  r_vec[Y] = r_target[Y] - pc->s.r[Y];
                  r_vec[Z] = r_target[Z] - pc->s.r[Z];

                  double r = sqrt(r_vec[X]*r_vec[X] + r_vec[Y]*r_vec[Y] + r_vec[Z]*r_vec[Z]);

                  if (r > rcut) continue;
                  if (r < 1.0e-10) continue;

                  /* Ewald real-space field magnitude */
                  double E_mag = psi_fft_pn_ewald_real_field(alpha, epsilon, r);

                  /* Add contribution: E = q · E_mag · r̂ */
                  double r_inv = 1.0 / r;
                  E_corr[X] += q_source * E_mag * r_vec[X] * r_inv;
                  E_corr[Y] += q_source * E_mag * r_vec[Y] * r_inv;
                  E_corr[Z] += q_source * E_mag * r_vec[Z] * r_inv;
                }
              }
            }
          }

          /* Add correction to particle's interpolated field */
          pc_target->Esub[X] += E_corr[X];
          pc_target->Esub[Y] += E_corr[Y];
          pc_target->Esub[Z] += E_corr[Z];
        }
      }
    }
  }

  return 0;
}

/*****************************************************************************
 *
 *  psi_fft_pn_ewald_correction_full
 *
 *  Complete Ewald real-space correction for FFT solver with subgrid particles.
 *
 *  The FFT solver with Ewald influence function computes only the long-range
 *  (reciprocal space) part:
 *    φ_FFT(k) = ρ̂(k) · (1/k²) · exp(-k²/4α²) / ε
 *
 *  This function adds the short-range (real-space) corrections:
 *
 *  1. POTENTIAL CORRECTION ON LATTICE NODES (for Nernst-Planck ion transport):
 *     For each node within rcut of a particle:
 *       Δψ(r_node) = q_particle · erfc(α|r|) / (ε|r|)
 *     where r = r_node - r_particle
 *
 *  2. FIELD CORRECTION ON PARTICLE (for subgrid particle force):
 *     The particle feels the real-space field from nearby lattice charges:
 *       E_real = Σ_nodes q_node · [-∇(erfc(αr)/(εr))]
 *     This is added to pc->Esub
 *
 *  3. REACTION FORCE ON FLUID (for momentum conservation - Newton III):
 *     F_fluid(node) = -F_particle contribution from that node
 *     Applied via hydro_f_local_add()
 *
 *  IMPORTANT: The β factor is NOT included here because the potential ψ
 *  in Ludwig is already dimensionless: ψ* = β·e·ψ_physical
 *  The field Esub is also dimensionless: E* = β·e·E_physical
 *
 *  Parameters:
 *    pn     - PN correction structure (contains FFT solver with Ewald params)
 *    cinfo  - colloid information
 *    hydro  - hydrodynamics structure for reaction forces (can be NULL)
 *
 *****************************************************************************/

int psi_fft_pn_ewald_correction_full(psi_fft_pn_t * pn, colloids_info_t * cinfo,
                                      hydro_t * hydro) {

  int ncell[3];
  int nlocal[3], offset[3];
  int ic, jc, kc;
  colloid_t * pc;
  psi_t * psi;
  double alpha;
  double rcut;
  double epsilon;
  double beta;
  double eunit;
  double kt;
  int irc;

  assert(pn);
  assert(cinfo);

  if (cinfo->nsubgrid == 0) return 0;

  psi = pn->psi;
  epsilon = pn->epsilon;
  beta = pn->beta;
  psi_unit_charge(psi, &eunit);
  kt = 1.0 / beta;

  /* Get Ewald parameters from FFT solver */
  if (pn->fft_solver == NULL ||
      pn->fft_solver->influence_type != PSI_FFT_INFLUENCE_EWALD) {
    /* Not using Ewald - nothing to do */
    pe_info(psi->pe, "WARNING: psi_fft_pn_ewald_correction_full called without Ewald solver\n");
    return 0;
  }

  alpha = pn->fft_solver->ewald_alpha;
  rcut = pn->fft_solver->ewald_rcut;
  irc = (int)ceil(rcut);

  cs_nlocal(cinfo->cs, nlocal);
  cs_nlocal_offset(cinfo->cs, offset);
  colloids_info_ncell(cinfo, ncell);

  /* ========================================================================
   * PART 1: Correct potential on lattice nodes near particles
   *         This is needed for correct Nernst-Planck ion transport
   *
   *         The FFT solver uses ψ* = β·e·ψ_physical (dimensionless potential)
   *         The FFT input is: ρ_real = ρ_elec · β (see copy_rho_to_real_kernel)
   *         So the correction must also include β:
   *           Δψ* = q_p · β · erfc(αr) / (ε·r)
   * ======================================================================== */

  for (ic = 0; ic <= ncell[X] + 1; ic++) {
    for (jc = 0; jc <= ncell[Y] + 1; jc++) {
      for (kc = 0; kc <= ncell[Z] + 1; kc++) {

        colloids_info_cell_list_head(cinfo, ic, jc, kc, &pc);

        for (; pc; pc = pc->next) {

          if (pc->s.bc != COLLOID_BC_SUBGRID) continue;

          double q_p = pc->s.q0 - pc->s.q1;

          /* Particle position in local coordinates */
          double r_p[3];
          r_p[X] = pc->s.r[X] - (double)offset[X];
          r_p[Y] = pc->s.r[Y] - (double)offset[Y];
          r_p[Z] = pc->s.r[Z] - (double)offset[Z];

          int i0 = (int)floor(r_p[X]);
          int j0 = (int)floor(r_p[Y]);
          int k0 = (int)floor(r_p[Z]);

          /* Loop over nearby nodes within cutoff */
          for (int di = -irc; di <= irc + 1; di++) {
            for (int dj = -irc; dj <= irc + 1; dj++) {
              for (int dk = -irc; dk <= irc + 1; dk++) {

                int node_i = i0 + di;
                int node_j = j0 + dj;
                int node_k = k0 + dk;

                if (node_i < 1 || node_i > nlocal[X]) continue;
                if (node_j < 1 || node_j > nlocal[Y]) continue;
                if (node_k < 1 || node_k > nlocal[Z]) continue;

                /* Vector from particle to node */
                double r_vec[3];
                r_vec[X] = (double)node_i - r_p[X];
                r_vec[Y] = (double)node_j - r_p[Y];
                r_vec[Z] = (double)node_k - r_p[Z];

                double r = sqrt(r_vec[X]*r_vec[X] + r_vec[Y]*r_vec[Y] + r_vec[Z]*r_vec[Z]);

                if (r > rcut || r < 1.0e-10) continue;

                /* Ewald real-space potential:
                 * phi_physical = erfc(ar) / (4*pi*eps*r)
                 * psi* = beta * eunit * phi_physical */
                PI_DOUBLE(pi);
                double phi_real = beta * eunit * erfc(alpha * r) / (4.0 * pi * epsilon * r);

                /* Add correction: Delta_psi* = q_p * phi_real */
                int index = cs_index(psi->cs, node_i, node_j, node_k);
                psi->psi->data[addr_rank0(psi->psi->nsites, index)] += q_p * phi_real;
              }
            }
          }
        }
      }
    }
  }

  /* ========================================================================
   * PART 2: Correct field on particles from nearby lattice charges
   *         Adds Ewald real-space E-field correction to particle Esub
   * ======================================================================== */

  PI_DOUBLE(pi);

  for (ic = 0; ic <= ncell[X] + 1; ic++) {
    for (jc = 0; jc <= ncell[Y] + 1; jc++) {
      for (kc = 0; kc <= ncell[Z] + 1; kc++) {

        colloids_info_cell_list_head(cinfo, ic, jc, kc, &pc);

        for (; pc; pc = pc->next) {

          if (pc->s.bc != COLLOID_BC_SUBGRID) continue;

          /* Particle position in local coordinates */
          double r_p2[3];
          r_p2[X] = pc->s.r[X] - (double)offset[X];
          r_p2[Y] = pc->s.r[Y] - (double)offset[Y];
          r_p2[Z] = pc->s.r[Z] - (double)offset[Z];

          int i02 = (int)floor(r_p2[X]);
          int j02 = (int)floor(r_p2[Y]);
          int k02 = (int)floor(r_p2[Z]);

          /* Accumulate field correction on particle from lattice charges */
          double E_corr_total[3] = {0.0, 0.0, 0.0};

          /* Loop over nearby lattice nodes within cutoff */
          for (int di = -irc; di <= irc + 1; di++) {
            for (int dj = -irc; dj <= irc + 1; dj++) {
              for (int dk = -irc; dk <= irc + 1; dk++) {

                int node_i = i02 + di;
                int node_j = j02 + dj;
                int node_k = k02 + dk;

                if (node_i < 1 || node_i > nlocal[X]) continue;
                if (node_j < 1 || node_j > nlocal[Y]) continue;
                if (node_k < 1 || node_k > nlocal[Z]) continue;

                int index = cs_index(psi->cs, node_i, node_j, node_k);

                /* Get charge density at lattice node */
                double rho0, rho1;
                psi_rho(psi, index, 0, &rho0);
                psi_rho(psi, index, 1, &rho1);
                double q_node = rho0 - rho1;

                if (fabs(q_node) < 1.0e-14) continue;

                /* Vector from node to particle */
                double r_vec[3];
                r_vec[X] = r_p2[X] - (double)node_i;
                r_vec[Y] = r_p2[Y] - (double)node_j;
                r_vec[Z] = r_p2[Z] - (double)node_k;

                double r = sqrt(r_vec[X]*r_vec[X] + r_vec[Y]*r_vec[Y] + r_vec[Z]*r_vec[Z]);

                if (r > rcut || r < 1.0e-10) continue;

                /* Ewald real-space field:
                 * E_physical = [erfc(ar)/r^2 + 2*a/sqrt(pi)*exp(-a^2r^2)/r] / (4*pi*eps)
                 * Esub = beta * eunit * E_physical */
                double ar = alpha * r;
                double ar2 = ar * ar;
                double r2 = r * r;
                double r_inv = 1.0 / r;

                double term1 = erfc(ar) / r2;
                double term2 = (2.0 * alpha / sqrt(pi)) * exp(-ar2) * r_inv;
                double E_mag = beta * eunit * (term1 + term2) / (4.0 * pi * epsilon);

                /* Field at particle from this node charge */
                E_corr_total[X] += q_node * E_mag * r_vec[X] * r_inv;
                E_corr_total[Y] += q_node * E_mag * r_vec[Y] * r_inv;
                E_corr_total[Z] += q_node * E_mag * r_vec[Z] * r_inv;
              }
            }
          }

          /* Add field correction to particle's Esub */
          pc->Esub[X] += E_corr_total[X];
          pc->Esub[Y] += E_corr_total[Y];
          pc->Esub[Z] += E_corr_total[Z];
        }
      }
    }
  }

  /* ========================================================================
   * PART 3: Calculate Ewald real-space force on ALL fluid nodes
   *         and apply via hydro_f_local_add
   *
   *         For each node with charge q1, compute force from:
   *         - All subgrid particles within cutoff
   *         - All other lattice nodes within cutoff
   *
   *         F = q1 * q2 * [erfc(ar)/r^2 + 2a/sqrt(pi)*exp(-a^2r^2)/r]
   *             * r_hat / (4*pi*eps)
   *
   *         This REPLACES the force from F = rho * E where E = -grad(psi)
   * ======================================================================== */

  if (hydro != NULL) {

    /* Sync hydro from device before modification */
    hydro_memcpy(hydro, tdpMemcpyDeviceToHost);

    /* Loop over all local fluid nodes */
    for (int i = 1; i <= nlocal[X]; i++) {
      for (int j = 1; j <= nlocal[Y]; j++) {
        for (int k = 1; k <= nlocal[Z]; k++) {

          int index1 = cs_index(psi->cs, i, j, k);

          /* Get charge at this node */
          double rho0_1, rho1_1;
          psi_rho(psi, index1, 0, &rho0_1);
          psi_rho(psi, index1, 1, &rho1_1);
          double q1 = rho0_1 - rho1_1;

          /* Accumulate real-space force on this node */
          double F_real[3] = {0.0, 0.0, 0.0};

          /* ----- Force from subgrid particles ----- */
          for (ic = 0; ic <= ncell[X] + 1; ic++) {
            for (jc = 0; jc <= ncell[Y] + 1; jc++) {
              for (kc = 0; kc <= ncell[Z] + 1; kc++) {

                colloids_info_cell_list_head(cinfo, ic, jc, kc, &pc);

                for (; pc; pc = pc->next) {

                  if (pc->s.bc != COLLOID_BC_SUBGRID) continue;

                  double q_p = pc->s.q0 - pc->s.q1;

                  /* Particle position in local coords */
                  double r_p[3];
                  r_p[X] = pc->s.r[X] - (double)offset[X];
                  r_p[Y] = pc->s.r[Y] - (double)offset[Y];
                  r_p[Z] = pc->s.r[Z] - (double)offset[Z];

                  /* Vector from particle to node */
                  double r_vec[3];
                  r_vec[X] = (double)i - r_p[X];
                  r_vec[Y] = (double)j - r_p[Y];
                  r_vec[Z] = (double)k - r_p[Z];

                  double r = sqrt(r_vec[X]*r_vec[X] + r_vec[Y]*r_vec[Y] + r_vec[Z]*r_vec[Z]);

                  if (r > rcut || r < 1.0e-10) continue;

                  /* Ewald real-space force magnitude */
                  double ar = alpha * r;
                  double ar2 = ar * ar;
                  double r2 = r * r;
                  double r_inv = 1.0 / r;

                  double term1 = erfc(ar) / r2;
                  double term2 = (2.0 * alpha / sqrt(pi)) * exp(-ar2) * r_inv;
                  double F_mag = q1 * q_p * (term1 + term2) / (4.0 * pi * epsilon);

                  F_real[X] += F_mag * r_vec[X] * r_inv;
                  F_real[Y] += F_mag * r_vec[Y] * r_inv;
                  F_real[Z] += F_mag * r_vec[Z] * r_inv;
                }
              }
            }
          }

          /* ----- Force from other lattice nodes within cutoff ----- */
          for (int di = -irc; di <= irc; di++) {
            for (int dj = -irc; dj <= irc; dj++) {
              for (int dk = -irc; dk <= irc; dk++) {

                if (di == 0 && dj == 0 && dk == 0) continue;  /* Skip self */

                int i2 = i + di;
                int j2 = j + dj;
                int k2 = k + dk;

                if (i2 < 1 || i2 > nlocal[X]) continue;
                if (j2 < 1 || j2 > nlocal[Y]) continue;
                if (k2 < 1 || k2 > nlocal[Z]) continue;

                int index2 = cs_index(psi->cs, i2, j2, k2);

                double rho0_2, rho1_2;
                psi_rho(psi, index2, 0, &rho0_2);
                psi_rho(psi, index2, 1, &rho1_2);
                double q2 = rho0_2 - rho1_2;

                if (fabs(q2) < 1.0e-14) continue;

                /* Vector from node2 to node1 */
                double r_vec[3];
                r_vec[X] = (double)(-di);
                r_vec[Y] = (double)(-dj);
                r_vec[Z] = (double)(-dk);

                double r = sqrt((double)(di*di + dj*dj + dk*dk));

                if (r > rcut) continue;

                /* Ewald real-space force magnitude */
                double ar = alpha * r;
                double ar2 = ar * ar;
                double r2 = r * r;
                double r_inv = 1.0 / r;

                double term1 = erfc(ar) / r2;
                double term2 = (2.0 * alpha / sqrt(pi)) * exp(-ar2) * r_inv;
                double F_mag = q1 * q2 * (term1 + term2) / (4.0 * pi * epsilon);

                F_real[X] += F_mag * r_vec[X] * r_inv;
                F_real[Y] += F_mag * r_vec[Y] * r_inv;
                F_real[Z] += F_mag * r_vec[Z] * r_inv;
              }
            }
          }

          /* Apply total real-space force to this fluid node */
          hydro_f_local_add(hydro, index1, F_real);
        }
      }
    }

    /* Sync hydro to device */
    hydro_memcpy(hydro, tdpMemcpyHostToDevice);
  }

  return 0;
}

/*****************************************************************************
 *
 *  psi_fft_pn_ewald_force_potential_full
 *
 *  Complete Ewald sum for potential and forces using traditional Ewald
 *  (Fourier sum over k-vectors + real-space), NOT FFT.
 *
 *  1. POTENTIAL on all lattice nodes:
 *     - Fourier: sum over k-vectors
 *     - Real-space: particle->node + node->node within cutoff
 *
 *  2. FORCE on all lattice nodes (applied to hydro):
 *     - Fourier: sum over k-vectors
 *     - Real-space: particle->node + node->node within cutoff
 *
 *  3. FIELD (Esub) and FORCE (fex) on subgrid particles:
 *     - Fourier: sum over k-vectors
 *     - Real-space: node->particle + particle->particle within cutoff
 *
 *  This replaces subgrid_update_Esub and psi_force calculations.
 *
 *****************************************************************************/

int psi_fft_pn_ewald_force_potential_full(psi_fft_pn_t * pn, colloids_info_t * cinfo,
                                           hydro_t * hydro) {

  int nlocal[3], offset[3], ntotal[3];
  int ncell[3];
  int ic, jc, kc;
  colloid_t * pc;
  psi_t * psi;
  double alpha, rcut, epsilon, beta, eunit;
  int irc;
  int nk[3], nkmax;
  double ltot[3];
  PI_DOUBLE(pi);

  assert(pn);
  assert(cinfo);

  psi = pn->psi;
  epsilon = pn->epsilon;
  beta = pn->beta;
  psi_unit_charge(psi, &eunit);

  /* Get Ewald parameters from FFT solver */
  if (pn->fft_solver == NULL ||
      pn->fft_solver->influence_type != PSI_FFT_INFLUENCE_EWALD) {
    pe_info(psi->pe, "WARNING: psi_fft_pn_ewald_force_potential_full needs Ewald solver\n");
    return 0;
  }

  alpha = pn->fft_solver->ewald_alpha;
  rcut = pn->fft_solver->ewald_rcut;
  irc = (int)ceil(rcut);

  cs_nlocal(cinfo->cs, nlocal);
  cs_nlocal_offset(cinfo->cs, offset);
  cs_ntotal(cinfo->cs, ntotal);
  cs_ltot(cinfo->cs, ltot);
  colloids_info_ncell(cinfo, ncell);

  /* Number of k-vectors: k_max such that exp(-k^2/4*alpha^2) is small */
  for (int d = 0; d < 3; d++) {
    nk[d] = (int)(alpha * ltot[d] / pi) + 1;
    if (nk[d] > ntotal[d]/2) nk[d] = ntotal[d]/2;
  }
  nkmax = nk[X];
  if (nk[Y] > nkmax) nkmax = nk[Y];
  if (nk[Z] > nkmax) nkmax = nk[Z];

  /* Allocate sin/cos tables */
  double * coskr_x = (double *)malloc((2*nkmax + 1) * sizeof(double));
  double * sinkr_x = (double *)malloc((2*nkmax + 1) * sizeof(double));
  double * coskr_y = (double *)malloc((2*nkmax + 1) * sizeof(double));
  double * sinkr_y = (double *)malloc((2*nkmax + 1) * sizeof(double));
  double * coskr_z = (double *)malloc((2*nkmax + 1) * sizeof(double));
  double * sinkr_z = (double *)malloc((2*nkmax + 1) * sizeof(double));

  /* Count k-vectors */
  int nktot = 0;
  for (int kz = 0; kz <= nk[Z]; kz++) {
    for (int ky = -nk[Y]; ky <= nk[Y]; ky++) {
      for (int kx = -nk[X]; kx <= nk[X]; kx++) {
        if (kx == 0 && ky == 0 && kz == 0) continue;
        if (kz == 0 && (ky < 0 || (ky == 0 && kx < 0))) continue;
        nktot++;
      }
    }
  }

  /* Allocate structure factor and k-vector arrays (using Kahan sums) */
  kahan_t * Sk_cos = (kahan_t *)malloc(nktot * sizeof(kahan_t));
  kahan_t * Sk_sin = (kahan_t *)malloc(nktot * sizeof(kahan_t));
  for (int n = 0; n < nktot; n++) {
    Sk_cos[n] = kahan_zero();
    Sk_sin[n] = kahan_zero();
  }
  double * kvec = (double *)malloc(3 * nktot * sizeof(double));
  double * Gk_arr = (double *)malloc(nktot * sizeof(double));

  /* Precompute k-vectors and Ewald Green function */
  int kn = 0;
  for (int kz = 0; kz <= nk[Z]; kz++) {
    for (int ky = -nk[Y]; ky <= nk[Y]; ky++) {
      for (int kx = -nk[X]; kx <= nk[X]; kx++) {
        if (kx == 0 && ky == 0 && kz == 0) continue;
        if (kz == 0 && (ky < 0 || (ky == 0 && kx < 0))) continue;

        double kv[3];
        kv[X] = 2.0 * pi * kx / ltot[X];
        kv[Y] = 2.0 * pi * ky / ltot[Y];
        kv[Z] = 2.0 * pi * kz / ltot[Z];
        double k2 = kv[X]*kv[X] + kv[Y]*kv[Y] + kv[Z]*kv[Z];

        kvec[3*kn + 0] = kv[X];
        kvec[3*kn + 1] = kv[Y];
        kvec[3*kn + 2] = kv[Z];

        double V = ltot[X] * ltot[Y] * ltot[Z];
        Gk_arr[kn] = 4.0 * pi * exp(-k2 / (4.0 * alpha * alpha)) / (k2 * V * epsilon);
        kn++;
      }
    }
  }

  /* Macro to compute sin/cos tables for position (rx, ry, rz) */
  #define COMPUTE_SINCOS_TABLES(rx, ry, rz) do { \
    double dx = 2.0 * pi * (rx) / ltot[X]; \
    double dy = 2.0 * pi * (ry) / ltot[Y]; \
    double dz = 2.0 * pi * (rz) / ltot[Z]; \
    coskr_x[nkmax] = 1.0; sinkr_x[nkmax] = 0.0; \
    coskr_y[nkmax] = 1.0; sinkr_y[nkmax] = 0.0; \
    coskr_z[nkmax] = 1.0; sinkr_z[nkmax] = 0.0; \
    coskr_x[nkmax+1] = cos(dx); sinkr_x[nkmax+1] = sin(dx); \
    coskr_y[nkmax+1] = cos(dy); sinkr_y[nkmax+1] = sin(dy); \
    coskr_z[nkmax+1] = cos(dz); sinkr_z[nkmax+1] = sin(dz); \
    coskr_x[nkmax-1] = coskr_x[nkmax+1]; sinkr_x[nkmax-1] = -sinkr_x[nkmax+1]; \
    coskr_y[nkmax-1] = coskr_y[nkmax+1]; sinkr_y[nkmax-1] = -sinkr_y[nkmax+1]; \
    coskr_z[nkmax-1] = coskr_z[nkmax+1]; sinkr_z[nkmax-1] = -sinkr_z[nkmax+1]; \
    for (int m = 2; m <= nkmax; m++) { \
      coskr_x[nkmax+m] = coskr_x[nkmax+m-1]*coskr_x[nkmax+1] - sinkr_x[nkmax+m-1]*sinkr_x[nkmax+1]; \
      sinkr_x[nkmax+m] = sinkr_x[nkmax+m-1]*coskr_x[nkmax+1] + coskr_x[nkmax+m-1]*sinkr_x[nkmax+1]; \
      coskr_y[nkmax+m] = coskr_y[nkmax+m-1]*coskr_y[nkmax+1] - sinkr_y[nkmax+m-1]*sinkr_y[nkmax+1]; \
      sinkr_y[nkmax+m] = sinkr_y[nkmax+m-1]*coskr_y[nkmax+1] + coskr_y[nkmax+m-1]*sinkr_y[nkmax+1]; \
      coskr_z[nkmax+m] = coskr_z[nkmax+m-1]*coskr_z[nkmax+1] - sinkr_z[nkmax+m-1]*sinkr_z[nkmax+1]; \
      sinkr_z[nkmax+m] = sinkr_z[nkmax+m-1]*coskr_z[nkmax+1] + coskr_z[nkmax+m-1]*sinkr_z[nkmax+1]; \
      coskr_x[nkmax-m] = coskr_x[nkmax+m]; sinkr_x[nkmax-m] = -sinkr_x[nkmax+m]; \
      coskr_y[nkmax-m] = coskr_y[nkmax+m]; sinkr_y[nkmax-m] = -sinkr_y[nkmax+m]; \
      coskr_z[nkmax-m] = coskr_z[nkmax+m]; sinkr_z[nkmax-m] = -sinkr_z[nkmax+m]; \
    } \
  } while(0)

  /* ========================================================================
   * Compute structure factors S(k) from all charges (particles + nodes)
   * ======================================================================== */

  /* Structure factors from subgrid particles */
  for (ic = 0; ic <= ncell[X] + 1; ic++) {
    for (jc = 0; jc <= ncell[Y] + 1; jc++) {
      for (kc = 0; kc <= ncell[Z] + 1; kc++) {
        colloids_info_cell_list_head(cinfo, ic, jc, kc, &pc);
        for (; pc; pc = pc->next) {
          if (pc->s.bc != COLLOID_BC_SUBGRID) continue;
          double q_p = pc->s.q0 - pc->s.q1;

          COMPUTE_SINCOS_TABLES(pc->s.r[X], pc->s.r[Y], pc->s.r[Z]);

          kn = 0;
          for (int kzz = 0; kzz <= nk[Z]; kzz++) {
            for (int kyy = -nk[Y]; kyy <= nk[Y]; kyy++) {
              for (int kxx = -nk[X]; kxx <= nk[X]; kxx++) {
                if (kxx == 0 && kyy == 0 && kzz == 0) continue;
                if (kzz == 0 && (kyy < 0 || (kyy == 0 && kxx < 0))) continue;

                double cos_kr = coskr_x[nkmax+kxx]*coskr_y[nkmax+kyy]*coskr_z[nkmax+kzz]
                              - coskr_x[nkmax+kxx]*sinkr_y[nkmax+kyy]*sinkr_z[nkmax+kzz]
                              - sinkr_x[nkmax+kxx]*coskr_y[nkmax+kyy]*sinkr_z[nkmax+kzz]
                              - sinkr_x[nkmax+kxx]*sinkr_y[nkmax+kyy]*coskr_z[nkmax+kzz];
                double sin_kr = sinkr_x[nkmax+kxx]*coskr_y[nkmax+kyy]*coskr_z[nkmax+kzz]
                              + coskr_x[nkmax+kxx]*sinkr_y[nkmax+kyy]*coskr_z[nkmax+kzz]
                              + coskr_x[nkmax+kxx]*coskr_y[nkmax+kyy]*sinkr_z[nkmax+kzz]
                              - sinkr_x[nkmax+kxx]*sinkr_y[nkmax+kyy]*sinkr_z[nkmax+kzz];

                kahan_add_double(&Sk_cos[kn], q_p * cos_kr);
                kahan_add_double(&Sk_sin[kn], q_p * sin_kr);
                kn++;
              }
            }
          }
        }
      }
    }
  }

  /* Structure factors from lattice nodes */
  for (int i = 1; i <= nlocal[X]; i++) {
    for (int j = 1; j <= nlocal[Y]; j++) {
      for (int k = 1; k <= nlocal[Z]; k++) {
        int index = cs_index(psi->cs, i, j, k);
        double rho0, rho1;
        psi_rho(psi, index, 0, &rho0);
        psi_rho(psi, index, 1, &rho1);
        double q_node = rho0 - rho1;
        if (fabs(q_node) < 1.0e-14) continue;

        double r_global[3];
        r_global[X] = (double)(offset[X] + i);
        r_global[Y] = (double)(offset[Y] + j);
        r_global[Z] = (double)(offset[Z] + k);

        COMPUTE_SINCOS_TABLES(r_global[X], r_global[Y], r_global[Z]);

        kn = 0;
        for (int kzz = 0; kzz <= nk[Z]; kzz++) {
          for (int kyy = -nk[Y]; kyy <= nk[Y]; kyy++) {
            for (int kxx = -nk[X]; kxx <= nk[X]; kxx++) {
              if (kxx == 0 && kyy == 0 && kzz == 0) continue;
              if (kzz == 0 && (kyy < 0 || (kyy == 0 && kxx < 0))) continue;

              double cos_kr = coskr_x[nkmax+kxx]*coskr_y[nkmax+kyy]*coskr_z[nkmax+kzz]
                            - coskr_x[nkmax+kxx]*sinkr_y[nkmax+kyy]*sinkr_z[nkmax+kzz]
                            - sinkr_x[nkmax+kxx]*coskr_y[nkmax+kyy]*sinkr_z[nkmax+kzz]
                            - sinkr_x[nkmax+kxx]*sinkr_y[nkmax+kyy]*coskr_z[nkmax+kzz];
              double sin_kr = sinkr_x[nkmax+kxx]*coskr_y[nkmax+kyy]*coskr_z[nkmax+kzz]
                            + coskr_x[nkmax+kxx]*sinkr_y[nkmax+kyy]*coskr_z[nkmax+kzz]
                            + coskr_x[nkmax+kxx]*coskr_y[nkmax+kyy]*sinkr_z[nkmax+kzz]
                            - sinkr_x[nkmax+kxx]*sinkr_y[nkmax+kyy]*sinkr_z[nkmax+kzz];

              kahan_add_double(&Sk_cos[kn], q_node * cos_kr);
              kahan_add_double(&Sk_sin[kn], q_node * sin_kr);
              kn++;
            }
          }
        }
      }
    }
  }

  /* MPI reduce structure factors using Kahan-aware reduction */
  MPI_Comm comm;
  cs_cart_comm(psi->cs, &comm);
  MPI_Datatype kahan_dt;
  MPI_Op kahan_op;
  kahan_mpi_datatype(&kahan_dt);
  kahan_mpi_op_sum(&kahan_op);
  MPI_Allreduce(MPI_IN_PLACE, Sk_cos, nktot, kahan_dt, kahan_op, comm);
  MPI_Allreduce(MPI_IN_PLACE, Sk_sin, nktot, kahan_dt, kahan_op, comm);
  MPI_Type_free(&kahan_dt);
  MPI_Op_free(&kahan_op);

  /* ========================================================================
   * PART 1: Potential on lattice nodes
   * ======================================================================== */

  for (int i = 1; i <= nlocal[X]; i++) {
    for (int j = 1; j <= nlocal[Y]; j++) {
      for (int k = 1; k <= nlocal[Z]; k++) {
        int index = cs_index(psi->cs, i, j, k);
        double r_global[3];
        r_global[X] = (double)(offset[X] + i);
        r_global[Y] = (double)(offset[Y] + j);
        r_global[Z] = (double)(offset[Z] + k);

        /* Fourier part (using Kahan summation) */
        kahan_t phi_fourier_k = kahan_zero();
        COMPUTE_SINCOS_TABLES(r_global[X], r_global[Y], r_global[Z]);

        kn = 0;
        for (int kzz = 0; kzz <= nk[Z]; kzz++) {
          for (int kyy = -nk[Y]; kyy <= nk[Y]; kyy++) {
            for (int kxx = -nk[X]; kxx <= nk[X]; kxx++) {
              if (kxx == 0 && kyy == 0 && kzz == 0) continue;
              if (kzz == 0 && (kyy < 0 || (kyy == 0 && kxx < 0))) continue;

              double cos_kr = coskr_x[nkmax+kxx]*coskr_y[nkmax+kyy]*coskr_z[nkmax+kzz]
                            - coskr_x[nkmax+kxx]*sinkr_y[nkmax+kyy]*sinkr_z[nkmax+kzz]
                            - sinkr_x[nkmax+kxx]*coskr_y[nkmax+kyy]*sinkr_z[nkmax+kzz]
                            - sinkr_x[nkmax+kxx]*sinkr_y[nkmax+kyy]*coskr_z[nkmax+kzz];
              double sin_kr = sinkr_x[nkmax+kxx]*coskr_y[nkmax+kyy]*coskr_z[nkmax+kzz]
                            + coskr_x[nkmax+kxx]*sinkr_y[nkmax+kyy]*coskr_z[nkmax+kzz]
                            + coskr_x[nkmax+kxx]*coskr_y[nkmax+kyy]*sinkr_z[nkmax+kzz]
                            - sinkr_x[nkmax+kxx]*sinkr_y[nkmax+kyy]*sinkr_z[nkmax+kzz];

              double factor = (kzz > 0) ? 2.0 : 1.0;
              double Sk_cos_val = kahan_sum(&Sk_cos[kn]);
              double Sk_sin_val = kahan_sum(&Sk_sin[kn]);
              kahan_add_double(&phi_fourier_k, factor * Gk_arr[kn] * (Sk_cos_val*cos_kr + Sk_sin_val*sin_kr));
              kn++;
            }
          }
        }
        double phi_fourier = kahan_sum(&phi_fourier_k);

        /* Real-space from particles */
        double phi_real = 0.0;
        double r_local[3] = {(double)i, (double)j, (double)k};

        for (ic = 0; ic <= ncell[X] + 1; ic++) {
          for (jc = 0; jc <= ncell[Y] + 1; jc++) {
            for (kc = 0; kc <= ncell[Z] + 1; kc++) {
              colloids_info_cell_list_head(cinfo, ic, jc, kc, &pc);
              for (; pc; pc = pc->next) {
                if (pc->s.bc != COLLOID_BC_SUBGRID) continue;
                double q_p = pc->s.q0 - pc->s.q1;
                double r_p[3];
                r_p[X] = pc->s.r[X] - (double)offset[X];
                r_p[Y] = pc->s.r[Y] - (double)offset[Y];
                r_p[Z] = pc->s.r[Z] - (double)offset[Z];
                double dr[3] = {r_local[X]-r_p[X], r_local[Y]-r_p[Y], r_local[Z]-r_p[Z]};
                double r = sqrt(dr[X]*dr[X] + dr[Y]*dr[Y] + dr[Z]*dr[Z]);
                if (r < rcut && r > 1.0e-10) {
                  phi_real += q_p * erfc(alpha*r) / (4.0*pi*epsilon*r);
                }
              }
            }
          }
        }

        /* Real-space from other nodes */
        for (int di = -irc; di <= irc; di++) {
          for (int dj = -irc; dj <= irc; dj++) {
            for (int dk = -irc; dk <= irc; dk++) {
              if (di == 0 && dj == 0 && dk == 0) continue;
              int i2 = i + di, j2 = j + dj, k2 = k + dk;
              if (i2 < 1 || i2 > nlocal[X]) continue;
              if (j2 < 1 || j2 > nlocal[Y]) continue;
              if (k2 < 1 || k2 > nlocal[Z]) continue;
              int index2 = cs_index(psi->cs, i2, j2, k2);
              double rho0_2, rho1_2;
              psi_rho(psi, index2, 0, &rho0_2);
              psi_rho(psi, index2, 1, &rho1_2);
              double q2 = rho0_2 - rho1_2;
              if (fabs(q2) < 1.0e-14) continue;
              double r = sqrt((double)(di*di + dj*dj + dk*dk));
              if (r < rcut) {
                phi_real += q2 * erfc(alpha*r) / (4.0*pi*epsilon*r);
              }
            }
          }
        }

        /* Store: psi* = beta * eunit * phi */
        psi->psi->data[addr_rank0(psi->psi->nsites, index)] = beta * eunit * (phi_fourier + phi_real);
      }
    }
  }

  /* ========================================================================
   * PART 2: Force on lattice nodes -> hydro
   * ======================================================================== */

  if (hydro != NULL) {
    // /* Sync hydro from device before modification */
    // hydro_memcpy(hydro, tdpMemcpyDeviceToHost);

    for (int i = 1; i <= nlocal[X]; i++) {
      for (int j = 1; j <= nlocal[Y]; j++) {
        for (int k = 1; k <= nlocal[Z]; k++) {
          int index = cs_index(psi->cs, i, j, k);
          double rho0, rho1;
          psi_rho(psi, index, 0, &rho0);
          psi_rho(psi, index, 1, &rho1);
          double q1 = rho0 - rho1;

          double r_global[3];
          r_global[X] = (double)(offset[X] + i);
          r_global[Y] = (double)(offset[Y] + j);
          r_global[Z] = (double)(offset[Z] + k);

          /* Fourier part: F = -q * grad(phi) = q * sum_k G(k)*k*Im[S*exp(-ikr)] (using Kahan) */
          kahan_t F_fourier_k[3] = {kahan_zero(), kahan_zero(), kahan_zero()};
          COMPUTE_SINCOS_TABLES(r_global[X], r_global[Y], r_global[Z]);

          kn = 0;
          for (int kzz = 0; kzz <= nk[Z]; kzz++) {
            for (int kyy = -nk[Y]; kyy <= nk[Y]; kyy++) {
              for (int kxx = -nk[X]; kxx <= nk[X]; kxx++) {
                if (kxx == 0 && kyy == 0 && kzz == 0) continue;
                if (kzz == 0 && (kyy < 0 || (kyy == 0 && kxx < 0))) continue;

                double cos_kr = coskr_x[nkmax+kxx]*coskr_y[nkmax+kyy]*coskr_z[nkmax+kzz]
                              - coskr_x[nkmax+kxx]*sinkr_y[nkmax+kyy]*sinkr_z[nkmax+kzz]
                              - sinkr_x[nkmax+kxx]*coskr_y[nkmax+kyy]*sinkr_z[nkmax+kzz]
                              - sinkr_x[nkmax+kxx]*sinkr_y[nkmax+kyy]*coskr_z[nkmax+kzz];
                double sin_kr = sinkr_x[nkmax+kxx]*coskr_y[nkmax+kyy]*coskr_z[nkmax+kzz]
                              + coskr_x[nkmax+kxx]*sinkr_y[nkmax+kyy]*coskr_z[nkmax+kzz]
                              + coskr_x[nkmax+kxx]*coskr_y[nkmax+kyy]*sinkr_z[nkmax+kzz]
                              - sinkr_x[nkmax+kxx]*sinkr_y[nkmax+kyy]*sinkr_z[nkmax+kzz];

                double factor = (kzz > 0) ? 2.0 : 1.0;
                double Sk_cos_val = kahan_sum(&Sk_cos[kn]);
                double Sk_sin_val = kahan_sum(&Sk_sin[kn]);
                double im_part = Sk_sin_val*cos_kr - Sk_cos_val*sin_kr;

                kahan_add_double(&F_fourier_k[X], factor * q1 * Gk_arr[kn] * kvec[3*kn+0] * im_part);
                kahan_add_double(&F_fourier_k[Y], factor * q1 * Gk_arr[kn] * kvec[3*kn+1] * im_part);
                kahan_add_double(&F_fourier_k[Z], factor * q1 * Gk_arr[kn] * kvec[3*kn+2] * im_part);
                kn++;
              }
            }
          }
          double F_total[3] = {kahan_sum(&F_fourier_k[X]), kahan_sum(&F_fourier_k[Y]), kahan_sum(&F_fourier_k[Z])};

          /* Real-space from particles */
          double r_local[3] = {(double)i, (double)j, (double)k};
          for (ic = 0; ic <= ncell[X] + 1; ic++) {
            for (jc = 0; jc <= ncell[Y] + 1; jc++) {
              for (kc = 0; kc <= ncell[Z] + 1; kc++) {
                colloids_info_cell_list_head(cinfo, ic, jc, kc, &pc);
                for (; pc; pc = pc->next) {
                  if (pc->s.bc != COLLOID_BC_SUBGRID) continue;
                  double q_p = pc->s.q0 - pc->s.q1;
                  double r_p[3];
                  r_p[X] = pc->s.r[X] - (double)offset[X];
                  r_p[Y] = pc->s.r[Y] - (double)offset[Y];
                  r_p[Z] = pc->s.r[Z] - (double)offset[Z];
                  double dr[3] = {r_local[X]-r_p[X], r_local[Y]-r_p[Y], r_local[Z]-r_p[Z]};
                  double r = sqrt(dr[X]*dr[X] + dr[Y]*dr[Y] + dr[Z]*dr[Z]);
                  if (r < rcut && r > 1.0e-10) {
                    double ar = alpha*r;
                    double r_inv = 1.0/r;
                    double term1 = erfc(ar)/(r*r);
                    double term2 = (2.0*alpha/sqrt(pi))*exp(-ar*ar)*r_inv;
                    double F_mag = q1*q_p*(term1+term2)/(4.0*pi*epsilon);
                    F_total[X] += F_mag*dr[X]*r_inv;
                    F_total[Y] += F_mag*dr[Y]*r_inv;
                    F_total[Z] += F_mag*dr[Z]*r_inv;
                  }
                }
              }
            }
          }

          /* Real-space from other nodes */
          for (int di = -irc; di <= irc; di++) {
            for (int dj = -irc; dj <= irc; dj++) {
              for (int dk = -irc; dk <= irc; dk++) {
                if (di == 0 && dj == 0 && dk == 0) continue;
                int i2 = i+di, j2 = j+dj, k2 = k+dk;
                if (i2 < 1 || i2 > nlocal[X]) continue;
                if (j2 < 1 || j2 > nlocal[Y]) continue;
                if (k2 < 1 || k2 > nlocal[Z]) continue;
                int index2 = cs_index(psi->cs, i2, j2, k2);
                double rho0_2, rho1_2;
                psi_rho(psi, index2, 0, &rho0_2);
                psi_rho(psi, index2, 1, &rho1_2);
                double q2 = rho0_2 - rho1_2;
                if (fabs(q2) < 1.0e-14) continue;
                double r = sqrt((double)(di*di + dj*dj + dk*dk));
                if (r < rcut) {
                  double ar = alpha*r;
                  double r_inv = 1.0/r;
                  double term1 = erfc(ar)/(r*r);
                  double term2 = (2.0*alpha/sqrt(pi))*exp(-ar*ar)*r_inv;
                  double F_mag = q1*q2*(term1+term2)/(4.0*pi*epsilon);
                  F_total[X] += F_mag*(double)(-di)*r_inv;
                  F_total[Y] += F_mag*(double)(-dj)*r_inv;
                  F_total[Z] += F_mag*(double)(-dk)*r_inv;
                }
              }
            }
          }

          hydro_f_local_add(hydro, index, F_total);
        }
      }
    }
  }

  /* ========================================================================
   * PART 3: Field (Esub) and Force (fex) on subgrid particles
   * ======================================================================== */

  for (ic = 0; ic <= ncell[X] + 1; ic++) {
    for (jc = 0; jc <= ncell[Y] + 1; jc++) {
      for (kc = 0; kc <= ncell[Z] + 1; kc++) {
        colloids_info_cell_list_head(cinfo, ic, jc, kc, &pc);
        for (; pc; pc = pc->next) {
          if (pc->s.bc != COLLOID_BC_SUBGRID) continue;
          double q_p = pc->s.q0 - pc->s.q1;

          /* Fourier part (using Kahan summation) */
          kahan_t E_fourier_k[3] = {kahan_zero(), kahan_zero(), kahan_zero()};
          COMPUTE_SINCOS_TABLES(pc->s.r[X], pc->s.r[Y], pc->s.r[Z]);

          kn = 0;
          for (int kzz = 0; kzz <= nk[Z]; kzz++) {
            for (int kyy = -nk[Y]; kyy <= nk[Y]; kyy++) {
              for (int kxx = -nk[X]; kxx <= nk[X]; kxx++) {
                if (kxx == 0 && kyy == 0 && kzz == 0) continue;
                if (kzz == 0 && (kyy < 0 || (kyy == 0 && kxx < 0))) continue;

                double cos_kr = coskr_x[nkmax+kxx]*coskr_y[nkmax+kyy]*coskr_z[nkmax+kzz]
                              - coskr_x[nkmax+kxx]*sinkr_y[nkmax+kyy]*sinkr_z[nkmax+kzz]
                              - sinkr_x[nkmax+kxx]*coskr_y[nkmax+kyy]*sinkr_z[nkmax+kzz]
                              - sinkr_x[nkmax+kxx]*sinkr_y[nkmax+kyy]*coskr_z[nkmax+kzz];
                double sin_kr = sinkr_x[nkmax+kxx]*coskr_y[nkmax+kyy]*coskr_z[nkmax+kzz]
                              + coskr_x[nkmax+kxx]*sinkr_y[nkmax+kyy]*coskr_z[nkmax+kzz]
                              + coskr_x[nkmax+kxx]*coskr_y[nkmax+kyy]*sinkr_z[nkmax+kzz]
                              - sinkr_x[nkmax+kxx]*sinkr_y[nkmax+kyy]*sinkr_z[nkmax+kzz];

                double factor = (kzz > 0) ? 2.0 : 1.0;
                double Sk_cos_val = kahan_sum(&Sk_cos[kn]);
                double Sk_sin_val = kahan_sum(&Sk_sin[kn]);
                double im_part = Sk_sin_val*cos_kr - Sk_cos_val*sin_kr;

                /* Esub = beta * eunit * E_physical */
                kahan_add_double(&E_fourier_k[X], factor * beta * eunit * Gk_arr[kn] * kvec[3*kn+0] * im_part);
                kahan_add_double(&E_fourier_k[Y], factor * beta * eunit * Gk_arr[kn] * kvec[3*kn+1] * im_part);
                kahan_add_double(&E_fourier_k[Z], factor * beta * eunit * Gk_arr[kn] * kvec[3*kn+2] * im_part);
                kn++;
              }
            }
          }
          double E_total[3] = {kahan_sum(&E_fourier_k[X]), kahan_sum(&E_fourier_k[Y]), kahan_sum(&E_fourier_k[Z])};

          /* Real-space from lattice nodes */
          double r_p[3];
          r_p[X] = pc->s.r[X] - (double)offset[X];
          r_p[Y] = pc->s.r[Y] - (double)offset[Y];
          r_p[Z] = pc->s.r[Z] - (double)offset[Z];
          int i0 = (int)floor(r_p[X]);
          int j0 = (int)floor(r_p[Y]);
          int k0 = (int)floor(r_p[Z]);

          for (int di = -irc; di <= irc+1; di++) {
            for (int dj = -irc; dj <= irc+1; dj++) {
              for (int dk = -irc; dk <= irc+1; dk++) {
                int ni = i0+di, nj = j0+dj, nk_idx = k0+dk;
                if (ni < 1 || ni > nlocal[X]) continue;
                if (nj < 1 || nj > nlocal[Y]) continue;
                if (nk_idx < 1 || nk_idx > nlocal[Z]) continue;
                int index = cs_index(psi->cs, ni, nj, nk_idx);
                double rho0, rho1;
                psi_rho(psi, index, 0, &rho0);
                psi_rho(psi, index, 1, &rho1);
                double q_node = rho0 - rho1;
                if (fabs(q_node) < 1.0e-14) continue;
                double dr[3] = {r_p[X]-(double)ni, r_p[Y]-(double)nj, r_p[Z]-(double)nk_idx};
                double r = sqrt(dr[X]*dr[X] + dr[Y]*dr[Y] + dr[Z]*dr[Z]);
                if (r < rcut && r > 1.0e-10) {
                  double ar = alpha*r;
                  double r_inv = 1.0/r;
                  double term1 = erfc(ar)/(r*r);
                  double term2 = (2.0*alpha/sqrt(pi))*exp(-ar*ar)*r_inv;
                  double E_mag = beta*eunit*(term1+term2)/(4.0*pi*epsilon);
                  E_total[X] += q_node*E_mag*dr[X]*r_inv;
                  E_total[Y] += q_node*E_mag*dr[Y]*r_inv;
                  E_total[Z] += q_node*E_mag*dr[Z]*r_inv;
                }
              }
            }
          }

          /* Real-space from other particles */
          for (int ic2 = 0; ic2 <= ncell[X]+1; ic2++) {
            for (int jc2 = 0; jc2 <= ncell[Y]+1; jc2++) {
              for (int kc2 = 0; kc2 <= ncell[Z]+1; kc2++) {
                colloid_t * pc2;
                colloids_info_cell_list_head(cinfo, ic2, jc2, kc2, &pc2);
                for (; pc2; pc2 = pc2->next) {
                  if (pc2 == pc) continue;
                  if (pc2->s.bc != COLLOID_BC_SUBGRID) continue;
                  double q_p2 = pc2->s.q0 - pc2->s.q1;
                  double r_p2[3];
                  r_p2[X] = pc2->s.r[X] - (double)offset[X];
                  r_p2[Y] = pc2->s.r[Y] - (double)offset[Y];
                  r_p2[Z] = pc2->s.r[Z] - (double)offset[Z];
                  double dr[3] = {r_p[X]-r_p2[X], r_p[Y]-r_p2[Y], r_p[Z]-r_p2[Z]};
                  double r = sqrt(dr[X]*dr[X] + dr[Y]*dr[Y] + dr[Z]*dr[Z]);
                  if (r < rcut && r > 1.0e-10) {
                    double ar = alpha*r;
                    double r_inv = 1.0/r;
                    double term1 = erfc(ar)/(r*r);
                    double term2 = (2.0*alpha/sqrt(pi))*exp(-ar*ar)*r_inv;
                    double E_mag = beta*eunit*(term1+term2)/(4.0*pi*epsilon);
                    E_total[X] += q_p2*E_mag*dr[X]*r_inv;
                    E_total[Y] += q_p2*E_mag*dr[Y]*r_inv;
                    E_total[Z] += q_p2*E_mag*dr[Z]*r_inv;
                  }
                }
              }
            }
          }

          /* Store field and force */
          pc->Esub[X] = E_total[X];
          pc->Esub[Y] = E_total[Y];
          pc->Esub[Z] = E_total[Z];

          double kt = 1.0/beta;
          pc->fex[X] = q_p * E_total[X] * kt / eunit;
          pc->fex[Y] = q_p * E_total[Y] * kt / eunit;
          pc->fex[Z] = q_p * E_total[Z] * kt / eunit;
        }
      }
    }
  }

  // if (hydro != NULL) {
  //   hydro_memcpy(hydro, tdpMemcpyHostToDevice);
  // }

  free(coskr_x); free(sinkr_x);
  free(coskr_y); free(sinkr_y);
  free(coskr_z); free(sinkr_z);
  free(Sk_cos); free(Sk_sin);
  free(kvec); free(Gk_arr);

  #undef COMPUTE_SINCOS_TABLES

  return 0;
}
