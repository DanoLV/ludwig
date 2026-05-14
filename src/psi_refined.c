/*****************************************************************************
 *
 *  psi_refined.c
 *
 *  Refined grid Poisson solver implementation.
 *
 *  Algorithm overview:
 *  1. Create refined grid with spacing dx_ref = 1.0/rfactor
 *  2. Transfer fluid charges: LB node charges go to coincident refined nodes
 *  3. Distribute particle charges using Peskin delta on refined grid
 *  4. Solve Poisson with FFT on refined grid
 *  5. Transfer potential back to LB nodes at coincident positions
 *  6. Compute electric field on refined grid
 *  7. Interpolate field to particles using Peskin on refined grid
 *
 *  Key insight: Peskin delta function on refined grid has better resolution
 *  for particle charge distribution and field interpolation.
 *
 *  CHANGE INIT - 20260121 Refined grid Poisson solver
 *
 *****************************************************************************/

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "pe.h"
#include "coords.h"
#include "psi_refined.h"

/* Only compile CUDA parts if using nvcc */
#ifdef __NVCC__
#include <cufft.h>
#include <cuda_runtime.h>
#endif

/* Peskin delta function support radius in grid units */
#define PESKIN_RANGE 2.0

/*****************************************************************************
 *
 *  psi_refined_d_peskin
 *
 *  Peskin delta function scaled for refined grid spacing.
 *  The input r is in refined grid units (already scaled by dx_ref).
 *
 *****************************************************************************/

double psi_refined_d_peskin(double r, double dx) {

  double rmod = fabs(r);
  double delta = 0.0;

  /* Peskin 4-point stencil, scaled by 1/dx for normalization */
  if (rmod <= 1.0) {
    delta = 0.125 * (3.0 - 2.0*rmod + sqrt(1.0 + 4.0*rmod - 4.0*rmod*rmod));
  }
  else if (rmod <= 2.0) {
    delta = 0.125 * (5.0 - 2.0*rmod - sqrt(-7.0 + 12.0*rmod - 4.0*rmod*rmod));
  }

  /* Scale by 1/dx for 1D, total 3D weight will be (1/dx)^3 */
  return delta / dx;
}

/*****************************************************************************
 *
 *  psi_refined_create
 *
 *  Create refined grid structure and allocate memory.
 *
 *****************************************************************************/

int psi_refined_create(pe_t * pe, cs_t * cs, psi_t * psi,
                        const psi_refined_options_t * opts,
                        psi_refined_t ** pref) {

  psi_refined_t * obj = NULL;
  int nlocal[3], ntotal[3];
  double ltot[3];

  assert(pe);
  assert(cs);
  assert(psi);
  assert(opts);
  assert(pref);

  obj = (psi_refined_t *) calloc(1, sizeof(psi_refined_t));
  if (obj == NULL) {
    pe_fatal(pe, "psi_refined_create: calloc failed\n");
    return -1;
  }

  obj->pe = pe;
  obj->cs = cs;
  obj->psi_lb = psi;
  obj->is_initialised = 0;

  /* Store refinement factor */
  obj->rfactor = opts->refinement_factor;
  if (obj->rfactor < 1) {
    pe_fatal(pe, "psi_refined: refinement_factor must be >= 1\n");
    free(obj);
    return -1;
  }

  obj->dx_ref = 1.0 / (double) obj->rfactor;

  /* Get LB grid dimensions */
  cs_nlocal(cs, nlocal);
  cs_ntotal(cs, ntotal);
  cs_ltot(cs, ltot);

  obj->nx_lb = nlocal[X];
  obj->ny_lb = nlocal[Y];
  obj->nz_lb = nlocal[Z];
  obj->ntotal_lb[X] = ntotal[X];
  obj->ntotal_lb[Y] = ntotal[Y];
  obj->ntotal_lb[Z] = ntotal[Z];

  /* Compute refined grid dimensions */
  obj->nx_ref = nlocal[X] * obj->rfactor;
  obj->ny_ref = nlocal[Y] * obj->rfactor;
  obj->nz_ref = nlocal[Z] * obj->rfactor;
  obj->ntotal_ref[X] = ntotal[X] * obj->rfactor;
  obj->ntotal_ref[Y] = ntotal[Y] * obj->rfactor;
  obj->ntotal_ref[Z] = ntotal[Z] * obj->rfactor;

  obj->lx = ltot[X];
  obj->ly = ltot[Y];
  obj->lz = ltot[Z];

  /* Physical parameters */
  obj->epsilon = opts->epsilon;
  obj->beta = opts->beta;
  obj->e = opts->e;

  /* Allocate host arrays for refined grid */
  int n_ref = obj->nx_ref * obj->ny_ref * obj->nz_ref;

  obj->rho_ref = (double *) calloc(n_ref, sizeof(double));
  obj->psi_ref = (double *) calloc(n_ref, sizeof(double));
  obj->efield_ref = (double *) calloc(3 * n_ref, sizeof(double));

  if (obj->rho_ref == NULL || obj->psi_ref == NULL || obj->efield_ref == NULL) {
    pe_fatal(pe, "psi_refined_create: calloc for host arrays failed\n");
    if (obj->rho_ref) free(obj->rho_ref);
    if (obj->psi_ref) free(obj->psi_ref);
    if (obj->efield_ref) free(obj->efield_ref);
    free(obj);
    return -1;
  }

#ifdef __NVCC__
  /* Allocate device arrays for FFT */
  int n_complex = obj->nx_ref * obj->ny_ref * (obj->nz_ref / 2 + 1);
  cudaError_t cuda_status;

  cuda_status = cudaMalloc(&obj->rho_real_d, n_ref * sizeof(double));
  if (cuda_status != cudaSuccess) {
    pe_fatal(pe, "psi_refined: cudaMalloc rho_real failed\n");
    free(obj->rho_ref);
    free(obj->psi_ref);
    free(obj->efield_ref);
    free(obj);
    return -1;
  }

  cuda_status = cudaMalloc(&obj->rho_complex_d, n_complex * sizeof(cufftDoubleComplex));
  if (cuda_status != cudaSuccess) {
    pe_fatal(pe, "psi_refined: cudaMalloc rho_complex failed\n");
    cudaFree(obj->rho_real_d);
    free(obj->rho_ref);
    free(obj->psi_ref);
    free(obj->efield_ref);
    free(obj);
    return -1;
  }

  cuda_status = cudaMalloc(&obj->psi_complex_d, n_complex * sizeof(cufftDoubleComplex));
  if (cuda_status != cudaSuccess) {
    pe_fatal(pe, "psi_refined: cudaMalloc psi_complex failed\n");
    cudaFree(obj->rho_real_d);
    cudaFree(obj->rho_complex_d);
    free(obj->rho_ref);
    free(obj->psi_ref);
    free(obj->efield_ref);
    free(obj);
    return -1;
  }

  cuda_status = cudaMalloc(&obj->psi_real_d, n_ref * sizeof(double));
  if (cuda_status != cudaSuccess) {
    pe_fatal(pe, "psi_refined: cudaMalloc psi_real failed\n");
    cudaFree(obj->rho_real_d);
    cudaFree(obj->rho_complex_d);
    cudaFree(obj->psi_complex_d);
    free(obj->rho_ref);
    free(obj->psi_ref);
    free(obj->efield_ref);
    free(obj);
    return -1;
  }

  /* Compute discrete Laplacian eigenvalues for refined grid */
  double * eigenval_h = (double *) malloc(n_complex * sizeof(double));
  if (eigenval_h == NULL) {
    pe_fatal(pe, "psi_refined: malloc eigenval failed\n");
    cudaFree(obj->rho_real_d);
    cudaFree(obj->rho_complex_d);
    cudaFree(obj->psi_complex_d);
    cudaFree(obj->psi_real_d);
    free(obj->rho_ref);
    free(obj->psi_ref);
    free(obj->efield_ref);
    free(obj);
    return -1;
  }

  int nx = obj->nx_ref;
  int ny = obj->ny_ref;
  int nz = obj->nz_ref;
  double dx = obj->dx_ref;

  /* Eigenvalues for 7-point Laplacian stencil:
   * laplacian[i] = (1/dx^2) * (psi[i+1] + psi[i-1] + psi[j+1] + psi[j-1]
   *                          + psi[k+1] + psi[k-1] - 6*psi[i])
   *
   * In Fourier space:
   * lambda(k) = (2/dx^2) * (cos(kx*dx) + cos(ky*dx) + cos(kz*dx) - 3)
   */
  for (int ix = 0; ix < nx; ix++) {
    for (int iy = 0; iy < ny; iy++) {
      for (int iz = 0; iz < nz / 2 + 1; iz++) {
        int idx = ix * ny * (nz / 2 + 1) + iy * (nz / 2 + 1) + iz;

        /* Wave numbers (integer indices) */
        double kx_idx = (ix <= nx/2) ? (double)ix : (double)(ix - nx);
        double ky_idx = (iy <= ny/2) ? (double)iy : (double)(iy - ny);
        double kz_idx = (double)iz;

        /* Wave vector components: k = 2*pi*n/L, where L = N*dx */
        double kx = 2.0 * M_PI * kx_idx / (nx * dx);
        double ky = 2.0 * M_PI * ky_idx / (ny * dx);
        double kz = 2.0 * M_PI * kz_idx / (nz * dx);

        /* Discrete Laplacian eigenvalue */
        double lambda = (2.0 / (dx * dx)) *
                        (cos(kx * dx) + cos(ky * dx) + cos(kz * dx) - 3.0);

        /* Store -lambda (positive for k != 0) */
        eigenval_h[idx] = -lambda;
      }
    }
  }

  cuda_status = cudaMalloc(&obj->laplacian_eigenval_d, n_complex * sizeof(double));
  if (cuda_status != cudaSuccess) {
    pe_fatal(pe, "psi_refined: cudaMalloc eigenval failed\n");
    free(eigenval_h);
    cudaFree(obj->rho_real_d);
    cudaFree(obj->rho_complex_d);
    cudaFree(obj->psi_complex_d);
    cudaFree(obj->psi_real_d);
    free(obj->rho_ref);
    free(obj->psi_ref);
    free(obj->efield_ref);
    free(obj);
    return -1;
  }

  cudaMemcpy(obj->laplacian_eigenval_d, eigenval_h, n_complex * sizeof(double),
             cudaMemcpyHostToDevice);
  free(eigenval_h);

  /* Create cuFFT plans */
  cufftHandle * plan_fwd = (cufftHandle *) malloc(sizeof(cufftHandle));
  cufftHandle * plan_bwd = (cufftHandle *) malloc(sizeof(cufftHandle));
  cufftResult cufft_status;

  cufft_status = cufftPlan3d(plan_fwd, nx, ny, nz, CUFFT_D2Z);
  if (cufft_status != CUFFT_SUCCESS) {
    pe_fatal(pe, "psi_refined: cufftPlan3d forward failed\n");
    cudaFree(obj->rho_real_d);
    cudaFree(obj->rho_complex_d);
    cudaFree(obj->psi_complex_d);
    cudaFree(obj->psi_real_d);
    cudaFree(obj->laplacian_eigenval_d);
    free(plan_fwd);
    free(plan_bwd);
    free(obj->rho_ref);
    free(obj->psi_ref);
    free(obj->efield_ref);
    free(obj);
    return -1;
  }

  cufft_status = cufftPlan3d(plan_bwd, nx, ny, nz, CUFFT_Z2D);
  if (cufft_status != CUFFT_SUCCESS) {
    pe_fatal(pe, "psi_refined: cufftPlan3d backward failed\n");
    cufftDestroy(*plan_fwd);
    cudaFree(obj->rho_real_d);
    cudaFree(obj->rho_complex_d);
    cudaFree(obj->psi_complex_d);
    cudaFree(obj->psi_real_d);
    cudaFree(obj->laplacian_eigenval_d);
    free(plan_fwd);
    free(plan_bwd);
    free(obj->rho_ref);
    free(obj->psi_ref);
    free(obj->efield_ref);
    free(obj);
    return -1;
  }

  obj->plan_forward = plan_fwd;
  obj->plan_backward = plan_bwd;

#endif /* __NVCC__ */

  obj->is_initialised = 1;
  *pref = obj;

  return 0;
}

/*****************************************************************************
 *
 *  psi_refined_free
 *
 *****************************************************************************/

int psi_refined_free(psi_refined_t ** pref) {

  psi_refined_t * obj = NULL;

  assert(pref);
  obj = *pref;

  if (obj) {
    /* Free host arrays */
    if (obj->rho_ref) free(obj->rho_ref);
    if (obj->psi_ref) free(obj->psi_ref);
    if (obj->efield_ref) free(obj->efield_ref);

#ifdef __NVCC__
    /* Free device arrays */
    if (obj->rho_real_d) cudaFree(obj->rho_real_d);
    if (obj->rho_complex_d) cudaFree(obj->rho_complex_d);
    if (obj->psi_complex_d) cudaFree(obj->psi_complex_d);
    if (obj->psi_real_d) cudaFree(obj->psi_real_d);
    if (obj->laplacian_eigenval_d) cudaFree(obj->laplacian_eigenval_d);

    /* Destroy cuFFT plans */
    if (obj->plan_forward) {
      cufftDestroy(*((cufftHandle *) obj->plan_forward));
      free(obj->plan_forward);
    }
    if (obj->plan_backward) {
      cufftDestroy(*((cufftHandle *) obj->plan_backward));
      free(obj->plan_backward);
    }
#endif

    free(obj);
    *pref = NULL;
  }

  return 0;
}

/*****************************************************************************
 *
 *  psi_refined_reset
 *
 *  Reset charge density and potential arrays to zero.
 *
 *****************************************************************************/

int psi_refined_reset(psi_refined_t * pref) {

  int n_ref;

  assert(pref);

  n_ref = pref->nx_ref * pref->ny_ref * pref->nz_ref;

  memset(pref->rho_ref, 0, n_ref * sizeof(double));
  memset(pref->psi_ref, 0, n_ref * sizeof(double));
  memset(pref->efield_ref, 0, 3 * n_ref * sizeof(double));

  return 0;
}

/*****************************************************************************
 *
 *  psi_refined_lb_to_ref_index
 *
 *  Convert LB grid indices (i_lb, j_lb, k_lb) starting from 1 to
 *  refined grid linear index.
 *
 *  LB node at (i_lb, j_lb, k_lb) corresponds to refined node at:
 *    i_ref = (i_lb - 1) * rfactor
 *    j_ref = (j_lb - 1) * rfactor
 *    k_ref = (k_lb - 1) * rfactor
 *
 *****************************************************************************/

int psi_refined_lb_to_ref_index(psi_refined_t * pref,
                                 int i_lb, int j_lb, int k_lb) {

  int i_ref, j_ref, k_ref;
  int idx;

  assert(pref);

  /* LB indices start at 1, refined indices start at 0 */
  i_ref = (i_lb - 1) * pref->rfactor;
  j_ref = (j_lb - 1) * pref->rfactor;
  k_ref = (k_lb - 1) * pref->rfactor;

  /* Linear index in row-major order */
  idx = i_ref * (pref->ny_ref * pref->nz_ref) + j_ref * pref->nz_ref + k_ref;

  return idx;
}

/*****************************************************************************
 *
 *  psi_refined_lb_to_ref_position
 *
 *  Convert position from LB units to refined grid units.
 *  Position in LB goes from 1 to N, in refined from 0 to N*rfactor-1.
 *
 *****************************************************************************/

void psi_refined_lb_to_ref_position(psi_refined_t * pref,
                                     const double r_lb[3],
                                     double r_ref[3]) {

  assert(pref);

  /* LB position r_lb in [1, N] maps to refined position in [0, N*rfactor)
   * r_ref = (r_lb - 1) * rfactor */
  r_ref[X] = (r_lb[X] - 1.0) * pref->rfactor;
  r_ref[Y] = (r_lb[Y] - 1.0) * pref->rfactor;
  r_ref[Z] = (r_lb[Z] - 1.0) * pref->rfactor;
}

/*****************************************************************************
 *
 *  psi_refined_transfer_charges_from_lb
 *
 *  Transfer fluid charge densities from LB grid to refined grid.
 *  Only the nodes that coincide (every rfactor-th node in refined grid)
 *  receive the LB charge values.
 *
 *  The charge density is scaled by rfactor^3 to maintain total charge
 *  since the refined cell volume is (1/rfactor)^3 of LB cell.
 *
 *****************************************************************************/

int psi_refined_transfer_charges_from_lb(psi_refined_t * pref) {

  psi_t * psi = NULL;
  int nlocal[3], nhalo;
  int ic, jc, kc;
  int i_ref, j_ref, k_ref, idx_ref;
  double rho0, rho1, rho_elec;

  assert(pref);
  assert(pref->psi_lb);

  psi = pref->psi_lb;

  cs_nlocal(pref->cs, nlocal);
  cs_nhalo(pref->cs, &nhalo);

  /* Loop over LB interior nodes */
  for (ic = 1; ic <= nlocal[X]; ic++) {
    for (jc = 1; jc <= nlocal[Y]; jc++) {
      for (kc = 1; kc <= nlocal[Z]; kc++) {

        int index_lb = cs_index(pref->cs, ic, jc, kc);

        /* Get charge densities from LB psi */
        psi_rho(psi, index_lb, 0, &rho0);
        psi_rho(psi, index_lb, 1, &rho1);

        /* Net charge density (rho_elec = rho0 - rho1) scaled by beta */
        rho_elec = (rho0 - rho1) * pref->beta;

        /* Refined grid indices (start at 0) */
        i_ref = (ic - 1) * pref->rfactor;
        j_ref = (jc - 1) * pref->rfactor;
        k_ref = (kc - 1) * pref->rfactor;

        /* Linear index in refined grid (row-major) */
        idx_ref = i_ref * (pref->ny_ref * pref->nz_ref) +
                  j_ref * pref->nz_ref + k_ref;

        /* Store charge with volume scaling:
         * Total charge Q = rho_lb * V_lb = rho_ref * V_ref
         * V_lb = 1, V_ref = dx_ref^3 = (1/rfactor)^3
         * So rho_ref = rho_lb * rfactor^3 */
        pref->rho_ref[idx_ref] = rho_elec * pref->rfactor * pref->rfactor * pref->rfactor;
      }
    }
  }

  return 0;
}

/*****************************************************************************
 *
 *  psi_refined_distribute_particle_charges
 *
 *  Distribute particle charges onto refined grid using Peskin delta.
 *  The Peskin function operates in refined grid units.
 *
 *****************************************************************************/

int psi_refined_distribute_particle_charges(psi_refined_t * pref,
                                             colloids_info_t * cinfo) {

  int ic, jc, kc;
  int ncell[3];
  int nlocal_lb[3], offset_lb[3];
  colloid_t * pc = NULL;

  assert(pref);
  assert(cinfo);

  if (cinfo->nsubgrid == 0) return 0;

  cs_nlocal(pref->cs, nlocal_lb);
  cs_nlocal_offset(pref->cs, offset_lb);
  colloids_info_ncell(cinfo, ncell);

  int rfactor = pref->rfactor;
  double dx_ref = pref->dx_ref;
  int nx_ref = pref->nx_ref;
  int ny_ref = pref->ny_ref;
  int nz_ref = pref->nz_ref;

  /* Peskin support radius in refined grid units */
  double range_ref = PESKIN_RANGE;

  /* Loop over colloid cells */
  for (ic = 0; ic <= ncell[X] + 1; ic++) {
    for (jc = 0; jc <= ncell[Y] + 1; jc++) {
      for (kc = 0; kc <= ncell[Z] + 1; kc++) {

        colloids_info_cell_list_head(cinfo, ic, jc, kc, &pc);

        for (; pc; pc = pc->next) {

          if (pc->s.bc != COLLOID_BC_SUBGRID) continue;

          /* Particle position in local LB coordinates */
          double r0_lb[3];
          r0_lb[X] = pc->s.r[X] - 1.0 * offset_lb[X];
          r0_lb[Y] = pc->s.r[Y] - 1.0 * offset_lb[Y];
          r0_lb[Z] = pc->s.r[Z] - 1.0 * offset_lb[Z];

          /* Convert to refined grid coordinates */
          double r0_ref[3];
          psi_refined_lb_to_ref_position(pref, r0_lb, r0_ref);

          /* Particle net charge scaled by beta */
          double q_net = (pc->s.q0 - pc->s.q1) * pref->beta;

          /* Find range of affected refined grid nodes */
          int i_min = (int) floor(r0_ref[X] - range_ref);
          int i_max = (int) ceil(r0_ref[X] + range_ref);
          int j_min = (int) floor(r0_ref[Y] - range_ref);
          int j_max = (int) ceil(r0_ref[Y] + range_ref);
          int k_min = (int) floor(r0_ref[Z] - range_ref);
          int k_max = (int) ceil(r0_ref[Z] + range_ref);

          /* Clamp to valid range (with periodic handling) */
          /* For now, assume we stay within local domain */

          /* Distribute charge */
          for (int i = i_min; i <= i_max; i++) {
            for (int j = j_min; j <= j_max; j++) {
              for (int k = k_min; k <= k_max; k++) {

                /* Handle periodic boundary (wrap indices) */
                int ii = i;
                int jj = j;
                int kk = k;

                /* Periodic wrap */
                if (ii < 0) ii += nx_ref;
                if (ii >= nx_ref) ii -= nx_ref;
                if (jj < 0) jj += ny_ref;
                if (jj >= ny_ref) jj -= ny_ref;
                if (kk < 0) kk += nz_ref;
                if (kk >= nz_ref) kk -= nz_ref;

                /* Skip if out of bounds (non-periodic case) */
                if (ii < 0 || ii >= nx_ref) continue;
                if (jj < 0 || jj >= ny_ref) continue;
                if (kk < 0 || kk >= nz_ref) continue;

                /* Distance in refined grid units */
                double dr_x = r0_ref[X] - (double)i;
                double dr_y = r0_ref[Y] - (double)j;
                double dr_z = r0_ref[Z] - (double)k;

                /* Peskin weight (includes 1/dx scaling for each dimension) */
                double weight = psi_refined_d_peskin(dr_x, 1.0) *
                                psi_refined_d_peskin(dr_y, 1.0) *
                                psi_refined_d_peskin(dr_z, 1.0);

                /* Add charge to refined grid node */
                int idx_ref = ii * (ny_ref * nz_ref) + jj * nz_ref + kk;
                pref->rho_ref[idx_ref] += q_net * weight;
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
 *  psi_refined_fft_solve (CUDA implementation)
 *
 *****************************************************************************/

#ifdef __NVCC__

/* CUDA kernel: divide by Laplacian eigenvalue */
__global__ void psi_refined_divide_kernel(cufftDoubleComplex * psi_hat,
                                           const cufftDoubleComplex * rho_hat,
                                           const double * lambda,
                                           double epsilon,
                                           int nx, int ny, int nz_complex) {

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int n_total = nx * ny * nz_complex;

  if (idx >= n_total) return;

  double lam = lambda[idx];

  /* Handle k=0 mode */
  if (lam < 1.0e-15) {
    psi_hat[idx].x = 0.0;
    psi_hat[idx].y = 0.0;
  }
  else {
    double factor = 1.0 / (epsilon * lam);
    psi_hat[idx].x = rho_hat[idx].x * factor;
    psi_hat[idx].y = rho_hat[idx].y * factor;
  }
}

int psi_refined_fft_solve(psi_refined_t * pref) {

  int nx = pref->nx_ref;
  int ny = pref->ny_ref;
  int nz = pref->nz_ref;
  int n_real = nx * ny * nz;
  int n_complex = nx * ny * (nz / 2 + 1);
  cufftResult cufft_status;

  assert(pref);
  assert(pref->is_initialised);

  /* Copy rho from host to device */
  cudaMemcpy(pref->rho_real_d, pref->rho_ref, n_real * sizeof(double),
             cudaMemcpyHostToDevice);

  /* Forward FFT */
  cufft_status = cufftExecD2Z(
      *((cufftHandle *) pref->plan_forward),
      (cufftDoubleReal *) pref->rho_real_d,
      (cufftDoubleComplex *) pref->rho_complex_d);

  if (cufft_status != CUFFT_SUCCESS) {
    pe_fatal(pref->pe, "psi_refined: forward FFT failed\n");
    return -1;
  }

  /* Divide by Laplacian eigenvalue */
  int threads = 256;
  int blocks = (n_complex + threads - 1) / threads;

  psi_refined_divide_kernel<<<blocks, threads>>>(
      (cufftDoubleComplex *) pref->psi_complex_d,
      (cufftDoubleComplex *) pref->rho_complex_d,
      (double *) pref->laplacian_eigenval_d,
      pref->epsilon,
      nx, ny, nz / 2 + 1);

  cudaDeviceSynchronize();

  /* Inverse FFT */
  cufft_status = cufftExecZ2D(
      *((cufftHandle *) pref->plan_backward),
      (cufftDoubleComplex *) pref->psi_complex_d,
      (cufftDoubleReal *) pref->psi_real_d);

  if (cufft_status != CUFFT_SUCCESS) {
    pe_fatal(pref->pe, "psi_refined: inverse FFT failed\n");
    return -1;
  }

  /* Copy result to host and normalize */
  cudaMemcpy(pref->psi_ref, pref->psi_real_d, n_real * sizeof(double),
             cudaMemcpyDeviceToHost);

  double norm = 1.0 / (double) n_real;
  for (int i = 0; i < n_real; i++) {
    pref->psi_ref[i] *= norm;
  }

  return 0;
}

#else /* Non-CUDA stub */

int psi_refined_fft_solve(psi_refined_t * pref) {
  pe_fatal(pref->pe, "psi_refined_fft_solve requires CUDA\n");
  return -1;
}

#endif /* __NVCC__ */

/*****************************************************************************
 *
 *  psi_refined_transfer_potential_to_lb
 *
 *  Transfer potential from refined grid to LB grid at coincident nodes.
 *
 *****************************************************************************/

int psi_refined_transfer_potential_to_lb(psi_refined_t * pref) {

  psi_t * psi = NULL;
  int nlocal[3];
  int ic, jc, kc;
  int i_ref, j_ref, k_ref, idx_ref;

  assert(pref);
  assert(pref->psi_lb);

  psi = pref->psi_lb;
  cs_nlocal(pref->cs, nlocal);

  /* Loop over LB interior nodes */
  for (ic = 1; ic <= nlocal[X]; ic++) {
    for (jc = 1; jc <= nlocal[Y]; jc++) {
      for (kc = 1; kc <= nlocal[Z]; kc++) {

        int index_lb = cs_index(pref->cs, ic, jc, kc);

        /* Corresponding refined grid node */
        i_ref = (ic - 1) * pref->rfactor;
        j_ref = (jc - 1) * pref->rfactor;
        k_ref = (kc - 1) * pref->rfactor;

        idx_ref = i_ref * (pref->ny_ref * pref->nz_ref) +
                  j_ref * pref->nz_ref + k_ref;

        /* Set potential on LB grid */
        psi_psi_set(psi, index_lb, pref->psi_ref[idx_ref]);
      }
    }
  }

  /* Update halos */
  psi_halo_psi(psi);

  return 0;
}

/*****************************************************************************
 *
 *  psi_refined_compute_efield
 *
 *  Compute electric field on refined grid using finite differences.
 *  E = -grad(psi)
 *
 *****************************************************************************/

int psi_refined_compute_efield(psi_refined_t * pref) {

  int nx = pref->nx_ref;
  int ny = pref->ny_ref;
  int nz = pref->nz_ref;
  double dx = pref->dx_ref;

  assert(pref);

  /* Central difference for interior, periodic at boundaries */
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      for (int k = 0; k < nz; k++) {

        int idx = i * (ny * nz) + j * nz + k;

        /* Neighbor indices with periodic wrap */
        int ip = (i + 1) % nx;
        int im = (i - 1 + nx) % nx;
        int jp = (j + 1) % ny;
        int jm = (j - 1 + ny) % ny;
        int kp = (k + 1) % nz;
        int km = (k - 1 + nz) % nz;

        int idx_ip = ip * (ny * nz) + j * nz + k;
        int idx_im = im * (ny * nz) + j * nz + k;
        int idx_jp = i * (ny * nz) + jp * nz + k;
        int idx_jm = i * (ny * nz) + jm * nz + k;
        int idx_kp = i * (ny * nz) + j * nz + kp;
        int idx_km = i * (ny * nz) + j * nz + km;

        /* E = -grad(psi) using central difference */
        pref->efield_ref[3*idx + X] = -(pref->psi_ref[idx_ip] - pref->psi_ref[idx_im]) / (2.0 * dx);
        pref->efield_ref[3*idx + Y] = -(pref->psi_ref[idx_jp] - pref->psi_ref[idx_jm]) / (2.0 * dx);
        pref->efield_ref[3*idx + Z] = -(pref->psi_ref[idx_kp] - pref->psi_ref[idx_km]) / (2.0 * dx);
      }
    }
  }

  return 0;
}

/*****************************************************************************
 *
 *  psi_refined_interpolate_efield_to_particles
 *
 *  Interpolate electric field from refined grid to particles using Peskin.
 *
 *****************************************************************************/

int psi_refined_interpolate_efield_to_particles(psi_refined_t * pref,
                                                  colloids_info_t * cinfo) {

  int ic, jc, kc;
  int ncell[3];
  int nlocal_lb[3], offset_lb[3];
  colloid_t * pc = NULL;

  assert(pref);
  assert(cinfo);

  if (cinfo->nsubgrid == 0) return 0;

  cs_nlocal(pref->cs, nlocal_lb);
  cs_nlocal_offset(pref->cs, offset_lb);
  colloids_info_ncell(cinfo, ncell);

  int nx_ref = pref->nx_ref;
  int ny_ref = pref->ny_ref;
  int nz_ref = pref->nz_ref;

  double range_ref = PESKIN_RANGE;

  /* Loop over colloid cells */
  for (ic = 0; ic <= ncell[X] + 1; ic++) {
    for (jc = 0; jc <= ncell[Y] + 1; jc++) {
      for (kc = 0; kc <= ncell[Z] + 1; kc++) {

        colloids_info_cell_list_head(cinfo, ic, jc, kc, &pc);

        for (; pc; pc = pc->next) {

          if (pc->s.bc != COLLOID_BC_SUBGRID) continue;

          /* Particle position in local LB coordinates */
          double r0_lb[3];
          r0_lb[X] = pc->s.r[X] - 1.0 * offset_lb[X];
          r0_lb[Y] = pc->s.r[Y] - 1.0 * offset_lb[Y];
          r0_lb[Z] = pc->s.r[Z] - 1.0 * offset_lb[Z];

          /* Convert to refined grid coordinates */
          double r0_ref[3];
          psi_refined_lb_to_ref_position(pref, r0_lb, r0_ref);

          /* Initialize interpolated field */
          double E_interp[3] = {0.0, 0.0, 0.0};

          /* Find range of nodes for interpolation */
          int i_min = (int) floor(r0_ref[X] - range_ref);
          int i_max = (int) ceil(r0_ref[X] + range_ref);
          int j_min = (int) floor(r0_ref[Y] - range_ref);
          int j_max = (int) ceil(r0_ref[Y] + range_ref);
          int k_min = (int) floor(r0_ref[Z] - range_ref);
          int k_max = (int) ceil(r0_ref[Z] + range_ref);

          /* Interpolate using Peskin */
          for (int i = i_min; i <= i_max; i++) {
            for (int j = j_min; j <= j_max; j++) {
              for (int k = k_min; k <= k_max; k++) {

                /* Handle periodic boundary */
                int ii = i;
                int jj = j;
                int kk = k;

                if (ii < 0) ii += nx_ref;
                if (ii >= nx_ref) ii -= nx_ref;
                if (jj < 0) jj += ny_ref;
                if (jj >= ny_ref) jj -= ny_ref;
                if (kk < 0) kk += nz_ref;
                if (kk >= nz_ref) kk -= nz_ref;

                if (ii < 0 || ii >= nx_ref) continue;
                if (jj < 0 || jj >= ny_ref) continue;
                if (kk < 0 || kk >= nz_ref) continue;

                /* Distance in refined grid units */
                double dr_x = r0_ref[X] - (double)i;
                double dr_y = r0_ref[Y] - (double)j;
                double dr_z = r0_ref[Z] - (double)k;

                /* Peskin weight */
                double weight = psi_refined_d_peskin(dr_x, 1.0) *
                                psi_refined_d_peskin(dr_y, 1.0) *
                                psi_refined_d_peskin(dr_z, 1.0);

                /* Get field at this node */
                int idx_ref = ii * (ny_ref * nz_ref) + jj * nz_ref + kk;

                E_interp[X] += weight * pref->efield_ref[3*idx_ref + X];
                E_interp[Y] += weight * pref->efield_ref[3*idx_ref + Y];
                E_interp[Z] += weight * pref->efield_ref[3*idx_ref + Z];
              }
            }
          }

          /* Store interpolated field in particle (Esub in colloid) */
          pc->Esub[X] = E_interp[X];
          pc->Esub[Y] = E_interp[Y];
          pc->Esub[Z] = E_interp[Z];
        }
      }
    }
  }

  return 0;
}

/*****************************************************************************
 *
 *  psi_refined_solve
 *
 *  Main solve routine combining all steps.
 *
 *****************************************************************************/

int psi_refined_solve(psi_refined_t * pref, colloids_info_t * cinfo,
                       int ntimestep) {

  assert(pref);
  (void) ntimestep;

  /* Step 0: Reset arrays */
  psi_refined_reset(pref);

  /* Step 1: Transfer fluid charges from LB to refined grid */
  psi_refined_transfer_charges_from_lb(pref);

  /* Step 2: Distribute particle charges on refined grid */
  if (cinfo) {
    psi_refined_distribute_particle_charges(pref, cinfo);
  }

  /* Step 3: Solve Poisson equation on refined grid */
  psi_refined_fft_solve(pref);

  /* Step 4: Transfer potential back to LB grid */
  psi_refined_transfer_potential_to_lb(pref);

  /* Step 5: Compute electric field on refined grid */
  psi_refined_compute_efield(pref);

  /* Step 6: Interpolate field to particles */
  if (cinfo) {
    psi_refined_interpolate_efield_to_particles(pref, cinfo);
  }

  return 0;
}

/*****************************************************************************
 *
 *  psi_refined_info
 *
 *****************************************************************************/

int psi_refined_info(psi_refined_t * pref) {

  assert(pref);

  pe_info(pref->pe, "\n");
  pe_info(pref->pe, "Refined Grid Poisson Solver\n");
  pe_info(pref->pe, "---------------------------\n");
  pe_info(pref->pe, "Refinement factor:  %d\n", pref->rfactor);
  pe_info(pref->pe, "LB grid:            %d x %d x %d\n",
          pref->nx_lb, pref->ny_lb, pref->nz_lb);
  pe_info(pref->pe, "Refined grid:       %d x %d x %d\n",
          pref->nx_ref, pref->ny_ref, pref->nz_ref);
  pe_info(pref->pe, "Refined spacing:    %g\n", pref->dx_ref);
  pe_info(pref->pe, "Permittivity:       %14.7e\n", pref->epsilon);
  pe_info(pref->pe, "\n");

  return 0;
}
