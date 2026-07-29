/*****************************************************************************
 *
 *  psi_fft.c
 *
 *  Poisson solver using cuFFT for GPU acceleration.
 *
 *  Solves the Poisson equation:
 *    nabla^2 psi = -rho_elec / epsilon
 *
 *  In Fourier space:
 *    -k^2 * psi_hat(k) = -rho_hat(k) / epsilon
 *    psi_hat(k) = rho_hat(k) / (epsilon * k^2)
 *
 *  Algorithm:
 *    1. Copy charge density rho from psi_t to device array
 *    2. Forward FFT: rho(r) -> rho_hat(k)
 *    3. Divide by epsilon * k^2 (k=0 mode set to zero for periodic BC)
 *    4. Inverse FFT: psi_hat(k) -> psi(r)
 *    5. Copy result back to psi_t
 *
 *  CHANGE INIT - 20260119 New FFT-based Poisson solver using cuFFT
 *
 *****************************************************************************/

#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>

#include "pe.h"
#include "coords.h"
#include "field.h"
#include "psi_fft.h"
#include "psi_fft_pn.h"
#include "target.h"

 /* Forward declarations: kernels defined in subgrid.c */
 /*CHANGE INIT - P3M Peskin deconvolution */
double d_peskin(double r);
/*CHANGE END - P3M Peskin deconvolution */
/*CHANGE INIT - 20260422 B-spline deconvolution kernels */
double d_bspline4(double r);
double d_bspline6(double r);
double d_kb4(double r);
double d_peskin6(double r);
/*CHANGE END - 20260422 B-spline deconvolution kernels */

/* Only compile if using CUDA */
#ifdef __NVCC__

#include <cufft.h>
#include <cuda_runtime.h>

/* Virtual table for psi_solver interface */
static psi_solver_vt_t vt_ = {
  (psi_solver_free_ft)psi_solver_fft_free,
  (psi_solver_solve_ft)psi_solver_fft_solve
};

/* Forward declarations of CUDA kernels */
__global__ void poisson_divide_kernel(cufftDoubleComplex* psi_hat,
                                       const cufftDoubleComplex* rho_hat,
                                       double epsilon,
                                       double kx_factor, double ky_factor, double kz_factor,
                                       int nx, int ny, int nz_complex);

__global__ void poisson_divide_discrete_kernel(cufftDoubleComplex* psi_hat,
                                                const cufftDoubleComplex* rho_hat,
                                                const double* laplacian_eigenval,
                                                double epsilon,
                                                int nx, int ny, int nz_complex);

__global__ void poisson_divide_hockney_kernel(cufftDoubleComplex* psi_hat,
                                               const cufftDoubleComplex* rho_hat,
                                               const double* G_opt,
                                               int nx, int ny, int nz_complex);

__global__ void copy_rho_to_real_kernel(double* rho_real,
                                         const double* rho_field,
                                         double beta,
                                         int nlocal_x, int nlocal_y, int nlocal_z,
                                         int nhalo);

__global__ void copy_psi_from_real_kernel(double* psi_field,
                                           const double* psi_real,
                                           double e0x, double e0y, double e0z,
                                           int nlocal_x, int nlocal_y, int nlocal_z,
                                           int nhalo);

/* Forward declaration for internal create function */
static int psi_solver_fft_create_internal(psi_t* psi, psi_fft_laplacian_t laplacian_type,
                                           psi_solver_fft_t** psolver);

/*****************************************************************************
 *
 *  psi_solver_fft_create
 *
 *  Create FFT-based Poisson solver with default (analytic) Laplacian.
 *
 *****************************************************************************/

int psi_solver_fft_create(psi_t* psi, psi_solver_fft_t** psolver) {
  return psi_solver_fft_create_internal(psi, PSI_FFT_LAPLACIAN_ANALYTIC, psolver);
}

/*****************************************************************************
 *
 *  psi_solver_fft_create_opt
 *
 *  Create FFT-based Poisson solver with specified Laplacian type.
 *
 *****************************************************************************/

int psi_solver_fft_create_opt(psi_t* psi, psi_fft_laplacian_t laplacian_type,
                               psi_solver_fft_t** psolver) {
  return psi_solver_fft_create_internal(psi, laplacian_type, psolver);
}

/*****************************************************************************
 *
 *  psi_solver_fft_create_internal
 *
 *  Create and initialize the FFT-based Poisson solver.
 *
 *  For discrete Laplacian, computes the eigenvalues:
 *    lambda(k) = sum_p wlaplacian[p] * exp(i * k . cv[p])
 *
 *  Since wlaplacian and cv are real, and cv are symmetric, lambda is real:
 *    lambda(k) = sum_p wlaplacian[p] * cos(k . cv[p])
 *
 *****************************************************************************/

static int psi_solver_fft_create_internal(psi_t* psi, psi_fft_laplacian_t laplacian_type,
                                           psi_solver_fft_t** psolver) {

  psi_solver_fft_t* solver = NULL;
  int nlocal[3];
  int ntotal[3];
  double ltot[3];
  cufftResult cufft_status;
  cudaError_t cuda_status;

  assert(psi);
  assert(psolver);

  solver = (psi_solver_fft_t*)calloc(1, sizeof(psi_solver_fft_t));
  if (solver == NULL) {
    pe_fatal(psi->pe, "psi_solver_fft_create: calloc failed\n");
    return -1;
  }

  solver->super.impl = &vt_;
  solver->psi = psi;
  solver->is_initialised = 0;
  solver->laplacian_type = laplacian_type;
  solver->laplacian_eigenval_d = NULL;
  solver->G_opt_d = NULL;
  solver->influence_type = PSI_FFT_INFLUENCE_SIMPLE;
  solver->shape_a = 1.0;
  solver->n_alias_max = 2;

  /* Get grid dimensions */
  cs_nlocal(psi->cs, nlocal);
  cs_ntotal(psi->cs, ntotal);
  cs_ltot(psi->cs, ltot);

  solver->nx = nlocal[X];
  solver->ny = nlocal[Y];
  solver->nz = nlocal[Z];
  solver->ntotal[X] = ntotal[X];
  solver->ntotal[Y] = ntotal[Y];
  solver->ntotal[Z] = ntotal[Z];
  solver->lx = ltot[X];
  solver->ly = ltot[Y];
  solver->lz = ltot[Z];

  /* Get permittivity */
  psi_epsilon(psi, &solver->epsilon);
  psi_beta(psi, &solver->beta);

  /* For single-process case, local = total */
  /* TODO: Add MPI support with distributed FFT */

  int nx = solver->nx;
  int ny = solver->ny;
  int nz = solver->nz;
  int n_real = nx * ny * nz;
  int n_complex = nx * ny * (nz / 2 + 1);

  /* Allocate device memory for FFT */
  cuda_status = cudaMalloc(&solver->rho_real_d, n_real * sizeof(cufftDoubleReal));
  if (cuda_status != cudaSuccess) {
    pe_fatal(psi->pe, "psi_solver_fft: cudaMalloc rho_real failed\n");
    free(solver);
    return -1;
  }

  cuda_status = cudaMalloc(&solver->rho_complex_d, n_complex * sizeof(cufftDoubleComplex));
  if (cuda_status != cudaSuccess) {
    pe_fatal(psi->pe, "psi_solver_fft: cudaMalloc rho_complex failed\n");
    cudaFree(solver->rho_real_d);
    free(solver);
    return -1;
  }

  cuda_status = cudaMalloc(&solver->psi_complex_d, n_complex * sizeof(cufftDoubleComplex));
  if (cuda_status != cudaSuccess) {
    pe_fatal(psi->pe, "psi_solver_fft: cudaMalloc psi_complex failed\n");
    cudaFree(solver->rho_real_d);
    cudaFree(solver->rho_complex_d);
    free(solver);
    return -1;
  }

  cuda_status = cudaMalloc(&solver->psi_real_d, n_real * sizeof(cufftDoubleReal));
  if (cuda_status != cudaSuccess) {
    pe_fatal(psi->pe, "psi_solver_fft: cudaMalloc psi_real failed\n");
    cudaFree(solver->rho_real_d);
    cudaFree(solver->rho_complex_d);
    cudaFree(solver->psi_complex_d);
    free(solver);
    return -1;
  }

  /* Compute and allocate discrete Laplacian eigenvalues if needed */
  if (laplacian_type == PSI_FFT_LAPLACIAN_DISCRETE) {
    stencil_t* s = psi->stencil;

    if (s == NULL) {
      pe_fatal(psi->pe, "psi_solver_fft: discrete Laplacian requires stencil\n");
      cudaFree(solver->rho_real_d);
      cudaFree(solver->rho_complex_d);
      cudaFree(solver->psi_complex_d);
      cudaFree(solver->psi_real_d);
      free(solver);
      return -1;
    }

    /* Allocate host array for eigenvalues */
    double* eigenval_h = (double*)malloc(n_complex * sizeof(double));
    if (eigenval_h == NULL) {
      pe_fatal(psi->pe, "psi_solver_fft: malloc eigenval failed\n");
      cudaFree(solver->rho_real_d);
      cudaFree(solver->rho_complex_d);
      cudaFree(solver->psi_complex_d);
      cudaFree(solver->psi_real_d);
      free(solver);
      return -1;
    }

    /* Compute eigenvalues: lambda(k) = sum_p wlaplacian[p] * cos(k . cv[p])
     * The wave vector k = (2*pi*ix/nx, 2*pi*iy/ny, 2*pi*iz/nz) in grid units
     * cv[p] are the stencil offsets (integers)
     * We store -lambda because Poisson: laplacian(psi) = -rho/epsilon
     * so psi_hat = rho_hat / (-lambda * epsilon) = rho_hat / (neg_lambda * epsilon)
     */
    for (int ix = 0; ix < nx; ix++) {
      for (int iy = 0; iy < ny; iy++) {
        for (int iz = 0; iz < nz / 2 + 1; iz++) {
          int idx = ix * ny * (nz / 2 + 1) + iy * (nz / 2 + 1) + iz;

          /* Wave numbers (integer indices) */
          double kx_idx = (ix <= nx / 2) ? (double)ix : (double)(ix - nx);
          double ky_idx = (iy <= ny / 2) ? (double)iy : (double)(iy - ny);
          double kz_idx = (double)iz;

          /* Wave vector components: k = 2*pi*n/N */
          double kx = 2.0 * M_PI * kx_idx / (double)nx;
          double ky = 2.0 * M_PI * ky_idx / (double)ny;
          double kz = 2.0 * M_PI * kz_idx / (double)nz;

          /* Compute eigenvalue: lambda = sum_p wlaplacian[p] * cos(k . cv[p])
           * For the discrete Laplacian stencil:
           *   wlaplacian[0] > 0 (central point)
           *   wlaplacian[p] < 0 for p != 0 (neighbors)
           *   sum_p wlaplacian[p] = 0
           * This gives lambda >= 0 for all k, with lambda = 0 only at k = 0.
           *
           * Poisson equation: nabla^2 psi = -rho/epsilon
           * In discrete form: sum_p wlaplacian[p] * psi_{i+cv[p]} = -rho_i/epsilon
           * In Fourier: lambda(k) * psi_hat(k) = -rho_hat(k)/epsilon
           * So: psi_hat = -rho_hat / (epsilon * lambda) = rho_hat / (epsilon * (-lambda))
           *
           * But since lambda > 0, we have: psi_hat = -rho_hat / (epsilon * lambda)
           * We store lambda directly (positive) and divide by it.
           */
          double lambda = 0.0;
          for (int p = 0; p < s->npoints; p++) {
            double kdotc = kx * s->cv[p][X] + ky * s->cv[p][Y] + kz * s->cv[p][Z];
            lambda += s->wlaplacian[p] * cos(kdotc);
          }

          /* Store lambda directly (positive for k != 0) */
          eigenval_h[idx] = lambda;
        }
      }
    }

    /* Copy to device */
    cuda_status = cudaMalloc(&solver->laplacian_eigenval_d, n_complex * sizeof(double));
    if (cuda_status != cudaSuccess) {
      pe_fatal(psi->pe, "psi_solver_fft: cudaMalloc eigenval failed\n");
      free(eigenval_h);
      cudaFree(solver->rho_real_d);
      cudaFree(solver->rho_complex_d);
      cudaFree(solver->psi_complex_d);
      cudaFree(solver->psi_real_d);
      free(solver);
      return -1;
    }

    cudaMemcpy(solver->laplacian_eigenval_d, eigenval_h, n_complex * sizeof(double),
               cudaMemcpyHostToDevice);
    free(eigenval_h);

    pe_info(psi->pe, "psi_solver_fft: using discrete Laplacian (stencil %d points)\n",
            s->npoints);
  }

  /* Create cuFFT plans for 3D R2C and C2R transforms */
  cufftHandle* plan_fwd = (cufftHandle*)malloc(sizeof(cufftHandle));
  cufftHandle* plan_bwd = (cufftHandle*)malloc(sizeof(cufftHandle));

  cufft_status = cufftPlan3d(plan_fwd, nx, ny, nz, CUFFT_D2Z);
  if (cufft_status != CUFFT_SUCCESS) {
    pe_fatal(psi->pe, "psi_solver_fft: cufftPlan3d forward failed: %d\n", cufft_status);
    if (solver->laplacian_eigenval_d) cudaFree(solver->laplacian_eigenval_d);
    cudaFree(solver->rho_real_d);
    cudaFree(solver->rho_complex_d);
    cudaFree(solver->psi_complex_d);
    cudaFree(solver->psi_real_d);
    free(plan_fwd);
    free(plan_bwd);
    free(solver);
    return -1;
  }

  cufft_status = cufftPlan3d(plan_bwd, nx, ny, nz, CUFFT_Z2D);
  if (cufft_status != CUFFT_SUCCESS) {
    pe_fatal(psi->pe, "psi_solver_fft: cufftPlan3d backward failed: %d\n", cufft_status);
    cufftDestroy(*plan_fwd);
    if (solver->laplacian_eigenval_d) cudaFree(solver->laplacian_eigenval_d);
    cudaFree(solver->rho_real_d);
    cudaFree(solver->rho_complex_d);
    cudaFree(solver->psi_complex_d);
    cudaFree(solver->psi_real_d);
    free(plan_fwd);
    free(plan_bwd);
    free(solver);
    return -1;
  }

  solver->plan_forward = plan_fwd;
  solver->plan_backward = plan_bwd;
  solver->is_initialised = 1;

  *psolver = solver;

  return 0;
}

/*****************************************************************************
 *
 *  psi_solver_fft_free
 *
 *****************************************************************************/

int psi_solver_fft_free(psi_solver_fft_t** psolver) {

  psi_solver_fft_t* solver = NULL;

  assert(psolver);
  solver = *psolver;

  if (solver) {
    if (solver->is_initialised) {
      /* Destroy cuFFT plans */
      if (solver->plan_forward) {
        cufftDestroy(*((cufftHandle*)solver->plan_forward));
        free(solver->plan_forward);
      }
      if (solver->plan_backward) {
        cufftDestroy(*((cufftHandle*)solver->plan_backward));
        free(solver->plan_backward);
      }

      /* Free device memory */
      if (solver->rho_real_d) cudaFree(solver->rho_real_d);
      if (solver->rho_complex_d) cudaFree(solver->rho_complex_d);
      if (solver->psi_complex_d) cudaFree(solver->psi_complex_d);
      if (solver->psi_real_d) cudaFree(solver->psi_real_d);
      if (solver->laplacian_eigenval_d) cudaFree(solver->laplacian_eigenval_d);
      if (solver->G_opt_d) cudaFree(solver->G_opt_d);
    }

    free(solver);
    *psolver = NULL;
  }

  return 0;
}

/*****************************************************************************
 *
 *  poisson_divide_kernel
 *
 *  CUDA kernel to compute psi_hat(k) = rho_hat(k) / (epsilon * k^2)
 *
 *  The k=0 mode is set to zero (removes arbitrary constant in potential).
 *
 *****************************************************************************/

__global__ void poisson_divide_kernel(cufftDoubleComplex* psi_hat,
                                       const cufftDoubleComplex* rho_hat,
                                       double epsilon,
                                       double kx_factor, double ky_factor, double kz_factor,
                                       int nx, int ny, int nz_complex) {

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int n_total = nx * ny * nz_complex;

  if (idx >= n_total) return;

  /* Convert linear index to 3D indices */
  int iz = idx % nz_complex;
  int iy = (idx / nz_complex) % ny;
  int ix = idx / (nz_complex * ny);

  /* Compute wave numbers */
  /* For real FFT, kz goes from 0 to nz/2 */
  /* kx and ky go from 0 to n/2, then -n/2+1 to -1 (wrapped) */
  double kx = (ix <= nx / 2) ? (double)ix : (double)(ix - nx);
  double ky = (iy <= ny / 2) ? (double)iy : (double)(iy - ny);
  double kz = (double)iz;  /* Always positive for R2C */

  /* Scale by 2*pi/L */
  kx *= kx_factor;
  ky *= ky_factor;
  kz *= kz_factor;

  double k_sq = kx * kx + ky * ky + kz * kz;

  /* Handle k=0 mode: set potential to zero */
  if (k_sq < 1.0e-15) {
    psi_hat[idx].x = 0.0;
    psi_hat[idx].y = 0.0;
  }
  else {
    /* psi_hat = rho_hat / (epsilon * k^2) */
    double factor = 1.0 / (epsilon * k_sq);
    psi_hat[idx].x = rho_hat[idx].x * factor;
    psi_hat[idx].y = rho_hat[idx].y * factor;
  }
}

/*****************************************************************************
 *
 *  poisson_divide_discrete_kernel
 *
 *  CUDA kernel using precomputed discrete Laplacian eigenvalues.
 *
 *  The discrete Laplacian stencil has eigenvalues lambda(k) >= 0, with
 *  lambda(0) = 0 (constant mode).
 *
 *  Poisson equation: (discrete Laplacian) psi = -rho/epsilon
 *  In Fourier: lambda(k) * psi_hat(k) = -rho_hat(k)/epsilon
 *  Solution: psi_hat(k) = -rho_hat(k) / (epsilon * lambda(k))
 *
 *  The k=0 mode (lambda ~ 0) is set to zero (arbitrary constant).
 *
 *****************************************************************************/

__global__ void poisson_divide_discrete_kernel(cufftDoubleComplex* psi_hat,
                                                const cufftDoubleComplex* rho_hat,
                                                const double* laplacian_eigenval,
                                                double epsilon,
                                                int nx, int ny, int nz_complex) {

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int n_total = nx * ny * nz_complex;

  if (idx >= n_total) return;

  double lambda = laplacian_eigenval[idx];

  /* Handle k=0 mode: lambda ~ 0, set potential to zero */
  if (lambda < 1.0e-15) {
    psi_hat[idx].x = 0.0;
    psi_hat[idx].y = 0.0;
  }
  else {
    /* psi_hat = -rho_hat / (epsilon * lambda) */
    // double factor = -1.0 / (epsilon * lambda);
    double factor = 1.0 / (epsilon * lambda);
    psi_hat[idx].x = rho_hat[idx].x * factor;
    psi_hat[idx].y = rho_hat[idx].y * factor;
  }
}

/*****************************************************************************
 *
 *  poisson_divide_hockney_kernel
 *
 *  CUDA kernel using precomputed optimal influence function Ĝ_opt(k) from
 *  Hockney & Eastwood Eq. 8-22.
 *
 *  This minimizes the error Q between mesh-calculated force and the
 *  reference Coulomb force for P3M algorithms.
 *
 *  psi_hat(k) = rho_hat(k) * Ĝ_opt(k)
 *
 *****************************************************************************/

__global__ void poisson_divide_hockney_kernel(cufftDoubleComplex* psi_hat,
                                               const cufftDoubleComplex* rho_hat,
                                               const double* G_opt,
                                               int nx, int ny, int nz_complex) {

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int n_total = nx * ny * nz_complex;

  if (idx >= n_total) return;

  double G = G_opt[idx];

  /* psi_hat = rho_hat * G_opt */
  psi_hat[idx].x = rho_hat[idx].x * G;
  psi_hat[idx].y = rho_hat[idx].y * G;
}

/*****************************************************************************
 *
 *  copy_rho_to_real_kernel
 *
 *  Copy charge density from Ludwig's field structure to contiguous array
 *  for FFT. Ludwig uses halo cells which we skip.
 *
 *****************************************************************************/

__global__ void copy_rho_to_real_kernel(double* rho_real,
                                         const double* rho_field,
                                         double beta,
                                         int nx, int ny, int nz,
                                         int nhalo,
                                         int nk,
                                         int nsites) {  /* total sites including halos */

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int n_total = nx * ny * nz;

  if (idx >= n_total) return;

  /* Convert linear index to 3D (ix, iy, iz) for FFT array
   * FFT uses row-major: idx = ix*(ny*nz) + iy*nz + iz */
  int iz = idx % nz;
  int iy = (idx / nz) % ny;
  int ix = idx / (nz * ny);

  /* Ludwig index calculation (Z fastest, column-major like Fortran):
   * str[Z] = 1
   * str[Y] = (nz + 2*nhalo)
   * str[X] = str[Y] * (ny + 2*nhalo)
   * index = str[X]*(nhalo + ic - 1) + str[Y]*(nhalo + jc - 1) + str[Z]*(nhalo + kc - 1)
   *
   * ic, jc, kc go from 1 to nx, ny, nz for interior cells
   * So ic = ix + 1, jc = iy + 1, kc = iz + 1
   */
  int str_z = 1;
  int str_y = (nz + 2 * nhalo);
  int str_x = str_y * (ny + 2 * nhalo);
  int ic = ix + 1;  /* Ludwig coords start at 1 */
  int jc = iy + 1;
  int kc = iz + 1;
  int ludwig_idx = str_x * (nhalo + ic - 1) + str_y * (nhalo + jc - 1) + str_z * (nhalo + kc - 1);

  /* Sum charge densities: rho_elec = sum_k z_k * rho_k
   * With ADDR_SOA: addr_rank1(nsites, nk, index, n) = nsites*n + index
   * rho[species 0] at nsites*0 + ludwig_idx
   * rho[species 1] at nsites*1 + ludwig_idx
   * scale with rho_real = rho_elec * eunit * beta;
   */
  double rho0 = rho_field[nsites * 0 + ludwig_idx];
  double rho1 = rho_field[nsites * 1 + ludwig_idx];
  double rho_elec = rho0 - rho1;

  /* CHANGE INIT - eunit fix: was missing eunit factor (harmless when eunit=1) */
  /* rho_real[idx] = rho_elec * beta; */
  rho_real[idx] = rho_elec * beta;  /* eunit=1 assumed; fix if eunit != 1 */
  /* CHANGE END */
}

/*****************************************************************************
 *
 *  copy_psi_from_real_kernel
 *
 *  Copy potential from FFT result back to Ludwig's field structure.
 *  cuFFT inverse transform is unnormalized, so we divide by N.
 *
 *  Also adds external electric field contribution: psi_ext = -E0 . r
 *  where r is the position in lattice units (ic, jc, kc go from 1 to n).
 *
 *****************************************************************************/

__global__ void copy_psi_from_real_kernel(double* psi_field,
                                           const double* psi_real,
                                           double e0x, double e0y, double e0z,
                                           int nx, int ny, int nz,
                                           int nhalo) {

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int n_total = nx * ny * nz;

  if (idx >= n_total) return;

  /* Convert linear index to 3D (ix, iy, iz) for FFT array
   * FFT uses row-major: idx = ix*(ny*nz) + iy*nz + iz */
  int iz = idx % nz;
  int iy = (idx / nz) % ny;
  int ix = idx / (nz * ny);

  /* Ludwig index calculation (Z fastest, column-major like Fortran):
   * str[Z] = 1
   * str[Y] = (nz + 2*nhalo)
   * str[X] = str[Y] * (ny + 2*nhalo)
   * index = str[X]*(nhalo + ic - 1) + str[Y]*(nhalo + jc - 1) + str[Z]*(nhalo + kc - 1)
   */
  int str_z = 1;
  int str_y = (nz + 2 * nhalo);
  int str_x = str_y * (ny + 2 * nhalo);
  int ic = ix + 1;  /* Ludwig coords start at 1 */
  int jc = iy + 1;
  int kc = iz + 1;
  int ludwig_idx = str_x * (nhalo + ic - 1) + str_y * (nhalo + jc - 1) + str_z * (nhalo + kc - 1);

  /* Normalize by N (cuFFT is unnormalized) */
  double norm = 1.0 / (double)(nx * ny * nz);

  /* FFT solution (Poisson equation) */
  double psi_poisson = psi_real[idx] * norm;

  /* Add external field contribution: psi_ext = -E0 . r
   * Position r uses Ludwig coordinates (ic, jc, kc) which go from 1 to n.
   * This gives E = -grad(psi) = E0 as expected. */
  double psi_ext = -(e0x * (double)ic + e0y * (double)jc + e0z * (double)kc);

  psi_field[ludwig_idx] = psi_poisson + psi_ext;
}

/*****************************************************************************
 *
 *  psi_solver_fft_solve
 *
 *  Main solve routine:
 *    1. Copy rho to device
 *    2. Forward FFT
 *    3. Divide by epsilon * k^2
 *    4. Inverse FFT
 *    5. Copy psi back
 *
 *****************************************************************************/

int psi_solver_fft_solve(psi_solver_fft_t* solver, int ntimestep) {

  psi_t* psi = NULL;
  cufftResult cufft_status;
  int nlocal[3];
  int nhalo;

  assert(solver);
  assert(solver->is_initialised);

  psi = solver->psi;

  cs_nlocal(psi->cs, nlocal);
  cs_nhalo(psi->cs, &nhalo);

  int nx = solver->nx;
  int ny = solver->ny;
  int nz = solver->nz;
  int n_real = nx * ny * nz;
  int n_complex = nx * ny * (nz / 2 + 1);

  /* Ensure rho data is on device (sync from host if needed) */
  field_memcpy(psi->rho, tdpMemcpyHostToDevice);

  /* Get the device data pointer.
   * In Ludwig, when ndevice > 0, target is allocated on device and
   * target->data points to device memory. We need to get that pointer.
   * The host structure field->target is a device pointer, so we copy
   * the data pointer from the device-side structure. */
  double* rho_data_d = NULL;
  double* psi_data_d = NULL;

  /* Get device pointer: field->target is device ptr, need to read target->data from device */
  size_t data_offset = offsetof(field_t, data);
  cudaMemcpy(&rho_data_d, (char*)(psi->rho->target) + data_offset, sizeof(double*),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(&psi_data_d, (char*)(psi->psi->target) + data_offset, sizeof(double*),
             cudaMemcpyDeviceToHost);

  /* Step 1: Copy charge density to contiguous array */
  int nk;
  psi_nk(psi, &nk);

  /* nsites = total sites including halos, needed for SOA addressing */
  int nsites = psi->nsites;

  int threads = 256;
  int blocks = (n_real + threads - 1) / threads;

  copy_rho_to_real_kernel << <blocks, threads >> > (
      (double*)solver->rho_real_d,
      rho_data_d,
      solver->beta,
      nx, ny, nz, nhalo, nk, nsites);

  cudaDeviceSynchronize();

  /* Step 2: Forward FFT (R2C) */
  cufft_status = cufftExecD2Z(
      *((cufftHandle*)solver->plan_forward),
      (cufftDoubleReal*)solver->rho_real_d,
      (cufftDoubleComplex*)solver->rho_complex_d);

  if (cufft_status != CUFFT_SUCCESS) {
    pe_fatal(psi->pe, "psi_solver_fft: forward FFT failed: %d\n", cufft_status);
    return -1;
  }

  /* Step 3: Divide by Laplacian eigenvalue (analytic, discrete, or Hockney optimal) */
  blocks = (n_complex + threads - 1) / threads;

  if ((solver->influence_type == PSI_FFT_INFLUENCE_HOCKNEY ||
    solver->influence_type == PSI_FFT_INFLUENCE_EWALD) && solver->G_opt_d != NULL) {
    /* Use precomputed influence function Ĝ(k) (Hockney optimal or Ewald) */
    poisson_divide_hockney_kernel << <blocks, threads >> > (
        (cufftDoubleComplex*)solver->psi_complex_d,
        (cufftDoubleComplex*)solver->rho_complex_d,
        (double*)solver->G_opt_d,
        nx, ny, nz / 2 + 1);
  }
  else if (solver->laplacian_type == PSI_FFT_LAPLACIAN_DISCRETE) {
    /* Use precomputed discrete Laplacian eigenvalues */
    poisson_divide_discrete_kernel << <blocks, threads >> > (
        (cufftDoubleComplex*)solver->psi_complex_d,
        (cufftDoubleComplex*)solver->rho_complex_d,
        (double*)solver->laplacian_eigenval_d,
        solver->epsilon,
        nx, ny, nz / 2 + 1);
  }
  else {
    /* Use analytic k^2 */
    double kx_factor = 2.0 * M_PI / solver->lx;
    double ky_factor = 2.0 * M_PI / solver->ly;
    double kz_factor = 2.0 * M_PI / solver->lz;

    poisson_divide_kernel << <blocks, threads >> > (
        (cufftDoubleComplex*)solver->psi_complex_d,
        (cufftDoubleComplex*)solver->rho_complex_d,
        solver->epsilon,
        kx_factor, ky_factor, kz_factor,
        nx, ny, nz / 2 + 1);
  }

  cudaDeviceSynchronize();

  /* Step 4: Inverse FFT (C2R) */
  cufft_status = cufftExecZ2D(
      *((cufftHandle*)solver->plan_backward),
      (cufftDoubleComplex*)solver->psi_complex_d,
      (cufftDoubleReal*)solver->psi_real_d);

  if (cufft_status != CUFFT_SUCCESS) {
    pe_fatal(psi->pe, "psi_solver_fft: inverse FFT failed: %d\n", cufft_status);
    return -1;
  }

  /* Step 5: Copy result back to Ludwig field and add external field */
  /* External field contribution: psi_ext = -E0 . r */
  double e0x = psi->e0[X];
  double e0y = psi->e0[Y];
  double e0z = psi->e0[Z];

  blocks = (n_real + threads - 1) / threads;

  copy_psi_from_real_kernel << <blocks, threads >> > (
      psi_data_d,
      (double*)solver->psi_real_d,
      e0x, e0y, e0z,
      nx, ny, nz, nhalo);

  cudaDeviceSynchronize();

  /* Sync back to host and update halos for psi field */
  field_memcpy(psi->psi, tdpMemcpyDeviceToHost);
  psi_halo_psi(psi);
  field_memcpy(psi->psi, tdpMemcpyHostToDevice);

  return 0;
}

/*****************************************************************************
 *
 *  psi_solver_fft_solve_with_subgrid
 *
 *  Solve Poisson equation including subgrid particle charges.
 *  This adds the particle charges to the mesh before solving.
 *
 *****************************************************************************/

int psi_solver_fft_solve_with_subgrid(psi_solver_fft_t* solver,
                                       colloids_info_t* cinfo,
                                       int ntimestep) {

  /* TODO: Implement subgrid charge addition similar to psi_petsc version */
  /* For now, just call the regular solve */

  assert(solver);
  (void)cinfo;  /* Suppress unused warning for now */

  return psi_solver_fft_solve(solver, ntimestep);
}

/*****************************************************************************
 *
 *  psi_solver_fft_info
 *
 *****************************************************************************/

int psi_solver_fft_info(psi_solver_fft_t* solver) {

  assert(solver);

  pe_info(solver->psi->pe, "\n");
  pe_info(solver->psi->pe, "Poisson solver (FFT/cuFFT)\n");
  pe_info(solver->psi->pe, "---------------------------\n");
  pe_info(solver->psi->pe, "Grid:         %d x %d x %d\n",
          solver->nx, solver->ny, solver->nz);
  pe_info(solver->psi->pe, "Permittivity: %14.7e\n", solver->epsilon);
  pe_info(solver->psi->pe, "Laplacian:    %s\n",
          solver->laplacian_type == PSI_FFT_LAPLACIAN_DISCRETE ?
          "discrete (stencil)" : "analytic (k^2)");
  pe_info(solver->psi->pe, "\n");

  return 0;
}

/*****************************************************************************
 *
 *  psi_solver_fft_create_hockney
 *
 *  Create FFT-based Poisson solver with Hockney optimal influence function.
 *  This uses Eq. 8-22 from Hockney & Eastwood to minimize mesh force errors.
 *
 *  Parameters:
 *    psi           - psi_t structure
 *    laplacian_type - type of Laplacian discretization
 *    shape_a       - S2 shape parameter (typically = H = 1.0)
 *    n_alias_max   - maximum alias index for sums (typically 2)
 *    psolver       - output solver
 *
 *****************************************************************************/

int psi_solver_fft_create_hockney(psi_t* psi, psi_fft_laplacian_t laplacian_type,
                                   double shape_a, int n_alias_max,
                                   psi_solver_fft_t** psolver) {
  int ierr = 0;

  assert(psi);
  assert(psolver);

  /* First create the standard solver */
  ierr = psi_solver_fft_create_internal(psi, laplacian_type, psolver);
  if (ierr != 0) return ierr;

  /* Then enable Hockney optimal influence function */
  ierr = psi_solver_fft_set_influence_hockney(*psolver, shape_a, n_alias_max);

  return ierr;
}

/*****************************************************************************
 *
 *  psi_solver_fft_set_influence_hockney
 *
 *  Enable Hockney optimal influence function Ĝ_opt(k) (Eq. 8-22).
 *  This precomputes Ĝ_opt for all k-vectors and stores them on the GPU.
 *
 *  The optimal influence function minimizes the error between the mesh force
 *  and the reference force for a given assignment scheme (U), reference (R),
 *  and differential operator (D).
 *
 *  Parameters:
 *    solver       - FFT solver (already created)
 *    shape_a      - S2 shape parameter (typically = H = 1.0, mesh spacing)
 *    n_alias_max  - maximum alias index for aliasing sums (typically 2)
 *
 *****************************************************************************/

int psi_solver_fft_set_influence_hockney(psi_solver_fft_t* solver,
                                          double shape_a, int n_alias_max) {
  int nx, ny, nz, nz_complex;
  int n_complex;
  double* G_opt_h = NULL;  /* Host array */
  double H;                  /* Mesh spacing (assumed uniform = 1.0) */
  double kx_factor, ky_factor, kz_factor;

  assert(solver);
  assert(solver->is_initialised);

  nx = solver->nx;
  ny = solver->ny;
  nz = solver->nz;
  nz_complex = nz / 2 + 1;
  n_complex = nx * ny * nz_complex;

  H = 1.0;  /* Lattice spacing in lattice units */

  /* Wavenumber factors: k_i = (2*pi/L_i) * n_i for n_i = 0, 1, ..., N-1 */
  kx_factor = 2.0 * M_PI / solver->lx;
  ky_factor = 2.0 * M_PI / solver->ly;
  kz_factor = 2.0 * M_PI / solver->lz;

  /* Allocate host array */
  G_opt_h = (double*)malloc(n_complex * sizeof(double));
  if (G_opt_h == NULL) {
    pe_fatal(solver->psi->pe, "psi_solver_fft_set_influence_hockney: malloc failed\n");
    return -1;
  }

  /* Compute Ĝ_opt(k) for all k-vectors */
  for (int ix = 0; ix < nx; ix++) {
    /* Wrap ix for negative frequencies: ix >= nx/2 corresponds to negative */
    int kx_idx = (ix <= nx / 2) ? ix : ix - nx;
    double kx = kx_factor * kx_idx;

    for (int iy = 0; iy < ny; iy++) {
      int ky_idx = (iy <= ny / 2) ? iy : iy - ny;
      double ky = ky_factor * ky_idx;

      for (int iz = 0; iz < nz_complex; iz++) {
        /* For R2C transform, iz goes from 0 to nz/2 (only positive kz stored) */
        double kz = kz_factor * iz;

        int idx = ix * (ny * nz_complex) + iy * nz_complex + iz;

        /* k = 0 mode: set G_opt = 0 (overall charge neutrality) */
        if (ix == 0 && iy == 0 && iz == 0) {
          G_opt_h[idx] = 0.0;
          continue;
        }

        /* Compute discrete Laplacian eigenvalue D̂²(k) = λ(k) */
        double cx = cos(kx * H);
        double cy = cos(ky * H);
        double cz = cos(kz * H);
        double lambda_k = 2.0 * (3.0 - cx - cy - cz);  /* = 6 - 2(cx + cy + cz) */

        /* Compute optimal influence function using Hockney Eq. 8-22 */
        G_opt_h[idx] = psi_fft_pn_optimal_influence_function(
            kx, ky, kz,
            H, shape_a,
            solver->epsilon,
            lambda_k,
            n_alias_max);
      }
    }
  }

  /* Allocate GPU array if not already allocated */
  if (solver->G_opt_d == NULL) {
    cudaMalloc(&solver->G_opt_d, n_complex * sizeof(double));
  }

  /* Copy to GPU */
  cudaMemcpy(solver->G_opt_d, G_opt_h, n_complex * sizeof(double),
             cudaMemcpyHostToDevice);

  /* Free host array */
  free(G_opt_h);

  /* Update solver parameters */
  solver->influence_type = PSI_FFT_INFLUENCE_HOCKNEY;
  solver->shape_a = shape_a;
  solver->n_alias_max = n_alias_max;

  pe_info(solver->psi->pe, "\n");
  pe_info(solver->psi->pe, "Hockney optimal influence function enabled:\n");
  pe_info(solver->psi->pe, "  S2 shape parameter a:  %14.7e\n", shape_a);
  pe_info(solver->psi->pe, "  Max alias index:       %d\n", n_alias_max);
  pe_info(solver->psi->pe, "\n");

  return 0;
}

/*****************************************************************************
 *
 *  psi_solver_fft_create_ewald
 *
 *  Create FFT-based Poisson solver with Ewald sum influence function.
 *  This uses the Ewald splitting to improve convergence.
 * Parameters:
 *   psi           - psi_t structure
 *   laplacian_type - type of Laplacian discretization
 *   alpha         - Ewald splitting parameter
 *   rcut          - Real-space cutoff for short-range correction
 *   psolver       - output solver
 *
 *****************************************************************************/

 /*CHANGE INIT - 20260422 deconv parameter in create_ewald */
int psi_solver_fft_create_ewald(psi_t* psi, psi_fft_laplacian_t laplacian_type,
                                   double alpha, double rcut,
                                   psi_fft_deconv_t deconv,
                                   psi_solver_fft_t** psolver) {
  int ierr = 0;

  assert(psi);
  assert(psolver);

  ierr = psi_solver_fft_create_internal(psi, laplacian_type, psolver);
  if (ierr != 0) return ierr;

  ierr = psi_solver_fft_set_influence_ewald(*psolver, alpha, rcut, deconv);

  return ierr;
}
/*CHANGE END - 20260422 deconv parameter in create_ewald */
/*****************************************************************************
 *
 *  psi_solver_fft_set_influence_ewald
 *
 *  Enable Ewald sum influence function.
 *
 *  The Ewald sum splits the Coulomb potential into:
 *    φ(r) = φ_real(r) + φ_recip(r) - φ_self
 *
 *  In Fourier space, the reciprocal part has the influence function:
 *    Ĝ_ewald(k) = (4π/k²) · exp(-k²/(4α²)) / ε
 *
 *  The Gaussian factor exp(-k²/4α²) makes the sum converge rapidly.
 *  The real-space correction adds: erfc(α·r)/(4πε·r) for nearby pairs.
 *
 *  Parameters:
 *    solver - FFT solver (already created)
 *    alpha  - Ewald splitting parameter
 *    rcut   - Real-space cutoff for short-range correction
 *
 *  Typical choice for alpha:
 *    alpha ≈ π/rcut  for balanced real/reciprocal computation
 *    or alpha ≈ 5/L for accuracy ~10^{-7}
 *
 *****************************************************************************/

 /*CHANGE INIT - 20260422 deconv parameter in set_influence_ewald */
int psi_solver_fft_set_influence_ewald(psi_solver_fft_t* solver,
                                        double alpha, double rcut,
                                        psi_fft_deconv_t deconv) {
  int nx, ny, nz, nz_complex;
  int n_complex;
  double* G_ewald_h = NULL;  /* Host array */
  double kx_factor, ky_factor, kz_factor;
  double alpha_sq_4;  /* 4·α² for Gaussian factor */

  assert(solver);
  assert(solver->is_initialised);
  assert(alpha > 0.0);
  assert(rcut > 0.0);

  nx = solver->nx;
  ny = solver->ny;
  nz = solver->nz;
  nz_complex = nz / 2 + 1;
  n_complex = nx * ny * nz_complex;

  alpha_sq_4 = 4.0 * alpha * alpha;

  /* Wavenumber factors: k_i = (2*pi/L_i) * n_i */
  kx_factor = 2.0 * M_PI / solver->lx;
  ky_factor = 2.0 * M_PI / solver->ly;
  kz_factor = 2.0 * M_PI / solver->lz;

  /* Allocate host array */
  G_ewald_h = (double*)malloc(n_complex * sizeof(double));
  if (G_ewald_h == NULL) {
    pe_fatal(solver->psi->pe, "psi_solver_fft_set_influence_ewald: malloc failed\n");
    return -1;
  }

  /* Compute Ĝ_ewald(k) for all k-vectors */
  for (int ix = 0; ix < nx; ix++) {
    /* Wrap ix for negative frequencies */
    int kx_idx = (ix <= nx / 2) ? ix : ix - nx;
    double kx = kx_factor * kx_idx;

    for (int iy = 0; iy < ny; iy++) {
      int ky_idx = (iy <= ny / 2) ? iy : iy - ny;
      double ky = ky_factor * ky_idx;

      for (int iz = 0; iz < nz_complex; iz++) {
        double kz = kz_factor * iz;

        int idx = ix * (ny * nz_complex) + iy * nz_complex + iz;

        /* k = 0 mode: set G_ewald = 0 (charge neutrality) */
        if (ix == 0 && iy == 0 && iz == 0) {
          G_ewald_h[idx] = 0.0;
          continue;
        }

        double k_sq = kx * kx + ky * ky + kz * kz;

        /* Ĝ_ewald(k) = (4π/k²) · exp(-k²/(4α²)) / ε */
        // double G_ewald = (4.0 * M_PI / k_sq) * exp(-k_sq / alpha_sq_4) / solver->epsilon;
        // double G_ewald = (1 / k_sq) * (1.0 - exp(-k_sq / alpha_sq_4)) / solver->epsilon;
        // double G_ewald = (1 / k_sq) * exp(-k_sq / alpha_sq_4) / solver->epsilon;

        /*CHANGE INIT - 20260422 deconvolution kernel selection
         * Ĝ(k) = exp(-k²/4α²) / (ε · k² · |P̂(k)|²)
         * P̂(k) = Πᵢ Σ_m w(m) cos(kᵢ·m)  separable in 3D */
        double px = 0.0, py = 0.0, pz = 0.0;
        if (deconv == PSI_FFT_DECONV_BSPLINE4) {
          for (int m = -2; m <= 2; m++) {
            double w = d_bspline4((double)m);
            px += w * cos(kx * m);
            py += w * cos(ky * m);
            pz += w * cos(kz * m);
          }
        }
        else if (deconv == PSI_FFT_DECONV_BSPLINE6) {
          for (int m = -3; m <= 3; m++) {
            double w = d_bspline6((double)m);
            px += w * cos(kx * m);
            py += w * cos(ky * m);
            pz += w * cos(kz * m);
          }
        }
        else if (deconv == PSI_FFT_DECONV_PESKIN) {
          /* PSI_FFT_DECONV_PESKIN  */
          for (int m = -2; m <= 2; m++) {
            double w = d_peskin((double)m);
            px += w * cos(kx * m);
            py += w * cos(ky * m);
            pz += w * cos(kz * m);
          }
        }
        else if (deconv == PSI_FFT_DECONV_KB4) {
          for (int m = -2; m <= 2; m++) {
            double w = d_kb4((double)m);
            px += w * cos(kx * m);
            py += w * cos(ky * m);
            pz += w * cos(kz * m);
          }
        }
        else if (deconv == PSI_FFT_DECONV_PESKIN6) {
          for (int m = -3; m <= 3; m++) {
            double w = d_peskin6((double)m);
            px += w * cos(kx * m);
            py += w * cos(ky * m);
            pz += w * cos(kz * m);
          }
        }
        else if (deconv == PSI_FFT_DECONV_NONE) {
          /* No deconvolution */
          px = 1.0;
          py = 1.0;
          pz = 1.0;
        }

        double P_hat = px * py * pz;
        double P_hat_sq = P_hat * P_hat;
        if (P_hat_sq < 1.0e-10) {
          G_ewald_h[idx] = 0.0;
          continue;
        }
        double G_ewald = exp(-k_sq / alpha_sq_4) / (solver->epsilon * k_sq * P_hat_sq);
        /*CHANGE END - 20260422 deconvolution kernel selection */
        G_ewald_h[idx] = G_ewald;
      }
    }
  }

  /* Allocate GPU array if not already allocated */
  if (solver->G_opt_d == NULL) {
    cudaMalloc(&solver->G_opt_d, n_complex * sizeof(double));
  }

  /* Copy to GPU */
  cudaMemcpy(solver->G_opt_d, G_ewald_h, n_complex * sizeof(double),
             cudaMemcpyHostToDevice);

  /* Free host array */
  free(G_ewald_h);

  /* Update solver parameters */
  solver->influence_type = PSI_FFT_INFLUENCE_EWALD;
  solver->ewald_alpha = alpha;
  solver->ewald_rcut = rcut;

  pe_info(solver->psi->pe, "\n");
  pe_info(solver->psi->pe, "Ewald sum influence function enabled:\n");
  pe_info(solver->psi->pe, "  Ewald alpha:           %14.7e\n", alpha);
  pe_info(solver->psi->pe, "  Real-space cutoff:     %14.7e\n", rcut);
  const char* deconv_name;
  switch (deconv) {
    case PSI_FFT_DECONV_NONE:     deconv_name = "None";             break;
    case PSI_FFT_DECONV_PESKIN:   deconv_name = "Peskin";           break;
    case PSI_FFT_DECONV_BSPLINE4: deconv_name = "B-spline order 4"; break;
    case PSI_FFT_DECONV_BSPLINE6: deconv_name = "B-spline order 6"; break;
    case PSI_FFT_DECONV_KB4:      deconv_name = "Kaiser-Bessel W=4"; break;
    case PSI_FFT_DECONV_PESKIN6:  deconv_name = "Peskin 6-point C3 (Bao 2016)"; break;
    default:                      deconv_name = "Unknown";          break;
  }
  pe_info(solver->psi->pe, "  Deconvolution kernel:  %s\n", deconv_name);
  pe_info(solver->psi->pe, "\n");

  return 0;
}
/*CHANGE END - 20260422 deconv parameter in set_influence_ewald */

/*CHANGE INIT - 20260710 spectral (FFT phase) shift of a scalar lattice field
 *
 *  psi_fft_shift_field
 *
 *  Translates a periodic scalar lattice field by a fractional displacement
 *  (ux, uy, uz) using the Fourier shift theorem: each mode k is multiplied
 *  by exp(-i k.u), which is an EXACT, unitary translation of the lattice
 *  field (no amplitude loss for any mode, unlike linear interpolation which
 *  kills the Nyquist mode completely at u = 0.5). For integer u this reduces
 *  to an exact index relabelling.
 *
 *  Nyquist handling: for even N the m = N/2 mode is shared between +N/2 and
 *  -N/2; keeping the output real requires replacing that dimension's phase
 *  factor by its real part cos(pi*u) (standard treatment for fractional
 *  spectral shifts of real fields; exact +/-1 for integer u).
 *
 *  field_data: HOST pointer to one scalar field in Ludwig halo layout
 *  (e.g. psi->psi->data, or &psi->rho->data[nsites*n] for species n).
 *  Only the interior is transformed; caller must re-halo afterwards.
 *
 *  Serial / single-domain only (checked): the FFT needs the full box.
 *
 *****************************************************************************/

/* CHANGE 20260714: the Nyquist factor per dimension is now an explicit
 * argument (nyqx/nyqy/nyqz). Policy set by the host wrappers:
 *   - integer u:     nyq = cos(pi*u) = +/-1  (exact relabelling)
 *   - fractional u:  nyq = 0                 (kill the plane: within the
 *     band-limited subspace |m| < N/2 the phase shift is exactly unitary
 *     for ANY u; the previous cos(pi*u) attenuation was u-dependent and
 *     non-unitary on round trips)
 *   - pure filtering: u = 0, nyq = 0 kills the plane, nyq = 1 no-op. */
__global__ void spectral_shift_phase_kernel(cufftDoubleComplex* fh,
                                            int nx, int ny, int nz,
                                            double ux, double uy, double uz,
                                            double nyqx, double nyqy,
                                            double nyqz) {
  int nzh = nz / 2 + 1;
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int n_total = nx * ny * nzh;
  if (idx >= n_total) return;

  int iz = idx % nzh;
  int iy = (idx / nzh) % ny;
  int ix = idx / (nzh * ny);

  /* Signed mode numbers */
  int mx = (ix <= nx / 2) ? ix : ix - nx;
  int my = (iy <= ny / 2) ? iy : iy - ny;
  int mz = iz;   /* R2C half-spectrum: 0 .. nz/2 */

  const double twopi = 2.0 * M_PI;

  double fr = 1.0, fi = 0.0;

  double th, cr, ci, tr;
  /* X */
  if (2 * mx == nx || 2 * mx == -nx) { cr = nyqx; ci = 0.0; }
  else { th = twopi * mx * ux / nx; cr = cos(th); ci = -sin(th); }
  tr = fr * cr - fi * ci; fi = fr * ci + fi * cr; fr = tr;
  /* Y */
  if (2 * my == ny || 2 * my == -ny) { cr = nyqy; ci = 0.0; }
  else { th = twopi * my * uy / ny; cr = cos(th); ci = -sin(th); }
  tr = fr * cr - fi * ci; fi = fr * ci + fi * cr; fr = tr;
  /* Z */
  if (2 * mz == nz) { cr = nyqz; ci = 0.0; }
  else { th = twopi * mz * uz / nz; cr = cos(th); ci = -sin(th); }
  tr = fr * cr - fi * ci; fi = fr * ci + fi * cr; fr = tr;

  double re = fh[idx].x, im = fh[idx].y;
  fh[idx].x = re * fr - im * fi;
  fh[idx].y = re * fi + im * fr;
}

/* Internal worker: spectral phase shift with explicit Nyquist factors. */
static int psi_fft_spectral_apply(psi_t* psi, double* field_data,
                                  double ux, double uy, double uz,
                                  double nyqx, double nyqy, double nyqz) {

  assert(psi);
  assert(field_data);

  int nlocal[3], ntotal[3], nhalo;
  cs_nlocal(psi->cs, nlocal);
  cs_ntotal(psi->cs, ntotal);
  cs_nhalo(psi->cs, &nhalo);

  if (nlocal[X] != ntotal[X] || nlocal[Y] != ntotal[Y] ||
      nlocal[Z] != ntotal[Z]) {
    pe_fatal(psi->pe, "psi_fft_shift_field: serial/single-domain only\n");
  }

  int nx = nlocal[X], ny = nlocal[Y], nz = nlocal[Z];
  int n = nx * ny * nz;
  int nzh = nz / 2 + 1;

  /* Cached plans and device buffers (single-threaded host code) */
  static int sh_nx = 0, sh_ny = 0, sh_nz = 0;
  static cufftHandle sh_fwd, sh_bwd;
  static double* d_real = NULL;
  static cufftDoubleComplex* d_cplx = NULL;
  static double* h_real = NULL;

  if (sh_nx != nx || sh_ny != ny || sh_nz != nz) {
    if (d_real) { cudaFree(d_real);  d_real = NULL; }
    if (d_cplx) { cudaFree(d_cplx);  d_cplx = NULL; }
    if (h_real) { free(h_real);      h_real = NULL; }
    if (sh_nx > 0) { cufftDestroy(sh_fwd); cufftDestroy(sh_bwd); }

    if (cufftPlan3d(&sh_fwd, nx, ny, nz, CUFFT_D2Z) != CUFFT_SUCCESS ||
        cufftPlan3d(&sh_bwd, nx, ny, nz, CUFFT_Z2D) != CUFFT_SUCCESS) {
      pe_fatal(psi->pe, "psi_fft_shift_field: cufftPlan3d failed\n");
    }
    cudaMalloc((void**)&d_real, (size_t)n * sizeof(double));
    cudaMalloc((void**)&d_cplx, (size_t)nx * ny * nzh * sizeof(cufftDoubleComplex));
    h_real = (double*)malloc((size_t)n * sizeof(double));
    sh_nx = nx; sh_ny = ny; sh_nz = nz;
  }

  /* Pack interior (Ludwig halo layout -> row-major x,y,z) on host */
  int str_z = 1;
  int str_y = (nz + 2 * nhalo);
  int str_x = str_y * (ny + 2 * nhalo);
  for (int ix = 0; ix < nx; ix++) {
    for (int iy = 0; iy < ny; iy++) {
      for (int iz = 0; iz < nz; iz++) {
        int lidx = str_x * (nhalo + ix) + str_y * (nhalo + iy) + str_z * (nhalo + iz);
        h_real[ix * (ny * nz) + iy * nz + iz] = field_data[lidx];
      }
    }
  }

  cudaMemcpy(d_real, h_real, (size_t)n * sizeof(double), cudaMemcpyHostToDevice);

  if (cufftExecD2Z(sh_fwd, d_real, d_cplx) != CUFFT_SUCCESS)
    pe_fatal(psi->pe, "psi_fft_shift_field: forward FFT failed\n");

  int nthreads = 256;
  int nblocks = (nx * ny * nzh + nthreads - 1) / nthreads;
  spectral_shift_phase_kernel<<<nblocks, nthreads>>>(d_cplx, nx, ny, nz,
                                                     ux, uy, uz,
                                                     nyqx, nyqy, nyqz);
  cudaDeviceSynchronize();

  if (cufftExecZ2D(sh_bwd, d_cplx, d_real) != CUFFT_SUCCESS)
    pe_fatal(psi->pe, "psi_fft_shift_field: inverse FFT failed\n");

  cudaMemcpy(h_real, d_real, (size_t)n * sizeof(double), cudaMemcpyDeviceToHost);

  /* Unpack with 1/N normalisation (cuFFT is unnormalised) */
  double norm = 1.0 / (double)n;
  for (int ix = 0; ix < nx; ix++) {
    for (int iy = 0; iy < ny; iy++) {
      for (int iz = 0; iz < nz; iz++) {
        int lidx = str_x * (nhalo + ix) + str_y * (nhalo + iy) + str_z * (nhalo + iz);
        field_data[lidx] = h_real[ix * (ny * nz) + iy * nz + iz] * norm;
      }
    }
  }

  return 0;
}

/* Nyquist factor for one axis: exact +/-1 phase for integer u; 0 (kill the
 * plane) for fractional u — see policy note above the phase kernel. */
static double psi_fft_nyq_factor(double u) {
  double frac = u - floor(u);
  return (frac == 0.0) ? cos(M_PI * u) : 0.0;
}

int psi_fft_shift_field(psi_t* psi, double* field_data,
                        double ux, double uy, double uz) {
  return psi_fft_spectral_apply(psi, field_data, ux, uy, uz,
                                psi_fft_nyq_factor(ux),
                                psi_fft_nyq_factor(uy),
                                psi_fft_nyq_factor(uz));
}

/* CHANGE 20260715: shift that LEAVES the Nyquist plane in place (factor 1).
 * For transporting a real DENSITY (rho) round-trip: the |m|<N/2 modes
 * translate exactly (unitary phase) and the Nyquist plane is preserved
 * intact rather than killed. The half-cell translation of the (-1)^i mode
 * has no real representation, so it stays put — but with factor 1 the round
 * trip T_{+u} . T_{-u} is the exact identity on the WHOLE field, so NP's
 * regenerated Nyquist content (~1e-9) is no longer destroyed u-dependently
 * each step (which was the source of the random-u force noise). */
int psi_fft_shift_field_keepnyq(psi_t* psi, double* field_data,
                                double ux, double uy, double uz) {
  double fx = (ux - floor(ux) == 0.0) ? cos(M_PI*ux) : 1.0;
  double fy = (uy - floor(uy) == 0.0) ? cos(M_PI*uy) : 1.0;
  double fz = (uz - floor(uz) == 0.0) ? cos(M_PI*uz) : 1.0;
  return psi_fft_spectral_apply(psi, field_data, ux, uy, uz, fx, fy, fz);
}

/*****************************************************************************
 *
 *  psi_fft_kill_nyquist
 *
 *  Band-limit a lattice field by zeroing the Nyquist plane(s) of the
 *  selected axes (killx/killy/killz nonzero). Projects the field onto the
 *  subspace |m| < N/2 in those axes, within which fractional spectral
 *  shifts are exactly unitary. Applied to the particle's scattered charge
 *  so that the Poisson solution has no content the -u un-shift can lose.
 *
 *****************************************************************************/
int psi_fft_kill_nyquist(psi_t* psi, double* field_data,
                         int killx, int killy, int killz) {
  return psi_fft_spectral_apply(psi, field_data, 0.0, 0.0, 0.0,
                                killx ? 0.0 : 1.0,
                                killy ? 0.0 : 1.0,
                                killz ? 0.0 : 1.0);
}

/*CHANGE INIT - 20260715 diagnostic: L2 norm of the Nyquist-x plane of a field.
 * Returns sqrt(sum over k with mx=nx/2 of |field_hat(k)|^2) / N, i.e. the
 * real-space RMS amplitude carried by the Nyquist-x plane. Used to trace
 * where u-dependent Nyquist content enters/leaves rho and psi. Standalone
 * FFT (own cached plan) so it does not disturb the shift buffers. */
double psi_fft_nyquist_x_norm(psi_t* psi, const double* field_data) {

  int nlocal[3], ntotal[3], nhalo;
  cs_nlocal(psi->cs, nlocal);
  cs_ntotal(psi->cs, ntotal);
  cs_nhalo(psi->cs, &nhalo);
  if (nlocal[X] != ntotal[X]) return -1.0;

  int nx = nlocal[X], ny = nlocal[Y], nz = nlocal[Z];
  int n = nx * ny * nz, nzh = nz / 2 + 1;

  static int dn_nx = 0, dn_ny = 0, dn_nz = 0;
  static cufftHandle dn_fwd;
  static double* dn_dreal = NULL;
  static cufftDoubleComplex* dn_dcplx = NULL;
  static double* dn_hreal = NULL;
  static cufftDoubleComplex* dn_hcplx = NULL;

  if (dn_nx != nx || dn_ny != ny || dn_nz != nz) {
    if (dn_dreal) cudaFree(dn_dreal);
    if (dn_dcplx) cudaFree(dn_dcplx);
    if (dn_hreal) free(dn_hreal);
    if (dn_hcplx) free(dn_hcplx);
    if (dn_nx > 0) cufftDestroy(dn_fwd);
    cufftPlan3d(&dn_fwd, nx, ny, nz, CUFFT_D2Z);
    cudaMalloc((void**)&dn_dreal, (size_t)n * sizeof(double));
    cudaMalloc((void**)&dn_dcplx, (size_t)nx * ny * nzh * sizeof(cufftDoubleComplex));
    dn_hreal = (double*)malloc((size_t)n * sizeof(double));
    dn_hcplx = (cufftDoubleComplex*)malloc((size_t)nx * ny * nzh * sizeof(cufftDoubleComplex));
    dn_nx = nx; dn_ny = ny; dn_nz = nz;
  }

  int str_z = 1, str_y = (nz + 2*nhalo), str_x = str_y * (ny + 2*nhalo);
  for (int ix = 0; ix < nx; ix++)
    for (int iy = 0; iy < ny; iy++)
      for (int iz = 0; iz < nz; iz++) {
        int lidx = str_x*(nhalo+ix) + str_y*(nhalo+iy) + str_z*(nhalo+iz);
        dn_hreal[ix*(ny*nz) + iy*nz + iz] = field_data[lidx];
      }
  cudaMemcpy(dn_dreal, dn_hreal, (size_t)n*sizeof(double), cudaMemcpyHostToDevice);
  cufftExecD2Z(dn_fwd, dn_dreal, dn_dcplx);
  cudaMemcpy(dn_hcplx, dn_dcplx,
             (size_t)nx*ny*nzh*sizeof(cufftDoubleComplex), cudaMemcpyDeviceToHost);

  /* Sum |.|^2 over the mx = nx/2 plane (Hermitian: the half-spectrum holds it) */
  double s = 0.0;
  int mxn = nx / 2;
  for (int iy = 0; iy < ny; iy++)
    for (int iz = 0; iz < nzh; iz++) {
      cufftDoubleComplex c = dn_hcplx[mxn*(ny*nzh) + iy*nzh + iz];
      s += c.x*c.x + c.y*c.y;
    }
  return sqrt(s) / (double)n;
}
/*CHANGE END - 20260715 */
/*CHANGE END - 20260710 spectral shift */

#else /* __NVCC__ not defined */

/*****************************************************************************
 *
 *  Stub implementations when CUDA is not available
 *
 *****************************************************************************/

int psi_solver_fft_create(psi_t* psi, psi_solver_fft_t** psolver) {
  (void)psi;
  (void)psolver;
  pe_fatal(psi->pe, "psi_solver_fft requires CUDA (nvcc)\n");
  return -1;
}

int psi_solver_fft_create_opt(psi_t* psi, psi_fft_laplacian_t laplacian_type,
                               psi_solver_fft_t** psolver) {
  (void)psi;
  (void)laplacian_type;
  (void)psolver;
  pe_fatal(psi->pe, "psi_solver_fft requires CUDA (nvcc)\n");
  return -1;
}

int psi_solver_fft_create_hockney(psi_t* psi, psi_fft_laplacian_t laplacian_type,
                                   double shape_a, int n_alias_max,
                                   psi_solver_fft_t** psolver) {
  (void)psi;
  (void)laplacian_type;
  (void)shape_a;
  (void)n_alias_max;
  (void)psolver;
  pe_fatal(psi->pe, "psi_solver_fft requires CUDA (nvcc)\n");
  return -1;
}

int psi_solver_fft_set_influence_hockney(psi_solver_fft_t* solver,
                                          double shape_a, int n_alias_max) {
  (void)solver;
  (void)shape_a;
  (void)n_alias_max;
  return -1;
}

int psi_solver_fft_set_influence_ewald(psi_solver_fft_t* solver,
                                        double alpha, double rcut,
                                        psi_fft_deconv_t deconv) {
  (void)solver;
  (void)alpha;
  (void)rcut;
  (void)deconv;
  return -1;
}

int psi_solver_fft_free(psi_solver_fft_t** psolver) {
  (void)psolver;
  return 0;
}

int psi_solver_fft_solve(psi_solver_fft_t* solver, int ntimestep) {
  (void)solver;
  (void)ntimestep;
  return -1;
}

/*CHANGE INIT - 20260710 spectral shift stub (no CUDA) */
int psi_fft_shift_field(psi_t* psi, double* field_data,
                        double ux, double uy, double uz) {
  (void)field_data; (void)ux; (void)uy; (void)uz;
  pe_fatal(psi->pe, "psi_fft_shift_field requires the CUDA/cuFFT build\n");
  return -1;
}

int psi_fft_kill_nyquist(psi_t* psi, double* field_data,
                         int killx, int killy, int killz) {
  (void)field_data; (void)killx; (void)killy; (void)killz;
  pe_fatal(psi->pe, "psi_fft_kill_nyquist requires the CUDA/cuFFT build\n");
  return -1;
}

double psi_fft_nyquist_x_norm(psi_t* psi, const double* field_data) {
  (void)field_data;
  pe_fatal(psi->pe, "psi_fft_nyquist_x_norm requires the CUDA/cuFFT build\n");
  return -1.0;
}

int psi_fft_shift_field_keepnyq(psi_t* psi, double* field_data,
                                double ux, double uy, double uz) {
  (void)field_data; (void)ux; (void)uy; (void)uz;
  pe_fatal(psi->pe, "psi_fft_shift_field_keepnyq requires the CUDA/cuFFT build\n");
  return -1;
}
/*CHANGE END - 20260710 */

int psi_solver_fft_solve_with_subgrid(psi_solver_fft_t* solver,
                                       colloids_info_t* cinfo,
                                       int ntimestep) {
  (void)solver;
  (void)cinfo;
  (void)ntimestep;
  return -1;
}

int psi_solver_fft_info(psi_solver_fft_t* solver) {
  (void)solver;
  return -1;
}

#endif /* __NVCC__ */
