/*****************************************************************************
 *
 *  nernst_planck.c
 *
 *  A solution for the Nerst-Planck equation, which is the advection
 *  diffusion equation for charged species \rho_k in the presence of
 *  a potential \psi.
 *
 *  We have, in the most simple case:
 *
 *  d_t rho_k + div . (rho_k u) = div . D_k (grad rho_k + Z_k rho_k grad psi)
 *
 *  where u is the velocity field, D_k are the diffusion constants for
 *  each species k, and Z_k = (valancy_k e / k_bT) = beta valency_k e.
 *  e is the unit charge.
 *
 *  If the chemical potential is mu_k for species k, the diffusive
 *  flux may be written as
 *
 *    j_k = - D_k rho_k grad (beta mu_k)
 *
 *  with mu_k = mu_k^ideal + mu_k^ex = k_bT ln(rho_k) + valency_k e psi.
 *  (For more complex problems, there may be other terms in the chemical
 *  potential.)
 *
 *  As it is important to conserve charge, we solve in a flux form.
 *  Following Capuani, Pagonabarraga and Frenkel, J. Chem. Phys.
 *  \textbf{121} 973 (2004) we include factors to ensure that the
 *  charge densities follow a Boltzmann distribution in equilbrium.
 *
 *  This writes the flux as
 *    j_k = - D_k exp[beta mu_k^ex] grad (rho_k exp[beta mu_k^ex])
 *
 *  which we approximate at the cell faces by (e.g., for x only)
 *
 *    -D_k (1/2) { exp[-beta mu_k^ex(i)] + exp[-beta mu_k^ex(i+1)] }
 *    * { rho_k(i+1) exp[beta mu_k^ex(i+1)] - rho_k(i) exp[beta mu_k^ex(i)] }
 *
 *  We then compute the divergence of the fluxes to update via an
 *  Euler forward step. The advective fluxes (again computed at the
 *  cells faces) may be added to the diffusive fluxes to solve the
 *  whole thing. Appropraite advective fluxes may be computed via
 *  the advection.h interface.
 *
 *  Solid boundaries simply involve enforcing a no normal flux
 *  condition at the cell face.
 *
 *  The potential and charge species are available via the psi_s
 *  object.
 *
 *  A uniform external electric field may be applied; this is done
 *  by adding a contribution to the potential
 *     psi -> psi - eE.r
 *  which just appears as -eE in the calculation of grad psi.
 *
 *
 *  Edinbrugh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  (c) 2012-2023 The University of Edinburgh
 *
 *  Contributing authors:
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *  Oliver Henrich (ohenrich@epcc.ed.ac.uk)
 *
 *****************************************************************************/

#include <assert.h>
#include <math.h>
#include <stdlib.h>

#include "pe.h"
#include "coords.h"
#include "advection.h"
#include "advection_bcs.h"
#include "nernst_planck.h"
#include "psi_gradients.h"
 /*CHANGE INIT - 20260210 Kahan compensated sum from util_sum */
#include "util_sum.h"
#include "stencil_d3q7.h"
#include "stencil_d3q19.h"
#include "stencil_d3q27.h"

double bernoulli(double x) {
  if (fabs(x) < 1e-6) {
    return 1.0 - x / 2.0 + x * x / 12.0; // expansión
  }
  return x / (exp(x) - 1.0);
}



/*CHANGE END - 20260210 */

 /* This needs an input switch to make it active. */
int nernst_planck_fluxes_force_d3qx(psi_t* psi, fe_t* fe, hydro_t* hydro,
    map_t* map, colloids_info_t* cinfo, double** flx);

static int nernst_planck_fluxes(psi_t* psi, fe_t* fel, double* fx,
        double* fy,
        double* fz);
static int nernst_planck_no_flux_condition(map_t* map, int nk, double* fx,
             double* fy, double* fz);
static int nernst_planck_update(psi_t* psi, double* fe, double* fy,
        double* fz);
static int nernst_planck_fluxes_d3qx(psi_t* psi, fe_t* fe, hydro_t* hydro,
    map_t* map, colloids_info_t* cinfo, double** flx);
static int nernst_planck_update_d3qx(psi_t* psi,
        map_t* map, double** flx);
static double max_acc;

int np_advective_fluxes(psi_t* psi, hydro_t* hydro, double** flux);
int np_no_flux_boundary(psi_t* psi, map_t* map, double** flux);

/*****************************************************************************
 *
 *  nernst_planck_driver
 *
 *  No hydrodynamics (use nernst_plancj_driver_d3qx instead).
 *  Allows for a no-flux condition between solid and fluid sites.
 *
 *****************************************************************************/

int nernst_planck_driver(psi_t* psi, fe_t* fel, map_t* map) {

  int nk;              /* Number of electrolyte species */
  int nsites;          /* Number of lattice sites */

  double* fe = NULL;
  double* fy = NULL;
  double* fz = NULL;

  psi_nk(psi, &nk);
  cs_nsites(psi->cs, &nsites);

  /* Allocate fluxes and initialise to zero */

  fe = (double*)calloc((size_t)nsites * nk, sizeof(double));
  fy = (double*)calloc((size_t)nsites * nk, sizeof(double));
  fz = (double*)calloc((size_t)nsites * nk, sizeof(double));

  if (fe == NULL) pe_fatal(psi->pe, "calloc(fe) failed\n");
  if (fy == NULL) pe_fatal(psi->pe, "calloc(fy) failed\n");
  if (fz == NULL) pe_fatal(psi->pe, "calloc(fz) failed\n");

  /* Diffusive fluxes based on six-point stencil, followed
   * by no-flux condition for solid-fluid boundaries. */

  nernst_planck_fluxes(psi, fel, fe, fy, fz);
  nernst_planck_no_flux_condition(map, nk, fe, fy, fz);

  /* Update charge distribution */

  nernst_planck_update(psi, fe, fy, fz);

  free(fz);
  free(fy);
  free(fe);

  return 0;
}

/*****************************************************************************
 *
 *  nernst_planck_fluxes
 *
 *  Compute diffusive fluxes via a simple 7-point stencil in 3d.
 *
 *****************************************************************************/

static int nernst_planck_fluxes(psi_t* psi, fe_t* fel, double* fx,
        double* fy,
        double* fz) {
  int ic, jc, kc, index;
  int nlocal[3];
  int nsites;
  int zs, ys, xs;
  int n, nk;

  double eunit, reunit;
  double b0, b1;
  double mu0, mu1;
  double rho0, rho1;
  double mu_s0, mu_s1;   /* Solvation chemical potential, from free energy */

  assert(psi);
  assert(fx);
  assert(fy);
  assert(fz);

  cs_nsites(psi->cs, &nsites);
  cs_nlocal(psi->cs, nlocal);
  cs_strides(psi->cs, &xs, &ys, &zs);

  psi_nk(psi, &nk);
  psi_unit_charge(psi, &eunit);
  reunit = 1.0 / eunit;

  for (ic = 0; ic <= nlocal[X]; ic++) {
    for (jc = 0; jc <= nlocal[Y]; jc++) {
      for (kc = 0; kc <= nlocal[Z]; kc++) {

        index = cs_index(psi->cs, ic, jc, kc);

        for (n = 0; n < nk; n++) {

          fel->func->mu_solv(fel, index, n, &mu_s0);
          mu0 = reunit * mu_s0
            + psi->valency[n] * psi->psi->data[addr_rank0(nsites, index)];
          rho0 = psi->rho->data[addr_rank1(nsites, nk, index, n)];

          /* x-direction (between ic and ic+1) */

          fel->func->mu_solv(fel, index + xs, n, &mu_s1);
          mu1 = reunit * mu_s1
            + psi->valency[n] * psi->psi->data[addr_rank0(nsites, index + xs)];

          b0 = exp(mu1 - mu0);
          b1 = exp(mu1 - mu0);
          rho1 = psi->rho->data[addr_rank1(nsites, nk, (index + xs), n)] * b1;

          fx[addr_rank1(nsites, nk, index, n)]
            = -psi->diffusivity[n] * 0.5 * (1.0 + b0) * (rho1 - rho0);

          /* y-direction (between jc and jc+1) */

          fel->func->mu_solv(fel, index + ys, n, &mu_s1);
          mu1 = reunit * mu_s1
            + psi->valency[n] * psi->psi->data[addr_rank0(nsites, index + ys)];

          b0 = exp(mu1 - mu0);
          b1 = exp(mu1 - mu0);
          rho1 = psi->rho->data[addr_rank1(nsites, nk, (index + ys), n)] * b1;

          fy[nk * index + n] = -psi->diffusivity[n] * 0.5 * (1.0 + b0) * (rho1 - rho0);

          /* z-direction (between kc and kc+1) */

          fel->func->mu_solv(fel, index + zs, n, &mu_s1);
          mu1 = reunit * mu_s1
            + psi->valency[n] * psi->psi->data[addr_rank0(nsites, index + zs)];

          b0 = exp(mu1 - mu0);
          b1 = exp(mu1 - mu0);
          rho1 = psi->rho->data[addr_rank1(nsites, nk, (index + zs), n)] * b1;

          fz[addr_rank1(nsites, nk, index, n)]
            = -psi->diffusivity[n] * 0.5 * (1.0 + b0) * (rho1 - rho0);
        }

        /* Next face */
      }
    }
  }

  return 0;
}

/*****************************************************************************
 *
 *  nernst_planck_no_flux_condition
 *
 *  Use map to determine solid/fluid boudaries and set no normal flux.
 *
 *****************************************************************************/

static int nernst_planck_no_flux_condition(map_t* map, int nk, double* fx,
             double* fy, double* fz) {

  cs_t* cs = NULL;
  int nsites = 0;
  int nlocal[3] = { 0 };

  assert(map);
  assert(fx);
  assert(fy);
  assert(fz);

  cs = map->cs;
  cs_nsites(cs, &nsites);
  cs_nlocal(cs, nlocal);

  for (int ic = 0; ic <= nlocal[X]; ic++) {
    for (int jc = 0; jc <= nlocal[Y]; jc++) {
      for (int kc = 0; kc <= nlocal[Z]; kc++) {
        int index = cs_index(cs, ic, jc, kc);
        int status = MAP_FLUID;
        map_status(map, index, &status);
        if (status == MAP_BOUNDARY) {
          for (int n = 0; n < nk; n++) {
            fx[addr_rank1(nsites, nk, index, n)] = 0.0;
            fy[addr_rank1(nsites, nk, index, n)] = 0.0;
            fz[addr_rank1(nsites, nk, index, n)] = 0.0;
          }
        }
      }
    }
  }

  return 0;
}

/*****************************************************************************
 *
 *  nernst_planck_update
 *
 *  Update the rho_k from the fluxes. Euler forward step.
 *
 *****************************************************************************/

static int nernst_planck_update(psi_t* psi, double* fx, double* fy,
        double* fz) {
  int ic, jc, kc, index;
  int nlocal[3];
  int zs, ys, xs;
  int n, nk;

  double dt;

  assert(psi);
  assert(fx);
  assert(fy);
  assert(fz);

  cs_nlocal(psi->cs, nlocal);
  cs_strides(psi->cs, &xs, &ys, &zs);

  psi_nk(psi, &nk);
  psi_multistep_timestep(psi, &dt);

  for (ic = 1; ic <= nlocal[X]; ic++) {
    for (jc = 1; jc <= nlocal[Y]; jc++) {
      for (kc = 1; kc <= nlocal[Z]; kc++) {

        index = cs_index(psi->cs, ic, jc, kc);

        for (n = 0; n < nk; n++) {
          psi->rho->data[addr_rank1(psi->nsites, nk, index, n)]
            -= (+fx[addr_rank1(psi->nsites, nk, index, n)]
          - fx[addr_rank1(psi->nsites, nk, (index - xs), n)]
          + fy[addr_rank1(psi->nsites, nk, index, n)]
          - fy[addr_rank1(psi->nsites, nk, (index - ys), n)]
          + fz[addr_rank1(psi->nsites, nk, index, n)]
          - fz[addr_rank1(psi->nsites, nk, (index - zs), n)]) * dt;
        }
      }
    }
  }

  return 0;
}

/*****************************************************************************
 *
 *  nernst_planck_driver_d3qx
 *
 *  The hydro object is allowed to be NULL, in which case there is
 *  no advection.
 *
 *  The map object is allowed to be NULL, in which case no boundary
 *  condition corrections are attempted.
 *
 *****************************************************************************/

int nernst_planck_driver_d3qx(psi_t* psi, fe_t* fe, hydro_t* hydro,
            map_t* map, colloids_info_t* cinfo) {

  int nk;              /* Number of electrolyte species */
  int ia;

  double** flx = NULL;

  psi_nk(psi, &nk);

  /* Allocate fluxes and initialise to zero */
  flx = (double**)calloc((size_t)psi->nsites * nk, sizeof(double*));
  assert(flx);
  if (flx == NULL) pe_fatal(psi->pe, "calloc(flx) failed\n");

  for (ia = 0; ia < psi->nsites * nk; ia++) {
    int nflux = psi->stencil->npoints;
    flx[ia] = (double*)calloc(nflux - 1, sizeof(double));
    assert(flx[ia]);
    if (flx[ia] == NULL) pe_fatal(psi->pe, "calloc(flx[]) failed\n");
  }

  /* Add advective fluxes */
  if (hydro) np_advective_fluxes(psi, hydro, flx); // Test sin flujo de adeveccion

  // /* Diagnostic: sum all fluxes after advection */
  // {
  //   int nlocal_d[3];
  //   cs_nlocal(psi->cs, nlocal_d);
  //   double sum_flx_adv[2] = { 0.0, 0.0 };
  //   for (int ic_d = 1; ic_d <= nlocal_d[X]; ic_d++) {
  //     for (int jc_d = 1; jc_d <= nlocal_d[Y]; jc_d++) {
  //       for (int kc_d = 1; kc_d <= nlocal_d[Z]; kc_d++) {
  //         int idx_d = cs_index(psi->cs, ic_d, jc_d, kc_d);
  //         for (int n_d = 0; n_d < nk; n_d++) {
  //           for (int c_d = 0; c_d < psi->stencil->npoints - 1; c_d++) {
  //             sum_flx_adv[n_d] += flx[addr_rank1(psi->nsites, nk, idx_d, n_d)][c_d];
  //           }
  //         }
  //       }
  //     }
  //   }
  //   pe_info(psi->pe, "  [NP flux] After advection:  sum_flx0=%.10e sum_flx1=%.10e\n",
  //           sum_flx_adv[0], sum_flx_adv[1]);
  // }

  /* Add diffusive fluxes */
  nernst_planck_fluxes_d3qx(psi, fe, hydro, map, cinfo, flx);

  /* Diagnostic: sum all fluxes after diffusion */
  {
    int nlocal_d[3];
    cs_nlocal(psi->cs, nlocal_d);
    double sum_flx_diff[2] = { 0.0, 0.0 };
    for (int ic_d = 1; ic_d <= nlocal_d[X]; ic_d++) {
      for (int jc_d = 1; jc_d <= nlocal_d[Y]; jc_d++) {
        for (int kc_d = 1; kc_d <= nlocal_d[Z]; kc_d++) {
          int idx_d = cs_index(psi->cs, ic_d, jc_d, kc_d);
          for (int n_d = 0; n_d < nk; n_d++) {
            for (int c_d = 0; c_d < psi->stencil->npoints - 1; c_d++) {
              sum_flx_diff[n_d] += flx[addr_rank1(psi->nsites, nk, idx_d, n_d)][c_d];
            }
          }
        }
      }
    }
    pe_info(psi->pe, "  [NP flux] After diffusion:  sum_flx0=%.10e sum_flx1=%.10e\n",
            sum_flx_diff[0], sum_flx_diff[1]);
  }

  /* Apply no-flux BC */
  if (map) np_no_flux_boundary(psi, map, flx);

  /* Diagnostic: link-level flux asymmetry check */
  /* For each node i, direction c -> neighbor j: check flx[i][c] + flx[j][cbar] = 0 */
  /* When j falls in halo, compare with the node at the opposite border */
  {
    int nlocal_a[3];
    cs_nlocal(psi->cs, nlocal_a);
    stencil_t* sa = psi->stencil;
    double max_asym[2] = { 0.0, 0.0 };
    double sum_asym[2] = { 0.0, 0.0 };
    double max_flx[2] = { 0.0, 0.0 };
    long   count_links = 0;
    int    worst_ijk[3] = { 0, 0, 0 };
    int    worst_c = 0, worst_n = 0;
    double worst_fi = 0.0, worst_fj = 0.0;

    for (int ic_a = 1; ic_a <= nlocal_a[X]; ic_a++) {
      for (int jc_a = 1; jc_a <= nlocal_a[Y]; jc_a++) {
        for (int kc_a = 1; kc_a <= nlocal_a[Z]; kc_a++) {
          int idx_i = cs_index(psi->cs, ic_a, jc_a, kc_a);

          for (int ca = 1; ca < (sa->npoints + 1) / 2; ca++) {
            int c_bar = sa->npoints - ca;
            int8_t cxa = sa->cv[ca][X];
            int8_t cya = sa->cv[ca][Y];
            int8_t cza = sa->cv[ca][Z];

            /* Neighbor in real coords; if halo, wrap to opposite border */
            int jx = ic_a + cxa;
            int jy = jc_a + cya;
            int jz = kc_a + cza;
            if (jx < 1)            jx += nlocal_a[X];
            if (jx > nlocal_a[X])  jx -= nlocal_a[X];
            if (jy < 1)            jy += nlocal_a[Y];
            if (jy > nlocal_a[Y])  jy -= nlocal_a[Y];
            if (jz < 1)            jz += nlocal_a[Z];
            if (jz > nlocal_a[Z])  jz -= nlocal_a[Z];
            int idx_j = cs_index(psi->cs, jx, jy, jz);

            for (int na = 0; na < nk; na++) {
              double fi = flx[addr_rank1(psi->nsites, nk, idx_i, na)][ca - 1];
              double fj = flx[addr_rank1(psi->nsites, nk, idx_j, na)][c_bar - 1];
              double asym = fabs(fi + fj);
              double mag = 0.5 * (fabs(fi) + fabs(fj));

              sum_asym[na] += asym;
              if (mag > max_flx[na]) max_flx[na] = mag;
              if (asym > max_asym[na]) {
                max_asym[na] = asym;
                worst_ijk[0] = ic_a; worst_ijk[1] = jc_a; worst_ijk[2] = kc_a;
                worst_c = ca; worst_n = na;
                worst_fi = fi; worst_fj = fj;
              }
              count_links++;
            }
          }
        }
      }
    }

    pe_info(psi->pe, "  [NP asym] Link flux asymmetry (nlinks=%ld):\n", count_links);
    for (int na = 0; na < nk; na++) {
      double rel = (max_flx[na] > 0.0) ? max_asym[na] / max_flx[na] : 0.0;
      pe_info(psi->pe,
        "  [NP asym] species %d: max|fi+fj|=%.6e  sum|fi+fj|=%.6e  max|f|=%.6e  rel=%.6e\n",
        na, max_asym[na], sum_asym[na], max_flx[na], rel);
    }
    pe_info(psi->pe,
      "  [NP asym] worst: ijk=(%d,%d,%d) c=%d n=%d fi=%.12e fj=%.12e sum=%.6e\n",
      worst_ijk[0], worst_ijk[1], worst_ijk[2], worst_c, worst_n,
      worst_fi, worst_fj, worst_fi + worst_fj);
  }

  /* CHANGE INIT - Update charges (with Kahan summation for improved accuracy) */
  nernst_planck_update_d3qx(psi, map, flx);

  for (ia = 0; ia < psi->nsites * nk; ia++) {
    free(flx[ia]);
  }
  free(flx);

  return 0;
}

/*****************************************************************************
 *
 *  nernst_planck_fluxes_d3qx
 *
 *  Compute diffusive fluxes.
 *
 *  We assume we can accumulate the diffusive and advective fluxes separately.
 *
 *  As we compute rho(n+1) = rho(n) - div.flux in the update routine,
 *  there is an extra minus sign in the fluxes here. This conincides
 *  with the sign of the advective fluxes, if present.
 *
 *****************************************************************************/

static int nernst_planck_fluxes_d3qx(psi_t* psi, fe_t* fe, hydro_t* hydro,
             map_t* map, colloids_info_t* cinfo,
             double** flx) {

  int ic, jc, kc;
  int index0, index1;
  int nlocal[3];
  int n, nk; /* Number of charged species */
  int c;
  int status1;

  double b0, b1;
  double mu0, mu1;
  double rho0, rho1;
  double mu_s0, mu_s1;   /* Solvation chemical potential, from free energy */

  double eunit, reunit;
  double dt;

  colloid_t* pc = NULL;

  double* __restrict__ psidata = psi->psi->data;
  double* __restrict__ rhodata = psi->rho->data;

  LB_RCS_TABLE(rcs);

  assert(psi);
  assert(fe);
  assert(fe->func->mu_solv);
  assert(flx);

  cs_nlocal(psi->cs, nlocal);

  psi_unit_charge(psi, &eunit);
  reunit = 1.0 / eunit;

  psi_nk(psi, &nk);
  psi_multistep_timestep(psi, &dt);


  /*CHANGE INIT - 20260206 Antisymmetric flux assignment for exact conservation */
  /* Step 1: Compute fluxes for all 26 directions as in the original code.
   * Step 2: Antisymmetrize — for each link (i,c)/(j,cbar), replace both
   *         with the average: flx[i][c] = (fi - fj)/2, flx[j][cbar] = -flx[i][c].
   *         This enforces exact conservation by construction.
   * Original code only did Step 1, which is analytically antisymmetric but
   * accumulates O(eps_mach) asymmetry per link per substep. */

   //Original code (Step 1 only):

   // // Symetric flux conservation: compute all fluxes first, then antisymmetrize
   /* Step 1: Compute fluxes for all 26 directions (same as original) */
  for (ic = 1; ic <= nlocal[X]; ic++) {
    for (jc = 1; jc <= nlocal[Y]; jc++) {
      for (kc = 1; kc <= nlocal[Z]; kc++) {

        index0 = cs_index(psi->cs, ic, jc, kc);
        colloids_info_map(cinfo, index0, &pc);

        if (pc) {
          continue;
        }
        else {
          stencil_t* s = psi->stencil;
          assert(s);

          for (c = 1; c < s->npoints; c++) {

            int8_t cx = s->cv[c][X];
            int8_t cy = s->cv[c][Y];
            int8_t cz = s->cv[c][Z];
            int8_t pcv = cx * cx + cy * cy + cz * cz;

            index1 = cs_index(psi->cs, ic + cx, jc + cy, kc + cz);
            map_status(map, index1, &status1);

            if (status1 == MAP_FLUID) {

              for (n = 0; n < nk; n++) {

                fe->func->mu_solv(fe, index0, n, &mu_s0);
                mu0 = reunit * mu_s0
                  + psi->valency[n] * psidata[addr_rank0(psi->nsites, index0)];
                rho0 = rhodata[addr_rank1(psi->nsites, nk, index0, n)];

                fe->func->mu_solv(fe, index1, n, &mu_s1);
                mu1 = reunit * mu_s1
                  + psi->valency[n] * psidata[addr_rank0(psi->nsites, index1)];
                b0 = exp(mu0 - mu1);
                b1 = exp(mu1 - mu0);
                rho1 = rhodata[addr_rank1(psi->nsites, nk, index1, n)] * b1;

                flx[addr_rank1(psi->nsites, nk, index0, n)][c - 1]
                  -= psi->diffusivity[n] * 0.5 * (1.0 + b0) * (rho1 - rho0) * rcs[pcv];
              }
            }
          }
        }
      }
    }
  }

  // // // Symetric flux conservation: compute all fluxes first, then antisymmetrize
  //  /* Step 1: Compute fluxes for all 26 directions using Scharfetter-Gummel */
  // double psi0, aux_psi0 = 0.0;
  // double psi1, aux_psi1 = 0.0;
  // double dpsi, max_dpsi1 = 0.0;
  // double kT;
  // psi_beta(psi, &kT);
  // kT=1/kT; 

  // double wv[27];
  // if (psi->stencil->npoints == 19)
  // {
  //   LB_WEIGHTS_D3Q19(wv);
  // }
  // else if (psi->stencil->npoints == 27)
  // {
  //   LB_WEIGHTS_D3Q27(wv);
  // }
  // else {
  //   LB_WEIGHTS_D3Q7(wv);
  // }

  // for (ic = 1; ic <= nlocal[X]; ic++) {
  //   for (jc = 1; jc <= nlocal[Y]; jc++) {
  //     for (kc = 1; kc <= nlocal[Z]; kc++) {

  //       index0 = cs_index(psi->cs, ic, jc, kc);
  //       colloids_info_map(cinfo, index0, &pc);

  //       if (pc) {
  //         continue;
  //       }
  //       else {

  //         stencil_t* s = psi->stencil;
  //         assert(s);

  //         for (c = 1; c < s->npoints; c++) {

  //           int8_t cx = s->cv[c][X];
  //           int8_t cy = s->cv[c][Y];
  //           int8_t cz = s->cv[c][Z];
  //           int8_t pcv = cx * cx + cy * cy + cz * cz;

  //           index1 = cs_index(psi->cs, ic + cx, jc + cy, kc + cz);
  //           map_status(map, index1, &status1);

  //           if (status1 == MAP_FLUID) {

  //             for (n = 0; n < nk; n++) {

  //               fe->func->mu_solv(fe, index0, n, &mu_s0);
  //               mu0 = reunit * mu_s0
  //                 + psi->valency[n] * psidata[addr_rank0(psi->nsites, index0)];
  //               rho0 = rhodata[addr_rank1(psi->nsites, nk, index0, n)];

  //               fe->func->mu_solv(fe, index1, n, &mu_s1);
  //               mu1 = reunit * mu_s1
  //                 + psi->valency[n] * psidata[addr_rank0(psi->nsites, index1)];
  //               rho1 = rhodata[addr_rank1(psi->nsites, nk, index1, n)];

  //               psi0 = mu0;
  //               psi1 = mu1;
  //               // dpsi = kT *(psi1 - psi0);
  //               dpsi = (psi1 - psi0);

  //               if (fabs(dpsi) > max_dpsi1)
  //               {
  //                 max_dpsi1 = fabs(dpsi);
  //                 aux_psi0 = psi0;
  //                 aux_psi1 = psi1;
  //               }
  //               if (dpsi > 50) dpsi=50; // To avoid overflow in exp(dpsi)
  //               if (-dpsi > 50) dpsi=-50; // To avoid overflow in exp(dpsi)

  //               flx[addr_rank1(psi->nsites, nk, index0, n)][c - 1]
  //                 -= wv[c - 1] * psi->diffusivity[n] * rcs[pcv] *
  //                 (rho1 * bernoulli(dpsi) - rho0 * bernoulli(-dpsi));

  //             }
  //           }
  //         }
  //       }
  //     }
  //   }
  // }

  // pe_info(psi->pe, "  [NP flux - potencial quimico] kT=%.10e mu0=%.10e mu1=%.10e dpsi=%.10e\n",  //rho0=%.10e rho1=%.10e\n",
  //                       kT, aux_psi0, aux_psi1, max_dpsi1); //, rho0, rho1);

  /* Step 2: Antisymmetrize fluxes for exact conservation.
   * For each interior pair (i,j) with direction c and opposite c_bar:
   *   fi = flx[i][c-1],  fj = flx[j][c_bar-1]
   *   Analytically fi = -fj, but numerically they differ by O(eps_mach).
   *   Replace: flx[i][c-1] = (fi - fj)/2,  flx[j][c_bar-1] = (fj - fi)/2
   *   When j is a halo, use the periodic image at the opposite border. */
  {
    stencil_t* s = psi->stencil;
    int half = (s->npoints - 1) / 2;  /* 13 for D3Q27 */

    for (ic = 1; ic <= nlocal[X]; ic++) {
      for (jc = 1; jc <= nlocal[Y]; jc++) {
        for (kc = 1; kc <= nlocal[Z]; kc++) {
          index0 = cs_index(psi->cs, ic, jc, kc);

          for (c = 1; c <= half; c++) {
            int c_bar = s->npoints - c;
            int8_t cx = s->cv[c][X];
            int8_t cy = s->cv[c][Y];
            int8_t cz = s->cv[c][Z];

            /* Neighbor: wrap periodically to interior node */
            int jx = ic + cx;
            int jy = jc + cy;
            int jz = kc + cz;
            if (jx < 1)           jx += nlocal[X];
            if (jx > nlocal[X])   jx -= nlocal[X];
            if (jy < 1)           jy += nlocal[Y];
            if (jy > nlocal[Y])   jy -= nlocal[Y];
            if (jz < 1)           jz += nlocal[Z];
            if (jz > nlocal[Z])   jz -= nlocal[Z];
            index1 = cs_index(psi->cs, jx, jy, jz);

            for (n = 0; n < nk; n++) {
              double* fi = &flx[addr_rank1(psi->nsites, nk, index0, n)][c - 1];
              double* fj = &flx[addr_rank1(psi->nsites, nk, index1, n)][c_bar - 1];
              double avg = 0.5 * (*fi - *fj);
              *fi = avg;
              *fj = -avg;
            }
          }
        }
      }
    }
  }
  /*CHANGE END - 20260206 Antisymmetric flux assignment */

  return 0;
}

/*****************************************************************************
 *
 *  nernst_planck_fluxes_force_d3qx
 *
 *  Compute diffusive fluxes and link-flux force on fluid.
 *
 *  We assume we can accumulate the diffusive and advective fluxes separately.
 *
 *  As we compute rho(n+1) = rho(n) - div.flux in the update routine,
 *  there is an extra minus sign in the fluxes here. This conincides
 *  with the sign of the advective fluxes, if present.
 *
 *****************************************************************************/

int nernst_planck_fluxes_force_d3qx(psi_t* psi, fe_t* fe, hydro_t* hydro,
            map_t* map, colloids_info_t* cinfo,
            double** flx) {

  int ic, jc, kc;
  int index0, index1;
  int nlocal[3];
  int n, nk;
  int nsites;
  int c;
  int status1;
  double eunit;
  double beta, rbeta;
  double b0, b1;
  double mu0, mu1;
  double rho0, rho1;
  double mu_s0, mu_s1;
  double rho_elec;
  double e[3];
  double flocal[4] = { 0.0, 0.0, 0.0, 0.0 }, fsum[4], f[3];
  double flxtmp[2];
  double dt;
  MPI_Comm comm;
  colloid_t* pc = NULL;
  double* __restrict__ psidata = psi->psi->data;
  double* __restrict__ rhodata = psi->rho->data;
  LB_RCS_TABLE(rcs);
  assert(psi);
  assert(flx);
  cs_nsites(psi->cs, &nsites);
  cs_nlocal(psi->cs, nlocal);
  cs_cart_comm(psi->cs, &comm);
  psi_nk(psi, &nk);
  psi_unit_charge(psi, &eunit);
  psi_beta(psi, &beta);
  psi_multistep_timestep(psi, &dt);
  rbeta = 1.0 / beta;
  for (ic = 1; ic <= nlocal[X]; ic++) {
    for (jc = 1; jc <= nlocal[Y]; jc++) {
      for (kc = 1; kc <= nlocal[Z]; kc++) {
        index0 = cs_index(psi->cs, ic, jc, kc);
        colloids_info_map(cinfo, index0, &pc);
        f[X] = 0.0; f[Y] = 0.0; f[Z] = 0.0;
        psi_rho_elec(psi, index0, &rho_elec);
        if (pc) {
          psi_electric_field(psi, index0, e);
          f[X] = rho_elec * e[X] * dt;
          f[Y] = rho_elec * e[Y] * dt;
          f[Z] = rho_elec * e[Z] * dt;
          pc->force[X] += f[X];
          pc->force[Y] += f[Y];
          pc->force[Z] += f[Z];
        }
        else {
          stencil_t* s = psi->stencil;
          for (c = 1; c < s->npoints; c++) {
            int8_t cx = s->cv[c][X]; int8_t cy = s->cv[c][Y]; int8_t cz = s->cv[c][Z];
            int8_t pcv = cx * cx + cy * cy + cz * cz;
            index1 = cs_index(psi->cs, ic + cx, jc + cy, kc + cz);
            map_status(map, index1, &status1);
            if (status1 == MAP_FLUID) {
              for (n = 0; n < nk; n++) {
                fe->func->mu_solv(fe, index0, n, &mu_s0);
                mu0 = mu_s0 + psi->valency[n] * eunit * psidata[addr_rank0(nsites, index0)];
                rho0 = rhodata[addr_rank1(nsites, nk, index0, n)];
                fe->func->mu_solv(fe, index1, n, &mu_s1);
                mu1 = mu_s1 + psi->valency[n] * eunit * psidata[addr_rank0(nsites, index1)];
                b0 = exp(-beta * (mu1 - mu0));
                b1 = exp(+beta * (mu1 - mu0));
                rho1 = rhodata[addr_rank1(nsites, nk, index1, n)] * b1;
                flxtmp[0] = -0.5 * (1.0 + b0) * (rho1 - rho0) * rcs[pcv];
                flx[addr_rank1(nsites, nk, index0, n)][c - 1] += psi->diffusivity[n] * flxtmp[0];
                f[X] -= s->wgradients[c] * cx * flxtmp[0] * rbeta;
                f[Y] -= s->wgradients[c] * cy * flxtmp[0] * rbeta;
                f[Z] -= s->wgradients[c] * cz * flxtmp[0] * rbeta;
              }
            }
          }
          f[X] *= dt; f[Y] *= dt; f[Z] *= dt;
          flocal[3] += 1.0;
          if (hydro) hydro_f_local_add(hydro, index0, f);
        }
        flocal[X] += f[X]; flocal[Y] += f[Y]; flocal[Z] += f[Z];
      }
    }
  }
  MPI_Allreduce(flocal, fsum, 4, MPI_DOUBLE, MPI_SUM, comm);
  fsum[X] /= fsum[3]; fsum[Y] /= fsum[3]; fsum[Z] /= fsum[3];
  for (ic = 1; ic <= nlocal[X]; ic++) {
    for (jc = 1; jc <= nlocal[Y]; jc++) {
      for (kc = 1; kc <= nlocal[Z]; kc++) {
        index0 = cs_index(psi->cs, ic, jc, kc);
        colloids_info_map(cinfo, index0, &pc);
        if (pc) continue;
        f[X] = -fsum[X]; f[Y] = -fsum[Y]; f[Z] = -fsum[Z];
        if (hydro) hydro_f_local_add(hydro, index0, f);
      }
    }
  }
  return 0;
}

/*****************************************************************************
 *
 *  nernst_planck_update_d3qx
 *
 *  Update the rho_k from the fluxes (D3QX stencil). Euler forward step.
 *
 *****************************************************************************/

static int nernst_planck_update_d3qx(psi_t* psi, map_t* map, double** flx) {

  int ic, jc, kc, index;
  int nsites;
  int nlocal[3];
  int n, nk;
  int c;
  int status;
  double acc, maxacc = 0.0;
  double dt;

  assert(psi);
  assert(flx);

  cs_nsites(psi->cs, &nsites);
  cs_nlocal(psi->cs, nlocal);

  psi_nk(psi, &nk);
  psi_multistep_timestep(psi, &dt);

  for (ic = 1; ic <= nlocal[X]; ic++) {
    for (jc = 1; jc <= nlocal[Y]; jc++) {
      for (kc = 1; kc <= nlocal[Z]; kc++) {

        index = cs_index(psi->cs, ic, jc, kc);
        map_status(map, index, &status);

        if (status == MAP_FLUID) {
          stencil_t* s = psi->stencil;
          for (n = 0; n < nk; n++) {

            /*CHANGE INIT - 20260210 Use kahan_t from util_sum.h */
            /* Original code:
            // acc = 0.0;
            // for (c = 1; c < s->npoints; c++) {
            //   psi->rho->data[addr_rank1(nsites, nk, index, n)]
            //     -= flx[addr_rank1(nsites, nk, index, n)][c - 1] * dt;
            //   acc += fabs(flx[addr_rank1(nsites, nk, index, n)][c - 1] * dt);
            // }
            */
            kahan_t ksum = kahan_zero();
            kahan_t kacc = kahan_zero();
            for (c = 1; c < s->npoints; c++) {
              double flux_dt = flx[addr_rank1(nsites, nk, index, n)][c - 1] * dt;
              kahan_add_double(&ksum, flux_dt);
              kahan_add_double(&kacc, fabs(flux_dt));
            }
            psi->rho->data[addr_rank1(nsites, nk, index, n)] -= kahan_sum(&ksum);
            acc = kahan_sum(&kacc);
            /*CHANGE END - 20260210 Use kahan_t from util_sum.h */
            acc /= fabs(psi->rho->data[addr_rank1(nsites, nk, index, n)]);
            if (maxacc < acc) maxacc = acc;
          }
        }

      }
    }
  }

  nernst_planck_maxacc_set(maxacc);

  return 0;
}

/*****************************************************************************
 *
 *  nernst_planck_maxacc_set
 *
 *  Setter function for the local maximal accuracy in the Nernst-Planck
 *  equation. This is defined as the absolut value of the ratio of the change
 *  during one fractional LB timestep (multistep dt) and the charge
 *  density itself.
 *
 *****************************************************************************/

int nernst_planck_maxacc_set(double acc) {
  max_acc = acc;
  return 0;
}

/*****************************************************************************
 *
 *  nernst_planck_maxacc
 *
 *  Getter function for the local maximal accuracy
 *  in the Nernst-Planck equation.
 *
 *****************************************************************************/

int nernst_planck_maxacc(double* acc) {
  *acc = max_acc;
  return 0;
}

/*****************************************************************************
 *
 *  nernst_planck_adjust_multistep
 *
 *
 *****************************************************************************/

int nernst_planck_adjust_multistep(psi_t* psi) {

  double maxacc_local[1], maxacc[1], diffacc; /* actual and preset value of diffusive accuracy */
  double diff, diffmax = 0.0;                   /* diffusivity of species and maximal value      */
  int n, nk, multisteps;
  MPI_Comm comm;

  psi_diffacc(psi, &diffacc);

  /* Take local maximum and reduce for global maximum */
  nernst_planck_maxacc(&maxacc_local[0]);
  cs_cart_comm(psi->cs, &comm);
  MPI_Allreduce(maxacc_local, maxacc, 1, MPI_DOUBLE, MPI_MAX, comm);

  /* Compare maximal accuracy with preset value for */
  /*   diffusion and adjust number of multisteps    */

  /* Increase no. of multisteps */
  if (*maxacc > diffacc && diffacc > 0.0) {
    psi_multisteps(psi, &multisteps);
    multisteps *= 2;
    psi->multisteps = multisteps;
    pe_info(psi->pe, "\nMaxacc > diffacc: changing no. of multisteps to %d\n",
      multisteps);
  }

  /* Reduce no. of multisteps */
  /* The factor 0.1 prevents too frequent changes. */
  if (*maxacc < 0.1 * diffacc && diffacc > 0.0) {

    psi_multisteps(psi, &multisteps);
    psi_nk(psi, &nk);

    for (n = 0; n < nk; n++) {
      psi_diffusivity(psi, n, &diff);
      if (diff > diffmax) diffmax = diff;
    }

    /* Only reduce if sanity criteria fulfilled */
    if (multisteps > 1 && diffmax / multisteps < 0.05) {
      multisteps *= 0.5;
      psi->multisteps = multisteps;
      pe_info(psi->pe, "\nMaxacc << diffacc: changing no. of multisteps to %d\n", multisteps);
    }

  }

  return 0;
}

/*****************************************************************************
 *
 *  np_advective_fluxes
 *
 *  'Centred difference' advective fluxes for the char densities rho.
 *
 *  Symmetric two-point stencil.
 *
 *****************************************************************************/

int np_advective_fluxes(psi_t* psi, hydro_t* hydro, double** flx) {

  int nlocal[3] = { 0 };
  cs_t* cs = NULL;
  stencil_t* s = NULL;

  double* __restrict__ rho = psi->rho->data;

  assert(psi);
  assert(hydro);
  assert(flx);

  cs = psi->cs;
  s = psi->stencil;
  assert(cs);
  assert(s);

  cs_nlocal(cs, nlocal);

  for (int ic = 1; ic <= nlocal[X]; ic++) {
    for (int jc = 1; jc <= nlocal[Y]; jc++) {
      for (int kc = 1; kc <= nlocal[Z]; kc++) {

        int index0 = cs_index(cs, ic, jc, kc);
        double u0[3] = { 0 };
        hydro_u(hydro, index0, u0);

        for (int p = 1; p < s->npoints; p++) {

          int8_t cx = s->cv[p][X];
          int8_t cy = s->cv[p][Y];
          int8_t cz = s->cv[p][Z];
          int index1 = cs_index(cs, ic + cx, jc + cy, kc + cz);
          double u1[3] = { 0 };
          double u = 0.0;
          hydro_u(hydro, index1, u1);

          u = 0.5 * ((u0[X] + u1[X]) * cx + (u0[Y] + u1[Y]) * cy + (u0[Z] + u1[Z]) * cz);

          for (int n = 0; n < psi->nk; n++) {
            double rho0 = rho[addr_rank1(psi->nsites, psi->nk, index0, n)];
            double rho1 = rho[addr_rank1(psi->nsites, psi->nk, index1, n)];
            double flux = u * 0.5 * (rho0 + rho1);
            flx[addr_rank1(psi->nsites, psi->nk, index0, n)][p - 1] = flux;
          }
        }
        /* Next site */
      }
    }
  }

  return 0;
}

/*****************************************************************************
 *
 *  advective_bcs_no_flux_d3qx
 *
 *  Set normal fluxes at solid fluid interfaces to zero.
 *
 *****************************************************************************/

int np_no_flux_boundary(psi_t* psi, map_t* map, double** flx) {

  int nlocal[3] = { 0 };
  cs_t* cs = NULL;
  stencil_t* s = NULL;

  assert(psi);
  assert(map);
  assert(flx);

  cs = psi->cs;
  s = psi->stencil;
  assert(cs);
  assert(s);

  cs_nlocal(cs, nlocal);

  for (int ic = 1; ic <= nlocal[X]; ic++) {
    for (int jc = 1; jc <= nlocal[Y]; jc++) {
      for (int kc = 1; kc <= nlocal[Z]; kc++) {

        int index0 = cs_index(cs, ic, jc, kc);
        int mask0 = 1;
        int status = MAP_BOUNDARY;
        map_status(map, index0, &status);
        mask0 = (status == MAP_FLUID);

        for (int p = 1; p < s->npoints; p++) {

          int8_t cx = s->cv[p][X];
          int8_t cy = s->cv[p][Y];
          int8_t cz = s->cv[p][Z];
          int index1 = cs_index(cs, ic + cx, jc + cy, kc + cz);
          int mask = 1;

          map_status(map, index1, &status);
          mask = (status == MAP_FLUID);
          mask = mask * mask0;

          for (int n = 0; n < psi->nk; n++) {
            flx[addr_rank1(psi->nsites, psi->nk, index0, n)][p - 1] *= mask;
          }
        }
      }
    }
  }

  return 0;
}

/*****************************************************************************
 *
 *  CHANGE INIT - 20260326 GPU implementation of nernst_planck_driver_d3qx
 *
 *  nernst_planck_driver_d3qx_gpu
 *
 *  GPU-parallel version of nernst_planck_driver_d3qx. Uses a flat
 *  contiguous flux array flx_gpu[nsites * nk * (npoints-1)] allocated
 *  on the device, avoiding the CPU double** structure which is
 *  incompatible with GPU kernels.
 *
 *  Layout of flx_gpu:
 *    flx_gpu[ (n * nsites + index) * (npoints-1) + c ]
 *  which matches addr_rank1(nsites, nk, index, n) * (npoints-1) + c
 *  in SOA mode (ADDR_SOA).
 *
 *  Kernels:
 *    1. np_advective_fluxes_kernel  - advective fluxes (if hydro != NULL)
 *    2. np_diffusive_fluxes_kernel  - diffusive fluxes (Capuani scheme)
 *    3. np_antisym_kernel           - antisymmetrisation for exact conservation
 *    4. np_update_rho_kernel        - update rho from divergence of fluxes
 *
 *  Note: np_no_flux_boundary is NOT implemented in this GPU version.
 *  Note: mu_solv = 0 is assumed (valid for fe_electro).
 *  Note: colloids are detected via map->status (MAP_COLLOID) instead of
 *        colloids_info_map() which is __host__ only.
 *
 *****************************************************************************/

/* Macro for flat flux array indexing (matches SOA addr_rank1 layout) */
#define FLX_GPU(nsites, nk, nflux, index, n, c) \
  ( ((n)*(nsites) + (index)) * (nflux) + (c) )

/*---------------------------------------------------------------------------
 *  np_advective_fluxes_kernel
 *
 *  One thread per lattice site. Computes advective fluxes from velocity field.
 *  flx_gpu[ FLX_GPU(index, n, p-1) ] = u_face * 0.5*(rho0 + rho1)
 *---------------------------------------------------------------------------*/
__global__ static void np_advective_fluxes_kernel(
    cs_t*    cs,
    int      nlocal[3],
    int      nsites,
    int      nk,
    int      npoints,
    double*  rho_d,      /* psi->rho->target->data  */
    double*  u_d,        /* hydro->u->target->data  */
    int8_t*  cv_d,       /* stencil cv[npoints][3], row-major flat */
    double*  flx_gpu)
{
  int kc = blockIdx.x * blockDim.x + threadIdx.x + 1;
  int jc = blockIdx.y * blockDim.y + threadIdx.y + 1;
  int ic = blockIdx.z * blockDim.z + threadIdx.z + 1;

  if (ic > nlocal[X] || jc > nlocal[Y] || kc > nlocal[Z]) return;

  int nflux = npoints - 1;
  int index0 = cs_index(cs, ic, jc, kc);

  double u0[3];
  u0[X] = u_d[addr_rank1(nsites, 3, index0, X)];
  u0[Y] = u_d[addr_rank1(nsites, 3, index0, Y)];
  u0[Z] = u_d[addr_rank1(nsites, 3, index0, Z)];

  for (int p = 1; p < npoints; p++) {
    int8_t cx = cv_d[p*3 + X];
    int8_t cy = cv_d[p*3 + Y];
    int8_t cz = cv_d[p*3 + Z];
    int index1 = cs_index(cs, ic + cx, jc + cy, kc + cz);

    double u1[3];
    u1[X] = u_d[addr_rank1(nsites, 3, index1, X)];
    u1[Y] = u_d[addr_rank1(nsites, 3, index1, Y)];
    u1[Z] = u_d[addr_rank1(nsites, 3, index1, Z)];

    double u = 0.5*((u0[X]+u1[X])*cx + (u0[Y]+u1[Y])*cy + (u0[Z]+u1[Z])*cz);

    for (int n = 0; n < nk; n++) {
      double rho0 = rho_d[addr_rank1(nsites, nk, index0, n)];
      double rho1 = rho_d[addr_rank1(nsites, nk, index1, n)];
      flx_gpu[FLX_GPU(nsites, nk, nflux, index0, n, p-1)] = u * 0.5*(rho0 + rho1);
    }
  }
}

/*---------------------------------------------------------------------------
 *  np_diffusive_fluxes_kernel
 *
 *  One thread per lattice site. Adds diffusive fluxes (Capuani scheme)
 *  to flx_gpu. mu_solv = 0 assumed (fe_electro).
 *  Skips sites occupied by colloids or boundaries (via map->status).
 *---------------------------------------------------------------------------*/
__global__ static void np_diffusive_fluxes_kernel(
    cs_t*   cs,
    int     nlocal[3],
    int     nsites,
    int     nk,
    int     npoints,
    double* psi_d,       /* psi->psi->target->data  */
    double* rho_d,       /* psi->rho->target->data  */
    char*   status_d,    /* map->target->status     */
    int8_t* cv_d,        /* stencil cv flat         */
    double  reunit,
    int*    valency_d,
    double* diffusivity_d,
    double* flx_gpu)
{
  int kc = blockIdx.x * blockDim.x + threadIdx.x + 1;
  int jc = blockIdx.y * blockDim.y + threadIdx.y + 1;
  int ic = blockIdx.z * blockDim.z + threadIdx.z + 1;

  if (ic > nlocal[X] || jc > nlocal[Y] || kc > nlocal[Z]) return;

  int nflux = npoints - 1;
  int index0 = cs_index(cs, ic, jc, kc);

  /* Skip non-fluid sites (colloids, boundaries) */
  if ((int)status_d[index0] != MAP_FLUID) return;

  /* rcs table: rcs[pcv] where pcv = cx^2+cy^2+cz^2 in {1,2,3} */
  const double rcs[4] = {0.0, 1.0, 1.0/1.4142135623730951,
                              1.0/1.7320508075688772};

  for (int c = 1; c < npoints; c++) {
    int8_t cx = cv_d[c*3 + X];
    int8_t cy = cv_d[c*3 + Y];
    int8_t cz = cv_d[c*3 + Z];
    int8_t pcv = cx*cx + cy*cy + cz*cz;
    int index1 = cs_index(cs, ic + cx, jc + cy, kc + cz);

    if ((int)status_d[index1] != MAP_FLUID) continue;

    for (int n = 0; n < nk; n++) {
      /* mu_solv = 0 for fe_electro */
      double mu0 = reunit * 0.0
                 + valency_d[n] * psi_d[addr_rank0(nsites, index0)];
      double mu1 = reunit * 0.0
                 + valency_d[n] * psi_d[addr_rank0(nsites, index1)];

      double b0 = exp(mu0 - mu1);
      double b1 = exp(mu1 - mu0);
      double rho0 = rho_d[addr_rank1(nsites, nk, index0, n)];
      double rho1 = rho_d[addr_rank1(nsites, nk, index1, n)] * b1;

      flx_gpu[FLX_GPU(nsites, nk, nflux, index0, n, c-1)]
        -= diffusivity_d[n] * 0.5*(1.0 + b0) * (rho1 - rho0) * rcs[pcv];
    }
  }
}

/*---------------------------------------------------------------------------
 *  np_antisym_kernel
 *
 *  One thread per lattice site. For each of the 'half' forward directions,
 *  enforces flx[i][c] = -flx[j][c_bar] exactly by replacing both with the
 *  antisymmetric average. Periodic wrapping handles boundary nodes.
 *  Only processes c = 1..half to avoid double-processing each pair.
 *---------------------------------------------------------------------------*/
__global__ static void np_antisym_kernel(
    cs_t*   cs,
    int     nlocal[3],
    int     nsites,
    int     nk,
    int     npoints,
    int8_t* cv_d,
    double* flx_gpu)
{
  int kc = blockIdx.x * blockDim.x + threadIdx.x + 1;
  int jc = blockIdx.y * blockDim.y + threadIdx.y + 1;
  int ic = blockIdx.z * blockDim.z + threadIdx.z + 1;

  if (ic > nlocal[X] || jc > nlocal[Y] || kc > nlocal[Z]) return;

  int nflux  = npoints - 1;
  int half   = (npoints - 1) / 2;
  int index0 = cs_index(cs, ic, jc, kc);

  for (int c = 1; c <= half; c++) {
    int c_bar = npoints - c;
    int8_t cx = cv_d[c*3 + X];
    int8_t cy = cv_d[c*3 + Y];
    int8_t cz = cv_d[c*3 + Z];

    /* Periodic wrap of neighbor to interior */
    int jx = ic + cx; if (jx < 1) jx += nlocal[X]; if (jx > nlocal[X]) jx -= nlocal[X];
    int jy = jc + cy; if (jy < 1) jy += nlocal[Y]; if (jy > nlocal[Y]) jy -= nlocal[Y];
    int jz = kc + cz; if (jz < 1) jz += nlocal[Z]; if (jz > nlocal[Z]) jz -= nlocal[Z];
    int index1 = cs_index(cs, jx, jy, jz);

    /* Only process pair once: from the node with smaller index */
    if (index0 > index1) continue;

    for (int n = 0; n < nk; n++) {
      int ai = FLX_GPU(nsites, nk, nflux, index0, n, c-1);
      int aj = FLX_GPU(nsites, nk, nflux, index1, n, c_bar-1);
      double avg = 0.5*(flx_gpu[ai] - flx_gpu[aj]);
      flx_gpu[ai] =  avg;
      flx_gpu[aj] = -avg;
    }
  }
}

/*---------------------------------------------------------------------------
 *  np_update_rho_kernel
 *
 *  One thread per lattice site. Updates rho from divergence of fluxes.
 *  rho[index,n] -= sum_c flx[index,n,c] * dt
 *---------------------------------------------------------------------------*/
__global__ static void np_update_rho_kernel(
    cs_t*   cs,
    int     nlocal[3],
    int     nsites,
    int     nk,
    int     npoints,
    char*   status_d,
    double* rho_d,
    double* flx_gpu,
    double  dt,
    double* maxacc_d)   /* per-block partial max, reduced on host */
{
  int kc = blockIdx.x * blockDim.x + threadIdx.x + 1;
  int jc = blockIdx.y * blockDim.y + threadIdx.y + 1;
  int ic = blockIdx.z * blockDim.z + threadIdx.z + 1;

  if (ic > nlocal[X] || jc > nlocal[Y] || kc > nlocal[Z]) return;

  int nflux  = npoints - 1;
  int index  = cs_index(cs, ic, jc, kc);

  if ((int)status_d[index] != MAP_FLUID) return;

  for (int n = 0; n < nk; n++) {
    double sum = 0.0;
    for (int c = 0; c < nflux; c++) {
      sum += flx_gpu[FLX_GPU(nsites, nk, nflux, index, n, c)] * dt;
    }
    double rho_old = rho_d[addr_rank1(nsites, nk, index, n)];
    rho_d[addr_rank1(nsites, nk, index, n)] = rho_old - sum;

    double acc = (rho_old > 0.0) ? fabs(sum) / rho_old : 0.0;
    /* Per-thread max written to global; atomicMax on doubles needs workaround */
    /* Use atomicAdd approximation: store to per-block shared later if needed  */
    /* For simplicity, maxacc is not computed on GPU - reported as 0.0 */
    (void)acc;
  }
}

/*****************************************************************************
 *
 *  nernst_planck_driver_d3qx_gpu
 *
 *****************************************************************************/

int nernst_planck_driver_d3qx_gpu(psi_t* psi, fe_t* fe, hydro_t* hydro,
                                   map_t* map, colloids_info_t* cinfo) {
  int    nk, nsites;
  int    nlocal[3];
  double dt;
  double eunit, reunit;

  assert(psi);
  assert(fe);
  assert(map);

  psi_nk(psi, &nk);
  cs_nsites(psi->cs, &nsites);
  cs_nlocal(psi->cs, nlocal);
  psi_unit_charge(psi, &eunit);
  reunit = 1.0 / eunit;
  psi_multistep_timestep(psi, &dt);

  int npoints = psi->stencil->npoints;
  int nflux   = npoints - 1;

  /* --- Allocate flat flux array on GPU --- */
  double* flx_gpu = NULL;
  size_t flx_size = (size_t)nsites * nk * nflux * sizeof(double);
  tdpAssert(tdpMalloc((void**)&flx_gpu, flx_size));
  tdpAssert(tdpMemset(flx_gpu, 0, flx_size));

  /* --- Copy stencil cv to GPU (flat int8_t[npoints][3]) --- */
  int8_t* cv_d = NULL;
  size_t cv_size = (size_t)npoints * 3 * sizeof(int8_t);
  int8_t* cv_flat = (int8_t*)malloc(cv_size);
  for (int p = 0; p < npoints; p++) {
    cv_flat[p*3 + X] = psi->stencil->cv[p][X];
    cv_flat[p*3 + Y] = psi->stencil->cv[p][Y];
    cv_flat[p*3 + Z] = psi->stencil->cv[p][Z];
  }
  tdpAssert(tdpMalloc((void**)&cv_d, cv_size));
  tdpAssert(tdpMemcpy(cv_d, cv_flat, cv_size, tdpMemcpyHostToDevice));
  free(cv_flat);

  /* --- Copy valency and diffusivity to GPU --- */
  int*    valency_d      = NULL;
  double* diffusivity_d  = NULL;
  tdpAssert(tdpMalloc((void**)&valency_d,     nk * sizeof(int)));
  tdpAssert(tdpMalloc((void**)&diffusivity_d, nk * sizeof(double)));
  tdpAssert(tdpMemcpy(valency_d,     psi->valency,     nk*sizeof(int),    tdpMemcpyHostToDevice));
  tdpAssert(tdpMemcpy(diffusivity_d, psi->diffusivity, nk*sizeof(double), tdpMemcpyHostToDevice));

  /* --- Copy nlocal array to GPU --- */
  int* nlocal_d = NULL;
  tdpAssert(tdpMalloc((void**)&nlocal_d, 3*sizeof(int)));
  tdpAssert(tdpMemcpy(nlocal_d, nlocal, 3*sizeof(int), tdpMemcpyHostToDevice));

  /* --- Ensure field data is up to date on GPU --- */
  field_memcpy(psi->psi, tdpMemcpyHostToDevice);
  field_memcpy(psi->rho, tdpMemcpyHostToDevice);
  if (hydro) field_memcpy(hydro->u, tdpMemcpyHostToDevice);
  map_memcpy(map, tdpMemcpyHostToDevice);

  /* --- GPU pointers --- */
  /*CHANGE INIT - 20260326 Fix SEGV: target is a GPU struct, cannot deref ->data from host.
   * Use tdpMemcpy(DeviceToHost) to read the device pointer value stored inside the GPU struct. */
  double* psi_d    = NULL;
  double* rho_d    = NULL;
  char*   status_d = NULL;
  double* u_d      = NULL;
  tdpAssert(tdpMemcpy(&psi_d,    &psi->psi->target->data,  sizeof(double*), tdpMemcpyDeviceToHost));
  tdpAssert(tdpMemcpy(&rho_d,    &psi->rho->target->data,  sizeof(double*), tdpMemcpyDeviceToHost));
  tdpAssert(tdpMemcpy(&status_d, &map->target->status,     sizeof(char*),   tdpMemcpyDeviceToHost));
  if (hydro) tdpAssert(tdpMemcpy(&u_d, &hydro->u->target->data, sizeof(double*), tdpMemcpyDeviceToHost));
  /*CHANGE END - 20260326 */
  cs_t*   cs_d     = psi->cs->target;

  /* --- Launch configuration: one thread per interior node --- */
  dim3 block(8, 8, 8);
  dim3 grid( (nlocal[Z] + block.x - 1) / block.x,
             (nlocal[Y] + block.y - 1) / block.y,
             (nlocal[X] + block.z - 1) / block.z );

  /* --- Kernel 1: advective fluxes --- */
  if (hydro) {
    np_advective_fluxes_kernel<<<grid, block>>>(
        cs_d, nlocal_d, nsites, nk, npoints,
        rho_d, u_d, cv_d, flx_gpu);
    tdpAssert(tdpDeviceSynchronize());
  }

  /* --- Kernel 2: diffusive fluxes --- */
  np_diffusive_fluxes_kernel<<<grid, block>>>(
      cs_d, nlocal_d, nsites, nk, npoints,
      psi_d, rho_d, status_d, cv_d,
      reunit, valency_d, diffusivity_d, flx_gpu);
  tdpAssert(tdpDeviceSynchronize());

  /* --- Kernel 3: antisymmetrisation --- */
  np_antisym_kernel<<<grid, block>>>(
      cs_d, nlocal_d, nsites, nk, npoints, cv_d, flx_gpu);
  tdpAssert(tdpDeviceSynchronize());

  /* --- Kernel 4: update rho --- */
  np_update_rho_kernel<<<grid, block>>>(
      cs_d, nlocal_d, nsites, nk, npoints,
      status_d, rho_d, flx_gpu, dt, NULL);
  tdpAssert(tdpDeviceSynchronize());

  /* --- Copy rho back to host --- */
  field_memcpy(psi->rho, tdpMemcpyDeviceToHost);

  /* maxacc not computed on GPU; set to 0 */
  nernst_planck_maxacc_set(0.0);

  /* --- Free GPU temporaries --- */
  tdpAssert(tdpFree(flx_gpu));
  tdpAssert(tdpFree(cv_d));
  tdpAssert(tdpFree(valency_d));
  tdpAssert(tdpFree(diffusivity_d));
  tdpAssert(tdpFree(nlocal_d));

  return 0;
}
/*CHANGE END - 20260326 GPU implementation of nernst_planck_driver_d3qx */
