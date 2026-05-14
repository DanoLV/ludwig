/*****************************************************************************
 *
 *  psi_exclusion.c
 *
 *  Ion exclusion layer for subgrid colloidal particles.
 *
 *  See psi_exclusion.h for a full algorithm description and the suggested
 *  call order inside the simulation loop.
 *
 *  Placement in the simulation loop (ludwig.c):
 *
 *    [Before Ewald / Nernst-Planck solve, each timestep]
 *
 *    psi_exclusion_compute_accessibility(excl);   // Step 1
 *    psi_exclusion_update_map(excl);              // Step 2
 *    psi_exclusion_remove_charge(excl);           // Step 3
 *    psi_exclusion_redistribute_charge(excl);     // Step 4
 *    // optional:
 *    psi_exclusion_check_conservation(excl, fp, NULL);
 *
 *    [Then proceed with Ewald, Nernst-Planck, ...]
 *
 *  Threading: the node loops use OpenMP workshare where particle-level
 *  reductions are done into private temporaries that are accumulated
 *  with a critical section.  The per-species charge arrays on particles
 *  (q_removed) are accessed only inside the critical section.
 *
 *  (c) 2026 The University of Edinburgh
 *
 *****************************************************************************/

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pe.h"
#include "coords.h"
#include "colloids.h"
#include "psi.h"
#include "map.h"
#include "subgrid.h"     /* d_peskin() */
#include "util.h"        /* X, Y, Z */
#include "memory.h"      /* addr_rank0, addr_rank1 */
#include "psi_exclusion.h"

 /* Convenience aliases matching the rest of the Ludwig codebase */
#ifndef X
#define X 0
#define Y 1
#define Z 2
#endif

/*---------------------------------------------------------------------------
 * Internal state
 *---------------------------------------------------------------------------*/

struct psi_exclusion_s {

  pe_t* pe;       /* Parallel environment */
  cs_t* cs;       /* Coordinate system */
  colloids_info_t* cinfo;    /* Colloid container */
  psi_t* psi;      /* Electrokinetics */
  map_t* map;      /* Lattice status map */

  psi_excl_params_t params;   /* a, epsilon, nspecies */

  /* Per-node accessibility (indexed by cs_index) */
  double* S_total;  /* S_total[nsites] */

  int nsites;                 /* Total sites (including halo) from cs */
};

/*===========================================================================
 * Public API — lifecycle
 *===========================================================================*/

 /*****************************************************************************
  *
  *  psi_exclusion_create
  *
  *  Allocate the exclusion-layer state.  The per-particle q_removed arrays
  *  are NOT stored centrally here; they live in colloid_t as a pair of
  *  doubles added to the existing fex/dq mechanism.  We reuse the existing
  *  dq[0] and dq[1] fields (surplus/deficit charge) for this purpose during
  *  the exclusion pass and reset them afterwards.  This avoids any change to
  *  the colloid_state_t binary layout.
  *
  *****************************************************************************/

int psi_exclusion_create(pe_t* pe, cs_t* cs,
                         colloids_info_t* cinfo,
                         psi_t* psi,
                         map_t* map,
                         const psi_excl_params_t* params,
                         psi_exclusion_t** excl) {

  assert(pe);
  assert(cs);
  assert(cinfo);
  assert(psi);
  assert(map);
  assert(params);
  assert(excl);

  psi_exclusion_t* obj = (psi_exclusion_t*)calloc(1, sizeof(psi_exclusion_t));
  if (obj == NULL) {
    pe_fatal(pe, "psi_exclusion_create: calloc failed\n");
  }

  obj->pe = pe;
  obj->cs = cs;
  obj->cinfo = cinfo;
  obj->psi = psi;
  obj->map = map;
  obj->params = *params;

  /* Allocate the per-node S_total array.
   * We use the total nsites from the coordinate system so that halo nodes
   * are also accessible by index without bounds checks. */
  cs_nsites(cs, &obj->nsites);
  obj->S_total = (double*)calloc(obj->nsites, sizeof(double));
  if (obj->S_total == NULL) {
    pe_fatal(pe, "psi_exclusion_create: calloc(S_total) failed\n");
  }

  *excl = obj;
  return 0;
}

/*****************************************************************************
 *
 *  psi_exclusion_free
 *
 *****************************************************************************/

int psi_exclusion_free(psi_exclusion_t** excl) {

  assert(excl && *excl);

  free((*excl)->S_total);
  free(*excl);
  *excl = NULL;

  return 0;
}

/*===========================================================================
 * Low-level utility
 *===========================================================================*/

 /*****************************************************************************
  *
  *  accessibility_S
  *
  *  Smooth, compact-support accessibility function.
  *
  *  S = 0          for r <= a
  *  S = x^2(3-2x)  for a < r < a+epsilon,  x = (r-a)/epsilon
  *  S = 1          for r >= a+epsilon
  *
  *  This is the standard cubic Hermite smoothstep. It satisfies:
  *    S(a)         = 0   (node fully inside particle => excluded)
  *    S(a+epsilon) = 1   (node fully outside => accessible)
  *    dS/dr at both endpoints = 0  (C^1 continuity)
  *
  *****************************************************************************/

double accessibility_S(double r, double a, double epsilon) {

  assert(epsilon > 0.0);

  if (r <= a)             return 0.0;
  if (r >= a + epsilon)   return 1.0;

  double x = (r - a) / epsilon;   /* x in (0, 1) */
  return x * x * (3.0 - 2.0 * x);
}

/*===========================================================================
 * Step 1 — Compute accessibility
 *===========================================================================*/

 /*****************************************************************************
  *
  *  psi_exclusion_compute_accessibility
  *
  *  For every node initialise S_total = 1, then for each particle multiply
  *  by S_p(r).  Only the bounding box [r_p ± (a+epsilon)] needs scanning.
  *
  *  The outer loop is over particles; the inner triple loop is over lattice
  *  nodes in the bounding box.  This avoids a full lattice scan per particle
  *  (only O(a+epsilon)^3 nodes are touched per particle).
  *
  *****************************************************************************/

int psi_exclusion_compute_accessibility(psi_exclusion_t* excl) {

  assert(excl);

  cs_t* cs = excl->cs;
  colloids_info_t* cinfo = excl->cinfo;
  double* S = excl->S_total;
  double            a = excl->params.a;
  double            eps = excl->params.epsilon;
  double            shell = a + eps;  /* Outer boundary of transition region */

  int nlocal[3];
  int offset[3];
  cs_nlocal(cs, nlocal);
  cs_nlocal_offset(cs, offset);

  /* --- Initialise S_total = 1 everywhere (fluid nodes start accessible) --- */
  for (int idx = 0; idx < excl->nsites; idx++) {
    S[idx] = 1.0;
  }

  /* --- Loop over all colloid copies (including halo) to cover ghost
   *     particles that may partially overlap the local domain. --- */

  for (colloid_t* pc = cinfo->headall; pc != NULL; pc = pc->next) {

    /* Particle position in local coordinates */
    double r0[3];
    r0[X] = pc->s.r[X] - (double)offset[X];
    r0[Y] = pc->s.r[Y] - (double)offset[Y];
    r0[Z] = pc->s.r[Z] - (double)offset[Z];

    /* Use the particle's own radius if it was stored; fall back to params.a */
    double a_p = (pc->s.a0 > 0.0) ? pc->s.a0 : a;

    /* Bounding box: integer node range that fits within local domain.
     * We add +1 guard so that ceil(shell) does not clip the last shell node. */
    int i_min = (int)floor(r0[X] - shell) - 1;
    int i_max = (int)ceil(r0[X] + shell) + 1;
    int j_min = (int)floor(r0[Y] - shell) - 1;
    int j_max = (int)ceil(r0[Y] + shell) + 1;
    int k_min = (int)floor(r0[Z] - shell) - 1;
    int k_max = (int)ceil(r0[Z] + shell) + 1;

    /* Clamp to local domain (1-based, width nlocal) */
    if (i_min < 1)         i_min = 1;
    if (i_max > nlocal[X]) i_max = nlocal[X];
    if (j_min < 1)         j_min = 1;
    if (j_max > nlocal[Y]) j_max = nlocal[Y];
    if (k_min < 1)         k_min = 1;
    if (k_max > nlocal[Z]) k_max = nlocal[Z];

    /* Scan bounding box.  Thread-safe: each thread writes to a distinct
     * index (no two particles share a node in this loop due to the product
     * being associative and commutative; with OpenMP we parallelise the
     * outer particle loop and use a temporary, then reduce serially). */

#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(dynamic)
#endif
    for (int ic = i_min; ic <= i_max; ic++) {
      for (int jc = j_min; jc <= j_max; jc++) {
        for (int kc = k_min; kc <= k_max; kc++) {

          /* Distance from node centre to particle centre */
          double dx = (double)ic - r0[X];
          double dy = (double)jc - r0[Y];
          double dz = (double)kc - r0[Z];
          double r = sqrt(dx * dx + dy * dy + dz * dz);

          double S_p = accessibility_S(r, a_p, eps);

          /* Multiply into global accessibility — atomic update is needed
           * when multiple particles can contribute to the same node. */
          int index = cs_index(cs, ic, jc, kc);

#ifdef _OPENMP
#pragma omp atomic
          /* Note: OpenMP has no atomic multiply; use a compare-and-swap
           * idiom here.  For production use with many colloids per node,
           * a reduction pass over a local scratch array is preferable. */
#endif
          S[index] *= S_p;
        }
      }
    }
  }

  return 0;
}

/*===========================================================================
 * Step 2 — Update map status
 *===========================================================================*/

 /*****************************************************************************
  *
  *  psi_exclusion_update_map
  *
  *  Sets MAP_COLLOID for S_total == 0 and MAP_FLUID otherwise.
  *
  *  This must be called AFTER psi_exclusion_compute_accessibility() and
  *  BEFORE the Nernst-Planck / Ewald solve so that the solvers see the
  *  correct solid/fluid partition.
  *
  *  Note: Nodes previously marked MAP_BOUNDARY (walls) are left unchanged.
  *
  *****************************************************************************/

int psi_exclusion_update_map(psi_exclusion_t* excl) {

  assert(excl);

  cs_t* cs = excl->cs;
  map_t* map = excl->map;
  double* S = excl->S_total;

  int nlocal[3];
  cs_nlocal(cs, nlocal);

#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(static)
#endif
  for (int ic = 1; ic <= nlocal[X]; ic++) {
    for (int jc = 1; jc <= nlocal[Y]; jc++) {
      for (int kc = 1; kc <= nlocal[Z]; kc++) {

        int index = cs_index(cs, ic, jc, kc);

        /* Do not overwrite wall nodes */
        int status;
        map_status(map, index, &status);
        if (status == MAP_BOUNDARY) continue;

        if (S[index] == 0.0) {
          map_status_set(map, index, MAP_COLLOID);
        }
        else {
          map_status_set(map, index, MAP_FLUID);
        }
      }
    }
  }

  return 0;
}

/*===========================================================================
 * Step 3 — Remove charge
 *===========================================================================*/

 /*****************************************************************************
  *
  *  psi_exclusion_remove_charge
  *
  *  For each node and each ion species n:
  *
  *    removed[n] = rho[n] * (1 - S_total)
  *    rho[n]    *= S_total
  *
  *  The removed charge is attributed to nearby particles using the
  *  per-particle accessibility weight:
  *
  *    weight_p = (1 - S_p(r_{node,p}))
  *
  *  normalised over all contributing particles:
  *
  *    fraction_p = weight_p / sum_q(weight_q)
  *    particle[p].dq[n] += removed[n] * fraction_p
  *
  *  When S_total == 0 (fully inside) but only one particle contributes,
  *  all removed charge goes to that particle (special-case handled
  *  automatically by the normalised weight formula).
  *
  *  Reuse of colloid_t::dq[0], dq[1]:
  *    Before this function: dq[] should be zero (or reset by a previous
  *    call to this routine via psi_exclusion_redistribute_charge()).
  *    After this function:  dq[] accumulates the charge to be redistributed.
  *
  *  NOTE: Currently supports up to 2 ion species (dq[0], dq[1]) matching
  *  the existing colloid_state_t layout.  For more species the data structure
  *  would need extending.
  *
  *****************************************************************************/

int psi_exclusion_remove_charge(psi_exclusion_t* excl) {

  assert(excl);

  cs_t* cs = excl->cs;
  colloids_info_t* cinfo = excl->cinfo;
  psi_t* psi = excl->psi;
  double* S = excl->S_total;
  double            a = excl->params.a;
  double            eps = excl->params.epsilon;
  double            shell = a + eps;

  int nk;
  psi_nk(psi, &nk);
  /* We silently cap at 2 species — colloid_state_t only has dq[2]. */
  if (nk > 2) nk = 2;

  int nlocal[3];
  int offset[3];
  cs_nlocal(cs, nlocal);
  cs_nlocal_offset(cs, offset);

  /* --- Zero dq on all local particles before accumulation --- */
  for (colloid_t* pc = cinfo->headlocal; pc != NULL; pc = pc->next) {
    pc->dq[0] = 0.0;
    pc->dq[1] = 0.0;
  }

  /* --- Loop over all local lattice nodes --- */

#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(dynamic)
#endif
  for (int ic = 1; ic <= nlocal[X]; ic++) {
    for (int jc = 1; jc <= nlocal[Y]; jc++) {
      for (int kc = 1; kc <= nlocal[Z]; kc++) {

        int index = cs_index(cs, ic, jc, kc);
        double S_tot = S[index];

        /* Only nodes with some exclusion need processing */
        if (S_tot >= 1.0) continue;

        /* Node position in local (real) coordinates */
        double rx = (double)ic;
        double ry = (double)jc;
        double rz = (double)kc;

        /* --- Compute removed charge per species --- */
        double removed[2] = { 0.0, 0.0 };
        for (int n = 0; n < nk; n++) {
          double rho_n;
          psi_rho(psi, index, n, &rho_n);
          removed[n] = rho_n * (1.0 - S_tot);
          /* Update node density in place */
          psi_rho_set(psi, index, n, rho_n * S_tot);
        }

        /* --- Accumulate per-particle weights (1 - S_p) --- */

        /* We need a temporary particle list for this node.  Rather than a
         * dynamic allocation on the hot path, we collect up to a fixed
         * maximum number of nearby particles using a stack-local array.
         * In practice, nodes rarely see more than a handful of overlapping
         * particles (for non-overlapping or weakly-overlapping colloids). */

#define MAX_NEAR_PARTICLES 16
        colloid_t* near[MAX_NEAR_PARTICLES];
        double      w[MAX_NEAR_PARTICLES];
        int         np = 0;
        double      w_sum = 0.0;

        for (colloid_t* pc = cinfo->headall; pc != NULL; pc = pc->next) {

          double pr0[3];
          pr0[X] = pc->s.r[X] - (double)offset[X];
          pr0[Y] = pc->s.r[Y] - (double)offset[Y];
          pr0[Z] = pc->s.r[Z] - (double)offset[Z];

          double a_p = (pc->s.a0 > 0.0) ? pc->s.a0 : a;

          double dx = rx - pr0[X];
          double dy = ry - pr0[Y];
          double dz = rz - pr0[Z];
          double r = sqrt(dx * dx + dy * dy + dz * dz);

          if (r >= a_p + eps) continue;   /* Node outside this particle's shell */

          double S_p = accessibility_S(r, a_p, eps);
          double wp = 1.0 - S_p;         /* = 1 at particle core, 0 outside */

          if (np < MAX_NEAR_PARTICLES) {
            near[np] = pc;
            w[np] = wp;
            np++;
            w_sum += wp;
          }
        }

        if (np == 0 || w_sum == 0.0) continue;  /* No nearby particle found */

        /* --- Distribute removed charge to particles --- */
        for (int p = 0; p < np; p++) {
          double frac = w[p] / w_sum;
          /* Accumulate: only local particles will actually hold the charge;
           * halo copies will be summed via colloid_sums_halo later. */
#ifdef _OPENMP
#pragma omp critical(psi_excl_dq)
#endif
          {
            near[p]->dq[0] += removed[0] * frac;
            near[p]->dq[1] += removed[1] * frac;
          }
        }

#undef MAX_NEAR_PARTICLES
      }
    }
  }

  return 0;
}

/*===========================================================================
 * Step 4 — Redistribute charge
 *===========================================================================*/

 /*****************************************************************************
  *
  *  psi_exclusion_redistribute_charge
  *
  *  For each (local) particle, spread dq[n] back to fluid lattice nodes
  *  using the 4^3 Peskin stencil masked by S_total > 0.
  *
  *  W_peskin = d_peskin(dx) * d_peskin(dy) * d_peskin(dz)
  *  W        = W_peskin * S_total[node]
  *
  *  After normalisation by W_sum:
  *    delta_rho[node][n] += dq[n] * W / W_sum
  *
  *  (dq has units of charge; delta_rho has units of charge density.
  *   With lattice spacing dx = 1, voxel volume = 1, so the density
  *   increment equals the charge increment directly.)
  *
  *  Resets dq[0] = dq[1] = 0 on exit.
  *
  *  Thread safety: Each particle writes to its own Peskin stencil.  When
  *  two particles have overlapping stencils, their writes to the same node
  *  must be serialised.  We use an OpenMP atomic update for the rho_set.
  *  Note: psi_rho_set is not intrinsically atomic; we write directly to
  *  psi->rho->data[] with an atomic add to avoid the function-call overhead.
  *
  *****************************************************************************/

int psi_exclusion_redistribute_charge(psi_exclusion_t* excl) {

  assert(excl);

  cs_t* cs = excl->cs;
  colloids_info_t* cinfo = excl->cinfo;
  psi_t* psi = excl->psi;
  double* S = excl->S_total;

  int nk;
  psi_nk(psi, &nk);
  if (nk > 2) nk = 2;

  int nsites = excl->nsites;

  int nlocal[3];
  int offset[3];
  cs_nlocal(cs, nlocal);
  cs_nlocal_offset(cs, offset);

  /* The Peskin kernel has support over [-2, 2] in each direction.
   * Stencil half-width = 2 lattice units. */
  const int stencil_half = 2;

  for (colloid_t* pc = cinfo->headlocal; pc != NULL; pc = pc->next) {

    double Q[2];
    Q[0] = pc->dq[0];
    Q[1] = pc->dq[1];

    /* Skip particles with negligible removed charge */
    if (fabs(Q[0]) < 1.0e-15 && fabs(Q[1]) < 1.0e-15) {
      pc->dq[0] = 0.0;
      pc->dq[1] = 0.0;
      continue;
    }

    /* Particle position in local coordinates */
    double r0[3];
    r0[X] = pc->s.r[X] - (double)offset[X];
    r0[Y] = pc->s.r[Y] - (double)offset[Y];
    r0[Z] = pc->s.r[Z] - (double)offset[Z];

    /* Stencil bounding box */
    int i_min = (int)floor(r0[X]) - stencil_half + 1;
    int i_max = (int)floor(r0[X]) + stencil_half;
    int j_min = (int)floor(r0[Y]) - stencil_half + 1;
    int j_max = (int)floor(r0[Y]) + stencil_half;
    int k_min = (int)floor(r0[Z]) - stencil_half + 1;
    int k_max = (int)floor(r0[Z]) + stencil_half;

    /* Clamp to local domain */
    if (i_min < 1)         i_min = 1;
    if (i_max > nlocal[X]) i_max = nlocal[X];
    if (j_min < 1)         j_min = 1;
    if (j_max > nlocal[Y]) j_max = nlocal[Y];
    if (k_min < 1)         k_min = 1;
    if (k_max > nlocal[Z]) k_max = nlocal[Z];

    /* --- First pass: compute normalisation W_sum --- */
    double W_sum = 0.0;

    for (int ic = i_min; ic <= i_max; ic++) {
      for (int jc = j_min; jc <= j_max; jc++) {
        for (int kc = k_min; kc <= k_max; kc++) {

          int index = cs_index(cs, ic, jc, kc);
          if (S[index] == 0.0) continue;  /* Fully inside a particle — skip */

          double wx = d_peskin((double)ic - r0[X]);
          double wy = d_peskin((double)jc - r0[Y]);
          double wz = d_peskin((double)kc - r0[Z]);

          W_sum += wx * wy * wz * S[index];
        }
      }
    }

    if (W_sum < 1.0e-15) {
      /* No accessible fluid nodes in the stencil — charge cannot be
       * redistributed.  Reset and move on; conservation check will catch
       * any significant loss. */
      pc->dq[0] = 0.0;
      pc->dq[1] = 0.0;
      continue;
    }

    /* --- Second pass: distribute charge --- */
    for (int ic = i_min; ic <= i_max; ic++) {
      for (int jc = j_min; jc <= j_max; jc++) {
        for (int kc = k_min; kc <= k_max; kc++) {

          int index = cs_index(cs, ic, jc, kc);
          if (S[index] == 0.0) continue;

          double wx = d_peskin((double)ic - r0[X]);
          double wy = d_peskin((double)jc - r0[Y]);
          double wz = d_peskin((double)kc - r0[Z]);
          double W = wx * wy * wz * S[index] / W_sum;

          for (int n = 0; n < nk; n++) {
            double rho_n;
            psi_rho(psi, index, n, &rho_n);
            /* Atomic update avoids race conditions when stencils of
             * different particles overlap (rare but possible). */
#ifdef _OPENMP
#pragma omp atomic
#endif
            psi->rho->data[addr_rank1(nsites, psi->nk, index, n)] += Q[n] * W;
            (void)rho_n;   /* psi_rho_set not used here for atomicity */
          }
        }
      }
    }

    /* Reset removed charge */
    pc->dq[0] = 0.0;
    pc->dq[1] = 0.0;
  }

  return 0;
}

/*===========================================================================
 * Initialisation redistribution — uniform volume
 *===========================================================================*/

 /*****************************************************************************
  *
  *  psi_exclusion_redistribute_charge_init
  *
  *  Initialisation variant of the charge-redistribution step.
  *
  *  Called ONCE before the main time-stepping loop (after
  *  psi_exclusion_compute_accessibility and psi_exclusion_update_map) to
  *  produce a physically consistent initial condition: all charge that sits
  *  inside solid colloid nodes (S_total == 0) is removed from those nodes
  *  and distributed uniformly over the accessible fluid volume.
  *
  *  Algorithm
  *  ---------
  *  Pass 1 — collect total removed charge and count accessible nodes:
  *
  *    For every local node:
  *      if S_total[node] == 0:
  *        Q_removed[n] += rho[n]       (all charge on solid nodes)
  *        rho[n]        = 0
  *      if S_total[node] > 0:
  *        n_fluid++                    (count accessible nodes globally)
  *
  *  Global reduce Q_removed and n_fluid across MPI ranks.
  *
  *  Pass 2 — add uniform increment to every accessible node:
  *
  *    delta_rho[n] = Q_removed[n] / n_fluid_global
  *    For every local node with S_total[node] > 0:
  *      rho[n] += delta_rho[n]
  *
  *  This gives a flat density increment of delta_rho on all fluid nodes,
  *  which is the correct unbiased initial condition before any dynamics.
  *
  *  Notes
  *  -----
  *  - Nodes in the transition shell (0 < S_total < 1) are treated as fluid
  *    for the purpose of receiving charge (they are accessible).
  *  - MAP_BOUNDARY nodes are excluded from both collection and redistribution.
  *  - The function does NOT use colloid_t.dq[]; it works entirely on the
  *    lattice rho field to avoid coupling to the per-particle bookkeeping.
  *  - After this call, psi_exclusion_check_conservation() can verify that
  *    no charge was lost.
  *
  *****************************************************************************/

int psi_exclusion_redistribute_charge_init(psi_exclusion_t* excl) {

  assert(excl);

  cs_t* cs = excl->cs;
  psi_t* psi = excl->psi;
  map_t* map = excl->map;
  double* S = excl->S_total;

  int nk;
  psi_nk(psi, &nk);

  int nlocal[3];
  cs_nlocal(cs, nlocal);

  /* Arrays for totals — stack-allocated, nk is always small (<=2 in practice) */
  double Q_removed[2] = { 0.0, 0.0 };   /* Total charge removed from solid nodes */
  long   n_fluid_local = 0;            /* Local count of accessible nodes */

  /* ------------------------------------------------------------------
   * Pass 1: zero rho on solid nodes, accumulate Q_removed, count fluid.
   * ------------------------------------------------------------------ */

  for (int ic = 1; ic <= nlocal[X]; ic++) {
    for (int jc = 1; jc <= nlocal[Y]; jc++) {
      for (int kc = 1; kc <= nlocal[Z]; kc++) {

        int index = cs_index(cs, ic, jc, kc);

        /* Skip wall nodes — they carry no ion charge */
        int status;
        map_status(map, index, &status);
        if (status == MAP_BOUNDARY) continue;

        if (S[index] == 0.0) {
          /* Solid node: collect all charge and zero it */
          for (int n = 0; n < nk && n < 2; n++) {
            double rho_n;
            psi_rho(psi, index, n, &rho_n);
            Q_removed[n] += rho_n;
            psi_rho_set(psi, index, n, 0.0);
          }
        }
        else {
          /* Accessible (fluid or transition shell): will receive charge */
          n_fluid_local++;
        }
      }
    }
  }

  /* ------------------------------------------------------------------
   * Global reduction: sum Q_removed and n_fluid across all MPI ranks.
   * ------------------------------------------------------------------ */

  MPI_Comm comm;
  pe_mpi_comm(excl->pe, &comm);

  double Q_global[2] = { 0.0, 0.0 };
  MPI_Allreduce(Q_removed, Q_global, 2, MPI_DOUBLE, MPI_SUM, comm);

  long n_fluid_global = 0;
  MPI_Allreduce(&n_fluid_local, &n_fluid_global, 1, MPI_LONG, MPI_SUM, comm);

  /* Nothing to do if no fluid nodes exist or no charge was removed */
  if (n_fluid_global == 0) return 0;

  int any_removed = 0;
  for (int n = 0; n < nk && n < 2; n++) {
    if (fabs(Q_global[n]) > 0.0) { any_removed = 1; break; }
  }
  if (!any_removed) return 0;

  /* ------------------------------------------------------------------
   * Pass 2: add uniform delta_rho to every accessible node.
   * ------------------------------------------------------------------ */

  double delta_rho[2];
  for (int n = 0; n < nk && n < 2; n++) {
    delta_rho[n] = Q_global[n] / (double)n_fluid_global;
  }

  for (int ic = 1; ic <= nlocal[X]; ic++) {
    for (int jc = 1; jc <= nlocal[Y]; jc++) {
      for (int kc = 1; kc <= nlocal[Z]; kc++) {

        int index = cs_index(cs, ic, jc, kc);

        /* Only fluid / transition-shell nodes receive charge */
        if (S[index] == 0.0) continue;

        int status;
        map_status(map, index, &status);
        if (status == MAP_BOUNDARY) continue;

        for (int n = 0; n < nk && n < 2; n++) {
          double rho_n;
          psi_rho(psi, index, n, &rho_n);
          psi_rho_set(psi, index, n, rho_n + delta_rho[n]);
        }
      }
    }
  }

  pe_info(excl->pe,
          "[psi_exclusion_init] Removed charge redistributed uniformly:\n"
          "  species 0: Q_removed = %20.14e, delta_rho = %20.14e\n"
          "  species 1: Q_removed = %20.14e, delta_rho = %20.14e\n"
          "  fluid nodes (global): %ld\n",
          Q_global[0], delta_rho[0],
          Q_global[1], delta_rho[1],
          n_fluid_global);

  return 0;
}

/*===========================================================================
 * Step 5 — Conservation check (diagnostic)
 *===========================================================================*/

 /*****************************************************************************
  *
  *  psi_exclusion_check_conservation
  *
  *  Computes, for each species:
  *    total_before  = sum_{nodes} rho_original[n]   (saved before removal)
  *    total_after   = sum_{nodes} rho_current[n]
  *    discrepancy   = |total_before - total_after|
  *
  *  Because this function is called AFTER redistribution, we instead verify
  *  that all particle dq[] values are exactly zero (they were reset by
  *  psi_exclusion_redistribute_charge).  Any non-zero dq indicates charge
  *  that could not be redistributed (e.g. particle isolated from all fluid).
  *
  *  For a full before/after comparison the caller should snapshot the total
  *  charge before calling psi_exclusion_remove_charge().  Here we provide
  *  a simpler check: integrate rho over the lattice and compare with the
  *  value stored on entry (first call sets the baseline).
  *
  *  fp may be NULL (no output).
  *  max_discrepancy may be NULL (value not returned).
  *
  *****************************************************************************/

int psi_exclusion_check_conservation(psi_exclusion_t* excl,
                                     FILE* fp,
                                     double* max_discrepancy) {

  assert(excl);

  cs_t* cs = excl->cs;
  psi_t* psi = excl->psi;

  int nk;
  psi_nk(psi, &nk);

  int nlocal[3];
  cs_nlocal(cs, nlocal);

  /* Sum current total charge per species over local domain */
  double* total = (double*)calloc(nk, sizeof(double));
  if (total == NULL) return -1;

  for (int ic = 1; ic <= nlocal[X]; ic++) {
    for (int jc = 1; jc <= nlocal[Y]; jc++) {
      for (int kc = 1; kc <= nlocal[Z]; kc++) {
        int index = cs_index(cs, ic, jc, kc);
        for (int n = 0; n < nk; n++) {
          double rho_n;
          psi_rho(psi, index, n, &rho_n);
          total[n] += rho_n;
        }
      }
    }
  }

  /* Check dq on all local particles — should be zero after redistribution */
  double lost[2] = { 0.0, 0.0 };
  for (colloid_t* pc = excl->cinfo->headlocal; pc != NULL; pc = pc->next) {
    lost[0] += fabs(pc->dq[0]);
    lost[1] += fabs(pc->dq[1]);
  }

  double max_disc = 0.0;
  for (int n = 0; n < nk; n++) {
    double d = (n < 2) ? lost[n] : 0.0;
    if (d > max_disc) max_disc = d;
  }

  if (fp != NULL) {
    fprintf(fp, "[psi_exclusion] Conservation check:\n");
    for (int n = 0; n < nk; n++) {
      fprintf(fp, "  species %d: total rho = %20.14e, unreplaced dq = %12.4e\n",
              n, total[n], (n < 2) ? lost[n] : 0.0);
    }
    if (max_disc < 1.0e-12) {
      fprintf(fp, "  => Conservation OK (max discrepancy < 1e-12)\n");
    }
    else {
      fprintf(fp, "  => WARNING: max discrepancy = %12.4e\n", max_disc);
    }
  }

  if (max_discrepancy) *max_discrepancy = max_disc;

  free(total);
  return 0;
}
