/*****************************************************************************
 *
 *  psi_exclusion.h
 *
 *  Ion exclusion layer for subgrid colloidal particles.
 *
 *  Implements a smooth accessibility function S(r) that enforces ion
 *  exclusion inside colloidal particles and redistributes removed charge
 *  back to the surrounding fluid nodes via the Peskin interpolation kernel.
 *
 *  Overview of the algorithm (call order in the simulation loop):
 *
 *  1. psi_exclusion_compute_accessibility()
 *       Computes S_total[node] = product of per-particle accessibility
 *       values S_p(r) = smoothstep on the shell [a, a+epsilon].
 *
 *  2. psi_exclusion_update_map()
 *       Sets MAP_COLLOID / MAP_FLUID for each node based on S_total.
 *
 *  3. psi_exclusion_remove_charge()
 *       removed = rho * (1 - S_total), rho *= S_total.
 *       Assigns removed charge to nearby particles weighted by (1 - S_p).
 *
 *  4. psi_exclusion_redistribute_charge()
 *       Redistributes particle.q_removed back to surrounding fluid nodes
 *       using the Peskin kernel masked by S_total > 0.
 *
 *  5. psi_exclusion_check_conservation() — optional diagnostic.
 *
 *  (c) 2026 The University of Edinburgh
 *
 *****************************************************************************/

#ifndef LUDWIG_PSI_EXCLUSION_H
#define LUDWIG_PSI_EXCLUSION_H

#include "pe.h"
#include "coords.h"
#include "colloids.h"
#include "psi.h"
#include "map.h"

/*---------------------------------------------------------------------------
 * Parameters
 *---------------------------------------------------------------------------*/

/**
 * psi_excl_params_t
 *
 * Parameters controlling the smooth exclusion layer.
 *
 *  a         : colloid hard-sphere radius (lattice units).  Read from the
 *              first particle (s.a0); override here for a monodisperse run.
 *  epsilon   : transition-layer thickness (lattice units).
 *              Recommended value: 0.5 * dx  (with dx = 1 in lattice units).
 *  nspecies  : number of ion species (must equal psi->nk).
 */
typedef struct psi_excl_params_s {
  double a;           /* Colloid radius (lattice units) */
  double epsilon;     /* Transition-layer thickness (lattice units) */
  int    nspecies;    /* Ion species count */
} psi_excl_params_t;

/*---------------------------------------------------------------------------
 * Opaque state handle
 *---------------------------------------------------------------------------*/

typedef struct psi_exclusion_s psi_exclusion_t;

/*---------------------------------------------------------------------------
 * Lifecycle
 *---------------------------------------------------------------------------*/

/**
 * psi_exclusion_create()
 *
 * Allocate the exclusion-layer state, including the per-node S_total array
 * and the per-particle q_removed arrays.
 *
 * Arguments:
 *   pe      - parallel environment
 *   cs      - coordinate system (provides nsites, nlocal, offset)
 *   cinfo   - colloid container (provides particle positions and count)
 *   psi     - electrokinetics object (provides nk, rho access)
 *   map     - lattice status map (will be modified in place)
 *   params  - parameters (a, epsilon, nspecies)
 *   excl    - OUTPUT: newly allocated object
 *
 * Returns 0 on success, non-zero on failure.
 */
int psi_exclusion_create(pe_t * pe, cs_t * cs,
                         colloids_info_t * cinfo,
                         psi_t * psi,
                         map_t * map,
                         const psi_excl_params_t * params,
                         psi_exclusion_t ** excl);

/**
 * psi_exclusion_free()
 *
 * Release all memory associated with the exclusion-layer object.
 */
int psi_exclusion_free(psi_exclusion_t ** excl);

/*---------------------------------------------------------------------------
 * Core algorithm steps
 *---------------------------------------------------------------------------*/

/**
 * psi_exclusion_compute_accessibility()
 *
 * Step 1.  Compute S_total[node] for every lattice node.
 *
 * S_total is initialised to 1 and then multiplied by S_p(r) for every
 * colloid whose exclusion shell overlaps the node:
 *
 *   S_total[node] *= S_p(r_{node,p})
 *
 * where S_p is the cubic smoothstep on [a, a+epsilon].
 *
 * Only nodes within distance (a + epsilon) of any colloid centre are
 * affected; all others retain S_total = 1.
 */
int psi_exclusion_compute_accessibility(psi_exclusion_t * excl);

/**
 * psi_exclusion_update_map()
 *
 * Step 2.  Update map status from S_total:
 *
 *   S_total == 0  =>  MAP_COLLOID
 *   S_total  > 0  =>  MAP_FLUID
 */
int psi_exclusion_update_map(psi_exclusion_t * excl);

/**
 * psi_exclusion_remove_charge()
 *
 * Step 3.  Remove charge from nodes that overlap with colloids.
 *
 * For each species n and each node:
 *   removed[n] = rho[n] * (1 - S_total)
 *   rho[n]    *= S_total
 *
 * The removed charge is accumulated to colloid.q_removed[n] weighted by
 * the per-particle accessibility (1 - S_p), normalised across all
 * contributing colloids.
 */
int psi_exclusion_remove_charge(psi_exclusion_t * excl);

/**
 * psi_exclusion_redistribute_charge()
 *
 * Step 4.  Spread q_removed[n] back to the fluid using the Peskin kernel
 * masked by S_total > 0.
 *
 * For each colloid p:
 *   Loop over the 4^3 Peskin stencil nodes around the colloid centre.
 *   W = d_peskin(dx) * d_peskin(dy) * d_peskin(dz) * S_total[node]
 *   After normalisation: delta_rho[node][n] += Q_removed[n] * W / W_sum
 *
 * Resets colloid.q_removed to zero on exit.
 */
int psi_exclusion_redistribute_charge(psi_exclusion_t * excl);

/**
 * psi_exclusion_redistribute_charge_init()
 *
 * Initialisation-only redistribution.  Call ONCE before the main loop,
 * after psi_exclusion_compute_accessibility() and psi_exclusion_update_map().
 *
 * Removes all charge from nodes where S_total == 0 (fully inside a colloid)
 * and redistributes it uniformly over all accessible fluid nodes
 * (S_total > 0, not MAP_BOUNDARY):
 *
 *   delta_rho[n] = Q_removed_global[n] / n_fluid_global
 *
 * This produces a flat, unbiased initial condition consistent with the
 * particle geometry, without any Peskin-stencil bias near particle surfaces.
 *
 * Does NOT use colloid_t.dq[].  Works directly on the lattice rho field.
 * Performs one MPI_Allreduce to synchronise Q_removed and n_fluid.
 */
int psi_exclusion_redistribute_charge_init(psi_exclusion_t * excl);

/**
 * psi_exclusion_check_conservation()
 *
 * Step 5 (diagnostic).  Verify that total charge is conserved.
 *
 * Prints to fp (may be NULL to suppress output).
 * Returns the maximum absolute discrepancy across all species.
 */
int psi_exclusion_check_conservation(psi_exclusion_t * excl,
                                     FILE * fp,
                                     double * max_discrepancy);

/*---------------------------------------------------------------------------
 * Low-level utility — exposed for unit tests
 *---------------------------------------------------------------------------*/

/**
 * accessibility_S()
 *
 * Smooth accessibility function with compact support:
 *
 *   r <= a         => 0
 *   a < r < a+eps  => x^2*(3 - 2x),  x = (r-a)/eps
 *   r >= a+eps     => 1
 *
 * S(a) = 0, S(a+eps) = 1, S'(a) = S'(a+eps) = 0 (C^1 smooth).
 */
double accessibility_S(double r, double a, double epsilon);

#endif /* LUDWIG_PSI_EXCLUSION_H */
