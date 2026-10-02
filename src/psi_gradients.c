/****************************************************************************
 *
 *  psi_gradients.c
 *
 *  Currently just routines for the electric field, aka, the gradient
 *  of the potential.
 *
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  Oliver Henrich (ohenrich@epcc.ed.ac.uk)
 *
 *  (c) 2014-2023 The University of Edinburgh
 *
 ****************************************************************************/

#include <assert.h>

#include "coords.h"
#include "psi_gradients.h"

/*****************************************************************************
 *
 *  psi_electric_field
 *
 *  Return the electric field associated with the current potential.
 *  E_a = - \nabla_a \psi
 *
 *****************************************************************************/

/*CHANGE INIT - 20260926 constant-potential walls (psi_petsc.c).
 * With walls (psi->wall_map set), a fluid node uses the DIFFERENCE form
 *   E_a = - sum_p wgradients[p] c_pa [psi(r+c_p) - psi(r)] (1 + chi(r+c_p))
 * with chi = 1 at MAP_BOUNDARY sites: the same factor 2 on wall links as
 * the wall Laplacian (Asta et al. 2019, eq. 15), because the physical wall
 * sits at HALF a spacing from the last fluid node. It is the average of the
 * two face fluxes that Laplacian uses, i.e. consistent with the discrete
 * Gauss law. The difference form is required: the plain sum below only
 * equals it because sum_p wgradients[p] c_pa = 0, which no longer holds once
 * a link is doubled. The plain sum next to a wall gives exactly 3/4 of the
 * field for a linear profile and does not converge with resolution.
 * E = 0 at wall sites (inside the conductor); the plain sum there mixed both
 * electrodes across the periodic boundary.
 * Without walls the original code below runs unchanged.
 * Still open: whether Nernst-Planck needs anything near the wall (its link
 * fluxes skip non-fluid neighbours and the static Boltzmann profile is
 * exact; the charging dynamics has not been checked). */
/*CHANGE END - 20260926 */
int psi_electric_field(psi_t * psi, int index, double e[3]) {

  int ijk[3] = {0};
  cs_t * cs = NULL;
  stencil_t * s = NULL;

  assert(psi);

  cs = psi->cs;
  s  = psi->stencil;
  assert(cs);
  assert(s);

  cs_index_to_ijk(cs, index, ijk);

  e[X] = 0;
  e[Y] = 0;
  e[Z] = 0;

  /*CHANGE INIT - 20260926 constant-potential walls (see comment above) */
  if (psi->wall_map) {
    int status = MAP_FLUID;
    map_status(psi->wall_map, index, &status);
    if (status == MAP_BOUNDARY) return 0;
    {
      double psi00 = psi->psi->data[addr_rank0(psi->nsites, index)];
      for (int p = 1; p < s->npoints; p++) {
        int8_t cx = s->cv[p][X];
        int8_t cy = s->cv[p][Y];
        int8_t cz = s->cv[p][Z];
        int index1 = cs_index(cs, ijk[X] + cx, ijk[Y] + cy, ijk[Z] + cz);
        double psi1 = psi->psi->data[addr_rank0(psi->nsites, index1)];
        int status1 = MAP_FLUID;
        double fac = 1.0;
        map_status(psi->wall_map, index1, &status1);
        if (status1 == MAP_BOUNDARY) fac = 2.0;
        e[X] -= s->wgradients[p]*cx*(psi1 - psi00)*fac;
        e[Y] -= s->wgradients[p]*cy*(psi1 - psi00)*fac;
        e[Z] -= s->wgradients[p]*cz*(psi1 - psi00)*fac;
      }
    }
    return 0;
  }
  /*CHANGE END - 20260926 */

  for (int p = 1; p < s->npoints; p++) {

    int8_t cx = s->cv[p][X];
    int8_t cy = s->cv[p][Y];
    int8_t cz = s->cv[p][Z];

    int index1 = cs_index(cs, ijk[X] + cx, ijk[Y] + cy, ijk[Z] + cz);
    double psi0 = psi->psi->data[addr_rank0(psi->nsites, index1)];

    /* E_a = -\nabla_a psi */
    e[X] -= s->wgradients[p]*cx*psi0;
    e[Y] -= s->wgradients[p]*cy*psi0;
    e[Z] -= s->wgradients[p]*cz*psi0;
  }

  return 0;
}
