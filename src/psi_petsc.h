/*****************************************************************************
 *
 *  psi_petsc.h
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  (c) 2012-2023 The University of Edinburgh
 *
 *  Contributing Authors:
 *    Oliver Henrich (ohenrich@epcc.ed.ac.uk)
 *    Kevin Stratford (kevin@epcc.ed.ac.uk)
 *
 *****************************************************************************/

#ifndef LUDWIG_PSI_SOLVER_PETSC_H
#define LUDWIG_PSI_SOLVER_PETSC_H

#include "psi_solver.h"

typedef struct psi_solver_petsc_s psi_solver_petsc_t;
typedef struct psi_solver_petsc_block_s psi_solver_petsc_block_t;

struct psi_solver_petsc_s {
  psi_solver_t super;                /* Superclass block */
  psi_t * psi;                       /* Retain a reference to psi_t */
  fe_t * fe;                         /* Free energy */
  var_epsilon_ft epsilon;            /* Variable dielectric model */
  psi_solver_petsc_block_t * block;  /* Opaque internal information. */
};

int psi_solver_petsc_create(psi_t * psi, psi_solver_petsc_t ** solver);
int psi_solver_petsc_free(psi_solver_petsc_t ** solver);
int psi_solver_petsc_solve(psi_solver_petsc_t * solver, int ntimestep);

int psi_solver_petsc_var_epsilon_create(psi_t * psi, var_epsilon_t epsilon,
					psi_solver_petsc_t ** solver);
int psi_solver_petsc_var_epsilon_solve(psi_solver_petsc_t * solver, int nt);

/*CHANGE INIT - Add subgrid particle charges to PETSc RHS */
#include "colloids.h"
int psi_solver_petsc_add_subgrid_charges(psi_solver_petsc_t * solver,
                                          colloids_info_t * cinfo);
int psi_solver_petsc_solve_with_subgrid(psi_solver_petsc_t * solver,
                                         colloids_info_t * cinfo,
                                         int ntimestep);
/*CHANGE END - Add subgrid particle charges to PETSc RHS */

/*CHANGE INIT - 20260926 constant-potential (Dirichlet) walls, following
 * Asta, Palaia, Trizac, Levesque & Rotenberg, arXiv:1907.04732, sec. II-C.
 * MAP_BOUNDARY sites are held at a prescribed potential; the Laplacian at a
 * fluid site doubles every link that reaches a wall site (their eq. 15), which
 * puts the wall at the mid-plane between the last fluid and first solid node
 * (second order). Without a call to psi_solver_petsc_wall_set() the solver is
 * the original fully periodic one. */
#include "map.h"
int psi_solver_petsc_wall_set(psi_solver_petsc_t * solver, map_t * map,
                              int axis, const double psi_wall[2]);
int psi_solver_petsc_wall_charge(const psi_solver_petsc_t * solver,
                                 double q[2], double * q_fluid);
/*CHANGE END - 20260926 */

#endif

