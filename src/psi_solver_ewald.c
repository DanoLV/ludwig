/*****************************************************************************
 *
 *  psi_solver_ewald.c
 *
 *  Thin wrapper that exposes ewald_charge_sum_full_gpu as a psi_solver_t
 *  so it can be assigned to ludwig->poisson and called via impl->solve.
 *
 *  CHANGE INIT - 20260512 new file
 *
 *****************************************************************************/

#include <assert.h>
#include <stdlib.h>

#include "psi_solver_ewald.h"

/* Forward declarations for vtable */
static int ewald_solver_free_vt (psi_solver_t ** psolver);
static int ewald_solver_solve_vt(psi_solver_t *  solver, int ntimestep);

static const psi_solver_vt_t vt_ = {
  .free  = ewald_solver_free_vt,
  .solve = ewald_solver_solve_vt,
};

/*****************************************************************************
 *
 *  psi_solver_ewald_create
 *
 *****************************************************************************/

int psi_solver_ewald_create(ewald_charge_t * ewald,
                            FILE * fp,
                            ewald_self_table_t * self_table,
                            psi_solver_ewald_t ** psolver) {

  assert(ewald);
  assert(psolver && *psolver == NULL);

  psi_solver_ewald_t * obj = (psi_solver_ewald_t *) calloc(1, sizeof(*obj));
  if (obj == NULL) return -1;

  obj->super.impl = &vt_;
  obj->ewald      = ewald;
  obj->fp         = fp;
  obj->self_table = self_table;

  *psolver = obj;
  return 0;
}

/*****************************************************************************
 *
 *  psi_solver_ewald_free
 *
 *****************************************************************************/

int psi_solver_ewald_free(psi_solver_ewald_t ** psolver) {

  assert(psolver && *psolver);

  free(*psolver);
  *psolver = NULL;
  return 0;
}

/*****************************************************************************
 *
 *  psi_solver_ewald_solve
 *
 *****************************************************************************/

int psi_solver_ewald_solve(psi_solver_ewald_t * solver, int ntimestep) {

  assert(solver);
  return ewald_charge_sum_full_gpu(solver->ewald, solver->fp, solver->self_table);
}

/* CHANGE INIT - 20260513 */

int psi_solver_ewald_set_fp(psi_solver_ewald_t * solver, FILE * fp) {

  assert(solver);

  solver->fp = fp;
  return 0;
}

/* CHANGE END - 20260513 */

/* vtable trampolines */

static int ewald_solver_free_vt(psi_solver_t ** psolver) {
  return psi_solver_ewald_free((psi_solver_ewald_t **) psolver);
}

static int ewald_solver_solve_vt(psi_solver_t * solver, int ntimestep) {
  return psi_solver_ewald_solve((psi_solver_ewald_t *) solver, ntimestep);
}

/* CHANGE INIT - Gaussian_Ewald */

static int ewald_gaussian_solver_free_vt (psi_solver_t ** psolver);
static int ewald_gaussian_solver_solve_vt(psi_solver_t *  solver, int ntimestep);

static const psi_solver_vt_t vt_gaussian_ = {
  .free  = ewald_gaussian_solver_free_vt,
  .solve = ewald_gaussian_solver_solve_vt,
};

int psi_solver_ewald_gaussian_create(ewald_charge_t * ewald,
                                     FILE * fp,
                                     double sigma,
                                     psi_solver_ewald_gaussian_t ** psolver) {

  assert(ewald);
  assert(psolver && *psolver == NULL);

  psi_solver_ewald_gaussian_t * obj =
    (psi_solver_ewald_gaussian_t *) calloc(1, sizeof(*obj));
  if (obj == NULL) return -1;

  obj->super.impl = &vt_gaussian_;
  obj->ewald       = ewald;
  obj->fp          = fp;
  obj->sigma       = sigma;

  *psolver = obj;
  return 0;
}

int psi_solver_ewald_gaussian_free(psi_solver_ewald_gaussian_t ** psolver) {

  assert(psolver && *psolver);

  free(*psolver);
  *psolver = NULL;
  return 0;
}

int psi_solver_ewald_gaussian_solve(psi_solver_ewald_gaussian_t * solver,
                                    int ntimestep) {

  assert(solver);
  return ewald_charge_sum_full_gaussian_gpu(solver->ewald, solver->fp,
                                            solver->sigma);
}

/* CHANGE INIT - 20260513 */

int psi_solver_ewald_gaussian_set_fp(psi_solver_ewald_gaussian_t * solver,
                                     FILE * fp) {
  assert(solver);

  solver->fp = fp;
  return 0;
}

/* CHANGE END - 20260513 */

static int ewald_gaussian_solver_free_vt(psi_solver_t ** psolver) {
  return psi_solver_ewald_gaussian_free((psi_solver_ewald_gaussian_t **) psolver);
}

static int ewald_gaussian_solver_solve_vt(psi_solver_t * solver, int ntimestep) {
  return psi_solver_ewald_gaussian_solve((psi_solver_ewald_gaussian_t *) solver,
                                         ntimestep);
}

/* CHANGE END - Gaussian_Ewald */

/* CHANGE INIT - Gaussian_Ewald_Dual */

static int ewald_gaussian_dual_solver_free_vt (psi_solver_t ** psolver);
static int ewald_gaussian_dual_solver_solve_vt(psi_solver_t *  solver, int ntimestep);

static const psi_solver_vt_t vt_gaussian_dual_ = {
  .free  = ewald_gaussian_dual_solver_free_vt,
  .solve = ewald_gaussian_dual_solver_solve_vt,
};

int psi_solver_ewald_gaussian_dual_create(ewald_charge_t * ewald,
                                          FILE * fp,
                                          double sigma_p,
                                          double sigma_f,
                                          psi_solver_ewald_gaussian_dual_t ** psolver) {

  assert(ewald);
  assert(psolver && *psolver == NULL);

  psi_solver_ewald_gaussian_dual_t * obj =
    (psi_solver_ewald_gaussian_dual_t *) calloc(1, sizeof(*obj));
  if (obj == NULL) return -1;

  obj->super.impl = &vt_gaussian_dual_;
  obj->ewald      = ewald;
  obj->fp         = fp;
  obj->sigma_p    = sigma_p;
  obj->sigma_f    = sigma_f;

  *psolver = obj;
  return 0;
}

int psi_solver_ewald_gaussian_dual_free(psi_solver_ewald_gaussian_dual_t ** psolver) {

  assert(psolver && *psolver);

  free(*psolver);
  *psolver = NULL;
  return 0;
}

int psi_solver_ewald_gaussian_dual_solve(psi_solver_ewald_gaussian_dual_t * solver,
                                         int ntimestep) {

  assert(solver);
  return ewald_charge_sum_full_gaussian_dual_gpu(solver->ewald, solver->fp,
                                                 solver->sigma_p, solver->sigma_f);
}

int psi_solver_ewald_gaussian_dual_set_fp(psi_solver_ewald_gaussian_dual_t * solver,
                                          FILE * fp) {
  assert(solver);

  solver->fp = fp;
  return 0;
}

static int ewald_gaussian_dual_solver_free_vt(psi_solver_t ** psolver) {
  return psi_solver_ewald_gaussian_dual_free(
      (psi_solver_ewald_gaussian_dual_t **) psolver);
}

static int ewald_gaussian_dual_solver_solve_vt(psi_solver_t * solver, int ntimestep) {
  return psi_solver_ewald_gaussian_dual_solve(
      (psi_solver_ewald_gaussian_dual_t *) solver, ntimestep);
}

/* CHANGE END - Gaussian_Ewald_Dual */
