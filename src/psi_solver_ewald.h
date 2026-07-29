/*****************************************************************************
 *
 *  psi_solver_ewald.h
 *
 *  Wrapper that exposes ewald_charge_sum_full_gpu as a psi_solver_t.
 *
 *  CHANGE INIT - 20260512 new file
 *
 *****************************************************************************/

#ifndef LUDWIG_PSI_SOLVER_EWALD_H
#define LUDWIG_PSI_SOLVER_EWALD_H

#include <stdio.h>
#include "psi_solver.h"
#include "ewald_charge.h"

typedef struct psi_solver_ewald_s psi_solver_ewald_t;

struct psi_solver_ewald_s {
  psi_solver_t        super;       /* superclass block — must be first */
  ewald_charge_t    * ewald;       /* Ewald solver object */
  FILE              * fp;          /* output file (may be NULL) */
  ewald_self_table_t* self_table;  /* self-field table (may be NULL) */
};

int psi_solver_ewald_create(ewald_charge_t * ewald,
                            FILE * fp,
                            ewald_self_table_t * self_table,
                            psi_solver_ewald_t ** psolver);

int psi_solver_ewald_free(psi_solver_ewald_t ** psolver);

int psi_solver_ewald_solve(psi_solver_ewald_t * solver, int ntimestep);

/* CHANGE INIT - 20260513 */
int psi_solver_ewald_set_fp(psi_solver_ewald_t * solver, FILE * fp);
/* CHANGE END - 20260513 */

/* CHANGE INIT - Gaussian_Ewald */
typedef struct psi_solver_ewald_gaussian_s psi_solver_ewald_gaussian_t;

struct psi_solver_ewald_gaussian_s {
  psi_solver_t    super;   /* superclass block — must be first */
  ewald_charge_t* ewald;   /* Ewald solver object */
  FILE          * fp;      /* output file (may be NULL) */
  double          sigma;   /* Gaussian width σ (lattice units) */
};

int psi_solver_ewald_gaussian_create(ewald_charge_t * ewald,
                                     FILE * fp,
                                     double sigma,
                                     psi_solver_ewald_gaussian_t ** psolver);

int psi_solver_ewald_gaussian_free(psi_solver_ewald_gaussian_t ** psolver);

int psi_solver_ewald_gaussian_solve(psi_solver_ewald_gaussian_t * solver,
                                    int ntimestep);

/* CHANGE INIT - 20260513 */
int psi_solver_ewald_gaussian_set_fp(psi_solver_ewald_gaussian_t * solver,
                                     FILE * fp);
/* CHANGE END - 20260513 */

/* CHANGE END - Gaussian_Ewald */

/* CHANGE INIT - Gaussian_Ewald_Dual */
typedef struct psi_solver_ewald_gaussian_dual_s psi_solver_ewald_gaussian_dual_t;

struct psi_solver_ewald_gaussian_dual_s {
  psi_solver_t    super;     /* superclass block — must be first */
  ewald_charge_t* ewald;     /* Ewald solver object */
  FILE          * fp;        /* output file (may be NULL) */
  double          sigma_p;   /* Particle Gaussian width σ_p */
  double          sigma_f;   /* Fluid node Gaussian width σ_f */
};

int psi_solver_ewald_gaussian_dual_create(ewald_charge_t * ewald,
                                          FILE * fp,
                                          double sigma_p,
                                          double sigma_f,
                                          psi_solver_ewald_gaussian_dual_t ** psolver);

int psi_solver_ewald_gaussian_dual_free(psi_solver_ewald_gaussian_dual_t ** psolver);

int psi_solver_ewald_gaussian_dual_solve(psi_solver_ewald_gaussian_dual_t * solver,
                                         int ntimestep);

int psi_solver_ewald_gaussian_dual_set_fp(psi_solver_ewald_gaussian_dual_t * solver,
                                          FILE * fp);
/* CHANGE END - Gaussian_Ewald_Dual */

#endif
