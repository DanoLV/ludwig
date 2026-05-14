/*****************************************************************************
 *
 *  psi_petsc.c
 *
 *  A solution of the Poisson equation for the potential and
 *  charge densities stored in the psi_t object.
 *
 *  This uses the PETSc library.
 *
 *  The Poisson equation with homogeneous permittivity looks like
 *
 *    nabla^2 \psi = - rho_elec / epsilon
 *
 *  where psi is the potential, rho_elec is the free charge density, and
 *  epsilon is a permittivity.
 *
 *  There is also a version for non-uniform dielectric.
 *
 *
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  (c) 2013-2023 The University of Edinburgh
 *
 *  Contributing Authors:
 *  Oliver Henrich  (now U. Strathclyde)
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *
 *****************************************************************************/

#include <assert.h>
#include <float.h>
#include <math.h>

#include "pe.h"
#include "coords.h"
#include "util.h"
#include "psi_petsc.h"

static psi_solver_vt_t vt_ = {
  (psi_solver_free_ft)psi_solver_petsc_free,
  (psi_solver_solve_ft)psi_solver_petsc_solve
};

static psi_solver_vt_t vart_ = {
  (psi_solver_free_ft)psi_solver_petsc_free,
  (psi_solver_solve_ft)psi_solver_petsc_var_epsilon_solve
};

int psi_solver_petsc_initialise(psi_t* psi, psi_solver_petsc_t* solver);
int psi_solver_petsc_matrix_set(psi_solver_petsc_t* solver);
int psi_solver_petsc_rhs_set(psi_solver_petsc_t* solver);
int psi_solver_petsc_psi_to_da(psi_solver_petsc_t* solver);
int psi_solver_petsc_da_to_psi(psi_solver_petsc_t* solver);

int psi_solver_petsc_var_epsilon_initialise(psi_t* psi, var_epsilon_t epsilon,
              psi_solver_petsc_t* solver);
int psi_solver_petsc_var_epsilon_matrix_set(psi_solver_petsc_t* solver);
int psi_solver_petsc_var_epsilon_rhs_set(psi_solver_petsc_t* solver);


/*****************************************************************************
 *
 *  psi_solver_petsc_create
 *
 *****************************************************************************/

int psi_solver_petsc_create(psi_t* psi, psi_solver_petsc_t** solver) {

  int ifail = -1;                     /* Check PETSC is available */
  int isInitialised = 0;

  PetscInitialised(&isInitialised);

  if (isInitialised) {
    psi_solver_petsc_t* petsc = NULL;

    petsc = (psi_solver_petsc_t*)calloc(1, sizeof(psi_solver_petsc_t));
    assert(petsc);

    if (petsc != NULL) {
      /* initialise ... */
      petsc->super.impl = &vt_;
      petsc->psi = psi;
      ifail = psi_solver_petsc_initialise(psi, petsc);
      if (ifail != 0) free(petsc);
      if (ifail == 0) *solver = petsc;
    }
  }

  return ifail;
}

/*****************************************************************************
 *
 *  psi_solver_petsc_var_epsilon_create
 *
 *****************************************************************************/

int psi_solver_petsc_var_epsilon_create(psi_t* psi, var_epsilon_t user,
          psi_solver_petsc_t** solver) {

  int ifail = -1;                     /* Check PETSC is available */
  int isInitialised = 0;

  PetscInitialised(&isInitialised);

  if (isInitialised) {
    psi_solver_petsc_t* petsc = NULL;

    petsc = (psi_solver_petsc_t*)calloc(1, sizeof(psi_solver_petsc_t));
    assert(petsc);

    if (petsc != NULL) {
      /* initialise ... */
      petsc->super.impl = &vart_;
      petsc->psi = psi;
      petsc->fe = user.fe;
      petsc->epsilon = user.epsilon;
      ifail = psi_solver_petsc_initialise(psi, petsc);
      if (ifail != 0) free(petsc);
      if (ifail == 0) *solver = petsc;
    }
  }

  return ifail;
}

/*****************************************************************************
 *
 *  psi_solver_petsc_free
 *
 *****************************************************************************/

int psi_solver_petsc_free(psi_solver_petsc_t** solver) {

  assert(solver && *solver);

  free(*solver);
  *solver = NULL;

  return 0;
}


#ifndef PETSC

/*****************************************************************************
 *
 *  psi_solver_petsc_initialise
 *
 *  There are two stub routines here for the case that PETSc is not
 *  avialable.
 *
 *****************************************************************************/

int psi_solver_petsc_initialise(psi_t* psi, psi_solver_petsc_t* solver) {

  /* No implementation */
  return -1;
}

/*****************************************************************************
 *
 *  psi_solver_petsc_solve
 *
 *****************************************************************************/

int psi_solver_petsc_solve(psi_solver_petsc_t* solver, int ntimestep) {

  /* No implementation */
  return -1;
}

/*****************************************************************************
 *
 *  psi_solver_petsc_var_epsilon_solve
 *
 *****************************************************************************/

int psi_solver_petsc_var_epsilon_solve(psi_solver_petsc_t* solver, int nt) {

  /* No implementation */
  return -1;
}

#else

#include "petscdmda.h"
#include "petscksp.h"

/* Here's the internal state. */

struct psi_solver_petsc_block_s {
  DM  da;         /* Domain management */
  Mat a;          /* System matrix */
  Vec x;          /* Unknown (potential) */
  Vec b;          /* Right-hand side */
  KSP ksp;        /* Krylov solver context */
};

/*****************************************************************************
 *
 *  psi_solver_petsc_initialise
 *
 *****************************************************************************/

int psi_solver_petsc_initialise(psi_t* psi, psi_solver_petsc_t* solver) {

  assert(psi);
  assert(solver);

  {
    size_t sz = sizeof(psi_solver_petsc_block_t);
    solver->block = (psi_solver_petsc_block_t*)calloc(1, sz);
    assert(solver->block);
    if (solver->block == NULL) return -1;
  }

  /* In order for the DMDA and the Cartesian MPI communicator
   * to share the same part of the domain decomposition it is
   *  necessary to renumber the process ranks of the default
   *  PETSc communicator. Default PETSc is column major decomposition. */

  {
    cs_t* cs = psi->cs;
    int coords[3] = { 0 };
    int cartsz[3] = { 0 };
    int ntotal[3] = { 0 };
    int nhalo = -1;
    int rank = -1;

    MPI_Comm comm = MPI_COMM_NULL;
    DMBoundaryType periodic = DM_BOUNDARY_PERIODIC;

    cs_cartsz(cs, cartsz);
    cs_cart_coords(cs, coords);
    cs_nhalo(cs, &nhalo);
    cs_ntotal(cs, ntotal);

    /* Set new rank according to PETSc ordering */
    /* Create communicator with new ranks according to PETSc ordering */
    /* Override default PETSc communicator */

    rank = coords[Z] * cartsz[Y] * cartsz[X] + coords[Y] * cartsz[X] + coords[X];
    MPI_Comm_split(PETSC_COMM_WORLD, 1, rank, &comm);
    PETSC_COMM_WORLD = comm;

    /* Create 3D distributed array (always periodic) */

    DMDACreate3d(PETSC_COMM_WORLD, periodic, periodic, periodic,
     DMDA_STENCIL_BOX, ntotal[X], ntotal[Y], ntotal[Z],
     cartsz[X], cartsz[Y], cartsz[Z], 1, nhalo,
     NULL, NULL, NULL, &solver->block->da);

    /*CHANGE INIT - 20260326 Use GPU vec/mat types when -vec_type cuda is set via .petscrc.
     * DMSetVecType/DMSetMatType must match the requested type BEFORE DMSetUp and
     * DMCreateGlobalVector, otherwise -vec_type cuda in .petscrc is ignored and
     * hypre PC fails with "HYPRE_MEMORY_DEVICE expects a device vector".
     * Original (CPU only):
     *   PetscCall(DMSetVecType(solver->block->da, VECSTANDARD));
     *   PetscCall(DMSetMatType(solver->block->da, MATMPIAIJ)); */
    {
      PetscBool use_cuda = PETSC_FALSE;
      char vtype[64];
      PetscOptionsGetString(NULL, NULL, "-vec_type", vtype, sizeof(vtype), &use_cuda);
      if (use_cuda && !strcmp(vtype, "cuda")) {
        PetscCall(DMSetVecType(solver->block->da, VECCUDA));
        PetscCall(DMSetMatType(solver->block->da, MATAIJCUSPARSE));
      } else {
        PetscCall(DMSetVecType(solver->block->da, VECSTANDARD));
        PetscCall(DMSetMatType(solver->block->da, MATMPIAIJ));
      }
    }
    /*CHANGE END - 20260326 */
    PetscCall(DMSetUp(solver->block->da));
  }

  /* Create global vectors and matrix */

  DMCreateMatrix(solver->block->da, &solver->block->a);
  DMCreateGlobalVector(solver->block->da, &solver->block->x);
  VecDuplicate(solver->block->x, &solver->block->b);

  /* Initialise solver context */

  {
    PetscReal abstol = psi->solver.abstol;
    PetscReal rtol = psi->solver.reltol;
    PetscInt  maxits = psi->solver.maxits;

    KSPCreate(PETSC_COMM_WORLD, &solver->block->ksp);
    KSPSetOperators(solver->block->ksp, solver->block->a, solver->block->a);
    KSPSetTolerances(solver->block->ksp, rtol, abstol, PETSC_DEFAULT, maxits);
    /*CHANGE INIT - 20260326 Allow .petscrc / command-line options for KSP solver */
    KSPSetFromOptions(solver->block->ksp);
    /*CHANGE END - 20260326 */
  }

  /* Not required in var-epsilon case, but no harm. */
  psi_solver_petsc_matrix_set(solver);

  KSPSetUp(solver->block->ksp);

  return 0;
}

/*****************************************************************************
 *
 *  psi_solver_petsc_matrix_set
 *
 *****************************************************************************/

int psi_solver_petsc_matrix_set(psi_solver_petsc_t* solver) {

  int xs, ys, zs;
  int xw, yw, zw;
  int xe, ye, ze;
  double epsilon;

  double v[27] = { 0 };        /* Accomodate largest current stencil */
  MatStencil col[27] = { 0 };  /* Ditto */

  stencil_t* s = solver->psi->stencil;

  assert(solver);
  assert(solver->psi->solver.nstencil <= 27);

  /* Obtain start and width ... */
  DMDAGetCorners(solver->block->da, &xs, &ys, &zs, &xw, &yw, &zw);

  xe = xs + xw;
  ye = ys + yw;
  ze = zs + zw;

  /* 3D-Laplacian with periodic BCs */
  /* Uniform dielectric constant */

  psi_epsilon(solver->psi, &epsilon);

  for (int k = zs; k < ze; k++) {
    for (int j = ys; j < ye; j++) {
      for (int i = xs; i < xe; i++) {

        /*CHANGE INIT - 20251119 CUDA C++ compatibility fix */
        /* Original: MatStencil row = {.i = i, .j = j, .k = k}; */
        /* CUDA nvcc requires explicit member assignment instead of designated initializers */
        MatStencil row;
        row.i = i;
        row.j = j;
        row.k = k;
        /*CHANGE END*/

        for (int p = 0; p < s->npoints; p++) {
          col[p].i = i + s->cv[p][X];
          col[p].j = j + s->cv[p][Y];
          col[p].k = k + s->cv[p][Z];
          v[p] = s->wlaplacian[p] * epsilon;
        }
        MatSetValuesStencil(solver->block->a, 1, &row, s->npoints, col, v,
                INSERT_VALUES);
      }
    }
  }

  /* Matrix assembly & halo swap */
  /* Retain the non-zero structure of the matrix */

  MatAssemblyBegin(solver->block->a, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(solver->block->a, MAT_FINAL_ASSEMBLY);
  MatSetOption(solver->block->a, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE);

  /* Set the matrix, and the nullspace */
  KSPSetOperators(solver->block->ksp, solver->block->a, solver->block->a);

  {
    MatNullSpace nullsp;
    MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, NULL, &nullsp);
    MatSetNullSpace(solver->block->a, nullsp);
    MatNullSpaceDestroy(&nullsp);
  }

  return 0;
}

/*****************************************************************************
 *
 *  psi_solver_petsc_solve
 *
 *****************************************************************************/

int psi_solver_petsc_solve(psi_solver_petsc_t* solver, int ntimestep) {

  /*CHANGE INIT - 20250112 Fix NULL pointer handling for test compatibility */
  /* Original: assert(solver); */
  /* The assert causes core dump in tests that check NULL handling.
   * Changed to return error code instead. */
  if (solver == NULL) return -1;
  /*CHANGE END - 20250112 */

  psi_solver_petsc_rhs_set(solver);
  psi_solver_petsc_psi_to_da(solver);

  /*CHANGE INIT - 20260326 Fix DIVERGED_INDEFINITE_MAT with aijcusparse:
   * MatSetNullSpace does not project the null space from the RHS/solution when
   * using GPU matrices. Explicit projection via MatNullSpaceRemove on b and x
   * ensures CG sees a consistent SPD system (no constant-mode contamination). */
  {
    MatNullSpace nullsp;
    MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, NULL, &nullsp);
    MatNullSpaceRemove(nullsp, solver->block->b);
    MatNullSpaceRemove(nullsp, solver->block->x);
    MatNullSpaceDestroy(&nullsp);
  }
  /*CHANGE END - 20260326 */

  KSPSetInitialGuessNonzero(solver->block->ksp, PETSC_TRUE);
  KSPSolve(solver->block->ksp, solver->block->b, solver->block->x);

  if (ntimestep % solver->psi->solver.nfreq == 0) {
    /* Report on progress of the solver.
     * Note the default Petsc residual is the preconditioned L2 norm. */
    pe_t* pe = solver->psi->pe;
    PetscInt  its = 0;
    PetscReal norm = 0.0;
    PetscCall(KSPGetIterationNumber(solver->block->ksp, &its));
    PetscCall(KSPGetResidualNorm(solver->block->ksp, &norm));
    pe_info(pe, "\n");
    pe_info(pe, "Krylov solver\n");
    pe_info(pe, "Norm of residual %g at %d iterations\n", norm, its);
  }

  psi_solver_petsc_da_to_psi(solver);

  return 0;
}

/*****************************************************************************
 *
 *  psi_solver_petsc_rhs_set
 *
 *****************************************************************************/

int psi_solver_petsc_rhs_set(psi_solver_petsc_t* solver) {

  cs_t* cs = NULL;
  int xs, ys, zs;
  int xw, yw, zw;
  int xe, ye, ze;
  int offset[3] = { 0 };
  double e0[3] = { 0 };
  double*** rho_3d = { 0 };

  assert(solver);

  cs = solver->psi->cs;
  cs_nlocal_offset(cs, offset);

  DMDAGetCorners(solver->block->da, &xs, &ys, &zs, &xw, &yw, &zw);
  DMDAVecGetArray(solver->block->da, solver->block->b, &rho_3d);

  xe = xs + xw;
  ye = ys + yw;
  ze = zs + zw;

  for (int k = zs; k < ze; k++) {
    int kc = k - offset[Z] + 1;
    for (int j = ys; j < ye; j++) {
      int jc = j - offset[Y] + 1;
      for (int i = xs; i < xe; i++) {

        int ic = i - offset[X] + 1;
        int index = cs_index(cs, ic, jc, kc);
        double rho_elec = 0.0;
        /* Non-dimensional potential in Poisson eqn requires e/kT */
        double eunit = solver->psi->e;
        double beta = solver->psi->beta;

        psi_rho_elec(solver->psi, index, &rho_elec);
        rho_3d[k][j][i] = rho_elec * eunit * beta;
      }
    }
  }

  /* Modify right hand side for external electric field */
  /* The system must be periodic, so no need to check. */

  e0[X] = solver->psi->e0[X];
  e0[Y] = solver->psi->e0[Y];
  e0[Z] = solver->psi->e0[Z];

  if (e0[X] || e0[Y] || e0[Z]) {

    int ntotal[3] = { 0 };
    int mpi_coords[3] = { 0 };
    int mpi_cartsz[3] = { 0 };
    double epsilon = 0.0;

    cs_ntotal(cs, ntotal);
    cs_cart_coords(cs, mpi_coords);
    cs_cartsz(cs, mpi_cartsz);

    psi_epsilon(solver->psi, &epsilon);

    if (e0[X] && mpi_coords[X] == 0) {
      for (int k = zs; k < ze; k++) {
        for (int j = ys; j < ye; j++) {
          rho_3d[k][j][0] += epsilon * e0[X] * ntotal[X];
        }
      }
    }

    if (e0[X] && mpi_coords[X] == mpi_cartsz[X] - 1) {
      for (int k = zs; k < ze; k++) {
        for (int j = ys; j < ye; j++) {
          rho_3d[k][j][xe - 1] -= epsilon * e0[X] * ntotal[X];
        }
      }
    }

    if (e0[Y] && mpi_coords[Y] == 0) {
      for (int k = zs; k < ze; k++) {
        for (int i = xs; i < xe; i++) {
          rho_3d[k][0][i] += epsilon * e0[Y] * ntotal[Y];
        }
      }
    }

    if (e0[Y] && mpi_coords[Y] == mpi_cartsz[Y] - 1) {
      for (int k = zs; k < ze; k++) {
        for (int i = xs; i < xe; i++) {
          rho_3d[k][ye - 1][i] -= epsilon * e0[Y] * ntotal[Y];
        }
      }
    }

    if (e0[Z] && mpi_coords[Z] == 0) {
      for (int j = ys; j < ye; j++) {
        for (int i = xs; i < xe; i++) {
          rho_3d[0][j][i] += epsilon * e0[Z] * ntotal[Z];
        }
      }
    }

    if (e0[Z] && mpi_coords[Z] == mpi_cartsz[Z] - 1) {
      for (int j = ys; j < ye; j++) {
        for (int i = xs; i < xe; i++) {
          rho_3d[ze - 1][j][i] -= epsilon * e0[Z] * ntotal[Z];
        }
      }
    }

  }

  DMDAVecRestoreArray(solver->block->da, solver->block->b, &rho_3d);

  return 0;
}

// /*CHANGE INIT - Add subgrid particle charges to PETSc RHS */
// /*****************************************************************************
//  *
//  *  psi_solver_petsc_add_subgrid_charges
//  *
//  *  Add subgrid particle charges to the PETSc RHS vector using a regularized
//  *  delta function (Gaussian or Peskin).
//  *
//  *  This modifies the RHS to include point charges at continuous positions:
//  *    rhs[i,j,k] += Q_particle * delta_regularized(r - r_particle)
//  *
//  *  Call this AFTER psi_solver_petsc_rhs_set() and BEFORE solving.
//  *
//  *****************************************************************************/

// #include "subgrid.h"  /* For d_peskin */
// #include <math.h>

// int psi_solver_petsc_add_subgrid_charges(psi_solver_petsc_t * solver,
//                                           colloids_info_t * cinfo) {

//   cs_t * cs = NULL;
//   int xs, ys, zs;
//   int xw, yw, zw;
//   int xe, ye, ze;
//   int offset[3] = {0};
//   int nlocal[3] = {0};
//   double *** rho_3d = NULL;
//   int ncell[3];
//   int ic, jc, kc;
//   colloid_t* pc = NULL;

//   assert(solver);
//   assert(cinfo);

//   if (cinfo->nsubgrid == 0) return 0;

//   cs = solver->psi->cs;
//   cs_nlocal(cs, nlocal);
//   cs_nlocal_offset(cs, offset);
//   colloids_info_ncell(cinfo, ncell);

//   DMDAGetCorners(solver->block->da, &xs, &ys, &zs, &xw, &yw, &zw);
//   DMDAVecGetArray(solver->block->da, solver->block->b, &rho_3d);

//   xe = xs + xw;
//   ye = ys + yw;
//   ze = zs + zw;

//   double eunit = solver->psi->e;
//   double beta  = solver->psi->beta;
//   int valency[2];
//   int nk;

//   psi_nk(solver->psi, &nk);
//   assert(nk == 2);  /* Currently only support 2 species */
//   psi_valency(solver->psi, 0, &valency[0]);
//   psi_valency(solver->psi, 1, &valency[1]);

//   /* Loop over all particles */
//   for (ic = 0; ic <= ncell[X] + 1; ic++) {
//     for (jc = 0; jc <= ncell[Y] + 1; jc++) {
//       for (kc = 0; kc <= ncell[Z] + 1; kc++) {

//         colloids_info_cell_list_head(cinfo, ic, jc, kc, &pc);

//         for (; pc; pc = pc->next) {

//           if (pc->s.bc != COLLOID_BC_SUBGRID) continue;

//           /* Net charge on particle accounting for valencies */
//           /* This matches psi_rho_elec: sum over species of valency[n] * q[n] */
//           double q_net = valency[0] * pc->s.q0 + valency[1] * pc->s.q1;

//           /*DEBUG - Print particle charge info */
//           static int debug_count = 0;
//           if (debug_count < 1) {
//             printf("DEBUG: Particle charge: q0=%.15e, q1=%.15e, q_net=%.15e\n",
//                    pc->s.q0, pc->s.q1, q_net);
//             printf("DEBUG: Valencies: v[0]=%d, v[1]=%d, eunit=%.15e, beta=%.15e\n",
//                    valency[0], valency[1], eunit, beta);
//             debug_count++;
//           }

//           /* Particle position in Ludwig coordinates (base-1 global) */
//           double r_ludwig[3];
//           r_ludwig[X] = pc->s.r[X];
//           r_ludwig[Y] = pc->s.r[Y];
//           r_ludwig[Z] = pc->s.r[Z];

//           /* Convert to PETSc coordinates (base-0 global) */
//           /* Ludwig: r=1.0 is first node, PETSc: i=0 is first node */
//           double r_petsc[3];
//           r_petsc[X] = r_ludwig[X] - 1.0;
//           r_petsc[Y] = r_ludwig[Y] - 1.0;
//           r_petsc[Z] = r_ludwig[Z] - 1.0;

//           /* Determine range in PETSc GLOBAL indices (base-0) */
//           int ig_min = imax(xs, (int)floor(r_petsc[X] - 2.0));
//           int ig_max = imin(xe - 1, (int)ceil(r_petsc[X] + 2.0));
//           int jg_min = imax(ys, (int)floor(r_petsc[Y] - 2.0));
//           int jg_max = imin(ye - 1, (int)ceil(r_petsc[Y] + 2.0));
//           int kg_min = imax(zs, (int)floor(r_petsc[Z] - 2.0));
//           int kg_max = imin(ze - 1, (int)ceil(r_petsc[Z] + 2.0));

//           /* Add charge contribution to nearby nodes */
//           for (int kg = kg_min; kg <= kg_max; kg++) {
//             for (int jg = jg_min; jg <= jg_max; jg++) {
//               for (int ig = ig_min; ig <= ig_max; ig++) {

//                 /* Distance from particle to node */
//                 /* PETSc node (ig,jg,kg) is at position (ig,jg,kg) in PETSc coords */
//                 double r[3];
//                 r[X] = r_petsc[X] - 1.0 * ig;
//                 r[Y] = r_petsc[Y] - 1.0 * jg;
//                 r[Z] = r_petsc[Z] - 1.0 * kg;

//                 /* Regularized delta function (Peskin) */
//                 extern double d_peskin(double);
//                 double delta = d_peskin(r[X]) * d_peskin(r[Y]) * d_peskin(r[Z]);

//                 /* Add charge contribution to RHS */
//                 /* Factor eunit*beta for non-dimensional potential */
//                 rho_3d[kg][jg][ig] += q_net * delta * eunit * beta;
//               }
//             }
//           }
//         }
//       }
//     }
//   }

//   DMDAVecRestoreArray(solver->block->da, solver->block->b, &rho_3d);

//   return 0;
// }

// /*****************************************************************************
//  *
//  *  psi_solver_petsc_solve_with_subgrid
//  *
//  *  Solve the Poisson-Boltzmann equation including subgrid particle charges
//  *  as regularized delta functions in the RHS.
//  *
//  *  This is a wrapper that:
//  *    1. Sets up RHS with ionic charge densities
//  *    2. Adds subgrid particle charges as regularized sources
//  *    3. Solves the Poisson equation
//  *    4. Copies solution back to psi
//  *
//  *  Usage:
//  *    psi_solver_petsc_solve_with_subgrid(solver, cinfo, timestep);
//  *
//  *****************************************************************************/

// int psi_solver_petsc_solve_with_subgrid(psi_solver_petsc_t * solver,
//                                          colloids_info_t * cinfo,
//                                          int ntimestep) {

//   assert(solver);
//   assert(cinfo);

//   /* 1. Build RHS with ionic charge densities (without particle distribution) */
//   psi_solver_petsc_rhs_set(solver);

//   /*DEBUG - Print RHS before adding subgrid charges */
//   if (ntimestep % 100 == 0) {
//     PetscScalar ***rho_3d_debug;
//     DMDAVecGetArray(solver->block->da, solver->block->b, &rho_3d_debug);
//     int xs, ys, zs, xw, yw, zw;
//     DMDAGetCorners(solver->block->da, &xs, &ys, &zs, &xw, &yw, &zw);
//     if (xs == 0 && ys == 0 && zs == 0) {
//       printf("DEBUG step %d: RHS[0,0,0] before subgrid = %.15e\n", ntimestep, rho_3d_debug[0][0][0]);
//     }
//     DMDAVecRestoreArray(solver->block->da, solver->block->b, &rho_3d_debug);
//   }

//   /* 2. Add subgrid particle charges as regularized delta functions */
//   psi_solver_petsc_add_subgrid_charges(solver, cinfo);

//   /*DEBUG - Print RHS after adding subgrid charges */
//   if (ntimestep % 100 == 0) {
//     PetscScalar ***rho_3d_debug;
//     DMDAVecGetArray(solver->block->da, solver->block->b, &rho_3d_debug);
//     int xs, ys, zs, xw, yw, zw;
//     DMDAGetCorners(solver->block->da, &xs, &ys, &zs, &xw, &yw, &zw);
//     if (xs == 0 && ys == 0 && zs == 0) {
//       printf("DEBUG step %d: RHS[0,0,0] after subgrid = %.15e\n", ntimestep, rho_3d_debug[0][0][0]);
//     }
//     DMDAVecRestoreArray(solver->block->da, solver->block->b, &rho_3d_debug);
//   }

//   /* 3. Set initial guess from current potential */
//   psi_solver_petsc_psi_to_da(solver);

//   /* 4. Solve the Poisson equation */
//   KSPSetInitialGuessNonzero(solver->block->ksp, PETSC_TRUE);
//   KSPSolve(solver->block->ksp, solver->block->b, solver->block->x);

//   /* 5. Report solver progress if needed */
//   if (ntimestep % solver->psi->solver.nfreq == 0) {
//     PetscReal rnorm = 0.0;
//     PetscInt nits = 0;
//     KSPGetIterationNumber(solver->block->ksp, &nits);
//     KSPGetResidualNorm(solver->block->ksp, &rnorm);
//     pe_info(solver->psi->pe, "PETSc (subgrid): %4d iterations; residual: %14.7e\n",
//             nits, rnorm);
//   }

//   /* 6. Copy solution back to psi */
//   psi_solver_petsc_da_to_psi(solver);

//   return 0;
// }
// /*CHANGE END - Add subgrid particle charges to PETSc RHS */

/*****************************************************************************
 *
 *  psi_solver_petsc_psi_to_da
 *
 *  Copy the potential from the psi_t represetation to the solution
 *  vector as an initial guess.
 *
 *****************************************************************************/

int psi_solver_petsc_psi_to_da(psi_solver_petsc_t* solver) {

  cs_t* cs = NULL;
  int xs, ys, zs;
  int xw, yw, zw;
  int xe, ye, ze;
  int offset[3] = { 0 };
  double*** psi_3d = NULL;

  assert(solver);

  cs = solver->psi->cs;
  cs_nlocal_offset(cs, offset);

  DMDAGetCorners(solver->block->da, &xs, &ys, &zs, &xw, &yw, &zw);
  DMDAVecGetArray(solver->block->da, solver->block->x, &psi_3d);

  xe = xs + xw;
  ye = ys + yw;
  ze = zs + zw;

  for (int k = zs; k < ze; k++) {
    int kc = k - offset[Z] + 1;
    for (int j = ys; j < ye; j++) {
      int jc = j - offset[Y] + 1;
      for (int i = xs; i < xe; i++) {
        int ic = i - offset[X] + 1;
        int index = cs_index(cs, ic, jc, kc);
        psi_3d[k][j][i] = solver->psi->psi->data[index];
      }
    }
  }

  DMDAVecRestoreArray(solver->block->da, solver->block->x, &psi_3d);

  return 0;
}

/*****************************************************************************
 *
 *  psi_solver_petsc_da_to_psi
 *
 *  Copy te Petsc solution back to the psi_t represetation.
 *
 *****************************************************************************/

int psi_solver_petsc_da_to_psi(psi_solver_petsc_t* solver) {

  cs_t* cs = NULL;
  int xs, ys, zs;
  int xw, yw, zw;
  int xe, ye, ze;
  int offset[3] = { 0 };
  double*** psi_3d = NULL;

  assert(solver);

  cs = solver->psi->cs;
  cs_nlocal_offset(cs, offset);

  DMDAGetCorners(solver->block->da, &xs, &ys, &zs, &xw, &yw, &zw);

  /*CHANGE INIT - 20260326 If the solution vector is on GPU (VECCUDA), bind it to CPU
   * so that DMDAVecGetArray can access it from the host. */
  VecBindToCPU(solver->block->x, PETSC_TRUE);
  /*CHANGE END - 20260326 */

  DMDAVecGetArray(solver->block->da, solver->block->x, &psi_3d);

  xe = xs + xw;
  ye = ys + yw;
  ze = zs + zw;

  for (int k = zs; k < ze; k++) {
    int kc = k - offset[Z] + 1;
    for (int j = ys; j < ye; j++) {
      int jc = j - offset[Y] + 1;
      for (int i = xs; i < xe; i++) {
        int ic = i - offset[X] + 1;
        int index = cs_index(cs, ic, jc, kc);
        solver->psi->psi->data[index] = psi_3d[k][j][i];
      }
    }
  }

  DMDAVecRestoreArray(solver->block->da, solver->block->x, &psi_3d);

  /*CHANGE INIT - 20260326 Restore GPU binding after host read. */
  VecBindToCPU(solver->block->x, PETSC_FALSE);
  /*CHANGE END - 20260326 */

  return 0;
}

/*****************************************************************************
 *
 *  psi_solver_petsc_var_epsilon_solve
 *
 *****************************************************************************/

int psi_solver_petsc_var_epsilon_solve(psi_solver_petsc_t* solver, int nt) {

  assert(solver);

  psi_solver_petsc_var_epsilon_matrix_set(solver);
  psi_solver_petsc_var_epsilon_rhs_set(solver);

  psi_solver_petsc_psi_to_da(solver);

  KSPSetInitialGuessNonzero(solver->block->ksp, PETSC_TRUE);
  KSPSolve(solver->block->ksp, solver->block->b, solver->block->x);

  if (nt % solver->psi->solver.nfreq == 0) {
    /* Report on progress of the solver.
     * Note the default Petsc residual is the preconditioned L2 norm. */
    pe_t* pe = solver->psi->pe;
    PetscInt  its = 0;
    PetscReal norm = 0.0;
    PetscCall(KSPGetIterationNumber(solver->block->ksp, &its));
    PetscCall(KSPGetResidualNorm(solver->block->ksp, &norm));
    pe_info(pe, "\n");
    pe_info(pe, "Krylov solver (with dielectric contrast)\n");
    pe_info(pe, "Norm of residual %g at %d iterations\n", norm, its);
  }

  psi_solver_petsc_da_to_psi(solver);

  return 0;
}

/*****************************************************************************
 *
 *  psi_solver_petsc_var_epsilon_matrix_set
 *
 *****************************************************************************/

int psi_solver_petsc_var_epsilon_matrix_set(psi_solver_petsc_t* solver) {

  cs_t* cs = NULL;
  int xs, ys, zs;
  int xw, yw, zw;
  int xe, ye, ze;
  int offset[3] = { 0 };

  double v[27] = { 0 };
  MatStencil col[27] = { 0 };
  stencil_t* s = solver->psi->stencil;

  assert(solver);

  cs = solver->psi->cs;
  cs_nlocal_offset(cs, offset);

  /* Get details of the distributed array data structure.
     The PETSc directives return global indices, but
     every process works only on its local copy. */

  DMDAGetCorners(solver->block->da, &xs, &ys, &zs, &xw, &yw, &zw);

  xe = xs + xw;
  ye = ys + yw;
  ze = zs + zw;

  /* 3D-operator with periodic BCs */

  for (int k = zs; k < ze; k++) {
    int kc = 1 + k - offset[Z];
    for (int j = ys; j < ye; j++) {
      int jc = 1 + j - offset[Y];
      for (int i = xs; i < xe; i++) {

        int ic = 1 + i - offset[X];
        int index = cs_index(cs, ic, jc, kc);
        double epsilon0 = 0.0;
        double gradeps[3] = { 0 };

        /*CHANGE INIT - 20251119 CUDA C++ compatibility fix */
        /* Original: MatStencil row = {.i = i, .j = j, .k = k}; */
        /* CUDA nvcc requires explicit member assignment instead of designated initializers */
        MatStencil row;
        row.i = i;
        row.j = j;
        row.k = k;
        /*CHANGE END*/

        solver->epsilon(solver->fe, index, &epsilon0);

        /* Local approx. to grad epsilon ... */
        for (int p = 1; p < s->npoints; p++) {
          int ic1 = ic + s->cv[p][X];
          int jc1 = jc + s->cv[p][Y];
          int kc1 = kc + s->cv[p][Z];
          int index1 = cs_index(cs, ic1, jc1, kc1);
          double epsilon1 = 0.0;
          solver->epsilon(solver->fe, index1, &epsilon1);
          gradeps[X] += s->wgradients[p] * s->cv[p][X] * epsilon1;
          gradeps[Y] += s->wgradients[p] * s->cv[p][Y] * epsilon1;
          gradeps[Z] += s->wgradients[p] * s->cv[p][Z] * epsilon1;
        }

        for (int p = 0; p < s->npoints; p++) {
          col[p].i = i + s->cv[p][X];
          col[p].j = j + s->cv[p][Y];
          col[p].k = k + s->cv[p][Z];

          /* Laplacian part of operator */
          v[p] = s->wlaplacian[p] * epsilon0;

          /* Addtional terms in generalised Poisson equation */
          v[p] += s->wgradients[p] * s->cv[p][X] * gradeps[X];
          v[p] += s->wgradients[p] * s->cv[p][Y] * gradeps[Y];
          v[p] += s->wgradients[p] * s->cv[p][Z] * gradeps[Z];
        }

        MatSetValuesStencil(solver->block->a, 1, &row, s->npoints, col, v,
                INSERT_VALUES);
      }
    }
  }

  /* Matrix assembly & halo swap */
  /* Retain the non-zero structure of the matrix */

  MatAssemblyBegin(solver->block->a, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(solver->block->a, MAT_FINAL_ASSEMBLY);
  MatSetOption(solver->block->a, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE);

  /* Set the matrix, preconditioner and nullspace */
  KSPSetOperators(solver->block->ksp, solver->block->a, solver->block->a);

  {
    MatNullSpace nullsp;
    MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, NULL, &nullsp);
    MatSetNullSpace(solver->block->a, nullsp);
    MatNullSpaceDestroy(&nullsp);
  }

  KSPSetFromOptions(solver->block->ksp);

  return 0;
}

/*****************************************************************************
 *
 *  psi_solver_petsc_var_epsilon_rhs_set
 *
 *  The only difference here (cf. uniform epsilon) is in the external
 *  field terms.
 *
 *****************************************************************************/

int psi_solver_petsc_var_epsilon_rhs_set(psi_solver_petsc_t* solver) {

  cs_t* cs = NULL;
  int xs, ys, zs;
  int xw, yw, zw;
  int xe, ye, ze;
  int offset[3] = { 0 };
  double e0[3] = { 0 };
  double*** rho_3d = { 0 };

  assert(solver);

  cs = solver->psi->cs;
  cs_nlocal_offset(cs, offset);

  DMDAGetCorners(solver->block->da, &xs, &ys, &zs, &xw, &yw, &zw);
  DMDAVecGetArray(solver->block->da, solver->block->b, &rho_3d);

  xe = xs + xw;
  ye = ys + yw;
  ze = zs + zw;

  for (int k = zs; k < ze; k++) {
    int kc = k - offset[Z] + 1;
    for (int j = ys; j < ye; j++) {
      int jc = j - offset[Y] + 1;
      for (int i = xs; i < xe; i++) {

        int ic = i - offset[X] + 1;
        int index = cs_index(cs, ic, jc, kc);
        double rho_elec = 0.0;
        /* Non-dimensional potential in Poisson eqn requires e/kT */
        double eunit = solver->psi->e;
        double beta = solver->psi->beta;

        psi_rho_elec(solver->psi, index, &rho_elec);
        rho_3d[k][j][i] = rho_elec * eunit * beta;
      }
    }
  }

  /* Modify right hand side for external electric field */
  /* The system must be periodic, so no need to check. */

  e0[X] = solver->psi->e0[X];
  e0[Y] = solver->psi->e0[Y];
  e0[Z] = solver->psi->e0[Z];

  if (e0[X] || e0[Y] || e0[Z]) {

    int ntotal[3] = { 0 };
    int mpi_coords[3] = { 0 };
    int mpi_cartsz[3] = { 0 };

    cs_ntotal(cs, ntotal);
    cs_cart_coords(cs, mpi_coords);
    cs_cartsz(cs, mpi_cartsz);

    if (e0[X] && mpi_coords[X] == 0) {
      for (int k = zs; k < ze; k++) {
        int kc = 1 + k - offset[Z];
        for (int j = ys; j < ye; j++) {
          int jc = 1 + j - offset[Y];
          int index = cs_index(cs, 1, jc, kc);
          double epsilon = 0.0;
          solver->epsilon(solver->fe, index, &epsilon);
          rho_3d[k][j][0] += epsilon * e0[X] * ntotal[X];
        }
      }
    }

    if (e0[X] && mpi_coords[X] == mpi_cartsz[X] - 1) {
      for (int k = zs; k < ze; k++) {
        int kc = 1 + k - offset[Z];
        for (int j = ys; j < ye; j++) {
          int jc = 1 + j - offset[Y];
          int ic = xe - offset[X];
          int index = cs_index(cs, ic, jc, kc);
          double epsilon = 0.0;
          solver->epsilon(solver->fe, index, &epsilon);
          rho_3d[k][j][xe - 1] -= epsilon * e0[X] * ntotal[X];
        }
      }
    }

    if (e0[Y] && mpi_coords[Y] == 0) {
      for (int k = zs; k < ze; k++) {
        int kc = 1 + k - offset[Z];
        for (int i = xs; i < xe; i++) {
          int ic = 1 + i - offset[X];
          int index = cs_index(cs, ic, 1, kc);
          double epsilon = 0.0;
          solver->epsilon(solver->fe, index, &epsilon);
          rho_3d[k][0][i] += epsilon * e0[Y] * ntotal[Y];
        }
      }
    }

    if (e0[Y] && mpi_coords[Y] == mpi_cartsz[Y] - 1) {
      for (int k = zs; k < ze; k++) {
        int kc = 1 + k - offset[Z];
        for (int i = xs; i < xe; i++) {
          int jc = ye - offset[Y];
          int ic = 1 + i - offset[X];
          int index = cs_index(cs, ic, jc, kc);
          double epsilon = 0.0;
          solver->epsilon(solver->fe, index, &epsilon);
          rho_3d[k][ye - 1][i] -= epsilon * e0[Y] * ntotal[Y];
        }
      }
    }

    if (e0[Z] && mpi_coords[Z] == 0) {
      for (int j = ys; j < ye; j++) {
        int jc = 1 + j - offset[Y];
        for (int i = xs; i < xe; i++) {
          int ic = 1 + i - offset[X];
          int index = cs_index(cs, ic, jc, 1);
          double epsilon = 0.0;
          solver->epsilon(solver->fe, index, &epsilon);
          rho_3d[0][j][i] += epsilon * e0[Z] * ntotal[Z];
        }
      }
    }

    if (e0[Z] && mpi_coords[Z] == mpi_cartsz[Z] - 1) {
      int kc = ze - offset[Z];
      for (int j = ys; j < ye; j++) {
        int jc = 1 + j - offset[Y];
        for (int i = xs; i < xe; i++) {
          int ic = 1 + i - offset[X];
          int index = cs_index(cs, ic, jc, kc);
          double epsilon = 0.0;
          solver->epsilon(solver->fe, index, &epsilon);
          rho_3d[ze - 1][j][i] -= epsilon * e0[Z] * ntotal[Z];
        }
      }
    }
  }

  DMDAVecRestoreArray(solver->block->da, solver->block->b, &rho_3d);

  return 0;
}

/*****************************************************************************
 *
 *  psi_solver_petsc_finalise
 *
 *****************************************************************************/

int psi_solver_petsc_finalise(psi_solver_petsc_t* solver) {

  assert(solver);

  PetscCall(KSPDestroy(&solver->block->ksp));
  PetscCall(VecDestroy(&solver->block->x));
  PetscCall(VecDestroy(&solver->block->b));
  PetscCall(MatDestroy(&solver->block->a));
  PetscCall(DMDestroy(&solver->block->da));

  free(solver->block);

  return 0;
}

#endif
