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
/*CHANGE INIT - 20260926 constant-potential walls */
int psi_solver_petsc_wall_charge_compute(psi_solver_petsc_t* solver);
/*CHANGE END - 20260926 */


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

/*CHANGE INIT - 20260926 constant-potential walls: stubs without PETSc */
int psi_solver_petsc_wall_set(psi_solver_petsc_t* solver, map_t* map,
                              int axis, const double psi_wall[2]) {
  return -1;
}

int psi_solver_petsc_wall_charge(const psi_solver_petsc_t* solver,
                                 double q[2], double* q_fluid) {
  return -1;
}
/*CHANGE END - 20260926 */

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
  /*CHANGE INIT - 20260926 constant-potential walls (all zero = periodic) */
  map_t* map;            /* non-NULL once psi_solver_petsc_wall_set() ran */
  int wall_axis;         /* X, Y or Z: which half a wall site belongs to */
  double wall_psi[2];    /* prescribed potential, lower / upper half */
  double wall_diag;      /* diagonal of the decoupled Dirichlet rows */
  double wall_q[2];      /* induced charge, lower / upper wall (last solve) */
  double fluid_q;        /* sum of rho_elec over fluid sites (last solve) */
  /*CHANGE END - 20260926 */
};

/*CHANGE INIT - 20260926 constant-potential walls: helpers.
 * Arguments are LOCAL lattice coordinates; halo values are allowed (the map
 * halo is refreshed in psi_solver_petsc_wall_set()). */

static int petsc_wall_site(const psi_solver_petsc_t* solver,
                           int ic, int jc, int kc) {
  int status = MAP_FLUID;
  int index = cs_index(solver->psi->cs, ic, jc, kc);
  map_status(solver->block->map, index, &status);
  return (status == MAP_BOUNDARY);
}

/* A wall site belongs to the lower (0) or upper (1) electrode according to
 * which half of the box it lies in along wall_axis. Halo coordinates are
 * wrapped back into 1..ntotal first. */

static int petsc_wall_side(const psi_solver_petsc_t* solver,
                           int ic, int jc, int kc) {
  int ntotal[3] = { 0 };
  int noffset[3] = { 0 };
  int lc[3] = { ic, jc, kc };
  int a = solver->block->wall_axis;
  cs_ntotal(solver->psi->cs, ntotal);
  cs_nlocal_offset(solver->psi->cs, noffset);
  {
    int n = ntotal[a];
    int g = noffset[a] + lc[a];
    g = ((g - 1) % n + n) % n + 1;
    return (2 * g <= n) ? 0 : 1;
  }
}

static double petsc_wall_value(const psi_solver_petsc_t* solver,
                               int ic, int jc, int kc) {
  return solver->block->wall_psi[petsc_wall_side(solver, ic, jc, kc)];
}
/*CHANGE END - 20260926 */

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

    /*CHANGE INIT - 20260921 DMDA stencil width = reach of the Poisson stencil,
     * not the lattice halo.
     * The width only sizes the matrix preallocation, (2w+1)^3 entries per row
     * for DMDA_STENCIL_BOX. The Poisson matrix is assembled from cv[p] with
     * components in {-1,0,1} for D3Q7/19/27, so w=1 is exact. Passing nhalo
     * (3 or 4 with Hann/Peskin6 kernels) preallocated 343-729 entries per row
     * for a 27-point operator, and the PETSc 3.25.4 / hypre 3.1.0 stack
     * rebuilt on 2026-08-04 aborts with "CUDA ERROR (code = 2, out of memory)
     * at memory.c:304" (hypre's own memory.c) for w >= 3 -- reproduced in a
     * standalone DMDA test with no Ludwig involved. Only global vectors are
     * used here (no DMGlobalToLocal), so no code relies on wider ghosts.
     * Original:
     *   DMDACreate3d(..., cartsz[X], cartsz[Y], cartsz[Z], 1, nhalo, ...); */
    (void) nhalo;
    DMDACreate3d(PETSC_COMM_WORLD, periodic, periodic, periodic,
     DMDA_STENCIL_BOX, ntotal[X], ntotal[Y], ntotal[Z],
     cartsz[X], cartsz[Y], cartsz[Z], 1, 1,
     NULL, NULL, NULL, &solver->block->da);
    /*CHANGE END - 20260921 */

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

  /*CHANGE INIT - 20260926 constant-potential walls.
   * Rebuilding (wall_set runs after initialise) must start from zero values:
   * the nonzero pattern is frozen (MAT_NEW_NONZERO_LOCATIONS false) and a row
   * that now writes fewer entries would otherwise keep the old ones. */
  int offset[3] = { 0 };
  cs_nlocal_offset(solver->psi->cs, offset);
  if (solver->block->map) MatZeroEntries(solver->block->a);
  /*CHANGE END - 20260926 */

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

        /*CHANGE INIT - 20260926 constant-potential walls (Asta et al. 2019,
         * eq. 15). Wall site: decoupled row D x = D psi_wall (rhs_set).
         * Fluid site: each link to a wall site is doubled, and its column is
         * ELIMINATED (moved to the right-hand side in rhs_set) so the matrix
         * stays symmetric -- CG + boomeramg need that. The diagonal can no
         * longer be read from wlaplacian[0]: it is minus the sum of the
         * scaled links, rebuilt per site. */
        if (solver->block->map) {
          int ic = i - offset[X] + 1;
          int jc = j - offset[Y] + 1;
          int kc = k - offset[Z] + 1;

          if (petsc_wall_site(solver, ic, jc, kc)) {
            v[0] = solver->block->wall_diag;
            MatSetValuesStencil(solver->block->a, 1, &row, 1, &row, v,
                                INSERT_VALUES);
          }
          else {
            double diag = 0.0;
            int np = 1;                    /* entry 0 is the diagonal */
            col[0] = row;
            for (int p = 1; p < s->npoints; p++) {
              int wall = petsc_wall_site(solver, ic + s->cv[p][X],
                                         jc + s->cv[p][Y], kc + s->cv[p][Z]);
              double coef = s->wlaplacian[p] * epsilon * (1.0 + wall);
              diag -= coef;
              if (!wall) {
                col[np].i = i + s->cv[p][X];
                col[np].j = j + s->cv[p][Y];
                col[np].k = k + s->cv[p][Z];
                v[np] = coef;
                np += 1;
              }
            }
            v[0] = diag;
            MatSetValuesStencil(solver->block->a, 1, &row, np, col, v,
                                INSERT_VALUES);
          }
          continue;
        }
        /*CHANGE END - 20260926 */

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

  /*CHANGE INIT - 20260926 constant-potential walls: a Dirichlet wall makes
   * the operator non-singular, so the constant-mode nullspace must go.
   * Original (unconditional):
   *   MatNullSpaceCreate(...); MatSetNullSpace(a, nullsp); ... */
  if (solver->block->map) {
    MatSetNullSpace(solver->block->a, NULL);
  }
  else
  /*CHANGE END - 20260926 */
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
  /*CHANGE INIT - 20260926 ... but not with walls: there is no nullspace then,
   * and projecting out the mean would shift the prescribed wall potential. */
  if (solver->block->map == NULL)
  /*CHANGE END - 20260926 */
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

  /*CHANGE INIT - 20260926 constant-potential walls: induced charge */
  if (solver->block->map) {
    psi_solver_petsc_wall_charge_compute(solver);
    if (ntimestep % solver->psi->solver.nfreq == 0) {
      pe_t* pe = solver->psi->pe;
      double* q = solver->block->wall_q;
      pe_info(pe, "Wall charge lower %22.15e upper %22.15e\n", q[0], q[1]);
      pe_info(pe, "Fluid charge      %22.15e total %22.15e\n",
              solver->block->fluid_q, q[0] + q[1] + solver->block->fluid_q);
    }
  }
  /*CHANGE END - 20260926 */

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
  double fluid_q_local = 0.0;   /*CHANGE 20260926 constant-potential walls */

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

        /*CHANGE INIT - 20260926 constant-potential walls.
         * Wall site: the row is D x = b, so b = D psi_wall.
         * Fluid site: move each eliminated wall column to the RHS,
         * b -= A_fw psi_wall with A_fw = 2 wlaplacian[p] epsilon. */
        if (solver->block->map) {
          if (petsc_wall_site(solver, ic, jc, kc)) {
            rho_3d[k][j][i] = solver->block->wall_diag
                            * petsc_wall_value(solver, ic, jc, kc);
            continue;
          }
          {
            stencil_t* s = solver->psi->stencil;
            double epsilon = 0.0;
            psi_epsilon(solver->psi, &epsilon);
            psi_rho_elec(solver->psi, index, &rho_elec);
            rho_3d[k][j][i] = rho_elec * eunit * beta;
            fluid_q_local += rho_elec;
            for (int p = 1; p < s->npoints; p++) {
              int icn = ic + s->cv[p][X];
              int jcn = jc + s->cv[p][Y];
              int kcn = kc + s->cv[p][Z];
              if (petsc_wall_site(solver, icn, jcn, kcn)) {
                rho_3d[k][j][i] -= 2.0 * s->wlaplacian[p] * epsilon
                                 * petsc_wall_value(solver, icn, jcn, kcn);
              }
            }
          }
          continue;
        }
        /*CHANGE END - 20260926 */

        psi_rho_elec(solver->psi, index, &rho_elec);
        rho_3d[k][j][i] = rho_elec * eunit * beta;
      }
    }
  }

  /*CHANGE INIT - 20260926 constant-potential walls: total fluid charge, for
   * the neutrality check wall_q[0] + wall_q[1] + fluid_q = 0. An external
   * field e0 is refused in psi_solver_petsc_wall_set(): the code below
   * assumes a periodic system. */
  if (solver->block->map) {
    MPI_Comm comm = MPI_COMM_NULL;
    cs_cart_comm(cs, &comm);
    MPI_Allreduce(&fluid_q_local, &solver->block->fluid_q, 1, MPI_DOUBLE,
                  MPI_SUM, comm);
  }
  /*CHANGE END - 20260926 */

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

/*CHANGE INIT - 20260926 constant-potential walls */
/*****************************************************************************
 *
 *  psi_solver_petsc_wall_set
 *
 *  Turn MAP_BOUNDARY sites into constant-potential walls (Asta et al. 2019,
 *  arXiv:1907.04732, sec. II-C). Sites in the lower half of the box along
 *  `axis` are held at psi_wall[0], the upper half at psi_wall[1], in the same
 *  units as the psi field. Must be called once the map is final (it is
 *  created after the solver). The matrix is rebuilt here, once: the walls do
 *  not move.
 *
 *  Returns 0 on success; -1 variable permittivity (not supported), -2 an
 *  external field e0 is set (rhs_set assumes periodicity for it), -3 the
 *  map has no MAP_BOUNDARY site.
 *
 *****************************************************************************/

int psi_solver_petsc_wall_set(psi_solver_petsc_t* solver, map_t* map,
                              int axis, const double psi_wall[2]) {

  assert(solver);
  assert(map);
  assert(axis == X || axis == Y || axis == Z);

  if (solver->super.impl != &vt_) return -1;
  if (solver->psi->e0[X] || solver->psi->e0[Y] || solver->psi->e0[Z]) {
    return -2;
  }

  {
    int nwall = 0;
    map_volume_allreduce(map, MAP_BOUNDARY, &nwall);
    if (nwall == 0) return -3;
  }

  /* Neighbour status is read in the halo */
  map_halo(map);

  {
    double epsilon = 0.0;
    stencil_t* s = solver->psi->stencil;
    psi_epsilon(solver->psi, &epsilon);

    solver->block->map = map;
    solver->psi->wall_map = map;   /* same wall rule in psi_electric_field */
    solver->block->wall_axis = axis;
    solver->block->wall_psi[0] = psi_wall[0];
    solver->block->wall_psi[1] = psi_wall[1];
    /* Same size as a fluid diagonal, so Jacobi/AMG see comparable rows */
    solver->block->wall_diag = epsilon * s->wlaplacian[0];
  }

  /* The field itself holds the prescribed value at wall sites: it is the
   * initial guess, and what any code reading psi sees before a solve. */
  {
    int nlocal[3] = { 0 };
    cs_nlocal(solver->psi->cs, nlocal);
    for (int ic = 1; ic <= nlocal[X]; ic++) {
      for (int jc = 1; jc <= nlocal[Y]; jc++) {
        for (int kc = 1; kc <= nlocal[Z]; kc++) {
          if (petsc_wall_site(solver, ic, jc, kc)) {
            int index = cs_index(solver->psi->cs, ic, jc, kc);
            psi_psi_set(solver->psi, index,
                        petsc_wall_value(solver, ic, jc, kc));
          }
        }
      }
    }
  }

  psi_solver_petsc_matrix_set(solver);
  KSPSetUp(solver->block->ksp);

  return 0;
}

/*****************************************************************************
 *
 *  psi_solver_petsc_wall_charge_compute
 *
 *  Induced charge on each wall, in the units of rho_elec (summed over
 *  sites), from Gauss's law through the wall mid-plane:
 *
 *    q = (epsilon / e beta) sum_{links wall w -> fluid f} 2 wlaplacian[p]
 *                                                     (psi_f - psi_w)
 *
 *  This is the paper's eq. 16 with the SYMMETRIC boundary rule (factor 2 on
 *  every link joining a wall and a fluid site, seen from the wall side).
 *  Their compact eq. 15 read literally at a wall site gives exactly zero,
 *  and the unmodified eq. 12 gives half the charge; both were checked
 *  numerically (2026-09-26). By construction q[0] + q[1] + fluid_q = 0 to
 *  solver tolerance. Wall-wall links are skipped: across the periodic
 *  boundary they join the two electrodes, which is not a physical link.
 *
 *****************************************************************************/

int psi_solver_petsc_wall_charge_compute(psi_solver_petsc_t* solver) {

  int xs, ys, zs, xw, yw, zw;
  int offset[3] = { 0 };
  double q_local[2] = { 0.0, 0.0 };
  double epsilon = 0.0;
  double eunit = solver->psi->e;
  double beta = solver->psi->beta;
  stencil_t* s = solver->psi->stencil;
  const double*** x3 = NULL;
  Vec xl;

  assert(solver);
  assert(solver->block->map);

  psi_epsilon(solver->psi, &epsilon);
  cs_nlocal_offset(solver->psi->cs, offset);
  DMDAGetCorners(solver->block->da, &xs, &ys, &zs, &xw, &yw, &zw);

  /* Ghosted copy of the solution (stencil width 1), read on the host as
   * in psi_solver_petsc_da_to_psi() */
  VecBindToCPU(solver->block->x, PETSC_TRUE);
  DMGetLocalVector(solver->block->da, &xl);
  VecBindToCPU(xl, PETSC_TRUE);
  DMGlobalToLocalBegin(solver->block->da, solver->block->x, INSERT_VALUES, xl);
  DMGlobalToLocalEnd(solver->block->da, solver->block->x, INSERT_VALUES, xl);
  DMDAVecGetArrayRead(solver->block->da, xl, &x3);

  for (int k = zs; k < zs + zw; k++) {
    int kc = k - offset[Z] + 1;
    for (int j = ys; j < ys + yw; j++) {
      int jc = j - offset[Y] + 1;
      for (int i = xs; i < xs + xw; i++) {
        int ic = i - offset[X] + 1;
        if (!petsc_wall_site(solver, ic, jc, kc)) continue;
        {
          int side = petsc_wall_side(solver, ic, jc, kc);
          for (int p = 1; p < s->npoints; p++) {
            int icn = ic + s->cv[p][X];
            int jcn = jc + s->cv[p][Y];
            int kcn = kc + s->cv[p][Z];
            if (petsc_wall_site(solver, icn, jcn, kcn)) continue;
            q_local[side] += 2.0 * s->wlaplacian[p]
              * (x3[k + s->cv[p][Z]][j + s->cv[p][Y]][i + s->cv[p][X]]
                 - x3[k][j][i]);
          }
        }
      }
    }
  }

  DMDAVecRestoreArrayRead(solver->block->da, xl, &x3);
  VecBindToCPU(xl, PETSC_FALSE);
  DMRestoreLocalVector(solver->block->da, &xl);
  VecBindToCPU(solver->block->x, PETSC_FALSE);

  q_local[0] *= epsilon / (eunit * beta);
  q_local[1] *= epsilon / (eunit * beta);

  {
    MPI_Comm comm = MPI_COMM_NULL;
    cs_cart_comm(solver->psi->cs, &comm);
    MPI_Allreduce(q_local, solver->block->wall_q, 2, MPI_DOUBLE, MPI_SUM,
                  comm);
  }

  return 0;
}

/*****************************************************************************
 *
 *  psi_solver_petsc_wall_charge
 *
 *  Charge on the lower/upper wall and in the fluid at the last solve.
 *
 *****************************************************************************/

int psi_solver_petsc_wall_charge(const psi_solver_petsc_t* solver,
                                 double q[2], double* q_fluid) {

  assert(solver);

  if (solver->block->map == NULL) return -1;

  q[0] = solver->block->wall_q[0];
  q[1] = solver->block->wall_q[1];
  if (q_fluid) *q_fluid = solver->block->fluid_q;

  return 0;
}
/*CHANGE END - 20260926 */

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
