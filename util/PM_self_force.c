/*****************************************************************************
 *
 *  PM_self_force.c
 *
 *  Compute the particle-mesh (PM) self-force on a single charged subgrid
 *  particle using the Peskin4 spread/gather kernel and an FFT Poisson solver.
 *
 *  For each position on a grid covering one voxel (between lattice nodes),
 *  the program:
 *    1. Places a particle with charge q at position r.
 *    2. Spreads charge to the mesh with the Peskin4 kernel.
 *    3. Solves Poisson's equation via FFT.
 *    4. Computes the electric field E = -grad(psi).
 *    5. Gathers E back to the particle with the same Peskin4 kernel.
 *    6. Computes the self-force F = q * E_gathered.
 *
 *  Output: tab-separated table to stdout:
 *    dx  dy  dz  Fx  Fy  Fz  |F|
 *  where (dx,dy,dz) is the fractional offset within the voxel [0,1).
 *
 *  Parameters (hard-coded, change as needed):
 *    N         system size (cubic)
 *    NSTEPS    grid points per axis within one voxel
 *    Q0, Q1    species charges on the particle
 *    EPSILON   permittivity
 *    BETA      1/(k_B T)
 *
 *  Compilation (requires the CUDA+MPI build active in config.mk):
 *    The FFT Poisson solver uses cuFFT, so a GPU build is needed.
 *    From the util/ directory, with config.mk pointing to the CUDA config:
 *
 *    nvcc -O2 -I../src -I../mpi_s -I../target \
 *         -I/usr/local/ompi/include \
 *         PM_self_force.c ../src/libludwig.a \
 *         -L/usr/local/ompi/lib -lmpi \
 *         -lcufft -lcuda -lm -o PM_self_force
 *
 *  Or build via the project Makefile after adding a rule for this file.
 *
 *  CHANGE INIT - 20260518 PM self-force utility using Peskin4 spread/gather
 *
 *****************************************************************************/

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../src/pe.h"
#include "../src/coords.h"
#include "../src/colloids.h"
#include "../src/psi.h"
#include "../src/psi_options.h"
#include "../src/psi_fft.h"
#include "../src/subgrid.h"
#include "../src/physics.h"

/* ------------------------------------------------------------------
 * User-adjustable parameters
 * ------------------------------------------------------------------ */

/* System size (cubic lattice N x N x N). Must be large enough that the
 * Peskin4 kernel (support radius 2) fits without boundary issues.
 * Minimum recommended: 8.  Must be even for FFT. */
#define N      16

/* Number of sample points per axis inside one voxel (0 to 1 exclusive).
 * NSTEPS=5 gives offsets 0, 0.2, 0.4, 0.6, 0.8.
 * NSTEPS=4 gives offsets 0, 0.25, 0.5, 0.75. */
#define NSTEPS  5

/* Particle charge species: q0 > 0 (cation), q1 = 0 (no counter-charge) */
#define Q0      1.0
#define Q1      0.0

/* Permittivity (lattice units) */
#define EPSILON 1.0

/* Boltzmann factor beta = 1/(k_B T) */
#define BETA    1.0

/* Unit charge */
#define EUNIT   1.0

/* Zero the potential field on the lattice (no library equivalent) */
static int sf_psi_zero(psi_t *obj) {
    int nsites = obj->nsites;
    for (int i = 0; i < nsites; i++) {
        psi_psi_set(obj, i, 0.0);
    }
    return 0;
}

/* Map a kernel name to the enum. Returns -1 if unknown. */
static int sf_kernel_from_name(const char *name, subgrid_kernel_t *kernel) {
    if      (strcmp(name, "peskin4")  == 0) *kernel = SUBGRID_KERNEL_PESKIN4;
    else if (strcmp(name, "peskin6")  == 0) *kernel = SUBGRID_KERNEL_PESKIN6;
    else if (strcmp(name, "bspline4") == 0) *kernel = SUBGRID_KERNEL_BSPLINE4;
    else if (strcmp(name, "bspline6") == 0) *kernel = SUBGRID_KERNEL_BSPLINE6;
    else if (strcmp(name, "kb4")      == 0) *kernel = SUBGRID_KERNEL_KB4;
    /*CHANGE INIT - 20260706 Hann kernel in PM self-force utility */
    else if (strcmp(name, "hann")     == 0) *kernel = SUBGRID_KERNEL_HANN;
    /*CHANGE END - 20260706 Hann kernel in PM self-force utility */
    else return -1;
    return 0;
}

/* ------------------------------------------------------------------
 * Main
 * ------------------------------------------------------------------ */

int main(int argc, char **argv) {

    pe_t             *pe       = NULL;
    cs_t             *cs       = NULL;
    psi_t            *psi      = NULL;
    colloids_info_t  *cinfo    = NULL;
    physics_t        *phys     = NULL;
    psi_solver_fft_t *fft      = NULL;
    FILE             *fp_diag  = NULL;

    int ntotal[3]    = {N, N, N};
    int periodic[3]  = {1, 1, 1};
    int decomp[3]    = {1, 1, 1};
    int ncell[3]     = {4, 4, 4};   /* Colloid cell list */

    /* Kernel selectable from the command line: argv[1] in
     * {peskin4, peskin6, bspline4, bspline6, kb4, hann}. Default: peskin4.
     * For hann, argv[2] (optional) is the kernel order n = full support
     * width in lattice units (default 4, cf. subgrid_hann_order input key). */
    subgrid_kernel_t kernel = SUBGRID_KERNEL_PESKIN4;
    const char *kernel_name = "peskin4";
    if (argc > 1) {
        if (sf_kernel_from_name(argv[1], &kernel) != 0) {
            fprintf(stderr, "Unknown kernel '%s' "
                    "(use peskin4|peskin6|bspline4|bspline6|kb4|hann)\n",
                    argv[1]);
            return 1;
        }
        kernel_name = argv[1];
    }
    /*CHANGE INIT - 20260706 Hann kernel in PM self-force utility */
    int hann_order = 4;
    if (kernel == SUBGRID_KERNEL_HANN) {
        if (argc > 2) {
            hann_order = atoi(argv[2]);
            if (hann_order <= 0) {
                fprintf(stderr, "hann order must be a positive integer "
                        "(got '%s')\n", argv[2]);
                return 1;
            }
        }
        subgrid_set_hann_order((double) hann_order);
    }
    /* Halo width from the kernel support radius (order-dependent for hann) */
    /* int nhalo = (kernel == SUBGRID_KERNEL_PESKIN6 ||
                    kernel == SUBGRID_KERNEL_BSPLINE6) ? 3 : 2; */
    int nhalo = subgrid_get_range(kernel);
    if (nhalo < 2) nhalo = 2;
    /*CHANGE END - 20260706 Hann kernel in PM self-force utility */

    /* Centre of the lattice (particle will be placed near here) */
    double r_centre[3] = {0.5*N + 1.0, 0.5*N + 1.0, 0.5*N + 1.0};

    /* -------------------------------------------------------------- */
    /* 1.  Parallel environment and coordinate system                  */
    /* -------------------------------------------------------------- */

    MPI_Init(&argc, &argv);

    pe_create(MPI_COMM_WORLD, PE_QUIET, &pe);

    cs_create(pe, &cs);
    cs_ntotal_set(cs, ntotal);
    cs_periodicity_set(cs, periodic);
    cs_decomposition_set(cs, decomp);
    cs_nhalo_set(cs, nhalo);   /* support radius: 2 (4-point) or 3 (6-point) */
    cs_init(cs);
    cs_commit(cs);

    /* -------------------------------------------------------------- */
    /* 2.  Electrokinetic potential (psi)                              */
    /* -------------------------------------------------------------- */

    psi_options_t opts = psi_options_default(nhalo);
    opts.nk            = 2;
    opts.e             = EUNIT;
    opts.beta          = BETA;
    opts.epsilon1      = EPSILON;
    opts.epsilon2      = EPSILON;
    opts.valency[0]    =  1;
    opts.valency[1]    = -1;
    opts.diffusivity[0] = 0.01;
    opts.diffusivity[1] = 0.01;

    psi_create(pe, cs, &opts, &psi);
    assert(psi);

    /* -------------------------------------------------------------- */
    /* 3.  FFT Poisson solver                                          */
    /* -------------------------------------------------------------- */

    psi_solver_fft_create(psi, &fft);
    assert(fft);

    /* -------------------------------------------------------------- */
    /* 4.  Physics (needed only for kt, which we set explicitly)       */
    /* -------------------------------------------------------------- */

    physics_create(pe, &phys);
    physics_kt_set(phys, 1.0 / BETA);

    /* -------------------------------------------------------------- */
    /* 5.  Colloid info (cell list)                                    */
    /* -------------------------------------------------------------- */

    colloids_info_create(pe, cs, ncell, &cinfo);
    assert(cinfo);

    /* -------------------------------------------------------------- */
    /* 6.  Output CSV and stdout header                                */
    /* -------------------------------------------------------------- */

    fp_diag = fopen("PM_self_force_diag.csv", "w");
    assert(fp_diag);
    fprintf(fp_diag, "dx;dy;dz;Fx;Fy;Fz;Fmod\n");

    /* fp_subgrid: discard internal diagnostics from subgrid_update_Esub */
    FILE *fp_subgrid = fopen("/dev/null", "w");
    assert(fp_subgrid);

    /*CHANGE INIT - 20260706 Hann kernel in PM self-force utility */
    if (kernel == SUBGRID_KERNEL_HANN) {
        printf("# PM self-force: %s (order %d) spread/gather + FFT Poisson\n",
               kernel_name, hann_order);
    }
    else
    /*CHANGE END - 20260706 Hann kernel in PM self-force utility */
    printf("# PM self-force: %s spread/gather + FFT Poisson\n", kernel_name);
    printf("# System: %d x %d x %d,  Q0=%.4g Q1=%.4g epsilon=%.4g beta=%.4g nhalo=%d\n",
           N, N, N, Q0, Q1, EPSILON, BETA, nhalo);
    printf("# Columns: dx dy dz  Fx Fy Fz  |F|\n");
    printf("%-12s %-12s %-12s   %-15s %-15s %-15s   %-15s\n",
           "dx", "dy", "dz", "Fx", "Fy", "Fz", "|F|");

    /* -------------------------------------------------------------- */
    /* 7.  Loop over fractional offsets within one voxel               */
    /* -------------------------------------------------------------- */

    for (int ix = 0; ix < NSTEPS; ix++) {
        for (int iy = 0; iy < NSTEPS; iy++) {
            for (int iz = 0; iz < NSTEPS; iz++) {

                double dx = (double)ix / NSTEPS;
                double dy = (double)iy / NSTEPS;
                double dz = (double)iz / NSTEPS;

                /* Particle position: centre of lattice + fractional offset */
                double r[3];
                r[0] = r_centre[0] + dx;
                r[1] = r_centre[1] + dy;
                r[2] = r_centre[2] + dz;

                /* -------------------------------------------------- */
                /* 7a. Reset lattice charge and potential               */
                /* -------------------------------------------------- */

                psi_rho_zero(psi);   /* from psi.h: zeroes all charge species */
                sf_psi_zero(psi);    /* local helper: zeroes potential */

                /* -------------------------------------------------- */
                /* 7b. Add (or replace) the single particle             */
                /* -------------------------------------------------- */

                /* Remove any previously added colloid by re-creating cinfo */
                colloids_info_free(cinfo);
                colloids_info_create(pe, cs, ncell, &cinfo);

                colloid_t *pc = NULL;
                colloids_info_add_local(cinfo, 1, r, &pc);
                if (pc) {
                    pc->s.bc  = COLLOID_BC_SUBGRID;
                    pc->s.q0  = Q0;
                    pc->s.q1  = Q1;
                    pc->Esub[0] = 0.0;
                    pc->Esub[1] = 0.0;
                    pc->Esub[2] = 0.0;
                    pc->fex[0]  = 0.0;
                    pc->fex[1]  = 0.0;
                    pc->fex[2]  = 0.0;
                }

                /* Set nsubgrid so that subgrid functions don't early-return */
                cinfo->nsubgrid = 1;

                /* -------------------------------------------------- */
                /* 7c. Spread charge from particle to lattice (Peskin4) */
                /* -------------------------------------------------- */

                distributed_charge_klein_t *charge = NULL;
                subgrid_charge_from_particles(cinfo, psi, &charge, kernel);
                subgrid_free_distributed_charge_t(&charge);

                /* -------------------------------------------------- */
                /* 7d. Halo exchange for charge, then solve Poisson    */
                /* -------------------------------------------------- */

                psi_halo_rho(psi);
                psi_solver_fft_solve(fft, 0);

                /* -------------------------------------------------- */
                /* 7e. Compute electric field E = -grad(psi)           */
                /* -------------------------------------------------- */

                psi_halo_psi(psi);
                psi_compute_electric_field(psi);

                /* -------------------------------------------------- */
                /* 7f. Gather E onto particle (Peskin4)                 */
                /* -------------------------------------------------- */

                subgrid_update_Esub(cinfo, psi, 0, fp_subgrid, pe, kernel);

                /* -------------------------------------------------- */
                /* 7g. Compute self-force F = q_net * E_gathered       */
                /* -------------------------------------------------- */

                double Fx = 0.0, Fy = 0.0, Fz = 0.0;
                if (pc) {
                    double q_net = pc->s.q0 - pc->s.q1;
                    Fx = q_net * pc->Esub[0];
                    Fy = q_net * pc->Esub[1];
                    Fz = q_net * pc->Esub[2];
                }

                double Fmag = sqrt(Fx*Fx + Fy*Fy + Fz*Fz);

                /* -------------------------------------------------- */
                /* 7h. Print row                                        */
                /* -------------------------------------------------- */

                printf("%-12.6f %-12.6f %-12.6f   %-15.8e %-15.8e %-15.8e   %-15.8e\n",
                       dx, dy, dz, Fx, Fy, Fz, Fmag);
                fprintf(fp_diag, "%.6f;%.6f;%.6f;%.15e;%.15e;%.15e;%.15e\n",
                        dx, dy, dz, Fx, Fy, Fz, Fmag);
            }
        }
    }

    /* -------------------------------------------------------------- */
    /* 8.  Cleanup                                                     */
    /* -------------------------------------------------------------- */

    fclose(fp_diag);
    fclose(fp_subgrid);

    colloids_info_free(cinfo);
    psi_solver_fft_free(&fft);
    psi_free(&psi);
    physics_free(phys);
    cs_free(cs);
    pe_free(pe);

    MPI_Finalize();

    return 0;
}

/* CHANGE END - 20260518 PM self-force utility using Peskin4 spread/gather */
