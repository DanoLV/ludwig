/*****************************************************************************
 *
 *  ewald_charge.c
 *
 *  The Ewald summation for point charges (Coulomb/Poisson).
 *  Adapted from ewald.c (magnetic dipoles).
 *
 *  This implementation handles BOTH:
 *    1. Point charges on colloids (discrete particles)
 *    2. Charge density on the LB lattice nodes (continuous field)
 *
 *  The total electrostatic potential is:
 *    phi(r) = sum_j q_j / (4*pi*epsilon*|r - r_j|)
 *           + integral rho(r') / (4*pi*epsilon*|r - r'|) d^3r'
 *
 *  For the lattice, we treat each node as a point charge with q = rho_elec * dV
 *  where dV = 1 in lattice units.
 *
 *  The Ewald sum splits this into:
 *    phi = phi_real + phi_fourier - phi_self
 *
 *  See Allen and Tildesley, Computer Simulation of Liquids, Chapter 5.5
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  Contributing authors:
 *  Based on ewald.c by Grace Kim and Kevin Stratford (kevin@epcc.ed.ac.uk)
 *
 *****************************************************************************/

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stddef.h>

#include "pe.h"
#include "coords.h"
#include "colloids.h"
#include "hydro.h"
#include "field.h"
#include "ewald_charge.h"
#include "timer.h"
#include "util.h"
#include "util_sum.h"

static int ewald_charge_on_ = 0;
static int nk_[3];
static int nkmax_;
static int nktot_;
// static int nk_limit_ = 3000;      /* Max nk allowed (0 = no limit) */
static int nk_limit_ = 000;      /* Max nk allowed (0 = no limit) */
static double ewald_rc_ = 0.0;
static double alpha_ = 0.0;       /* Ewald splitting parameter */
static double eps_reg_ = 0.0;    /* Real-space regularization parameter */
static double kmax_;
static double rpi_;
static double epsilon_ = 1.0;     /* Permittivity */
static double epsilon_prime_ = 1.0; /* Dielectric constant of surrounding medium (Deserno & Holm convention) */
static double beta_ = 1.0;        /* 1/kT */
static double eunit_ = 1.0;       /* Unit charge */

static double ereal_ = 0.0;
static double efourier_ = 0.0;

static double* sinx_;   /* The term S(k) for each k: sum_i q_i * sin(k.r_i) */
static double* cosx_;   /* The term C(k) for each k: sum_i q_i * cos(k.r_i) */
static double* sinkr_;  /* Table for sin(kr) values */
static double* coskr_;  /* Table for cos(kr) values */

struct ewald_charge_s {
  pe_t* pe;                 /* Parallel environment */
  cs_t* cs;                 /* Coordinate system */
  colloids_info_t* cinfo;   /* Colloid information (can be NULL) */
  psi_t* psi;               /* Charge density field (can be NULL) */
  map_t* map;               /* Map for solid/fluid (can be NULL) */
  hydro_t* hydro;           /* Hydrodynamics for reaction forces (can be NULL) */
  ewald_charge_source_t sources;  /* Which sources to include */
  double epsilon_prime;      /* Dielectric constant of surrounding medium (Deserno & Holm convention) */
  double kappa_debye_stored; /* Stored kappa for self-field helper (set by self_build functions) */

  /* CHANGE INIT - PoissonVerification - Arrays for face sampling, Poisson check, Maxwell stress force */
  double* E_face;          /* [Nnodes * 6*N*N * 3] E field at face quadrature points
                            * Layout per node: face order -X,+X,-Y,+Y,-Z,+Z; N*N points each; 3 components */
  double* div_E;           /* [Nnodes] numerical div(E) at each node via Gauss theorem */
  double* force_divstress; /* [Nnodes * 3] force from Maxwell stress tensor divergence */
  double* poisson_error;   /* [Nnodes] |div(E) - rho_eff/epsilon| at each node */
  int* node_flags;      /* [Nnodes] EWALD_NODE_NORMAL or EWALD_NODE_NEAR_PARTICLE (diagnostic) */
  /* CHANGE END - PoissonVerification */
};

static int ewald_charge_sum_sin_cos_terms(ewald_charge_t* ewald);
static int ewald_charge_sum_sin_cos_colloids(ewald_charge_t* ewald);
static int ewald_charge_sum_sin_cos_lattice(ewald_charge_t* ewald);
static int ewald_charge_get_number_fourier_terms(ewald_charge_t* ewald);
static int ewald_charge_set_kr_table(ewald_charge_t* ewald, double r[3]);

/* CHANGE INIT - PoissonVerification - Forward declarations for face-sampling helper functions */
static int ewald_flag_near_particle_nodes(ewald_charge_t* ewald, int nlocal[3], int noffset[3]);
static int ewald_efield_at_r(ewald_charge_t* ewald, const double r[3], double E[3]);
static int ewald_face_sample_cpu(ewald_charge_t* ewald, int nlocal[3], int noffset[3]);
static int ewald_calc_div_and_force(ewald_charge_t* ewald, int nlocal[3], int noffset[3]);
static int ewald_write_poisson_verification(ewald_charge_t* ewald, int nlocal[3], int noffset[3], int step);
static int ewald_write_force_divstress(ewald_charge_t* ewald, int nlocal[3], int noffset[3], int step);
static int ewald_write_poisson_stencil(ewald_charge_t* ewald, int nlocal[3], int noffset[3], int step);
static int ewald_write_force_fluid(ewald_charge_t* ewald, int nlocal[3], int noffset[3], int step);
/* CHANGE END - PoissonVerification */

/*****************************************************************************
 *
 *  ewald_charge_create
 *
 *  Create Ewald sum for point charges and/or lattice charge density.
 *
 *  Parameters:
 *    epsilon   - permittivity of the medium
 *    rc        - real space cutoff
 *    alpha     - Ewald splitting parameter (if 0, use default 5/(2*rc))
 *    cinfo     - colloid information (can be NULL if no colloids)
 *    psi       - charge density field (can be NULL if no lattice charges)
 *    map       - map for solid/fluid nodes (can be NULL)
 *    sources   - which charge sources to include
 *
 *****************************************************************************/

int ewald_charge_create(pe_t* pe, cs_t* cs, double rc,
                        double alpha, double epsilon_prime,
                        colloids_info_t* cinfo, psi_t* psi,
                        map_t* map, hydro_t* hydro, ewald_charge_source_t sources,
                        ewald_charge_t** pewald) {
  int nk;
  double ltot[3];
  PI_DOUBLE(pi);
  ewald_charge_t* ewald = NULL;

  assert(pe);
  assert(cs);
  assert(pewald);

  /* Validate: at least one source must be specified */
  // if (sources & EWALD_SOURCE_COLLOIDS) assert(cinfo);
  if (sources & EWALD_SOURCE_LATTICE) assert(psi);

  ewald = (ewald_charge_t*)calloc(1, sizeof(ewald_charge_t));
  assert(ewald);
  if (ewald == NULL) pe_fatal(pe, "calloc(ewald_charge) failed");

  ewald->pe = pe;
  ewald->cs = cs;
  ewald->cinfo = cinfo;
  ewald->psi = psi;
  ewald->map = map;
  ewald->hydro = hydro;
  ewald->sources = sources;

  cs_ltot(cs, ltot);

  /* Get parameters from psi */
  psi_epsilon(psi, &epsilon_);
  psi_beta(psi, &beta_);
  psi_unit_charge(psi, &eunit_);

  /* Set constants */
  rpi_ = 1.0 / sqrt(pi);
  epsilon_prime_ = epsilon_prime;
  ewald_rc_ = rc;
  ewald_charge_on_ = 1;

  /* Ewald parameter: use provided alpha or default */
  if (alpha > 0.0) {
    alpha_ = alpha;
  }
  else {
    // alpha_ = 5.0/(2.0*ewald_rc_);  /* Default as in dipole version */
    // alpha_ = 5.0/(1.0*ewald_rc_);  /* mas terminos de fourier*/
    alpha_ = 5.0 * 0.7 / ewald_rc_;  /* mas terminos de fourier*/
  }

  /* Number of k-vectors needed for convergence */
  nk = ceil(alpha_ * alpha_ * ewald_rc_ * ltot[X] / pi);

  /* Limit maximum nk to avoid excessive computation */
  /* nk_limit_ = 0 means no limit */
  if (nk_limit_ > 0 && nk > nk_limit_) {
    pe_info(pe, "Ewald charge: limiting nk from %d to %d\n", nk, nk_limit_);
    nk = nk_limit_;
  }

  nk_[X] = nk;
  nk_[Y] = nk;
  nk_[Z] = nk;
  kmax_ = pow(2.0 * pi * nk / ltot[X], 2);
  nkmax_ = nk + 1;
  nktot_ = ewald_charge_get_number_fourier_terms(ewald);
  assert(nktot_ > 0);

  sinx_ = (double*)malloc(nktot_ * sizeof(double));
  cosx_ = (double*)malloc(nktot_ * sizeof(double));

  if (sinx_ == NULL) pe_fatal(pe, "Ewald charge sum malloc(sinx_) failed\n");
  if (cosx_ == NULL) pe_fatal(pe, "Ewald charge sum malloc(cosx_) failed\n");

  sinkr_ = (double*)malloc(3 * nkmax_ * sizeof(double));
  coskr_ = (double*)malloc(3 * nkmax_ * sizeof(double));

  if (sinkr_ == NULL) pe_fatal(pe, "Ewald charge sum malloc(sinkr_) failed\n");
  if (coskr_ == NULL) pe_fatal(pe, "Ewald charge sum malloc(coskr_) failed\n");

  /* Default: vacuum boundary conditions (epsilon' = 1) as in Deserno & Holm */
  ewald->epsilon_prime = epsilon_prime_;

  // /* CHANGE INIT - PoissonVerification - Allocate face-sampling arrays */
  // {
  //   int nlocal[3];
  //   cs_nlocal(cs, nlocal);
  //   int ntotal_nodes = nlocal[X] * nlocal[Y] * nlocal[Z];
  //   int ns = EWALD_NSAMPLE_FACE;
  //   int nsample_per_node = 6 * ns * ns;

  //   ewald->E_face = (double*)calloc(ntotal_nodes * nsample_per_node * 3, sizeof(double));
  //   ewald->div_E = (double*)calloc(ntotal_nodes, sizeof(double));
  //   ewald->force_divstress = (double*)calloc(ntotal_nodes * 3, sizeof(double));
  //   ewald->poisson_error = (double*)calloc(ntotal_nodes, sizeof(double));
  //   ewald->node_flags = (int*)calloc(ntotal_nodes, sizeof(int));

  //   if (!ewald->E_face || !ewald->div_E || !ewald->force_divstress ||
  //       !ewald->poisson_error || !ewald->node_flags) {
  //     pe_fatal(pe, "calloc failed for Poisson verification arrays\n");
  //   }
  // }
  // /* CHANGE END - PoissonVerification */

  *pewald = ewald;

  return 0;
}

/*****************************************************************************
 *
 *  ewald_charge_set_eps_reg
 *
 *  Set the real-space regularization parameter eps_reg.
 *  Replaces 1/r -> 1/sqrt(r^2 + eps_reg^2) in all real-space GPU kernels.
 *  Default is 0.0 (no regularization).
 *
 *****************************************************************************/

int ewald_charge_set_eps_reg(double eps_reg) {
  eps_reg_ = eps_reg;
  return 0;
}

/*****************************************************************************
 *
 *  ewald_charge_set_epsilon_prime
 *
 *  Set the dielectric constant of the surrounding medium (epsilon').
 *  Uses Deserno & Holm convention (J. Chem. Phys. 1998).
 *
 *  The dipole correction is: E^(d) = 2*pi / ((1 + 2*epsilon')*V) * |M|^2
 *
 *  Special cases:
 *    epsilon' = 1      : vacuum boundary conditions (default)
 *    epsilon' -> infty : metallic/tinfoil boundary conditions (no correction)
 *
 *****************************************************************************/

int ewald_charge_set_epsilon_prime(ewald_charge_t* ewald, double epsilon_prime) {

  assert(ewald);
  assert(epsilon_prime >= 1.0);  /* epsilon' >= 1 for physical media */

  ewald->epsilon_prime = epsilon_prime;

  return 0;
}

/*****************************************************************************
 *
 *  ewald_charge_set_cinfo
 *
 *  CHANGE INIT - 20260513 setter for colloid info pointer
 *
 *****************************************************************************/

int ewald_charge_set_cinfo(ewald_charge_t * ewald, colloids_info_t * cinfo) {

  assert(ewald);

  ewald->cinfo = cinfo;
  return 0;
}

int ewald_charge_set_map(ewald_charge_t * ewald, map_t * map) {

  assert(ewald);

  ewald->map = map;
  return 0;
}

int ewald_charge_set_hydro(ewald_charge_t * ewald, hydro_t * hydro) {

  assert(ewald);

  ewald->hydro = hydro;
  return 0;
}

/* CHANGE END - 20260513 */

/*****************************************************************************
 *
 *  ewald_charge_get_epsilon_prime
 *
 *****************************************************************************/

int ewald_charge_get_epsilon_prime(ewald_charge_t* ewald, double* epsilon_prime) {

  assert(ewald);
  assert(epsilon_prime);

  *epsilon_prime = ewald->epsilon_prime;

  return 0;
}

/*****************************************************************************
 *
 *  ewald_charge_free
 *
 *****************************************************************************/

int ewald_charge_free(ewald_charge_t* ewald) {

  assert(ewald);

  free(sinx_);
  free(cosx_);
  free(sinkr_);
  free(coskr_);

  /* CHANGE INIT - PoissonVerification - Free face-sampling arrays */
  free(ewald->E_face);
  free(ewald->div_E);
  free(ewald->force_divstress);
  free(ewald->poisson_error);
  free(ewald->node_flags);
  /* CHANGE END - PoissonVerification */

  free(ewald);

  ewald_charge_on_ = 0;

  return 0;
}

/*****************************************************************************
 *
 *  ewald_charge_info
 *
 *****************************************************************************/

int ewald_charge_info(ewald_charge_t* ewald) {

  int ncolloid = 0;
  double eself;
  double sum_q2_coll = 0.0;
  double sum_q2_lattice = 0.0;

  assert(ewald);

  ewald_charge_self_energy(ewald, &eself);

  /* Compute sum of q^2 for colloids */
  if ((ewald->sources & EWALD_SOURCE_COLLOIDS) && ewald->cinfo) {
    colloid_t* pc = NULL;
    colloids_info_ntotal(ewald->cinfo, &ncolloid);
    colloids_info_local_head(ewald->cinfo, &pc);
    for (; pc; pc = pc->nextlocal) {
      double q = pc->s.q0 - pc->s.q1;
      sum_q2_coll += q * q;
    }
    MPI_Comm comm = MPI_COMM_NULL;
    cs_cart_comm(ewald->cs, &comm);
    MPI_Allreduce(MPI_IN_PLACE, &sum_q2_coll, 1, MPI_DOUBLE, MPI_SUM, comm);
  }

  /* Compute sum of q^2 for lattice nodes */
  if ((ewald->sources & EWALD_SOURCE_LATTICE) && ewald->psi) {
    int nlocal[3];
    cs_nlocal(ewald->cs, nlocal);

    for (int ic = 1; ic <= nlocal[X]; ic++) {
      for (int jc = 1; jc <= nlocal[Y]; jc++) {
        for (int kc = 1; kc <= nlocal[Z]; kc++) {
          int index = cs_index(ewald->cs, ic, jc, kc);

          /* Skip solid nodes if map is provided */
          if (ewald->map) {
            int status;
            map_status(ewald->map, index, &status);
            if (status != MAP_FLUID) continue;
          }

          double rho_elec;
          psi_rho_elec(ewald->psi, index, &rho_elec);
          /* In lattice units, q = rho_elec * dV = rho_elec * 1 */
          sum_q2_lattice += rho_elec * rho_elec;
        }
      }
    }
    MPI_Comm comm = MPI_COMM_NULL;
    cs_cart_comm(ewald->cs, &comm);
    MPI_Allreduce(MPI_IN_PLACE, &sum_q2_lattice, 1, MPI_DOUBLE, MPI_SUM, comm);
  }

  pe_info(ewald->pe, "\n");
  pe_info(ewald->pe, "Ewald sum for charges (Coulomb + Lattice)\n");
  pe_info(ewald->pe, "-----------------------------------------\n");
  pe_info(ewald->pe, "Sources included:                        ");
  if (ewald->sources & EWALD_SOURCE_COLLOIDS) pe_info(ewald->pe, " COLLOIDS");
  if (ewald->sources & EWALD_SOURCE_LATTICE) pe_info(ewald->pe, " LATTICE");
  pe_info(ewald->pe, "\n");
  pe_info(ewald->pe, "Number of colloids:                       %d\n", ncolloid);
  pe_info(ewald->pe, "Real space cut off:                      %14.7e\n", ewald_rc_);
  pe_info(ewald->pe, "Permittivity epsilon:                    %14.7e\n", epsilon_);
  pe_info(ewald->pe, "Ewald parameter alpha:                   %14.7e\n", alpha_);
  pe_info(ewald->pe, "Sum of q^2 (colloids):                   %14.7e\n", sum_q2_coll);
  pe_info(ewald->pe, "Sum of q^2 (lattice):                    %14.7e\n", sum_q2_lattice);
  pe_info(ewald->pe, "Self energy (constant):                  %14.7e\n", eself);
  pe_info(ewald->pe, "Maximum square wavevector:               %14.7e\n", kmax_);
  pe_info(ewald->pe, "Max nk limit (0=unlimited):               %d\n", nk_limit_);
  pe_info(ewald->pe, "Max. term retained in Fourier space sum:  %d\n", nkmax_);
  pe_info(ewald->pe, "Total terms kept in Fourier space sum:    %d\n\n", nktot_);

  return 0;
}

/*****************************************************************************
 *
 *  ewald_charge_check_electroneutrality
 *
 *  Compute the total charge of the system (colloids + fluid lattice nodes)
 *  and report it via pe_info. Call before the main loop to verify that
 *  initial conditions are charge-neutral.
 *
 *  Q_total = sum_p (q0_p - q1_p)  +  sum_r (rho0(r) - rho1(r))
 *
 *  Solid nodes are excluded from the lattice sum, consistent with the
 *  Ewald sum itself. Uses Kahan compensated summation.
 *
 *****************************************************************************/

int ewald_charge_check_electroneutrality(ewald_charge_t* ewald) {

  assert(ewald);

  int nlocal[3];
  cs_nlocal(ewald->cs, nlocal);

  MPI_Comm comm = MPI_COMM_NULL;
  cs_cart_comm(ewald->cs, &comm);

  kahan_t Q_colloids_k = kahan_zero();
  kahan_t Q_lattice_k = kahan_zero();

  /* Colloid contribution */
  if ((ewald->sources & EWALD_SOURCE_COLLOIDS) && ewald->cinfo) {
    colloid_t* pc = NULL;
    colloids_info_local_head(ewald->cinfo, &pc);
    for (; pc; pc = pc->nextlocal) {
      kahan_add_double(&Q_colloids_k, pc->s.q0 - pc->s.q1);
    }
  }

  /* Lattice (fluid nodes) contribution */
  if ((ewald->sources & EWALD_SOURCE_LATTICE) && ewald->psi) {
    for (int ic = 1; ic <= nlocal[X]; ic++) {
      for (int jc = 1; jc <= nlocal[Y]; jc++) {
        for (int kc = 1; kc <= nlocal[Z]; kc++) {
          int index = cs_index(ewald->cs, ic, jc, kc);

          if (ewald->map) {
            int status;
            map_status(ewald->map, index, &status);
            if (status != MAP_FLUID) continue;
          }

          double rho_elec;
          psi_rho_elec(ewald->psi, index, &rho_elec);
          kahan_add_double(&Q_lattice_k, rho_elec);
        }
      }
    }
  }

  /* MPI global reduce */
  {
    kahan_t buf[2] = { Q_colloids_k, Q_lattice_k };
    MPI_Datatype kahan_dt; MPI_Op kahan_op;
    kahan_mpi_datatype(&kahan_dt);
    kahan_mpi_op_sum(&kahan_op);
    MPI_Allreduce(MPI_IN_PLACE, buf, 2, kahan_dt, kahan_op, comm);
    MPI_Type_free(&kahan_dt);
    MPI_Op_free(&kahan_op);
    Q_colloids_k = buf[0];
    Q_lattice_k = buf[1];
  }

  double Q_colloids = kahan_sum(&Q_colloids_k);
  double Q_lattice = kahan_sum(&Q_lattice_k);
  double Q_total = Q_colloids + Q_lattice;

  pe_info(ewald->pe, "\n");
  pe_info(ewald->pe, "Charge neutrality check\n");
  pe_info(ewald->pe, "-----------------------\n");
  pe_info(ewald->pe, "  Q_colloids   = %+.15e\n", Q_colloids);
  pe_info(ewald->pe, "  Q_lattice    = %+.15e\n", Q_lattice);
  pe_info(ewald->pe, "  Q_total      = %+.15e\n", Q_total);
  if (fabs(Q_total) < 1.0e-10) {
    pe_info(ewald->pe, "  Status: NEUTRAL (|Q| < 1e-10)\n\n");
  }
  else {
    pe_info(ewald->pe, "  Status: NOT NEUTRAL -- charge imbalance detected!\n\n");
  }

  return 0;
}

/*****************************************************************************
 *
 *  ewald_charge_kappa
 *
 *  Return the value of the Ewald parameter alpha.
 *
 *****************************************************************************/

int ewald_charge_kappa(ewald_charge_t* ewald, double* kappa) {

  assert(ewald);
  assert(kappa);

  *kappa = alpha_;

  return 0;
}

/*****************************************************************************
 *
 *  ewald_charge_sum
 *
 *  Main routine: accumulate forces on all charged particles.
 *
 *  Includes:
 *    - Fourier space contribution (long-range)
 *    - Real space contribution (short-range)
 *    - External electric field contribution
 *
 *****************************************************************************/

int ewald_charge_sum(ewald_charge_t* ewald) {

  if (ewald == NULL) return 0;

  TIMER_start(TIMER_EWALD_TOTAL);

  ewald_charge_fourier_space_sum(ewald);
  ewald_charge_real_space_sum(ewald);
  ewald_charge_external_field(ewald);

  TIMER_stop(TIMER_EWALD_TOTAL);

  return 0;
}

/*****************************************************************************
 *
 *  ewald_charge_sum_full
 *
 *  Complete Ewald summation for potential and forces.
 *
 *  This function computes EVERYTHING using traditional Ewald summation:
 *
 *  1. POTENTIAL on all lattice nodes (stored in psi):
 *     - Fourier: sum over k-vectors with Gaussian screening
 *     - Real-space: particle->node + node->node within cutoff
 *
 *  2. FORCE on all lattice nodes (applied to hydro via hydro_f_local_add):
 *     - Fourier: sum over k-vectors
 *     - Real-space: particle->node + node->node within cutoff
 *
 *  3. FIELD (Esub) and FORCE (fex) on subgrid particles:
 *     - Fourier: sum over k-vectors
 *     - Real-space: node->particle + particle->particle within cutoff
 *
 *  Uses Kahan summation for numerical stability.
 *
 *****************************************************************************/

int ewald_charge_sum_full(ewald_charge_t* ewald, FILE* fp) {

  int nlocal[3], noffset[3], ncell[3];
  double ltot[3];
  double fkx, fky, fkz;
  double r4alpha_sq;
  double b0;
  int irc;
  PI_DOUBLE(pi);

  if (ewald == NULL) return 0;

  TIMER_start(TIMER_EWALD_TOTAL);

  cs_nlocal(ewald->cs, nlocal);
  cs_nlocal_offset(ewald->cs, noffset);
  cs_ltot(ewald->cs, ltot);
  if (ewald->cinfo) colloids_info_ncell(ewald->cinfo, ncell);

  irc = (int)ceil(ewald_rc_);

  fkx = 2.0 * pi / ltot[X];
  fky = 2.0 * pi / ltot[Y];
  fkz = 2.0 * pi / ltot[Z];
  r4alpha_sq = 1.0 / (4.0 * alpha_ * alpha_);

  /* Fourier prefactor: 1/(V*epsilon) */
  b0 = 1.0 / (ltot[X] * ltot[Y] * ltot[Z] * epsilon_);

  int ntotal_nodes = nlocal[X] * nlocal[Y] * nlocal[Z];

  /* Kahan accumulators for total forces (momentum conservation check) */
  kahan_t F_fluid_total[3] = { kahan_zero(), kahan_zero(), kahan_zero() };
  kahan_t F_particle_total[3] = { kahan_zero(), kahan_zero(), kahan_zero() };

  pe_info(ewald->pe, "\n");
  pe_info(ewald->pe, "Ewald charge sum full: computing potential, forces, and fields...\n");
  pe_info(ewald->pe, "  Local nodes: %d x %d x %d = %d\n", nlocal[X], nlocal[Y], nlocal[Z], ntotal_nodes);
  pe_info(ewald->pe, "  Fourier terms: %d\n", nktot_);
  pe_info(ewald->pe, "  Real-space cutoff: %.2f (irc=%d)\n", ewald_rc_, irc);

  /* ========================================================================
   * Precompute structure factors S(k) and C(k) using Kahan summation
   * ======================================================================== */

  pe_info(ewald->pe, "  [1/5] Computing structure factors S(k), C(k)...\n");

  /* Allocate Kahan arrays for structure factors */
  kahan_t* Sk_sin = (kahan_t*)malloc(nktot_ * sizeof(kahan_t));
  kahan_t* Sk_cos = (kahan_t*)malloc(nktot_ * sizeof(kahan_t));
  for (int n = 0; n < nktot_; n++) {
    Sk_sin[n] = kahan_zero();
    Sk_cos[n] = kahan_zero();
  }

  /* Precompute k-vectors and Green function */
  double* kvec = (double*)malloc(3 * nktot_ * sizeof(double));
  double* Gk_arr = (double*)malloc(nktot_ * sizeof(double));

  int kn = 0;
  for (int kz = 0; kz <= nk_[Z]; kz++) {
    for (int ky = -nk_[Y]; ky <= nk_[Y]; ky++) {
      for (int kx = -nk_[X]; kx <= nk_[X]; kx++) {
        double k[3], ksq;
        k[X] = fkx * kx;
        k[Y] = fky * ky;
        k[Z] = fkz * kz;
        ksq = k[X] * k[X] + k[Y] * k[Y] + k[Z] * k[Z];
        if (ksq <= 0.0 || ksq > kmax_) continue;

        kvec[3 * kn + X] = k[X];
        kvec[3 * kn + Y] = k[Y];
        kvec[3 * kn + Z] = k[Z];
        Gk_arr[kn] = b0 * exp(-r4alpha_sq * ksq) / ksq;
        kn++;
      }
    }
  }

  /* Structure factors from colloids */
  if ((ewald->sources & EWALD_SOURCE_COLLOIDS) && ewald->cinfo) {
    colloid_t* pc;
    colloids_info_local_head(ewald->cinfo, &pc);
    for (; pc; pc = pc->nextlocal) {
      double q = pc->s.q0 - pc->s.q1;
      ewald_charge_set_kr_table(ewald, pc->s.r);

      kn = 0;
      for (int kz = 0; kz <= nk_[Z]; kz++) {
        for (int ky = -nk_[Y]; ky <= nk_[Y]; ky++) {
          for (int kx = -nk_[X]; kx <= nk_[X]; kx++) {
            double ksq = (fkx * kx) * (fkx * kx) + (fky * ky) * (fky * ky) + (fkz * kz) * (fkz * kz);
            if (ksq <= 0.0 || ksq > kmax_) continue;

            double skr[3], ckr[3];
            skr[X] = sinkr_[3 * abs(kx) + X]; if (kx < 0) skr[X] = -skr[X];
            skr[Y] = sinkr_[3 * abs(ky) + Y]; if (ky < 0) skr[Y] = -skr[Y];
            skr[Z] = sinkr_[3 * kz + Z];
            ckr[X] = coskr_[3 * abs(kx) + X];
            ckr[Y] = coskr_[3 * abs(ky) + Y];
            ckr[Z] = coskr_[3 * kz + Z];

            double sinkr = skr[X] * ckr[Y] * ckr[Z] + ckr[X] * skr[Y] * ckr[Z]
              + ckr[X] * ckr[Y] * skr[Z] - skr[X] * skr[Y] * skr[Z];
            double coskr = ckr[X] * ckr[Y] * ckr[Z] - ckr[X] * skr[Y] * skr[Z]
              - skr[X] * ckr[Y] * skr[Z] - skr[X] * skr[Y] * ckr[Z];

            kahan_add_double(&Sk_sin[kn], q * sinkr);
            kahan_add_double(&Sk_cos[kn], q * coskr);
            kn++;
          }
        }
      }
    }
  }

  /* Structure factors from lattice nodes */
  if ((ewald->sources & EWALD_SOURCE_LATTICE) && ewald->psi) {
    for (int ic = 1; ic <= nlocal[X]; ic++) {
      for (int jc = 1; jc <= nlocal[Y]; jc++) {
        for (int kc = 1; kc <= nlocal[Z]; kc++) {
          int index = cs_index(ewald->cs, ic, jc, kc);

          if (ewald->map) {
            int status;
            map_status(ewald->map, index, &status);
            if (status != MAP_FLUID) continue;
          }

          double rho_elec;
          psi_rho_elec(ewald->psi, index, &rho_elec);
          double q = rho_elec;

          double r[3];
          r[X] = (double)(noffset[X] + ic);
          r[Y] = (double)(noffset[Y] + jc);
          r[Z] = (double)(noffset[Z] + kc);
          ewald_charge_set_kr_table(ewald, r);

          kn = 0;
          for (int kz = 0; kz <= nk_[Z]; kz++) {
            for (int ky = -nk_[Y]; ky <= nk_[Y]; ky++) {
              for (int kx = -nk_[X]; kx <= nk_[X]; kx++) {
                double ksq = (fkx * kx) * (fkx * kx) + (fky * ky) * (fky * ky) + (fkz * kz) * (fkz * kz);
                if (ksq <= 0.0 || ksq > kmax_) continue;

                double skr[3], ckr[3];
                skr[X] = sinkr_[3 * abs(kx) + X]; if (kx < 0) skr[X] = -skr[X];
                skr[Y] = sinkr_[3 * abs(ky) + Y]; if (ky < 0) skr[Y] = -skr[Y];
                skr[Z] = sinkr_[3 * kz + Z];
                ckr[X] = coskr_[3 * abs(kx) + X];
                ckr[Y] = coskr_[3 * abs(ky) + Y];
                ckr[Z] = coskr_[3 * kz + Z];

                double sinkr = skr[X] * ckr[Y] * ckr[Z] + ckr[X] * skr[Y] * ckr[Z]
                  + ckr[X] * ckr[Y] * skr[Z] - skr[X] * skr[Y] * skr[Z];
                double coskr = ckr[X] * ckr[Y] * ckr[Z] - ckr[X] * skr[Y] * skr[Z]
                  - skr[X] * ckr[Y] * skr[Z] - skr[X] * skr[Y] * ckr[Z];

                kahan_add_double(&Sk_sin[kn], q * sinkr);
                kahan_add_double(&Sk_cos[kn], q * coskr);
                kn++;
              }
            }
          }
        }
      }
    }
  }

  /* MPI reduce structure factors */
  pe_info(ewald->pe, "  [2/5] MPI reducing structure factors...\n");
  MPI_Comm comm;
  cs_cart_comm(ewald->cs, &comm);
  MPI_Datatype kahan_dt;
  MPI_Op kahan_op;
  kahan_mpi_datatype(&kahan_dt);
  kahan_mpi_op_sum(&kahan_op);
  MPI_Allreduce(MPI_IN_PLACE, Sk_sin, nktot_, kahan_dt, kahan_op, comm);
  MPI_Allreduce(MPI_IN_PLACE, Sk_cos, nktot_, kahan_dt, kahan_op, comm);
  MPI_Type_free(&kahan_dt);
  MPI_Op_free(&kahan_op);

  /* ========================================================================
   * Compute dipole moment M = sum_i q_i * r_i for dipole correction
   * The dipole correction field is: E_dipole = -(4*pi / (1+2*epsilon')*V) * M
   * ======================================================================== */

  kahan_t M_dipole[3] = { kahan_zero(), kahan_zero(), kahan_zero() };
  kahan_t Q_total_k = kahan_zero();  /* Total charge for diagnostics */

  /* Dipole from colloids */
  if ((ewald->sources & EWALD_SOURCE_COLLOIDS) && ewald->cinfo) {
    colloid_t* pc;
    colloids_info_local_head(ewald->cinfo, &pc);
    for (; pc; pc = pc->nextlocal) {
      double q = pc->s.q0 - pc->s.q1;
      kahan_add_double(&M_dipole[X], q * pc->s.r[X]);
      kahan_add_double(&M_dipole[Y], q * pc->s.r[Y]);
      kahan_add_double(&M_dipole[Z], q * pc->s.r[Z]);
      kahan_add_double(&Q_total_k, q);
    }
  }

  /* Dipole from lattice nodes */
  if ((ewald->sources & EWALD_SOURCE_LATTICE) && ewald->psi) {
    for (int ic = 1; ic <= nlocal[X]; ic++) {
      for (int jc = 1; jc <= nlocal[Y]; jc++) {
        for (int kc = 1; kc <= nlocal[Z]; kc++) {
          int index = cs_index(ewald->cs, ic, jc, kc);

          if (ewald->map) {
            int status;
            map_status(ewald->map, index, &status);
            if (status != MAP_FLUID) continue;
          }

          double rho_elec;
          psi_rho_elec(ewald->psi, index, &rho_elec);
          double q = rho_elec;

          double rx = (double)(noffset[X] + ic);
          double ry = (double)(noffset[Y] + jc);
          double rz = (double)(noffset[Z] + kc);

          kahan_add_double(&M_dipole[X], q * rx);
          kahan_add_double(&M_dipole[Y], q * ry);
          kahan_add_double(&M_dipole[Z], q * rz);
          kahan_add_double(&Q_total_k, q);
        }
      }
    }
  }

  /* MPI reduce dipole moment and total charge */
  {
    kahan_t M_reduce[4] = { M_dipole[X], M_dipole[Y], M_dipole[Z], Q_total_k };
    kahan_mpi_datatype(&kahan_dt);
    kahan_mpi_op_sum(&kahan_op);
    MPI_Allreduce(MPI_IN_PLACE, M_reduce, 4, kahan_dt, kahan_op, comm);
    MPI_Type_free(&kahan_dt);
    MPI_Op_free(&kahan_op);
    M_dipole[X] = M_reduce[0];
    M_dipole[Y] = M_reduce[1];
    M_dipole[Z] = M_reduce[2];
    Q_total_k = M_reduce[3];
  }

  double M[3] = { kahan_sum(&M_dipole[X]), kahan_sum(&M_dipole[Y]), kahan_sum(&M_dipole[Z]) };
  double Q_total = kahan_sum(&Q_total_k);
  double V = ltot[X] * ltot[Y] * ltot[Z];

  /* Dipole correction using Deserno & Holm convention (J. Chem. Phys. 1998)
   * Energy:    E^(d) = 2*pi / ((1 + 2*epsilon')*V) * |M|^2
   * Force:     F^(d)_i = -4*pi*q_i / ((1 + 2*epsilon')*V) * M
   * Field:     E_dipole = -4*pi / ((1 + 2*epsilon')*V) * M
   *
   * epsilon' = 1     : vacuum boundary conditions
   * epsilon' = infty : metallic (tinfoil) boundary conditions (correction vanishes)
   */
  double dipole_prefactor = 4.0 * pi / ((1.0 + 2.0 * ewald->epsilon_prime) * V);
  double E_dipole[3] = { -dipole_prefactor * M[X],
                        -dipole_prefactor * M[Y],
                        -dipole_prefactor * M[Z] };

  pe_info(ewald->pe, "  Dipole moment M = (%1.15e, %1.15e, %1.15e)\n", M[X], M[Y], M[Z]);
  pe_info(ewald->pe, "  Total charge Q = %1.15e\n", Q_total);
  pe_info(ewald->pe, "  Dipole epsilon' = %14.7e (Deserno & Holm convention)\n", ewald->epsilon_prime);
  pe_info(ewald->pe, "  Dipole field E_dip = (%1.15e, %1.15e, %1.15e)\n", E_dipole[X], E_dipole[Y], E_dipole[Z]);

  /* ========================================================================
   * PART 1: Potential on lattice nodes
   * ======================================================================== */

  pe_info(ewald->pe, "  [3/5] Computing potential on lattice nodes...\n");
  if (ewald->psi) {
    for (int ic = 1; ic <= nlocal[X]; ic++) {
      for (int jc = 1; jc <= nlocal[Y]; jc++) {
        for (int kc = 1; kc <= nlocal[Z]; kc++) {
          int index = cs_index(ewald->cs, ic, jc, kc);

          if (ewald->map) {
            int status;
            map_status(ewald->map, index, &status);
            if (status != MAP_FLUID) continue;
          }

          double r[3];
          r[X] = (double)(noffset[X] + ic);
          r[Y] = (double)(noffset[Y] + jc);
          r[Z] = (double)(noffset[Z] + kc);

          /* Fourier part (Kahan) */
          kahan_t phi_fourier_k = kahan_zero();
          ewald_charge_set_kr_table(ewald, r);

          kn = 0;
          for (int kz = 0; kz <= nk_[Z]; kz++) {
            for (int ky = -nk_[Y]; ky <= nk_[Y]; ky++) {
              for (int kx = -nk_[X]; kx <= nk_[X]; kx++) {
                double ksq = (fkx * kx) * (fkx * kx) + (fky * ky) * (fky * ky) + (fkz * kz) * (fkz * kz);
                if (ksq <= 0.0 || ksq > kmax_) continue;

                double skr[3], ckr[3];
                skr[X] = sinkr_[3 * abs(kx) + X]; if (kx < 0) skr[X] = -skr[X];
                skr[Y] = sinkr_[3 * abs(ky) + Y]; if (ky < 0) skr[Y] = -skr[Y];
                skr[Z] = sinkr_[3 * kz + Z];
                ckr[X] = coskr_[3 * abs(kx) + X];
                ckr[Y] = coskr_[3 * abs(ky) + Y];
                ckr[Z] = coskr_[3 * kz + Z];

                double sinkr = skr[X] * ckr[Y] * ckr[Z] + ckr[X] * skr[Y] * ckr[Z]
                  + ckr[X] * ckr[Y] * skr[Z] - skr[X] * skr[Y] * skr[Z];
                double coskr = ckr[X] * ckr[Y] * ckr[Z] - ckr[X] * skr[Y] * skr[Z]
                  - skr[X] * ckr[Y] * skr[Z] - skr[X] * skr[Y] * ckr[Z];

                double factor = (kz > 0) ? 2.0 : 1.0;
                double Sk_sin_val = kahan_sum(&Sk_sin[kn]);
                double Sk_cos_val = kahan_sum(&Sk_cos[kn]);
                kahan_add_double(&phi_fourier_k, factor * Gk_arr[kn] * (Sk_cos_val * coskr + Sk_sin_val * sinkr));
                kn++;
              }
            }
          }
          double phi_fourier = kahan_sum(&phi_fourier_k);

          /* Real-space from particles */
          double phi_real = 0.0;
          double r_local[3] = { (double)ic, (double)jc, (double)kc };

          if ((ewald->sources & EWALD_SOURCE_COLLOIDS) && ewald->cinfo) {
            for (int icc = 0; icc <= ncell[X] + 1; icc++) {
              for (int jcc = 0; jcc <= ncell[Y] + 1; jcc++) {
                for (int kcc = 0; kcc <= ncell[Z] + 1; kcc++) {
                  colloid_t* pc;
                  colloids_info_cell_list_head(ewald->cinfo, icc, jcc, kcc, &pc);
                  for (; pc; pc = pc->next) {
                    double q_p = pc->s.q0 - pc->s.q1;
                    double r_p[3];
                    r_p[X] = pc->s.r[X] - (double)noffset[X];
                    r_p[Y] = pc->s.r[Y] - (double)noffset[Y];
                    r_p[Z] = pc->s.r[Z] - (double)noffset[Z];
                    double dr[3] = { r_local[X] - r_p[X], r_local[Y] - r_p[Y], r_local[Z] - r_p[Z] };
                    double dist = sqrt(dr[X] * dr[X] + dr[Y] * dr[Y] + dr[Z] * dr[Z]);
                    if (dist < ewald_rc_ && dist > 1.0e-10) {
                      phi_real += q_p * erfc(alpha_ * dist) / (4.0 * pi * epsilon_ * dist);
                    }
                  }
                }
              }
            }
          }

          /* Real-space from other nodes */
          if (ewald->sources & EWALD_SOURCE_LATTICE) {
            for (int di = -irc; di <= irc; di++) {
              for (int dj = -irc; dj <= irc; dj++) {
                for (int dk = -irc; dk <= irc; dk++) {
                  if (di == 0 && dj == 0 && dk == 0) continue;
                  int i2 = ic + di, j2 = jc + dj, k2 = kc + dk;
                  if (i2 < 1 || i2 > nlocal[X]) continue;
                  if (j2 < 1 || j2 > nlocal[Y]) continue;
                  if (k2 < 1 || k2 > nlocal[Z]) continue;
                  int index2 = cs_index(ewald->cs, i2, j2, k2);
                  double rho_elec;
                  psi_rho_elec(ewald->psi, index2, &rho_elec);
                  double q2 = rho_elec;
                  if (fabs(q2) < 1.0e-14) continue;
                  double dist = sqrt((double)(di * di + dj * dj + dk * dk));
                  if (dist < ewald_rc_) {
                    phi_real += q2 * erfc(alpha_ * dist) / (4.0 * pi * epsilon_ * dist);
                  }
                }
              }
            }
          }

          /* Dipole correction to potential: phi_dipole = (4*pi / (1+2*epsilon')*V) * M . r
           * This is equivalent to -E_dipole . r */
          double phi_dipole = dipole_prefactor * (M[X] * r[X] + M[Y] * r[Y] + M[Z] * r[Z]);
          //  double phi_dipole = 0.0;  /* For now, we can set this to zero and only apply the dipole correction to forces and fields */

          /* Store: psi* = beta * eunit * phi */
          ewald->psi->psi->data[addr_rank0(ewald->psi->psi->nsites, index)] =
            beta_ * eunit_ * (phi_fourier + phi_real + phi_dipole);
        }
      }
    }
  }

  /* ========================================================================
   * PART 2: Force on lattice nodes -> hydro
   * ======================================================================== */

  pe_info(ewald->pe, "  [4/5] Computing force on lattice nodes...\n");
  if (ewald->hydro != NULL) {
    // /* Sync hydro from device before modification */
    // hydro_memcpy(ewald->hydro, tdpMemcpyDeviceToHost);

    for (int ic = 1; ic <= nlocal[X]; ic++) {
      for (int jc = 1; jc <= nlocal[Y]; jc++) {
        for (int kc = 1; kc <= nlocal[Z]; kc++) {
          int index = cs_index(ewald->cs, ic, jc, kc);

          if (ewald->map) {
            int status;
            map_status(ewald->map, index, &status);
            if (status != MAP_FLUID) continue;
          }

          double rho_elec;
          psi_rho_elec(ewald->psi, index, &rho_elec);
          double q1 = rho_elec;

          double r[3];
          r[X] = (double)(noffset[X] + ic);
          r[Y] = (double)(noffset[Y] + jc);
          r[Z] = (double)(noffset[Z] + kc);

          /* Fourier part (Kahan) */
          kahan_t F_fourier_k[3] = { kahan_zero(), kahan_zero(), kahan_zero() };
          ewald_charge_set_kr_table(ewald, r);

          kn = 0;
          for (int kz = 0; kz <= nk_[Z]; kz++) {
            for (int ky = -nk_[Y]; ky <= nk_[Y]; ky++) {
              for (int kx = -nk_[X]; kx <= nk_[X]; kx++) {
                double ksq = (fkx * kx) * (fkx * kx) + (fky * ky) * (fky * ky) + (fkz * kz) * (fkz * kz);
                if (ksq <= 0.0 || ksq > kmax_) continue;

                double skr[3], ckr[3];
                skr[X] = sinkr_[3 * abs(kx) + X]; if (kx < 0) skr[X] = -skr[X];
                skr[Y] = sinkr_[3 * abs(ky) + Y]; if (ky < 0) skr[Y] = -skr[Y];
                skr[Z] = sinkr_[3 * kz + Z];
                ckr[X] = coskr_[3 * abs(kx) + X];
                ckr[Y] = coskr_[3 * abs(ky) + Y];
                ckr[Z] = coskr_[3 * kz + Z];

                double sinkr = skr[X] * ckr[Y] * ckr[Z] + ckr[X] * skr[Y] * ckr[Z]
                  + ckr[X] * ckr[Y] * skr[Z] - skr[X] * skr[Y] * skr[Z];
                double coskr = ckr[X] * ckr[Y] * ckr[Z] - ckr[X] * skr[Y] * skr[Z]
                  - skr[X] * ckr[Y] * skr[Z] - skr[X] * skr[Y] * ckr[Z];

                double factor = (kz > 0) ? 2.0 : 1.0;
                double Sk_sin_val = kahan_sum(&Sk_sin[kn]);
                double Sk_cos_val = kahan_sum(&Sk_cos[kn]);
                double im_part = Sk_sin_val * coskr - Sk_cos_val * sinkr;

                kahan_add_double(&F_fourier_k[X], factor * q1 * Gk_arr[kn] * kvec[3 * kn + X] * im_part);
                kahan_add_double(&F_fourier_k[Y], factor * q1 * Gk_arr[kn] * kvec[3 * kn + Y] * im_part);
                kahan_add_double(&F_fourier_k[Z], factor * q1 * Gk_arr[kn] * kvec[3 * kn + Z] * im_part);
                kn++;
              }
            }
          }
          double F_total[3] = { kahan_sum(&F_fourier_k[X]), kahan_sum(&F_fourier_k[Y]), kahan_sum(&F_fourier_k[Z]) };

          /* Real-space from particles */
          double r_local[3] = { (double)ic, (double)jc, (double)kc };
          if ((ewald->sources & EWALD_SOURCE_COLLOIDS) && ewald->cinfo) {
            for (int icc = 0; icc <= ncell[X] + 1; icc++) {
              for (int jcc = 0; jcc <= ncell[Y] + 1; jcc++) {
                for (int kcc = 0; kcc <= ncell[Z] + 1; kcc++) {
                  colloid_t* pc;
                  colloids_info_cell_list_head(ewald->cinfo, icc, jcc, kcc, &pc);
                  for (; pc; pc = pc->next) {
                    double q_p = pc->s.q0 - pc->s.q1;
                    double r_p[3];
                    r_p[X] = pc->s.r[X] - (double)noffset[X];
                    r_p[Y] = pc->s.r[Y] - (double)noffset[Y];
                    r_p[Z] = pc->s.r[Z] - (double)noffset[Z];
                    double dr[3] = { r_local[X] - r_p[X], r_local[Y] - r_p[Y], r_local[Z] - r_p[Z] };
                    double dist = sqrt(dr[X] * dr[X] + dr[Y] * dr[Y] + dr[Z] * dr[Z]);
                    if (dist < ewald_rc_ && dist > 1.0e-10) {
                      double ar = alpha_ * dist;
                      double r_inv = 1.0 / dist;
                      double term1 = erfc(ar) / (dist * dist);
                      double term2 = (2.0 * alpha_ / sqrt(pi)) * exp(-ar * ar) * r_inv;
                      double F_mag = q1 * q_p * (term1 + term2) / (4.0 * pi * epsilon_);
                      F_total[X] += F_mag * dr[X] * r_inv;
                      F_total[Y] += F_mag * dr[Y] * r_inv;
                      F_total[Z] += F_mag * dr[Z] * r_inv;
                    }
                  }
                }
              }
            }
          }

          /* Real-space from other nodes */
          if (ewald->sources & EWALD_SOURCE_LATTICE) {
            for (int di = -irc; di <= irc; di++) {
              for (int dj = -irc; dj <= irc; dj++) {
                for (int dk = -irc; dk <= irc; dk++) {
                  if (di == 0 && dj == 0 && dk == 0) continue;
                  int i2 = ic + di, j2 = jc + dj, k2 = kc + dk;
                  if (i2 < 1 || i2 > nlocal[X]) continue;
                  if (j2 < 1 || j2 > nlocal[Y]) continue;
                  if (k2 < 1 || k2 > nlocal[Z]) continue;
                  int index2 = cs_index(ewald->cs, i2, j2, k2);
                  double rho_elec;
                  psi_rho_elec(ewald->psi, index2, &rho_elec);
                  double q2 = rho_elec;
                  if (fabs(q2) < 1.0e-14) continue;
                  double dist = sqrt((double)(di * di + dj * dj + dk * dk));
                  if (dist < ewald_rc_) {
                    double ar = alpha_ * dist;
                    double r_inv = 1.0 / dist;
                    double term1 = erfc(ar) / (dist * dist);
                    double term2 = (2.0 * alpha_ / sqrt(pi)) * exp(-ar * ar) * r_inv;
                    double F_mag = q1 * q2 * (term1 + term2) / (4.0 * pi * epsilon_);
                    F_total[X] += F_mag * (double)(-di) * r_inv;
                    F_total[Y] += F_mag * (double)(-dj) * r_inv;
                    F_total[Z] += F_mag * (double)(-dk) * r_inv;
                  }
                }
              }
            }
          }

          /* Dipole correction to force: F_dipole = q * E_dipole */
          F_total[X] += q1 * E_dipole[X];
          F_total[Y] += q1 * E_dipole[Y];
          F_total[Z] += q1 * E_dipole[Z];

          hydro_f_local_add(ewald->hydro, index, F_total);

          /* Accumulate total force on fluid (for momentum conservation check) */
          kahan_add_double(&F_fluid_total[X], F_total[X]);
          kahan_add_double(&F_fluid_total[Y], F_total[Y]);
          kahan_add_double(&F_fluid_total[Z], F_total[Z]);
        }
      }
    }

    // /* Sync hydro back to device */
    // hydro_memcpy(ewald->hydro, tdpMemcpyHostToDevice);
  }

  /* ========================================================================
   * PART 3: Field (Esub) and Force (fex) on subgrid particles
   * ======================================================================== */

  pe_info(ewald->pe, "  [5/5] Computing field and force on particles...\n");
  if ((ewald->sources & EWALD_SOURCE_COLLOIDS) && ewald->cinfo) {
    for (int icc = 0; icc <= ncell[X] + 1; icc++) {
      for (int jcc = 0; jcc <= ncell[Y] + 1; jcc++) {
        for (int kcc = 0; kcc <= ncell[Z] + 1; kcc++) {
          colloid_t* pc;
          colloids_info_cell_list_head(ewald->cinfo, icc, jcc, kcc, &pc);
          for (; pc; pc = pc->next) {
            double q_p = pc->s.q0 - pc->s.q1;

            /* Fourier part (Kahan) */
            kahan_t E_fourier_k[3] = { kahan_zero(), kahan_zero(), kahan_zero() };
            ewald_charge_set_kr_table(ewald, pc->s.r);

            kn = 0;
            for (int kz = 0; kz <= nk_[Z]; kz++) {
              for (int ky = -nk_[Y]; ky <= nk_[Y]; ky++) {
                for (int kx = -nk_[X]; kx <= nk_[X]; kx++) {
                  double ksq = (fkx * kx) * (fkx * kx) + (fky * ky) * (fky * ky) + (fkz * kz) * (fkz * kz);
                  if (ksq <= 0.0 || ksq > kmax_) continue;

                  double skr[3], ckr[3];
                  skr[X] = sinkr_[3 * abs(kx) + X]; if (kx < 0) skr[X] = -skr[X];
                  skr[Y] = sinkr_[3 * abs(ky) + Y]; if (ky < 0) skr[Y] = -skr[Y];
                  skr[Z] = sinkr_[3 * kz + Z];
                  ckr[X] = coskr_[3 * abs(kx) + X];
                  ckr[Y] = coskr_[3 * abs(ky) + Y];
                  ckr[Z] = coskr_[3 * kz + Z];

                  double sinkr = skr[X] * ckr[Y] * ckr[Z] + ckr[X] * skr[Y] * ckr[Z]
                    + ckr[X] * ckr[Y] * skr[Z] - skr[X] * skr[Y] * skr[Z];
                  double coskr = ckr[X] * ckr[Y] * ckr[Z] - ckr[X] * skr[Y] * skr[Z]
                    - skr[X] * ckr[Y] * skr[Z] - skr[X] * skr[Y] * ckr[Z];

                  double factor = (kz > 0) ? 2.0 : 1.0;
                  double Sk_sin_val = kahan_sum(&Sk_sin[kn]);
                  double Sk_cos_val = kahan_sum(&Sk_cos[kn]);
                  double im_part = Sk_sin_val * coskr - Sk_cos_val * sinkr;

                  /* Esub = beta * eunit * E_physical */
                  kahan_add_double(&E_fourier_k[X], factor * beta_ * eunit_ * Gk_arr[kn] * kvec[3 * kn + X] * im_part);
                  kahan_add_double(&E_fourier_k[Y], factor * beta_ * eunit_ * Gk_arr[kn] * kvec[3 * kn + Y] * im_part);
                  kahan_add_double(&E_fourier_k[Z], factor * beta_ * eunit_ * Gk_arr[kn] * kvec[3 * kn + Z] * im_part);
                  kn++;
                }
              }
            }
            double E_total[3] = { kahan_sum(&E_fourier_k[X]), kahan_sum(&E_fourier_k[Y]), kahan_sum(&E_fourier_k[Z]) };

            /* Real-space from lattice nodes */
            double r_p[3];
            r_p[X] = pc->s.r[X] - (double)noffset[X];
            r_p[Y] = pc->s.r[Y] - (double)noffset[Y];
            r_p[Z] = pc->s.r[Z] - (double)noffset[Z];
            int i0 = (int)floor(r_p[X]);
            int j0 = (int)floor(r_p[Y]);
            int k0 = (int)floor(r_p[Z]);

            if (ewald->sources & EWALD_SOURCE_LATTICE) {
              for (int di = -irc; di <= irc + 1; di++) {
                for (int dj = -irc; dj <= irc + 1; dj++) {
                  for (int dk = -irc; dk <= irc + 1; dk++) {
                    int ni = i0 + di, nj = j0 + dj, nk_idx = k0 + dk;
                    if (ni < 1 || ni > nlocal[X]) continue;
                    if (nj < 1 || nj > nlocal[Y]) continue;
                    if (nk_idx < 1 || nk_idx > nlocal[Z]) continue;
                    int index = cs_index(ewald->cs, ni, nj, nk_idx);
                    double rho_elec;
                    psi_rho_elec(ewald->psi, index, &rho_elec);
                    double q_node = rho_elec;
                    double dr[3] = { r_p[X] - (double)ni, r_p[Y] - (double)nj, r_p[Z] - (double)nk_idx };
                    double dist = sqrt(dr[X] * dr[X] + dr[Y] * dr[Y] + dr[Z] * dr[Z]);
                    if (dist < ewald_rc_ && dist > 1.0e-10) {
                      double ar = alpha_ * dist;
                      double r_inv = 1.0 / dist;
                      double term1 = erfc(ar) / (dist * dist);
                      double term2 = (2.0 * alpha_ / sqrt(pi)) * exp(-ar * ar) * r_inv;
                      double E_mag = beta_ * eunit_ * (term1 + term2) / (4.0 * pi * epsilon_);
                      E_total[X] += q_node * E_mag * dr[X] * r_inv;
                      E_total[Y] += q_node * E_mag * dr[Y] * r_inv;
                      E_total[Z] += q_node * E_mag * dr[Z] * r_inv;
                    }
                  }
                }
              }
            }

            /* Real-space from other particles */
            for (int ic2 = 0; ic2 <= ncell[X] + 1; ic2++) {
              for (int jc2 = 0; jc2 <= ncell[Y] + 1; jc2++) {
                for (int kc2 = 0; kc2 <= ncell[Z] + 1; kc2++) {
                  colloid_t* pc2;
                  colloids_info_cell_list_head(ewald->cinfo, ic2, jc2, kc2, &pc2);
                  for (; pc2; pc2 = pc2->next) {
                    if (pc2 == pc) continue;
                    double q_p2 = pc2->s.q0 - pc2->s.q1;
                    double r_p2[3];
                    r_p2[X] = pc2->s.r[X] - (double)noffset[X];
                    r_p2[Y] = pc2->s.r[Y] - (double)noffset[Y];
                    r_p2[Z] = pc2->s.r[Z] - (double)noffset[Z];
                    double dr[3] = { r_p[X] - r_p2[X], r_p[Y] - r_p2[Y], r_p[Z] - r_p2[Z] };
                    double dist = sqrt(dr[X] * dr[X] + dr[Y] * dr[Y] + dr[Z] * dr[Z]);
                    if (dist < ewald_rc_ && dist > 1.0e-10) {
                      double ar = alpha_ * dist;
                      double r_inv = 1.0 / dist;
                      double term1 = erfc(ar) / (dist * dist);
                      double term2 = (2.0 * alpha_ / sqrt(pi)) * exp(-ar * ar) * r_inv;
                      double E_mag = beta_ * eunit_ * (term1 + term2) / (4.0 * pi * epsilon_);
                      E_total[X] += q_p2 * E_mag * dr[X] * r_inv;
                      E_total[Y] += q_p2 * E_mag * dr[Y] * r_inv;
                      E_total[Z] += q_p2 * E_mag * dr[Z] * r_inv;
                    }
                  }
                }
              }
            }

            /* Add dipole correction to field: Esub_dipole = beta * eunit * E_dipole */
            E_total[X] += beta_ * eunit_ * E_dipole[X];
            E_total[Y] += beta_ * eunit_ * E_dipole[Y];
            E_total[Z] += beta_ * eunit_ * E_dipole[Z];

            /* Store field and force */
            pc->Esub[X] = E_total[X];
            pc->Esub[Y] = E_total[Y];
            pc->Esub[Z] = E_total[Z];

            double kt = 1.0 / beta_;
            pc->fex[X] = q_p * E_total[X] * kt / eunit_;
            pc->fex[Y] = q_p * E_total[Y] * kt / eunit_;
            pc->fex[Z] = q_p * E_total[Z] * kt / eunit_;

            /* Accumulate total force on particles (for momentum conservation check) */
            kahan_add_double(&F_particle_total[X], pc->fex[X]);
            kahan_add_double(&F_particle_total[Y], pc->fex[Y]);
            kahan_add_double(&F_particle_total[Z], pc->fex[Z]);
          }
        }
      }
    }
  }

  /* External field contribution */
  ewald_charge_external_field(ewald);

  /* Cleanup */
  free(Sk_sin);
  free(Sk_cos);
  free(kvec);
  free(Gk_arr);

  /* MPI reduce total forces for momentum conservation check */
  {
    MPI_Comm comm;
    cs_cart_comm(ewald->cs, &comm);
    MPI_Datatype kahan_dt;
    MPI_Op kahan_op;
    kahan_mpi_datatype(&kahan_dt);
    kahan_mpi_op_sum(&kahan_op);
    MPI_Allreduce(MPI_IN_PLACE, F_fluid_total, 3, kahan_dt, kahan_op, comm);
    MPI_Allreduce(MPI_IN_PLACE, F_particle_total, 3, kahan_dt, kahan_op, comm);
    MPI_Type_free(&kahan_dt);
    MPI_Op_free(&kahan_op);
  }

  /* Print momentum conservation check */
  double F_fluid[3] = { kahan_sum(&F_fluid_total[X]), kahan_sum(&F_fluid_total[Y]), kahan_sum(&F_fluid_total[Z]) };
  double F_particle[3] = { kahan_sum(&F_particle_total[X]), kahan_sum(&F_particle_total[Y]), kahan_sum(&F_particle_total[Z]) };
  double F_diff[3] = { F_fluid[X] + F_particle[X], F_fluid[Y] + F_particle[Y], F_fluid[Z] + F_particle[Z] };

  fprintf(fp, "%1.15e,%1.15e, %1.15e, %1.15e,%1.15e, %1.15e, %1.15e,%1.15e, %1.15e, %1.15e,\n",
                sqrt(F_diff[X] * F_diff[X] + F_diff[Y] * F_diff[Y] + F_diff[Z] * F_diff[Z]),
                F_diff[X], F_diff[Y], F_diff[Z],
                F_fluid[X], F_fluid[Y], F_fluid[Z],
                F_particle[X], F_particle[Y], F_particle[Z]);

  pe_info(ewald->pe, "\n");
  pe_info(ewald->pe, "  Momentum conservation check:\n");
  pe_info(ewald->pe, "    F_fluid    = (%14.7e, %14.7e, %14.7e)\n", F_fluid[X], F_fluid[Y], F_fluid[Z]);
  pe_info(ewald->pe, "    F_particle = (%14.7e, %14.7e, %14.7e)\n", F_particle[X], F_particle[Y], F_particle[Z]);
  pe_info(ewald->pe, "    F_total    = (%14.7e, %14.7e, %14.7e)\n", F_diff[X], F_diff[Y], F_diff[Z]);
  pe_info(ewald->pe, "    |F_total|  = %14.7e\n", sqrt(F_diff[X] * F_diff[X] + F_diff[Y] * F_diff[Y] + F_diff[Z] * F_diff[Z]));

  pe_info(ewald->pe, "Ewald charge sum full: complete.\n\n");

  TIMER_stop(TIMER_EWALD_TOTAL);

  return 0;
}

/*****************************************************************************
 *
 *  ewald_charge_external_field
 *
 *  Add force from external electric field E0 on each charged colloid.
 *
 *  F_ext = q * E0
 *
 *****************************************************************************/

int ewald_charge_external_field(ewald_charge_t* ewald) {

  double e0[3];
  int ncell[3];

  if (ewald == NULL) return 0;
  if (ewald->psi == NULL) return 0;
  if (!(ewald->sources & EWALD_SOURCE_COLLOIDS)) return 0;
  if (ewald->cinfo == NULL) return 0;

  /* Get external field from psi structure */
  e0[X] = ewald->psi->e0[X];
  e0[Y] = ewald->psi->e0[Y];
  e0[Z] = ewald->psi->e0[Z];

  /* Skip if no external field */
  if (fabs(e0[X]) < 1.0e-14 && fabs(e0[Y]) < 1.0e-14 && fabs(e0[Z]) < 1.0e-14) {
    return 0;
  }

  colloids_info_ncell(ewald->cinfo, ncell);

  for (int ic = 1; ic <= ncell[X]; ic++) {
    for (int jc = 1; jc <= ncell[Y]; jc++) {
      for (int kc = 1; kc <= ncell[Z]; kc++) {

        colloid_t* pc;
        colloids_info_cell_list_head(ewald->cinfo, ic, jc, kc, &pc);

        for (; pc; pc = pc->next) {

          double q = pc->s.q0 - pc->s.q1;
          if (fabs(q) < 1.0e-12) continue;

          /* F = q * E0 */
          // pc->force[X] += q * e0[X];
          // pc->force[Y] += q * e0[Y];
          // pc->force[Z] += q * e0[Z];
          pc->Esub[X] += e0[X]; // Transformo en campo electrico
          pc->Esub[Y] += e0[Y]; // Transformo en campo electrico
          pc->Esub[Z] += e0[Z]; // Transformo en campo electrico
        }
      }
    }
  }

  return 0;
}

/*****************************************************************************
 *
 *  ewald_charge_real_space_energy
 *
 *  Compute contribution to the energy from a single pair interaction.
 *
 *  E_real = q1 * q2 * erfc(alpha*r) / (4*pi*epsilon*r)
 *
 *****************************************************************************/

int ewald_charge_real_space_energy(ewald_charge_t* ewald, double q1, double q2,
                                   const double r12[3], double* ereal) {
  double r;
  PI_DOUBLE(pi);

  assert(ewald);
  assert(ereal);

  *ereal = 0.0;

  r = sqrt(r12[X] * r12[X] + r12[Y] * r12[Y] + r12[Z] * r12[Z]);

  if (r < ewald_rc_ && r > 0.0) {
    /* Coulomb energy with erfc screening */
    *ereal = q1 * q2 * erfc(alpha_ * r) / (4.0 * pi * epsilon_ * r);
  }

  return 0;
}

/*****************************************************************************
 *
 *  ewald_charge_fourier_space_energy
 *
 *  Fourier-space part of the Ewald summation for the energy.
 *
 *  E_fourier = (1/2V*epsilon) * sum_{k!=0} (4*pi/k^2) * exp(-k^2/4*alpha^2)
 *              * [S(k)^2 + C(k)^2]
 *
 *  where S(k) = sum_i q_i * sin(k.r_i)  (includes both colloids and lattice)
 *        C(k) = sum_i q_i * cos(k.r_i)
 *
 *****************************************************************************/

int ewald_charge_fourier_space_energy(ewald_charge_t* ewald, double* ef) {

  double e = 0.0;
  double k[3], ksq;
  double fkx, fky, fkz;
  double b0, b;
  double r4alpha_sq;
  double ltot[3];
  int kx, ky, kz, kn = 0;
  PI_DOUBLE(pi);

  assert(ewald);

  cs_ltot(ewald->cs, ltot);
  ewald_charge_sum_sin_cos_terms(ewald);

  fkx = 2.0 * pi / ltot[X];
  fky = 2.0 * pi / ltot[Y];
  fkz = 2.0 * pi / ltot[Z];

  /* Prefactor: (4*pi) / (2 * V * epsilon) = 2*pi / (V * epsilon) */
  b0 = 2.0 * pi / (ltot[X] * ltot[Y] * ltot[Z] * epsilon_);
  r4alpha_sq = 1.0 / (4.0 * alpha_ * alpha_);

  /* Sum over k to get the energy. */

  for (kz = 0; kz <= nk_[Z]; kz++) {
    for (ky = -nk_[Y]; ky <= nk_[Y]; ky++) {
      for (kx = -nk_[X]; kx <= nk_[X]; kx++) {

        k[X] = fkx * kx;
        k[Y] = fky * ky;
        k[Z] = fkz * kz;
        ksq = k[X] * k[X] + k[Y] * k[Y] + k[Z] * k[Z];

        if (ksq <= 0.0 || ksq > kmax_) continue;

        b = b0 * exp(-r4alpha_sq * ksq) / ksq;

        if (kz == 0) {
          e += 0.5 * b * (sinx_[kn] * sinx_[kn] + cosx_[kn] * cosx_[kn]);
        }
        else {
          e += b * (sinx_[kn] * sinx_[kn] + cosx_[kn] * cosx_[kn]);
        }
        kn++;
      }
    }
  }

  *ef = e;

  return 0;
}

/*****************************************************************************
 *
 *  ewald_charge_sum_sin_cos_terms
 *
 *  Compute S(k) and C(k) from all charge sources (colloids + lattice).
 *
 *****************************************************************************/

static int ewald_charge_sum_sin_cos_terms(ewald_charge_t* ewald) {

  assert(ewald);

  /* Initialize S(k) and C(k) to zero */
  for (int kn = 0; kn < nktot_; kn++) {
    sinx_[kn] = 0.0;
    cosx_[kn] = 0.0;
  }

  /* Add colloid contributions */
  if (ewald->sources & EWALD_SOURCE_COLLOIDS) {
    ewald_charge_sum_sin_cos_colloids(ewald);
  }

  /* Add lattice contributions */
  if (ewald->sources & EWALD_SOURCE_LATTICE) {
    ewald_charge_sum_sin_cos_lattice(ewald);
  }

  /* MPI reduction to get global S(k) and C(k) */
  {
    double* subsin = NULL;
    double* subcos = NULL;
    MPI_Comm comm = MPI_COMM_NULL;

    cs_cart_comm(ewald->cs, &comm);

    subsin = (double*)calloc(nktot_, sizeof(double));
    subcos = (double*)calloc(nktot_, sizeof(double));
    assert(subsin);
    assert(subcos);
    if (subsin == NULL) pe_fatal(ewald->pe, "calloc(subsin) failed\n");
    if (subcos == NULL) pe_fatal(ewald->pe, "calloc(subcos) failed\n");

    for (int kn = 0; kn < nktot_; kn++) {
      subsin[kn] = sinx_[kn];
      subcos[kn] = cosx_[kn];
    }

    MPI_Allreduce(subsin, sinx_, nktot_, MPI_DOUBLE, MPI_SUM, comm);
    MPI_Allreduce(subcos, cosx_, nktot_, MPI_DOUBLE, MPI_SUM, comm);

    free(subsin);
    free(subcos);
  }

  return 0;
}

/*****************************************************************************
 *
 *  ewald_charge_sum_sin_cos_colloids
 *
 *  Contribution to S(k) and C(k) from colloid point charges.
 *
 *  S(k) += sum_i q_i * sin(k.r_i)
 *  C(k) += sum_i q_i * cos(k.r_i)
 *
 *****************************************************************************/

static int ewald_charge_sum_sin_cos_colloids(ewald_charge_t* ewald) {

  double fkx, fky, fkz;
  double ltot[3];
  PI_DOUBLE(pi);
  colloid_t* pc = NULL;

  assert(ewald);
  if (ewald->cinfo == NULL) return 0;

  cs_ltot(ewald->cs, ltot);

  fkx = 2.0 * pi / ltot[X];
  fky = 2.0 * pi / ltot[Y];
  fkz = 2.0 * pi / ltot[Z];

  colloids_info_local_head(ewald->cinfo, &pc);

  for (; pc; pc = pc->nextlocal) {

    int kn = 0;
    double q = pc->s.q0 - pc->s.q1;  /* Net charge of particle */

    // if (fabs(q) < 1.0e-12) continue;  /* Skip uncharged particles */

    ewald_charge_set_kr_table(ewald, pc->s.r);

    /* Loop over wavevectors */
    for (int kz = 0; kz <= nk_[Z]; kz++) {
      for (int ky = -nk_[Y]; ky <= nk_[Y]; ky++) {
        for (int kx = -nk_[X]; kx <= nk_[X]; kx++) {
          double k[3], ksq;
          double sinkr, coskr;
          double skr[3], ckr[3];

          k[X] = fkx * kx;
          k[Y] = fky * ky;
          k[Z] = fkz * kz;
          ksq = k[X] * k[X] + k[Y] * k[Y] + k[Z] * k[Z];

          if (ksq <= 0.0 || ksq > kmax_) continue;
          assert(kn < nktot_);

          skr[X] = sinkr_[3 * abs(kx) + X];
          skr[Y] = sinkr_[3 * abs(ky) + Y];
          skr[Z] = sinkr_[3 * kz + Z];
          ckr[X] = coskr_[3 * abs(kx) + X];
          ckr[Y] = coskr_[3 * abs(ky) + Y];
          ckr[Z] = coskr_[3 * kz + Z];

          if (kx < 0) skr[X] = -skr[X];
          if (ky < 0) skr[Y] = -skr[Y];

          /* sin(k.r) and cos(k.r) using addition formulas */
          sinkr = skr[X] * ckr[Y] * ckr[Z] + ckr[X] * skr[Y] * ckr[Z]
            + ckr[X] * ckr[Y] * skr[Z] - skr[X] * skr[Y] * skr[Z];

          coskr = ckr[X] * ckr[Y] * ckr[Z] - ckr[X] * skr[Y] * skr[Z]
            - skr[X] * ckr[Y] * skr[Z] - skr[X] * skr[Y] * ckr[Z];

          sinx_[kn] += q * sinkr;
          cosx_[kn] += q * coskr;

          kn++;
        }
      }
    }
  }

  return 0;
}

/*****************************************************************************
 *
 *  ewald_charge_sum_sin_cos_lattice
 *
 *  Contribution to S(k) and C(k) from lattice charge density.
 *
 *  For each lattice node at position r with charge density rho_elec:
 *    q_node = rho_elec * dV  (dV = 1 in lattice units)
 *
 *  S(k) += sum_nodes q_node * sin(k.r_node)
 *  C(k) += sum_nodes q_node * cos(k.r_node)
 *
 *****************************************************************************/

static int ewald_charge_sum_sin_cos_lattice(ewald_charge_t* ewald) {

  double fkx, fky, fkz;
  double ltot[3];
  int nlocal[3], noffset[3];
  PI_DOUBLE(pi);

  assert(ewald);
  if (ewald->psi == NULL) return 0;

  cs_ltot(ewald->cs, ltot);
  cs_nlocal(ewald->cs, nlocal);
  cs_nlocal_offset(ewald->cs, noffset);

  fkx = 2.0 * pi / ltot[X];
  fky = 2.0 * pi / ltot[Y];
  fkz = 2.0 * pi / ltot[Z];

  /* Loop over local lattice nodes */
  for (int ic = 1; ic <= nlocal[X]; ic++) {
    for (int jc = 1; jc <= nlocal[Y]; jc++) {
      for (int kc = 1; kc <= nlocal[Z]; kc++) {

        int index = cs_index(ewald->cs, ic, jc, kc);

        /* Skip solid nodes if map is provided */
        if (ewald->map) {
          int status;
          map_status(ewald->map, index, &status);
          if (status != MAP_FLUID) continue;
        }

        double rho_elec;
        psi_rho_elec(ewald->psi, index, &rho_elec);

        /* q = rho_elec * dV = rho_elec (dV = 1 in lattice units) */
        double q = rho_elec;
        // if (fabs(q) < 1.0e-14) continue;  /* Skip nodes with negligible charge */

        /* Position of this node (global coordinates) */
        double r[3];
        r[X] = (double)(noffset[X] + ic);
        r[Y] = (double)(noffset[Y] + jc);
        r[Z] = (double)(noffset[Z] + kc);

        ewald_charge_set_kr_table(ewald, r);

        /* Loop over wavevectors */
        int kn = 0;
        for (int kzi = 0; kzi <= nk_[Z]; kzi++) {
          for (int kyi = -nk_[Y]; kyi <= nk_[Y]; kyi++) {
            for (int kxi = -nk_[X]; kxi <= nk_[X]; kxi++) {
              double k[3], ksq;
              double sinkr, coskr;
              double skr[3], ckr[3];

              k[X] = fkx * kxi;
              k[Y] = fky * kyi;
              k[Z] = fkz * kzi;
              ksq = k[X] * k[X] + k[Y] * k[Y] + k[Z] * k[Z];

              if (ksq <= 0.0 || ksq > kmax_) continue;
              assert(kn < nktot_);

              skr[X] = sinkr_[3 * abs(kxi) + X];
              skr[Y] = sinkr_[3 * abs(kyi) + Y];
              skr[Z] = sinkr_[3 * kzi + Z];
              ckr[X] = coskr_[3 * abs(kxi) + X];
              ckr[Y] = coskr_[3 * abs(kyi) + Y];
              ckr[Z] = coskr_[3 * kzi + Z];

              if (kxi < 0) skr[X] = -skr[X];
              if (kyi < 0) skr[Y] = -skr[Y];

              sinkr = skr[X] * ckr[Y] * ckr[Z] + ckr[X] * skr[Y] * ckr[Z]
                + ckr[X] * ckr[Y] * skr[Z] - skr[X] * skr[Y] * skr[Z];

              coskr = ckr[X] * ckr[Y] * ckr[Z] - ckr[X] * skr[Y] * skr[Z]
                - skr[X] * ckr[Y] * skr[Z] - skr[X] * skr[Y] * ckr[Z];

              sinx_[kn] += q * sinkr;
              cosx_[kn] += q * coskr;

              kn++;
            }
          }
        }
      }
    }
  }

  return 0;
}

/*****************************************************************************
 *
 *  ewald_charge_self_energy
 *
 *  Return the value of the self energy term.
 *
 *  E_self = -alpha / (sqrt(pi) * 4*pi*epsilon) * sum_i q_i^2
 *
 *  Includes both colloid and lattice contributions.
 *
 *****************************************************************************/

int ewald_charge_self_energy(ewald_charge_t* ewald, double* eself) {

  double sum_q2 = 0.0;
  PI_DOUBLE(pi);

  assert(ewald);

  /* Colloid contribution */
  if ((ewald->sources & EWALD_SOURCE_COLLOIDS) && ewald->cinfo) {
    colloid_t* pc = NULL;
    colloids_info_local_head(ewald->cinfo, &pc);
    for (; pc; pc = pc->nextlocal) {
      double q = pc->s.q0 - pc->s.q1;
      sum_q2 += q * q;
    }
  }

  /* Lattice contribution */
  if ((ewald->sources & EWALD_SOURCE_LATTICE) && ewald->psi) {
    int nlocal[3];
    cs_nlocal(ewald->cs, nlocal);

    for (int ic = 1; ic <= nlocal[X]; ic++) {
      for (int jc = 1; jc <= nlocal[Y]; jc++) {
        for (int kc = 1; kc <= nlocal[Z]; kc++) {
          int index = cs_index(ewald->cs, ic, jc, kc);

          if (ewald->map) {
            int status;
            map_status(ewald->map, index, &status);
            if (status != MAP_FLUID) continue;
          }

          double rho_elec;
          psi_rho_elec(ewald->psi, index, &rho_elec);
          double q = rho_elec;  /* q = rho_elec * dV, dV = 1 */
          sum_q2 += q * q;
        }
      }
    }
  }

  /* MPI reduction */
  {
    MPI_Comm comm = MPI_COMM_NULL;
    cs_cart_comm(ewald->cs, &comm);
    MPI_Allreduce(MPI_IN_PLACE, &sum_q2, 1, MPI_DOUBLE, MPI_SUM, comm);
  }

  *eself = -alpha_ * rpi_ * sum_q2 / (4.0 * pi * epsilon_);

  return 0;
}

/*****************************************************************************
 *
 *  ewald_charge_total_energy
 *
 *  Return the contributions to the energy.
 *
 *****************************************************************************/

int ewald_charge_total_energy(ewald_charge_t* ewald, double* ereal,
                              double* efour, double* eself) {

  if (ewald) {
    *ereal = ereal_;
    *efour = efourier_;
    ewald_charge_self_energy(ewald, eself);
  }
  else {
    *ereal = 0.0;
    *efour = 0.0;
    *eself = 0.0;
  }

  return 0;
}

/*****************************************************************************
 *
 *  ewald_charge_real_space_sum
 *
 *  Look for interactions in real space and accumulate the force
 *  on each particle involved.
 *
 *  This computes:
 *    - Colloid-colloid interactions
 *    - Colloid-lattice interactions
 *    - (Lattice-lattice is handled by FFT solver, not here)
 *
 *  Force: F = q1 * q2 * r / (4*pi*epsilon*r^3)
 *           * [erfc(alpha*r) + 2*alpha*r/sqrt(pi) * exp(-alpha^2*r^2)]
 *
 *****************************************************************************/

int ewald_charge_real_space_sum(ewald_charge_t* ewald) {

  int ncell[3];
  PI_DOUBLE(pi);

  TIMER_start(TIMER_EWALD_REAL_SPACE);

  assert(ewald);

  ereal_ = 0.0;

  /* Colloid-colloid real space interactions */
  if ((ewald->sources & EWALD_SOURCE_COLLOIDS) && ewald->cinfo) {

    colloids_info_ncell(ewald->cinfo, ncell);

    for (int ic = 1; ic <= ncell[X]; ic++) {
      for (int jc = 1; jc <= ncell[Y]; jc++) {
        for (int kc = 1; kc <= ncell[Z]; kc++) {

          colloid_t* p_c1;
          colloids_info_cell_list_head(ewald->cinfo, ic, jc, kc, &p_c1);

          for (; p_c1; p_c1 = p_c1->next) {

            double q1 = p_c1->s.q0 - p_c1->s.q1;
            // if (fabs(q1) < 1.0e-12) continue;

            for (int dx = -1; dx <= +1; dx++) {
              for (int dy = -1; dy <= +1; dy++) {
                for (int dz = -1; dz <= +1; dz++) {

                  int id = ic + dx;
                  int jd = jc + dy;
                  int kd = kc + dz;

                  colloid_t* p_c2;
                  colloids_info_cell_list_head(ewald->cinfo, id, jd, kd, &p_c2);

                  for (; p_c2; p_c2 = p_c2->next) {

                    double q2 = p_c2->s.q0 - p_c2->s.q1;
                    if (fabs(q2) < 1.0e-12) continue;

                    if (p_c1->s.index < p_c2->s.index) {
                      double r12[3], r;

                      cs_minimum_distance(ewald->cs, p_c2->s.r, p_c1->s.r, r12);
                      r = sqrt(r12[X] * r12[X] + r12[Y] * r12[Y] + r12[Z] * r12[Z]);

                      if (r < ewald_rc_ && r > 0.0) {
                        double rr = 1.0 / r;
                        double alpha_r = alpha_ * r;
                        // double prefac = q1*q2 / (4.0*pi*epsilon_);
                        double prefac = eunit_ * beta_ / (4.0 * pi * epsilon_); // Transformo en campo electrico
                        double erfc_term, exp_term;
                        double f_mag;

                        erfc_term = erfc(alpha_r);
                        ereal_ += prefac * erfc_term * rr;

                        exp_term = 2.0 * alpha_ * rpi_ * exp(-alpha_r * alpha_r);
                        f_mag = prefac * rr * rr * (erfc_term + exp_term * r);

                        for (int i = 0; i < 3; i++) {
                          double f_i = f_mag * r12[i] * rr;
                          // p_c1->force[i] -= f_i; 
                          // p_c2->force[i] += f_i; 
                          p_c1->Esub[i] -= f_i * q2; // Transformo en campo electrico
                          p_c2->Esub[i] += f_i * q1; // Transformo en campo electrico
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  /* Colloid-lattice real space interactions */
  if ((ewald->sources & EWALD_SOURCE_COLLOIDS) &&
      (ewald->sources & EWALD_SOURCE_LATTICE) &&
      ewald->cinfo && ewald->psi) {

    int nlocal[3], noffset[3];
    cs_nlocal(ewald->cs, nlocal);
    cs_nlocal_offset(ewald->cs, noffset);

    colloid_t* pc = NULL;
    colloids_info_local_head(ewald->cinfo, &pc);

    for (; pc; pc = pc->nextlocal) {

      double q_coll = pc->s.q0 - pc->s.q1;
      // if (fabs(q_coll) < 1.0e-12) continue;

      /* Search nearby lattice nodes within cutoff */
      int ic_min = (int)floor(pc->s.r[X] - noffset[X] - ewald_rc_);
      int ic_max = (int)ceil(pc->s.r[X] - noffset[X] + ewald_rc_);
      int jc_min = (int)floor(pc->s.r[Y] - noffset[Y] - ewald_rc_);
      int jc_max = (int)ceil(pc->s.r[Y] - noffset[Y] + ewald_rc_);
      int kc_min = (int)floor(pc->s.r[Z] - noffset[Z] - ewald_rc_);
      int kc_max = (int)ceil(pc->s.r[Z] - noffset[Z] + ewald_rc_);

      /* Clamp to local domain */
      if (ic_min < 1) ic_min = 1;
      if (ic_max > nlocal[X]) ic_max = nlocal[X];
      if (jc_min < 1) jc_min = 1;
      if (jc_max > nlocal[Y]) jc_max = nlocal[Y];
      if (kc_min < 1) kc_min = 1;
      if (kc_max > nlocal[Z]) kc_max = nlocal[Z];

      for (int ic = ic_min; ic <= ic_max; ic++) {
        for (int jc = jc_min; jc <= jc_max; jc++) {
          for (int kc = kc_min; kc <= kc_max; kc++) {

            int index = cs_index(ewald->cs, ic, jc, kc);

            if (ewald->map) {
              int status;
              map_status(ewald->map, index, &status);
              if (status != MAP_FLUID) continue;
            }

            double rho_elec;
            psi_rho_elec(ewald->psi, index, &rho_elec);
            double q_lat = rho_elec;

            // if (fabs(q_lat) < 1.0e-14) continue;

            /* Position of lattice node */
            double r_lat[3];
            r_lat[X] = (double)(noffset[X] + ic);
            r_lat[Y] = (double)(noffset[Y] + jc);
            r_lat[Z] = (double)(noffset[Z] + kc);

            /* Distance from colloid to lattice node */
            double r12[3], r;
            cs_minimum_distance(ewald->cs, r_lat, pc->s.r, r12);
            r = sqrt(r12[X] * r12[X] + r12[Y] * r12[Y] + r12[Z] * r12[Z]);

            if (r < ewald_rc_ && r > 0.0) {
              double rr = 1.0 / r;
              double alpha_r = alpha_ * r;
              // double prefac = q_coll * q_lat / (4.0*pi*epsilon_);
              double prefac = q_lat / (4.0 * pi * epsilon_);   // Transformo en campo electrico
              double erfc_term, exp_term;
              double f_mag;

              erfc_term = erfc(alpha_r);
              ereal_ += prefac * erfc_term * rr;

              exp_term = 2.0 * alpha_ * rpi_ * exp(-alpha_r * alpha_r);
              f_mag = prefac * rr * rr * (erfc_term + exp_term * r);

              /* Force on colloid from lattice charge */
              double force_on_particle[3];
              for (int i = 0; i < 3; i++) {
                double f_i = f_mag * r12[i] * rr;
                pc->Esub[i] += f_i; // Campo electrico en la particula
                /* F_particle = q_coll * E, so reaction force on fluid = -F_particle */
                force_on_particle[i] = q_coll * f_i;
              }

              /* Apply reaction force to fluid node (Newton III) */
              if (ewald->hydro) {
                double f_reaction[3];
                f_reaction[X] = -force_on_particle[X];
                f_reaction[Y] = -force_on_particle[Y];
                f_reaction[Z] = -force_on_particle[Z];
                hydro_f_local_add(ewald->hydro, index, f_reaction);
              }
            }
          }
        }
      }
    }
  }

  hydro_memcpy(ewald->hydro, tdpMemcpyHostToDevice);

  TIMER_stop(TIMER_EWALD_REAL_SPACE);

  return 0;
}

/*****************************************************************************
 *
 *  ewald_charge_fourier_space_sum
 *
 *  Accumulate the force on each colloid particle arising from
 *  the Fourier space part of the Ewald sum.
 *
 *  Force on particle i:
 *    F_i = q_i * (4*pi)/(V*epsilon) * sum_{k!=0} (k/k^2) * exp(-k^2/4*alpha^2)
 *          * [S(k)*sin(k.r_i) - C(k)*cos(k.r_i)]
 *
 *  Note: S(k) and C(k) include contributions from BOTH colloids and lattice.
 *
 *****************************************************************************/

int ewald_charge_fourier_space_sum(ewald_charge_t* ewald) {

  double k[3], ksq;
  double b0, b;
  double fkx, fky, fkz;
  double r4alpha_sq;
  double ltot[3];

  int kx, ky, kz, kn = 0;
  int ncell[3];
  PI_DOUBLE(pi);

  TIMER_start(TIMER_EWALD_FOURIER_SPACE);

  assert(ewald);

  cs_ltot(ewald->cs, ltot);
  ewald_charge_sum_sin_cos_terms(ewald);

  fkx = 2.0 * pi / ltot[X];
  fky = 2.0 * pi / ltot[Y];
  fkz = 2.0 * pi / ltot[Z];
  r4alpha_sq = 1.0 / (4.0 * alpha_ * alpha_);

  /* Fourier prefactor: 1/(V*epsilon) - no 4*pi here as it's in real space */
  b0 = 1.0 / (ltot[X] * ltot[Y] * ltot[Z] * epsilon_);

  /* Forces on colloids */
  if ((ewald->sources & EWALD_SOURCE_COLLOIDS) && ewald->cinfo) {

    colloids_info_ncell(ewald->cinfo, ncell);

    for (int ic = 1; ic <= ncell[X]; ic++) {
      for (int jc = 1; jc <= ncell[Y]; jc++) {
        for (int kc = 1; kc <= ncell[Z]; kc++) {

          colloid_t* p_colloid;
          colloids_info_cell_list_head(ewald->cinfo, ic, jc, kc, &p_colloid);

          for (; p_colloid; p_colloid = p_colloid->next) {

            double q = p_colloid->s.q0 - p_colloid->s.q1;
            // if (fabs(q) < 1.0e-12) continue;

            double f[3] = { 0.0, 0.0, 0.0 };

            ewald_charge_set_kr_table(ewald, p_colloid->s.r);

            efourier_ = 0.0;

            kn = 0;
            for (kz = 0; kz <= nk_[Z]; kz++) {
              for (ky = -nk_[Y]; ky <= nk_[Y]; ky++) {
                for (kx = -nk_[X]; kx <= nk_[X]; kx++) {

                  double coskr, sinkr, ckr[3], skr[3];

                  k[X] = fkx * kx;
                  k[Y] = fky * ky;
                  k[Z] = fkz * kz;
                  ksq = k[X] * k[X] + k[Y] * k[Y] + k[Z] * k[Z];

                  if (ksq <= 0.0 || ksq > kmax_) continue;

                  b = b0 * exp(-r4alpha_sq * ksq) / ksq;

                  if (kz > 0) b *= 2.0;
                  efourier_ += 0.5 * b * (sinx_[kn] * sinx_[kn] + cosx_[kn] * cosx_[kn]);

                  skr[X] = sinkr_[3 * abs(kx) + X];
                  skr[Y] = sinkr_[3 * abs(ky) + Y];
                  skr[Z] = sinkr_[3 * kz + Z];
                  ckr[X] = coskr_[3 * abs(kx) + X];
                  ckr[Y] = coskr_[3 * abs(ky) + Y];
                  ckr[Z] = coskr_[3 * kz + Z];

                  if (kx < 0) skr[X] = -skr[X];
                  if (ky < 0) skr[Y] = -skr[Y];

                  sinkr = skr[X] * ckr[Y] * ckr[Z] + ckr[X] * skr[Y] * ckr[Z]
                    + ckr[X] * ckr[Y] * skr[Z] - skr[X] * skr[Y] * skr[Z];

                  coskr = ckr[X] * ckr[Y] * ckr[Z] - ckr[X] * skr[Y] * skr[Z]
                    - skr[X] * ckr[Y] * skr[Z] - skr[X] * skr[Y] * ckr[Z];

                  for (int i = 0; i < 3; i++) {
                    f[i] += b * k[i] * (sinx_[kn] * sinkr - cosx_[kn] * coskr);
                  }

                  kn++;
                }
              }
            }

            for (int i = 0; i < 3; i++) {
              // p_colloid->force[i] += q * f[i];
              p_colloid->Esub[i] += f[i]; // Transformo en campo electrico
            }
          }
        }
      }
    }
  }

  TIMER_stop(TIMER_EWALD_FOURIER_SPACE);

  return 0;
}

/*****************************************************************************
 *
 *  ewald_charge_get_number_fourier_terms
 *
 *****************************************************************************/

static int ewald_charge_get_number_fourier_terms(ewald_charge_t* ewald) {

  int kx, ky, kz, kn = 0;
  double k[3], ksq;
  double fkx, fky, fkz;
  double ltot[3];
  PI_DOUBLE(pi);

  assert(ewald);

  cs_ltot(ewald->cs, ltot);

  fkx = 2.0 * pi / ltot[X];
  fky = 2.0 * pi / ltot[Y];
  fkz = 2.0 * pi / ltot[Z];

  for (kz = 0; kz <= nk_[Z]; kz++) {
    for (ky = -nk_[Y]; ky <= nk_[Y]; ky++) {
      for (kx = -nk_[X]; kx <= nk_[X]; kx++) {

        k[0] = fkx * kx;
        k[1] = fky * ky;
        k[2] = fkz * kz;
        ksq = k[0] * k[0] + k[1] * k[1] + k[2] * k[2];

        if (ksq <= 0.0 || ksq > kmax_) continue;
        kn++;
      }
    }
  }

  return kn;
}

/*****************************************************************************
 *
 *  ewald_charge_set_kr_table
 *
 *****************************************************************************/

static int ewald_charge_set_kr_table(ewald_charge_t* ewald, double r[3]) {

  int i, k;
  double c2[3];
  double ltot[3];
  PI_DOUBLE(pi);

  assert(ewald);

  cs_ltot(ewald->cs, ltot);

  for (i = 0; i < 3; i++) {
    sinkr_[3 * 0 + i] = 0.0;
    coskr_[3 * 0 + i] = 1.0;
    sinkr_[3 * 1 + i] = sin(2.0 * pi * r[i] / ltot[i]);
    coskr_[3 * 1 + i] = cos(2.0 * pi * r[i] / ltot[i]);
    c2[i] = 2.0 * coskr_[3 * 1 + i];
  }

  for (k = 2; k < nkmax_; k++) {
    for (i = 0; i < 3; i++) {
      sinkr_[3 * k + i] = c2[i] * sinkr_[3 * (k - 1) + i] - sinkr_[3 * (k - 2) + i];
      coskr_[3 * k + i] = c2[i] * coskr_[3 * (k - 1) + i] - coskr_[3 * (k - 2) + i];
    }
  }

  return 0;
}

/*****************************************************************************
 *
 *  ewald_charge_potential_at_colloid
 *
 *  Compute the electrostatic potential at a colloid position.
 *
 *****************************************************************************/

int ewald_charge_potential_at_colloid(ewald_charge_t* ewald, colloid_t* pc,
                                      double* phi) {
  double r_target[3];

  assert(ewald);
  assert(pc);
  assert(phi);

  /* Must precompute S(k) and C(k) before calling potential_at_node */
  ewald_charge_sum_sin_cos_terms(ewald);

  r_target[X] = pc->s.r[X];
  r_target[Y] = pc->s.r[Y];
  r_target[Z] = pc->s.r[Z];

  return ewald_charge_potential_at_node(ewald,
           (int)r_target[X], (int)r_target[Y], (int)r_target[Z], phi);
}

/*****************************************************************************
 *
 *  ewald_charge_potential_at_node
 *
 *  Compute the electrostatic potential at a lattice node position.
 *  This includes contributions from all other charges (excluding self).
 *
 *****************************************************************************/

int ewald_charge_potential_at_node(ewald_charge_t* ewald, int ic, int jc, int kc,
                                   double* phi) {
  double phi_real = 0.0;
  double phi_fourier = 0.0;
  double ltot[3];
  int noffset[3];
  PI_DOUBLE(pi);

  assert(ewald);
  assert(phi);

  cs_ltot(ewald->cs, ltot);
  cs_nlocal_offset(ewald->cs, noffset);

  /* Target position */
  double r_target[3];
  r_target[X] = (double)(noffset[X] + ic);
  r_target[Y] = (double)(noffset[Y] + jc);
  r_target[Z] = (double)(noffset[Z] + kc);

  /* Real space contribution from colloids */
  if ((ewald->sources & EWALD_SOURCE_COLLOIDS) && ewald->cinfo) {
    colloid_t* pc = NULL;
    colloids_info_local_head(ewald->cinfo, &pc);

    for (; pc; pc = pc->nextlocal) {
      double q = pc->s.q0 - pc->s.q1;
      // if (fabs(q) < 1.0e-12) continue;

      double r12[3], r;
      cs_minimum_distance(ewald->cs, pc->s.r, r_target, r12);
      r = sqrt(r12[X] * r12[X] + r12[Y] * r12[Y] + r12[Z] * r12[Z]);

      if (r < ewald_rc_ && r > 0.0) {
        phi_real += q * erfc(alpha_ * r) / (4.0 * pi * epsilon_ * r);
        // phi_real += q*eunit_*erfc(alpha_*r) / (4.0*pi*epsilon_*beta_*r);
      }
    }
  }

  /* Fourier space contribution */
  /* NOTE: sinx_[], cosx_[] must be precomputed by caller via
   * ewald_charge_sum_sin_cos_terms() before calling this function */
  {
    double fkx = 2.0 * pi / ltot[X];
    double fky = 2.0 * pi / ltot[Y];
    double fkz = 2.0 * pi / ltot[Z];
    /* Fourier prefactor: 1/(V*epsilon) - no 4*pi here as it's in real space */
    double b0 = 1.0 / (ltot[X] * ltot[Y] * ltot[Z] * epsilon_);
    double r4alpha_sq = 1.0 / (4.0 * alpha_ * alpha_);

    ewald_charge_set_kr_table(ewald, r_target);

    int kn = 0;
    for (int kzi = 0; kzi <= nk_[Z]; kzi++) {
      for (int kyi = -nk_[Y]; kyi <= nk_[Y]; kyi++) {
        for (int kxi = -nk_[X]; kxi <= nk_[X]; kxi++) {
          double k[3], ksq;
          double b;
          double sinkr, coskr, skr[3], ckr[3];

          k[X] = fkx * kxi;
          k[Y] = fky * kyi;
          k[Z] = fkz * kzi;
          ksq = k[X] * k[X] + k[Y] * k[Y] + k[Z] * k[Z];

          if (ksq <= 0.0 || ksq > kmax_) continue;

          b = b0 * exp(-r4alpha_sq * ksq) / ksq;
          if (kzi > 0) b *= 2.0;

          skr[X] = sinkr_[3 * abs(kxi) + X];
          skr[Y] = sinkr_[3 * abs(kyi) + Y];
          skr[Z] = sinkr_[3 * kzi + Z];
          ckr[X] = coskr_[3 * abs(kxi) + X];
          ckr[Y] = coskr_[3 * abs(kyi) + Y];
          ckr[Z] = coskr_[3 * kzi + Z];

          if (kxi < 0) skr[X] = -skr[X];
          if (kyi < 0) skr[Y] = -skr[Y];

          sinkr = skr[X] * ckr[Y] * ckr[Z] + ckr[X] * skr[Y] * ckr[Z]
            + ckr[X] * ckr[Y] * skr[Z] - skr[X] * skr[Y] * skr[Z];

          coskr = ckr[X] * ckr[Y] * ckr[Z] - ckr[X] * skr[Y] * skr[Z]
            - skr[X] * ckr[Y] * skr[Z] - skr[X] * skr[Y] * ckr[Z];

          /* Potential: phi = sum_k b * [S(k)*cos(k.r) + C(k)*sin(k.r)] */
          phi_fourier += b * (sinx_[kn] * sinkr + cosx_[kn] * coskr);

          kn++;
        }
      }
    }
  }

  /* External field contribution: phi_ext = -E0 . r */
  double phi_ext = 0.0;
  if (ewald->psi) {
    phi_ext = -(ewald->psi->e0[X] * r_target[X]
              + ewald->psi->e0[Y] * r_target[Y]
              + ewald->psi->e0[Z] * r_target[Z]);
  }

  *phi = (phi_real + phi_fourier + phi_ext) * eunit_ * beta_;

  return 0;
}

/*****************************************************************************
 *
 *  ewald_charge_compute_psi
 *
 *  Compute and store the potential on the entire lattice using Ewald sum.
 *  This can be used as an alternative to the FFT-based solver.
 *
 *  WARNING: This is O(N * n_k) per node, so it's very expensive for
 *  large systems. Use only for validation or small systems.
 *
 *****************************************************************************/

int ewald_charge_compute_psi(ewald_charge_t* ewald) {

  int nlocal[3];
  int ntotal_nodes;
  int node_count = 0;

  assert(ewald);
  assert(ewald->psi);

  cs_nlocal(ewald->cs, nlocal);
  ntotal_nodes = nlocal[X] * nlocal[Y] * nlocal[Z];

  pe_info(ewald->pe, "\n");
  pe_info(ewald->pe, "Ewald charge: computing psi on lattice...\n");
  pe_info(ewald->pe, "  Local nodes: %d x %d x %d = %d\n",
          nlocal[X], nlocal[Y], nlocal[Z], ntotal_nodes);
  pe_info(ewald->pe, "  Fourier terms per node: %d\n", nktot_);
  pe_info(ewald->pe, "  Fourier terms limit per node: %d\n", nk_limit_);
  pe_info(ewald->pe, "  WARNING: This is O(N * n_k) and may be slow!\n");

  /* Precompute S(k) and C(k) once for ALL nodes */
  pe_info(ewald->pe, "  Computing S(k) and C(k)...\n");
  ewald_charge_sum_sin_cos_terms(ewald);
  pe_info(ewald->pe, "  Done. Now computing potential at each node...\n");

  for (int ic = 1; ic <= nlocal[X]; ic++) {

    /* Progress indicator every 10% */
    if (ic % (nlocal[X] / 10 + 1) == 0) {
      pe_info(ewald->pe, "  Progress: %d%%\n", (100 * ic) / nlocal[X]);
    }

    for (int jc = 1; jc <= nlocal[Y]; jc++) {
      for (int kc = 1; kc <= nlocal[Z]; kc++) {

        int index = cs_index(ewald->cs, ic, jc, kc);

        /* Skip solid nodes */
        if (ewald->map) {
          int status;
          map_status(ewald->map, index, &status);
          if (status != MAP_FLUID) continue;
        }

        double phi;
        ewald_charge_potential_at_node(ewald, ic, jc, kc, &phi);

        /* Store in psi field */
        psi_psi_set(ewald->psi, index, phi);
        node_count++;
      }
    }
  }

  pe_info(ewald->pe, "  Computed potential at %d nodes.\n", node_count);

  /* Update halos */
  pe_info(ewald->pe, "  Updating halos...\n");
  psi_halo_psi(ewald->psi);

  pe_info(ewald->pe, "Ewald charge: psi computation complete.\n\n");

  return 0;
}

/*****************************************************************************
 *
 *  GPU Implementation of Ewald Sum
 *
 *  This section contains CUDA kernels and the GPU version of
 *  ewald_charge_sum_full for accelerated computation.
 *
 *****************************************************************************/

#ifdef __NVCC__

#include <cuda_runtime.h>
#include <cufft.h>
#include "target.h"

 /* GPU constants - copied to device at initialization */
__constant__ double d_alpha;
__constant__ double d_epsilon;
__constant__ double d_beta;
__constant__ double d_eunit;
__constant__ double d_rpi;
__constant__ double d_ewald_rc;
__constant__ double d_eps_reg;   /* Real-space regularization: 1/r -> 1/sqrt(r^2+eps_reg^2) */
__constant__ double d_fkx, d_fky, d_fkz;
__constant__ double d_r4alpha_sq;
__constant__ double d_b0;
__constant__ double d_kmax;
__constant__ int d_nk[3];
__constant__ double d_E_dipole[3];        /* Dipole correction field */
__constant__ double d_dipole_prefactor;   /* 4*pi / (3*V*epsilon) */
__constant__ double d_M_dipole[3];        /* Dipole moment */
__constant__ int    d_nlocal[3];
__constant__ int    d_noffset[3];
__constant__ double d_ltot[3];   /* total system size, for minimum image in real-space */

/*****************************************************************************
 *
 *  ewald_peskin_1d
 *
 *  GPU device function: Peskin delta kernel for one dimension.
 *  Same formula as d_peskin() in subgrid.c.
 *  Support: |r| <= 2; returns 0 outside.
 *
 *****************************************************************************/

__device__ __forceinline__ double ewald_peskin_1d(double r) {
  double rmod = fabs(r);
  if (rmod <= 1.0) {
    return 0.125 * (3.0 - 2.0 * rmod + sqrt(1.0 + 4.0 * rmod - 4.0 * rmod * rmod));
  }
  else if (rmod <= 2.0) {
    return 0.125 * (5.0 - 2.0 * rmod - sqrt(-7.0 + 12.0 * rmod - 4.0 * rmod * rmod));
  }
  return 0.0;
}

/*****************************************************************************
 *
 *  ewald_structure_factor_lattice_kernel
 *
 *  Compute structure factor contribution from lattice nodes.
 *  Each thread handles one node, atomically adds to Sk_sin/Sk_cos.
 *
 *****************************************************************************/

__global__ void ewald_structure_factor_lattice_kernel(
    const double* __restrict__ rho_data,
    double* __restrict__ Sk_sin,
    double* __restrict__ Sk_cos,
    const double* __restrict__ kvec,
    int nktot,
    int nsites,
    int nhalo) {

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int ntotal = d_nlocal[0] * d_nlocal[1] * d_nlocal[2];

  if (idx >= ntotal) return;

  /* Convert linear index to (ic, jc, kc) in Ludwig local coords (1-based) */
  int kc = idx % d_nlocal[2] + 1;
  int jc = (idx / d_nlocal[2]) % d_nlocal[1] + 1;
  int ic = idx / (d_nlocal[2] * d_nlocal[1]) + 1;

  /* Ludwig index calculation */
  int str_z = 1;
  int str_y = (d_nlocal[2] + 2 * nhalo);
  int str_x = str_y * (d_nlocal[1] + 2 * nhalo);
  int ludwig_idx = str_x * (nhalo + ic - 1) + str_y * (nhalo + jc - 1) + str_z * (nhalo + kc - 1);

  /* Get charge: rho_elec = rho[0] - rho[1] (SOA layout) */
  double rho0 = rho_data[nsites * 0 + ludwig_idx];
  double rho1 = rho_data[nsites * 1 + ludwig_idx];
  double q = rho0 - rho1;

  if (fabs(q) < 1.0e-14) return;

  /* Global position */
  double rx = (double)(d_noffset[0] + ic);
  double ry = (double)(d_noffset[1] + jc);
  double rz = (double)(d_noffset[2] + kc);

  /* Loop over k-vectors and add contributions */
  for (int kn = 0; kn < nktot; kn++) {
    double kx = kvec[3 * kn + 0];
    double ky = kvec[3 * kn + 1];
    double kz = kvec[3 * kn + 2];

    double kr = kx * rx + ky * ry + kz * rz;
    double sinkr = sin(kr);
    double coskr = cos(kr);

    atomicAdd(&Sk_sin[kn], q * sinkr);
    atomicAdd(&Sk_cos[kn], q * coskr);
  }
}

/*****************************************************************************
 *
 *  ewald_structure_factor_particle_kernel
 *
 *  Compute structure factor contribution from colloid particles.
 *  Each thread handles one particle, loops over k-vectors.
 *
 *****************************************************************************/

__global__ void ewald_structure_factor_particle_kernel(
    const double* __restrict__ particle_r,
    const double* __restrict__ particle_q,
    double* __restrict__ Sk_sin,
    double* __restrict__ Sk_cos,
    const double* __restrict__ kvec,
    int nparticles,
    int nktot) {

  int p = blockIdx.x * blockDim.x + threadIdx.x;

  if (p >= nparticles) return;

  double q = particle_q[p];
  if (fabs(q) < 1.0e-14) return;

  /* Particle position (global coordinates) */
  double rx = particle_r[3 * p + 0];
  double ry = particle_r[3 * p + 1];
  double rz = particle_r[3 * p + 2];

  /* Loop over k-vectors and add contributions */
  for (int kn = 0; kn < nktot; kn++) {
    double kx = kvec[3 * kn + 0];
    double ky = kvec[3 * kn + 1];
    double kz = kvec[3 * kn + 2];

    double kr = kx * rx + ky * ry + kz * rz;
    double sinkr = sin(kr);
    double coskr = cos(kr);

    atomicAdd(&Sk_sin[kn], q * sinkr);
    atomicAdd(&Sk_cos[kn], q * coskr);
  }
}

/*****************************************************************************
 *
 *  ewald_potential_fourier_kernel
 *
 *  Compute Fourier part of potential at each lattice node.
 *  Each thread handles one node.
 *
 *****************************************************************************/

__global__ void ewald_potential_fourier_kernel(
    double* __restrict__ psi_data,
    const double* __restrict__ Sk_sin,
    const double* __restrict__ Sk_cos,
    const double* __restrict__ kvec,
    const double* __restrict__ Gk,
    const int* __restrict__ kz_arr,
    int nktot,
    int nsites,
    int nhalo) {

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int ntotal = d_nlocal[0] * d_nlocal[1] * d_nlocal[2];

  if (idx >= ntotal) return;

  /* Convert to Ludwig local coords (1-based) */
  int kc = idx % d_nlocal[2] + 1;
  int jc = (idx / d_nlocal[2]) % d_nlocal[1] + 1;
  int ic = idx / (d_nlocal[2] * d_nlocal[1]) + 1;

  /* Ludwig index */
  int str_z = 1;
  int str_y = (d_nlocal[2] + 2 * nhalo);
  int str_x = str_y * (d_nlocal[1] + 2 * nhalo);
  int ludwig_idx = str_x * (nhalo + ic - 1) + str_y * (nhalo + jc - 1) + str_z * (nhalo + kc - 1);

  /* Global position */
  double rx = (double)(d_noffset[0] + ic);
  double ry = (double)(d_noffset[1] + jc);
  double rz = (double)(d_noffset[2] + kc);

  /* Sum Fourier contributions */
  double phi_fourier = 0.0;

  for (int kn = 0; kn < nktot; kn++) {
    double kx = kvec[3 * kn + 0];
    double ky = kvec[3 * kn + 1];
    double kz_val = kvec[3 * kn + 2];

    double kr = kx * rx + ky * ry + kz_val * rz;
    double sinkr = sin(kr);
    double coskr = cos(kr);

    double factor = (kz_arr[kn] > 0) ? 2.0 : 1.0;
    phi_fourier += factor * Gk[kn] * (Sk_cos[kn] * coskr + Sk_sin[kn] * sinkr);
  }

  /* Dipole correction: phi_dipole = (4*pi / 3*V*epsilon) * M . r = -E_dipole . r */
  double phi_dipole = d_dipole_prefactor * (d_M_dipole[0] * rx + d_M_dipole[1] * ry + d_M_dipole[2] * rz);

  /* Store: psi* = beta * eunit * phi */
  psi_data[ludwig_idx] = d_beta * d_eunit * (phi_fourier + phi_dipole);
}

/*****************************************************************************
 *
 *  ewald_potential_real_particle_kernel
 *
 *  Add real-space contribution from particles to lattice nodes.
 *  Each thread handles one lattice node, loops over nearby particles.
 *
 *****************************************************************************/

__global__ void ewald_potential_real_particle_kernel(
    double* __restrict__ psi_data,
    const double* __restrict__ particle_r,
    const double* __restrict__ particle_q,
    int nparticles,
    int nsites,
    int nhalo,
    int irc) {

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int ntotal = d_nlocal[0] * d_nlocal[1] * d_nlocal[2];

  if (idx >= ntotal) return;

  /* Convert to Ludwig local coords (1-based) */
  int kc = idx % d_nlocal[2] + 1;
  int jc = (idx / d_nlocal[2]) % d_nlocal[1] + 1;
  int ic = idx / (d_nlocal[2] * d_nlocal[1]) + 1;

  /* Ludwig index */
  int str_z = 1;
  int str_y = (d_nlocal[2] + 2 * nhalo);
  int str_x = str_y * (d_nlocal[1] + 2 * nhalo);
  int ludwig_idx = str_x * (nhalo + ic - 1) + str_y * (nhalo + jc - 1) + str_z * (nhalo + kc - 1);

  /* Local position in lattice units */
  double r_local[3] = { (double)ic, (double)jc, (double)kc };

  double phi_real = 0.0;

  /* Loop over all particles */
  for (int p = 0; p < nparticles; p++) {
    double q_p = particle_q[p];
    if (fabs(q_p) < 1.0e-14) continue;

    /* Particle position in local coords */
    double r_p[3];
    r_p[0] = particle_r[3 * p + 0] - (double)d_noffset[0];
    r_p[1] = particle_r[3 * p + 1] - (double)d_noffset[1];
    r_p[2] = particle_r[3 * p + 2] - (double)d_noffset[2];

    double dr[3] = { r_local[0] - r_p[0], r_local[1] - r_p[1], r_local[2] - r_p[2] };
    double dist = sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]);

    if (dist < d_ewald_rc && dist > 1.0e-10) {
      double r_reg = sqrt(dist * dist + d_eps_reg * d_eps_reg);
      /* phi_real += q_p * erfc(d_alpha*dist) / (4.0*M_PI*d_epsilon*dist); */
      phi_real += q_p * erfc(d_alpha * dist) / (4.0 * M_PI * d_epsilon * r_reg);
    }
  }

  /* Add to existing potential (Fourier part already stored) */
  psi_data[ludwig_idx] += d_beta * d_eunit * phi_real;
}

/*****************************************************************************
 *
 *  ewald_force_fourier_kernel
 *
 *  Compute Fourier part of force at each lattice node.
 *  Each thread handles one node.
 *
 *****************************************************************************/

__global__ void ewald_force_fourier_kernel(
    double* __restrict__ force_data,
    const double* __restrict__ rho_data,
    const double* __restrict__ Sk_sin,
    const double* __restrict__ Sk_cos,
    const double* __restrict__ kvec,
    const double* __restrict__ Gk,
    const int* __restrict__ kz_arr,
    int nktot,
    int nsites,
    int nhalo) {

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int ntotal = d_nlocal[0] * d_nlocal[1] * d_nlocal[2];

  if (idx >= ntotal) return;

  /* Convert to Ludwig local coords (1-based) */
  int kc = idx % d_nlocal[2] + 1;
  int jc = (idx / d_nlocal[2]) % d_nlocal[1] + 1;
  int ic = idx / (d_nlocal[2] * d_nlocal[1]) + 1;

  /* Ludwig index */
  int str_z = 1;
  int str_y = (d_nlocal[2] + 2 * nhalo);
  int str_x = str_y * (d_nlocal[1] + 2 * nhalo);
  int ludwig_idx = str_x * (nhalo + ic - 1) + str_y * (nhalo + jc - 1) + str_z * (nhalo + kc - 1);

  /* Get charge at this node */
  double rho0 = rho_data[nsites * 0 + ludwig_idx];
  double rho1 = rho_data[nsites * 1 + ludwig_idx];
  double q1 = rho0 - rho1;

  /* Global position */
  double rx = (double)(d_noffset[0] + ic);
  double ry = (double)(d_noffset[1] + jc);
  double rz = (double)(d_noffset[2] + kc);

  /* Sum Fourier contributions to force */
  double F_fourier[3] = { 0.0, 0.0, 0.0 };

  for (int kn = 0; kn < nktot; kn++) {
    double kx = kvec[3 * kn + 0];
    double ky = kvec[3 * kn + 1];
    double kz_val = kvec[3 * kn + 2];

    double kr = kx * rx + ky * ry + kz_val * rz;
    double sinkr = sin(kr);
    double coskr = cos(kr);

    double factor = (kz_arr[kn] > 0) ? 2.0 : 1.0;
    double im_part = Sk_sin[kn] * coskr - Sk_cos[kn] * sinkr;

    F_fourier[0] += factor * q1 * Gk[kn] * kx * im_part;
    F_fourier[1] += factor * q1 * Gk[kn] * ky * im_part;
    F_fourier[2] += factor * q1 * Gk[kn] * kz_val * im_part;
  }

  /* Dipole correction to force: F_dipole = q * E_dipole */
  F_fourier[0] += q1 * d_E_dipole[0];
  F_fourier[1] += q1 * d_E_dipole[1];
  F_fourier[2] += q1 * d_E_dipole[2];

  /* Store force (hydro force layout: force_data[3*nsites]) */
  /* Using atomic add since multiple kernels may update */
  atomicAdd(&force_data[3 * ludwig_idx + 0], F_fourier[0]);
  atomicAdd(&force_data[3 * ludwig_idx + 1], F_fourier[1]);
  atomicAdd(&force_data[3 * ludwig_idx + 2], F_fourier[2]);
}

/*****************************************************************************
 *
 *  ewald_force_real_particle_kernel
 *
 *  Add real-space force contribution from particles to lattice nodes.
 *
 *****************************************************************************/

__global__ void ewald_force_real_particle_kernel(
    double* __restrict__ force_data,
    const double* __restrict__ rho_data,
    const double* __restrict__ particle_r,
    const double* __restrict__ particle_q,
    int nparticles,
    int nsites,
    int nhalo,
    int irc) {

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int ntotal = d_nlocal[0] * d_nlocal[1] * d_nlocal[2];

  if (idx >= ntotal) return;

  /* Convert to Ludwig local coords (1-based) */
  int kc = idx % d_nlocal[2] + 1;
  int jc = (idx / d_nlocal[2]) % d_nlocal[1] + 1;
  int ic = idx / (d_nlocal[2] * d_nlocal[1]) + 1;

  /* Ludwig index */
  int str_z = 1;
  int str_y = (d_nlocal[2] + 2 * nhalo);
  int str_x = str_y * (d_nlocal[1] + 2 * nhalo);
  int ludwig_idx = str_x * (nhalo + ic - 1) + str_y * (nhalo + jc - 1) + str_z * (nhalo + kc - 1);

  /* Get charge at this node */
  double rho0 = rho_data[nsites * 0 + ludwig_idx];
  double rho1 = rho_data[nsites * 1 + ludwig_idx];
  double q1 = rho0 - rho1;

  /* Local position */
  double r_local[3] = { (double)ic, (double)jc, (double)kc };

  double F_real[3] = { 0.0, 0.0, 0.0 };

  /* Loop over all particles */
  for (int p = 0; p < nparticles; p++) {
    double q_p = particle_q[p];
    if (fabs(q_p) < 1.0e-14) continue;

    /* Particle position in local coords */
    double r_p[3];
    r_p[0] = particle_r[3 * p + 0] - (double)d_noffset[0];
    r_p[1] = particle_r[3 * p + 1] - (double)d_noffset[1];
    r_p[2] = particle_r[3 * p + 2] - (double)d_noffset[2];

    double dr[3] = { r_local[0] - r_p[0], r_local[1] - r_p[1], r_local[2] - r_p[2] };
    double dist = sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]);

    if (dist < d_ewald_rc && dist > 1.0e-10) {
      double ar = d_alpha * dist;
      double r_reg = sqrt(dist * dist + d_eps_reg * d_eps_reg);
      double r_inv = 1.0 / r_reg;
      /* double r_inv = 1.0/dist; */
      /* double term1 = erfc(ar)/(dist*dist); */
      double term1 = erfc(ar) / (r_reg * r_reg);
      double term2 = (2.0 * d_alpha * d_rpi) * exp(-ar * ar) * r_inv;
      double F_mag = q1 * q_p * (term1 + term2) / (4.0 * M_PI * d_epsilon);

      F_real[0] += F_mag * dr[0] * r_inv;
      F_real[1] += F_mag * dr[1] * r_inv;
      F_real[2] += F_mag * dr[2] * r_inv;
    }
  }

  /* Add to force */
  atomicAdd(&force_data[3 * ludwig_idx + 0], F_real[0]);
  atomicAdd(&force_data[3 * ludwig_idx + 1], F_real[1]);
  atomicAdd(&force_data[3 * ludwig_idx + 2], F_real[2]);
}

/*CHANGE INIT - 20251203 Electric field output */
/*****************************************************************************
 *
 *  ewald_efield_fourier_kernel
 *
 *  Compute Fourier part of electric field at each lattice node.
 *  This is the same as ewald_force_fourier_kernel but without multiplying by q.
 *  Each thread handles one node.
 *
 *****************************************************************************/

__global__ void ewald_efield_fourier_kernel(
    double* __restrict__ efield_data,
    const double* __restrict__ Sk_sin,
    const double* __restrict__ Sk_cos,
    const double* __restrict__ kvec,
    const double* __restrict__ Gk,
    const int* __restrict__ kz_arr,
    int nktot,
    int nsites,
    int nhalo) {

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int ntotal = d_nlocal[0] * d_nlocal[1] * d_nlocal[2];

  if (idx >= ntotal) return;

  /* Convert to Ludwig local coords (1-based) */
  int kc = idx % d_nlocal[2] + 1;
  int jc = (idx / d_nlocal[2]) % d_nlocal[1] + 1;
  int ic = idx / (d_nlocal[2] * d_nlocal[1]) + 1;

  /* Ludwig index */
  int str_z = 1;
  int str_y = (d_nlocal[2] + 2 * nhalo);
  int str_x = str_y * (d_nlocal[1] + 2 * nhalo);
  int ludwig_idx = str_x * (nhalo + ic - 1) + str_y * (nhalo + jc - 1) + str_z * (nhalo + kc - 1);

  /* Global position */
  double rx = (double)(d_noffset[0] + ic);
  double ry = (double)(d_noffset[1] + jc);
  double rz = (double)(d_noffset[2] + kc);

  /* Sum Fourier contributions to electric field (without charge) */
  double E_fourier[3] = { 0.0, 0.0, 0.0 };

  for (int kn = 0; kn < nktot; kn++) {
    double kx = kvec[3 * kn + 0];
    double ky = kvec[3 * kn + 1];
    double kz_val = kvec[3 * kn + 2];

    double kr = kx * rx + ky * ry + kz_val * rz;
    double sinkr = sin(kr);
    double coskr = cos(kr);

    double factor = (kz_arr[kn] > 0) ? 2.0 : 1.0;
    /* E = -grad(phi), phi = Gk*(Sk_cos*cos + Sk_sin*sin) with S(k)=sum q_j e^{+ikr_j}
     * -grad(phi)/dk = Gk*k*(Sk_cos*sin - Sk_sin*cos) = -im_part */
    double im_part = Sk_sin[kn] * coskr - Sk_cos[kn] * sinkr;

    E_fourier[0] -= factor * Gk[kn] * kx * im_part;
    E_fourier[1] -= factor * Gk[kn] * ky * im_part;
    E_fourier[2] -= factor * Gk[kn] * kz_val * im_part;
  }

  /* Dipole correction to field */
  E_fourier[0] += d_E_dipole[0];
  E_fourier[1] += d_E_dipole[1];
  E_fourier[2] += d_E_dipole[2];

  /* Store electric field */
  efield_data[3 * ludwig_idx + 0] = E_fourier[0];
  efield_data[3 * ludwig_idx + 1] = E_fourier[1];
  efield_data[3 * ludwig_idx + 2] = E_fourier[2];
}

/*****************************************************************************
 *
 *  ewald_efield_real_particle_kernel
 *
 *  Add real-space electric field contribution from particles to lattice nodes.
 *  This is the same as ewald_force_real_particle_kernel but without multiplying by q_node.
 *
 *****************************************************************************/

__global__ void ewald_efield_real_particle_kernel(
    double* __restrict__ efield_data,
    const double* __restrict__ particle_r,
    const double* __restrict__ particle_q,
    int nparticles,
    int nsites,
    int nhalo,
    int irc) {

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int ntotal = d_nlocal[0] * d_nlocal[1] * d_nlocal[2];

  if (idx >= ntotal) return;

  /* Convert to Ludwig local coords (1-based) */
  int kc = idx % d_nlocal[2] + 1;
  int jc = (idx / d_nlocal[2]) % d_nlocal[1] + 1;
  int ic = idx / (d_nlocal[2] * d_nlocal[1]) + 1;

  /* Ludwig index */
  int str_z = 1;
  int str_y = (d_nlocal[2] + 2 * nhalo);
  int str_x = str_y * (d_nlocal[1] + 2 * nhalo);
  int ludwig_idx = str_x * (nhalo + ic - 1) + str_y * (nhalo + jc - 1) + str_z * (nhalo + kc - 1);

  /* Local position */
  double r_local[3] = { (double)ic, (double)jc, (double)kc };

  double E_real[3] = { 0.0, 0.0, 0.0 };

  /* Loop over all particles */
  for (int p = 0; p < nparticles; p++) {
    double q_p = particle_q[p];
    if (fabs(q_p) < 1.0e-14) continue;

    /* Particle position in local coords */
    double r_p[3];
    r_p[0] = particle_r[3 * p + 0] - (double)d_noffset[0];
    r_p[1] = particle_r[3 * p + 1] - (double)d_noffset[1];
    r_p[2] = particle_r[3 * p + 2] - (double)d_noffset[2];

    double dr[3] = { r_local[0] - r_p[0], r_local[1] - r_p[1], r_local[2] - r_p[2] };
    double dist = sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]);

    if (dist < d_ewald_rc && dist > 1.0e-10) {
      double ar = d_alpha * dist;
      double r_reg = sqrt(dist * dist + d_eps_reg * d_eps_reg);
      double r_inv = 1.0 / r_reg;
      /* double r_inv = 1.0/dist; */
      /* double term1 = erfc(ar)/(dist*dist); */
      double term1 = erfc(ar) / (r_reg * r_reg);
      double term2 = (2.0 * d_alpha * d_rpi) * exp(-ar * ar) * r_inv;
      /* E_mag = F_mag / q_node, so we compute field without q_node */
      double E_mag = q_p * (term1 + term2) / (4.0 * M_PI * d_epsilon);

      E_real[0] += E_mag * dr[0] * r_inv;
      E_real[1] += E_mag * dr[1] * r_inv;
      E_real[2] += E_mag * dr[2] * r_inv;
    }
  }

  /* Add to electric field */
  atomicAdd(&efield_data[3 * ludwig_idx + 0], E_real[0]);
  atomicAdd(&efield_data[3 * ludwig_idx + 1], E_real[1]);
  atomicAdd(&efield_data[3 * ludwig_idx + 2], E_real[2]);
}

/*****************************************************************************
 *
 *  ewald_efield_real_lattice_kernel
 *
 *  Add real-space electric field contribution from lattice nodes to other lattice nodes.
 *  This is the lattice-lattice correction term.
 *
 *****************************************************************************/

__global__ void ewald_efield_real_lattice_kernel(
    double* __restrict__ efield_data,
    const double* __restrict__ rho_data,
    int nsites,
    int nhalo,
    int irc) {

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int ntotal = d_nlocal[0] * d_nlocal[1] * d_nlocal[2];

  if (idx >= ntotal) return;

  /* Convert to Ludwig local coords (1-based) */
  int kc = idx % d_nlocal[2] + 1;
  int jc = (idx / d_nlocal[2]) % d_nlocal[1] + 1;
  int ic = idx / (d_nlocal[2] * d_nlocal[1]) + 1;

  /* Ludwig index */
  int str_z = 1;
  int str_y = (d_nlocal[2] + 2 * nhalo);
  int str_x = str_y * (d_nlocal[1] + 2 * nhalo);
  int ludwig_idx = str_x * (nhalo + ic - 1) + str_y * (nhalo + jc - 1) + str_z * (nhalo + kc - 1);

  /* Local position */
  double r_local[3] = { (double)ic, (double)jc, (double)kc };

  double E_real[3] = { 0.0, 0.0, 0.0 };

  /* Loop over nearby lattice nodes */
  for (int di = -irc; di <= irc; di++) {
    int i2 = ic + di;
    /* Allow access to halos for periodic boundaries */
    if (i2 < 1 - nhalo || i2 > d_nlocal[0] + nhalo) continue;

    for (int dj = -irc; dj <= irc; dj++) {
      int j2 = jc + dj;
      if (j2 < 1 - nhalo || j2 > d_nlocal[1] + nhalo) continue;

      for (int dk = -irc; dk <= irc; dk++) {
        int k2 = kc + dk;
        if (k2 < 1 - nhalo || k2 > d_nlocal[2] + nhalo) continue;

        /* Skip self-interaction */
        if (di == 0 && dj == 0 && dk == 0) continue;

        /* Neighbor index */
        int ludwig_idx2 = str_x * (nhalo + i2 - 1) + str_y * (nhalo + j2 - 1) + str_z * (nhalo + k2 - 1);

        /* Get charge density difference at neighbor node */
        double rho0_2 = rho_data[nsites * 0 + ludwig_idx2];
        double rho1_2 = rho_data[nsites * 1 + ludwig_idx2];
        double q_node = rho0_2 - rho1_2;

        /* Distance vector: from source (neighbor) to receptor (current node) */
        double dr[3] = { -(double)di, -(double)dj, -(double)dk };
        double dist = sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]);

        if (dist < d_ewald_rc && dist > 1.0e-10) {
          double ar = d_alpha * dist;
          double r_reg = sqrt(dist * dist + d_eps_reg * d_eps_reg);
          double r_inv = 1.0 / r_reg;
          /* double r_inv = 1.0/dist; */
          /* double term1 = erfc(ar)/(dist*dist); */
          double term1 = erfc(ar) / (r_reg * r_reg);
          double term2 = (2.0 * d_alpha * d_rpi) * exp(-ar * ar) * r_inv;
          /* Same formula as particle kernel (no d_beta * d_eunit) */
          double E_mag = q_node * (term1 + term2) / (4.0 * M_PI * d_epsilon);

          E_real[0] += E_mag * dr[0] * r_inv;
          E_real[1] += E_mag * dr[1] * r_inv;
          E_real[2] += E_mag * dr[2] * r_inv;
        }
      }
    }
  }

  /* Add to electric field */
  atomicAdd(&efield_data[3 * ludwig_idx + 0], E_real[0]);
  atomicAdd(&efield_data[3 * ludwig_idx + 1], E_real[1]);
  atomicAdd(&efield_data[3 * ludwig_idx + 2], E_real[2]);
}
/*CHANGE END - 20251203 */

/*****************************************************************************
 *
 *  ewald_potential_real_lattice_kernel
 *
 *  Add real-space potential contribution from lattice nodes to other lattice nodes.
 *  This is the lattice-lattice correction term, analogous to ewald_efield_real_lattice_kernel.
 *
 *****************************************************************************/

__global__ void ewald_potential_real_lattice_kernel(
    double* __restrict__ psi_data,
    const double* __restrict__ rho_data,
    int nsites,
    int nhalo,
    int irc) {

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int ntotal = d_nlocal[0] * d_nlocal[1] * d_nlocal[2];

  if (idx >= ntotal) return;

  /* Convert to Ludwig local coords (1-based) */
  int kc = idx % d_nlocal[2] + 1;
  int jc = (idx / d_nlocal[2]) % d_nlocal[1] + 1;
  int ic = idx / (d_nlocal[2] * d_nlocal[1]) + 1;

  /* Ludwig index */
  int str_z = 1;
  int str_y = (d_nlocal[2] + 2 * nhalo);
  int str_x = str_y * (d_nlocal[1] + 2 * nhalo);
  int ludwig_idx = str_x * (nhalo + ic - 1) + str_y * (nhalo + jc - 1) + str_z * (nhalo + kc - 1);

  double phi_real = 0.0;

  /* Loop over nearby lattice nodes */
  for (int di = -irc; di <= irc; di++) {
    int i2 = ic + di;
    if (i2 < 1 - nhalo || i2 > d_nlocal[0] + nhalo) continue;

    for (int dj = -irc; dj <= irc; dj++) {
      int j2 = jc + dj;
      if (j2 < 1 - nhalo || j2 > d_nlocal[1] + nhalo) continue;

      for (int dk = -irc; dk <= irc; dk++) {
        int k2 = kc + dk;
        if (k2 < 1 - nhalo || k2 > d_nlocal[2] + nhalo) continue;

        /* Skip self-interaction */
        if (di == 0 && dj == 0 && dk == 0) continue;

        /* Neighbor index */
        int ludwig_idx2 = str_x * (nhalo + i2 - 1) + str_y * (nhalo + j2 - 1) + str_z * (nhalo + k2 - 1);

        /* Get net charge at neighbor node */
        double rho0_2 = rho_data[nsites * 0 + ludwig_idx2];
        double rho1_2 = rho_data[nsites * 1 + ludwig_idx2];
        double q_node = rho0_2 - rho1_2;

        if (fabs(q_node) < 1.0e-14) continue;

        double dist = sqrt((double)(di * di + dj * dj + dk * dk));

        if (dist < d_ewald_rc && dist > 1.0e-10) {
          double r_reg = sqrt(dist * dist + d_eps_reg * d_eps_reg);
          /* phi_real += q_node * erfc(d_alpha*dist) / (4.0*M_PI*d_epsilon*dist); */
          phi_real += q_node * erfc(d_alpha * dist) / (4.0 * M_PI * d_epsilon * r_reg);
        }
      }
    }
  }

  /* Add to potential: beta * eunit * phi_real (same scaling as particle kernel) */
  atomicAdd(&psi_data[ludwig_idx], d_beta * d_eunit * phi_real);
}

/*****************************************************************************
 *
 *  ewald_particle_field_fourier_kernel
 *
 *  Compute Fourier part of field at each particle position.
 *  Each thread handles one particle.
 *
 *****************************************************************************/

__global__ void ewald_particle_field_fourier_kernel(
    double* __restrict__ Esub_data,
    double* __restrict__ fex_data,
    const double* __restrict__ particle_r,
    const double* __restrict__ particle_q,
    const double* __restrict__ Sk_sin,
    const double* __restrict__ Sk_cos,
    const double* __restrict__ kvec,
    const double* __restrict__ Gk,
    const int* __restrict__ kz_arr,
    int nparticles,
    int nktot) {

  int p = blockIdx.x * blockDim.x + threadIdx.x;
  if (p >= nparticles) return;

  double q_p = particle_q[p];
  double rx = particle_r[3 * p + 0];
  double ry = particle_r[3 * p + 1];
  double rz = particle_r[3 * p + 2];

  /* Sum Fourier contributions to E field */
  double E_fourier[3] = { 0.0, 0.0, 0.0 };

  for (int kn = 0; kn < nktot; kn++) {
    double kx = kvec[3 * kn + 0];
    double ky = kvec[3 * kn + 1];
    double kz_val = kvec[3 * kn + 2];

    double kr = kx * rx + ky * ry + kz_val * rz;
    double sinkr = sin(kr);
    double coskr = cos(kr);

    double factor = (kz_arr[kn] > 0) ? 2.0 : 1.0;
    double im_part = Sk_sin[kn] * coskr - Sk_cos[kn] * sinkr;

    /* Esub = beta * eunit * E_physical */
    E_fourier[0] += factor * d_beta * d_eunit * Gk[kn] * kx * im_part;
    E_fourier[1] += factor * d_beta * d_eunit * Gk[kn] * ky * im_part;
    E_fourier[2] += factor * d_beta * d_eunit * Gk[kn] * kz_val * im_part;
  }

  /* Add dipole correction: Esub_dipole = beta * eunit * E_dipole */
  E_fourier[0] += d_beta * d_eunit * d_E_dipole[0];
  E_fourier[1] += d_beta * d_eunit * d_E_dipole[1];
  E_fourier[2] += d_beta * d_eunit * d_E_dipole[2];

  /* Store Esub (will add real-space contribution later) */
  Esub_data[3 * p + 0] = E_fourier[0];
  Esub_data[3 * p + 1] = E_fourier[1];
  Esub_data[3 * p + 2] = E_fourier[2];

  /* Compute fex = q * Esub * kt / eunit */
  double kt = 1.0 / d_beta;
  fex_data[3 * p + 0] = q_p * E_fourier[0] * kt / d_eunit;
  fex_data[3 * p + 1] = q_p * E_fourier[1] * kt / d_eunit;
  fex_data[3 * p + 2] = q_p * E_fourier[2] * kt / d_eunit;
}

/*****************************************************************************
 *
 *  ewald_particle_field_real_kernel
 *
 *  Add real-space contribution to particle field from lattice and particles.
 *  Each thread handles one particle.
 *
 *****************************************************************************/

__global__ void ewald_particle_field_real_kernel(
    double* __restrict__ Esub_data,
    double* __restrict__ fex_data,
    const double* __restrict__ particle_r,
    const double* __restrict__ particle_q,
    const double* __restrict__ rho_data,
    int nparticles,
    int nsites,
    int nhalo,
    int irc) {

  int p = blockIdx.x * blockDim.x + threadIdx.x;
  if (p >= nparticles) return;

  double q_p = particle_q[p];
  double r_p[3];
  r_p[0] = particle_r[3 * p + 0] - (double)d_noffset[0];
  r_p[1] = particle_r[3 * p + 1] - (double)d_noffset[1];
  r_p[2] = particle_r[3 * p + 2] - (double)d_noffset[2];

  int i0 = (int)floor(r_p[0]);
  int j0 = (int)floor(r_p[1]);
  int k0 = (int)floor(r_p[2]);

  double E_real[3] = { 0.0, 0.0, 0.0 };

  /* Real-space from lattice nodes */
  int str_z = 1;
  int str_y = (d_nlocal[2] + 2 * nhalo);
  int str_x = str_y * (d_nlocal[1] + 2 * nhalo);

  for (int di = -irc; di <= irc + 1; di++) {
    for (int dj = -irc; dj <= irc + 1; dj++) {
      for (int dk = -irc; dk <= irc + 1; dk++) {
        int ni = i0 + di;
        int nj = j0 + dj;
        int nk_idx = k0 + dk;

        if (ni < 1 || ni > d_nlocal[0]) continue;
        if (nj < 1 || nj > d_nlocal[1]) continue;
        if (nk_idx < 1 || nk_idx > d_nlocal[2]) continue;

        int ludwig_idx = str_x * (nhalo + ni - 1) + str_y * (nhalo + nj - 1) + str_z * (nhalo + nk_idx - 1);

        double rho0 = rho_data[nsites * 0 + ludwig_idx];
        double rho1 = rho_data[nsites * 1 + ludwig_idx];
        double q_node = rho0 - rho1;

        double dr[3] = { r_p[0] - (double)ni, r_p[1] - (double)nj, r_p[2] - (double)nk_idx };
        double dist = sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]);

        if (dist < d_ewald_rc && dist > 1.0e-10) {
          double ar = d_alpha * dist;
          double r_reg = sqrt(dist * dist + d_eps_reg * d_eps_reg);
          double r_inv = 1.0 / r_reg;
          /* double r_inv = 1.0/dist; */
          /* double term1 = erfc(ar)/(dist*dist); */
          double term1 = erfc(ar) / (r_reg * r_reg);
          double term2 = (2.0 * d_alpha * d_rpi) * exp(-ar * ar) * r_inv;
          double E_mag = d_beta * d_eunit * (term1 + term2) / (4.0 * M_PI * d_epsilon);

          E_real[0] += q_node * E_mag * dr[0] * r_inv;
          E_real[1] += q_node * E_mag * dr[1] * r_inv;
          E_real[2] += q_node * E_mag * dr[2] * r_inv;
        }
      }
    }
  }

  /* Real-space from other particles */
  for (int p2 = 0; p2 < nparticles; p2++) {
    if (p2 == p) continue;

    double q_p2 = particle_q[p2];
    if (fabs(q_p2) < 1.0e-14) continue;

    double r_p2[3];
    r_p2[0] = particle_r[3 * p2 + 0] - (double)d_noffset[0];
    r_p2[1] = particle_r[3 * p2 + 1] - (double)d_noffset[1];
    r_p2[2] = particle_r[3 * p2 + 2] - (double)d_noffset[2];

    double dr[3] = { r_p[0] - r_p2[0], r_p[1] - r_p2[1], r_p[2] - r_p2[2] };
    double dist = sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]);

    if (dist < d_ewald_rc && dist > 1.0e-10) {
      double ar = d_alpha * dist;
      double r_reg = sqrt(dist * dist + d_eps_reg * d_eps_reg);
      double r_inv = 1.0 / r_reg;
      /* double r_inv = 1.0/dist; */
      /* double term1 = erfc(ar)/(dist*dist); */
      double term1 = erfc(ar) / (r_reg * r_reg);
      double term2 = (2.0 * d_alpha * d_rpi) * exp(-ar * ar) * r_inv;
      double E_mag = d_beta * d_eunit * (term1 + term2) / (4.0 * M_PI * d_epsilon);

      E_real[0] += q_p2 * E_mag * dr[0] * r_inv;
      E_real[1] += q_p2 * E_mag * dr[1] * r_inv;
      E_real[2] += q_p2 * E_mag * dr[2] * r_inv;
    }
  }

  /* Add real-space to Fourier contribution already stored */
  Esub_data[3 * p + 0] += E_real[0];
  Esub_data[3 * p + 1] += E_real[1];
  Esub_data[3 * p + 2] += E_real[2];

  /* Update fex with real-space contribution */
  double kt = 1.0 / d_beta;
  fex_data[3 * p + 0] += q_p * E_real[0] * kt / d_eunit;
  fex_data[3 * p + 1] += q_p * E_real[1] * kt / d_eunit;
  fex_data[3 * p + 2] += q_p * E_real[2] * kt / d_eunit;
}

/*****************************************************************************
 *
 *  ewald_sum_force_kernel
 *
 *  Sum total force for momentum conservation check (reduction).
 *
 *****************************************************************************/

__global__ void ewald_sum_force_kernel(
    const double* __restrict__ force_data,
    double* __restrict__ F_total,
    int nsites,
    int nhalo) {

  extern __shared__ double sdata[];

  int tid = threadIdx.x;
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int ntotal = d_nlocal[0] * d_nlocal[1] * d_nlocal[2];

  /* Each thread loads one element to shared memory */
  double fx = 0.0, fy = 0.0, fz = 0.0;

  if (idx < ntotal) {
    int kc = idx % d_nlocal[2] + 1;
    int jc = (idx / d_nlocal[2]) % d_nlocal[1] + 1;
    int ic = idx / (d_nlocal[2] * d_nlocal[1]) + 1;

    int str_z = 1;
    int str_y = (d_nlocal[2] + 2 * nhalo);
    int str_x = str_y * (d_nlocal[1] + 2 * nhalo);
    int ludwig_idx = str_x * (nhalo + ic - 1) + str_y * (nhalo + jc - 1) + str_z * (nhalo + kc - 1);

    fx = force_data[3 * ludwig_idx + 0];
    fy = force_data[3 * ludwig_idx + 1];
    fz = force_data[3 * ludwig_idx + 2];
  }

  sdata[3 * tid + 0] = fx;
  sdata[3 * tid + 1] = fy;
  sdata[3 * tid + 2] = fz;
  __syncthreads();

  /* Reduction in shared memory */
  for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1) {
    if (tid < s) {
      sdata[3 * tid + 0] += sdata[3 * (tid + s) + 0];
      sdata[3 * tid + 1] += sdata[3 * (tid + s) + 1];
      sdata[3 * tid + 2] += sdata[3 * (tid + s) + 2];
    }
    __syncthreads();
  }

  /* Write result for this block */
  if (tid == 0) {
    atomicAdd(&F_total[0], sdata[0]);
    atomicAdd(&F_total[1], sdata[1]);
    atomicAdd(&F_total[2], sdata[2]);
  }
}

/*CHANGE INIT - 20260309 GPU kernels for self-field correction table */

/*****************************************************************************
 *
 *  self_fill_rho_kernel  (K1a)
 *
 *  Fill rho_d[L^3] with the Debye-Hückel charge density
 *    rho(r) = exp(-kappa * |r - r_p|) / |r - r_p|
 *  using minimum-image PBC, and accumulate per-block partial sums
 *  into partial_sum_d[] for subsequent mean subtraction.
 *
 *  Node positions follow the Ludwig convention: node (ax,ay,az) in flat
 *  index ax*L*L+ay*L+az sits at integer position (ax+1, ay+1, az+1),
 *  matching ewald_structure_factor_lattice_kernel (rx = noffset + ic,
 *  noffset=0, ic=1..L) and ewald_particle_field_real_kernel (dr = r_p - ni).
 *
 *  One thread per node.  Node index layout: idx = ax*L*L + ay*L + az.
 *
 *****************************************************************************/

__global__ void self_fill_rho_kernel(
    double* __restrict__ rho,
    double* __restrict__ partial_sum,
    double xp, double yp, double zp,
    double kappa, int L) {

  extern __shared__ double sdata[];

  int tid = threadIdx.x;
  int idx = blockIdx.x * blockDim.x + tid;
  int Ltot = L * L * L;

  double val = 0.0;
  if (idx < Ltot) {
    int az = idx % L;
    int ay = (idx / L) % L;
    int ax = idx / (L * L);

    /* Node at integer position (ax+1, ay+1, az+1) — same as Ludwig ic=ax+1 */
    double dx = (double)(ax + 1) - xp;
    double dy = (double)(ay + 1) - yp;
    double dz = (double)(az + 1) - zp;

    /* Minimum-image wrap */
    double Ld = (double)L;
    dx -= Ld * floor(dx / Ld + 0.5);
    dy -= Ld * floor(dy / Ld + 0.5);
    dz -= Ld * floor(dz / Ld + 0.5);

    double r = sqrt(dx * dx + dy * dy + dz * dz);
    val = (r > 1.0e-10) ? exp(-kappa * r) / r : 0.0;
    rho[idx] = val;
  }

  /* Shared-memory reduction for partial sum */
  sdata[tid] = val;
  __syncthreads();
  for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1) {
    if (tid < s) sdata[tid] += sdata[tid + s];
    __syncthreads();
  }
  if (tid == 0) partial_sum[blockIdx.x] = sdata[0];
}

/*****************************************************************************
 *
 *  self_subtract_mean_kernel  (K1b)
 *
 *  Subtract the mean charge density from every node so that
 *  sum rho = 0 (removes the k=0 Ewald divergence).
 *
 *****************************************************************************/

 /* norm: divide rho by norm_sum first, then subtract 1/Ltot so that
  * sum(rho) = 0 and sum(rho_before_mean_subtraction) = 1.
  * This ensures the cloud integrates to unit charge — same units as
  * q_colloid = 1 in psi->rho (density per site, summing to q_total). */
__global__ void self_subtract_mean_kernel(
    double* __restrict__ rho, double norm_sum, int Ltot) {

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < Ltot) rho[idx] = rho[idx] / norm_sum - 1.0 / (double)Ltot;
}

/*CHANGE END - 20260309 GPU kernels for self-field correction table */

/*CHANGE INIT - 20260212 New ewald_charge_sum_FFT_full_gpu using cuFFT for Fourier-space*/

/*****************************************************************************
 *
 *  ewald_copy_rho_to_fft_kernel (K1)
 *
 *  Copy rho_elec = rho[0] - rho[1] from Ludwig SOA layout (with halo)
 *  to a contiguous FFT array (no halo, row-major).
 *
 *****************************************************************************/

__global__ void ewald_copy_rho_to_fft_kernel(
    double* __restrict__ rho_real,
    const double* __restrict__ rho_data,
    int nsites,
    int nhalo,
    int nx, int ny, int nz) {

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int n_total = nx * ny * nz;
  if (idx >= n_total) return;

  /* FFT row-major: idx = ix*(ny*nz) + iy*nz + iz */
  int iz = idx % nz;
  int iy = (idx / nz) % ny;
  int ix = idx / (nz * ny);

  /* Ludwig index (Z fastest, 1-based local coords) */
  int str_z = 1;
  int str_y = (nz + 2 * nhalo);
  int str_x = str_y * (ny + 2 * nhalo);
  int ic = ix + 1;
  int jc = iy + 1;
  int kc = iz + 1;
  int ludwig_idx = str_x * (nhalo + ic - 1) + str_y * (nhalo + jc - 1)
    + str_z * (nhalo + kc - 1);

  /* rho_elec = rho[species 0] - rho[species 1] (SOA layout) */
  double rho0 = rho_data[nsites * 0 + ludwig_idx];
  double rho1 = rho_data[nsites * 1 + ludwig_idx];

  rho_real[idx] = rho0 - rho1;
}

/*****************************************************************************
 *
 *  ewald_green_multiply_potential_kernel (K2)
 *
 *  Multiply rho_hat_total(k) by G_ewald(k) to get psi_hat(k).
 *  Simple element-wise multiplication in Fourier space.
 *
 *****************************************************************************/

__global__ void ewald_green_multiply_potential_kernel(
    cufftDoubleComplex* __restrict__ psi_hat,
    const cufftDoubleComplex* __restrict__ rho_hat,
    const double* __restrict__ G_ewald,
    int n_complex) {

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx >= n_complex) return;

  double G = G_ewald[idx];
  psi_hat[idx].x = rho_hat[idx].x * G;
  psi_hat[idx].y = rho_hat[idx].y * G;
}

/*****************************************************************************
 *
 *  ewald_efield_multiply_kernel (K3)
 *
 *  Compute E_hat_alpha(k) = +i * k_alpha * psi_hat(k) for one component.
 *  Convention: F = q * grad(psi), matching the direct Ewald kernel.
 *  In cuFFT convention (rho_hat = conj(S)), +i gives the correct sign.
 *
 *  +i * k_alpha * (a + bi) = -k_alpha*b + i*k_alpha*a
 *
 *****************************************************************************/

__global__ void ewald_efield_multiply_kernel(
    cufftDoubleComplex* __restrict__ E_hat,
    const cufftDoubleComplex* __restrict__ psi_hat,
    int component,
    double kx_factor, double ky_factor, double kz_factor,
    int nx, int ny, int nz_complex) {

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int n_total = nx * ny * nz_complex;
  if (idx >= n_total) return;

  int iz = idx % nz_complex;
  int iy = (idx / nz_complex) % ny;
  int ix = idx / (nz_complex * ny);

  double k_alpha;
  if (component == 0) {
    int kx_idx = (ix <= nx / 2) ? ix : ix - nx;
    k_alpha = kx_factor * kx_idx;
  }
  else if (component == 1) {
    int ky_idx = (iy <= ny / 2) ? iy : iy - ny;
    k_alpha = ky_factor * ky_idx;
  }
  else {
    k_alpha = kz_factor * iz;
  }

  /* E_hat = +i * k_alpha * psi_hat = (-k_alpha*b, +k_alpha*a) */
  double a = psi_hat[idx].x;
  double b = psi_hat[idx].y;
  E_hat[idx].x = -k_alpha * b;
  E_hat[idx].y = k_alpha * a;
}

/*****************************************************************************
 *
 *  ewald_copy_psi_from_fft_kernel (K4)
 *
 *  Copy inverse FFT result to Ludwig layout, normalize by 1/N,
 *  add dipole correction. Store as beta * eunit * phi.
 *
 *****************************************************************************/

__global__ void ewald_copy_psi_from_fft_kernel(
    double* __restrict__ psi_data,
    const double* __restrict__ psi_real,
    int nx, int ny, int nz,
    int nhalo,
    double norm) {

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int n_total = nx * ny * nz;
  if (idx >= n_total) return;

  int iz = idx % nz;
  int iy = (idx / nz) % ny;
  int ix = idx / (nz * ny);

  int str_z = 1;
  int str_y = (nz + 2 * nhalo);
  int str_x = str_y * (ny + 2 * nhalo);
  int ic = ix + 1;
  int jc = iy + 1;
  int kc = iz + 1;
  int ludwig_idx = str_x * (nhalo + ic - 1) + str_y * (nhalo + jc - 1)
    + str_z * (nhalo + kc - 1);

  /* Global position for dipole correction */
  double rx = (double)(d_noffset[0] + ic);
  double ry = (double)(d_noffset[1] + jc);
  double rz = (double)(d_noffset[2] + kc);

  double phi_fourier = psi_real[idx] * norm;

  /* Dipole correction: phi_dipole = dipole_prefactor * (M . r) */
  double phi_dipole = d_dipole_prefactor
    * (d_M_dipole[0] * rx + d_M_dipole[1] * ry + d_M_dipole[2] * rz);

  psi_data[ludwig_idx] = d_beta * d_eunit * (phi_fourier + phi_dipole);
}

/*****************************************************************************
 *
 *  ewald_copy_force_from_fft_kernel (K5)
 *
 *  After 3 IFFTs for Ex, Ey, Ez: multiply by q(r), add dipole correction,
 *  store to force array.
 *
 *****************************************************************************/

__global__ void ewald_copy_force_from_fft_kernel(
    double* __restrict__ force_data,
    const double* __restrict__ Ex_real,
    const double* __restrict__ Ey_real,
    const double* __restrict__ Ez_real,
    const double* __restrict__ rho_data,
    int nsites,
    int nx, int ny, int nz,
    int nhalo,
    double norm) {

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int n_total = nx * ny * nz;
  if (idx >= n_total) return;

  int iz = idx % nz;
  int iy = (idx / nz) % ny;
  int ix = idx / (nz * ny);

  int str_z = 1;
  int str_y = (nz + 2 * nhalo);
  int str_x = str_y * (ny + 2 * nhalo);
  int ic = ix + 1;
  int jc = iy + 1;
  int kc = iz + 1;
  int ludwig_idx = str_x * (nhalo + ic - 1) + str_y * (nhalo + jc - 1)
    + str_z * (nhalo + kc - 1);

  /* Get charge at this node */
  double rho0 = rho_data[nsites * 0 + ludwig_idx];
  double rho1 = rho_data[nsites * 1 + ludwig_idx];
  double q = rho0 - rho1;

  /* E field from IFFT (normalized) */
  double Ex = Ex_real[idx] * norm;
  double Ey = Ey_real[idx] * norm;
  double Ez = Ez_real[idx] * norm;

  /* Force: F = q * (E_fourier + E_dipole) */
  force_data[3 * ludwig_idx + 0] = q * (Ex + d_E_dipole[0]);
  force_data[3 * ludwig_idx + 1] = q * (Ey + d_E_dipole[1]);
  force_data[3 * ludwig_idx + 2] = q * (Ez + d_E_dipole[2]);
}

/*****************************************************************************
 *
 *  ewald_particle_sfactor_fft_kernel (K6)
 *
 *  Add particle contributions to the cuFFT complex array.
 *  Each thread handles one particle, loops over all (ix,iy,iz).
 *  Uses cuFFT sign convention: exp(-ik.r).
 *
 *****************************************************************************/

__global__ void ewald_particle_sfactor_fft_kernel(
    cufftDoubleComplex* __restrict__ rho_hat,
    const double* __restrict__ particle_r,
    const double* __restrict__ particle_q,
    int nparticles,
    double kx_factor, double ky_factor, double kz_factor,
    int nx, int ny, int nz_complex) {

  int p = blockIdx.x * blockDim.x + threadIdx.x;
  if (p >= nparticles) return;

  double q = particle_q[p];
  if (fabs(q) < 1.0e-14) return;

  double rx = particle_r[3 * p + 0];
  double ry = particle_r[3 * p + 1];
  double rz = particle_r[3 * p + 2];

  for (int ix = 0; ix < nx; ix++) {
    int kx_idx = (ix <= nx / 2) ? ix : ix - nx;
    double kx = kx_factor * kx_idx;
    for (int iy = 0; iy < ny; iy++) {
      int ky_idx = (iy <= ny / 2) ? iy : iy - ny;
      double ky = ky_factor * ky_idx;
      for (int iz = 0; iz < nz_complex; iz++) {
        double kz = kz_factor * iz;
        int idx = ix * ny * nz_complex + iy * nz_complex + iz;

        double kr = kx * rx + ky * ry + kz * rz;

        /* cuFFT convention: exp(-ik.r) = cos(kr) - i*sin(kr) */
        atomicAdd(&rho_hat[idx].x, q * cos(kr));
        atomicAdd(&rho_hat[idx].y, -q * sin(kr));
      }
    }
  }
}

/*****************************************************************************
 *
 *  ewald_particle_field_fft_kernel (K7)
 *
 *  Compute Fourier-space electric field at particle positions using
 *  the cuFFT complex array. Each thread handles one particle.
 *  Handles Hermitian symmetry: factor 2 for 0 < iz < nz/2.
 *
 *****************************************************************************/

__global__ void ewald_particle_field_fft_kernel(
    double* __restrict__ Esub_data,
    double* __restrict__ fex_data,
    const double* __restrict__ particle_r,
    const double* __restrict__ particle_q,
    const cufftDoubleComplex* __restrict__ rho_hat_total,
    const double* __restrict__ G_ewald,
    int nparticles,
    double kx_factor, double ky_factor, double kz_factor,
    int nx, int ny, int nz_complex,
    double norm) {

  int p = blockIdx.x * blockDim.x + threadIdx.x;
  if (p >= nparticles) return;

  double q_p = particle_q[p];
  double rx = particle_r[3 * p + 0];
  double ry = particle_r[3 * p + 1];
  double rz = particle_r[3 * p + 2];

  double E[3] = { 0.0, 0.0, 0.0 };

  for (int ix = 0; ix < nx; ix++) {
    int kx_idx = (ix <= nx / 2) ? ix : ix - nx;
    double kx = kx_factor * kx_idx;
    for (int iy = 0; iy < ny; iy++) {
      int ky_idx = (iy <= ny / 2) ? iy : iy - ny;
      double ky = ky_factor * ky_idx;
      for (int iz = 0; iz < nz_complex; iz++) {
        double kz = kz_factor * iz;
        int idx = ix * ny * nz_complex + iy * nz_complex + iz;

        double G = G_ewald[idx];
        double a = rho_hat_total[idx].x;
        double b = rho_hat_total[idx].y;

        /* psi_hat = G * rho_hat */
        double psi_re = G * a;
        double psi_im = G * b;

        /* E_hat_alpha = +i * k_alpha * psi_hat  (convention: F = q * grad psi)
         *             = (-k_alpha*psi_im, +k_alpha*psi_re)
         *
         * E_alpha(r_p) = Re[E_hat * exp(+ik.r_p)] * herm_factor * norm
         *   Re[(-k_a*psi_im + i*k_a*psi_re)(cos(kr) + i*sin(kr))]
         *   = -k_a*(psi_im*cos(kr) + psi_re*sin(kr))
         */
        double kr = kx * rx + ky * ry + kz * rz;
        double cos_kr = cos(kr);
        double sin_kr = sin(kr);

        /* Hermitian symmetry: kz=0 and kz=nz/2 count once, rest twice */
        double herm_factor = (iz > 0 && iz < nz_complex - 1) ? 2.0 : 1.0;

        double common = herm_factor * norm;
        double E_contrib = -(psi_im * cos_kr + psi_re * sin_kr);

        E[0] += common * kx * E_contrib;
        E[1] += common * ky * E_contrib;
        E[2] += common * kz * E_contrib;
      }
    }
  }

  /* Add dipole correction */
  E[0] += d_E_dipole[0];
  E[1] += d_E_dipole[1];
  E[2] += d_E_dipole[2];

  /* Store Esub = beta * eunit * E_physical */
  Esub_data[3 * p + 0] = d_beta * d_eunit * E[0];
  Esub_data[3 * p + 1] = d_beta * d_eunit * E[1];
  Esub_data[3 * p + 2] = d_beta * d_eunit * E[2];

  /* fex = q * Esub * kt / eunit = q * E_physical */
  double kt = 1.0 / d_beta;
  fex_data[3 * p + 0] = q_p * Esub_data[3 * p + 0] * kt / d_eunit;
  fex_data[3 * p + 1] = q_p * Esub_data[3 * p + 1] * kt / d_eunit;
  fex_data[3 * p + 2] = q_p * Esub_data[3 * p + 2] * kt / d_eunit;
}

/*****************************************************************************
 *
 *  ewald_fft_fill_efield_fourier_kernel (K8)
 *
 *  Fill the Fourier-space electric field array from the IFFT results.
 *  The IFFT of (+ik_alpha * psi_hat) gives +d_alpha(phi) = -E_physical.
 *  So E_physical_alpha = -(Ex_real * norm).
 *  Also adds dipole correction to give total Fourier efield.
 *
 *****************************************************************************/

__global__ void ewald_fft_fill_efield_fourier_kernel(
    double* __restrict__ efield_data,
    const double* __restrict__ Ex_real,
    const double* __restrict__ Ey_real,
    const double* __restrict__ Ez_real,
    int nx, int ny, int nz,
    int nhalo,
    double norm) {

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int n_total = nx * ny * nz;
  if (idx >= n_total) return;

  int iz = idx % nz;
  int iy = (idx / nz) % ny;
  int ix = idx / (nz * ny);

  int str_z = 1;
  int str_y = (nz + 2 * nhalo);
  int str_x = str_y * (ny + 2 * nhalo);
  int ic = ix + 1;
  int jc = iy + 1;
  int kc = iz + 1;
  int ludwig_idx = str_x * (nhalo + ic - 1) + str_y * (nhalo + jc - 1)
    + str_z * (nhalo + kc - 1);

  /* E_physical = -grad(phi); IFFT[+ik*psi_hat] gives +grad(phi), so negate */
  double Ex = -(Ex_real[idx] * norm);
  double Ey = -(Ey_real[idx] * norm);
  double Ez = -(Ez_real[idx] * norm);

  /* Add dipole correction */
  Ex += d_E_dipole[0];
  Ey += d_E_dipole[1];
  Ez += d_E_dipole[2];

  efield_data[3 * ludwig_idx + 0] = Ex;
  efield_data[3 * ludwig_idx + 1] = Ey;
  efield_data[3 * ludwig_idx + 2] = Ez;
}

/*CHANGE END - 20260212*/

/*****************************************************************************
 *
 *  ewald_charge_sum_full_gpu
 *
 *  GPU-accelerated version of ewald_charge_sum_full.
 *  self_table: if non-NULL, the self-field correction E_self is subtracted
 *              from pc->Esub and pc->fex is recomputed accordingly.
 *
 *****************************************************************************/

int ewald_charge_sum_full_gpu(ewald_charge_t* ewald, FILE* fp,
                              ewald_self_table_t* self_table) {

  int nlocal[3], noffset[3];
  double ltot[3];
  double fkx, fky, fkz;
  double r4alpha_sq, b0;
  int irc;
  PI_DOUBLE(pi);

  if (ewald == NULL) return 0;
  assert(fp);

  TIMER_start(TIMER_EWALD_TOTAL);

  cs_nlocal(ewald->cs, nlocal);
  cs_nlocal_offset(ewald->cs, noffset);
  cs_ltot(ewald->cs, ltot);

  int nhalo;
  cs_nhalo(ewald->cs, &nhalo);

  irc = (int)ceil(ewald_rc_);

  fkx = 2.0 * pi / ltot[X];
  fky = 2.0 * pi / ltot[Y];
  fkz = 2.0 * pi / ltot[Z];
  r4alpha_sq = 1.0 / (4.0 * alpha_ * alpha_);
  b0 = 1.0 / (ltot[X] * ltot[Y] * ltot[Z] * epsilon_);

  int ntotal_nodes = nlocal[X] * nlocal[Y] * nlocal[Z];
  int nsites = ewald->psi->nsites;

  pe_info(ewald->pe, "\n");
  pe_info(ewald->pe, "Ewald charge sum full (GPU): computing potential, forces, and fields...\n");
  pe_info(ewald->pe, "  Local nodes: %d x %d x %d = %d\n", nlocal[X], nlocal[Y], nlocal[Z], ntotal_nodes);
  pe_info(ewald->pe, "  Fourier terms: %d\n", nktot_);
  pe_info(ewald->pe, "  Real-space cutoff: %.2f (irc=%d)\n", ewald_rc_, irc);

  /* Copy constants to device */
  cudaMemcpyToSymbol(d_alpha, &alpha_, sizeof(double));
  cudaMemcpyToSymbol(d_eps_reg, &eps_reg_, sizeof(double));
  cudaMemcpyToSymbol(d_epsilon, &epsilon_, sizeof(double));
  cudaMemcpyToSymbol(d_beta, &beta_, sizeof(double));
  cudaMemcpyToSymbol(d_eunit, &eunit_, sizeof(double));
  cudaMemcpyToSymbol(d_rpi, &rpi_, sizeof(double));
  cudaMemcpyToSymbol(d_ewald_rc, &ewald_rc_, sizeof(double));
  cudaMemcpyToSymbol(d_fkx, &fkx, sizeof(double));
  cudaMemcpyToSymbol(d_fky, &fky, sizeof(double));
  cudaMemcpyToSymbol(d_fkz, &fkz, sizeof(double));
  cudaMemcpyToSymbol(d_r4alpha_sq, &r4alpha_sq, sizeof(double));
  cudaMemcpyToSymbol(d_b0, &b0, sizeof(double));
  cudaMemcpyToSymbol(d_kmax, &kmax_, sizeof(double));
  cudaMemcpyToSymbol(d_nk, nk_, 3 * sizeof(int));
  cudaMemcpyToSymbol(d_nlocal, nlocal, 3 * sizeof(int));
  cudaMemcpyToSymbol(d_noffset, noffset, 3 * sizeof(int));
  cudaMemcpyToSymbol(d_ltot, ltot, 3 * sizeof(double));

  /* ========================================================================
   * Precompute k-vectors and Green function on CPU, then copy to GPU
   * ======================================================================== */

  pe_info(ewald->pe, "  [1/6] Precomputing k-vectors and Green function...\n");

  double* kvec_h = (double*)malloc(3 * nktot_ * sizeof(double));
  double* Gk_h = (double*)malloc(nktot_ * sizeof(double));
  int* kz_arr_h = (int*)malloc(nktot_ * sizeof(int));

  // Se itera sobre todos los vectores k = (kx,ky,kz) en el espacio recíproco
  // que cumplen 0 < |k|² ≤ k_max
  int kn = 0;
  for (int kz = 0; kz <= nk_[Z]; kz++) {
    for (int ky = -nk_[Y]; ky <= nk_[Y]; ky++) {
      for (int kx = -nk_[X]; kx <= nk_[X]; kx++) {
        double k[3], ksq;
        k[X] = fkx * kx;
        k[Y] = fky * ky;
        k[Z] = fkz * kz;
        ksq = k[X] * k[X] + k[Y] * k[Y] + k[Z] * k[Z];
        if (ksq <= 0.0 || ksq > kmax_) continue;

        kvec_h[3 * kn + X] = k[X];
        kvec_h[3 * kn + Y] = k[Y];
        kvec_h[3 * kn + Z] = k[Z];
        Gk_h[kn] = b0 * exp(-r4alpha_sq * ksq) / ksq;
        kz_arr_h[kn] = kz;
        kn++;
      }
    }
  }
  int nk_actual = kn;

  /* Allocate device arrays */
  double* kvec_d, * Gk_d, * Sk_sin_d, * Sk_cos_d;
  int* kz_arr_d;

  cudaMalloc(&kvec_d, 3 * nk_actual * sizeof(double));
  cudaMalloc(&Gk_d, nk_actual * sizeof(double));
  cudaMalloc(&kz_arr_d, nk_actual * sizeof(int));
  cudaMalloc(&Sk_sin_d, nk_actual * sizeof(double));
  cudaMalloc(&Sk_cos_d, nk_actual * sizeof(double));

  cudaMemcpy(kvec_d, kvec_h, 3 * nk_actual * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(Gk_d, Gk_h, nk_actual * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(kz_arr_d, kz_arr_h, nk_actual * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemset(Sk_sin_d, 0, nk_actual * sizeof(double));
  cudaMemset(Sk_cos_d, 0, nk_actual * sizeof(double));

  /* Get device pointers for rho and psi fields */
  double* rho_data_d = NULL;
  double* psi_data_d = NULL;

  size_t data_offset = offsetof(field_t, data);
  cudaMemcpy(&rho_data_d, (char*)(ewald->psi->rho->target) + data_offset,
             sizeof(double*), cudaMemcpyDeviceToHost);
  cudaMemcpy(&psi_data_d, (char*)(ewald->psi->psi->target) + data_offset,
             sizeof(double*), cudaMemcpyDeviceToHost);

  /* Ensure data is on device */
  field_memcpy(ewald->psi->rho, tdpMemcpyHostToDevice);

  /* ========================================================================
   * Collect particle data and copy to GPU
   * ======================================================================== */

  int nparticles = 0;
  double* particle_r_h = NULL, * particle_q_h = NULL;
  double* particle_r_d = NULL, * particle_q_d = NULL;
  double* Esub_d = NULL, * fex_d = NULL;

  if ((ewald->sources & EWALD_SOURCE_COLLOIDS) && ewald->cinfo) {
    /* Count particles */
    colloid_t* pc;
    colloids_info_local_head(ewald->cinfo, &pc);
    for (; pc; pc = pc->nextlocal) nparticles++;

    if (nparticles > 0) {
      particle_r_h = (double*)malloc(3 * nparticles * sizeof(double));
      particle_q_h = (double*)malloc(nparticles * sizeof(double));

      int p = 0;
      colloids_info_local_head(ewald->cinfo, &pc);
      for (; pc; pc = pc->nextlocal) {
        particle_r_h[3 * p + 0] = pc->s.r[X];
        particle_r_h[3 * p + 1] = pc->s.r[Y];
        particle_r_h[3 * p + 2] = pc->s.r[Z];
        particle_q_h[p] = pc->s.q0 - pc->s.q1;
        p++;
      }

      cudaMalloc(&particle_r_d, 3 * nparticles * sizeof(double));
      cudaMalloc(&particle_q_d, nparticles * sizeof(double));
      cudaMalloc(&Esub_d, 3 * nparticles * sizeof(double));
      cudaMalloc(&fex_d, 3 * nparticles * sizeof(double));

      cudaMemcpy(particle_r_d, particle_r_h, 3 * nparticles * sizeof(double), cudaMemcpyHostToDevice);
      cudaMemcpy(particle_q_d, particle_q_h, nparticles * sizeof(double), cudaMemcpyHostToDevice);
    }
  }

  pe_info(ewald->pe, "  Particles: %d\n", nparticles);

  /* ========================================================================
   * PART 1: Compute structure factors on GPU
   * ======================================================================== */

  pe_info(ewald->pe, "  [2/6] Computing structure factors S(k), C(k) on GPU...\n");

  int threads = 256;
  int blocks = (ntotal_nodes + threads - 1) / threads;

  /* Structure factor contribution from lattice */
  if (ewald->sources & EWALD_SOURCE_LATTICE) {
    ewald_structure_factor_lattice_kernel << <blocks, threads >> > (
        rho_data_d, Sk_sin_d, Sk_cos_d, kvec_d, nk_actual, nsites, nhalo);
    cudaDeviceSynchronize();
  }

  /* Structure factor contribution from particles */
  if ((ewald->sources & EWALD_SOURCE_COLLOIDS) && nparticles > 0) {
    int p_blocks = (nparticles + threads - 1) / threads;
    ewald_structure_factor_particle_kernel << <p_blocks, threads >> > (
        particle_r_d, particle_q_d, Sk_sin_d, Sk_cos_d, kvec_d, nparticles, nk_actual);
    cudaDeviceSynchronize();
  }

  /* MPI reduce structure factors */
  pe_info(ewald->pe, "  [3/6] MPI reducing structure factors...\n");

  double* Sk_sin_h = (double*)malloc(nk_actual * sizeof(double));
  double* Sk_cos_h = (double*)malloc(nk_actual * sizeof(double));
  cudaMemcpy(Sk_sin_h, Sk_sin_d, nk_actual * sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(Sk_cos_h, Sk_cos_d, nk_actual * sizeof(double), cudaMemcpyDeviceToHost);

  MPI_Comm comm;
  cs_cart_comm(ewald->cs, &comm);
  MPI_Allreduce(MPI_IN_PLACE, Sk_sin_h, nk_actual, MPI_DOUBLE, MPI_SUM, comm);
  MPI_Allreduce(MPI_IN_PLACE, Sk_cos_h, nk_actual, MPI_DOUBLE, MPI_SUM, comm);

  cudaMemcpy(Sk_sin_d, Sk_sin_h, nk_actual * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(Sk_cos_d, Sk_cos_h, nk_actual * sizeof(double), cudaMemcpyHostToDevice);

  /* CHANGE INIT - SyncSinCos_GPU - Copy Sk to global sinx_/cosx_ so ewald_efield_at_r (CPU)
   * can use the Fourier contribution. In the GPU path these are never populated otherwise. */
  assert(nk_actual == nktot_);
  memcpy(sinx_, Sk_sin_h, nk_actual * sizeof(double));
  memcpy(cosx_, Sk_cos_h, nk_actual * sizeof(double));
  /* CHANGE END - SyncSinCos_GPU */

  /* ========================================================================
   * Compute dipole moment M = sum_i q_i * r_i for dipole correction
   * The dipole correction field is: E_dipole = -(4*pi / (1 + 2*epsilon_prime)*V) * M
   * ======================================================================== */

  kahan_t M_dipole[3] = { kahan_zero(), kahan_zero(), kahan_zero() };
  kahan_t Q_total_k = kahan_zero();

  /* Dipole from colloids */
  if (nparticles > 0) {
    for (int p = 0; p < nparticles; p++) {
      double q = particle_q_h[p];
      kahan_add_double(&M_dipole[X], q * particle_r_h[3 * p + 0]);
      kahan_add_double(&M_dipole[Y], q * particle_r_h[3 * p + 1]);
      kahan_add_double(&M_dipole[Z], q * particle_r_h[3 * p + 2]);
      kahan_add_double(&Q_total_k, q);
    }
  }

  /* Dipole from lattice nodes (need to read rho from device) */
  if ((ewald->sources & EWALD_SOURCE_LATTICE) && ewald->psi) {
    /* Sync rho to host for dipole calculation */
    field_memcpy(ewald->psi->rho, tdpMemcpyDeviceToHost);

    /* Accumulate each species separately to avoid catastrophic cancellation
     * when rho[0] and rho[1] are nearly equal (q_node = v0*rho0 + v1*rho1 ~ 0) */
    int nk_ewald;
    psi_nk(ewald->psi, &nk_ewald);
    kahan_t* Q_species_k = (kahan_t*)calloc(nk_ewald, sizeof(kahan_t));
    kahan_t* Mx_species_k = (kahan_t*)calloc(nk_ewald, sizeof(kahan_t));
    kahan_t* My_species_k = (kahan_t*)calloc(nk_ewald, sizeof(kahan_t));
    kahan_t* Mz_species_k = (kahan_t*)calloc(nk_ewald, sizeof(kahan_t));
    for (int s = 0; s < nk_ewald; s++) {
      Q_species_k[s] = kahan_zero();
      Mx_species_k[s] = kahan_zero();
      My_species_k[s] = kahan_zero();
      Mz_species_k[s] = kahan_zero();
    }

    for (int ic = 1; ic <= nlocal[X]; ic++) {
      for (int jc = 1; jc <= nlocal[Y]; jc++) {
        for (int kc = 1; kc <= nlocal[Z]; kc++) {
          int index = cs_index(ewald->cs, ic, jc, kc);
          double rx = (double)(noffset[X] + ic);
          double ry = (double)(noffset[Y] + jc);
          double rz = (double)(noffset[Z] + kc);
          for (int s = 0; s < nk_ewald; s++) {
            int val;
            double rho_s;
            psi_valency(ewald->psi, s, &val);
            psi_rho(ewald->psi, index, s, &rho_s);
            double q_s = val * rho_s;
            kahan_add_double(&Q_species_k[s], q_s);
            kahan_add_double(&Mx_species_k[s], q_s * rx);
            kahan_add_double(&My_species_k[s], q_s * ry);
            kahan_add_double(&Mz_species_k[s], q_s * rz);
          }
        }
      }
    }

    for (int s = 0; s < nk_ewald; s++) {
      kahan_add_double(&Q_total_k, kahan_sum(&Q_species_k[s]));
      kahan_add_double(&M_dipole[X], kahan_sum(&Mx_species_k[s]));
      kahan_add_double(&M_dipole[Y], kahan_sum(&My_species_k[s]));
      kahan_add_double(&M_dipole[Z], kahan_sum(&Mz_species_k[s]));
    }
    free(Q_species_k); free(Mx_species_k); free(My_species_k); free(Mz_species_k);
  }

  /* MPI reduce dipole moment */
  {
    kahan_t M_reduce[4] = { M_dipole[X], M_dipole[Y], M_dipole[Z], Q_total_k };
    MPI_Datatype kahan_dt;
    MPI_Op kahan_op;
    kahan_mpi_datatype(&kahan_dt);
    kahan_mpi_op_sum(&kahan_op);
    MPI_Allreduce(MPI_IN_PLACE, M_reduce, 4, kahan_dt, kahan_op, comm);
    MPI_Type_free(&kahan_dt);
    MPI_Op_free(&kahan_op);
    M_dipole[X] = M_reduce[0];
    M_dipole[Y] = M_reduce[1];
    M_dipole[Z] = M_reduce[2];
    Q_total_k = M_reduce[3];
  }

  double M[3] = { kahan_sum(&M_dipole[X]), kahan_sum(&M_dipole[Y]), kahan_sum(&M_dipole[Z]) };
  double Q_total = kahan_sum(&Q_total_k);
  double V = ltot[X] * ltot[Y] * ltot[Z];

  /* Dipole correction using Deserno & Holm convention (J. Chem. Phys. 1998)
   * Energy:    E^(d) = 2*pi / ((1 + 2*epsilon')*V) * |M|^2
   * Force:     F^(d)_i = -4*pi*q_i / ((1 + 2*epsilon')*V) * M
   * Field:     E_dipole = -4*pi / ((1 + 2*epsilon')*V) * M
   *
   * epsilon' = 1     : vacuum boundary conditions
   * epsilon' = infty : metallic (tinfoil) boundary conditions (correction vanishes)
   */
  double dipole_prefactor = 4.0 * pi / ((1.0 + 2.0 * ewald->epsilon_prime) * V);
  double E_dipole[3] = { -dipole_prefactor * M[X],
                        -dipole_prefactor * M[Y],
                        -dipole_prefactor * M[Z] };

  pe_info(ewald->pe, "  Dipole moment M = (%14.7e, %14.7e, %14.7e)\n", M[X], M[Y], M[Z]);
  pe_info(ewald->pe, "  Total charge Q = %14.7e\n", Q_total);
  pe_info(ewald->pe, "  Dipole epsilon' = %14.7e (Deserno & Holm convention)\n", ewald->epsilon_prime);
  pe_info(ewald->pe, "  Dipole field E_dip = (%14.7e, %14.7e, %14.7e)\n", E_dipole[X], E_dipole[Y], E_dipole[Z]);

  /* Copy dipole field to device constants */
  cudaMemcpyToSymbol(d_E_dipole, E_dipole, 3 * sizeof(double));
  cudaMemcpyToSymbol(d_M_dipole, M, 3 * sizeof(double));
  double dipole_pref_dev = dipole_prefactor;
  cudaMemcpyToSymbol(d_dipole_prefactor, &dipole_pref_dev, sizeof(double));

  /* ========================================================================
   * PART 2: Compute potential on GPU
   * ======================================================================== */

  pe_info(ewald->pe, "  [4/6] Computing potential on lattice nodes (GPU)...\n");

  /* Fourier part */
  ewald_potential_fourier_kernel << <blocks, threads >> > (
      psi_data_d, Sk_sin_d, Sk_cos_d, kvec_d, Gk_d, kz_arr_d,
      nk_actual, nsites, nhalo);
  cudaDeviceSynchronize();

  /* Real-space from particles */
  if (nparticles > 0) {
    ewald_potential_real_particle_kernel << <blocks, threads >> > (
        psi_data_d, particle_r_d, particle_q_d, nparticles, nsites, nhalo, irc);
    cudaDeviceSynchronize();
  }

  /* Real-space from lattice nodes */
  ewald_potential_real_lattice_kernel << <blocks, threads >> > (
      psi_data_d, rho_data_d, nsites, nhalo, irc);
  cudaDeviceSynchronize();

  /* Sync psi back to host */
  field_memcpy(ewald->psi->psi, tdpMemcpyDeviceToHost);

  /*CHANGE INIT - 20251203 Electric field output */
  /* Compute diagnostic psi_fourier and psi_real on separate GPU arrays */
  if (ewald->psi && (ewald->psi->psi_fourier || ewald->psi->psi_real)) {
    double* psi_fourier_d = NULL;
    double* psi_real_d = NULL;

    /* Allocate temporary arrays (scalar fields: 1 value per site) */
    cudaMalloc(&psi_fourier_d, nsites * sizeof(double));
    cudaMemset(psi_fourier_d, 0, nsites * sizeof(double));
    cudaMalloc(&psi_real_d, nsites * sizeof(double));
    cudaMemset(psi_real_d, 0, nsites * sizeof(double));

    /* Fill Fourier component */
    ewald_potential_fourier_kernel << <blocks, threads >> > (
        psi_fourier_d, Sk_sin_d, Sk_cos_d, kvec_d, Gk_d, kz_arr_d,
        nk_actual, nsites, nhalo);
    cudaDeviceSynchronize();

    /* Fill real component (real kernel adds to array, which starts at zero) */
    if (nparticles > 0) {
      ewald_potential_real_particle_kernel << <blocks, threads >> > (
          psi_real_d, particle_r_d, particle_q_d, nparticles, nsites, nhalo, irc);
      cudaDeviceSynchronize();
    }
    ewald_potential_real_lattice_kernel << <blocks, threads >> > (
        psi_real_d, rho_data_d, nsites, nhalo, irc);
    cudaDeviceSynchronize();

    /* Copy back to host and store in psi fields */
    double* psi_fourier_h = (double*)malloc(nsites * sizeof(double));
    double* psi_real_h = (double*)malloc(nsites * sizeof(double));
    cudaMemcpy(psi_fourier_h, psi_fourier_d, nsites * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(psi_real_h, psi_real_d, nsites * sizeof(double), cudaMemcpyDeviceToHost);

    for (int ic = 1; ic <= nlocal[X]; ic++) {
      for (int jc = 1; jc <= nlocal[Y]; jc++) {
        for (int kc = 1; kc <= nlocal[Z]; kc++) {
          int index = cs_index(ewald->cs, ic, jc, kc);
          if (ewald->psi->psi_fourier) {
            field_scalar_set(ewald->psi->psi_fourier, index, psi_fourier_h[index]);
          }
          if (ewald->psi->psi_real) {
            field_scalar_set(ewald->psi->psi_real, index, psi_real_h[index]);
          }
        }
      }
    }

    free(psi_fourier_h);
    free(psi_real_h);
    cudaFree(psi_fourier_d);
    cudaFree(psi_real_d);

    /* Sync diagnostic fields to device */
    if (ewald->psi->psi_fourier) {
      field_memcpy(ewald->psi->psi_fourier, tdpMemcpyHostToDevice);
    }
    if (ewald->psi->psi_real) {
      field_memcpy(ewald->psi->psi_real, tdpMemcpyHostToDevice);
    }
  }
  /*CHANGE END - 20251203 */

  /* ========================================================================
   * PART 3: Compute force on fluid nodes
   * ======================================================================== */

  pe_info(ewald->pe, "  [5/6] Computing force on lattice nodes (GPU)...\n");

  /* Kahan accumulators for momentum conservation check (outside if block) */
  kahan_t F_fluid_total[3] = { kahan_zero(), kahan_zero(), kahan_zero() };

  if (ewald->hydro != NULL) {
    /* Sync hydro from device */
    hydro_memcpy(ewald->hydro, tdpMemcpyDeviceToHost);

    /*CHANGE INIT - 20251203 Electric field output */
    /* Allocate temporary arrays on device for electric field and force */
    double* efield_d = NULL;
    double* efield_fourier_d = NULL;
    double* efield_real_d = NULL;
    double* force_d;
    cudaMalloc(&force_d, 3 * nsites * sizeof(double));
    cudaMemset(force_d, 0, 3 * nsites * sizeof(double));

    /* Allocate efield arrays if needed */
    if (ewald->psi && ewald->psi->efield) {
      cudaMalloc(&efield_d, 3 * nsites * sizeof(double));
      cudaMemset(efield_d, 0, 3 * nsites * sizeof(double));

      /* Allocate separate arrays for Fourier and real components (for diagnostics) */
      cudaMalloc(&efield_fourier_d, 3 * nsites * sizeof(double));
      cudaMemset(efield_fourier_d, 0, 3 * nsites * sizeof(double));
      cudaMalloc(&efield_real_d, 3 * nsites * sizeof(double));
      cudaMemset(efield_real_d, 0, 3 * nsites * sizeof(double));

      /* Compute Fourier part of electric field (without charge) */
      ewald_efield_fourier_kernel << <blocks, threads >> > (
          efield_d, Sk_sin_d, Sk_cos_d, kvec_d, Gk_d, kz_arr_d,
          nk_actual, nsites, nhalo);
      cudaDeviceSynchronize();

      /* Copy Fourier part to diagnostic array before adding real part */
      cudaMemcpy(efield_fourier_d, efield_d, 3 * nsites * sizeof(double), cudaMemcpyDeviceToDevice);

      /* Real-space contribution from particles to electric field */
      if (nparticles > 0) {
        ewald_efield_real_particle_kernel << <blocks, threads >> > (
            efield_d, particle_r_d, particle_q_d,
            nparticles, nsites, nhalo, irc);
        cudaDeviceSynchronize();
        /* Also add to real-space diagnostic array */
        ewald_efield_real_particle_kernel << <blocks, threads >> > (
            efield_real_d, particle_r_d, particle_q_d,
            nparticles, nsites, nhalo, irc);
        cudaDeviceSynchronize();
      }

      /* Real-space contribution from lattice nodes to electric field */
      ewald_efield_real_lattice_kernel << <blocks, threads >> > (
          efield_d, rho_data_d, nsites, nhalo, irc);
      cudaDeviceSynchronize();
      /* Also add to real-space diagnostic array */
      ewald_efield_real_lattice_kernel << <blocks, threads >> > (
          efield_real_d, rho_data_d, nsites, nhalo, irc);
      cudaDeviceSynchronize();
    }
    /*CHANGE END - 20251203 */

    /* Fourier part of force */
    ewald_force_fourier_kernel << <blocks, threads >> > (
        force_d, rho_data_d, Sk_sin_d, Sk_cos_d, kvec_d, Gk_d, kz_arr_d,
        nk_actual, nsites, nhalo);
    cudaDeviceSynchronize();

    /* Real-space from particles */
    if (nparticles > 0) {
      ewald_force_real_particle_kernel << <blocks, threads >> > (
          force_d, rho_data_d, particle_r_d, particle_q_d,
          nparticles, nsites, nhalo, irc);
      cudaDeviceSynchronize();
    }

    /* Copy force back and add to hydro */
    double* force_h = (double*)malloc(3 * nsites * sizeof(double));
    cudaMemcpy(force_h, force_d, 3 * nsites * sizeof(double), cudaMemcpyDeviceToHost);

    /*CHANGE INIT - 20251203 Electric field output */
    /* Copy electric field back to host if computed */
    double* efield_h = NULL;
    double* efield_fourier_h = NULL;
    double* efield_real_h = NULL;
    if (efield_d != NULL) {
      efield_h = (double*)malloc(3 * nsites * sizeof(double));
      cudaMemcpy(efield_h, efield_d, 3 * nsites * sizeof(double), cudaMemcpyDeviceToHost);
      pe_info(ewald->pe, "  Electric field computed and copied to host\n");

      /* Copy diagnostic components */
      if (efield_fourier_d != NULL) {
        efield_fourier_h = (double*)malloc(3 * nsites * sizeof(double));
        cudaMemcpy(efield_fourier_h, efield_fourier_d, 3 * nsites * sizeof(double), cudaMemcpyDeviceToHost);
      }
      if (efield_real_d != NULL) {
        efield_real_h = (double*)malloc(3 * nsites * sizeof(double));
        cudaMemcpy(efield_real_h, efield_real_d, 3 * nsites * sizeof(double), cudaMemcpyDeviceToHost);
      }
    }
    /*CHANGE END - 20251203 */

    /*CHANGE INIT - 20251203 Electric field output */
    /* Debug: sum up electric field values to check they are not zero */
    double e_sum[3] = { 0.0, 0.0, 0.0 };
    int e_count = 0;
    /*CHANGE END - 20251203 */

    for (int ic = 1; ic <= nlocal[X]; ic++) {
      for (int jc = 1; jc <= nlocal[Y]; jc++) {
        for (int kc = 1; kc <= nlocal[Z]; kc++) {
          int index = cs_index(ewald->cs, ic, jc, kc);
          double f[3] = { force_h[3 * index + 0], force_h[3 * index + 1], force_h[3 * index + 2] };
          hydro_f_local_add(ewald->hydro, index, f);

          // /* CHANGE INIT - ForceFluidOutput - Store Ewald force density in psi->force_fluid for output */
          // if (ewald->psi && ewald->psi->force_fluid) {
          //   field_vector_set(ewald->psi->force_fluid, index, f);
          // }
          // /* CHANGE END - ForceFluidOutput */

          /*CHANGE INIT - 20251203 Electric field output */
          /* Store electric field computed on GPU */
          if (efield_h != NULL) {
            double e_field[3] = { efield_h[3 * index + 0], efield_h[3 * index + 1], efield_h[3 * index + 2] };
            field_vector_set(ewald->psi->efield, index, e_field);
            /* Debug: accumulate for checking */
            e_sum[X] += fabs(e_field[X]);
            e_sum[Y] += fabs(e_field[Y]);
            e_sum[Z] += fabs(e_field[Z]);
            e_count++;

            /* Store diagnostic components */
            if (efield_fourier_h != NULL && ewald->psi->efield_fourier) {
              double e_fourier[3] = { efield_fourier_h[3 * index + 0], efield_fourier_h[3 * index + 1], efield_fourier_h[3 * index + 2] };
              field_vector_set(ewald->psi->efield_fourier, index, e_fourier);
            }
            if (efield_real_h != NULL && ewald->psi->efield_real) {
              double e_real[3] = { efield_real_h[3 * index + 0], efield_real_h[3 * index + 1], efield_real_h[3 * index + 2] };
              field_vector_set(ewald->psi->efield_real, index, e_real);
            }
          }
          /*CHANGE END - 20251203 */

          kahan_add_double(&F_fluid_total[X], f[X]);
          kahan_add_double(&F_fluid_total[Y], f[Y]);
          kahan_add_double(&F_fluid_total[Z], f[Z]);
        }
      }
    }

    /*CHANGE INIT - 20251203 Electric field output */
    /* Debug: print sum of electric field */
    if (efield_h != NULL && e_count > 0) {
      pe_info(ewald->pe, "  Electric field sum |E|: (%14.7e, %14.7e, %14.7e) over %d sites\n",
              e_sum[X] / e_count, e_sum[Y] / e_count, e_sum[Z] / e_count, e_count);
    }
    /*CHANGE END - 20251203 */

    /*CHANGE INIT - 20251203 Electric field output */
    /* Free temporary arrays */
    if (efield_h != NULL) free(efield_h);
    if (efield_fourier_h != NULL) free(efield_fourier_h);
    if (efield_real_h != NULL) free(efield_real_h);
    if (efield_d != NULL) cudaFree(efield_d);
    if (efield_fourier_d != NULL) cudaFree(efield_fourier_d);
    if (efield_real_d != NULL) cudaFree(efield_real_d);

    /* Sync efield back to device if it was computed */
    if (ewald->psi && ewald->psi->efield) {
      field_memcpy(ewald->psi->efield, tdpMemcpyHostToDevice);
    }
    /* Sync diagnostic fields */
    if (ewald->psi && ewald->psi->efield_fourier) {
      field_memcpy(ewald->psi->efield_fourier, tdpMemcpyHostToDevice);
    }
    if (ewald->psi && ewald->psi->efield_real) {
      field_memcpy(ewald->psi->efield_real, tdpMemcpyHostToDevice);
    }
    /*CHANGE END - 20251203 */
    // /* CHANGE INIT - ForceFluidOutput - Sync force_fluid back to device */
    // if (ewald->psi && ewald->psi->force_fluid) {
    //   field_memcpy(ewald->psi->force_fluid, tdpMemcpyHostToDevice);
    // }
    // /* CHANGE END - ForceFluidOutput */
    free(force_h);
    cudaFree(force_d);

    /* Sync hydro back to device */
    hydro_memcpy(ewald->hydro, tdpMemcpyHostToDevice);
  }

  /* MPI reduce F_fluid_total (outside if block) */
  {
    MPI_Datatype kahan_dt;
    MPI_Op kahan_op;
    kahan_mpi_datatype(&kahan_dt);
    kahan_mpi_op_sum(&kahan_op);
    MPI_Allreduce(MPI_IN_PLACE, F_fluid_total, 3, kahan_dt, kahan_op, comm);
    MPI_Type_free(&kahan_dt);
    MPI_Op_free(&kahan_op);
  }

  /* ========================================================================
   * PART 4: Compute field and force on particles
   * ======================================================================== */

  pe_info(ewald->pe, "  [6/6] Computing field and force on particles (GPU)...\n");

  kahan_t F_particle_total[3] = { kahan_zero(), kahan_zero(), kahan_zero() };

  if (nparticles > 0) {
    int p_blocks = (nparticles + threads - 1) / threads;

    /* Fourier part */
    ewald_particle_field_fourier_kernel << <p_blocks, threads >> > (
        Esub_d, fex_d, particle_r_d, particle_q_d,
        Sk_sin_d, Sk_cos_d, kvec_d, Gk_d, kz_arr_d,
        nparticles, nk_actual);
    cudaDeviceSynchronize();

    /* Real-space from lattice and other particles */
    ewald_particle_field_real_kernel << <p_blocks, threads >> > (
        Esub_d, fex_d, particle_r_d, particle_q_d, rho_data_d,
        nparticles, nsites, nhalo, irc);
    cudaDeviceSynchronize();

    /* Copy results back to host and store in colloid structures */
    double* Esub_h = (double*)malloc(3 * nparticles * sizeof(double));
    double* fex_h = (double*)malloc(3 * nparticles * sizeof(double));
    cudaMemcpy(Esub_h, Esub_d, 3 * nparticles * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(fex_h, fex_d, 3 * nparticles * sizeof(double), cudaMemcpyDeviceToHost);

    int p = 0;
    /*CHANGE INIT - 20260310 Self-field correction */
    double kt = 1.0 / beta_;
    /*CHANGE END - 20260310*/
    colloid_t* pc;
    colloids_info_local_head(ewald->cinfo, &pc);
    for (; pc; pc = pc->nextlocal) {
      /*CHANGE INIT - 20260310 Self-field correction */
      double Esub_x = Esub_h[3 * p + 0];
      double Esub_y = Esub_h[3 * p + 1];
      double Esub_z = Esub_h[3 * p + 2];

      /* Subtract self-field artefact if table provided.
       * The table stores E_self for a +1 Debye cloud (normalised to unit charge).
       * The real counter-cloud has charge -q_colloid, so the artefact is
       *   E_artefact = -q_p * E_self_table
       * Subtracting it: Esub_corrected = Esub - E_artefact = Esub + q_p * E_self */
      if (self_table) {
        double xf = pc->s.r[X] - floor(pc->s.r[X]);
        double yf = pc->s.r[Y] - floor(pc->s.r[Y]);
        double zf = pc->s.r[Z] - floor(pc->s.r[Z]);
        double Eself[3];
        ewald_charge_self_field_interpolate(self_table, xf, yf, zf, Eself);
        double q_p = pc->s.q0 - pc->s.q1;
        Esub_x += q_p * Eself[0];
        Esub_y += q_p * Eself[1];
        Esub_z += q_p * Eself[2];

        double q_pc = pc->s.q0 - pc->s.q1;
        pc->Esub[X] = Esub_x;
        pc->Esub[Y] = Esub_y;
        pc->Esub[Z] = Esub_z;
        pc->fex[X] = q_pc * Esub_x * kt / eunit_;
        pc->fex[Y] = q_pc * Esub_y * kt / eunit_;
        pc->fex[Z] = q_pc * Esub_z * kt / eunit_;

      }
      else {
        /*CHANGE END - 20260310*/
        /* Original (no self-field correction):*/
        // pc->Esub[X] = Esub_h[3*p + 0];
        // pc->Esub[Y] = Esub_h[3*p + 1];
        // pc->Esub[Z] = Esub_h[3*p + 2];
        pc->Esub[X] = Esub_x;
        pc->Esub[Y] = Esub_y;
        pc->Esub[Z] = Esub_z;
        pc->fex[X] = fex_h[3 * p + 0];
        pc->fex[Y] = fex_h[3 * p + 1];
        pc->fex[Z] = fex_h[3 * p + 2];
      }

      kahan_add_double(&F_particle_total[X], pc->fex[X]);
      kahan_add_double(&F_particle_total[Y], pc->fex[Y]);
      kahan_add_double(&F_particle_total[Z], pc->fex[Z]);
      p++;
    }

    free(Esub_h);
    free(fex_h);
  }

  /* MPI reduce particle forces */
  {
    MPI_Datatype kahan_dt;
    MPI_Op kahan_op;
    kahan_mpi_datatype(&kahan_dt);
    kahan_mpi_op_sum(&kahan_op);
    MPI_Allreduce(MPI_IN_PLACE, F_particle_total, 3, kahan_dt, kahan_op, comm);
    MPI_Type_free(&kahan_dt);
    MPI_Op_free(&kahan_op);
  }

  double F_fluid[3] = { kahan_sum(&F_fluid_total[X]), kahan_sum(&F_fluid_total[Y]), kahan_sum(&F_fluid_total[Z]) };
  double F_particle[3] = { kahan_sum(&F_particle_total[X]), kahan_sum(&F_particle_total[Y]), kahan_sum(&F_particle_total[Z]) };
  double F_diff[3] = { F_fluid[X] + F_particle[X], F_fluid[Y] + F_particle[Y], F_fluid[Z] + F_particle[Z] };

  /* Write to file (same format as ewald_charge_sum_full) */
  fprintf(fp, "%1.15e,%1.15e, %1.15e, %1.15e,%1.15e, %1.15e, %1.15e,%1.15e, %1.15e, %1.15e,\n",
          sqrt(F_diff[X] * F_diff[X] + F_diff[Y] * F_diff[Y] + F_diff[Z] * F_diff[Z]),
          F_diff[X], F_diff[Y], F_diff[Z],
          F_fluid[X], F_fluid[Y], F_fluid[Z],
          F_particle[X], F_particle[Y], F_particle[Z]);

  /* Print momentum conservation check */
  pe_info(ewald->pe, "\n");
  pe_info(ewald->pe, "  Momentum conservation check:\n");
  pe_info(ewald->pe, "    F_fluid    = (%14.7e, %14.7e, %14.7e)\n", F_fluid[X], F_fluid[Y], F_fluid[Z]);
  pe_info(ewald->pe, "    F_particle = (%14.7e, %14.7e, %14.7e)\n", F_particle[X], F_particle[Y], F_particle[Z]);
  pe_info(ewald->pe, "    F_total    = (%14.7e, %14.7e, %14.7e)\n", F_diff[X], F_diff[Y], F_diff[Z]);
  pe_info(ewald->pe, "    |F_total|  = %14.7e\n", sqrt(F_diff[X] * F_diff[X] + F_diff[Y] * F_diff[Y] + F_diff[Z] * F_diff[Z]));

  /* External field contribution */
  ewald_charge_external_field(ewald);

  /* ========================================================================
   * Cleanup
   * ======================================================================== */

  free(kvec_h);
  free(Gk_h);
  free(kz_arr_h);
  free(Sk_sin_h);
  free(Sk_cos_h);

  cudaFree(kvec_d);
  cudaFree(Gk_d);
  cudaFree(kz_arr_d);
  cudaFree(Sk_sin_d);
  cudaFree(Sk_cos_d);

  if (nparticles > 0) {
    free(particle_r_h);
    free(particle_q_h);
    cudaFree(particle_r_d);
    cudaFree(particle_q_d);
    cudaFree(Esub_d);
    cudaFree(fex_d);
  }

  pe_info(ewald->pe, "Ewald charge sum full (GPU): complete.\n\n");

  TIMER_stop(TIMER_EWALD_TOTAL);

  return 0;
}

/*****************************************************************************
 *
 *  ewald_face_E_gpu_kernel  (replaces ewald_face_div_gpu_kernel)
 *
 *  CHANGE INIT - FaceDiv_GPU - Deduplicated face-sampling for div(E)
 *
 *  Each unique face is evaluated exactly once by one thread.
 *  Unique faces for a window of size ni×nj×nk_win:
 *    ⊥X:  (ni+1) × nj    × nk_win  faces, flat id = fi*nj*nk_win + fj*nk_win + fk
 *    ⊥Y:  ni    × (nj+1) × nk_win  faces, offset by nfX
 *    ⊥Z:  ni    × nj     × (nk_win+1) faces, offset by nfX+nfY
 *
 *  For a ⊥X face at window x-index fi (fi=0..ni), the face plane is at
 *    x = (d_noffset[0] + i0 + fi) - 0.5    (i.e. the -X face of node i0+fi)
 *  and covers nodes i0+fi-1 (its +X face) and i0+fi (its -X face).
 *
 *  flux stored in face_flux_d[face_id] = integrated E·n̂ over the face
 *  (positive outward for the node in the +direction; negate for the other node).
 *
 *  A second kernel (ewald_node_div_gpu_kernel) accumulates fluxes per node.
 *
 *****************************************************************************/

 /* Helper: evaluate E at point r[], same physics as before. */
__device__ static void eval_E_at_r(
    const double* __restrict__ Sk_sin,
    const double* __restrict__ Sk_cos,
    const double* __restrict__ kvec,
    const double* __restrict__ Gk,
    const int* __restrict__ kz_arr,
    const double* __restrict__ rho_data,
    const double* __restrict__ part_r,
    const double* __restrict__ part_q,
    int nktot, int nsites, int nhalo,
    int nk_species,
    const int* __restrict__ valency,
    int np,
    double irc, int irc_i,
    const double r[3],
    double E[3])
{
  E[0] = 0.0; E[1] = 0.0; E[2] = 0.0;

  /* Ludwig strides */
  int str_z = 1;
  int str_y = d_nlocal[2] + 2 * nhalo;
  int str_x = str_y * (d_nlocal[1] + 2 * nhalo);

  /* --- Fourier contribution --- */
  for (int kn = 0; kn < nktot; kn++) {
    double kx = kvec[3 * kn + 0];
    double ky = kvec[3 * kn + 1];
    double kz = kvec[3 * kn + 2];
    double kr = kx * r[0] + ky * r[1] + kz * r[2];
    double cos_kr = cos(kr);
    double sin_kr = sin(kr);
    double factor = (kz_arr[kn] > 0) ? 2.0 : 1.0;
    double impart = Sk_sin[kn] * cos_kr - Sk_cos[kn] * sin_kr;
    E[0] -= factor * Gk[kn] * kx * impart;
    E[1] -= factor * Gk[kn] * ky * impart;
    E[2] -= factor * Gk[kn] * kz * impart;
  }

  /* Dipole correction */
  E[0] += d_E_dipole[0];
  E[1] += d_E_dipole[1];
  E[2] += d_E_dipole[2];

  /* --- Real-space: lattice nodes --- */
  int rx0 = (int)floor(r[0]);
  int ry0 = (int)floor(r[1]);
  int rz0 = (int)floor(r[2]);

  for (int dx = -irc_i; dx <= irc_i; dx++) {
    for (int dy = -irc_i; dy <= irc_i; dy++) {
      for (int dz = -irc_i; dz <= irc_i; dz++) {
        int ix_g = rx0 + dx;
        int iy_g = ry0 + dy;
        int iz_g = rz0 + dz;

        double drx = r[0] - ix_g;
        double dry = r[1] - iy_g;
        double drz = r[2] - iz_g;

        if (drx > 0.5 * d_ltot[0]) drx -= d_ltot[0];
        if (drx < -0.5 * d_ltot[0]) drx += d_ltot[0];
        if (dry > 0.5 * d_ltot[1]) dry -= d_ltot[1];
        if (dry < -0.5 * d_ltot[1]) dry += d_ltot[1];
        if (drz > 0.5 * d_ltot[2]) drz -= d_ltot[2];
        if (drz < -0.5 * d_ltot[2]) drz += d_ltot[2];

        double r2 = drx * drx + dry * dry + drz * drz;
        if (r2 <= 0.0 || r2 > irc * irc) continue;

        double dist = sqrt(r2);

        int ix_l = ix_g - d_noffset[0];
        int iy_l = iy_g - d_noffset[1];
        int iz_l = iz_g - d_noffset[2];
        /* CHANGE INIT - EvalE_Halo - include halo nodes so border face points get full real-space sum */
        // Original: skipped nodes outside [0, nlocal), missing halo contributions at domain borders
        // if (ix_l < 0 || ix_l >= d_nlocal[0]) continue;
        // if (iy_l < 0 || iy_l >= d_nlocal[1]) continue;
        // if (iz_l < 0 || iz_l >= d_nlocal[2]) continue;
        if (ix_l < -nhalo || ix_l >= d_nlocal[0] + nhalo) continue;
        if (iy_l < -nhalo || iy_l >= d_nlocal[1] + nhalo) continue;
        if (iz_l < -nhalo || iz_l >= d_nlocal[2] + nhalo) continue;
        /* CHANGE END - EvalE_Halo */

        int n_idx = str_x * (nhalo + ix_l) + str_y * (nhalo + iy_l) + str_z * (nhalo + iz_l);

        double rho_node = 0.0;
        for (int n = 0; n < nk_species; n++) {
          int irho = nsites * n + n_idx;
          rho_node += d_eunit * valency[n] * rho_data[irho];
        }
        if (rho_node == 0.0) continue;

        double ar = d_alpha * dist;
        double erfcval = erfc(ar);
        double expval = exp(-ar * ar);
        double coeff = rho_node / (4.0 * M_PI * d_epsilon) *
          (erfcval / (r2 * dist) + 2.0 * d_alpha * d_rpi * expval / r2);
        E[0] += coeff * drx;
        E[1] += coeff * dry;
        E[2] += coeff * drz;
      }
    }
  }

  /* --- Real-space: particles --- */
  for (int p = 0; p < np; p++) {
    double prx = r[0] - part_r[3 * p + 0];
    double pry = r[1] - part_r[3 * p + 1];
    double prz = r[2] - part_r[3 * p + 2];

    if (prx > 0.5 * d_ltot[0]) prx -= d_ltot[0];
    if (prx < -0.5 * d_ltot[0]) prx += d_ltot[0];
    if (pry > 0.5 * d_ltot[1]) pry -= d_ltot[1];
    if (pry < -0.5 * d_ltot[1]) pry += d_ltot[1];
    if (prz > 0.5 * d_ltot[2]) prz -= d_ltot[2];
    if (prz < -0.5 * d_ltot[2]) prz += d_ltot[2];

    double r2 = prx * prx + pry * pry + prz * prz;
    if (r2 <= 0.0 || r2 > irc * irc) continue;

    double dist = sqrt(r2);
    double qp = part_q[p];
    double ar = d_alpha * dist;
    double coeff = qp / (4.0 * M_PI * d_epsilon) *
      (erfc(ar) / (r2 * dist) + 2.0 * d_alpha * d_rpi * exp(-ar * ar) / r2);
    E[0] += coeff * prx;
    E[1] += coeff * pry;
    E[2] += coeff * prz;
  }
}

/* -----------------------------------------------------------------------
 * ewald_face_E_gpu_kernel
 *
 * One thread per unique face.  Integrates E·n̂ over ns×ns quadrature points
 * and stores the result in face_flux_d[face_id].
 *
 * Face indexing (all 0-based window indices):
 *   ⊥X faces: fi ∈ [0, ni],  fj ∈ [0, nj-1],  fk ∈ [0, nk_win-1]
 *              id = fi*nj*nk_win + fj*nk_win + fk
 *   ⊥Y faces: fi ∈ [0, ni-1], fj ∈ [0, nj],   fk ∈ [0, nk_win-1]
 *              id = nfX + fi*(nj+1)*nk_win + fj*nk_win + fk
 *   ⊥Z faces: fi ∈ [0, ni-1], fj ∈ [0, nj-1], fk ∈ [0, nk_win]
 *              id = nfX + nfY + fi*nj*(nk_win+1) + fj*(nk_win+1) + fk
 *
 * A ⊥X face at fi has its plane at x = (noffset[X] + i0 + fi) - 0.5.
 * It is the +X face of node (i0+fi-1, j0+fj, k0+fk) [outward normal +X]
 * and the -X face of node (i0+fi,   j0+fj, k0+fk) [outward normal -X].
 * ewald_node_div_gpu_kernel applies the correct sign per node.
 * ----------------------------------------------------------------------- */
__global__ void ewald_face_E_gpu_kernel(
    double* __restrict__ face_flux_d,   /* [nfX + nfY + nfZ] integrated fluxes */
    const double* __restrict__ Sk_sin,
    const double* __restrict__ Sk_cos,
    const double* __restrict__ kvec,
    const double* __restrict__ Gk,
    const int* __restrict__ kz_arr,
    const double* __restrict__ rho_data,
    const double* __restrict__ part_r,
    const double* __restrict__ part_q,
    int nktot, int nsites, int nhalo,
    int nk_species,
    const int* __restrict__ valency,
    int np,
    int ns,
    int i0, int j0, int k0,
    int ni, int nj, int nk_win)
{
  /* Total unique faces */
  int nfX = (ni + 1) * nj * nk_win;
  int nfY = ni * (nj + 1) * nk_win;
  /* nfZ = ni * nj * (nk_win+1); */

  int face_id = blockIdx.x * blockDim.x + threadIdx.x;
  int nfaces_total = nfX + nfY + ni * nj * (nk_win + 1);
  if (face_id >= nfaces_total) return;

  double irc = d_ewald_rc;
  int    irc_i = (int)ceil(irc) + 1;
  double face_area_elt = 1.0 / (ns * ns);

  /* Decode face_id → axis, fi, fj, fk, face plane position, tangent axes */
  int ax;       /* normal axis: 0=X, 1=Y, 2=Z */
  int fi, fj, fk;

  if (face_id < nfX) {
    /* ⊥X face */
    ax = 0;
    int tmp = face_id;
    fi = tmp / (nj * nk_win);
    tmp -= fi * nj * nk_win;
    fj = tmp / nk_win;
    fk = tmp % nk_win;
  }
  else if (face_id < nfX + nfY) {
    /* ⊥Y face */
    ax = 1;
    int tmp = face_id - nfX;
    fi = tmp / ((nj + 1) * nk_win);
    tmp -= fi * (nj + 1) * nk_win;
    fj = tmp / nk_win;
    fk = tmp % nk_win;
  }
  else {
    /* ⊥Z face */
    ax = 2;
    int tmp = face_id - nfX - nfY;
    fi = tmp / (nj * (nk_win + 1));
    tmp -= fi * nj * (nk_win + 1);
    fj = tmp / (nk_win + 1);
    fk = tmp % (nk_win + 1);
  }

  /* Tangent axes */
  int t1 = (ax + 1) % 3;
  int t2 = (ax + 2) % 3;

  /* Face center in physical [X,Y,Z] coordinates.
   * The face is perpendicular to axis ax.  Along ax: plane at (noffset[ax] + i0/j0/k0 + f_ax) - 0.5
   * Along the two tangent axes: center of the corresponding node span.
   * fi,fj,fk always mean: fi→X extent, fj→Y extent, fk→Z extent regardless of ax. */
  double fc[3];
  fc[0] = (double)(d_noffset[0] + i0 + fi) - (ax == 0 ? 0.5 : 0.0);
  fc[1] = (double)(d_noffset[1] + j0 + fj) - (ax == 1 ? 0.5 : 0.0);
  fc[2] = (double)(d_noffset[2] + k0 + fk) - (ax == 2 ? 0.5 : 0.0);

  /* Integrate E·n̂ (with n̂ = +axis direction) over ns×ns quadrature points.
   * The sign relative to each adjacent node is applied in ewald_node_div_gpu_kernel. */
  double flux = 0.0;
  for (int s1 = 0; s1 < ns; s1++) {
    for (int s2 = 0; s2 < ns; s2++) {
      double r[3];
      r[ax] = fc[ax];
      r[t1] = fc[t1] - 0.5 + (s1 + 0.5) / ns;
      r[t2] = fc[t2] - 0.5 + (s2 + 0.5) / ns;

      double E[3];
      eval_E_at_r(Sk_sin, Sk_cos, kvec, Gk, kz_arr,
                  rho_data, part_r, part_q,
                  nktot, nsites, nhalo, nk_species, valency, np,
                  irc, irc_i, r, E);

      /* n̂ = +axis: store E[ax] (caller multiplies by ±1 per node) */
      flux += E[ax] * face_area_elt;
    }
  }

  face_flux_d[face_id] = flux;
}

/* -----------------------------------------------------------------------
 * ewald_node_div_gpu_kernel
 *
 * One thread per window node. Accumulates div(E) from the 6 adjacent face
 * fluxes stored by ewald_face_E_gpu_kernel.
 *
 * For node (wi, wj, wk) (0-based window indices):
 *   -X face: ⊥X face fi=wi,   fj=wj, fk=wk  →  flux × (-1)  [outward normal -X]
 *   +X face: ⊥X face fi=wi+1, fj=wj, fk=wk  →  flux × (+1)  [outward normal +X]
 *   -Y face: ⊥Y face fi=wi, fj=wj,   fk=wk  →  flux × (-1)
 *   +Y face: ⊥Y face fi=wi, fj=wj+1, fk=wk  →  flux × (+1)
 *   -Z face: ⊥Z face fi=wi, fj=wj, fk=wk    →  flux × (-1)
 *   +Z face: ⊥Z face fi=wi, fj=wj, fk=wk+1  →  flux × (+1)
 * ----------------------------------------------------------------------- */
__global__ void ewald_node_div_gpu_kernel(
    double* __restrict__ div_E_d,         /* [nwin] output */
    const double* __restrict__ face_flux_d,
    int ni, int nj, int nk_win)
{
  int win_idx = blockIdx.x * blockDim.x + threadIdx.x;
  int nwin = ni * nj * nk_win;
  if (win_idx >= nwin) return;

  int wi = win_idx / (nj * nk_win);
  int wj = (win_idx / nk_win) % nj;
  int wk = win_idx % nk_win;

  int nfX = (ni + 1) * nj * nk_win;
  int nfY = ni * (nj + 1) * nk_win;

  /* ⊥X face ids */
  int id_mX = (wi)*nj * nk_win + wj * nk_win + wk;
  int id_pX = (wi + 1) * nj * nk_win + wj * nk_win + wk;

  /* ⊥Y face ids */
  int id_mY = nfX + wi * (nj + 1) * nk_win + (wj)*nk_win + wk;
  int id_pY = nfX + wi * (nj + 1) * nk_win + (wj + 1) * nk_win + wk;

  /* ⊥Z face ids */
  int id_mZ = nfX + nfY + wi * nj * (nk_win + 1) + wj * (nk_win + 1) + (wk);
  int id_pZ = nfX + nfY + wi * nj * (nk_win + 1) + wj * (nk_win + 1) + (wk + 1);

  /* flux stored with n̂ = +axis; outward normals: -X→-1, +X→+1, etc. */
  div_E_d[win_idx] = -face_flux_d[id_mX]
    + face_flux_d[id_pX]
    - face_flux_d[id_mY]
    + face_flux_d[id_pY]
    - face_flux_d[id_mZ]
    + face_flux_d[id_pZ];
}

/* -----------------------------------------------------------------------
 * ewald_node_rho_eff_gpu_kernel
 *
 * One thread per window node. Computes rho_eff (ions + particle charge)
 * at the node center, same convention as ewald_write_poisson_verification:
 *   node center x0 = noffset[X] + (i0 + wi)   (0-based, no +1)
 * A particle at position r_p is assigned to this node if |x0 - r_p| ≤ 0.5.
 * ----------------------------------------------------------------------- */
__global__ void ewald_node_rho_eff_gpu_kernel(
    double* __restrict__ rho_eff_d,
    const double* __restrict__ rho_data,
    const double* __restrict__ part_r,
    const double* __restrict__ part_q,
    int nsites, int nhalo,
    int nk_species,
    const int* __restrict__ valency,
    int np,
    int i0, int j0, int k0,
    int ni, int nj, int nk_win)
{
  int win_idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (win_idx >= ni * nj * nk_win) return;

  int wi = win_idx / (nj * nk_win);
  int wj = (win_idx / nk_win) % nj;
  int wk = win_idx % nk_win;

  int i = i0 + wi;
  int j = j0 + wj;
  int k = k0 + wk;

  double cx = (double)(d_noffset[0] + i);
  double cy = (double)(d_noffset[1] + j);
  double cz = (double)(d_noffset[2] + k);

  int str_z = 1;
  int str_y = d_nlocal[2] + 2 * nhalo;
  int str_x = str_y * (d_nlocal[1] + 2 * nhalo);
  int node_idx = str_x * (nhalo + i) + str_y * (nhalo + j) + str_z * (nhalo + k);

  double rho_eff = 0.0;
  for (int n = 0; n < nk_species; n++) {
    int irho = nsites * n + node_idx;
    rho_eff += d_eunit * valency[n] * rho_data[irho];
  }
  for (int p = 0; p < np; p++) {
    double dx = cx - part_r[3 * p + 0];
    double dy = cy - part_r[3 * p + 1];
    double dz = cz - part_r[3 * p + 2];
    if (dx > 0.5 * d_ltot[0]) dx -= d_ltot[0];
    if (dx < -0.5 * d_ltot[0]) dx += d_ltot[0];
    if (dy > 0.5 * d_ltot[1]) dy -= d_ltot[1];
    if (dy < -0.5 * d_ltot[1]) dy += d_ltot[1];
    if (dz > 0.5 * d_ltot[2]) dz -= d_ltot[2];
    if (dz < -0.5 * d_ltot[2]) dz += d_ltot[2];
    if (fabs(dx) <= 0.5 && fabs(dy) <= 0.5 && fabs(dz) <= 0.5) {
      int share = ((fabs(dx) == 0.5) ? 2 : 1)
        * ((fabs(dy) == 0.5) ? 2 : 1)
        * ((fabs(dz) == 0.5) ? 2 : 1);
      rho_eff += part_q[p] / share;
    }
  }
  if (rho_eff_d != NULL) rho_eff_d[win_idx] = rho_eff;
}
/* CHANGE END - FaceDiv_GPU */

/*****************************************************************************
 *
 *  ewald_charge_sum_full_gpu_poisson_force
 *
 *  Same as ewald_charge_sum_full_gpu, plus:
 *    [N+1] Flag nodes near particle faces (diagnostic)
 *    [N+2] Sample E field on cube faces (CPU, after sinx_/cosx_ are ready)
 *    [N+3] Compute div(E), Poisson error, Maxwell stress force
 *    [N+4] Write results to files (same timestep as efield output)
 *
 *  The face sampling uses the precomputed Fourier structure factors sinx_/cosx_
 *  from ewald_charge_sum_sin_cos_terms(), called inside this function.
 *
 *****************************************************************************/

 /*CHANGE INIT - PoissonForceFunction - Poisson verification via GPU face-sampling (replaces CPU path) */

int ewald_charge_sum_full_gpu_poisson_force(ewald_charge_t* ewald, FILE* fp,
                                             ewald_self_table_t* self_table,
                                             int step) {

  /* This function assumes ewald_charge_sum_full_gpu() was already called this
   * timestep, so sinx_/cosx_[] are valid.  It does NOT repeat the Ewald sum.
   * Computes div(E) via GPU face-sampling (ewald_face_div_gpu_kernel), then
   * writes poisson_divergence and poisson_stencil output files.
   */

  if (ewald == NULL) return 0;

  PI_DOUBLE(pi);

  int nlocal[3], noffset[3];
  double ltot[3];
  cs_nlocal(ewald->cs, nlocal);
  cs_nlocal_offset(ewald->cs, noffset);
  cs_ltot(ewald->cs, ltot);

  int nhalo;
  cs_nhalo(ewald->cs, &nhalo);
  int nsites = ewald->psi->nsites;

  /* [1] Flag nodes where particle is close to a cube face (diagnostic only) */
  ewald_flag_near_particle_nodes(ewald, nlocal, noffset);

  // /* -----------------------------------------------------------------------
  //  * [2] Rebuild k-vector arrays (same as ewald_charge_sum_full_gpu)
  //  * --------------------------------------------------------------------- */
  // double fkx = 2.0 * pi / ltot[X];
  // double fky = 2.0 * pi / ltot[Y];
  // double fkz = 2.0 * pi / ltot[Z];
  // double r4alpha_sq = 1.0 / (4.0 * alpha_ * alpha_);
  // double b0 = 1.0 / (ltot[X] * ltot[Y] * ltot[Z] * epsilon_);

  // double* kvec_h   = (double*)malloc(3 * nktot_ * sizeof(double));
  // double* Gk_h     = (double*)malloc(nktot_ * sizeof(double));
  // int*    kz_arr_h = (int*)malloc(nktot_ * sizeof(int));

  // int kn = 0;
  // for (int kz = 0; kz <= nk_[Z]; kz++) {
  //   for (int ky = -nk_[Y]; ky <= nk_[Y]; ky++) {
  //     for (int kx = -nk_[X]; kx <= nk_[X]; kx++) {
  //       double k[3], ksq;
  //       k[X] = fkx * kx;
  //       k[Y] = fky * ky;
  //       k[Z] = fkz * kz;
  //       ksq = k[X]*k[X] + k[Y]*k[Y] + k[Z]*k[Z];
  //       if (ksq <= 0.0 || ksq > kmax_) continue;
  //       kvec_h[3*kn + X] = k[X];
  //       kvec_h[3*kn + Y] = k[Y];
  //       kvec_h[3*kn + Z] = k[Z];
  //       Gk_h[kn]     = b0 * exp(-r4alpha_sq * ksq) / ksq;
  //       kz_arr_h[kn] = kz;
  //       kn++;
  //     }
  //   }
  // }
  // int nk_actual = kn;
  // assert(nk_actual == nktot_);

  // double* kvec_d, *Gk_d;
  // int*    kz_arr_d;
  // double* Sk_sin_d, *Sk_cos_d;

  // cudaMalloc(&kvec_d,   3 * nk_actual * sizeof(double));
  // cudaMalloc(&Gk_d,     nk_actual * sizeof(double));
  // cudaMalloc(&kz_arr_d, nk_actual * sizeof(int));
  // cudaMalloc(&Sk_sin_d, nk_actual * sizeof(double));
  // cudaMalloc(&Sk_cos_d, nk_actual * sizeof(double));

  // cudaMemcpy(kvec_d,   kvec_h,   3 * nk_actual * sizeof(double), cudaMemcpyHostToDevice);
  // cudaMemcpy(Gk_d,     Gk_h,     nk_actual * sizeof(double),     cudaMemcpyHostToDevice);
  // cudaMemcpy(kz_arr_d, kz_arr_h, nk_actual * sizeof(int),        cudaMemcpyHostToDevice);
  // /* sinx_/cosx_ already populated from last ewald_charge_sum_full_gpu call */
  // cudaMemcpy(Sk_sin_d, sinx_, nk_actual * sizeof(double), cudaMemcpyHostToDevice);
  // cudaMemcpy(Sk_cos_d, cosx_, nk_actual * sizeof(double), cudaMemcpyHostToDevice);

  // free(kvec_h);
  // free(Gk_h);
  // free(kz_arr_h);

  // /* -----------------------------------------------------------------------
  //  * [3] Copy device constants
  //  * --------------------------------------------------------------------- */
  // cudaMemcpyToSymbol(d_alpha,    &alpha_,    sizeof(double));
  // cudaMemcpyToSymbol(d_epsilon,  &epsilon_,  sizeof(double));
  // cudaMemcpyToSymbol(d_eunit,    &eunit_,    sizeof(double));
  // cudaMemcpyToSymbol(d_rpi,      &rpi_,      sizeof(double));
  // cudaMemcpyToSymbol(d_ewald_rc, &ewald_rc_, sizeof(double));
  // cudaMemcpyToSymbol(d_nlocal,   nlocal,     3 * sizeof(int));
  // cudaMemcpyToSymbol(d_noffset,  noffset,    3 * sizeof(int));
  // cudaMemcpyToSymbol(d_ltot,     ltot,       3 * sizeof(double));

  // /* Dipole correction (copy from last Ewald sum — stored in d_E_dipole already,
  //  * but re-copy to be safe; use zero if dipole correction not needed here) */
  // double E_dipole_zero[3] = {0.0, 0.0, 0.0};
  // cudaMemcpyToSymbol(d_E_dipole, E_dipole_zero, 3 * sizeof(double));

  // /* -----------------------------------------------------------------------
  //  * [4] rho field on device
  //  * --------------------------------------------------------------------- */
  // field_memcpy(ewald->psi->rho, tdpMemcpyHostToDevice);

  // size_t data_offset = offsetof(field_t, data);
  // double* rho_data_d = NULL;
  // cudaMemcpy(&rho_data_d, (char*)(ewald->psi->rho->target) + data_offset,
  //            sizeof(double*), cudaMemcpyDeviceToHost);

  // int nk_species = 0;
  // psi_nk(ewald->psi, &nk_species);

  // /* Build valency array on host and copy to device */
  // int* valency_h = (int*)malloc(nk_species * sizeof(int));
  // for (int n = 0; n < nk_species; n++) {
  //   psi_valency(ewald->psi, n, &valency_h[n]);
  // }
  // int* valency_d;
  // cudaMalloc(&valency_d, nk_species * sizeof(int));
  // cudaMemcpy(valency_d, valency_h, nk_species * sizeof(int), cudaMemcpyHostToDevice);
  // free(valency_h);

  // /* -----------------------------------------------------------------------
  //  * [5] Particle data
  //  * --------------------------------------------------------------------- */
  // int nparticles = 0;
  // double* particle_r_h = NULL, *particle_q_h = NULL;
  // double* particle_r_d = NULL, *particle_q_d = NULL;

  // if (ewald->cinfo) {
  //   colloid_t* pc;
  //   colloids_info_local_head(ewald->cinfo, &pc);
  //   for (; pc; pc = pc->nextlocal) nparticles++;

  //   if (nparticles > 0) {
  //     particle_r_h = (double*)malloc(3 * nparticles * sizeof(double));
  //     particle_q_h = (double*)malloc(nparticles * sizeof(double));

  //     int p = 0;
  //     colloids_info_local_head(ewald->cinfo, &pc);
  //     for (; pc; pc = pc->nextlocal, p++) {
  //       particle_r_h[3*p + 0] = pc->s.r[X];
  //       particle_r_h[3*p + 1] = pc->s.r[Y];
  //       particle_r_h[3*p + 2] = pc->s.r[Z];
  //       particle_q_h[p] = pc->s.q0 - pc->s.q1;
  //     }

  //     cudaMalloc(&particle_r_d, 3 * nparticles * sizeof(double));
  //     cudaMalloc(&particle_q_d, nparticles * sizeof(double));
  //     cudaMemcpy(particle_r_d, particle_r_h, 3 * nparticles * sizeof(double), cudaMemcpyHostToDevice);
  //     cudaMemcpy(particle_q_d, particle_q_h, nparticles * sizeof(double), cudaMemcpyHostToDevice);
  //   }
  // }

  // /* -----------------------------------------------------------------------
  //  * [6] Compute window bounds and launch GPU kernel
  //  * --------------------------------------------------------------------- */
  // int i0, i1, j0, j1, k0, k1;
  // EWALD_WINDOW_BOUNDS(nlocal[X], i0, i1);
  // EWALD_WINDOW_BOUNDS(nlocal[Y], j0, j1);
  // EWALD_WINDOW_BOUNDS(nlocal[Z], k0, k1);

  // int ni     = i1 - i0 + 1;
  // int nj     = j1 - j0 + 1;
  // int nk_win = k1 - k0 + 1;
  // int nwin   = ni * nj * nk_win;

  // /* Unique face counts: each interior face computed once instead of twice */
  // int nfX = (ni+1) * nj * nk_win;
  // int nfY = ni * (nj+1) * nk_win;
  // int nfZ = ni * nj * (nk_win+1);
  // int nfaces = nfX + nfY + nfZ;

  // double* face_flux_d = NULL;
  // double* div_E_d     = NULL;
  // cudaMalloc(&face_flux_d, nfaces * sizeof(double));
  // cudaMalloc(&div_E_d,     nwin   * sizeof(double));

  // int threads = 128;

  // pe_info(ewald->pe, "  [Poisson] GPU face-sampling (NSAMPLE_FACE=%d, window %dx%dx%d, %d unique faces)...\n",
  //         EWALD_NSAMPLE_FACE, ni, nj, nk_win, nfaces);

  // /* [6a] Kernel 1: compute flux on each unique face (1 thread per face) */
  // int blocks_face = (nfaces + threads - 1) / threads;
  // ewald_face_E_gpu_kernel<<<blocks_face, threads>>>(
  //     face_flux_d,
  //     Sk_sin_d, Sk_cos_d,
  //     kvec_d, Gk_d, kz_arr_d,
  //     rho_data_d,
  //     particle_r_d, particle_q_d,
  //     nk_actual, nsites, nhalo,
  //     nk_species,
  //     valency_d,
  //     nparticles,
  //     EWALD_NSAMPLE_FACE,
  //     i0, j0, k0,
  //     ni, nj, nk_win);
  // cudaDeviceSynchronize();

  // /* [6b] Kernel 2: accumulate 6 face fluxes per node → div(E) */
  // int blocks_node = (nwin + threads - 1) / threads;
  // ewald_node_div_gpu_kernel<<<blocks_node, threads>>>(
  //     div_E_d,
  //     face_flux_d,
  //     ni, nj, nk_win);
  // cudaDeviceSynchronize();

  // /* -----------------------------------------------------------------------
  //  * [7] Copy div_E back to host and scatter into ewald->div_E[]
  //  * --------------------------------------------------------------------- */
  // double* div_E_h = (double*)malloc(nwin * sizeof(double));
  // cudaMemcpy(div_E_h, div_E_d, nwin * sizeof(double), cudaMemcpyDeviceToHost);

  // for (int wi = 0; wi < ni; wi++) {
  // for (int wj = 0; wj < nj; wj++) {
  // for (int wk = 0; wk < nk_win; wk++) {
  //   int win_idx = wi * nj * nk_win + wj * nk_win + wk;
  //   int i = i0 + wi;
  //   int j = j0 + wj;
  //   int k = k0 + wk;
  //   int idx = i * nlocal[Y] * nlocal[Z] + j * nlocal[Z] + k;
  //   ewald->div_E[idx] = div_E_h[win_idx];
  // }}}

  // free(div_E_h);

  /* -----------------------------------------------------------------------
   * [8] Write output files
   * --------------------------------------------------------------------- */
   // ewald_write_poisson_verification(ewald, nlocal, noffset, step);

  // ewald_write_poisson_stencil(ewald, nlocal, noffset, step);
  /* CHANGE INIT - ForceFluidOutput_PoissonForce - Write force_fluid to file at same frequency as efield */
  // ewald_write_force_fluid(ewald, nlocal, noffset, step);
  /* CHANGE END - ForceFluidOutput_PoissonForce */

  pe_info(ewald->pe, "  [Poisson] Done (step %d).\n", step);

  /* -----------------------------------------------------------------------
   * [9] Cleanup
   * --------------------------------------------------------------------- */
   // cudaFree(face_flux_d);
   // cudaFree(div_E_d);
   // cudaFree(kvec_d);
   // cudaFree(Gk_d);
   // cudaFree(kz_arr_d);
   // cudaFree(Sk_sin_d);
   // cudaFree(Sk_cos_d);
   // cudaFree(valency_d);
   // if (nparticles > 0) {
   //   free(particle_r_h);
   //   free(particle_q_h);
   //   cudaFree(particle_r_d);
   //   cudaFree(particle_q_d);
   // }

  return 0;
}

/*CHANGE END - PoissonForceFunction */

/*****************************************************************************
 *
 *  ewald_charge_sum_full_gpu_fourier
 *
 *  Same as ewald_charge_sum_full_gpu but computes only the Fourier-space
 *  (reciprocal-space) contributions.  Real-space and dipole corrections are
 *  omitted:
 *    - No real-space potential from particles or lattice nodes.
 *    - No real-space electric field from particles or lattice nodes.
 *    - No real-space force on fluid nodes from particles.
 *    - No real-space force on particles from lattice / other particles.
 *    - No dipole correction (E_dipole = 0).
 *
 *****************************************************************************/

int ewald_charge_sum_full_gpu_fourier(ewald_charge_t* ewald, FILE* fp,
                                      ewald_self_table_t* self_table) {

  int nlocal[3], noffset[3];
  double ltot[3];
  double fkx, fky, fkz;
  double r4alpha_sq, b0;
  PI_DOUBLE(pi);

  if (ewald == NULL) return 0;
  assert(fp);

  TIMER_start(TIMER_EWALD_TOTAL);

  cs_nlocal(ewald->cs, nlocal);
  cs_nlocal_offset(ewald->cs, noffset);
  cs_ltot(ewald->cs, ltot);

  int nhalo;
  cs_nhalo(ewald->cs, &nhalo);

  fkx = 2.0 * pi / ltot[X];
  fky = 2.0 * pi / ltot[Y];
  fkz = 2.0 * pi / ltot[Z];
  r4alpha_sq = 1.0 / (4.0 * alpha_ * alpha_);
  b0 = 1.0 / (ltot[X] * ltot[Y] * ltot[Z] * epsilon_);

  int ntotal_nodes = nlocal[X] * nlocal[Y] * nlocal[Z];
  int nsites = ewald->psi->nsites;

  pe_info(ewald->pe, "\n");
  pe_info(ewald->pe, "Ewald charge sum full Fourier-only (GPU): computing potential, forces, and fields...\n");
  pe_info(ewald->pe, "  Local nodes: %d x %d x %d = %d\n", nlocal[X], nlocal[Y], nlocal[Z], ntotal_nodes);
  pe_info(ewald->pe, "  Fourier terms: %d\n", nktot_);

  /* Copy constants to device */
  cudaMemcpyToSymbol(d_alpha, &alpha_, sizeof(double));
  cudaMemcpyToSymbol(d_eps_reg, &eps_reg_, sizeof(double));
  cudaMemcpyToSymbol(d_epsilon, &epsilon_, sizeof(double));
  cudaMemcpyToSymbol(d_beta, &beta_, sizeof(double));
  cudaMemcpyToSymbol(d_eunit, &eunit_, sizeof(double));
  cudaMemcpyToSymbol(d_rpi, &rpi_, sizeof(double));
  cudaMemcpyToSymbol(d_ewald_rc, &ewald_rc_, sizeof(double));
  cudaMemcpyToSymbol(d_fkx, &fkx, sizeof(double));
  cudaMemcpyToSymbol(d_fky, &fky, sizeof(double));
  cudaMemcpyToSymbol(d_fkz, &fkz, sizeof(double));
  cudaMemcpyToSymbol(d_r4alpha_sq, &r4alpha_sq, sizeof(double));
  cudaMemcpyToSymbol(d_b0, &b0, sizeof(double));
  cudaMemcpyToSymbol(d_kmax, &kmax_, sizeof(double));
  cudaMemcpyToSymbol(d_nk, nk_, 3 * sizeof(int));
  cudaMemcpyToSymbol(d_nlocal, nlocal, 3 * sizeof(int));
  cudaMemcpyToSymbol(d_noffset, noffset, 3 * sizeof(int));

  /* No dipole correction: set dipole field to zero on device */
  double E_dipole_zero[3] = { 0.0, 0.0, 0.0 };
  double M_zero[3] = { 0.0, 0.0, 0.0 };
  double dipole_pref_zero = 0.0;
  cudaMemcpyToSymbol(d_E_dipole, E_dipole_zero, 3 * sizeof(double));
  cudaMemcpyToSymbol(d_M_dipole, M_zero, 3 * sizeof(double));
  cudaMemcpyToSymbol(d_dipole_prefactor, &dipole_pref_zero, sizeof(double));

  /* ========================================================================
   * Precompute k-vectors and Green function on CPU, then copy to GPU
   * ======================================================================== */

  pe_info(ewald->pe, "  [1/5] Precomputing k-vectors and Green function...\n");

  double* kvec_h = (double*)malloc(3 * nktot_ * sizeof(double));
  double* Gk_h = (double*)malloc(nktot_ * sizeof(double));
  int* kz_arr_h = (int*)malloc(nktot_ * sizeof(int));

  int kn = 0;
  for (int kz = 0; kz <= nk_[Z]; kz++) {
    for (int ky = -nk_[Y]; ky <= nk_[Y]; ky++) {
      for (int kx = -nk_[X]; kx <= nk_[X]; kx++) {
        double k[3], ksq;
        k[X] = fkx * kx;
        k[Y] = fky * ky;
        k[Z] = fkz * kz;
        ksq = k[X] * k[X] + k[Y] * k[Y] + k[Z] * k[Z];
        if (ksq <= 0.0 || ksq > kmax_) continue;

        kvec_h[3 * kn + X] = k[X];
        kvec_h[3 * kn + Y] = k[Y];
        kvec_h[3 * kn + Z] = k[Z];
        Gk_h[kn] = b0 * exp(-r4alpha_sq * ksq) / ksq;
        kz_arr_h[kn] = kz;
        kn++;
      }
    }
  }
  int nk_actual = kn;

  /* Allocate device arrays */
  double* kvec_d, * Gk_d, * Sk_sin_d, * Sk_cos_d;
  int* kz_arr_d;

  cudaMalloc(&kvec_d, 3 * nk_actual * sizeof(double));
  cudaMalloc(&Gk_d, nk_actual * sizeof(double));
  cudaMalloc(&kz_arr_d, nk_actual * sizeof(int));
  cudaMalloc(&Sk_sin_d, nk_actual * sizeof(double));
  cudaMalloc(&Sk_cos_d, nk_actual * sizeof(double));

  cudaMemcpy(kvec_d, kvec_h, 3 * nk_actual * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(Gk_d, Gk_h, nk_actual * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(kz_arr_d, kz_arr_h, nk_actual * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemset(Sk_sin_d, 0, nk_actual * sizeof(double));
  cudaMemset(Sk_cos_d, 0, nk_actual * sizeof(double));

  /* Get device pointers for rho and psi fields */
  double* rho_data_d = NULL;
  double* psi_data_d = NULL;

  size_t data_offset = offsetof(field_t, data);
  cudaMemcpy(&rho_data_d, (char*)(ewald->psi->rho->target) + data_offset,
             sizeof(double*), cudaMemcpyDeviceToHost);
  cudaMemcpy(&psi_data_d, (char*)(ewald->psi->psi->target) + data_offset,
             sizeof(double*), cudaMemcpyDeviceToHost);

  /* Ensure rho is on device */
  field_memcpy(ewald->psi->rho, tdpMemcpyHostToDevice);

  /* ========================================================================
   * Collect particle data and copy to GPU
   * ======================================================================== */

  int nparticles = 0;
  double* particle_r_h = NULL, * particle_q_h = NULL;
  double* particle_r_d = NULL, * particle_q_d = NULL;
  double* Esub_d = NULL, * fex_d = NULL;

  if ((ewald->sources & EWALD_SOURCE_COLLOIDS) && ewald->cinfo) {
    colloid_t* pc;
    colloids_info_local_head(ewald->cinfo, &pc);
    for (; pc; pc = pc->nextlocal) nparticles++;

    if (nparticles > 0) {
      particle_r_h = (double*)malloc(3 * nparticles * sizeof(double));
      particle_q_h = (double*)malloc(nparticles * sizeof(double));

      int p = 0;
      colloids_info_local_head(ewald->cinfo, &pc);
      for (; pc; pc = pc->nextlocal) {
        particle_r_h[3 * p + 0] = pc->s.r[X];
        particle_r_h[3 * p + 1] = pc->s.r[Y];
        particle_r_h[3 * p + 2] = pc->s.r[Z];
        particle_q_h[p] = pc->s.q0 - pc->s.q1;
        p++;
      }

      cudaMalloc(&particle_r_d, 3 * nparticles * sizeof(double));
      cudaMalloc(&particle_q_d, nparticles * sizeof(double));
      cudaMalloc(&Esub_d, 3 * nparticles * sizeof(double));
      cudaMalloc(&fex_d, 3 * nparticles * sizeof(double));

      cudaMemcpy(particle_r_d, particle_r_h, 3 * nparticles * sizeof(double), cudaMemcpyHostToDevice);
      cudaMemcpy(particle_q_d, particle_q_h, nparticles * sizeof(double), cudaMemcpyHostToDevice);
    }
  }

  pe_info(ewald->pe, "  Particles: %d\n", nparticles);

  /* ========================================================================
   * PART 1: Compute structure factors on GPU (Fourier-space only)
   * ======================================================================== */

  pe_info(ewald->pe, "  [2/5] Computing structure factors S(k), C(k) on GPU...\n");

  int threads = 256;
  int blocks = (ntotal_nodes + threads - 1) / threads;

  if (ewald->sources & EWALD_SOURCE_LATTICE) {
    ewald_structure_factor_lattice_kernel << <blocks, threads >> > (
        rho_data_d, Sk_sin_d, Sk_cos_d, kvec_d, nk_actual, nsites, nhalo);
    cudaDeviceSynchronize();
  }

  if ((ewald->sources & EWALD_SOURCE_COLLOIDS) && nparticles > 0) {
    int p_blocks = (nparticles + threads - 1) / threads;
    ewald_structure_factor_particle_kernel << <p_blocks, threads >> > (
        particle_r_d, particle_q_d, Sk_sin_d, Sk_cos_d, kvec_d, nparticles, nk_actual);
    cudaDeviceSynchronize();
  }

  /* MPI reduce structure factors */
  pe_info(ewald->pe, "  [3/5] MPI reducing structure factors...\n");

  MPI_Comm comm;
  cs_cart_comm(ewald->cs, &comm);

  double* Sk_sin_h = (double*)malloc(nk_actual * sizeof(double));
  double* Sk_cos_h = (double*)malloc(nk_actual * sizeof(double));
  cudaMemcpy(Sk_sin_h, Sk_sin_d, nk_actual * sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(Sk_cos_h, Sk_cos_d, nk_actual * sizeof(double), cudaMemcpyDeviceToHost);
  MPI_Allreduce(MPI_IN_PLACE, Sk_sin_h, nk_actual, MPI_DOUBLE, MPI_SUM, comm);
  MPI_Allreduce(MPI_IN_PLACE, Sk_cos_h, nk_actual, MPI_DOUBLE, MPI_SUM, comm);
  cudaMemcpy(Sk_sin_d, Sk_sin_h, nk_actual * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(Sk_cos_d, Sk_cos_h, nk_actual * sizeof(double), cudaMemcpyHostToDevice);

  /* ========================================================================
   * PART 2: Compute potential on GPU (Fourier part only)
   * ======================================================================== */

  pe_info(ewald->pe, "  [4/5] Computing potential on lattice nodes (Fourier-only, GPU)...\n");

  ewald_potential_fourier_kernel << <blocks, threads >> > (
      psi_data_d, Sk_sin_d, Sk_cos_d, kvec_d, Gk_d, kz_arr_d,
      nk_actual, nsites, nhalo);
  cudaDeviceSynchronize();

  field_memcpy(ewald->psi->psi, tdpMemcpyDeviceToHost);

  if (ewald->psi && ewald->psi->psi_fourier) {
    double* psi_fourier_d = NULL;
    cudaMalloc(&psi_fourier_d, nsites * sizeof(double));
    cudaMemset(psi_fourier_d, 0, nsites * sizeof(double));

    ewald_potential_fourier_kernel << <blocks, threads >> > (
        psi_fourier_d, Sk_sin_d, Sk_cos_d, kvec_d, Gk_d, kz_arr_d,
        nk_actual, nsites, nhalo);
    cudaDeviceSynchronize();

    double* psi_fourier_h = (double*)malloc(nsites * sizeof(double));
    cudaMemcpy(psi_fourier_h, psi_fourier_d, nsites * sizeof(double), cudaMemcpyDeviceToHost);

    for (int ic = 1; ic <= nlocal[X]; ic++) {
      for (int jc = 1; jc <= nlocal[Y]; jc++) {
        for (int kc = 1; kc <= nlocal[Z]; kc++) {
          int index = cs_index(ewald->cs, ic, jc, kc);
          field_scalar_set(ewald->psi->psi_fourier, index, psi_fourier_h[index]);
        }
      }
    }
    free(psi_fourier_h);
    cudaFree(psi_fourier_d);
    field_memcpy(ewald->psi->psi_fourier, tdpMemcpyHostToDevice);
  }

  /* ========================================================================
   * PART 3: Compute force on fluid nodes (Fourier part only)
   * ======================================================================== */

  pe_info(ewald->pe, "  [5/5] Computing force on lattice nodes (Fourier-only, GPU)...\n");

  kahan_t F_fluid_total[3] = { kahan_zero(), kahan_zero(), kahan_zero() };

  if (ewald->hydro != NULL) {
    hydro_memcpy(ewald->hydro, tdpMemcpyDeviceToHost);

    double* efield_d = NULL;
    double* efield_fourier_d = NULL;
    double* force_d;
    cudaMalloc(&force_d, 3 * nsites * sizeof(double));
    cudaMemset(force_d, 0, 3 * nsites * sizeof(double));

    if (ewald->psi && ewald->psi->efield) {
      cudaMalloc(&efield_d, 3 * nsites * sizeof(double));
      cudaMemset(efield_d, 0, 3 * nsites * sizeof(double));
      cudaMalloc(&efield_fourier_d, 3 * nsites * sizeof(double));
      cudaMemset(efield_fourier_d, 0, 3 * nsites * sizeof(double));

      ewald_efield_fourier_kernel << <blocks, threads >> > (
          efield_d, Sk_sin_d, Sk_cos_d, kvec_d, Gk_d, kz_arr_d,
          nk_actual, nsites, nhalo);
      cudaDeviceSynchronize();
      cudaMemcpy(efield_fourier_d, efield_d, 3 * nsites * sizeof(double), cudaMemcpyDeviceToDevice);
    }

    ewald_force_fourier_kernel << <blocks, threads >> > (
        force_d, rho_data_d, Sk_sin_d, Sk_cos_d, kvec_d, Gk_d, kz_arr_d,
        nk_actual, nsites, nhalo);
    cudaDeviceSynchronize();

    double* force_h = (double*)malloc(3 * nsites * sizeof(double));
    cudaMemcpy(force_h, force_d, 3 * nsites * sizeof(double), cudaMemcpyDeviceToHost);

    double* efield_h = NULL;
    double* efield_fourier_h = NULL;
    if (efield_d != NULL) {
      efield_h = (double*)malloc(3 * nsites * sizeof(double));
      cudaMemcpy(efield_h, efield_d, 3 * nsites * sizeof(double), cudaMemcpyDeviceToHost);
      if (efield_fourier_d != NULL) {
        efield_fourier_h = (double*)malloc(3 * nsites * sizeof(double));
        cudaMemcpy(efield_fourier_h, efield_fourier_d, 3 * nsites * sizeof(double), cudaMemcpyDeviceToHost);
      }
    }

    double e_sum[3] = { 0.0, 0.0, 0.0 };
    int e_count = 0;

    for (int ic = 1; ic <= nlocal[X]; ic++) {
      for (int jc = 1; jc <= nlocal[Y]; jc++) {
        for (int kc = 1; kc <= nlocal[Z]; kc++) {
          int index = cs_index(ewald->cs, ic, jc, kc);
          double f[3] = { force_h[3 * index + 0], force_h[3 * index + 1], force_h[3 * index + 2] };
          hydro_f_local_add(ewald->hydro, index, f);

          if (efield_h != NULL) {
            double e_field[3] = { efield_h[3 * index + 0], efield_h[3 * index + 1], efield_h[3 * index + 2] };
            field_vector_set(ewald->psi->efield, index, e_field);
            e_sum[X] += fabs(e_field[X]);
            e_sum[Y] += fabs(e_field[Y]);
            e_sum[Z] += fabs(e_field[Z]);
            e_count++;

            if (efield_fourier_h != NULL && ewald->psi->efield_fourier) {
              double e_fourier[3] = { efield_fourier_h[3 * index + 0], efield_fourier_h[3 * index + 1], efield_fourier_h[3 * index + 2] };
              field_vector_set(ewald->psi->efield_fourier, index, e_fourier);
            }
          }

          kahan_add_double(&F_fluid_total[X], f[X]);
          kahan_add_double(&F_fluid_total[Y], f[Y]);
          kahan_add_double(&F_fluid_total[Z], f[Z]);
        }
      }
    }

    if (efield_h != NULL && e_count > 0) {
      pe_info(ewald->pe, "  Electric field sum |E|: (%14.7e, %14.7e, %14.7e) over %d sites\n",
              e_sum[X] / e_count, e_sum[Y] / e_count, e_sum[Z] / e_count, e_count);
    }

    if (efield_h != NULL) free(efield_h);
    if (efield_fourier_h != NULL) free(efield_fourier_h);
    if (efield_d != NULL) cudaFree(efield_d);
    if (efield_fourier_d != NULL) cudaFree(efield_fourier_d);

    if (ewald->psi && ewald->psi->efield) {
      field_memcpy(ewald->psi->efield, tdpMemcpyHostToDevice);
    }
    if (ewald->psi && ewald->psi->efield_fourier) {
      field_memcpy(ewald->psi->efield_fourier, tdpMemcpyHostToDevice);
    }

    free(force_h);
    cudaFree(force_d);
    hydro_memcpy(ewald->hydro, tdpMemcpyHostToDevice);
  }

  /* MPI reduce F_fluid_total */
  {
    MPI_Datatype kahan_dt;
    MPI_Op kahan_op;
    kahan_mpi_datatype(&kahan_dt);
    kahan_mpi_op_sum(&kahan_op);
    MPI_Allreduce(MPI_IN_PLACE, F_fluid_total, 3, kahan_dt, kahan_op, comm);
    MPI_Type_free(&kahan_dt);
    MPI_Op_free(&kahan_op);
  }

  /* ========================================================================
   * PART 4: Compute field and force on particles (Fourier part only)
   * ======================================================================== */

  kahan_t F_particle_total[3] = { kahan_zero(), kahan_zero(), kahan_zero() };

  if (nparticles > 0) {
    int p_blocks = (nparticles + threads - 1) / threads;

    ewald_particle_field_fourier_kernel << <p_blocks, threads >> > (
        Esub_d, fex_d, particle_r_d, particle_q_d,
        Sk_sin_d, Sk_cos_d, kvec_d, Gk_d, kz_arr_d,
        nparticles, nk_actual);
    cudaDeviceSynchronize();

    double* Esub_h = (double*)malloc(3 * nparticles * sizeof(double));
    double* fex_h = (double*)malloc(3 * nparticles * sizeof(double));
    cudaMemcpy(Esub_h, Esub_d, 3 * nparticles * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(fex_h, fex_d, 3 * nparticles * sizeof(double), cudaMemcpyDeviceToHost);

    int p = 0;
    double kt = 1.0 / beta_;
    colloid_t* pc;
    colloids_info_local_head(ewald->cinfo, &pc);
    for (; pc; pc = pc->nextlocal) {
      double Esub_x = Esub_h[3 * p + 0];
      double Esub_y = Esub_h[3 * p + 1];
      double Esub_z = Esub_h[3 * p + 2];

      if (self_table) {
        double xf = pc->s.r[X] - floor(pc->s.r[X]);
        double yf = pc->s.r[Y] - floor(pc->s.r[Y]);
        double zf = pc->s.r[Z] - floor(pc->s.r[Z]);
        double Eself[3];
        ewald_charge_self_field_interpolate(self_table, xf, yf, zf, Eself);
        double q_p = pc->s.q0 - pc->s.q1;
        Esub_x += q_p * Eself[0];
        Esub_y += q_p * Eself[1];
        Esub_z += q_p * Eself[2];

        double q_pc = pc->s.q0 - pc->s.q1;
        pc->Esub[X] = Esub_x;
        pc->Esub[Y] = Esub_y;
        pc->Esub[Z] = Esub_z;
        pc->fex[X] = q_pc * Esub_x * kt / eunit_;
        pc->fex[Y] = q_pc * Esub_y * kt / eunit_;
        pc->fex[Z] = q_pc * Esub_z * kt / eunit_;
      }
      else {
        pc->Esub[X] = Esub_x;
        pc->Esub[Y] = Esub_y;
        pc->Esub[Z] = Esub_z;
        pc->fex[X] = fex_h[3 * p + 0];
        pc->fex[Y] = fex_h[3 * p + 1];
        pc->fex[Z] = fex_h[3 * p + 2];
      }

      kahan_add_double(&F_particle_total[X], pc->fex[X]);
      kahan_add_double(&F_particle_total[Y], pc->fex[Y]);
      kahan_add_double(&F_particle_total[Z], pc->fex[Z]);
      p++;
    }

    free(Esub_h);
    free(fex_h);
  }

  /* MPI reduce particle forces */
  {
    MPI_Datatype kahan_dt;
    MPI_Op kahan_op;
    kahan_mpi_datatype(&kahan_dt);
    kahan_mpi_op_sum(&kahan_op);
    MPI_Allreduce(MPI_IN_PLACE, F_particle_total, 3, kahan_dt, kahan_op, comm);
    MPI_Type_free(&kahan_dt);
    MPI_Op_free(&kahan_op);
  }

  double F_fluid[3] = { kahan_sum(&F_fluid_total[X]),    kahan_sum(&F_fluid_total[Y]),    kahan_sum(&F_fluid_total[Z]) };
  double F_particle[3] = { kahan_sum(&F_particle_total[X]), kahan_sum(&F_particle_total[Y]), kahan_sum(&F_particle_total[Z]) };
  double F_diff[3] = { F_fluid[X] + F_particle[X], F_fluid[Y] + F_particle[Y], F_fluid[Z] + F_particle[Z] };

  fprintf(fp, "%1.15e,%1.15e, %1.15e, %1.15e,%1.15e, %1.15e, %1.15e,%1.15e, %1.15e, %1.15e,\n",
          sqrt(F_diff[X] * F_diff[X] + F_diff[Y] * F_diff[Y] + F_diff[Z] * F_diff[Z]),
          F_diff[X], F_diff[Y], F_diff[Z],
          F_fluid[X], F_fluid[Y], F_fluid[Z],
          F_particle[X], F_particle[Y], F_particle[Z]);

  pe_info(ewald->pe, "\n");
  pe_info(ewald->pe, "  Momentum conservation check:\n");
  pe_info(ewald->pe, "    F_fluid    = (%14.7e, %14.7e, %14.7e)\n", F_fluid[X], F_fluid[Y], F_fluid[Z]);
  pe_info(ewald->pe, "    F_particle = (%14.7e, %14.7e, %14.7e)\n", F_particle[X], F_particle[Y], F_particle[Z]);
  pe_info(ewald->pe, "    F_total    = (%14.7e, %14.7e, %14.7e)\n", F_diff[X], F_diff[Y], F_diff[Z]);
  pe_info(ewald->pe, "    |F_total|  = %14.7e\n", sqrt(F_diff[X] * F_diff[X] + F_diff[Y] * F_diff[Y] + F_diff[Z] * F_diff[Z]));

  ewald_charge_external_field(ewald);

  /* ========================================================================
   * Cleanup
   * ======================================================================== */

  free(kvec_h);
  free(Gk_h);
  free(kz_arr_h);
  free(Sk_sin_h);
  free(Sk_cos_h);

  cudaFree(kvec_d);
  cudaFree(Gk_d);
  cudaFree(kz_arr_d);
  cudaFree(Sk_sin_d);
  cudaFree(Sk_cos_d);

  if (nparticles > 0) {
    free(particle_r_h);
    free(particle_q_h);
    cudaFree(particle_r_d);
    cudaFree(particle_q_d);
    cudaFree(Esub_d);
    cudaFree(fex_d);
  }

  pe_info(ewald->pe, "Ewald charge sum full Fourier-only (GPU): complete.\n\n");

  TIMER_stop(TIMER_EWALD_TOTAL);

  return 0;
}

/*CHANGE INIT - 20260212 New ewald_charge_sum_FFT_full_gpu using cuFFT for Fourier-space*/

/*****************************************************************************
 *
 *  ewald_fft_efield_fourier_kernel
 *
 *  Compute Fourier part of electric field at each lattice node using the
 *  IFFT result of (-ik * psi_hat). This gives E = -grad(phi).
 *  Equivalent to ewald_efield_fourier_kernel but using IFFT output.
 *
 *  The IFFT of (+ik * psi_hat) = +grad(phi), so negate to get E = -grad(phi).
 *  Add dipole correction.
 *
 *****************************************************************************/

__global__ void ewald_fft_efield_fourier_kernel(
    double* __restrict__ efield_data,
    const double* __restrict__ Ex_real,
    const double* __restrict__ Ey_real,
    const double* __restrict__ Ez_real,
    int nx, int ny, int nz,
    int nhalo,
    double norm) {

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int n_total = nx * ny * nz;
  if (idx >= n_total) return;

  int iz = idx % nz;
  int iy = (idx / nz) % ny;
  int ix = idx / (nz * ny);

  int str_z = 1;
  int str_y = (nz + 2 * nhalo);
  int str_x = str_y * (ny + 2 * nhalo);
  int ic = ix + 1;
  int jc = iy + 1;
  int kc = iz + 1;
  int ludwig_idx = str_x * (nhalo + ic - 1) + str_y * (nhalo + jc - 1)
    + str_z * (nhalo + kc - 1);

  /* IFFT[+ik*psi_hat] = +grad(phi); negate to get E = -grad(phi) */
  double Ex = -(Ex_real[idx] * norm) + d_E_dipole[0];
  double Ey = -(Ey_real[idx] * norm) + d_E_dipole[1];
  double Ez = -(Ez_real[idx] * norm) + d_E_dipole[2];

  efield_data[3 * ludwig_idx + 0] = Ex;
  efield_data[3 * ludwig_idx + 1] = Ey;
  efield_data[3 * ludwig_idx + 2] = Ez;
}

/*****************************************************************************
 *
 *  ewald_charge_sum_FFT_full_gpu
 *
 *  Same as ewald_charge_sum_full_gpu but uses cuFFT for the Fourier-space
 *  lattice structure factor and potential/force computation.
 *
 *  Replaces:
 *    ewald_structure_factor_lattice_kernel  -> cuFFT forward
 *    ewald_potential_fourier_kernel         -> K2 (G*rho_hat) + IFFT
 *    ewald_efield_fourier_kernel            -> K3 (-ik*psi_hat) + IFFT
 *    ewald_force_fourier_kernel             -> K3 (+ik*psi_hat) + IFFT
 *    ewald_particle_field_fourier_kernel    -> K7 (direct sum on cuFFT array)
 *
 *  Everything else (real-space, dipole, MPI, particle handling) is identical
 *  to ewald_charge_sum_full_gpu.
 *
 *****************************************************************************/

int ewald_charge_sum_FFT_full_gpu(ewald_charge_t* ewald, FILE* fp) {

  int nlocal[3], noffset[3];
  double ltot[3];
  double fkx, fky, fkz;
  double r4alpha_sq, b0;
  int irc;
  PI_DOUBLE(pi);

  if (ewald == NULL) return 0;
  assert(fp);

  TIMER_start(TIMER_EWALD_TOTAL);

  cs_nlocal(ewald->cs, nlocal);
  cs_nlocal_offset(ewald->cs, noffset);
  cs_ltot(ewald->cs, ltot);

  int nhalo;
  cs_nhalo(ewald->cs, &nhalo);

  irc = (int)ceil(ewald_rc_);

  fkx = 2.0 * pi / ltot[X];
  fky = 2.0 * pi / ltot[Y];
  fkz = 2.0 * pi / ltot[Z];
  r4alpha_sq = 1.0 / (4.0 * alpha_ * alpha_);
  b0 = 1.0 / (ltot[X] * ltot[Y] * ltot[Z] * epsilon_);

  int nx = nlocal[X];
  int ny = nlocal[Y];
  int nz = nlocal[Z];
  int n_real = nx * ny * nz;
  int nz_complex = nz / 2 + 1;
  int n_complex = nx * ny * nz_complex;
  int nsites = ewald->psi->nsites;
  double norm = 1.0 / (double)n_real;

  pe_info(ewald->pe, "\n");
  pe_info(ewald->pe, "Ewald charge sum FFT (GPU): computing potential, forces, and fields...\n");
  pe_info(ewald->pe, "  Local nodes: %d x %d x %d = %d\n", nx, ny, nz, n_real);
  pe_info(ewald->pe, "  Real-space cutoff: %.2f (irc=%d)\n", ewald_rc_, irc);

  /* Copy constants to device */
  cudaMemcpyToSymbol(d_alpha, &alpha_, sizeof(double));
  cudaMemcpyToSymbol(d_eps_reg, &eps_reg_, sizeof(double));
  cudaMemcpyToSymbol(d_epsilon, &epsilon_, sizeof(double));
  cudaMemcpyToSymbol(d_beta, &beta_, sizeof(double));
  cudaMemcpyToSymbol(d_eunit, &eunit_, sizeof(double));
  cudaMemcpyToSymbol(d_rpi, &rpi_, sizeof(double));
  cudaMemcpyToSymbol(d_ewald_rc, &ewald_rc_, sizeof(double));
  cudaMemcpyToSymbol(d_fkx, &fkx, sizeof(double));
  cudaMemcpyToSymbol(d_fky, &fky, sizeof(double));
  cudaMemcpyToSymbol(d_fkz, &fkz, sizeof(double));
  cudaMemcpyToSymbol(d_r4alpha_sq, &r4alpha_sq, sizeof(double));
  cudaMemcpyToSymbol(d_b0, &b0, sizeof(double));
  cudaMemcpyToSymbol(d_kmax, &kmax_, sizeof(double));
  cudaMemcpyToSymbol(d_nk, nk_, 3 * sizeof(int));
  cudaMemcpyToSymbol(d_nlocal, nlocal, 3 * sizeof(int));
  cudaMemcpyToSymbol(d_noffset, noffset, 3 * sizeof(int));

  /* ========================================================================
   * [1] Precompute G_ewald(k) in cuFFT layout on CPU, copy to GPU.
   *     G(k) = (1/k^2) * exp(-k^2/(4*alpha^2)) / epsilon = b0/k^2 * exp(-r4alpha_sq*k^2)
   *     G(k=0) = 0.  Same formula as Gk_h in ewald_charge_sum_full_gpu.
   * ======================================================================== */

  pe_info(ewald->pe, "  [1/6] Precomputing Ewald Green function G(k)...\n");

  double* G_ewald_h = (double*)malloc(n_complex * sizeof(double));

  for (int ix = 0; ix < nx; ix++) {
    int kx_idx = (ix <= nx / 2) ? ix : ix - nx;
    double kx = fkx * kx_idx;
    for (int iy = 0; iy < ny; iy++) {
      int ky_idx = (iy <= ny / 2) ? iy : iy - ny;
      double ky = fky * ky_idx;
      for (int iz = 0; iz < nz_complex; iz++) {
        double kz = fkz * iz;
        int idx = ix * ny * nz_complex + iy * nz_complex + iz;

        if (ix == 0 && iy == 0 && iz == 0) {
          G_ewald_h[idx] = 0.0;
          continue;
        }

        double ksq = kx * kx + ky * ky + kz * kz;
        G_ewald_h[idx] = b0 * exp(-r4alpha_sq * ksq) / ksq;
      }
    }
  }

  double* G_ewald_d = NULL;
  cudaMalloc(&G_ewald_d, n_complex * sizeof(double));
  cudaMemcpy(G_ewald_d, G_ewald_h, n_complex * sizeof(double), cudaMemcpyHostToDevice);
  free(G_ewald_h);

  /* ========================================================================
   * [2] Get device pointers for rho and psi fields, ensure rho is on device
   * ======================================================================== */

  double* rho_data_d = NULL;
  double* psi_data_d = NULL;

  size_t data_offset = offsetof(field_t, data);
  cudaMemcpy(&rho_data_d, (char*)(ewald->psi->rho->target) + data_offset,
             sizeof(double*), cudaMemcpyDeviceToHost);
  cudaMemcpy(&psi_data_d, (char*)(ewald->psi->psi->target) + data_offset,
             sizeof(double*), cudaMemcpyDeviceToHost);

  field_memcpy(ewald->psi->rho, tdpMemcpyHostToDevice);

  /* ========================================================================
   * [3] Collect particle data and copy to GPU (same as ewald_charge_sum_full_gpu)
   * ======================================================================== */

  int nparticles = 0;
  double* particle_r_h = NULL, * particle_q_h = NULL;
  double* particle_r_d = NULL, * particle_q_d = NULL;
  double* Esub_d = NULL, * fex_d = NULL;

  if ((ewald->sources & EWALD_SOURCE_COLLOIDS) && ewald->cinfo) {
    colloid_t* pc;
    colloids_info_local_head(ewald->cinfo, &pc);
    for (; pc; pc = pc->nextlocal) nparticles++;

    if (nparticles > 0) {
      particle_r_h = (double*)malloc(3 * nparticles * sizeof(double));
      particle_q_h = (double*)malloc(nparticles * sizeof(double));

      int p = 0;
      colloids_info_local_head(ewald->cinfo, &pc);
      for (; pc; pc = pc->nextlocal) {
        particle_r_h[3 * p + 0] = pc->s.r[X];
        particle_r_h[3 * p + 1] = pc->s.r[Y];
        particle_r_h[3 * p + 2] = pc->s.r[Z];
        particle_q_h[p] = pc->s.q0 - pc->s.q1;
        p++;
      }

      cudaMalloc(&particle_r_d, 3 * nparticles * sizeof(double));
      cudaMalloc(&particle_q_d, nparticles * sizeof(double));
      cudaMalloc(&Esub_d, 3 * nparticles * sizeof(double));
      cudaMalloc(&fex_d, 3 * nparticles * sizeof(double));

      cudaMemcpy(particle_r_d, particle_r_h, 3 * nparticles * sizeof(double), cudaMemcpyHostToDevice);
      cudaMemcpy(particle_q_d, particle_q_h, nparticles * sizeof(double), cudaMemcpyHostToDevice);
    }
  }

  pe_info(ewald->pe, "  Particles: %d\n", nparticles);

  /* MPI comm */
  MPI_Comm comm;
  cs_cart_comm(ewald->cs, &comm);

  /* ========================================================================
   * [4] Compute structure factor S(k) using FFT
   *
   *  In ewald_charge_sum_full_gpu, S(k) = Sk_cos[kn] + i*Sk_sin[kn]
   *  where Sk_cos = Re[S] = sum q*cos(k.r), Sk_sin = Im[S] = sum q*sin(k.r).
   *
   *  cuFFT R2C convention: rho_hat[k] = sum_r rho(r)*exp(-ik.r)
   *    rho_hat.x = sum q*cos(k.r)  = Sk_cos
   *    rho_hat.y = -sum q*sin(k.r) = -Sk_sin
   *
   *  So: Sk_cos = rho_hat.x, Sk_sin = -rho_hat.y
   *
   *  Step 1: K1 copies rho_lattice to contiguous FFT array
   *  Step 2: FFT forward (D2Z)
   *  Step 3: MPI reduce (sum over domains)
   *  Step 4: K6 adds particle structure factors
   *
   *  Result: rho_complex_d[k].x = Sk_cos_total, rho_complex_d[k].y = -Sk_sin_total
   * ======================================================================== */

  pe_info(ewald->pe, "  [2/6] Computing structure factors S(k) via FFT...\n");

  int threads = 256;
  int blocks_r = (n_real + threads - 1) / threads;
  int blocks_c = (n_complex + threads - 1) / threads;

  double* rho_real_d = NULL;
  cufftDoubleComplex* rho_complex_d = NULL;

  cudaMalloc(&rho_real_d, n_real * sizeof(double));
  cudaMalloc(&rho_complex_d, n_complex * sizeof(cufftDoubleComplex));
  cudaMemset(rho_complex_d, 0, n_complex * sizeof(cufftDoubleComplex));

  cufftHandle plan_fwd, plan_bwd;
  cufftPlan3d(&plan_fwd, nx, ny, nz, CUFFT_D2Z);
  cufftPlan3d(&plan_bwd, nx, ny, nz, CUFFT_Z2D);

  /* K1: copy rho_lattice to FFT array */
  if (ewald->sources & EWALD_SOURCE_LATTICE) {
    ewald_copy_rho_to_fft_kernel << <blocks_r, threads >> > (
        rho_real_d, rho_data_d, nsites, nhalo, nx, ny, nz);
    cudaDeviceSynchronize();
  }

  /* FFT forward: rho_hat_lattice = FFT(rho_lattice) */
  cufftExecD2Z(plan_fwd, (cufftDoubleReal*)rho_real_d,
               (cufftDoubleComplex*)rho_complex_d);
  cudaDeviceSynchronize();

  /* MPI reduce rho_hat (sum over MPI domains) */
  {
    int npe;
    MPI_Comm_size(comm, &npe);
    if (npe > 1) {
      pe_info(ewald->pe, "  [3/6] MPI reducing structure factors...\n");
      double* complex_h = (double*)malloc(2 * n_complex * sizeof(double));
      cudaMemcpy(complex_h, rho_complex_d,
                 n_complex * sizeof(cufftDoubleComplex), cudaMemcpyDeviceToHost);
      MPI_Allreduce(MPI_IN_PLACE, complex_h, 2 * n_complex, MPI_DOUBLE, MPI_SUM, comm);
      cudaMemcpy(rho_complex_d, complex_h,
                 n_complex * sizeof(cufftDoubleComplex), cudaMemcpyHostToDevice);
      free(complex_h);
    }
  }

  /* K6: add particle structure factors.
   * cuFFT convention: exp(-ik.r) = cos(kr) - i*sin(kr)
   * So rho_hat[k].x += q*cos(kr), rho_hat[k].y += -q*sin(kr)
   * This matches: Sk_cos = rho_hat.x, Sk_sin = -rho_hat.y */
  if ((ewald->sources & EWALD_SOURCE_COLLOIDS) && nparticles > 0) {
    /* Need all particles across MPI for structure factor */
    int nparticles_local = nparticles;
    int nparticles_all = nparticles;
    double* all_particle_r_d = particle_r_d;
    double* all_particle_q_d = particle_q_d;
    double* all_particle_r_h = particle_r_h;
    double* all_particle_q_h = particle_q_h;

    int npe;
    MPI_Comm_size(comm, &npe);
    if (npe > 1) {
      int* counts = (int*)malloc(npe * sizeof(int));
      MPI_Allgather(&nparticles_local, 1, MPI_INT, counts, 1, MPI_INT, comm);
      nparticles_all = 0;
      for (int i = 0; i < npe; i++) nparticles_all += counts[i];

      if (nparticles_all > 0) {
        all_particle_r_h = (double*)malloc(3 * nparticles_all * sizeof(double));
        all_particle_q_h = (double*)malloc(nparticles_all * sizeof(double));

        int* displs_r = (int*)malloc(npe * sizeof(int));
        int* counts_r = (int*)malloc(npe * sizeof(int));
        displs_r[0] = 0; counts_r[0] = counts[0] * 3;
        for (int i = 1; i < npe; i++) {
          counts_r[i] = counts[i] * 3;
          displs_r[i] = displs_r[i - 1] + counts_r[i - 1];
        }
        MPI_Allgatherv(particle_r_h, nparticles_local * 3, MPI_DOUBLE,
                       all_particle_r_h, counts_r, displs_r, MPI_DOUBLE, comm);

        int* displs_q = (int*)malloc(npe * sizeof(int));
        displs_q[0] = 0;
        for (int i = 1; i < npe; i++) displs_q[i] = displs_q[i - 1] + counts[i - 1];
        MPI_Allgatherv(particle_q_h, nparticles_local, MPI_DOUBLE,
                       all_particle_q_h, counts, displs_q, MPI_DOUBLE, comm);

        cudaMalloc(&all_particle_r_d, 3 * nparticles_all * sizeof(double));
        cudaMalloc(&all_particle_q_d, nparticles_all * sizeof(double));
        cudaMemcpy(all_particle_r_d, all_particle_r_h, 3 * nparticles_all * sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy(all_particle_q_d, all_particle_q_h, nparticles_all * sizeof(double), cudaMemcpyHostToDevice);

        free(displs_r); free(counts_r); free(displs_q);
      }
      free(counts);
    }

    if (nparticles_all > 0) {
      int p_blocks = (nparticles_all + threads - 1) / threads;
      ewald_particle_sfactor_fft_kernel << <p_blocks, threads >> > (
          rho_complex_d, all_particle_r_d, all_particle_q_d, nparticles_all,
          fkx, fky, fkz, nx, ny, nz_complex);
      cudaDeviceSynchronize();
    }

    if (npe > 1 && all_particle_r_d != particle_r_d) {
      if (all_particle_r_h != particle_r_h) free(all_particle_r_h);
      if (all_particle_q_h != particle_q_h) free(all_particle_q_h);
      cudaFree(all_particle_r_d);
      cudaFree(all_particle_q_d);
    }
  }

  /* rho_complex_d now contains rho_hat_total:
   *   rho_hat.x = Sk_cos_total = sum_all q*cos(k.r)
   *   rho_hat.y = -Sk_sin_total = -sum_all q*sin(k.r)
   */

   /* ========================================================================
    * [5] Compute dipole moment (same as ewald_charge_sum_full_gpu)
    * ======================================================================== */

  kahan_t M_dipole[3] = { kahan_zero(), kahan_zero(), kahan_zero() };
  kahan_t Q_total_k = kahan_zero();

  if (nparticles > 0) {
    for (int p = 0; p < nparticles; p++) {
      double q = particle_q_h[p];
      kahan_add_double(&M_dipole[X], q * particle_r_h[3 * p + 0]);
      kahan_add_double(&M_dipole[Y], q * particle_r_h[3 * p + 1]);
      kahan_add_double(&M_dipole[Z], q * particle_r_h[3 * p + 2]);
      kahan_add_double(&Q_total_k, q);
    }
  }

  if ((ewald->sources & EWALD_SOURCE_LATTICE) && ewald->psi) {
    field_memcpy(ewald->psi->rho, tdpMemcpyDeviceToHost);

    /* Accumulate each species separately to avoid catastrophic cancellation
     * when rho[0] and rho[1] are nearly equal (q_node = v0*rho0 + v1*rho1 ~ 0) */
    int nk_ewald;
    psi_nk(ewald->psi, &nk_ewald);
    kahan_t* Q_species_k = (kahan_t*)calloc(nk_ewald, sizeof(kahan_t));
    kahan_t* M_species_k[3];
    for (int d = 0; d < 3; d++) {
      M_species_k[d] = (kahan_t*)calloc(nk_ewald, sizeof(kahan_t));
    }
    for (int s = 0; s < nk_ewald; s++) Q_species_k[s] = kahan_zero();
    for (int d = 0; d < 3; d++)
      for (int s = 0; s < nk_ewald; s++) M_species_k[d][s] = kahan_zero();

    for (int ic = 1; ic <= nlocal[X]; ic++) {
      for (int jc = 1; jc <= nlocal[Y]; jc++) {
        for (int kc = 1; kc <= nlocal[Z]; kc++) {
          int index = cs_index(ewald->cs, ic, jc, kc);

          double rx = (double)(noffset[X] + ic);
          double ry = (double)(noffset[Y] + jc);
          double rz = (double)(noffset[Z] + kc);

          for (int s = 0; s < nk_ewald; s++) {
            int val;
            double rho_s;
            psi_valency(ewald->psi, s, &val);
            psi_rho(ewald->psi, index, s, &rho_s);
            double q_s = val * rho_s;
            kahan_add_double(&Q_species_k[s], q_s);
            kahan_add_double(&M_species_k[X][s], q_s * rx);
            kahan_add_double(&M_species_k[Y][s], q_s * ry);
            kahan_add_double(&M_species_k[Z][s], q_s * rz);
          }
        }
      }
    }

    /* Combine species after full summation */
    for (int s = 0; s < nk_ewald; s++) {
      kahan_add_double(&Q_total_k, kahan_sum(&Q_species_k[s]));
      kahan_add_double(&M_dipole[X], kahan_sum(&M_species_k[X][s]));
      kahan_add_double(&M_dipole[Y], kahan_sum(&M_species_k[Y][s]));
      kahan_add_double(&M_dipole[Z], kahan_sum(&M_species_k[Z][s]));
    }

    for (int d = 0; d < 3; d++) free(M_species_k[d]);
    free(Q_species_k);
  }

  {
    kahan_t M_reduce[4] = { M_dipole[X], M_dipole[Y], M_dipole[Z], Q_total_k };
    MPI_Datatype kahan_dt;
    MPI_Op kahan_op;
    kahan_mpi_datatype(&kahan_dt);
    kahan_mpi_op_sum(&kahan_op);
    MPI_Allreduce(MPI_IN_PLACE, M_reduce, 4, kahan_dt, kahan_op, comm);
    MPI_Type_free(&kahan_dt);
    MPI_Op_free(&kahan_op);
    M_dipole[X] = M_reduce[0];
    M_dipole[Y] = M_reduce[1];
    M_dipole[Z] = M_reduce[2];
    Q_total_k = M_reduce[3];
  }

  double M[3] = { kahan_sum(&M_dipole[X]), kahan_sum(&M_dipole[Y]), kahan_sum(&M_dipole[Z]) };
  double Q_total = kahan_sum(&Q_total_k);
  double V = ltot[X] * ltot[Y] * ltot[Z];

  double dipole_prefactor = 4.0 * pi / ((1.0 + 2.0 * ewald->epsilon_prime) * V);
  double E_dipole[3] = { -dipole_prefactor * M[X],
                        -dipole_prefactor * M[Y],
                        -dipole_prefactor * M[Z] };

  pe_info(ewald->pe, "  Dipole moment M = (%14.7e, %14.7e, %14.7e)\n", M[X], M[Y], M[Z]);
  pe_info(ewald->pe, "  Total charge Q = %14.7e\n", Q_total);
  pe_info(ewald->pe, "  Dipole epsilon' = %14.7e\n", ewald->epsilon_prime);
  pe_info(ewald->pe, "  Dipole field E_dip = (%14.7e, %14.7e, %14.7e)\n",
          E_dipole[X], E_dipole[Y], E_dipole[Z]);

  cudaMemcpyToSymbol(d_E_dipole, E_dipole, 3 * sizeof(double));
  cudaMemcpyToSymbol(d_M_dipole, M, 3 * sizeof(double));
  double dipole_pref_dev = dipole_prefactor;
  cudaMemcpyToSymbol(d_dipole_prefactor, &dipole_pref_dev, sizeof(double));

  /* ========================================================================
   * [6a] Allocate shared FFT work arrays
   * ======================================================================== */

  cufftDoubleComplex* psi_complex_d = NULL;
  cufftDoubleComplex* efield_complex_d = NULL;
  double* psi_real_d = NULL;
  double* Ex_real_d = NULL;
  double* Ey_real_d = NULL;
  double* Ez_real_d = NULL;

  cudaMalloc(&psi_complex_d, n_complex * sizeof(cufftDoubleComplex));
  cudaMalloc(&efield_complex_d, n_complex * sizeof(cufftDoubleComplex));
  cudaMalloc(&psi_real_d, n_real * sizeof(double));
  cudaMalloc(&Ex_real_d, n_real * sizeof(double));
  cudaMalloc(&Ey_real_d, n_real * sizeof(double));
  cudaMalloc(&Ez_real_d, n_real * sizeof(double));

  /* K2: psi_hat = G * rho_hat_total */
  ewald_green_multiply_potential_kernel << <blocks_c, threads >> > (
      psi_complex_d, rho_complex_d, G_ewald_d, n_complex);
  cudaDeviceSynchronize();

  /* ========================================================================
   * [6b] Potential on lattice nodes
   *
   *  phi(r) = IFFT[psi_hat] / N + dipole correction
   *  Equivalent to ewald_potential_fourier_kernel which computes:
   *    phi(r) = 2 * sum_{kz>0} Gk*(Sk_cos*cos + Sk_sin*sin)
   *           + sum_{kz=0}    Gk*(Sk_cos*cos + Sk_sin*sin)
   *  The IFFT gives exactly this sum (with factor 1/N for normalization).
   * ======================================================================== */

  pe_info(ewald->pe, "  [4/6] Computing potential on lattice nodes...\n");

  /* IFFT: psi(r) = IFFT[psi_hat] (unnormalized) */
  cufftExecZ2D(plan_bwd, psi_complex_d, psi_real_d);
  cudaDeviceSynchronize();

  /* K4: copy to Ludwig layout, normalize by 1/N, add dipole */
  ewald_copy_psi_from_fft_kernel << <blocks_r, threads >> > (
      psi_data_d, psi_real_d, nx, ny, nz, nhalo, norm);
  cudaDeviceSynchronize();

  /* Real-space from particles (unchanged from direct version) */
  if (nparticles > 0) {
    int blocks_nodes = (n_real + threads - 1) / threads;
    ewald_potential_real_particle_kernel << <blocks_nodes, threads >> > (
        psi_data_d, particle_r_d, particle_q_d, nparticles, nsites, nhalo, irc);
    cudaDeviceSynchronize();
  }

  field_memcpy(ewald->psi->psi, tdpMemcpyDeviceToHost);

  /* ========================================================================
   * [6c] Force and electric field on lattice nodes
   *
   *  grad(phi)_alpha(r) = IFFT[+i*k_alpha * psi_hat] / N
   *  F_alpha(r) = q(r) * (grad(phi)_alpha + E_dipole_alpha)   [= K5]
   *  E_alpha(r) = -grad(phi)_alpha + E_dipole_alpha            [= K8]
   *
   *  This is equivalent to ewald_force_fourier_kernel which uses:
   *    F_alpha = q * 2*sum Gk * k_alpha * im_part  (im_part = Sk_sin*cos - Sk_cos*sin)
   *  and ewald_efield_fourier_kernel which uses:
   *    E_alpha = -2*sum Gk * k_alpha * im_part + E_dipole
   * ======================================================================== */

  pe_info(ewald->pe, "  [5/6] Computing force on lattice nodes...\n");

  kahan_t F_fluid_total[3] = { kahan_zero(), kahan_zero(), kahan_zero() };

  if (ewald->hydro != NULL) {
    hydro_memcpy(ewald->hydro, tdpMemcpyDeviceToHost);

    double* efield_d = NULL;
    double* efield_fourier_d = NULL;
    double* efield_real_d = NULL;
    double* force_d = NULL;

    cudaMalloc(&force_d, 3 * nsites * sizeof(double));
    cudaMemset(force_d, 0, 3 * nsites * sizeof(double));

    if (ewald->psi && ewald->psi->efield) {
      cudaMalloc(&efield_d, 3 * nsites * sizeof(double));
      cudaMemset(efield_d, 0, 3 * nsites * sizeof(double));
      cudaMalloc(&efield_fourier_d, 3 * nsites * sizeof(double));
      cudaMemset(efield_fourier_d, 0, 3 * nsites * sizeof(double));
      cudaMalloc(&efield_real_d, 3 * nsites * sizeof(double));
      cudaMemset(efield_real_d, 0, 3 * nsites * sizeof(double));
    }

    /* Compute grad(phi)_x, _y, _z via FFT and store results */

    /* grad(phi)_x = IFFT[+ik_x * psi_hat] / N */
    ewald_efield_multiply_kernel << <blocks_c, threads >> > (
        efield_complex_d, psi_complex_d, 0, fkx, fky, fkz, nx, ny, nz_complex);
    cudaDeviceSynchronize();
    cufftExecZ2D(plan_bwd, efield_complex_d, Ex_real_d);
    cudaDeviceSynchronize();

    /* grad(phi)_y */
    ewald_efield_multiply_kernel << <blocks_c, threads >> > (
        efield_complex_d, psi_complex_d, 1, fkx, fky, fkz, nx, ny, nz_complex);
    cudaDeviceSynchronize();
    cufftExecZ2D(plan_bwd, efield_complex_d, Ey_real_d);
    cudaDeviceSynchronize();

    /* grad(phi)_z */
    ewald_efield_multiply_kernel << <blocks_c, threads >> > (
        efield_complex_d, psi_complex_d, 2, fkx, fky, fkz, nx, ny, nz_complex);
    cudaDeviceSynchronize();
    cufftExecZ2D(plan_bwd, efield_complex_d, Ez_real_d);
    cudaDeviceSynchronize();

    /* K5: F = q * (grad(phi) + E_dipole) — adds Fourier force to force_d */
    ewald_copy_force_from_fft_kernel << <blocks_r, threads >> > (
        force_d, Ex_real_d, Ey_real_d, Ez_real_d, rho_data_d,
        nsites, nx, ny, nz, nhalo, norm);
    cudaDeviceSynchronize();

    /* Real-space force from particles (same as ewald_charge_sum_full_gpu) */
    if (nparticles > 0) {
      int blocks_nodes = (n_real + threads - 1) / threads;
      ewald_force_real_particle_kernel << <blocks_nodes, threads >> > (
          force_d, rho_data_d, particle_r_d, particle_q_d,
          nparticles, nsites, nhalo, irc);
      cudaDeviceSynchronize();
    }

    /* Electric field output (if requested) */
    if (efield_fourier_d != NULL) {
      /* K8: E_fourier = -grad(phi) + E_dipole = -(IFFT[+ik*psi_hat]/N) + E_dipole */
      ewald_fft_efield_fourier_kernel << <blocks_r, threads >> > (
          efield_fourier_d, Ex_real_d, Ey_real_d, Ez_real_d, nx, ny, nz, nhalo, norm);
      cudaDeviceSynchronize();

      /* Copy Fourier efield to total efield (real part added below) */
      cudaMemcpy(efield_d, efield_fourier_d, 3 * nsites * sizeof(double), cudaMemcpyDeviceToDevice);

      /* Real-space efield from particles */
      if (nparticles > 0) {
        ewald_efield_real_particle_kernel << <blocks_r, threads >> > (
            efield_d, particle_r_d, particle_q_d, nparticles, nsites, nhalo, irc);
        cudaDeviceSynchronize();
        ewald_efield_real_particle_kernel << <blocks_r, threads >> > (
            efield_real_d, particle_r_d, particle_q_d, nparticles, nsites, nhalo, irc);
        cudaDeviceSynchronize();
      }

      /* Real-space efield from lattice nodes */
      ewald_efield_real_lattice_kernel << <blocks_r, threads >> > (
          efield_d, rho_data_d, nsites, nhalo, irc);
      cudaDeviceSynchronize();
      ewald_efield_real_lattice_kernel << <blocks_r, threads >> > (
          efield_real_d, rho_data_d, nsites, nhalo, irc);
      cudaDeviceSynchronize();
    }

    /* Copy force and efield back to host */
    double* force_h = (double*)malloc(3 * nsites * sizeof(double));
    cudaMemcpy(force_h, force_d, 3 * nsites * sizeof(double), cudaMemcpyDeviceToHost);

    double* efield_h = NULL;
    double* efield_fourier_h = NULL;
    double* efield_real_h = NULL;

    if (efield_d != NULL) {
      efield_h = (double*)malloc(3 * nsites * sizeof(double));
      cudaMemcpy(efield_h, efield_d, 3 * nsites * sizeof(double), cudaMemcpyDeviceToHost);
      pe_info(ewald->pe, "  Electric field (FFT) computed and copied to host\n");

      if (efield_fourier_d != NULL) {
        efield_fourier_h = (double*)malloc(3 * nsites * sizeof(double));
        cudaMemcpy(efield_fourier_h, efield_fourier_d, 3 * nsites * sizeof(double), cudaMemcpyDeviceToHost);
      }
      if (efield_real_d != NULL) {
        efield_real_h = (double*)malloc(3 * nsites * sizeof(double));
        cudaMemcpy(efield_real_h, efield_real_d, 3 * nsites * sizeof(double), cudaMemcpyDeviceToHost);
      }
    }

    double e_sum[3] = { 0.0, 0.0, 0.0 };
    int e_count = 0;

    for (int ic = 1; ic <= nlocal[X]; ic++) {
      for (int jc = 1; jc <= nlocal[Y]; jc++) {
        for (int kc = 1; kc <= nlocal[Z]; kc++) {
          int index = cs_index(ewald->cs, ic, jc, kc);
          double f[3] = { force_h[3 * index + 0], force_h[3 * index + 1], force_h[3 * index + 2] };
          hydro_f_local_add(ewald->hydro, index, f);

          if (efield_h != NULL) {
            double e_field[3] = { efield_h[3 * index + 0], efield_h[3 * index + 1], efield_h[3 * index + 2] };
            field_vector_set(ewald->psi->efield, index, e_field);
            e_sum[X] += fabs(e_field[X]);
            e_sum[Y] += fabs(e_field[Y]);
            e_sum[Z] += fabs(e_field[Z]);
            e_count++;

            if (efield_fourier_h != NULL && ewald->psi->efield_fourier) {
              double e_fourier[3] = { efield_fourier_h[3 * index + 0],
                                     efield_fourier_h[3 * index + 1],
                                     efield_fourier_h[3 * index + 2] };
              field_vector_set(ewald->psi->efield_fourier, index, e_fourier);
            }
            if (efield_real_h != NULL && ewald->psi->efield_real) {
              double e_real[3] = { efield_real_h[3 * index + 0],
                                  efield_real_h[3 * index + 1],
                                  efield_real_h[3 * index + 2] };
              field_vector_set(ewald->psi->efield_real, index, e_real);
            }
          }

          kahan_add_double(&F_fluid_total[X], f[X]);
          kahan_add_double(&F_fluid_total[Y], f[Y]);
          kahan_add_double(&F_fluid_total[Z], f[Z]);
        }
      }
    }

    if (efield_h != NULL && e_count > 0) {
      pe_info(ewald->pe, "  Electric field sum |E|: (%14.7e, %14.7e, %14.7e) over %d sites\n",
              e_sum[X] / e_count, e_sum[Y] / e_count, e_sum[Z] / e_count, e_count);
    }

    if (efield_h != NULL) free(efield_h);
    if (efield_fourier_h != NULL) free(efield_fourier_h);
    if (efield_real_h != NULL) free(efield_real_h);
    if (efield_d != NULL) cudaFree(efield_d);
    if (efield_fourier_d != NULL) cudaFree(efield_fourier_d);
    if (efield_real_d != NULL) cudaFree(efield_real_d);

    if (ewald->psi && ewald->psi->efield) {
      field_memcpy(ewald->psi->efield, tdpMemcpyHostToDevice);
    }
    if (ewald->psi && ewald->psi->efield_fourier) {
      field_memcpy(ewald->psi->efield_fourier, tdpMemcpyHostToDevice);
    }
    if (ewald->psi && ewald->psi->efield_real) {
      field_memcpy(ewald->psi->efield_real, tdpMemcpyHostToDevice);
    }

    free(force_h);
    cudaFree(force_d);

    hydro_memcpy(ewald->hydro, tdpMemcpyHostToDevice);
  }

  /* MPI reduce F_fluid_total */
  {
    MPI_Datatype kahan_dt;
    MPI_Op kahan_op;
    kahan_mpi_datatype(&kahan_dt);
    kahan_mpi_op_sum(&kahan_op);
    MPI_Allreduce(MPI_IN_PLACE, F_fluid_total, 3, kahan_dt, kahan_op, comm);
    MPI_Type_free(&kahan_dt);
    MPI_Op_free(&kahan_op);
  }

  /* ========================================================================
   * [7] Field and force on particles
   *
   *  Fourier part uses K7: direct sum over rho_complex_d (same as
   *  ewald_particle_field_fourier_kernel but using cuFFT array instead of
   *  Sk_sin/Sk_cos arrays).
   *
   *  The K7 kernel uses rho_complex_d where:
   *    rho_hat.x = Sk_cos, rho_hat.y = -Sk_sin
   *
   *  Comparing with ewald_particle_field_fourier_kernel:
   *    im_part = Sk_sin*cos(kr) - Sk_cos*sin(kr)
   *            = -rho_hat.y*cos(kr) - rho_hat.x*sin(kr)
   *    E_alpha += factor * beta*eunit * Gk * k_alpha * im_part
   *
   *  K7 formula:
   *    psi_re = G*rho_hat.x = G*Sk_cos
   *    psi_im = G*rho_hat.y = -G*Sk_sin
   *    E_contrib = -(psi_im*cos + psi_re*sin)
   *              = -((-G*Sk_sin)*cos + G*Sk_cos*sin)
   *              = G*(Sk_sin*cos - Sk_cos*sin) = G*im_part
   *    E_alpha += herm*norm * k_alpha * E_contrib = herm*norm * k_alpha * G*im_part
   *
   *  In ewald_particle_field_fourier_kernel:
   *    E_alpha += factor * beta*eunit * Gk * k_alpha * im_part
   *  where Gk = b0*exp(-r4a*k^2)/k^2, norm = 1/N = 1/(Lx*Ly*Lz).
   *  Note: b0 = 1/(V*epsilon), so Gk*norm = b0/k^2*exp(...)/N = G_ewald.
   *  Also factor beta*eunit is applied inside K7 for Esub.
   *
   *  Real-space from ewald_particle_field_real_kernel (unchanged).
   * ======================================================================== */

  pe_info(ewald->pe, "  [6/6] Computing field and force on particles...\n");

  kahan_t F_particle_total[3] = { kahan_zero(), kahan_zero(), kahan_zero() };

  if (nparticles > 0) {
    int p_blocks = (nparticles + threads - 1) / threads;

    /* Fourier part: K7 (equivalent to ewald_particle_field_fourier_kernel) */
    ewald_particle_field_fft_kernel << <p_blocks, threads >> > (
        Esub_d, fex_d, particle_r_d, particle_q_d,
        rho_complex_d, G_ewald_d, nparticles,
        fkx, fky, fkz, nx, ny, nz_complex, norm);
    cudaDeviceSynchronize();

    /* Real-space from lattice and other particles (same as direct version) */
    ewald_particle_field_real_kernel << <p_blocks, threads >> > (
        Esub_d, fex_d, particle_r_d, particle_q_d, rho_data_d,
        nparticles, nsites, nhalo, irc);
    cudaDeviceSynchronize();

    /* Copy results back to host and store in colloid structures */
    double* Esub_h = (double*)malloc(3 * nparticles * sizeof(double));
    double* fex_h = (double*)malloc(3 * nparticles * sizeof(double));
    cudaMemcpy(Esub_h, Esub_d, 3 * nparticles * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(fex_h, fex_d, 3 * nparticles * sizeof(double), cudaMemcpyDeviceToHost);

    int p = 0;
    colloid_t* pc;
    colloids_info_local_head(ewald->cinfo, &pc);
    for (; pc; pc = pc->nextlocal) {
      pc->Esub[X] = Esub_h[3 * p + 0];
      pc->Esub[Y] = Esub_h[3 * p + 1];
      pc->Esub[Z] = Esub_h[3 * p + 2];
      pc->fex[X] = fex_h[3 * p + 0];
      pc->fex[Y] = fex_h[3 * p + 1];
      pc->fex[Z] = fex_h[3 * p + 2];

      kahan_add_double(&F_particle_total[X], pc->fex[X]);
      kahan_add_double(&F_particle_total[Y], pc->fex[Y]);
      kahan_add_double(&F_particle_total[Z], pc->fex[Z]);
      p++;
    }

    free(Esub_h);
    free(fex_h);
  }

  /* MPI reduce particle forces */
  {
    MPI_Datatype kahan_dt;
    MPI_Op kahan_op;
    kahan_mpi_datatype(&kahan_dt);
    kahan_mpi_op_sum(&kahan_op);
    MPI_Allreduce(MPI_IN_PLACE, F_particle_total, 3, kahan_dt, kahan_op, comm);
    MPI_Type_free(&kahan_dt);
    MPI_Op_free(&kahan_op);
  }

  double F_fluid[3] = { kahan_sum(&F_fluid_total[X]),    kahan_sum(&F_fluid_total[Y]),    kahan_sum(&F_fluid_total[Z]) };
  double F_particle[3] = { kahan_sum(&F_particle_total[X]), kahan_sum(&F_particle_total[Y]), kahan_sum(&F_particle_total[Z]) };
  double F_diff[3] = { F_fluid[X] + F_particle[X], F_fluid[Y] + F_particle[Y], F_fluid[Z] + F_particle[Z] };

  fprintf(fp, "%1.15e,%1.15e, %1.15e, %1.15e,%1.15e, %1.15e, %1.15e,%1.15e, %1.15e, %1.15e,\n",
          sqrt(F_diff[X] * F_diff[X] + F_diff[Y] * F_diff[Y] + F_diff[Z] * F_diff[Z]),
          F_diff[X], F_diff[Y], F_diff[Z],
          F_fluid[X], F_fluid[Y], F_fluid[Z],
          F_particle[X], F_particle[Y], F_particle[Z]);

  pe_info(ewald->pe, "\n");
  pe_info(ewald->pe, "  Momentum conservation check:\n");
  pe_info(ewald->pe, "    F_fluid    = (%14.7e, %14.7e, %14.7e)\n", F_fluid[X], F_fluid[Y], F_fluid[Z]);
  pe_info(ewald->pe, "    F_particle = (%14.7e, %14.7e, %14.7e)\n", F_particle[X], F_particle[Y], F_particle[Z]);
  pe_info(ewald->pe, "    F_total    = (%14.7e, %14.7e, %14.7e)\n", F_diff[X], F_diff[Y], F_diff[Z]);
  pe_info(ewald->pe, "    |F_total|  = %14.7e\n", sqrt(F_diff[X] * F_diff[X] + F_diff[Y] * F_diff[Y] + F_diff[Z] * F_diff[Z]));

  ewald_charge_external_field(ewald);

  /* ========================================================================
   * Cleanup
   * ======================================================================== */

  cufftDestroy(plan_fwd);
  cufftDestroy(plan_bwd);

  cudaFree(rho_real_d);
  cudaFree(rho_complex_d);
  cudaFree(psi_complex_d);
  cudaFree(efield_complex_d);
  cudaFree(psi_real_d);
  cudaFree(Ex_real_d);
  cudaFree(Ey_real_d);
  cudaFree(Ez_real_d);
  cudaFree(G_ewald_d);

  if (nparticles > 0) {
    free(particle_r_h);
    free(particle_q_h);
    cudaFree(particle_r_d);
    cudaFree(particle_q_d);
    cudaFree(Esub_d);
    cudaFree(fex_d);
  }

  pe_info(ewald->pe, "Ewald charge sum FFT (GPU): complete.\n\n");

  TIMER_stop(TIMER_EWALD_TOTAL);

  return 0;
}

/*CHANGE END - 20260212*/

/*CHANGE INIT - 20260302 Peskin-spread FFT Ewald version */

/*****************************************************************************
 *
 *  ewald_peskin_spread_kernel
 *
 *  Spread one charge source onto the FFT grid using the Peskin delta kernel.
 *  Each thread handles one source (lattice node or particle).
 *
 *  r_src: global 1-based coordinates.
 *  rho_real: FFT array (0-based, no halo, size nx*ny*nz).
 *
 *****************************************************************************/

__global__ void ewald_peskin_spread_kernel(
    double* __restrict__ rho_real,
    const double* __restrict__ r_src,
    const double* __restrict__ q_src,
    int nsrc,
    int nx, int ny, int nz) {

  int s = blockIdx.x * blockDim.x + threadIdx.x;
  if (s >= nsrc) return;

  double q = q_src[s];
  if (fabs(q) < 1.0e-14) return;

  /* Convert global 1-based to local 1-based (subtract noffset) */
  double r0x = r_src[3 * s + 0] - (double)d_noffset[0];
  double r0y = r_src[3 * s + 1] - (double)d_noffset[1];
  double r0z = r_src[3 * s + 2] - (double)d_noffset[2];

  /* Peskin support |r|<=2: window floor(r0)-1 .. floor(r0)+2 (1-based local) */
  int i_min = (int)floor(r0x) - 1;
  int i_max = (int)floor(r0x) + 2;
  int j_min = (int)floor(r0y) - 1;
  int j_max = (int)floor(r0y) + 2;
  int k_min = (int)floor(r0z) - 1;
  int k_max = (int)floor(r0z) + 2;

  if (i_min < 1) i_min = 1;  if (i_max > nx) i_max = nx;
  if (j_min < 1) j_min = 1;  if (j_max > ny) j_max = ny;
  if (k_min < 1) k_min = 1;  if (k_max > nz) k_max = nz;

  for (int i = i_min; i <= i_max; i++) {
    double wx = ewald_peskin_1d(r0x - (double)i);
    if (wx == 0.0) continue;
    for (int j = j_min; j <= j_max; j++) {
      double wy = ewald_peskin_1d(r0y - (double)j);
      if (wy == 0.0) continue;
      for (int k = k_min; k <= k_max; k++) {
        double wz = ewald_peskin_1d(r0z - (double)k);
        if (wz == 0.0) continue;
        int fft_idx = (i - 1) * ny * nz + (j - 1) * nz + (k - 1);
        atomicAdd(&rho_real[fft_idx], q * wx * wy * wz);
      }
    }
  }
}

/*****************************************************************************
 *
 *  ewald_peskin_interpolate_field_kernel
 *
 *  Interpolate E = -grad(phi) to particle position using the Peskin kernel.
 *  Ex_real, Ey_real, Ez_real contain IFFT[+ik*psi_hat] (= +grad(phi)).
 *  Negation gives E = -grad(phi). Dipole correction added.
 *
 *  Sets Esub = beta*eunit*E_physical and fex = q*E_physical.
 *
 *****************************************************************************/

__global__ void ewald_peskin_interpolate_field_kernel(
    double* __restrict__ Esub_data,
    double* __restrict__ fex_data,
    const double* __restrict__ Ex_real,
    const double* __restrict__ Ey_real,
    const double* __restrict__ Ez_real,
    const double* __restrict__ particle_r,
    const double* __restrict__ particle_q,
    int nparticles,
    int nx, int ny, int nz,
    double norm) {

  int p = blockIdx.x * blockDim.x + threadIdx.x;
  if (p >= nparticles) return;

  double q_p = particle_q[p];

  double r0x = particle_r[3 * p + 0] - (double)d_noffset[0];
  double r0y = particle_r[3 * p + 1] - (double)d_noffset[1];
  double r0z = particle_r[3 * p + 2] - (double)d_noffset[2];

  int i_min = (int)floor(r0x) - 1;
  int i_max = (int)floor(r0x) + 2;
  int j_min = (int)floor(r0y) - 1;
  int j_max = (int)floor(r0y) + 2;
  int k_min = (int)floor(r0z) - 1;
  int k_max = (int)floor(r0z) + 2;

  if (i_min < 1) i_min = 1;  if (i_max > nx) i_max = nx;
  if (j_min < 1) j_min = 1;  if (j_max > ny) j_max = ny;
  if (k_min < 1) k_min = 1;  if (k_max > nz) k_max = nz;

  double Ex = 0.0, Ey = 0.0, Ez = 0.0;

  for (int i = i_min; i <= i_max; i++) {
    double wx = ewald_peskin_1d(r0x - (double)i);
    if (wx == 0.0) continue;
    for (int j = j_min; j <= j_max; j++) {
      double wy = ewald_peskin_1d(r0y - (double)j);
      if (wy == 0.0) continue;
      for (int k = k_min; k <= k_max; k++) {
        double wz = ewald_peskin_1d(r0z - (double)k);
        if (wz == 0.0) continue;
        int fft_idx = (i - 1) * ny * nz + (j - 1) * nz + (k - 1);
        double w = wx * wy * wz;
        /* negate: IFFT[+ik*psi_hat]/N = +grad(phi), so E = -(...) */
        Ex -= w * Ex_real[fft_idx] * norm;
        Ey -= w * Ey_real[fft_idx] * norm;
        Ez -= w * Ez_real[fft_idx] * norm;
      }
    }
  }

  /* Dipole correction */
  Ex += d_E_dipole[0];
  Ey += d_E_dipole[1];
  Ez += d_E_dipole[2];

  /* Esub = beta * eunit * E_physical */
  Esub_data[3 * p + 0] = d_beta * d_eunit * Ex;
  Esub_data[3 * p + 1] = d_beta * d_eunit * Ey;
  Esub_data[3 * p + 2] = d_beta * d_eunit * Ez;

  /* fex = q * E_physical */
  double kt = 1.0 / d_beta;
  fex_data[3 * p + 0] = q_p * Esub_data[3 * p + 0] * kt / d_eunit;
  fex_data[3 * p + 1] = q_p * Esub_data[3 * p + 1] * kt / d_eunit;
  fex_data[3 * p + 2] = q_p * Esub_data[3 * p + 2] * kt / d_eunit;
}

/*****************************************************************************
 *
 *  ewald_charge_sum_FFT_full_gpu_peskin
 *
 *  Same as ewald_charge_sum_FFT_full_gpu but ALL charges (lattice nodes AND
 *  particles) are spread to the FFT grid with the Peskin delta kernel before
 *  the forward FFT. No continuous-position structure factors (no K6).
 *
 *  The charge density on the FFT grid is:
 *    rho_fft(r_n) = sum_{nodos} q_n * peskin(r_n - r_nodo)
 *                 + sum_{particulas} q_p * peskin(r_n - r_p)
 *
 *  The Fourier force on each particle is obtained by interpolating the
 *  FFT-computed E field back to r_p with the same Peskin kernel.
 *  Real-space (erfc) corrections are unchanged.
 *
 *****************************************************************************/

int ewald_charge_sum_FFT_full_gpu_peskin(ewald_charge_t* ewald, FILE* fp) {

  int nlocal[3], noffset[3];
  double ltot[3];
  double fkx, fky, fkz;
  double r4alpha_sq, b0;
  int irc;
  PI_DOUBLE(pi);

  if (ewald == NULL) return 0;
  assert(fp);

  TIMER_start(TIMER_EWALD_TOTAL);

  cs_nlocal(ewald->cs, nlocal);
  cs_nlocal_offset(ewald->cs, noffset);
  cs_ltot(ewald->cs, ltot);

  int nhalo;
  cs_nhalo(ewald->cs, &nhalo);

  irc = (int)ceil(ewald_rc_);

  fkx = 2.0 * pi / ltot[X];
  fky = 2.0 * pi / ltot[Y];
  fkz = 2.0 * pi / ltot[Z];
  r4alpha_sq = 1.0 / (4.0 * alpha_ * alpha_);
  b0 = 1.0 / (ltot[X] * ltot[Y] * ltot[Z] * epsilon_);

  int nx = nlocal[X];
  int ny = nlocal[Y];
  int nz = nlocal[Z];
  int n_real = nx * ny * nz;
  int nz_complex = nz / 2 + 1;
  int n_complex = nx * ny * nz_complex;
  int nsites = ewald->psi->nsites;
  double norm = 1.0 / (double)n_real;

  pe_info(ewald->pe, "\n");
  pe_info(ewald->pe, "Ewald charge sum FFT-Peskin (GPU): computing potential, forces, and fields...\n");
  pe_info(ewald->pe, "  Local nodes: %d x %d x %d = %d\n", nx, ny, nz, n_real);
  pe_info(ewald->pe, "  Real-space cutoff: %.2f (irc=%d)\n", ewald_rc_, irc);

  cudaMemcpyToSymbol(d_alpha, &alpha_, sizeof(double));
  cudaMemcpyToSymbol(d_eps_reg, &eps_reg_, sizeof(double));
  cudaMemcpyToSymbol(d_epsilon, &epsilon_, sizeof(double));
  cudaMemcpyToSymbol(d_beta, &beta_, sizeof(double));
  cudaMemcpyToSymbol(d_eunit, &eunit_, sizeof(double));
  cudaMemcpyToSymbol(d_rpi, &rpi_, sizeof(double));
  cudaMemcpyToSymbol(d_ewald_rc, &ewald_rc_, sizeof(double));
  cudaMemcpyToSymbol(d_fkx, &fkx, sizeof(double));
  cudaMemcpyToSymbol(d_fky, &fky, sizeof(double));
  cudaMemcpyToSymbol(d_fkz, &fkz, sizeof(double));
  cudaMemcpyToSymbol(d_r4alpha_sq, &r4alpha_sq, sizeof(double));
  cudaMemcpyToSymbol(d_b0, &b0, sizeof(double));
  cudaMemcpyToSymbol(d_kmax, &kmax_, sizeof(double));
  cudaMemcpyToSymbol(d_nk, nk_, 3 * sizeof(int));
  cudaMemcpyToSymbol(d_nlocal, nlocal, 3 * sizeof(int));
  cudaMemcpyToSymbol(d_noffset, noffset, 3 * sizeof(int));

  /* [1] G_ewald(k) */
  pe_info(ewald->pe, "  [1/6] Precomputing Ewald Green function G(k)...\n");

  double* G_ewald_h = (double*)malloc(n_complex * sizeof(double));
  for (int ix = 0; ix < nx; ix++) {
    int kx_idx = (ix <= nx / 2) ? ix : ix - nx;
    double kx = fkx * kx_idx;
    for (int iy = 0; iy < ny; iy++) {
      int ky_idx = (iy <= ny / 2) ? iy : iy - ny;
      double ky = fky * ky_idx;
      for (int iz = 0; iz < nz_complex; iz++) {
        double kz = fkz * iz;
        int idx = ix * ny * nz_complex + iy * nz_complex + iz;
        if (ix == 0 && iy == 0 && iz == 0) { G_ewald_h[idx] = 0.0; continue; }
        double ksq = kx * kx + ky * ky + kz * kz;
        G_ewald_h[idx] = b0 * exp(-r4alpha_sq * ksq) / ksq;
      }
    }
  }
  double* G_ewald_d = NULL;
  cudaMalloc(&G_ewald_d, n_complex * sizeof(double));
  cudaMemcpy(G_ewald_d, G_ewald_h, n_complex * sizeof(double), cudaMemcpyHostToDevice);
  free(G_ewald_h);

  /* [2] Device pointers for rho and psi */
  double* rho_data_d = NULL;
  double* psi_data_d = NULL;
  size_t data_offset = offsetof(field_t, data);
  cudaMemcpy(&rho_data_d, (char*)(ewald->psi->rho->target) + data_offset,
             sizeof(double*), cudaMemcpyDeviceToHost);
  cudaMemcpy(&psi_data_d, (char*)(ewald->psi->psi->target) + data_offset,
             sizeof(double*), cudaMemcpyDeviceToHost);
  field_memcpy(ewald->psi->rho, tdpMemcpyHostToDevice);

  /* [3] Collect particle data */
  int nparticles = 0;
  double* particle_r_h = NULL, * particle_q_h = NULL;
  double* particle_r_d = NULL, * particle_q_d = NULL;
  double* Esub_d = NULL, * fex_d = NULL;

  if ((ewald->sources & EWALD_SOURCE_COLLOIDS) && ewald->cinfo) {
    colloid_t* pc;
    colloids_info_local_head(ewald->cinfo, &pc);
    for (; pc; pc = pc->nextlocal) nparticles++;
    if (nparticles > 0) {
      particle_r_h = (double*)malloc(3 * nparticles * sizeof(double));
      particle_q_h = (double*)malloc(nparticles * sizeof(double));
      int p = 0;
      colloids_info_local_head(ewald->cinfo, &pc);
      for (; pc; pc = pc->nextlocal) {
        particle_r_h[3 * p + 0] = pc->s.r[X];
        particle_r_h[3 * p + 1] = pc->s.r[Y];
        particle_r_h[3 * p + 2] = pc->s.r[Z];
        particle_q_h[p] = pc->s.q0 - pc->s.q1;
        p++;
      }
      cudaMalloc(&particle_r_d, 3 * nparticles * sizeof(double));
      cudaMalloc(&particle_q_d, nparticles * sizeof(double));
      cudaMalloc(&Esub_d, 3 * nparticles * sizeof(double));
      cudaMalloc(&fex_d, 3 * nparticles * sizeof(double));
      cudaMemcpy(particle_r_d, particle_r_h, 3 * nparticles * sizeof(double), cudaMemcpyHostToDevice);
      cudaMemcpy(particle_q_d, particle_q_h, nparticles * sizeof(double), cudaMemcpyHostToDevice);
    }
  }
  pe_info(ewald->pe, "  Particles: %d\n", nparticles);

  MPI_Comm comm;
  cs_cart_comm(ewald->cs, &comm);

  /* ========================================================================
   * [4] Build rho_fft by spreading ALL charges with Peskin kernel.
   *
   *  Lattice nodes: position (noffset+ic, noffset+jc, noffset+kc) global 1-based.
   *  Particles: position r_p global 1-based (continuous).
   *  Both use ewald_peskin_spread_kernel.
   * ======================================================================== */

  pe_info(ewald->pe, "  [2/6] Spreading charges onto FFT grid via Peskin kernel...\n");

  int threads = 256;
  int blocks_r = (n_real + threads - 1) / threads;
  int blocks_c = (n_complex + threads - 1) / threads;

  double* rho_real_d = NULL;
  cufftDoubleComplex* rho_complex_d = NULL;
  cudaMalloc(&rho_real_d, n_real * sizeof(double));
  cudaMemset(rho_real_d, 0, n_real * sizeof(double));
  cudaMalloc(&rho_complex_d, n_complex * sizeof(cufftDoubleComplex));
  cudaMemset(rho_complex_d, 0, n_complex * sizeof(cufftDoubleComplex));

  cufftHandle plan_fwd, plan_bwd;
  cufftPlan3d(&plan_fwd, nx, ny, nz, CUFFT_D2Z);
  cufftPlan3d(&plan_bwd, nx, ny, nz, CUFFT_Z2D);

  /* Spread lattice nodes */
  if (ewald->sources & EWALD_SOURCE_LATTICE) {
    field_memcpy(ewald->psi->rho, tdpMemcpyDeviceToHost);

    int nlattice = 0;
    double* lat_r_h = (double*)malloc(3 * n_real * sizeof(double));
    double* lat_q_h = (double*)malloc(n_real * sizeof(double));

    for (int ic = 1; ic <= nlocal[X]; ic++) {
      for (int jc = 1; jc <= nlocal[Y]; jc++) {
        for (int kc = 1; kc <= nlocal[Z]; kc++) {
          int index = cs_index(ewald->cs, ic, jc, kc);
          double rho_elec;
          psi_rho_elec(ewald->psi, index, &rho_elec);
          if (fabs(rho_elec) < 1.0e-14) continue;
          lat_r_h[3 * nlattice + 0] = (double)(noffset[X] + ic);
          lat_r_h[3 * nlattice + 1] = (double)(noffset[Y] + jc);
          lat_r_h[3 * nlattice + 2] = (double)(noffset[Z] + kc);
          lat_q_h[nlattice] = rho_elec;
          nlattice++;
        }
      }
    }
    pe_info(ewald->pe, "  Lattice sources: %d non-zero nodes\n", nlattice);

    if (nlattice > 0) {
      double* lat_r_d = NULL, * lat_q_d = NULL;
      cudaMalloc(&lat_r_d, 3 * nlattice * sizeof(double));
      cudaMalloc(&lat_q_d, nlattice * sizeof(double));
      cudaMemcpy(lat_r_d, lat_r_h, 3 * nlattice * sizeof(double), cudaMemcpyHostToDevice);
      cudaMemcpy(lat_q_d, lat_q_h, nlattice * sizeof(double), cudaMemcpyHostToDevice);
      int lat_blocks = (nlattice + threads - 1) / threads;
      ewald_peskin_spread_kernel << <lat_blocks, threads >> > (
          rho_real_d, lat_r_d, lat_q_d, nlattice, nx, ny, nz);
      cudaDeviceSynchronize();
      cudaFree(lat_r_d);
      cudaFree(lat_q_d);
    }
    free(lat_r_h);
    free(lat_q_h);
  }

  /* Spread particles (need all MPI ranks) */
  if ((ewald->sources & EWALD_SOURCE_COLLOIDS) && nparticles > 0) {
    int nparticles_local = nparticles;
    int nparticles_all = nparticles;
    double* all_particle_r_d = particle_r_d;
    double* all_particle_q_d = particle_q_d;
    double* all_particle_r_h = particle_r_h;
    double* all_particle_q_h = particle_q_h;

    int npe;
    MPI_Comm_size(comm, &npe);
    if (npe > 1) {
      int* counts = (int*)malloc(npe * sizeof(int));
      MPI_Allgather(&nparticles_local, 1, MPI_INT, counts, 1, MPI_INT, comm);
      nparticles_all = 0;
      for (int i = 0; i < npe; i++) nparticles_all += counts[i];
      if (nparticles_all > 0) {
        all_particle_r_h = (double*)malloc(3 * nparticles_all * sizeof(double));
        all_particle_q_h = (double*)malloc(nparticles_all * sizeof(double));
        int* displs_r = (int*)malloc(npe * sizeof(int));
        int* counts_r = (int*)malloc(npe * sizeof(int));
        displs_r[0] = 0; counts_r[0] = counts[0] * 3;
        for (int i = 1; i < npe; i++) {
          counts_r[i] = counts[i] * 3;
          displs_r[i] = displs_r[i - 1] + counts_r[i - 1];
        }
        MPI_Allgatherv(particle_r_h, nparticles_local * 3, MPI_DOUBLE,
                       all_particle_r_h, counts_r, displs_r, MPI_DOUBLE, comm);
        int* displs_q = (int*)malloc(npe * sizeof(int));
        displs_q[0] = 0;
        for (int i = 1; i < npe; i++) displs_q[i] = displs_q[i - 1] + counts[i - 1];
        MPI_Allgatherv(particle_q_h, nparticles_local, MPI_DOUBLE,
                       all_particle_q_h, counts, displs_q, MPI_DOUBLE, comm);
        cudaMalloc(&all_particle_r_d, 3 * nparticles_all * sizeof(double));
        cudaMalloc(&all_particle_q_d, nparticles_all * sizeof(double));
        cudaMemcpy(all_particle_r_d, all_particle_r_h, 3 * nparticles_all * sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy(all_particle_q_d, all_particle_q_h, nparticles_all * sizeof(double), cudaMemcpyHostToDevice);
        free(displs_r); free(counts_r); free(displs_q);
      }
      free(counts);
    }
    if (nparticles_all > 0) {
      int p_blocks = (nparticles_all + threads - 1) / threads;
      ewald_peskin_spread_kernel << <p_blocks, threads >> > (
          rho_real_d, all_particle_r_d, all_particle_q_d,
          nparticles_all, nx, ny, nz);
      cudaDeviceSynchronize();
    }
    if (npe > 1 && all_particle_r_d != particle_r_d) {
      if (all_particle_r_h != particle_r_h) free(all_particle_r_h);
      if (all_particle_q_h != particle_q_h) free(all_particle_q_h);
      cudaFree(all_particle_r_d);
      cudaFree(all_particle_q_d);
    }
  }

  /* FFT forward */
  cufftExecD2Z(plan_fwd, (cufftDoubleReal*)rho_real_d,
               (cufftDoubleComplex*)rho_complex_d);
  cudaDeviceSynchronize();

  /* MPI reduce rho_hat */
  {
    int npe;
    MPI_Comm_size(comm, &npe);
    if (npe > 1) {
      pe_info(ewald->pe, "  [3/6] MPI reducing structure factors...\n");
      double* complex_h = (double*)malloc(2 * n_complex * sizeof(double));
      cudaMemcpy(complex_h, rho_complex_d,
                 n_complex * sizeof(cufftDoubleComplex), cudaMemcpyDeviceToHost);
      MPI_Allreduce(MPI_IN_PLACE, complex_h, 2 * n_complex, MPI_DOUBLE, MPI_SUM, comm);
      cudaMemcpy(rho_complex_d, complex_h,
                 n_complex * sizeof(cufftDoubleComplex), cudaMemcpyHostToDevice);
      free(complex_h);
    }
  }

  /* [5] Dipole moment */
  kahan_t M_dipole[3] = { kahan_zero(), kahan_zero(), kahan_zero() };
  kahan_t Q_total_k = kahan_zero();

  if (nparticles > 0) {
    for (int p = 0; p < nparticles; p++) {
      double q = particle_q_h[p];
      kahan_add_double(&M_dipole[X], q * particle_r_h[3 * p + 0]);
      kahan_add_double(&M_dipole[Y], q * particle_r_h[3 * p + 1]);
      kahan_add_double(&M_dipole[Z], q * particle_r_h[3 * p + 2]);
      kahan_add_double(&Q_total_k, q);
    }
  }
  if ((ewald->sources & EWALD_SOURCE_LATTICE) && ewald->psi) {
    field_memcpy(ewald->psi->rho, tdpMemcpyDeviceToHost);
    for (int ic = 1; ic <= nlocal[X]; ic++) {
      for (int jc = 1; jc <= nlocal[Y]; jc++) {
        for (int kc = 1; kc <= nlocal[Z]; kc++) {
          int index = cs_index(ewald->cs, ic, jc, kc);
          double rho_elec;
          psi_rho_elec(ewald->psi, index, &rho_elec);
          kahan_add_double(&M_dipole[X], rho_elec * (double)(noffset[X] + ic));
          kahan_add_double(&M_dipole[Y], rho_elec * (double)(noffset[Y] + jc));
          kahan_add_double(&M_dipole[Z], rho_elec * (double)(noffset[Z] + kc));
          kahan_add_double(&Q_total_k, rho_elec);
        }
      }
    }
  }
  {
    kahan_t M_reduce[4] = { M_dipole[X], M_dipole[Y], M_dipole[Z], Q_total_k };
    MPI_Datatype kahan_dt; MPI_Op kahan_op;
    kahan_mpi_datatype(&kahan_dt); kahan_mpi_op_sum(&kahan_op);
    MPI_Allreduce(MPI_IN_PLACE, M_reduce, 4, kahan_dt, kahan_op, comm);
    MPI_Type_free(&kahan_dt); MPI_Op_free(&kahan_op);
    M_dipole[X] = M_reduce[0]; M_dipole[Y] = M_reduce[1];
    M_dipole[Z] = M_reduce[2]; Q_total_k = M_reduce[3];
  }
  double M[3] = { kahan_sum(&M_dipole[X]), kahan_sum(&M_dipole[Y]), kahan_sum(&M_dipole[Z]) };
  double Q_total = kahan_sum(&Q_total_k);
  double V = ltot[X] * ltot[Y] * ltot[Z];
  double dipole_prefactor = 4.0 * pi / ((1.0 + 2.0 * ewald->epsilon_prime) * V);
  double E_dipole[3] = { -dipole_prefactor * M[X], -dipole_prefactor * M[Y], -dipole_prefactor * M[Z] };

  pe_info(ewald->pe, "  Dipole moment M = (%14.7e, %14.7e, %14.7e)\n", M[X], M[Y], M[Z]);
  pe_info(ewald->pe, "  Total charge Q = %14.7e\n", Q_total);
  pe_info(ewald->pe, "  Dipole epsilon' = %14.7e\n", ewald->epsilon_prime);
  pe_info(ewald->pe, "  Dipole field E_dip = (%14.7e, %14.7e, %14.7e)\n",
          E_dipole[X], E_dipole[Y], E_dipole[Z]);

  cudaMemcpyToSymbol(d_E_dipole, E_dipole, 3 * sizeof(double));
  cudaMemcpyToSymbol(d_M_dipole, M, 3 * sizeof(double));
  double dipole_pref_dev = dipole_prefactor;
  cudaMemcpyToSymbol(d_dipole_prefactor, &dipole_pref_dev, sizeof(double));

  /* [6] FFT work arrays, potential, force, efield on lattice nodes */
  cufftDoubleComplex* psi_complex_d = NULL;
  cufftDoubleComplex* efield_complex_d = NULL;
  double* psi_real_d = NULL;
  double* Ex_real_d = NULL, * Ey_real_d = NULL, * Ez_real_d = NULL;

  cudaMalloc(&psi_complex_d, n_complex * sizeof(cufftDoubleComplex));
  cudaMalloc(&efield_complex_d, n_complex * sizeof(cufftDoubleComplex));
  cudaMalloc(&psi_real_d, n_real * sizeof(double));
  cudaMalloc(&Ex_real_d, n_real * sizeof(double));
  cudaMalloc(&Ey_real_d, n_real * sizeof(double));
  cudaMalloc(&Ez_real_d, n_real * sizeof(double));

  ewald_green_multiply_potential_kernel << <blocks_c, threads >> > (
      psi_complex_d, rho_complex_d, G_ewald_d, n_complex);
  cudaDeviceSynchronize();

  /* Potential */
  pe_info(ewald->pe, "  [4/6] Computing potential on lattice nodes...\n");
  cufftExecZ2D(plan_bwd, psi_complex_d, psi_real_d);
  cudaDeviceSynchronize();
  ewald_copy_psi_from_fft_kernel << <blocks_r, threads >> > (
      psi_data_d, psi_real_d, nx, ny, nz, nhalo, norm);
  cudaDeviceSynchronize();
  if (nparticles > 0) {
    int blocks_nodes = (n_real + threads - 1) / threads;
    ewald_potential_real_particle_kernel << <blocks_nodes, threads >> > (
        psi_data_d, particle_r_d, particle_q_d, nparticles, nsites, nhalo, irc);
    cudaDeviceSynchronize();
  }
  field_memcpy(ewald->psi->psi, tdpMemcpyDeviceToHost);

  /* Force and efield on lattice nodes */
  pe_info(ewald->pe, "  [5/6] Computing force on lattice nodes...\n");
  kahan_t F_fluid_total[3] = { kahan_zero(), kahan_zero(), kahan_zero() };

  if (ewald->hydro != NULL) {
    hydro_memcpy(ewald->hydro, tdpMemcpyDeviceToHost);

    double* efield_d = NULL, * efield_fourier_d = NULL, * efield_real_d = NULL;
    double* force_d = NULL;
    cudaMalloc(&force_d, 3 * nsites * sizeof(double));
    cudaMemset(force_d, 0, 3 * nsites * sizeof(double));

    if (ewald->psi && ewald->psi->efield) {
      cudaMalloc(&efield_d, 3 * nsites * sizeof(double)); cudaMemset(efield_d, 0, 3 * nsites * sizeof(double));
      cudaMalloc(&efield_fourier_d, 3 * nsites * sizeof(double)); cudaMemset(efield_fourier_d, 0, 3 * nsites * sizeof(double));
      cudaMalloc(&efield_real_d, 3 * nsites * sizeof(double)); cudaMemset(efield_real_d, 0, 3 * nsites * sizeof(double));
    }

    ewald_efield_multiply_kernel << <blocks_c, threads >> > (efield_complex_d, psi_complex_d, 0, fkx, fky, fkz, nx, ny, nz_complex);
    cudaDeviceSynchronize(); cufftExecZ2D(plan_bwd, efield_complex_d, Ex_real_d); cudaDeviceSynchronize();

    ewald_efield_multiply_kernel << <blocks_c, threads >> > (efield_complex_d, psi_complex_d, 1, fkx, fky, fkz, nx, ny, nz_complex);
    cudaDeviceSynchronize(); cufftExecZ2D(plan_bwd, efield_complex_d, Ey_real_d); cudaDeviceSynchronize();

    ewald_efield_multiply_kernel << <blocks_c, threads >> > (efield_complex_d, psi_complex_d, 2, fkx, fky, fkz, nx, ny, nz_complex);
    cudaDeviceSynchronize(); cufftExecZ2D(plan_bwd, efield_complex_d, Ez_real_d); cudaDeviceSynchronize();

    ewald_copy_force_from_fft_kernel << <blocks_r, threads >> > (
        force_d, Ex_real_d, Ey_real_d, Ez_real_d, rho_data_d,
        nsites, nx, ny, nz, nhalo, norm);
    cudaDeviceSynchronize();

    if (nparticles > 0) {
      int blocks_nodes = (n_real + threads - 1) / threads;
      ewald_force_real_particle_kernel << <blocks_nodes, threads >> > (
          force_d, rho_data_d, particle_r_d, particle_q_d,
          nparticles, nsites, nhalo, irc);
      cudaDeviceSynchronize();
    }

    if (efield_fourier_d != NULL) {
      ewald_fft_efield_fourier_kernel << <blocks_r, threads >> > (
          efield_fourier_d, Ex_real_d, Ey_real_d, Ez_real_d, nx, ny, nz, nhalo, norm);
      cudaDeviceSynchronize();
      cudaMemcpy(efield_d, efield_fourier_d, 3 * nsites * sizeof(double), cudaMemcpyDeviceToDevice);
      if (nparticles > 0) {
        ewald_efield_real_particle_kernel << <blocks_r, threads >> > (efield_d, particle_r_d, particle_q_d, nparticles, nsites, nhalo, irc); cudaDeviceSynchronize();
        ewald_efield_real_particle_kernel << <blocks_r, threads >> > (efield_real_d, particle_r_d, particle_q_d, nparticles, nsites, nhalo, irc); cudaDeviceSynchronize();
      }
      ewald_efield_real_lattice_kernel << <blocks_r, threads >> > (efield_d, rho_data_d, nsites, nhalo, irc); cudaDeviceSynchronize();
      ewald_efield_real_lattice_kernel << <blocks_r, threads >> > (efield_real_d, rho_data_d, nsites, nhalo, irc); cudaDeviceSynchronize();
    }

    double* force_h = (double*)malloc(3 * nsites * sizeof(double));
    cudaMemcpy(force_h, force_d, 3 * nsites * sizeof(double), cudaMemcpyDeviceToHost);

    double* efield_h = NULL, * efield_fourier_h = NULL, * efield_real_h = NULL;
    if (efield_d != NULL) {
      efield_h = (double*)malloc(3 * nsites * sizeof(double));
      cudaMemcpy(efield_h, efield_d, 3 * nsites * sizeof(double), cudaMemcpyDeviceToHost);
      pe_info(ewald->pe, "  Electric field (FFT-Peskin) computed and copied to host\n");
      if (efield_fourier_d != NULL) { efield_fourier_h = (double*)malloc(3 * nsites * sizeof(double)); cudaMemcpy(efield_fourier_h, efield_fourier_d, 3 * nsites * sizeof(double), cudaMemcpyDeviceToHost); }
      if (efield_real_d != NULL) { efield_real_h = (double*)malloc(3 * nsites * sizeof(double)); cudaMemcpy(efield_real_h, efield_real_d, 3 * nsites * sizeof(double), cudaMemcpyDeviceToHost); }
    }

    double e_sum[3] = { 0.0,0.0,0.0 }; int e_count = 0;
    for (int ic = 1; ic <= nlocal[X]; ic++) {
      for (int jc = 1; jc <= nlocal[Y]; jc++) {
        for (int kc = 1; kc <= nlocal[Z]; kc++) {
          int index = cs_index(ewald->cs, ic, jc, kc);
          double f[3] = { force_h[3 * index + 0], force_h[3 * index + 1], force_h[3 * index + 2] };
          hydro_f_local_add(ewald->hydro, index, f);
          if (efield_h != NULL) {
            double e_field[3] = { efield_h[3 * index + 0], efield_h[3 * index + 1], efield_h[3 * index + 2] };
            field_vector_set(ewald->psi->efield, index, e_field);
            e_sum[X] += fabs(e_field[X]); e_sum[Y] += fabs(e_field[Y]); e_sum[Z] += fabs(e_field[Z]); e_count++;
            if (efield_fourier_h != NULL && ewald->psi->efield_fourier) {
              double ef[3] = { efield_fourier_h[3 * index + 0], efield_fourier_h[3 * index + 1], efield_fourier_h[3 * index + 2] };
              field_vector_set(ewald->psi->efield_fourier, index, ef);
            }
            if (efield_real_h != NULL && ewald->psi->efield_real) {
              double er[3] = { efield_real_h[3 * index + 0], efield_real_h[3 * index + 1], efield_real_h[3 * index + 2] };
              field_vector_set(ewald->psi->efield_real, index, er);
            }
          }
          kahan_add_double(&F_fluid_total[X], f[X]);
          kahan_add_double(&F_fluid_total[Y], f[Y]);
          kahan_add_double(&F_fluid_total[Z], f[Z]);
        }
      }
    }
    if (efield_h != NULL && e_count > 0)
      pe_info(ewald->pe, "  Electric field sum |E|: (%14.7e, %14.7e, %14.7e) over %d sites\n",
              e_sum[X] / e_count, e_sum[Y] / e_count, e_sum[Z] / e_count, e_count);

    if (efield_h != NULL) free(efield_h);
    if (efield_fourier_h != NULL) free(efield_fourier_h);
    if (efield_real_h != NULL) free(efield_real_h);
    if (efield_d != NULL) cudaFree(efield_d);
    if (efield_fourier_d != NULL) cudaFree(efield_fourier_d);
    if (efield_real_d != NULL) cudaFree(efield_real_d);
    if (ewald->psi && ewald->psi->efield) field_memcpy(ewald->psi->efield, tdpMemcpyHostToDevice);
    if (ewald->psi && ewald->psi->efield_fourier) field_memcpy(ewald->psi->efield_fourier, tdpMemcpyHostToDevice);
    if (ewald->psi && ewald->psi->efield_real) field_memcpy(ewald->psi->efield_real, tdpMemcpyHostToDevice);
    free(force_h); cudaFree(force_d);
    hydro_memcpy(ewald->hydro, tdpMemcpyHostToDevice);
  }

  {
    MPI_Datatype kahan_dt; MPI_Op kahan_op;
    kahan_mpi_datatype(&kahan_dt); kahan_mpi_op_sum(&kahan_op);
    MPI_Allreduce(MPI_IN_PLACE, F_fluid_total, 3, kahan_dt, kahan_op, comm);
    MPI_Type_free(&kahan_dt); MPI_Op_free(&kahan_op);
  }

  /* ========================================================================
   * [7] Force on particles: Peskin interpolation of FFT field + real-space
   * ======================================================================== */

  pe_info(ewald->pe, "  [6/6] Computing field and force on particles (Peskin interpolation)...\n");
  kahan_t F_particle_total[3] = { kahan_zero(), kahan_zero(), kahan_zero() };

  if (nparticles > 0) {
    int p_blocks = (nparticles + threads - 1) / threads;

    ewald_peskin_interpolate_field_kernel << <p_blocks, threads >> > (
        Esub_d, fex_d, Ex_real_d, Ey_real_d, Ez_real_d,
        particle_r_d, particle_q_d, nparticles, nx, ny, nz, norm);
    cudaDeviceSynchronize();

    ewald_particle_field_real_kernel << <p_blocks, threads >> > (
        Esub_d, fex_d, particle_r_d, particle_q_d, rho_data_d,
        nparticles, nsites, nhalo, irc);
    cudaDeviceSynchronize();

    double* Esub_h = (double*)malloc(3 * nparticles * sizeof(double));
    double* fex_h = (double*)malloc(3 * nparticles * sizeof(double));
    cudaMemcpy(Esub_h, Esub_d, 3 * nparticles * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(fex_h, fex_d, 3 * nparticles * sizeof(double), cudaMemcpyDeviceToHost);

    int p = 0;
    colloid_t* pc;
    colloids_info_local_head(ewald->cinfo, &pc);
    for (; pc; pc = pc->nextlocal) {
      pc->Esub[X] = Esub_h[3 * p + 0]; pc->Esub[Y] = Esub_h[3 * p + 1]; pc->Esub[Z] = Esub_h[3 * p + 2];
      pc->fex[X] = fex_h[3 * p + 0];  pc->fex[Y] = fex_h[3 * p + 1];  pc->fex[Z] = fex_h[3 * p + 2];
      kahan_add_double(&F_particle_total[X], pc->fex[X]);
      kahan_add_double(&F_particle_total[Y], pc->fex[Y]);
      kahan_add_double(&F_particle_total[Z], pc->fex[Z]);
      p++;
    }
    free(Esub_h); free(fex_h);
  }

  {
    MPI_Datatype kahan_dt; MPI_Op kahan_op;
    kahan_mpi_datatype(&kahan_dt); kahan_mpi_op_sum(&kahan_op);
    MPI_Allreduce(MPI_IN_PLACE, F_particle_total, 3, kahan_dt, kahan_op, comm);
    MPI_Type_free(&kahan_dt); MPI_Op_free(&kahan_op);
  }

  double F_fluid[3] = { kahan_sum(&F_fluid_total[X]),    kahan_sum(&F_fluid_total[Y]),    kahan_sum(&F_fluid_total[Z]) };
  double F_particle[3] = { kahan_sum(&F_particle_total[X]), kahan_sum(&F_particle_total[Y]), kahan_sum(&F_particle_total[Z]) };
  double F_diff[3] = { F_fluid[X] + F_particle[X], F_fluid[Y] + F_particle[Y], F_fluid[Z] + F_particle[Z] };

  fprintf(fp, "%1.15e,%1.15e, %1.15e, %1.15e,%1.15e, %1.15e, %1.15e,%1.15e, %1.15e, %1.15e,\n",
          sqrt(F_diff[X] * F_diff[X] + F_diff[Y] * F_diff[Y] + F_diff[Z] * F_diff[Z]),
          F_diff[X], F_diff[Y], F_diff[Z], F_fluid[X], F_fluid[Y], F_fluid[Z],
          F_particle[X], F_particle[Y], F_particle[Z]);

  pe_info(ewald->pe, "\n");
  pe_info(ewald->pe, "  Momentum conservation check:\n");
  pe_info(ewald->pe, "    F_fluid    = (%14.7e, %14.7e, %14.7e)\n", F_fluid[X], F_fluid[Y], F_fluid[Z]);
  pe_info(ewald->pe, "    F_particle = (%14.7e, %14.7e, %14.7e)\n", F_particle[X], F_particle[Y], F_particle[Z]);
  pe_info(ewald->pe, "    F_total    = (%14.7e, %14.7e, %14.7e)\n", F_diff[X], F_diff[Y], F_diff[Z]);
  pe_info(ewald->pe, "    |F_total|  = %14.7e\n", sqrt(F_diff[X] * F_diff[X] + F_diff[Y] * F_diff[Y] + F_diff[Z] * F_diff[Z]));

  ewald_charge_external_field(ewald);

  /* Cleanup */
  cufftDestroy(plan_fwd); cufftDestroy(plan_bwd);
  cudaFree(rho_real_d); cudaFree(rho_complex_d);
  cudaFree(psi_complex_d); cudaFree(efield_complex_d);
  cudaFree(psi_real_d); cudaFree(Ex_real_d); cudaFree(Ey_real_d); cudaFree(Ez_real_d);
  cudaFree(G_ewald_d);
  if (nparticles > 0) {
    free(particle_r_h); free(particle_q_h);
    cudaFree(particle_r_d); cudaFree(particle_q_d);
    cudaFree(Esub_d); cudaFree(fex_d);
  }

  pe_info(ewald->pe, "Ewald charge sum FFT-Peskin (GPU): complete.\n\n");
  TIMER_stop(TIMER_EWALD_TOTAL);
  return 0;
}

/*CHANGE END - 20260302 Peskin-spread FFT Ewald version */

#else /* __NVCC__ not defined */

 /*****************************************************************************
  *
  *  ewald_charge_sum_full_gpu (stub for non-CUDA builds)
  *
  *****************************************************************************/

int ewald_charge_sum_full_gpu(ewald_charge_t* ewald, FILE* fp,
                              ewald_self_table_t* self_table) {

  (void)self_table; /* unused in CPU stub */

  if (ewald && ewald->pe) {
    pe_info(ewald->pe, "ewald_charge_sum_full_gpu: CUDA not available, using CPU version.\n");
  }

  return ewald_charge_sum_full(ewald, fp);
}

/*CHANGE INIT - 20260212 New ewald_charge_sum_FFT_full_gpu using cuFFT for Fourier-space*/
int ewald_charge_sum_FFT_full_gpu(ewald_charge_t* ewald, FILE* fp) {

  if (ewald && ewald->pe) {
    pe_info(ewald->pe, "ewald_charge_sum_FFT_full_gpu: CUDA not available, using CPU version.\n");
  }

  return ewald_charge_sum_full(ewald, fp);
}
/*CHANGE END - 20260212*/

#endif /* __NVCC__ */

/* ========================================================================
 * Self-field correction for lattice-discretisation artefact
 * ========================================================================
 *
 * When a Debye cloud is discretised onto a cubic lattice it loses its
 * perfect spherical symmetry.  The resulting spurious electric field at
 * the particle position (E_self) must be subtracted from the Ewald field
 * so that the net force on the particle is not biased by the grid.
 *
 * Strategy
 * ---------
 *  1. Sweep fractional positions (xf, yf, zf) over one unit cell on a
 *     regular Nx x Ny x Nz sub-grid.
 *  2. For each fractional position place the particle at r_p = (xf, yf, zf)
 *     in a periodic box of side L.
 *  3. Fill every lattice node with the Debye-Hückel charge density
 *        rho(r) = A * exp(-kappa * |r - r_p|) / |r - r_p|
 *     normalised so that the total charge integrates to zero (subtract
 *     the mean so the Ewald k=0 term vanishes).
 *  4. Compute E at r_p from this charge cloud using the *same* Ewald
 *     k-space formula used in the production run.  Because the cloud is
 *     symmetric in the continuum, any non-zero result is purely due to
 *     discretisation.
 *  5. Store E_self[ix][iy][iz][3] for later trilinear interpolation.
 *
 * Only the Fourier (k-space) part of the Ewald sum is used here.  The
 * real-space contribution from the discretised cloud would require a
 * real-space cutoff loop; because the cloud is smooth and converges
 * rapidly in k-space the Fourier part already captures the artefact.
 * ======================================================================== */

struct ewald_self_table_s {
  int    nx, ny, nz;          /* Sub-grid dimensions (e.g. 16 x 16 x 16) */
  double (*table)[3];        /* Flat array [nx*ny*nz][3]: E_self values   */
};

/*****************************************************************************
 *
 *  ewald_charge_sum_full_gpu_self_build
 *
 *  Build the self-field lookup table.
 *
 *  Parameters:
 *    ewald       - the Ewald object (must have been created already)
 *    kappa_debye - inverse Debye screening length (1/lambda_D, lattice units)
 *    nx, ny, nz  - lookup-table grid dimensions (typically 16)
 *    L           - side length of the scratch periodic box (lattice units).
 *                  Should be large enough that the Debye cloud is negligible
 *                  at the box boundaries; L ~ 8*lambda_D is usually safe.
 *    ptable      - output pointer; call ewald_charge_self_table_free() when done.
 *
 *****************************************************************************/

int ewald_charge_sum_full_gpu_self_build(ewald_charge_t* ewald,
                                          double kappa_debye,
                                          int    nx, int ny, int nz,
                                          int    L,
                                          ewald_self_table_t** ptable) {

  PI_DOUBLE(pi);

  assert(ewald);
  assert(ptable);
  assert(kappa_debye > 0.0);
  assert(nx > 0 && ny > 0 && nz > 0);
  assert(L > 0);

  /* -----------------------------------------------------------------------
   * Use EXACTLY the same Ewald parameters as ewald_charge_sum_full_gpu.
   * The scratch box is cubic with side L; nlocal = [L,L,L], noffset = [0,0,0].
   * The k-vectors are built from the scratch-box cell lengths ltot=[L,L,L].
   * ----------------------------------------------------------------------- */

  double ltot_self[3] = { (double)L, (double)L, (double)L };
  int    nlocal_self[3] = { L, L, L };
  int    noffset_self[3] = { 0, 0, 0 };
  int    nhalo_self = 0;

  double fkx = 2.0 * pi / ltot_self[X];
  double fky = 2.0 * pi / ltot_self[Y];
  double fkz = 2.0 * pi / ltot_self[Z];
  double r4alpha_sq = 1.0 / (4.0 * alpha_ * alpha_);
  double b0 = 1.0 / (ltot_self[X] * ltot_self[Y] * ltot_self[Z] * epsilon_);
  int    irc = (int)ceil(ewald_rc_);

  /* Number of k-vectors: same formula as ewald_charge_create */
  int nk_self = (int)ceil(alpha_ * alpha_ * ewald_rc_ * ltot_self[X] / pi);
  if (nk_self < 1) nk_self = 1;
  double kmax_self = pow(2.0 * pi * nk_self / ltot_self[X], 2.0);

  const int Ltot = L * L * L;
  /* nsites for SOA layout with nhalo=0: just L*L*L */
  const int nsites_self = Ltot;

  /* -----------------------------------------------------------------------
   * Build the cache filename from all parameters that define the table.
   * Format: efield_self_nx%d_ny%d_nz%d_L%d_kappa%.6g_alpha%.6g_rc%.6g_eps%.6g_beta%.6g_eunit%.6g.bin
   * ----------------------------------------------------------------------- */
  char self_fname[512];
  snprintf(self_fname, sizeof(self_fname),
           "efield_self_nx%d_ny%d_nz%d_L%d_kappa%.6g_alpha%.6g_rc%.6g_eps%.6g_beta%.6g_eunit%.6g.bin",
           nx, ny, nz, L,
           kappa_debye, alpha_, ewald_rc_, epsilon_, beta_, eunit_);

  pe_info(ewald->pe, "\n");
  pe_info(ewald->pe, "=== Ewald self-field correction table (%dx%dx%d) ===\n", nx, ny, nz);
  pe_info(ewald->pe, "  Cache file: %s\n", self_fname);
  pe_info(ewald->pe, "  Scratch box L = %d,  kappa_debye = %g\n", L, kappa_debye);
  pe_info(ewald->pe, "  alpha = %g,  rc = %g,  nk_self = %d\n", alpha_, ewald_rc_, nk_self);
  pe_info(ewald->pe, "  epsilon = %g,  beta = %g,  eunit = %g\n", epsilon_, beta_, eunit_);

  /* -----------------------------------------------------------------------
   * Allocate lookup table (host).
   * ----------------------------------------------------------------------- */

  ewald_self_table_t* table = (ewald_self_table_t*)malloc(sizeof(ewald_self_table_t));
  if (!table) pe_fatal(ewald->pe, "ewald_charge_sum_full_gpu_self_build: malloc table\n");

  table->nx = nx;
  table->ny = ny;
  table->nz = nz;
  const int ntab = nx * ny * nz;
  table->table = (double (*)[3]) calloc(ntab, 3 * sizeof(double));
  if (!table->table) {
    free(table);
    pe_fatal(ewald->pe, "ewald_charge_sum_full_gpu_self_build: calloc table->table\n");
  }

  /* -----------------------------------------------------------------------
   * Try to load from cache file.  If successful, skip computation entirely.
   * Binary format: ntab * 3 doubles in row-major order (ix fast outermost).
   * ----------------------------------------------------------------------- */
  {
    FILE* fc = fopen(self_fname, "rb");
    if (fc) {
      size_t nread = fread(table->table, 3 * sizeof(double), (size_t)ntab, fc);
      fclose(fc);
      if ((int)nread == ntab) {
        pe_info(ewald->pe, "  Loaded self-field table from cache: %s\n", self_fname);
        pe_info(ewald->pe, "=== Self-field correction table complete (from cache) ===\n\n");
        *ptable = table;
        return 0;
      }
      else {
        pe_info(ewald->pe, "  Cache file incomplete (%zu/%d entries), recomputing.\n",
                nread, ntab);
      }
    }
    else {
      pe_info(ewald->pe, "  Cache file not found, computing table.\n");
    }
  }

  /* -----------------------------------------------------------------------
   * Pre-compute k-vectors, G(k), and kz flags on host.
   * Identical loop to ewald_charge_sum_full_gpu (lines 3935-3954).
   * ----------------------------------------------------------------------- */

  int nk_upper = (2 * nk_self + 1) * (2 * nk_self + 1) * (nk_self + 1);
  double* kvec_h = (double*)malloc(3 * nk_upper * sizeof(double));
  double* Gk_h = (double*)malloc(nk_upper * sizeof(double));
  int* kzarr_h = (int*)malloc(nk_upper * sizeof(int));
  if (!kvec_h || !Gk_h || !kzarr_h)
    pe_fatal(ewald->pe, "ewald_charge_sum_full_gpu_self_build: malloc k-arrays\n");

  int nk_actual = 0;
  for (int kz = 0; kz <= nk_self; kz++) {
    for (int ky = -nk_self; ky <= nk_self; ky++) {
      for (int kx = -nk_self; kx <= nk_self; kx++) {
        double kxv = fkx * kx, kyv = fky * ky, kzv = fkz * kz;
        double ksq = kxv * kxv + kyv * kyv + kzv * kzv;
        if (ksq <= 0.0 || ksq > kmax_self) continue;
        kvec_h[3 * nk_actual + X] = kxv;
        kvec_h[3 * nk_actual + Y] = kyv;
        kvec_h[3 * nk_actual + Z] = kzv;
        Gk_h[nk_actual] = b0 * exp(-r4alpha_sq * ksq) / ksq;
        kzarr_h[nk_actual] = kz;
        nk_actual++;
      }
    }
  }
  pe_info(ewald->pe, "  nk_actual = %d\n", nk_actual);

#ifdef __NVCC__

  /* -----------------------------------------------------------------------
   * GPU path.
   *
   * The scratch-box rho is stored in SOA layout matching psi->rho:
   *   rho_self_d[nsites_self * 0 + idx] = rho (especie 0, la nube de Debye)
   *   rho_self_d[nsites_self * 1 + idx] = 0   (especie 1, vacía)
   *
   * With nhalo=0 the Ludwig index simplifies to:
   *   ludwig_idx = ic*(L*L) + jc*L + kc   (0-based: ic,jc,kc in 0..L-1)
   *
   * We then call EXACTLY the same kernels used in ewald_charge_sum_full_gpu:
   *   ewald_structure_factor_lattice_kernel
   *   ewald_particle_field_fourier_kernel
   *   ewald_particle_field_real_kernel
   *
   * The "particle" is the single probe at (xp,yp,zp) with charge q=1
   * (we want E/q, i.e. the field per unit charge; the actual particle charge
   *  cancels when we subtract Eself from Esub in production).
   * ----------------------------------------------------------------------- */

  const int threads = 256;
  const int blocks = (Ltot + threads - 1) / threads;

  /* Override device constants with scratch-box parameters.
   * This is safe because self_build is called once, before any simulation
   * step, and constants are reset at the start of each ewald_charge_sum_full_gpu call. */
  cudaMemcpyToSymbol(d_alpha, &alpha_, sizeof(double));
  cudaMemcpyToSymbol(d_eps_reg, &eps_reg_, sizeof(double));
  cudaMemcpyToSymbol(d_epsilon, &epsilon_, sizeof(double));
  cudaMemcpyToSymbol(d_beta, &beta_, sizeof(double));
  cudaMemcpyToSymbol(d_eunit, &eunit_, sizeof(double));
  cudaMemcpyToSymbol(d_rpi, &rpi_, sizeof(double));
  cudaMemcpyToSymbol(d_ewald_rc, &ewald_rc_, sizeof(double));
  cudaMemcpyToSymbol(d_fkx, &fkx, sizeof(double));
  cudaMemcpyToSymbol(d_fky, &fky, sizeof(double));
  cudaMemcpyToSymbol(d_fkz, &fkz, sizeof(double));
  cudaMemcpyToSymbol(d_r4alpha_sq, &r4alpha_sq, sizeof(double));
  cudaMemcpyToSymbol(d_b0, &b0, sizeof(double));
  cudaMemcpyToSymbol(d_kmax, &kmax_self, sizeof(double));
  cudaMemcpyToSymbol(d_nlocal, nlocal_self, 3 * sizeof(int));
  cudaMemcpyToSymbol(d_noffset, noffset_self, 3 * sizeof(int));

  /* Dipole correction = 0 for a charge-neutral cloud centred in the box */
  double E_dipole_zero[3] = { 0.0, 0.0, 0.0 };
  double M_dipole_zero[3] = { 0.0, 0.0, 0.0 };
  double dipole_pref_zero = 0.0;
  cudaMemcpyToSymbol(d_E_dipole, E_dipole_zero, 3 * sizeof(double));
  cudaMemcpyToSymbol(d_M_dipole, M_dipole_zero, 3 * sizeof(double));
  cudaMemcpyToSymbol(d_dipole_prefactor, &dipole_pref_zero, sizeof(double));

  /* Allocate scratch-box rho in SOA layout: 2 species * nsites_self */
  double* rho_self_d = NULL;
  cudaMalloc(&rho_self_d, 2 * nsites_self * sizeof(double));
  cudaMemset(rho_self_d, 0, 2 * nsites_self * sizeof(double));

  /* k-space arrays on device */
  double* kvec_d, * Gk_d, * Sk_sin_d, * Sk_cos_d;
  int* kzarr_d;
  cudaMalloc(&kvec_d, 3 * nk_actual * sizeof(double));
  cudaMalloc(&Gk_d, nk_actual * sizeof(double));
  cudaMalloc(&kzarr_d, nk_actual * sizeof(int));
  cudaMalloc(&Sk_sin_d, nk_actual * sizeof(double));
  cudaMalloc(&Sk_cos_d, nk_actual * sizeof(double));

  cudaMemcpy(kvec_d, kvec_h, 3 * nk_actual * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(Gk_d, Gk_h, nk_actual * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(kzarr_d, kzarr_h, nk_actual * sizeof(int), cudaMemcpyHostToDevice);

  /* Probe particle: single particle with charge q=1 */
  const int nparticles_self = 1;
  double* particle_r_d, * particle_q_d, * Esub_d, * fex_d;
  cudaMalloc(&particle_r_d, 3 * sizeof(double));
  cudaMalloc(&particle_q_d, sizeof(double));
  cudaMalloc(&Esub_d, 3 * sizeof(double));
  cudaMalloc(&fex_d, 3 * sizeof(double));
  double q_probe = 1.0;
  cudaMemcpy(particle_q_d, &q_probe, sizeof(double), cudaMemcpyHostToDevice);

  /* Host buffer for partial sums (mean subtraction of rho) */
  const int blocks_fill = (Ltot + threads - 1) / threads;
  double* partial_sum_h = (double*)malloc(blocks_fill * sizeof(double));
  double* partial_sum_d = NULL;
  if (!partial_sum_h) pe_fatal(ewald->pe, "ewald_charge_sum_full_gpu_self_build: malloc partial_sum_h\n");
  cudaMalloc(&partial_sum_d, blocks_fill * sizeof(double));

  const int p_blocks = 1;   /* single particle */
  int completed = 0;
  const int ntab_total = ntab;

  /* -----------------------------------------------------------------------
   * Main loop: one tab entry per (ix,iy,iz).
   * ----------------------------------------------------------------------- */

  for (int ix = 0; ix < nx; ix++) {
    double xf = ix / (double)nx;
    for (int iy = 0; iy < ny; iy++) {
      double yf = iy / (double)ny;
      for (int iz = 0; iz < nz; iz++) {
        double zf = iz / (double)nz;
        /* Particle sits in the central voxel of the scratch box (node L/2).
         * positions 1..L (same convention as ewald_structure_factor_lattice_kernel
         * with noffset=0, ic=1..L).  xf in [0,1) is the fractional offset
         * Particle placed in the central voxel: node L/2 + fractional offset xf.
         * Nodes are integers 1..L, so the central node is L/2 (integer for even L).
         * xf in [0,1) is the sub-voxel offset from that node. */
        int Lhalf = L / 2;
        double xp = (double)Lhalf + xf, yp = (double)Lhalf + yf, zp = (double)Lhalf + zf;

        /* -- Step 1: fill rho_self with Debye cloud (SOA species-0 slot) -- */
        /* K1a: fill + partial sums */
        self_fill_rho_kernel << <blocks_fill, threads, threads * sizeof(double) >> > (
            rho_self_d,          /* writes to species-0 slot directly */
            partial_sum_d,
            xp, yp, zp, kappa_debye, L);
        cudaDeviceSynchronize();

        /* Compute normalisation sum on host */
        cudaMemcpy(partial_sum_h, partial_sum_d,
                   blocks_fill * sizeof(double), cudaMemcpyDeviceToHost);
        double rho_sum = 0.0;
        for (int b = 0; b < blocks_fill; b++) rho_sum += partial_sum_h[b];
        /* rho_sum = Σ exp(-κr)/r.  K1b normalises: rho[i] = rho[i]/rho_sum - 1/Ltot
         * so that Σrho = 0 and the un-subtracted cloud integrates to 1 (unit charge),
         * matching the units of psi->rho in production (q_node sums to q_colloid). */

         /* K1b: normalise and subtract mean */
        self_subtract_mean_kernel << <blocks_fill, threads >> > (
            rho_self_d, rho_sum, Ltot);
        cudaDeviceSynchronize();
        /* Species-1 stays zero: q = rho[0] - rho[1] = rho[0] */

        /* -- Step 2: structure factors (same kernel as production) -- */
        cudaMemset(Sk_sin_d, 0, nk_actual * sizeof(double));
        cudaMemset(Sk_cos_d, 0, nk_actual * sizeof(double));

        ewald_structure_factor_lattice_kernel << <blocks, threads >> > (
            rho_self_d, Sk_sin_d, Sk_cos_d, kvec_d,
            nk_actual, nsites_self, nhalo_self);
        cudaDeviceSynchronize();

        /* -- Step 3: particle position on device -- */
        double r_probe_h[3] = { xp, yp, zp };
        cudaMemcpy(particle_r_d, r_probe_h, 3 * sizeof(double), cudaMemcpyHostToDevice);

        cudaMemset(Esub_d, 0, 3 * sizeof(double));
        cudaMemset(fex_d, 0, 3 * sizeof(double));

        /* -- Step 4: Fourier field at particle (same kernel as production) -- */
        ewald_particle_field_fourier_kernel << <p_blocks, threads >> > (
            Esub_d, fex_d,
            particle_r_d, particle_q_d,
            Sk_sin_d, Sk_cos_d, kvec_d, Gk_d, kzarr_d,
            nparticles_self, nk_actual);
        cudaDeviceSynchronize();

        /* -- Step 5: real-space field at particle (same kernel as production) -- */
        ewald_particle_field_real_kernel << <p_blocks, threads >> > (
            Esub_d, fex_d,
            particle_r_d, particle_q_d,
            rho_self_d,
            nparticles_self, nsites_self, nhalo_self, irc);
        cudaDeviceSynchronize();

        /* -- Step 6: read E_self from device -- */
        double Eself_h[3];
        cudaMemcpy(Eself_h, Esub_d, 3 * sizeof(double), cudaMemcpyDeviceToHost);

        int tab_idx = ix * ny * nz + iy * nz + iz;
        table->table[tab_idx][0] = Eself_h[0];
        table->table[tab_idx][1] = Eself_h[1];
        table->table[tab_idx][2] = Eself_h[2];

        /* Progress: print every ~5% */
        completed++;
        if (completed % (ntab_total / 20 + 1) == 0 || completed == ntab_total) {
          int pct = (int)(100.0 * completed / ntab_total);
          pe_info(ewald->pe,
                  "  Self-field table: %4d / %4d (%3d%%)  "
                  "E_self = (%10.3e, %10.3e, %10.3e)\n",
                  completed, ntab_total, pct,
                  Eself_h[0], Eself_h[1], Eself_h[2]);
        }

      } /* iz */
    } /* iy */
  } /* ix */

  /* Free device buffers */
  cudaFree(rho_self_d);
  cudaFree(kvec_d);
  cudaFree(Gk_d);
  cudaFree(kzarr_d);
  cudaFree(Sk_sin_d);
  cudaFree(Sk_cos_d);
  cudaFree(particle_r_d);
  cudaFree(particle_q_d);
  cudaFree(Esub_d);
  cudaFree(fex_d);
  cudaFree(partial_sum_d);
  free(partial_sum_h);

#else /* __NVCC__ not defined: CPU fallback */

  /* -----------------------------------------------------------------------
   * CPU path: mirrors the GPU path exactly using the same mathematical
   * formulas, without CUDA calls.
   * ----------------------------------------------------------------------- */

  double* rho_scratch = (double*)calloc(Ltot, sizeof(double));
  if (!rho_scratch)
    pe_fatal(ewald->pe, "ewald_charge_sum_full_gpu_self_build: calloc rho_scratch\n");

  int completed = 0;
  const int ntab_total = ntab;

  for (int ix = 0; ix < nx; ix++) {
    double xf = ix / (double)nx;
    for (int iy = 0; iy < ny; iy++) {
      double yf = iy / (double)ny;
      for (int iz = 0; iz < nz; iz++) {
        double zf = iz / (double)nz;
        /* Particle in central voxel: node L/2 + fractional offset (CPU path). */
        int Lhalf = L / 2;
        double xp = (double)Lhalf + xf, yp = (double)Lhalf + yf, zp = (double)Lhalf + zf;

        /* Step 1: Debye cloud + mean subtraction */
        double rho_sum = 0.0;
        for (int ax = 0; ax < L; ax++) {
          for (int ay = 0; ay < L; ay++) {
            for (int az = 0; az < L; az++) {
              /* Node at integer position (ax+1, ay+1, az+1) */
              double dx = (double)(ax + 1) - xp;
              double dy = (double)(ay + 1) - yp;
              double dz = (double)(az + 1) - zp;
              double Ld = (double)L;
              dx -= Ld * floor(dx / Ld + 0.5);
              dy -= Ld * floor(dy / Ld + 0.5);
              dz -= Ld * floor(dz / Ld + 0.5);
              double r = sqrt(dx * dx + dy * dy + dz * dz);
              double v = (r > 1.0e-10) ? exp(-kappa_debye * r) / r : 0.0;
              rho_scratch[ax * L * L + ay * L + az] = v;
              rho_sum += v;
            }
          }
        }
        /* Normalise so that Σρ_before_mean = 1 (unit charge, matching psi->rho units),
         * then subtract mean so Σρ_after = 0 (removes k=0 Ewald divergence). */
        for (int n = 0; n < Ltot; n++)
          rho_scratch[n] = rho_scratch[n] / rho_sum - 1.0 / (double)Ltot;

        /* Step 2: Fourier part (same formula as ewald_particle_field_fourier_kernel) */
        double Eself[3] = { 0.0, 0.0, 0.0 };
        for (int kn = 0; kn < nk_actual; kn++) {
          double kxv = kvec_h[3 * kn + X];
          double kyv = kvec_h[3 * kn + Y];
          double kzv = kvec_h[3 * kn + Z];
          double factor = (kzarr_h[kn] > 0) ? 2.0 : 1.0;

          /* Structure factors over all nodes (integer positions ax+1, same as GPU) */
          double S_sin = 0.0, S_cos = 0.0;
          for (int ax = 0; ax < L; ax++) {
            double rx = (double)(ax + 1) + noffset_self[X];
            for (int ay = 0; ay < L; ay++) {
              double ry = (double)(ay + 1) + noffset_self[Y];
              for (int az = 0; az < L; az++) {
                double rz = (double)(az + 1) + noffset_self[Z];
                double q = rho_scratch[ax * L * L + ay * L + az];
                if (fabs(q) < 1.0e-20) continue;
                double kr = kxv * rx + kyv * ry + kzv * rz;
                S_sin += q * sin(kr);
                S_cos += q * cos(kr);
              }
            }
          }

          /* Field at particle (ewald_particle_field_fourier_kernel formula) */
          double kr_p = kxv * xp + kyv * yp + kzv * zp;
          double im_part = S_sin * cos(kr_p) - S_cos * sin(kr_p);
          double coeff = factor * beta_ * eunit_ * Gk_h[kn] * im_part;
          Eself[0] += coeff * kxv;
          Eself[1] += coeff * kyv;
          Eself[2] += coeff * kzv;
        }
        /* Dipole correction = 0 (charge-neutral cloud) */

        /* Step 3: real-space part (same formula as ewald_particle_field_real_kernel) */
        int i0 = (int)floor(xp - noffset_self[X]);
        int j0 = (int)floor(yp - noffset_self[Y]);
        int k0 = (int)floor(zp - noffset_self[Z]);
        double Ld = (double)L;

        for (int di = -irc; di <= irc + 1; di++) {
          for (int dj = -irc; dj <= irc + 1; dj++) {
            for (int dk = -irc; dk <= irc + 1; dk++) {
              /* ni_phys = i0+di is the 1-based physical node index.
               * rho_scratch uses 0-based indexing: ax=0 → position ax+1=1.
               * So 0-based index = ni_phys - 1, wrapped periodically. */
              int ni = (((i0 + di) - 1) % L + L) % L;
              int nj = (((j0 + dj) - 1) % L + L) % L;
              int nk = (((k0 + dk) - 1) % L + L) % L;

              double q = rho_scratch[ni * L * L + nj * L + nk];
              if (fabs(q) < 1.0e-20) continue;

              /* Node at integer position (i0+di), consistent with ax+1 convention */
              double drx = xp - (double)(i0 + di);
              double dry = yp - (double)(j0 + dj);
              double drz = zp - (double)(k0 + dk);
              drx -= Ld * floor(drx / Ld + 0.5);
              dry -= Ld * floor(dry / Ld + 0.5);
              drz -= Ld * floor(drz / Ld + 0.5);

              double dist = sqrt(drx * drx + dry * dry + drz * drz);
              if (dist < ewald_rc_ && dist > 1.0e-10) {
                double ar = alpha_ * dist;
                double r_inv = 1.0 / dist;
                double term1 = erfc(ar) / (dist * dist);
                double term2 = (2.0 * alpha_ * rpi_) * exp(-ar * ar) * r_inv;
                double E_mag = beta_ * eunit_ * (term1 + term2) / (4.0 * pi * epsilon_);
                Eself[0] += q * E_mag * drx * r_inv;
                Eself[1] += q * E_mag * dry * r_inv;
                Eself[2] += q * E_mag * drz * r_inv;
              }
            }
          }
        }

        int tab_idx = ix * ny * nz + iy * nz + iz;
        table->table[tab_idx][0] = Eself[0];
        table->table[tab_idx][1] = Eself[1];
        table->table[tab_idx][2] = Eself[2];

        /* Progress: print every ~5% */
        completed++;
        if (completed % (ntab_total / 20 + 1) == 0 || completed == ntab_total) {
          int pct = (int)(100.0 * completed / ntab_total);
          pe_info(ewald->pe,
                  "  Self-field table: %4d / %4d (%3d%%)  "
                  "E_self = (%10.3e, %10.3e, %10.3e)\n",
                  completed, ntab_total, pct,
                  Eself[0], Eself[1], Eself[2]);
        }

      } /* iz */
    } /* iy */
  } /* ix */

  free(rho_scratch);

#endif /* __NVCC__ */

  free(kvec_h);
  free(Gk_h);
  free(kzarr_h);

  /* -----------------------------------------------------------------------
   * Save computed table: binary cache + human-readable text.
   * Text format: one line per entry "ix iy iz xf yf zf Ex Ey Ez"
   * ----------------------------------------------------------------------- */
  {
    FILE* fc = fopen(self_fname, "wb");
    if (fc) {
      size_t nwritten = fwrite(table->table, 3 * sizeof(double), (size_t)ntab, fc);
      fclose(fc);
      if ((int)nwritten == ntab) {
        pe_info(ewald->pe, "  Saved self-field table to cache: %s\n", self_fname);
      }
      else {
        pe_info(ewald->pe, "  Warning: only %zu/%d entries written to cache.\n",
                nwritten, ntab);
      }
    }
    else {
      pe_info(ewald->pe, "  Warning: could not open cache file for writing: %s\n",
              self_fname);
    }
  }

  /* Text file: same base name with .txt extension */
  {
    char txt_fname[512 + 4];
    snprintf(txt_fname, sizeof(txt_fname), "%s.txt", self_fname);
    FILE* ft = fopen(txt_fname, "w");
    if (ft) {
      fprintf(ft, "# Ewald self-field correction table\n");
      fprintf(ft, "# nx=%d ny=%d nz=%d  L=%d  kappa=%.6g\n", nx, ny, nz, L, kappa_debye);
      fprintf(ft, "# alpha=%.6g  rc=%.6g  epsilon=%.6g  beta=%.6g  eunit=%.6g\n",
              alpha_, ewald_rc_, epsilon_, beta_, eunit_);
      fprintf(ft, "# Columns: ix iy iz  xf yf zf  Ex Ey Ez\n");
      for (int ix = 0; ix < nx; ix++) {
        double xf = ix / (double)nx;
        for (int iy = 0; iy < ny; iy++) {
          double yf = iy / (double)ny;
          for (int iz = 0; iz < nz; iz++) {
            double zf = iz / (double)nz;
            int tab_idx = ix * ny * nz + iy * nz + iz;
            fprintf(ft, "%3d %3d %3d  %10.6f %10.6f %10.6f  %14.6e %14.6e %14.6e\n",
                    ix, iy, iz, xf, yf, zf,
                    table->table[tab_idx][0],
                    table->table[tab_idx][1],
                    table->table[tab_idx][2]);
          }
        }
      }
      fclose(ft);
      pe_info(ewald->pe, "  Saved self-field table (text): %s\n", txt_fname);
    }
    else {
      pe_info(ewald->pe, "  Warning: could not open text file for writing: %s\n", txt_fname);
    }
  }

  *ptable = table;

  pe_info(ewald->pe, "=== Self-field correction table complete ===\n\n");

  return 0;
}

/*****************************************************************************
 *
 *  ewald_charge_sum_full_gpu_self_build_symmetric
 *
 *  Build the self-field lookup table exploiting the Oh cubic symmetry of the
 *  discretised Debye cloud.  Only the irreducible simplex
 *
 *      ix >= iy >= iz,   ix in [0, n/2]
 *
 *  is computed (roughly 1/48 of all table entries for large n).  The rest are
 *  filled analytically from the symmetry relations:
 *
 *  Permutation of axes (P is any permutation of {x,y,z}):
 *      E_a( P(ix,iy,iz) ) = E_{P(a)}(ix,iy,iz)
 *
 *  Inversion of axis a (for the field component along that axis):
 *      E_a( ..., N-i_a, ... ) = -E_a( ..., i_a, ... )
 *      E_b( ..., N-i_a, ... ) =  E_b( ..., i_a, ... )   (b != a)
 *
 *  These two together generate the full Oh group restricted to [0,N)^3.
 *
 *  Requires n even (so that n/2 is an integer) and nx=ny=nz=n.
 *
 *****************************************************************************/

 /* Helper: compute one table entry (ix,iy,iz) using the same GPU/CPU path
  * as ewald_charge_sum_full_gpu_self_build.  Writes result into E[3].
  * All arrays (kvec_h, Gk_h, kzarr_h, nk_actual, L, etc.) come from the
  * caller's scope via the shared ewald object and pre-built k-arrays.
  *
  * We factor the per-point computation into a static helper so that
  * ewald_charge_sum_full_gpu_self_build_symmetric can share the GPU setup
  * code without duplicating it.
  */

  /* Forward-declare the helper type used in the symmetric builder */
typedef struct {
  double* kvec_h;
  double* Gk_h;
  int* kzarr_h;
  int      nk_actual;
  int      L;
  int      nsites_self;
  int      nhalo_self;
  int      irc;
  int      noffset_self[3];
#ifdef __NVCC__
  double* rho_self_d;
  double* kvec_d;
  double* Gk_d;
  int* kzarr_d;
  double* Sk_sin_d;
  double* Sk_cos_d;
  double* particle_r_d;
  double* particle_q_d;
  double* Esub_d;
  double* fex_d;
  double* partial_sum_d;
  double* partial_sum_h;
  int      blocks_fill;
  int      threads;
  int      blocks;
  int      p_blocks;
#else
  double* rho_scratch;
#endif
} ewald_self_ctx_t;

/* Compute E_self for fractional position (xf,yf,zf) within the unit cell.
 * Central node is at Lhalf = L/2; particle at (Lhalf+xf, Lhalf+yf, Lhalf+zf). */
static void ewald_self_compute_one(ewald_charge_t* ewald,
                                   ewald_self_ctx_t* ctx,
                                   double xf, double yf, double zf,
                                   double E[3]) {
  PI_DOUBLE(pi);
  int    L = ctx->L;
  int    Lhalf = L / 2;
  double xp = (double)Lhalf + xf;
  double yp = (double)Lhalf + yf;
  double zp = (double)Lhalf + zf;

#ifdef __NVCC__
  int threads = ctx->threads;
  int blocks_fill = ctx->blocks_fill;
  int blocks = ctx->blocks;
  int p_blocks = ctx->p_blocks;
  int Ltot = L * L * L;
  int nk_actual = ctx->nk_actual;

  self_fill_rho_kernel << <blocks_fill, threads, threads * sizeof(double) >> > (
      ctx->rho_self_d, ctx->partial_sum_d, xp, yp, zp, ewald->kappa_debye_stored, L);
  cudaDeviceSynchronize();

  cudaMemcpy(ctx->partial_sum_h, ctx->partial_sum_d,
             blocks_fill * sizeof(double), cudaMemcpyDeviceToHost);
  double rho_sum = 0.0;
  for (int b = 0; b < blocks_fill; b++) rho_sum += ctx->partial_sum_h[b];

  self_subtract_mean_kernel << <blocks_fill, threads >> > (
      ctx->rho_self_d, rho_sum, Ltot);
  cudaDeviceSynchronize();

  cudaMemset(ctx->Sk_sin_d, 0, nk_actual * sizeof(double));
  cudaMemset(ctx->Sk_cos_d, 0, nk_actual * sizeof(double));
  ewald_structure_factor_lattice_kernel << <blocks, threads >> > (
      ctx->rho_self_d, ctx->Sk_sin_d, ctx->Sk_cos_d, ctx->kvec_d,
      nk_actual, ctx->nsites_self, ctx->nhalo_self);
  cudaDeviceSynchronize();

  double r_probe_h[3] = { xp, yp, zp };
  cudaMemcpy(ctx->particle_r_d, r_probe_h, 3 * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemset(ctx->Esub_d, 0, 3 * sizeof(double));
  cudaMemset(ctx->fex_d, 0, 3 * sizeof(double));

  ewald_particle_field_fourier_kernel << <p_blocks, threads >> > (
      ctx->Esub_d, ctx->fex_d,
      ctx->particle_r_d, ctx->particle_q_d,
      ctx->Sk_sin_d, ctx->Sk_cos_d, ctx->kvec_d, ctx->Gk_d, ctx->kzarr_d,
      1, nk_actual);
  cudaDeviceSynchronize();

  ewald_particle_field_real_kernel << <p_blocks, threads >> > (
      ctx->Esub_d, ctx->fex_d,
      ctx->particle_r_d, ctx->particle_q_d,
      ctx->rho_self_d,
      1, ctx->nsites_self, ctx->nhalo_self, ctx->irc);
  cudaDeviceSynchronize();

  cudaMemcpy(E, ctx->Esub_d, 3 * sizeof(double), cudaMemcpyDeviceToHost);

#else /* CPU path */
  int    L3 = L * L * L;
  double Ld = (double)L;
  double* rho = ctx->rho_scratch;
  double rho_sum = 0.0;
  int    noffset_self[3] = { 0,0,0 };

  for (int ax = 0; ax < L; ax++) {
    for (int ay = 0; ay < L; ay++) {
      for (int az = 0; az < L; az++) {
        double dx = (double)(ax + 1) - xp;
        double dy = (double)(ay + 1) - yp;
        double dz = (double)(az + 1) - zp;
        dx -= Ld * floor(dx / Ld + 0.5);
        dy -= Ld * floor(dy / Ld + 0.5);
        dz -= Ld * floor(dz / Ld + 0.5);
        double r = sqrt(dx * dx + dy * dy + dz * dz);
        double v = (r > 1.0e-10) ? exp(-ewald->kappa_debye_stored * r) / r : 0.0;
        rho[ax * L * L + ay * L + az] = v;
        rho_sum += v;
      }
    }
  }
  for (int n = 0; n < L3; n++)
    rho[n] = rho[n] / rho_sum - 1.0 / (double)L3;

  E[0] = E[1] = E[2] = 0.0;
  double* kvec_h = ctx->kvec_h;
  double* Gk_h = ctx->Gk_h;
  int* kzarr_h = ctx->kzarr_h;
  int      nk_actual = ctx->nk_actual;

  for (int kn = 0; kn < nk_actual; kn++) {
    double kxv = kvec_h[3 * kn + X], kyv = kvec_h[3 * kn + Y], kzv = kvec_h[3 * kn + Z];
    double factor = (kzarr_h[kn] > 0) ? 2.0 : 1.0;
    double S_sin = 0.0, S_cos = 0.0;
    for (int ax = 0; ax < L; ax++) {
      double rx = (double)(ax + 1);
      for (int ay = 0; ay < L; ay++) {
        double ry = (double)(ay + 1);
        for (int az = 0; az < L; az++) {
          double rz = (double)(az + 1);
          double q = rho[ax * L * L + ay * L + az];
          if (fabs(q) < 1.0e-20) continue;
          double kr = kxv * rx + kyv * ry + kzv * rz;
          S_sin += q * sin(kr);
          S_cos += q * cos(kr);
        }
      }
    }
    double kr_p = kxv * xp + kyv * yp + kzv * zp;
    double im_part = S_sin * cos(kr_p) - S_cos * sin(kr_p);
    double coeff = factor * beta_ * eunit_ * Gk_h[kn] * im_part;
    E[0] += coeff * kxv;
    E[1] += coeff * kyv;
    E[2] += coeff * kzv;
  }

  int i0 = (int)floor(xp);
  int j0 = (int)floor(yp);
  int k0 = (int)floor(zp);
  int irc = ctx->irc;

  for (int di = -irc; di <= irc + 1; di++) {
    for (int dj = -irc; dj <= irc + 1; dj++) {
      for (int dk = -irc; dk <= irc + 1; dk++) {
        int ni = (((i0 + di) - 1) % L + L) % L;
        int nj = (((j0 + dj) - 1) % L + L) % L;
        int nk = (((k0 + dk) - 1) % L + L) % L;
        double q = rho[ni * L * L + nj * L + nk];
        if (fabs(q) < 1.0e-20) continue;
        double drx = xp - (double)(i0 + di);
        double dry = yp - (double)(j0 + dj);
        double drz = zp - (double)(k0 + dk);
        drx -= Ld * floor(drx / Ld + 0.5);
        dry -= Ld * floor(dry / Ld + 0.5);
        drz -= Ld * floor(drz / Ld + 0.5);
        double dist = sqrt(drx * drx + dry * dry + drz * drz);
        if (dist < ewald_rc_ && dist > 1.0e-10) {
          double ar = alpha_ * dist;
          double r_inv = 1.0 / dist;
          double term1 = erfc(ar) / (dist * dist);
          double term2 = (2.0 * alpha_ * rpi_) * exp(-ar * ar) * r_inv;
          double E_mag = beta_ * eunit_ * (term1 + term2) / (4.0 * pi * epsilon_);
          E[0] += q * E_mag * drx * r_inv;
          E[1] += q * E_mag * dry * r_inv;
          E[2] += q * E_mag * drz * r_inv;
        }
      }
    }
  }
#endif
}

int ewald_charge_sum_full_gpu_self_build_symmetric(ewald_charge_t* ewald,
                                                   double kappa_debye,
                                                   int    n,
                                                   int    L,
                                                   ewald_self_table_t** ptable) {
  PI_DOUBLE(pi);

  assert(ewald);
  assert(ptable);
  assert(kappa_debye > 0.0);
  assert(n > 0 && n % 2 == 0 && "n must be even");
  assert(L > 0 && L % 2 == 0 && "L must be even");

  /* -----------------------------------------------------------------------
   * Same Ewald parameters and k-vector setup as ewald_charge_sum_full_gpu_self_build.
   * ----------------------------------------------------------------------- */
  double ltot_self[3] = { (double)L, (double)L, (double)L };
  int    nlocal_self[3] = { L, L, L };
  int    noffset_self[3] = { 0, 0, 0 };
  int    nhalo_self = 0;
  int    nsites_self = L * L * L;

  double fkx = 2.0 * pi / ltot_self[X];
  double fky = 2.0 * pi / ltot_self[Y];
  double fkz = 2.0 * pi / ltot_self[Z];
  double r4alpha_sq = 1.0 / (4.0 * alpha_ * alpha_);
  double b0 = 1.0 / (ltot_self[X] * ltot_self[Y] * ltot_self[Z] * epsilon_);
  int    irc = (int)ceil(ewald_rc_);

  int nk_self = (int)ceil(alpha_ * alpha_ * ewald_rc_ * ltot_self[X] / pi);
  if (nk_self < 1) nk_self = 1;
  double kmax_self = pow(2.0 * pi * nk_self / ltot_self[X], 2.0);

  /* Build filename */
  char self_fname[512];
  snprintf(self_fname, sizeof(self_fname),
           "efield_self_sym_n%d_L%d_kappa%.6g_alpha%.6g_rc%.6g_eps%.6g_beta%.6g_eunit%.6g.bin",
           n, L, kappa_debye, alpha_, ewald_rc_, epsilon_, beta_, eunit_);

  pe_info(ewald->pe, "\n");
  pe_info(ewald->pe, "=== Ewald self-field correction table SYMMETRIC (%dx%dx%d) ===\n", n, n, n);
  pe_info(ewald->pe, "  Cache file: %s\n", self_fname);
  pe_info(ewald->pe, "  Scratch box L = %d,  kappa_debye = %g\n", L, kappa_debye);
  pe_info(ewald->pe, "  alpha = %g,  rc = %g\n", alpha_, ewald_rc_);

  /* Allocate table */
  const int ntab = n * n * n;
  ewald_self_table_t* table = (ewald_self_table_t*)malloc(sizeof(ewald_self_table_t));
  if (!table) pe_fatal(ewald->pe, "ewald_charge_sum_full_gpu_self_build_symmetric: malloc\n");
  table->nx = table->ny = table->nz = n;
  table->table = (double (*)[3]) calloc(ntab, 3 * sizeof(double));
  if (!table->table) { free(table); pe_fatal(ewald->pe, "calloc table->table\n"); }

  /* Try cache */
  {
    FILE* fc = fopen(self_fname, "rb");
    if (fc) {
      size_t nread = fread(table->table, 3 * sizeof(double), (size_t)ntab, fc);
      fclose(fc);
      if ((int)nread == ntab) {
        pe_info(ewald->pe, "  Loaded from cache: %s\n", self_fname);
        pe_info(ewald->pe, "=== Self-field (symmetric) table complete (from cache) ===\n\n");
        *ptable = table;
        return 0;
      }
      pe_info(ewald->pe, "  Cache incomplete, recomputing.\n");
    }
    else {
      pe_info(ewald->pe, "  Cache not found, computing.\n");
    }
  }

  /* k-vectors */
  int nk_upper = (2 * nk_self + 1) * (2 * nk_self + 1) * (nk_self + 1);
  double* kvec_h = (double*)malloc(3 * nk_upper * sizeof(double));
  double* Gk_h = (double*)malloc(nk_upper * sizeof(double));
  int* kzarr_h = (int*)malloc(nk_upper * sizeof(int));
  if (!kvec_h || !Gk_h || !kzarr_h)
    pe_fatal(ewald->pe, "malloc k-arrays\n");

  int nk_actual = 0;
  for (int kz = 0; kz <= nk_self; kz++) {
    for (int ky = -nk_self; ky <= nk_self; ky++) {
      for (int kx = -nk_self; kx <= nk_self; kx++) {
        double kxv = fkx * kx, kyv = fky * ky, kzv = fkz * kz;
        double ksq = kxv * kxv + kyv * kyv + kzv * kzv;
        if (ksq <= 0.0 || ksq > kmax_self) continue;
        kvec_h[3 * nk_actual + X] = kxv;
        kvec_h[3 * nk_actual + Y] = kyv;
        kvec_h[3 * nk_actual + Z] = kzv;
        Gk_h[nk_actual] = b0 * exp(-r4alpha_sq * ksq) / ksq;
        kzarr_h[nk_actual] = kz;
        nk_actual++;
      }
    }
  }
  pe_info(ewald->pe, "  nk_actual = %d\n", nk_actual);

  /* Build context */
  ewald_self_ctx_t ctx;
  memset(&ctx, 0, sizeof(ctx));
  ctx.kvec_h = kvec_h;
  ctx.Gk_h = Gk_h;
  ctx.kzarr_h = kzarr_h;
  ctx.nk_actual = nk_actual;
  ctx.L = L;
  ctx.nsites_self = nsites_self;
  ctx.nhalo_self = nhalo_self;
  ctx.irc = irc;
  ctx.noffset_self[0] = ctx.noffset_self[1] = ctx.noffset_self[2] = 0;

#ifdef __NVCC__
  const int threads = 256;
  const int blocks = (nsites_self + threads - 1) / threads;
  const int blocks_fill = blocks;
  ctx.threads = threads;
  ctx.blocks = blocks;
  ctx.blocks_fill = blocks_fill;
  ctx.p_blocks = 1;

  cudaMemcpyToSymbol(d_alpha, &alpha_, sizeof(double));
  cudaMemcpyToSymbol(d_eps_reg, &eps_reg_, sizeof(double));
  cudaMemcpyToSymbol(d_epsilon, &epsilon_, sizeof(double));
  cudaMemcpyToSymbol(d_beta, &beta_, sizeof(double));
  cudaMemcpyToSymbol(d_eunit, &eunit_, sizeof(double));
  cudaMemcpyToSymbol(d_rpi, &rpi_, sizeof(double));
  cudaMemcpyToSymbol(d_ewald_rc, &ewald_rc_, sizeof(double));
  cudaMemcpyToSymbol(d_fkx, &fkx, sizeof(double));
  cudaMemcpyToSymbol(d_fky, &fky, sizeof(double));
  cudaMemcpyToSymbol(d_fkz, &fkz, sizeof(double));
  cudaMemcpyToSymbol(d_r4alpha_sq, &r4alpha_sq, sizeof(double));
  cudaMemcpyToSymbol(d_b0, &b0, sizeof(double));
  cudaMemcpyToSymbol(d_kmax, &kmax_self, sizeof(double));
  cudaMemcpyToSymbol(d_nlocal, nlocal_self, 3 * sizeof(int));
  cudaMemcpyToSymbol(d_noffset, noffset_self, 3 * sizeof(int));
  double E_dip0[3] = { 0,0,0 }, M_dip0[3] = { 0,0,0 }, dpref0 = 0.0;
  cudaMemcpyToSymbol(d_E_dipole, E_dip0, 3 * sizeof(double));
  cudaMemcpyToSymbol(d_M_dipole, M_dip0, 3 * sizeof(double));
  cudaMemcpyToSymbol(d_dipole_prefactor, &dpref0, sizeof(double));

  cudaMalloc(&ctx.rho_self_d, 2 * nsites_self * sizeof(double));
  cudaMemset(ctx.rho_self_d, 0, 2 * nsites_self * sizeof(double));
  cudaMalloc(&ctx.kvec_d, 3 * nk_actual * sizeof(double));
  cudaMalloc(&ctx.Gk_d, nk_actual * sizeof(double));
  cudaMalloc(&ctx.kzarr_d, nk_actual * sizeof(int));
  cudaMalloc(&ctx.Sk_sin_d, nk_actual * sizeof(double));
  cudaMalloc(&ctx.Sk_cos_d, nk_actual * sizeof(double));
  cudaMemcpy(ctx.kvec_d, kvec_h, 3 * nk_actual * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(ctx.Gk_d, Gk_h, nk_actual * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(ctx.kzarr_d, kzarr_h, nk_actual * sizeof(int), cudaMemcpyHostToDevice);
  cudaMalloc(&ctx.particle_r_d, 3 * sizeof(double));
  cudaMalloc(&ctx.particle_q_d, sizeof(double));
  cudaMalloc(&ctx.Esub_d, 3 * sizeof(double));
  cudaMalloc(&ctx.fex_d, 3 * sizeof(double));
  double q_probe = 1.0;
  cudaMemcpy(ctx.particle_q_d, &q_probe, sizeof(double), cudaMemcpyHostToDevice);
  ctx.partial_sum_h = (double*)malloc(blocks_fill * sizeof(double));
  cudaMalloc(&ctx.partial_sum_d, blocks_fill * sizeof(double));

  /* Store kappa_debye in ewald object for use by helper (GPU reads it via kernel arg) */
  ewald->kappa_debye_stored = kappa_debye;

#else
  ctx.rho_scratch = (double*)calloc(nsites_self, sizeof(double));
  ewald->kappa_debye_stored = kappa_debye;
#endif

  /* -----------------------------------------------------------------------
   * Main loop: compute only the irreducible simplex ix >= iy >= iz in [0,n/2].
   * Then expand by symmetry to fill the full [0,n)^3 table.
   * ----------------------------------------------------------------------- */

  int Nhalf = n / 2;
  int completed = 0;
  /* Count points in simplex for progress reporting */
  int nsimplex = 0;
  for (int ix = 0; ix <= Nhalf; ix++)
    for (int iy = 0; iy <= ix; iy++)
      for (int iz = 0; iz <= iy; iz++)
        nsimplex++;

  pe_info(ewald->pe, "  Irreducible simplex: %d points (out of %d total)\n",
          nsimplex, ntab);

  /* Temporary storage for irreducible points before symmetry expansion */
  /* We store computed E directly into table at (ix,iy,iz) and expand afterwards */

  for (int ix = 0; ix <= Nhalf; ix++) {
    double xf = ix / (double)n;
    for (int iy = 0; iy <= ix; iy++) {
      double yf = iy / (double)n;
      for (int iz = 0; iz <= iy; iz++) {
        double zf = iz / (double)n;

        double E[3];
        ewald_self_compute_one(ewald, &ctx, xf, yf, zf, E);

        /* Store in irreducible slot */
        table->table[ix * n * n + iy * n + iz][0] = E[0];
        table->table[ix * n * n + iy * n + iz][1] = E[1];
        table->table[ix * n * n + iy * n + iz][2] = E[2];

        completed++;
        if (completed % (nsimplex / 20 + 1) == 0 || completed == nsimplex) {
          int pct = (int)(100.0 * completed / nsimplex);
          pe_info(ewald->pe,
                  "  Self-field (sym): %4d / %4d (%3d%%)  ix=%d iy=%d iz=%d"
                  "  E=(%10.3e,%10.3e,%10.3e)\n",
                  completed, nsimplex, pct, ix, iy, iz, E[0], E[1], E[2]);
        }
      }
    }
  }

  /* -----------------------------------------------------------------------
   * Symmetry expansion.
   *
   * For a point (ix,iy,iz) in the simplex, the 48 images are generated by
   * all permutations of (ix,iy,iz) and all sign flips of each coordinate
   * (with the inversion x -> n-x mod n, since xf = ix/n and the field is
   * antiperiodic in each component along its own axis).
   *
   * Given E[3] computed at (ix,iy,iz):
   *   - Permuting (ix,iy,iz) -> (ia,ib,ic) also permutes E: E' = (E[p0],E[p1],E[p2])
   *   - Flipping axis a: index a -> (n - i_a) % n, and E[a] -> -E[a], others unchanged
   * ----------------------------------------------------------------------- */

   /* All 6 permutations of (x,y,z) */
  static const int perms[6][3] = {
    {0,1,2}, {0,2,1}, {1,0,2}, {1,2,0}, {2,0,1}, {2,1,0}
  };

  for (int ix = 0; ix <= Nhalf; ix++) {
    for (int iy = 0; iy <= ix; iy++) {
      for (int iz = 0; iz <= iy; iz++) {
        /* Read the computed value */
        double Ex = table->table[ix * n * n + iy * n + iz][0];
        double Ey = table->table[ix * n * n + iy * n + iz][1];
        double Ez = table->table[ix * n * n + iy * n + iz][2];
        double Ebase[3] = { Ex, Ey, Ez };
        int    Ibase[3] = { ix, iy, iz };

        /* Loop over all 6 permutations */
        for (int p = 0; p < 6; p++) {
          int   Ip[3] = { Ibase[perms[p][0]], Ibase[perms[p][1]], Ibase[perms[p][2]] };
          double Ep[3] = { Ebase[perms[p][0]], Ebase[perms[p][1]], Ebase[perms[p][2]] };

          /* Loop over all 8 sign combinations (sx,sy,sz in {0,1}) */
          for (int sx = 0; sx < 2; sx++) {
            for (int sy = 0; sy < 2; sy++) {
              for (int sz = 0; sz < 2; sz++) {
                /* Flipped indices (periodic wrap: 0 stays 0, i -> n-i) */
                int jx = (sx && Ip[0] > 0) ? n - Ip[0] : Ip[0];
                int jy = (sy && Ip[1] > 0) ? n - Ip[1] : Ip[1];
                int jz = (sz && Ip[2] > 0) ? n - Ip[2] : Ip[2];

                /* Flipped field (flip sign of component along flipped axis) */
                double Fx = (sx) ? -Ep[0] : Ep[0];
                double Fy = (sy) ? -Ep[1] : Ep[1];
                double Fz = (sz) ? -Ep[2] : Ep[2];

                int tidx = jx * n * n + jy * n + jz;
                table->table[tidx][0] = Fx;
                table->table[tidx][1] = Fy;
                table->table[tidx][2] = Fz;
              }
            }
          }
        }
      }
    }
  }

  /* Cleanup context */
#ifdef __NVCC__
  cudaFree(ctx.rho_self_d);
  cudaFree(ctx.kvec_d);
  cudaFree(ctx.Gk_d);
  cudaFree(ctx.kzarr_d);
  cudaFree(ctx.Sk_sin_d);
  cudaFree(ctx.Sk_cos_d);
  cudaFree(ctx.particle_r_d);
  cudaFree(ctx.particle_q_d);
  cudaFree(ctx.Esub_d);
  cudaFree(ctx.fex_d);
  cudaFree(ctx.partial_sum_d);
  free(ctx.partial_sum_h);
#else
  free(ctx.rho_scratch);
#endif

  free(kvec_h);
  free(Gk_h);
  free(kzarr_h);

  /* Save binary cache */
  {
    FILE* fc = fopen(self_fname, "wb");
    if (fc) {
      size_t nw = fwrite(table->table, 3 * sizeof(double), (size_t)ntab, fc);
      fclose(fc);
      if ((int)nw == ntab)
        pe_info(ewald->pe, "  Saved to cache: %s\n", self_fname);
      else
        pe_info(ewald->pe, "  Warning: only %zu/%d written.\n", nw, ntab);
    }
    else {
      pe_info(ewald->pe, "  Warning: could not write cache: %s\n", self_fname);
    }
  }

  /* Save text file */
  {
    char txt_fname[520];
    snprintf(txt_fname, sizeof(txt_fname), "%s.txt", self_fname);
    FILE* ft = fopen(txt_fname, "w");
    if (ft) {
      fprintf(ft, "# Ewald self-field correction table (symmetric)\n");
      fprintf(ft, "# n=%d  L=%d  kappa=%.6g  alpha=%.6g  rc=%.6g\n",
              n, L, kappa_debye, alpha_, ewald_rc_);
      fprintf(ft, "# epsilon=%.6g  beta=%.6g  eunit=%.6g\n", epsilon_, beta_, eunit_);
      fprintf(ft, "# Columns: ix iy iz  xf yf zf  Ex Ey Ez\n");
      for (int ix = 0; ix < n; ix++) {
        double xf = ix / (double)n;
        for (int iy = 0; iy < n; iy++) {
          double yf = iy / (double)n;
          for (int iz = 0; iz < n; iz++) {
            double zf = iz / (double)n;
            int tidx = ix * n * n + iy * n + iz;
            fprintf(ft, "%3d %3d %3d  %10.6f %10.6f %10.6f  %14.6e %14.6e %14.6e\n",
                    ix, iy, iz, xf, yf, zf,
                    table->table[tidx][0],
                    table->table[tidx][1],
                    table->table[tidx][2]);
          }
        }
      }
      fclose(ft);
      pe_info(ewald->pe, "  Saved text: %s\n", txt_fname);
    }
  }

  *ptable = table;
  pe_info(ewald->pe, "=== Self-field (symmetric) table complete ===\n\n");
  return 0;
}

/*****************************************************************************
 *
 *  ewald_charge_self_table_free
 *
 *  Free the self-field lookup table.
 *
 *****************************************************************************/

int ewald_charge_self_table_free(ewald_self_table_t** ptable) {

  assert(ptable);

  if (*ptable == NULL) return 0;

  free((*ptable)->table);
  free(*ptable);
  *ptable = NULL;

  return 0;
}

/*****************************************************************************
 *
 *  ewald_charge_self_field_interpolate
 *
 *  Return the self-field correction at fractional position (xf, yf, zf)
 *  using trilinear interpolation of the lookup table.
 *
 *  xf, yf, zf must be in [0, 1).  The result Eself[3] is the spurious
 *  electric field that should be subtracted from the Ewald field:
 *
 *    F_corrected = q * (E_ewald - Eself)
 *
 *  (Same units as pc->Esub, i.e. beta * eunit * E_physical.)
 *
 *****************************************************************************/

int ewald_charge_self_field_interpolate(const ewald_self_table_t* table,
                                        double xf, double yf, double zf,
                                        double Eself[3]) {

  assert(table);
  assert(Eself);

  int nx = table->nx;
  int ny = table->ny;
  int nz = table->nz;

  /* Map fractional coordinates to table indices and weights.
   * Table entries sit at xf_table = ix/nx (ix = 0..nx-1), so the
   * continuous index is simply xf * nx.  We wrap periodically. */

  double tx = xf * (double)nx;   /* continuous index, 0-based */
  double ty = yf * (double)ny;
  double tz = zf * (double)nz;

  /* Wrap to [0, N) */
  tx = tx - (double)nx * floor(tx / (double)nx);
  ty = ty - (double)ny * floor(ty / (double)ny);
  tz = tz - (double)nz * floor(tz / (double)nz);

  int i0 = (int)floor(tx);
  int j0 = (int)floor(ty);
  int k0 = (int)floor(tz);

  double wx = tx - (double)i0;   /* fractional weight in [0,1) */
  double wy = ty - (double)j0;
  double wz = tz - (double)k0;

  /* Periodic neighbours */
  int i1 = (i0 + 1) % nx;
  int j1 = (j0 + 1) % ny;
  int k1 = (k0 + 1) % nz;

  /* Eight corner indices */
#define TAB(a,b,c) table->table[(a)*ny*nz + (b)*nz + (c)]

  for (int d = 0; d < 3; d++) {
    Eself[d] =
      (1.0 - wx) * (1.0 - wy) * (1.0 - wz) * TAB(i0, j0, k0)[d] +
      wx * (1.0 - wy) * (1.0 - wz) * TAB(i1, j0, k0)[d] +
      (1.0 - wx) * wy * (1.0 - wz) * TAB(i0, j1, k0)[d] +
      wx * wy * (1.0 - wz) * TAB(i1, j1, k0)[d] +
      (1.0 - wx) * (1.0 - wy) * wz * TAB(i0, j0, k1)[d] +
      wx * (1.0 - wy) * wz * TAB(i1, j0, k1)[d] +
      (1.0 - wx) * wy * wz * TAB(i0, j1, k1)[d] +
      wx * wy * wz * TAB(i1, j1, k1)[d];
  }

#undef TAB

  return 0;
}
/*CHANGE END - 20260309 Self-field correction for lattice-discretisation artefact */

/* CHANGE INIT - Gaussian_Ewald - Gaussian Ewald sum implementation
 *
 * All kernels and the main function are copies of the standard Ewald versions
 * with suffix _gaussian.  The physics changes are:
 *
 *   Real-space:  erfc(alpha*r)/r  →  [erf(r/sigma) - erf(eta*r)] / r
 *                The potential is finite at r=0; no eps_reg regularization needed.
 *
 *   Structure factor:  q * exp(i k·r)  →  q * fk * exp(i k·r)
 *                      where fk = exp(-k² sigma² / 4)  (Gaussian form factor)
 *                      Applied to BOTH lattice nodes and particles.
 *
 *   Green function:  G(k) = (4π/εk²) exp(-k²/4α²)  →  (4π/εk²) exp(-k²/4η²)
 *                   Only alpha→eta substitution; sigma factor is in rho(k), NOT here,
 *                   to avoid double-counting.
 *
 *   sigma: physical Gaussian width (controls short-range smoothing)
 *   eta:   Ewald splitting parameter (algorithmic, controls real/Fourier balance)
 *          Analogous to alpha in standard Ewald.
 *
 * Reference: GAUSSIAN_EWALD_MODIFICATIONS.md in claude/
 */

#ifdef __NVCC__

 /* -------------------------------------------------------------------------
  * ewald_structure_factor_lattice_kernel_gaussian
  *
  * Copy of ewald_structure_factor_lattice_kernel with Gaussian form factor.
  * Each charge q at a node contributes  q * fk  to the structure factor,
  * where fk = exp(-k² sigma² / 4) is precomputed in gaussian_fk_d[kn].
  * ------------------------------------------------------------------------- */

__global__ void ewald_structure_factor_lattice_kernel_gaussian(
    const double* __restrict__ rho_data,
    double* __restrict__ Sk_sin,
    double* __restrict__ Sk_cos,
    const double* __restrict__ kvec,
    const double* __restrict__ gaussian_fk,   /* exp(-k²σ²/4) per k-vector */
    int nktot,
    int nsites,
    int nhalo) {

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int ntotal = d_nlocal[0] * d_nlocal[1] * d_nlocal[2];

  if (idx >= ntotal) return;

  int kc = idx % d_nlocal[2] + 1;
  int jc = (idx / d_nlocal[2]) % d_nlocal[1] + 1;
  int ic = idx / (d_nlocal[2] * d_nlocal[1]) + 1;

  int str_z = 1;
  int str_y = (d_nlocal[2] + 2 * nhalo);
  int str_x = str_y * (d_nlocal[1] + 2 * nhalo);
  int ludwig_idx = str_x * (nhalo + ic - 1) + str_y * (nhalo + jc - 1) + str_z * (nhalo + kc - 1);

  double rho0 = rho_data[nsites * 0 + ludwig_idx];
  double rho1 = rho_data[nsites * 1 + ludwig_idx];
  double q = rho0 - rho1;

  if (fabs(q) < 1.0e-14) return;

  double rx = (double)(d_noffset[0] + ic);
  double ry = (double)(d_noffset[1] + jc);
  double rz = (double)(d_noffset[2] + kc);

  for (int kn = 0; kn < nktot; kn++) {
    double kx = kvec[3 * kn + 0];
    double ky = kvec[3 * kn + 1];
    double kz = kvec[3 * kn + 2];

    double kr = kx * rx + ky * ry + kz * rz;
    double sinkr = sin(kr);
    double coskr = cos(kr);

    /* CHANGE: multiply charge by Gaussian form factor fk = exp(-k²σ²/4) */
    double q_eff = q * gaussian_fk[kn];
    atomicAdd(&Sk_sin[kn], q_eff * sinkr);
    atomicAdd(&Sk_cos[kn], q_eff * coskr);
  }
}

/* -------------------------------------------------------------------------
 * ewald_structure_factor_particle_kernel_gaussian
 *
 * Copy of ewald_structure_factor_particle_kernel with Gaussian form factor.
 * ------------------------------------------------------------------------- */

__global__ void ewald_structure_factor_particle_kernel_gaussian(
    const double* __restrict__ particle_r,
    const double* __restrict__ particle_q,
    double* __restrict__ Sk_sin,
    double* __restrict__ Sk_cos,
    const double* __restrict__ kvec,
    const double* __restrict__ gaussian_fk,   /* exp(-k²σ²/4) per k-vector */
    int nparticles,
    int nktot) {

  int p = blockIdx.x * blockDim.x + threadIdx.x;

  if (p >= nparticles) return;

  double q = particle_q[p];
  if (fabs(q) < 1.0e-14) return;

  double rx = particle_r[3 * p + 0];
  double ry = particle_r[3 * p + 1];
  double rz = particle_r[3 * p + 2];

  for (int kn = 0; kn < nktot; kn++) {
    double kx = kvec[3 * kn + 0];
    double ky = kvec[3 * kn + 1];
    double kz = kvec[3 * kn + 2];

    double kr = kx * rx + ky * ry + kz * rz;
    double sinkr = sin(kr);
    double coskr = cos(kr);

    /* CHANGE: multiply charge by Gaussian form factor fk = exp(-k²σ²/4) */
    double q_eff = q * gaussian_fk[kn];
    atomicAdd(&Sk_sin[kn], q_eff * sinkr);
    atomicAdd(&Sk_cos[kn], q_eff * coskr);
  }
}

/* -------------------------------------------------------------------------
 * ewald_potential_real_particle_kernel_gaussian
 *
 * Copy of ewald_potential_real_particle_kernel with Gaussian real-space potential.
 * Uses [erf(r/sigma) - erf(eta*r)] / r  instead of  erfc(alpha*r) / r.
 * Potential is finite at r=0; no eps_reg needed.
 * d_gauss_sigma and d_gauss_eta are device constants set in the main function.
 * ------------------------------------------------------------------------- */

__device__ __constant__ double d_gauss_sigma;   /* Gaussian width σ */
__device__ __constant__ double d_gauss_eta;     /* Ewald splitting η */

/* CHANGE INIT - Gaussian_Ewald_Dual - per-species widths and precomputed σ_eff for pair interactions */
__device__ __constant__ double d_gauss_sigma_p;     /* Particle Gaussian width σ_p */
__device__ __constant__ double d_gauss_sigma_f;     /* Fluid node Gaussian width σ_f */
__device__ __constant__ double d_gauss_sigma_eff_pp; /* sqrt(2)·σ_p — particle-particle pair */
__device__ __constant__ double d_gauss_sigma_eff_ff; /* sqrt(2)·σ_f — fluid-fluid pair */
__device__ __constant__ double d_gauss_sigma_eff_pf; /* sqrt(σ_p²+σ_f²) — particle-fluid pair */
/* CHANGE END - Gaussian_Ewald_Dual */

__global__ void ewald_potential_real_particle_kernel_gaussian(
    double* __restrict__ psi_data,
    const double* __restrict__ particle_r,
    const double* __restrict__ particle_q,
    int nparticles,
    int nsites,
    int nhalo,
    int irc) {

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int ntotal = d_nlocal[0] * d_nlocal[1] * d_nlocal[2];

  if (idx >= ntotal) return;

  int kc = idx % d_nlocal[2] + 1;
  int jc = (idx / d_nlocal[2]) % d_nlocal[1] + 1;
  int ic = idx / (d_nlocal[2] * d_nlocal[1]) + 1;

  int str_z = 1;
  int str_y = (d_nlocal[2] + 2 * nhalo);
  int str_x = str_y * (d_nlocal[1] + 2 * nhalo);
  int ludwig_idx = str_x * (nhalo + ic - 1) + str_y * (nhalo + jc - 1) + str_z * (nhalo + kc - 1);

  double r_local[3] = { (double)ic, (double)jc, (double)kc };
  double phi_real = 0.0;

  for (int p = 0; p < nparticles; p++) {
    double q_p = particle_q[p];
    if (fabs(q_p) < 1.0e-14) continue;

    double r_p[3];
    r_p[0] = particle_r[3 * p + 0] - (double)d_noffset[0];
    r_p[1] = particle_r[3 * p + 1] - (double)d_noffset[1];
    r_p[2] = particle_r[3 * p + 2] - (double)d_noffset[2];

    double dr[3] = { r_local[0] - r_p[0], r_local[1] - r_p[1], r_local[2] - r_p[2] };
    double dist = sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]);

    if (dist < d_ewald_rc) {
      double phi_diff;
      if (dist < 1.0e-6) {
        /* Taylor: [erf(r/σ) - erf(ηr)]/r → (2/√π)(1/σ - η) at r=0 */
        phi_diff = (2.0 / sqrt(M_PI)) * (1.0 / d_gauss_sigma - d_gauss_eta);
      }
      else {
        /* CHANGE: Gaussian real-space potential */
        phi_diff = erf(dist / d_gauss_sigma) / dist - erf(d_gauss_eta * dist) / dist;
      }
      phi_real += q_p * phi_diff / (4.0 * M_PI * d_epsilon);
    }
  }

  psi_data[ludwig_idx] += d_beta * d_eunit * phi_real;
}

/* -------------------------------------------------------------------------
 * ewald_potential_real_lattice_kernel_gaussian
 *
 * Copy of ewald_potential_real_lattice_kernel with Gaussian real-space potential.
 * ------------------------------------------------------------------------- */

__global__ void ewald_potential_real_lattice_kernel_gaussian(
    double* __restrict__ psi_data,
    const double* __restrict__ rho_data,
    int nsites,
    int nhalo,
    int irc) {

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int ntotal = d_nlocal[0] * d_nlocal[1] * d_nlocal[2];

  if (idx >= ntotal) return;

  int kc = idx % d_nlocal[2] + 1;
  int jc = (idx / d_nlocal[2]) % d_nlocal[1] + 1;
  int ic = idx / (d_nlocal[2] * d_nlocal[1]) + 1;

  int str_z = 1;
  int str_y = (d_nlocal[2] + 2 * nhalo);
  int str_x = str_y * (d_nlocal[1] + 2 * nhalo);
  int ludwig_idx = str_x * (nhalo + ic - 1) + str_y * (nhalo + jc - 1) + str_z * (nhalo + kc - 1);

  double phi_real = 0.0;

  for (int di = -irc; di <= irc; di++) {
    int i2 = ic + di;
    if (i2 < 1 - nhalo || i2 > d_nlocal[0] + nhalo) continue;

    for (int dj = -irc; dj <= irc; dj++) {
      int j2 = jc + dj;
      if (j2 < 1 - nhalo || j2 > d_nlocal[1] + nhalo) continue;

      for (int dk = -irc; dk <= irc; dk++) {
        int k2 = kc + dk;
        if (k2 < 1 - nhalo || k2 > d_nlocal[2] + nhalo) continue;

        if (di == 0 && dj == 0 && dk == 0) continue;

        int ludwig_idx2 = str_x * (nhalo + i2 - 1) + str_y * (nhalo + j2 - 1) + str_z * (nhalo + k2 - 1);

        double rho0_2 = rho_data[nsites * 0 + ludwig_idx2];
        double rho1_2 = rho_data[nsites * 1 + ludwig_idx2];
        double q_node = rho0_2 - rho1_2;

        if (fabs(q_node) < 1.0e-14) continue;

        double dist = sqrt((double)(di * di + dj * dj + dk * dk));

        if (dist < d_ewald_rc) {
          double phi_diff;
          if (dist < 1.0e-6) {
            phi_diff = (2.0 / sqrt(M_PI)) * (1.0 / d_gauss_sigma - d_gauss_eta);
          }
          else {
            /* CHANGE: Gaussian real-space potential */
            phi_diff = erf(dist / d_gauss_sigma) / dist - erf(d_gauss_eta * dist) / dist;
          }
          phi_real += q_node * phi_diff / (4.0 * M_PI * d_epsilon);
        }
      }
    }
  }

  atomicAdd(&psi_data[ludwig_idx], d_beta * d_eunit * phi_real);
}

/* -------------------------------------------------------------------------
 * ewald_force_real_particle_kernel_gaussian
 *
 * Copy of ewald_force_real_particle_kernel with Gaussian real-space force.
 * d/dr [erf(r/σ)/r - erf(ηr)/r] computed directly to avoid cancellation.
 * ------------------------------------------------------------------------- */

__global__ void ewald_force_real_particle_kernel_gaussian(
    double* __restrict__ force_data,
    const double* __restrict__ rho_data,
    const double* __restrict__ particle_r,
    const double* __restrict__ particle_q,
    int nparticles,
    int nsites,
    int nhalo,
    int irc) {

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int ntotal = d_nlocal[0] * d_nlocal[1] * d_nlocal[2];

  if (idx >= ntotal) return;

  int kc = idx % d_nlocal[2] + 1;
  int jc = (idx / d_nlocal[2]) % d_nlocal[1] + 1;
  int ic = idx / (d_nlocal[2] * d_nlocal[1]) + 1;

  int str_z = 1;
  int str_y = (d_nlocal[2] + 2 * nhalo);
  int str_x = str_y * (d_nlocal[1] + 2 * nhalo);
  int ludwig_idx = str_x * (nhalo + ic - 1) + str_y * (nhalo + jc - 1) + str_z * (nhalo + kc - 1);

  double rho0 = rho_data[nsites * 0 + ludwig_idx];
  double rho1 = rho_data[nsites * 1 + ludwig_idx];
  double q1 = rho0 - rho1;

  double r_local[3] = { (double)ic, (double)jc, (double)kc };
  double F_real[3] = { 0.0, 0.0, 0.0 };

  for (int p = 0; p < nparticles; p++) {
    double q_p = particle_q[p];
    if (fabs(q_p) < 1.0e-14) continue;

    double r_p[3];
    r_p[0] = particle_r[3 * p + 0] - (double)d_noffset[0];
    r_p[1] = particle_r[3 * p + 1] - (double)d_noffset[1];
    r_p[2] = particle_r[3 * p + 2] - (double)d_noffset[2];

    double dr[3] = { r_local[0] - r_p[0], r_local[1] - r_p[1], r_local[2] - r_p[2] };
    double dist = sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]);

    if (dist < d_ewald_rc) {
      double dphi_dr;
      if (dist < 1.0e-5) {
        /* Taylor: d/dr[erf(r/σ)/r - erf(ηr)/r] ≈ -(4/3√π)(1/σ³ - η³)·r */
        double inv_s3 = 1.0 / (d_gauss_sigma * d_gauss_sigma * d_gauss_sigma);
        double eta3 = d_gauss_eta * d_gauss_eta * d_gauss_eta;
        dphi_dr = -(4.0 / (3.0 * sqrt(M_PI))) * (inv_s3 - eta3) * dist;
      }
      else {
        double r2 = dist * dist;
        double r_inv = 1.0 / dist;
        double u_sig = dist / d_gauss_sigma;
        double u_eta = d_gauss_eta * dist;
        /* CHANGE: d/dr[erf(r/σ)/r] - d/dr[erf(ηr)/r] */
        double dg_sig = (2.0 / (sqrt(M_PI) * d_gauss_sigma)) * exp(-u_sig * u_sig) * r_inv
          - erf(u_sig) / r2;
        double dg_eta = (2.0 * d_gauss_eta / sqrt(M_PI)) * exp(-u_eta * u_eta) * r_inv
          - erf(u_eta) / r2;
        dphi_dr = dg_sig - dg_eta;
      }
      double r_inv = (dist > 1.0e-10) ? 1.0 / dist : 0.0;
      double F_mag = q1 * q_p * (-dphi_dr) / (4.0 * M_PI * d_epsilon);

      F_real[0] += F_mag * dr[0] * r_inv;
      F_real[1] += F_mag * dr[1] * r_inv;
      F_real[2] += F_mag * dr[2] * r_inv;
    }
  }

  atomicAdd(&force_data[3 * ludwig_idx + 0], F_real[0]);
  atomicAdd(&force_data[3 * ludwig_idx + 1], F_real[1]);
  atomicAdd(&force_data[3 * ludwig_idx + 2], F_real[2]);
}

/* -------------------------------------------------------------------------
 * ewald_efield_real_particle_kernel_gaussian
 *
 * Copy of ewald_efield_real_particle_kernel with Gaussian real-space field.
 * Same as force kernel but without multiplying by q1 (field, not force).
 * ------------------------------------------------------------------------- */

__global__ void ewald_efield_real_particle_kernel_gaussian(
    double* __restrict__ efield_data,
    const double* __restrict__ particle_r,
    const double* __restrict__ particle_q,
    int nparticles,
    int nsites,
    int nhalo,
    int irc) {

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int ntotal = d_nlocal[0] * d_nlocal[1] * d_nlocal[2];

  if (idx >= ntotal) return;

  int kc = idx % d_nlocal[2] + 1;
  int jc = (idx / d_nlocal[2]) % d_nlocal[1] + 1;
  int ic = idx / (d_nlocal[2] * d_nlocal[1]) + 1;

  int str_z = 1;
  int str_y = (d_nlocal[2] + 2 * nhalo);
  int str_x = str_y * (d_nlocal[1] + 2 * nhalo);
  int ludwig_idx = str_x * (nhalo + ic - 1) + str_y * (nhalo + jc - 1) + str_z * (nhalo + kc - 1);

  double r_local[3] = { (double)ic, (double)jc, (double)kc };
  double E_real[3] = { 0.0, 0.0, 0.0 };

  for (int p = 0; p < nparticles; p++) {
    double q_p = particle_q[p];
    if (fabs(q_p) < 1.0e-14) continue;

    double r_p[3];
    r_p[0] = particle_r[3 * p + 0] - (double)d_noffset[0];
    r_p[1] = particle_r[3 * p + 1] - (double)d_noffset[1];
    r_p[2] = particle_r[3 * p + 2] - (double)d_noffset[2];

    double dr[3] = { r_local[0] - r_p[0], r_local[1] - r_p[1], r_local[2] - r_p[2] };
    double dist = sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]);

    if (dist < d_ewald_rc) {
      double dphi_dr;
      if (dist < 1.0e-5) {
        double inv_s3 = 1.0 / (d_gauss_sigma * d_gauss_sigma * d_gauss_sigma);
        double eta3 = d_gauss_eta * d_gauss_eta * d_gauss_eta;
        dphi_dr = -(4.0 / (3.0 * sqrt(M_PI))) * (inv_s3 - eta3) * dist;
      }
      else {
        double r2 = dist * dist;
        double r_inv = 1.0 / dist;
        double u_sig = dist / d_gauss_sigma;
        double u_eta = d_gauss_eta * dist;
        double dg_sig = (2.0 / (sqrt(M_PI) * d_gauss_sigma)) * exp(-u_sig * u_sig) * r_inv
          - erf(u_sig) / r2;
        double dg_eta = (2.0 * d_gauss_eta / sqrt(M_PI)) * exp(-u_eta * u_eta) * r_inv
          - erf(u_eta) / r2;
        dphi_dr = dg_sig - dg_eta;
      }
      double r_inv = (dist > 1.0e-10) ? 1.0 / dist : 0.0;
      /* CHANGE: E_mag = -q_p * dphi_dr / (4πε) — field without q_node */
      double E_mag = q_p * (-dphi_dr) / (4.0 * M_PI * d_epsilon);

      E_real[0] += E_mag * dr[0] * r_inv;
      E_real[1] += E_mag * dr[1] * r_inv;
      E_real[2] += E_mag * dr[2] * r_inv;
    }
  }

  atomicAdd(&efield_data[3 * ludwig_idx + 0], E_real[0]);
  atomicAdd(&efield_data[3 * ludwig_idx + 1], E_real[1]);
  atomicAdd(&efield_data[3 * ludwig_idx + 2], E_real[2]);
}

/* -------------------------------------------------------------------------
 * ewald_efield_real_lattice_kernel_gaussian
 *
 * Gaussian version of ewald_efield_real_lattice_kernel.
 * Uses [erf(r/σ) - erf(ηr)]/r instead of erfc(αr)/r for real-space field.
 * d/dr of that potential gives the field magnitude.
 * ------------------------------------------------------------------------- */

 /* CHANGE INIT - EfieldLatticGaussian - Gaussian real-space efield nodo->nodo */

__global__ void ewald_efield_real_lattice_kernel_gaussian(
    double* __restrict__ efield_data,
    const double* __restrict__ rho_data,
    int nsites,
    int nhalo,
    int irc) {

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int ntotal = d_nlocal[0] * d_nlocal[1] * d_nlocal[2];

  if (idx >= ntotal) return;

  int kc = idx % d_nlocal[2] + 1;
  int jc = (idx / d_nlocal[2]) % d_nlocal[1] + 1;
  int ic = idx / (d_nlocal[2] * d_nlocal[1]) + 1;

  int str_z = 1;
  int str_y = (d_nlocal[2] + 2 * nhalo);
  int str_x = str_y * (d_nlocal[1] + 2 * nhalo);
  int ludwig_idx = str_x * (nhalo + ic - 1) + str_y * (nhalo + jc - 1) + str_z * (nhalo + kc - 1);

  double E_real[3] = { 0.0, 0.0, 0.0 };

  for (int di = -irc; di <= irc; di++) {
    int i2 = ic + di;
    if (i2 < 1 - nhalo || i2 > d_nlocal[0] + nhalo) continue;

    for (int dj = -irc; dj <= irc; dj++) {
      int j2 = jc + dj;
      if (j2 < 1 - nhalo || j2 > d_nlocal[1] + nhalo) continue;

      for (int dk = -irc; dk <= irc; dk++) {
        int k2 = kc + dk;
        if (k2 < 1 - nhalo || k2 > d_nlocal[2] + nhalo) continue;

        if (di == 0 && dj == 0 && dk == 0) continue;

        int ludwig_idx2 = str_x * (nhalo + i2 - 1) + str_y * (nhalo + j2 - 1) + str_z * (nhalo + k2 - 1);

        double rho0_2 = rho_data[nsites * 0 + ludwig_idx2];
        double rho1_2 = rho_data[nsites * 1 + ludwig_idx2];
        double q_node = rho0_2 - rho1_2;

        if (fabs(q_node) < 1.0e-14) continue;

        /* Vector from source node to receptor node */
        double dr[3] = { -(double)di, -(double)dj, -(double)dk };
        double dist = sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]);

        if (dist < d_ewald_rc) {
          double dphi_dr;
          if (dist < 1.0e-5) {
            double inv_s3 = 1.0 / (d_gauss_sigma * d_gauss_sigma * d_gauss_sigma);
            double eta3 = d_gauss_eta * d_gauss_eta * d_gauss_eta;
            dphi_dr = -(4.0 / (3.0 * sqrt(M_PI))) * (inv_s3 - eta3) * dist;
          }
          else {
            double r2 = dist * dist;
            double r_inv = 1.0 / dist;
            double u_sig = dist / d_gauss_sigma;
            double u_eta = d_gauss_eta * dist;
            double dg_sig = (2.0 / (sqrt(M_PI) * d_gauss_sigma)) * exp(-u_sig * u_sig) * r_inv
              - erf(u_sig) / r2;
            double dg_eta = (2.0 * d_gauss_eta / sqrt(M_PI)) * exp(-u_eta * u_eta) * r_inv
              - erf(u_eta) / r2;
            dphi_dr = dg_sig - dg_eta;
          }
          double r_inv = (dist > 1.0e-10) ? 1.0 / dist : 0.0;
          double E_mag = q_node * (-dphi_dr) / (4.0 * M_PI * d_epsilon);

          E_real[0] += E_mag * dr[0] * r_inv;
          E_real[1] += E_mag * dr[1] * r_inv;
          E_real[2] += E_mag * dr[2] * r_inv;
        }
      }
    }
  }

  atomicAdd(&efield_data[3 * ludwig_idx + 0], E_real[0]);
  atomicAdd(&efield_data[3 * ludwig_idx + 1], E_real[1]);
  atomicAdd(&efield_data[3 * ludwig_idx + 2], E_real[2]);
}

/* CHANGE END - EfieldLatticGaussian */

/* -------------------------------------------------------------------------
 * ewald_particle_field_real_kernel_gaussian
 *
 * Copy of ewald_particle_field_real_kernel with Gaussian real-space field.
 * Computes the real-space part of the electric field at particle positions
 * due to (a) all lattice nodes and (b) all other particles.
 *
 * Standard:  E_mag = beta*eunit * (erfc(αr)/r² + 2α/√π·exp(-α²r²)/r) / (4πε)
 * Gaussian:  E_mag = beta*eunit * (-dphi_dr) / (4πε)
 *            where phi = erf(r/σ)/r - erf(ηr)/r
 * No eps_reg regularisation: potential is finite at r=0.
 * Taylor at r<1e-5: -dphi_dr ≈ (4/3√π)(1/σ³ - η³)·r  (from force derivation)
 * ------------------------------------------------------------------------- */

__global__ void ewald_particle_field_real_kernel_gaussian(
    double* __restrict__ Esub_data,
    double* __restrict__ fex_data,
    const double* __restrict__ particle_r,
    const double* __restrict__ particle_q,
    const double* __restrict__ rho_data,
    int nparticles,
    int nsites,
    int nhalo,
    int irc) {

  int p = blockIdx.x * blockDim.x + threadIdx.x;
  if (p >= nparticles) return;

  double q_p = particle_q[p];
  double r_p[3];
  r_p[0] = particle_r[3 * p + 0] - (double)d_noffset[0];
  r_p[1] = particle_r[3 * p + 1] - (double)d_noffset[1];
  r_p[2] = particle_r[3 * p + 2] - (double)d_noffset[2];

  int i0 = (int)floor(r_p[0]);
  int j0 = (int)floor(r_p[1]);
  int k0 = (int)floor(r_p[2]);

  double E_real[3] = { 0.0, 0.0, 0.0 };

  int str_z = 1;
  int str_y = (d_nlocal[2] + 2 * nhalo);
  int str_x = str_y * (d_nlocal[1] + 2 * nhalo);

  /* Real-space field from lattice nodes */
  for (int di = -irc; di <= irc + 1; di++) {
    for (int dj = -irc; dj <= irc + 1; dj++) {
      for (int dk = -irc; dk <= irc + 1; dk++) {
        int ni = i0 + di;
        int nj = j0 + dj;
        int nk_idx = k0 + dk;

        if (ni < 1 || ni > d_nlocal[0]) continue;
        if (nj < 1 || nj > d_nlocal[1]) continue;
        if (nk_idx < 1 || nk_idx > d_nlocal[2]) continue;

        int ludwig_idx = str_x * (nhalo + ni - 1) + str_y * (nhalo + nj - 1) + str_z * (nhalo + nk_idx - 1);

        double rho0 = rho_data[nsites * 0 + ludwig_idx];
        double rho1 = rho_data[nsites * 1 + ludwig_idx];
        double q_node = rho0 - rho1;

        if (fabs(q_node) < 1.0e-14) continue;

        double dr[3] = { r_p[0] - (double)ni, r_p[1] - (double)nj, r_p[2] - (double)nk_idx };
        double dist = sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]);

        if (dist < d_ewald_rc) {
          double dphi_dr;
          if (dist < 1.0e-5) {
            /* Taylor: -dphi_dr ≈ (4/3√π)(1/σ³ - η³)·r */
            double inv_s3 = 1.0 / (d_gauss_sigma * d_gauss_sigma * d_gauss_sigma);
            double eta3 = d_gauss_eta * d_gauss_eta * d_gauss_eta;
            dphi_dr = -(4.0 / (3.0 * sqrt(M_PI))) * (inv_s3 - eta3) * dist;
          }
          else {
            double r2 = dist * dist;
            double r_inv = 1.0 / dist;
            double u_sig = dist / d_gauss_sigma;
            double u_eta = d_gauss_eta * dist;
            /* CHANGE: d/dr[erf(r/σ)/r] - d/dr[erf(ηr)/r] */
            double dg_sig = (2.0 / (sqrt(M_PI) * d_gauss_sigma)) * exp(-u_sig * u_sig) * r_inv
              - erf(u_sig) / r2;
            double dg_eta = (2.0 * d_gauss_eta / sqrt(M_PI)) * exp(-u_eta * u_eta) * r_inv
              - erf(u_eta) / r2;
            dphi_dr = dg_sig - dg_eta;
          }
          double r_inv = (dist > 1.0e-10) ? 1.0 / dist : 0.0;
          /* CHANGE: E_mag = beta * eunit * (-dphi_dr) / (4πε) */
          double E_mag = d_beta * d_eunit * (-dphi_dr) / (4.0 * M_PI * d_epsilon);

          E_real[0] += q_node * E_mag * dr[0] * r_inv;
          E_real[1] += q_node * E_mag * dr[1] * r_inv;
          E_real[2] += q_node * E_mag * dr[2] * r_inv;
        }
      }
    }
  }

  /* Real-space field from other particles */
  for (int p2 = 0; p2 < nparticles; p2++) {
    if (p2 == p) continue;

    double q_p2 = particle_q[p2];
    if (fabs(q_p2) < 1.0e-14) continue;

    double r_p2[3];
    r_p2[0] = particle_r[3 * p2 + 0] - (double)d_noffset[0];
    r_p2[1] = particle_r[3 * p2 + 1] - (double)d_noffset[1];
    r_p2[2] = particle_r[3 * p2 + 2] - (double)d_noffset[2];

    double dr[3] = { r_p[0] - r_p2[0], r_p[1] - r_p2[1], r_p[2] - r_p2[2] };
    double dist = sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]);

    if (dist < d_ewald_rc) {
      double dphi_dr;
      if (dist < 1.0e-5) {
        double inv_s3 = 1.0 / (d_gauss_sigma * d_gauss_sigma * d_gauss_sigma);
        double eta3 = d_gauss_eta * d_gauss_eta * d_gauss_eta;
        dphi_dr = -(4.0 / (3.0 * sqrt(M_PI))) * (inv_s3 - eta3) * dist;
      }
      else {
        double r2 = dist * dist;
        double r_inv = 1.0 / dist;
        double u_sig = dist / d_gauss_sigma;
        double u_eta = d_gauss_eta * dist;
        double dg_sig = (2.0 / (sqrt(M_PI) * d_gauss_sigma)) * exp(-u_sig * u_sig) * r_inv
          - erf(u_sig) / r2;
        double dg_eta = (2.0 * d_gauss_eta / sqrt(M_PI)) * exp(-u_eta * u_eta) * r_inv
          - erf(u_eta) / r2;
        dphi_dr = dg_sig - dg_eta;
      }
      double r_inv = (dist > 1.0e-10) ? 1.0 / dist : 0.0;
      double E_mag = d_beta * d_eunit * (-dphi_dr) / (4.0 * M_PI * d_epsilon);

      E_real[0] += q_p2 * E_mag * dr[0] * r_inv;
      E_real[1] += q_p2 * E_mag * dr[1] * r_inv;
      E_real[2] += q_p2 * E_mag * dr[2] * r_inv;
    }
  }

  /* Add real-space to Fourier contribution already stored */
  Esub_data[3 * p + 0] += E_real[0];
  Esub_data[3 * p + 1] += E_real[1];
  Esub_data[3 * p + 2] += E_real[2];

  /* Update fex with real-space contribution */
  double kt = 1.0 / d_beta;
  fex_data[3 * p + 0] += q_p * E_real[0] * kt / d_eunit;
  fex_data[3 * p + 1] += q_p * E_real[1] * kt / d_eunit;
  fex_data[3 * p + 2] += q_p * E_real[2] * kt / d_eunit;
}

/* =========================================================================
 * CHANGE INIT - Gaussian_Ewald_Dual — Dual-sigma kernels (Paso 3)
 *
 * Each kernel below is the analogue of its "_gaussian" counterpart but uses
 * the appropriate σ_eff for the pair type:
 *   *_pp suffix → particle-particle pair → d_gauss_sigma_eff_pp
 *   *_ff suffix → fluid-fluid pair       → d_gauss_sigma_eff_ff
 *   *_pf suffix → particle-fluid pair    → d_gauss_sigma_eff_pf
 *
 * The Ewald splitting parameter η (d_gauss_eta) is shared across all pairs.
 * ========================================================================= */

/*CHANGE INIT - 20260611 Macro for dphi_dr with sigma->0 limit guard
 * Used by all _dual real-space kernels. When sigma<1e-14 (point charge limit):
 *   dg_sig = d/dr[1/r] = -1/r^2  (Coulomb)
 * otherwise use the standard Gaussian derivative formula.
 * DPHI_DR_DUAL(dphi_dr_, sigma_, dist_, r2_, r_inv_) sets dphi_dr_ given
 * the four precomputed values. The dist<1e-5 Taylor branch is also guarded.
 */
#define DPHI_DR_DUAL(dphi_dr_, sigma_, dist_, r2_, r_inv_)                     \
  do {                                                                          \
    if ((sigma_) < 1.0e-14) {                                                  \
      double u_eta_ = d_gauss_eta * (dist_);                                   \
      double dg_eta_ = (2.0 * d_gauss_eta / sqrt(M_PI))                       \
                         * exp(-u_eta_ * u_eta_) * (r_inv_)                   \
                       - erf(u_eta_) / (r2_);                                  \
      (dphi_dr_) = -1.0 / (r2_) - dg_eta_;                                    \
    } else if ((dist_) < 1.0e-5) {                                             \
      double inv_s3_ = 1.0 / ((sigma_) * (sigma_) * (sigma_));                \
      double eta3_   = d_gauss_eta * d_gauss_eta * d_gauss_eta;               \
      (dphi_dr_) = -(4.0 / (3.0 * sqrt(M_PI))) * (inv_s3_ - eta3_) * (dist_);\
    } else {                                                                    \
      double u_sig_ = (dist_) / (sigma_);                                      \
      double u_eta_ = d_gauss_eta * (dist_);                                   \
      double dg_sig_ = (2.0 / (sqrt(M_PI) * (sigma_)))                        \
                         * exp(-u_sig_ * u_sig_) * (r_inv_)                   \
                       - erf(u_sig_) / (r2_);                                  \
      double dg_eta_ = (2.0 * d_gauss_eta / sqrt(M_PI))                       \
                         * exp(-u_eta_ * u_eta_) * (r_inv_)                   \
                       - erf(u_eta_) / (r2_);                                  \
      (dphi_dr_) = dg_sig_ - dg_eta_;                                          \
    }                                                                           \
  } while (0)
/*CHANGE END - 20260611 */

/* φ(node) from particles — particle-fluid pair → σ_eff_pf */
__global__ void ewald_potential_real_particle_kernel_dual_pf(
    double* __restrict__ psi_data,
    const double* __restrict__ particle_r,
    const double* __restrict__ particle_q,
    int nparticles,
    int nsites,
    int nhalo,
    int irc) {

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int ntotal = d_nlocal[0] * d_nlocal[1] * d_nlocal[2];
  if (idx >= ntotal) return;

  int kc = idx % d_nlocal[2] + 1;
  int jc = (idx / d_nlocal[2]) % d_nlocal[1] + 1;
  int ic = idx / (d_nlocal[2] * d_nlocal[1]) + 1;

  int str_z = 1;
  int str_y = (d_nlocal[2] + 2 * nhalo);
  int str_x = str_y * (d_nlocal[1] + 2 * nhalo);
  int ludwig_idx = str_x * (nhalo + ic - 1) + str_y * (nhalo + jc - 1) + str_z * (nhalo + kc - 1);

  double r_local[3] = { (double)ic, (double)jc, (double)kc };
  double phi_real = 0.0;
  double sigma = d_gauss_sigma_eff_pf;

  for (int p = 0; p < nparticles; p++) {
    double q_p = particle_q[p];
    if (fabs(q_p) < 1.0e-14) continue;

    double r_p[3];
    r_p[0] = particle_r[3 * p + 0] - (double)d_noffset[0];
    r_p[1] = particle_r[3 * p + 1] - (double)d_noffset[1];
    r_p[2] = particle_r[3 * p + 2] - (double)d_noffset[2];

    double dr[3] = { r_local[0] - r_p[0], r_local[1] - r_p[1], r_local[2] - r_p[2] };
    double dist = sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]);

    if (dist < d_ewald_rc) {
      double phi_diff;
      /*CHANGE INIT - 20260611 sigma->0 guard for _pf potential */
      if (sigma < 1.0e-14) {
        double u_eta = d_gauss_eta * dist;
        phi_diff = 1.0 / dist - erf(u_eta) / dist;
      } else if (dist < 1.0e-6) {
        /*CHANGE END - 20260611 */
        phi_diff = (2.0 / sqrt(M_PI)) * (1.0 / sigma - d_gauss_eta);
      } else {
        phi_diff = erf(dist / sigma) / dist - erf(d_gauss_eta * dist) / dist;
      }
      phi_real += q_p * phi_diff / (4.0 * M_PI * d_epsilon);
    }
  }

  psi_data[ludwig_idx] += d_beta * d_eunit * phi_real;
}

/* φ(node) from other fluid nodes — fluid-fluid pair → σ_eff_ff */
__global__ void ewald_potential_real_lattice_kernel_dual_ff(
    double* __restrict__ psi_data,
    const double* __restrict__ rho_data,
    int nsites,
    int nhalo,
    int irc) {

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int ntotal = d_nlocal[0] * d_nlocal[1] * d_nlocal[2];
  if (idx >= ntotal) return;

  int kc = idx % d_nlocal[2] + 1;
  int jc = (idx / d_nlocal[2]) % d_nlocal[1] + 1;
  int ic = idx / (d_nlocal[2] * d_nlocal[1]) + 1;

  int str_z = 1;
  int str_y = (d_nlocal[2] + 2 * nhalo);
  int str_x = str_y * (d_nlocal[1] + 2 * nhalo);
  int ludwig_idx = str_x * (nhalo + ic - 1) + str_y * (nhalo + jc - 1) + str_z * (nhalo + kc - 1);

  double phi_real = 0.0;
  double sigma = d_gauss_sigma_eff_ff;

  for (int di = -irc; di <= irc; di++) {
    int i2 = ic + di;
    if (i2 < 1 - nhalo || i2 > d_nlocal[0] + nhalo) continue;
    for (int dj = -irc; dj <= irc; dj++) {
      int j2 = jc + dj;
      if (j2 < 1 - nhalo || j2 > d_nlocal[1] + nhalo) continue;
      for (int dk = -irc; dk <= irc; dk++) {
        int k2 = kc + dk;
        if (k2 < 1 - nhalo || k2 > d_nlocal[2] + nhalo) continue;
        if (di == 0 && dj == 0 && dk == 0) continue;

        int ludwig_idx2 = str_x * (nhalo + i2 - 1) + str_y * (nhalo + j2 - 1) + str_z * (nhalo + k2 - 1);
        double rho0_2 = rho_data[nsites * 0 + ludwig_idx2];
        double rho1_2 = rho_data[nsites * 1 + ludwig_idx2];
        double q_node = rho0_2 - rho1_2;
        if (fabs(q_node) < 1.0e-14) continue;

        double dist = sqrt((double)(di * di + dj * dj + dk * dk));
        if (dist < d_ewald_rc) {
          double phi_diff;
          /*CHANGE INIT - 20260611 Handle sigma_ff=0 limit: erf(r/0)/r -> 1/r */
          if (sigma < 1.0e-14) {
            phi_diff = 1.0 / dist - erf(d_gauss_eta * dist) / dist;
          } else if (dist < 1.0e-6) {
            /*CHANGE END - 20260611 */
            phi_diff = (2.0 / sqrt(M_PI)) * (1.0 / sigma - d_gauss_eta);
          } else {
            phi_diff = erf(dist / sigma) / dist - erf(d_gauss_eta * dist) / dist;
          }
          phi_real += q_node * phi_diff / (4.0 * M_PI * d_epsilon);
        }
      }
    }
  }

  psi_data[ludwig_idx] += d_beta * d_eunit * phi_real;
}

/* F(node) from particles — particle-fluid pair → σ_eff_pf */
__global__ void ewald_force_real_particle_kernel_dual_pf(
    double* __restrict__ force_data,
    const double* __restrict__ rho_data,
    const double* __restrict__ particle_r,
    const double* __restrict__ particle_q,
    int nparticles,
    int nsites,
    int nhalo,
    int irc) {

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int ntotal = d_nlocal[0] * d_nlocal[1] * d_nlocal[2];
  if (idx >= ntotal) return;

  int kc = idx % d_nlocal[2] + 1;
  int jc = (idx / d_nlocal[2]) % d_nlocal[1] + 1;
  int ic = idx / (d_nlocal[2] * d_nlocal[1]) + 1;

  int str_z = 1;
  int str_y = (d_nlocal[2] + 2 * nhalo);
  int str_x = str_y * (d_nlocal[1] + 2 * nhalo);
  int ludwig_idx = str_x * (nhalo + ic - 1) + str_y * (nhalo + jc - 1) + str_z * (nhalo + kc - 1);

  double rho0 = rho_data[nsites * 0 + ludwig_idx];
  double rho1 = rho_data[nsites * 1 + ludwig_idx];
  double q1 = rho0 - rho1;

  double r_local[3] = { (double)ic, (double)jc, (double)kc };
  double F_real[3] = { 0.0, 0.0, 0.0 };
  double sigma = d_gauss_sigma_eff_pf;

  for (int p = 0; p < nparticles; p++) {
    double q_p = particle_q[p];
    if (fabs(q_p) < 1.0e-14) continue;

    double r_p[3];
    r_p[0] = particle_r[3 * p + 0] - (double)d_noffset[0];
    r_p[1] = particle_r[3 * p + 1] - (double)d_noffset[1];
    r_p[2] = particle_r[3 * p + 2] - (double)d_noffset[2];

    double dr[3] = { r_local[0] - r_p[0], r_local[1] - r_p[1], r_local[2] - r_p[2] };
    double dist = sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]);

    if (dist < d_ewald_rc) {
      double r2 = dist * dist;
      double r_inv = (dist > 1.0e-10) ? 1.0 / dist : 0.0;
      double dphi_dr;
      DPHI_DR_DUAL(dphi_dr, sigma, dist, r2, r_inv);
      double F_mag = q1 * q_p * (-dphi_dr) / (4.0 * M_PI * d_epsilon);
      F_real[0] += F_mag * dr[0] * r_inv;
      F_real[1] += F_mag * dr[1] * r_inv;
      F_real[2] += F_mag * dr[2] * r_inv;
    }
  }

  atomicAdd(&force_data[3 * ludwig_idx + 0], F_real[0]);
  atomicAdd(&force_data[3 * ludwig_idx + 1], F_real[1]);
  atomicAdd(&force_data[3 * ludwig_idx + 2], F_real[2]);
}

/* E(node) from particles — particle-fluid pair → σ_eff_pf */
__global__ void ewald_efield_real_particle_kernel_dual_pf(
    double* __restrict__ efield_data,
    const double* __restrict__ particle_r,
    const double* __restrict__ particle_q,
    int nparticles,
    int nsites,
    int nhalo,
    int irc) {

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int ntotal = d_nlocal[0] * d_nlocal[1] * d_nlocal[2];
  if (idx >= ntotal) return;

  int kc = idx % d_nlocal[2] + 1;
  int jc = (idx / d_nlocal[2]) % d_nlocal[1] + 1;
  int ic = idx / (d_nlocal[2] * d_nlocal[1]) + 1;

  int str_z = 1;
  int str_y = (d_nlocal[2] + 2 * nhalo);
  int str_x = str_y * (d_nlocal[1] + 2 * nhalo);
  int ludwig_idx = str_x * (nhalo + ic - 1) + str_y * (nhalo + jc - 1) + str_z * (nhalo + kc - 1);

  double r_local[3] = { (double)ic, (double)jc, (double)kc };
  double E_real[3] = { 0.0, 0.0, 0.0 };
  double sigma = d_gauss_sigma_eff_pf;

  for (int p = 0; p < nparticles; p++) {
    double q_p = particle_q[p];
    if (fabs(q_p) < 1.0e-14) continue;

    double r_p[3];
    r_p[0] = particle_r[3 * p + 0] - (double)d_noffset[0];
    r_p[1] = particle_r[3 * p + 1] - (double)d_noffset[1];
    r_p[2] = particle_r[3 * p + 2] - (double)d_noffset[2];

    double dr[3] = { r_local[0] - r_p[0], r_local[1] - r_p[1], r_local[2] - r_p[2] };
    double dist = sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]);

    if (dist < d_ewald_rc) {
      double r2 = dist * dist;
      double r_inv = (dist > 1.0e-10) ? 1.0 / dist : 0.0;
      double dphi_dr;
      DPHI_DR_DUAL(dphi_dr, sigma, dist, r2, r_inv);
      /* Match efield_real_particle_kernel_gaussian: no beta·eunit factor here.
       * The field stored is the physical E divided by 4πε. */
      double E_mag = q_p * (-dphi_dr) / (4.0 * M_PI * d_epsilon);
      E_real[0] += E_mag * dr[0] * r_inv;
      E_real[1] += E_mag * dr[1] * r_inv;
      E_real[2] += E_mag * dr[2] * r_inv;
    }
  }

  atomicAdd(&efield_data[3 * ludwig_idx + 0], E_real[0]);
  atomicAdd(&efield_data[3 * ludwig_idx + 1], E_real[1]);
  atomicAdd(&efield_data[3 * ludwig_idx + 2], E_real[2]);
}

/* E(node) from other fluid nodes — fluid-fluid pair → σ_eff_ff */
__global__ void ewald_efield_real_lattice_kernel_dual_ff(
    double* __restrict__ efield_data,
    const double* __restrict__ rho_data,
    int nsites,
    int nhalo,
    int irc) {

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int ntotal = d_nlocal[0] * d_nlocal[1] * d_nlocal[2];
  if (idx >= ntotal) return;

  int kc = idx % d_nlocal[2] + 1;
  int jc = (idx / d_nlocal[2]) % d_nlocal[1] + 1;
  int ic = idx / (d_nlocal[2] * d_nlocal[1]) + 1;

  int str_z = 1;
  int str_y = (d_nlocal[2] + 2 * nhalo);
  int str_x = str_y * (d_nlocal[1] + 2 * nhalo);
  int ludwig_idx = str_x * (nhalo + ic - 1) + str_y * (nhalo + jc - 1) + str_z * (nhalo + kc - 1);

  double E_real[3] = { 0.0, 0.0, 0.0 };
  double sigma = d_gauss_sigma_eff_ff;

  for (int di = -irc; di <= irc; di++) {
    int i2 = ic + di;
    if (i2 < 1 - nhalo || i2 > d_nlocal[0] + nhalo) continue;
    for (int dj = -irc; dj <= irc; dj++) {
      int j2 = jc + dj;
      if (j2 < 1 - nhalo || j2 > d_nlocal[1] + nhalo) continue;
      for (int dk = -irc; dk <= irc; dk++) {
        int k2 = kc + dk;
        if (k2 < 1 - nhalo || k2 > d_nlocal[2] + nhalo) continue;
        if (di == 0 && dj == 0 && dk == 0) continue;

        int ludwig_idx2 = str_x * (nhalo + i2 - 1) + str_y * (nhalo + j2 - 1) + str_z * (nhalo + k2 - 1);
        double rho0_2 = rho_data[nsites * 0 + ludwig_idx2];
        double rho1_2 = rho_data[nsites * 1 + ludwig_idx2];
        double q_node = rho0_2 - rho1_2;
        if (fabs(q_node) < 1.0e-14) continue;

        double dr[3] = { -(double)di, -(double)dj, -(double)dk };
        double dist = sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]);

        if (dist < d_ewald_rc) {
          double r2 = dist * dist;
          double r_inv = (dist > 1.0e-10) ? 1.0 / dist : 0.0;
          double dphi_dr;
          DPHI_DR_DUAL(dphi_dr, sigma, dist, r2, r_inv);
          double E_mag = q_node * (-dphi_dr) / (4.0 * M_PI * d_epsilon);
          E_real[0] += E_mag * dr[0] * r_inv;
          E_real[1] += E_mag * dr[1] * r_inv;
          E_real[2] += E_mag * dr[2] * r_inv;
        }
      }
    }
  }

  atomicAdd(&efield_data[3 * ludwig_idx + 0], E_real[0]);
  atomicAdd(&efield_data[3 * ludwig_idx + 1], E_real[1]);
  atomicAdd(&efield_data[3 * ludwig_idx + 2], E_real[2]);
}

/*CHANGE INIT - 20260611 Force on fluid nodes from lattice nodes — fluid-fluid pair → σ_eff_ff */
/* F(node) from lattice nodes — fluid-fluid pair → σ_eff_ff
 *
 * Analogous to ewald_force_real_particle_kernel_dual_pf but for node-node pairs.
 * F_mag = q_receiver * q_source * (-dphi_dr) / (4π·ε)
 * Uses sigma = σ_eff_ff = √2 · σ_f.
 * When σ_f = 0: dg_sig = -1/r² (Coulomb limit, same as efield_real_lattice_dual_ff).
 */
__global__ void ewald_force_real_lattice_kernel_dual_ff(
    double* __restrict__ force_data,
    const double* __restrict__ rho_data,
    int nsites,
    int nhalo,
    int irc) {

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int ntotal = d_nlocal[0] * d_nlocal[1] * d_nlocal[2];
  if (idx >= ntotal) return;

  int kc = idx % d_nlocal[2] + 1;
  int jc = (idx / d_nlocal[2]) % d_nlocal[1] + 1;
  int ic = idx / (d_nlocal[2] * d_nlocal[1]) + 1;

  int str_z = 1;
  int str_y = (d_nlocal[2] + 2 * nhalo);
  int str_x = str_y * (d_nlocal[1] + 2 * nhalo);
  int ludwig_idx = str_x * (nhalo + ic - 1) + str_y * (nhalo + jc - 1) + str_z * (nhalo + kc - 1);

  double rho0 = rho_data[nsites * 0 + ludwig_idx];
  double rho1 = rho_data[nsites * 1 + ludwig_idx];
  double q1 = rho0 - rho1;
  if (fabs(q1) < 1.0e-14) return;

  double F_real[3] = { 0.0, 0.0, 0.0 };
  double sigma = d_gauss_sigma_eff_ff;

  for (int di = -irc; di <= irc; di++) {
    int i2 = ic + di;
    if (i2 < 1 - nhalo || i2 > d_nlocal[0] + nhalo) continue;
    for (int dj = -irc; dj <= irc; dj++) {
      int j2 = jc + dj;
      if (j2 < 1 - nhalo || j2 > d_nlocal[1] + nhalo) continue;
      for (int dk = -irc; dk <= irc; dk++) {
        int k2 = kc + dk;
        if (k2 < 1 - nhalo || k2 > d_nlocal[2] + nhalo) continue;
        if (di == 0 && dj == 0 && dk == 0) continue;

        int ludwig_idx2 = str_x * (nhalo + i2 - 1) + str_y * (nhalo + j2 - 1) + str_z * (nhalo + k2 - 1);
        double rho0_2 = rho_data[nsites * 0 + ludwig_idx2];
        double rho1_2 = rho_data[nsites * 1 + ludwig_idx2];
        double q_node = rho0_2 - rho1_2;
        if (fabs(q_node) < 1.0e-14) continue;

        double dr[3] = { -(double)di, -(double)dj, -(double)dk };
        double dist = sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]);

        if (dist < d_ewald_rc) {
          double r2 = dist * dist;
          double r_inv = (dist > 1.0e-10) ? 1.0 / dist : 0.0;
          double dphi_dr;
          DPHI_DR_DUAL(dphi_dr, sigma, dist, r2, r_inv);
          double F_mag = q1 * q_node * (-dphi_dr) / (4.0 * M_PI * d_epsilon);
          F_real[0] += F_mag * dr[0] * r_inv;
          F_real[1] += F_mag * dr[1] * r_inv;
          F_real[2] += F_mag * dr[2] * r_inv;
        }
      }
    }
  }

  atomicAdd(&force_data[3 * ludwig_idx + 0], F_real[0]);
  atomicAdd(&force_data[3 * ludwig_idx + 1], F_real[1]);
  atomicAdd(&force_data[3 * ludwig_idx + 2], F_real[2]);
}
/*CHANGE END - 20260611 */

/* E on particles — two contributions:
 *   from fluid nodes  → particle-fluid pair → σ_eff_pf
 *   from other parts  → particle-particle pair → σ_eff_pp
 */
__global__ void ewald_particle_field_real_kernel_dual(
    double* __restrict__ Esub_data,
    double* __restrict__ fex_data,
    const double* __restrict__ particle_r,
    const double* __restrict__ particle_q,
    const double* __restrict__ rho_data,
    int nparticles,
    int nsites,
    int nhalo,
    int irc) {

  int p = blockIdx.x * blockDim.x + threadIdx.x;
  if (p >= nparticles) return;

  double q_p = particle_q[p];
  double r_p[3];
  r_p[0] = particle_r[3 * p + 0] - (double)d_noffset[0];
  r_p[1] = particle_r[3 * p + 1] - (double)d_noffset[1];
  r_p[2] = particle_r[3 * p + 2] - (double)d_noffset[2];

  int i0 = (int)floor(r_p[0]);
  int j0 = (int)floor(r_p[1]);
  int k0 = (int)floor(r_p[2]);

  double E_real[3] = { 0.0, 0.0, 0.0 };

  int str_z = 1;
  int str_y = (d_nlocal[2] + 2 * nhalo);
  int str_x = str_y * (d_nlocal[1] + 2 * nhalo);

  /* Field from lattice nodes — σ_eff_pf */
  {
    double sigma = d_gauss_sigma_eff_pf;
    for (int di = -irc; di <= irc + 1; di++) {
      for (int dj = -irc; dj <= irc + 1; dj++) {
        for (int dk = -irc; dk <= irc + 1; dk++) {
          int ni = i0 + di;
          int nj = j0 + dj;
          int nk_idx = k0 + dk;
          if (ni < 1 || ni > d_nlocal[0]) continue;
          if (nj < 1 || nj > d_nlocal[1]) continue;
          if (nk_idx < 1 || nk_idx > d_nlocal[2]) continue;

          int ludwig_idx = str_x * (nhalo + ni - 1) + str_y * (nhalo + nj - 1) + str_z * (nhalo + nk_idx - 1);
          double rho0 = rho_data[nsites * 0 + ludwig_idx];
          double rho1 = rho_data[nsites * 1 + ludwig_idx];
          double q_node = rho0 - rho1;
          if (fabs(q_node) < 1.0e-14) continue;

          double dr[3] = { r_p[0] - (double)ni, r_p[1] - (double)nj, r_p[2] - (double)nk_idx };
          double dist = sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]);

          if (dist < d_ewald_rc) {
            double r2 = dist * dist;
            double r_inv = (dist > 1.0e-10) ? 1.0 / dist : 0.0;
            double dphi_dr;
            DPHI_DR_DUAL(dphi_dr, sigma, dist, r2, r_inv);
            double E_mag = d_beta * d_eunit * (-dphi_dr) / (4.0 * M_PI * d_epsilon);
            E_real[0] += q_node * E_mag * dr[0] * r_inv;
            E_real[1] += q_node * E_mag * dr[1] * r_inv;
            E_real[2] += q_node * E_mag * dr[2] * r_inv;
          }
        }
      }
    }
  }

  /* Field from other particles — σ_eff_pp */
  {
    double sigma = d_gauss_sigma_eff_pp;
    for (int p2 = 0; p2 < nparticles; p2++) {
      if (p2 == p) continue;
      double q_p2 = particle_q[p2];
      if (fabs(q_p2) < 1.0e-14) continue;

      double r_p2[3];
      r_p2[0] = particle_r[3 * p2 + 0] - (double)d_noffset[0];
      r_p2[1] = particle_r[3 * p2 + 1] - (double)d_noffset[1];
      r_p2[2] = particle_r[3 * p2 + 2] - (double)d_noffset[2];

      double dr[3] = { r_p[0] - r_p2[0], r_p[1] - r_p2[1], r_p[2] - r_p2[2] };
      double dist = sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]);

      if (dist < d_ewald_rc) {
        double r2 = dist * dist;
        double r_inv = (dist > 1.0e-10) ? 1.0 / dist : 0.0;
        double dphi_dr;
        DPHI_DR_DUAL(dphi_dr, sigma, dist, r2, r_inv);
        double E_mag = d_beta * d_eunit * (-dphi_dr) / (4.0 * M_PI * d_epsilon);
        E_real[0] += q_p2 * E_mag * dr[0] * r_inv;
        E_real[1] += q_p2 * E_mag * dr[1] * r_inv;
        E_real[2] += q_p2 * E_mag * dr[2] * r_inv;
      }
    }
  }

  Esub_data[3 * p + 0] += E_real[0];
  Esub_data[3 * p + 1] += E_real[1];
  Esub_data[3 * p + 2] += E_real[2];

  double kt = 1.0 / d_beta;
  fex_data[3 * p + 0] += q_p * E_real[0] * kt / d_eunit;
  fex_data[3 * p + 1] += q_p * E_real[1] * kt / d_eunit;
  fex_data[3 * p + 2] += q_p * E_real[2] * kt / d_eunit;
}

/* CHANGE END - Gaussian_Ewald_Dual (Paso 3) */

/* =========================================================================
 * CHANGE INIT - Gaussian_Ewald_Dual (Paso 3b) — Fourier kernels with receiver form factor
 *
 * The structure factor S(k) is built with f_p(k) for particle charges and
 * f_f(k) for lattice charges (Paso 3 already does this).  But when computing
 * the field/potential AT a node or AT a particle, we must also multiply by
 * the RECEIVER form factor: f_f(k) for a node receiver, f_p(k) for a particle
 * receiver. Without this, momentum conservation breaks whenever σ_p ≠ σ_f.
 *
 * In the single-sigma case (σ_p = σ_f = σ), f_p = f_f = exp(-k²σ²/4); the
 * existing kernels still produce correct physics because S(k) already
 * contains f² and Gk·S(k) gives the right per-receiver field.  But with
 * different sigmas the asymmetry is real.
 * ========================================================================= */

__global__ void ewald_potential_fourier_kernel_dual(
    double* __restrict__ psi_data,
    const double* __restrict__ Sk_sin,
    const double* __restrict__ Sk_cos,
    const double* __restrict__ kvec,
    const double* __restrict__ Gk,
    const double* __restrict__ fk_recv,   /* receiver form factor — f_f for nodes */
    const int* __restrict__ kz_arr,
    int nktot,
    int nsites,
    int nhalo) {

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int ntotal = d_nlocal[0] * d_nlocal[1] * d_nlocal[2];
  if (idx >= ntotal) return;

  int kc = idx % d_nlocal[2] + 1;
  int jc = (idx / d_nlocal[2]) % d_nlocal[1] + 1;
  int ic = idx / (d_nlocal[2] * d_nlocal[1]) + 1;

  int str_z = 1;
  int str_y = (d_nlocal[2] + 2 * nhalo);
  int str_x = str_y * (d_nlocal[1] + 2 * nhalo);
  int ludwig_idx = str_x * (nhalo + ic - 1) + str_y * (nhalo + jc - 1) + str_z * (nhalo + kc - 1);

  double rx = (double)(d_noffset[0] + ic);
  double ry = (double)(d_noffset[1] + jc);
  double rz = (double)(d_noffset[2] + kc);

  double phi_fourier = 0.0;

  for (int kn = 0; kn < nktot; kn++) {
    double kx = kvec[3 * kn + 0];
    double ky = kvec[3 * kn + 1];
    double kz_val = kvec[3 * kn + 2];

    double kr = kx * rx + ky * ry + kz_val * rz;
    double sinkr = sin(kr);
    double coskr = cos(kr);

    double factor = (kz_arr[kn] > 0) ? 2.0 : 1.0;
    phi_fourier += factor * Gk[kn] * fk_recv[kn]
      * (Sk_cos[kn] * coskr + Sk_sin[kn] * sinkr);
  }

  double phi_dipole = d_dipole_prefactor * (d_M_dipole[0] * rx + d_M_dipole[1] * ry + d_M_dipole[2] * rz);

  psi_data[ludwig_idx] = d_beta * d_eunit * (phi_fourier + phi_dipole);
}

__global__ void ewald_efield_fourier_kernel_dual(
    double* __restrict__ efield_data,
    const double* __restrict__ Sk_sin,
    const double* __restrict__ Sk_cos,
    const double* __restrict__ kvec,
    const double* __restrict__ Gk,
    const double* __restrict__ fk_recv,   /* f_f for nodes */
    const int* __restrict__ kz_arr,
    int nktot,
    int nsites,
    int nhalo) {

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int ntotal = d_nlocal[0] * d_nlocal[1] * d_nlocal[2];
  if (idx >= ntotal) return;

  int kc = idx % d_nlocal[2] + 1;
  int jc = (idx / d_nlocal[2]) % d_nlocal[1] + 1;
  int ic = idx / (d_nlocal[2] * d_nlocal[1]) + 1;

  int str_z = 1;
  int str_y = (d_nlocal[2] + 2 * nhalo);
  int str_x = str_y * (d_nlocal[1] + 2 * nhalo);
  int ludwig_idx = str_x * (nhalo + ic - 1) + str_y * (nhalo + jc - 1) + str_z * (nhalo + kc - 1);

  double rx = (double)(d_noffset[0] + ic);
  double ry = (double)(d_noffset[1] + jc);
  double rz = (double)(d_noffset[2] + kc);

  double E_fourier[3] = { 0.0, 0.0, 0.0 };

  for (int kn = 0; kn < nktot; kn++) {
    double kx = kvec[3 * kn + 0];
    double ky = kvec[3 * kn + 1];
    double kz_val = kvec[3 * kn + 2];

    double kr = kx * rx + ky * ry + kz_val * rz;
    double sinkr = sin(kr);
    double coskr = cos(kr);

    double factor = (kz_arr[kn] > 0) ? 2.0 : 1.0;
    double im_part = Sk_sin[kn] * coskr - Sk_cos[kn] * sinkr;
    double w = factor * Gk[kn] * fk_recv[kn];

    E_fourier[0] -= w * kx * im_part;
    E_fourier[1] -= w * ky * im_part;
    E_fourier[2] -= w * kz_val * im_part;
  }

  E_fourier[0] += d_E_dipole[0];
  E_fourier[1] += d_E_dipole[1];
  E_fourier[2] += d_E_dipole[2];

  efield_data[3 * ludwig_idx + 0] = E_fourier[0];
  efield_data[3 * ludwig_idx + 1] = E_fourier[1];
  efield_data[3 * ludwig_idx + 2] = E_fourier[2];
}

__global__ void ewald_force_fourier_kernel_dual(
    double* __restrict__ force_data,
    const double* __restrict__ rho_data,
    const double* __restrict__ Sk_sin,
    const double* __restrict__ Sk_cos,
    const double* __restrict__ kvec,
    const double* __restrict__ Gk,
    const double* __restrict__ fk_recv,   /* f_f for nodes */
    const int* __restrict__ kz_arr,
    int nktot,
    int nsites,
    int nhalo) {

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int ntotal = d_nlocal[0] * d_nlocal[1] * d_nlocal[2];
  if (idx >= ntotal) return;

  int kc = idx % d_nlocal[2] + 1;
  int jc = (idx / d_nlocal[2]) % d_nlocal[1] + 1;
  int ic = idx / (d_nlocal[2] * d_nlocal[1]) + 1;

  int str_z = 1;
  int str_y = (d_nlocal[2] + 2 * nhalo);
  int str_x = str_y * (d_nlocal[1] + 2 * nhalo);
  int ludwig_idx = str_x * (nhalo + ic - 1) + str_y * (nhalo + jc - 1) + str_z * (nhalo + kc - 1);

  double rho0 = rho_data[nsites * 0 + ludwig_idx];
  double rho1 = rho_data[nsites * 1 + ludwig_idx];
  double q1 = rho0 - rho1;

  double rx = (double)(d_noffset[0] + ic);
  double ry = (double)(d_noffset[1] + jc);
  double rz = (double)(d_noffset[2] + kc);

  double F_fourier[3] = { 0.0, 0.0, 0.0 };

  for (int kn = 0; kn < nktot; kn++) {
    double kx = kvec[3 * kn + 0];
    double ky = kvec[3 * kn + 1];
    double kz_val = kvec[3 * kn + 2];

    double kr = kx * rx + ky * ry + kz_val * rz;
    double sinkr = sin(kr);
    double coskr = cos(kr);

    double factor = (kz_arr[kn] > 0) ? 2.0 : 1.0;
    double im_part = Sk_sin[kn] * coskr - Sk_cos[kn] * sinkr;
    double w = factor * q1 * Gk[kn] * fk_recv[kn];

    F_fourier[0] += w * kx * im_part;
    F_fourier[1] += w * ky * im_part;
    F_fourier[2] += w * kz_val * im_part;
  }

  F_fourier[0] += q1 * d_E_dipole[0];
  F_fourier[1] += q1 * d_E_dipole[1];
  F_fourier[2] += q1 * d_E_dipole[2];

  atomicAdd(&force_data[3 * ludwig_idx + 0], F_fourier[0]);
  atomicAdd(&force_data[3 * ludwig_idx + 1], F_fourier[1]);
  atomicAdd(&force_data[3 * ludwig_idx + 2], F_fourier[2]);
}

__global__ void ewald_particle_field_fourier_kernel_dual(
    double* __restrict__ Esub_data,
    double* __restrict__ fex_data,
    const double* __restrict__ particle_r,
    const double* __restrict__ particle_q,
    const double* __restrict__ Sk_sin,
    const double* __restrict__ Sk_cos,
    const double* __restrict__ kvec,
    const double* __restrict__ Gk,
    const double* __restrict__ fk_recv,   /* f_p for particles */
    const int* __restrict__ kz_arr,
    int nparticles,
    int nktot) {

  int p = blockIdx.x * blockDim.x + threadIdx.x;
  if (p >= nparticles) return;

  double q_p = particle_q[p];
  double rx = particle_r[3 * p + 0];
  double ry = particle_r[3 * p + 1];
  double rz = particle_r[3 * p + 2];

  double E_fourier[3] = { 0.0, 0.0, 0.0 };

  for (int kn = 0; kn < nktot; kn++) {
    double kx = kvec[3 * kn + 0];
    double ky = kvec[3 * kn + 1];
    double kz_val = kvec[3 * kn + 2];

    double kr = kx * rx + ky * ry + kz_val * rz;
    double sinkr = sin(kr);
    double coskr = cos(kr);

    double factor = (kz_arr[kn] > 0) ? 2.0 : 1.0;
    double im_part = Sk_sin[kn] * coskr - Sk_cos[kn] * sinkr;
    double w = factor * d_beta * d_eunit * Gk[kn] * fk_recv[kn];

    E_fourier[0] += w * kx * im_part;
    E_fourier[1] += w * ky * im_part;
    E_fourier[2] += w * kz_val * im_part;
  }

  E_fourier[0] += d_beta * d_eunit * d_E_dipole[0];
  E_fourier[1] += d_beta * d_eunit * d_E_dipole[1];
  E_fourier[2] += d_beta * d_eunit * d_E_dipole[2];

  Esub_data[3 * p + 0] = E_fourier[0];
  Esub_data[3 * p + 1] = E_fourier[1];
  Esub_data[3 * p + 2] = E_fourier[2];

  double kt = 1.0 / d_beta;
  fex_data[3 * p + 0] = q_p * E_fourier[0] * kt / d_eunit;
  fex_data[3 * p + 1] = q_p * E_fourier[1] * kt / d_eunit;
  fex_data[3 * p + 2] = q_p * E_fourier[2] * kt / d_eunit;
}

/* CHANGE END - Gaussian_Ewald_Dual (Paso 3b) */

#endif /* __NVCC__ (kernels) */

/*****************************************************************************
 *
 *  ewald_charge_sum_full_gaussian_gpu
 *
 *  Gaussian Ewald sum: all charges (particles and lattice nodes) are
 *  treated as Gaussian distributions g_sigma(r) = exp(-r²/σ²)/(π^(3/2)σ³).
 *
 *  The physical interaction is  erf(r/σ)/r  instead of  1/r.
 *  Ewald splitting:  erf(r/σ)/r = [erf(r/σ)-erf(ηr)]/r + erf(ηr)/r
 *
 *  Parameters:
 *    sigma  — Gaussian width (physical, units of lattice spacing)
 *    eta    — Ewald splitting parameter (algorithmic, replaces alpha)
 *
 *  Differences vs ewald_charge_sum_full_gpu:
 *    - Precomputes gaussian_fk[kn] = exp(-k²σ²/4) and passes to kernels
 *    - Green function uses alpha (same as standard, no sigma factor in G(k))
 *    - Structure factor kernels multiply q by gaussian_fk[kn]
 *    - Real-space kernels use [erf(r/σ)-erf(αr)]/r with Taylor at r→0
 *    - No eps_reg needed (potential is finite at r=0)
 *
 *****************************************************************************/

int ewald_charge_sum_full_gaussian_gpu(ewald_charge_t* ewald, FILE* fp,
                                       double sigma) {

#ifndef __NVCC__
  pe_info(ewald->pe, "ewald_charge_sum_full_gaussian_gpu: CUDA not available.\n");
  return -1;
#else

  int nlocal[3], noffset[3];
  double ltot[3];
  double fkx, fky, fkz;
  double r4alpha_sq, b0;
  int irc;
  PI_DOUBLE(pi);

  if (ewald == NULL) return 0;
  assert(fp);

  TIMER_start(TIMER_EWALD_TOTAL);

  cs_nlocal(ewald->cs, nlocal);
  cs_nlocal_offset(ewald->cs, noffset);
  cs_ltot(ewald->cs, ltot);

  int nhalo;
  cs_nhalo(ewald->cs, &nhalo);

  irc = (int)ceil(ewald_rc_);

  fkx = 2.0 * pi / ltot[X];
  fky = 2.0 * pi / ltot[Y];
  fkz = 2.0 * pi / ltot[Z];
  r4alpha_sq = 1.0 / (4.0 * alpha_ * alpha_);
  b0 = 1.0 / (ltot[X] * ltot[Y] * ltot[Z] * epsilon_);

  int ntotal_nodes = nlocal[X] * nlocal[Y] * nlocal[Z];
  int nsites = ewald->psi->nsites;

  pe_info(ewald->pe, "\n");
  pe_info(ewald->pe, "Ewald charge sum full GAUSSIAN (GPU): sigma=%.4f alpha=%.4f\n", sigma, alpha_);
  pe_info(ewald->pe, "  Local nodes: %d x %d x %d = %d\n", nlocal[X], nlocal[Y], nlocal[Z], ntotal_nodes);
  pe_info(ewald->pe, "  Fourier terms: %d\n", nktot_);
  pe_info(ewald->pe, "  Real-space cutoff: %.2f (irc=%d)\n", ewald_rc_, irc);

  /* Copy standard constants to device */
  cudaMemcpyToSymbol(d_alpha, &alpha_, sizeof(double));  /* kept for Fourier kernels that use d_alpha */
  cudaMemcpyToSymbol(d_eps_reg, &eps_reg_, sizeof(double));
  cudaMemcpyToSymbol(d_epsilon, &epsilon_, sizeof(double));
  cudaMemcpyToSymbol(d_beta, &beta_, sizeof(double));
  cudaMemcpyToSymbol(d_eunit, &eunit_, sizeof(double));
  cudaMemcpyToSymbol(d_rpi, &rpi_, sizeof(double));
  cudaMemcpyToSymbol(d_ewald_rc, &ewald_rc_, sizeof(double));
  cudaMemcpyToSymbol(d_fkx, &fkx, sizeof(double));
  cudaMemcpyToSymbol(d_fky, &fky, sizeof(double));
  cudaMemcpyToSymbol(d_fkz, &fkz, sizeof(double));
  cudaMemcpyToSymbol(d_r4alpha_sq, &r4alpha_sq, sizeof(double));
  cudaMemcpyToSymbol(d_b0, &b0, sizeof(double));
  cudaMemcpyToSymbol(d_kmax, &kmax_, sizeof(double));
  cudaMemcpyToSymbol(d_nk, nk_, 3 * sizeof(int));
  cudaMemcpyToSymbol(d_nlocal, nlocal, 3 * sizeof(int));
  cudaMemcpyToSymbol(d_noffset, noffset, 3 * sizeof(int));

  /* copy Gaussian-specific constants */
  cudaMemcpyToSymbol(d_gauss_sigma, &sigma, sizeof(double));
  cudaMemcpyToSymbol(d_gauss_eta, &alpha_, sizeof(double));

  /* ========================================================================
   * Precompute k-vectors, Green function, and Gaussian form factors on CPU
   * ======================================================================== */

  pe_info(ewald->pe, "  [1/6] Precomputing k-vectors, Green function, gaussian_fk...\n");

  double* kvec_h = (double*)malloc(3 * nktot_ * sizeof(double));
  double* Gk_h = (double*)malloc(nktot_ * sizeof(double));
  double* gaussian_fk_h = (double*)malloc(nktot_ * sizeof(double));
  int* kz_arr_h = (int*)malloc(nktot_ * sizeof(int));

  int kn = 0;
  for (int kz = 0; kz <= nk_[Z]; kz++) {
    for (int ky = -nk_[Y]; ky <= nk_[Y]; ky++) {
      for (int kx = -nk_[X]; kx <= nk_[X]; kx++) {
        double k[3], ksq;
        k[X] = fkx * kx;
        k[Y] = fky * ky;
        k[Z] = fkz * kz;
        ksq = k[X] * k[X] + k[Y] * k[Y] + k[Z] * k[Z];
        if (ksq <= 0.0 || ksq > kmax_) continue;

        kvec_h[3 * kn + X] = k[X];
        kvec_h[3 * kn + Y] = k[Y];
        kvec_h[3 * kn + Z] = k[Z];

        /* G(k) uses alpha (r4alpha_sq), sigma factor is NOT here */
        Gk_h[kn] = b0 * exp(-r4alpha_sq * ksq) / ksq;

        /* CHANGE: Gaussian form factor fk = exp(-k²σ²/4) */
        gaussian_fk_h[kn] = exp(-ksq * sigma * sigma / 4.0);

        kz_arr_h[kn] = kz;
        kn++;
      }
    }
  }
  int nk_actual = kn;

  /* Allocate and copy to device */
  double* kvec_d, * Gk_d, * gaussian_fk_d, * Sk_sin_d, * Sk_cos_d;
  int* kz_arr_d;

  cudaMalloc(&kvec_d, 3 * nk_actual * sizeof(double));
  cudaMalloc(&Gk_d, nk_actual * sizeof(double));
  cudaMalloc(&gaussian_fk_d, nk_actual * sizeof(double));
  cudaMalloc(&kz_arr_d, nk_actual * sizeof(int));
  cudaMalloc(&Sk_sin_d, nk_actual * sizeof(double));
  cudaMalloc(&Sk_cos_d, nk_actual * sizeof(double));

  cudaMemcpy(kvec_d, kvec_h, 3 * nk_actual * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(Gk_d, Gk_h, nk_actual * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(gaussian_fk_d, gaussian_fk_h, nk_actual * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(kz_arr_d, kz_arr_h, nk_actual * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemset(Sk_sin_d, 0, nk_actual * sizeof(double));
  cudaMemset(Sk_cos_d, 0, nk_actual * sizeof(double));

  /* Get device pointers for rho and psi fields */
  double* rho_data_d = NULL;
  double* psi_data_d = NULL;

  size_t data_offset = offsetof(field_t, data);
  cudaMemcpy(&rho_data_d, (char*)(ewald->psi->rho->target) + data_offset,
             sizeof(double*), cudaMemcpyDeviceToHost);
  cudaMemcpy(&psi_data_d, (char*)(ewald->psi->psi->target) + data_offset,
             sizeof(double*), cudaMemcpyDeviceToHost);

  field_memcpy(ewald->psi->rho, tdpMemcpyHostToDevice);

  /* ========================================================================
   * Collect particle data
   * ======================================================================== */

  int nparticles = 0;
  double* particle_r_h = NULL, * particle_q_h = NULL;
  double* particle_r_d = NULL, * particle_q_d = NULL;
  double* Esub_d = NULL, * fex_d = NULL;

  if ((ewald->sources & EWALD_SOURCE_COLLOIDS) && ewald->cinfo) {
    colloid_t* pc;
    colloids_info_local_head(ewald->cinfo, &pc);
    for (; pc; pc = pc->nextlocal) nparticles++;

    if (nparticles > 0) {
      particle_r_h = (double*)malloc(3 * nparticles * sizeof(double));
      particle_q_h = (double*)malloc(nparticles * sizeof(double));

      int p = 0;
      colloids_info_local_head(ewald->cinfo, &pc);
      for (; pc; pc = pc->nextlocal) {
        particle_r_h[3 * p + 0] = pc->s.r[X];
        particle_r_h[3 * p + 1] = pc->s.r[Y];
        particle_r_h[3 * p + 2] = pc->s.r[Z];
        particle_q_h[p] = pc->s.q0 - pc->s.q1;
        p++;
      }

      cudaMalloc(&particle_r_d, 3 * nparticles * sizeof(double));
      cudaMalloc(&particle_q_d, nparticles * sizeof(double));
      cudaMalloc(&Esub_d, 3 * nparticles * sizeof(double));
      cudaMalloc(&fex_d, 3 * nparticles * sizeof(double));

      cudaMemcpy(particle_r_d, particle_r_h, 3 * nparticles * sizeof(double), cudaMemcpyHostToDevice);
      cudaMemcpy(particle_q_d, particle_q_h, nparticles * sizeof(double), cudaMemcpyHostToDevice);
    }
  }

  pe_info(ewald->pe, "  Particles: %d\n", nparticles);

  /* ========================================================================
   * PART 1: Compute structure factors (with Gaussian form factor)
   * ======================================================================== */

  pe_info(ewald->pe, "  [2/6] Computing Gaussian structure factors S(k), C(k) on GPU...\n");

  int threads = 256;
  int blocks = (ntotal_nodes + threads - 1) / threads;

  if (ewald->sources & EWALD_SOURCE_LATTICE) {
    /* CHANGE: _gaussian kernel — multiplies q by gaussian_fk[kn] */
    ewald_structure_factor_lattice_kernel_gaussian << <blocks, threads >> > (
        rho_data_d, Sk_sin_d, Sk_cos_d, kvec_d, gaussian_fk_d,
        nk_actual, nsites, nhalo);
    cudaDeviceSynchronize();
  }

  if ((ewald->sources & EWALD_SOURCE_COLLOIDS) && nparticles > 0) {
    int p_blocks = (nparticles + threads - 1) / threads;
    /* CHANGE: _gaussian kernel */
    ewald_structure_factor_particle_kernel_gaussian << <p_blocks, threads >> > (
        particle_r_d, particle_q_d, Sk_sin_d, Sk_cos_d, kvec_d, gaussian_fk_d,
        nparticles, nk_actual);
    cudaDeviceSynchronize();
  }

  /* MPI reduce structure factors (identical to standard version) */
  pe_info(ewald->pe, "  [3/6] MPI reducing structure factors...\n");

  double* Sk_sin_h = (double*)malloc(nk_actual * sizeof(double));
  double* Sk_cos_h = (double*)malloc(nk_actual * sizeof(double));
  cudaMemcpy(Sk_sin_h, Sk_sin_d, nk_actual * sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(Sk_cos_h, Sk_cos_d, nk_actual * sizeof(double), cudaMemcpyDeviceToHost);

  MPI_Comm comm;
  cs_cart_comm(ewald->cs, &comm);
  MPI_Allreduce(MPI_IN_PLACE, Sk_sin_h, nk_actual, MPI_DOUBLE, MPI_SUM, comm);
  MPI_Allreduce(MPI_IN_PLACE, Sk_cos_h, nk_actual, MPI_DOUBLE, MPI_SUM, comm);

  cudaMemcpy(Sk_sin_d, Sk_sin_h, nk_actual * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(Sk_cos_d, Sk_cos_h, nk_actual * sizeof(double), cudaMemcpyHostToDevice);

  // /* ========================================================================
  //  * Dipole moment (identical to standard version — geometry unchanged)
  //  * ======================================================================== */

  // kahan_t M_dipole[3] = { kahan_zero(), kahan_zero(), kahan_zero() };
  // kahan_t Q_total_k = kahan_zero();

  // if (nparticles > 0) {
  //   for (int p = 0; p < nparticles; p++) {
  //     double q = particle_q_h[p];
  //     kahan_add_double(&M_dipole[X], q * particle_r_h[3 * p + 0]);
  //     kahan_add_double(&M_dipole[Y], q * particle_r_h[3 * p + 1]);
  //     kahan_add_double(&M_dipole[Z], q * particle_r_h[3 * p + 2]);
  //     kahan_add_double(&Q_total_k, q);
  //   }
  // }

  // if ((ewald->sources & EWALD_SOURCE_LATTICE) && ewald->psi) {
  //   field_memcpy(ewald->psi->rho, tdpMemcpyDeviceToHost);

  //   int nk_ewald;
  //   psi_nk(ewald->psi, &nk_ewald);
  //   kahan_t* Q_species_k = (kahan_t*)calloc(nk_ewald, sizeof(kahan_t));
  //   kahan_t* Mx_species_k = (kahan_t*)calloc(nk_ewald, sizeof(kahan_t));
  //   kahan_t* My_species_k = (kahan_t*)calloc(nk_ewald, sizeof(kahan_t));
  //   kahan_t* Mz_species_k = (kahan_t*)calloc(nk_ewald, sizeof(kahan_t));
  //   for (int s = 0; s < nk_ewald; s++) {
  //     Q_species_k[s] = kahan_zero();
  //     Mx_species_k[s] = kahan_zero();
  //     My_species_k[s] = kahan_zero();
  //     Mz_species_k[s] = kahan_zero();
  //   }

  //   for (int ic = 1; ic <= nlocal[X]; ic++) {
  //     for (int jc = 1; jc <= nlocal[Y]; jc++) {
  //       for (int kc = 1; kc <= nlocal[Z]; kc++) {
  //         int index = cs_index(ewald->cs, ic, jc, kc);
  //         double rx = (double)(noffset[X] + ic);
  //         double ry = (double)(noffset[Y] + jc);
  //         double rz = (double)(noffset[Z] + kc);
  //         for (int s = 0; s < nk_ewald; s++) {
  //           int val;
  //           double rho_s;
  //           psi_valency(ewald->psi, s, &val);
  //           psi_rho(ewald->psi, index, s, &rho_s);
  //           double q_s = val * rho_s;
  //           kahan_add_double(&Q_species_k[s], q_s);
  //           kahan_add_double(&Mx_species_k[s], q_s * rx);
  //           kahan_add_double(&My_species_k[s], q_s * ry);
  //           kahan_add_double(&Mz_species_k[s], q_s * rz);
  //         }
  //       }
  //     }
  //   }

  //   for (int s = 0; s < nk_ewald; s++) {
  //     kahan_add_double(&Q_total_k, kahan_sum(&Q_species_k[s]));
  //     kahan_add_double(&M_dipole[X], kahan_sum(&Mx_species_k[s]));
  //     kahan_add_double(&M_dipole[Y], kahan_sum(&My_species_k[s]));
  //     kahan_add_double(&M_dipole[Z], kahan_sum(&Mz_species_k[s]));
  //   }
  //   free(Q_species_k); free(Mx_species_k); free(My_species_k); free(Mz_species_k);
  // }

  // {
  //   kahan_t M_reduce[4] = { M_dipole[X], M_dipole[Y], M_dipole[Z], Q_total_k };
  //   MPI_Datatype kahan_dt;
  //   MPI_Op kahan_op;
  //   kahan_mpi_datatype(&kahan_dt);
  //   kahan_mpi_op_sum(&kahan_op);
  //   MPI_Allreduce(MPI_IN_PLACE, M_reduce, 4, kahan_dt, kahan_op, comm);
  //   MPI_Type_free(&kahan_dt);
  //   MPI_Op_free(&kahan_op);
  //   M_dipole[X] = M_reduce[0];
  //   M_dipole[Y] = M_reduce[1];
  //   M_dipole[Z] = M_reduce[2];
  //   Q_total_k = M_reduce[3];
  // }

  // double M[3] = { kahan_sum(&M_dipole[X]), kahan_sum(&M_dipole[Y]), kahan_sum(&M_dipole[Z]) };
  // double Q_total = kahan_sum(&Q_total_k);
  // double V = ltot[X] * ltot[Y] * ltot[Z];

  // double dipole_prefactor = 4.0 * pi / ((1.0 + 2.0 * ewald->epsilon_prime) * V);
  // double E_dipole[3] = { -dipole_prefactor * M[X],
  //                        -dipole_prefactor * M[Y],
  //                        -dipole_prefactor * M[Z] };

  // pe_info(ewald->pe, "  Dipole moment M = (%14.7e, %14.7e, %14.7e)\n", M[X], M[Y], M[Z]);
  // pe_info(ewald->pe, "  Total charge Q = %14.7e\n", Q_total);

  // cudaMemcpyToSymbol(d_E_dipole, E_dipole, 3 * sizeof(double));
  // cudaMemcpyToSymbol(d_M_dipole, M, 3 * sizeof(double));
  // cudaMemcpyToSymbol(d_dipole_prefactor, &dipole_prefactor, sizeof(double));

  /* Zero dipole device constants so the Fourier kernel applies no dipole correction */
  {
    double zero3[3] = { 0.0, 0.0, 0.0 };
    double zero1 = 0.0;
    cudaMemcpyToSymbol(d_E_dipole, zero3, 3 * sizeof(double));
    cudaMemcpyToSymbol(d_M_dipole, zero3, 3 * sizeof(double));
    cudaMemcpyToSymbol(d_dipole_prefactor, &zero1, sizeof(double));
  }

  /* ========================================================================
   * PART 2: Potential on lattice nodes
   * ======================================================================== */

  pe_info(ewald->pe, "  [4/6] Computing Gaussian potential on lattice nodes (GPU)...\n");

  /* Fourier part — unchanged kernel (Gk already uses eta) */
  ewald_potential_fourier_kernel << <blocks, threads >> > (
      psi_data_d, Sk_sin_d, Sk_cos_d, kvec_d, Gk_d, kz_arr_d,
      nk_actual, nsites, nhalo);
  cudaDeviceSynchronize();

  /* Real-space from particles — CHANGE: Gaussian kernel */
  if (nparticles > 0) {
    ewald_potential_real_particle_kernel_gaussian << <blocks, threads >> > (
        psi_data_d, particle_r_d, particle_q_d, nparticles, nsites, nhalo, irc);
    cudaDeviceSynchronize();
  }

  /* Real-space from lattice nodes — CHANGE: Gaussian kernel */
  ewald_potential_real_lattice_kernel_gaussian << <blocks, threads >> > (
      psi_data_d, rho_data_d, nsites, nhalo, irc);
  cudaDeviceSynchronize();

  field_memcpy(ewald->psi->psi, tdpMemcpyDeviceToHost);

  /* ========================================================================
   * PART 3: Force on fluid nodes
   * ======================================================================== */

  pe_info(ewald->pe, "  [5/6] Computing Gaussian force on lattice nodes (GPU)...\n");

  kahan_t F_fluid_total[3] = { kahan_zero(), kahan_zero(), kahan_zero() };

  if (ewald->hydro != NULL) {
    hydro_memcpy(ewald->hydro, tdpMemcpyDeviceToHost);

    double* efield_d = NULL;
    double* force_d;
    cudaMalloc(&force_d, 3 * nsites * sizeof(double));
    cudaMemset(force_d, 0, 3 * nsites * sizeof(double));

    if (ewald->psi && ewald->psi->efield) {
      cudaMalloc(&efield_d, 3 * nsites * sizeof(double));
      cudaMemset(efield_d, 0, 3 * nsites * sizeof(double));

      /* Fourier efield — unchanged kernel */
      ewald_efield_fourier_kernel << <blocks, threads >> > (
          efield_d, Sk_sin_d, Sk_cos_d, kvec_d, Gk_d, kz_arr_d,
          nk_actual, nsites, nhalo);
      cudaDeviceSynchronize();

      /* Real efield from particles — CHANGE: Gaussian kernel */
      if (nparticles > 0) {
        ewald_efield_real_particle_kernel_gaussian << <blocks, threads >> > (
            efield_d, particle_r_d, particle_q_d,
            nparticles, nsites, nhalo, irc);
        cudaDeviceSynchronize();
      }

      /* CHANGE INIT - EfieldLatticGaussian - Gaussian efield nodo->nodo */
      // Original: ewald_efield_real_lattice_kernel (uses erfc, incorrect for gaussian charges)
      ewald_efield_real_lattice_kernel_gaussian << <blocks, threads >> > (
          efield_d, rho_data_d, nsites, nhalo, irc);
      cudaDeviceSynchronize();
      /* CHANGE END - EfieldLatticGaussian */
    }

    /* Fourier force — unchanged kernel */
    ewald_force_fourier_kernel << <blocks, threads >> > (
        force_d, rho_data_d, Sk_sin_d, Sk_cos_d, kvec_d, Gk_d, kz_arr_d,
        nk_actual, nsites, nhalo);
    cudaDeviceSynchronize();

    /* Real-space force from particles — CHANGE: Gaussian kernel */
    if (nparticles > 0) {
      ewald_force_real_particle_kernel_gaussian << <blocks, threads >> > (
          force_d, rho_data_d, particle_r_d, particle_q_d,
          nparticles, nsites, nhalo, irc);
      cudaDeviceSynchronize();
    }

    double* force_h = (double*)malloc(3 * nsites * sizeof(double));
    cudaMemcpy(force_h, force_d, 3 * nsites * sizeof(double), cudaMemcpyDeviceToHost);

    double* efield_h = NULL;
    if (efield_d != NULL) {
      efield_h = (double*)malloc(3 * nsites * sizeof(double));
      cudaMemcpy(efield_h, efield_d, 3 * nsites * sizeof(double), cudaMemcpyDeviceToHost);
    }

    for (int ic = 1; ic <= nlocal[X]; ic++) {
      for (int jc = 1; jc <= nlocal[Y]; jc++) {
        for (int kc = 1; kc <= nlocal[Z]; kc++) {
          int index = cs_index(ewald->cs, ic, jc, kc);
          double f[3] = { force_h[3 * index + 0], force_h[3 * index + 1], force_h[3 * index + 2] };
          hydro_f_local_add(ewald->hydro, index, f);

          if (efield_h != NULL) {
            double e_field[3] = { efield_h[3 * index + 0], efield_h[3 * index + 1], efield_h[3 * index + 2] };
            field_vector_set(ewald->psi->efield, index, e_field);
          }

          kahan_add_double(&F_fluid_total[X], f[X]);
          kahan_add_double(&F_fluid_total[Y], f[Y]);
          kahan_add_double(&F_fluid_total[Z], f[Z]);
        }
      }
    }

    if (efield_h) free(efield_h);
    if (efield_d) cudaFree(efield_d);
    free(force_h);
    cudaFree(force_d);

    if (ewald->psi && ewald->psi->efield)
      field_memcpy(ewald->psi->efield, tdpMemcpyHostToDevice);

    hydro_memcpy(ewald->hydro, tdpMemcpyHostToDevice);
  }

  {
    MPI_Datatype kahan_dt;
    MPI_Op kahan_op;
    kahan_mpi_datatype(&kahan_dt);
    kahan_mpi_op_sum(&kahan_op);
    MPI_Allreduce(MPI_IN_PLACE, F_fluid_total, 3, kahan_dt, kahan_op, comm);
    MPI_Type_free(&kahan_dt);
    MPI_Op_free(&kahan_op);
  }

  /* ========================================================================
   * PART 4: Field and force on particles
   * ======================================================================== */

  pe_info(ewald->pe, "  [6/6] Computing Gaussian field and force on particles (GPU)...\n");

  kahan_t F_particle_total[3] = { kahan_zero(), kahan_zero(), kahan_zero() };

  if (nparticles > 0) {
    int p_blocks = (nparticles + threads - 1) / threads;

    /* Fourier part — unchanged kernel */
    ewald_particle_field_fourier_kernel << <p_blocks, threads >> > (
        Esub_d, fex_d, particle_r_d, particle_q_d,
        Sk_sin_d, Sk_cos_d, kvec_d, Gk_d, kz_arr_d,
        nparticles, nk_actual);
    cudaDeviceSynchronize();

    /* Real-space part — CHANGE: Gaussian kernel */
    ewald_particle_field_real_kernel_gaussian << <p_blocks, threads >> > (
        Esub_d, fex_d, particle_r_d, particle_q_d, rho_data_d,
        nparticles, nsites, nhalo, irc);
    cudaDeviceSynchronize();

    double* Esub_h = (double*)malloc(3 * nparticles * sizeof(double));
    double* fex_h = (double*)malloc(3 * nparticles * sizeof(double));
    cudaMemcpy(Esub_h, Esub_d, 3 * nparticles * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(fex_h, fex_d, 3 * nparticles * sizeof(double), cudaMemcpyDeviceToHost);

    int p = 0;
    double kt = 1.0 / beta_;
    colloid_t* pc;
    colloids_info_local_head(ewald->cinfo, &pc);
    for (; pc; pc = pc->nextlocal) {
      pc->Esub[X] = Esub_h[3 * p + 0];
      pc->Esub[Y] = Esub_h[3 * p + 1];
      pc->Esub[Z] = Esub_h[3 * p + 2];
      pc->fex[X] = fex_h[3 * p + 0];
      pc->fex[Y] = fex_h[3 * p + 1];
      pc->fex[Z] = fex_h[3 * p + 2];
      kahan_add_double(&F_particle_total[X], pc->fex[X]);
      kahan_add_double(&F_particle_total[Y], pc->fex[Y]);
      kahan_add_double(&F_particle_total[Z], pc->fex[Z]);
      p++;
    }
    (void)kt;  /* suppress unused warning; kt used in standard version for Esub→fex */

    free(Esub_h);
    free(fex_h);
  }

  {
    MPI_Datatype kahan_dt;
    MPI_Op kahan_op;
    kahan_mpi_datatype(&kahan_dt);
    kahan_mpi_op_sum(&kahan_op);
    MPI_Allreduce(MPI_IN_PLACE, F_particle_total, 3, kahan_dt, kahan_op, comm);
    MPI_Type_free(&kahan_dt);
    MPI_Op_free(&kahan_op);
  }

  double F_fluid[3] = { kahan_sum(&F_fluid_total[X]),    kahan_sum(&F_fluid_total[Y]),    kahan_sum(&F_fluid_total[Z]) };
  double F_particle[3] = { kahan_sum(&F_particle_total[X]), kahan_sum(&F_particle_total[Y]), kahan_sum(&F_particle_total[Z]) };
  double F_diff[3] = { F_fluid[X] + F_particle[X], F_fluid[Y] + F_particle[Y], F_fluid[Z] + F_particle[Z] };

  fprintf(fp, "%1.15e,%1.15e, %1.15e, %1.15e,%1.15e, %1.15e, %1.15e,%1.15e, %1.15e, %1.15e,\n",
          sqrt(F_diff[X] * F_diff[X] + F_diff[Y] * F_diff[Y] + F_diff[Z] * F_diff[Z]),
          F_diff[X], F_diff[Y], F_diff[Z],
          F_fluid[X], F_fluid[Y], F_fluid[Z],
          F_particle[X], F_particle[Y], F_particle[Z]);

  pe_info(ewald->pe, "\n");
  pe_info(ewald->pe, "  Momentum conservation check (Gaussian Ewald):\n");
  pe_info(ewald->pe, "    F_fluid    = (%14.7e, %14.7e, %14.7e)\n", F_fluid[X], F_fluid[Y], F_fluid[Z]);
  pe_info(ewald->pe, "    F_particle = (%14.7e, %14.7e, %14.7e)\n", F_particle[X], F_particle[Y], F_particle[Z]);
  pe_info(ewald->pe, "    F_total    = (%14.7e, %14.7e, %14.7e)\n", F_diff[X], F_diff[Y], F_diff[Z]);
  pe_info(ewald->pe, "    |F_total|  = %14.7e\n",
          sqrt(F_diff[X] * F_diff[X] + F_diff[Y] * F_diff[Y] + F_diff[Z] * F_diff[Z]));

  ewald_charge_external_field(ewald);

  /* ========================================================================
   * Cleanup
   * ======================================================================== */

  free(kvec_h);
  free(Gk_h);
  free(gaussian_fk_h);
  free(kz_arr_h);
  free(Sk_sin_h);
  free(Sk_cos_h);

  cudaFree(kvec_d);
  cudaFree(Gk_d);
  cudaFree(gaussian_fk_d);
  cudaFree(kz_arr_d);
  cudaFree(Sk_sin_d);
  cudaFree(Sk_cos_d);

  if (nparticles > 0) {
    free(particle_r_h);
    free(particle_q_h);
    cudaFree(particle_r_d);
    cudaFree(particle_q_d);
    cudaFree(Esub_d);
    cudaFree(fex_d);
  }

  pe_info(ewald->pe, "Ewald charge sum full GAUSSIAN (GPU): complete.\n\n");

  TIMER_stop(TIMER_EWALD_TOTAL);

  return 0;

#endif /* __NVCC__ */
}

/*****************************************************************************
 *
 *  ewald_charge_sum_full_gaussian_dual_gpu
 *
 *  CHANGE INIT - Gaussian_Ewald_Dual
 *
 *  Same as ewald_charge_sum_full_gaussian_gpu but with independent Gaussian
 *  widths for particles (σ_p) and fluid nodes (σ_f).
 *
 *  Pair interaction widths (convolution of two Gaussians):
 *    σ_eff_pp = sqrt(2) · σ_p     (particle-particle)
 *    σ_eff_ff = sqrt(2) · σ_f     (fluid-fluid)
 *    σ_eff_pf = sqrt(σ_p² + σ_f²) (particle-fluid)
 *
 *  Fourier-space form factors (per-charge, used in structure factor):
 *    f_p(k) = exp(-k² σ_p² / 4)   for particle charges
 *    f_f(k) = exp(-k² σ_f² / 4)   for fluid node charges
 *  The product fp·ff naturally gives the pair convolution exp(-k² σ_eff_pf²/2).
 *
 *  Real-space cutoff uses ewald_rc_ for all pairs.  Each pair uses its own
 *  σ_eff in the smoothing term [erf(r/σ_eff) - erf(η r)] / r.
 *
 *  Status: STUB — kernels and main loop pending (Pasos 2-4).
 *
 *****************************************************************************/

int ewald_charge_sum_full_gaussian_dual_gpu(ewald_charge_t* ewald, FILE* fp,
                                            double sigma_p, double sigma_f) {

#ifndef __NVCC__
  pe_info(ewald->pe, "ewald_charge_sum_full_gaussian_dual_gpu: CUDA not available.\n");
  return -1;
#else

  int nlocal[3], noffset[3];
  double ltot[3];
  double fkx, fky, fkz;
  double r4alpha_sq, b0;
  int irc;
  PI_DOUBLE(pi);

  if (ewald == NULL) return 0;
  assert(fp);

  TIMER_start(TIMER_EWALD_TOTAL);

  cs_nlocal(ewald->cs, nlocal);
  cs_nlocal_offset(ewald->cs, noffset);
  cs_ltot(ewald->cs, ltot);

  int nhalo;
  cs_nhalo(ewald->cs, &nhalo);

  irc = (int)ceil(ewald_rc_);

  fkx = 2.0 * pi / ltot[X];
  fky = 2.0 * pi / ltot[Y];
  fkz = 2.0 * pi / ltot[Z];
  r4alpha_sq = 1.0 / (4.0 * alpha_ * alpha_);
  b0 = 1.0 / (ltot[X] * ltot[Y] * ltot[Z] * epsilon_);

  int ntotal_nodes = nlocal[X] * nlocal[Y] * nlocal[Z];
  int nsites = ewald->psi->nsites;

  /* Pair-effective widths */
  double sigma_eff_pp = sqrt(2.0) * sigma_p;
  double sigma_eff_ff = sqrt(2.0) * sigma_f;
  double sigma_eff_pf = sqrt(sigma_p * sigma_p + sigma_f * sigma_f);

  pe_info(ewald->pe, "\n");
  pe_info(ewald->pe, "Ewald charge sum full GAUSSIAN DUAL (GPU):\n");
  pe_info(ewald->pe, "  sigma_p=%.4f sigma_f=%.4f alpha=%.4f\n",
          sigma_p, sigma_f, alpha_);
  pe_info(ewald->pe, "  sigma_eff: pp=%.4f ff=%.4f pf=%.4f\n",
          sigma_eff_pp, sigma_eff_ff, sigma_eff_pf);
  pe_info(ewald->pe, "  Local nodes: %d x %d x %d = %d\n",
          nlocal[X], nlocal[Y], nlocal[Z], ntotal_nodes);
  pe_info(ewald->pe, "  Fourier terms: %d\n", nktot_);
  pe_info(ewald->pe, "  Real-space cutoff: %.2f (irc=%d)\n", ewald_rc_, irc);

  /* Copy standard constants to device */
  cudaMemcpyToSymbol(d_alpha, &alpha_, sizeof(double));
  cudaMemcpyToSymbol(d_eps_reg, &eps_reg_, sizeof(double));
  cudaMemcpyToSymbol(d_epsilon, &epsilon_, sizeof(double));
  cudaMemcpyToSymbol(d_beta, &beta_, sizeof(double));
  cudaMemcpyToSymbol(d_eunit, &eunit_, sizeof(double));
  cudaMemcpyToSymbol(d_rpi, &rpi_, sizeof(double));
  cudaMemcpyToSymbol(d_ewald_rc, &ewald_rc_, sizeof(double));
  cudaMemcpyToSymbol(d_fkx, &fkx, sizeof(double));
  cudaMemcpyToSymbol(d_fky, &fky, sizeof(double));
  cudaMemcpyToSymbol(d_fkz, &fkz, sizeof(double));
  cudaMemcpyToSymbol(d_r4alpha_sq, &r4alpha_sq, sizeof(double));
  cudaMemcpyToSymbol(d_b0, &b0, sizeof(double));
  cudaMemcpyToSymbol(d_kmax, &kmax_, sizeof(double));
  cudaMemcpyToSymbol(d_nk, nk_, 3 * sizeof(int));
  cudaMemcpyToSymbol(d_nlocal, nlocal, 3 * sizeof(int));
  cudaMemcpyToSymbol(d_noffset, noffset, 3 * sizeof(int));

  /* Dual-sigma device constants */
  cudaMemcpyToSymbol(d_gauss_sigma_p,      &sigma_p,      sizeof(double));
  cudaMemcpyToSymbol(d_gauss_sigma_f,      &sigma_f,      sizeof(double));
  cudaMemcpyToSymbol(d_gauss_sigma_eff_pp, &sigma_eff_pp, sizeof(double));
  cudaMemcpyToSymbol(d_gauss_sigma_eff_ff, &sigma_eff_ff, sizeof(double));
  cudaMemcpyToSymbol(d_gauss_sigma_eff_pf, &sigma_eff_pf, sizeof(double));
  cudaMemcpyToSymbol(d_gauss_eta,          &alpha_,       sizeof(double));

  /* ========================================================================
   * Precompute k-vectors, Green function, and two form-factor arrays
   * ======================================================================== */

  pe_info(ewald->pe, "  [1/6] Precomputing k-vectors, G(k), gaussian_fk_p, gaussian_fk_f...\n");

  double* kvec_h          = (double*)malloc(3 * nktot_ * sizeof(double));
  double* Gk_h            = (double*)malloc(    nktot_ * sizeof(double));
  double* gaussian_fk_p_h = (double*)malloc(    nktot_ * sizeof(double));
  double* gaussian_fk_f_h = (double*)malloc(    nktot_ * sizeof(double));
  int*    kz_arr_h        = (int*)   malloc(    nktot_ * sizeof(int));

  int kn = 0;
  for (int kz = 0; kz <= nk_[Z]; kz++) {
    for (int ky = -nk_[Y]; ky <= nk_[Y]; ky++) {
      for (int kx = -nk_[X]; kx <= nk_[X]; kx++) {
        double k[3], ksq;
        k[X] = fkx * kx;
        k[Y] = fky * ky;
        k[Z] = fkz * kz;
        ksq = k[X] * k[X] + k[Y] * k[Y] + k[Z] * k[Z];
        if (ksq <= 0.0 || ksq > kmax_) continue;

        kvec_h[3 * kn + X] = k[X];
        kvec_h[3 * kn + Y] = k[Y];
        kvec_h[3 * kn + Z] = k[Z];

        Gk_h[kn] = b0 * exp(-r4alpha_sq * ksq) / ksq;

        /* Per-species form factors */
        gaussian_fk_p_h[kn] = exp(-ksq * sigma_p * sigma_p / 4.0);
        gaussian_fk_f_h[kn] = exp(-ksq * sigma_f * sigma_f / 4.0);

        kz_arr_h[kn] = kz;
        kn++;
      }
    }
  }
  int nk_actual = kn;

  double* kvec_d, * Gk_d, * gaussian_fk_p_d, * gaussian_fk_f_d, * Sk_sin_d, * Sk_cos_d;
  int* kz_arr_d;

  cudaMalloc(&kvec_d,          3 * nk_actual * sizeof(double));
  cudaMalloc(&Gk_d,                nk_actual * sizeof(double));
  cudaMalloc(&gaussian_fk_p_d,     nk_actual * sizeof(double));
  cudaMalloc(&gaussian_fk_f_d,     nk_actual * sizeof(double));
  cudaMalloc(&kz_arr_d,            nk_actual * sizeof(int));
  cudaMalloc(&Sk_sin_d,            nk_actual * sizeof(double));
  cudaMalloc(&Sk_cos_d,            nk_actual * sizeof(double));

  cudaMemcpy(kvec_d, kvec_h, 3 * nk_actual * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(Gk_d,   Gk_h,       nk_actual * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(gaussian_fk_p_d, gaussian_fk_p_h, nk_actual * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(gaussian_fk_f_d, gaussian_fk_f_h, nk_actual * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(kz_arr_d, kz_arr_h, nk_actual * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemset(Sk_sin_d, 0, nk_actual * sizeof(double));
  cudaMemset(Sk_cos_d, 0, nk_actual * sizeof(double));

  /* Device pointers for rho/psi */
  double* rho_data_d = NULL;
  double* psi_data_d = NULL;
  size_t data_offset = offsetof(field_t, data);
  cudaMemcpy(&rho_data_d, (char*)(ewald->psi->rho->target) + data_offset,
             sizeof(double*), cudaMemcpyDeviceToHost);
  cudaMemcpy(&psi_data_d, (char*)(ewald->psi->psi->target) + data_offset,
             sizeof(double*), cudaMemcpyDeviceToHost);

  field_memcpy(ewald->psi->rho, tdpMemcpyHostToDevice);

  /* ========================================================================
   * Collect particle data
   * ======================================================================== */

  int nparticles = 0;
  double* particle_r_h = NULL, * particle_q_h = NULL;
  double* particle_r_d = NULL, * particle_q_d = NULL;
  double* Esub_d = NULL, * fex_d = NULL;

  if ((ewald->sources & EWALD_SOURCE_COLLOIDS) && ewald->cinfo) {
    colloid_t* pc;
    colloids_info_local_head(ewald->cinfo, &pc);
    for (; pc; pc = pc->nextlocal) nparticles++;

    if (nparticles > 0) {
      particle_r_h = (double*)malloc(3 * nparticles * sizeof(double));
      particle_q_h = (double*)malloc(    nparticles * sizeof(double));

      int p = 0;
      colloids_info_local_head(ewald->cinfo, &pc);
      for (; pc; pc = pc->nextlocal) {
        particle_r_h[3 * p + 0] = pc->s.r[X];
        particle_r_h[3 * p + 1] = pc->s.r[Y];
        particle_r_h[3 * p + 2] = pc->s.r[Z];
        particle_q_h[p] = pc->s.q0 - pc->s.q1;
        p++;
      }

      cudaMalloc(&particle_r_d, 3 * nparticles * sizeof(double));
      cudaMalloc(&particle_q_d,     nparticles * sizeof(double));
      cudaMalloc(&Esub_d,       3 * nparticles * sizeof(double));
      cudaMalloc(&fex_d,        3 * nparticles * sizeof(double));

      cudaMemcpy(particle_r_d, particle_r_h, 3 * nparticles * sizeof(double), cudaMemcpyHostToDevice);
      cudaMemcpy(particle_q_d, particle_q_h,     nparticles * sizeof(double), cudaMemcpyHostToDevice);
    }
  }

  pe_info(ewald->pe, "  Particles: %d\n", nparticles);

  /* ========================================================================
   * PART 1: Structure factors with per-species form factors
   *   Lattice nodes use gaussian_fk_f (σ_f), particles use gaussian_fk_p (σ_p)
   * ======================================================================== */

  pe_info(ewald->pe, "  [2/6] Computing dual-sigma structure factors S(k), C(k)...\n");

  int threads = 256;
  int blocks = (ntotal_nodes + threads - 1) / threads;

  if (ewald->sources & EWALD_SOURCE_LATTICE) {
    ewald_structure_factor_lattice_kernel_gaussian << <blocks, threads >> > (
        rho_data_d, Sk_sin_d, Sk_cos_d, kvec_d, gaussian_fk_f_d,
        nk_actual, nsites, nhalo);
    cudaDeviceSynchronize();
  }

  if ((ewald->sources & EWALD_SOURCE_COLLOIDS) && nparticles > 0) {
    int p_blocks = (nparticles + threads - 1) / threads;
    ewald_structure_factor_particle_kernel_gaussian << <p_blocks, threads >> > (
        particle_r_d, particle_q_d, Sk_sin_d, Sk_cos_d, kvec_d, gaussian_fk_p_d,
        nparticles, nk_actual);
    cudaDeviceSynchronize();
  }

  /* MPI reduce structure factors */
  pe_info(ewald->pe, "  [3/6] MPI reducing structure factors...\n");

  double* Sk_sin_h = (double*)malloc(nk_actual * sizeof(double));
  double* Sk_cos_h = (double*)malloc(nk_actual * sizeof(double));
  cudaMemcpy(Sk_sin_h, Sk_sin_d, nk_actual * sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(Sk_cos_h, Sk_cos_d, nk_actual * sizeof(double), cudaMemcpyDeviceToHost);

  MPI_Comm comm;
  cs_cart_comm(ewald->cs, &comm);
  MPI_Allreduce(MPI_IN_PLACE, Sk_sin_h, nk_actual, MPI_DOUBLE, MPI_SUM, comm);
  MPI_Allreduce(MPI_IN_PLACE, Sk_cos_h, nk_actual, MPI_DOUBLE, MPI_SUM, comm);

  cudaMemcpy(Sk_sin_d, Sk_sin_h, nk_actual * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(Sk_cos_d, Sk_cos_h, nk_actual * sizeof(double), cudaMemcpyHostToDevice);

  /* Zero dipole device constants (no dipole correction in this variant) */
  {
    double zero3[3] = { 0.0, 0.0, 0.0 };
    double zero1 = 0.0;
    cudaMemcpyToSymbol(d_E_dipole, zero3, 3 * sizeof(double));
    cudaMemcpyToSymbol(d_M_dipole, zero3, 3 * sizeof(double));
    cudaMemcpyToSymbol(d_dipole_prefactor, &zero1, sizeof(double));
  }

  /* ========================================================================
   * PART 2: Potential on lattice nodes
   *   Real-space contributions use σ_eff_pf (particle→node) and σ_eff_ff (node→node)
   * ======================================================================== */

  pe_info(ewald->pe, "  [4/6] Computing dual-sigma potential on lattice nodes...\n");

  /* Fourier part — dual kernel: multiplies by f_f (node receiver) */
  ewald_potential_fourier_kernel_dual << <blocks, threads >> > (
      psi_data_d, Sk_sin_d, Sk_cos_d, kvec_d, Gk_d, gaussian_fk_f_d, kz_arr_d,
      nk_actual, nsites, nhalo);
  cudaDeviceSynchronize();

  if (nparticles > 0) {
    ewald_potential_real_particle_kernel_dual_pf << <blocks, threads >> > (
        psi_data_d, particle_r_d, particle_q_d, nparticles, nsites, nhalo, irc);
    cudaDeviceSynchronize();
  }

  ewald_potential_real_lattice_kernel_dual_ff << <blocks, threads >> > (
      psi_data_d, rho_data_d, nsites, nhalo, irc);
  cudaDeviceSynchronize();

  field_memcpy(ewald->psi->psi, tdpMemcpyDeviceToHost);

  /* ========================================================================
   * PART 3: Force and E-field on fluid nodes
   * ======================================================================== */

  pe_info(ewald->pe, "  [5/6] Computing dual-sigma force on lattice nodes...\n");

  kahan_t F_fluid_total[3] = { kahan_zero(), kahan_zero(), kahan_zero() };

  if (ewald->hydro != NULL) {
    hydro_memcpy(ewald->hydro, tdpMemcpyDeviceToHost);

    double* efield_d = NULL;
    double* force_d;
    cudaMalloc(&force_d, 3 * nsites * sizeof(double));
    cudaMemset(force_d, 0, 3 * nsites * sizeof(double));

    if (ewald->psi && ewald->psi->efield) {
      cudaMalloc(&efield_d, 3 * nsites * sizeof(double));
      cudaMemset(efield_d, 0, 3 * nsites * sizeof(double));

      /* Fourier efield — dual kernel: multiplies by f_f (node receiver) */
      ewald_efield_fourier_kernel_dual << <blocks, threads >> > (
          efield_d, Sk_sin_d, Sk_cos_d, kvec_d, Gk_d, gaussian_fk_f_d, kz_arr_d,
          nk_actual, nsites, nhalo);
      cudaDeviceSynchronize();

      if (nparticles > 0) {
        ewald_efield_real_particle_kernel_dual_pf << <blocks, threads >> > (
            efield_d, particle_r_d, particle_q_d,
            nparticles, nsites, nhalo, irc);
        cudaDeviceSynchronize();
      }

      ewald_efield_real_lattice_kernel_dual_ff << <blocks, threads >> > (
          efield_d, rho_data_d, nsites, nhalo, irc);
      cudaDeviceSynchronize();
    }

    /* Fourier force — dual kernel: multiplies by f_f (node receiver) */
    ewald_force_fourier_kernel_dual << <blocks, threads >> > (
        force_d, rho_data_d, Sk_sin_d, Sk_cos_d, kvec_d, Gk_d, gaussian_fk_f_d, kz_arr_d,
        nk_actual, nsites, nhalo);
    cudaDeviceSynchronize();

    if (nparticles > 0) {
      ewald_force_real_particle_kernel_dual_pf << <blocks, threads >> > (
          force_d, rho_data_d, particle_r_d, particle_q_d,
          nparticles, nsites, nhalo, irc);
      cudaDeviceSynchronize();
    }

    /*CHANGE INIT - 20260611 Real-space force from lattice nodes (fluid-fluid pair, σ_eff_ff) */
    ewald_force_real_lattice_kernel_dual_ff << <blocks, threads >> > (
        force_d, rho_data_d, nsites, nhalo, irc);
    cudaDeviceSynchronize();
    /*CHANGE END - 20260611 */

    double* force_h = (double*)malloc(3 * nsites * sizeof(double));
    cudaMemcpy(force_h, force_d, 3 * nsites * sizeof(double), cudaMemcpyDeviceToHost);

    double* efield_h = NULL;
    if (efield_d != NULL) {
      efield_h = (double*)malloc(3 * nsites * sizeof(double));
      cudaMemcpy(efield_h, efield_d, 3 * nsites * sizeof(double), cudaMemcpyDeviceToHost);
    }

    for (int ic = 1; ic <= nlocal[X]; ic++) {
      for (int jc = 1; jc <= nlocal[Y]; jc++) {
        for (int kc = 1; kc <= nlocal[Z]; kc++) {
          int index = cs_index(ewald->cs, ic, jc, kc);
          double f[3] = { force_h[3 * index + 0], force_h[3 * index + 1], force_h[3 * index + 2] };
          hydro_f_local_add(ewald->hydro, index, f);

          if (efield_h != NULL) {
            double e_field[3] = { efield_h[3 * index + 0], efield_h[3 * index + 1], efield_h[3 * index + 2] };
            field_vector_set(ewald->psi->efield, index, e_field);
          }

          kahan_add_double(&F_fluid_total[X], f[X]);
          kahan_add_double(&F_fluid_total[Y], f[Y]);
          kahan_add_double(&F_fluid_total[Z], f[Z]);
        }
      }
    }

    if (efield_h) free(efield_h);
    if (efield_d) cudaFree(efield_d);
    free(force_h);
    cudaFree(force_d);

    if (ewald->psi && ewald->psi->efield)
      field_memcpy(ewald->psi->efield, tdpMemcpyHostToDevice);

    hydro_memcpy(ewald->hydro, tdpMemcpyHostToDevice);
  }

  {
    MPI_Datatype kahan_dt;
    MPI_Op kahan_op;
    kahan_mpi_datatype(&kahan_dt);
    kahan_mpi_op_sum(&kahan_op);
    MPI_Allreduce(MPI_IN_PLACE, F_fluid_total, 3, kahan_dt, kahan_op, comm);
    MPI_Type_free(&kahan_dt);
    MPI_Op_free(&kahan_op);
  }

  /* ========================================================================
   * PART 4: Field and force on particles (uses σ_eff_pf for nodes, σ_eff_pp for parts)
   * ======================================================================== */

  pe_info(ewald->pe, "  [6/6] Computing dual-sigma field and force on particles...\n");

  kahan_t F_particle_total[3] = { kahan_zero(), kahan_zero(), kahan_zero() };

  if (nparticles > 0) {
    int p_blocks = (nparticles + threads - 1) / threads;

    /* Fourier part — dual kernel: multiplies by f_p (particle receiver) */
    ewald_particle_field_fourier_kernel_dual << <p_blocks, threads >> > (
        Esub_d, fex_d, particle_r_d, particle_q_d,
        Sk_sin_d, Sk_cos_d, kvec_d, Gk_d, gaussian_fk_p_d, kz_arr_d,
        nparticles, nk_actual);
    cudaDeviceSynchronize();

    ewald_particle_field_real_kernel_dual << <p_blocks, threads >> > (
        Esub_d, fex_d, particle_r_d, particle_q_d, rho_data_d,
        nparticles, nsites, nhalo, irc);
    cudaDeviceSynchronize();

    double* Esub_h = (double*)malloc(3 * nparticles * sizeof(double));
    double* fex_h  = (double*)malloc(3 * nparticles * sizeof(double));
    cudaMemcpy(Esub_h, Esub_d, 3 * nparticles * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(fex_h,  fex_d,  3 * nparticles * sizeof(double), cudaMemcpyDeviceToHost);

    int p = 0;
    colloid_t* pc;
    colloids_info_local_head(ewald->cinfo, &pc);
    for (; pc; pc = pc->nextlocal) {
      pc->Esub[X] = Esub_h[3 * p + 0];
      pc->Esub[Y] = Esub_h[3 * p + 1];
      pc->Esub[Z] = Esub_h[3 * p + 2];
      pc->fex[X]  = fex_h[3 * p + 0];
      pc->fex[Y]  = fex_h[3 * p + 1];
      pc->fex[Z]  = fex_h[3 * p + 2];
      kahan_add_double(&F_particle_total[X], pc->fex[X]);
      kahan_add_double(&F_particle_total[Y], pc->fex[Y]);
      kahan_add_double(&F_particle_total[Z], pc->fex[Z]);
      p++;
    }

    free(Esub_h);
    free(fex_h);
  }

  {
    MPI_Datatype kahan_dt;
    MPI_Op kahan_op;
    kahan_mpi_datatype(&kahan_dt);
    kahan_mpi_op_sum(&kahan_op);
    MPI_Allreduce(MPI_IN_PLACE, F_particle_total, 3, kahan_dt, kahan_op, comm);
    MPI_Type_free(&kahan_dt);
    MPI_Op_free(&kahan_op);
  }

  double F_fluid[3]    = { kahan_sum(&F_fluid_total[X]),    kahan_sum(&F_fluid_total[Y]),    kahan_sum(&F_fluid_total[Z]) };
  double F_particle[3] = { kahan_sum(&F_particle_total[X]), kahan_sum(&F_particle_total[Y]), kahan_sum(&F_particle_total[Z]) };
  double F_diff[3]     = { F_fluid[X] + F_particle[X], F_fluid[Y] + F_particle[Y], F_fluid[Z] + F_particle[Z] };

  fprintf(fp, "%1.15e,%1.15e, %1.15e, %1.15e,%1.15e, %1.15e, %1.15e,%1.15e, %1.15e, %1.15e,\n",
          sqrt(F_diff[X] * F_diff[X] + F_diff[Y] * F_diff[Y] + F_diff[Z] * F_diff[Z]),
          F_diff[X], F_diff[Y], F_diff[Z],
          F_fluid[X], F_fluid[Y], F_fluid[Z],
          F_particle[X], F_particle[Y], F_particle[Z]);

  pe_info(ewald->pe, "\n");
  pe_info(ewald->pe, "  Momentum conservation check (Gaussian Dual Ewald):\n");
  pe_info(ewald->pe, "    F_fluid    = (%14.7e, %14.7e, %14.7e)\n", F_fluid[X], F_fluid[Y], F_fluid[Z]);
  pe_info(ewald->pe, "    F_particle = (%14.7e, %14.7e, %14.7e)\n", F_particle[X], F_particle[Y], F_particle[Z]);
  pe_info(ewald->pe, "    F_total    = (%14.7e, %14.7e, %14.7e)\n", F_diff[X], F_diff[Y], F_diff[Z]);
  pe_info(ewald->pe, "    |F_total|  = %14.7e\n",
          sqrt(F_diff[X] * F_diff[X] + F_diff[Y] * F_diff[Y] + F_diff[Z] * F_diff[Z]));

  ewald_charge_external_field(ewald);

  /* Cleanup */
  free(kvec_h);
  free(Gk_h);
  free(gaussian_fk_p_h);
  free(gaussian_fk_f_h);
  free(kz_arr_h);
  free(Sk_sin_h);
  free(Sk_cos_h);

  cudaFree(kvec_d);
  cudaFree(Gk_d);
  cudaFree(gaussian_fk_p_d);
  cudaFree(gaussian_fk_f_d);
  cudaFree(kz_arr_d);
  cudaFree(Sk_sin_d);
  cudaFree(Sk_cos_d);

  if (nparticles > 0) {
    free(particle_r_h);
    free(particle_q_h);
    cudaFree(particle_r_d);
    cudaFree(particle_q_d);
    cudaFree(Esub_d);
    cudaFree(fex_d);
  }

  pe_info(ewald->pe, "Ewald charge sum full GAUSSIAN DUAL (GPU): complete.\n\n");

  TIMER_stop(TIMER_EWALD_TOTAL);

  return 0;

#endif /* __NVCC__ */
}

/* CHANGE END - Gaussian_Ewald_Dual */

/*****************************************************************************
 *
 *  ewald_flag_near_particle_nodes
 *
 *  Marks nodes where a particle is close to a cube face (diagnostic only).
 *  EWALD_NODE_NEAR_PARTICLE means the quadrature on that face may be
 *  less accurate, but the force is still computed (not excluded).
 *
 *  A particle at fractional local position (rpx, rpy, rpz) has distance
 *  to the closest face of cube centered at node (i,j,k):
 *    dist_to_face = 0.5 - max(|rpx-i|, |rpy-j|, |rpz-k|)  (in each axis)
 *  If any axis gives dist < EWALD_FACE_SAFE_RADIUS → flag the node.
 *
 *****************************************************************************/

 /*CHANGE INIT - PoissonVerification_FlagNodes - Diagnose nodes near particle faces */

static int ewald_flag_near_particle_nodes(ewald_charge_t* ewald,
                                           int nlocal[3], int noffset[3]) {
  int nnodes = nlocal[X] * nlocal[Y] * nlocal[Z];
  memset(ewald->node_flags, 0, nnodes * sizeof(int));

  if (ewald->cinfo == NULL) return 0;

  colloid_t* pc;
  colloids_info_local_head(ewald->cinfo, &pc);
  for (; pc; pc = pc->nextlocal) {
    /* Particle position in local coords */
    double rpx = pc->s.r[X] - noffset[X];
    double rpy = pc->s.r[Y] - noffset[Y];
    double rpz = pc->s.r[Z] - noffset[Z];

    /* Check only neighboring nodes (particle is within 1.5 of any node it can affect) */
    int ic = (int)round(rpx);
    int jc = (int)round(rpy);
    int kc = (int)round(rpz);

    for (int di = -1; di <= 1; di++) {
      for (int dj = -1; dj <= 1; dj++) {
        for (int dk = -1; dk <= 1; dk++) {
          int i = ic + di, j = jc + dj, k = kc + dk;
          if (i < 0 || i >= nlocal[X]) continue;
          if (j < 0 || j >= nlocal[Y]) continue;
          if (k < 0 || k >= nlocal[Z]) continue;

          /* Distance from particle to each face of cube centered at (i,j,k)
           * Positive value: particle inside the cube on that axis
           * Negative value: particle outside */
          double dx_face = 0.5 - fabs(rpx - i);
          double dy_face = 0.5 - fabs(rpy - j);
          double dz_face = 0.5 - fabs(rpz - k);
          double min_dist = fmin(fmin(dx_face, dy_face), dz_face);

          /* Particle inside voxel (min_dist > 0) and close to a face */
          if (min_dist >= 0.0 && min_dist < EWALD_FACE_SAFE_RADIUS) {
            int idx = i * nlocal[Y] * nlocal[Z] + j * nlocal[Z] + k;
            ewald->node_flags[idx] = EWALD_NODE_NEAR_PARTICLE;
          }
        }      
}
    }
  }
  return 0;
}

/*CHANGE END - PoissonVerification_FlagNodes */

/*****************************************************************************
 *
 *  ewald_efield_at_r
 *
 *  Compute the electric field E at an arbitrary continuous position r[3]
 *  using the Ewald sum (real-space + Fourier-space contributions).
 *
 *  Requires that sinx_[], cosx_[], sinkr_[], coskr_[] have already been
 *  computed by ewald_charge_sum_sin_cos_terms() in the main Ewald call.
 *
 *  Real-space: direct Coulomb sum over lattice nodes and particles within rc.
 *  Fourier-space: sum over k-vectors using precomputed S(k) = cosx_+i*sinx_.
 *
 *****************************************************************************/

 /*CHANGE INIT - PoissonVerification_EfieldAtPoint - E field at arbitrary continuous position */

static int ewald_efield_at_r(ewald_charge_t* ewald, const double r[3], double E[3]) {

  int nlocal[3], noffset[3];
  double ltot[3];
  int irc;
  PI_DOUBLE(pi);

  cs_nlocal(ewald->cs, nlocal);
  cs_nlocal_offset(ewald->cs, noffset);
  cs_ltot(ewald->cs, ltot);

  double fkx = 2.0 * pi / ltot[X];
  double fky = 2.0 * pi / ltot[Y];
  double fkz = 2.0 * pi / ltot[Z];
  double b0 = 1.0 / (ltot[X] * ltot[Y] * ltot[Z] * epsilon_);
  double r4alpha_sq = 1.0 / (4.0 * alpha_ * alpha_);

  E[X] = 0.0; E[Y] = 0.0; E[Z] = 0.0;

  /* -----------------------------------------------------------------------
   * [1] Real-space contribution from lattice nodes within rc
   * CHANGE INIT - RealSpace_ExactCenter - Use floor(r) as loop origin instead of round(r).
   * round(r) was correct only when r coincides with a node center (integer coords).
   * For face-sampling points r lies on a cell face (non-integer), so round(r) shifts
   * the search box and misses nodes within rc on the opposite side.
   * floor(r) + irc+1 guarantees all integer nodes within rc from the exact point r are visited. */
  irc = (int)ceil(ewald_rc_) + 1;

  if ((ewald->sources & EWALD_SOURCE_LATTICE) && ewald->psi) {
    int rx0 = (int)floor(r[X]);
    int ry0 = (int)floor(r[Y]);
    int rz0 = (int)floor(r[Z]);

    for (int dx = -irc; dx <= irc; dx++) {
      for (int dy = -irc; dy <= irc; dy++) {
        for (int dz = -irc; dz <= irc; dz++) {
          int ix_g = rx0 + dx;   /* global lattice index */
          int iy_g = ry0 + dy;
          int iz_g = rz0 + dz;

          /* Minimum image displacement */
          double drx = r[X] - ix_g;
          double dry = r[Y] - iy_g;
          double drz = r[Z] - iz_g;

          /* Apply minimum image for periodic BC */
          if (drx > 0.5 * ltot[X]) drx -= ltot[X];
          if (drx < -0.5 * ltot[X]) drx += ltot[X];
          if (dry > 0.5 * ltot[Y]) dry -= ltot[Y];
          if (dry < -0.5 * ltot[Y]) dry += ltot[Y];
          if (drz > 0.5 * ltot[Z]) drz -= ltot[Z];
          if (drz < -0.5 * ltot[Z]) drz += ltot[Z];

          double r2 = drx * drx + dry * dry + drz * drz;
          if (r2 <= 0.0 || r2 > ewald_rc_ * ewald_rc_) continue;

          double dist = sqrt(r2);

          /* Map global index to local (cs_index uses 1-based with halo) */
          int ix_l = ix_g - noffset[X];
          int iy_l = iy_g - noffset[Y];
          int iz_l = iz_g - noffset[Z];
          if (ix_l < 0 || ix_l >= nlocal[X]) continue;
          if (iy_l < 0 || iy_l >= nlocal[Y]) continue;
          if (iz_l < 0 || iz_l >= nlocal[Z]) continue;

          /* Get charge at this node (sum over ion species with sign) */
          int nk;
          psi_nk(ewald->psi, &nk);
          double rho_node = 0.0;
          int node_index = cs_index(ewald->cs, ix_l + 1, iy_l + 1, iz_l + 1);
          for (int n = 0; n < nk; n++) {
            double rho_n;
            int valency;
            psi_valency(ewald->psi, n, &valency);
            psi_rho(ewald->psi, node_index, n, &rho_n);
            /* CHANGE INIT - RhoUnits_EfieldAtR - psi_rho returns NUMBER density; multiply by eunit_ to get charge density [e/Δx³].
             * Colloid charges (pc->s.q0 etc.) are already in charge units, so this makes ions consistent. */
             // Original (wrong): rho_node += valency * rho_n;
            rho_node += eunit_ * valency * rho_n;
            /* CHANGE END - RhoUnits_EfieldAtR */
          }
          if (rho_node == 0.0) continue;

          /* Real-space Ewald field: E = q * (erfc(alpha*r)/r^2 + 2*alpha/sqrt(pi)*exp(-alpha^2*r^2)/r) * dr/r */
          double ar = alpha_ * dist;
          double erfcval = erfc(ar);
          double expval = exp(-ar * ar);
          double coeff = rho_node / (4.0 * pi * epsilon_) *
            (erfcval / (r2 * dist) + 2.0 * alpha_ * rpi_ * expval / r2);

          E[X] += coeff * drx;
          E[Y] += coeff * dry;
          E[Z] += coeff * drz;
        }
      }
    }
  }
  /* CHANGE END - RealSpace_ExactCenter */

  /* Real-space from particles */
  if ((ewald->sources & EWALD_SOURCE_COLLOIDS) && ewald->cinfo) {
    colloid_t* pc;
    colloids_info_local_head(ewald->cinfo, &pc);
    for (; pc; pc = pc->nextlocal) {
      double drx = r[X] - pc->s.r[X];
      double dry = r[Y] - pc->s.r[Y];
      double drz = r[Z] - pc->s.r[Z];

      if (drx > 0.5 * ltot[X]) drx -= ltot[X];
      if (drx < -0.5 * ltot[X]) drx += ltot[X];
      if (dry > 0.5 * ltot[Y]) dry -= ltot[Y];
      if (dry < -0.5 * ltot[Y]) dry += ltot[Y];
      if (drz > 0.5 * ltot[Z]) drz -= ltot[Z];
      if (drz < -0.5 * ltot[Z]) drz += ltot[Z];

      double r2 = drx * drx + dry * dry + drz * drz;
      if (r2 <= 0.0 || r2 > ewald_rc_ * ewald_rc_) continue;

      double dist = sqrt(r2);
      double q_p = pc->s.q0 - pc->s.q1;
      double ar = alpha_ * dist;
      double coeff = q_p / (4.0 * pi * epsilon_) *
        (erfc(ar) / (r2 * dist) + 2.0 * alpha_ * rpi_ * exp(-ar * ar) / r2);

      E[X] += coeff * drx;
      E[Y] += coeff * dry;
      E[Z] += coeff * drz;
    }
  }

  /* -----------------------------------------------------------------------
   * [2] Fourier-space contribution using precomputed S(k) = cosx_ + i*sinx_
   *
   * E_fourier(r) = (1/V) * sum_k [ -i*k * G(k) * S(k) * exp(i*k.r) ]
   *             = (1/V) * sum_k [ k * G(k) * (sinx_[kn]*cos(k.r) - cosx_[kn]*sin(k.r)) ]
   *
   * IMPORTANT: loop order and kn counting MUST match ewald_charge_sum_full_gpu:
   *   - loop order: kz (0..nk), ky (-nk..nk), kx (-nk..nk)
   *   - kn incremented only for valid k (ksq > 0 && ksq <= kmax_)
   * ----------------------------------------------------------------------- */
  {
    int kn = 0;
    for (int kz = 0; kz <= nk_[Z]; kz++) {
      for (int ky = -nk_[Y]; ky <= nk_[Y]; ky++) {
        for (int kx = -nk_[X]; kx <= nk_[X]; kx++) {

          double kvx = fkx * kx;
          double kvy = fky * ky;
          double kvz = fkz * kz;
          double k2 = kvx * kvx + kvy * kvy + kvz * kvz;

          if (k2 <= 0.0 || k2 > kmax_) continue;   /* NO kn++ here — matches storage order */

          /* Factor 2 for kz>0: the loop covers only half k-space (kz>=0).
           * Negative-kz contributes identically (S(-k)=S*(k) → same real part),
           * so we double all kz>0 terms to recover the full sum.
           * This MUST match the `factor = (kz > 0) ? 2.0 : 1.0` in ewald_charge_sum_full_gpu. */
          double factor = (kz > 0) ? 2.0 : 1.0;
          double Gk = factor * b0 * exp(-r4alpha_sq * k2) / k2;
          double kr = kvx * r[X] + kvy * r[Y] + kvz * r[Z];
          double cos_kr = cos(kr);
          double sin_kr = sin(kr);

          /* E = -grad(phi_F): phi_F = Gk*(Sk_cos*cos + Sk_sin*sin)
           * -grad = Gk*k*(Sk_cos*sin - Sk_sin*cos) = -Gk*k*impart
           * CHANGE INIT - EfieldAtR_Sign - correct sign of Fourier E field */
           // Original (wrong): E[X] += Gk * kvx * impart;
          double impart = sinx_[kn] * cos_kr - cosx_[kn] * sin_kr;
          E[X] -= Gk * kvx * impart;
          E[Y] -= Gk * kvy * impart;
          E[Z] -= Gk * kvz * impart;
          /* CHANGE END - EfieldAtR_Sign */

          kn++;
        }
      }
    }
  }

  return 0;
}

/*CHANGE END - PoissonVerification_EfieldAtPoint */

/*****************************************************************************
 *
 *  ewald_face_sample_cpu
 *
 *  For each lattice node, compute E at EWALD_NSAMPLE_FACE² quadrature points
 *  on each of the 6 faces of the unit cube centered at the node.
 *
 *  E_face layout per node (nsample_per_node = 6*N*N entries, each 3 doubles):
 *    face 0 (-X): samples [0     .. N²-1 ]
 *    face 1 (+X): samples [N²    .. 2N²-1]
 *    face 2 (-Y): samples [2N²   .. 3N²-1]
 *    face 3 (+Y): samples [3N²   .. 4N²-1]
 *    face 4 (-Z): samples [4N²   .. 5N²-1]
 *    face 5 (+Z): samples [5N²   .. 6N²-1]
 *
 *  Quadrature points on each face use a uniform (sy+0.5)/N grid in the two
 *  tangential directions, centered on the face.  This is the midpoint rule
 *  with equal weights 1/N² per point (face area = 1 in LB units).
 *
 *  Called AFTER ewald_charge_sum_sin_cos_terms() so sinx_/cosx_ are valid.
 *
 *****************************************************************************/

 /*CHANGE INIT - PoissonVerification_FaceSampling - Sample E field on cube faces (CPU) */

static int ewald_face_sample_cpu(ewald_charge_t* ewald,
                                  int nlocal[3], int noffset[3]) {
  int ns = EWALD_NSAMPLE_FACE;
  int npts = ns * ns;               /* points per face */
  int nsample_per_node = 6 * npts;
  int ntotal = nlocal[X] * nlocal[Y] * nlocal[Z];

  /* CHANGE INIT - FaceSample_Window - Restrict sampling to ±WINDOW nodes around domain centre */
  int i0, i1, j0, j1, k0, k1;
  EWALD_WINDOW_BOUNDS(nlocal[X], i0, i1);
  EWALD_WINDOW_BOUNDS(nlocal[Y], j0, j1);
  EWALD_WINDOW_BOUNDS(nlocal[Z], k0, k1);
  int nwindow = (i1 - i0 + 1) * (j1 - j0 + 1) * (k1 - k0 + 1);
  pe_info(ewald->pe, "  [Poisson] Face sampling: window [%d-%d, %d-%d, %d-%d], %d nodes (of %d)\n",
          i0, i1, j0, j1, k0, k1, nwindow, ntotal);
  /* CHANGE END - FaceSample_Window */

  for (int i = i0; i <= i1; i++) {
    for (int j = j0; j <= j1; j++) {
      for (int k = k0; k <= k1; k++) {

        int flat_idx = i * nlocal[Y] * nlocal[Z] + j * nlocal[Z] + k;

        /* Global coords of node center */
        double x0 = i + noffset[X];
        double y0 = j + noffset[Y];
        double z0 = k + noffset[Z];

        int base = flat_idx * nsample_per_node * 3;

        /* Face 0 — -X: x = x0-0.5, normal=(-1,0,0) */
        for (int sy = 0; sy < ns; sy++) {
          for (int sz = 0; sz < ns; sz++) {
            int s = sy * ns + sz;
            double r[3];
            r[X] = x0 - 0.5;
            r[Y] = y0 - 0.5 + (sy + 0.5) / ns;
            r[Z] = z0 - 0.5 + (sz + 0.5) / ns;
            double E[3];
            ewald_efield_at_r(ewald, r, E);
            ewald->E_face[base + (0 * npts + s) * 3 + X] = E[X];
            ewald->E_face[base + (0 * npts + s) * 3 + Y] = E[Y];
            ewald->E_face[base + (0 * npts + s) * 3 + Z] = E[Z];
          }
        }

        /* Face 1 — +X: x = x0+0.5, normal=(+1,0,0) */
        for (int sy = 0; sy < ns; sy++) {
          for (int sz = 0; sz < ns; sz++) {
            int s = sy * ns + sz;
            double r[3];
            r[X] = x0 + 0.5;
            r[Y] = y0 - 0.5 + (sy + 0.5) / ns;
            r[Z] = z0 - 0.5 + (sz + 0.5) / ns;
            double E[3];
            ewald_efield_at_r(ewald, r, E);
            ewald->E_face[base + (1 * npts + s) * 3 + X] = E[X];
            ewald->E_face[base + (1 * npts + s) * 3 + Y] = E[Y];
            ewald->E_face[base + (1 * npts + s) * 3 + Z] = E[Z];
          }
        }

        /* Face 2 — -Y: y = y0-0.5, normal=(0,-1,0) */
        for (int sx = 0; sx < ns; sx++) {
          for (int sz = 0; sz < ns; sz++) {
            int s = sx * ns + sz;
            double r[3];
            r[X] = x0 - 0.5 + (sx + 0.5) / ns;
            r[Y] = y0 - 0.5;
            r[Z] = z0 - 0.5 + (sz + 0.5) / ns;
            double E[3];
            ewald_efield_at_r(ewald, r, E);
            ewald->E_face[base + (2 * npts + s) * 3 + X] = E[X];
            ewald->E_face[base + (2 * npts + s) * 3 + Y] = E[Y];
            ewald->E_face[base + (2 * npts + s) * 3 + Z] = E[Z];
          }
        }

        /* Face 3 — +Y: y = y0+0.5, normal=(0,+1,0) */
        for (int sx = 0; sx < ns; sx++) {
          for (int sz = 0; sz < ns; sz++) {
            int s = sx * ns + sz;
            double r[3];
            r[X] = x0 - 0.5 + (sx + 0.5) / ns;
            r[Y] = y0 + 0.5;
            r[Z] = z0 - 0.5 + (sz + 0.5) / ns;
            double E[3];
            ewald_efield_at_r(ewald, r, E);
            ewald->E_face[base + (3 * npts + s) * 3 + X] = E[X];
            ewald->E_face[base + (3 * npts + s) * 3 + Y] = E[Y];
            ewald->E_face[base + (3 * npts + s) * 3 + Z] = E[Z];
          }
        }

        /* Face 4 — -Z: z = z0-0.5, normal=(0,0,-1) */
        for (int sx = 0; sx < ns; sx++) {
          for (int sy = 0; sy < ns; sy++) {
            int s = sx * ns + sy;
            double r[3];
            r[X] = x0 - 0.5 + (sx + 0.5) / ns;
            r[Y] = y0 - 0.5 + (sy + 0.5) / ns;
            r[Z] = z0 - 0.5;
            double E[3];
            ewald_efield_at_r(ewald, r, E);
            ewald->E_face[base + (4 * npts + s) * 3 + X] = E[X];
            ewald->E_face[base + (4 * npts + s) * 3 + Y] = E[Y];
            ewald->E_face[base + (4 * npts + s) * 3 + Z] = E[Z];
          }        
}

        /* Face 5 — +Z: z = z0+0.5, normal=(0,0,+1) */
        for (int sx = 0; sx < ns; sx++) {
          for (int sy = 0; sy < ns; sy++) {
            int s = sx * ns + sy;
            double r[3];
            r[X] = x0 - 0.5 + (sx + 0.5) / ns;
            r[Y] = y0 - 0.5 + (sy + 0.5) / ns;
            r[Z] = z0 + 0.5;
            double E[3];
            ewald_efield_at_r(ewald, r, E);
            ewald->E_face[base + (5 * npts + s) * 3 + X] = E[X];
            ewald->E_face[base + (5 * npts + s) * 3 + Y] = E[Y];
            ewald->E_face[base + (5 * npts + s) * 3 + Z] = E[Z];
          }
        }

      }
    }
  }

  pe_info(ewald->pe, "  [Poisson] Face sampling done.\n");

  return 0;
}

/*CHANGE END - PoissonVerification_FaceSampling */

/*****************************************************************************
 *
 *  ewald_calc_div_and_force
 *
 *  From the face-sampled E field, compute at each node:
 *    div_E[idx]           = numerical div(E) via Gauss theorem
 *    poisson_error[idx]   = |div(E) - rho_eff/epsilon|
 *    force_divstress[idx] = integral of Maxwell stress tensor T over faces
 *                           = force on all matter inside the voxel
 *
 *  Face normals (outward):
 *    face 0 (-X): n = (-1,0,0)   face 1 (+X): n = (+1,0,0)
 *    face 2 (-Y): n = (0,-1,0)   face 3 (+Y): n = (0,+1,0)
 *    face 4 (-Z): n = (0,0,-1)   face 5 (+Z): n = (0,0,+1)
 *
 *  face_area_elt = Δx² / N² = 1/N² (Δx=1 in LB units)
 *
 *****************************************************************************/

 /*CHANGE INIT - PoissonVerification_CalcDivAndForce - div(E), Poisson error, Maxwell stress force */

static int ewald_calc_div_and_force(ewald_charge_t* ewald,
                                     int nlocal[3], int noffset[3]) {
  int ns = EWALD_NSAMPLE_FACE;
  int npts = ns * ns;
  int nsample_per_node = 6 * npts;
  /* face_area_elt: each sample point represents area 1/N² (Δx=1) */
  double face_area_elt = 1.0 / (double)npts;

  pe_info(ewald->pe, "  [Poisson] Computing div(E) and Maxwell stress force...\n");

  /* CHANGE INIT - CalcDiv_Window - Same window as face_sample so we only read valid E_face data */
  int i0, i1, j0, j1, k0, k1;
  EWALD_WINDOW_BOUNDS(nlocal[X], i0, i1);
  EWALD_WINDOW_BOUNDS(nlocal[Y], j0, j1);
  EWALD_WINDOW_BOUNDS(nlocal[Z], k0, k1);
  /* CHANGE END - CalcDiv_Window */

  for (int i = i0; i <= i1; i++) {
    for (int j = j0; j <= j1; j++) {
      for (int k = k0; k <= k1; k++) {

        int idx = i * nlocal[Y] * nlocal[Z] + j * nlocal[Z] + k;
        int base = idx * nsample_per_node * 3;

        /* --- Poisson: numerical div(E) via Gauss theorem ---
         * Phi_face = sum_s (E · n_face) * face_area_elt
         * Outward normals: -X→-1, +X→+1, -Y→-1, +Y→+1, -Z→-1, +Z→+1 */
        double Phi_xm = 0.0, Phi_xp = 0.0;
        double Phi_ym = 0.0, Phi_yp = 0.0;
        double Phi_zm = 0.0, Phi_zp = 0.0;

        for (int s = 0; s < npts; s++) {
          Phi_xm += (-ewald->E_face[base + (0 * npts + s) * 3 + X]) * face_area_elt;
          Phi_xp += (ewald->E_face[base + (1 * npts + s) * 3 + X]) * face_area_elt;
          Phi_ym += (-ewald->E_face[base + (2 * npts + s) * 3 + Y]) * face_area_elt;
          Phi_yp += (ewald->E_face[base + (3 * npts + s) * 3 + Y]) * face_area_elt;
          Phi_zm += (-ewald->E_face[base + (4 * npts + s) * 3 + Z]) * face_area_elt;
          Phi_zp += (ewald->E_face[base + (5 * npts + s) * 3 + Z]) * face_area_elt;
        }

        double div_e = Phi_xm + Phi_xp + Phi_ym + Phi_yp + Phi_zm + Phi_zp;
        ewald->div_E[idx] = div_e;

        /* Effective charge in this voxel: ions + particle charges inside cube */
        double rho_eff = 0.0;
        if (ewald->psi) {
          int nk_species;
          psi_nk(ewald->psi, &nk_species);
          int node_index = cs_index(ewald->cs, i + 1, j + 1, k + 1);
          for (int n = 0; n < nk_species; n++) {
            int valency;
            double rho_n;
            psi_valency(ewald->psi, n, &valency);
            psi_rho(ewald->psi, node_index, n, &rho_n);
            /* CHANGE INIT - RhoUnits_CalcDiv - psi_rho returns NUMBER density; multiply by eunit_ to get charge density.
             * rho_eff must be in charge units [e/Δx³] so that rho_scaled = beta*eunit*rho_eff/eps
             * matches div(E*) in LB units where E* = beta*eunit*E_phys. */
             // Original (wrong): rho_eff += valency * rho_n;
            rho_eff += eunit_ * valency * rho_n;
            /* CHANGE END - RhoUnits_CalcDiv */
          }
        }
        if (ewald->cinfo) {
          colloid_t* pc;
          colloids_info_local_head(ewald->cinfo, &pc);
          for (; pc; pc = pc->nextlocal) {
            double dx = pc->s.r[X] - (i + noffset[X]);
            double dy = pc->s.r[Y] - (j + noffset[Y]);
            double dz = pc->s.r[Z] - (k + noffset[Z]);
            /* CHANGE INIT - ColoidAssignment - Distribute charge over all voxels containing the colloid
             * A colloid exactly on a face/edge/corner is shared by 2/4/8 voxels equally.
             * Count how many voxels claim this particle, then divide q by that count. */
             // Original: if (fabs(dx) < 0.5 && fabs(dy) < 0.5 && fabs(dz) < 0.5)
            if (fabs(dx) <= 0.5 && fabs(dy) <= 0.5 && fabs(dz) <= 0.5) {
              int share = ((fabs(dx) == 0.5) ? 2 : 1)
                * ((fabs(dy) == 0.5) ? 2 : 1)
                * ((fabs(dz) == 0.5) ? 2 : 1);
              rho_eff += (pc->s.q0 - pc->s.q1) / share;  /* q already in charge units [e] */
            }
            /* CHANGE END - ColoidAssignment */
          }
        }

        /* div_e here is from face-sampling of ewald_efield_at_r → physical units [E_phys].
         * Gauss's law (physical): div(E_phys) = rho_charge/eps.
         * rho_eff is in charge units [e/Δx³], so compare with rho_eff/epsilon_.
         * CHANGE INIT - DivEUnits_CalcDiv - corrected unit comparison for face-sampled div_E */
         // Original (wrong, mixed physical vs LB units): ewald->poisson_error[idx] = fabs(div_e - beta_ * eunit_ * rho_eff / epsilon_);
        ewald->poisson_error[idx] = fabs(div_e - rho_eff / epsilon_);
        /* CHANGE END - DivEUnits_CalcDiv */

    /* CHANGE INIT - PoissonOnly - Skip Maxwell stress when EWALD_POISSON_ONLY is defined */
#ifndef EWALD_POISSON_ONLY
    /* --- Maxwell stress tensor force ---
     * T_ij = ε*(E_i*E_j - δ_ij*E²/2)
     * F_i = ∮ T_ij * n_j dS  (always computed, also for near-particle nodes)
     * If particle inside voxel: result includes force on fluid + particle */
        double Fx = 0.0, Fy = 0.0, Fz = 0.0;
        double Ex, Ey, Ez, Esq;

        for (int s = 0; s < npts; s++) {

          /* Face +X: n=(+1,0,0) → F_i += T_i1 * dA */
          Ex = ewald->E_face[base + (1 * npts + s) * 3 + X];
          Ey = ewald->E_face[base + (1 * npts + s) * 3 + Y];
          Ez = ewald->E_face[base + (1 * npts + s) * 3 + Z];
          Esq = Ex * Ex + Ey * Ey + Ez * Ez;
          Fx += epsilon_ * (Ex * Ex - 0.5 * Esq) * face_area_elt;
          Fy += epsilon_ * (Ey * Ex) * face_area_elt;
          Fz += epsilon_ * (Ez * Ex) * face_area_elt;

          /* Face -X: n=(-1,0,0) → F_i -= T_i1 * dA */
          Ex = ewald->E_face[base + (0 * npts + s) * 3 + X];
          Ey = ewald->E_face[base + (0 * npts + s) * 3 + Y];
          Ez = ewald->E_face[base + (0 * npts + s) * 3 + Z];
          Esq = Ex * Ex + Ey * Ey + Ez * Ez;
          Fx -= epsilon_ * (Ex * Ex - 0.5 * Esq) * face_area_elt;
          Fy -= epsilon_ * (Ey * Ex) * face_area_elt;
          Fz -= epsilon_ * (Ez * Ex) * face_area_elt;

          /* Face +Y: n=(0,+1,0) → F_i += T_i2 * dA */
          Ex = ewald->E_face[base + (3 * npts + s) * 3 + X];
          Ey = ewald->E_face[base + (3 * npts + s) * 3 + Y];
          Ez = ewald->E_face[base + (3 * npts + s) * 3 + Z];
          Esq = Ex * Ex + Ey * Ey + Ez * Ez;
          Fx += epsilon_ * (Ex * Ey) * face_area_elt;
          Fy += epsilon_ * (Ey * Ey - 0.5 * Esq) * face_area_elt;
          Fz += epsilon_ * (Ez * Ey) * face_area_elt;

          /* Face -Y: n=(0,-1,0) → F_i -= T_i2 * dA */
          Ex = ewald->E_face[base + (2 * npts + s) * 3 + X];
          Ey = ewald->E_face[base + (2 * npts + s) * 3 + Y];
          Ez = ewald->E_face[base + (2 * npts + s) * 3 + Z];
          Esq = Ex * Ex + Ey * Ey + Ez * Ez;
          Fx -= epsilon_ * (Ex * Ey) * face_area_elt;
          Fy -= epsilon_ * (Ey * Ey - 0.5 * Esq) * face_area_elt;
          Fz -= epsilon_ * (Ez * Ey) * face_area_elt;

          /* Face +Z: n=(0,0,+1) → F_i += T_i3 * dA */
          Ex = ewald->E_face[base + (5 * npts + s) * 3 + X];
          Ey = ewald->E_face[base + (5 * npts + s) * 3 + Y];
          Ez = ewald->E_face[base + (5 * npts + s) * 3 + Z];
          Esq = Ex * Ex + Ey * Ey + Ez * Ez;
          Fx += epsilon_ * (Ex * Ez) * face_area_elt;
          Fy += epsilon_ * (Ey * Ez) * face_area_elt;
          Fz += epsilon_ * (Ez * Ez - 0.5 * Esq) * face_area_elt;

          /* Face -Z: n=(0,0,-1) → F_i -= T_i3 * dA */
          Ex = ewald->E_face[base + (4 * npts + s) * 3 + X];
          Ey = ewald->E_face[base + (4 * npts + s) * 3 + Y];
          Ez = ewald->E_face[base + (4 * npts + s) * 3 + Z];
          Esq = Ex * Ex + Ey * Ey + Ez * Ez;
          Fx -= epsilon_ * (Ex * Ez) * face_area_elt;
          Fy -= epsilon_ * (Ey * Ez) * face_area_elt;
          Fz -= epsilon_ * (Ez * Ez - 0.5 * Esq) * face_area_elt;
        }

        ewald->force_divstress[idx * 3 + X] = Fx;
        ewald->force_divstress[idx * 3 + Y] = Fy;
        ewald->force_divstress[idx * 3 + Z] = Fz;
#endif /* EWALD_POISSON_ONLY */
        /* CHANGE END - PoissonOnly */

      }
    }
  }

  pe_info(ewald->pe, "  [Poisson] div(E) and force: done.\n");

  return 0;
}

/*CHANGE END - PoissonVerification_CalcDivAndForce */

/*****************************************************************************
 *
 *  ewald_write_poisson_verification
 *  ewald_write_force_divstress
 *
 *  Write div(E), Poisson error, and Maxwell stress force to ASCII files.
 *  Called at the same output frequency as efield files.
 *  Columns: ix iy iz  value(s)  node_flag
 *
 *****************************************************************************/

 /*CHANGE INIT - PoissonVerification_Output - Write Poisson and force data to files */

static int ewald_write_poisson_verification(ewald_charge_t* ewald,
                                             int nlocal[3], int noffset[3],
                                             int step) {
  char filename[512];
  snprintf(filename, sizeof(filename), "poisson_divergence_%06d.dat", step);
  FILE* f = fopen(filename, "w");
  if (f == NULL) {
    pe_info(ewald->pe, "  [Poisson] WARNING: could not open %s\n", filename);
    return -1;
  }

  /* CHANGE INIT - WritePoissonWindow - Only write nodes that were sampled */
  int i0, i1, j0, j1, k0, k1;
  EWALD_WINDOW_BOUNDS(nlocal[X], i0, i1);
  EWALD_WINDOW_BOUNDS(nlocal[Y], j0, j1);
  EWALD_WINDOW_BOUNDS(nlocal[Z], k0, k1);
  /* CHANGE END - WritePoissonWindow */

  fprintf(f, "ix iy iz div_E rho_over_eps err_divE rel_err_divE node_flag\n");
  for (int i = i0; i <= i1; i++) {
    for (int j = j0; j <= j1; j++) {
      for (int k = k0; k <= k1; k++) {
        int idx = i * nlocal[Y] * nlocal[Z] + j * nlocal[Z] + k;

        /* Recompute rho_eff for output (same logic as ewald_calc_div_and_force) */
        double rho_eff = 0.0;
        if (ewald->psi) {
          int nk_species;
          psi_nk(ewald->psi, &nk_species);
          int node_index = cs_index(ewald->cs, i + 1, j + 1, k + 1);
          for (int n = 0; n < nk_species; n++) {
            int valency; double rho_n;
            psi_valency(ewald->psi, n, &valency);
            psi_rho(ewald->psi, node_index, n, &rho_n);
            /* CHANGE INIT - RhoUnits_WritePoissonVer - psi_rho returns NUMBER density; multiply by eunit_ to get charge density */
            // Original (wrong): rho_eff += valency * rho_n;
            rho_eff += eunit_ * valency * rho_n;
            /* CHANGE END - RhoUnits_WritePoissonVer */
          }
        }
        if (ewald->cinfo) {
          colloid_t* pc;
          colloids_info_local_head(ewald->cinfo, &pc);
          for (; pc; pc = pc->nextlocal) {
            double dx = pc->s.r[X] - (i + noffset[X]);
            double dy = pc->s.r[Y] - (j + noffset[Y]);
            double dz = pc->s.r[Z] - (k + noffset[Z]);
            /* CHANGE INIT - ColoidAssignment - Distribute charge over all voxels containing the colloid */
            // Original: if (fabs(dx) < 0.5 && fabs(dy) < 0.5 && fabs(dz) < 0.5)
            if (fabs(dx) <= 0.5 && fabs(dy) <= 0.5 && fabs(dz) <= 0.5) {
              int share = ((fabs(dx) == 0.5) ? 2 : 1)
                * ((fabs(dy) == 0.5) ? 2 : 1)
                * ((fabs(dz) == 0.5) ? 2 : 1);
              rho_eff += (pc->s.q0 - pc->s.q1) / share;
            }
            /* CHANGE END - ColoidAssignment */
          }
        }

        /* CHANGE INIT - WritePoissonVer_RelErr - add relative error */
        double rho_over_eps = rho_eff / epsilon_;
        double abs_err = fabs(ewald->div_E[idx] - rho_over_eps);
        double rel_err = (fabs(rho_over_eps) > 0.0) ? abs_err / fabs(rho_over_eps) : 0.0;
        fprintf(f, "%d %d %d  %.12e  %.12e  %.12e  %.12e  %d\n",
                i + noffset[X], j + noffset[Y], k + noffset[Z],
                ewald->div_E[idx],
                rho_over_eps,
                abs_err,
                rel_err,
                ewald->node_flags[idx]);
        /* CHANGE END - WritePoissonVer_RelErr */
      }    
}
  }

  fclose(f);
  pe_info(ewald->pe, "  [Poisson] Written: %s\n", filename);
  return 0;
}

static int ewald_write_force_divstress(ewald_charge_t* ewald,
                                        int nlocal[3], int noffset[3],
                                        int step) {
  char filename[512];
  snprintf(filename, sizeof(filename), "force_divstress_%06d.dat", step);
  FILE* f = fopen(filename, "w");
  if (f == NULL) {
    pe_info(ewald->pe, "  [Force] WARNING: could not open %s\n", filename);
    return -1;
  }

  /* CHANGE INIT - WriteForceWindow - Only write nodes that were sampled */
  int i0f, i1f, j0f, j1f, k0f, k1f;
  EWALD_WINDOW_BOUNDS(nlocal[X], i0f, i1f);
  EWALD_WINDOW_BOUNDS(nlocal[Y], j0f, j1f);
  EWALD_WINDOW_BOUNDS(nlocal[Z], k0f, k1f);
  /* CHANGE END - WriteForceWindow */

  fprintf(f, "# ix iy iz  Fx  Fy  Fz  node_flag\n");
  for (int i = i0f; i <= i1f; i++) {
    for (int j = j0f; j <= j1f; j++) {
      for (int k = k0f; k <= k1f; k++) {
        int idx = i * nlocal[Y] * nlocal[Z] + j * nlocal[Z] + k;
        fprintf(f, "%d %d %d  %.12e  %.12e  %.12e  %d\n",
                i + noffset[X], j + noffset[Y], k + noffset[Z],
                ewald->force_divstress[idx * 3 + X],
                ewald->force_divstress[idx * 3 + Y],
                ewald->force_divstress[idx * 3 + Z],
                ewald->node_flags[idx]);
      }    
}
  }

  fclose(f);
  pe_info(ewald->pe, "  [Force]   Written: %s\n", filename);
  return 0;
}

/* CHANGE INIT - PoissonStencil - Poisson check using discrete stencil
 *
 *  Two independent verifications written to poisson_stencil_NNNNNN.dat:
 *
 *  (A) div(E) via wgradients applied to psi->efield:
 *        div_E = sum_p  wgradients[p] * (cx*Ex + cy*Ey + cz*Ez)[p]
 *      E is stored in LB units: E* = beta*eunit*E_phys.
 *      Poisson requires: div_E* = beta*eunit*rho/eps = rho_scaled.
 *
 *  (B) Laplacian of ψ* via wlaplacian applied to psi->psi:
 *        lap_psi* = sum_p  wlaplacian[p] * psi*[p]
 *      psi* = beta*eunit*phi_phys, so Poisson becomes:
 *        lap_psi* = -beta*eunit*rho/eps = -rho_scaled
 *      Equivalent to (A) if div·grad = lap (true for D3Q7, approximate for D3Q27).
 *
 *  Columns: ix iy iz  div_E  lap_psi  rho_scaled  err_divE  err_lapPsi  node_flag
 *    rho_scaled = beta*eunit*rho_eff/epsilon
 *    err_divE   = |div_E  - rho_scaled|
 *    err_lapPsi = |lap_psi + rho_scaled|   (lap=-rho/eps → lap+rho/eps=0)
 */

static int ewald_write_poisson_stencil(ewald_charge_t* ewald,
                                        int nlocal[3], int noffset[3],
                                        int step) {
  if (!ewald->psi || !ewald->psi->stencil) {
    pe_info(ewald->pe, "  [PoissonStencil] psi or stencil not available, skipping.\n");
    return 0;
  }

  cs_t* cs = ewald->cs;
  stencil_t* s = ewald->psi->stencil;
  field_t* ef = ewald->psi->efield;   /* may be NULL */
  int      nsites = ewald->psi->nsites;

  char filename[512];
  snprintf(filename, sizeof(filename), "poisson_stencil_%06d.dat", step);
  FILE* f = fopen(filename, "w");
  if (f == NULL) {
    pe_info(ewald->pe, "  [PoissonStencil] WARNING: could not open %s\n", filename);
    return -1;
  }

  int i0, i1, j0, j1, k0, k1;
  EWALD_WINDOW_BOUNDS(nlocal[X], i0, i1);
  EWALD_WINDOW_BOUNDS(nlocal[Y], j0, j1);
  EWALD_WINDOW_BOUNDS(nlocal[Z], k0, k1);

  fprintf(f, "ix iy iz div_E lap_psi rho_scaled err_divE err_lapPsi rel_err_divE rel_err_lapPsi node_flag\n");

  for (int i = i0; i <= i1; i++) {
    for (int j = j0; j <= j1; j++) {
      for (int k = k0; k <= k1; k++) {

        int ic = i + 1, jc = j + 1, kc = k + 1;   /* 1-based for cs_index */
        int index0 = cs_index(cs, ic, jc, kc);
        int flat_idx = i * nlocal[Y] * nlocal[Z] + j * nlocal[Z] + k;

        /* (A) div(E) using wgradients on stored efield.
         * psi->efield stores E_phys (physical units, no beta*eunit factor).
         * Multiply by beta_*eunit_ to convert to LB units, matching lap_psi and rho_scaled.
         * CHANGE INIT - DivE_LBunits_Stencil - scale div_E to LB units (beta*eunit*div(E_phys)) */
        double div_e = 0.0;
        if (ef) {
          for (int p = 1; p < s->npoints; p++) {
            int8_t cx = s->cv[p][X];
            int8_t cy = s->cv[p][Y];
            int8_t cz = s->cv[p][Z];
            int index1 = cs_index(cs, ic + cx, jc + cy, kc + cz);
            double E[3];
            field_vector(ef, index1, E);
            div_e += s->wgradients[p] * (cx * E[X] + cy * E[Y] + cz * E[Z]);
          }
          div_e *= beta_ * eunit_;   /* convert E_phys → E* = beta*eunit*E_phys */
        }
        /* CHANGE END - DivE_LBunits_Stencil */

        /* (B) Laplacian de ψ*.
         * El stencil D3Q27 tiene wlaplacian[p] = -6*wv[p] < 0 para vecinos →
         * Σ_p wlaplacian[p]*f[p] = -∇²f.  Se niega para obtener ∇²ψ*.
         * Poisson: ∇²ψ* = -rho_scaled  →  lap_psi debe tener signo opuesto a rho_scaled.
         * Residuo: |lap_psi + rho_scaled| = 0 exacto. */
        double lap_psi = 0.0;
        for (int p = 0; p < s->npoints; p++) {
          int8_t cx = s->cv[p][X];
          int8_t cy = s->cv[p][Y];
          int8_t cz = s->cv[p][Z];
          int index1 = cs_index(cs, ic + cx, jc + cy, kc + cz);
          double psi_val = ewald->psi->psi->data[addr_rank0(nsites, index1)];
          lap_psi += s->wlaplacian[p] * psi_val;
        }
        lap_psi = -lap_psi;   /* stencil da -∇²ψ; negar para obtener ∇²ψ* */

        /* rho_eff (same logic as ewald_write_poisson_verification) */
        double rho_eff = 0.0;
        if (ewald->psi) {
          int nk_species;
          psi_nk(ewald->psi, &nk_species);
          for (int n = 0; n < nk_species; n++) {
            int valency; double rho_n;
            psi_valency(ewald->psi, n, &valency);
            psi_rho(ewald->psi, index0, n, &rho_n);
            /* CHANGE INIT - RhoUnits_PoissonStencil - psi_rho returns NUMBER density; multiply by eunit_ to get charge density.
             * rho_scaled = beta*eunit*rho_eff/eps must match lap(psi*) and div(E*) in LB units. */
             // Original (wrong): rho_eff += valency * rho_n;
            rho_eff += eunit_ * valency * rho_n;
            /* CHANGE END - RhoUnits_PoissonStencil */
          }
        }
        if (ewald->cinfo) {
          colloid_t* pc;
          colloids_info_local_head(ewald->cinfo, &pc);
          for (; pc; pc = pc->nextlocal) {
            double dx = pc->s.r[X] - (i + noffset[X]);
            double dy = pc->s.r[Y] - (j + noffset[Y]);
            double dz = pc->s.r[Z] - (k + noffset[Z]);
            if (fabs(dx) <= 0.5 && fabs(dy) <= 0.5 && fabs(dz) <= 0.5) {
              int share = ((fabs(dx) == 0.5) ? 2 : 1)
                * ((fabs(dy) == 0.5) ? 2 : 1)
                * ((fabs(dz) == 0.5) ? 2 : 1);
              rho_eff += (pc->s.q0 - pc->s.q1) / share;  /* colloid q already in charge units [e] */
            }
          }
        }

        /* rho_scaled = beta*eunit * rho_eff/epsilon (unidades LB).
         * lap_psi = -∇²ψ* = +rho_scaled  →  residuo: |lap_psi - rho_scaled| = 0 exacto.
         * div_e   = +∇·E* = +rho_scaled   →  residuo: |div_e   - rho_scaled| = 0 exacto. */
         /* CHANGE INIT - PoissonStencil_RelErr - add relative errors */
        double rho_scaled = beta_ * eunit_ * rho_eff / epsilon_;
        double abs_ref = fabs(rho_scaled);
        double err_divE = fabs(div_e - rho_scaled);   /* Gauss: ∇·E* = +rho_scaled */
        double err_lapPsi = fabs(lap_psi + rho_scaled);   /* Poisson: ∇²ψ* = -rho_scaled */
        double rel_divE = (abs_ref > 0.0) ? err_divE / abs_ref : 0.0;
        double rel_lapPsi = (abs_ref > 0.0) ? err_lapPsi / abs_ref : 0.0;
        fprintf(f, "%d %d %d  %.12e  %.12e  %.12e  %.12e  %.12e  %.12e  %.12e  %d\n",
                i + noffset[X], j + noffset[Y], k + noffset[Z],
                div_e,
                lap_psi,
                rho_scaled,
                err_divE, err_lapPsi,
                rel_divE, rel_lapPsi,
                ewald->node_flags[flat_idx]);
        /* CHANGE END - PoissonStencil_RelErr */
      }
    }
  }

  fclose(f);
  pe_info(ewald->pe, "  [PoissonStencil] Written: %s\n", filename);
  return 0;
}

// static int ewald_write_force_fluid(ewald_charge_t* ewald,
//                                    int nlocal[3], int noffset[3],
//                                    int step) {
//   if (!ewald->psi || !ewald->psi->force_fluid) return 0;

//   cs_t* cs = ewald->cs;

//   char filename[512];
//   snprintf(filename, sizeof(filename), "force_fluid_%06d.dat", step);
//   FILE* f = fopen(filename, "w");
//   if (f == NULL) {
//     pe_info(ewald->pe, "  [ForceFluid] WARNING: could not open %s\n", filename);
//     return -1;
//   }

//   int i0, i1, j0, j1, k0, k1;
//   EWALD_WINDOW_BOUNDS(nlocal[X], i0, i1);
//   EWALD_WINDOW_BOUNDS(nlocal[Y], j0, j1);
//   EWALD_WINDOW_BOUNDS(nlocal[Z], k0, k1);

//   fprintf(f, "ix iy iz  Fx  Fy  Fz\n");
//   for (int i = i0; i <= i1; i++) {
//     for (int j = j0; j <= j1; j++) {
//       for (int k = k0; k <= k1; k++) {
//         int index = cs_index(cs, i + 1, j + 1, k + 1);
//         double F[3];
//         field_vector(ewald->psi->force_fluid, index, F);
//         fprintf(f, "%d %d %d  %.12e  %.12e  %.12e\n",
//                 i + noffset[X], j + noffset[Y], k + noffset[Z],
//                 F[X], F[Y], F[Z]);
//       }
//     }
//   }

//   fclose(f);
//   pe_info(ewald->pe, "  [ForceFluid] Written: %s\n", filename);
//   return 0;
// }