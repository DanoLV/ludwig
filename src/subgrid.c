/*****************************************************************************
 *
 *  subgrid.c
 *
 *  Routines for point-like particles.
 *
 *  See Nash et al. (2007).
 *
 *  Overview.
 *
 *  Two-way coupling between sub-grid particles and the fluid is implemented
 *  in roughly two phases.
 *
 *  (1) Force on each particle from the immediately surrounding fluid.
 *  (2) Influence of the particles on local fluid nodes;
 *
 *  (1) subgrid_update() is responsible for the particle update and setting
 *      any position increment dr. Schematically:
 *
 *      -> subgrid_interpolation()
 *         accumulates contributions to the force on the particle from
 *         local fluid nodes to local particle copies (fsub[3]);
 *      -> COLLOID_SUM_SUBGRID
 *         ensures all copies agree on the net force per particle fsub[3].
 *      -> all copies update v and dr = v.dt from fsub[3] and must agree.
 *      -> Actual position updates must be deferred until the start of
 *         the next time step and solloids_info_position_update().
 *
 *  (2) subgrid_force_from_particle() is responsible for computing
 *      the force on local fluid nodes from particles. Schematically;
 *
 *      -> On entry, fex[3] may contain pair interaction and other
 *         "external" forces on the particle;
 *      -> subgrid_wall_lubrication()
 *          detect and compute particle/wall lubrication forces,
 *          and accumulate to fex[3] once per particle (i.e. local copies
 *          only involved);
 *      -> COLLOID_SUM_FORCE_EXT_ONLY
 *         => all copies agree on fext[3], text[3]
 *
 *      -> All particle copies contribute \delta(r - R_i) fext[3] to local
 *         fluid nodes only via hydro_f_local_add() at position r.
 *      -> This force may then enter the fluid collision stage.
 *
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  (c) 2010-2022 The University of Edinburgh
 *
 *  Contributing authors:
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *
 *****************************************************************************/

#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <string.h>

#include "pe.h"
#include "coords.h"
#include "physics.h"
#include "colloids.h"
#include "colloid_sums.h"
#include "util.h"
#include "subgrid.h"
 /*CHANGE INIT - Subgrid charge */
#include "psi_colloid.h"
#include "psi_gradients.h"
#include "psi_petsc.h"
#include "field.h"
#include "target.h"
/*CHANGE END - Subgrid charge */

// static double d_peskin(double);
static int subgrid_interpolation(colloids_info_t* cinfo, hydro_t* hydro);
// CHANGE INIT - 20260117 cambio de soporte de peskin
// static const double drange_ = 1.0; /* Max. range of interpolation - 1 */
static const double drange_ = 1.0; /* Max. range of interpolation - 1 */
// CHANGE END - 20260117 cambio de soporte de peskin

/*CHANGE INIT - Poisson-Boltzmann kappa (inverse Debye length) */
static double kappa_pb_ = 0.0;  /* Stored kappa value, computed once */
static int kappa_initialized_ = 0;  /* Flag to check if kappa has been computed */
/*CHANGE END - Poisson-Boltzmann kappa */

/*CHANGE INIT - 20260117 Short-range corrections infrastructure (Level 0)
 * See docs/SHORT_RANGE_CORRECTIONS_ANALYSIS.md for full documentation.
 * This implements P3M-style short-range corrections for electric interactions
 * of subgrid particles to eliminate spurious velocities and improve accuracy.
 */

 /* Short-range correction parameters */
 /* CHANGE 20260118: Reduced kappa_E from 2.0 to 0.5 to avoid over-correction.
  * With kappa_E = 2.0, the correction was ~23x larger than the mesh field.
  * With kappa_E = 0.5, the correction is ~5-6x the mesh field, which is
  * still large but more reasonable for testing. This needs calibration.
  * OLD VALUE: 2.0
  */
#define SHORTRANGE_KAPPA_E  0.5    /* Ewald splitting parameter (numerical, NOT Debye-Huckel) */
#define SHORTRANGE_RC       1.5    /* Cutoff radius in lattice units */

  /* Short-range correction level selector */
typedef enum {
	SHORTRANGE_LEVEL_NONE = 0,      /* No corrections */
	SHORTRANGE_LEVEL_1_SELF = 1,    /* Self-field only (Alternative A) */
	SHORTRANGE_LEVEL_2_PP = 2,      /* + Particle-Particle (Alternative B) */
	SHORTRANGE_LEVEL_3_PN = 3,      /* + Particle-Node (Alternative F) */
	SHORTRANGE_LEVEL_4_P3M = 4,     /* + Reference function (Alternative D) */
	SHORTRANGE_LEVEL_5_CONSERV = 5  /* + Mass conservation (Alternative E) */
} shortrange_level_t;

/* Current correction level (can be set from input) */
static shortrange_level_t shortrange_level_ = SHORTRANGE_LEVEL_NONE;

/*****************************************************************************
 *
 *  shortrange_erf_kernel
 *
 *  Computes the short-range correction kernel: erf(kappa*r) / r^2
 *  This kernel represents the part of the interaction that the mesh
 *  does NOT capture well and must be corrected directly.
 *
 *  Note: kappa_E is a NUMERICAL splitting parameter, NOT the physical
 *  Debye-Huckel screening length. See Section 1.5 of the analysis document.
 *
 *****************************************************************************/
static double shortrange_erf_kernel(double r, double kappa_E) {
	if (r < 1.0e-10) return 0.0;  /* Avoid singularity at r=0 */
	return erf(kappa_E * r) / (r * r);
}

/*****************************************************************************
 *
 *  shortrange_coulomb_prefactor
 *
 *  Computes the Coulomb field prefactor: 1 / (4*pi*epsilon)
 *  This is used for computing electric fields from charges.
 *
 *  Note: For the electric FIELD E = q/(4*pi*epsilon*r^2), we need 1/(4*pi*epsilon).
 *  The factor e^2/(4*pi*epsilon) would be for FORCE or ENERGY, not field.
 *
 *  CHANGE 20260118: Fixed prefactor - removed e^2 factor
 *  OLD CODE: return (eunit * eunit) / (4.0 * pi * epsilon);
 *
 *****************************************************************************/
static double shortrange_coulomb_prefactor(psi_t* psi) {
	double epsilon;
	PI_DOUBLE(pi);
	psi_epsilon(psi, &epsilon);
	return 1.0 / (4.0 * pi * epsilon);
}

/*****************************************************************************
 *
 *  shortrange_force
 *
 *  Computes the short-range correction force vector between two charges.
 *  F_corr = prefactor * q1 * q2 * erf(kappa*r)/r^2 * r_hat
 *
 *  Arguments:
 *    dist      - distance between charges
 *    q1, q2    - charges (in units of elementary charge)
 *    prefactor - Coulomb prefactor e^2/(4*pi*epsilon)
 *    kappa_E   - Ewald splitting parameter
 *    r_vec[3]  - separation vector (from 1 to 2)
 *    F[3]      - OUTPUT: force vector on charge 1
 *
 *****************************************************************************/
static void shortrange_force(double dist, double q1, double q2,
							  double prefactor, double kappa_E,
							  const double r_vec[3], double F[3]) {
	double kernel = shortrange_erf_kernel(dist, kappa_E);
	double mag = prefactor * q1 * q2 * kernel / dist;  /* Force magnitude / r */
	F[X] = mag * r_vec[X];
	F[Y] = mag * r_vec[Y];
	F[Z] = mag * r_vec[Z];
}

/*CHANGE END - 20260117 Short-range corrections infrastructure (Level 0) */

/*****************************************************************************
 *
 *  subgrid_force_from_particles()
 *
 *  For each particle, accumulate the force on the relevant surrounding
 *  lattice nodes. Only nodes in the local domain are involved.
 *
 *  If there are no subgrid particles, hydro is allowed to be NULL.
 *
 *****************************************************************************/
int subgrid_force_from_particles(colloids_info_t* cinfo, hydro_t* hydro,
				 wall_t* wall) {

	int ic, jc, kc;
	int i, j, k, i_min, i_max, j_min, j_max, k_min, k_max;
	int index;
	int nlocal[3], offset[3];
	int ncell[3];

	double r[3], r0[3], force[3];
	double dr;
	colloid_t* p_colloid = NULL;  /* Subgrid colloid */
	colloid_t* presolved = NULL;  /* Resolved colloid occupuing node */

	assert(cinfo);
	assert(wall);

	if (cinfo->nsubgrid == 0) return 0;

	cs_nlocal(cinfo->cs, nlocal);
	cs_nlocal_offset(cinfo->cs, offset);
	colloids_info_ncell(cinfo, ncell);

	/* Add any wall lubrication corrections before communication to
	 * find total external force on each particle */

	subgrid_wall_lubrication(cinfo, wall);
	colloid_sums_halo(cinfo, COLLOID_SUM_FORCE_EXT_ONLY);

	/* While there is no device implementation, must copy back-and forth
	 * the force. */

	assert(hydro);
	hydro_memcpy(hydro, tdpMemcpyDeviceToHost);

	/* Loop through all cells (including the halo cells) */

	for (ic = 0; ic <= ncell[X] + 1; ic++) {
		for (jc = 0; jc <= ncell[Y] + 1; jc++) {
			for (kc = 0; kc <= ncell[Z] + 1; kc++) {

				colloids_info_cell_list_head(cinfo, ic, jc, kc, &p_colloid);

				for (; p_colloid; p_colloid = p_colloid->next) {

					if (p_colloid->s.bc != COLLOID_BC_SUBGRID) continue;

					/* Need to translate the colloid position to "local"
					 * coordinates, so that the correct range of lattice
					 * nodes is found */

					r0[X] = p_colloid->s.r[X] - 1.0 * offset[X];
					r0[Y] = p_colloid->s.r[Y] - 1.0 * offset[Y];
					r0[Z] = p_colloid->s.r[Z] - 1.0 * offset[Z];

					/* Work out which local lattice sites are involved
					 * and loop around */

					i_min = imax(1, (int)floor(r0[X] - drange_));
					i_max = imin(nlocal[X], (int)ceil(r0[X] + drange_));
					j_min = imax(1, (int)floor(r0[Y] - drange_));
					j_max = imin(nlocal[Y], (int)ceil(r0[Y] + drange_));
					k_min = imax(1, (int)floor(r0[Z] - drange_));
					k_max = imin(nlocal[Z], (int)ceil(r0[Z] + drange_));

					for (i = i_min; i <= i_max; i++) {
						for (j = j_min; j <= j_max; j++) {
							for (k = k_min; k <= k_max; k++) {

								index = cs_index(cinfo->cs, i, j, k);

								/* Separation between r0 and the coordinate position of
								 * this site */

								r[X] = r0[X] - 1.0 * i;
								r[Y] = r0[Y] - 1.0 * j;
								r[Z] = r0[Z] - 1.0 * k;

								dr = d_peskin(r[X]) * d_peskin(r[Y]) * d_peskin(r[Z]);

								force[X] = p_colloid->fex[X] * dr;
								force[Y] = p_colloid->fex[Y] * dr;
								force[Z] = p_colloid->fex[Z] * dr;

								colloids_info_map(cinfo, index, &presolved);

								if (presolved == NULL) {
									hydro_f_local_add(hydro, index, force);
								}
								else {
									double rd[3] = { 0 };
									double torque[3] = { 0 };
									presolved->force[X] += force[X];
									presolved->force[Y] += force[Y];
									presolved->force[Z] += force[Z];
									rd[X] = 1.0 * i - (presolved->s.r[X] - 1.0 * offset[X]);
									rd[Y] = 1.0 * j - (presolved->s.r[Y] - 1.0 * offset[Y]);
									rd[Z] = 1.0 * k - (presolved->s.r[Z] - 1.0 * offset[Z]);
									cross_product(rd, force, torque);
									presolved->torque[X] += torque[X];
									presolved->torque[Y] += torque[Y];
									presolved->torque[Z] += torque[Z];
								}

							}
						}
					}

					/* Next colloid */
				}

				/* Next cell */
			}
		}
	}

	hydro_memcpy(hydro, tdpMemcpyHostToDevice);

	return 0;
}
/*CHANGE INIT - Subgrid charge */
/*****************************************************************************
 *
 *  binary_search_charge_index()
 *
 *  Binary search to find position to insert/find cs_index in sorted array
 *
 *****************************************************************************/
static int binary_search_charge_index(distributed_charge_klein_t* charge, int cs_index) {
	int left = 0;
	int right = charge->count - 1;

	while (left <= right) {
		int mid = left + (right - left) / 2;
		if (charge->entries[mid]->cs_index == cs_index) {
			return mid;  /* Found */
		}
		if (charge->entries[mid]->cs_index < cs_index) {
			left = mid + 1;
		}
		else {
			right = mid - 1;
		}
	}
	return left;  /* Position to insert */
}

/*****************************************************************************
 *
 *  add_charge_to_array()
 *
 *  Add or accumulate charge at given cs_index
 *
 *****************************************************************************/
void add_charge_to_array(distributed_charge_klein_t** charge_ptr, int cs_index,
								double q0_dr, double q1_dr, psi_t* obj) {
	distributed_charge_klein_t* charge = *charge_ptr;

	/* Initialize if needed */
	if (charge == NULL) {
		charge = (distributed_charge_klein_t*)malloc(sizeof(distributed_charge_klein_t));
		charge->entries = (distributed_charge_klein_entry_t**)malloc(16 * sizeof(distributed_charge_klein_entry_t*));
		charge->count = 0;
		charge->capacity = 16;
		*charge_ptr = charge;
	}

	/* Binary search for index */
	int pos = binary_search_charge_index(charge, cs_index);

	/* Check if index already exists */
	if (pos < charge->count&& charge->entries[pos]->cs_index == cs_index) {
		/* Accumulate to existing entry */
		klein_add_double(charge->entries[pos]->rho0_sum, q0_dr);
		klein_add_double(charge->entries[pos]->rho1_sum, q1_dr);
	}
	else {
		/* Need to insert new entry */
		if (charge->count >= charge->capacity) {
			charge->capacity *= 2;
			charge->entries = (distributed_charge_klein_entry_t**)realloc(charge->entries,
							   charge->capacity * sizeof(distributed_charge_klein_entry_t*));
		}

		/* Shift elements to make space */
		for (int i = charge->count; i > pos; i--) {
			charge->entries[i] = charge->entries[i - 1];
		}

		/* Create new entry */
		distributed_charge_klein_entry_t* entry = (distributed_charge_klein_entry_t*)malloc(sizeof(distributed_charge_klein_entry_t));
		entry->cs_index = cs_index;

		/* Store original values */
		psi_rho(obj, cs_index, 0, &entry->rho0_original);
		psi_rho(obj, cs_index, 1, &entry->rho1_original);

		/* Allocate and initialize Klein sums */
		entry->rho0_sum = (klein_t*)malloc(sizeof(klein_t));
		entry->rho1_sum = (klein_t*)malloc(sizeof(klein_t));
		*entry->rho0_sum = klein_zero();
		*entry->rho1_sum = klein_zero();

		/* Add first contribution */
		klein_add_double(entry->rho0_sum, q0_dr);
		klein_add_double(entry->rho1_sum, q1_dr);

		charge->entries[pos] = entry;
		charge->count++;
	}
}

/*****************************************************************************
 *
 *  subgrid_charge_from_particles()
 *
 *  For each particle, accumulate the charge on the relevant surrounding
 *  lattice nodes. Only nodes in the local domain are involved.
 *
 *****************************************************************************/
int subgrid_charge_from_particles(colloids_info_t* cinfo, psi_t* obj, distributed_charge_klein_t** charge)
{
	int ic, jc, kc;
	int i, j, k, i_min, i_max, j_min, j_max, k_min, k_max;
	int index;
	int nlocal[3], offset[3];
	int ncell[3];
	double r[3], r0[3];
	double dr, q0_dr, q1_dr;
	colloid_t* p_colloid = NULL;

	assert(cinfo);
	assert(obj);
	assert(charge);

	if (cinfo->nsubgrid == 0) return 0;

	cs_nlocal(cinfo->cs, nlocal);
	cs_nlocal_offset(cinfo->cs, offset);
	colloids_info_ncell(cinfo, ncell);

	/* Loop through all cells (including the halo cells) */
	for (ic = 0; ic <= ncell[X] + 1; ic++) {
		for (jc = 0; jc <= ncell[Y] + 1; jc++) {
			for (kc = 0; kc <= ncell[Z] + 1; kc++) {

				colloids_info_cell_list_head(cinfo, ic, jc, kc, &p_colloid);

				for (; p_colloid; p_colloid = p_colloid->next) {

					if (p_colloid->s.bc != COLLOID_BC_SUBGRID) continue;

					r0[X] = p_colloid->s.r[X] - 1.0 * offset[X];
					r0[Y] = p_colloid->s.r[Y] - 1.0 * offset[Y];
					r0[Z] = p_colloid->s.r[Z] - 1.0 * offset[Z];

					// Peskin - 2 neigbours
					i_min = imax(0, (int)floor(r0[X] - drange_));
					i_max = imin(nlocal[X] + 1, (int)ceil(r0[X] + drange_));
					j_min = imax(0, (int)floor(r0[Y] - drange_));
					j_max = imin(nlocal[Y] + 1, (int)ceil(r0[Y] + drange_));
					k_min = imax(0, (int)floor(r0[Z] - drange_));
					k_max = imin(nlocal[Z] + 1, (int)ceil(r0[Z] + drange_));

					// // first neigbours
					// i_min = imax(0, (int)floor(r0[X]));
					// i_max = imin(nlocal[X] + 1, (int)ceil(r0[X]));
					// j_min = imax(0, (int)floor(r0[Y]));
					// j_max = imin(nlocal[Y] + 1, (int)ceil(r0[Y]));
					// k_min = imax(0, (int)floor(r0[Z]));
					// k_max = imin(nlocal[Z] + 1, (int)ceil(r0[Z]));

					for (i = i_min; i <= i_max; i++) {
						for (j = j_min; j <= j_max; j++) {
							for (k = k_min; k <= k_max; k++) {

								index = cs_index(cinfo->cs, i, j, k);

								r[X] = r0[X] - 1.0 * i;
								r[Y] = r0[Y] - 1.0 * j;
								r[Z] = r0[Z] - 1.0 * k;

								/*CHANGE INIT - Use Poisson-Boltzmann weight instead of Peskin */
								dr = d_peskin(r[X]) * d_peskin(r[Y]) * d_peskin(r[Z]);
								// dr = d_idw(r0, i, j, k, i_min, i_max, j_min, j_max, k_min, k_max);
								// dr = d_isdw(r0, i, j, k, i_min, i_max, j_min, j_max, k_min, k_max);
								// dr = d_trilinear(r[X]) * d_trilinear(r[Y]) * d_trilinear(r[Z]);
								// dr = d_pb(r0, i, j, k, i_min, i_max, j_min, j_max, k_min, k_max, subgrid_get_kappa());
								/*CHANGE END - Use Poisson-Boltzmann weight instead of Peskin */

								q0_dr = p_colloid->s.q0 * dr;
								q1_dr = p_colloid->s.q1 * dr;

								/* Add to Klein sum array with binary search */
								add_charge_to_array(charge, index, q0_dr, q1_dr, obj);
							}
						}
					}
				}
			}
		}
	}

	/* Now apply all accumulated charges to psi */
	if (*charge != NULL) {
		for (int i = 0; i < (*charge)->count; i++) {
			distributed_charge_klein_entry_t* entry = (*charge)->entries[i];
			double new_rho0 = entry->rho0_original + klein_sum(entry->rho0_sum);
			double new_rho1 = entry->rho1_original + klein_sum(entry->rho1_sum);
			psi_rho_set(obj, entry->cs_index, 0, new_rho0);
			psi_rho_set(obj, entry->cs_index, 1, new_rho1);
		}
	}

	return 0;
}

/*****************************************************************************
 *
 *  subgrid_charge_from_particles_compenzate()
 *
 *  For each particle, accumulate the charge on the relevant surrounding
 *  lattice nodes. Only nodes in the local domain are involved.
 *
 *****************************************************************************/
int subgrid_charge_from_particles_compenzate(colloids_info_t* cinfo, psi_t* obj, distributed_charge_klein_t** charge)
{
	int ic, jc, kc;
	int i, j, k, i_min, i_max, j_min, j_max, k_min, k_max;
	int index;
	int nlocal[3], offset[3];
	int ncell[3];
	double r[3], r0[3];
	double dr, q0_dr, q1_dr;
	colloid_t* p_colloid = NULL;

	assert(cinfo);
	assert(obj);
	assert(charge);

	if (cinfo->nsubgrid == 0) return 0;

	cs_nlocal(cinfo->cs, nlocal);
	cs_nlocal_offset(cinfo->cs, offset);
	colloids_info_ncell(cinfo, ncell);

	/* Loop through all cells (including the halo cells) */
	for (ic = 0; ic <= ncell[X] + 1; ic++) {
		for (jc = 0; jc <= ncell[Y] + 1; jc++) {
			for (kc = 0; kc <= ncell[Z] + 1; kc++) {

				colloids_info_cell_list_head(cinfo, ic, jc, kc, &p_colloid);

				for (; p_colloid; p_colloid = p_colloid->next) {

					if (p_colloid->s.bc != COLLOID_BC_SUBGRID) continue;

					r0[X] = p_colloid->s.r[X] - 1.0 * offset[X];
					r0[Y] = p_colloid->s.r[Y] - 1.0 * offset[Y];
					r0[Z] = p_colloid->s.r[Z] - 1.0 * offset[Z];

					i_min = imax(0, (int)floor(r0[X] - drange_));
					i_max = imin(nlocal[X] + 1, (int)ceil(r0[X] + drange_));
					j_min = imax(0, (int)floor(r0[Y] - drange_));
					j_max = imin(nlocal[Y] + 1, (int)ceil(r0[Y] + drange_));
					k_min = imax(0, (int)floor(r0[Z] - drange_));
					k_max = imin(nlocal[Z] + 1, (int)ceil(r0[Z] + drange_));

					for (i = i_min; i <= i_max; i++) {
						for (j = j_min; j <= j_max; j++) {
							for (k = k_min; k <= k_max; k++) {

								index = cs_index(cinfo->cs, i, j, k);

								r[X] = r0[X] - 1.0 * i;
								r[Y] = r0[Y] - 1.0 * j;
								r[Z] = r0[Z] - 1.0 * k;

								dr = d_peskin(r[X]) * d_peskin(r[Y]) * d_peskin(r[Z]);
								// dr = d_isdw(r0, i, j, k, i_min, i_max, j_min, j_max, k_min, k_max);

								q0_dr = p_colloid->s.q0 * dr;
								q1_dr = p_colloid->s.q1 * dr;

								/* Add to Klein sum array with binary search */
								add_charge_to_array(charge, index, q0_dr, q1_dr, obj);
							}
						}
					}
				}
			}
		}
	}

	/* Now apply all accumulated charges to psi */
	if (*charge != NULL) {
		for (int i = 0; i < (*charge)->count; i++) {
			distributed_charge_klein_entry_t* entry = (*charge)->entries[i];
			double new_rho0 = -entry->rho0_original + klein_sum(entry->rho1_sum);
			double new_rho1 = -entry->rho1_original + klein_sum(entry->rho0_sum);
			psi_rho_set(obj, entry->cs_index, 0, new_rho0);
			psi_rho_set(obj, entry->cs_index, 1, new_rho1);
		}
	}

	return 0;
}

/*****************************************************************************
 *
 *  subgrid_charge_from_grid()
 *
 *  Copy original fluid charge from psi into charge (saving originals),
 *  then redistribute each node's charge via the Peskin kernel so that
 *  the field in psi becomes the Peskin-smoothed version.
 *
 *  The charge structure stores the pre-smoothing rho values so that
 *  subgrid_charge_from_particles_restore() can recover them.
 *
 *****************************************************************************/
/*CHANGE INIT - subgrid_charge_from_grid */
int subgrid_charge_from_grid(colloids_info_t* cinfo, psi_t* obj,
                              distributed_charge_klein_t** charge)
{
	int i, j, k;
	int i2, j2, k2;
	int i_min, i_max, j_min, j_max, k_min, k_max;
	int index, index2;
	int nlocal[3];
	double rho0, rho1, dr;
	int drange_i = (int)ceil(drange_);

	/* Temporary list of charged source nodes with their coordinates */
	typedef struct { int idx, si, sj, sk; double rho0, rho1; } src_t;
	int src_cap = 64, src_n = 0;
	src_t* srcs = (src_t*)malloc(src_cap * sizeof(src_t));

	assert(cinfo);
	assert(obj);
	assert(charge);
	assert(srcs);

	cs_nlocal(cinfo->cs, nlocal);

	/* Pass 1: collect charged nodes, save originals via charge, zero psi */
	for (i = 1; i <= nlocal[X]; i++) {
		for (j = 1; j <= nlocal[Y]; j++) {
			for (k = 1; k <= nlocal[Z]; k++) {

				index = cs_index(cinfo->cs, i, j, k);

				psi_rho(obj, index, 0, &rho0);
				psi_rho(obj, index, 1, &rho1);

				if (rho0 == 0.0 && rho1 == 0.0) continue;

				add_charge_to_array(charge, index, 0.0, 0.0, obj);

				psi_rho_set(obj, index, 0, 0.0);
				psi_rho_set(obj, index, 1, 0.0);

				if (src_n >= src_cap) {
					src_cap *= 2;
					srcs = (src_t*)realloc(srcs, src_cap * sizeof(src_t));
				}
				srcs[src_n].idx  = index;
				srcs[src_n].si   = i;
				srcs[src_n].sj   = j;
				srcs[src_n].sk   = k;
				srcs[src_n].rho0 = rho0;
				srcs[src_n].rho1 = rho1;
				src_n++;
			}
		}
	}

	/* Pass 2: Peskin-spread each source node's charge to neighbours */
	for (int e = 0; e < src_n; e++) {

		int si = srcs[e].si, sj = srcs[e].sj, sk = srcs[e].sk;
		rho0 = srcs[e].rho0;
		rho1 = srcs[e].rho1;

		i_min = imax(1, si - drange_i);
		i_max = imin(nlocal[X], si + drange_i);
		j_min = imax(1, sj - drange_i);
		j_max = imin(nlocal[Y], sj + drange_i);
		k_min = imax(1, sk - drange_i);
		k_max = imin(nlocal[Z], sk + drange_i);

		for (i2 = i_min; i2 <= i_max; i2++) {
			for (j2 = j_min; j2 <= j_max; j2++) {
				for (k2 = k_min; k2 <= k_max; k2++) {

					dr = d_peskin((double)(si - i2))
					   * d_peskin((double)(sj - j2))
					   * d_peskin((double)(sk - k2));
					if (dr == 0.0) continue;

					index2 = cs_index(cinfo->cs, i2, j2, k2);
					add_charge_to_array(charge, index2, rho0 * dr, rho1 * dr, obj);
				}
			}
		}
	}

	free(srcs);

	if (*charge == NULL) return 0;

	/* Pass 3: write Peskin-smoothed charges into psi */
	for (int e = 0; e < (*charge)->count; e++) {
		distributed_charge_klein_entry_t* entry = (*charge)->entries[e];
		psi_rho_set(obj, entry->cs_index, 0, klein_sum(entry->rho0_sum));
		psi_rho_set(obj, entry->cs_index, 1, klein_sum(entry->rho1_sum));
	}

	return 0;
}
/*CHANGE END - subgrid_charge_from_grid */

/*****************************************************************************
 *
 *  subgrid_charge_from_particles_restore()
 *
 *  Restore original charge values from before particles were added
 *
 *****************************************************************************/
int subgrid_charge_from_particles_restore(colloids_info_t* cinfo, psi_t* obj, distributed_charge_klein_t** charge)
{
	assert(cinfo);
	assert(obj);
	assert(charge);

	if (*charge == NULL) return 0;

	/* Restore original values */
	for (int i = 0; i < (*charge)->count; i++) {
		distributed_charge_klein_entry_t* entry = (*charge)->entries[i];
		psi_rho_set(obj, entry->cs_index, 0, entry->rho0_original);
		psi_rho_set(obj, entry->cs_index, 1, entry->rho1_original);
	}

	return 0;
}

/*****************************************************************************
 *
 *  subgrid_free_distributed_charge_t()
 *
 *  Free memory allocated for distributed charge structure
 *
 *****************************************************************************/
void subgrid_free_distributed_charge_t(distributed_charge_klein_t** charge)
{
	if (charge == NULL || *charge == NULL) return;

	distributed_charge_klein_t* c = *charge;

	/* Free each entry */
	for (int i = 0; i < c->count; i++) {
		if (c->entries[i] != NULL) {
			if (c->entries[i]->rho0_sum != NULL) free(c->entries[i]->rho0_sum);
			if (c->entries[i]->rho1_sum != NULL) free(c->entries[i]->rho1_sum);
			free(c->entries[i]);
		}
	}

	/* Free array and structure */
	if (c->entries != NULL) free(c->entries);
	free(c);
	*charge = NULL;
}

/*****************************************************************************
 *
 *  subgrid_charge_from_particles_substract()
 *
 *  For each particle, accumulate the charge on the relevant surrounding
 *  lattice nodes. Only nodes in the local domain are involved.
 *
 *****************************************************************************/
int subgrid_charge_from_particles_substract(colloids_info_t* cinfo, psi_t* obj, distributed_charge_klein_t** charge)
{
	int ic, jc, kc;
	int i, j, k, i_min, i_max, j_min, j_max, k_min, k_max;
	int index;
	int nlocal[3], offset[3];
	int ncell[3];
	double rho0, rho1;
	double r[3], r0[3];
	double dr;
	colloid_t* p_colloid = NULL;  /* Subgrid colloid */

	assert(cinfo);
	assert(obj);

	if (cinfo->nsubgrid == 0) return 0;

	/* Substrac if already have the sum o q0 and q1*/
	if (*charge != NULL) {
		for (int i = 0; i < (*charge)->count; i++) {
			distributed_charge_klein_entry_t* entry = (*charge)->entries[i];
			double new_rho0;
			double new_rho1;
			//Get the current rho
			psi_rho(obj, entry->cs_index, 0, &new_rho0);
			psi_rho(obj, entry->cs_index, 1, &new_rho1);

			//Subtract the accumulated sum
			new_rho0 -= klein_sum(entry->rho0_sum);
			new_rho1 -= klein_sum(entry->rho1_sum);

			//Set the new rho
			psi_rho_set(obj, entry->cs_index, 0, new_rho0);
			psi_rho_set(obj, entry->cs_index, 1, new_rho1);
		}

		return 0;
	}

	cs_nlocal(cinfo->cs, nlocal);
	cs_nlocal_offset(cinfo->cs, offset);
	colloids_info_ncell(cinfo, ncell);

	/* Loop through all cells (including the halo cells) */

	for (ic = 0; ic <= ncell[X] + 1; ic++) {
		for (jc = 0; jc <= ncell[Y] + 1; jc++) {
			for (kc = 0; kc <= ncell[Z] + 1; kc++) {

				colloids_info_cell_list_head(cinfo, ic, jc, kc, &p_colloid);

				for (; p_colloid; p_colloid = p_colloid->next) {

					if (p_colloid->s.bc != COLLOID_BC_SUBGRID) continue;

					/* Need to translate the colloid position to "local"
					 * coordinates, so that the correct range of lattice
					 * nodes is found */

					r0[X] = p_colloid->s.r[X] - 1.0 * offset[X];
					r0[Y] = p_colloid->s.r[Y] - 1.0 * offset[Y];
					r0[Z] = p_colloid->s.r[Z] - 1.0 * offset[Z];

					// Work out which lattice sites are involved including Halo

					// Peskin - 2 neigbours
					// i_min = imax(0, (int)floor(r0[X] - drange_));
					// i_max = imin(nlocal[X] + 1, (int)ceil(r0[X] + drange_));
					// j_min = imax(0, (int)floor(r0[Y] - drange_));
					// j_max = imin(nlocal[Y] + 1, (int)ceil(r0[Y] + drange_));
					// k_min = imax(0, (int)floor(r0[Z] - drange_));
					// k_max = imin(nlocal[Z] + 1, (int)ceil(r0[Z] + drange_));

					// // first neigbours
					i_min = imax(0, (int)floor(r0[X]));
					i_max = imin(nlocal[X] + 1, (int)ceil(r0[X]));
					j_min = imax(0, (int)floor(r0[Y]));
					j_max = imin(nlocal[Y] + 1, (int)ceil(r0[Y]));
					k_min = imax(0, (int)floor(r0[Z]));
					k_max = imin(nlocal[Z] + 1, (int)ceil(r0[Z]));

					for (i = i_min; i <= i_max; i++) {
						for (j = j_min; j <= j_max; j++) {
							for (k = k_min; k <= k_max; k++) {

								index = cs_index(cinfo->cs, i, j, k);

								/* Separation between r0 and the coordinate position of
								 * this site */

								r[X] = r0[X] - 1.0 * i;
								r[Y] = r0[Y] - 1.0 * j;
								r[Z] = r0[Z] - 1.0 * k;

								// dr = d_peskin(r[X]) * d_peskin(r[Y]) * d_peskin(r[Z]);
								// dr = d_idw(r0, i, j, k, i_min, i_max, j_min, j_max, k_min, k_max);
								dr = d_isdw(r0, i, j, k, i_min, i_max, j_min, j_max, k_min, k_max);

								psi_rho(obj, index, 0, &rho0);
								rho0 = rho0 - p_colloid->s.q0 * dr;

								psi_rho(obj, index, 1, &rho1);
								rho1 = rho1 - p_colloid->s.q1 * dr;

								psi_rho_set(obj, index, 0, rho0);
								psi_rho_set(obj, index, 1, rho1);

							}
						}
					}
					/* Next colloid */
				}
				/* Next cell */
			}
		}
	}
	return 0;
}

/*****************************************************************************
 *
 *  binary_search_force_index()
 *
 *  Binary search to find position to insert/find cs_index in sorted array
 *
 *****************************************************************************/
static int binary_search_force_index(distributed_force_klein_t* force, int cs_index) {
	int left = 0;
	int right = force->count - 1;

	while (left <= right) {
		int mid = left + (right - left) / 2;
		if (force->entries[mid]->cs_index == cs_index) {
			return mid;  /* Found */
		}
		if (force->entries[mid]->cs_index < cs_index) {
			left = mid + 1;
		}
		else {
			right = mid - 1;
		}
	}
	return left;  /* Position to insert */
}

/*****************************************************************************
 *
 *  add_force_to_array()
 *
 *  Add or accumulate force at given cs_index
 *
 *****************************************************************************/
static void add_force_to_array(distributed_force_klein_t** force_ptr, int cs_index,
								double force_add[3]) {
	distributed_force_klein_t* force = *force_ptr;

	/* Initialize if needed */
	if (force == NULL) {
		force = (distributed_force_klein_t*)malloc(sizeof(distributed_force_klein_t));
		force->entries = (distributed_force_klein_entry_t**)malloc(16 * sizeof(distributed_force_klein_entry_t*));
		force->count = 0;
		force->capacity = 16;
		*force_ptr = force;
	}

	/* Binary search for index */
	int pos = binary_search_force_index(force, cs_index);

	/* Check if index already exists */
	if (pos < force->count&& force->entries[pos]->cs_index == cs_index) {
		/* Accumulate to existing entry */
		klein_add_double(force->entries[pos]->force[0], force_add[0]);
		klein_add_double(force->entries[pos]->force[1], force_add[1]);
		klein_add_double(force->entries[pos]->force[2], force_add[2]);
	}
	else {
		/* Need to insert new entry */
		if (force->count >= force->capacity) {
			force->capacity *= 2;
			force->entries = (distributed_force_klein_entry_t**)realloc(force->entries,
							   force->capacity * sizeof(distributed_force_klein_entry_t*));
		}

		/* Shift elements to make space */
		for (int i = force->count; i > pos; i--) {
			force->entries[i] = force->entries[i - 1];
		}

		/* Create new entry */
		distributed_force_klein_entry_t* entry = (distributed_force_klein_entry_t*)malloc(sizeof(distributed_force_klein_entry_t));
		entry->cs_index = cs_index;

		/* Allocate and initialize Klein sums */
		// entry->force = (klein_t**)malloc(sizeof(klein_t*));
		entry->force = (klein_t**)malloc(3 * sizeof(klein_t*));
		for (int i = 0; i < 3; i++) {
			/*CHANGE INIT - 20251119 CUDA C++ compatibility fix */
		/* Original: entry->force[i] = malloc(sizeof(klein_t)); */
		/* C++ requires explicit cast from void* */
			entry->force[i] = (klein_t*)malloc(sizeof(klein_t));
			/*CHANGE END*/
			*(entry->force[i]) = klein_zero();
		}

		/* Add first contribution */
		klein_add_double(entry->force[0], force_add[0]);
		klein_add_double(entry->force[1], force_add[1]);
		klein_add_double(entry->force[2], force_add[2]);

		force->entries[pos] = entry;
		force->count++;
	}
}

/*****************************************************************************
 *
 *  subgrid_free_distributed_force_t()
 *
 *  Free memory allocated for distributed force structure
 *
 *****************************************************************************/
void subgrid_free_distributed_force_t(distributed_force_klein_t** force)
{
	if (force == NULL || *force == NULL) return;

	distributed_force_klein_t* c = *force;

	/* Free each entry */
	for (int i = 0; i < c->count; i++) {
		if (c->entries[i] != NULL) {
			if (c->entries[i]->force != NULL) free(c->entries[i]->force);
			free(c->entries[i]);
		}
	}

	/* Free array and structure */
	if (c->entries != NULL) free(c->entries);
	free(c);
	*force = NULL;
}

/*****************************************************************************
 *
 *  subgrid_update_forces_electrokinetics
 *
 *  Accumulate single particle force contributions from electric fields.
 *  Sum goes to fex as this force is calculated ouside of ludwig_colloids_update
 *  in order to apply the force in the same step as the applied field
 *****************************************************************************/
int subgrid_update_forces_electrokinetics(colloids_info_t* cinfo,
										  map_t* map,
										  physics_t* phys,
										  psi_t* psi,
										  hydro_t* hydro) {

	int i, j, k, ic, jc, kc, i_min, i_max, j_min, j_max, k_min, k_max;
	int ncell[3];
	int index;
	int nlocal[3], offset[3];
	int nsfluid;
	double kt, eunit, reunit, dr;
	double r[3], r0[3];
	double e[3];           /* electric field */
	// double E_field[3] = { 0.0, 0.0, 0.0 };      /* force on particle from this lattice site */
	double force[3] = { 0.0, 0.0, 0.0 };      /* force on particle from this lattice site */
	double flocal[4] = { 0.0, 0.0, 0.0, 0.0 }; /* cumulative forces and fluid node count */
	double fsum[4] = { 0.0, 0.0, 0.0, 0.0 }; /* global sum of forces and fluid node count */
	klein_t flocal_k[3];
	flocal_k[X] = klein_zero();
	flocal_k[Y] = klein_zero();
	flocal_k[Z] = klein_zero();
	distributed_force_klein_t* force_k_indexed = NULL;

	colloid_t* pc;
	MPI_Comm comm;

	assert(cinfo);
	assert(map);
	assert(psi);

	if (cinfo->nsubgrid == 0) return 0;

	cs_nlocal(cinfo->cs, nlocal);
	cs_nlocal_offset(cinfo->cs, offset);
	cs_cart_comm(cinfo->cs, &comm);
	colloids_info_ncell(cinfo, ncell);
	physics_kt(phys, &kt);
	psi_unit_charge(psi, &eunit);
	reunit = 1.0 / eunit;

	/* While there is no device implementation, must copy back-and forth
	 * the force. */
	assert(hydro);
	hydro_memcpy(hydro, tdpMemcpyDeviceToHost);

    // // INIT VERSION - Ewald											
	// /* Calculate electric forces on particles and accumulate total force */
	// for (ic = 0; ic <= ncell[X] + 1; ic++) {
	// 	for (jc = 0; jc <= ncell[Y] + 1; jc++) {
	// 		for (kc = 0; kc <= ncell[Z] + 1; kc++) {

	// 			colloids_info_cell_list_head(cinfo, ic, jc, kc, &pc);

	// 			for (; pc; pc = pc->next) {


	// 				if (pc->s.bc != COLLOID_BC_SUBGRID) continue;

	// 				klein_t force_k[3];
	// 				force_k[X] = klein_zero();
	// 				force_k[Y] = klein_zero();
	// 				force_k[Z] = klein_zero();

	// 				// Si uso Ewald	se asigna directamente la fuerza calculada en Ewald, que ya incluye el término de carga
	// 				force[X] = pc->fex[X];
	// 				force[Y] = pc->fex[Y];
	// 				force[Z] = pc->fex[Z];

	// 				/* Translate colloid position to local coordinates */
	// 				r0[X] = pc->s.r[X] - 1.0 * offset[X];
	// 				r0[Y] = pc->s.r[Y] - 1.0 * offset[Y];
	// 				r0[Z] = pc->s.r[Z] - 1.0 * offset[Z];

	// 				/* Work out which local lattice sites are involved */
	// 				subgrid_get_lattice_index(r0, nlocal, &i_min, &i_max, &j_min, &j_max, &k_min, &k_max);

	// 				for (i = i_min; i <= i_max; i++) {
	// 					for (j = j_min; j <= j_max; j++) {
	// 						for (k = k_min; k <= k_max; k++) {

	// 							double force_aux[3] = { 0.0, 0.0, 0.0 };      /* force on particle from this lattice site */

	// 							index = cs_index(cinfo->cs, i, j, k);

	// 							/* Separation between r0 and the lattice site */
	// 							r[X] = r0[X] - 1.0 * i;
	// 							r[Y] = r0[Y] - 1.0 * j;
	// 							r[Z] = r0[Z] - 1.0 * k;

	// 							dr = d_peskin(r[X]) * d_peskin(r[Y]) * d_peskin(r[Z]);

	// 							/* Force on particle from electric field at this site index*/
	// 							force_aux[X] = force[X] * dr;
	// 							force_aux[Y] = force[Y] * dr;
	// 							force_aux[Z] = force[Z] * dr;

	// 							/* Add to Klein sum array with binary search */
	// 							add_force_to_array(&force_k_indexed, index, force_aux);

	// 						}
	// 					}
	// 				}
	// 			}
	// 		}
	// 	}
	// }
	// // // END VERSION - Ewald

	// INIT VERSION - Peskin 
	/* Calculate electric forces on particles and accumulate total force */
	for (ic = 0; ic <= ncell[X] + 1; ic++) {
		for (jc = 0; jc <= ncell[Y] + 1; jc++) {
			for (kc = 0; kc <= ncell[Z] + 1; kc++) {

				colloids_info_cell_list_head(cinfo, ic, jc, kc, &pc);

				for (; pc; pc = pc->next) {


					if (pc->s.bc != COLLOID_BC_SUBGRID) continue;

					// Force goes to fex as this force is calculated ouside of ludwig_colloids_update
					force[X] = kt * reunit * pc->Esub[X] * (pc->s.q0 - pc->s.q1);
					force[Y] = kt * reunit * pc->Esub[Y] * (pc->s.q0 - pc->s.q1);
					force[Z] = kt * reunit * pc->Esub[Z] * (pc->s.q0 - pc->s.q1);
					// pc->fex[X] = kt * reunit * pc->Esub[X] * (pc->s.q0 - pc->s.q1);
					// pc->fex[Y] = kt * reunit * pc->Esub[Y] * (pc->s.q0 - pc->s.q1);
					// pc->fex[Z] = kt * reunit * pc->Esub[Z] * (pc->s.q0 - pc->s.q1);


					/* Translate colloid position to local coordinates */
					r0[X] = pc->s.r[X] - 1.0 * offset[X];
					r0[Y] = pc->s.r[Y] - 1.0 * offset[Y];
					r0[Z] = pc->s.r[Z] - 1.0 * offset[Z];

					/* Work out which local lattice sites are involved */
					subgrid_get_lattice_index(r0, nlocal, &i_min, &i_max, &j_min, &j_max, &k_min, &k_max);

					for (i = i_min; i <= i_max; i++) {
						for (j = j_min; j <= j_max; j++) {
							for (k = k_min; k <= k_max; k++) {

								double force_aux[3] = { 0.0, 0.0, 0.0 };  /* force on particle from this lattice site */

								index = cs_index(cinfo->cs, i, j, k);

								/* Separation between r0 and the lattice site */
								r[X] = r0[X] - 1.0 * i;
								r[Y] = r0[Y] - 1.0 * j;
								r[Z] = r0[Z] - 1.0 * k;

								dr = d_peskin(r[X]) * d_peskin(r[Y]) * d_peskin(r[Z]);

								/* Force on particle from electric field at this site index*/
								force_aux[X] = force[X] * dr;
								force_aux[Y] = force[Y] * dr;
								force_aux[Z] = force[Z] * dr;

								/* Add to Klein sum array with binary search */
								add_force_to_array(&force_k_indexed, index, force_aux);

							}
						}
					}

					pc->fex[X] += force[X];
					pc->fex[Y] += force[Y];
					pc->fex[Z] += force[Z];

				}
			}
		}
	}
	// END VERSION - Peskin
	
	colloid_sums_halo(cinfo, COLLOID_SUM_FORCE_EXT_ONLY);

	/* Now apply all accumulated external electric forces to hydro by index */
	if (force_k_indexed != NULL) {
		for (int i = 0; i < (force_k_indexed)->count; i++) {

			distributed_force_klein_entry_t* entry = (force_k_indexed)->entries[i];

			double new_force[3];
			new_force[X] = klein_sum(entry->force[X]);
			new_force[Y] = klein_sum(entry->force[Y]);
			new_force[Z] = klein_sum(entry->force[Z]);

			hydro_f_local_add(hydro, entry->cs_index, new_force);

		}

		// Free memory
		subgrid_free_distributed_force_t(&force_k_indexed);

	}

	hydro_memcpy(hydro, tdpMemcpyHostToDevice);

	return 0;

}

/*****************************************************************************
 *
 *  subgrid_update_forces_electrokinetics_theory
 *
 *  Calculate electric forces on particles using direct Coulomb interaction
 *  with Debye-Hückel screening from charge densities at nodes.
 *
 *  For Poisson-Boltzmann, calculates F = q_particle * sum_over_nodes(E_yukawa)
 *  where E_yukawa includes Debye screening:
 *    E = (prefactor) * (rho0 - rho1) * exp(-κr) * (1/r² + κ/r) * r_ij
 *  with prefactor = eunit * kt / (4*pi*epsilon) and κ = 1/λ_Debye
 *
 *  The force on particle is calculated WITHOUT Peskin weighting.
 *  Only the reaction force to fluid uses Peskin distribution.
 *
 *  Sum goes to fex as this force is calculated outside of ludwig_colloids_update.
 *****************************************************************************/
int subgrid_update_forces_electrokinetics_theory(colloids_info_t* cinfo,
												  map_t* map,
												  physics_t* phys,
												  psi_t* psi,
												  hydro_t* hydro) {

	int i, j, k, ic, jc, kc, i_min, i_max, j_min, j_max, k_min, k_max;
	int ncell[3];
	int index;
	int nlocal[3], offset[3];
	double kt, eunit, epsilon, dr_peskin;
	double r_particle_local[3], r_particle_global[3];
	double r_node_global[3], r_ij[3];
	double dist;
	double rho0, rho1, rho_net;
	double prefactor;
	double kappa;              /* Inverse Debye length */
	double lambda_debye;       /* Debye length */
	double rho_bulk;           /* Bulk ionic strength for Debye calculation */
	double lb;                 /* Bjerrum length */
	PI_DOUBLE(pi);
	distributed_force_klein_t* force_k_indexed = NULL;

	colloid_t* pc;
	MPI_Comm comm;

	assert(cinfo);
	assert(map);
	assert(psi);

	if (cinfo->nsubgrid == 0) return 0;

	cs_nlocal(cinfo->cs, nlocal);
	cs_nlocal_offset(cinfo->cs, offset);
	cs_cart_comm(cinfo->cs, &comm);
	colloids_info_ncell(cinfo, ncell);

	physics_kt(phys, &kt);
	psi_unit_charge(psi, &eunit);
	psi_epsilon(psi, &epsilon);

	/* Coulomb prefactor: (eunit * kt) / (4*pi*epsilon) */
	prefactor = (eunit * kt) / (4.0 * pi * epsilon);

	/* Calculate Debye screening parameter κ */
	/* Bjerrum length: l_B = e^2 / (4*pi*epsilon*kT) */
	lb = (eunit * eunit) / (4.0 * pi * epsilon * kt);

	/* Get bulk ionic density for Debye length calculation */
	/* Using average of rho0 and rho1 at a bulk location as estimate */
	/* For symmetric electrolyte: λ_D = 1/sqrt(8π × l_B × rho_b) */
	/* This is a simplification - ideally should get from psi_options */
	psi_rho(psi, cs_index(cinfo->cs, nlocal[X] / 2, nlocal[Y] / 2, nlocal[Z] / 2), 0, &rho0);
	psi_rho(psi, cs_index(cinfo->cs, nlocal[X] / 2, nlocal[Y] / 2, nlocal[Z] / 2), 1, &rho1);
	rho_bulk = 0.5 * (rho0 + rho1);

	if (rho_bulk > 1e-15) {
		lambda_debye = 1.0 / sqrt(8.0 * pi * lb * rho_bulk);
		kappa = 1.0 / lambda_debye;
	}
	else {
		/* No screening if no bulk electrolyte */
		kappa = 0.0;
	}

	/* While there is no device implementation, must copy back-and forth */
	assert(hydro);
	hydro_memcpy(hydro, tdpMemcpyDeviceToHost);

	/* Calculate electric forces using direct Coulomb interaction from charge at nodes */
	for (ic = 0; ic <= ncell[X] + 1; ic++) {
		for (jc = 0; jc <= ncell[Y] + 1; jc++) {
			for (kc = 0; kc <= ncell[Z] + 1; kc++) {

				colloids_info_cell_list_head(cinfo, ic, jc, kc, &pc);

				for (; pc; pc = pc->next) {

					if (pc->s.bc != COLLOID_BC_SUBGRID) continue;

					klein_t force_particle_k[3];
					force_particle_k[X] = klein_zero();
					force_particle_k[Y] = klein_zero();
					force_particle_k[Z] = klein_zero();

					/* Net charge on particle */
					double q_particle = (pc->s.q0 - pc->s.q1);

					/* Particle position in global coordinates */
					r_particle_global[X] = pc->s.r[X];
					r_particle_global[Y] = pc->s.r[Y];
					r_particle_global[Z] = pc->s.r[Z];

					/* Work out which local lattice sites are involved (Peskin range) */
					subgrid_get_lattice_index(r_particle_local, nlocal,
											  &i_min, &i_max, &j_min, &j_max, &k_min, &k_max);

					/* Loop over nodes in Peskin range to calculate total force on particle */
					for (i = i_min; i <= i_max; i++) {
						for (j = j_min; j <= j_max; j++) {
							for (k = k_min; k <= k_max; k++) {

								double E_coulomb[3] = { 0.0, 0.0, 0.0 };
								double force_from_node[3] = { 0.0, 0.0, 0.0 };

								index = cs_index(cinfo->cs, i, j, k);

								/* Get charge density at this node */
								psi_rho(psi, index, 0, &rho0);
								psi_rho(psi, index, 1, &rho1);
								rho_net = rho0 - rho1;  /* Net charge density at node */

								if (fabs(rho_net) < 1e-15) continue; /* Skip if no charge */

								/* Node position in global coordinates */
								r_node_global[X] = 1.0 * (i + offset[X]);
								r_node_global[Y] = 1.0 * (j + offset[Y]);
								r_node_global[Z] = 1.0 * (k + offset[Z]);

								/* Vector from node to particle */
								r_ij[X] = r_particle_global[X] - r_node_global[X];
								r_ij[Y] = r_particle_global[Y] - r_node_global[Y];
								r_ij[Z] = r_particle_global[Z] - r_node_global[Z];

								dist = sqrt(r_ij[X] * r_ij[X] + r_ij[Y] * r_ij[Y] + r_ij[Z] * r_ij[Z]);

								/* Skip if too close to avoid singularity */
								if (dist < 0.1) continue;

								/* Electric field with Debye-Hückel screening (Yukawa potential) */
								/* E = prefactor * q * exp(-κr) * (1/r² + κ/r) * r_unit */
								/* For κ=0, this reduces to standard Coulomb */
								double exp_kr = (kappa > 0.0) ? exp(-kappa * dist) : 1.0;
								double screening_factor = exp_kr * (1.0 / dist + kappa) / dist;

								E_coulomb[X] = prefactor * rho_net * r_ij[X] * screening_factor;
								E_coulomb[Y] = prefactor * rho_net * r_ij[Y] * screening_factor;
								E_coulomb[Z] = prefactor * rho_net * r_ij[Z] * screening_factor;

								/* Force on particle from this node: F = q_particle * E_coulomb (NO PESKIN) */
								force_from_node[X] = q_particle * E_coulomb[X];
								force_from_node[Y] = q_particle * E_coulomb[Y];
								force_from_node[Z] = q_particle * E_coulomb[Z];

								/* Add to total particle force using Klein summation */
								klein_add_double(&force_particle_k[X], force_from_node[X]);
								klein_add_double(&force_particle_k[Y], force_from_node[Y]);
								klein_add_double(&force_particle_k[Z], force_from_node[Z]);

							}
						}
					}

					/* Total force on particle */
					double force_particle[3];
					force_particle[X] = klein_sum(&force_particle_k[X]);
					force_particle[Y] = klein_sum(&force_particle_k[Y]);
					force_particle[Z] = klein_sum(&force_particle_k[Z]);

					/* Apply force to particle */
					pc->fex[X] += force_particle[X];
					pc->fex[Y] += force_particle[Y];
					pc->fex[Z] += force_particle[Z];

					/* Now distribute the reaction force to fluid nodes using Peskin */
					for (i = i_min; i <= i_max; i++) {
						for (j = j_min; j <= j_max; j++) {
							for (k = k_min; k <= k_max; k++) {

								index = cs_index(cinfo->cs, i, j, k);

								/* Peskin weight for this node */
								double r_rel[3];
								r_rel[X] = r_particle_local[X] - 1.0 * i;
								r_rel[Y] = r_particle_local[Y] - 1.0 * j;
								r_rel[Z] = r_particle_local[Z] - 1.0 * k;

								dr_peskin = d_peskin(r_rel[X]) * d_peskin(r_rel[Y]) * d_peskin(r_rel[Z]);

								/* Distribute reaction force to fluid with Peskin weight */
								double force_hydro[3];
								force_hydro[X] = -force_particle[X] * dr_peskin;
								force_hydro[Y] = -force_particle[Y] * dr_peskin;
								force_hydro[Z] = -force_particle[Z] * dr_peskin;

								add_force_to_array(&force_k_indexed, index, force_hydro);

							}
						}
					}

				}
			}
		}
	}

	colloid_sums_halo(cinfo, COLLOID_SUM_FORCE_EXT_ONLY);

	/* Apply accumulated electric forces to hydro by index */
	if (force_k_indexed != NULL) {
		for (int i = 0; i < (force_k_indexed)->count; i++) {

			distributed_force_klein_entry_t* entry = (force_k_indexed)->entries[i];

			double new_force[3];
			new_force[X] = klein_sum(entry->force[X]);
			new_force[Y] = klein_sum(entry->force[Y]);
			new_force[Z] = klein_sum(entry->force[Z]);

			hydro_f_local_add(hydro, entry->cs_index, new_force);

		}

		subgrid_free_distributed_force_t(&force_k_indexed);
	}

	hydro_memcpy(hydro, tdpMemcpyHostToDevice);

	return 0;
}

/*****************************************************************************
 *
 *  subgrid_update_Esub
 *
 *  Calculate electric field on subgrid particles
 *****************************************************************************/
int subgrid_update_Esub(colloids_info_t* cinfo,
						psi_t* psi,
						int step,
						FILE* fp,
						pe_t* pe) {

	int i, j, k, ic, jc, kc, i_min, i_max, j_min, j_max, k_min, k_max;
	int ncell[3];
	int index;
	int nlocal[3], offset[3];
	double dr;
	double r[3], r0[3];
	double phi_node = 0.0;
	double e[3];           						/* Total electric field */
	double E_field[3] = { 0.0, 0.0, 0.0 };      /* Electric field on particle from this lattice site */
	double E_self[3] = { 0.0, 0.0, 0.0 };      	/* Self-field */

	/* Variables for PB method calculation */
	double r_particle_global[3];
	double r_node_global[3], r_ij[3];
	double dist;
	double rho0, rho1, rho_net;
	double kt, eunit, epsilon;
	double prefactor;
	double kappa;
	PI_DOUBLE(pi);

	colloid_t* pc;
	MPI_Comm comm;

	// FILE* fp;
	int rank;

	assert(cinfo);
	assert(psi);

	if (cinfo->nsubgrid == 0) return 0;

	cs_nlocal(cinfo->cs, nlocal);
	cs_nlocal_offset(cinfo->cs, offset);
	cs_cart_comm(cinfo->cs, &comm);
	colloids_info_ncell(cinfo, ncell);
	MPI_Comm_rank(comm, &rank);

	/* Get parameters for PB calculation */
	psi_unit_charge(psi, &eunit);
	psi_epsilon(psi, &epsilon);
	psi_beta(psi, &kt);
	kt = 1.0 / kt;  /* Convert beta to kT */

	/* Coulomb prefactor: (eunit * kt) / (4*pi*epsilon) */
	prefactor = (eunit * kt) / (4.0 * pi * epsilon);

	/* Use stored kappa value */
	kappa = subgrid_get_kappa();

	// /* Initialize Esub*/
	// for (ic = 0; ic <= ncell[X] + 1; ic++) {
	// 	for (jc = 0; jc <= ncell[Y] + 1; jc++) {
	// 		for (kc = 0; kc <= ncell[Z] + 1; kc++) {

	// 			colloids_info_cell_list_head(cinfo, ic, jc, kc, &pc);

	// 			for (; pc; pc = pc->next) {

	// 				if (pc->s.bc != COLLOID_BC_SUBGRID) continue;

	// 				pc->Esub[X] = 0.0;
	// 				pc->Esub[Y] = 0.0;
	// 				pc->Esub[Z] = 0.0;

	// 			}
	// 		}
	// 	}
	// }

	// /* Second pass: Calculate electric field on particles from Eflied on nodes*/
	// for (ic = 0; ic <= ncell[X] + 1; ic++) {
	// 	for (jc = 0; jc <= ncell[Y] + 1; jc++) {
	// 		for (kc = 0; kc <= ncell[Z] + 1; kc++) {

	// 			colloids_info_cell_list_head(cinfo, ic, jc, kc, &pc);

	// 			for (; pc; pc = pc->next) {

	// 				if (pc->s.bc != COLLOID_BC_SUBGRID) continue;

	// 				/* METHOD 1: Current method using Peskin weights on local nodes */
	// 				klein_t E_field_k[3];
	// 				E_field_k[X] = klein_zero();
	// 				E_field_k[Y] = klein_zero();
	// 				E_field_k[Z] = klein_zero();

	// 				/* Translate colloid position to local coordinates */
	// 				r0[X] = pc->s.r[X] - 1.0 * offset[X];
	// 				r0[Y] = pc->s.r[Y] - 1.0 * offset[Y];
	// 				r0[Z] = pc->s.r[Z] - 1.0 * offset[Z];

	// 				/* Work out which local lattice sites are involved */
	// 				subgrid_get_lattice_index(r0, nlocal, &i_min, &i_max, &j_min, &j_max, &k_min, &k_max);
	// 				// subgrid_get_lattice_index_fn(r0, nlocal, &i_min, &i_max, &j_min, &j_max, &k_min, &k_max);

	// 				for (i = i_min; i <= i_max; i++) {
	// 					for (j = j_min; j <= j_max; j++) {
	// 						for (k = k_min; k <= k_max; k++) {

	// 							index = cs_index(cinfo->cs, i, j, k);

	// 							/* Separation between r0 and the lattice site */
	// 							r[X] = r0[X] - 1.0 * i;
	// 							r[Y] = r0[Y] - 1.0 * j;
	// 							r[Z] = r0[Z] - 1.0 * k;

	// 							/*CHANGE INIT - Use Poisson-Boltzmann weight instead of Peskin */
	// 							dr = d_peskin(r[X]) * d_peskin(r[Y]) * d_peskin(r[Z]);
	// 							// dr = d_idw(r0, i, j, k, i_min, i_max, j_min, j_max, k_min, k_max);
	// 							// dr = d_isdw(r0, i, j, k, i_min, i_max, j_min, j_max, k_min, k_max);
	// 							// dr = d_trilinear(r[X]) * d_trilinear(r[Y]) * d_trilinear(r[Z]);
	// 							// dr = d_pb(r0, i, j, k, i_min, i_max, j_min, j_max, k_min, k_max, subgrid_get_kappa());
	// 							/*CHANGE END - Use Poisson-Boltzmann weight instead of Peskin */

	// 							/* Electric field at this lattice site */
	// 							psi_electric_field(psi, index, e);

	// 							// pe_info(pe, "campo en el indice %d : Ex=%.20f   Ey=%.20f   Ez=%.20f\n", index, e[X], e[Y], e[Z]);

	// 							/* Field on particle from electric field at this site */
	// 							E_field[X] = e[X] * dr;
	// 							E_field[Y] = e[Y] * dr;
	// 							E_field[Z] = e[Z] * dr;

	// 							/* Add to particle field */
	// 							klein_add_double(&E_field_k[X], E_field[X]);
	// 							klein_add_double(&E_field_k[Y], E_field[Y]);
	// 							klein_add_double(&E_field_k[Z], E_field[Z]);

	// 						}
	// 					}
	// 				}

	// 				pc->Esub[X] = klein_sum(&E_field_k[X]);
	// 				pc->Esub[Y] = klein_sum(&E_field_k[Y]);
	// 				pc->Esub[Z] = klein_sum(&E_field_k[Z]);

	// 				// /* METHOD 2: PB method - Direct Coulomb sum over ALL nodes */
	// 				// klein_t E_PB_k[3];
	// 				// E_PB_k[X] = klein_zero();
	// 				// E_PB_k[Y] = klein_zero();
	// 				// E_PB_k[Z] = klein_zero();

	// 				// /* Particle position in global coordinates */
	// 				// r_particle_global[X] = pc->s.r[X];
	// 				// r_particle_global[Y] = pc->s.r[Y];
	// 				// r_particle_global[Z] = pc->s.r[Z];

	// 				// /* Loop over ALL local lattice sites */
	// 				// for (i = 1; i <= nlocal[X]; i++) {
	// 				// 	for (j = 1; j <= nlocal[Y]; j++) {
	// 				// 		for (k = 1; k <= nlocal[Z]; k++) {

	// 				// 			double E_coulomb[3] = { 0.0, 0.0, 0.0 };

	// 				// 			index = cs_index(cinfo->cs, i, j, k);

	// 				// 			/* Get charge density at this node */
	// 				// 			psi_rho(psi, index, 0, &rho0);
	// 				// 			psi_rho(psi, index, 1, &rho1);
	// 				// 			rho_net = rho0 - rho1;  /* Net charge density at node */

	// 				// 			/* Node position in global coordinates */
	// 				// 			r_node_global[X] = 1.0 * (i + offset[X]);
	// 				// 			r_node_global[Y] = 1.0 * (j + offset[Y]);
	// 				// 			r_node_global[Z] = 1.0 * (k + offset[Z]);

	// 				// 			/* Vector from node to particle */
	// 				// 			r_ij[X] = r_particle_global[X] - r_node_global[X];
	// 				// 			r_ij[Y] = r_particle_global[Y] - r_node_global[Y];
	// 				// 			r_ij[Z] = r_particle_global[Z] - r_node_global[Z];

	// 				// 			dist = sqrt(r_ij[X] * r_ij[X] + r_ij[Y] * r_ij[Y] + r_ij[Z] * r_ij[Z]);

	// 				// 			/* Skip if too close to avoid singularity */
	// 				// 			if (dist < 0.000001) continue;

	// 				// 			/* Electric field with Debye-Hückel screening (Yukawa potential) */
	// 				// 			/* E = prefactor * q * exp(-κr) * (1/r² + κ/r) * r_unit */
	// 				// 			double exp_kr = (kappa > 0.0) ? exp(-kappa * dist) : 1.0;
	// 				// 			double screening_factor = exp_kr * (1.0 / dist + kappa) / dist;

	// 				// 			E_coulomb[X] = prefactor * rho_net * r_ij[X] * screening_factor;
	// 				// 			E_coulomb[Y] = prefactor * rho_net * r_ij[Y] * screening_factor;
	// 				// 			E_coulomb[Z] = prefactor * rho_net * r_ij[Z] * screening_factor;

	// 				// 			/* Accumulate field contribution using Klein summation */
	// 				// 			klein_add_double(&E_PB_k[X], E_coulomb[X]);
	// 				// 			klein_add_double(&E_PB_k[Y], E_coulomb[Y]);
	// 				// 			klein_add_double(&E_PB_k[Z], E_coulomb[Z]);

	// 				// 		}
	// 				// 	}
	// 				// }

	// 				// /* Store PB field in temporary array for output */
	// 				// pc->fc0[X] = klein_sum(&E_PB_k[X]);
	// 				// pc->fc0[Y] = klein_sum(&E_PB_k[Y]);
	// 				// pc->fc0[Z] = klein_sum(&E_PB_k[Z]);

	// 				/* CHANGE INIT - Moved fprintf to after colloid_sums_halo to avoid writing halo particles */
	// 				// // /* Write data for each particle */
	// 				// char string[256];
	// 				// double Emod = sqrt(pc->Esub[X] * pc->Esub[X] + pc->Esub[Y] * pc->Esub[Y] + pc->Esub[Z] * pc->Esub[Z]);
	// 				// // sprintf(string, "%d;%d;%.15e;%.15e; %.15e;%.15e;%.15e; %.15e;%.15e\n", step, pc->s.index,  Emod, pc->Esub[X], pc->Esub[Y], pc->Esub[Z], E_self[X], E_self[Y], E_self[Z]);
	// 				// sprintf(string, "%d;%d;%.15e;%.15e; %.15e;%.15e\n", step, pc->s.index, Emod, pc->Esub[X], pc->Esub[Y], pc->Esub[Z]);
	// 				// for (i = 0;i < strlen(string);i++)if (string[i] == '.')string[i] = ',';
	// 				// fprintf(fp, "%s", string);
	// 				/* CHANGE END */

	// 				// /* NUEVO: Calcular el autocampo de la partícula */
	// 				// if (subgrid_compute_self_field_single_particle(pc, cinfo, psi, E_self) == 0) {

	// 				// 	/* Write data for each particle */
	// 				// 	char string[256];
	// 				// 	double Emod = sqrt(pc->Esub[X]*pc->Esub[X] + pc->Esub[Y]*pc->Esub[Y] + pc->Esub[Z]*pc->Esub[Z]);
	// 				// 	sprintf(string, "%d;%d;%.15e;%.15e; %.15e;%.15e;%.15e; %.15e;%.15e\n", step, pc->s.index,  Emod, pc->Esub[X], pc->Esub[Y], pc->Esub[Z], E_self[X], E_self[Y], E_self[Z]);
	// 				// 	for (i = 0;i < strlen(string);i++)if (string[i] == '.')string[i] = ',';
	// 				// 	fprintf(fp, "%s", string);

	// 				// 	/* Restar el autocampo del campo total */
	// 				// 	pc->Esub[X] -= E_self[X];
	// 				// 	pc->Esub[Y] -= E_self[Y];
	// 				// 	pc->Esub[Z] -= E_self[Z];

	// 				// }

	// 			}
	// 		}
	// 	}
	// }

	// INIT VERSION - Peskin all
	{
		/* Second pass: Calculate electric field on particles from Phi on nodes*/
		for (ic = 0; ic <= ncell[X] + 1; ic++) {
			for (jc = 0; jc <= ncell[Y] + 1; jc++) {
				for (kc = 0; kc <= ncell[Z] + 1; kc++) {

					colloids_info_cell_list_head(cinfo, ic, jc, kc, &pc);

					for (; pc; pc = pc->next) {

						if (pc->s.bc != COLLOID_BC_SUBGRID) continue;

						/* METHOD 1: Current method using Peskin weights on local nodes */
						klein_t E_field_k[3];
						E_field_k[X] = klein_zero();
						E_field_k[Y] = klein_zero();
						E_field_k[Z] = klein_zero();

						/* Translate colloid position to local coordinates */
						r0[X] = pc->s.r[X] - 1.0 * offset[X];
						r0[Y] = pc->s.r[Y] - 1.0 * offset[Y];
						r0[Z] = pc->s.r[Z] - 1.0 * offset[Z];

						/* Work out which local lattice sites are involved */
						subgrid_get_lattice_index(r0, nlocal, &i_min, &i_max, &j_min, &j_max, &k_min, &k_max);

						for (i = i_min; i <= i_max; i++) {
							for (j = j_min; j <= j_max; j++) {
								for (k = k_min; k <= k_max; k++) {

									index = cs_index(cinfo->cs, i, j, k);

									/* Separation between r0 and the lattice site */
									r[X] = r0[X] - 1.0 * i;
									r[Y] = r0[Y] - 1.0 * j;
									r[Z] = r0[Z] - 1.0 * k;

									dr = d_peskin(r[X]) * d_peskin(r[Y]) * d_peskin(r[Z]);

									/* Electric field at this lattice site */
									psi_psi(psi, index, &phi_node);

									// pe_info(pe, "campo en el indice %d : Ex=%.20f   Ey=%.20f   Ez=%.20f\n", index, e[X], e[Y], e[Z]);

									/* Electric field at this lattice site */
									psi_electric_field(psi, index, e);

									// pe_info(pe, "campo en el indice %d : Ex=%.20f   Ey=%.20f   Ez=%.20f\n", index, e[X], e[Y], e[Z]);

									/* Field on particle from electric field at this site */
									E_field[X] = e[X] * dr;
									E_field[Y] = e[Y] * dr;
									E_field[Z] = e[Z] * dr;
									// E_field[X] = phi_node * d_peskin_derivative(r[X]) * d_peskin(r[Y]) * d_peskin(r[Z]);
									// E_field[Y] = phi_node * d_peskin(r[X]) * d_peskin_derivative(r[Y]) * d_peskin(r[Z]);
									// E_field[Z] = phi_node * d_peskin(r[X]) * d_peskin(r[Y]) * d_peskin_derivative(r[Z]);
									
									/* Add to particle field */
									klein_add_double(&E_field_k[X], E_field[X]);
									klein_add_double(&E_field_k[Y], E_field[Y]);
									klein_add_double(&E_field_k[Z], E_field[Z]);

								}
							}
						}

						pc->Esub[X] = klein_sum(&E_field_k[X]);
						pc->Esub[Y] = klein_sum(&E_field_k[Y]);
						pc->Esub[Z] = klein_sum(&E_field_k[Z]);

					}
				}
			}
		}
	}
	// END VERSION - Peskin all

	colloid_sums_halo(cinfo, COLLOID_SUM_ELECTRIC_FIELD);

	/* CHANGE INIT - Write particle data after colloid_sums_halo, only for local (non-halo) particles */
	/* Loop only over non-halo cells to write data */
	for (ic = 1; ic <= ncell[X]; ic++) {
		for (jc = 1; jc <= ncell[Y]; jc++) {
			for (kc = 1; kc <= ncell[Z]; kc++) {

				colloids_info_cell_list_head(cinfo, ic, jc, kc, &pc);

				for (; pc; pc = pc->next) {

					if (pc->s.bc != COLLOID_BC_SUBGRID) continue;

					/* Write data for each particle - including both methods */
					char string[512];
					double Emod = sqrt(pc->Esub[X] * pc->Esub[X] + pc->Esub[Y] * pc->Esub[Y] + pc->Esub[Z] * pc->Esub[Z]);
					// double Emod_PB = sqrt(pc->fc0[X] * pc->fc0[X] + pc->fc0[Y] * pc->fc0[Y] + pc->fc0[Z] * pc->fc0[Z]);
					/* Format: step;index;Emod_Peskin;Ex_Peskin;Ey_Peskin;Ez_Peskin;Emod_PB;Ex_PB;Ey_PB;Ez_PB */
					// sprintf(string, "%d;%d;%.15e;%.15e;%.15e;%.15e;%.15e;%.15e;%.15e;%.15e\n",
					sprintf(string, "%d;%d;%.15e;%.15e;%.15e;%.15e;%.15e;%.15e;%.15e;%.15e\n",
							step, pc->s.index,
							Emod, pc->Esub[X], pc->Esub[Y], pc->Esub[Z]);
					// Emod_PB, pc->fc0[X], pc->fc0[Y], pc->fc0[Z]);
					for (i = 0; i < strlen(string); i++) if (string[i] == '.') string[i] = ',';
					fprintf(fp, "%s", string);
				}
			}
		}
	}
	/* CHANGE END */

	return 0;

}

/*CHANGE INIT - Poisson-Boltzmann force calculation */
/*****************************************************************************
 *
 *  subgrid_force_poisson_boltzmann
 *
 *  Calculate force on subgrid particles by summing the screened Coulomb
 *  (Debye-Hückel) electric field contribution from ALL lattice nodes
 *  at the particle position.
 *
 *  For each node with charge density rho_el, compute the field at particle:
 *  E = (prefactor) * rho_el * exp(-κr) * (1/r² + κ/r) * r̂
 *
 *  Then: F = q_particle * E_total
 *
 *  Sum goes to fex as this force is calculated outside of ludwig_colloids_update.
 * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 * !!!!    Necesita repensarse para incluir imagens periodicas en todas las direcciones posibles
 * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 *****************************************************************************/
int subgrid_force_poisson_boltzmann(colloids_info_t* cinfo,
									 map_t* map,
									 physics_t* phys,
									 psi_t* psi,
									 hydro_t* hydro,
									 int step,
									 FILE* fp) {

	int i, j, k, ic, jc, kc;
	int ncell[3];
	int index;
	int nlocal[3], offset[3];
	double kt, eunit, epsilon;
	double r_particle_global[3];
	double r_node_global[3], r_ij[3];
	double dist;
	double rho0, rho1, rho_net;
	double prefactor;
	double kappa;              /* Inverse Debye length from subgrid_get_kappa() */
	PI_DOUBLE(pi);
	distributed_force_klein_t* force_k_indexed = NULL;

	colloid_t* pc;
	MPI_Comm comm;
	int rank;  /* MPI rank for file writing */

	assert(cinfo);
	assert(map);
	assert(psi);

	if (cinfo->nsubgrid == 0) return 0;

	int n, ntotal[3];
	cs_ntotal(cinfo->cs, ntotal);
	cs_nlocal(cinfo->cs, nlocal);
	cs_nlocal_offset(cinfo->cs, offset);
	cs_cart_comm(cinfo->cs, &comm);
	colloids_info_ncell(cinfo, ncell);
	MPI_Comm_rank(comm, &rank);

	physics_kt(phys, &kt);
	psi_unit_charge(psi, &eunit);
	psi_epsilon(psi, &epsilon);

	/* Coulomb prefactor: (eunit * kt) / (4*pi*epsilon) */
	prefactor = (eunit * kt) / (4.0 * pi * epsilon);

	/* Use stored kappa value from subgrid_compute_kappa */
	/* This ensures consistency with subgrid_add_psi_grad_from_particles */
	kappa = subgrid_get_kappa();

	/* While there is no device implementation, must copy back-and forth */
	assert(hydro);
	hydro_memcpy(hydro, tdpMemcpyDeviceToHost);

	/* Calculate electric forces using direct Coulomb from ALL nodes in system */
	for (ic = 0; ic <= ncell[X] + 1; ic++) {
		for (jc = 0; jc <= ncell[Y] + 1; jc++) {
			for (kc = 0; kc <= ncell[Z] + 1; kc++) {

				colloids_info_cell_list_head(cinfo, ic, jc, kc, &pc);

				for (; pc; pc = pc->next) {

					if (pc->s.bc != COLLOID_BC_SUBGRID) continue;

					klein_t E_total_k[3];
					E_total_k[X] = klein_zero();
					E_total_k[Y] = klein_zero();
					E_total_k[Z] = klein_zero();

					/* Net charge on particle */
					double q_particle = (pc->s.q0 - pc->s.q1);

					/* Particle position in global coordinates */
					r_particle_global[X] = pc->s.r[X];
					r_particle_global[Y] = pc->s.r[Y];
					r_particle_global[Z] = pc->s.r[Z];

					/* Loop over ALL local lattice sites */
					for (i = 1; i <= nlocal[X]; i++) {
						for (j = 1; j <= nlocal[Y]; j++) {
							for (k = 1; k <= nlocal[Z]; k++) {

								double E_coulomb[3] = { 0.0, 0.0, 0.0 };

								index = cs_index(cinfo->cs, i, j, k);

								/* Skip solid sites */
								int status;
								map_status(map, index, &status);
								if (status != MAP_FLUID) continue;

								/* Get charge density at this node */
								psi_rho(psi, index, 0, &rho0);
								psi_rho(psi, index, 1, &rho1);
								rho_net = rho0 - rho1;  /* Net charge density at node */

								// if (fabs(rho_net) < 1e-17) continue; /* Skip if no charge */

								/* Node position in global coordinates */
								r_node_global[X] = 1.0 * (i + offset[X]);
								r_node_global[Y] = 1.0 * (j + offset[Y]);
								r_node_global[Z] = 1.0 * (k + offset[Z]);

								/* Vector from node to particle */
								r_ij[X] = r_particle_global[X] - r_node_global[X];
								r_ij[Y] = r_particle_global[Y] - r_node_global[Y];
								r_ij[Z] = r_particle_global[Z] - r_node_global[Z];

								double dist2 = r_ij[X] * r_ij[X] + r_ij[Y] * r_ij[Y] + r_ij[Z] * r_ij[Z];

								dist = sqrt(r_ij[X] * r_ij[X] + r_ij[Y] * r_ij[Y] + r_ij[Z] * r_ij[Z]);

								/* Skip if too close to avoid singularity */
								if (dist < 0.000001) continue;

								/* Electric field with Debye-Hückel screening (Yukawa potential) */
								/* E = prefactor * q * exp(-κr) * (1/r² + κ/r) * r_unit */
								/* For κ=0, this reduces to standard Coulomb */
								// double exp_kr = (kappa > 0.0) ? exp(-kappa * dist) : 1.0;
								// double screening_factor = exp_kr * (1.0 / dist + kappa) / dist;

								// E_coulomb[X] = prefactor * rho_net * r_ij[X] * screening_factor;
								// E_coulomb[Y] = prefactor * rho_net * r_ij[Y] * screening_factor;
								// E_coulomb[Z] = prefactor * rho_net * r_ij[Z] * screening_factor;

								// // Fuerza con Coulomb porque esta fuera de equilibrio el sistema
								double factor = (1.0 / dist2);
								for (int n = 1; n < 10; n++)
								{
									factor += 2.0 / (2 * (n * dist) * (n * dist) + dist2);
								}

								E_coulomb[X] = prefactor * rho_net * r_ij[X] * factor;
								E_coulomb[Y] = prefactor * rho_net * r_ij[Y] * factor;
								E_coulomb[Z] = prefactor * rho_net * r_ij[Z] * factor;

								/* Accumulate field contribution using Klein summation */
								klein_add_double(&E_total_k[X], E_coulomb[X]);
								klein_add_double(&E_total_k[Y], E_coulomb[Y]);
								klein_add_double(&E_total_k[Z], E_coulomb[Z]);

							}
						}
					}

					/* Store total field on particle (Esub) */
					pc->Esub[X] = klein_sum(&E_total_k[X]);
					pc->Esub[Y] = klein_sum(&E_total_k[Y]);
					pc->Esub[Z] = klein_sum(&E_total_k[Z]);

					/* Total force on particle: F = q * E */
					double force_particle[3];
					force_particle[X] = q_particle * pc->Esub[X];
					force_particle[Y] = q_particle * pc->Esub[Y];
					force_particle[Z] = q_particle * pc->Esub[Z];

					/* Apply force to particle */
					pc->fex[X] += force_particle[X];
					pc->fex[Y] += force_particle[Y];
					pc->fex[Z] += force_particle[Z];

					// /* Distribute reaction force to ALL fluid nodes using Yukawa weight */
					// /* This ensures momentum conservation by using the same weight */
					// /* that was used to calculate the force on the particle */
					// for (i = 1; i <= nlocal[X]; i++) {
					// 	for (j = 1; j <= nlocal[Y]; j++) {
					// 		for (k = 1; k <= nlocal[Z]; k++) {

					// 			index = cs_index(cinfo->cs, i, j, k);

					// 			/* Skip solid sites */
					// 			int status_hydro;
					// 			map_status(map, index, &status_hydro);
					// 			if (status_hydro != MAP_FLUID) continue;

					// 			/* Get charge density at this node */
					// 			double rho0_hydro, rho1_hydro, rho_net_hydro;
					// 			psi_rho(psi, index, 0, &rho0_hydro);
					// 			psi_rho(psi, index, 1, &rho1_hydro);
					// 			rho_net_hydro = rho0_hydro - rho1_hydro;

					// 			/* Node position in global coordinates */
					// 			double r_node_hydro[3];
					// 			r_node_hydro[X] = 1.0 * (i + offset[X]);
					// 			r_node_hydro[Y] = 1.0 * (j + offset[Y]);
					// 			r_node_hydro[Z] = 1.0 * (k + offset[Z]);

					// 			/* Vector from node to particle */
					// 			double r_ij_hydro[3];
					// 			r_ij_hydro[X] = r_particle_global[X] - r_node_hydro[X];
					// 			r_ij_hydro[Y] = r_particle_global[Y] - r_node_hydro[Y];
					// 			r_ij_hydro[Z] = r_particle_global[Z] - r_node_hydro[Z];

					// 			double dist_hydro = sqrt(r_ij_hydro[X]*r_ij_hydro[X] +
					// 			                        r_ij_hydro[Y]*r_ij_hydro[Y] +
					// 			                        r_ij_hydro[Z]*r_ij_hydro[Z]);

					// 			/* Skip if too close to avoid singularity */
					// 			if (dist_hydro < 0.0001) continue;

					// 			/* Calculate reaction force on this node */
					// 			/* F_node = -q_particle * (field contribution from this node) */
					// 			/* field = prefactor * rho_net * exp(-κr) * (1/r² + κ/r) * r̂ */
					// 			double exp_kr_hydro = (kappa > 0.0) ? exp(-kappa * dist_hydro) : 1.0;
					// 			double screening_factor_hydro = exp_kr_hydro * (1.0/dist_hydro + kappa) / dist_hydro;

					// 			double force_hydro[3];
					// 			force_hydro[X] = -q_particle * prefactor * rho_net_hydro * r_ij_hydro[X] * screening_factor_hydro;
					// 			force_hydro[Y] = -q_particle * prefactor * rho_net_hydro * r_ij_hydro[Y] * screening_factor_hydro;
					// 			force_hydro[Z] = -q_particle * prefactor * rho_net_hydro * r_ij_hydro[Z] * screening_factor_hydro;

					// 			add_force_to_array(&force_k_indexed, index, force_hydro);

					// 		}
					// 	}
					// }

				}
			}
		}
	}

	colloid_sums_halo(cinfo, COLLOID_SUM_FORCE_EXT_ONLY);

	/*CHANGE INIT - Write electric field data to file */
	/* Write particle data after colloid_sums_halo, only for local (non-halo) particles */
	/* Loop only over non-halo cells to write data */
	if (fp != NULL) {
		for (ic = 1; ic <= ncell[X]; ic++) {
			for (jc = 1; jc <= ncell[Y]; jc++) {
				for (kc = 1; kc <= ncell[Z]; kc++) {

					colloids_info_cell_list_head(cinfo, ic, jc, kc, &pc);

					for (; pc; pc = pc->next) {

						if (pc->s.bc != COLLOID_BC_SUBGRID) continue;

						/* Write data for each particle */
						char string[256];
						double Emod = sqrt(pc->Esub[X] * pc->Esub[X] +
										   pc->Esub[Y] * pc->Esub[Y] +
										   pc->Esub[Z] * pc->Esub[Z]);
						sprintf(string, "%d;%d;%.15e;%.15e;%.15e;%.15e\n",
								step, pc->s.index, Emod,
								pc->Esub[X], pc->Esub[Y], pc->Esub[Z]);
						for (i = 0; i < strlen(string); i++)
							if (string[i] == '.') string[i] = ',';
						fprintf(fp, "%s", string);
					}
				}
			}
		}
	}
	/*CHANGE END - Write electric field data to file */

	// /* Apply accumulated electric forces to hydro by index */
	// /* This distributes reaction forces using Yukawa/Debye-Hückel weights */
	// if (force_k_indexed != NULL) {
	// 	for (int idx = 0; idx < (force_k_indexed)->count; idx++) {

	// 		distributed_force_klein_entry_t* entry = (force_k_indexed)->entries[idx];

	// 		double new_force[3];
	// 		new_force[X] = klein_sum(entry->force[X]);
	// 		new_force[Y] = klein_sum(entry->force[Y]);
	// 		new_force[Z] = klein_sum(entry->force[Z]);

	// 		hydro_f_local_add(hydro, entry->cs_index, new_force);

	// 	}

	// 	subgrid_free_distributed_force_t(&force_k_indexed);
	// }

	// /*CHANGE INIT - Momentum conservation correction */
	// /* Calculate total force on all particles and apply correction to fluid */
	// {
	// 	int nsfluid;
	// 	double flocal[4] = { 0.0, 0.0, 0.0, 0.0 };  /* Local force sum and fluid count */
	// 	double fsum[4] = { 0.0, 0.0, 0.0, 0.0 };    /* Global force sum and fluid count */
	// 	double force_correction[3];

	// 	/* Sum forces on all local particles */
	// 	for (ic = 0; ic <= ncell[X] + 1; ic++) {
	// 		for (jc = 0; jc <= ncell[Y] + 1; jc++) {
	// 			for (kc = 0; kc <= ncell[Z] + 1; kc++) {

	// 				colloids_info_cell_list_head(cinfo, ic, jc, kc, &pc);

	// 				for (; pc; pc = pc->next) {
	// 					if (pc->s.bc != COLLOID_BC_SUBGRID) continue;

	// 					flocal[X] += pc->fex[X];
	// 					flocal[Y] += pc->fex[Y];
	// 					flocal[Z] += pc->fex[Z];
	// 				}
	// 			}
	// 		}
	// 	}

	// 	/* Count local fluid nodes */
	// 	map_volume_local(map, MAP_FLUID, &nsfluid);
	// 	flocal[3] = (double)nsfluid;

	// 	/* Sum across all MPI ranks */
	// 	MPI_Allreduce(flocal, fsum, 4, MPI_DOUBLE, MPI_SUM, comm);

	// 	/* Calculate average force per fluid node */
	// 	if (fsum[3] > 0.0) {
	// 		force_correction[X] = -fsum[X] / fsum[3];
	// 		force_correction[Y] = -fsum[Y] / fsum[3];
	// 		force_correction[Z] = -fsum[Z] / fsum[3];

	// 		/* Apply correction force to all fluid nodes */
	// 		for (i = 1; i <= nlocal[X]; i++) {
	// 			for (j = 1; j <= nlocal[Y]; j++) {
	// 				for (k = 1; k <= nlocal[Z]; k++) {

	// 					index = cs_index(cinfo->cs, i, j, k);

	// 					/* Check if this is a fluid node */
	// 					int status;
	// 					map_status(map, index, &status);
	// 					if (status != MAP_FLUID) continue;

	// 					/* Apply negative of average force to conserve momentum */
	// 					hydro_f_local_add(hydro, index, force_correction);
	// 				}
	// 			}
	// 		}
	// 	}
	// }
	// /*CHANGE END - Momentum conservation correction */

	hydro_memcpy(hydro, tdpMemcpyHostToDevice);

	return 0;

}
/*CHANGE END - Poisson-Boltzmann force calculation */

/*****************************************************************************
 *
 *  subgrid_compute_self_field_single_particle
 *
 *  Calculates the self-field (autocampo) for a single particle using PETSc
 *  Poisson solver with Peskin charge distribution
 *
 *****************************************************************************/
int subgrid_compute_self_field_single_particle(colloid_t* pc,
											   colloids_info_t* cinfo,
											   psi_t* psi_global,
											   double E_self[3]) {

	int i, j, k, ia, ja, ka, i_min, i_max, j_min, j_max, k_min, k_max;
	int nlocal[3], offset[3];
	double r_i[3], r_j[3], r0[3], r_ij[3];
	double dr_i, dr_j, q_i;
	double dist, dist3;
	double epsilon, eunit, beta, kt;
	double prefactor;
	PI_DOUBLE(pi);
	klein_t E_k[3];

	cs_t* cs = cinfo->cs;

	assert(pc);
	assert(cinfo);
	assert(psi_global);
	assert(E_self);

	cs_nlocal(cs, nlocal);
	cs_nlocal_offset(cs, offset);

	psi_epsilon(psi_global, &epsilon);
	psi_unit_charge(psi_global, &eunit);
	psi_beta(psi_global, &beta);
	kt = 1.0 / beta;

	/* Prefactor for Coulomb force: (eunit * kt) / (4 * pi * epsilon) */
	prefactor = eunit * kt / (4.0 * pi * epsilon);

	/* Initialize Klein sums */
	E_k[X] = klein_zero();
	E_k[Y] = klein_zero();
	E_k[Z] = klein_zero();

	/* Particle position in local coordinates */
	r0[X] = pc->s.r[X] - 1.0 * offset[X];
	r0[Y] = pc->s.r[Y] - 1.0 * offset[Y];
	r0[Z] = pc->s.r[Z] - 1.0 * offset[Z];

	/* Range of lattice sites where charge is distributed */
	i_min = imax(1, (int)floor(r0[X] - drange_));
	i_max = imin(nlocal[X], (int)ceil(r0[X] + drange_));
	j_min = imax(1, (int)floor(r0[Y] - drange_));
	j_max = imin(nlocal[Y], (int)ceil(r0[Y] + drange_));
	k_min = imax(1, (int)floor(r0[Z] - drange_));
	k_max = imin(nlocal[Z], (int)ceil(r0[Z] + drange_));

	/* Double loop: for each pair of lattice sites (i,j) with distributed charge */
	/* compute the Coulomb interaction from site i acting on site j */
	/* then interpolate back to particle center */

	for (i = i_min; i <= i_max; i++) {
		for (j = j_min; j <= j_max; j++) {
			for (k = k_min; k <= k_max; k++) {

				/* Position of source site i */
				r_i[X] = r0[X] - 1.0 * i;
				r_i[Y] = r0[Y] - 1.0 * j;
				r_i[Z] = r0[Z] - 1.0 * k;

				/* Peskin weight and charge at site i */
				dr_i = d_peskin(r_i[X]) * d_peskin(r_i[Y]) * d_peskin(r_i[Z]);
				q_i = (pc->s.q0 - pc->s.q1) * dr_i;

				if (fabs(q_i) < 1e-15) continue; /* Skip if no charge */

				/* For each source site, compute field at particle center */
				/* from ALL charged sites (including itself if dist > 0) */

				for (ia = i_min; ia <= i_max; ia++) {
					for (ja = j_min; ja <= j_max; ja++) {
						for (ka = k_min; ka <= k_max; ka++) {

							/* Position of field evaluation site j relative to particle */
							r_j[X] = r0[X] - 1.0 * ia;
							r_j[Y] = r0[Y] - 1.0 * ja;
							r_j[Z] = r0[Z] - 1.0 * ka;

							/* Peskin weight at site j for interpolation back */
							dr_j = d_peskin(r_j[X]) * d_peskin(r_j[Y]) * d_peskin(r_j[Z]);

							/* Distance vector from site i to site j */
							r_ij[X] = (1.0 * ia - 1.0 * i);
							r_ij[Y] = (1.0 * ja - 1.0 * j);
							r_ij[Z] = (1.0 * ka - 1.0 * k);

							dist = sqrt(r_ij[X] * r_ij[X] + r_ij[Y] * r_ij[Y] + r_ij[Z] * r_ij[Z]);

							/* Skip self-interaction at same site */
							if (dist < 0.1) continue;

							dist3 = dist * dist * dist;

							/* Coulomb field from charge at i evaluated at j */
							/* E = prefactor * q_i * r_ij / |r_ij|^3 */
							/* Then interpolate with weight dr_j */

							double E_contrib[3];
							E_contrib[X] = prefactor * q_i * r_ij[X] / dist3 * dr_j;
							E_contrib[Y] = prefactor * q_i * r_ij[Y] / dist3 * dr_j;
							E_contrib[Z] = prefactor * q_i * r_ij[Z] / dist3 * dr_j;

							klein_add_double(&E_k[X], E_contrib[X]);
							klein_add_double(&E_k[Y], E_contrib[Y]);
							klein_add_double(&E_k[Z], E_contrib[Z]);
						}
					}
				}
			}
		}
	}

	/* Extract final self-field values */
	E_self[X] = klein_sum(&E_k[X]);
	E_self[Y] = klein_sum(&E_k[Y]);
	E_self[Z] = klein_sum(&E_k[Z]);

	return 0;

	return 0;
}

/*****************************************************************************
 *
 *  subgrid_get_lattice_index
 *
 *  Get indexes for neigbour lattice sites
 *
 *****************************************************************************/
void subgrid_get_lattice_index(double r0[3], int nlocal[3], int* i_min, int* i_max, int* j_min, int* j_max, int* k_min, int* k_max)
{
	// Rangos de Peskin - 4x4x4
	*i_min = imax(1, (int)floor(r0[X] - drange_));
	*i_max = imin(nlocal[X], (int)ceil(r0[X] + drange_));
	*j_min = imax(1, (int)floor(r0[Y] - drange_));
	*j_max = imin(nlocal[Y], (int)ceil(r0[Y] + drange_));
	*k_min = imax(1, (int)floor(r0[Z] - drange_));
	*k_max = imin(nlocal[Z], (int)ceil(r0[Z] + drange_));

}

/*****************************************************************************
 *
 *  subgrid_get_lattice_index_fn
 *
 *  Get indexes for firts neigbour lattice sites
 *
 *****************************************************************************/
void subgrid_get_lattice_index_fn(double r0[3], int nlocal[3], int* i_min, int* i_max, int* j_min, int* j_max, int* k_min, int* k_max)
{
	// Rangos de Peskin - 4x4x4
	*i_min = imax(1, (int)floor(r0[X]));
	*i_max = imin(nlocal[X], (int)ceil(r0[X]));
	*j_min = imax(1, (int)floor(r0[Y]));
	*j_max = imin(nlocal[Y], (int)ceil(r0[Y]));
	*k_min = imax(1, (int)floor(r0[Z]));
	*k_max = imin(nlocal[Z], (int)ceil(r0[Z]));

}

/*CHANGE INIT - 20260117 Short-range corrections Level 1: Self-field
 * See docs/SHORT_RANGE_CORRECTIONS_ANALYSIS.md Section 7, Level 1.
 *
 * This corrects the self-field: the electric field that a particle
 * generates on itself through its distributed charge on the mesh.
 * Without this correction, the particle "feels" its own charge.
 */

 /*****************************************************************************
  *
  *  subgrid_shortrange_set_level
  *
  *  Set the short-range correction level.
  *  0 = none, 1 = self-field, 2 = +PP, 3 = +PN, etc.
  *
  *****************************************************************************/
void subgrid_shortrange_set_level(int level) {
	if (level >= SHORTRANGE_LEVEL_NONE && level <= SHORTRANGE_LEVEL_5_CONSERV) {
		shortrange_level_ = (shortrange_level_t)level;
	}
}

/*****************************************************************************
 *
 *  subgrid_shortrange_get_level
 *
 *  Get the current short-range correction level.
 *
 *****************************************************************************/
int subgrid_shortrange_get_level(void) {
	return (int)shortrange_level_;
}

/*****************************************************************************
 *
 *  subgrid_shortrange_level1
 *
 *  Level 1 correction: Self-field only (Alternative A)
 *
 *  DISABLED 20260118: This approach does not conserve momentum because:
 *  - The mesh field (used by Nernst-Planck for fluid forces) is not modified
 *  - Only the particle force is corrected
 *  - The fluid sees the original field with self-field, particle sees corrected
 *  - This creates a force imbalance
 *
 *  For momentum conservation, use Level 3 (PN correction) instead.
 *
 *****************************************************************************/
int subgrid_shortrange_level1(colloids_info_t* cinfo, psi_t* psi, hydro_t* hydro) {

	/* DISABLED - does not conserve momentum
	 * The self-field correction cannot be applied without also modifying
	 * the mesh field that Nernst-Planck uses. Since modifying the mesh
	 * is expensive, we skip Level 1 and use Level 3 (PN) instead.
	 */
	(void)cinfo;
	(void)psi;
	(void)hydro;

	return 0;
}

/*****************************************************************************
 *
 *  subgrid_shortrange_level3_pn
 *
 *  Level 3 correction: Particle-Node (PN) short-range forces
 *
 *  CHANGE 20260118: New approach for momentum conservation.
 *
 *  Instead of adding force to fex (which gets distributed by Peskin later),
 *  we now:
 *  1. Calculate the electric field from each node charge to the particle
 *  2. Add that field to pc->Esub (so particle force comes through normal mechanism)
 *  3. Apply the reaction force directly to the fluid node (Newton III)
 *
 *  This ensures momentum conservation because:
 *  - Particle feels F = q_p * E_correction via Esub
 *  - Fluid node feels -F directly via hydro_f_local_add
 *
 *  The short-range kernel is erf(kappa*r)/r^2, which captures the part
 *  of the Coulomb interaction that the mesh doesn't resolve well.
 *
 *  Arguments:
 *    cinfo - colloid info structure
 *    psi   - electrokinetics psi structure (for charge density)
 *    hydro - hydrodynamics structure (for applying forces to fluid)
 *
 *****************************************************************************/
int subgrid_shortrange_level3_pn(colloids_info_t* cinfo, psi_t* psi, hydro_t* hydro) {

	double prefactor;
	double kappa_E = SHORTRANGE_KAPPA_E;
	double kt, eunit;
	int ncell[3];
	int nlocal[3], offset[3];
	int ic, jc, kc;
	int di, dj, dk;
	int i0, j0, k0;
	int index;
	colloid_t* pc;

	assert(cinfo);
	assert(psi);
	assert(hydro);

	if (cinfo->nsubgrid == 0) return 0;

	prefactor = shortrange_coulomb_prefactor(psi);

	/* Get kT and eunit for field/force conversion */
	psi_unit_charge(psi, &eunit);
	psi_beta(psi, &kt);
	kt = 1.0 / kt;  /* Convert beta to kT */

	cs_nlocal(cinfo->cs, nlocal);
	cs_nlocal_offset(cinfo->cs, offset);
	colloids_info_ncell(cinfo, ncell);

	/* Loop over all particles */
	for (ic = 0; ic <= ncell[X] + 1; ic++) {
		for (jc = 0; jc <= ncell[Y] + 1; jc++) {
			for (kc = 0; kc <= ncell[Z] + 1; kc++) {

				colloids_info_cell_list_head(cinfo, ic, jc, kc, &pc);

				for (; pc; pc = pc->next) {

					if (pc->s.bc != COLLOID_BC_SUBGRID) continue;

					/* Net charge of this particle */
					double q_p = pc->s.q0 - pc->s.q1;

					if (fabs(q_p) < 1.0e-12) continue;  /* Skip uncharged particles */

					/* Klein sums for field correction at particle */
					klein_t E_corr_k[3];
					E_corr_k[X] = klein_zero();
					E_corr_k[Y] = klein_zero();
					E_corr_k[Z] = klein_zero();

					/* Particle position in local coordinates */
					double r0[3];
					r0[X] = pc->s.r[X] - 1.0 * offset[X];
					r0[Y] = pc->s.r[Y] - 1.0 * offset[Y];
					r0[Z] = pc->s.r[Z] - 1.0 * offset[Z];

					/* Base index for short-range support */
					i0 = (int)floor(r0[X]);
					j0 = (int)floor(r0[Y]);
					k0 = (int)floor(r0[Z]);

					/* Loop over nearby nodes within cutoff */
					/* Using same range as Peskin for simplicity (-1 to +2) */
					for (di = -1; di <= 2; di++) {
						for (dj = -1; dj <= 2; dj++) {
							for (dk = -1; dk <= 2; dk++) {

								int node_i = i0 + di;
								int node_j = j0 + dj;
								int node_k = k0 + dk;

								/* Check bounds */
								if (node_i < 1 || node_i > nlocal[X]) continue;
								if (node_j < 1 || node_j > nlocal[Y]) continue;
								if (node_k < 1 || node_k > nlocal[Z]) continue;

								/* Separation vector: node -> particle (for field direction) */
								double r_np[3];
								r_np[X] = r0[X] - (double)node_i;
								r_np[Y] = r0[Y] - (double)node_j;
								r_np[Z] = r0[Z] - (double)node_k;

								double dist = sqrt(r_np[X] * r_np[X] +
												   r_np[Y] * r_np[Y] +
												   r_np[Z] * r_np[Z]);

								/* Skip if too close (singularity) or beyond cutoff */
								if (dist < 0.01) continue;
								if (dist > SHORTRANGE_RC) continue;

								/* Get charge density at this node */
								index = cs_index(cinfo->cs, node_i, node_j, node_k);
								double rho0, rho1;
								psi_rho(psi, index, 0, &rho0);
								psi_rho(psi, index, 1, &rho1);
								double q_node = rho0 - rho1;  /* Net charge at node */

								/* Electric field at particle position from node charge:
								 * E = prefactor * q_node * erf(kappa*r)/r^2 * r_hat
								 *
								 * r_np/dist is unit vector from node to particle
								 * Field points away from positive charge (outward)
								 */
								double kernel = shortrange_erf_kernel(dist, kappa_E);
								double E_mag = prefactor * q_node * kernel;

								double E_from_node[3];
								E_from_node[X] = E_mag * r_np[X] / dist;
								E_from_node[Y] = E_mag * r_np[Y] / dist;
								E_from_node[Z] = E_mag * r_np[Z] / dist;

								/* Accumulate field correction at particle */
								klein_add_double(&E_corr_k[X], E_from_node[X]);
								klein_add_double(&E_corr_k[Y], E_from_node[Y]);
								klein_add_double(&E_corr_k[Z], E_from_node[Z]);

								/* Force on particle from this field: F_p = q_p * E
								 * Newton III reaction on fluid node: F_node = -F_p = -q_p * E
								 *
								 * Convert to lattice force units for fluid
								 */
								double F_on_node[3];
								F_on_node[X] = -q_p * E_from_node[X] * kt / eunit;
								F_on_node[Y] = -q_p * E_from_node[Y] * kt / eunit;
								F_on_node[Z] = -q_p * E_from_node[Z] * kt / eunit;

								/* Apply reaction force to fluid node */
								hydro_f_local_add(hydro, index, F_on_node);
							}
						}
					}

					/* Add field correction to Esub
					 * The particle force will come from: F = q_p * (Esub_mesh + E_correction)
					 * through the normal mechanism in subgrid_on_velocity
					 */
					pc->Esub[X] += klein_sum(&E_corr_k[X]);
					pc->Esub[Y] += klein_sum(&E_corr_k[Y]);
					pc->Esub[Z] += klein_sum(&E_corr_k[Z]);
				}
			}
		}
	}

	return 0;
}

/*****************************************************************************
 *
 *  subgrid_shortrange_corrections
 *
 *  Main entry point for short-range corrections.
 *  Calls the appropriate level based on shortrange_level_.
 *
 *  CHANGE 20260118: Level 1 disabled (doesn't conserve momentum).
 *  Level 3 (PN) now implements particle-node short-range forces with
 *  proper Newton III pairs to conserve momentum.
 *
 *****************************************************************************/
int subgrid_shortrange_corrections(colloids_info_t* cinfo, psi_t* psi,
									hydro_t* hydro) {

	switch (shortrange_level_) {
	case SHORTRANGE_LEVEL_NONE:
		return 0;

	case SHORTRANGE_LEVEL_1_SELF:
		/* Level 1 disabled - doesn't conserve momentum */
		return subgrid_shortrange_level1(cinfo, psi, hydro);  /* Returns 0, does nothing */

	case SHORTRANGE_LEVEL_2_PP:
		/* Level 2: PP only (to be implemented) */
		return 0;  /* Not implemented yet */

	case SHORTRANGE_LEVEL_3_PN:
		/* Level 3: Particle-Node short-range forces */
		return subgrid_shortrange_level3_pn(cinfo, psi, hydro);

	case SHORTRANGE_LEVEL_4_P3M:
		/* Level 4: P3M with reference function (to be implemented) */
		return subgrid_shortrange_level3_pn(cinfo, psi, hydro);  /* Fallback to Level 3 */

	case SHORTRANGE_LEVEL_5_CONSERV:
		/* Level 5: Full conservation (to be implemented) */
		return subgrid_shortrange_level3_pn(cinfo, psi, hydro);  /* Fallback to Level 3 */

	default:
		return -1;
	}
}

/*CHANGE END - 20260117 Short-range corrections Level 1 */

/*CHANGE END - Subgrid charge */


/*****************************************************************************
 *
 *  subgrid_update
 *
 *  This function is responsible for update of position for
 *  sub-gridscale particles. It takes the place of BBL for
 *  fully resolved particles.
 *
 *****************************************************************************/
int subgrid_update(colloids_info_t* cinfo, hydro_t* hydro, int noise_flag) {
	int ia;
	int ic, jc, kc;
	int ncell[3];
	double drag, reta;
	double eta;
	PI_DOUBLE(pi);
	colloid_t* p_colloid;
	physics_t* phys = NULL;

	double ran[2];    /* Random numbers for fluctuation dissipation correction */
	double frand[3];  /* Random force */
	double kt;        /* Temperature */

	assert(cinfo);
	assert(hydro);

	if (cinfo->nsubgrid == 0) return 0;

	colloids_info_ncell(cinfo, ncell);

	subgrid_interpolation(cinfo, hydro);
	colloid_sums_halo(cinfo, COLLOID_SUM_SUBGRID);

	/* Loop through all cells (including the halo cells) */

	physics_ref(&phys);
	physics_eta_shear(phys, &eta);
	physics_kt(phys, &kt);
	reta = 1.0 / (6.0 * pi * eta);

	for (ic = 0; ic <= ncell[X] + 1; ic++) {
		for (jc = 0; jc <= ncell[Y] + 1; jc++) {
			for (kc = 0; kc <= ncell[Z] + 1; kc++) {

				colloids_info_cell_list_head(cinfo, ic, jc, kc, &p_colloid);

				for (; p_colloid; p_colloid = p_colloid->next) {

					if (p_colloid->s.bc != COLLOID_BC_SUBGRID) continue;

					drag = reta * (1.0 / p_colloid->s.ah - 1.0 / p_colloid->s.al);

					if (noise_flag == 0) {
						frand[X] = 0.0; frand[Y] = 0.0; frand[Z] = 0.0;
					}
					else {
						for (ia = 0; ia < 3; ia++) {
							while (1) {
								/* To keep the random correction smaller than 3 sigma.
								 * Otherwise, a large thermal fluctuation may cause a
								 * numerical problem. */
								util_ranlcg_reap_gaussian(&p_colloid->s.rng, ran);
								if (fabs(ran[0]) < 3.0) {
									frand[ia] = sqrt(2.0 * kt * drag) * ran[0];
									break;
								}
								if (fabs(ran[1]) < 3.0) {
									frand[ia] = sqrt(2.0 * kt * drag) * ran[1];
									break;
								}
							}
						}
					}

					for (ia = 0; ia < 3; ia++) {
						p_colloid->s.v[ia] = p_colloid->fsub[ia] + drag * p_colloid->fex[ia]
							+ frand[ia];
						p_colloid->s.dr[ia] = p_colloid->s.v[ia];
					}
				}
				/* Next cell */
			}
		}
	}

	return 0;
}

/*****************************************************************************
 *
 *  subgrid_interpolation
 *
 *  Interpolate (delta function method) the lattice velocity field
 *  to the position of the particles.
 *
 *****************************************************************************/

static int subgrid_interpolation(colloids_info_t* cinfo, hydro_t* hydro) {

	int ic, jc, kc;
	int i, j, k, i_min, i_max, j_min, j_max, k_min, k_max;
	int index;
	int nlocal[3], offset[3];
	int ncell[3];

	double r0[3], r[3], u[3];
	double dr;
	colloid_t* p_colloid;

	assert(cinfo);
	assert(hydro);

	cs_nlocal(cinfo->cs, nlocal);
	cs_nlocal_offset(cinfo->cs, offset);
	colloids_info_ncell(cinfo, ncell);

	/* While there is no subgrid device implementation,
	   need to recover the current velocity. */

	hydro_memcpy(hydro, tdpMemcpyDeviceToHost);

	//CHANGE INIT - Subgrid charge
	/* Loop through all cells (including the halo cells) and set
	 * the velocity at each particle to zero for this step. */

	 // for (ic = 0; ic <= ncell[X] + 1; ic++) {
	 // 	for (jc = 0; jc <= ncell[Y] + 1; jc++) {
	 // 		for (kc = 0; kc <= ncell[Z] + 1; kc++) {

	 // 			colloids_info_cell_list_head(cinfo, ic, jc, kc, &p_colloid);

	 // 			for (; p_colloid; p_colloid = p_colloid->next) {

	 // 				if (p_colloid->s.bc != COLLOID_BC_SUBGRID) continue;

	 // 				p_colloid->fsub[X] = 0.0;
	 // 				p_colloid->fsub[Y] = 0.0;
	 // 				p_colloid->fsub[Z] = 0.0;
	 // 			}
	 // 		}
	 // 	}
	 // }
	 //CHANGE END - Subgrid charge

	 /* And add up the contributions to the velocity from the lattice. */
	for (ic = 0; ic <= ncell[X] + 1; ic++) {
		for (jc = 0; jc <= ncell[Y] + 1; jc++) {
			for (kc = 0; kc <= ncell[Z] + 1; kc++) {

				colloids_info_cell_list_head(cinfo, ic, jc, kc, &p_colloid);

				for (; p_colloid; p_colloid = p_colloid->next) {

					if (p_colloid->s.bc != COLLOID_BC_SUBGRID) continue;

					//CHANGE INIT - Subgrid charge
					klein_t force_k[3];
					force_k[X] = klein_zero();
					force_k[Y] = klein_zero();
					force_k[Z] = klein_zero();
					//CHANGE END - Subgrid charge

					/* Need to translate the colloid position to "local"
					* coordinates, so that the correct range of lattice
					* nodes is found */

					r0[X] = p_colloid->s.r[X] - 1.0 * offset[X];
					r0[Y] = p_colloid->s.r[Y] - 1.0 * offset[Y];
					r0[Z] = p_colloid->s.r[Z] - 1.0 * offset[Z];

					/* Work out which local lattice sites are involved
					 * and loop around */

					i_min = imax(1, (int)floor(r0[X] - drange_));
					i_max = imin(nlocal[X], (int)ceil(r0[X] + drange_));
					j_min = imax(1, (int)floor(r0[Y] - drange_));
					j_max = imin(nlocal[Y], (int)ceil(r0[Y] + drange_));
					k_min = imax(1, (int)floor(r0[Z] - drange_));
					k_max = imin(nlocal[Z], (int)ceil(r0[Z] + drange_));

					for (i = i_min; i <= i_max; i++) {
						for (j = j_min; j <= j_max; j++) {
							for (k = k_min; k <= k_max; k++) {

								index = cs_index(cinfo->cs, i, j, k);

								/* Separation between r0 and the coordinate position of
								 * this site */

								r[X] = r0[X] - 1.0 * i;
								r[Y] = r0[Y] - 1.0 * j;
								r[Z] = r0[Z] - 1.0 * k;

								dr = d_peskin(r[X]) * d_peskin(r[Y]) * d_peskin(r[Z]);
								hydro_u(hydro, index, u);

								//CHANGE INIT - Subgrid charge
								// p_colloid->fsub[X] += u[X] * dr;
								// p_colloid->fsub[Y] += u[Y] * dr;
								// p_colloid->fsub[Z] += u[Z] * dr;

								/* Add to particle force */
								klein_add_double(&force_k[X], u[X] * dr);
								klein_add_double(&force_k[Y], u[Y] * dr);
								klein_add_double(&force_k[Z], u[Z] * dr);
								//CHANGE END - Subgrid charge

							}
						}
					}

					//CHANGE INIT - Subgrid charge
					p_colloid->fsub[X] = klein_sum(&force_k[X]);
					p_colloid->fsub[Y] = klein_sum(&force_k[Y]);
					p_colloid->fsub[Z] = klein_sum(&force_k[Z]);
					//CHANGE END - Subgrid charge

		/* Next colloid */
				}
				/* Next cell */
			}
		}
	}

	return 0;
}

/*****************************************************************************
 *
 *  subgrid_wall_lubrication
 *
 *  Accumulate lubrication corrections to the external force on each particle.
 *
 *****************************************************************************/

int subgrid_wall_lubrication(colloids_info_t* cinfo, wall_t* wall) {

	double drag[3];
	colloid_t* pc = NULL;

	double f[3] = { 0.0, 0.0, 0.0 };

	assert(cinfo);
	assert(wall);

	colloids_info_local_head(cinfo, &pc);

	for (; pc; pc = pc->nextlocal) {
		if (pc->s.bc != COLLOID_BC_SUBGRID) continue;
		wall_lubr_sphere(wall, pc->s.ah, pc->s.r, drag);
		pc->fex[X] += drag[X] * pc->s.v[X];
		pc->fex[Y] += drag[Y] * pc->s.v[Y];
		pc->fex[Z] += drag[Z] * pc->s.v[Z];
		f[X] -= drag[X] * pc->s.v[X];
		f[Y] -= drag[Y] * pc->s.v[Y];
		f[Z] -= drag[Z] * pc->s.v[Z];
	}

	wall_momentum_add(wall, f);

	return 0;
}

/*****************************************************************************
 *
 *  d_peskin
 *
 *  Approximation to \delta(r) according to Peskin.
 *
 *****************************************************************************/

 // static double d_peskin(double r) {
double d_peskin(double r) {

	double rmod;
	double delta = 0.0;

	/*CHANGE INIT - 20260418 Scale Peskin kernel with drange_ to support variable stencil width.
	 * drange_=1 (original, support=2): s=1, no scaling, identical to original.
	 * drange_=2 (support=4): r mapped to r/2, result divided by 2 to preserve partition-of-unity. */
	// rmod = fabs(r);
	double s = (drange_+1.0)/2.0;  /* Scale factor for mapping r to r/s */
	rmod = fabs(r / s);
	/*CHANGE END - 20260418 Scale Peskin kernel with drange_ */

	if (rmod <= 1.0) {
		delta = 0.125 * (3.0 - 2.0 * rmod + sqrt(1.0 + 4.0 * rmod - 4.0 * rmod * rmod));
	}
	else if (rmod <= 2.0) {
		delta = 0.125 * (5.0 - 2.0 * rmod - sqrt(-7.0 + 12.0 * rmod - 4.0 * rmod * rmod));
	}

	/*CHANGE INIT - 20260418 Scale Peskin kernel with drange_ */
	// return delta;
	return delta / s;
	/*CHANGE END - 20260418 Scale Peskin kernel with drange_ */
}

/*****************************************************************************
 *
 *  d_peskin_derivative
 *
 *  Derivative of the Peskin delta function.
 *
 * dx += d_peskin_derivative(dx) * d_peskin(dy) * d_peskin(dz)
 * dy += d_peskin(dx) * d_peskin_derivative(dy) * d_peskin(dz)
 * dz += d_peskin(dx) * d_peskin(dy) * d_peskin_derivative(dz)
 *****************************************************************************/
double d_peskin_derivative(double x) {

	double r = fabs(x/drange_);
	double sign;
	double val = 0.0;

	if (x > 0.0) sign = 1.0;
	else if (x < 0.0) sign = -1.0;
	else sign = 0.0;

	if (r <= 1.0) {
		double tmp = sqrt(1.0 + 4.0 * r - 4.0 * r * r);
		val = 0.125 * (-2.0 + (4.0 - 8.0 * r) / (2.0 * tmp));
	}
	else if (r <= 2.0) {
		double tmp = sqrt(-7.0 + 12.0 * r - 4.0 * r * r);
		val = 0.125 * (-2.0 - (12.0 - 8.0 * r) / (2.0 * tmp));
	}

	return val * sign;
}

/*****************************************************************************
 *
 *  d_trilinear
 *
 *  Trilinear interpolation weight function for nearest neighbors.
 *  Uses linear interpolation in each dimension, which when combined
 *  gives trilinear interpolation over the 8 nearest lattice nodes.
 *
 *  For a given 1D distance r from a lattice node:
 *    - If |r| <= 1.0: weight = 1.0 - |r| (linear decay)
 *    - If |r| > 1.0:  weight = 0.0 (no contribution)
 *
 *  When used as: dr = d_trilinear(r[X]) * d_trilinear(r[Y]) * d_trilinear(r[Z])
 *  this produces trilinear interpolation weights for the 8 nearest neighbors.
 *
 *  Example: If particle is at (i+0.3, j+0.7, k+0.2):
 *    - Node (i,j,k):     weight = (1-0.3) * (1-0.7) * (1-0.2) = 0.168
 *    - Node (i+1,j,k):   weight = 0.3 * (1-0.7) * (1-0.2) = 0.072
 *    - Node (i,j+1,k):   weight = (1-0.3) * 0.7 * (1-0.2) = 0.392
 *    - etc. for all 8 neighbors
 *
 *  Sum of all 8 weights = 1.0 (exact conservation)
 *
 *****************************************************************************/

double d_trilinear(double r) {

	double rmod;
	double weight = 0.0;

	rmod = fabs(r);

	if (rmod <= 1.0) {
		weight = 1.0 - rmod;
	}

	return weight;
}

/*****************************************************************************
 *
 *  d_idw
 *
 *  Linear inverse distance weighting interpolation function.
 *  Distributes quantities from a particle position to neighboring nodes
 *  proportional to 1/distance, normalized by the sum of all 1/distance.
 *
 *  This function must be called for ALL neighboring nodes to compute
 *  the proper normalization. It takes:
 *    - r0[3]: particle position in local coordinates
 *    - node_i, node_j, node_k: indices of the current node
 *    - i_min, i_max, j_min, j_max, k_min, k_max: range of all nodes
 *
 *  Returns the weight for this specific node.
 *
 *****************************************************************************/

double d_idw(double r0[3], int node_i, int node_j, int node_k,
					int i_min, int i_max, int j_min, int j_max,
					int k_min, int k_max) {

	int i, j, k;
	double r_current[3], r_node[3];
	double dist_current, dist_node;
	double inv_dist_sum = 0.0;
	double inv_dist_current;
	double weight = 0.0;
	const double epsilon = 1e-15;  /* Small value to avoid division by zero */

	/* Position of current node relative to particle */
	r_current[X] = r0[X] - 1.0 * node_i;
	r_current[Y] = r0[Y] - 1.0 * node_j;
	r_current[Z] = r0[Z] - 1.0 * node_k;

	/* Distance from particle to current node */
	dist_current = sqrt(r_current[X] * r_current[X] +
						r_current[Y] * r_current[Y] +
						r_current[Z] * r_current[Z]);

	/* Compute sum of inverse distances to all neighboring nodes */
	for (i = i_min; i <= i_max; i++) {
		for (j = j_min; j <= j_max; j++) {
			for (k = k_min; k <= k_max; k++) {

				r_node[X] = r0[X] - 1.0 * i;
				r_node[Y] = r0[Y] - 1.0 * j;
				r_node[Z] = r0[Z] - 1.0 * k;

				dist_node = sqrt(r_node[X] * r_node[X] +
								r_node[Y] * r_node[Y] +
								r_node[Z] * r_node[Z]);

				/* Add epsilon to avoid division by zero */
				inv_dist_sum += 1.0 / (dist_node + epsilon);

			}
		}
	}

	/* Calculate normalized weight for current node */
	inv_dist_current = 1.0 / (dist_current + epsilon);

	if (inv_dist_sum > epsilon) {
		weight = inv_dist_current / inv_dist_sum;
	}

	return weight;
}
/*****************************************************************************
 *
 *  d_isdw
 *
 *  Linear square inverse distance weighting interpolation function.
 *  Distributes quantities from a particle position to neighboring nodes
 *  proportional to 1/distance, normalized by the sum of all 1/distance.
 *
 *  This function must be called for ALL neighboring nodes to compute
 *  the proper normalization. It takes:
 *    - r0[3]: particle position in local coordinates
 *    - node_i, node_j, node_k: indices of the current node
 *    - i_min, i_max, j_min, j_max, k_min, k_max: range of all nodes
 *
 *  Returns the weight for this specific node.
 *
 *****************************************************************************/

double d_isdw(double r0[3], int node_i, int node_j, int node_k,
					int i_min, int i_max, int j_min, int j_max,
					int k_min, int k_max) {

	int i, j, k;
	double r_current[3], r_node[3];
	double dist_current, dist_node;
	double inv_dist_sum = 0.0;
	double inv_dist_current;
	double weight = 0.0;
	const double epsilon = 1e-16;  /* Small value to avoid division by zero */

	/* Position of current node relative to particle */
	r_current[X] = r0[X] - 1.0 * node_i;
	r_current[Y] = r0[Y] - 1.0 * node_j;
	r_current[Z] = r0[Z] - 1.0 * node_k;

	/* Distance from particle to current node */
	dist_current = sqrt(r_current[X] * r_current[X] +
						r_current[Y] * r_current[Y] +
						r_current[Z] * r_current[Z]);

	/* Compute sum of inverse distances to all neighboring nodes */
	for (i = i_min; i <= i_max; i++) {
		for (j = j_min; j <= j_max; j++) {
			for (k = k_min; k <= k_max; k++) {

				r_node[X] = r0[X] - 1.0 * i;
				r_node[Y] = r0[Y] - 1.0 * j;
				r_node[Z] = r0[Z] - 1.0 * k;

				dist_node = sqrt(r_node[X] * r_node[X] +
								r_node[Y] * r_node[Y] +
								r_node[Z] * r_node[Z]);

				/* Add epsilon to avoid division by zero */
				inv_dist_sum += 1.0 / (dist_node + epsilon) * (dist_node + epsilon);

			}
		}
	}

	/* Calculate normalized weight for current node */
	inv_dist_current = 1.0 / (dist_current + epsilon) * (dist_current + epsilon);

	if (inv_dist_sum > epsilon) {
		weight = inv_dist_current / inv_dist_sum;
	}

	return weight;
}

/*CHANGE INIT - Poisson-Boltzmann weight */
/*****************************************************************************
 *
 *  d_pb()
 *
 *  Poisson-Boltzmann (Debye-Hückel) weight function for charge distribution.
 *
 *  The weight is defined as the screened Coulomb field contribution from
 *  this node divided by the sum of all field contributions:
 *
 *    w_i = E_i / sum_j(E_j)
 *
 *  where E_i = exp(-κr_i) * (1/r_i + κ) / r_i  (Yukawa/Debye-Hückel field)
 *
 *  For κ=0, this reduces to Coulomb: E_i = 1/r_i²
 *
 *  Parameters:
 *    r0      - particle position in local coordinates
 *    node_i,j,k - indices of current node
 *    i_min..k_max - range of nodes for distribution
 *    kappa   - inverse Debye length (1/λ_D), set to 0 for pure Coulomb
 *
 *****************************************************************************/
double d_pb(double r0[3], int node_i, int node_j, int node_k,
			int i_min, int i_max, int j_min, int j_max,
			int k_min, int k_max, double kappa) {

	int i, j, k;
	double r_current[3], r_node[3];
	double dist_current, dist_node;
	double E_sum = 0.0;
	double E_current;
	double weight = 0.0;
	const double epsilon = 1e-16;  /* Small value to avoid division by zero */

	/* Position of current node relative to particle */
	r_current[X] = r0[X] - 1.0 * node_i;
	r_current[Y] = r0[Y] - 1.0 * node_j;
	r_current[Z] = r0[Z] - 1.0 * node_k;

	/* Distance from particle to current node */
	dist_current = sqrt(r_current[X] * r_current[X] +
						r_current[Y] * r_current[Y] +
						r_current[Z] * r_current[Z]);

	/* Compute sum of Debye-Hückel field contributions from all nodes */
	for (i = i_min; i <= i_max; i++) {
		for (j = j_min; j <= j_max; j++) {
			for (k = k_min; k <= k_max; k++) {

				r_node[X] = r0[X] - 1.0 * i;
				r_node[Y] = r0[Y] - 1.0 * j;
				r_node[Z] = r0[Z] - 1.0 * k;

				dist_node = sqrt(r_node[X] * r_node[X] +
								 r_node[Y] * r_node[Y] +
								 r_node[Z] * r_node[Z]);

				/* Debye-Hückel field: E = exp(-κr) * (1/r + κ) / r */
				/* For κ=0: E = 1/r² (Coulomb) */
				if (dist_node > epsilon) {
					double exp_kr = (kappa > 0.0) ? exp(-kappa * dist_node) : 1.0;
					double E_node = exp_kr * (1.0 / dist_node + kappa) / dist_node;
					E_sum += E_node;
				}

			}
		}
	}

	/* Calculate Debye-Hückel field for current node */
	if (dist_current > epsilon) {
		double exp_kr = (kappa > 0.0) ? exp(-kappa * dist_current) : 1.0;
		E_current = exp_kr * (1.0 / dist_current + kappa) / dist_current;
	}
	else {
		/* If particle is exactly at node, give it full weight */
		E_current = 1.0 / epsilon;
	}

	/* Calculate normalized weight */
	if (E_sum > epsilon) {
		weight = E_current / E_sum;
	}

	return weight;
}

/*****************************************************************************
 *
 *  subgrid_add_psi_grad_from_particles
 *
 *  Add the electric potential from subgrid particles (calculated using
 *  Poisson-Boltzmann/Debye-Hückel solution) to the fluid electric potential.
 *  This avoids distributing particle charges to the grid.
 *
 *  The Yukawa/Debye-Hückel screened potential from a charged particle is:
 *    φ(r) = (q * e / 4πε) * exp(-κr) / r
 *
 *  This function:
 *  1. First pass: Calculate particle potential at all nodes and compute mean
 *  2. Use MPI_Allreduce to get global mean across all processes
 *  3. Second pass: Add particle potential (with mean removed) to fluid potential
 *  4. This ensures the particle potential has the same reference (zero mean)
 *     as the PETSc solution for the fluid
 *  5. The electric field E = -∇φ will then be calculated from the modified
 *     potential by the existing psi_electric_field() function
 *
 *****************************************************************************/

int subgrid_add_psi_grad_from_particles(colloids_info_t* cinfo, psi_t* psi,
										 physics_t* phys) {
	int ic, jc, kc;
	int i, j, k;
	int index;
	int nlocal[3], offset[3];
	int ncell[3];
	int ntotal[3];
	double eunit, epsilon;
	double kappa;
	colloid_t* pc = NULL;
	MPI_Comm comm;
	PI_DOUBLE(pi);

	assert(cinfo);
	assert(psi);
	assert(phys);

	if (cinfo->nsubgrid == 0) return 0;

	cs_nlocal(psi->cs, nlocal);
	cs_nlocal_offset(psi->cs, offset);
	cs_ntotal(psi->cs, ntotal);
	colloids_info_ncell(cinfo, ncell);
	cs_cart_comm(psi->cs, &comm);

	/* Get physical parameters */
	psi_unit_charge(psi, &eunit);
	psi_epsilon(psi, &epsilon);

	/* Compute or get stored kappa */
	subgrid_compute_kappa(psi, &kappa);

	/* Coulomb prefactor for potential: eunit / (4*pi*epsilon) */
	double prefactor = eunit / (4.0 * pi * epsilon);

	/* First pass: Calculate particle potential at all nodes and compute mean */
	double phi_sum_local = 0.0;
	int nsites_local = nlocal[X] * nlocal[Y] * nlocal[Z];

	for (i = 1; i <= nlocal[X]; i++) {
		for (j = 1; j <= nlocal[Y]; j++) {
			for (k = 1; k <= nlocal[Z]; k++) {

				/* Global position of this lattice node */
				double r_node[3];
				r_node[X] = 1.0 * (i + offset[X]);
				r_node[Y] = 1.0 * (j + offset[Y]);
				r_node[Z] = 1.0 * (k + offset[Z]);

				/* Electric potential contribution from all subgrid particles */
				double phi_particles = 0.0;

				/* Loop over all particles in all cells */
				for (ic = 0; ic <= ncell[X] + 1; ic++) {
					for (jc = 0; jc <= ncell[Y] + 1; jc++) {
						for (kc = 0; kc <= ncell[Z] + 1; kc++) {

							colloids_info_cell_list_head(cinfo, ic, jc, kc, &pc);

							for (; pc; pc = pc->next) {

								if (pc->s.bc != COLLOID_BC_SUBGRID) continue;

								/* Net charge on particle (accounting for valencies) */
								double q_net = pc->s.q0 - pc->s.q1;

								if (fabs(q_net) < 1.0e-14) continue; /* Skip neutral particles */

								/* Particle position */
								double r_particle[3];
								r_particle[X] = pc->s.r[X];
								r_particle[Y] = pc->s.r[Y];
								r_particle[Z] = pc->s.r[Z];

								/* Vector from particle to node */
								double dr[3];
								cs_minimum_distance(psi->cs, r_node, r_particle, dr);

								double r = sqrt(dr[X] * dr[X] + dr[Y] * dr[Y] + dr[Z] * dr[Z]);

								/* Avoid singularity at particle center */
								if (r < 0.1) continue; /* Skip nodes too close to particle */

								/* Yukawa/Debye-Hückel screened potential:
								 * φ = prefactor * q * exp(-κr) / r */
								double exp_kr = (kappa > 0.0) ? exp(-kappa * r) : 1.0;
								double phi = prefactor * q_net * exp_kr / r;

								/* Accumulate potential contribution from this particle */
								phi_particles += phi;
							}
						}
					}
				}

				/* Accumulate sum for mean calculation */
				phi_sum_local += phi_particles;
			}
		}
	}

	/* Compute global mean of particle potential */
	double phi_sum_global = 0.0;
	int nsites_global = ntotal[X] * ntotal[Y] * ntotal[Z];

	MPI_Allreduce(&phi_sum_local, &phi_sum_global, 1, MPI_DOUBLE, MPI_SUM, comm);
	double phi_mean = phi_sum_global / (double)nsites_global;

	/* Second pass: Add particle potential (with mean removed) to fluid potential */
	for (i = 1; i <= nlocal[X]; i++) {
		for (j = 1; j <= nlocal[Y]; j++) {
			for (k = 1; k <= nlocal[Z]; k++) {

				index = cs_index(psi->cs, i, j, k);

				/* Global position of this lattice node */
				double r_node[3];
				r_node[X] = 1.0 * (i + offset[X]);
				r_node[Y] = 1.0 * (j + offset[Y]);
				r_node[Z] = 1.0 * (k + offset[Z]);

				/* Electric potential contribution from all subgrid particles */
				double phi_particles = 0.0;

				/* Loop over all particles in all cells */
				for (ic = 0; ic <= ncell[X] + 1; ic++) {
					for (jc = 0; jc <= ncell[Y] + 1; jc++) {
						for (kc = 0; kc <= ncell[Z] + 1; kc++) {

							colloids_info_cell_list_head(cinfo, ic, jc, kc, &pc);

							for (; pc; pc = pc->next) {

								if (pc->s.bc != COLLOID_BC_SUBGRID) continue;

								/* Net charge on particle (accounting for valencies) */
								double q_net = pc->s.q0 - pc->s.q1;

								if (fabs(q_net) < 1.0e-14) continue; /* Skip neutral particles */

								/* Particle position */
								double r_particle[3];
								r_particle[X] = pc->s.r[X];
								r_particle[Y] = pc->s.r[Y];
								r_particle[Z] = pc->s.r[Z];

								/* Vector from particle to node */
								double dr[3];
								cs_minimum_distance(psi->cs, r_node, r_particle, dr);

								double r = sqrt(dr[X] * dr[X] + dr[Y] * dr[Y] + dr[Z] * dr[Z]);

								/* Avoid singularity at particle center */
								if (r < 0.1) continue; /* Skip nodes too close to particle */

								/* Yukawa/Debye-Hückel screened potential:
								 * φ = prefactor * q * exp(-κr) / r */
								double exp_kr = (kappa > 0.0) ? exp(-kappa * r) : 1.0;
								double phi = prefactor * q_net * exp_kr / r;

								/* Accumulate potential contribution from this particle */
								phi_particles += phi;
							}
						}
					}
				}

				/* Add particle potential (with mean removed) to fluid potential */
				/* This ensures both have the same reference (zero mean) */
				double phi_current = 0.0;
				psi_psi(psi, index, &phi_current);
				psi_psi_set(psi, index, phi_current + phi_particles - phi_mean);
			}
		}
	}

	return 0;
}

/*****************************************************************************
 *
 *  subgrid_compute_kappa()
 *
 *  Compute the inverse Debye length (kappa) from the psi object.
 *  This should be called once at the beginning of the simulation.
 *
 *  kappa = 1/λ_D where λ_D = 1/sqrt(8π × l_B × rho_bulk)
 *  and l_B = e²/(4πεkT) is the Bjerrum length.
 *
 *  Returns 0 on success, stores kappa in output parameter and static variable.
 *
 *****************************************************************************/
int subgrid_compute_kappa(psi_t* psi, double* kappa) {

	double lambda_debye;       /* Debye length */
	double rho_bulk;           /* Bulk ionic strength */
	double lb;                 /* Bjerrum length */
	double kt, eunit, epsilon, beta;
	double rho0_bulk, rho1_bulk;
	int nlocal[3];
	PI_DOUBLE(pi);

	assert(psi);
	assert(kappa);

	cs_nlocal(psi->cs, nlocal);

	psi_unit_charge(psi, &eunit);
	psi_epsilon(psi, &epsilon);
	psi_beta(psi, &beta);
	kt = 1.0 / beta;  /* beta = 1/kT, so kT = 1/beta */

	/* Bjerrum length: l_B = e^2 / (4*pi*epsilon*kT) */
	lb = (eunit * eunit) / (4.0 * pi * epsilon * kt);

	/* Get bulk ionic density from center of domain */
	psi_rho(psi, cs_index(psi->cs, nlocal[X] / 2, nlocal[Y] / 2, nlocal[Z] / 2), 0, &rho0_bulk);
	psi_rho(psi, cs_index(psi->cs, nlocal[X] / 2, nlocal[Y] / 2, nlocal[Z] / 2), 1, &rho1_bulk);
	rho_bulk = 0.5 * (rho0_bulk + rho1_bulk);

	/* Calculate Debye length and kappa */
	if (rho_bulk > 1e-15) {
		lambda_debye = 1.0 / sqrt(8.0 * pi * lb * rho_bulk);
		*kappa = 1.0 / lambda_debye;
	}
	else {
		/* No screening if no bulk electrolyte (pure Coulomb) */
		*kappa = 0.0;
	}

	/* Store in static variable */
	kappa_pb_ = *kappa;
	kappa_initialized_ = 1;

	return 0;
}

/*****************************************************************************
 *
 *  subgrid_get_kappa()
 *
 *  Return the stored kappa value. Should be called after subgrid_compute_kappa()
 *  or subgrid_set_kappa() has been called.
 *
 *****************************************************************************/
double subgrid_get_kappa(void) {

	return kappa_pb_;
}

/*****************************************************************************
 *
 *  subgrid_set_kappa()
 *
 *  Manually set the kappa value. Useful for testing or when kappa is known
 *  from external sources.
 *
 *****************************************************************************/
void subgrid_set_kappa(double kappa) {

	kappa_pb_ = kappa;
	kappa_initialized_ = 1;
}

/*CHANGE INIT - Peskin scatter kernel on GPU */

/* Flat structure to pass particle data to GPU kernel */
typedef struct {
	double rx, ry, rz;  /* local-domain position (offset already subtracted) */
	double q0, q1;
} peskin_particle_t;

#ifdef __NVCC__
#include <cuda_runtime.h>

__device__ static double d_peskin_gpu(double r) {
	double rmod = fabs(r);
	if (rmod <= 1.0)
		return 0.125 * (3.0 - 2.0 * rmod + sqrt(1.0 + 4.0 * rmod - 4.0 * rmod * rmod));
	if (rmod <= 2.0)
		return 0.125 * (5.0 - 2.0 * rmod - sqrt(-7.0 + 12.0 * rmod - 4.0 * rmod * rmod));
	return 0.0;
}

/*****************************************************************************
 *
 *  peskin_scatter_nodes_kernel
 *
 *  Each thread = one interior LB node from psi_src.
 *  Scatters its charge to the Peskin neighbourhood in rho_dst via atomicAdd.
 *  SOA layout: rho->data[nsites*species + index]
 *  Ludwig index stride (Z fastest):
 *    str_z=1, str_y=nz+2*nhalo, str_x=str_y*(ny+2*nhalo)
 *    index(i,j,k) = str_x*(nhalo+i-1) + str_y*(nhalo+j-1) + (nhalo+k-1)
 *
 *****************************************************************************/
__global__ void peskin_scatter_nodes_kernel(
	const double* __restrict__ rho_src,
	double* rho_dst,
	int nx, int ny, int nz, int nhalo, int nsites)
{
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if (tid >= nx * ny * nz) return;

	int kc = tid % nz + 1;
	int jc = (tid / nz) % ny + 1;
	int ic = tid / (nz * ny) + 1;

	int str_z = 1;
	int str_y = nz + 2 * nhalo;
	int str_x = str_y * (ny + 2 * nhalo);

	int src_idx = str_x * (nhalo + ic - 1) + str_y * (nhalo + jc - 1) + (nhalo + kc - 1);
	double rho0 = rho_src[nsites * 0 + src_idx];
	double rho1 = rho_src[nsites * 1 + src_idx];

	int i_min = max(1, ic - 1);  int i_max = min(nx, ic + 1);
	int j_min = max(1, jc - 1);  int j_max = min(ny, jc + 1);
	int k_min = max(1, kc - 1);  int k_max = min(nz, kc + 1);

	for (int i = i_min; i <= i_max; i++) {
		double wx = d_peskin_gpu((double)(ic - i));
		for (int j = j_min; j <= j_max; j++) {
			double wy = d_peskin_gpu((double)(jc - j));
			for (int k = k_min; k <= k_max; k++) {
				double dr = wx * wy * d_peskin_gpu((double)(kc - k));
				int dst = str_x * (nhalo + i - 1) + str_y * (nhalo + j - 1) + (nhalo + k - 1);
				atomicAdd(&rho_dst[nsites * 0 + dst], rho0 * dr);
				atomicAdd(&rho_dst[nsites * 1 + dst], rho1 * dr);
			}
		}
	}
}

/*****************************************************************************
 *
 *  peskin_scatter_particles_kernel
 *
 *  Each thread = one subgrid particle.
 *  Scatters q0/q1 to the Peskin neighbourhood in rho_dst via atomicAdd.
 *
 *****************************************************************************/
__global__ void peskin_scatter_particles_kernel(
	const peskin_particle_t* __restrict__ particles,
	int                                    nparticles,
	double* rho_dst,
	int nx, int ny, int nz, int nhalo, int nsites)
{
	int pid = blockIdx.x * blockDim.x + threadIdx.x;
	if (pid >= nparticles) return;

	double rx = particles[pid].rx;
	double ry = particles[pid].ry;
	double rz = particles[pid].rz;
	double q0 = particles[pid].q0;
	double q1 = particles[pid].q1;

	int str_z = 1;
	int str_y = nz + 2 * nhalo;
	int str_x = str_y * (ny + 2 * nhalo);

	int i_min = max(0, (int)floor(rx - 1.0));
	int i_max = min(nx + 1, (int)ceil(rx + 1.0));
	int j_min = max(0, (int)floor(ry - 1.0));
	int j_max = min(ny + 1, (int)ceil(ry + 1.0));
	int k_min = max(0, (int)floor(rz - 1.0));
	int k_max = min(nz + 1, (int)ceil(rz + 1.0));

	for (int i = i_min; i <= i_max; i++) {
		double wx = d_peskin_gpu(rx - i);
		for (int j = j_min; j <= j_max; j++) {
			double wy = d_peskin_gpu(ry - j);
			for (int k = k_min; k <= k_max; k++) {
				double dr = wx * wy * d_peskin_gpu(rz - k);
				int dst = str_x * (nhalo + i - 1) + str_y * (nhalo + j - 1) + (nhalo + k - 1);
				atomicAdd(&rho_dst[nsites * 0 + dst], q0 * dr);
				atomicAdd(&rho_dst[nsites * 1 + dst], q1 * dr);
			}
		}
	}
}

#endif /* __NVCC__ */

/*****************************************************************************
 *
 *  subgrid_peskin_scatter_rho_buf
 *
 *  Peskin scatter of fluid nodes + subgrid particles into a plain double[]
 *  buffer (rho_smooth, already zeroed by caller, size = nsites * nk).
 *  No second psi_t needed — allocates a temporary GPU buffer internally.
 *
 *  On entry:  psi_src->rho host data is valid; rho_smooth is zeroed.
 *  On exit:   rho_smooth contains the Peskin-smoothed charge densities.
 *
 *****************************************************************************/
/*CHANGE INIT - subgrid_peskin_scatter_rho_buf */
int subgrid_peskin_scatter_rho_buf(colloids_info_t* cinfo,
                                    psi_t* psi_src,
                                    double* rho_smooth, int ndata)
{
	assert(cinfo);
	assert(psi_src);
	assert(rho_smooth);

#ifdef __NVCC__

	int nlocal[3], offset[3], nhalo;
	cs_nlocal(psi_src->cs, nlocal);
	cs_nlocal_offset(psi_src->cs, offset);
	cs_nhalo(psi_src->cs, &nhalo);

	int nx = nlocal[X], ny = nlocal[Y], nz = nlocal[Z];
	int nsites = psi_src->nsites;

	/* Upload psi_src->rho to device */
	field_memcpy(psi_src->rho, tdpMemcpyHostToDevice);

	/* Get device pointer for source rho */
	size_t data_offset = offsetof(field_t, data);
	double* rho_src_d = NULL;
	cudaMemcpy(&rho_src_d, (char*)(psi_src->rho->target) + data_offset,
	           sizeof(double*), cudaMemcpyDeviceToHost);

	/* Allocate zeroed destination buffer on GPU */
	double* rho_dst_d = NULL;
	cudaMalloc(&rho_dst_d, (size_t)ndata * sizeof(double));
	cudaMemset(rho_dst_d, 0, (size_t)ndata * sizeof(double));

	int threads = 256;

	/* Kernel 1: scatter LB-node charges */
	int n_interior = nx * ny * nz;
	peskin_scatter_nodes_kernel<<<(n_interior + threads - 1) / threads, threads>>>(
		rho_src_d, rho_dst_d, nx, ny, nz, nhalo, nsites);

	/* Build flat particle array and scatter on GPU */
	int ncell[3];
	colloids_info_ncell(cinfo, ncell);

	int npart = 0;
	for (int ic = 0; ic <= ncell[X] + 1; ic++)
		for (int jc = 0; jc <= ncell[Y] + 1; jc++)
			for (int kc = 0; kc <= ncell[Z] + 1; kc++) {
				colloid_t* p = NULL;
				colloids_info_cell_list_head(cinfo, ic, jc, kc, &p);
				for (; p; p = p->next)
					if (p->s.bc == COLLOID_BC_SUBGRID) npart++;
			}

	if (npart > 0) {
		peskin_particle_t* parts_h =
			(peskin_particle_t*)malloc(npart * sizeof(peskin_particle_t));
		int idx = 0;
		for (int ic = 0; ic <= ncell[X] + 1; ic++)
			for (int jc = 0; jc <= ncell[Y] + 1; jc++)
				for (int kc = 0; kc <= ncell[Z] + 1; kc++) {
					colloid_t* p = NULL;
					colloids_info_cell_list_head(cinfo, ic, jc, kc, &p);
					for (; p; p = p->next) {
						if (p->s.bc != COLLOID_BC_SUBGRID) continue;
						parts_h[idx].rx = p->s.r[X] - 1.0 * offset[X];
						parts_h[idx].ry = p->s.r[Y] - 1.0 * offset[Y];
						parts_h[idx].rz = p->s.r[Z] - 1.0 * offset[Z];
						parts_h[idx].q0 = p->s.q0;
						parts_h[idx].q1 = p->s.q1;
						idx++;
					}
				}

		peskin_particle_t* parts_d = NULL;
		cudaMalloc(&parts_d, npart * sizeof(peskin_particle_t));
		cudaMemcpy(parts_d, parts_h, npart * sizeof(peskin_particle_t),
		           cudaMemcpyHostToDevice);

		/* Kernel 2: scatter particle charges */
		peskin_scatter_particles_kernel<<<(npart + threads - 1) / threads, threads>>>(
			parts_d, npart, rho_dst_d, nx, ny, nz, nhalo, nsites);

		cudaFree(parts_d);
		free(parts_h);
	}

	cudaDeviceSynchronize();

	/* Download result into host buffer */
	cudaMemcpy(rho_smooth, rho_dst_d, (size_t)ndata * sizeof(double),
	           cudaMemcpyDeviceToHost);
	cudaFree(rho_dst_d);

#else
	/* CPU fallback */
	int nlocal[3], offset[3];
	cs_nlocal(psi_src->cs, nlocal);
	cs_nlocal_offset(psi_src->cs, offset);
	int nsites = psi_src->nsites;
	int nk     = psi_src->nk;

	/* Scatter LB-node charges */
	for (int i = 1; i <= nlocal[X]; i++) {
		for (int j = 1; j <= nlocal[Y]; j++) {
			for (int k = 1; k <= nlocal[Z]; k++) {
				int src = cs_index(psi_src->cs, i, j, k);
				double rho0, rho1;
				psi_rho(psi_src, src, 0, &rho0);
				psi_rho(psi_src, src, 1, &rho1);
				int i_min = imax(1, i - 1), i_max = imin(nlocal[X], i + 1);
				int j_min = imax(1, j - 1), j_max = imin(nlocal[Y], j + 1);
				int k_min = imax(1, k - 1), k_max = imin(nlocal[Z], k + 1);
				for (int ii = i_min; ii <= i_max; ii++)
					for (int jj = j_min; jj <= j_max; jj++)
						for (int kk = k_min; kk <= k_max; kk++) {
							double dr = d_peskin(i - ii) * d_peskin(j - jj) * d_peskin(k - kk);
							int dst = cs_index(psi_src->cs, ii, jj, kk);
							rho_smooth[addr_rank1(nsites, nk, dst, 0)] += rho0 * dr;
							rho_smooth[addr_rank1(nsites, nk, dst, 1)] += rho1 * dr;
						}
			}
		}
	}

	/* Scatter subgrid particle charges */
	if (cinfo->nsubgrid > 0) {
		int ncell[3];
		colloids_info_ncell(cinfo, ncell);
		for (int ic = 0; ic <= ncell[X] + 1; ic++)
			for (int jc = 0; jc <= ncell[Y] + 1; jc++)
				for (int kc = 0; kc <= ncell[Z] + 1; kc++) {
					colloid_t* p = NULL;
					colloids_info_cell_list_head(cinfo, ic, jc, kc, &p);
					for (; p; p = p->next) {
						if (p->s.bc != COLLOID_BC_SUBGRID) continue;
						double r0[3] = { p->s.r[X] - 1.0*offset[X],
						                 p->s.r[Y] - 1.0*offset[Y],
						                 p->s.r[Z] - 1.0*offset[Z] };
						int i_min = imax(1, (int)floor(r0[X]-1.0));
						int i_max = imin(nlocal[X], (int)ceil(r0[X]+1.0));
						int j_min = imax(1, (int)floor(r0[Y]-1.0));
						int j_max = imin(nlocal[Y], (int)ceil(r0[Y]+1.0));
						int k_min = imax(1, (int)floor(r0[Z]-1.0));
						int k_max = imin(nlocal[Z], (int)ceil(r0[Z]+1.0));
						for (int i = i_min; i <= i_max; i++)
							for (int j = j_min; j <= j_max; j++)
								for (int k = k_min; k <= k_max; k++) {
									double dr = d_peskin(r0[X]-i) * d_peskin(r0[Y]-j) * d_peskin(r0[Z]-k);
									int dst = cs_index(psi_src->cs, i, j, k);
									rho_smooth[addr_rank1(nsites, nk, dst, 0)] += p->s.q0 * dr;
									rho_smooth[addr_rank1(nsites, nk, dst, 1)] += p->s.q1 * dr;
								}
					}
				}
	}
#endif

	return 0;
}
/*CHANGE END - subgrid_peskin_scatter_rho_buf */

/*****************************************************************************
 *
 *  subgrid_peskin_scatter_rho_gpu
 *
 *  Distributes onto psi_dst->rho:
 *    1) Fluid charge from psi_src->rho (LB nodes) via Peskin scatter on GPU
 *    2) Subgrid particle charges via Peskin scatter on GPU
 *
 *  Replaces both the LB-node loop and subgrid_charge_from_particles.
 *  psi_dst->rho must be zeroed by the caller before this function.
 *
 *  On entry:  psi_src->rho is valid on the host; psi_dst->rho is zeroed.
 *  On exit:   psi_dst->rho is updated on the host.
 *
 *****************************************************************************/
int subgrid_peskin_scatter_rho_gpu(colloids_info_t* cinfo,
								   psi_t* psi_src, psi_t* psi_dst)
{
	assert(cinfo);
	assert(psi_src);
	assert(psi_dst);

#ifdef __NVCC__

	int nlocal[3], offset[3], nhalo;
	cs_nlocal(psi_src->cs, nlocal);
	cs_nlocal_offset(psi_src->cs, offset);
	cs_nhalo(psi_src->cs, &nhalo);

	int nx = nlocal[X], ny = nlocal[Y], nz = nlocal[Z];
	int nsites = psi_src->nsites;

	/* --- Upload rho fields to device --- */
	field_memcpy(psi_src->rho, tdpMemcpyHostToDevice);
	field_memcpy(psi_dst->rho, tdpMemcpyHostToDevice);

	size_t data_offset = offsetof(field_t, data);
	double* rho_src_d = NULL;
	double* rho_dst_d = NULL;
	cudaMemcpy(&rho_src_d, (char*)(psi_src->rho->target) + data_offset,
			   sizeof(double*), cudaMemcpyDeviceToHost);
	cudaMemcpy(&rho_dst_d, (char*)(psi_dst->rho->target) + data_offset,
			   sizeof(double*), cudaMemcpyDeviceToHost);

	/* --- Kernel 1: scatter LB-node charges --- */
	int n_interior = nx * ny * nz;
	int threads = 256;
	peskin_scatter_nodes_kernel << <(n_interior + threads - 1) / threads, threads >> > (
		rho_src_d, rho_dst_d, nx, ny, nz, nhalo, nsites);

	/* --- Build flat particle array on CPU, upload, scatter on GPU --- */
	int ncell[3];
	colloids_info_ncell(cinfo, ncell);

	/* Count subgrid particles first to alloc */
	int npart = 0;
	for (int ic = 0; ic <= ncell[X] + 1; ic++)
		for (int jc = 0; jc <= ncell[Y] + 1; jc++)
			for (int kc = 0; kc <= ncell[Z] + 1; kc++) {
				colloid_t* p = NULL;
				colloids_info_cell_list_head(cinfo, ic, jc, kc, &p);
				for (; p; p = p->next)
					if (p->s.bc == COLLOID_BC_SUBGRID) npart++;
			}

	if (npart > 0) {
		peskin_particle_t* parts_h =
			(peskin_particle_t*)malloc(npart * sizeof(peskin_particle_t));

		int idx = 0;
		for (int ic = 0; ic <= ncell[X] + 1; ic++)
			for (int jc = 0; jc <= ncell[Y] + 1; jc++)
				for (int kc = 0; kc <= ncell[Z] + 1; kc++) {
					colloid_t* p = NULL;
					colloids_info_cell_list_head(cinfo, ic, jc, kc, &p);
					for (; p; p = p->next) {
						if (p->s.bc != COLLOID_BC_SUBGRID) continue;
						parts_h[idx].rx = p->s.r[X] - 1.0 * offset[X];
						parts_h[idx].ry = p->s.r[Y] - 1.0 * offset[Y];
						parts_h[idx].rz = p->s.r[Z] - 1.0 * offset[Z];
						parts_h[idx].q0 = p->s.q0;
						parts_h[idx].q1 = p->s.q1;
						idx++;
					}
				}

		peskin_particle_t* parts_d = NULL;
		cudaMalloc(&parts_d, npart * sizeof(peskin_particle_t));
		cudaMemcpy(parts_d, parts_h, npart * sizeof(peskin_particle_t),
				   cudaMemcpyHostToDevice);

		/* --- Kernel 2: scatter particle charges --- */
		peskin_scatter_particles_kernel << <(npart + threads - 1) / threads, threads >> > (
			parts_d, npart, rho_dst_d, nx, ny, nz, nhalo, nsites);

		cudaFree(parts_d);
		free(parts_h);
	}

	cudaDeviceSynchronize();

	/* Download result back to host */
	field_memcpy(psi_dst->rho, tdpMemcpyDeviceToHost);

#else
	/* CPU fallback (compiled without CUDA) */
	int nlocal[3], offset[3];
	cs_nlocal(psi_src->cs, nlocal);
	cs_nlocal_offset(psi_src->cs, offset);
	const double range = 1.0;

	/* LB nodes */
	for (int i = 1; i <= nlocal[X]; i++) {
		for (int j = 1; j <= nlocal[Y]; j++) {
			for (int k = 1; k <= nlocal[Z]; k++) {
				int src_idx = cs_index(psi_src->cs, i, j, k);
				double rho0, rho1;
				psi_rho(psi_src, src_idx, 0, &rho0);
				psi_rho(psi_src, src_idx, 1, &rho1);
				int i_min = imax(1, (int)floor(i - range)), i_max = imin(nlocal[X], (int)ceil(i + range));
				int j_min = imax(1, (int)floor(j - range)), j_max = imin(nlocal[Y], (int)ceil(j + range));
				int k_min = imax(1, (int)floor(k - range)), k_max = imin(nlocal[Z], (int)ceil(k + range));
				for (int ii = i_min; ii <= i_max; ii++)
					for (int jj = j_min; jj <= j_max; jj++)
						for (int kk = k_min; kk <= k_max; kk++) {
							double dr = d_peskin(i - ii) * d_peskin(j - jj) * d_peskin(k - kk);
							int dst_idx = cs_index(psi_dst->cs, ii, jj, kk);
							double r0, r1;
							psi_rho(psi_dst, dst_idx, 0, &r0);
							psi_rho(psi_dst, dst_idx, 1, &r1);
							psi_rho_set(psi_dst, dst_idx, 0, r0 + rho0 * dr);
							psi_rho_set(psi_dst, dst_idx, 1, r1 + rho1 * dr);
						}
			}
		}
	}

	/* Subgrid particles */
	if (cinfo->nsubgrid > 0) {
		int ncell[3];
		colloids_info_ncell(cinfo, ncell);
		for (int ic = 0; ic <= ncell[X] + 1; ic++)
			for (int jc = 0; jc <= ncell[Y] + 1; jc++)
				for (int kc = 0; kc <= ncell[Z] + 1; kc++) {
					colloid_t* p = NULL;
					colloids_info_cell_list_head(cinfo, ic, jc, kc, &p);
					for (; p; p = p->next) {
						if (p->s.bc != COLLOID_BC_SUBGRID) continue;
						double r0[3] = { p->s.r[X] - 1.0 * offset[X],
										 p->s.r[Y] - 1.0 * offset[Y],
										 p->s.r[Z] - 1.0 * offset[Z] };
						int i_min = imax(1, (int)floor(r0[X] - range)), i_max = imin(nlocal[X], (int)ceil(r0[X] + range));
						int j_min = imax(1, (int)floor(r0[Y] - range)), j_max = imin(nlocal[Y], (int)ceil(r0[Y] + range));
						int k_min = imax(1, (int)floor(r0[Z] - range)), k_max = imin(nlocal[Z], (int)ceil(r0[Z] + range));
						for (int i = i_min; i <= i_max; i++)
							for (int j = j_min; j <= j_max; j++)
								for (int k = k_min; k <= k_max; k++) {
									double dr = d_peskin(r0[X] - i) * d_peskin(r0[Y] - j) * d_peskin(r0[Z] - k);
									int dst_idx = cs_index(psi_dst->cs, i, j, k);
									double rh0, rh1;
									psi_rho(psi_dst, dst_idx, 0, &rh0);
									psi_rho(psi_dst, dst_idx, 1, &rh1);
									psi_rho_set(psi_dst, dst_idx, 0, rh0 + p->s.q0 * dr);
									psi_rho_set(psi_dst, dst_idx, 1, rh1 + p->s.q1 * dr);
								}
					}
				}
	}
#endif /* __NVCC__ */

	return 0;
}
/*CHANGE END - Peskin scatter kernel on GPU */
/*CHANGE END - Poisson-Boltzmann weight */