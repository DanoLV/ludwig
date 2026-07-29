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

/*CHANGE INIT - 20260424 Kaiser-Bessel globals */
static double subgrid_kb4_beta_ = 6.0;
static double subgrid_kb4_norm_fluid_ = 1.0;
static double i0_series(double x);
/*CHANGE END - 20260424 Kaiser-Bessel globals */

/*CHANGE INIT - 20260630 Hann spread/gather kernel order (configurable) */
/* Order n = full support width (in lattice units). w_i = 0 for |r| >= n/2.
 * Default n = 4 (4-point support, |r| <= 2). Set from input via
 * subgrid_set_hann_order(). Must be a positive even integer for the support
 * to align with an integer number of nodes. */
static double subgrid_hann_order_ = 4.0;
/*CHANGE END - 20260630 */

/*CHANGE INIT - 20260425 forward declare for interlacing offset variant */
int subgrid_charge_from_grid_offset(colloids_info_t* cinfo, psi_t* obj,
									 distributed_charge_klein_t** charge,
									 subgrid_kernel_t kernel, double mesh_offset);
/*CHANGE END - 20260425 */

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
 /*CHANGE INIT - 20260422 kernel parameter for subgrid_charge_from_particles */
int subgrid_charge_from_particles(colloids_info_t* cinfo, psi_t* obj, distributed_charge_klein_t** charge,
								   subgrid_kernel_t kernel)
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

	// range for used kernel, e.g. 2 for B-spline4, 3 for B-spline6, 2 for Peskin4
	int krange = subgrid_get_range(kernel);

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

					subgrid_get_lattice_index_range(r0, krange, nlocal, &i_min, &i_max, &j_min, &j_max, &k_min, &k_max);

					double kb_w[4][4][4] = { {{0}} };
					double kb_w_sum = 1.0;
					if (kernel == SUBGRID_KERNEL_KB4) {
						kb_w_sum = 0.0;
						for (i = i_min; i <= i_max; i++)
							for (j = j_min; j <= j_max; j++)
								for (k = k_min; k <= k_max; k++) {
									double wx = d_kb4(r0[X] - i);
									double wy = d_kb4(r0[Y] - j);
									double wz = d_kb4(r0[Z] - k);
									kb_w[i - i_min][j - j_min][k - k_min] = wx * wy * wz;
									kb_w_sum += wx * wy * wz;
								}
					}

					for (i = i_min; i <= i_max; i++) {
						for (j = j_min; j <= j_max; j++) {
							for (k = k_min; k <= k_max; k++) {

								index = cs_index(cinfo->cs, i, j, k);

								r[X] = r0[X] - 1.0 * i;
								r[Y] = r0[Y] - 1.0 * j;
								r[Z] = r0[Z] - 1.0 * k;

								if (kernel == SUBGRID_KERNEL_BSPLINE6) {
									dr = d_bspline6(r[X]) * d_bspline6(r[Y]) * d_bspline6(r[Z]);
								}
								else if (kernel == SUBGRID_KERNEL_BSPLINE4) {
									dr = d_bspline4(r[X]) * d_bspline4(r[Y]) * d_bspline4(r[Z]);
								}
								else if (kernel == SUBGRID_KERNEL_PESKIN4) {
									dr = d_peskin(r[X]) * d_peskin(r[Y]) * d_peskin(r[Z]);
								}
								else if (kernel == SUBGRID_KERNEL_KB4) {
									dr = kb_w[i - i_min][j - j_min][k - k_min] / kb_w_sum;
								}
								else if (kernel == SUBGRID_KERNEL_PESKIN6) {
									dr = d_peskin6(r[X]) * d_peskin6(r[Y]) * d_peskin6(r[Z]);
								}
								/*CHANGE INIT - 20260630 Hann kernel */
								else if (kernel == SUBGRID_KERNEL_HANN) {
									dr = d_hann(r[X]) * d_hann(r[Y]) * d_hann(r[Z]);
								}
								/*CHANGE END - 20260630 */

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
/*CHANGE END - 20260422 kernel parameter for subgrid_charge_from_particles */

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
 /*CHANGE INIT - 20260422 kernel parameter for subgrid_charge_from_grid */
int subgrid_charge_from_grid(colloids_info_t* cinfo, psi_t* obj,
							  distributed_charge_klein_t** charge,
							  subgrid_kernel_t kernel)
{
	return subgrid_charge_from_grid_offset(cinfo, obj, charge, kernel, 0.0);
}

/*CHANGE INIT - 20260427 subgrid_scatter_fluid_offset: scatter fluid rho into auxiliary buffer */
/*****************************************************************************
 *
 *  subgrid_scatter_fluid_offset
 *
 *  Scatters fluid charge from psi->rho into rho_buf using the kernel with
 *  mesh_offset applied to the X axis. Writes directly to rho_buf (no
 *  add_charge_to_array), preserving charge exactly via periodic wrap.
 *
 *  On entry:  psi->rho contains the fluid-only charge (rho_saved state).
 *             rho_buf is zeroed by caller, size = nsites * nk.
 *  On exit:   rho_buf contains the offset-smoothed fluid charge.
 *             psi->rho is unchanged.
 *
 *****************************************************************************/
int subgrid_scatter_fluid_offset(colloids_info_t* cinfo, psi_t* obj,
								  subgrid_kernel_t kernel, double mesh_offset,
								  double* rho_buf, int ndata)
{
	int i, j, k, i2, j2, k2;
	int i_min, i_max, j_min, j_max, k_min, k_max;
	int nlocal[3];

	assert(cinfo); assert(obj); assert(rho_buf);
	cs_nlocal(cinfo->cs, nlocal);
	int krange = subgrid_get_range(kernel);

	/* Fractional part of offset determines if border treatment is needed */
	double frac_x = mesh_offset - floor(mesh_offset);

	/* Build X source list: interior only, unless frac_x != 0 */
	int ix_excl_lo = 0, ix_excl_hi = -1;
	int ix_halo_lo = 1, ix_halo_hi = 0;
	if (frac_x > 0.0) {
		ix_excl_lo = nlocal[X] - krange + 1; ix_excl_hi = nlocal[X];
		ix_halo_lo = 1 - krange;             ix_halo_hi = 0;
	}
	else if (frac_x < 0.0) {
		ix_excl_lo = 1;               ix_excl_hi = krange;
		ix_halo_lo = nlocal[X] + 1;  ix_halo_hi = nlocal[X] + krange;
	}

	int list_cap = nlocal[X] + 2 * krange + 4;
	int* i_list = (int*)malloc(list_cap * sizeof(int)); int i_list_n = 0;
	for (i = 1; i <= nlocal[X]; i++) { if (i < ix_excl_lo || i > ix_excl_hi) i_list[i_list_n++] = i; }
	for (i = ix_halo_lo; i <= ix_halo_hi; i++) i_list[i_list_n++] = i;

	for (int ii = 0; ii < i_list_n; ii++) {
		i = i_list[ii];
		int iw = i; if (iw < 1) iw += nlocal[X]; else if (iw > nlocal[X]) iw -= nlocal[X];
		for (j = 1; j <= nlocal[Y]; j++) {
			for (k = 1; k <= nlocal[Z]; k++) {

				int index_src = cs_index(cinfo->cs, iw, j, k);
				double rho0, rho1;
				psi_rho(obj, index_src, 0, &rho0);
				psi_rho(obj, index_src, 1, &rho1);
				if (rho0 == 0.0 && rho1 == 0.0) continue;

				double r0[3] = { (double)i + mesh_offset, (double)j, (double)k };
				subgrid_get_lattice_index_range_halo(r0, krange, nlocal,
													 &i_min, &i_max, &j_min, &j_max, &k_min, &k_max);

				for (i2 = i_min; i2 <= i_max; i2++) {
					for (j2 = j_min; j2 <= j_max; j2++) {
						for (k2 = k_min; k2 <= k_max; k2++) {

							double dx = r0[X] - (double)i2;
							double dy = r0[Y] - (double)j2;
							double dz = r0[Z] - (double)k2;
							double dr;
							if (kernel == SUBGRID_KERNEL_BSPLINE6) dr = d_bspline6(dx) * d_bspline6(dy) * d_bspline6(dz);
							else if (kernel == SUBGRID_KERNEL_BSPLINE4) dr = d_bspline4(dx) * d_bspline4(dy) * d_bspline4(dz);
							else if (kernel == SUBGRID_KERNEL_PESKIN4)  dr = d_peskin(dx) * d_peskin(dy) * d_peskin(dz);
							else if (kernel == SUBGRID_KERNEL_KB4)      dr = d_kb4(dx) * d_kb4(dy) * d_kb4(dz) / subgrid_kb4_norm_fluid();
							else if (kernel == SUBGRID_KERNEL_PESKIN6)  dr = d_peskin6(dx) * d_peskin6(dy) * d_peskin6(dz);
							/*CHANGE INIT - 20260710 trilinear: exact relabelling at integer offsets */
							else if (kernel == SUBGRID_KERNEL_TRILINEAR) dr = d_trilinear(dx) * d_trilinear(dy) * d_trilinear(dz);
							/*CHANGE END - 20260710 */
							/*CHANGE INIT - 20260630 Hann kernel */
							else if (kernel == SUBGRID_KERNEL_HANN)     dr = d_hann(dx) * d_hann(dy) * d_hann(dz);
							/*CHANGE END - 20260630 */
							else dr = 0.0;
							if (dr == 0.0) continue;

							/* Periodic wrap for destination */
							int idw = i2, jdw = j2, kdw = k2;
							if (idw < 1) idw += nlocal[X]; else if (idw > nlocal[X]) idw -= nlocal[X];
							if (jdw < 1) jdw += nlocal[Y]; else if (jdw > nlocal[Y]) jdw -= nlocal[Y];
							if (kdw < 1) kdw += nlocal[Z]; else if (kdw > nlocal[Z]) kdw -= nlocal[Z];
							int index_dst = cs_index(cinfo->cs, idw, jdw, kdw);

							rho_buf[addr_rank1(obj->nsites, obj->nk, index_dst, 0)] += rho0 * dr;
							rho_buf[addr_rank1(obj->nsites, obj->nk, index_dst, 1)] += rho1 * dr;
						}
					}
				}
			}
		}
	}
	free(i_list);

	/* DIAG: print net charge in rho_buf per X column (first call only) */
	{
		static int sfo_diag = 0; sfo_diag++;
		if (sfo_diag <= 1) {
			double q_tot = 0.0;
			for (int ix = 1; ix <= nlocal[X]; ix++) {
				double q_ix = 0.0;
				for (int jx = 1; jx <= nlocal[Y]; jx++)
					for (int kx = 1; kx <= nlocal[Z]; kx++) {
						int idx = cs_index(obj->cs, ix, jx, kx);
						double r0v = rho_buf[addr_rank1(obj->nsites, obj->nk, idx, 0)];
						double r1v = rho_buf[addr_rank1(obj->nsites, obj->nk, idx, 1)];
						q_ix += r0v - r1v;
					}
				printf("[sfo_col] call=%d offset=%.3f ix=%d Q_net=%.6e\n", sfo_diag, mesh_offset, ix, q_ix);
				q_tot += q_ix;
			}
			printf("[sfo_col] call=%d offset=%.3f total_Q_net=%.6e\n", sfo_diag, mesh_offset, q_tot);

		}
	}

	return 0;
}
/*CHANGE END - 20260427 subgrid_scatter_fluid_offset */

/*CHANGE INIT - 20260425 subgrid_charge_from_grid_offset for interlacing */
int subgrid_charge_from_grid_offset(colloids_info_t* cinfo, psi_t* obj,
							  distributed_charge_klein_t** charge,
							  subgrid_kernel_t kernel, double mesh_offset)
{
	int i, j, k;
	int i2, j2, k2;
	int i_min, i_max, j_min, j_max, k_min, k_max;
	int index, index2;
	int nlocal[3];
	double rho0, rho1, dr;

	/* Temporary list of charged source nodes with their coordinates */
	typedef struct { int idx, si, sj, sk; double rho0, rho1; } src_t;
	int src_cap = 64, src_n = 0;
	src_t* srcs = (src_t*)malloc(src_cap * sizeof(src_t));

	assert(cinfo);
	assert(obj);
	assert(charge);
	assert(srcs);

	cs_nlocal(cinfo->cs, nlocal);

	/* Integer half-support radius matching the chosen kernel */
	// range for used kernel, e.g. 2 for B-spline4, 3 for B-spline6, 2 for Peskin4
	int krange = subgrid_get_range(kernel);

	/* Pass 1: collect source nodes. With mesh_offset != 0, the kernel support
	 * shifts so that border interior nodes scatter charge outside the domain.
	 * To handle periodic BC correctly, we replace those border interior nodes
	 * with the corresponding halo nodes (which hold the same charge values
	 * after psi_halo_rho). This ensures charge is redistributed correctly
	 * across the periodic boundary without double-counting.
	 *
	 * Convention: offset > 0 shifts the source to the right (+X), so:
	 *   - exclude interior nodes i = nlocal-krange+1..nlocal (right border)
	 *   - include halo nodes     i = 1-krange..0             (left halo)
	 *   (halo nodes hold the same charge as the right-border interior nodes
	 *    but scattering from their actual halo position gives correct wrapping)
	 * For offset < 0, the exclusion/inclusion is mirrored.
	 * For offset == 0, use standard interior range (no change). */

	 /* Build source index lists for each axis. With a non-integer offset the
	  * kernel support from border nodes crosses the periodic boundary
	  * asymmetrically. Replace those border nodes with the corresponding halo
	  * nodes so the scatter lands correctly without double-counting.
	  * With an integer offset (frac == 0) no replacement is needed. */
	  /* Exclusion/inclusion only on X axis (offset applied only in X).
	   * With a non-integer offset the kernel support from X-border nodes
	   * crosses the periodic boundary: replace those nodes with halo nodes. */
	double frac_x = mesh_offset - floor(mesh_offset);
	int ix_excl_lo = 0, ix_excl_hi = -1; /* empty by default */
	int ix_halo_lo = 1, ix_halo_hi = 0;
	if (frac_x > 0.0) {
		ix_excl_lo = nlocal[X] - krange + 1; ix_excl_hi = nlocal[X];
		ix_halo_lo = 1 - krange;             ix_halo_hi = 0;
	}
	else if (frac_x < 0.0) {
		ix_excl_lo = 1;               ix_excl_hi = krange;
		ix_halo_lo = nlocal[X] + 1;  ix_halo_hi = nlocal[X] + krange;
	}

	int list_cap = nlocal[X] + 2 * krange + 4;
	int* i_list = (int*)malloc(list_cap * sizeof(int)); int i_list_n = 0;

	for (i = 1; i <= nlocal[X]; i++) { if (i < ix_excl_lo || i > ix_excl_hi) i_list[i_list_n++] = i; }
	for (i = ix_halo_lo; i <= ix_halo_hi; i++) i_list[i_list_n++] = i;

	for (int ii = 0; ii < i_list_n; ii++) {
		i = i_list[ii];
		int iw = i; if (iw < 1) iw += nlocal[X]; else if (iw > nlocal[X]) iw -= nlocal[X];
		for (j = 1; j <= nlocal[Y]; j++) {
			for (k = 1; k <= nlocal[Z]; k++) {

				index = cs_index(cinfo->cs, iw, j, k);

				psi_rho(obj, index, 0, &rho0);
				psi_rho(obj, index, 1, &rho1);

				if (rho0 == 0.0 && rho1 == 0.0) continue;

				/* Zero the interior node (halo i nodes are read-only) */
				if (iw >= 1 && iw <= nlocal[X]) {
					add_charge_to_array(charge, index, 0.0, 0.0, obj);
					psi_rho_set(obj, index, 0, 0.0);
					psi_rho_set(obj, index, 1, 0.0);
				}

				if (src_n >= src_cap) {
					src_cap *= 2;
					srcs = (src_t*)realloc(srcs, src_cap * sizeof(src_t));
				}
				srcs[src_n].idx = index;
				srcs[src_n].si = i;  /* unwarped: carries offset info */
				srcs[src_n].sj = j;
				srcs[src_n].sk = k;
				srcs[src_n].rho0 = rho0;
				srcs[src_n].rho1 = rho1;
				src_n++;
			}  /* k */
		}  /* j */
	}  /* ii */
	free(i_list);

	/* Diagnostic: total charge collected in Pass 1 */
	{
		static int p1_call = 0; p1_call++;
		if (p1_call <= 2) {
			double q0_p1 = 0.0, q1_p1 = 0.0;
			for (int e = 0; e < src_n; e++) { q0_p1 += srcs[e].rho0; q1_p1 += srcs[e].rho1; }
			printf("[scatter_p1] call=%d offset=%.3f src_n=%d Q0=%.6e Q1=%.6e\n",
				   p1_call, mesh_offset, src_n, q0_p1, q1_p1);
		}
	}

	/* Pass 2: spread each source node's charge to neighbours */

	static int p2_call = 0; p2_call++;
	double p2_q0_deposited = 0.0, p2_q1_deposited = 0.0;
	double p2_dr_min = 1e10;
	int p2_bad_si = -1, p2_bad_sj = -1, p2_bad_sk = -1;
	double p2_bad_sum = 0.0;

	for (int e = 0; e < src_n; e++) {

		int si = srcs[e].si, sj = srcs[e].sj, sk = srcs[e].sk;
		/* Apply mesh_offset: source node appears shifted by -mesh_offset on Grid B */
		double r0[3] = { (double)si + mesh_offset, (double)sj, (double)sk };
		rho0 = srcs[e].rho0;
		rho1 = srcs[e].rho1;

		subgrid_get_lattice_index_range_halo(r0, krange, nlocal, &i_min, &i_max, &j_min, &j_max, &k_min, &k_max);

		/* Track kernel weight sum for this source */
		double dr_sum = 0.0;

		for (i2 = i_min; i2 <= i_max; i2++) {
			for (j2 = j_min; j2 <= j_max; j2++) {
				for (k2 = k_min; k2 <= k_max; k2++) {

					double dx = r0[X] - (double)i2;
					double dy = r0[Y] - (double)j2;
					double dz = r0[Z] - (double)k2;

					if (kernel == SUBGRID_KERNEL_BSPLINE6) {
						dr = d_bspline6(dx) * d_bspline6(dy) * d_bspline6(dz);
					}
					else if (kernel == SUBGRID_KERNEL_BSPLINE4) {
						dr = d_bspline4(dx) * d_bspline4(dy) * d_bspline4(dz);
					}
					else if (kernel == SUBGRID_KERNEL_PESKIN4) {
						dr = d_peskin(dx) * d_peskin(dy) * d_peskin(dz);
					}
					else if (kernel == SUBGRID_KERNEL_KB4) {
						dr = d_kb4(dx) * d_kb4(dy) * d_kb4(dz)
							/ subgrid_kb4_norm_fluid_;
					}
					else if (kernel == SUBGRID_KERNEL_PESKIN6) {
						dr = d_peskin6(dx) * d_peskin6(dy) * d_peskin6(dz);
					}
					/*CHANGE INIT - 20260630 Hann kernel */
					else if (kernel == SUBGRID_KERNEL_HANN) {
						dr = d_hann(dx) * d_hann(dy) * d_hann(dz);
					}
					/*CHANGE END - 20260630 */
					// if (kernel == SUBGRID_KERNEL_BSPLINE6) {
					// 	dr = d_bspline6(r0[X] - (double)i2) 
					// 		* d_bspline6(r0[Y] - (double)j2)
					// 		* d_bspline6(r0[Z] - (double)k2);
					// }
					// else if (kernel == SUBGRID_KERNEL_BSPLINE4) {
					// 	dr = d_bspline4(r0[X] - (double)i2)
					// 		* d_bspline4(r0[Y] - (double)j2)
					// 		* d_bspline4(r0[Z] - (double)k2);
					// }
					// else if (kernel == SUBGRID_KERNEL_PESKIN4) {
					// 	dr = d_peskin(r0[X] - (double)i2)
					// 		* d_peskin(r0[Y] - (double)j2)
					// 		* d_peskin(r0[Z] - (double)k2);
					// }
					if (dr == 0.0) continue;

					/* Wrap halo indices to interior for periodic BC:
					 * with mesh_offset != 0, kernel support can reach halo nodes
					 * which get overwritten by psi_halo_rho — map them to interior. */
					int iw = i2, jw = j2, kw = k2;
					if (iw < 1) iw += nlocal[X]; else if (iw > nlocal[X]) iw -= nlocal[X];
					if (jw < 1) jw += nlocal[Y]; else if (jw > nlocal[Y]) jw -= nlocal[Y];
					if (kw < 1) kw += nlocal[Z]; else if (kw > nlocal[Z]) kw -= nlocal[Z];
					index2 = cs_index(cinfo->cs, iw, jw, kw);
					add_charge_to_array(charge, index2, rho0 * dr, rho1 * dr, obj);
					dr_sum += dr;
				}
			}
		}
		/* Track worst (furthest from 1) kernel sum */
		double diff = fabs(dr_sum - 1.0);
		if (diff > fabs(p2_dr_min - 1.0)) {
			p2_dr_min = dr_sum;
			p2_bad_si = si; p2_bad_sj = sj; p2_bad_sk = sk;
			p2_bad_sum = dr_sum;
		}
	}

	if (p2_call <= 2)
		printf("[scatter_p2] call=%d offset=%.3f worst_sum=%.8f at si=(%d,%d,%d)\n",
			   p2_call, mesh_offset, p2_bad_sum, p2_bad_si, p2_bad_sj, p2_bad_sk);

	free(srcs);

	if (*charge == NULL) return 0;

	/* Pass 3: write smoothed charges into psi */
	double q0_total = 0.0, q1_total = 0.0;
	int halo_count = 0;
	/* Accumulate per-ix charge for diagnostics */
	double* q0_per_ix = (double*)calloc(nlocal[X] + 2, sizeof(double));
	double* q1_per_ix = (double*)calloc(nlocal[X] + 2, sizeof(double));
	for (int e = 0; e < (*charge)->count; e++) {
		distributed_charge_klein_entry_t* entry = (*charge)->entries[e];
		double v0 = klein_sum(entry->rho0_sum);
		double v1 = klein_sum(entry->rho1_sum);
		q0_total += v0;
		q1_total += v1;
		int ijk[3];
		cs_index_to_ijk(obj->cs, entry->cs_index, ijk);
		if (ijk[X] < 1 || ijk[X] > nlocal[X] ||
			ijk[Y] < 1 || ijk[Y] > nlocal[Y] ||
			ijk[Z] < 1 || ijk[Z] > nlocal[Z]) halo_count++;
		else {
			q0_per_ix[ijk[X]] += v0;
			q1_per_ix[ijk[X]] += v1;
		}
		psi_rho_set(obj, entry->cs_index, 0, v0);
		psi_rho_set(obj, entry->cs_index, 1, v1);
	}
	/* Always print first call, then only if anomalies */
	static int scatter_call_count = 0;
	scatter_call_count++;
	if (scatter_call_count <= 2) {
		printf("[scatter_fluid] call=%d offset=%.3f src=%d dst=%d halo_dst=%d Q0=%.6e Q1=%.6e\n",
			   scatter_call_count, mesh_offset, src_n, (*charge)->count, halo_count, q0_total, q1_total);
		for (int ix = 1; ix <= nlocal[X]; ix++) {
			if (q0_per_ix[ix] != 0.0 || q1_per_ix[ix] != 0.0)
				printf("[scatter_fluid]   ix=%d Q0=%.6e Q1=%.6e\n", ix, q0_per_ix[ix], q1_per_ix[ix]);
		}
	}
	free(q0_per_ix);
	free(q1_per_ix);

	return 0;
}
/*CHANGE END - 20260422 kernel parameter for subgrid_charge_from_grid */
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
 /*CHANGE INIT - 20260422 kernel parameter for subgrid_update_forces_electrokinetics */
int subgrid_update_forces_electrokinetics(colloids_info_t* cinfo,
										  map_t* map,
										  physics_t* phys,
										  psi_t* psi,
										  hydro_t* hydro,
										  subgrid_kernel_t kernel) {

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

	/* Half-support radius matching the chosen kernel */
	// double krange;
	// if (kernel == SUBGRID_KERNEL_BSPLINE6) 	krange = 2.0;
	// else if (kernel == SUBGRID_KERNEL_BSPLINE4) 	krange = 1.0;
	// else if (kernel == SUBGRID_KERNEL_PESKIN4)  	krange = 1.0;
	// else                                        	krange = drange_;
	int krange = subgrid_get_range(kernel);

	/* Scatter electric force from particles onto fluid nodes */
	for (ic = 0; ic <= ncell[X] + 1; ic++) {
		for (jc = 0; jc <= ncell[Y] + 1; jc++) {
			for (kc = 0; kc <= ncell[Z] + 1; kc++) {

				colloids_info_cell_list_head(cinfo, ic, jc, kc, &pc);

				for (; pc; pc = pc->next) {

					if (pc->s.bc != COLLOID_BC_SUBGRID) continue;

					force[X] = kt * reunit * pc->Esub[X] * (pc->s.q0 - pc->s.q1);
					force[Y] = kt * reunit * pc->Esub[Y] * (pc->s.q0 - pc->s.q1);
					force[Z] = kt * reunit * pc->Esub[Z] * (pc->s.q0 - pc->s.q1);

					r0[X] = pc->s.r[X] - 1.0 * offset[X];
					r0[Y] = pc->s.r[Y] - 1.0 * offset[Y];
					r0[Z] = pc->s.r[Z] - 1.0 * offset[Z];

					subgrid_get_lattice_index_range(r0, krange, nlocal, &i_min, &i_max, &j_min, &j_max, &k_min, &k_max);

					double kb_w[4][4][4] = { {{0}} };
					double kb_w_sum = 1.0;
					if (kernel == SUBGRID_KERNEL_KB4) {
						kb_w_sum = 0.0;
						for (i = i_min; i <= i_max; i++)
							for (j = j_min; j <= j_max; j++)
								for (k = k_min; k <= k_max; k++) {
									double wx = d_kb4(r0[X] - i);
									double wy = d_kb4(r0[Y] - j);
									double wz = d_kb4(r0[Z] - k);
									kb_w[i - i_min][j - j_min][k - k_min] = wx * wy * wz;
									kb_w_sum += wx * wy * wz;
								}
					}

					for (i = i_min; i <= i_max; i++) {
						for (j = j_min; j <= j_max; j++) {
							for (k = k_min; k <= k_max; k++) {

								double force_aux[3] = { 0.0, 0.0, 0.0 };

								index = cs_index(cinfo->cs, i, j, k);

								r[X] = r0[X] - 1.0 * i;
								r[Y] = r0[Y] - 1.0 * j;
								r[Z] = r0[Z] - 1.0 * k;

								if (kernel == SUBGRID_KERNEL_BSPLINE6) {
									dr = d_bspline6(r[X]) * d_bspline6(r[Y]) * d_bspline6(r[Z]);
								}
								else if (kernel == SUBGRID_KERNEL_BSPLINE4) {
									dr = d_bspline4(r[X]) * d_bspline4(r[Y]) * d_bspline4(r[Z]);
								}
								else if (kernel == SUBGRID_KERNEL_PESKIN4) {
									dr = d_peskin(r[X]) * d_peskin(r[Y]) * d_peskin(r[Z]);
								}
								else if (kernel == SUBGRID_KERNEL_KB4) {
									dr = kb_w[i - i_min][j - j_min][k - k_min] / kb_w_sum;
								}
								else if (kernel == SUBGRID_KERNEL_PESKIN6) {
									dr = d_peskin6(r[X]) * d_peskin6(r[Y]) * d_peskin6(r[Z]);
								}
								/*CHANGE INIT - 20260630 Hann kernel */
								else if (kernel == SUBGRID_KERNEL_HANN) {
									dr = d_hann(r[X]) * d_hann(r[Y]) * d_hann(r[Z]);
								}
								/*CHANGE END - 20260630 */

								force_aux[X] = force[X] * dr;
								force_aux[Y] = force[Y] * dr;
								force_aux[Z] = force[Z] * dr;

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
	/*CHANGE END - 20260422 kernel parameter for subgrid_update_forces_electrokinetics */

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
 *  subgrid_update_forces_electrokinetics_ewald
 *
 *  Accumulate single particle force contributions from electric fields.
 *  Sum goes to fex as this force is calculated ouside of ludwig_colloids_update
 *  in order to apply the force in the same step as the applied field
 *****************************************************************************/
int subgrid_update_forces_electrokinetics_ewald(colloids_info_t* cinfo,
												map_t* map,
												physics_t* phys,
												psi_t* psi,
										  		hydro_t* hydro,
										  		subgrid_kernel_t kernel) {

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

	// INIT VERSION - Ewald											
	/* Calculate electric forces on particles and accumulate total force */
	for (ic = 0; ic <= ncell[X] + 1; ic++) {
		for (jc = 0; jc <= ncell[Y] + 1; jc++) {
			for (kc = 0; kc <= ncell[Z] + 1; kc++) {

				colloids_info_cell_list_head(cinfo, ic, jc, kc, &pc);

				for (; pc; pc = pc->next) {


					if (pc->s.bc != COLLOID_BC_SUBGRID) continue;

					klein_t force_k[3];
					force_k[X] = klein_zero();
					force_k[Y] = klein_zero();
					force_k[Z] = klein_zero();

					// Si uso Ewald	se asigna directamente la fuerza calculada en Ewald, que ya incluye el término de carga
					force[X] = pc->fex[X];
					force[Y] = pc->fex[Y];
					force[Z] = pc->fex[Z];

					/* Translate colloid position to local coordinates */
					r0[X] = pc->s.r[X] - 1.0 * offset[X];
					r0[Y] = pc->s.r[Y] - 1.0 * offset[Y];
					r0[Z] = pc->s.r[Z] - 1.0 * offset[Z];

					/* Work out which local lattice sites are involved */
					subgrid_get_lattice_index(r0, nlocal, &i_min, &i_max, &j_min, &j_max, &k_min, &k_max);

					double kb_w[4][4][4] = { {{0}} };
					double kb_w_sum = 1.0;
					if (kernel == SUBGRID_KERNEL_KB4) {
						kb_w_sum = 0.0;
						for (i = i_min; i <= i_max; i++)
							for (j = j_min; j <= j_max; j++)
								for (k = k_min; k <= k_max; k++) {
									double wx = d_kb4(r0[X] - i);
									double wy = d_kb4(r0[Y] - j);
									double wz = d_kb4(r0[Z] - k);
									kb_w[i - i_min][j - j_min][k - k_min] = wx * wy * wz;
									kb_w_sum += wx * wy * wz;
								}
					}

					for (i = i_min; i <= i_max; i++) {
						for (j = j_min; j <= j_max; j++) {
							for (k = k_min; k <= k_max; k++) {

								double force_aux[3] = { 0.0, 0.0, 0.0 };      /* force on particle from this lattice site */

								index = cs_index(cinfo->cs, i, j, k);

								/* Separation between r0 and the lattice site */
								r[X] = r0[X] - 1.0 * i;
								r[Y] = r0[Y] - 1.0 * j;
								r[Z] = r0[Z] - 1.0 * k;

								// dr = d_peskin(r[X]) * d_peskin(r[Y]) * d_peskin(r[Z]);

								if (kernel == SUBGRID_KERNEL_BSPLINE6) {
									dr = d_bspline6(r[X]) * d_bspline6(r[Y]) * d_bspline6(r[Z]);
								}
								else if (kernel == SUBGRID_KERNEL_BSPLINE4) {
									dr = d_bspline4(r[X]) * d_bspline4(r[Y]) * d_bspline4(r[Z]);
								}
								else if (kernel == SUBGRID_KERNEL_PESKIN4) {
									dr = d_peskin(r[X]) * d_peskin(r[Y]) * d_peskin(r[Z]);
								}
								else if (kernel == SUBGRID_KERNEL_KB4) {
									dr = kb_w[i - i_min][j - j_min][k - k_min] / kb_w_sum;
								}
								else if (kernel == SUBGRID_KERNEL_PESKIN6) {
									dr = d_peskin6(r[X]) * d_peskin6(r[Y]) * d_peskin6(r[Z]);
								}
								/*CHANGE INIT - 20260630 Hann kernel */
								else if (kernel == SUBGRID_KERNEL_HANN) {
									dr = d_hann(r[X]) * d_hann(r[Y]) * d_hann(r[Z]);
								}
								/*CHANGE END - 20260630 */

								/* Force on particle from electric field at this site index*/
								force_aux[X] = force[X] * dr;
								force_aux[Y] = force[Y] * dr;
								force_aux[Z] = force[Z] * dr;

								/* Add to Klein sum array with binary search */
								add_force_to_array(&force_k_indexed, index, force_aux);

							}
						}
					}
				}
			}
		}
	}
	// // END VERSION - Ewald

	colloid_sums_halo(cinfo, COLLOID_SUM_FORCE_EXT_ONLY);

	/* Apply accumulated forces to hydro */
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
 *  subgrid_print_Esub
 *
 *  Print electric field on subgrid particles
 *****************************************************************************/
int subgrid_print_Esub(colloids_info_t* cinfo,
						int step,
						FILE* fp) {

	int i, j, k, ic, jc, kc;
	int ncell[3];

	colloid_t* pc;
	// MPI_Comm comm;


	assert(cinfo);

	if (cinfo->nsubgrid == 0) return 0;

	// cs_cart_comm(cinfo->cs, &comm);
	colloids_info_ncell(cinfo, ncell);

	// colloid_sums_halo(cinfo, COLLOID_SUM_ELECTRIC_FIELD);

	/* Loop only over non-halo cells to write data */
	for (ic = 1; ic <= ncell[X]; ic++) {
		for (jc = 1; jc <= ncell[Y]; jc++) {
			for (kc = 1; kc <= ncell[Z]; kc++) {

				colloids_info_cell_list_head(cinfo, ic, jc, kc, &pc);

				for (; pc; pc = pc->next) {

					if (pc->s.bc != COLLOID_BC_SUBGRID) continue;

					char string[512];
					double Emod = sqrt(pc->Esub[X] * pc->Esub[X] + pc->Esub[Y] * pc->Esub[Y] + pc->Esub[Z] * pc->Esub[Z]);

					sprintf(string, "%d;%d;%.15e;%.15e;%.15e;%.15e\n",
							step, pc->s.index,
							Emod, pc->Esub[X], pc->Esub[Y], pc->Esub[Z]);
					for (i = 0; i < (int)strlen(string); i++) if (string[i] == '.') string[i] = ',';
					fprintf(fp, "%s", string);
				}
			}
		}
	}
	fflush(fp);

	return 0;

}

/*****************************************************************************
 *
 *  subgrid_update_Esub
 *
 *  Calculate electric field on subgrid particles
 *****************************************************************************/
 /*CHANGE INIT - 20260422 kernel parameter for subgrid_update_Esub */
int subgrid_update_Esub(colloids_info_t* cinfo,
						psi_t* psi,
						int step,
						FILE* fp,
						pe_t* pe,
						subgrid_kernel_t kernel) {

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

	/* Half-support radius for gather loop, matching the kernel used for scatter */
	// double krange;
	// if (kernel == SUBGRID_KERNEL_BSPLINE6) krange = 2.0;
	// else if (kernel == SUBGRID_KERNEL_BSPLINE4) krange = 1.0;
	// else if (kernel == SUBGRID_KERNEL_PESKIN4)  krange = 1.0;
	// else                                         krange = drange_;
	int krange = subgrid_get_range(kernel);

	/* Gather electric field onto particles using the chosen interpolation kernel */
	for (ic = 0; ic <= ncell[X] + 1; ic++) {
		for (jc = 0; jc <= ncell[Y] + 1; jc++) {
			for (kc = 0; kc <= ncell[Z] + 1; kc++) {

				colloids_info_cell_list_head(cinfo, ic, jc, kc, &pc);

				for (; pc; pc = pc->next) {

					if (pc->s.bc != COLLOID_BC_SUBGRID) continue;

					klein_t E_field_k[3];
					E_field_k[X] = klein_zero();
					E_field_k[Y] = klein_zero();
					E_field_k[Z] = klein_zero();

					r0[X] = pc->s.r[X] - 1.0 * offset[X];
					r0[Y] = pc->s.r[Y] - 1.0 * offset[Y];
					r0[Z] = pc->s.r[Z] - 1.0 * offset[Z];

					subgrid_get_lattice_index_range(r0, krange, nlocal, &i_min, &i_max, &j_min, &j_max, &k_min, &k_max);

					double kb_w[4][4][4] = { {{0}} };
					double kb_w_sum = 1.0;
					if (kernel == SUBGRID_KERNEL_KB4) {
						kb_w_sum = 0.0;
						for (i = i_min; i <= i_max; i++)
							for (j = j_min; j <= j_max; j++)
								for (k = k_min; k <= k_max; k++) {
									double wx = d_kb4(r0[X] - i);
									double wy = d_kb4(r0[Y] - j);
									double wz = d_kb4(r0[Z] - k);
									kb_w[i - i_min][j - j_min][k - k_min] = wx * wy * wz;
									kb_w_sum += wx * wy * wz;
								}
					}

					for (i = i_min; i <= i_max; i++) {
						for (j = j_min; j <= j_max; j++) {
							for (k = k_min; k <= k_max; k++) {

								index = cs_index(cinfo->cs, i, j, k);

								r[X] = r0[X] - 1.0 * i;
								r[Y] = r0[Y] - 1.0 * j;
								r[Z] = r0[Z] - 1.0 * k;

								if (kernel == SUBGRID_KERNEL_BSPLINE6) {
									dr = d_bspline6(r[X]) * d_bspline6(r[Y]) * d_bspline6(r[Z]);
								}
								else if (kernel == SUBGRID_KERNEL_BSPLINE4) {
									dr = d_bspline4(r[X]) * d_bspline4(r[Y]) * d_bspline4(r[Z]);
								}
								else if (kernel == SUBGRID_KERNEL_PESKIN4) {
									dr = d_peskin(r[X]) * d_peskin(r[Y]) * d_peskin(r[Z]);
								}
								else if (kernel == SUBGRID_KERNEL_KB4) {
									dr = kb_w[i - i_min][j - j_min][k - k_min] / kb_w_sum;
								}
								else if (kernel == SUBGRID_KERNEL_PESKIN6) {
									dr = d_peskin6(r[X]) * d_peskin6(r[Y]) * d_peskin6(r[Z]);
								}
								/*CHANGE INIT - 20260630 Hann kernel */
								else if (kernel == SUBGRID_KERNEL_HANN) {
									dr = d_hann(r[X]) * d_hann(r[Y]) * d_hann(r[Z]);
								}
								/*CHANGE END - 20260630 */


								psi_electric_field(psi, index, e);

								E_field[X] = e[X] * dr;
								E_field[Y] = e[Y] * dr;
								E_field[Z] = e[Z] * dr;

								// double psi_val = 0.0;
								// psi_psi(psi, index, &psi_val);
								// E_field[X] = - psi_val *d_peskin_derivative(r[X])* d_peskin(r[Y])* d_peskin(r[Z]);
								// E_field[Y] = - psi_val *d_peskin_derivative(r[Y])* d_peskin(r[X])* d_peskin(r[Z]);
								// E_field[Z] = - psi_val *d_peskin_derivative(r[Z])* d_peskin(r[X])* d_peskin(r[Y]);


								klein_add_double(&E_field_k[X], E_field[X]);
								klein_add_double(&E_field_k[Y], E_field[Y]);
								klein_add_double(&E_field_k[Z], E_field[Z]);
							}
						}
					}

					pc->Esub[X] = klein_sum(&E_field_k[X]);
					pc->Esub[Y] = klein_sum(&E_field_k[Y]);
					pc->Esub[Z] = klein_sum(&E_field_k[Z]);

					// /* DIAG: print particle position, gathered E-field, and psi at nodes */
					// {
					// 	static int esub_diag = 0; esub_diag++;
					// 	if (esub_diag <= 4) {
					// 		printf("[Esub_DIAG] call=%d r=(%.4f,%.4f,%.4f) r0=(%.4f,%.4f,%.4f) "
					// 			   "Esub=(%+.6e,%+.6e,%+.6e) range=[%d..%d]\n",
					// 			   esub_diag, pc->s.r[X], pc->s.r[Y], pc->s.r[Z],
					// 			   r0[X], r0[Y], r0[Z],
					// 			   pc->Esub[X], pc->Esub[Y], pc->Esub[Z],
					// 			   i_min, i_max);
					// 		int jm = (int)round(r0[Y]), km = (int)round(r0[Z]);
					// 		if (jm < 1) jm = 1; if (jm > nlocal[Y]) jm = nlocal[Y];
					// 		if (km < 1) km = 1; if (km > nlocal[Z]) km = nlocal[Z];
					// 		for (int ii_d = i_min; ii_d <= i_max; ii_d++) {
					// 			int id = cs_index(cinfo->cs, ii_d, jm, km);
					// 			double e_d[3]; psi_electric_field(psi, id, e_d);
					// 			double dr_d = d_peskin(r0[X] - ii_d);
					// 			printf("[Esub_DIAG]   i=%d psi_psi=%.6e Ex=%.6e dr=%.6e\n",
					// 				   ii_d,
					// 				   psi->psi->data[addr_rank0(psi->nsites, id)],
					// 				   e_d[X], dr_d);
					// 		}
					// 	}
					// }
				}
			}
		}
	}
	/*CHANGE END - 20260422 kernel parameter for subgrid_update_Esub */

	colloid_sums_halo(cinfo, COLLOID_SUM_ELECTRIC_FIELD);

	/* CHANGE INIT - Write particle data after colloid_sums_halo, only for local (non-halo) particles */
	/* Loop only over non-halo cells to write data */
	for (ic = 1; ic <= ncell[X]; ic++) {
		for (jc = 1; jc <= ncell[Y]; jc++) {
			for (kc = 1; kc <= ncell[Z]; kc++) {

				colloids_info_cell_list_head(cinfo, ic, jc, kc, &pc);

				for (; pc; pc = pc->next) {

					if (pc->s.bc != COLLOID_BC_SUBGRID) continue;

					/*CHANGE INIT - 20260518 fix sprintf: remove dangling PB format specifiers */
					char string[512];
					double Emod = sqrt(pc->Esub[X] * pc->Esub[X] + pc->Esub[Y] * pc->Esub[Y] + pc->Esub[Z] * pc->Esub[Z]);
					/* Format: step;index;|E|;Ex;Ey;Ez */
					sprintf(string, "%d;%d;%.15e;%.15e;%.15e;%.15e\n",
							step, pc->s.index,
							Emod, pc->Esub[X], pc->Esub[Y], pc->Esub[Z]);
					for (i = 0; i < (int)strlen(string); i++) if (string[i] == '.') string[i] = ',';
					fprintf(fp, "%s", string);
					/*CHANGE END - 20260518 */
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
int subgrid_get_range(subgrid_kernel_t kernel)
{
	double krange;
	if (kernel == SUBGRID_KERNEL_BSPLINE6) 		return 2;
	else if (kernel == SUBGRID_KERNEL_BSPLINE4) return 1;
	else if (kernel == SUBGRID_KERNEL_PESKIN4)  return 1;
	else if (kernel == SUBGRID_KERNEL_KB4)      return 1;
	else if (kernel == SUBGRID_KERNEL_PESKIN6)  return 2;
	/*CHANGE INIT - 20260630 Hann kernel range = ceil(n/2) */
	else if (kernel == SUBGRID_KERNEL_HANN)     return (int) ceil(0.5 * subgrid_hann_order_);
	/*CHANGE END - 20260630 */
	/*CHANGE INIT - 20260710 trilinear (interpolating) kernel */
	else if (kernel == SUBGRID_KERNEL_TRILINEAR) return 1;
	/*CHANGE END - 20260710 */
	else                                        return drange_;

}
/*****************************************************************************
 *
 *  subgrid_get_lattice_index
 *
 *  Get indexes for neigbour lattice sites
 *
 *****************************************************************************/
void subgrid_get_lattice_index_range(double r0[3], int range, int nlocal[3], int* i_min, int* i_max, int* j_min, int* j_max, int* k_min, int* k_max)
{
	*i_min = imax(1, (int)floor(r0[X] - range));
	*i_max = imin(nlocal[X], (int)ceil(r0[X] + range));
	*j_min = imax(1, (int)floor(r0[Y] - range));
	*j_max = imin(nlocal[Y], (int)ceil(r0[Y] + range));
	*k_min = imax(1, (int)floor(r0[Z] - range));
	*k_max = imin(nlocal[Z], (int)ceil(r0[Z] + range));

}
/*****************************************************************************
 *
 *  subgrid_get_lattice_index
 *
 *  Get indexes for neigbour lattice sites
 *
 *****************************************************************************/
void subgrid_get_lattice_index_range_halo(double r0[3], int range, int nlocal[3], int* i_min, int* i_max, int* j_min, int* j_max, int* k_min, int* k_max)
{
	/* Allow sources in either halo ([-range, 0] or [nlocal+1, nlocal+range])
	 * to reach all their kernel neighbours. Destinations outside [1,nlocal]
	 * are wrapped by the caller. */
	*i_min = imax(1 - 2 * range, (int)floor(r0[X] - range));
	*i_max = imin(nlocal[X] + 2 * range, (int)ceil(r0[X] + range));
	*j_min = imax(1 - 2 * range, (int)floor(r0[Y] - range));
	*j_max = imin(nlocal[Y] + 2 * range, (int)ceil(r0[Y] + range));
	*k_min = imax(1 - 2 * range, (int)floor(r0[Z] - range));
	*k_max = imin(nlocal[Z] + 2 * range, (int)ceil(r0[Z] + range));

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
	rmod = fabs(r);
	// double s = (drange_ + 1.0) / 2.0;  /* Scale factor for mapping r to r/s */
	// rmod = fabs(r / s);
	/*CHANGE END - 20260418 Scale Peskin kernel with drange_ */

	if (rmod <= 1.0) {
		delta = 0.125 * (3.0 - 2.0 * rmod + sqrt(1.0 + 4.0 * rmod - 4.0 * rmod * rmod));
	}
	else if (rmod <= 2.0) {
		delta = 0.125 * (5.0 - 2.0 * rmod - sqrt(-7.0 + 12.0 * rmod - 4.0 * rmod * rmod));
	}

	/*CHANGE INIT - 20260418 Scale Peskin kernel with drange_ */
	return delta;
	// return delta / s;
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

	double r = fabs(x / drange_);
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

/*CHANGE INIT - 20260422 B-spline order-4 kernel */
/*****************************************************************************
 *
 *  d_bspline4
 *
 *  Cubic B-spline (order 4, degree 3) delta approximation.
 *  Support: [-2, 2] — same as Peskin with drange_=1.
 *
 *  Standard cubic B-spline (uniform knots, sum=1 on integer grid):
 *    |t| in [0,1): 2/3 - t^2 + t^3/2
 *    |t| in [1,2): (2 - |t|)^3 / 6
 *    |t| >= 2   : 0
 *
 *****************************************************************************/
double d_bspline4(double r) {

	double t = fabs(r);
	double val = 0.0;

	if (t < 1.0) {
		val = 2.0 / 3.0 - t * t + 0.5 * t * t * t;
	}
	else if (t < 2.0) {
		double u = 2.0 - t;
		val = u * u * u / 6.0;
	}

	return val;
}
/*CHANGE END - 20260422 B-spline order-4 kernel */

/*CHANGE INIT - 20260422 B-spline order-6 kernel */
/*****************************************************************************
 *
 *  d_bspline6
 *
 *  Quintic B-spline (order 6, degree 5) delta approximation.
 *  Support: [-3, 3] — 6-point stencil per dimension.
 *
 *  Derived from the standard de Boor recurrence; satisfies partition of
 *  unity on any integer-spaced grid. Coefficients verified numerically.
 *
 *  Let t = |r|:
 *    t in [0,1): 11/20 - t^2/2 + t^4/4 - t^5/12
 *    t in [1,2): Horner in u = t-1:
 *                13/60 + u*(-5/12 + u*(1/6 + u*(1/6 + u*(-1/6 + u*(1/24)))))
 *    t in [2,3): (3-t)^5 / 120
 *    t >= 3   : 0
 *
 *****************************************************************************/
double d_bspline6(double r) {

	double t = fabs(r);
	double val = 0.0;

	if (t < 1.0) {
		double t2 = t * t;
		val = 11.0 / 20.0 - t2 / 2.0 + t2 * t2 / 4.0 - t2 * t2 * t / 12.0;
	}
	else if (t < 2.0) {
		double u = t - 1.0;
		val = 13.0 / 60.0 + u * (-5.0 / 12.0 + u * (1.0 / 6.0 + u * (1.0 / 6.0 + u * (-1.0 / 6.0 + u * (1.0 / 24.0)))));
	}
	else if (t < 3.0) {
		double u = 3.0 - t;
		double u2 = u * u;
		val = u2 * u2 * u / 120.0;
	}

	return val;
}
/*CHANGE END - 20260422 B-spline order-6 kernel */

/*CHANGE INIT - 20260424 Peskin 6-point kernel (Bao et al. 2016) */
/*****************************************************************************
 *
 *  d_peskin6
 *
 *  New 6-point C^3 immersed-boundary kernel (Bao et al. 2016, Eq. 2.16-2.23).
 *  Support: [-3, 3].  Satisfies partition of unity on integer grid.
 *
 *  K = 59/60 - sqrt(29)/20 ≈ 0.7141.
 *  For r in [0,1): beta(r), gamma(r) defined, then phi computed at 6 nodes.
 *  Argument s = x_node - x_particle (signed distance).
 *
 *****************************************************************************/
/* Six phi formulas of the Bao et al. 2016 6-point kernel, evaluated at the
 * normalised fractional offset r in [0,1). segment in {0..5} selects the
 * weight for one of the 6 support nodes. */
static double peskin6_phi_segment(int segment, double r) {
	/* K = 59/60 - sqrt(29)/20 (Bao et al. 2016 Eq. 2.15) */
	static const double K = 0.71407520893979593;
	double beta = 9.0 / 4.0 - 1.5 * (K + r * r) + (22.0 / 3.0 - 7.0 * K) * r - (7.0 / 3.0) * r * r * r;
	double t1 = (3.0 * K - 1.0) * r + r * r * r;
	double t2 = (4.0 - 3.0 * K) * r - r * r * r;
	double gamma_r = -11.0 / 32.0 * r * r + 3.0 / 32.0 * (2.0 * K + r * r) * r * r
		+ t1 * t1 / 72.0 + t2 * t2 / 18.0;
	double phi_m3 = (-beta + sqrt(beta * beta - 112.0 * gamma_r)) / 56.0;

	switch (segment) {
	case 0: return phi_m3;
	case 1: return -3.0 * phi_m3 - 1.0 / 16.0 + (K + r * r) / 8.0
		+ (3.0 * K - 1.0) * r / 12.0 + r * r * r / 12.0;
	case 2: return  2.0 * phi_m3 + 1.0 / 4.0 + (4.0 - 3.0 * K) * r / 6.0 - r * r * r / 6.0;
	case 3: return  2.0 * phi_m3 + 5.0 / 8.0 - (K + r * r) / 4.0;
	case 4: return -3.0 * phi_m3 + 1.0 / 4.0 - (4.0 - 3.0 * K) * r / 6.0 + r * r * r / 6.0;
	case 5: return  phi_m3 - 1.0 / 16.0 + (K + r * r) / 8.0
		- (3.0 * K - 1.0) * r / 12.0 - r * r * r / 12.0;
	default: return 0.0;
	}
}

/*CHANGE INIT - 20260630 Fix d_peskin6 node-segment mapping.
 * Previous version derived (r, segment) per integer band of s with a
 * piecewise r definition, which broke partition of unity (Sum != 1) and
 * symmetry (d6(s) != d6(-s)) at integer offsets (particle exactly on a
 * node). That injected ~2.7% charge error wherever the particle / fluid
 * ions sit on integer nodes, producing the noisy Debye layer seen with
 * Peskin6. Correct mapping: particle at floor+off (off in [0,1)); the node
 * at floor+m (m in {-2..3}) gets phi_segment(3 - m, off). This gives
 * Sum = 1 and the correct first moment (centre of mass = off) for ALL
 * offsets, including integers, and is symmetric d6(s)=d6(-s). */
double d_peskin6(double s) {
	if (fabs(s) >= 3.0) return 0.0;

	double off = s - floor(s);          /* fractional offset in [0,1) */
	int    m   = (int)lround(off - s);  /* node index relative to floor */
	int    segment = 3 - m;             /* node floor+m -> segment 3-m */
	if (segment < 0 || segment > 5) return 0.0;

	return peskin6_phi_segment(segment, off);
}
/*CHANGE END - 20260630 Fix d_peskin6 node-segment mapping */
/*CHANGE END - 20260424 Peskin 6-point kernel */

/*CHANGE INIT - 20260424 Kaiser-Bessel order-4 kernel */
/*****************************************************************************
 *
 *  Kaiser-Bessel window kernel, support W=4 ([-2,2] per dimension).
 *
 *  phi(r) = I0(beta * sqrt(1 - (r/2)^2)) / (2 * I0(beta))   for |r| < 2
 *           0                                                  otherwise
 *
 *  I0 computed via Cephes polynomial approximation (Numerical Recipes).
 *  Two-region Chebyshev fit, error < 1e-7 for all x.
 *
 *****************************************************************************/
void  subgrid_set_kb4_beta(double beta) {
	subgrid_kb4_beta_ = beta;
	double i0b = i0_series(beta);
	double s = 0.0;
	for (int m = -2; m <= 2; m++) {
		double t = fabs((double)m);
		if (t < 2.0) s += i0_series(beta * sqrt(1.0 - (m / 2.0) * (m / 2.0))) / (2.0 * i0b);
	}
	subgrid_kb4_norm_fluid_ = s * s * s;
}

double subgrid_kb4_norm_fluid(void) {
	return subgrid_kb4_norm_fluid_;
}

static double i0_series(double x) {
	double ax = fabs(x);
	double y, ans;
	if (ax < 3.75) {
		y = x / 3.75;
		y *= y;
		ans = 1.0 + y * (3.5156229 + y * (3.0899424 + y * (1.2067492
			+ y * (0.2659732 + y * (0.0360768 + y * 0.0045813)))));
	}
	else {
		y = 3.75 / ax;
		ans = (exp(ax) / sqrt(ax)) * (0.39894228 + y * (0.01328592
			+ y * (0.00225319 + y * (-0.00157565 + y * (0.00916281
				+ y * (-0.02057706 + y * (0.02635537 + y * (-0.01647633
					+ y * 0.00392377))))))));
	}
	return ans;
}

double d_kb4(double r) {
	double t = fabs(r);
	if (t >= 2.0) return 0.0;
	double arg = sqrt(1.0 - (r / 2.0) * (r / 2.0));
	return i0_series(subgrid_kb4_beta_ * arg) / (2.0 * i0_series(subgrid_kb4_beta_));
}
/*CHANGE END - 20260424 Kaiser-Bessel order-4 kernel */

/*CHANGE INIT - 20260630 Hann (raised-cosine) spread/gather kernel */
/*****************************************************************************
 *
 *  subgrid_set_hann_order
 *
 *  Set the Hann kernel order n (full support width, lattice units). A larger
 *  n spreads the charge/force over more nodes, reducing the "snap-to-grid"
 *  effect at the cost of a wider stencil. n must be a positive integer;
 *  partition of unity holds for any integer n (odd n makes the active-node
 *  count vary between n and n-1 with the sub-grid offset, even n is uniform).
 *
 *****************************************************************************/
void subgrid_set_hann_order(double n) {
	/* n must be a positive integer for exact partition of unity (odd n is fine;
	 * it only makes the active-node count vary with the offset). Callers that
	 * read from input validate and report this; the assert guards direct use. */
	assert(n > 0.0);
	assert(n == floor(n));
	subgrid_hann_order_ = n;
}

/*****************************************************************************
 *
 *  d_hann
 *
 *  Raised-cosine (Hann) window used as a regularised delta for spread/gather:
 *
 *      w(r) = (1/n) [1 + cos(2*pi*r/n)],   |r| <= n/2
 *           = 0                            otherwise
 *
 *  where r = x_node - x_particle (signed distance in lattice units) and n is
 *  the order (support width). For integer n the window satisfies partition of
 *  unity exactly: sum_i w(i - c) = 1 for any offset c. This holds because the
 *  Fourier transform of the truncated raised cosine has exact zeros at all
 *  non-zero integer frequencies, so by Poisson summation the node sum is
 *  independent of c. (Verified numerically to ~1e-16 for n in {2,3,4,5,6,8}.)
 *  Hence no per-particle renormalisation is required, unlike KB4.
 *
 *  Reference: raised-cosine interpolation window; reduces grid-snapping of
 *  the particle<->mesh coupling (see ek_particle_coupling note).
 *
 *****************************************************************************/
double d_hann(double r) {
	double n = subgrid_hann_order_;
	double t = fabs(r);
	if (t >= 0.5 * n) return 0.0;
	PI_DOUBLE(pi);
	return (1.0 / n) * (1.0 + cos(2.0 * pi * r / n));
}
/*CHANGE END - 20260630 Hann spread/gather kernel */

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

/*CHANGE INIT - 20260422 subgrid_get_drange */
double subgrid_get_drange(void) {
	return drange_;
}
/*CHANGE END - 20260422 subgrid_get_drange */

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
	peskin_scatter_nodes_kernel << <(n_interior + threads - 1) / threads, threads >> > (
		rho_src_d, rho_dst_d, nx, ny, nz, nhalo, nsites);

	/* Build flat particle array and scatter on GPU */
	int ncell[3];
	colloids_info_ncell(cinfo, ncell);

	cudaDeviceSynchronize();

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
		peskin_scatter_particles_kernel << <(npart + threads - 1) / threads, threads >> > (
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
	int nk = psi_src->nk;

	/* Scatter LB-node charges */
	for (int i = 1; i <= nlocal[X]; i++) {
		for (int j = 1; j <= nlocal[Y]; j++) {
			for (int k = 1; k <= nlocal[Z]; k++) {
				int src = cs_index(psi_src->cs, i, j, k);
				double rho0, rho1;
				psi_rho(psi_src, src, 0, &rho0);
				psi_rho(psi_src, src, 1, &rho1);
				int i_min = imax(0, i - 1), i_max = imin(nlocal[X] + 1, i + 1);
				int j_min = imax(0, j - 1), j_max = imin(nlocal[Y] + 1, j + 1);
				int k_min = imax(0, k - 1), k_max = imin(nlocal[Z] + 1, k + 1);
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
						double r0[3] = { p->s.r[X] - 1.0 * offset[X],
										 p->s.r[Y] - 1.0 * offset[Y],
										 p->s.r[Z] - 1.0 * offset[Z] };
						int i_min = imax(0, (int)floor(r0[X] - 1.0));
						int i_max = imin(nlocal[X] + 1, (int)ceil(r0[X] + 1.0));
						int j_min = imax(0, (int)floor(r0[Y] - 1.0));
						int j_max = imin(nlocal[Y] + 1, (int)ceil(r0[Y] + 1.0));
						int k_min = imax(0, (int)floor(r0[Z] - 1.0));
						int k_max = imin(nlocal[Z] + 1, (int)ceil(r0[Z] + 1.0));
						for (int i = i_min; i <= i_max; i++)
							for (int j = j_min; j <= j_max; j++)
								for (int k = k_min; k <= k_max; k++) {
									double dr = d_peskin(r0[X] - i) * d_peskin(r0[Y] - j) * d_peskin(r0[Z] - k);
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

/*CHANGE INIT - 20260625 PM short-range correction table */

/*****************************************************************************
 *
 *  pm_sr_table_t
 *
 *  Lookup table for the short-range correction to the PM (FFT/PETSc) solver.
 *
 *  Stores delta_phi(r), delta_E(r), delta_F(r) = phi_ref - phi_PM_p, etc.,
 *  where the reference is a Gaussian charge distribution of width sigma and
 *  phi_PM_p is estimated from the same Gaussian smeared onto the lattice.
 *
 *  Two build strategies:
 *    pm_sr_table_build_radial  -- phi_PM_p averaged over all directions (1D)
 *    pm_sr_table_build_3d      -- phi_PM_p evaluated on a 3D grid (dx,dy,dz)
 *
 *  Parameters common to both:
 *    sigma    - Gaussian width (lattice units)
 *    r_cut    - correction cutoff radius (lattice units)
 *    n_voxel  - sub-voxel points per axis for averaging within each cell
 *    epsilon  - dielectric permittivity
 *
 *****************************************************************************/

#define PM_SR_NR  1024   /* radial bins for 1D table  */

struct pm_sr_table_s {
  int    nr;       /* radial bins (1D) or points per axis (3D: nr^3 total) */
  int    is_3d;    /* 0 = radial, 1 = 3D                                   */
  double r_cut;
  double sigma;
  double dr;       /* bin width: r_cut/nr (1D) or 2*r_cut/nr per axis (3D) */
  double* phi;     /* delta_phi table                                       */
  double* E;       /* delta_E: radial magnitude (1D only)                  */
  double* F;       /* delta_F = delta_E (1D only)                          */
  /* CHANGE 20260724: 3D VECTOR field correction. delta_E is a vector
   * (E_ref*r_hat - E_pm_vec), stored as 3 planes so the cubic anisotropy of
   * the kernel (which the radial table cannot represent) is captured. The
   * old scalar "E_pm_x * r/dx" reconstruction was singular at dx=0 and wrong
   * off-axis, so the 3D table never actually differed from the radial one. */
  double* Ex;      /* delta_E x-component, nr^3 (3D only) */
  double* Ey;      /* delta_E y-component, nr^3 (3D only) */
  double* Ez;      /* delta_E z-component, nr^3 (3D only) */
};

/* -------------------------------------------------------------------------
 * Internal: 1D kernel weight for the spread/gather kernel used to scatter
 * the particle charge. Mirrors the selection in subgrid_charge_from_particles
 * and subgrid_update_forces_electrokinetics.
 * ------------------------------------------------------------------------- */
static double pm_sr_kernel_weight_1d(double dr, subgrid_kernel_t kernel) {
  switch (kernel) {
  case SUBGRID_KERNEL_BSPLINE6: return d_bspline6(dr);
  case SUBGRID_KERNEL_BSPLINE4: return d_bspline4(dr);
  case SUBGRID_KERNEL_PESKIN6:  return d_peskin6(dr);
  /*CHANGE INIT - 20260630 Hann kernel */
  case SUBGRID_KERNEL_HANN:     return d_hann(dr);
  /*CHANGE END - 20260630 */
  case SUBGRID_KERNEL_KB4:      return d_kb4(dr);
  case SUBGRID_KERNEL_PESKIN4:
  default:                      return d_peskin(dr);
  }
}

/* -------------------------------------------------------------------------
 * Internal: analytic Gaussian phi and E (unit charge, epsilon given)
 * ------------------------------------------------------------------------- */
static double pm_sr_phi_gauss(double r, double sigma, double epsilon) {
  PI_DOUBLE(pi);
  if (r < 1.0e-12) return 0.0;
  return erf(r / (sigma * sqrt(2.0))) / (4.0 * pi * epsilon * r);
}

static double pm_sr_E_gauss(double r, double sigma, double epsilon) {
  PI_DOUBLE(pi);
  if (r < 1.0e-12) return 0.0;
  double r2   = r * r;
  double s2   = sigma * sqrt(2.0);
  double erfv = erf(r / s2);
  double gaus = sqrt(2.0 / pi) / sigma * exp(-r2 / (2.0 * sigma * sigma));
  return (erfv / r2 - gaus / r) / (4.0 * pi * epsilon);
}

/* -------------------------------------------------------------------------
 * Internal: phi_PM_p(r) along direction r_hat=(1,0,0).
 * Charge q=1 is smeared with trilinear weights onto 8 surrounding nodes
 * from fractional position (xf, yf, zf). The Coulomb potential of each
 * smeared node charge is evaluated at (xf + r, yf, zf).
 * Average over n_voxel^3 fractional positions within the unit cell.
 * ------------------------------------------------------------------------- */
static double pm_sr_phi_PM_radial(double r, double sigma, double epsilon,
                                   int n_voxel, subgrid_kernel_t kernel) {
  PI_DOUBLE(pi);
  int krange = subgrid_get_range(kernel);
  double sum = 0.0;
  for (int ix = 0; ix < n_voxel; ix++) {
    double xf = (ix + 0.5) / n_voxel;
    for (int iy = 0; iy < n_voxel; iy++) {
      double yf = (iy + 0.5) / n_voxel;
      for (int iz = 0; iz < n_voxel; iz++) {
        double zf = (iz + 0.5) / n_voxel;
        /* evaluation point at distance r along x from particle (xf,yf,zf) */
        double ex = xf + r, ey = yf, ez = zf;
        /* Smear unit charge with the real kernel over its full support */
        double phi = 0.0;
        for (int di = -krange; di <= krange + 1; di++)
          for (int dj = -krange; dj <= krange + 1; dj++)
            for (int dk = -krange; dk <= krange + 1; dk++) {
              double w = pm_sr_kernel_weight_1d(xf - di, kernel)
                       * pm_sr_kernel_weight_1d(yf - dj, kernel)
                       * pm_sr_kernel_weight_1d(zf - dk, kernel);
              if (w == 0.0) continue;
              double dist = sqrt((ex-di)*(ex-di) + (ey-dj)*(ey-dj) +
                                 (ez-dk)*(ez-dk));
              if (dist < 1.0e-12) continue;
              phi += w / (4.0*pi*epsilon*dist);
            }
        sum += phi;
      }
    }
  }
  return sum / (n_voxel * n_voxel * n_voxel);
}

/* Same but returns radial E component (x-axis) */
static double pm_sr_E_PM_radial(double r, double sigma, double epsilon,
                                 int n_voxel, subgrid_kernel_t kernel) {
  PI_DOUBLE(pi);
  int krange = subgrid_get_range(kernel);
  double sum = 0.0;
  for (int ix = 0; ix < n_voxel; ix++) {
    double xf = (ix + 0.5) / n_voxel;
    for (int iy = 0; iy < n_voxel; iy++) {
      double yf = (iy + 0.5) / n_voxel;
      for (int iz = 0; iz < n_voxel; iz++) {
        double zf = (iz + 0.5) / n_voxel;
        double ex = xf + r, ey = yf, ez = zf;
        double E = 0.0;
        for (int di = -krange; di <= krange + 1; di++)
          for (int dj = -krange; dj <= krange + 1; dj++)
            for (int dk = -krange; dk <= krange + 1; dk++) {
              double w = pm_sr_kernel_weight_1d(xf - di, kernel)
                       * pm_sr_kernel_weight_1d(yf - dj, kernel)
                       * pm_sr_kernel_weight_1d(zf - dk, kernel);
              if (w == 0.0) continue;
              double dx = ex - di, dy = ey - dj, dz = ez - dk;
              double d2 = dx*dx + dy*dy + dz*dz;
              if (d2 < 1.0e-24) continue;
              double d3 = sqrt(d2) * d2;
              E += w * dx / (4.0*pi*epsilon*d3);
            }
        sum += E;
      }
    }
  }
  return sum / (n_voxel * n_voxel * n_voxel);
}

/*****************************************************************************
 *
 *  pm_sr_table_build_radial
 *
 *  Build 1D radial table: phi_PM_p averaged over fractional positions
 *  (all along the same direction). Fast, approximate.
 *
 *****************************************************************************/
int pm_sr_table_build_radial(double sigma, double r_cut, int n_voxel,
                              double epsilon, subgrid_kernel_t kernel,
                              pm_sr_table_t** ptable) {
  assert(ptable && sigma > 0.0 && r_cut > 0.0 && n_voxel > 0 && epsilon > 0.0);

  pm_sr_table_t* t = (pm_sr_table_t*)calloc(1, sizeof(pm_sr_table_t));
  if (!t) return -1;

  t->nr    = PM_SR_NR;
  t->is_3d = 0;
  t->r_cut = r_cut;
  t->sigma = sigma;
  t->dr    = r_cut / PM_SR_NR;
  t->phi   = (double*)malloc(PM_SR_NR * sizeof(double));
  t->E     = (double*)malloc(PM_SR_NR * sizeof(double));
  t->F     = (double*)malloc(PM_SR_NR * sizeof(double));
  if (!t->phi || !t->E || !t->F) { pm_sr_table_free(&t); return -1; }

  for (int i = 0; i < PM_SR_NR; i++) {
    double r    = (i + 0.5) * t->dr;
    double dphi = pm_sr_phi_gauss(r, sigma, epsilon)
                - pm_sr_phi_PM_radial(r, sigma, epsilon, n_voxel, kernel);
    double dE   = pm_sr_E_gauss(r, sigma, epsilon)
                - pm_sr_E_PM_radial(r, sigma, epsilon, n_voxel, kernel);
    t->phi[i] = dphi;
    t->E[i]   = dE;
    t->F[i]   = dE;   /* per unit q_node */
  }

  *ptable = t;
  return 0;
}

/*****************************************************************************
 *
 *  pm_sr_table_build_3d
 *
 *  Build 3D table: phi_PM_p evaluated for each (dx, dy, dz) on a grid of
 *  nr^3 points in [-r_cut, r_cut]^3. Exact per direction, larger memory.
 *
 *  Layout: index = ix*nr*nr + iy*nr + iz
 *  Coordinate: dx = -r_cut + (ix + 0.5)*dr, same for dy, dz.
 *
 *****************************************************************************/
int pm_sr_table_build_3d(double sigma, double r_cut, int nr, int n_voxel,
                          double epsilon, subgrid_kernel_t kernel,
                          pm_sr_table_t** ptable) {
  assert(ptable && sigma > 0.0 && r_cut > 0.0 && nr > 0
         && n_voxel > 0 && epsilon > 0.0);
  int krange = subgrid_get_range(kernel);

  pm_sr_table_t* t = (pm_sr_table_t*)calloc(1, sizeof(pm_sr_table_t));
  if (!t) return -1;

  t->nr    = nr;
  t->is_3d = 1;
  t->r_cut = r_cut;
  t->sigma = sigma;
  t->dr    = 2.0 * r_cut / nr;   /* axis spacing */
  int ntot = nr * nr * nr;
  t->phi   = (double*)malloc(ntot * sizeof(double));
  t->Ex    = (double*)malloc(ntot * sizeof(double));
  t->Ey    = (double*)malloc(ntot * sizeof(double));
  t->Ez    = (double*)malloc(ntot * sizeof(double));
  if (!t->phi || !t->Ex || !t->Ey || !t->Ez) { pm_sr_table_free(&t); return -1; }

  PI_DOUBLE(pi);

  for (int ix = 0; ix < nr; ix++) {
    double dx = -r_cut + (ix + 0.5) * t->dr;
    for (int iy = 0; iy < nr; iy++) {
      double dy = -r_cut + (iy + 0.5) * t->dr;
      for (int iz = 0; iz < nr; iz++) {
        double dz = -r_cut + (iz + 0.5) * t->dr;
        double r  = sqrt(dx*dx + dy*dy + dz*dz);
        int idx   = ix*nr*nr + iy*nr + iz;

        /* Reference: Gaussian analytic (radial, exact) */
        double phi_ref = pm_sr_phi_gauss(r, sigma, epsilon);
        double E_ref   = pm_sr_E_gauss(r, sigma, epsilon);

        /* CHANGE 20260724b: NO sub-voxel averaging. The table index (dx,dy,dz)
         * IS the exact particle->node displacement used by the application
         * (dx = r0 - i, carrying the particle's fractional sub-cell position).
         * So the PM field is computed for that SINGLE displacement, not
         * averaged over offsets — this makes the table offset-aware.
         *
         * Geometry: the EVALUATION NODE is at the origin. The PARTICLE is at
         * (dx,dy,dz) relative to it (dx = particle - node). The particle's
         * charge is spread onto the integer nodes (m,n,l) NEAR THE PARTICLE
         * with weight kernel(particle - spread_node) = kernel(dx - m). Each
         * fragment w at node (m,n,l) contributes its Coulomb field measured at
         * the origin: direction (origin - node) = -(m,n,l). This offset-aware
         * PM field is what the radial table destroyed by averaging over xf. */
        double phi_pm = 0.0, E_pm_x = 0.0, E_pm_y = 0.0, E_pm_z = 0.0;
        int mlo_x = (int)floor(dx) - krange, mhi_x = (int)floor(dx) + krange + 1;
        int mlo_y = (int)floor(dy) - krange, mhi_y = (int)floor(dy) + krange + 1;
        int mlo_z = (int)floor(dz) - krange, mhi_z = (int)floor(dz) + krange + 1;
        for (int m = mlo_x; m <= mhi_x; m++)
          for (int n = mlo_y; n <= mhi_y; n++)
            for (int l = mlo_z; l <= mhi_z; l++) {
              double w = pm_sr_kernel_weight_1d(dx - m, kernel)
                       * pm_sr_kernel_weight_1d(dy - n, kernel)
                       * pm_sr_kernel_weight_1d(dz - l, kernel);
              if (w == 0.0) continue;
              /* field of fragment at (m,n,l) measured at the origin node */
              double vx = -(double)m, vy = -(double)n, vz = -(double)l;
              double d2 = vx*vx + vy*vy + vz*vz;
              if (d2 < 1.0e-24) continue;   /* fragment exactly on the node */
              double d3 = sqrt(d2) * d2;
              double c = w / (4.0*pi*epsilon);
              phi_pm += w / (4.0*pi*epsilon*sqrt(d2));
              E_pm_x += c * vx / d3;
              E_pm_y += c * vy / d3;
              E_pm_z += c * vz / d3;
            }

        /* Reference field at the origin from the particle at (dx,dy,dz):
         * points from particle to node = -(dx,dy,dz)/r. Both E_ref_vec and
         * E_pm_vec are the field AT THE NODE, so the correction stored is
         * delta_E = E_ref_vec - E_pm_vec, consistent frame. */
        double rhx = (r > 1.0e-12) ? -dx / r : 0.0;
        double rhy = (r > 1.0e-12) ? -dy / r : 0.0;
        double rhz = (r > 1.0e-12) ? -dz / r : 0.0;

        t->phi[idx] = phi_ref - phi_pm;
        t->Ex[idx]  = E_ref * rhx - E_pm_x;
        t->Ey[idx]  = E_ref * rhy - E_pm_y;
        t->Ez[idx]  = E_ref * rhz - E_pm_z;
      }
    }
  }

  *ptable = t;
  return 0;
}

/*****************************************************************************
 *
 *  pm_sr_table_free
 *
 *****************************************************************************/
void pm_sr_table_free(pm_sr_table_t** ptable) {
  if (!ptable || !*ptable) return;
  free((*ptable)->phi);
  free((*ptable)->E);   /* NULL for 3D (calloc), free(NULL) is safe */
  free((*ptable)->F);
  free((*ptable)->Ex);  /* NULL for 1D */
  free((*ptable)->Ey);
  free((*ptable)->Ez);
  free(*ptable);
  *ptable = NULL;
}

/*****************************************************************************
 *
 *  pm_sr_table_interpolate
 *
 *  Returns delta_phi, delta_E (radial), delta_F (radial, per unit q_node)
 *  at distance r for the 1D table.
 *  Returns 1 if r < r_cut, 0 otherwise (no correction).
 *
 *****************************************************************************/
int pm_sr_table_interpolate(const pm_sr_table_t* t, double r,
                             double* delta_phi, double* delta_E,
                             double* delta_F) {
  assert(t && !t->is_3d);
  *delta_phi = *delta_E = *delta_F = 0.0;
  if (r >= t->r_cut || r < 0.0) return 0;

  double pos  = r / t->dr;
  int    i0   = (int)floor(pos);
  int    i1   = i0 + 1;
  double frac = pos - i0;
  if (i1 >= t->nr) i1 = t->nr - 1;

  *delta_phi = (1.0-frac)*t->phi[i0] + frac*t->phi[i1];
  *delta_E   = (1.0-frac)*t->E[i0]   + frac*t->E[i1];
  *delta_F   = (1.0-frac)*t->F[i0]   + frac*t->F[i1];
  return 1;
}

/*****************************************************************************
 *
 *  pm_sr_table_interpolate_3d
 *
 *  Trilinear interpolation in the 3D table at displacement (dx, dy, dz).
 *  Returns delta_phi (scalar) and the VECTOR field correction
 *  (delta_Ex, delta_Ey, delta_Ez). Returns 1 if inside r_cut, 0 otherwise.
 *
 *****************************************************************************/
int pm_sr_table_interpolate_3d(const pm_sr_table_t* t,
                                double dx, double dy, double dz,
                                double* delta_phi,
                                double* delta_Ex, double* delta_Ey,
                                double* delta_Ez) {
  assert(t && t->is_3d);
  *delta_phi = *delta_Ex = *delta_Ey = *delta_Ez = 0.0;

  double r = sqrt(dx*dx + dy*dy + dz*dz);
  if (r >= t->r_cut || r < 1.0e-12) return 0;

  int    nr  = t->nr;
  double dr  = t->dr;
  double off = t->r_cut;   /* origin of axis: index 0 corresponds to -r_cut */

  double px = (dx + off) / dr - 0.5;
  double py = (dy + off) / dr - 0.5;
  double pz = (dz + off) / dr - 0.5;

  int ix = (int)floor(px),  iy = (int)floor(py),  iz = (int)floor(pz);
  double fx = px - ix,      fy = py - iy,          fz = pz - iz;

  /* clamp */
  int ix1 = ix+1, iy1 = iy+1, iz1 = iz+1;
  if (ix  <  0) { ix  = 0; fx = 0.0; }
  if (iy  <  0) { iy  = 0; fy = 0.0; }
  if (iz  <  0) { iz  = 0; fz = 0.0; }
  if (ix1 >= nr) ix1 = nr-1;
  if (iy1 >= nr) iy1 = nr-1;
  if (iz1 >= nr) iz1 = nr-1;

#define IDX3(a,b,c) ((a)*nr*nr + (b)*nr + (c))
#define TRILIN(FIELD) ( \
    (1-fx)*(1-fy)*(1-fz)*t->FIELD[IDX3(ix ,iy ,iz )] + \
    (1-fx)*(1-fy)*   fz *t->FIELD[IDX3(ix ,iy ,iz1)] + \
    (1-fx)*   fy *(1-fz)*t->FIELD[IDX3(ix ,iy1,iz )] + \
    (1-fx)*   fy *   fz *t->FIELD[IDX3(ix ,iy1,iz1)] + \
       fx *(1-fy)*(1-fz)*t->FIELD[IDX3(ix1,iy ,iz )] + \
       fx *(1-fy)*   fz *t->FIELD[IDX3(ix1,iy ,iz1)] + \
       fx *   fy *(1-fz)*t->FIELD[IDX3(ix1,iy1,iz )] + \
       fx *   fy *   fz *t->FIELD[IDX3(ix1,iy1,iz1)] )

  *delta_phi = TRILIN(phi);
  *delta_Ex  = TRILIN(Ex);
  *delta_Ey  = TRILIN(Ey);
  *delta_Ez  = TRILIN(Ez);
#undef TRILIN
#undef IDX3
  return 1;
}

/*****************************************************************************
 *
 *  pm_sr_table_print
 *
 *  Print the table to fp (stdout if NULL).
 *  1D: columns r  delta_phi  delta_E  delta_F
 *  3D: columns dx  dy  dz  r  delta_phi  delta_E  delta_F
 *
 *****************************************************************************/
void pm_sr_table_print(const pm_sr_table_t* t, FILE* fp) {
  assert(t);
  FILE* out = fp ? fp : stdout;
  fprintf(out, "# PM short-range correction table  is_3d=%d\n", t->is_3d);
  fprintf(out, "# sigma=%.6g  r_cut=%.6g  nr=%d  dr=%.6g\n",
          t->sigma, t->r_cut, t->nr, t->dr);

  if (!t->is_3d) {
    fprintf(out, "# r  delta_phi  delta_E  delta_F\n");
    for (int i = 0; i < t->nr; i++) {
      double r = (i + 0.5) * t->dr;
      fprintf(out, "%.8e  %.8e  %.8e  %.8e\n",
              r, t->phi[i], t->E[i], t->F[i]);
    }
  } else {
    int nr = t->nr;
    fprintf(out, "# dx  dy  dz  r  delta_phi  delta_Ex  delta_Ey  delta_Ez\n");
    for (int ix = 0; ix < nr; ix++) {
      double dx = -t->r_cut + (ix+0.5)*t->dr;
      for (int iy = 0; iy < nr; iy++) {
        double dy = -t->r_cut + (iy+0.5)*t->dr;
        for (int iz = 0; iz < nr; iz++) {
          double dz = -t->r_cut + (iz+0.5)*t->dr;
          double r  = sqrt(dx*dx + dy*dy + dz*dz);
          int idx   = ix*nr*nr + iy*nr + iz;
          fprintf(out, "%.8e  %.8e  %.8e  %.8e  %.8e  %.8e  %.8e  %.8e\n",
                  dx, dy, dz, r, t->phi[idx],
                  t->Ex[idx], t->Ey[idx], t->Ez[idx]);
        }
      }
    }
  }
}

/*CHANGE END - 20260625 PM short-range correction table */

/*CHANGE INIT - 20260625 pm_sr_correct_phi / pm_sr_correct_efield */

/*****************************************************************************
 *
 *  pm_sr_correct_phi
 *
 *  Step 1 of option C:
 *    - If phi_pm_buf != NULL, copy phi_PM into it before modifying psi->psi
 *    - Add q_p * delta_phi(r) to psi->psi for each node within r_cut
 *
 *****************************************************************************/
int pm_sr_correct_phi(colloids_info_t* cinfo, psi_t* psi,
                       const pm_sr_table_t* table, double* phi_pm_buf) {

  int ncell[3], nlocal[3], offset[3];

  assert(cinfo && psi && table);
  if (cinfo->nsubgrid == 0) return 0;

  cs_nlocal(cinfo->cs, nlocal);
  cs_nlocal_offset(cinfo->cs, offset);
  colloids_info_ncell(cinfo, ncell);

  /* Optionally save phi_PM before correction (caller may pass NULL) */
  if (phi_pm_buf)
    memcpy(phi_pm_buf, psi->psi->data, psi->nsites * sizeof(double));

  double eunit;
  psi_unit_charge(psi, &eunit);
  int range = (int)ceil(table->r_cut);

  for (int ic = 1; ic <= ncell[X]; ic++)
    for (int jc = 1; jc <= ncell[Y]; jc++)
      for (int kc = 1; kc <= ncell[Z]; kc++) {
        colloid_t* pc = NULL;
        colloids_info_cell_list_head(cinfo, ic, jc, kc, &pc);
        for (; pc; pc = pc->next) {
          if (pc->s.bc != COLLOID_BC_SUBGRID) continue;

          double q_p = (pc->s.q0 - pc->s.q1) * eunit;
          double r0[3] = { pc->s.r[X] - offset[X],
                           pc->s.r[Y] - offset[Y],
                           pc->s.r[Z] - offset[Z] };
          int i_min, i_max, j_min, j_max, k_min, k_max;
          subgrid_get_lattice_index_range(r0, range, nlocal,
                                          &i_min, &i_max,
                                          &j_min, &j_max,
                                          &k_min, &k_max);

          for (int i = i_min; i <= i_max; i++)
            for (int j = j_min; j <= j_max; j++)
              for (int k = k_min; k <= k_max; k++) {
                double dx = r0[X] - i, dy = r0[Y] - j, dz = r0[Z] - k;
                double r  = sqrt(dx*dx + dy*dy + dz*dz);
                double dphi = 0.0, dE = 0.0, dF = 0.0;
                double dEx = 0.0, dEy = 0.0, dEz = 0.0;
                int hit = table->is_3d
                  ? pm_sr_table_interpolate_3d(table, dx, dy, dz, &dphi, &dEx, &dEy, &dEz)
                  : pm_sr_table_interpolate(table, r, &dphi, &dE, &dF);
                if (!hit) continue;
                int idx = cs_index(cinfo->cs, i, j, k);
                double phi_old;
                psi_psi(psi, idx, &phi_old);
                psi_psi_set(psi, idx, phi_old + q_p * dphi);
              }
        }
      }

  return 0;
}

/*****************************************************************************
 *
 *  pm_sr_apply_force_correction
 *
 *  Apply the short-range force/field correction AFTER the uncorrected forces
 *  have been computed (psi_force_gradmu, subgrid_update_Esub,
 *  subgrid_update_forces_electrokinetics).
 *
 *  Two passes per particle to conserve momentum exactly:
 *
 *    Pass 1: sum the total correction force on the particle from all nodes j
 *            within r_cut (q_node = local ionic charge, q_p = particle charge,
 *            delta_F(r) the tabulated theory-PM force per unit of both charges,
 *            r_hat = (dx,dy,dz)/r):
 *
 *              dF_total += q_p * q_node * delta_F(r) * r_hat
 *
 *            and correct the diagnostic fields at the same time:
 *              psi->efield[j] += q_p   * delta_F(r) * r_hat * kt   (= dF/q_node)
 *              pc->Esub       += q_node* delta_F(r) * r_hat        (= dF/q_p)
 *
 *    Pass 2: distribute -dF_total onto the fluid with the SAME Peskin kernel
 *            used elsewhere, over the kernel support of the particle. Since
 *            the kernel is a partition of unity, the fluid receives exactly
 *            -dF_total, equal and opposite to the particle force.
 *
 *  Field corrections feed no force calculation (diagnostic only), so they do
 *  not affect momentum conservation.
 *
 *****************************************************************************/
int pm_sr_apply_force_correction(colloids_info_t* cinfo, map_t* map,
                                  psi_t* psi, hydro_t* hydro,
                                  const pm_sr_table_t* table,
                                  subgrid_kernel_t kernel) {

  int ncell[3], nlocal[3], offset[3];

  assert(cinfo && psi && hydro && table);
  if (cinfo->nsubgrid == 0) return 0;

  cs_nlocal(cinfo->cs, nlocal);
  cs_nlocal_offset(cinfo->cs, offset);
  colloids_info_ncell(cinfo, ncell);

  double eunit, beta;
  psi_unit_charge(psi, &eunit);
  psi_beta(psi, &beta);
  double kt     = 1.0 / beta;
  double reunit = 1.0 / eunit;
  int range  = (int)ceil(table->r_cut);
  /* Reaction force distributed with the SAME kernel used for charge spread
   * (partition of unity), matching subgrid_update_forces_electrokinetics. */
  int krange = subgrid_get_range(kernel);

  hydro_memcpy(hydro, tdpMemcpyDeviceToHost);

  for (int ic = 1; ic <= ncell[X]; ic++)
    for (int jc = 1; jc <= ncell[Y]; jc++)
      for (int kc = 1; kc <= ncell[Z]; kc++) {
        colloid_t* pc = NULL;
        colloids_info_cell_list_head(cinfo, ic, jc, kc, &pc);
        for (; pc; pc = pc->next) {
          if (pc->s.bc != COLLOID_BC_SUBGRID) continue;

          /* Particle charge in solver units (same as the Poisson source) */
          double q_p_solver = pc->s.q0 - pc->s.q1;
          /* Force prefactor: identical to subgrid_update_forces_electrokinetics
           * which does force = kt * reunit * Esub * (q0-q1). Here the field
           * correction delta_E is already in solver-potential units, and the
           * node charge below is also in solver units, so:
           *   dF = kt * reunit * (q_node_solver * delta_E) * (q0-q1) */
          double fpref = kt * reunit * q_p_solver;

          double r0[3] = { pc->s.r[X] - offset[X],
                           pc->s.r[Y] - offset[Y],
                           pc->s.r[Z] - offset[Z] };

          /* ---- Pass 1: total correction force on the particle ---- */
          double dF_total[3] = { 0.0, 0.0, 0.0 };

          int i_min, i_max, j_min, j_max, k_min, k_max;
          subgrid_get_lattice_index_range(r0, range, nlocal,
                                          &i_min, &i_max,
                                          &j_min, &j_max,
                                          &k_min, &k_max);

          for (int i = i_min; i <= i_max; i++)
            for (int j = j_min; j <= j_max; j++)
              for (int k = k_min; k <= k_max; k++) {
                double dx = r0[X] - i, dy = r0[Y] - j, dz = r0[Z] - k;
                double r  = sqrt(dx*dx + dy*dy + dz*dz);
                if (r < 1.0e-12) continue;

                /* Correction field VECTOR (dEcx,dEcy,dEcz). For the 3D table
                 * it comes straight from the interpolation (captures kernel
                 * anisotropy). For the radial table it is the scalar dE
                 * projected onto r_hat. */
                double dphi = 0.0, dEcx = 0.0, dEcy = 0.0, dEcz = 0.0;
                int hit;
                if (table->is_3d) {
                  hit = pm_sr_table_interpolate_3d(table, dx, dy, dz,
                                                   &dphi, &dEcx, &dEcy, &dEcz);
                } else {
                  double dE = 0.0, dF = 0.0;
                  hit = pm_sr_table_interpolate(table, r, &dphi, &dE, &dF);
                  double rx = dx/r, ry = dy/r, rz = dz/r;
                  dEcx = dE * rx; dEcy = dE * ry; dEcz = dE * rz;
                }
                if (!hit) continue;

                int idx = cs_index(cinfo->cs, i, j, k);

                /* Node charge in solver units (Poisson source = rho0 - rho1) */
                double rho0, rho1;
                psi_rho(psi, idx, 0, &rho0);
                psi_rho(psi, idx, 1, &rho1);
                double q_node = rho0 - rho1;

                /* Correction field at the particle from this node (vector) */
                double dEsx = q_node * dEcx;
                double dEsy = q_node * dEcy;
                double dEsz = q_node * dEcz;

                /* Pair force particle <-> node j (same prefactor as base force) */
                double dF_j[3];
                dF_j[X] = fpref * dEsx;
                dF_j[Y] = fpref * dEsy;
                dF_j[Z] = fpref * dEsz;

                /* Accumulate force on the particle */
                dF_total[X] += dF_j[X];
                dF_total[Y] += dF_j[Y];
                dF_total[Z] += dF_j[Z];

                /* Newton III: counter-force on the ion at node j (directly) */
                double f_node[3] = { -dF_j[X], -dF_j[Y], -dF_j[Z] };
                hydro_f_local_add(hydro, idx, f_node);

                /* diagnostic field at node (stored as kt*E): field from the
                 * particle charge q_p_solver seen at this node */
                double eb = kt * q_p_solver;
                psi->efield->data[addr_rank1(psi->efield->nsites,
                                             psi->efield->nf, idx, X)] += eb * dEcx;
                psi->efield->data[addr_rank1(psi->efield->nsites,
                                             psi->efield->nf, idx, Y)] += eb * dEcy;
                psi->efield->data[addr_rank1(psi->efield->nsites,
                                             psi->efield->nf, idx, Z)] += eb * dEcz;

                /* diagnostic effective field at particle (solver units) */
                pc->Esub[X] += dEsx;
                pc->Esub[Y] += dEsy;
                pc->Esub[Z] += dEsz;
              }

          /* Force on the particle */
          pc->fex[X] += dF_total[X];
          pc->fex[Y] += dF_total[Y];
          pc->fex[Z] += dF_total[Z];

          /* ---- Pass 2: communicate the particle force to the fluid with
           * Peskin, SAME sign as the particle (as subgrid_update_forces_
           * electrokinetics does for the electrostatic force). The Peskin
           * kernel is a partition of unity, so the sum over the support
           * equals +dF_total, exactly cancelling the -dF_j counter-forces
           * deposited node-by-node in Pass 1. Net fluid momentum from this
           * correction is therefore zero. ---- */
          subgrid_get_lattice_index_range(r0, krange, nlocal,
                                          &i_min, &i_max,
                                          &j_min, &j_max,
                                          &k_min, &k_max);

          for (int i = i_min; i <= i_max; i++)
            for (int j = j_min; j <= j_max; j++)
              for (int k = k_min; k <= k_max; k++) {
                double dr = pm_sr_kernel_weight_1d(r0[X]-i, kernel)
                          * pm_sr_kernel_weight_1d(r0[Y]-j, kernel)
                          * pm_sr_kernel_weight_1d(r0[Z]-k, kernel);
                if (dr == 0.0) continue;
                int idx = cs_index(cinfo->cs, i, j, k);
                double f_comm[3] = { dF_total[X]*dr,
                                     dF_total[Y]*dr,
                                     dF_total[Z]*dr };
                hydro_f_local_add(hydro, idx, f_comm);
              }
        }
      }

  hydro_memcpy(hydro, tdpMemcpyHostToDevice);

  /* NOTE: do NOT call colloid_sums_halo here. subgrid_update_forces_electro-
   * kinetics already synchronised pc->fex before this function ran; calling it
   * again would re-sum the full fex (not just our correction) and break the
   * action-reaction balance. The correction force on the particle and the
   * Peskin-distributed reaction on the fluid are computed locally and balance
   * by construction (serial / single-domain). */

  return 0;
}

/*CHANGE END - 20260625 pm_sr_correct_phi / pm_sr_apply_force_correction */