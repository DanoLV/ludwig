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
/*CHANGE END - Subgrid charge */

// static double d_peskin(double);
static int subgrid_interpolation(colloids_info_t* cinfo, hydro_t* hydro);
static const double drange_ = 1.0; /* Max. range of interpolation - 1 */

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
static void add_charge_to_array(distributed_charge_klein_t** charge_ptr, int cs_index,
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

								/* Separation between r0 and the coordinate position of
								 * this site */

								r[X] = r0[X] - 1.0 * i;
								r[Y] = r0[Y] - 1.0 * j;
								r[Z] = r0[Z] - 1.0 * k;

								dr = d_peskin(r[X]) * d_peskin(r[Y]) * d_peskin(r[Z]);

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
			entry->force[i] = malloc(sizeof(klein_t));
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
	double E_field[3] = { 0.0, 0.0, 0.0 };      /* force on particle from this lattice site */
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

					// Force goes to fex as this force is calculated ouside of ludwig_colloids_update
					force[X] = kt * reunit * pc->Esub[X] * (pc->s.q0 - pc->s.q1);
					force[Y] = kt * reunit * pc->Esub[Y] * (pc->s.q0 - pc->s.q1);
					force[Z] = kt * reunit * pc->Esub[Z] * (pc->s.q0 - pc->s.q1);

					/* Translate colloid position to local coordinates */
					r0[X] = pc->s.r[X] - 1.0 * offset[X];
					r0[Y] = pc->s.r[Y] - 1.0 * offset[Y];
					r0[Z] = pc->s.r[Z] - 1.0 * offset[Z];

					/* Work out which local lattice sites are involved */
					subgrid_get_lattice_index(r0, nlocal, &i_min, &i_max, &j_min, &j_max, &k_min, &k_max);

					for (i = i_min; i <= i_max; i++) {
						for (j = j_min; j <= j_max; j++) {
							for (k = k_min; k <= k_max; k++) {

								double force_aux[3] = { 0.0, 0.0, 0.0 };      /* force on particle from this lattice site */

								index = cs_index(cinfo->cs, i, j, k);

								/* Electric field at this lattice site */
								psi_electric_field(psi, index, e);

								/* Separation between r0 and the lattice site */
								r[X] = r0[X] - 1.0 * i;
								r[Y] = r0[Y] - 1.0 * j;
								r[Z] = r0[Z] - 1.0 * k;

								/* Peskin delta function weight */
								dr = d_peskin(r[X]) * d_peskin(r[Y]) * d_peskin(r[Z]);

								/* Force on particle from electric field at this site index*/
								force_aux[X] = force[X] * dr;
								force_aux[Y] = force[Y] * dr;
								force_aux[Z] = force[Z] * dr;

								/* Add to particle force */
								klein_add_double(&force_k[X], force_aux[X]);
								klein_add_double(&force_k[Y], force_aux[Y]);
								klein_add_double(&force_k[Z], force_aux[Z]);

								/* Add to Klein sum array with binary search */
								add_force_to_array(&force_k_indexed, index, force_aux);

								/* Accumulate contribution to total force on system */
								klein_add_double(&flocal_k[X], force_aux[X]);
								klein_add_double(&flocal_k[Y], force_aux[Y]);
								klein_add_double(&flocal_k[Z], force_aux[Z]);

							}
						}
					}
					// Sum goes to fex as this force is calculated ouside of ludwig_colloids_update
					pc->fex[X] += klein_sum(&force_k[X]);
					pc->fex[Y] += klein_sum(&force_k[Y]);
					pc->fex[Z] += klein_sum(&force_k[Z]);

				}
			}
		}
	}

	colloid_sums_halo(cinfo, COLLOID_SUM_FORCE_EXT_ONLY);

	/* Now apply all accumulated electric forces to hydro by index */
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

	//----------------------------------------------------------------
	// Momentum conservation correction
	//----------------------------------------------------------------

	/* Count fluid nodes */
	map_volume_local(map, MAP_FLUID, &nsfluid);

	flocal[3] = (double)nsfluid;
	flocal[X] = klein_sum(&flocal_k[X]);
	flocal[Y] = klein_sum(&flocal_k[Y]);
	flocal[Z] = klein_sum(&flocal_k[Z]);

	/* Sum across all MPI ranks */
	MPI_Allreduce(flocal, fsum, 4, MPI_DOUBLE, MPI_SUM, comm);

	/* Calculate average force per fluid node */
	if (fsum[3] > 0.0) {
		fsum[X] /= fsum[3];
		fsum[Y] /= fsum[3];
		fsum[Z] /= fsum[3];
	}

	/* Second pass: Apply correction force to all fluid nodes to conserve momentum */
	if (hydro && fsum[3] > 0.0) {
		for (ic = 1; ic <= nlocal[X]; ic++) {
			for (jc = 1; jc <= nlocal[Y]; jc++) {
				for (kc = 1; kc <= nlocal[Z]; kc++) {

					index = cs_index(cinfo->cs, ic, jc, kc);

					/* Check if this is a fluid node */
					colloids_info_map(cinfo, index, &pc);
					if (pc) continue;  // Skip colloid nodes. If there are solid nodes, this must be corrected.

					/* Apply negative of average force to conserve momentum */
					force[X] = -fsum[X];
					force[Y] = -fsum[Y];
					force[Z] = -fsum[Z];

					hydro_f_local_add(hydro, index, force);
				}
			}
		}
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
						FILE* fp) {

	int i, j, k, ic, jc, kc, i_min, i_max, j_min, j_max, k_min, k_max;
	int ncell[3];
	int index;
	int nlocal[3], offset[3];
	double dr;
	double r[3], r0[3];
	double e[3];           						/* Total electric field */
	double E_field[3] = { 0.0, 0.0, 0.0 };      /* Electric field on particle from this lattice site */
	double E_self[3] = { 0.0, 0.0, 0.0 };      	/* Self-field */

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

	/* Loop through all cells (including the halo cells) and set
	 * the field at each particle to zero for this step. */

	for (ic = 0; ic <= ncell[X] + 1; ic++) {
		for (jc = 0; jc <= ncell[Y] + 1; jc++) {
			for (kc = 0; kc <= ncell[Z] + 1; kc++) {

				colloids_info_cell_list_head(cinfo, ic, jc, kc, &pc);

				for (; pc; pc = pc->next) {

					if (pc->s.bc != COLLOID_BC_SUBGRID) continue;

					pc->Esub[X] = 0.0;
					pc->Esub[Y] = 0.0;
					pc->Esub[Z] = 0.0;
				}
			}
		}
	}

	/* First pass: Calculate electric field on particles*/
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

								/* Peskin delta function weight */
								dr = d_peskin(r[X]) * d_peskin(r[Y]) * d_peskin(r[Z]);

								/* Electric field at this lattice site */
								psi_electric_field(psi, index, e);

								/* Field on particle from electric field at this site */
								E_field[X] = e[X] * dr;
								E_field[Y] = e[Y] * dr;
								E_field[Z] = e[Z] * dr;

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

					// /* Write data for each particle */
					char string[256];
					double Emod = sqrt(pc->Esub[X]*pc->Esub[X] + pc->Esub[Y]*pc->Esub[Y] + pc->Esub[Z]*pc->Esub[Z]);
					sprintf(string, "%d;%d;%.15e;%.15e; %.15e;%.15e;%.15e; %.15e;%.15e\n", step, Emod, pc->s.index, pc->Esub[X], pc->Esub[Y], pc->Esub[Z], E_self[X], E_self[Y], E_self[Z]);
					for (i = 0;i < strlen(string);i++)if (string[i] == '.')string[i] = ',';
					fprintf(fp, "%s", string);

					// /* NUEVO: Calcular el autocampo de la partícula */
					// if (subgrid_compute_self_field_single_particle(pc, cinfo, psi, E_self) == 0) {

					// 	/* Write data for each particle */
					// 	char string[256];
					// 	double Emod = sqrt(pc->Esub[X]*pc->Esub[X] + pc->Esub[Y]*pc->Esub[Y] + pc->Esub[Z]*pc->Esub[Z]);
					// 	sprintf(string, "%d;%d;%.15e;%.15e; %.15e;%.15e;%.15e; %.15e;%.15e\n", step, Emod, pc->s.index, pc->Esub[X], pc->Esub[Y], pc->Esub[Z], E_self[X], E_self[Y], E_self[Z]);
					// 	for (i = 0;i < strlen(string);i++)if (string[i] == '.')string[i] = ',';
					// 	fprintf(fp, "%s", string);

					// 	/* Restar el autocampo del campo total */
					// 	pc->Esub[X] -= E_self[X];
					// 	pc->Esub[Y] -= E_self[Y];
					// 	pc->Esub[Z] -= E_self[Z];
					// }

				}
			}
		}
	}

	colloid_sums_halo(cinfo, COLLOID_SUM_ELECTRIC_FIELD);

	return 0;

}

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
	double dr_i, dr_j, q_i, q_j;
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
				// q_i = (pc->s.q0 - pc->s.q1) * dr_i;

				// prueba loca
				int index = cs_index(cinfo->cs, i, j, k);
				double rho0,rho1;
				psi_rho(psi_global, index, 0, &rho0);
				psi_rho(psi_global, index, 0, &rho1);
				q_i = ( (rho0 - rho1) + (pc->s.q0 - pc->s.q1) ) * dr_i;
				// prueba loca

				if (fabs(q_i) < 1e-14) continue; /* Skip if no charge */

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
	*i_min = imax(1, (int)floor(r0[X] - drange_));
	*i_max = imin(nlocal[X], (int)ceil(r0[X] + drange_));
	*j_min = imax(1, (int)floor(r0[Y] - drange_));
	*j_max = imin(nlocal[Y], (int)ceil(r0[Y] + drange_));
	*k_min = imax(1, (int)floor(r0[Z] - drange_));
	*k_max = imin(nlocal[Z], (int)ceil(r0[Z] + drange_));
}

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

	/* Loop through all cells (including the halo cells) and set
	 * the velocity at each particle to zero for this step. */

	for (ic = 0; ic <= ncell[X] + 1; ic++) {
		for (jc = 0; jc <= ncell[Y] + 1; jc++) {
			for (kc = 0; kc <= ncell[Z] + 1; kc++) {

				colloids_info_cell_list_head(cinfo, ic, jc, kc, &p_colloid);

				for (; p_colloid; p_colloid = p_colloid->next) {

					if (p_colloid->s.bc != COLLOID_BC_SUBGRID) continue;

					p_colloid->fsub[X] = 0.0;
					p_colloid->fsub[Y] = 0.0;
					p_colloid->fsub[Z] = 0.0;
				}
			}
		}
	}

	/* And add up the contributions to the velocity from the lattice. */

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
								hydro_u(hydro, index, u);

								p_colloid->fsub[X] += u[X] * dr;
								p_colloid->fsub[Y] += u[Y] * dr;
								p_colloid->fsub[Z] += u[Z] * dr;
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

	rmod = fabs(r);

	if (rmod <= 1.0) {
		delta = 0.125 * (3.0 - 2.0 * rmod + sqrt(1.0 + 4.0 * rmod - 4.0 * rmod * rmod));
	}
	else if (rmod <= 2.0) {
		delta = 0.125 * (5.0 - 2.0 * rmod - sqrt(-7.0 + 12.0 * rmod - 4.0 * rmod * rmod));
	}

	return delta;
}