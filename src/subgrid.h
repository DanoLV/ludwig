/*****************************************************************************
 *
 *  subgrid.h
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *  (c) 2010-2021 The University of Edinburgh
 *
 *****************************************************************************/

#ifndef LUDWIG_SUBGRID_H
#define LUDWIG_SUBGRID_H

#include "colloids.h"
#include "hydro.h"
#include "wall.h"
 /*CHANGE INIT - Subgrid charge */
#include "psi_colloid.h"
#include "util_sum.h"

typedef struct distributed_charge_klein_entry_s {
  int cs_index;                /* Lattice site index */
  double rho0_original;        /* Original value of rho0 */
  double rho1_original;        /* Original value of rho1 */
  klein_t* rho0_sum;           /* Pointer to Klein sum for rho0 */
  klein_t* rho1_sum;           /* Pointer to Klein sum for rho1 */
} distributed_charge_klein_entry_t;

typedef struct distributed_charge_klein_s {
  distributed_charge_klein_entry_t** entries;  /* Array of pointers to entries */
  int count;                   /* Number of entries */
  int capacity;                /* Allocated capacity */
} distributed_charge_klein_t;

typedef struct distributed_force_klein_entry_s {
  int cs_index;                /* Lattice site index */
  klein_t** force;           /* Pointer to Klein sum for force */
} distributed_force_klein_entry_t;

typedef struct distributed_force_klein_s {
  distributed_force_klein_entry_t** entries;  /* Array of pointers to entries */
  int count;                   /* Number of entries */
  int capacity;                /* Allocated capacity */
} distributed_force_klein_t;

double d_peskin(double r);
int subgrid_charge_from_particles(colloids_info_t* cinfo, psi_t* obj, distributed_charge_klein_t** charge);
int subgrid_charge_from_particles_substract(colloids_info_t* cinfo, psi_t* obj, distributed_charge_klein_t** charge);
int subgrid_charge_from_particles_restore(colloids_info_t* cinfo, psi_t* obj, distributed_charge_klein_t** charge);
void subgrid_free_distributed_charge_t(distributed_charge_klein_t** charge);
void subgrid_get_lattice_index(double r0[3], int nlocal[3], int* i_min, int* i_max, int* j_min, int* j_max, int* k_min, int* k_max);
int subgrid_update_forces_electrokinetics(colloids_info_t* cinfo, map_t* map, physics_t* phys, psi_t* psi, hydro_t* hydro);
int subgrid_update_Esub(colloids_info_t* cinfo, psi_t* psi, int step, FILE* fp);
int subgrid_compute_self_field_single_particle(colloid_t* pc, 
                                                colloids_info_t* cinfo,
                                                psi_t* psi_global,
                                                double E_self[3]);
/*CHANGE END - Subgrid charge */

int subgrid_update(colloids_info_t* cinfo, hydro_t* hydro, int noise_flag);
int subgrid_force_from_particles(colloids_info_t* cinfo, hydro_t* hydro, wall_t* wall);
int subgrid_wall_lubrication(colloids_info_t* cinfo, wall_t* wall);

#endif
