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
/*CHANGE END - Subgrid charge */

int subgrid_update(colloids_info_t * cinfo, hydro_t * hydro, int noise_flag);
int subgrid_force_from_particles(colloids_info_t * cinfo, hydro_t * hydro,
				 wall_t * wall);
int subgrid_wall_lubrication(colloids_info_t * cinfo, wall_t * wall);
/*CHANGE INIT - Subgrid charge */
double d_peskin(double r);
int subgrid_charge_from_particles(colloids_info_t* cinfo, psi_t* obj);
int subgrid_charge_from_particles_substract(colloids_info_t* cinfo, psi_t* obj);
void subgrid_get_lattice_index(double r0[3], int nlocal[3], int* i_min, int* i_max, int* j_min, int* j_max, int* k_min, int* k_max);
/*CHANGE END - Subgrid charge */
#endif
