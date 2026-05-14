/*****************************************************************************
 *
 *  psi_force.h
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  (c) 2013-2023 The University of Edinburgh
 *
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *    
 *****************************************************************************/

#ifndef LUDWIG_PSI_FORCE_H
#define LUDWIG_PSI_FORCE_H

#include "psi.h"
#include "free_energy.h"
#include "hydro.h"
#include "colloids.h"
#include "map.h"
/*CHANGE INIT - 20260422 kernel parameter in psi_force_gradmu */
#include "subgrid.h"

int psi_force_gradmu(psi_t * psi,  fe_t * fe, field_t * phi, hydro_t * hydro,
		     map_t * map, colloids_info_t * cinfo, subgrid_kernel_t kernel);
int psi_force_gradmu_offset(psi_t * psi,  fe_t * fe, field_t * phi, hydro_t * hydro,
		     map_t * map, colloids_info_t * cinfo, subgrid_kernel_t kernel,
		     double mesh_offset);
/*CHANGE END - 20260422 kernel parameter in psi_force_gradmu */
int psi_force_divstress(psi_t * psi, fe_t * fe, hydro_t * hydro, 
			colloids_info_t * cinfo);

#endif
