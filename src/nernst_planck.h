/*****************************************************************************
 *
 *  nernst_planck.h
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  (c) 2012-2022 The University of Edinburgh
 *
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *
 *****************************************************************************/

#ifndef LUDWIG_NERNST_PLANCK_H
#define LUDWIG_NERNST_PLANCK_H

#include "psi.h"
#include "free_energy.h"
#include "hydro.h"
#include "map.h"
#include "colloids.h"

int nernst_planck_driver(psi_t * psi, fe_t * fe, map_t * map);
int nernst_planck_driver_d3qx(psi_t * psi, fe_t * fe, hydro_t * hydro,
			      map_t * map, colloids_info_t * cinfo);
/*CHANGE INIT - 20260326 GPU driver for Nernst-Planck d3qx */
int nernst_planck_driver_d3qx_gpu(psi_t * psi, fe_t * fe, hydro_t * hydro,
                                   map_t * map, colloids_info_t * cinfo);
/*CHANGE END - 20260326 */
int nernst_planck_adjust_multistep(psi_t * psi);

/*CHANGE INIT - 20261112 Lattice (FLLA) advection of the ions.
 * DeVincentis & Thomas, Comput. Math. Appl. 144, 1 (2023): advective fluxes
 * from the LB distribution (their eq. 21) + linear reconstruction with the
 * superbee flux limiter (their eqs. 13-19). Selected with the input key
 *   electrokinetics_advection  centred | flla | flla_mass   (default centred)
 * GPU driver only. */
#include "lb_data.h"
int nernst_planck_advection_scheme_set(pe_t * pe, const char * name);
int nernst_planck_lattice_flux_compute(lb_t * lb, psi_t * psi);
/*CHANGE END - 20261112 */

int nernst_planck_maxacc(double * acc);
int nernst_planck_maxacc_set(double acc);

#endif
