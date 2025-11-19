/****************************************************************************
 *
 *  stats_velocity.h
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *  (c) 2011-2020 The University of Edinburgh
 *
 ****************************************************************************/

#ifndef STATS_VELOCITY_H
#define STATS_VELOCITY_H

#include "hydro.h"
#include "map.h"
/*CHANGE INIT - Total force calculation */
#include "colloids.h"
/*CHANGE END - Total force calculation */

typedef struct stats_vel_s stats_vel_t;

struct stats_vel_s {
  int print_vol_flux;
};

stats_vel_t stats_vel_default(void);

int stats_velocity_minmax(stats_vel_t * stat, hydro_t * hydro, map_t * map);

/*CHANGE INIT - Total force calculation */
int stats_total_force(hydro_t * hydro, map_t * map, colloids_info_t * cinfo,
                      double ffluid[3], double fcoll[3], double fsubgrid[3],
                      double ftotal[3]);
int stats_total_force_write(hydro_t * hydro, map_t * map, colloids_info_t * cinfo,
                             int timestep, const char * filename);
/*CHANGE END - Total force calculation */

#endif
