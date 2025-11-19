/****************************************************************************
 *
 *  stats_velocity.c
 *
 *  Basic statistics for the velocity field.
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  (c) 2011-2020 The University of Edinburgh
 *
 *  Contributing authors:
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *
 ****************************************************************************/

#include <assert.h>
#include <float.h>

#include "pe.h"
#include "coords.h"
#include "hydro.h"
#include "stats_velocity.h"
#include "util.h"
/*CHANGE INIT - Total force calculation */
#include "util_sum.h"
#include "colloids.h"
/*CHANGE END - Total force calculation */

/****************************************************************************
 *
 *  stats_vel_default
 *
 *  Always print minmax at moment.
 *
 ****************************************************************************/

stats_vel_t stats_vel_default(void) {

  stats_vel_t stat = {0};

  return stat;
}

/****************************************************************************
 *
 *  stats_velocity_minmax
 *
 *  The volume flux of is of interest for porous media calculations
 *  of permeability. Note that with the body force density f, the
 *  volume flux is the same as the total momemtum  plus 0.5*f per
 *  lattice site. So for complex porous media, the total momentum
 *  can actually look quite wrong (e.g., have the opposite sign to
 *  the flow).
 *
 ****************************************************************************/

int stats_velocity_minmax(stats_vel_t * stat, hydro_t * hydro, map_t * map) {

  int ic, jc, kc, ia, index;
  int nlocal[3];
  int status;

  double umin[3];
  double umax[3];
  double utmp[3];
  double usum_local[3], usum[3];

  MPI_Comm comm;

  assert(stat);
  assert(hydro);
  assert(map);

  cs_nlocal(hydro->cs, nlocal);
  pe_mpi_comm(hydro->pe, &comm);

  for (ia = 0; ia < 3; ia++) {
    umin[ia] = FLT_MAX;
    umax[ia] = FLT_MIN;
    usum_local[ia] = 0.0;
  }

  for (ic = 1; ic <= nlocal[X]; ic++) {
    for (jc = 1; jc <= nlocal[Y]; jc++) {
      for (kc = 1; kc <= nlocal[Z]; kc++) {

        index = cs_index(hydro->cs, ic, jc, kc);
	map_status(map, index, &status);

	if (status == MAP_FLUID) {

	  hydro_u(hydro, index, utmp);

	  for (ia = 0; ia < 3; ia++) {
	    umin[ia] = dmin(umin[ia], utmp[ia]);
	    umax[ia] = dmax(umax[ia], utmp[ia]);
	    usum_local[ia] += utmp[ia];
	  }
	}
      }
    }
  }

  MPI_Reduce(umin, utmp, 3, MPI_DOUBLE, MPI_MIN, 0, comm);

  for (ia = 0; ia < 3; ia++) {
    umin[ia] = utmp[ia];
  }

  MPI_Reduce(umax, utmp, 3, MPI_DOUBLE, MPI_MAX, 0, comm);

  for (ia = 0; ia < 3; ia++) {
    umax[ia] = utmp[ia];
  }

  MPI_Reduce(usum_local, usum, 3, MPI_DOUBLE, MPI_SUM, 0, comm);

  pe_info(hydro->pe, "\n");
  pe_info(hydro->pe, "Velocity - x y z\n");
  pe_info(hydro->pe, "[minimum ] %14.7e %14.7e %14.7e\n", umin[X], umin[Y], umin[Z]);
  pe_info(hydro->pe, "[maximum ] %14.7e %14.7e %14.7e\n", umax[X], umax[Y], umax[Z]);

  if (stat->print_vol_flux) {
    pe_info(hydro->pe, "[vol flux] %14.7e %14.7e %14.7e\n", usum[X], usum[Y], usum[Z]);
  }

  return 0;
}

/*CHANGE INIT - Total force calculation */
/****************************************************************************
 *
 *  stats_total_force
 *
 *  Calculate total force on the system using Klein summation:
 *  - Force on fluid nodes (from hydro->force)
 *  - Force on resolved colloids (force + fex)
 *  - Force on subgrid colloids (fex)
 *
 *  Returns individual components in ffluid, fcoll, fsubgrid, and ftotal.
 *
 ****************************************************************************/

int stats_total_force(hydro_t * hydro, map_t * map, colloids_info_t * cinfo,
                      double ffluid[3], double fcoll[3], double fsubgrid[3],
                      double ftotal[3]) {

  int ic, jc, kc, ia, index;
  int nlocal[3];
  int status;

  klein_t ffluid_local[3];
  klein_t fcoll_local[3];
  klein_t fsubgrid_local[3];
  double ftmp[3];

  colloid_t * pc = NULL;
  MPI_Comm comm;
  MPI_Datatype klein_mpi_type;
  MPI_Op klein_mpi_sum;

  assert(hydro);
  assert(map);
  assert(cinfo);
  assert(ffluid);
  assert(fcoll);
  assert(fsubgrid);
  assert(ftotal);

  cs_nlocal(hydro->cs, nlocal);
  pe_mpi_comm(hydro->pe, &comm);

  /* Initialize Klein sums */
  for (ia = 0; ia < 3; ia++) {
    ffluid_local[ia] = klein_zero();
    fcoll_local[ia] = klein_zero();
    fsubgrid_local[ia] = klein_zero();
  }

  /* Copy force data from device to host */
  hydro_memcpy(hydro, tdpMemcpyDeviceToHost);

  /* Sum forces on fluid nodes using Klein summation */
  for (ic = 1; ic <= nlocal[X]; ic++) {
    for (jc = 1; jc <= nlocal[Y]; jc++) {
      for (kc = 1; kc <= nlocal[Z]; kc++) {

        index = cs_index(hydro->cs, ic, jc, kc);
        map_status(map, index, &status);

        if (status == MAP_FLUID) {
          hydro_f_local(hydro, index, ftmp);
          for (ia = 0; ia < 3; ia++) {
            klein_add_double(&ffluid_local[ia], ftmp[ia]);
          }
        }
      }
    }
  }

  /* Sum forces on colloids using Klein summation */
  colloids_info_local_head(cinfo, &pc);

  for (; pc; pc = pc->nextlocal) {
    if (pc->s.bc == COLLOID_BC_SUBGRID) {
      /* Subgrid particles: only external force (fex) */
      for (ia = 0; ia < 3; ia++) {
        klein_add_double(&fsubgrid_local[ia], pc->fex[ia]);
      }
    }
    else {
      /* Resolved colloids: force + fex */
      for (ia = 0; ia < 3; ia++) {
        klein_add_double(&fcoll_local[ia], pc->force[ia]);
        klein_add_double(&fcoll_local[ia], pc->fex[ia]);
      }
    }
  }

  /* Set up MPI operations for Klein type */
  klein_mpi_datatype(&klein_mpi_type);
  klein_mpi_op_sum(&klein_mpi_sum);

  /* Reduce across all MPI ranks using Klein MPI sum */
  klein_t ffluid_klein[3], fcoll_klein[3], fsubgrid_klein[3];

  MPI_Reduce(ffluid_local, ffluid_klein, 3, klein_mpi_type, klein_mpi_sum, 0, comm);
  MPI_Reduce(fcoll_local, fcoll_klein, 3, klein_mpi_type, klein_mpi_sum, 0, comm);
  MPI_Reduce(fsubgrid_local, fsubgrid_klein, 3, klein_mpi_type, klein_mpi_sum, 0, comm);

  /* Extract sums */
  for (ia = 0; ia < 3; ia++) {
    ffluid[ia] = klein_sum(&ffluid_klein[ia]);
    fcoll[ia] = klein_sum(&fcoll_klein[ia]);
    fsubgrid[ia] = klein_sum(&fsubgrid_klein[ia]);
    // ftotal[ia] = ffluid[ia] + fcoll[ia] + fsubgrid[ia];
  }

  /* Free MPI operations */
  MPI_Type_free(&klein_mpi_type);
  MPI_Op_free(&klein_mpi_sum);

  return 0;
}

/****************************************************************************
 *
 *  stats_total_force_write
 *
 *  Write total force to file
 *
 ****************************************************************************/

int stats_total_force_write(hydro_t * hydro, map_t * map, colloids_info_t * cinfo,
                             int timestep, const char * filename) {

  double ffluid[3], fcoll[3], fsubgrid[3], ftotal[3];
  FILE * fp = NULL;
  int rank;
  MPI_Comm comm;

  assert(hydro);
  assert(map);
  assert(cinfo);
  assert(filename);

  pe_mpi_comm(hydro->pe, &comm);
  MPI_Comm_rank(comm, &rank);

  /* Calculate forces */
  stats_total_force(hydro, map, cinfo, ffluid, fcoll, fsubgrid, ftotal);

  /* Only rank 0 writes to file */
  if (rank == 0) {
    fp = fopen(filename, "a");
    if (fp == NULL) {
      pe_fatal(hydro->pe, "Failed to open force output file: %s\n", filename);
      return -1;
    }

    /* Write header if file is empty (first write) */
    fseek(fp, 0, SEEK_END);
    if (ftell(fp) == 0) {
      fprintf(fp, "# Total force output\n");
      fprintf(fp, "timestep;");// ftotal_x ftotal_y ftotal_z ");
      fprintf(fp, "ffluid_x;ffluid_y;ffluid_z;");
      // fprintf(fp, "fcoll_x;fcoll_y;fcoll_z;");
      fprintf(fp, "fsubgrid_x;fsubgrid_y;fsubgrid_z\n");
    }

    /* Write force data */
    // fprintf(fp, "%d;%.15e;%.15e;%.15e;", timestep, ftotal[X], ftotal[Y], ftotal[Z]);
    fprintf(fp, "%d;", timestep);
    fprintf(fp, "%.15e;%.15e;%.15e;", ffluid[X], ffluid[Y], ffluid[Z]);
    // fprintf(fp, "%.15e;%.15e;%.15e;", fcoll[X], fcoll[Y], fcoll[Z]);
    fprintf(fp, "%.15e;%.15e;%.15e\n", fsubgrid[X], fsubgrid[Y], fsubgrid[Z]);

    fclose(fp);
  }

  return 0;
}
/*CHANGE END - Total force calculation */
