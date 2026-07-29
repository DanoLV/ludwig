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
double d_peskin_derivative(double r);
double d_trilinear(double r);
/*CHANGE INIT - 20260422 B-spline order-4 kernel and kernel selection type */
double d_bspline4(double r);
double d_bspline6(double r);
double d_kb4(double r);
void subgrid_set_kb4_beta(double beta);
double subgrid_kb4_norm_fluid(void);
double d_peskin6(double r);
/*CHANGE INIT - 20260630 Hann spread/gather kernel */
double d_hann(double r);
void subgrid_set_hann_order(double n);
/*CHANGE END - 20260630 */

typedef enum {
  SUBGRID_KERNEL_PESKIN4 = 0,
  SUBGRID_KERNEL_BSPLINE4 = 1,
  SUBGRID_KERNEL_BSPLINE6 = 2,
  SUBGRID_KERNEL_KB4 = 3,
  SUBGRID_KERNEL_PESKIN6 = 4,
  /*CHANGE INIT - 20260630 Hann kernel id */
  SUBGRID_KERNEL_HANN = 5,
  /*CHANGE END - 20260630 */
  /*CHANGE INIT - 20260710 interpolating kernel for shifted-grid fluid ops.
   * Trilinear is a delta at integer offsets, so fluid scatter/gather with an
   * integer mesh shift reduces to an exact node relabelling (required for the
   * integer-shift identity test of the interlacing machinery). */
  SUBGRID_KERNEL_TRILINEAR = 6
  /*CHANGE END - 20260710 */
} subgrid_kernel_t;

int subgrid_charge_from_grid(colloids_info_t* cinfo, psi_t* obj, distributed_charge_klein_t** charge,
                              subgrid_kernel_t kernel);
int subgrid_charge_from_particles(colloids_info_t* cinfo, psi_t* obj, distributed_charge_klein_t** charge,
                                   subgrid_kernel_t kernel);
/*CHANGE END - 20260422 B-spline order-4 kernel and kernel selection type */
int subgrid_charge_from_particles_compenzate(colloids_info_t* cinfo, psi_t* obj, distributed_charge_klein_t** charge);
int subgrid_charge_from_particles_substract(colloids_info_t* cinfo, psi_t* obj, distributed_charge_klein_t** charge);
int subgrid_charge_from_particles_restore(colloids_info_t* cinfo, psi_t* obj, distributed_charge_klein_t** charge);
void subgrid_free_distributed_charge_t(distributed_charge_klein_t** charge);
void add_charge_to_array(distributed_charge_klein_t** charge_ptr, int cs_index,
                         double q0_dr, double q1_dr, psi_t* obj);
int subgrid_peskin_scatter_rho_buf(colloids_info_t* cinfo, psi_t* psi_src, double* rho_smooth, int ndata);
/*CHANGE INIT - 20260427 fluid scatter with offset into auxiliary buffer */
int subgrid_scatter_fluid_offset(colloids_info_t* cinfo, psi_t* obj,
                                  subgrid_kernel_t kernel, double mesh_offset,
                                  double* rho_buf, int ndata);
/*CHANGE END - 20260427 */
int subgrid_peskin_scatter_rho_gpu(colloids_info_t* cinfo, psi_t* psi_src, psi_t* psi_dst);
int subgrid_get_range(subgrid_kernel_t kernel);
void subgrid_get_lattice_index_range(double r0[3], int range, int nlocal[3], int* i_min, int* i_max, int* j_min, int* j_max, int* k_min, int* k_max);
void subgrid_get_lattice_index_range_halo(double r0[3], int range, int nlocal[3], int* i_min, int* i_max, int* j_min, int* j_max, int* k_min, int* k_max);
void subgrid_get_lattice_index(double r0[3], int nlocal[3], int* i_min, int* i_max, int* j_min, int* j_max, int* k_min, int* k_max);
void subgrid_get_lattice_index_fn(double r0[3], int nlocal[3], int* i_min, int* i_max, int* j_min, int* j_max, int* k_min, int* k_max);
/*CHANGE INIT - 20260422 kernel parameter for subgrid_update_forces_electrokinetics */
int subgrid_update_forces_electrokinetics(colloids_info_t* cinfo, map_t* map, physics_t* phys, psi_t* psi, hydro_t* hydro, subgrid_kernel_t kernel);
int subgrid_update_forces_electrokinetics_ewald(colloids_info_t* cinfo, map_t* map, physics_t* phys, psi_t* psi, hydro_t* hydro,subgrid_kernel_t kernel);
/*CHANGE END - 20260422 kernel parameter for subgrid_update_forces_electrokinetics */
int subgrid_update_forces_electrokinetics_theory(colloids_info_t* cinfo, map_t* map, physics_t* phys, psi_t* psi, hydro_t* hydro);
/*CHANGE INIT - 20260422 kernel parameter for subgrid_update_Esub */
int subgrid_print_Esub(colloids_info_t* cinfo,int step,FILE* fp);
int subgrid_update_Esub(colloids_info_t* cinfo, psi_t* psi, int step, FILE* fp, pe_t* pe,
                         subgrid_kernel_t kernel);
/*CHANGE END - 20260422 kernel parameter for subgrid_update_Esub */
/*CHANGE INIT - Poisson-Boltzmann force calculation */
int subgrid_force_poisson_boltzmann(colloids_info_t* cinfo, map_t* map,
                                     physics_t* phys, psi_t* psi, hydro_t* hydro,
                                     int step, FILE* fp);
/*CHANGE END - Poisson-Boltzmann force calculation */
int subgrid_compute_self_field_single_particle(colloid_t* pc,
                                                colloids_info_t* cinfo,
                                                psi_t* psi_global,
                                                double E_self[3]);
double d_idw(double r0[3], int node_i, int node_j, int node_k,
                    int i_min, int i_max, int j_min, int j_max,
                    int k_min, int k_max);
double d_isdw(double r0[3], int node_i, int node_j, int node_k,
                    int i_min, int i_max, int j_min, int j_max,
                    int k_min, int k_max);
/*CHANGE INIT - Poisson-Boltzmann weight */
double d_pb(double r0[3], int node_i, int node_j, int node_k,
            int i_min, int i_max, int j_min, int j_max,
            int k_min, int k_max, double kappa);
int subgrid_compute_kappa(psi_t* psi, double* kappa);
double subgrid_get_kappa(void);
void subgrid_set_kappa(double kappa);
/*CHANGE INIT - 20260422 subgrid_get_drange */
double subgrid_get_drange(void);
/*CHANGE END - 20260422 subgrid_get_drange */
int subgrid_add_psi_grad_from_particles(colloids_info_t* cinfo, psi_t* psi,
                                         physics_t* phys);
/*CHANGE END - Poisson-Boltzmann weight */

/*CHANGE INIT - 20260117 Short-range corrections declarations
 * See docs/SHORT_RANGE_CORRECTIONS_ANALYSIS.md for full documentation.
 * CHANGE 20260118: level1 disabled (no momentum conservation)
 *                  level3_pn implements PN correction with Newton III
 */
void subgrid_shortrange_set_level(int level);
int  subgrid_shortrange_get_level(void);
int  subgrid_shortrange_level1(colloids_info_t* cinfo, psi_t* psi, hydro_t* hydro);
int  subgrid_shortrange_level3_pn(colloids_info_t* cinfo, psi_t* psi, hydro_t* hydro);
int  subgrid_shortrange_corrections(colloids_info_t* cinfo, psi_t* psi,
                                     hydro_t* hydro);
/*CHANGE END - 20260117 Short-range corrections declarations */

/*CHANGE INIT - 20260625 PM short-range correction table */
typedef struct pm_sr_table_s pm_sr_table_t;

int  pm_sr_table_build_radial(double sigma, double r_cut, int n_voxel,
                               double epsilon, subgrid_kernel_t kernel,
                               pm_sr_table_t** ptable);
int  pm_sr_table_build_3d(double sigma, double r_cut, int nr, int n_voxel,
                           double epsilon, subgrid_kernel_t kernel,
                           pm_sr_table_t** ptable);
void pm_sr_table_free(pm_sr_table_t** ptable);
int  pm_sr_table_interpolate(const pm_sr_table_t* t, double r,
                              double* delta_phi, double* delta_E,
                              double* delta_F);
int  pm_sr_table_interpolate_3d(const pm_sr_table_t* t,
                                 double dx, double dy, double dz,
                                 double* delta_phi, double* delta_Ex,
                                 double* delta_Ey, double* delta_Ez);
void pm_sr_table_print(const pm_sr_table_t* t, FILE* fp);

int  pm_sr_correct_phi(colloids_info_t* cinfo, psi_t* psi,
                        const pm_sr_table_t* table, double* phi_pm_buf);
int  pm_sr_apply_force_correction(colloids_info_t* cinfo, map_t* map,
                                   psi_t* psi, hydro_t* hydro,
                                   const pm_sr_table_t* table,
                                   subgrid_kernel_t kernel);
/*CHANGE END - 20260625 PM short-range correction table */

/*CHANGE END - Subgrid charge */

int subgrid_update(colloids_info_t* cinfo, hydro_t* hydro, int noise_flag);
int subgrid_force_from_particles(colloids_info_t* cinfo, hydro_t* hydro, wall_t* wall);
int subgrid_wall_lubrication(colloids_info_t* cinfo, wall_t* wall);

#endif
