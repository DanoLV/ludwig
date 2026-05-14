/*****************************************************************************
 *
 *  ewald_charge.h
 *
 *  Ewald summation for charged particles and lattice charge density.
 *  Adapted from ewald.h (magnetic dipoles) for point charges.
 *
 *  This implementation handles BOTH:
 *    1. Point charges on colloids (discrete particles)
 *    2. Charge density on the LB lattice nodes (continuous field)
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  Contributing authors:
 *  Based on ewald.c by Grace Kim and Kevin Stratford
 *
 *****************************************************************************/

#ifndef LUDWIG_EWALD_CHARGE_H
#define LUDWIG_EWALD_CHARGE_H

typedef struct ewald_charge_s ewald_charge_t;

#include "pe.h"
#include "coords.h"
#include "colloids.h"
#include "psi.h"
#include "map.h"
#include "hydro.h"

/* Source flags: what charge sources to include in the Ewald sum */
typedef enum {
  EWALD_SOURCE_COLLOIDS = 1,    /* Include colloid point charges */
  EWALD_SOURCE_LATTICE  = 2,    /* Include lattice charge density */
  EWALD_SOURCE_ALL      = 3     /* Both colloids and lattice */
} ewald_charge_source_t;

int ewald_charge_create(pe_t * pe, cs_t * cs, double rc,
                        double alpha, double epsilon_prime,
                        colloids_info_t * cinfo, psi_t * psi,
                        map_t * map, hydro_t * hydro, ewald_charge_source_t sources,
                        ewald_charge_t ** pewald);
int ewald_charge_set_eps_reg(double eps_reg);
int ewald_charge_free(ewald_charge_t * ewald);
int ewald_charge_info(ewald_charge_t * ewald);
int ewald_charge_check_electroneutrality(ewald_charge_t * ewald);
int ewald_charge_kappa(ewald_charge_t * ewald, double * kappa);
int ewald_charge_sum(ewald_charge_t * ewald);
int ewald_charge_sum_full(ewald_charge_t * ewald, FILE * fp);
int ewald_charge_real_space_sum(ewald_charge_t * ewald);
int ewald_charge_fourier_space_sum(ewald_charge_t * ewald);
int ewald_charge_external_field(ewald_charge_t * ewald);

int ewald_charge_total_energy(ewald_charge_t * ewald, double * ereal,
                              double * efourier, double * eself);
int ewald_charge_fourier_space_energy(ewald_charge_t * ewald, double * ef);
int ewald_charge_self_energy(ewald_charge_t * ewald, double * es);
int ewald_charge_real_space_energy(ewald_charge_t * ewald, double q1, double q2,
                                   const double r12[3], double * er);

/* Get potential at a point (for field calculation) */
int ewald_charge_potential_at_colloid(ewald_charge_t * ewald, colloid_t * pc,
                                      double * phi);

/* Get potential at a lattice node */
int ewald_charge_potential_at_node(ewald_charge_t * ewald, int ic, int jc, int kc,
                                   double * phi);

/* Compute and store potential/field on entire lattice */
int ewald_charge_compute_psi(ewald_charge_t * ewald);

/*CHANGE INIT - 20260309 Self-field correction for lattice-discretisation artefact */
/* Lookup-table dimensions for the self-field correction */
#define EWALD_SELF_NX 16
#define EWALD_SELF_NY 16
#define EWALD_SELF_NZ 16

/*
 * Opaque handle for the self-field lookup table.
 * Allocated by ewald_charge_sum_full_gpu_self_build(), freed by
 * ewald_charge_self_table_free().
 */
typedef struct ewald_self_table_s ewald_self_table_t;

/* GPU-accelerated version of ewald_charge_sum_full */
/* self_table may be NULL (no correction applied) */
int ewald_charge_sum_full_gpu(ewald_charge_t * ewald, FILE * fp,
                              ewald_self_table_t * self_table);

/* Fourier-space only version: no real-space or dipole contributions */
int ewald_charge_sum_full_gpu_fourier(ewald_charge_t * ewald, FILE * fp,
                                      ewald_self_table_t * self_table);

/* Build the self-field lookup table (call once after ewald_charge_create). */
int ewald_charge_sum_full_gpu_self_build(ewald_charge_t * ewald,
                                         double kappa_debye,
                                         int    nx, int ny, int nz,
                                         int    L,
                                         ewald_self_table_t ** ptable);

/* Symmetric version: requires n even.  Computes only the irreducible simplex
 * ix >= iy >= iz in [0, n/2] and expands the rest by Oh cubic symmetry
 * (~1/48 of the compute work for large n). */
int ewald_charge_sum_full_gpu_self_build_symmetric(ewald_charge_t * ewald,
                                                   double kappa_debye,
                                                   int    n,
                                                   int    L,
                                                   ewald_self_table_t ** ptable);

/* Free the lookup table. */
int ewald_charge_self_table_free(ewald_self_table_t ** ptable);

/* Evaluate the self-field at fractional position (xf,yf,zf) in [0,1)^3
 * using trilinear interpolation.  Writes result to Eself[3]. */
int ewald_charge_self_field_interpolate(const ewald_self_table_t * table,
                                        double xf, double yf, double zf,
                                        double Eself[3]);
/*CHANGE END - 20260309*/

/* CHANGE INIT - PoissonVerification - Compile-time parameters for face sampling */
/* Number of quadrature points per face per direction (total per face: N×N) */
#define EWALD_NSAMPLE_FACE     32
/* Stride for node sampling: 1=all nodes, N>1=sample 1 in N nodes (debug) */
#define EWALD_POISSON_STRIDE   1
/* Window mode: sample only ±EWALD_POISSON_WINDOW nodes around domain centre.
 * Set to 0 to disable (use full domain with STRIDE). */
#define EWALD_POISSON_WINDOW   0

/* Helper macro: compute window bounds [i0,i1] given domain size N and centre c=N/2.
 * Used consistently in face_sample, calc_div, write_poisson, write_force. */
#if EWALD_POISSON_WINDOW > 0
#define EWALD_WINDOW_BOUNDS(N, i0, i1) \
  do { int _c = (N)/2, _w = EWALD_POISSON_WINDOW; \
       (i0) = (_c - _w) < 0       ? 0     : (_c - _w); \
       (i1) = (_c + _w) >= (N)    ? (N)-1 : (_c + _w); } while(0)
#else
#define EWALD_WINDOW_BOUNDS(N, i0, i1)  do { (i0) = 0; (i1) = (N)-1; } while(0)
#endif
/* Distance threshold: particle closer than this to a cube face → flag for diagnostic */
#define EWALD_FACE_SAFE_RADIUS 0.5
/* Node flag values */
#define EWALD_NODE_NORMAL        0
#define EWALD_NODE_NEAR_PARTICLE 1
/* CHANGE END - PoissonVerification */
#define EWALD_POISSON_ONLY 
/* Compile with -DEWALD_POISSON_ONLY to skip Maxwell stress (faster for Poisson check) */
/* CHANGE END - PoissonVerification */

/* CHANGE INIT - Gaussian_Ewald - Gaussian Ewald sum: particles and nodes as Gaussian charge distributions */
/* sigma: Gaussian width (physical); splitting parameter is alpha from ewald_charge_create */
int ewald_charge_sum_full_gaussian_gpu(ewald_charge_t * ewald, FILE * fp,
                                       double sigma);
/* CHANGE END - Gaussian_Ewald */

/*CHANGE INIT - 20260212 New ewald_charge_sum_FFT_full_gpu using cuFFT for Fourier-space*/
/* GPU-accelerated FFT version: hybrid FFT lattice + direct particle loops */
int ewald_charge_sum_FFT_full_gpu(ewald_charge_t * ewald, FILE * fp);
/* GPU-accelerated FFT version: all charges (lattice + particles) spread to
 * grid with Peskin kernel before FFT. No continuous-position structure factors. */
int ewald_charge_sum_FFT_full_gpu_peskin(ewald_charge_t * ewald, FILE * fp);
/*CHANGE END - 20260212*/

/* CHANGE INIT - PoissonVerification - New function: Poisson verification + Maxwell stress force */
/* Assumes ewald_charge_sum_full_gpu() was already called this timestep (sinx_/cosx_ valid).
 * step: simulation timestep, used for output filename. */
int ewald_charge_sum_full_gpu_poisson_force(ewald_charge_t * ewald, FILE * fp,
                                            ewald_self_table_t * self_table,
                                            int step);
/* CHANGE END - PoissonVerification */

/* Set/get dielectric constant of surrounding medium (Deserno & Holm convention)
 * epsilon' = 1     : vacuum boundary conditions (default)
 * epsilon' = infty : metallic (tinfoil) boundary conditions */
int ewald_charge_set_epsilon_prime(ewald_charge_t * ewald, double epsilon_prime);
/* CHANGE INIT - 20260513 */
int ewald_charge_set_cinfo(ewald_charge_t * ewald, colloids_info_t * cinfo);
int ewald_charge_set_map  (ewald_charge_t * ewald, map_t           * map);
int ewald_charge_set_hydro(ewald_charge_t * ewald, hydro_t         * hydro);
/* CHANGE END - 20260513 */
int ewald_charge_get_epsilon_prime(ewald_charge_t * ewald, double * epsilon_prime);

#endif
