/*****************************************************************************
 *
 *  psi_force.c
 *
 *  Compute the force on the fluid originating with charge.
 *
 *  Edinburgh Soft Matter and Statisitical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  (c) 2013-2023 The University of Edinburgh
 *
 *  Contributing authors:
 *    Kevin Stratford (kevin@epcc.ed.ac.uk)
 *    Ignacio Pagonabarraga
 *    Oliver Henrich
 *
 *****************************************************************************/

#include <assert.h>
#include <math.h>

#include "pe.h"
#include "coords.h"
#include "util.h"
#include "physics.h"
#include "fe_electro.h"
#include "fe_electro_symmetric.h"
#include "psi_force.h"
#include "psi_gradients.h"
#include "subgrid.h" //CHANGE

int psi_force_gradmu_e(psi_t* psi, fe_t* fe, hydro_t* hydro,
           colloids_info_t* cinfo);
int psi_force_gradmu_es(psi_t* psi, fe_t* fe, field_t* phi, hydro_t* hydro,
      colloids_info_t* cinfo);
/*CHANGE INIT - 20260422 psi_force_gradmu_e_ewald with kernel parameter */
int psi_force_gradmu_e_ewald(psi_t* psi, fe_t* fe, hydro_t* hydro,
           colloids_info_t* cinfo, subgrid_kernel_t kernel);
int psi_force_gradmu_e_ewald_offset(psi_t* psi, fe_t* fe, hydro_t* hydro,
           colloids_info_t* cinfo, subgrid_kernel_t kernel, double mesh_offset);
/*CHANGE END - 20260422 */

/*****************************************************************************
 *
 *  psi_force_gradmu
 *
 *  This routine computes the force on the fluid via the gradient
 *  of the chemical potential.
 *
 *****************************************************************************/

 /*CHANGE INIT - 20260422 kernel parameter in psi_force_gradmu */
int psi_force_gradmu(psi_t* psi, fe_t* fe, field_t* phi,
         hydro_t* hydro,
         map_t* map, colloids_info_t* cinfo,
         subgrid_kernel_t kernel) {

  assert(fe);

  switch (fe->id) {
  case FE_ELECTRO:
    psi_force_gradmu_e(psi, fe, hydro, cinfo);
    // psi_force_gradmu_e_ewald(psi, fe, hydro, cinfo, kernel);
    break;
  case FE_ELECTRO_EWALD:
    psi_force_gradmu_e_ewald(psi, fe, hydro, cinfo, kernel);
    break;
  case FE_ELECTRO_SYMMETRIC:
    psi_force_gradmu_es(psi, fe, phi, hydro, cinfo);
    break;
  default:
    pe_fatal(psi->pe, "Wrong free energy\n");
  }

  return 0;
}
/*CHANGE END - 20260422 kernel parameter in psi_force_gradmu */

/*CHANGE INIT - 20260426 psi_force_gradmu_offset for interlacing */
int psi_force_gradmu_offset(psi_t* psi, fe_t* fe, field_t* phi,
         hydro_t* hydro,
         map_t* map, colloids_info_t* cinfo,
         subgrid_kernel_t kernel, double mesh_offset) {

  assert(fe);

  switch (fe->id) {
  case FE_ELECTRO:
  case FE_ELECTRO_EWALD:
    psi_force_gradmu_e_ewald_offset(psi, fe, hydro, cinfo, kernel, mesh_offset);
    break;
  case FE_ELECTRO_SYMMETRIC:
    psi_force_gradmu_es(psi, fe, phi, hydro, cinfo);
    break;
  default:
    pe_fatal(psi->pe, "Wrong free energy\n");
  }

  return 0;
}
/*CHANGE END - 20260426 psi_force_gradmu_offset for interlacing */

/*****************************************************************************
 *
 *  psi_force_gradmu_e
 *
 *  The first of two versions, this one for FE_ELECTRO.
 *  There is some repetition of code which could be rationalised.
 *
 *  If hydro is NULL, there is no force on the fluid, but there
 *  can be a force on the colloids.
 *
 *****************************************************************************/
  // INIT VERSION - Peskin subgrid only
 int psi_force_gradmu_e(psi_t* psi, fe_t* fe, hydro_t* hydro,
            colloids_info_t* cinfo) {

     int ic, jc, kc;
     int ia;
     int nlocal[3];
     int index;
     int xs, ys, zs;       /* Coordinate strides */
     double rho_elec;      /* Species and electric charge density */
     double e[3];          /* Total electric field */
     double kt, eunit, reunit;
     double force[3];
     /* Cummulative forces for momentum correction */
     double flocal[4] = { 0.0, 0.0, 0.0, 0.0 };
     // CHANGE INIT - Subgrid charge
     klein_t flocal_k[3];
     flocal_k[X] = klein_zero();
     flocal_k[Y] = klein_zero();
     flocal_k[Z] = klein_zero();
     // CHANGE END - Subgrid charge
     double fsum[4];


     physics_t* phys = NULL;
     MPI_Comm comm;

     colloid_t* pc = NULL;

     assert(fe);
     assert(psi);
     assert(cinfo);

     cs_nlocal(psi->cs, nlocal);
     cs_strides(psi->cs, &xs, &ys, &zs);
     cs_cart_comm(psi->cs, &comm);

     physics_ref(&phys);
     physics_kt(phys, &kt);
     psi_unit_charge(psi, &eunit);
     reunit = 1.0 / eunit;

     for (ic = 1; ic <= nlocal[X]; ic++) {
       for (jc = 1; jc <= nlocal[Y]; jc++) {
         for (kc = 1; kc <= nlocal[Z]; kc++) {

           index = cs_index(psi->cs, ic, jc, kc);
           colloids_info_map(cinfo, index, &pc);

           /* Contribution from ionic electrostatic part
                    Note: The sum over the ionic species and the
                          gradient of the electrostatic potential
                          are implicitly calculated */

           psi_rho_elec(psi, index, &rho_elec);
           // CHANGE INIT - Subgrid charge
           psi_electric_field(psi, index, e);
           // CHANGE
           psi->efield->data[addr_rank1(psi->efield->nsites, psi->efield->nf, index, X)] = kt*e[X]; //kt*e[X];
           psi->efield->data[addr_rank1(psi->efield->nsites, psi->efield->nf, index, Y)] = kt*e[Y];
           psi->efield->data[addr_rank1(psi->efield->nsites, psi->efield->nf, index, Z)] = kt*e[Z];

          //  double phi_node = 0.0;
          //  double r[3];
          //  double dr;
          //  double r0[3] = { ic, jc, kc }; // Position of the lattice site
          //  int i, j, k, i_min, i_max, j_min, j_max, k_min, k_max;
          //  double E_field[3] = { 0.0, 0.0, 0.0 };      /* Electric field on particle from this lattice site */
          //  klein_t E_field_k[3];
          //  E_field_k[X] = klein_zero();
          //  E_field_k[Y] = klein_zero();
          //  E_field_k[Z] = klein_zero();

          //  subgrid_get_lattice_index(r0, nlocal, &i_min, &i_max, &j_min, &j_max, &k_min, &k_max);

          //  for (i = i_min; i <= i_max; i++) {
          //    for (j = j_min; j <= j_max; j++) {
          //      for (k = k_min; k <= k_max; k++) {

          //        index = cs_index(psi->cs, i, j, k);

          //        /* Separation between r0 and the lattice site */
          //        r[X] = r0[X] - 1.0 * i;
          //        r[Y] = r0[Y] - 1.0 * j;
          //        r[Z] = r0[Z] - 1.0 * k;

          //        dr = d_peskin(r[X]) * d_peskin(r[Y]) * d_peskin(r[Z]);

          //        /* Electric field at this lattice site */
          //        psi_psi(psi, index, &phi_node);

          //        // pe_info(pe, "campo en el indice %d : Ex=%.20f   Ey=%.20f   Ez=%.20f\n", index, e[X], e[Y], e[Z]);

          //        /* Field on particle from electric field at this site */
          //        E_field[X] = phi_node * d_peskin_derivative(r[X]) * d_peskin(r[Y]) * d_peskin(r[Z]);
          //        E_field[Y] = phi_node * d_peskin(r[X]) * d_peskin_derivative(r[Y]) * d_peskin(r[Z]);
          //        E_field[Z] = phi_node * d_peskin(r[X]) * d_peskin(r[Y]) * d_peskin_derivative(r[Z]);

          //        /* Add to particle field */
          //        klein_add_double(&E_field_k[X], E_field[X]);
          //        klein_add_double(&E_field_k[Y], E_field[Y]);
          //        klein_add_double(&E_field_k[Z], E_field[Z]);

          //      }
          //    }
          //  }

          //  e[X] = klein_sum(&E_field_k[X]);
          //  e[Y] = klein_sum(&E_field_k[Y]);
          //  e[Z] = klein_sum(&E_field_k[Z]);

          //  psi->efield->data[addr_rank1(psi->efield->nsites, psi->efield->nf, index, X)] = e[X] ;
          //  psi->efield->data[addr_rank1(psi->efield->nsites, psi->efield->nf, index, Y)] = e[Y];
          //  psi->efield->data[addr_rank1(psi->efield->nsites, psi->efield->nf, index, Z)] = e[Z];

           //CHANGE END - Subgrid charge

           for (ia = 0; ia < 3; ia++) {
             e[ia] *= kt * reunit;
             force[ia] = rho_elec * e[ia];
           }

           /* If solid, accumulate contribution to colloid;
              otherwise to fluid node */

           if (pc) {
             pc->force[X] += force[X];
             pc->force[Y] += force[Y];
             pc->force[Z] += force[Z];
           }
           else {
             if (hydro) hydro_f_local_add(hydro, index, force);
             flocal[3] += 1.0;
           }

           /* Accumulate contribution to total force on system */

           // // CHANGE INIT - Subgrid charge
           // // flocal[X] += force[X];
           // // flocal[Y] += force[Y];
           // // flocal[Z] += force[Z];
           // klein_add_double(&flocal_k[X], force[X]);
           // klein_add_double(&flocal_k[Y], force[Y]);
           // klein_add_double(&flocal_k[Z], force[Z]);
           // // CHANGE END - Subgrid charge

         }
       }
     }

     // // CHANGE INIT - Subgrid charge
     // flocal[X] = klein_sum(&flocal_k[X]);
     // flocal[Y] = klein_sum(&flocal_k[Y]);
     // flocal[Z] = klein_sum(&flocal_k[Z]);
     // // CHANGE END - Subgrid charge

     // MPI_Allreduce(flocal, fsum, 4, MPI_DOUBLE, MPI_SUM, comm);

     // fsum[X] /= fsum[3];
     // fsum[Y] /= fsum[3];
     // fsum[Z] /= fsum[3];

     // /* Now actually compute the force on the fluid with the correction
     //    (based on number of fluid nodes) and store */

     // for (ic = 1; ic <= nlocal[X]; ic++) {
     //   for (jc = 1; jc <= nlocal[Y]; jc++) {
     //     for (kc = 1; kc <= nlocal[Z]; kc++) {

     //       index = cs_index(psi->cs, ic, jc, kc);

     //       colloids_info_map(cinfo, index, &pc);
     //       if (pc) continue;

     //       force[X] = -fsum[X];
     //       force[Y] = -fsum[Y];
     //       force[Z] = -fsum[Z];

     //       if (hydro) hydro_f_local_add(hydro, index, force);
     //     }
     //   }
     // }

     return 0;
   }
   // END VERSION - Peskin subgrid only

//    // INIT VERSION - Peskin-all
// int psi_force_gradmu_e(psi_t* psi, fe_t* fe, hydro_t* hydro,
//            colloids_info_t* cinfo) {

//   int ic, jc, kc;
//   int ia;
//   int nlocal[3];
//   int index;
//   int xs, ys, zs;       /* Coordinate strides */
//   double rho_elec;      /* Species and electric charge density */
//   double e[3];          /* Total electric field */
//   double kt, eunit, reunit;
//   double force[3];
//   /* Cummulative forces for momentum correction */
//   double flocal[4] = { 0.0, 0.0, 0.0, 0.0 };
//   // CHANGE INIT - Subgrid charge
//   klein_t flocal_k[3];
//   flocal_k[X] = klein_zero();
//   flocal_k[Y] = klein_zero();
//   flocal_k[Z] = klein_zero();
//   // CHANGE END - Subgrid charge
//   double fsum[4];


//   physics_t* phys = NULL;
//   MPI_Comm comm;

//   colloid_t* pc = NULL;

//   assert(fe);
//   assert(psi);
//   assert(cinfo);

//   cs_nlocal(psi->cs, nlocal);
//   cs_strides(psi->cs, &xs, &ys, &zs);
//   cs_cart_comm(psi->cs, &comm);

//   physics_ref(&phys);
//   physics_kt(phys, &kt);
//   psi_unit_charge(psi, &eunit);
//   reunit = 1.0 / eunit;

//   for (ic = 1; ic <= nlocal[X]; ic++) {
//     for (jc = 1; jc <= nlocal[Y]; jc++) {
//       for (kc = 1; kc <= nlocal[Z]; kc++) {
//   // for (ic = 0; ic <= nlocal[X]+1; ic++) {
//   //   for (jc = 0; jc <= nlocal[Y]+1; jc++) {
//   //     for (kc = 0; kc <= nlocal[Z]+1; kc++) {

//         /*CHANGE INIT - psi_force_gradmu_e fix: preserve node index */
//         int index0 = cs_index(psi->cs, ic, jc, kc);

//         psi_rho_elec(psi, index0, &rho_elec);

//         double phi_node = 0.0;
//         double r[3];
//         double r0[3] = { ic, jc, kc };
//         int i, j, k, i_min, i_max, j_min, j_max, k_min, k_max;
//         klein_t E_field_k[3];
//         E_field_k[X] = klein_zero();
//         E_field_k[Y] = klein_zero();
//         E_field_k[Z] = klein_zero();

//         subgrid_get_lattice_index(r0, nlocal, &i_min, &i_max, &j_min, &j_max, &k_min, &k_max);

//         // for (i = i_min; i <= i_max; i++) {
//         //   for (j = j_min; j <= j_max; j++) {
//         //     for (k = k_min; k <= k_max; k++) {
//         for (i = i_min-1; i <= i_max+1; i++) {
//           for (j = j_min-1; j <= j_max+1; j++) {
//             for (k = k_min-1; k <= k_max+1; k++) {
//               index = cs_index(psi->cs, i, j, k);

//               r[X] = r0[X] - 1.0 * i;
//               r[Y] = r0[Y] - 1.0 * j;
//               r[Z] = r0[Z] - 1.0 * k;

//               // psi_psi(psi, index, &phi_node);

//               // klein_add_double(&E_field_k[X], phi_node * d_peskin_derivative(r[X]) * d_peskin(r[Y]) * d_peskin(r[Z]));
//               // klein_add_double(&E_field_k[Y], phi_node * d_peskin(r[X]) * d_peskin_derivative(r[Y]) * d_peskin(r[Z]));
//               // klein_add_double(&E_field_k[Z], phi_node * d_peskin(r[X]) * d_peskin(r[Y]) * d_peskin_derivative(r[Z]));

//               // pe_info(pe, "campo en el indice %d : Ex=%.20f   Ey=%.20f   Ez=%.20f\n", index, e[X], e[Y], e[Z]);

//               /* Electric field at this lattice site */
//               psi_electric_field(psi, index, e);

//               /* Field on particle from electric field at this site */
//               double dr = d_peskin(r[X]) * d_peskin(r[Y]) * d_peskin(r[Z]);
//               // Usando B-SPLINE
//               // double dr = d_bspline4(r[X]) * d_bspline4(r[Y]) * d_bspline4(r[Z]); 

//               /* E = -grad phi via Peskin derivative kernel */
//               klein_add_double(&E_field_k[X], e[X] * dr);
//               klein_add_double(&E_field_k[Y], e[Y] * dr);
//               klein_add_double(&E_field_k[Z], e[Z] * dr);
//             }
//           }
//         }

//         e[X] = klein_sum(&E_field_k[X]);
//         e[Y] = klein_sum(&E_field_k[Y]);
//         e[Z] = klein_sum(&E_field_k[Z]);

//         psi->efield->data[addr_rank1(psi->efield->nsites, psi->efield->nf, index0, X)] = kt * e[X];
//         psi->efield->data[addr_rank1(psi->efield->nsites, psi->efield->nf, index0, Y)] = kt * e[Y];
//         psi->efield->data[addr_rank1(psi->efield->nsites, psi->efield->nf, index0, Z)] = kt * e[Z];
//         /*CHANGE END - psi_force_gradmu_e fix */

//         for (ia = 0; ia < 3; ia++) {
//           e[ia] *= kt * reunit;
//           force[ia] = rho_elec * e[ia];
//         }

//         /* fluid node */
//         if (hydro) hydro_f_local_add(hydro, index0, force);

//       }
//     }
//   }

//   return 0;
// }
// // END VERSION - Peskin-all

/*CHANGE INIT - 20260422 psi_force_gradmu_e_ewald */
/*****************************************************************************
 *
 *  psi_force_gradmu_e_ewald
 *
 *  Version of psi_force_gradmu_e with selectable interpolation kernel.
 *  The gather loop ranges are computed from the kernel support instead of
 *  the fixed drange_ used by subgrid_get_lattice_index.
 *
 *****************************************************************************/
int psi_force_gradmu_e_ewald(psi_t* psi, fe_t* fe, hydro_t* hydro,
                              colloids_info_t* cinfo, subgrid_kernel_t kernel) {

  int ic, jc, kc;
  int ia;
  int nlocal[3];
  int index;
  int xs, ys, zs;
  double rho_elec;
  double e[3];
  double kt, eunit, reunit;
  double force[3];

  physics_t* phys = NULL;

  assert(fe);
  assert(psi);
  assert(cinfo);

  cs_nlocal(psi->cs, nlocal);
  cs_strides(psi->cs, &xs, &ys, &zs);

  physics_ref(&phys);
  physics_kt(phys, &kt);
  psi_unit_charge(psi, &eunit);
  reunit = 1.0 / eunit;

  // range for used kernel, e.g. 2 for B-spline4, 3 for B-spline6, 2 for Peskin4
  int krange = subgrid_get_range(kernel);

  for (ic = 1; ic <= nlocal[X]; ic++) {
    for (jc = 1; jc <= nlocal[Y]; jc++) {
      for (kc = 1; kc <= nlocal[Z]; kc++) {

        int index0 = cs_index(psi->cs, ic, jc, kc);

        psi_rho_elec(psi, index0, &rho_elec);

        double r[3];
        double r0[3] = { ic, jc, kc };
        int i, j, k, i_min, i_max, j_min, j_max, k_min, k_max;
        klein_t E_field_k[3];
        E_field_k[X] = klein_zero();
        E_field_k[Y] = klein_zero();
        E_field_k[Z] = klein_zero();

        subgrid_get_lattice_index_range_halo(r0, krange, nlocal, &i_min, &i_max, &j_min, &j_max, &k_min, &k_max);

        for (i = i_min; i <= i_max; i++) {
          for (j = j_min; j <= j_max; j++) {
            for (k = k_min; k <= k_max; k++) {

              index = cs_index(psi->cs, i, j, k);

              r[X] = r0[X] - 1.0 * i;
              r[Y] = r0[Y] - 1.0 * j;
              r[Z] = r0[Z] - 1.0 * k;

              double dr;
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
                dr = d_kb4(r[X]) * d_kb4(r[Y]) * d_kb4(r[Z]) / subgrid_kb4_norm_fluid();
              }
              else if (kernel == SUBGRID_KERNEL_PESKIN6) {
                dr = d_peskin6(r[X]) * d_peskin6(r[Y]) * d_peskin6(r[Z]);
              }

              psi_electric_field(psi, index, e);

              klein_add_double(&E_field_k[X], e[X] * dr);
              klein_add_double(&E_field_k[Y], e[Y] * dr);
              klein_add_double(&E_field_k[Z], e[Z] * dr);
            }
          }
        }

        e[X] = klein_sum(&E_field_k[X]);
        e[Y] = klein_sum(&E_field_k[Y]);
        e[Z] = klein_sum(&E_field_k[Z]);

        psi->efield->data[addr_rank1(psi->efield->nsites, psi->efield->nf, index0, X)] = kt * e[X];
        psi->efield->data[addr_rank1(psi->efield->nsites, psi->efield->nf, index0, Y)] = kt * e[Y];
        psi->efield->data[addr_rank1(psi->efield->nsites, psi->efield->nf, index0, Z)] = kt * e[Z];

        for (ia = 0; ia < 3; ia++) {
          e[ia] *= kt * reunit;
          force[ia] = rho_elec * e[ia];
        }

        if (hydro) hydro_f_local_add(hydro, index0, force);
      }
    }
  }

  return 0;
}
/*CHANGE END - 20260422 psi_force_gradmu_e_ewald */

/*CHANGE INIT - 20260426 psi_force_gradmu_e_ewald_offset */
/*****************************************************************************
 *
 *  psi_force_gradmu_e_ewald_offset
 *
 *  Fluid gather with mesh_offset for interlacing. Adjoint of the scatter:
 *  each fluid source node (si,sj,sk) with original charge rho_saved acts as
 *  if it is at position r_src = {si + offset, sj, sk}. The interpolated
 *  E-field at that position is gathered from the surrounding lattice nodes
 *  using the kernel weights w(r_src - j). The force F = rho_saved * E_interp
 *  is deposited on the source node (si,sj,sk).
 *
 *  rho_saved: fluid-only charge density before scatter (nsites * nk array)
 *
 *****************************************************************************/
int psi_force_gradmu_e_ewald_offset(psi_t* psi, fe_t* fe, hydro_t* hydro,
                              colloids_info_t* cinfo, subgrid_kernel_t kernel,
                              double mesh_offset) {

  int si, sj, sk;
  int ia;
  int nlocal[3];
  int xs, ys, zs;
  double e[3];
  double kt, eunit, reunit;
  double force[3];

  physics_t* phys = NULL;

  assert(fe);
  assert(psi);
  assert(cinfo);

  cs_nlocal(psi->cs, nlocal);
  cs_strides(psi->cs, &xs, &ys, &zs);

  physics_ref(&phys);
  physics_kt(phys, &kt);
  psi_unit_charge(psi, &eunit);
  reunit = 1.0 / eunit;

  int krange = subgrid_get_range(kernel);

  /* Mirror the scatter source selection — only X axis has offset. */
  double frac_x = mesh_offset - floor(mesh_offset);
  int ix_excl_lo = 0, ix_excl_hi = -1;
  int ix_halo_lo = 1, ix_halo_hi = 0;
  if (frac_x > 0.0) {
    ix_excl_lo = nlocal[X] - krange + 1; ix_excl_hi = nlocal[X];
    ix_halo_lo = 1 - krange;             ix_halo_hi = 0;
  } else if (frac_x < 0.0) {
    ix_excl_lo = 1;               ix_excl_hi = krange;
    ix_halo_lo = nlocal[X] + 1;  ix_halo_hi = nlocal[X] + krange;
  }

  int list_cap = nlocal[X] + 2 * krange + 4;
  int *i_list = (int*)malloc(list_cap*sizeof(int)); int i_list_n = 0;
  for (si = 1; si <= nlocal[X]; si++) { if (si < ix_excl_lo || si > ix_excl_hi) i_list[i_list_n++] = si; }
  for (si = ix_halo_lo; si <= ix_halo_hi; si++) i_list[i_list_n++] = si;

  for (int ii = 0; ii < i_list_n; ii++) {
  si = i_list[ii];
  int si_dest = si; if (si_dest < 1) si_dest += nlocal[X]; else if (si_dest > nlocal[X]) si_dest -= nlocal[X];
  for (sj = 1; sj <= nlocal[Y]; sj++) {
  for (sk = 1; sk <= nlocal[Z]; sk++) {

      int index0 = cs_index(psi->cs, si_dest, sj, sk);

        double rho_elec;
        psi_rho_elec(psi, index0, &rho_elec);
        if (rho_elec == 0.0) continue;

        /* Source position: X shifted, Y/Z not shifted — same as scatter */
        double r_src[3] = { (double)si + mesh_offset, (double)sj, (double)sk };

        int i, j, k, i_min, i_max, j_min, j_max, k_min, k_max;
        klein_t E_field_k[3];
        E_field_k[X] = klein_zero();
        E_field_k[Y] = klein_zero();
        E_field_k[Z] = klein_zero();

        subgrid_get_lattice_index_range_halo(r_src, krange, nlocal,
                                             &i_min, &i_max, &j_min, &j_max, &k_min, &k_max);

        for (i = i_min; i <= i_max; i++) {
          for (j = j_min; j <= j_max; j++) {
            for (k = k_min; k <= k_max; k++) {

              /* Wrap halo indices to interior for periodic BC */
              int iw = i, jw = j, kw = k;
              if (iw < 1) iw += nlocal[X]; else if (iw > nlocal[X]) iw -= nlocal[X];
              if (jw < 1) jw += nlocal[Y]; else if (jw > nlocal[Y]) jw -= nlocal[Y];
              if (kw < 1) kw += nlocal[Z]; else if (kw > nlocal[Z]) kw -= nlocal[Z];
              int index = cs_index(psi->cs, iw, jw, kw);

              double rx = r_src[X] - (double)i;
              double ry = r_src[Y] - (double)j;
              double rz = r_src[Z] - (double)k;

              double dr;
              if (kernel == SUBGRID_KERNEL_BSPLINE6) {
                dr = d_bspline6(rx) * d_bspline6(ry) * d_bspline6(rz);
              }
              else if (kernel == SUBGRID_KERNEL_BSPLINE4) {
                dr = d_bspline4(rx) * d_bspline4(ry) * d_bspline4(rz);
              }
              else if (kernel == SUBGRID_KERNEL_PESKIN4) {
                dr = d_peskin(rx) * d_peskin(ry) * d_peskin(rz);
              }
              else if (kernel == SUBGRID_KERNEL_KB4) {
                dr = d_kb4(rx) * d_kb4(ry) * d_kb4(rz) / subgrid_kb4_norm_fluid();
              }
              else if (kernel == SUBGRID_KERNEL_PESKIN6) {
                dr = d_peskin6(rx) * d_peskin6(ry) * d_peskin6(rz);
              }
              /*CHANGE INIT - 20260710 trilinear: pointwise E at integer offsets */
              else if (kernel == SUBGRID_KERNEL_TRILINEAR) {
                dr = d_trilinear(rx) * d_trilinear(ry) * d_trilinear(rz);
              }
              /*CHANGE END - 20260710 */
              else { dr = 0.0; }

              if (dr == 0.0) continue;

              psi_electric_field(psi, index, e);

              klein_add_double(&E_field_k[X], e[X] * dr);
              klein_add_double(&E_field_k[Y], e[Y] * dr);
              klein_add_double(&E_field_k[Z], e[Z] * dr);
            }
          }
        }

        e[X] = klein_sum(&E_field_k[X]);
        e[Y] = klein_sum(&E_field_k[Y]);
        e[Z] = klein_sum(&E_field_k[Z]);

        /* Store smoothed E-field at source node */
        psi->efield->data[addr_rank1(psi->efield->nsites, psi->efield->nf, index0, X)] = kt * e[X];
        psi->efield->data[addr_rank1(psi->efield->nsites, psi->efield->nf, index0, Y)] = kt * e[Y];
        psi->efield->data[addr_rank1(psi->efield->nsites, psi->efield->nf, index0, Z)] = kt * e[Z];

        for (ia = 0; ia < 3; ia++) {
          e[ia] *= kt * reunit;
          force[ia] = rho_elec * e[ia];
        }

        if (hydro) hydro_f_local_add(hydro, index0, force);
  } /* sk */
  } /* sj */
  } /* ii / si */

  free(i_list);
  return 0;
}
/*CHANGE END - 20260426 psi_force_gradmu_e_ewald_offset */

/*****************************************************************************
 *
 *  psi_force_gradmu_es
 *
 *  This for FE_ELECTRO_SYMMETRIC, including the solvation and
 *  composition-dependent terms.
 *
 *  Note: The ionic solvation free energy difference is in units
 *        of kt and must be dressed for the force calculation.
 *
 *****************************************************************************/

int psi_force_gradmu_es(psi_t* psi, fe_t* fe, field_t* phi, hydro_t* hydro,
      colloids_info_t* cinfo) {

  int ic, jc, kc;
  int in, nk;
  int ia;
  int nlocal[3];
  int index;
  int xs, ys, zs;       /* Coordinate strides */
  double rho, rho_elec; /* Species and electric charge density */
  double e[3];          /* Total electric field */
  double muphim1, muphip1, musm1, musp1;
  double phi0;          /* Compositional order parameter */
  double kt, eunit, reunit;
  double force[3];
  /* Cummulative forces for momentum correction */
  double flocal[4] = { 0.0, 0.0, 0.0, 0.0 };
  double fsum[4];

  physics_t* phys = NULL;
  MPI_Comm comm;

  colloid_t* pc = NULL;

  assert(fe);
  assert(psi);
  assert(phi);
  assert(cinfo);

  cs_nlocal(psi->cs, nlocal);
  cs_strides(psi->cs, &xs, &ys, &zs);
  cs_cart_comm(psi->cs, &comm);

  physics_ref(&phys);
  physics_kt(phys, &kt);
  psi_unit_charge(psi, &eunit);
  reunit = 1.0 / eunit;

  psi_nk(psi, &nk);
  assert(nk == 2); /* This routine is not completely general */

  for (ic = 1; ic <= nlocal[X]; ic++) {
    for (jc = 1; jc <= nlocal[Y]; jc++) {
      for (kc = 1; kc <= nlocal[Z]; kc++) {

        index = cs_index(psi->cs, ic, jc, kc);
        colloids_info_map(cinfo, index, &pc);

        /* X-direction */

        field_scalar(phi, index, &phi0);

        /* Contribution from composition part */
        fe->func->mu(fe, index - xs, &muphim1);
        fe->func->mu(fe, index + xs, &muphip1);

        force[X] = -phi0 * 0.5 * (muphip1 - muphim1);

        /* Contribution from ionic solvation part */
        for (in = 0; in < nk; in++) {
          psi_rho(psi, index, in, &rho);
          fe->func->mu_solv(fe, index - xs, in, &musm1);
          fe->func->mu_solv(fe, index + xs, in, &musp1);
          force[X] -= rho * 0.5 * (musp1 - musm1);
        }

        /* Y-direction */

        /* Contribution from composition part */
        fe->func->mu(fe, index - ys, &muphim1);
        fe->func->mu(fe, index + ys, &muphip1);
        force[Y] = -phi0 * 0.5 * (muphip1 - muphim1);

        /* Contribution from ionic solvation part */

        for (in = 0; in < nk; in++) {
          psi_rho(psi, index, in, &rho);
          fe->func->mu_solv(fe, index - ys, in, &musm1);
          fe->func->mu_solv(fe, index + ys, in, &musp1);
          force[Y] -= rho * 0.5 * (musp1 - musm1);
        }

        /* Z-direction */
        /* Contribution from composition part */

        fe->func->mu(fe, index - zs, &muphim1);
        fe->func->mu(fe, index + zs, &muphip1);
        force[Z] = -phi0 * 0.5 * (muphip1 - muphim1);

        /* Contribution from ionic solvation part */
        for (in = 0; in < nk; in++) {
          psi_rho(psi, index, in, &rho);
          fe->func->mu_solv(fe, index - zs, in, &musm1);
          fe->func->mu_solv(fe, index + zs, in, &musp1);
          force[Z] -= rho * 0.5 * (musp1 - musm1);
        }

        /* Contribution from ionic electrostatic part
                 Note: The sum over the ionic species and the
                       gradient of the electrostatic potential
                       are implicitly calculated */

        psi_rho_elec(psi, index, &rho_elec);
        psi_electric_field(psi, index, e);

        for (ia = 0; ia < 3; ia++) {
          e[ia] *= kt * reunit;
          force[ia] += rho_elec * e[ia];
        }

        /* If solid, accumulate contribution to colloid;
           otherwise to fluid node */

        if (pc) {

          pc->force[X] += force[X];
          pc->force[Y] += force[Y];
          pc->force[Z] += force[Z];

        }
        else {
          if (hydro) hydro_f_local_add(hydro, index, force);
          flocal[3] += 1.0;
        }

        /* Accumulate contribution to total force on system */

        flocal[X] += force[X];
        flocal[Y] += force[Y];
        flocal[Z] += force[Z];

      }
    }
  }

  MPI_Allreduce(flocal, fsum, 4, MPI_DOUBLE, MPI_SUM, comm);

  fsum[X] /= fsum[3];
  fsum[Y] /= fsum[3];
  fsum[Z] /= fsum[3];

  /* Now actually compute the force on the fluid with the correction
     (based on number of fluid nodes) and store */

  for (ic = 1; ic <= nlocal[X]; ic++) {
    for (jc = 1; jc <= nlocal[Y]; jc++) {
      for (kc = 1; kc <= nlocal[Z]; kc++) {

        index = cs_index(psi->cs, ic, jc, kc);

        colloids_info_map(cinfo, index, &pc);
        if (pc) continue;

        force[X] = -fsum[X];
        force[Y] = -fsum[Y];
        force[Z] = -fsum[Z];

        if (hydro) hydro_f_local_add(hydro, index, force);
      }
    }
  }

  return 0;
}

/*****************************************************************************
 *
 *  psi_force_divstress
 *
 *  A routine for force via divergence of stress, allowing
 *  the stress to be computed inside the colloids.
 *
 *  The stress is to include the full electric field.
 *
 *****************************************************************************/

int psi_force_divstress(psi_t* psi, fe_t* fe, hydro_t* hydro,
      colloids_info_t* cinfo) {

  int nlocal[3] = { 0 };
  cs_t* cs = NULL;
  stencil_t* s = NULL;

  assert(psi);
  assert(cinfo);

  cs = psi->cs;
  s = psi->stencil;
  assert(cs);
  assert(s);

  cs_nlocal(cs, nlocal);

  for (int ic = 1; ic <= nlocal[X]; ic++) {
    for (int jc = 1; jc <= nlocal[Y]; jc++) {
      for (int kc = 1; kc <= nlocal[Z]; kc++) {

        int index = cs_index(cs, ic, jc, kc);
        double force[3] = { 0 };
        colloid_t* pc = NULL;

        colloids_info_map(cinfo, index, &pc);

        /* Calculate divergence based on the stencil */
        for (int p = 1; p < s->npoints; p++) {

          int8_t cx = s->cv[p][X];
          int8_t cy = s->cv[p][Y];
          int8_t cz = s->cv[p][Z];
          int index1 = cs_index(cs, ic + cx, jc + cy, kc + cz);
          double pth[3][3] = { 0 };

          fe->func->stress(fe, index1, pth);

          for (int ia = 0; ia < 3; ia++) {
            for (int ib = 0; ib < 3; ib++) {
              force[ia] -= s->wgradients[p] * pth[ia][ib] * s->cv[p][ib];
            }
          }
        }

        /* Store the force on the colloid or on the lattice */

        if (pc) {
          pc->force[X] += force[X];
          pc->force[Y] += force[Y];
          pc->force[Z] += force[Z];
        }
        else {
          hydro_f_local_add(hydro, index, force);
        }

      }
    }
  }

  return 0;
}

// CHANGE INIT - Subgrid charge
/*****************************************************************************
 *
 *  psi_force_gradmu_e_subgrid
 *
 *  Calculate electric force on fluid when there is only subgrid particles
 *
 *****************************************************************************/

int psi_force_gradmu_e_subgrid(psi_t* psi, fe_t* fe, hydro_t* hydro,
           colloids_info_t* cinfo) {

  int ic, jc, kc;
  int ia;
  int nlocal[3];
  int index;
  // int xs, ys, zs;          /* Coordinate strides */
  double rho_elec;            /* Species and electric charge density */
  double e[3];                /* Total electric field */
  double kt, eunit, reunit;
  double force[3];


  physics_t* phys = NULL;
  MPI_Comm comm;

  colloid_t* pc = NULL;

  assert(fe);
  assert(psi);
  assert(cinfo);

  cs_nlocal(psi->cs, nlocal);
  // cs_strides(psi->cs, &xs, &ys, &zs);
  // cs_cart_comm(psi->cs, &comm);

  physics_ref(&phys);
  physics_kt(phys, &kt);
  psi_unit_charge(psi, &eunit);
  reunit = 1.0 / eunit;

  for (ic = 1; ic <= nlocal[X]; ic++) {
    for (jc = 1; jc <= nlocal[Y]; jc++) {
      for (kc = 1; kc <= nlocal[Z]; kc++) {

        index = cs_index(psi->cs, ic, jc, kc);

        psi_rho_elec(psi, index, &rho_elec);
        psi_electric_field(psi, index, e);

        for (ia = 0; ia < 3; ia++) {
          e[ia] *= kt * reunit;
          force[ia] = rho_elec * e[ia];
        }

        if (hydro) hydro_f_local_add(hydro, index, force);

      }
    }
  }

  return 0;
}
// CHANGE END - Subgrid charge