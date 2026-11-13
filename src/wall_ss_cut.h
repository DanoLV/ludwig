/*****************************************************************************
 *
 *  wall_ss_cut.h
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  (c) 2014-2022 The University of Edinburgh
 *  Juho Lintuvuori (juho.lintuvuori@u-psud.fr)
 *
 *****************************************************************************/

#ifndef LUDWIG_WALL_SS_CUT_H
#define LUDWIG_WALL_SS_CUT_H

#include "pe.h"
#include "coords.h"
#include "wall.h"
#include "colloids.h"
#include "interaction.h"

typedef struct wall_ss_cut_options_s wall_ss_cut_options_t;
typedef struct wall_ss_cut_s wall_ss_cut_t;

struct wall_ss_cut_options_s {
  double epsilon;       /* Energy scale */
  double sigma;         /* Sigma (length) */
  double nu;            /* exponent */
  double hc;            /* Surface-surface cut off */
  /*CHANGE INIT - 20261018 wall planes given explicitly (map-based walls).
   * The original potential only acts on the boundary walls of the wall_t
   * object (isboundary); constant-potential walls built with
   * porous_media_init wall_x|y|z are map sites and wall_t knows nothing of
   * them. Original struct ended at 'hc'. */
  int    use_planes;    /* 1: use the planes below instead of isboundary */
  int    axis;          /* Axis normal to the planes (0 = x, 1 = y, 2 = z) */
  double plane_lo;      /* Coordinate of the lower wall plane */
  double plane_hi;      /* Coordinate of the upper wall plane */
  double radius;        /* Sphere radius: h = distance to plane - radius */
  /*CHANGE END - 20261018 */
  /*CHANGE INIT - 20261021 alternative wall model: purely repulsive 10-4-3
   * (Steele) ion-wall potential of Hoang Ngoc Minh et al., J. Chem. Phys.
   * 158, 104103 (2023), Eqs. (7)-(8). model 0 = original soft sphere
   * (epsilon, sigma, nu, hc); model 1 = Steele, using epsilon = eps_w and
   * sigma = sigma_w (nu and hc are ignored; the cut-off is z* sigma_w). */
  int    model;         /* 0 soft sphere (default), 1 Steele repulsive */
  /*CHANGE END - 20261021 */
};

int wall_ss_cut_create(pe_t * pe, cs_t * cs, wall_t * wall,
		       const wall_ss_cut_options_t * opts,
		       wall_ss_cut_t ** pobj);
int wall_ss_cut_free(wall_ss_cut_t * obj);
int wall_ss_cut_info(wall_ss_cut_t * obj);
int wall_ss_cut_register(wall_ss_cut_t * obj, interact_t * parent);
int wall_ss_cut_compute(colloids_info_t * cinfo, void * self);
int wall_ss_cut_stats(void * self, double * stats);
int wall_ss_cut_single(wall_ss_cut_t * obj, double h, double * f, double * v);

#endif
