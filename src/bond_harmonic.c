//CHANGE3
/*****************************************************************************
 *
 *  bond_harmonic.c
 *
 *  Finite extensible elastic bond
 *
 *   V(r) = (1/2) k (r - r0)^2    
 *
 *  where r is the separation of bonded pairs.
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  (c) 2014-2020 The University of Edinburgh
 *
 *  Contributing authors:
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *  Kai Qi (kai.qi@epfl.ch)
 *
 *****************************************************************************/

#include <assert.h>
#include <math.h>
#include <stdlib.h>

#include "pe.h"
#include "coords.h"
#include "colloids.h"
#include "bond_harmonic.h"
#include "util.h"

struct bond_harmonic_s {
  pe_t * pe;           /* Parallel environment */
  cs_t * cs;           /* Coordinate system */
  double k[NUM_INT_PART_TYPES][NUM_INT_PART_TYPES];            /* spring constant */
  double r0[NUM_INT_PART_TYPES][NUM_INT_PART_TYPES];           /* natrual length */
  double vlocal;       /* Bond potential contribution */
  double rminlocal;    /* Minimum bond extent */
  double rmaxlocal;    /* Maximum bond extension */
  double bondlocal;    /* Number of bonds computed (double) */
};

/*****************************************************************************
 *
 *  bond_harmonic_create
 *
 *****************************************************************************/

int bond_harmonic_create(pe_t * pe, cs_t * cs, bond_harmonic_t ** pobj) {

  bond_harmonic_t * obj = NULL;

  assert(pe);
  assert(cs);
  assert(pobj);

  obj = (bond_harmonic_t *) calloc(1, sizeof(bond_harmonic_t));
  assert(obj);
  if (obj == NULL) pe_fatal(pe, "calloc(bond_harmonic_t) failed\n");

  obj->pe = pe;
  obj->cs = cs;

  *pobj = obj;

  return 0;
}

/*****************************************************************************
 *
 *  bond_harmonic_free
 *
 *****************************************************************************/

int bond_harmonic_free(bond_harmonic_t * obj) {

  assert(obj);

  free(obj);

  return 0;
}

/*****************************************************************************
 *
 *  bond_harmonic_param_set
 *
 *****************************************************************************/

int bond_harmonic_param_set(bond_harmonic_t * obj, double k[][NUM_INT_PART_TYPES], double r0[][NUM_INT_PART_TYPES]) {

  assert(obj);

  for (int i=0;i<NUM_INT_PART_TYPES;i++)
    for (int j=0;j<NUM_INT_PART_TYPES;j++) {
        obj->k[i][j] = k[i][j];
        obj->r0[i][j] = r0[i][j];
    }

  return 0;
}

/*****************************************************************************
 *
 *  bond_harmonic_info
 *
 *****************************************************************************/

int bond_harmonic_info(bond_harmonic_t * obj) {

  assert(obj);

  pe_info(obj->pe, "Harmonic bond\n");

  pe_info(obj->pe, "Spring constant:          ");
  for(int i=0;i<NUM_INT_PART_TYPES;i++) {
    for(int j=0;j<NUM_INT_PART_TYPES;j++) 
        pe_info(obj->pe, "%14.7e  ", obj->k[i][j]);
    pe_info(obj->pe, "\n                          ");
  }
  pe_info(obj->pe, "\n");

  pe_info(obj->pe, "Equilibrium separation:   ");
  for(int i=0;i<NUM_INT_PART_TYPES;i++) {
    for(int j=0;j<NUM_INT_PART_TYPES;j++) 
        pe_info(obj->pe, "%14.7e  ", obj->r0[i][j]);
    pe_info(obj->pe, "\n                          ");
  }
  pe_info(obj->pe, "\n");

  return 0;
}

/*****************************************************************************
 *
 *  bond_harmonic_register
 *
 *****************************************************************************/

int bond_harmonic_register(bond_harmonic_t * obj, interact_t * parent) {

  assert(obj);
  assert(parent);

  interact_potential_add(parent, INTERACT_BOND_HARMONIC, obj, bond_harmonic_compute);
  interact_statistic_add(parent, INTERACT_BOND_HARMONIC, obj, bond_harmonic_stats);

  double r0max=0;
  for(int i=0;i<NUM_INT_PART_TYPES;i++)
    for(int j=0;j<NUM_INT_PART_TYPES;j++)
        r0max=dmax(r0max,obj->r0[i][j]);
  interact_rc_set(parent, INTERACT_BOND_HARMONIC, 4.0*r0max);

  return 0;
}

/*****************************************************************************
 *
 *  bond_harmonic_compute
 *
 *****************************************************************************/

int bond_harmonic_compute(colloids_info_t * cinfo, void * self) {

  bond_harmonic_t * obj = (bond_harmonic_t *) self;

  int n;
  double r12[3];
  double r2min,r2max;
  double r2,f;
  double r2md;
  int it0,it1;

  colloid_t * pc = NULL;

  assert(cinfo);
  assert(obj);

  colloids_info_local_head(cinfo, &pc);

  r2min=0.0;
  for(int i=0;i<NUM_INT_PART_TYPES;i++)
    for(int j=0;j<NUM_INT_PART_TYPES;j++) 
        r2min=dmax(r2min,4*obj->r0[i][j]*obj->r0[i][j]);

  r2max = 0.0;

  obj->vlocal = 0;
  obj->bondlocal = 0.0;

  for (; pc; pc = pc->nextlocal) {
   
    if (pc->s.nbonds == 0) continue;

    for (n = 0; n < pc->s.nbonds; n++) {
      assert(pc->bonded[n]);
      if (pc->s.index > pc->bonded[n]->s.index) continue;

      /* Compute force arising on each particle from single bond */

      cs_minimum_distance(obj->cs, pc->s.r, pc->bonded[n]->s.r, r12);
      r2 = r12[X]*r12[X] + r12[Y]*r12[Y] + r12[Z]*r12[Z];
      r2md = sqrt(r2);

      it0=pc->s.inter_type;
      it1=pc->bonded[n]->s.inter_type;

      if (r2 < r2min) r2min = r2;
      if (r2 > r2max) r2max = r2;
      if (r2 > 4*obj->r0[it0][it1]*obj->r0[it0][it1]) pe_fatal(obj->pe, "Broken harmonic bond\n");

      obj->vlocal += 0.5*obj->k[it0][it1]*(r2md-obj->r0[it0][it1])*(r2md-obj->r0[it0][it1]);
      obj->bondlocal += 1.0;
      f = -obj->k[it0][it1]*(r2md-obj->r0[it0][it1])/r2md;

      pc->force[X] -= f*r12[X];
      pc->force[Y] -= f*r12[Y];
      pc->force[Z] -= f*r12[Z];

      pc->bonded[n]->force[X] += f*r12[X];
      pc->bonded[n]->force[Y] += f*r12[Y];
      pc->bonded[n]->force[Z] += f*r12[Z];
    }
  }

  obj->rminlocal = sqrt(r2min);
  obj->rmaxlocal = sqrt(r2max);

  return 0;
}

/*****************************************************************************
 *
 *  bond_harmonic_stats
 *
 *****************************************************************************/

int bond_harmonic_stats(void * self, double * stats) {

  bond_harmonic_t * obj = (bond_harmonic_t *) self;

  assert(obj);
  assert(stats);

  stats[INTERACT_STAT_VLOCAL]    = obj->vlocal;
  stats[INTERACT_STAT_RMINLOCAL] = obj->rminlocal;
  stats[INTERACT_STAT_RMAXLOCAL] = obj->rmaxlocal;

  return 0;
}

/*****************************************************************************
 *
 *  bond_harmonic_single
 *
 *  For a single bond, compute v and |f| given r.
 *
 *****************************************************************************/

int bond_harmonic_single(bond_harmonic_t * obj, double r, double * v, double * f,int it0, int it1) {

  assert(obj);
  assert(r < 2*obj->r0[it0][it1]);

  *v = 0.5*obj->k[it0][it1]*(r-obj->r0[it0][it1])*(r-obj->r0[it0][it1]);
  *f = -obj->k[it0][it1]*(r-obj->r0[it0][it1]);

  return 0;
}
