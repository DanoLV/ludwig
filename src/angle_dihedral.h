//CHANGE3
/*****************************************************************************
 *
 *  angle_dihedral.h
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  (c) 2017 The University of Edinburgh
 *
 *  Contributing authors:
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *  Kai Qi (kai.qi@epfl.ch)
 *
 *****************************************************************************/

#ifndef LUDWIG_ANGLE_DIHEDRAL_H
#define LUDWIG_ANGLE_DIHEDRAL_H

#include "pe.h"
#include "coords.h"
#include "colloids.h"
#include "interaction.h"

typedef struct angle_dihedral_s angle_dihedral_t;

int angle_dihedral_create(pe_t * pe, cs_t * cs, angle_dihedral_t ** pobj);
int angle_dihedral_free(angle_dihedral_t * obj);
int angle_dihedral_info(angle_dihedral_t * obj);
int angle_dihedral_param_set(angle_dihedral_t * obj, double kappa[][NUM_INT_PART_TYPES][NUM_INT_PART_TYPES][NUM_INT_PART_TYPES], int mu[][NUM_INT_PART_TYPES][NUM_INT_PART_TYPES][NUM_INT_PART_TYPES], double phi0[][NUM_INT_PART_TYPES][NUM_INT_PART_TYPES][NUM_INT_PART_TYPES]);
int angle_dihedral_register(angle_dihedral_t * obj, interact_t * parent);
int angle_dihedral_compute(colloids_info_t * cinfo, void * self);
int angle_dihedral_stats(void * self, double * stats);

#endif
