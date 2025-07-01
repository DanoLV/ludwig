#!/bin/bash
#------------------------------------------------------------------------------------
# Run microgel_poly_cross_density 
#------------------------------------------------------------------------------------

grid_size=50
number_monomers_surface=300
bond_length=0.8
distance_threshold=0.65
input_radius=0.05
hydrodynamic_radius=0.05
subgrid_particle_offset=0.5
density=8.0
charge=1.0
permittivity=1.0e3
crosslink_density=1.0
# sa=1.0 #1.0e-99
# saf=1.0 #1.0e-99
# drmax=1.0

./microgel_poly_cross_density \
                --gsize $grid_size \
                --nmon $number_monomers_surface \
                --lbond $bond_length \
                --distt $distance_threshold \
                --irad $input_radius \
                --hrad $hydrodynamic_radius \
                --offset $subgrid_particle_offset \
                --idensity $density \
                --charge $charge \
                --crosslink_density $crosslink_density
                # --permittivity $permittivity \
                # --sa $sa \
                # --saf $saf \
                # --drmax $drmax