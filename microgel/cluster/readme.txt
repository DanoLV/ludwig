Overview
----------------------------------
This Ludwig branch addresses the problem of assembling micropolymers
using subgrid particles and studying their response to electric fields
and charges.
New code to compute electrical forces over subgrid particles was incorporated.

This ludwig version incorporates bond harmonic, angle harmonic and angle dihedral
interactions coded by Kai Qi (kai.qi@epfl.ch)

Contributing authors:
Daniel La Valle

---------------------------------------------------------------------------------------------------------
util/microgel_poly_cross_density.c
----------------------------------
The utility microgel_poly_cross_density generates an input colloid file
to simulate a microgel distributed within a spherical volume, with
crosslinks created between polymer chains according to the specified
crosslinking density.

---------------------------------------------------------------------------------------------------------
Instructions to Run

1. Compile Ludwig and copy the executable to the working directory
   using the script ludwig_redo.sh.

2. Compile extract_colloids from /utils and copy the executable to the working directory

3. Compile microgel_poly_cross_density using the script
   microgel_redo.sh.

4. Run microgel_poly_cross_density using the script microgel_run.sh,
   setting the parameters to the desired values.

5. Run Ludwig using the script run.sh with the appropriate
   command-line parameters.
    Note:
    - Configurations files included for a 50_50_50 size.
    - Rename "config.cds.init.001-001-1 monomer" to config.cds.init.001-001 to run with 1 monomer.
    - Rename "config.cds.init.001-001-microgel" to config.cds.init.001-001 to run with 1 microgel.


Notes on Output
- If the colloid input file contains only one monomer,
  the resulting plot will display the velocity in the x-direction.
- If the colloid input file contains multiple monomers,
  the plot will show density, average bond length, and moments of inertia.