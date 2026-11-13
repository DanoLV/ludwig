#!/bin/bash
#------------------------------------------------------------------------------------
# Run ludwig microgel simulation with MPI support
#------------------------------------------------------------------------------------
# Run './runbg.sh --help' for detailed usage information
#------------------------------------------------------------------------------------

# Simulation parameters
while getopts "n:i:p:b:l:d:y:v:e:f:s:q:" flag
do
    case "${flag}" in
        p) paso=${OPTARG};;         # Delta steps for output
        l) zize=${OPTARG};;         # Size of the grid
        f) electric_e0=${OPTARG};;  # electric_e0
        s) single=${OPTARG};;       # Single monomer

    esac
done

./plot_electric_field_peskin.py -f ../psi-000030000.001-001 -s 32 32 32 -m plane -p xy -c x -o p30000-field_xy_Ex_peskin.png
./plot_electric_field.py -f ../psi-000030000.001-001 -s 32 32 32 -m plane -p xy -c y -o p30000-field_xy_Ey.png
./plot_electric_field.py -f ../psi-000030000.001-001 -s 32 32 32 -m plane -p xy -c z -o p30000-field_xy_Ez.png

./plot_electric_field.py -f ../psi-000030000.001-001 -s 32 32 32 -m plane -p xy -c psi -o p30000-potencial_xy.png
./plot_electric_field.py -f ../psi-000030000.001-001 -s 32 32 32 -m plane3d -p xy -c psi -o p30000-potencial_3d_xy.png

./plot_electric_field_peskin.py -f ../psi-000030000.001-001 -s 32 32 32 -m plane -p xy -c x -o p30000-field_xy_Ex_peskin.png
./plot_electric_field_peskin.py -f ../psi-000030000.001-001 -s 32 32 32 -m plane -p xy -c y -o p30000-field_xy_Ey_peskin.png
./plot_electric_field_peskin.py -f ../psi-000030000.001-001 -s 32 32 32 -m plane -p xy -c z -o p30000-field_xy_Ez_peskin.png