#!/bin/bash
#------------------------------------------------------------------------------------
# Run ludwig microgel simulation plots
#------------------------------------------------------------------------------------
# Input parameters:
#   -n  Nsteps       : Number of steps to calculate
#   -i  Ninicio      : Initial step number
#   -p  paso         : Step interval for output
#   -l  ladox        : Frame size in the X direction
#   -y  ladoyz       : Frame size in the Y and Z directions
#   -d  del          : Flag to run script that deletes files
#   -v  viscosidad   : Viscosity value
#   -e  energy       : Free energy model to use: 'fe_electro' or 'none'
#   -f  electric_e0  : Excternal electric field: Ex_Ey_Ez
#   -q  fluid_only   : Plot only fluid average velocity (y/n)
#------------------------------------------------------------------------------------
# Example:
#./run.sh -n 1000 -i 0000 -p 500 -l 26 -y 26 -v 1e-03 -e fe_electro -f 0.001_0.0_0.0 -s y -d y
# ./coloideacsv.sh -n 1100000 -i 1000000 -p 500
# ./calculosvelfluid.py -nciclos 1000000 -npaso 1000 -o datosfluid.csv
# ./calculosvel.py -nciclos 1100000 -npaso 500 -o datos.csv
#------------------------------------------------------------------------------------

# Simulation parameters
while getopts "n:i:p:b:l:d:y:v:e:f:s:q:" flag
do
    case "${flag}" in
        n) Nsteps=${OPTARG};;       # Number of steps to calculate
        i) Ninicio=${OPTARG};;      # Initial step number
        p) paso=${OPTARG};;         # Delta steps for output
        v) viscosidad=${OPTARG};;   # Viscosity
        e) energy=${OPTARG};;       # free_energy: fe_electro/none
        f) electric_e0=${OPTARG};;  # electric_e0
        s) single=${OPTARG};;       # Single monomer
        q) fluid_only=${OPTARG};;   # Plot only fluid average velocity
    esac
done

# Run Plot
echo "Inicia plot:"

echo "Nsteps=$Nsteps"
echo "Ninicio=$Ninicio"
echo "paso=$paso"
echo "single=$single"
echo "fluid_only=$fluid_only"

# Total steps of simulation
NT=$((Nsteps + Ninicio))

# Postprocesing - convert data to .cvs files
./coloideacsv.sh -n $NT -i $Ninicio -p $paso

# Ensure graficos subdirectory exists
mkdir -p graficos

# Check if we need to plot fluid velocities
if [ "$fluid_only" == "y" ]; then
    # Plot only fluid velocities (no colloids needed)
    # Run calculation in root dir (where vel-* files are), save output to graficos
    ./calculosvelfluidonly.py -nciclos $NT -ninicio $Ninicio -npaso $paso -o graficos/datosfluid.csv

    # Generate plot in graficos directory
    cd graficos
    ./plotvel.py -i datosfluid.csv --fluid-only
    cd ..

elif [ "$single" == "y" ]; then
    # Calculates and plots velocity for a single subgrid monomer
    # Run calculation in root dir (where colloids-*.csv are), save output to graficos
    ./calculosvelfluid.py -nciclos $NT -ninicio $Ninicio -npaso $paso -o graficos/datosfluid.csv

    # Generate plot in graficos directory
    cd graficos
    ./plotvel.py -i datosfluid.csv
    cd ..

else
    # Calculates and plots density, medium bond length and inertia moments for a microgel
    # Run calculation in root dir (where colloids-*.csv are), save output to graficos
    ./calculos.py -nciclos $NT -npaso $paso -o graficos/datos.csv

    # Generate plot in graficos directory
    cd graficos
    ./plot.py
    cd ..
fi
