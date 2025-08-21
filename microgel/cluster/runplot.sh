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
#------------------------------------------------------------------------------------
# Example:
#./run.sh -n 1000 -i 0000 -p 500 -l 26 -y 26 -v 1e-03 -e fe_electro -f 0.001_0.0_0.0 -s y -d y
# ./coloideacsv.sh -n 1100000 -i 1000000 -p 500
# ./calculosvel.py -nciclos 1100000 -npaso 500 -o datos.csv
#------------------------------------------------------------------------------------

# Simulation parameters
while getopts "n:i:p:b:l:d:y:v:e:f:s:" flag
do
    case "${flag}" in
        n) Nsteps=${OPTARG};;       # Number of steps to calculate
        i) Ninicio=${OPTARG};;      # Initial step number
        p) paso=${OPTARG};;         # Delta steps for output
        v) viscosidad=${OPTARG};;   # Viscosity
        e) energy=${OPTARG};;       # free_energy: fe_electro/none
        f) electric_e0=${OPTARG};;  # electric_e0 
        s) single=${OPTARG};;       # Single monomer
    esac
done

module purge

# Load modules
module load OpenMPI/5.0.7-GCC-14.2.0
module load Python/3.11.5-GCCcore-13.2.0
source ~/venvs/microgel-env/bin/activate

# Run Plot
echo "Inicia plot:"

echo "Nsteps=$Nsteps"
echo "Ninicio=$Ninicio"
echo "paso=$paso"

# Total steps of simulation
NT=$((Nsteps + Ninicio))

# Postprocesing - convert data to .cvs files
./coloideacsv.sh -n $NT -i $Ninicio -p $paso

if [ "$single" == "y" ]; then
# Calculates and plots velocity for a single subgrid monomer
./calculosvel.py -nciclos $NT -npaso $paso -o datos.csv
./plotvel.py
else 
# Calculates and plots density, medium bond length and inertia moments for a microgel
./calculos.py -nciclos $NT -npaso $paso -o datos.csv
./plot.py
fi

module purge

# Printing date and time
date

echo "Job finished successfully."