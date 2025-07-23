#!/bin/bash
#------------------------------------------------------------------------------------
# Run ludwig microgel simulation
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
#./run.sh -n 3000 -i 0000 -p 500 -l 50 -y 50 -v 5e-02 -e fe_electro -f 0.10_0.0_0.0 -s y -d y
#------------------------------------------------------------------------------------

clear

while getopts "n:i:p:l:d:y:v:e:f:s:" flag
do
    case "${flag}" in
        n) Nsteps=${OPTARG};;       # Number of steps to calculate
        i) Ninicio=${OPTARG};;      # Initiual step number
        p) paso=${OPTARG};;         # Delta steps for output
        l) ladox=${OPTARG};;        # Frame size in X direction
        y) ladoyz=${OPTARG};;       # Frame size in Y and Z directions
        d) del=${OPTARG};;          # Run script to delete files
        v) viscosidad=${OPTARG};;   # Viscosity
        e) energy=${OPTARG};;       # free_energy: fe_electro/none
        f) electric_e0=${OPTARG};;  # electric_e0 
        s) single=${OPTARG};;       # Single monomer
    esac
done

# Change parameters in input file
sed -i -e "/N_start/c\N_start $Ninicio" input
sed -i -e "/N_cycles/c\N_cycles $Nsteps" input
sed -i -e "/^viscosity /c\viscosity $viscosidad" input
sed -i -e "/^viscosity_bulk/c\viscosity_bulk $viscosidad" input
sed -i -e "/^free_energy/c\free_energy $energy" input
sed -i -e "/^electric_e0/c\electric_e0 $electric_e0" input
sed -i -e "/colloid_io_freq/c\colloid_io_freq $paso" input
sed -i -e "/size/c\size $ladox\_$ladoyz\_$ladoyz" input

# Total steps of simulation
NT=$((Nsteps + Ninicio))

#Delete files from previus runs
if [ "$del" == "y" ]; then
    ./del.sh
fi

# Run Ludwig
./Ludwig.exe

# Postprocesing - convert data to .cvs files
cp config.cds.init.001-001 config.cds00000000.001-001
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