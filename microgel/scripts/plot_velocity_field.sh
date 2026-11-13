#!/bin/bash
#------------------------------------------------------------------------------------
# Run multiple plots to compare teoretical Electric field to simulations
#------------------------------------------------------------------------------------

while getopts "p:x:y:z:" flag
do
    case "${flag}" in
        p) paso=${OPTARG};;        # Delta steps for output
        x) Lx=${OPTARG};;          # grid size X
        y) Ly=${OPTARG};;          # grid size Y
        z) Lz=${OPTARG};;          # grid size Z
    esac
done

for d in ./*/; do
    [ -d "$d" ] || continue

    # Ejecutar el script python dentro de la carpeta
    (
        # cd "$d" || exit
        echo "directorio: $d"
        # ls -l
        ../scripts/plot_velocity_field.py \
                 -d "$d" \
                 -t "$paso" \
                 -s $Lx $Ly $Lz \
                 -o "$d/plots/velocity-p_$paso"
    )
done