#!/bin/bash
#-------------------------------------------------------------------------------------------------------------
# Convert config colloid files by batch to .csv deleting vtk files
#-------------------------------------------------------------------------------------------------------------
# Parameters
#---------------------------------------------
while getopts "n:i:p:o:" flag
do
    case "${flag}" in
        n) Ntotal=${OPTARG};;
        i) Ninicio=${OPTARG};;
        p) paso=${OPTARG};;
        o) output_dir=${OPTARG};;
    esac
done

#---------------------------------------------
# file names  
#---------------------------------------------
dir=${output_dir:-colloid_data}
file='config.cds'
fileext='.001-001'
filen=()
for ((i=Ninicio;i<=Ntotal;i+=paso))
do
    filen+=($file$(printf "%0${9}d" "$i")$fileext)
done

#---------------------------------------------
# Convert files
#---------------------------------------------
for i in ${!filen[@]};
do
./scripts/extract_colloids ${filen[$i]}
done
