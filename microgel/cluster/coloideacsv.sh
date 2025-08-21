#!/bin/bash
#-------------------------------------------------------------------------------------------------------------
# Convert config colloid files by batch
#-------------------------------------------------------------------------------------------------------------
# Parameters
#---------------------------------------------
while getopts "n:i:p:" flag
do
    case "${flag}" in
        n) Ntotal=${OPTARG};;
        i) Ninicio=${OPTARG};;
        p) paso=${OPTARG};;
    esac
done

#---------------------------------------------
# file names  
#---------------------------------------------
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
./extract_colloids ${filen[$i]}
done
