#!/bin/bash
#------------------------------------------------------------------------------------
# Compile ludwig and copy program to current dir
#------------------------------------------------------------------------------------

# Printing date and time
date

# Load modules
module load OpenMPI/5.0.7-GCC-14.2.0  #OpenMPI/4.1.6-GCC-13.2.0  # OpenMPI/4.1.5-GCC-12.3.0    # OpenMPI/5.0.7-GCC-14.2.0

echo "Running on " `hostname`

# (cd ~/Sim/ludwig/ && make clean)
# (cd ~/Sim/ludwig/ && make -j 10)
# cp ../src/Ludwig.exe Ludwig.exe

(cd .. && make clean)
(cd .. && make -j 10)
cp ../src/Ludwig.exe Ludwig.exe
cp ../util/extract_colloids extract_colloids

module purge

echo "Build finished successfully."