#!/bin/bash
#------------------------------------------------------------------------------------
# Compile ludwig and copy program to current dir
#------------------------------------------------------------------------------------

# (cd ~/Sim/ludwig/ && make clean)
# (cd ~/Sim/ludwig/ && make -j 10)
# cp ../src/Ludwig.exe Ludwig.exe

# (cd ../.. && make clean)
# (cd ../.. && make -j 10)
# cp ../../src/Ludwig.exe Ludwig.exe
# # cp ../../util/extract_colloids extract_colloids
# # cp ../../util/microgel_poly_cross_density microgel_poly_cross_density

export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

(cd ../../src && make clean)
(cd ../../src && make -j 10)
cp ../../src/Ludwig.exe Ludwig.exe