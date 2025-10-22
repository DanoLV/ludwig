#!/bin/bash
#------------------------------------------------------------------------------------
# Compile microgel_poly_cross_density from ../util and copy program to current dir
#------------------------------------------------------------------------------------
(cd ../util/ && make -j10 microgel_poly_cross_density)
cp ../util/microgel_poly_cross_density microgel_poly_cross_density