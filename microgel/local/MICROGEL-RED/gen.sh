#!/bin/bash
# Red esferica entrecruzada (microgel): colocacion por rechazo con separacion
# minima 2.3 y enlace de TODO par a distancia [2.5, 3.5], que es lo que da la
# conectividad alta (z=4.72) que el generador de cadena+entrecruzamientos no
# alcanza (z=1.92). Radio 8 y densidad de saturacion para igualar el tamano
# del microgel anterior (Rg 6.67 vs 6.41; diametro 16.0 vs 15.9).
export LD_LIBRARY_PATH=/usr/local/petsc-cuda-hypre/lib:/usr/local/ompi/lib:$LD_LIBRARY_PATH
/home/danolv/Sim/ludwig/util/porous_solid_cube \
   --side 48 --sphere-radius 8.0 \
   --lbond 3.0 --delta 0.5 --min-distance 2.3 --density 0.15 \
   --max_links 6 --seed 20260924 --irad 0.1 --hrad 0.1
