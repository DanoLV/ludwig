#!/bin/bash
# Cola Hann cos^4: espera a las corridas chicas (dx = 2, 1) y corre las dos dx = 1/2 de a una, en paralelo con la cola
# del Peskin-6 (dos corridas grandes juntas entran en los 8 GB; cinco arrancando juntas no).
cd "$(dirname "$(readlink -f "$0")")"
for d in C-electrolito-hann{4p4-LD3-ref0.5,8p4-LD6-ref1}{,-q0.25}; do
  until grep -aqs BATERIA_LISTA $d/bateria.log; do sleep 30; done
done
for d in C-electrolito-hann16p4-LD12-ref2 C-electrolito-hann16p4-LD12-ref2-q0.25; do
  (cd $d && bash correr.sh 1 > bateria.log 2>&1)
done
echo COLA_LISTA
