#!/bin/bash
# Cola Peskin-6: espera a las corridas chicas (dx = 2, 1) y corre las dos dx = 1/2 de a una (8 GB de GPU no alcanzan
# para ambas junto con otras).
cd "$(dirname "$(readlink -f "$0")")"
for d in C-electrolito-peskin6s{1-LD3-ref0.5,2-LD6-ref1}{,-q0.25}; do
  until grep -aqs BATERIA_LISTA $d/bateria.log; do sleep 30; done
done
for d in C-electrolito-peskin6s4-LD12-ref2 C-electrolito-peskin6s4-LD12-ref2-q0.25; do
  (cd $d && bash correr.sh 1 > bateria.log 2>&1)
done
echo COLA_LISTA
