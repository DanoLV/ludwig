#!/bin/bash
# Test de linealidad: 4 baterias de 3 o 1 posiciones, 4 corridas en paralelo.
cd "$(dirname "$(readlink -f "$0")")"
for b in C-lineal-hann8-LD4-q0.25 C-lineal-hann8-LD4-q2 C-lineal-hann8-LD8-q0.25 C-lineal-hann8-LD8-q2; do
  bash $b/correr.sh 3 > $b/bateria.log 2>&1
done
echo LISTO > lineal_listo.txt
