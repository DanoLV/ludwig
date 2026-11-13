#!/bin/bash
# Corre las 6 baterias C-electrolito-* una tras otra, 3 posiciones en paralelo. Reanudable.
cd "$(dirname "$(readlink -f "$0")")"
for b in C-electrolito-hann8-LD6 C-electrolito-hann6-LD6 C-electrolito-hann8-LD4 C-electrolito-hann6-LD4 C-electrolito-hann8-LD8 C-electrolito-hann6-LD8; do
  bash $b/correr.sh 3 > $b/bateria.log 2>&1
done
echo TODAS_LISTAS > electrolito_listo.txt
