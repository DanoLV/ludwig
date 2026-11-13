#!/bin/bash
# Runs the four two-potential electrolyte batteries (lambda_D = 8, 10, 12, 14) in parallel.
cd "$(dirname "$(readlink -f "$0")")"
printf "%s\n" C-electrolito-dospotenciales-hann8-LD{14,12,10,8} | xargs -P 4 -I{} bash {}/correr.sh
echo TODAS_LISTAS > dospot_electrolito_listo.txt
