#!/bin/bash
# Actualiza estado_equilibrio.txt cada 30 min hasta que Ludwig termine.
cd "$(dirname "$(readlink -f "$0")")"
S=/home/danolv/Sim/ludwig/microgel/local/scripts
while :; do
  { date '+%F %T'; python3 $S/mgel_equilibrio.py . ; } > estado_equilibrio.txt.tmp 2>&1 \
    && mv estado_equilibrio.txt.tmp estado_equilibrio.txt
  grep -q "finished normally" run.log 2>/dev/null && break
  sleep 1800
done
python3 $S/mgel_serie.py . > /dev/null 2>&1   # serie.json final
