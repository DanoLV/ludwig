#!/bin/bash
# SAL-SERIE (2026-11-16): corre los casos de crear_serie_sal.py, de a NPAR en paralelo (default 3),
# dejando >= 2100 MiB libres de GPU. NO se corrió al crearlo (pedido del usuario: solo armar la serie).
# Costo estimado: casos de dx = 1 ~18 min c/u (0.035 s/paso x 3e4); DXH-* ~4.7 h c/u (8x nodos, 6e4 pasos).
# Uso: NPAR=3 bash microgel/scripts/INGO-1PART-LBTERM-DIAG/correr_serie_sal.sh [caso ...]   (sin casos: todos)
BASE=/home/danolv/Sim/ludwig/microgel/local/INGO-1PART-LBTERM-DIAG/SAL-SERIE
LIB=/usr/local/petsc-cuda-hypre/lib:/usr/local/ompi/lib
NPAR=${NPAR:-3}; RESERVA=2100
CASOS=${@:-$(ls $BASE)}
for c in $CASOS; do
  d=$BASE/$c
  if grep -q "finished normally" $d/run.log 2>/dev/null; then echo "$c: ya corrida, se salta"; continue; fi
  while [ $(jobs -r | wc -l) -ge $NPAR ] || \
        [ $(nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits | head -1) -lt $((RESERVA + 1000)) ]; do
    sleep 30
  done
  ( cd $d && t0=$(date +%s) && LD_LIBRARY_PATH=$LIB:$LD_LIBRARY_PATH ./Ludwig.exe > run.log 2>&1
    echo "$c: exit $? en $(( $(date +%s) - t0 )) s; $(grep -c 'finished normally' run.log) finished normally" ) &
  sleep 20
done
wait
