#!/bin/bash
# Prueba del ejecutable 45dc94bc (2026-11-11): hidro_si (debe reproducir hidro_si exacto) y hidro_no, en paralelo.
export LD_LIBRARY_PATH=/usr/local/petsc-cuda-hypre/lib:/usr/local/ompi/lib:$LD_LIBRARY_PATH
cd "$(dirname "$(readlink -f "$0")")"
for d in hidro_si_exe45dc94bc hidro_no_exe45dc94bc; do
  ( cd $d && stdbuf -o0 ./Ludwig.exe > run.log 2>&1; grep -aq "finished normally" run.log && echo "[OK] $d" || echo "[FALLO] $d" ) &
done
wait; echo LISTO
