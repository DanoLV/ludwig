#!/bin/bash
# Corre las continuaciones cont<NEND>/<variante>/<corrida> (crear_cont_visc.py, 2026-11-16) en paralelo dejando
# >= 2100 MiB libres en la GPU (regla del usuario); ~904 MiB por corrida 24^3. No borra nada.
# [OK] exige "finished normally" y "Re-starting simulation at step N0".
# Uso: scripts/correr_cont_visc.sh NEND   (corre todas las carpetas de cont<NEND> que no hayan terminado)
export LD_LIBRARY_PATH=/usr/local/petsc-cuda-hypre/lib:/usr/local/ompi/lib:$LD_LIBRARY_PATH
cd /home/danolv/Sim/ludwig/microgel/local/TEST-HIDRO-dx2/cont$1 || exit 1
RESERVA=2100; MEM=904
libre() { nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits | head -1; }
una() {
  n0=$(awk '/^N_start/ {print $2}' $1/input)
  (cd $1 && stdbuf -o0 ./Ludwig.exe > run.log 2>&1)
  grep -aq "finished normally" $1/run.log && grep -aq "Re-starting simulation at step $n0" $1/run.log \
    && echo "$(date '+%F %T') [OK] $1" || echo "$(date '+%F %T') [FALLO] $1"
}
for d in */hidro_si*; do
  if grep -aqs "finished normally" $d/run.log; then echo "[SKIP] $d"; continue; fi
  while [ "$(libre)" -lt $((RESERVA + MEM)) ]; do sleep 20; done
  echo "$(date '+%F %T') lanzo $d (libre antes $(libre) MiB)"
  una $d &
  sleep 30
done
wait
echo "$(date '+%F %T') CONT_LISTO"
