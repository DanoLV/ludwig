#!/bin/bash
# Corre las 18 corridas d = 8.2 ... 9.8 (dx = 2 primero, son cortas; despues dx = 1), crear_2part_d8a10.py, 2026-11-11.
# Lanza en paralelo las que entren dejando siempre >= RESERVA MiB libres en la GPU (criterio del usuario, 2026-11-10);
# cada corrida usa <= 1014 MiB (48^3, medido). Espera 60 s entre lanzamientos. No borra nada; log en run.log de cada una.
export LD_LIBRARY_PATH=/usr/local/petsc-cuda-hypre/lib:/usr/local/ompi/lib:$LD_LIBRARY_PATH
cd "$(dirname "$(readlink -f "$0")")"
RESERVA=2100; MEM=1014
libre() { nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits | head -1; }
una() {
  (cd $1 && stdbuf -o0 ./Ludwig.exe > run.log 2>&1)
  grep -aq "finished normally" $1/run.log && echo "$(date '+%F %T') [OK] $1" || echo "$(date '+%F %T') [FALLO] $1"
}
for d in SLIDES-LD8-hann04-corr_no-dx2-d8a10/r_* SLIDES-LD8-hann08-corr_no-N10000-d8a10/pos_*; do
  if grep -aqs "finished normally" $d/run.log; then echo "[SKIP] $d"; continue; fi
  if ls $d/psi-0* >/dev/null 2>&1; then echo "[NO TOCO, salida parcial] $d"; continue; fi
  while [ "$(libre)" -lt $((RESERVA + MEM)) ]; do sleep 60; done
  echo "$(date '+%F %T') lanzo $d (libre antes $(libre) MiB)"
  una $d &
  sleep 60
done
wait
echo "$(date '+%F %T') BATERIA_LISTA"
