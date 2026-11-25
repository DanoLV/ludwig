#!/bin/bash
# Corre advfix/* (crear_advfix.py, 2026-11-16) en paralelo dejando >= 2100 MiB libres en la GPU
# (regla del usuario); cada corrida 24^3 usa ~904 MiB (medido 2026-11-10). No borra nada.
export LD_LIBRARY_PATH=/usr/local/petsc-cuda-hypre/lib:/usr/local/ompi/lib:$LD_LIBRARY_PATH
cd /home/danolv/Sim/ludwig/microgel/local/TEST-HIDRO-dx2/advfix
RESERVA=2100; MEM=904
libre() { nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits | head -1; }
una() {
  (cd $1 && stdbuf -o0 ./Ludwig.exe > run.log 2>&1)
  grep -aq "finished normally" $1/run.log && echo "$(date '+%F %T') [OK] $1" || echo "$(date '+%F %T') [FALLO] $1"
}
for d in hidro_si_eta0.02 hidro_si hidro_si_eta0.1 hidro_si_eta0.2 hidro_si_eta0.5 hidro_si_eta1.0 hidro_si_eta2.0; do
  if grep -aqs "finished normally" $d/run.log; then echo "[SKIP] $d"; continue; fi
  while [ "$(libre)" -lt $((RESERVA + MEM)) ]; do sleep 30; done
  echo "$(date '+%F %T') lanzo $d (libre antes $(libre) MiB)"
  una $d &
  sleep 30
done
wait
echo "$(date '+%F %T') ADVFIX_LISTO"
