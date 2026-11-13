#!/bin/bash
# Barrido de viscosidad (2026-11-10). Lanza en PARALELO las viscosidades que entren en la GPU dejando siempre libre
# RESERVA MiB (pedido del usuario: "dejar una memoria libre como la que hay ahora"; medido ~2100-2300 MiB con las tres
# continuaciones de 96^3 corriendo). Cada corrida 24^3 usa ~904 MiB (medido): lanza una si libre >= RESERVA + 904,
# luego espera 60 s a que asigne su memoria antes de volver a medir. No borra nada; log de cada una en run.log.
export LD_LIBRARY_PATH=/usr/local/petsc-cuda-hypre/lib:/usr/local/ompi/lib:$LD_LIBRARY_PATH
cd "$(dirname "$(readlink -f "$0")")"
RESERVA=2100; MEM=904
libre() { nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits | head -1; }
una() {
  (cd $1 && stdbuf -o0 ./Ludwig.exe > run.log 2>&1)
  grep -aq "finished normally" $1/run.log && echo "$(date '+%F %T') [OK] $1" || echo "$(date '+%F %T') [FALLO] $1"
}
for d in hidro_si_eta0.02 hidro_si_eta0.1 hidro_si_eta0.2 hidro_si_eta1.0 hidro_si_eta2.0; do
  if grep -aqs "finished normally" $d/run.log; then echo "[SKIP] $d"; continue; fi
  while [ "$(libre)" -lt $((RESERVA + MEM)) ]; do sleep 60; done
  echo "$(date '+%F %T') lanzo $d (libre antes $(libre) MiB)"
  una $d &
  sleep 60
done
wait
echo "$(date '+%F %T') BARRIDO_LISTO"
