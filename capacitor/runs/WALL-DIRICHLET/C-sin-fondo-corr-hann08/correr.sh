#!/bin/bash
# Corre todas las posiciones z_* (1 paso cada una), N en paralelo. Reanudable.
export LD_LIBRARY_PATH=/usr/local/petsc-cuda-hypre/lib:/usr/local/ompi/lib:$LD_LIBRARY_PATH
cd "$(dirname "$(readlink -f "$0")")"
N=${1:-3}
uno() {
  cd "$1" || return 1
  if grep -aqs "finished normally" run.log && [ -s proceced_data/particle_Esub.csv ]; then echo "[SKIP] $1"; return 0; fi
  rm -f psi-* qsi-* efield-* dist-* rho-* vel-* config.cds0* proceced_data/*.csv
  stdbuf -o0 ./Ludwig.exe > run.log 2>&1
  grep -aq "finished normally" run.log && echo "[OK] $1" || echo "[FALLO] $1"
}
export -f uno
ls -d z_* | xargs -P "$N" -I{} bash -c 'uno {}'
echo BATERIA_LISTA
