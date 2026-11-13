#!/bin/bash
# Corre todas las posiciones d_* (2026-11-09), N en paralelo. Saltea las terminadas (no borra nada).
export LD_LIBRARY_PATH=/usr/local/petsc-cuda-hypre/lib:/usr/local/ompi/lib:$LD_LIBRARY_PATH
cd "$(dirname "$(readlink -f "$0")")"
N=${1:-3}
uno() {
  cd "$1" || return 1
  if grep -aqs "finished normally" run.log; then echo "[SKIP] $1"; return 0; fi
  if ls psi-0* >/dev/null 2>&1; then echo "[NO TOCO, salida parcial] $1"; return 1; fi
  stdbuf -o0 ./Ludwig.exe > run.log 2>&1
  grep -aq "finished normally" run.log && echo "[OK] $1" || echo "[FALLO] $1"
}
export -f uno
ls -d d_* | xargs -P "$N" -I{} bash -c 'uno {}'
echo BATERIA_LISTA
