#!/bin/bash
# Resumes every r_* from its step N0 (N in parallel). Deletes nothing; log in run_cont.log.
export LD_LIBRARY_PATH=/usr/local/petsc-cuda-hypre/lib:/usr/local/ompi/lib:$LD_LIBRARY_PATH
cd "$(dirname "$(readlink -f "$0")")"
N=${1:-3}
uno() {
  cd "$1" || return 1
  if grep -aqs "finished normally" run_cont.log; then echo "[SKIP] $1"; return 0; fi
  n0=$(awk '/^N_start/ {print $2}' input)
  stdbuf -o0 ./Ludwig.exe > run_cont.log 2>&1
  grep -aq "finished normally" run_cont.log && grep -aq "Re-starting simulation at step $n0" run_cont.log \
    && echo "[OK] $1" || echo "[FALLO] $1"
}
export -f uno
ls -d r_* | xargs -P "$N" -I{} bash -c 'uno {}'
echo BATERIA_LISTA
