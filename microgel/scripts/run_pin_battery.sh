#!/bin/bash
# Corre baterias de pin-to-grid (dirs con x*), 4 en paralelo. Reanudable.
# Uso: run_pin_battery.sh dir1 [dir2 ...]
export LD_LIBRARY_PATH=/usr/local/petsc-cuda-hypre/lib:/usr/local/ompi/lib:$LD_LIBRARY_PATH
export LC_ALL=C
run_one() { d="$1"; cd "$d" || return 1
  grep -aqs "finished normally" logs/output.txt && { echo "[SKIP] ${d#*local/}"; return 0; }
  rm -rf proceced_data logs; rm -f config.cds0* dist-0* vel-0* psi-0* efield-0* qsi-0* rho-0*
  mkdir -p proceced_data colloid_data logs
  ./Ludwig.exe > run.log 2>&1; rc=$?; cp run.log logs/output.txt
  grep -aq "finished normally" run.log && echo "[OK] ${d#*local/}" || echo "[FALLO rc=$rc] ${d#*local/}"; }
export -f run_one
for g in "$@"; do ls -d "$g"/x*; done | xargs -P 4 -I{} bash -c 'run_one "$@"' _ {}
echo PIN_BATERIA_LISTA
