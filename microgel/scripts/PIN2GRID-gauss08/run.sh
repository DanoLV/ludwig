#!/bin/bash
# Pin-to-grid con kernel gaussiano (support 8, sigma 1.44612), 4 en paralelo. Reanudable:
# saltea las posiciones terminadas y rehace desde cero las que quedaron a medias.
export LD_LIBRARY_PATH=/usr/local/petsc-cuda-hypre/lib:/usr/local/ompi/lib:$LD_LIBRARY_PATH
export LC_ALL=C
run_one() { d="$1"; cd "$d" || return 1
  grep -aqs "finished normally" logs/output.txt && { echo "[SKIP] ${d##*/}"; return 0; }
  rm -rf proceced_data logs; rm -f config.cds0* dist-0* vel-0* psi-0* efield-0* qsi-0* rho-0*
  mkdir -p proceced_data colloid_data logs
  ./Ludwig.exe > run.log 2>&1; rc=$?; cp run.log logs/output.txt
  grep -aq "finished normally" run.log && echo "[OK] ${d##*/}" || echo "[FALLO rc=$rc] ${d##*/}"; }
export -f run_one
ls -d "$(dirname "$(readlink -f "$0")")"/x* | xargs -P 4 -I{} bash -c 'run_one "$@"' _ {}
echo PIN_GAUSS_LISTO
