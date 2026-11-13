#!/bin/bash
# Corre baterias de fuerza entre 2 particulas (dirs con pos_*), N en paralelo. Reanudable:
# saltea las posiciones terminadas y rehace desde cero las que quedaron a medias.
# Uso: run_force_battery.sh NPAR dir1 [dir2 ...]   (las posiciones se intercalan entre dirs)
export LD_LIBRARY_PATH=/usr/local/petsc-cuda-hypre/lib:/usr/local/ompi/lib:$LD_LIBRARY_PATH
export LC_ALL=C
export EXTRACT=/home/danolv/Sim/ludwig/util/extract_colloids   # exportada: xargs la necesita
NPAR=$1; shift
run_one() {
  d="$1"; cd "$d" || return 1
  if [ -f colloid_data/colloids-00003000.csv ] && grep -aqs "finished normally" logs/output.txt; then
    echo "[SKIP] ${d#*local/}"; return 0; fi
  rm -rf proceced_data colloid_data logs
  rm -f config.cds0* vel-0* psi-0* dist-0* efield-0* qsi-0* rho-0* colloids-0*
  mkdir -p proceced_data colloid_data logs
  ./Ludwig.exe > run.log 2>&1; rc=$?; cp run.log logs/output.txt
  grep -aq "out of memory" run.log && { echo "[OOM] ${d#*local/}"; return 1; }
  grep -aq "finished normally" run.log || { echo "[FALLO rc=$rc] ${d#*local/}"; return 1; }
  "$EXTRACT" config.cds00003000.001-001 > /dev/null 2>&1
  mv -f colloids-00003000.csv colloid_data/ 2>/dev/null
  [ -f colloid_data/colloids-00003000.csv ] || { echo "[SIN CSV] ${d#*local/}"; return 1; }
  echo "[OK] ${d#*local/}"
}
export -f run_one
# intercalar: misma distancia de todos los dirs antes de pasar a la siguiente
python3 - "$@" << 'PY' | xargs -P "$NPAR" -I{} bash -c 'run_one "$@"' _ {}
import sys, glob, re, os
groups = [sorted(glob.glob(os.path.join(os.path.abspath(g), "pos_*")),
                 key=lambda p: float(re.search(r"pos_([0-9.]+)_", p).group(1))) for g in sys.argv[1:]]
for row in zip(*groups):
    for p in row: print(p)
PY
echo BATERIA_LISTA
