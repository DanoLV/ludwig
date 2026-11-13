#!/bin/bash
# Reproduce todas las pruebas de las paredes a potencial constante y regenera
# los graficos. Usa src/Ludwig.exe (enlazado en cada carpeta) y el .petscrc de
# GPU (CG + BoomerAMG). Duracion: A, segundos; regresion, ~3 min;
# B, ~10 min con la GPU libre; C, ~5 min por kernel.
#
# Uso: ./correr_tests.sh            (todo)
#      ./correr_tests.sh A B        (solo esas)
set -u
export LD_LIBRARY_PATH=/usr/local/petsc-cuda-hypre/lib:/usr/local/ompi/lib:${LD_LIBRARY_PATH:-}
cd "$(dirname "$(readlink -f "$0")")"
PY=/home/danolv/Sim/ludwig/.venv/bin/python3
tests=${*:-"A B C regresion"}

corre() {   # carpeta
  (cd "$1" && rm -f psi-* qsi-* efield-* dist-* rho-* vel-* config.cds0* \
   && stdbuf -o0 ./Ludwig.exe > run.log 2>&1)
}

for t in $tests; do
  case $t in
    A) echo "== A: capacitor sin iones"; corre A-sin-iones
       grep -E "Wall charge|Fluid charge|finished normally" A-sin-iones/run.log ;;
    B) echo "== B: Debye-Hückel (60000 pasos)"; corre B-debye-huckel
       $PY scripts/analizar_B.py B-debye-huckel 60000 ;;
    C) echo "== C: particula subgrid entre paredes a tierra (5 baterias: kernel, fondo, correccion)"
       for spec in carga-imagen-hann08:8 carga-imagen-hann06:6 sin-fondo-hann08:8 corr-hann08:8 sin-fondo-corr-hann08:8; do
         b=C-${spec%%:*}; n=${spec##*:}
         (cd $b && rm -f z_*/run.log && ./correr.sh 3 | tail -1)
         $PY scripts/analizar_C.py $b $n | tail -4
       done
       $PY scripts/comparar_C.py ;;
    regresion) echo "== regresion periodica (debe dar diferencia 0)"
       corre regresion-periodica
       $PY - <<'PYEOF'
import numpy as np
S='/home/danolv/Sim/ludwig/microgel/local/EFIELD-1PART-PETSC-LD8-hann08/corr_yes/'
for f in ('psi-000000001.001-001','qsi-000000001.001-001','efield-000000001.001-001'):
    a=np.loadtxt('regresion-periodica/'+f); b=np.loadtxt(S+f)
    print('  %-26s max|nuevo-original| = %.3e' % (f, abs(a-b).max()))
PYEOF
       ;;
  esac
done

echo "== pruebas numericas y graficos"
$PY scripts/verificar_carga_pared.py | tail -4
$PY scripts/convergencia_campo.py
$PY scripts/graficos_paredes.py
