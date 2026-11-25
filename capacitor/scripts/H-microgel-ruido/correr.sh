#!/bin/bash
# H-microgel-ruido: corre las 4 pruebas cortas en paralelo (reanudable: salta las terminadas).
# Uso: capacitor/local/WALL-DIRICHLET/H-microgel-ruido/correr.sh
export LD_LIBRARY_PATH=/usr/local/petsc-cuda-hypre/lib:/usr/local/ompi/lib:$LD_LIBRARY_PATH
cd "$(dirname "$(readlink -f "$0")")" || exit 1
for d in neutro-steeleG neutro-steele1e-5 superficie-steeleG superficie-steele1e-5; do
  if grep -aq "finished normally" "$d/logs/output.txt" 2>/dev/null; then echo "[YA] $d"; continue; fi
  ( cd "$d" && stdbuf -o0 ./Ludwig.exe > logs/output.txt 2>&1
    grep -aq "finished normally" logs/output.txt && echo "[OK] $d" || echo "[FALLO] $d" ) &
done
wait
echo "H_LISTAS"
