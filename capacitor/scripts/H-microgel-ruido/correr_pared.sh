#!/bin/bash
# Prueba de pared (4 corridas de 1 monómero) en paralelo.
export LD_LIBRARY_PATH=/usr/local/petsc-cuda-hypre/lib:/usr/local/ompi/lib:$LD_LIBRARY_PATH
cd "$(dirname "$(readlink -f "$0")")" || exit 1
for d in pared-q0-steeleG-d5 pared-q1-steeleG-d5 pared-q0-steele1e-5-d5 pared-q1-steele1e-5-d5; do
  if grep -aq "finished normally" "$d/logs/output.txt" 2>/dev/null; then echo "[YA] $d"; continue; fi
  ( cd "$d" && stdbuf -o0 ./Ludwig.exe > logs/output.txt 2>&1
    grep -aq "finished normally" logs/output.txt && echo "[OK] $d" || echo "[FALLO] $d" ) &
done
wait
echo "PARED_LISTAS"
