#!/bin/bash
# Corrida larga H-microgel-ruido: neutro-largo y superficie-largo en paralelo (reanudable).
export LD_LIBRARY_PATH=/usr/local/petsc-cuda-hypre/lib:/usr/local/ompi/lib:$LD_LIBRARY_PATH
cd "$(dirname "$(readlink -f "$0")")" || exit 1
for d in neutro-largo superficie-largo; do
  if grep -aq "finished normally" "$d/logs/output.txt" 2>/dev/null; then echo "[YA] $d"; continue; fi
  ( cd "$d" && stdbuf -o0 ./Ludwig.exe > logs/output.txt 2>&1
    grep -aq "finished normally" logs/output.txt && echo "[OK] $d" || echo "[FALLO] $d" ) &
done
wait
echo "LARGO_LISTAS"
