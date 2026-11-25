#!/bin/bash
# TEST-NOFLUX-ADV (2026-11-16): corre <ver>/{pared,sinpared} en serie (corridas cortas, 16^3).
# Uso: bash capacitor/scripts/TEST-NOFLUX-ADV/correr.sh <ver> [caso ...]
BASE=/home/danolv/Sim/ludwig/capacitor/local/WALL-DIRICHLET/TEST-NOFLUX-ADV
LIB=/usr/local/petsc-cuda-hypre/lib:/usr/local/ompi/lib
CASOS=${@:2}; for caso in ${CASOS:-pared sinpared}; do
  d=$BASE/$1/$caso
  cd "$d" || exit 1
  t0=$(date +%s)
  LD_LIBRARY_PATH=$LIB:$LD_LIBRARY_PATH ./Ludwig.exe > run.log 2>&1
  echo "$1/$caso: exit $? en $(( $(date +%s) - t0 )) s; $(grep -c 'finished normally' run.log) 'finished normally'"
done
