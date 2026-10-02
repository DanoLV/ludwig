#!/bin/bash
# Plantilla de lanzamiento para corridas de microgel: corre Ludwig.exe y, en
# paralelo, va convirtiendo cada config.cds recien volcado a colloids-*.csv
# en colloid_data/ (via csv_watch.sh / extract_colloids), para poder abrir
# la corrida en Paraview mientras avanza -- igual que hace runbg.sh con
# runplot.sh/coloideacsv.sh para las corridas de una sola particula.
export LD_LIBRARY_PATH=/usr/local/petsc-cuda-hypre/lib:/usr/local/ompi/lib:$LD_LIBRARY_PATH
cd "$(dirname "$(readlink -f "$0")")"
mkdir -p logs colloid_data

CSV_WATCH=/home/danolv/Sim/ludwig/microgel/local/scripts/csv_watch.sh

./Ludwig.exe > run.log 2>&1 &
pid=$!

if [ -x "$CSV_WATCH" ]; then
  bash "$CSV_WATCH" "$(pwd)" --pid "$pid" --poll 60 > logs/csv_watch.log 2>&1 &
fi

wait "$pid"
rc=$?
cp run.log logs/output.txt
echo "exit=$rc terminado=$(grep -ac 'finished normally' run.log)"
