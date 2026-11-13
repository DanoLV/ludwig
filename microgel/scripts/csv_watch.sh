#!/bin/bash
# Convierte cada config.cds recien volcado a un colloids-*.csv para Paraview,
# a medida que la simulacion corre -- el mismo paso que runbg.sh hace via
# runplot.sh/coloideacsv.sh, pero sin el resto de ese flujo (que asume un
# input armado por runbg.sh). Pensado para correr en paralelo a un run.sh
# que ya lanzo Ludwig.exe por su cuenta.
#
# Uso:
#   csv_watch.sh RUNDIR [--pid PID] [--poll SEGUNDOS]
#
# Sin --pid corre hasta que no aparezcan volcados nuevos durante --poll*3
# segundos seguidos (util para "engancharse" a una corrida ya lanzada, o
# para pasar sobre una que ya termino y dejo volcados sin convertir).

set -u

EXTRACT=/home/danolv/Sim/ludwig/util/extract_colloids

rundir=""
pid=""
poll=30

while [ $# -gt 0 ]; do
  case "$1" in
    --pid) pid="$2"; shift 2;;
    --poll) poll="$2"; shift 2;;
    *) rundir="$1"; shift;;
  esac
done

if [ -z "$rundir" ]; then
  echo "Uso: $0 RUNDIR [--pid PID] [--poll SEGUNDOS]" >&2
  exit 1
fi
if [ ! -x "$EXTRACT" ]; then
  echo "No se encuentra extract_colloids en $EXTRACT (make -C util extract_colloids)" >&2
  exit 1
fi

mkdir -p "$rundir/colloid_data"
cd "$rundir" || exit 1

idle=0
while :; do
  did_work=0
  for f in config.cds[0-9]*.001-001; do
    [ -e "$f" ] || continue
    step=$(echo "$f" | sed -E 's/^config\.cds0*([0-9]+)\.001-001$/\1/')
    [ -z "$step" ] && continue
    step8=$(printf "%08d" "$step")
    csv="colloid_data/colloids-${step8}.csv"
    [ -e "$csv" ] && continue
    "$EXTRACT" "$f" >/dev/null 2>&1
    if [ -e "colloids-${step8}.csv" ]; then
      mv -f "colloids-${step8}.csv" "$csv"
      did_work=1
    fi
    rm -f coll-"${step8}".vtk coll"${step8}".vtk 2>/dev/null
  done
  rm -f coll*.vtk 2>/dev/null

  if [ -n "$pid" ]; then
    kill -0 "$pid" 2>/dev/null || { echo "[csv_watch] PID $pid termino, ultima pasada."; sleep 2; continue_after_death=1; }
    if [ "${continue_after_death:-0}" = "1" ]; then
      # una pasada mas por si el ultimo volcado quedo sin convertir, y salir
      for f in config.cds[0-9]*.001-001; do
        [ -e "$f" ] || continue
        step=$(echo "$f" | sed -E 's/^config\.cds0*([0-9]+)\.001-001$/\1/')
        step8=$(printf "%08d" "$step")
        csv="colloid_data/colloids-${step8}.csv"
        [ -e "$csv" ] && continue
        "$EXTRACT" "$f" >/dev/null 2>&1
        [ -e "colloids-${step8}.csv" ] && mv -f "colloids-${step8}.csv" "$csv"
      done
      rm -f coll*.vtk 2>/dev/null
      echo "[csv_watch] terminado."
      exit 0
    fi
  else
    if [ "$did_work" = "0" ]; then
      idle=$((idle + 1))
      [ "$idle" -ge 3 ] && { echo "[csv_watch] sin novedades, salgo."; exit 0; }
    else
      idle=0
    fi
  fi

  sleep "$poll"
done
