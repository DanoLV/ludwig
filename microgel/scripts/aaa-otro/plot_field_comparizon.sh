#!/bin/bash
#------------------------------------------------------------------------------------
# Run multiple plots to compare teoretical Electric field to simulations
#------------------------------------------------------------------------------------

while getopts "e:k:p:t:" flag
do
    case "${flag}" in
        e) epsilon=${OPTARG};;     # electrostatic epsilo * kt
        k) kappa=${OPTARG};;       # inverse of Debye length
        p) paso=${OPTARG};;        # Delta steps for output
        t) kt=${OPTARG};;          # thermal energy
        
    esac
done

# Check mandatory options
if [[ -z "$epsilon" || -z "$kappa" || -z "$paso" ]]; then
        echo 'Missing mandatory input line parameters' >&2
        exit 1
fi

for d in */; do
    [ -d "$d" ] || continue

    echo "Procesando carpeta: $d"

    # Extraer parámetros usando regex en bash
    if [[ $d =~ pos_([0-9.]+)_([0-9.]+)_([0-9.]+) ]]; then
        p1="${BASH_REMATCH[1]}"
        p2="${BASH_REMATCH[2]}"
        p3="${BASH_REMATCH[3]}"
    else
        echo "   ⚠ No se pudo extraer parámetros de $d"
        continue
    fi

    echo "   Parámetros: p1=$p1  p2=$p2  p3=$p3"

    p1a=$(echo "$p1 - 0.5" | bc)
    p2a=$(echo "$p2 - 0.5" | bc)
    p3a=$(echo "$p3 - 0.5" | bc)

    echo "   Parámetros llamada: p1=$p1  p2=$p2  p3=$p3"

    # Ejecutar el script python dentro de la carpeta
    (
        cd "$d" || exit
        ./graficos/compare_field_theory.py \
                 -f ./psi-$(printf "%09d" "$paso").001-001 \
                 --charge-pos $p1a $p2a $p3a --size 32 32 32 \
                 --mode line --start 0 $p2a $p3a --end 32 $p2a $p3a \
                 --num-points 5000 --show-nodes \
                 --kappa $kappa \
                 --epsilon $epsilon \
                 --kt $kt \
                 --output ./graficos/Campo_eje_x-pos_${p1}_${p2}_${p3}-kappa_$kappa-eps_$epsilon-p_$paso.png
    )
done

# ./compare_field_theory.py -f ../psi-000015000.001-001 --charge-pos 16.10 16 16 --size 32 32 32 --mode line --start 0.5 16 16 --end 31.5 16 16 --num-points 5000 --show-nodes --kappa 0.200761489786388 --epsilon 0.1 --output Campo_eje_x.png