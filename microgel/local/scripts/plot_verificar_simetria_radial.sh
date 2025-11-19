#!/bin/bash
#------------------------------------------------------------------------------------
# Run multiple plots to compare teoretical Electric field to simulations
#------------------------------------------------------------------------------------

while getopts "d:e:k:l:o:p:q:r:s:f:t:" flag
do
    case "${flag}" in
        d) dim=${OPTARG};;         # system dimension
        e) epsilon=${OPTARG};;     # electrostatic epsilo * kt
        k) kappa=${OPTARG};;       # inverse of Debye length
        l) L=${OPTARG};;           # system size
        o) base_file=${OPTARG};;   # base file name (default psi)
        p) paso=${OPTARG};;        # Delta steps for output
        q) log=${OPTARG};;         # log plot
        r) radial=${OPTARG};;      # radial flag
        s) start=${OPTARG};;       # start position
        f) end=${OPTARG};;         # end position
        t) kt=${OPTARG};;          # thermal energy
    esac
done

# Check mandatory options
if [[ -z "$epsilon" || -z "$kappa" || -z "$paso" ]]; then
        echo 'Missing mandatory input line parameters' >&2
        exit 1
fi

if [[ "$dim" -eq 2 ]]; then
    Lz=4
elif [[ "$dim" -eq 3 ]]; then
    Lz=$L
else
    echo "Error: La dimensión '$dim' no es válida. Debe ser 2 o 3." >&2
    exit 1
fi

for d in ../*/; do
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

    if [[ -z $radial ]]; then
        if [[ $start =~ ([0-9.]+)_([0-9.]+)_([0-9.]+) ]]; then
            sx="${BASH_REMATCH[1]}"
            sy="${BASH_REMATCH[2]}"
            sz="${BASH_REMATCH[3]}"
        else
            echo "   ⚠ No se pudo extraer parámetros de start $s"
            continue
        fi

        if [[ $end =~ ([0-9.]+)_([0-9.]+)_([0-9.]+) ]]; then
            fx="${BASH_REMATCH[1]}"
            fy="${BASH_REMATCH[2]}"
            fz="${BASH_REMATCH[3]}"
        else
            echo "   ⚠ No se pudo extraer parámetros de end$e"
            continue
        fi
    fi

    echo "   Parámetros: p1=$p1  p2=$p2  p3=$p3"

    p1a=$(echo "$p1 - 0.5" | bc)
    p2a=$(echo "$p2 - 0.5" | bc)
    p3a=$(echo "$p3 - 0.5" | bc)

    echo "   Parámetros llamada: p1=$p1  p2=$p2  p3=$p3"
    echo "   L=$L  Lz=$Lz"
    # Ejecutar el script python dentro de la carpeta
    (
        cd "$d" || exit
        ../scripts/verificar_simetria_psi.py -f ./psi-$(printf "%09d" "$paso").001-001 \
                 -c ./config.cds$(printf "%08d" "$paso").001-001 \
                 -s $L $L $Lz --compare-theory --charge 1.0 \
                 --kappa $kappa \
                 --epsilon $epsilon \
                 --kt $kt \
                 --radial-only \
                 --shell-thickness 0.15 \
                 -o ./plots/simetria_DH-pos_${p1}_${p2}_${p3}-kappa_$kappa-eps_$epsilon-p_$paso-

        ../scripts/verificar_simetria_psi.py -f ./psi-$(printf "%09d" "$paso").001-001 \
                 -c ./config.cds$(printf "%08d" "$paso").001-001 \
                 -s $L $L $Lz --compare-theory --charge 1.0 \
                 --kappa $kappa \
                 --epsilon $epsilon \
                 --kt $kt \
                 --radial-only \
                 --log-scale xy \
                 --shell-thickness 0.15 \
                 -o ./plots/simetria_DH-pos_${p1}_${p2}_${p3}-kappa_$kappa-eps_$epsilon-p_$paso-log-
    )
done