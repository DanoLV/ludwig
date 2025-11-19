#!/bin/bash
#------------------------------------------------------------------------------------
# Run multiple plots to compare teoretical Electric field to simulations
#------------------------------------------------------------------------------------

while getopts "e:k:o:p:s:f:t:x:y:z:" flag
do
    case "${flag}" in
        e) epsilon=${OPTARG};;     # electrostatic epsilo * kt
        k) kappa=${OPTARG};;       # inverse of Debye length
        o) base_file=${OPTARG};;   # base file name (default psi)
        p) paso=${OPTARG};;        # Delta steps for output
        s) start=${OPTARG};;       # start position
        f) end=${OPTARG};;         # end position
        t) kt=${OPTARG};;          # thermal energy
        x) Lx=${OPTARG};;          # grid size X
        y) Ly=${OPTARG};;          # grid size Y
        z) Lz=${OPTARG};;          # grid size Z
    esac
done

# Check mandatory options
if [[ -z "$epsilon" || -z "$kappa" || -z "$paso" ]]; then
        echo 'Missing mandatory input line parameters' >&2
        exit 1
fi

# Set default grid size if not provided
if [[ -z "$Lx" || -z "$Ly" || -z "$Lz" ]]; then
        echo 'Grid size not provided. Trying to extract from directory name...'
        # Extract from directory name (e.g., Lx_16_Ly_16_Lz_4)
        current_dir=$(basename "$(pwd)")
        echo "Current directory: $current_dir"
        if [[ $current_dir =~ Lx_([0-9]+)_Ly_([0-9]+)_Lz_([0-9]+) ]]; then
                Lx="${BASH_REMATCH[1]}"
                Ly="${BASH_REMATCH[2]}"
                Lz="${BASH_REMATCH[3]}"
                echo "Extracted grid size: $Lx x $Ly x $Lz"
        else
                echo 'Could not extract grid size from directory name. Using default 32 32 32' >&2
                Lx=32
                Ly=32
                Lz=32
        fi
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

    echo "   Parámetros: p1=$p1  p2=$p2  p3=$p3"

    p1a=$(echo "$p1 - 0.5" | bc)
    p2a=$(echo "$p2 - 0.5" | bc)
    p3a=$(echo "$p3 - 0.5" | bc)

    echo "   Parámetros llamada: p1=$p1  p2=$p2  p3=$p3"

    # Ejecutar el script python dentro de la carpeta
    (
        cd "$d" || exit
        ../scripts/compare_field_theory.py \
                 -f ./psi-$(printf "%09d" "$paso").001-001 \
                 --charge-pos $p1a $p2a $p3a --size $Lx $Ly $Lz \
                 --mode line --start $sx $sy $sz --end $fx $fy $fz \
                 --num-points 5000 --show-nodes \
                 --kappa $kappa \
                 --epsilon $epsilon \
                 --kt $kt \
                 --output ./plots/${base_file}-pos_${p1}_${p2}_${p3}-kappa_$kappa-eps_$epsilon-p_$paso.png
    )
done