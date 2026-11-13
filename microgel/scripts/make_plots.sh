#!/bin/bash
#------------------------------------------------------------------------------------
# Run multiple plots to compare teoretical Electric field to simulations
#------------------------------------------------------------------------------------
# ./make_plots.sh -d 3 -l 32 -e 10000 -k 0.200761489786388 -o psi_field -p 2000 -t 0.00001
while getopts "d:e:l:o:p:t:x:y:z:" flag
do
    case "${flag}" in
        d) dim=${OPTARG};;         # system dimension
        e) epsilon=${OPTARG};;     # electrostatic epsilo * kt
        # k) kappa=${OPTARG};;       # inverse of Debye length
        l) L=${OPTARG};;           # system size
        o) base_file=${OPTARG};;   # base file name (default psi)
        p) paso=${OPTARG};;        # Delta steps for output
        # r) rhoel=${OPTARG};;      # charge density
        t) kt=${OPTARG};;          # thermal energy
        x) Lx=${OPTARG};;          # grid size X
        y) Ly=${OPTARG};;          # grid size Y
        z) Lz=${OPTARG};;          # grid size Z
    esac
done

shell_thick=0.2
parent_dir=$(dirname "$PWD")
parent_name=$(basename "$parent_dir")

if [[ -z $base_file ]]; then
    base_file="single-peskin-Efield"
fi

if [[ -z $Lx || -z $Ly || -z $Lz ]]; then

    if [[ $parent_name =~ Lx_([0-9]+)_Ly_([0-9]+)_Lz_([0-9]+) ]]; then
        Lx="${BASH_REMATCH[1]}"
        Ly="${BASH_REMATCH[2]}"
        Lz="${BASH_REMATCH[3]}"
    else
        echo "   ⚠ No se pudo extraer parámetros de L de $d"
        [translate:exit] 1
    fi

fi

if [[ -z $dim ]]; then

    # Si Lx, Ly y Lz son enteros
    if [[ $Lx -eq $Ly && $Lx -eq $Lz ]]; then
        dim=3
    elif [[ $Lx -eq $Ly || $Lx -eq $Lz || $Ly -eq $Lz ]]; then
        dim=2
    else
        echo "   ⚠ Dimension invalida. Debe ser 2 o 3."
        exit 1
    fi

fi

lista=($Lx $Ly $Lz)
L="${lista[0]}"
for val in "${lista[@]}"; do
    if (( $(echo "$val > $L" | bc -l) )); then
        L="$val"
    fi
done

medio=$(echo "$L / 2 + 0.5" | bc)

Lxm=$(echo "$Lx / 2 + 0.5" | bc)
Lym=$(echo "$Ly / 2 + 0.5" | bc)
Lzm=$(echo "$Lz / 2 + 0.5" | bc)


for d in ../*/; do

    parent_name=$(basename "$d")
    # echo  "parent_name = $parent_name"

    if [[ "$parent_name" == "scripts" ]]; then
        continue
    fi

    # Extraer parámetros usando regex en bash
    if [[ $parent_name =~ rhoel_([0-9]+([.][0-9]+)?([eE][+-]?[0-9]+)?) ]]; then
        rhoel="${BASH_REMATCH[1]}"
        # echo "rhoel: $rhoel"
        rhoel=$(printf "%f" "$rhoel")
    else
        echo "   ⚠ No se pudo extraer parámetros de $d"
        continue
    fi

     # Extraer parámetros usando regex en bash
    if [[ $parent_name =~ kT_([0-9]+([.][0-9]+)?([eE][+-]?[0-9]+)?) ]]; then
        kt="${BASH_REMATCH[1]}"
        # echo "kt: $kt"
        kt=$(printf "%f" "$kt")
    else
        echo "   ⚠ No se pudo extraer parámetros kt de $d"
        continue
    fi

     # Extraer parámetros usando regex en bash
    if [[ $parent_name =~ eps_([0-9]+([.][0-9]+)?([eE][+-]?[0-9]+)?) ]]; then
        epsilon="${BASH_REMATCH[1]}"
        # echo "epsilon: $epsilon"
        epsilon=$(printf "%f" "$epsilon")
    else
        echo "   ⚠ No se pudo extraer parámetros epsilon de $d"
        continue
    fi


    # l_B = e^2 * 4 pi epsilon k_BT
    l_B=$(echo "1/(4 * (4*a(1)) * $epsilon * $kt)" | bc -l)
    # echo "l_B: $l_B"

    rho_sub=$(echo "1 / ($Lx * $Ly * $Lz)" | bc -l)
    # echo "rho_sub: $rho_sub"

    sumk_rho=$(echo "2 * $rhoel + 2*$rho_sub" | bc -l)   # Aquí quité espacio al inicio y el paréntesis extra
    # echo "sumk_rho: $sumk_rho"

    # l_D = ( 4 pi l_B e sum_k (rho_k Z_k²) )^{-1/2}
    l_D=$(echo "1/(sqrt((4 * (4*a(1)) * $l_B * $sumk_rho)))" | bc -l)
    # echo "l_D: $l_D"

    # kappa = 1/l_D 
    kappa=$(echo "sqrt((4 * (4*a(1)) * $l_B * $sumk_rho))" | bc -l)
    # echo "kappa: $kappa"

    echo "Calculated kappa: $kappa"

    if (( $(echo "$kappa < 0" | bc -l) )); then 
        echo "   ⚠ Valor de kappa menor que 0"
        continue
    fi

    cd $d
    # campo electrico sobre la particula en diferentes posiciones segun carpeta
    ../scripts/analizar_campo_electrico.py $d --no-show
    # Linea sobre x medio
    ../scripts/plot_field_comparizon.sh -x "$Lx" -y "$Ly" -z "$Lz" -e "$epsilon" -k "$kappa" -o "${base_file}_x" -p "$paso" -s "0.0_${Lym}_${Lzm}" -f "${Lx}_${Lym}_${Lzm}" -t "$kt"
    # Linea sobre y medio
    ../scripts/plot_field_comparizon.sh -x "$Lx" -y "$Ly" -z "$Lz" -e "$epsilon" -k "$kappa"  -o "${base_file}_y" -p "$paso" -s "${Lxm}_0.0_${Lzm}" -f "${Lxm}_${Ly}_${Lzm}" -t "$kt"
     # Linea sobre diagonal en plano xy
    ../scripts/plot_field_comparizon.sh -x "$Lx" -y "$Ly" -z "$Lz" -e "$epsilon" -k "$kappa" -o "${base_file}_xy_diag" -p "$paso" -s "0.0_0.0_${Lzm}" -f "${Lx}_${Ly}_${Lzm}" -t "$kt"
    # Linea sobre diagonal del cubo
    if [[ "$dim" -eq 3 ]]; then
        ../scripts/plot_field_comparizon.sh -x "$Lx" -y "$Ly" -z "$Lz" -e "$epsilon" -k "$kappa" -o "${base_file}_xyz_diag" -p "$paso" -s "0.0_0.0_0.0" -f "${Lx}_${Ly}_${Lz}" -t "$kt"
    fi

    # Simetria radial del modulo del campo electrico
    ../scripts/plot_verificar_simetria_radial.sh -a "$shell_thick" -e "$epsilon" -k "$kappa" -o "${base_file}_radial" -p "$paso" -t "$kt" -d "$dim" -l "$L" -r y

    # Campo de velocidades
    ../scripts/plot_velocity_field.sh -p "$paso" -x "$Lx" -y "$Ly" -z "$Lz"

done

../scripts/plot_campo_electrico_rhoel.py log 1e-9
../scripts/plot_campo_electrico_rhoel.py linear