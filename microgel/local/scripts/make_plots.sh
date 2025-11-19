#!/bin/bash
#------------------------------------------------------------------------------------
# Run multiple plots to compare teoretical Electric field to simulations
#------------------------------------------------------------------------------------
# ./make_plots.sh -d 3 -l 32 -e 10000 -k 0.200761489786388 -o psi_field -p 2000 -t 0.00001
while getopts "d:e:k:l:o:p:r:t:x:y:z:" flag
do
    case "${flag}" in
        d) dim=${OPTARG};;         # system dimension
        e) epsilon=${OPTARG};;     # electrostatic epsilo * kt
        k) kappa=${OPTARG};;       # inverse of Debye length
        l) L=${OPTARG};;           # system size
        o) base_file=${OPTARG};;   # base file name (default psi)
        p) paso=${OPTARG};;        # Delta steps for output
        r) rhoel=${OPTARG};;      # charge density
        t) kt=${OPTARG};;          # thermal energy
        x) Lx=${OPTARG};;          # grid size X
        y) Ly=${OPTARG};;          # grid size Y
        z) Lz=${OPTARG};;          # grid size Z
    esac
done

medio=$(echo "$L / 2 + 0.5" | bc)

Lxm=$(echo "$Lx / 2 + 0.5" | bc)
Lym=$(echo "$Ly / 2 + 0.5" | bc)
Lzm=$(echo "$Lz / 2 + 0.5" | bc)

parent_dir=$(dirname "$PWD")
parent_name=$(basename "$parent_dir")

if [[ -z $kappa ]]; then 
# Extraer parámetros usando regex en bash
echo "paso1"
    if [[ -z $rhoel ]]; then    
    echo "paso2"
        if [[ $parent_name =~ rhoel_([0-9.]+) ]]; then
        echo "paso3"
        rhoel="${BASH_REMATCH[1]}"
        echo "rhoel: $rhoel"
        else
            echo "   ⚠ No se pudo extraer parámetros de $d"
            continue
        fi
    fi

    # # l_B = e^2 * 4 pi epsilon k_BT
    # # pi=$(echo "scale=20; 4*a(1)" | bc -l)
    # l_B=$(echo "4 * (4*a(1)) * $epsiolon * $kt" | bc -l)
    # echo "l_B: $l_B"
    # # l_D = ( 4 pi l_B e sum_k (rho_k Z_k²) )^{-1/2}
    # # charge densities
    # rho_sub=$(echo "1 / ($Lx * $Ly * $Lz)" | bc -l)
    # echo "rho_sub: $rho_sub"
    # sumk_rho=$(echo " 2 * $rhoel + $rho_sub)" | bc -l)
    # echo "sumk_rho: $sumk_rho"

    # l_D=$(echo "sqrt((4 * (4*a(1)) * $l_B * $sumk_rho))" | bc -l)
    # echo "l_D: $l_D"


    # kappa=$(echo "scale=10;1 / $l_D" | bc -l)
    # echo "kappa: $kappa"

    # l_B = e^2 * 4 pi epsilon k_BT
    l_B=$(echo "1/(4 * (4*a(1)) * $epsilon * $kt)" | bc -l)
    echo "l_B: $l_B"

    rho_sub=$(echo "1 / ($Lx * $Ly * $Lz)" | bc -l)
    # echo "rho_sub: $rho_sub"

    sumk_rho=$(echo "2 * $rhoel + $rho_sub" | bc -l)   # Aquí quité espacio al inicio y el paréntesis extra
    echo "sumk_rho: $sumk_rho"

    # l_D = ( 4 pi l_B e sum_k (rho_k Z_k²) )^{-1/2}
    l_D=$(echo "1/(sqrt((4 * (4*a(1)) * $l_B * $sumk_rho)))" | bc -l)
    echo "l_D: $l_D"

    # kappa = 1/l_D 
    kappa=$(echo "sqrt((4 * (4*a(1)) * $l_B * $sumk_rho))" | bc -l)
    echo "kappa: $kappa"

fi

echo "Calculated kappa: $kappa"

# campo electrico sobre la particula en diferentes posiciones segun carpeta
# ./analizar_campo_electrico.py
# Linea sobre x medio
./plot_field_comparizon.sh -x $Lx -y $Ly -z $Lz -e $epsilon -k $kappa -o ${base_file}_x -p $paso -s 0.0_${Lym}_${Lzm} -f ${L}_${Lym}_${Lzm} -t $kt
# Linea sobre y medio
./plot_field_comparizon.sh -x $Lx -y $Ly -z $Lz -e $epsilon -k $kappa -o ${base_file}_y -p $paso -s ${Lxm}_0.0_${Lzm} -f ${Lxm}_${Ly}_${Lzm} -t $kt
# Linea sobre diagonal en plano xy
./plot_field_comparizon.sh -x $Lx -y $Ly -z $Lz -e $epsilon -k $kappa -o ${base_file}_xy_diag -p $paso -s 0.0_0.0_${Lzm} -f ${Lx}_${Ly}_${Lzm} -t $kt
# Linea sobre diagonal del cubo
if [[ "$dim" -eq 3 ]]; then
    ./plot_field_comparizon.sh -x $Lx -y $Ly -z $Lz -e $epsilon -k $kappa -o ${base_file}_xyz_diag -p $paso -s 0.0_0.0_0.0 -f ${Lx}_${Ly}_${Lz} -t $kt
fi

# Simetria radial del modulo del campo electrico
./plot_verificar_simetria_radial.sh -e $epsilon -k $kappa -o ${base_file}_radial -p $paso -t $kt -d $dim -l $L -r y