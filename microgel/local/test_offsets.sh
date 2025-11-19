#!/bin/bash
# Script para probar diferentes offsets y encontrar el correcto

PSI_FILE="../psi-000030000.001-001"
CONFIG_POS_X=16.20
CONFIG_POS_Y=16.50
CONFIG_POS_Z=16.50
SIZE="32 32 32"

echo "=========================================="
echo "PROBANDO DIFERENTES OFFSETS"
echo "=========================================="
echo ""
echo "Posición en config.init: $CONFIG_POS_X, $CONFIG_POS_Y, $CONFIG_POS_Z"
echo ""

# Calcular posiciones con diferentes offsets
NO_OFFSET_X=$CONFIG_POS_X
NO_OFFSET_Y=$CONFIG_POS_Y
NO_OFFSET_Z=$CONFIG_POS_Z

OFFSET_05_X=$(python3 -c "print($CONFIG_POS_X - 0.5)")
OFFSET_05_Y=$(python3 -c "print($CONFIG_POS_Y - 0.5)")
OFFSET_05_Z=$(python3 -c "print($CONFIG_POS_Z - 0.5)")

OFFSET_10_X=$(python3 -c "print($CONFIG_POS_X - 1.0)")
OFFSET_10_Y=$(python3 -c "print($CONFIG_POS_Y - 1.0)")
OFFSET_10_Z=$(python3 -c "print($CONFIG_POS_Z - 1.0)")

cd /home/bater/Sim/ludwig/microgel/local/linea/single-r_0.8-v_0.02-kT_0.0005-ah_1.605000000000000e-02-lb_fluctuation_0-fe_electro-q_1.0-e_0.0_0.0_0.0-p200-N_50k_16.00_16.50_16.50/graficos

echo "----------------------------------------"
echo "TEST 1: SIN OFFSET (usar directamente)"
echo "----------------------------------------"
echo "Comando:"
echo "./compare_field_theory.py -f $PSI_FILE \\"
echo "    --charge-pos $NO_OFFSET_X $NO_OFFSET_Y $NO_OFFSET_Z \\"
echo "    --size $SIZE \\"
echo "    --mode line \\"
echo "    --start 0.5 $NO_OFFSET_Y $NO_OFFSET_Z \\"
echo "    --end 31.5 $NO_OFFSET_Y $NO_OFFSET_Z \\"
echo "    --num-points 100 --output campo_no_offset.png"
echo ""

./compare_field_theory.py -f $PSI_FILE \
    --charge-pos $NO_OFFSET_X $NO_OFFSET_Y $NO_OFFSET_Z \
    --size $SIZE \
    --mode line \
    --start 0.5 $NO_OFFSET_Y $NO_OFFSET_Z \
    --end 31.5 $NO_OFFSET_Y $NO_OFFSET_Z \
    --num-points 100 --output campo_no_offset.png 2>&1 | tail -5

echo ""
echo "----------------------------------------"
echo "TEST 2: OFFSET -0.5"
echo "----------------------------------------"
echo "Comando:"
echo "./compare_field_theory.py -f $PSI_FILE \\"
echo "    --charge-pos $OFFSET_05_X $OFFSET_05_Y $OFFSET_05_Z \\"
echo "    --size $SIZE \\"
echo "    --mode line \\"
echo "    --start 0.5 $OFFSET_05_Y $OFFSET_05_Z \\"
echo "    --end 31.5 $OFFSET_05_Y $OFFSET_05_Z \\"
echo "    --num-points 100 --output campo_offset_05.png"
echo ""

./compare_field_theory.py -f $PSI_FILE \
    --charge-pos $OFFSET_05_X $OFFSET_05_Y $OFFSET_05_Z \
    --size $SIZE \
    --mode line \
    --start 0.5 $OFFSET_05_Y $OFFSET_05_Z \
    --end 31.5 $OFFSET_05_Y $OFFSET_05_Z \
    --num-points 100 --output campo_offset_05.png 2>&1 | tail -5

echo ""
echo "----------------------------------------"
echo "TEST 3: OFFSET -1.0"
echo "----------------------------------------"
echo "Comando:"
echo "./compare_field_theory.py -f $PSI_FILE \\"
echo "    --charge-pos $OFFSET_10_X $OFFSET_10_Y $OFFSET_10_Z \\"
echo "    --size $SIZE \\"
echo "    --mode line \\"
echo "    --start 0.5 $OFFSET_10_Y $OFFSET_10_Z \\"
echo "    --end 31.5 $OFFSET_10_Y $OFFSET_10_Z \\"
echo "    --num-points 100 --output campo_offset_10.png"
echo ""

./compare_field_theory.py -f $PSI_FILE \
    --charge-pos $OFFSET_10_X $OFFSET_10_Y $OFFSET_10_Z \
    --size $SIZE \
    --mode line \
    --start 0.5 $OFFSET_10_Y $OFFSET_10_Z \
    --end 31.5 $OFFSET_10_Y $OFFSET_10_Z \
    --num-points 100 --output campo_offset_10.png 2>&1 | tail -5

echo ""
echo "=========================================="
echo "PRUEBAS COMPLETADAS"
echo "=========================================="
echo ""
echo "Archivos generados:"
echo "  - campo_no_offset.png  (offset = 0.0)"
echo "  - campo_offset_05.png  (offset = -0.5)"
echo "  - campo_offset_10.png  (offset = -1.0)"
echo ""
echo "Compara los tres gráficos y mira cuál tiene los campos"
echo "simulado y teórico mejor alineados."
echo ""
