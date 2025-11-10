#!/bin/bash
#------------------------------------------------------------------------------------
# Ejemplo de uso de runmulti.sh
#------------------------------------------------------------------------------------
# Este script muestra diferentes ejemplos de cómo usar runmulti.sh
#------------------------------------------------------------------------------------

echo "=========================================="
echo "Ejemplos de uso de runmulti.sh"
echo "=========================================="
echo ""
echo "NOTA: Estos ejemplos están comentados. Descomenta el que quieras ejecutar."
echo ""

# Ejemplo 1: Variar la carga de la partícula
# Ejecuta 4 simulaciones con diferentes valores de carga
# Carpetas resultantes: test-charge_0.5, test-charge_1.0, test-charge_1.5, test-charge_2.0
#
# ./runmulti.sh \
#   --charge 0.5,1.0,1.5,2.0 \
#   --base-dir test-charge \
#   -i 0 -n 10000 -s 500 \
#   -x 32 -y 32 -v 0.5 \
#   -t fe_electro -e 0.0_0.0_0.1

# Ejemplo 2: Variar la posición X de la partícula
# Ejecuta 5 simulaciones con diferentes posiciones en X
# Las posiciones Y y Z se mantienen según el archivo config original
#
# ./runmulti.sh \
#   --position-x 16.0,16.25,16.5,16.75,17.0 \
#   --base-dir test-position-x \
#   -i 0 -n 10000 -s 500 \
#   -x 32 -y 32 -v 0.5

# Ejemplo 3: Variar el parámetro AL
# Ejecuta 4 simulaciones con diferentes valores de AL
#
# ./runmulti.sh \
#   --al 0.4,0.5,0.6,0.7 \
#   --base-dir test-al \
#   -i 0 -n 10000 -s 500 \
#   -x 32 -y 32 -v 0.5 \
#   -t fe_electro

# Ejemplo 4: Variar carga Y posición (producto cartesiano)
# Ejecuta 4 simulaciones: 2 cargas × 2 posiciones = 4 combinaciones
# Carpetas: test-combined_1.0_16.0, test-combined_1.0_17.0,
#           test-combined_2.0_16.0, test-combined_2.0_17.0
#
# ./runmulti.sh \
#   --charge 1.0,2.0 \
#   --position-x 16.0,17.0 \
#   --base-dir test-combined \
#   -i 0 -n 10000 -s 500 \
#   -x 32 -y 32 -v 0.5 \
#   -t fe_electro -e 0.0_0.0_0.1

# Ejemplo 5: Simulación con MPI
# Ejecuta con 4 procesos MPI en grid 2x2x1, variando la carga
#
# ./runmulti.sh \
#   --charge 1.0,1.5,2.0 \
#   --base-dir test-mpi-charge \
#   -i 0 -n 50000 -s 1000 \
#   -x 64 -y 64 -v 0.5 \
#   -m 4 -g 2_2_1 -c 10 \
#   -t fe_electro -e 0.0_0.0_0.1

# Ejemplo 6: Variar posición completa (usando --position)
# Útil para posiciones específicas en una diagonal u otras trayectorias
#
# ./runmulti.sh \
#   --position 16.0_16.0_16.0,16.5_16.5_16.5,17.0_17.0_17.0 \
#   --base-dir test-diagonal \
#   -i 0 -n 10000 -s 500 \
#   -x 32 -y 32 -v 0.5 \
#   -t fe_electro

# Ejemplo 7: Escaneo en eje X (similar a tus carpetas existentes)
# Replica el patrón que veo en tus simulaciones "diagonal"
#
# ./runmulti.sh \
#   --position-x 16.00,16.05,16.10,16.20,16.25,16.30,16.40,16.45,16.50 \
#   --base-dir single-r_0.8-v_0.02-kT_0.0005-lb_fluctuation_0-fe_electro-q_1.0-e_0.0_0.0_0.0-p250-N_15k-pos \
#   -i 0 -n 15000 -s 250 \
#   -x 32 -y 32 -v 0.02 \
#   -k 0.0005 -w 0 \
#   -t fe_electro -e 0.0_0.0_0.0 \
#   -u y

# Ejemplo 8: Variar parámetros de input (temperatura, densidad, campo)
# Estos parámetros se pasan directamente a runbg.sh
# Necesitas un wrapper script si quieres variarlos
#
# #!/bin/bash
# # Variar temperatura
# for temp in 0.0005 0.001 0.002; do
#   ./runmulti.sh \
#     --charge 1.0 \
#     --base-dir test-temp-${temp} \
#     -i 0 -n 10000 -s 500 \
#     -x 32 -y 32 -v 0.5 \
#     -k $temp \
#     -t fe_electro -e 0.0_0.0_0.1
# done

# Ejemplo 9: Ejemplo completo con múltiples parámetros de runbg.sh
# Demuestra que TODOS los parámetros de runbg.sh funcionan
#
# ./runmulti.sh \
#   --charge 1.0,1.5 \
#   --position-x 16.0,17.0 \
#   --base-dir full-example \
#   -i 0 \
#   -n 10000 \
#   -s 500 \
#   -x 32 \
#   -y 32 \
#   -v 0.5 \
#   -r 0.8 \
#   -k 0.0005 \
#   -t fe_electro \
#   -e 0.0_0.0_0.1 \
#   -w 0 \
#   -q M10 \
#   -z petsc \
#   -u y

# Ejemplo 10: Ejecución en PARALELO (todas a la vez)
# ADVERTENCIA: Asegúrate de tener suficientes recursos (CPU, memoria)
#
# ./runmulti.sh \
#   --charge 1.0,1.5,2.0 \
#   --base-dir test-parallel \
#   --parallel \
#   -i 0 -n 10000 -s 500 \
#   -x 32 -y 32 -v 0.5

# Ejemplo 11: Ejecución en paralelo LIMITADO (recomendado)
# Solo 3 simulaciones a la vez, más seguro para recursos
#
# ./runmulti.sh \
#   --charge 0.5,1.0,1.5,2.0,2.5,3.0 \
#   --base-dir test-parallel-limited \
#   --parallel --max-parallel 3 \
#   -i 0 -n 10000 -s 500 \
#   -x 32 -y 32 -v 0.5

# Ejemplo 12: Test rápido (solo 50 pasos para verificar que funciona)
# Descomenta este para hacer una prueba rápida
#
# ./runmulti.sh \
#   --charge 1.0,1.5 \
#   --base-dir quick-test \
#   -i 0 -n 50 -s 10 \
#   -x 32 -y 32 -v 0.5

echo "Para ejecutar algún ejemplo, edita este archivo y descomenta el ejemplo deseado."
echo ""
echo "IMPORTANTE: Todos los parámetros de runbg.sh están disponibles:"
echo "  -i, -n, -s, -x, -y, -v, -r, -k, -t, -e, -m, -g, -c, -w, -q, -z,"
echo "  -a, -b, -j, -u, -l, -d, -f, -p"
echo ""
echo "NUEVAS OPCIONES DE PARALELIZACIÓN:"
echo "  --parallel              : Ejecutar todas las simulaciones simultáneamente"
echo "  --max-parallel N        : Limitar a N simulaciones concurrentes"
echo ""
echo "EJEMPLOS:"
echo "  Secuencial (default):   ./runmulti.sh --charge 1.0,2.0 --base-dir test ..."
echo "  Paralelo ilimitado:     ./runmulti.sh --charge 1.0,2.0 --base-dir test --parallel ..."
echo "  Paralelo limitado:      ./runmulti.sh --charge 1.0,2.0,3.0 --base-dir test --parallel --max-parallel 2 ..."
echo ""
echo "Ver README_runmulti.md para documentación completa."
echo ""
