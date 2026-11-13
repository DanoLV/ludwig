#!/bin/bash
#
# Ejemplos de uso del script de distribución de carga
#
# Instrucciones: Copia y pega los comandos que necesites, ajustando los parámetros

# =============================================================================
# CONFIGURACIÓN - AJUSTA ESTOS VALORES SEGÚN TU SIMULACIÓN
# =============================================================================

# Tamaño de la malla (obtener con: grep "size" input)
NX=32
NY=32
NZ=32

# Archivo qsi a analizar
QSI_FILE="qsi-000010250.001-001"

# Posición central
POS_CENTRO=$((NY / 2))

# =============================================================================
# EJEMPLO 1: ESTADÍSTICAS Y CARGA TOTAL
# =============================================================================

echo "=== 1. Calcular carga total y estadísticas ==="
# ./plot_charge_distribution.py -f $QSI_FILE -s $NX $NY $NZ -m stats

# =============================================================================
# EJEMPLO 2: PLANOS 2D
# =============================================================================

echo ""
echo "=== 2. Gráficos en planos 2D ==="

# 2a. Carga neta en plano XY central
echo "2a. Carga neta en plano XY"
# ./plot_charge_distribution.py -f $QSI_FILE -s $NX $NY $NZ -m plane -p xy -c net -o carga_neta_xy.png

# 2b. Especie 0 (cationes) en plano XZ
echo "2b. Cationes en plano XZ"
# ./plot_charge_distribution.py -f $QSI_FILE -s $NX $NY $NZ -m plane -p xz -c species0 -o cationes_xz.png

# 2c. Especie 1 (aniones) en plano YZ
echo "2c. Aniones en plano YZ"
# ./plot_charge_distribution.py -f $QSI_FILE -s $NX $NY $NZ -m plane -p yz -c species1 -o aniones_yz.png

# 2d. Densidad total en plano XY
echo "2d. Densidad total en plano XY"
# ./plot_charge_distribution.py -f $QSI_FILE -s $NX $NY $NZ -m plane -p xy -c total -o densidad_total_xy.png

# =============================================================================
# EJEMPLO 3: SUPERFICIES 3D
# =============================================================================

echo ""
echo "=== 3. Superficies 3D ==="

# 3a. Superficie 3D de carga neta
echo "3a. Carga neta en 3D"
# ./plot_charge_distribution.py -f $QSI_FILE -s $NX $NY $NZ -m plane3d -p xy -c net -o carga_neta_3d.png

# 3b. Superficie 3D de cationes con colormap viridis
echo "3b. Cationes en 3D con viridis"
# ./plot_charge_distribution.py -f $QSI_FILE -s $NX $NY $NZ -m plane3d -p xy -c species0 --colormap viridis -o cationes_3d.png

# 3c. Superficie 3D de aniones con colormap plasma
echo "3c. Aniones en 3D con plasma"
# ./plot_charge_distribution.py -f $QSI_FILE -s $NX $NY $NZ -m plane3d -p xz -c species1 --colormap plasma -o aniones_3d.png

# 3d. Carga neta con vista personalizada (seismic colormap)
echo "3d. Carga neta con colormap seismic"
# ./plot_charge_distribution.py -f $QSI_FILE -s $NX $NY $NZ -m plane3d -p xy -c net --colormap seismic --elevation 45 --azimuth -30 -o carga_neta_3d_seismic.png

# 3e. Vista superior de carga neta
echo "3e. Vista superior de carga neta"
# ./plot_charge_distribution.py -f $QSI_FILE -s $NX $NY $NZ -m plane3d -p xy -c net --elevation 90 --azimuth 0 -o carga_neta_vista_superior.png

# =============================================================================
# EJEMPLO 4: LÍNEAS 1D
# =============================================================================

echo ""
echo "=== 4. Perfiles de línea 1D ==="

# 4a. Todas las especies a lo largo del eje X
echo "4a. Todas las especies a lo largo del eje X"
# ./plot_charge_distribution.py -f $QSI_FILE -s $NX $NY $NZ -m line --start 0 $POS_CENTRO $POS_CENTRO --end $((NX-1)) $POS_CENTRO $POS_CENTRO -c all -o perfil_x_todas.png

# 4b. Solo carga neta a lo largo de la diagonal
echo "4b. Carga neta en diagonal"
# ./plot_charge_distribution.py -f $QSI_FILE -s $NX $NY $NZ -m line --start 0 0 0 --end $((NX-1)) $((NY-1)) $((NZ-1)) -c net -o perfil_diagonal_neta.png

# 4c. Especie 0 a lo largo del eje Z
echo "4c. Cationes a lo largo del eje Z"
# ./plot_charge_distribution.py -f $QSI_FILE -s $NX $NY $NZ -m line --start $POS_CENTRO $POS_CENTRO 0 --end $POS_CENTRO $POS_CENTRO $((NZ-1)) -c species0 -o perfil_z_cationes.png

# 4d. Especie 1 a lo largo del eje Y
echo "4d. Aniones a lo largo del eje Y"
# ./plot_charge_distribution.py -f $QSI_FILE -s $NX $NY $NZ -m line --start $POS_CENTRO 0 $POS_CENTRO --end $POS_CENTRO $((NY-1)) $POS_CENTRO -c species1 -o perfil_y_aniones.png

# 4e. Densidad total a lo largo del eje X con más puntos
echo "4e. Densidad total a lo largo del eje X"
# ./plot_charge_distribution.py -f $QSI_FILE -s $NX $NY $NZ -m line --start 0 $POS_CENTRO $POS_CENTRO --end $((NX-1)) $POS_CENTRO $POS_CENTRO -c total --num-points 200 -o perfil_x_total.png

# 4f. Perfil radial desde el centro (si hay partícula centrada)
echo "4f. Perfil radial desde el centro"
# ./plot_charge_distribution.py -f $QSI_FILE -s $NX $NY $NZ -m line --start $POS_CENTRO $POS_CENTRO $POS_CENTRO --end $((NX-1)) $POS_CENTRO $POS_CENTRO -c all -o perfil_radial.png

# =============================================================================
# EJEMPLO 5: COMPARACIÓN DE ESPECIES
# =============================================================================

echo ""
echo "=== 5. Comparar cationes y aniones ==="

# Generar gráficos lado a lado
# Cationes en XY
# ./plot_charge_distribution.py -f $QSI_FILE -s $NX $NY $NZ -m plane -p xy -c species0 -o comparacion_cationes_xy.png

# Aniones en XY
# ./plot_charge_distribution.py -f $QSI_FILE -s $NX $NY $NZ -m plane -p xy -c species1 -o comparacion_aniones_xy.png

# Carga neta en XY
# ./plot_charge_distribution.py -f $QSI_FILE -s $NX $NY $NZ -m plane -p xy -c net -o comparacion_neta_xy.png

# =============================================================================
# EJEMPLO 6: ANÁLISIS DE CAPA DOBLE ELÉCTRICA
# =============================================================================

echo ""
echo "=== 6. Análisis de capa doble eléctrica ==="

# 6a. Vista 3D de carga neta alrededor de partícula
echo "6a. Capa doble en 3D"
# ./plot_charge_distribution.py -f $QSI_FILE -s $NX $NY $NZ -m plane3d -p xy -c net --colormap seismic -o capa_doble_3d.png

# 6b. Perfil radial de carga (asumiendo partícula centrada)
echo "6b. Perfil radial de carga"
# ./plot_charge_distribution.py -f $QSI_FILE -s $NX $NY $NZ -m line --start $POS_CENTRO $POS_CENTRO $POS_CENTRO --end $((NX-1)) $POS_CENTRO $POS_CENTRO -c all -o capa_doble_perfil.png

# 6c. Corte vertical de carga neta
echo "6c. Corte vertical XZ de carga neta"
# ./plot_charge_distribution.py -f $QSI_FILE -s $NX $NY $NZ -m plane -p xz -c net -o capa_doble_xz.png

# =============================================================================
# EJEMPLO 7: ANÁLISIS TEMPORAL
# =============================================================================

echo ""
echo "=== 7. Análisis temporal (múltiples timesteps) ==="

# 7a. Verificar conservación de carga en el tiempo
echo "7a. Verificar conservación de carga"
# for file in qsi-*.001-001; do
#     echo "=== $file ==="
#     ./plot_charge_distribution.py -f $file -s $NX $NY $NZ -m stats | grep "Carga total"
#     echo ""
# done

# 7b. Generar serie de gráficos 3D
echo "7b. Serie de gráficos 3D"
# mkdir -p series_temporal_3d
# for file in qsi-*.001-001; do
#     timestep=$(echo $file | grep -o '[0-9]\{9\}')
#     ./plot_charge_distribution.py -f $file -s $NX $NY $NZ -m plane3d -p xy -c net -o series_temporal_3d/carga_$timestep.png
# done

# =============================================================================
# EJEMPLO 8: COMBINACIÓN CON CAMPO ELÉCTRICO
# =============================================================================

echo ""
echo "=== 8. Análisis combinado carga + campo eléctrico ==="

# 8a. Distribución de carga
echo "8a. Distribución de carga"
# ./plot_charge_distribution.py -f qsi-000010250.001-001 -s $NX $NY $NZ -m plane -p xy -c net -o carga_xy.png

# 8b. Campo eléctrico resultante (requiere archivo psi correspondiente)
echo "8b. Campo eléctrico resultante"
# ./plot_electric_field.py -f psi-000010250.001-001 -s $NX $NY $NZ -m plane -p xy -c magnitude -v -o campo_xy.png

# 8c. Potencial eléctrico
echo "8c. Potencial eléctrico"
# ./plot_electric_field.py -f psi-000010250.001-001 -s $NX $NY $NZ -m plane -p xy -c psi -o potencial_xy.png

# =============================================================================
# EJEMPLO 9: DIFERENTES PLANOS
# =============================================================================

echo ""
echo "=== 9. Analizar los tres planos ortogonales ==="

# Carga neta en los tres planos
# ./plot_charge_distribution.py -f $QSI_FILE -s $NX $NY $NZ -m plane -p xy -c net -o carga_plano_xy.png
# ./plot_charge_distribution.py -f $QSI_FILE -s $NX $NY $NZ -m plane -p xz -c net -o carga_plano_xz.png
# ./plot_charge_distribution.py -f $QSI_FILE -s $NX $NY $NZ -m plane -p yz -c net -o carga_plano_yz.png

# =============================================================================
# EJEMPLO 10: EXPORTAR DATOS PARA ANÁLISIS EXTERNO
# =============================================================================

echo ""
echo "=== 10. Información del script ==="

# Ver ayuda completa
echo "Ver ayuda:"
echo "  ./plot_charge_distribution.py --help"

# Ver estadísticas rápidas
echo "Estadísticas rápidas:"
echo "  ./plot_charge_distribution.py -f $QSI_FILE -s $NX $NY $NZ -m stats"

# Listar archivos qsi
echo "Listar archivos qsi:"
echo "  ls -lh qsi-*.001-001"

# Contar archivos
echo "Contar archivos:"
echo "  ls qsi-*.001-001 | wc -l"

echo ""
echo "Para usar estos ejemplos, descomenta (quita el #) las líneas que necesites."
echo "Recuerda ajustar las variables NX, NY, NZ y QSI_FILE al inicio del script."
