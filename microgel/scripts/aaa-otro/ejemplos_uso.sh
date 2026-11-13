#!/bin/bash
#
# Ejemplos de uso de los scripts de graficación del campo eléctrico
#
# Instrucciones: Copia y pega los comandos que necesites, ajustando los parámetros

# =============================================================================
# CONFIGURACIÓN - AJUSTA ESTOS VALORES SEGÚN TU SIMULACIÓN
# =============================================================================

# Tamaño de la malla (obtener con: grep "size" input)
NX=32
NY=32
NZ=32

# Archivo psi a analizar (o usa * para patrones)
PSI_FILE="psi-000010050.001-001"

# =============================================================================
# EJEMPLOS: SCRIPT INDIVIDUAL (plot_electric_field.py)
# =============================================================================

echo "=== Ejemplos de uso del script individual ==="

# 1. Plano XY central con magnitud del campo y vectores
echo "1. Plano XY central con vectores"
# ./plot_electric_field.py -f $PSI_FILE -s $NX $NY $NZ -m plane -p xy -v -o campo_xy.png

# 2. Plano XZ mostrando componente Ez en posición y=10
echo "2. Plano XZ con componente Ez"
# ./plot_electric_field.py -f $PSI_FILE -s $NX $NY $NZ -m plane -p xz -pos 10 -c z -o campo_xz_Ez.png

# 3. Plano YZ mostrando el potencial eléctrico
echo "3. Potencial en plano YZ"
# ./plot_electric_field.py -f $PSI_FILE -s $NX $NY $NZ -m plane -p yz -c psi -o potencial_yz.png

# 4. Todas las componentes a lo largo del eje X (centro YZ)
echo "4. Todas las componentes a lo largo del eje X"
POS_CENTRO=$((NY / 2))
# ./plot_electric_field.py -f $PSI_FILE -s $NX $NY $NZ -m line --start 0 $POS_CENTRO $POS_CENTRO --end $((NX-1)) $POS_CENTRO $POS_CENTRO -c all -o campo_eje_x_all.png

# 5. Solo componente Ez a lo largo de la diagonal
echo "5. Solo Ez a lo largo de la diagonal"
# ./plot_electric_field.py -f $PSI_FILE -s $NX $NY $NZ -m line --start 0 0 0 --end $((NX-1)) $((NY-1)) $((NZ-1)) -c z -o campo_diagonal_Ez.png

# 6. Magnitud del campo a lo largo del eje Z (centro XY)
echo "6. Magnitud a lo largo del eje Z"
# ./plot_electric_field.py -f $PSI_FILE -s $NX $NY $NZ -m line --start $POS_CENTRO $POS_CENTRO 0 --end $POS_CENTRO $POS_CENTRO $((NZ-1)) -c magnitude -o campo_eje_z_mag.png

# 6b. Solo potencial a lo largo del eje Y
echo "6b. Potencial a lo largo del eje Y"
# ./plot_electric_field.py -f $PSI_FILE -s $NX $NY $NZ -m line --start $POS_CENTRO 0 $POS_CENTRO --end $POS_CENTRO $((NY-1)) $POS_CENTRO -c psi -o potencial_eje_y.png

# 7. Superficie 3D del potencial en plano XY
echo "7. Superficie 3D del potencial"
# ./plot_electric_field.py -f $PSI_FILE -s $NX $NY $NZ -m plane3d -p xy -c psi -o potencial_3d_xy.png

# 8. Superficie 3D de la magnitud del campo con colormap plasma
echo "8. Superficie 3D del campo con plasma colormap"
# ./plot_electric_field.py -f $PSI_FILE -s $NX $NY $NZ -m plane3d -p xz -c magnitude --colormap plasma -o campo_3d_plasma.png

# 9. Superficie 3D con vista personalizada
echo "9. Superficie 3D con vista personalizada"
# ./plot_electric_field.py -f $PSI_FILE -s $NX $NY $NZ -m plane3d -p xy -c psi --elevation 45 --azimuth -30 -o potencial_3d_vista.png

# =============================================================================
# EJEMPLOS: PROCESAMIENTO POR LOTES (batch_plot_electric_field.py)
# =============================================================================

echo ""
echo "=== Ejemplos de procesamiento por lotes ==="

# 7. Procesar todos los archivos psi en plano XY
echo "7. Procesar todos los archivos (plano XY)"
# ./batch_plot_electric_field.py -d . -s $NX $NY $NZ -m plane -p xy -v --output-dir graficos_xy/

# 8. Procesar archivos con timestep específico (ej: 00001XXXX)
echo "8. Procesar archivos con patrón específico"
# ./batch_plot_electric_field.py -d . -s $NX $NY $NZ --pattern "psi-00001*.001-001" -m plane -p xy --output-dir graficos_t1/

# 9. Procesar cada 5 archivos (útil para datasets grandes)
echo "9. Procesar cada 5 archivos"
# ./batch_plot_electric_field.py -d . -s $NX $NY $NZ -m plane -p xz --skip 5 --output-dir graficos_xz_skip/

# 10. Procesar solo los primeros 10 archivos
echo "10. Procesar primeros 10 archivos"
# ./batch_plot_electric_field.py -d . -s $NX $NY $NZ -m plane -p xy --max-files 10 --output-dir graficos_primeros10/

# 11. Todas las componentes a lo largo del eje X para todos los archivos
echo "11. Líneas (todas componentes) para todos los archivos"
# ./batch_plot_electric_field.py -d . -s $NX $NY $NZ -m line --start 0 $POS_CENTRO $POS_CENTRO --end $((NX-1)) $POS_CENTRO $POS_CENTRO -c all --output-dir graficos_lineas_all/

# 11b. Solo componente Ez a lo largo del eje X para todos los archivos
echo "11b. Solo Ez en líneas para todos los archivos"
# ./batch_plot_electric_field.py -d . -s $NX $NY $NZ -m line --start 0 $POS_CENTRO $POS_CENTRO --end $((NX-1)) $POS_CENTRO $POS_CENTRO -c z --output-dir graficos_lineas_Ez/

# 12. Superficies 3D del potencial para todos los archivos
echo "12. Superficies 3D del potencial"
# ./batch_plot_electric_field.py -d . -s $NX $NY $NZ -m plane3d -p xy -c psi --output-dir graficos_3d_psi/

# 13. Superficies 3D de la magnitud del campo con colormap personalizado
echo "13. Superficies 3D con colormap plasma"
# ./batch_plot_electric_field.py -d . -s $NX $NY $NZ -m plane3d -p xz -c magnitude --colormap plasma --output-dir graficos_3d_plasma/

# =============================================================================
# EJEMPLOS: COMPARAR DIFERENTES COMPONENTES
# =============================================================================

echo ""
echo "=== Comparar diferentes componentes del campo ==="

# 12. Generar gráficos de todas las componentes + magnitud
echo "12. Todas las componentes en el mismo plano"
# ./plot_electric_field.py -f $PSI_FILE -s $NX $NY $NZ -m plane -p xy -c magnitude -o E_magnitud.png
# ./plot_electric_field.py -f $PSI_FILE -s $NX $NY $NZ -m plane -p xy -c x -o E_componente_x.png
# ./plot_electric_field.py -f $PSI_FILE -s $NX $NY $NZ -m plane -p xy -c y -o E_componente_y.png
# ./plot_electric_field.py -f $PSI_FILE -s $NX $NY $NZ -m plane -p xy -c z -o E_componente_z.png
# ./plot_electric_field.py -f $PSI_FILE -s $NX $NY $NZ -m plane -p xy -c psi -o potencial.png

# =============================================================================
# EJEMPLOS: ANÁLISIS EN DIFERENTES PLANOS
# =============================================================================

echo ""
echo "=== Análisis en los tres planos principales ==="

# 13. Generar gráficos en los tres planos
echo "13. Tres planos ortogonales"
# ./plot_electric_field.py -f $PSI_FILE -s $NX $NY $NZ -m plane -p xy -v -o plano_xy.png
# ./plot_electric_field.py -f $PSI_FILE -s $NX $NY $NZ -m plane -p xz -v -o plano_xz.png
# ./plot_electric_field.py -f $PSI_FILE -s $NX $NY $NZ -m plane -p yz -v -o plano_yz.png

# =============================================================================
# EJEMPLOS: RESTAR CAMPO EXTERNO
# =============================================================================

echo ""
echo "=== Restar campo eléctrico externo ==="

# 14a. Restar campo externo en dirección Z (ej: E_ext = 1.0 en Z)
echo "14a. Restar campo externo Ez = 1.0"
# ./plot_electric_field.py -f $PSI_FILE -s $NX $NY $NZ -m plane -p xy -c magnitude --external-field 0 0 1.0 -o campo_sin_externo_xy.png

# 14b. Comparar campo total vs campo sin externo (línea)
echo "14b. Comparación con y sin campo externo"
# Campo total
# ./plot_electric_field.py -f $PSI_FILE -s $NX $NY $NZ -m line --start 0 $POS_CENTRO $POS_CENTRO --end $((NX-1)) $POS_CENTRO $POS_CENTRO -c all -o campo_total_linea.png
# Campo sin externo
# ./plot_electric_field.py -f $PSI_FILE -s $NX $NY $NZ -m line --start 0 $POS_CENTRO $POS_CENTRO --end $((NX-1)) $POS_CENTRO $POS_CENTRO -c all --external-field 0 0 1.0 -o campo_cargas_linea.png

# 14c. Visualizar solo polarización de partícula (3D)
echo "14c. Polarización en 3D (sin campo externo)"
# ./plot_electric_field.py -f $PSI_FILE -s $NX $NY $NZ -m plane3d -p xy --external-field 0 0 1.0 -o polarizacion_3d.png

# 14d. Restar campo externo en dirección X (electroforesis horizontal)
echo "14d. Campo sin externo en dirección X"
# ./plot_electric_field.py -f $PSI_FILE -s $NX $NY $NZ -m plane -p xz -c magnitude --external-field 0.005 0 0 -o campo_sin_ext_x.png

# 14e. Visualizar capa doble sin campo aplicado
echo "14e. Capa doble eléctrica sin campo externo"
# ./plot_electric_field.py -f $PSI_FILE -s $NX $NY $NZ -m plane3d -p xz -c magnitude --external-field 0 0 1.0 --colormap RdBu_r -o capa_doble_sin_ext.png

# =============================================================================
# EJEMPLOS AVANZADOS
# =============================================================================

echo ""
echo "=== Ejemplos avanzados ==="

# 14. Vectores densos (stride=1) en región específica
echo "14. Vectores densos"
# ./plot_electric_field.py -f $PSI_FILE -s $NX $NY $NZ -m plane -p xy --vector-stride 1 -o campo_vectores_densos.png

# 15. Análisis de múltiples posiciones en el mismo plano
echo "15. Múltiples posiciones"
for POS in 8 16 24; do
    # ./plot_electric_field.py -f $PSI_FILE -s $NX $NY $NZ -m plane -p xy -pos $POS -o plano_xy_z${POS}.png
    echo "  - Posición z=$POS"
done

# 16. Crear video con todos los timesteps (requiere ffmpeg)
echo "16. Crear video (requiere procesar archivos primero)"
# ./batch_plot_electric_field.py -d . -s $NX $NY $NZ -m plane -p xy --output-dir video_frames/
# cd video_frames && ffmpeg -framerate 10 -pattern_type glob -i 'E_field_*.png' -c:v libx264 -pix_fmt yuv420p ../campo_electrico.mp4 && cd ..

# =============================================================================
# INFORMACIÓN ÚTIL
# =============================================================================

echo ""
echo "=== Comandos útiles ==="

# Ver tamaño de la malla desde el archivo input
echo "Ver tamaño de malla:"
echo "  grep 'size' input"

# Listar archivos psi disponibles
echo "Listar archivos psi:"
echo "  ls -lh psi-*.001-001"

# Contar archivos psi
echo "Contar archivos:"
echo "  ls psi-*.001-001 | wc -l"

# Ver ayuda detallada
echo "Ver ayuda:"
echo "  ./plot_electric_field.py --help"
echo "  ./batch_plot_electric_field.py --help"

echo ""
echo "Para usar estos ejemplos, descomenta (quita el #) las líneas que necesites."
echo "Recuerda ajustar las variables NX, NY, NZ al inicio del script según tu simulación."
