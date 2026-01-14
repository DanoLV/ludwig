#!/bin/bash

# Script para VISUALIZAR cómo se organizarían las carpetas de simulaciones por parámetro
# (DRY RUN - no mueve nada, solo muestra qué haría)
# Uso: ./organize_by_parameter_dryrun.sh <directorio_base> <parametro> [prefijo]

# Verificar argumentos
if [ $# -lt 2 ]; then
    echo "Uso: $0 <directorio_base> <parametro> [prefijo]"
    echo "Ejemplo: $0 /home/bater/Sim/ludwig/microgel/local/single-peskin-Lx_32_Ly_32_Lz_32-D3Q19-stencil_19-weight_PB_Fn-rhoel rhoel"
    echo "Ejemplo con prefijo: $0 /home/bater/Sim/ludwig/microgel/local/single-peskin-Lx_32_Ly_32_Lz_32-D3Q19-stencil_19-weight_PB_Sn rhoel PB_Sn"
    exit 1
fi

BASE_DIR="$1"
PARAM="$2"
PREFIX="${3:-}"  # Prefijo opcional (vacío por defecto)

# Verificar que el directorio existe
if [ ! -d "$BASE_DIR" ]; then
    echo "Error: El directorio $BASE_DIR no existe"
    exit 1
fi

cd "$BASE_DIR" || exit 1

echo "====================================="
echo "DRY RUN - No se moverá nada"
echo "====================================="
echo "Directorio base: $BASE_DIR"
echo "Parámetro: $PARAM"
if [ -n "$PREFIX" ]; then
    echo "Prefijo: $PREFIX"
fi
echo "----------------------------------------"

# Contador de carpetas
would_move=0
would_skip=0

# Iterar sobre todas las carpetas que contienen el parámetro en su nombre
for dir in */; do
    # Eliminar la barra final
    dir="${dir%/}"

    # Saltar la carpeta scripts y subcarpetas ya organizadas
    if [ "$dir" = "scripts" ] || [ ! -d "$dir" ]; then
        continue
    fi

    # Verificar si el nombre contiene el parámetro
    if [[ "$dir" =~ ${PARAM}_ ]]; then
        # Extraer el valor del parámetro
        if [[ "$dir" =~ ${PARAM}_([0-9]+\.[0-9]+E[+-][0-9]+) ]]; then
            param_value="${BASH_REMATCH[1]}"

            # Crear nombre de subcarpeta con o sin prefijo
            if [ -n "$PREFIX" ]; then
                subdir="${PREFIX}-${PARAM}_${param_value}"
            else
                subdir="${PARAM}_${param_value}"
            fi

            echo "📁 $dir"
            echo "   → Se movería a: $subdir/"

            if [ "$dir" != "$subdir" ] && [ ! -d "$subdir/$dir" ]; then
                ((would_move++))
            else
                echo "   ⚠ Ya existe en destino, se saltaría"
                ((would_skip++))
            fi
            echo ""
        else
            echo "⚠ $dir - No se pudo extraer el valor de $PARAM"
            ((would_skip++))
        fi
    fi
done

echo "----------------------------------------"
echo "Resumen (DRY RUN):"
echo "  Se moverían: $would_move carpetas"
echo "  Se saltarían: $would_skip carpetas"
echo ""
echo "Para ejecutar realmente, usa: organize_by_parameter.sh"
