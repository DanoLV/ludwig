#!/bin/bash

# Script para organizar carpetas de simulaciones por parámetro
# Uso: ./organize_by_parameter.sh <directorio_base> <parametro> [prefijo]
# Ejemplo: ./organize_by_parameter.sh /path/to/simulations rhoel PB_Sn

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

echo "Organizando carpetas en: $BASE_DIR"
echo "Parámetro: $PARAM"
if [ -n "$PREFIX" ]; then
    echo "Prefijo: $PREFIX"
fi
echo "----------------------------------------"

# Contador de carpetas movidas
moved_count=0
skipped_count=0

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
        # Extraer el valor del parámetro usando regex
        # Busca el patrón: parametro_valor donde valor puede contener puntos, E+, E-, números
        if [[ "$dir" =~ ${PARAM}_([0-9]+\.[0-9]+E[+-][0-9]+) ]]; then
            param_value="${BASH_REMATCH[1]}"

            # Crear nombre de subcarpeta con o sin prefijo
            if [ -n "$PREFIX" ]; then
                subdir="${PREFIX}-${PARAM}_${param_value}"
            else
                subdir="${PARAM}_${param_value}"
            fi

            echo "Procesando: $dir"
            echo "  Valor de $PARAM: $param_value"
            echo "  Subcarpeta destino: $subdir"

            # Crear subcarpeta si no existe
            if [ ! -d "$subdir" ]; then
                mkdir -p "$subdir"
                echo "  ✓ Subcarpeta creada: $subdir"
            fi

            # Mover la carpeta a la subcarpeta
            if [ "$dir" != "$subdir" ] && [ ! -d "$subdir/$dir" ]; then
                mv "$dir" "$subdir/"
                echo "  ✓ Movido a: $subdir/"
                ((moved_count++))
            else
                echo "  ⚠ Ya existe en destino, saltando"
                ((skipped_count++))
            fi
            echo ""
        else
            echo "⚠ No se pudo extraer el valor de $PARAM de: $dir"
            ((skipped_count++))
        fi
    fi
done

echo "----------------------------------------"
echo "Resumen:"
echo "  Carpetas movidas: $moved_count"
echo "  Carpetas saltadas: $skipped_count"
echo "Organización completada."
