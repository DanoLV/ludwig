#!/usr/bin/env python3
"""
Script para procesar múltiples archivos psi y generar gráficos del campo eléctrico.

Este script busca todos los archivos psi en un directorio y genera gráficos
para cada uno de ellos de forma automática.

Autor: Script generado para análisis de simulaciones Ludwig
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
import os
import glob
from pathlib import Path

# Importar clases del otro script
import importlib.util

def load_plot_module(script_path):
    """Carga el módulo plot_electric_field.py dinámicamente"""
    spec = importlib.util.spec_from_file_location("plot_electric_field", script_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def find_psi_files(directory, pattern='psi-*.001-001'):
    """
    Busca archivos psi en un directorio

    Args:
        directory: Directorio donde buscar
        pattern: Patrón de búsqueda (por defecto: psi-*.001-001)

    Returns:
        Lista de archivos encontrados, ordenados por nombre
    """
    search_path = os.path.join(directory, pattern)
    files = sorted(glob.glob(search_path))
    return files


def extract_timestep(filename):
    """
    Extrae el timestep del nombre del archivo psi

    Args:
        filename: Nombre del archivo (ej: psi-000010050.001-001)

    Returns:
        Timestep como entero
    """
    basename = os.path.basename(filename)
    # Formato: psi-NNNNNNNN.001-001
    parts = basename.split('-')
    if len(parts) >= 2:
        timestep_str = parts[1].split('.')[0]
        try:
            return int(timestep_str)
        except ValueError:
            return 0
    return 0


def main():
    parser = argparse.ArgumentParser(
        description='Procesar múltiples archivos psi y generar gráficos del campo eléctrico',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Ejemplos de uso:

1. Procesar todos los archivos psi en el directorio actual (plano XY):
   ./batch_plot_electric_field.py -d . -s 32 32 32 -m plane -p xy

2. Procesar archivos psi y generar líneas a lo largo del eje X:
   ./batch_plot_electric_field.py -d . -s 32 32 32 -m line --start 0 16 16 --end 31 16 16

3. Procesar con vectores del campo:
   ./batch_plot_electric_field.py -d . -s 32 32 32 -m plane -p xz -v -c magnitude

4. Especificar patrón de búsqueda personalizado:
   ./batch_plot_electric_field.py -d . -s 32 32 32 --pattern "psi-00001*.001-001" -m plane -p xy

5. Guardar en directorio específico:
   ./batch_plot_electric_field.py -d . -s 32 32 32 -m plane -p xy --output-dir graficos/
        """
    )

    # Argumentos obligatorios
    parser.add_argument('-d', '--directory', default='.',
                       help='Directorio donde buscar archivos psi (por defecto: directorio actual)')
    parser.add_argument('-s', '--size', nargs=3, type=int, required=True,
                       metavar=('NX', 'NY', 'NZ'),
                       help='Tamaño de la malla (nx ny nz)')
    parser.add_argument('-m', '--mode', choices=['plane', 'line', 'plane3d'], required=True,
                       help='Modo de graficación: plane (plano 2D), line (línea 1D), o plane3d (superficie 3D)')

    # Argumentos para búsqueda de archivos
    parser.add_argument('--pattern', default='psi-*.001-001',
                       help='Patrón de búsqueda de archivos (por defecto: psi-*.001-001)')
    parser.add_argument('--output-dir', default='electric_field_plots',
                       help='Directorio de salida (por defecto: electric_field_plots)')

    # Argumentos para modo plano
    parser.add_argument('-p', '--plane', choices=['xy', 'xz', 'yz'],
                       help='Plano a graficar (obligatorio para mode=plane)')
    parser.add_argument('-pos', '--position', type=int,
                       help='Posición del plano (índice). Por defecto: centro')
    parser.add_argument('-v', '--vectors', action='store_true',
                       help='Mostrar vectores del campo (solo para mode=plane)')
    parser.add_argument('--vector-stride', type=int, default=2,
                       help='Espaciado entre vectores (por defecto: 2)')

    # Argumentos para modo línea
    parser.add_argument('--start', nargs=3, type=float, metavar=('X', 'Y', 'Z'),
                       help='Punto inicial de la línea (obligatorio para mode=line)')
    parser.add_argument('--end', nargs=3, type=float, metavar=('X', 'Y', 'Z'),
                       help='Punto final de la línea (obligatorio para mode=line)')
    parser.add_argument('--num-points', type=int, default=100,
                       help='Número de puntos a lo largo de la línea (por defecto: 100)')

    # Argumentos compartidos
    parser.add_argument('-c', '--component',
                       choices=['magnitude', 'x', 'y', 'z', 'psi', 'all'],
                       default='magnitude',
                       help='Componente a graficar. Plane/Plane3D: magnitude/x/y/z/psi. '
                            'Line: all/magnitude/x/y/z/psi. (default: magnitude)')

    # Argumentos para modo 3D
    parser.add_argument('--colormap', default='viridis',
                       help='Mapa de colores para modo plane3d (default: viridis)')
    parser.add_argument('--elevation', type=float, default=30,
                       help='Ángulo de elevación para modo plane3d (default: 30)')
    parser.add_argument('--azimuth', type=float, default=-60,
                       help='Ángulo azimutal para modo plane3d (default: -60)')

    # Argumentos adicionales
    parser.add_argument('--max-files', type=int,
                       help='Procesar solo los primeros N archivos')
    parser.add_argument('--skip', type=int, default=1,
                       help='Procesar cada N archivos (por defecto: 1, procesar todos)')

    args = parser.parse_args()

    # Validar argumentos según el modo
    if args.mode in ['plane', 'plane3d'] and args.plane is None:
        parser.error(f"mode={args.mode} requiere especificar --plane")
    if args.mode == 'line' and (args.start is None or args.end is None):
        parser.error("mode=line requiere especificar --start y --end")

    # Cargar módulo de plotting
    script_dir = os.path.dirname(os.path.abspath(__file__))
    plot_module_path = os.path.join(script_dir, 'plot_electric_field.py')

    if not os.path.exists(plot_module_path):
        print(f"Error: No se encuentra el archivo plot_electric_field.py en {script_dir}")
        print("Asegúrate de que ambos scripts estén en el mismo directorio.")
        sys.exit(1)

    plot_module = load_plot_module(plot_module_path)
    PsiReader = plot_module.PsiReader
    ElectricFieldPlotter = plot_module.ElectricFieldPlotter

    # Buscar archivos psi
    print(f"Buscando archivos psi en {args.directory} con patrón {args.pattern}...")
    psi_files = find_psi_files(args.directory, args.pattern)

    if not psi_files:
        print(f"No se encontraron archivos psi en {args.directory}")
        sys.exit(1)

    # Aplicar filtros
    if args.skip > 1:
        psi_files = psi_files[::args.skip]

    if args.max_files:
        psi_files = psi_files[:args.max_files]

    print(f"Se encontraron {len(psi_files)} archivos a procesar.")

    # Crear directorio de salida
    os.makedirs(args.output_dir, exist_ok=True)
    print(f"Los gráficos se guardarán en: {args.output_dir}/")

    # Procesar cada archivo
    for i, psi_file in enumerate(psi_files):
        print(f"\n[{i+1}/{len(psi_files)}] Procesando {os.path.basename(psi_file)}...")

        # Leer archivo psi
        try:
            reader = PsiReader(psi_file, tuple(args.size))
            psi = reader.read()

            # Calcular campo eléctrico
            Ex, Ey, Ez = reader.compute_electric_field()

            # Crear graficador
            plotter = ElectricFieldPlotter(psi, Ex, Ey, Ez, tuple(args.size))

            # Generar nombre de archivo de salida
            timestep = extract_timestep(psi_file)
            basename = os.path.basename(psi_file)

            if args.mode == 'plane':
                output_filename = f"E_field_{args.plane}_{args.component}_{basename}.png"
                output_path = os.path.join(args.output_dir, output_filename)

                plotter.plot_plane(
                    plane=args.plane,
                    position=args.position,
                    component=args.component,
                    show_vectors=args.vectors,
                    vector_stride=args.vector_stride,
                    output=output_path
                )
            elif args.mode == 'plane3d':
                output_filename = f"E_field_3d_{args.plane}_{args.component}_{basename}.png"
                output_path = os.path.join(args.output_dir, output_filename)

                plotter.plot_plane_3d(
                    plane=args.plane,
                    position=args.position,
                    component=args.component,
                    colormap=args.colormap,
                    elevation=args.elevation,
                    azimuth=args.azimuth,
                    output=output_path
                )
            else:  # mode == 'line'
                output_filename = f"E_field_line_{args.component}_{basename}.png"
                output_path = os.path.join(args.output_dir, output_filename)

                plotter.plot_line(
                    start=args.start,
                    end=args.end,
                    num_points=args.num_points,
                    components=args.component,
                    output=output_path
                )

            print(f"   ✓ Guardado: {output_filename}")

        except Exception as e:
            print(f"   ✗ Error al procesar {psi_file}: {str(e)}")
            continue

    print(f"\n¡Proceso completado! Se generaron gráficos en {args.output_dir}/")


if __name__ == '__main__':
    main()
