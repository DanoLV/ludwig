#! /usr/bin/python3
"""
Script para combinar plots de velocidad agrupados por un parámetro específico
Lee los archivos datosfluid.csv y genera plots combinados

Uso:
    ./plot_velocities_grouped.py -d <directorio_padre> --group-by <parámetro> [opciones]

Parámetros disponibles para agrupar:
    - stencil: Agrupa por tipo de stencil (7, 19, 27, etc.)
    - rhoel: Agrupa por valor de rhoel
    - pos: Agrupa por posición
    - q: Agrupa por carga
    - eps: Agrupa por epsilon
    - etc.

Ejemplo:
    ./plot_velocities_grouped.py -d /path/to/parent/folder --group-by stencil
    ./plot_velocities_grouped.py -d /path/to/parent/folder --group-by rhoel --plot-type particle
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
import os
import re
from pathlib import Path
from collections import defaultdict


def limites_y(axis, x_min, x_max):
    """Ajusta los límites Y automáticamente para el rango de X especificado"""
    axis.set_xlim(x_min, x_max)

    ys_in_range = []
    for line in axis.get_lines():
        xd = np.asarray(line.get_xdata())
        yd = np.asarray(line.get_ydata())
        mask = (xd >= x_min) & (xd <= x_max) & np.isfinite(yd)
        if np.any(mask):
            ys_in_range.append(yd[mask])

    if ys_in_range:
        all_y = np.hstack(ys_in_range)
        y_min, y_max = all_y.min(), all_y.max()
        if y_min == y_max:
            pad = abs(y_min) * 0.01 if y_min != 0 else 1e-6
        else:
            pad = 0.05 * (y_max - y_min)
        axis.set_ylim(y_min - pad, y_max + pad)

    return 0


def parse_directory_name(dirname):
    """
    Extrae parámetros del nombre del directorio
    Ejemplo: single-peskin-Lx_32_Ly_32_Lz_32-D3Q19-stencil_var-n_2000-s_10-L_32-q_1.0-pos_16.20_16.50_16.50-e_0.0_0.0_0.0-kT_0.00001-rho_0.8-eta_0.05-eps_1.0e4-rhoel_0.000E+00-stencil_7
    """
    params = {}

    # Patrones comunes
    patterns = {
        'stencil': r'stencil_(\d+)',
        'rhoel': r'rhoel_([\d.E+-]+)',
        'n': r'-n_(\d+)',
        's': r'-s_(\d+)',
        'L': r'-L_(\d+)',
        'q': r'-q_([\d.]+)',
        'pos': r'pos_([\d.]+_[\d.]+_[\d.]+)',
        'e': r'-e_([\d.]+_[\d.]+_[\d.]+)',
        'kT': r'kT_([\d.E+-]+)',
        'rho': r'rho_([\d.]+)',
        'eta': r'eta_([\d.]+)',
        'eps': r'eps_([\d.E+-]+)',
        'Lx': r'Lx_(\d+)',
        'Ly': r'Ly_(\d+)',
        'Lz': r'Lz_(\d+)',
    }

    for param_name, pattern in patterns.items():
        match = re.search(pattern, dirname)
        if match:
            params[param_name] = match.group(1)

    return params


def find_data_files(parent_dir, depth=3):
    """
    Busca recursivamente archivos datosfluid.csv
    Retorna lista de (path_to_csv, params_dict)
    """
    parent_path = Path(parent_dir)
    data_files = []

    # Buscar archivos datosfluid.csv recursivamente
    for csv_file in parent_path.rglob('datosfluid.csv'):
        # Extraer parámetros de todos los directorios en la ruta
        params = {}
        for parent in csv_file.parents:
            if parent == parent_path:
                break
            dir_params = parse_directory_name(parent.name)
            params.update(dir_params)

        data_files.append((csv_file, params))

    return data_files


def read_velocity_data(data_file):
    """Lee el archivo datosfluid.csv y retorna los datos"""
    if not os.path.exists(data_file):
        return None

    try:
        data = np.genfromtxt(data_file, delimiter=',', skip_header=1)
        return data
    except Exception as e:
        print(f"Error leyendo {data_file}: {e}")
        return None


def group_data_by_parameter(data_files, group_params):
    """
    Agrupa archivos de datos por uno o varios parámetros
    group_params puede ser una lista de parámetros o un solo parámetro
    """
    # Asegurar que group_params sea una lista
    if isinstance(group_params, str):
        group_params = [group_params]

    grouped = defaultdict(list)

    for csv_file, params in data_files:
        # Verificar que todos los parámetros existan
        missing_params = [p for p in group_params if p not in params]
        if missing_params:
            print(f"Advertencia: {csv_file} no tiene parámetro(s) {missing_params}")
            continue

        # Crear clave compuesta si hay múltiples parámetros
        if len(group_params) == 1:
            key = params[group_params[0]]
        else:
            # Clave compuesta: "param1=value1, param2=value2"
            key = ", ".join([f"{p}={params[p]}" for p in group_params])

        grouped[key].append((csv_file, params))

    return dict(grouped)


# Set command line parameters
argParser = argparse.ArgumentParser(
    description="Combina plots de velocidad agrupados por parámetro"
)
argParser.add_argument(
    "-d", "--dir",
    required=True,
    help="Directorio padre que contiene los datos"
)
argParser.add_argument(
    "--group-by",
    required=True,
    nargs='+',
    help="Parámetro(s) para agrupar. Puede ser uno o varios separados por espacio (ej: stencil, stencil rhoel, stencil rhoel pos)"
)
argParser.add_argument(
    "--plot-type",
    choices=['both', 'particle', 'fluid'],
    default='particle',
    help="Tipo de plot: both (partícula y fluido), particle (solo partícula), fluid (solo fluido)"
)
argParser.add_argument(
    "--out_dir",
    help="Directorio de salida para los plots (se usa si no se especifica -o)",
    default="."
)
argParser.add_argument(
    "-o", "--output",
    help="Ruta completa del archivo de salida (ejemplo: ./plots/mi_grafico.png). Si se especifica, se ignora --out_dir"
)
argParser.add_argument(
    "--component",
    choices=['all', 'x', 'y', 'z'],
    default='all',
    help="Componente de velocidad a plotear (all, x, y, z)"
)

# Parse command line arguments
try:
    args = argParser.parse_args()
except:
    sys.exit("No se pudieron leer los parámetros de línea de comandos")

# Validar directorio padre
parent_dir = args.dir
if not os.path.exists(parent_dir):
    sys.exit(f"El directorio {parent_dir} no existe")

# Buscar archivos de datos
print(f"Buscando archivos datosfluid.csv en {parent_dir}...")
data_files = find_data_files(parent_dir)

if not data_files:
    sys.exit(f"No se encontraron archivos datosfluid.csv en {parent_dir}")

print(f"Encontrados {len(data_files)} archivos de datos")

# Agrupar por parámetro(s)
grouped_data = group_data_by_parameter(data_files, args.group_by)

if not grouped_data:
    params_str = "+".join(args.group_by)
    sys.exit(f"No se pudieron agrupar datos por parámetro(s) '{params_str}'")

# Mostrar información de agrupación
params_str = " + ".join(args.group_by)
print(f"\nDatos agrupados por '{params_str}':")
for key, files in sorted(grouped_data.items()):
    print(f"  {key}: {len(files)} archivo(s)")

# Archivo de salida
if args.output:
    # Si se especificó -o, usar esa ruta completa
    output_file = args.output
    # Crear el directorio si no existe
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Creado directorio: {output_dir}")
else:
    # Si no se especificó -o, usar out_dir con nombre automático
    params_str = "_".join(args.group_by)
    output_file = os.path.join(args.out_dir, f"velocidades_by_{params_str}.png")

# Plot definitions
plt.rc('axes', labelsize=22)
plt.rc('ytick', labelsize=18)
plt.rc('xtick', labelsize=18)
plt.rc('legend', fontsize=14)

# Colores para diferentes grupos
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
          '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

# Determinar número de paneles
if args.plot_type == 'both':
    num_panels = 2
else:
    num_panels = 1

# Crear figura
if num_panels == 2:
    figure, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 16), sharex=True)
else:
    figure, axis = plt.subplots(1, figsize=(14, 8))

# Plotear datos agrupados
for idx, (group_value, file_list) in enumerate(sorted(grouped_data.items())):
    color = colors[idx % len(colors)]

    # Para cada grupo, podría haber múltiples archivos
    # Vamos a plotear todos juntos con el mismo color pero diferente alpha
    for file_idx, (csv_file, params) in enumerate(file_list):
        data = read_velocity_data(csv_file)

        if data is None:
            continue

        # Extraer columnas
        x1 = data[:, 0]   # cycle
        vx = data[:, 4]   # vx (particle)
        vy = data[:, 5]   # vy (particle)
        vz = data[:, 6]   # vz (particle)
        vfx = data[:, 7]  # <vfx> (fluid)
        vfy = data[:, 8]  # <vfy> (fluid)
        vfz = data[:, 9]  # <vfz> (fluid)

        # Crear etiqueta
        if len(file_list) > 1:
            # Si hay múltiples archivos, agregar info adicional
            # Excluir los parámetros usados para agrupar
            other_params = {k: v for k, v in params.items() if k not in args.group_by}
            extra_info = ", ".join([f"{k}={v}" for k, v in sorted(other_params.items())[:2]])
            label_base = f"{group_value} ({extra_info})"
        else:
            label_base = group_value

        # Plotear según tipo
        if args.plot_type in ['particle', 'both']:
            target_ax = ax1 if args.plot_type == 'both' else axis

            if args.component in ['all', 'x']:
                target_ax.plot(x1, vx, color=color, linewidth=2.0,
                             label=f'{label_base} - Vx', alpha=0.8)
            if args.component in ['all', 'y']:
                target_ax.plot(x1, vy, color=color, linewidth=2.0, linestyle='--',
                             label=f'{label_base} - Vy', alpha=0.8)
            if args.component in ['all', 'z']:
                target_ax.plot(x1, vz, color=color, linewidth=2.0, linestyle='-.',
                             label=f'{label_base} - Vz', alpha=0.8)

        if args.plot_type in ['fluid', 'both']:
            target_ax = ax2 if args.plot_type == 'both' else axis

            if args.component in ['all', 'x']:
                target_ax.plot(x1, vfx, color=color, linewidth=2.0,
                             label=f'{label_base} - Vfx', alpha=0.8)
            if args.component in ['all', 'y']:
                target_ax.plot(x1, vfy, color=color, linewidth=2.0, linestyle='--',
                             label=f'{label_base} - Vfy', alpha=0.8)
            if args.component in ['all', 'z']:
                target_ax.plot(x1, vfz, color=color, linewidth=2.0, linestyle='-.',
                             label=f'{label_base} - Vfz', alpha=0.8)

# Configurar ejes y límites
if args.plot_type == 'both':
    # Configurar ambos paneles
    axes_list = [ax1, ax2]
    titles = ['Particle Velocity', 'Fluid Velocity (Average)']
    ylabels = ['Particle Velocity', 'Fluid Velocity']
else:
    axes_list = [axis]
    if args.plot_type == 'particle':
        titles = ['Particle Velocity']
        ylabels = ['Particle Velocity']
    else:
        titles = ['Fluid Velocity (Average)']
        ylabels = ['Fluid Velocity']

for ax, title, ylabel in zip(axes_list, titles, ylabels):
    if len(ax.get_lines()) > 0:
        x_data = ax.get_lines()[0].get_xdata()
        x_min_set = x_data.min() - (x_data.max() - x_data.min()) / 50
        x_max_set = x_data.max()
        limites_y(ax, x_min_set, x_max_set)

    ax.set_ylabel(ylabel, fontsize=22)
    ax.set_xlabel('Step', fontsize=22)
    params_str = " + ".join(args.group_by)
    ax.set_title(f'{title} - Grouped by {params_str}',
                 fontsize=24, fontweight='bold', pad=15)
    ax.legend(loc='best', ncol=1, fontsize=12)
    ax.grid(True, which='both', linewidth=1.5, linestyle='--', alpha=0.5)
    ax.tick_params(width=2)
    ax.tick_params(axis='x', labelbottom=True)
    for spine in ['top', 'bottom', 'left', 'right']:
        ax.spines[spine].set_linewidth(2)

plt.tight_layout()
plt.savefig(output_file, dpi=300, bbox_inches='tight')
print(f"\nPlot guardado en: {output_file}")
