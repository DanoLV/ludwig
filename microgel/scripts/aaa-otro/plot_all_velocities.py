#! /usr/bin/python3
"""
Script para combinar plots de velocidad de múltiples subcarpetas stencil
Lee los archivos datosfluid.csv de cada subcarpeta y genera un plot combinado

Uso:
    ./plot_all_velocities.py -d <directorio_padre> [--out_dir <directorio_salida>] [--fluid-only]

Ejemplo:
    ./plot_all_velocities.py -d /path/to/parent/folder --out_dir ./plots

El script busca automáticamente todas las subcarpetas que contengan 'stencil_' en su nombre
y combina los plots de velocidad de todas ellas en un solo gráfico.
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
import os
from pathlib import Path


def limites_y(axis, x_min, x_max):
    """Ajusta los límites Y automáticamente para el rango de X especificado"""
    axis.set_xlim(x_min, x_max)

    ys_in_range = []
    for line in axis.get_lines():
        xd = np.asarray(line.get_xdata())
        yd = np.asarray(line.get_ydata())
        # Seleccionar solo puntos con x dentro del rango y valores finitos en y
        mask = (xd >= x_min) & (xd <= x_max) & np.isfinite(yd)
        if np.any(mask):
            ys_in_range.append(yd[mask])

    if ys_in_range:
        all_y = np.hstack(ys_in_range)
        y_min, y_max = all_y.min(), all_y.max()
        if y_min == y_max:
            # Si todos los valores son iguales, expandimos un poco para visualizar
            pad = abs(y_min) * 0.01 if y_min != 0 else 1e-6
        else:
            pad = 0.05 * (y_max - y_min)   # 5% padding
        axis.set_ylim(y_min - pad, y_max + pad)

    return 0


def extract_stencil_number(dirname):
    """Extrae el número de stencil del nombre de directorio"""
    if 'stencil_' in dirname:
        try:
            return int(dirname.split('stencil_')[-1])
        except:
            return None
    return None


def find_stencil_subdirs(parent_dir):
    """Encuentra todas las subcarpetas que contienen 'stencil_' en su nombre"""
    parent_path = Path(parent_dir)
    stencil_dirs = []

    if not parent_path.exists():
        print(f"Error: El directorio {parent_dir} no existe")
        return []

    # Buscar directorios que contengan 'stencil_' en su nombre
    for item in parent_path.iterdir():
        if item.is_dir() and 'stencil_' in item.name:
            stencil_num = extract_stencil_number(item.name)
            if stencil_num is not None:
                stencil_dirs.append((stencil_num, item))

    # Ordenar por número de stencil
    stencil_dirs.sort(key=lambda x: x[0])

    return stencil_dirs


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


# Set command line parameters
argParser = argparse.ArgumentParser(
    description="Combina plots de velocidad de múltiples subcarpetas stencil"
)
argParser.add_argument(
    "-d", "--dir",
    required=True,
    help="Directorio padre que contiene las subcarpetas stencil"
)
argParser.add_argument(
    "--out_dir",
    help="Directorio de salida para los plots",
    default="."
)
argParser.add_argument(
    "--fluid-only",
    action="store_true",
    help="Plotear solo velocidad promedio del fluido"
)

# Parse command line arguments
try:
    args = argParser.parse_args()
except:
    sys.exit("No se pudieron leer los parámetros de línea de comandos")

# Directorio padre
parent_dir = args.dir
if not os.path.exists(parent_dir):
    sys.exit(f"El directorio {parent_dir} no existe")

# Buscar subcarpetas stencil
stencil_dirs = find_stencil_subdirs(parent_dir)
if not stencil_dirs:
    sys.exit(f"No se encontraron subcarpetas con 'stencil_' en {parent_dir}")

print(f"Encontradas {len(stencil_dirs)} subcarpetas stencil:")
for num, path in stencil_dirs:
    print(f"  - Stencil {num}: {path.name}")

# Archivo de salida
output_file = os.path.join(args.out_dir, "velocidades_combined.png")

# Plot definitions
plt.rc('axes', labelsize=22)
plt.rc('ytick', labelsize=18)
plt.rc('xtick', labelsize=18)
plt.rc('legend', fontsize=14)

# Colores y estilos para diferentes stencils
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
          '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
linestyles = ['-', '--', '-.', ':']

if args.fluid_only:
    #-------------------------------------------------------------
    # VELOCIDAD DEL FLUIDO SOLAMENTE - UN SOLO PANEL
    #-------------------------------------------------------------
    figure, axis = plt.subplots(1, figsize=(14, 8))

    for idx, (stencil_num, stencil_path) in enumerate(stencil_dirs):
        data_file = stencil_path / "proceced_data" / "datosfluid.csv"
        data = read_velocity_data(data_file)

        if data is None:
            print(f"Saltando stencil {stencil_num}: no se pudo leer el archivo de datos")
            continue

        # Extraer columnas de datos
        x1 = data[:, 0]  # cycle
        vfx = data[:, 7]  # <vfx>
        vfy = data[:, 8]  # <vfy>
        vfz = data[:, 9]  # <vfz>

        # Seleccionar color y estilo
        color = colors[idx % len(colors)]

        # Plot
        axis.plot(x1, vfx, color=color, linewidth=2.0,
                 label=f'Stencil {stencil_num} - Vfx', alpha=0.8)
        axis.plot(x1, vfy, color=color, linewidth=2.0,
                 linestyle='--', label=f'Stencil {stencil_num} - Vfy', alpha=0.8)
        axis.plot(x1, vfz, color=color, linewidth=2.0,
                 linestyle='-.', label=f'Stencil {stencil_num} - Vfz', alpha=0.8)

    # Ajustar límites
    if len(axis.get_lines()) > 0:
        x_data = axis.get_lines()[0].get_xdata()
        x_min_set = x_data.min() - (x_data.max() - x_data.min()) / 50
        x_max_set = x_data.max()
        limites_y(axis, x_min_set, x_max_set)

    axis.set_ylabel('Fluid Velocity', fontsize=22)
    axis.set_xlabel('Cycle', fontsize=22)
    axis.set_title('Fluid Velocity Comparison - All Stencils',
                   fontsize=24, fontweight='bold', pad=15)
    axis.legend(loc='best', ncol=2)
    axis.grid(True, which='both', linewidth=1.5, linestyle='--', alpha=0.5)
    axis.tick_params(width=2)
    for spine in ['top', 'bottom', 'left', 'right']:
        axis.spines[spine].set_linewidth(2)

else:
    #-------------------------------------------------------------
    # DOS PANELES: VELOCIDAD DE PARTÍCULA Y FLUIDO
    #-------------------------------------------------------------
    figure, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 16), sharex=True)

    for idx, (stencil_num, stencil_path) in enumerate(stencil_dirs):
        data_file = stencil_path / "proceced_data" / "datosfluid.csv"
        data = read_velocity_data(data_file)

        if data is None:
            print(f"Saltando stencil {stencil_num}: no se pudo leer el archivo de datos")
            continue

        # Extraer columnas de datos
        x1 = data[:, 0]   # cycle
        vx = data[:, 4]   # vx (particle)
        vy = data[:, 5]   # vy (particle)
        vz = data[:, 6]   # vz (particle)
        vfx = data[:, 7]  # <vfx> (fluid)
        vfy = data[:, 8]  # <vfy> (fluid)
        vfz = data[:, 9]  # <vfz> (fluid)

        # Seleccionar color
        color = colors[idx % len(colors)]

        #-------------------------------------------------------------
        # Panel 1: VELOCIDAD DE LA PARTÍCULA
        #-------------------------------------------------------------
        ax1.plot(x1, vx, color=color, linewidth=2.0,
                label=f'Stencil {stencil_num} - Vx', alpha=0.8)
        ax1.plot(x1, vy, color=color, linewidth=2.0, linestyle='--',
                label=f'Stencil {stencil_num} - Vy', alpha=0.8)
        ax1.plot(x1, vz, color=color, linewidth=2.0, linestyle='-.',
                label=f'Stencil {stencil_num} - Vz', alpha=0.8)

        #-------------------------------------------------------------
        # Panel 2: VELOCIDAD DEL FLUIDO
        #-------------------------------------------------------------
        ax2.plot(x1, vfx, color=color, linewidth=2.0,
                label=f'Stencil {stencil_num} - Vfx', alpha=0.8)
        ax2.plot(x1, vfy, color=color, linewidth=2.0, linestyle='--',
                label=f'Stencil {stencil_num} - Vfy', alpha=0.8)
        ax2.plot(x1, vfz, color=color, linewidth=2.0, linestyle='-.',
                label=f'Stencil {stencil_num} - Vfz', alpha=0.8)

    # Ajustar límites para ambos paneles
    if len(ax1.get_lines()) > 0:
        x_data = ax1.get_lines()[0].get_xdata()
        x_min_set = x_data.min() - (x_data.max() - x_data.min()) / 50
        x_max_set = x_data.max()

        limites_y(ax1, x_min_set, x_max_set)
        limites_y(ax2, x_min_set, x_max_set)

    # Configurar Panel 1
    ax1.set_ylabel('Particle Velocity', fontsize=22)
    ax1.set_title('Particle Velocity Comparison - All Stencils',
                  fontsize=24, fontweight='bold', pad=15)
    ax1.set_xlabel('Step', fontsize=22)
    ax1.legend(loc='best', ncol=2, fontsize=12)
    ax1.grid(True, which='both', linewidth=1.5, linestyle='--', alpha=0.5)
    ax1.tick_params(width=2)
    ax1.tick_params(axis='x', labelbottom=True)
    for spine in ['top', 'bottom', 'left', 'right']:
        ax1.spines[spine].set_linewidth(2)

    # Configurar Panel 2
    ax2.set_ylabel('Fluid Velocity', fontsize=22)
    ax2.set_xlabel('Step', fontsize=22)
    ax2.set_title('Fluid Velocity (Average) Comparison - All Stencils',
                  fontsize=24, fontweight='bold', pad=15)
    ax2.legend(loc='best', ncol=2, fontsize=12)
    ax2.grid(True, which='both', linewidth=1.5, linestyle='--', alpha=0.5)
    ax2.tick_params(width=2)
    for spine in ['top', 'bottom', 'left', 'right']:
        ax2.spines[spine].set_linewidth(2)

    # Ajustar espaciado entre subplots
    plt.tight_layout()

# Guardar figura
plt.savefig(output_file, dpi=300, bbox_inches='tight')
print(f"\nPlot guardado en: {output_file}")
