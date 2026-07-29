#!/usr/bin/env python3
"""
Script para graficar las componentes del campo E_sub sobre la partícula en
función de su posición, recorriendo las subcarpetas de un directorio padre.

Cada subcarpeta corresponde a una simulación con la partícula en una posición
distinta (codificada en el nombre como pos_X_Y_Z). De cada una se lee el
archivo proceced_data/particle_Esub.csv (separador ';', decimal ',') y se
extraen Esub_X, Esub_Y, Esub_Z (y Emod) para el step indicado.

Uso:
    python plot_Esub_vs_pos.py --step 1000 --output Esub_vs_pos.png
    python plot_Esub_vs_pos.py --step 500 -o salida.png --parent-dir /ruta/al/directorio
"""

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

font_size = 16

DEFAULT_PARENT = ('/home/bater/Sim/ludwig/microgel/local/'
                  'Ewald-gaussian_dual_fix-sigma_f_0.0-sigma_p_0.5-alpha_1.0-rc_6.0')

CSV_RELPATH = Path('proceced_data') / 'particle_Esub.csv'

_NUM = r'[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?'
POS_RE = re.compile(rf'pos_({_NUM})_({_NUM})_({_NUM})')


def parse_position(folder_name):
    """Extrae la posición (x, y, z) del nombre de la subcarpeta (pos_X_Y_Z)."""
    m = POS_RE.search(folder_name)
    if m is None:
        return None
    return np.array([float(m.group(1)), float(m.group(2)), float(m.group(3))])


def read_esub_at_step(csv_file, step):
    """
    Lee particle_Esub.csv y devuelve (Emod, Esub_X, Esub_Y, Esub_Z) para el step dado.

    Formato del archivo: '# Step;Index;Emod;Esub_X;Esub_Y;Esub_Z'
    con ';' como separador y ',' como separador decimal.

    Returns:
        numpy array [Emod, Ex, Ey, Ez], o None si el step no está en el archivo.
    """
    with open(csv_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split(';')
            if len(parts) < 6:
                continue
            try:
                row_step = int(parts[0])
            except ValueError:
                continue
            if row_step != step:
                continue
            vals = [float(p.replace(',', '.')) for p in parts[2:6]]
            return np.array(vals)
    return None


def collect_data(parent_dir, step):
    """
    Recorre las subcarpetas del directorio padre y junta posición + E_sub.

    Returns:
        positions: array (N, 3) con la posición de la partícula de cada subcarpeta
        fields:    array (N, 4) con [Emod, Esub_X, Esub_Y, Esub_Z]
        names:     lista de nombres de subcarpeta usados
    """
    parent = Path(parent_dir)
    if not parent.is_dir():
        print(f"ERROR: el directorio padre no existe: {parent}")
        sys.exit(1)

    positions = []
    fields = []
    names = []

    for sub in sorted(p for p in parent.iterdir() if p.is_dir()):
        pos = parse_position(sub.name)
        if pos is None:
            print(f"  Aviso: no se pudo extraer pos_X_Y_Z de '{sub.name}', se omite")
            continue

        csv_file = sub / CSV_RELPATH
        if not csv_file.is_file():
            print(f"  Aviso: no existe {csv_file}, se omite")
            continue

        e_vals = read_esub_at_step(csv_file, step)
        if e_vals is None:
            print(f"  Aviso: step {step} no encontrado en {csv_file}, se omite")
            continue

        positions.append(pos)
        fields.append(e_vals)
        names.append(sub.name)
        print(f"  {sub.name}: pos=({pos[0]:.3f},{pos[1]:.3f},{pos[2]:.3f}) "
              f"Esub=({e_vals[1]:.4e},{e_vals[2]:.4e},{e_vals[3]:.4e})")

    if not positions:
        print(f"ERROR: no se encontraron datos para el step {step} en {parent}")
        sys.exit(1)

    return np.array(positions), np.array(fields), names


def detect_varying_axis(positions):
    """Devuelve el índice (0,1,2) de la coordenada que más varía entre subcarpetas."""
    spans = positions.max(axis=0) - positions.min(axis=0)
    return int(np.argmax(spans))


def plot_esub(positions, fields, step, output_file, axis=None, show_emod=False):
    """
    Grafica Esub_X, Esub_Y, Esub_Z (y opcionalmente Emod) vs la posición de la
    partícula a lo largo del eje que varía.
    """
    axis_names = ['x', 'y', 'z']
    if axis is None:
        axis = detect_varying_axis(positions)
    coord = positions[:, axis]

    order = np.argsort(coord)
    coord = coord[order]
    fields = fields[order]

    fig, ax = plt.subplots(figsize=(10, 7))

    ax.plot(coord, fields[:, 1], 'o-', color='tab:blue',  label=r'$E_{sub,x}$')
    ax.plot(coord, fields[:, 2], 's-', color='tab:green', label=r'$E_{sub,y}$')
    ax.plot(coord, fields[:, 3], '^-', color='tab:red',   label=r'$E_{sub,z}$')
    if show_emod:
        ax.plot(coord, fields[:, 0], 'd--', color='black', label=r'$|E_{sub}|$')

    ax.axhline(0.0, color='gray', linewidth=0.8, linestyle=':')

    ax.set_xlabel(f'Position {axis_names[axis]}', fontsize=font_size)
    ax.set_ylabel('Electric field components', fontsize=font_size)
    ax.set_title(f'Electric field components vs position {axis_names[axis]} '
                 f'(step {step})', fontsize=font_size)
    ax.tick_params(labelsize=font_size - 2)
    ax.legend(fontsize=font_size - 2)
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(output_file, dpi=150)
    print(f"\nGráfico guardado: {output_file} ({len(coord)} puntos, eje {axis_names[axis]})")


def main():
    parser = argparse.ArgumentParser(
        description='Grafica E_sub (x,y,z) vs posición de la partícula para un step dado, '
                    'recorriendo las subcarpetas de un directorio padre.')
    parser.add_argument('--step', '-s', type=int, required=True,
                        help='Número de step a extraer de particle_Esub.csv')
    parser.add_argument('--output', '-o', type=str, required=True,
                        help='Nombre del archivo de salida (png)')
    parser.add_argument('--parent-dir', '-d', type=str, default=DEFAULT_PARENT,
                        help=f'Directorio padre con las subcarpetas (default: {DEFAULT_PARENT})')
    parser.add_argument('--axis', type=str, default=None, choices=['x', 'y', 'z'],
                        help='Eje de posición para el plot (default: el que más varía)')
    parser.add_argument('--show-emod', action='store_true',
                        help='Incluir también |E_sub| (Emod) en el gráfico')
    args = parser.parse_args()

    print(f"Directorio padre: {args.parent_dir}")
    print(f"Step: {args.step}\n")

    positions, fields, _ = collect_data(args.parent_dir, args.step)

    axis = {'x': 0, 'y': 1, 'z': 2}[args.axis] if args.axis else None
    plot_esub(positions, fields, args.step, args.output,
              axis=axis, show_emod=args.show_emod)


if __name__ == '__main__':
    main()
