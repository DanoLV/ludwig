#!/usr/bin/env python3
"""
Grafica planos de potencial químico μ_i centrados en la partícula.

El plano se define por su vector normal. El corte pasa por el nodo de la malla
más cercano a la posición de la partícula proyectado sobre ese eje.

Los archivos mui-NNNNNNNNN.001-001 tienen el mismo formato que qsi:
  N filas, 2 columnas (μ₊, μ₋), orden x-major de Ludwig.

La posición de la partícula se lee de config.cdsNNNNNNNN.001-001, línea 36.

Uso:
    python plot_mui_plane.py -n 530 -L 32 [-d directorio]
    python plot_mui_plane.py -n 530 -L 32 --normal 0,1,0   # plano XZ
    python plot_mui_plane.py -n 530 -L 32 --normal 1,0,0 --normal 0,1,0 --normal 0,0,1
    python plot_mui_plane.py -n 530 -L 32 --all-planes
    python plot_mui_plane.py --n-start 100 --n-end 500 --n-step 100 -L 32 --all-planes

Planos predefinidos por nombre de --normal: x/y/z o nx,ny,nz
"""

import numpy as np
import argparse
import sys
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm


# ---------------------------------------------------------------------------
# Lectores
# ---------------------------------------------------------------------------

class MuiReader:
    """Lee archivos mui de Ludwig: N filas, 2 columnas (μ₊, μ₋)"""

    def __init__(self, filename, grid_size):
        self.filename = filename
        self.nx, self.ny, self.nz = grid_size

    def read(self):
        data = np.loadtxt(self.filename)
        n_points = self.nx * self.ny * self.nz
        if data.shape[0] != n_points:
            raise ValueError(
                f"El archivo mui tiene {data.shape[0]} filas, "
                f"pero se esperaban {n_points}")
        mu0 = data[:, 0].reshape((self.nx, self.ny, self.nz))
        mu1 = data[:, 1].reshape((self.nx, self.ny, self.nz))
        return mu0, mu1


def read_particle_position(cds_file):
    """Lee posición desde config.cds*.001-001 (línea 36, 3 valores x y z)"""
    with open(cds_file, 'r') as f:
        lines = f.readlines()
    parts = lines[35].strip().split()
    if len(parts) < 3:
        raise ValueError(f"Línea 36 de {cds_file} no tiene 3 valores: '{lines[35].strip()}'")
    return np.array([float(parts[0]), float(parts[1]), float(parts[2])])


# ---------------------------------------------------------------------------
# Extracción del plano
# ---------------------------------------------------------------------------

def parse_normal(normal_str):
    """Convierte string 'x','y','z' o 'nx,ny,nz' a vector unitario numpy."""
    aliases = {
        'x': np.array([1., 0., 0.]),
        'y': np.array([0., 1., 0.]),
        'z': np.array([0., 0., 1.]),
    }
    s = normal_str.strip().lower()
    if s in aliases:
        v = aliases[s]
    else:
        try:
            parts = [float(x) for x in s.split(',')]
            if len(parts) != 3:
                raise ValueError
            v = np.array(parts)
        except ValueError:
            raise ValueError(f"Normal '{normal_str}' no reconocida. "
                             f"Use 'x','y','z' o 'nx,ny,nz'.")
    norm = np.linalg.norm(v)
    if norm < 1e-10:
        raise ValueError(f"Normal '{normal_str}' es el vector cero.")
    return v / norm


def normal_label(n):
    """Etiqueta legible para el vector normal."""
    for name, v in [('x', [1,0,0]), ('y', [0,1,0]), ('z', [0,0,1])]:
        if np.allclose(np.abs(n), v, atol=1e-6):
            return name
    return f"({n[0]:.2f},{n[1]:.2f},{n[2]:.2f})"


def extract_plane(field_3d, normal, particle_pos, grid_size):
    """
    Extrae un plano del campo 3D perpendicular a 'normal' que pasa
    por el nodo más cercano a particle_pos.

    Ludwig usa índices 1-based: el nodo (i,j,k) está en posición (i+1, j+1, k+1).
    particle_pos está en coordenadas 1-based.

    Devuelve:
      plane_data : array 2D con los valores del plano
      ax1_coords, ax2_coords : coordenadas de los ejes del plano (relativas a partícula)
      ax1_label, ax2_label   : nombres de los ejes
      slice_coord            : coordenada del corte a lo largo de normal (relativa a partícula)
    """
    nx, ny, nz = grid_size
    # Posición de la partícula en índices 0-based de la malla
    px, py, pz = particle_pos[0] - 1, particle_pos[1] - 1, particle_pos[2] - 1

    n = normal / np.linalg.norm(normal)
    # Eje dominante del vector normal
    dominant = int(np.argmax(np.abs(n)))

    if dominant == 0:  # normal ~ X → plano YZ
        slice_idx = int(round(px)) % nx
        plane_data = field_3d[slice_idx, :, :]           # (ny, nz)
        c1 = np.arange(ny) + 1 - (py + 1)               # relativo a partícula
        c2 = np.arange(nz) + 1 - (pz + 1)
        ax1_label, ax2_label = 'y', 'z'
        slice_coord = (slice_idx + 1) - (px + 1)
    elif dominant == 1:  # normal ~ Y → plano XZ
        slice_idx = int(round(py)) % ny
        plane_data = field_3d[:, slice_idx, :]           # (nx, nz)
        c1 = np.arange(nx) + 1 - (px + 1)
        c2 = np.arange(nz) + 1 - (pz + 1)
        ax1_label, ax2_label = 'x', 'z'
        slice_coord = (slice_idx + 1) - (py + 1)
    else:               # normal ~ Z → plano XY
        slice_idx = int(round(pz)) % nz
        plane_data = field_3d[:, :, slice_idx]           # (nx, ny)
        c1 = np.arange(nx) + 1 - (px + 1)
        c2 = np.arange(ny) + 1 - (py + 1)
        ax1_label, ax2_label = 'x', 'y'
        slice_coord = (slice_idx + 1) - (pz + 1)

    return plane_data, c1, c2, ax1_label, ax2_label, slice_coord


# ---------------------------------------------------------------------------
# Gráfico de un plano
# ---------------------------------------------------------------------------

def plot_plane(mu0_plane, mu1_plane, c1, c2, ax1_label, ax2_label,
               normal, slice_coord, nstep, particle_pos, output_file):
    """
    Genera figura con 2 subplots (μ₊ y μ₋) para un plano dado.
    Colormap divergente centrado en 0.
    """
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    planes  = [mu0_plane,   mu1_plane]
    titles  = [r'$\mu_+$ (cationes, $z=+1$)', r'$\mu_-$ (aniones, $z=-1$)']
    cmaps   = ['RdBu_r', 'RdBu_r']

    C1, C2 = np.meshgrid(c2, c1)   # meshgrid: eje horizontal=c2, vertical=c1

    for ax, plane, title, cmap in zip(axes, planes, titles, cmaps):
        valid = plane[np.isfinite(plane)]
        if len(valid) == 0:
            ax.set_title(f"{title}\n(sin datos)")
            continue
        vmax = np.nanmax(np.abs(valid))
        vmin = -vmax if vmax > 0 else -1.0
        norm = TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax) if vmax > 0 else None

        im = ax.pcolormesh(C2, C1, plane, cmap=cmap,
                           norm=norm, shading='auto')
        plt.colorbar(im, ax=ax, label=r'$\mu_i$')

        # Marca la posición de la partícula
        ax.axhline(0, color='k', linewidth=0.5, linestyle='--', alpha=0.5)
        ax.axvline(0, color='k', linewidth=0.5, linestyle='--', alpha=0.5)
        ax.plot(0, 0, 'k+', markersize=12, markeredgewidth=2,
                label=f'partícula ({particle_pos[0]:.1f},{particle_pos[1]:.1f},{particle_pos[2]:.1f})')

        ax.set_xlabel(f'{ax2_label} - {ax2_label}_p  (lu)', fontsize=11)
        ax.set_ylabel(f'{ax1_label} - {ax1_label}_p  (lu)', fontsize=11)
        nl = normal_label(normal)
        ax.set_title(f"{title}\nPlano n={nl}, corte={slice_coord:+.1f} lu, paso={nstep}",
                     fontsize=10)
        ax.set_aspect('equal')
        ax.legend(fontsize=8, loc='upper right')

    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Gráfico: {output_file}")


# ---------------------------------------------------------------------------
# Procesado de un paso
# ---------------------------------------------------------------------------

def process_step_plot(nstep, grid_size, base_dir, normals, output_prefix):
    """
    Lee mui y cds para el paso dado y genera los gráficos de plano.
    Devuelve False si los archivos no existen.
    """
    nstep_str9 = f"{nstep:09d}"
    nstep_str8 = f"{nstep:08d}"

    mui_file = os.path.join(base_dir, f"mui-{nstep_str9}.001-001")
    cds_file = os.path.join(base_dir, f"config.cds{nstep_str8}.001-001")

    if not os.path.exists(mui_file):
        print(f"  AVISO: paso {nstep} omitido — no se encontró {mui_file}")
        return False
    if not os.path.exists(cds_file):
        print(f"  AVISO: paso {nstep} omitido — no se encontró {cds_file}")
        return False

    mu0, mu1 = MuiReader(mui_file, grid_size).read()
    particle_pos = read_particle_position(cds_file)
    print(f"  Partícula: {particle_pos}")

    for normal in normals:
        mu0_plane, c1, c2, ax1_lbl, ax2_lbl, sc = extract_plane(
            mu0, normal, particle_pos, grid_size)
        mu1_plane, *_ = extract_plane(
            mu1, normal, particle_pos, grid_size)

        nl = normal_label(normal)
        out = os.path.join(base_dir,
                           f"{output_prefix}-n{nstep_str9}-plane_{nl}.png")
        plot_plane(mu0_plane, mu1_plane, c1, c2, ax1_lbl, ax2_lbl,
                   normal, sc, nstep, particle_pos, out)

    return True


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Grafica planos de μ_i centrados en la partícula',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Ejemplos:
  Plano XY (normal z), paso 530:
    python plot_mui_plane.py -n 530 -L 32

  Tres planos ortogonales, paso 530:
    python plot_mui_plane.py -n 530 -L 32 --all-planes

  Normal arbitraria:
    python plot_mui_plane.py -n 530 -L 32 --normal 1,1,0

  Serie temporal, todos los planos:
    python plot_mui_plane.py --n-start 100 --n-end 500 --n-step 100 -L 32 --all-planes
""")

    parser.add_argument('-n', '--nstep', type=int, default=None,
                        help='Paso de tiempo único')
    parser.add_argument('--n-start', type=int, default=None)
    parser.add_argument('--n-end',   type=int, default=None)
    parser.add_argument('--n-step',  type=int, default=None)

    parser.add_argument('-L', '--grid-size', type=int, required=True,
                        help='Tamaño de la malla cúbica L×L×L')
    parser.add_argument('-d', '--directory', default='.',
                        help='Directorio con los archivos (default: .)')
    parser.add_argument('--normal', action='append', default=None,
                        metavar='NX,NY,NZ',
                        help="Vector normal al plano. Se puede repetir. "
                             "Use 'x','y','z' o 'nx,ny,nz'. (default: z)")
    parser.add_argument('--all-planes', action='store_true', default=False,
                        help='Generar los tres planos ortogonales (x, y, z)')
    parser.add_argument('--output-prefix', default='mui_plane',
                        help='Prefijo del nombre de los archivos de salida (default: mui_plane)')

    args = parser.parse_args()

    # Lista de pasos
    if args.nstep is not None:
        steps = [args.nstep]
    elif args.n_start is not None and args.n_end is not None and args.n_step is not None:
        steps = list(range(args.n_start, args.n_end + 1, args.n_step))
    else:
        print("ERROR: Indica -n PASO o bien --n-start, --n-end y --n-step")
        sys.exit(1)

    # Lista de normales
    if args.all_planes:
        normals = [np.array([1.,0.,0.]), np.array([0.,1.,0.]), np.array([0.,0.,1.])]
    elif args.normal:
        try:
            normals = [parse_normal(s) for s in args.normal]
        except ValueError as e:
            print(f"ERROR: {e}")
            sys.exit(1)
    else:
        normals = [np.array([0., 0., 1.])]   # default: plano XY

    grid_size = (args.grid_size, args.grid_size, args.grid_size)

    normal_names = [normal_label(n) for n in normals]
    print(f"L={args.grid_size}, planos normales: {normal_names}")
    print(f"Pasos: {steps[0]} → {steps[-1]} ({len(steps)} pasos)\n")

    for nstep in steps:
        print(f"Paso {nstep:9d}:")
        process_step_plot(nstep, grid_size, args.directory, normals, args.output_prefix)

    print("\nListo.")


if __name__ == '__main__':
    main()
