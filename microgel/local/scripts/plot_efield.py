#!/usr/bin/env python3
"""
Script para graficar el campo eléctrico sobre una línea y sobre un plano
leyendo directamente los archivos efield de Ludwig.

Soporta archivos: efield, efield_real, efield_fourier

Ejemplos de uso:

  # Campo total a lo largo del eje X en y=16.5, z=16.5
  python plot_efield.py -d DIR -n 10000 -L 32 --mode line \
      --start 0.5 16.5 16.5 --end 31.5 16.5 16.5

  # Campo en plano XY (z fijo) coloreado por magnitud
  python plot_efield.py -d DIR -n 10000 -L 32 --mode plane \
      --plane xy --position 16

  # Línea a lo largo de X con comparación teórica (Debye-Hückel)
  python plot_efield.py -d DIR -n 10000 -L 32 --mode line \
      --start 0.5 16.5 16.5 --end 31.5 16.5 16.5 \
      --charge-pos 16.1 16.5 16.5 --charge 1.0 --kappa 0.1 --epsilon 1e4 --kt 1e-5

  # Plano XZ con quiver (flechas)
  python plot_efield.py -d DIR -n 10000 -L 32 --mode plane \
      --plane xz --position 16 --quiver

  # Solo componente real de Ewald
  python plot_efield.py -d DIR -n 10000 -L 32 --mode line \
      --start 0.5 16.5 16.5 --end 31.5 16.5 16.5 --field-type real
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
import sys
from scipy.interpolate import RegularGridInterpolator

font_size = 13


# ─────────────────────────────────────────────────────────────────────────────
# Lectura de archivos
# ─────────────────────────────────────────────────────────────────────────────

class EFieldReader:
    """Lee archivos efield de Ludwig (3 columnas: Ex Ey Ez)"""

    def __init__(self, filename, grid_size):
        self.filename = filename
        self.nx, self.ny, self.nz = grid_size

    def read(self):
        data = np.loadtxt(self.filename)
        n = self.nx * self.ny * self.nz

        if data.ndim == 1:
            if len(data) != n * 3:
                raise ValueError(f"Archivo {self.filename}: {len(data)} valores, "
                                 f"esperados {n*3}")
            data = data.reshape((n, 3))
        elif data.shape != (n, 3):
            raise ValueError(f"Archivo {self.filename}: forma {data.shape}, "
                             f"esperada ({n}, 3)")

        Ex = data[:, 0].reshape((self.nx, self.ny, self.nz))
        Ey = data[:, 1].reshape((self.nx, self.ny, self.nz))
        Ez = data[:, 2].reshape((self.nx, self.ny, self.nz))
        return Ex, Ey, Ez


def read_colloid_pos(filename):
    """Lee posición del coloide desde colloids-*.csv"""
    with open(filename) as f:
        for line in f:
            if line.startswith('#') or 'id' in line.lower():
                continue
            parts = line.strip().split(',')
            if len(parts) >= 4:
                return np.array([float(parts[1]), float(parts[2]), float(parts[3])])
    raise ValueError(f"No se pudo leer posición del coloide en {filename}")


# ─────────────────────────────────────────────────────────────────────────────
# Campo teórico
# ─────────────────────────────────────────────────────────────────────────────

def theory_field_components(points, charge_pos, scale_factor, kappa=0.0):
    """
    Calcula las componentes Ex, Ey, Ez del campo teórico (Coulomb o Debye-Hückel)
    en un array de puntos (N, 3).

    Retorna (Ex, Ey, Ez) arrays de forma (N,).
    """
    xc, yc, zc = charge_pos
    dx = points[:, 0] - xc
    dy = points[:, 1] - yc
    dz = points[:, 2] - zc

    r2 = dx**2 + dy**2 + dz**2
    r = np.sqrt(r2)
    r = np.where(r < 0.01, 0.01, r)
    r2 = r**2

    if kappa > 0:
        exp_f = np.exp(-kappa * r)
        mag = scale_factor * exp_f * (kappa / r + 1.0 / r2)
    else:
        mag = scale_factor / (r2 * r)   # = scale_factor / r³

    return mag * dx / r, mag * dy / r, mag * dz / r


def fit_scale_factor(Ex, Ey, Ez, charge_pos, grid_size, kappa=0.0):
    """Ajuste mínimos cuadrados del factor de escala entre simulación y patrón 1/r²"""
    nx, ny, nz = grid_size
    xc, yc, zc = charge_pos
    x = np.arange(0.5, nx)
    y = np.arange(0.5, ny)
    z = np.arange(0.5, nz)
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

    r2 = (X - xc)**2 + (Y - yc)**2 + (Z - zc)**2
    r = np.sqrt(np.where(r2 < 0.0001, 0.0001, r2))
    r2 = r**2

    if kappa > 0:
        exp_f = np.exp(-kappa * r)
        mag_unit = exp_f * (kappa / r + 1.0 / r2)
    else:
        mag_unit = 1.0 / (r2 * r)

    dx = (X - xc); dy = (Y - yc); dz = (Z - zc)
    Exu = mag_unit * dx / r
    Eyu = mag_unit * dy / r
    Ezu = mag_unit * dz / r

    mask = r2 > 4.0
    num = (np.sum(Ex[mask] * Exu[mask]) + np.sum(Ey[mask] * Eyu[mask]) +
           np.sum(Ez[mask] * Ezu[mask]))
    den = (np.sum(Exu[mask]**2) + np.sum(Eyu[mask]**2) + np.sum(Ezu[mask]**2))
    return num / den


# ─────────────────────────────────────────────────────────────────────────────
# Extracción sobre una línea
# ─────────────────────────────────────────────────────────────────────────────

def extract_line(Ex, Ey, Ez, start, end, grid_size, num_points=200, show_nodes=False):
    """
    Interpola (o extrae en nodos) Ex, Ey, Ez a lo largo de una línea.

    El eje de la línea se parametriza como coordenada física a lo largo de la
    dirección de la línea:
      - Si la línea es paralela a un eje (X, Y o Z), el eje de la gráfica es
        la coordenada física de ese eje (igual que la del archivo efield).
      - Si la línea es oblicua, el eje es la distancia acumulada desde start.

    Retorna (coord_along_line, Ex_line, Ey_line, Ez_line, points_3d).
    """
    nx, ny, nz = grid_size
    x = np.arange(0.5, nx)
    y = np.arange(0.5, ny)
    z = np.arange(0.5, nz)

    interp_Ex = RegularGridInterpolator((x, y, z), Ex, bounds_error=False, fill_value=0)
    interp_Ey = RegularGridInterpolator((x, y, z), Ey, bounds_error=False, fill_value=0)
    interp_Ez = RegularGridInterpolator((x, y, z), Ez, bounds_error=False, fill_value=0)

    start = np.array(start, dtype=float)
    end   = np.array(end,   dtype=float)
    line_vec = end - start

    # Detectar qué componente varía (para usar coordenada física en el eje)
    varying = np.where(np.abs(line_vec) > 1e-10)[0]
    if len(varying) == 1:
        # Línea paralela a un eje -> eje de gráfica = coordenada de ese eje
        axis = varying[0]
        use_coord = True
    else:
        use_coord = False

    if show_nodes:
        line_len2 = np.dot(line_vec, line_vec)
        tol = 0.01
        nodes = []
        coords = []
        for i, xi in enumerate(x):
            for j, yj in enumerate(y):
                for k, zk in enumerate(z):
                    node = np.array([xi, yj, zk])
                    t = np.dot(node - start, line_vec) / line_len2
                    if 0 <= t <= 1:
                        proj = start + t * line_vec
                        if np.sum((node - proj)**2) < tol**2:
                            nodes.append(node)
                            coords.append(node[axis] if use_coord
                                          else np.linalg.norm(proj - start))
        if not nodes:
            print("ADVERTENCIA: No se encontraron nodos en la línea.")
            return None, None, None, None, None
        pts = np.array(nodes)
        coord = np.array(coords)
        # Ordenar por coordenada
        order = np.argsort(coord)
        pts = pts[order]
        coord = coord[order]
    else:
        pts = np.array([np.linspace(start[i], end[i], num_points) for i in range(3)]).T
        if use_coord:
            coord = pts[:, axis]
        else:
            coord = np.linalg.norm(pts - start, axis=1)

    Ex_l = interp_Ex(pts)
    Ey_l = interp_Ey(pts)
    Ez_l = interp_Ez(pts)
    return coord, Ex_l, Ey_l, Ez_l, pts


# ─────────────────────────────────────────────────────────────────────────────
# Gráfico sobre línea
# ─────────────────────────────────────────────────────────────────────────────

def plot_line(Ex, Ey, Ez, start, end, grid_size, num_points=200,
              show_nodes=False, charge_pos=None, kappa=0.0, field_label='E',
              output=None, title_extra=''):
    """
    Genera una figura con 4 paneles (Ex, Ey, Ez, |E|) a lo largo de una línea.

    El eje X de la gráfica es la coordenada física del eje que varía (misma
    escala que el archivo efield), NO la distancia desde start. Esto garantiza
    que la posición del coloide en el eje X coincide con su coordenada real.

    Opcionalmente superpone la curva teórica si se proporciona charge_pos.
    """
    coord_s, Ex_l, Ey_l, Ez_l, pts = extract_line(
        Ex, Ey, Ez, start, end, grid_size, num_points, show_nodes)
    if coord_s is None:
        return

    E_mag = np.sqrt(Ex_l**2 + Ey_l**2 + Ez_l**2)

    # Determinar etiqueta del eje
    start_a = np.array(start, dtype=float)
    end_a   = np.array(end,   dtype=float)
    line_vec = end_a - start_a
    varying = np.where(np.abs(line_vec) > 1e-10)[0]
    if len(varying) == 1:
        xlabel = ['X [lu]', 'Y [lu]', 'Z [lu]'][varying[0]]
    else:
        xlabel = 'Distancia desde start [lu]'

    # Teórico
    has_theory = charge_pos is not None
    if has_theory:
        scale = fit_scale_factor(Ex, Ey, Ez, charge_pos, grid_size, kappa=kappa)
        print(f"  Factor de escala ajustado: {scale:.4e}")
        # Puntos de alta resolución para teoría, usando la misma parametrización
        pts_th = np.array([np.linspace(start_a[i], end_a[i], 500) for i in range(3)]).T
        if len(varying) == 1:
            coord_th = pts_th[:, varying[0]]
        else:
            coord_th = np.linalg.norm(pts_th - start_a, axis=1)
        Ex_th, Ey_th, Ez_th = theory_field_components(pts_th, charge_pos, scale, kappa)
        E_mag_th = np.sqrt(Ex_th**2 + Ey_th**2 + Ez_th**2)

    fig, axes = plt.subplots(2, 2, figsize=(14, 9))

    sim_style = dict(marker='o', linestyle='none', markersize=4, alpha=0.8) if show_nodes \
               else dict(linestyle='-', linewidth=1.8)

    data_pairs = [
        (axes[0, 0], Ex_l, 'Ex'),
        (axes[0, 1], Ey_l, 'Ey'),
        (axes[1, 0], Ez_l, 'Ez'),
        (axes[1, 1], E_mag, '|E|'),
    ]
    if has_theory:
        theory_data = [Ex_th, Ey_th, Ez_th, E_mag_th]

    for idx, (ax, vals, comp) in enumerate(data_pairs):
        ax.plot(coord_s, vals, color='steelblue', label=f'{field_label} simulado', **sim_style)
        if has_theory:
            ax.plot(coord_th, theory_data[idx], 'r--', linewidth=1.5,
                    label='Teórico (DH)' if kappa > 0 else 'Teórico (Coulomb)')
            # Marcar posición del coloide
            if len(varying) == 1:
                ax.axvline(charge_pos[varying[0]], color='green', linestyle=':',
                           linewidth=1.2, alpha=0.7, label=f'coloide x={charge_pos[varying[0]]:.2f}')
        ax.set_xlabel(xlabel, fontsize=font_size)
        ax.set_ylabel(comp, fontsize=font_size)
        ax.set_title(comp, fontsize=font_size)
        ax.legend(fontsize=font_size - 2)
        ax.grid(True, alpha=0.3)

    nx, ny, nz = grid_size
    suptitle = (f'Campo eléctrico {field_label} a lo largo de la línea\n'
                f'Desde ({start_a[0]:.1f},{start_a[1]:.1f},{start_a[2]:.1f}) '
                f'hasta ({end_a[0]:.1f},{end_a[1]:.1f},{end_a[2]:.1f})  |  '
                f'L=({nx},{ny},{nz}){title_extra}')
    plt.suptitle(suptitle, fontsize=font_size)
    plt.tight_layout()

    if output:
        plt.savefig(output, dpi=200, bbox_inches='tight')
        print(f"Guardado: {output}")
    else:
        plt.show()
    plt.close()


# ─────────────────────────────────────────────────────────────────────────────
# Gráfico sobre un plano
# ─────────────────────────────────────────────────────────────────────────────

def plot_plane(Ex, Ey, Ez, plane, position_coord, grid_size,
               quiver=False, quiver_step=2, field_label='E',
               charge_pos=None, output=None, title_extra=''):
    """
    Genera una figura con la magnitud |E| en un plano de corte.

    position_coord: coordenada física del plano de corte (misma escala que el
    archivo efield). Se busca el nodo más cercano. Ej: z=16.0 -> nodo z=15.5
    o z=16.5, el más próximo.

    Opcionalmente superpone flechas (quiver) y marca la posición del coloide.
    """
    nx, ny, nz = grid_size
    plane = plane.lower()

    # Coordenadas físicas de los nodos
    xn = np.arange(0.5, nx)
    yn = np.arange(0.5, ny)
    zn = np.arange(0.5, nz)

    if plane == 'xy':
        idx = int(np.argmin(np.abs(zn - position_coord)))
        coord_val = zn[idx]
        E_mag = np.sqrt(Ex[:, :, idx]**2 + Ey[:, :, idx]**2 + Ez[:, :, idx]**2)
        C1, C2 = Ex[:, :, idx], Ey[:, :, idx]
        ax1c, ax2c = xn, yn
        xl, yl = 'X [lu]', 'Y [lu]'
        plane_label = f'Z={coord_val:.1f} (índice {idx})'
        cp1 = charge_pos[0] if charge_pos is not None else None
        cp2 = charge_pos[1] if charge_pos is not None else None
    elif plane == 'xz':
        idx = int(np.argmin(np.abs(yn - position_coord)))
        coord_val = yn[idx]
        E_mag = np.sqrt(Ex[:, idx, :]**2 + Ey[:, idx, :]**2 + Ez[:, idx, :]**2)
        C1, C2 = Ex[:, idx, :], Ez[:, idx, :]
        ax1c, ax2c = xn, zn
        xl, yl = 'X [lu]', 'Z [lu]'
        plane_label = f'Y={coord_val:.1f} (índice {idx})'
        cp1 = charge_pos[0] if charge_pos is not None else None
        cp2 = charge_pos[2] if charge_pos is not None else None
    elif plane == 'yz':
        idx = int(np.argmin(np.abs(xn - position_coord)))
        coord_val = xn[idx]
        E_mag = np.sqrt(Ex[idx, :, :]**2 + Ey[idx, :, :]**2 + Ez[idx, :, :]**2)
        C1, C2 = Ey[idx, :, :], Ez[idx, :, :]
        ax1c, ax2c = yn, zn
        xl, yl = 'Y [lu]', 'Z [lu]'
        plane_label = f'X={coord_val:.1f} (índice {idx})'
        cp1 = charge_pos[1] if charge_pos is not None else None
        cp2 = charge_pos[2] if charge_pos is not None else None
    else:
        raise ValueError(f"Plano no reconocido: {plane}")

    A1, A2 = np.meshgrid(ax1c, ax2c, indexing='ij')

    fig, ax = plt.subplots(1, 1, figsize=(7, 6))

    im = ax.contourf(A1, A2, E_mag, levels=30, cmap='plasma')
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('|E| [lu]', fontsize=font_size)

    if quiver:
        s = quiver_step
        ax.quiver(A1[::s, ::s], A2[::s, ::s],
                  C1[::s, ::s], C2[::s, ::s],
                  color='white', alpha=0.6, scale=None, width=0.003)

    # Marcar posición del coloide en el plano
    if charge_pos is not None and cp1 is not None and cp2 is not None:
        ax.plot(cp1, cp2, 'w+', markersize=12, markeredgewidth=2, label='coloide')
        ax.legend(fontsize=font_size - 2, loc='upper right')

    ax.set_xlabel(xl, fontsize=font_size)
    ax.set_ylabel(yl, fontsize=font_size)
    ax.set_aspect('equal')
    ax.set_title(
        f'|{field_label}| en plano {plane.upper()}, {plane_label}\n'
        f'L=({nx},{ny},{nz}){title_extra}',
        fontsize=font_size)

    plt.tight_layout()

    if output:
        plt.savefig(output, dpi=200, bbox_inches='tight')
        print(f"Guardado: {output}")
    else:
        plt.show()
    plt.close()


# ─────────────────────────────────────────────────────────────────────────────
# main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Graficar campo eléctrico sobre línea o plano desde archivos efield de Ludwig',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)

    # Archivos / malla
    parser.add_argument('-d', '--directory', default='.',
                        help='Directorio con los archivos efield (default: .)')
    parser.add_argument('-n', '--nstep', type=int, required=True,
                        help='Número de paso temporal')
    parser.add_argument('-L', '--grid-size', type=int, nargs='+', required=True,
                        metavar='N',
                        help='Tamaño de malla: un valor (LxLxL) o tres (NX NY NZ)')
    parser.add_argument('--field-type', choices=['total', 'real', 'fourier'],
                        default='total',
                        help='Tipo de archivo efield a leer (default: total)')
    parser.add_argument('--colloid-file', default=None,
                        help='Archivo colloids-*.csv para leer posición del coloide')

    # Modo
    parser.add_argument('--mode', choices=['line', 'plane'], required=True,
                        help='Modo de visualización')

    # Línea
    parser.add_argument('--start', nargs=3, type=float, metavar=('X', 'Y', 'Z'),
                        help='Punto inicial de la línea')
    parser.add_argument('--end',   nargs=3, type=float, metavar=('X', 'Y', 'Z'),
                        help='Punto final de la línea')
    parser.add_argument('--num-points', type=int, default=200,
                        help='Puntos de interpolación en la línea (default: 200)')
    parser.add_argument('--show-nodes', action='store_true',
                        help='Mostrar solo nodos de la malla (sin interpolación)')

    # Plano
    parser.add_argument('--plane', choices=['xy', 'xz', 'yz'],
                        help='Plano de corte')
    parser.add_argument('--position', type=float, default=None,
                        help='Coordenada física del plano de corte (ej: 16.0). '
                             'Se selecciona el nodo más cercano. (default: centro)')
    parser.add_argument('--quiver', action='store_true',
                        help='Superponer flechas de campo en el plano')
    parser.add_argument('--quiver-step', type=int, default=2,
                        help='Paso de submuestreo para quiver (default: 2)')

    # Teoría (solo modo línea)
    parser.add_argument('--charge-pos', nargs=3, type=float, metavar=('X', 'Y', 'Z'),
                        help='Posición de la carga para comparación teórica')
    parser.add_argument('--charge', type=float, default=1.0,
                        help='Carga q (default: 1.0)')
    parser.add_argument('--kappa', type=float, default=0.0,
                        help='Parámetro de Debye κ (default: 0 = Coulomb)')
    parser.add_argument('--epsilon', type=float, default=1.0,
                        help='Permitividad ε (informativo, default: 1.0)')
    parser.add_argument('--kt', type=float, default=1.0,
                        help='Energía térmica kT (informativo, default: 1.0)')

    # Salida
    parser.add_argument('-o', '--output', default=None,
                        help='Archivo de salida (default: mostrar en pantalla)')

    args = parser.parse_args()

    # Tamaño de malla
    if len(args.grid_size) == 1:
        L = args.grid_size[0]
        grid_size = (L, L, L)
    elif len(args.grid_size) == 3:
        grid_size = tuple(args.grid_size)
    else:
        parser.error('--grid-size acepta 1 o 3 valores')

    # Nombre del archivo efield
    nstep_str = f"{args.nstep:09d}"
    prefix = {'total': 'efield', 'real': 'efield_real', 'fourier': 'efield_fourier'}[args.field_type]
    efield_file = os.path.join(args.directory, f"{prefix}-{nstep_str}.001-001")

    if not os.path.exists(efield_file):
        print(f"ERROR: No se encuentra {efield_file}")
        sys.exit(1)

    print(f"Leyendo {efield_file}...")
    reader = EFieldReader(efield_file, grid_size)
    Ex, Ey, Ez = reader.read()
    E_mag_full = np.sqrt(Ex**2 + Ey**2 + Ez**2)
    print(f"  Malla: {grid_size[0]}×{grid_size[1]}×{grid_size[2]}")
    print(f"  |E| rango: [{E_mag_full.min():.4e}, {E_mag_full.max():.4e}]")

    # Posición del coloide (opcional)
    charge_pos = args.charge_pos
    if charge_pos is None and args.colloid_file is not None:
        charge_pos = read_colloid_pos(args.colloid_file).tolist()
        print(f"  Posición del coloide: {charge_pos}")

    # Info extra para títulos
    title_extra = ''
    if args.field_type != 'total':
        title_extra += f'  [{args.field_type}]'
    if args.kappa > 0:
        title_extra += f'  κ={args.kappa:.3f} (λ_D={1/args.kappa:.2f})'
    if args.epsilon != 1.0:
        title_extra += f'  ε={args.epsilon:.2e}'
    if args.kt != 1.0:
        title_extra += f'  kT={args.kt:.2e}'

    field_label = {'total': 'E', 'real': 'E_real', 'fourier': 'E_fourier'}[args.field_type]

    if args.mode == 'line':
        if args.start is None or args.end is None:
            parser.error('--mode line requiere --start y --end')
        print(f"\nGraficando línea desde {args.start} hasta {args.end}...")
        plot_line(Ex, Ey, Ez,
                  args.start, args.end, grid_size,
                  num_points=args.num_points,
                  show_nodes=args.show_nodes,
                  charge_pos=charge_pos,
                  kappa=args.kappa,
                  field_label=field_label,
                  output=args.output,
                  title_extra=title_extra)

    else:  # plane
        if args.plane is None:
            parser.error('--mode plane requiere --plane')
        pos = args.position
        if pos is None:
            # Centro geométrico: coordenada del nodo central
            L = grid_size[{'xy': 2, 'xz': 1, 'yz': 0}[args.plane]]
            pos = L / 2.0   # coordenada física central
        print(f"\nGraficando plano {args.plane.upper()} en coordenada {pos}...")
        plot_plane(Ex, Ey, Ez,
                   args.plane, pos, grid_size,
                   quiver=args.quiver,
                   quiver_step=args.quiver_step,
                   field_label=field_label,
                   charge_pos=charge_pos,
                   output=args.output,
                   title_extra=title_extra)

    print("¡Listo!")


if __name__ == '__main__':
    main()
