#!/usr/bin/env python3
"""
plot_qsi_subgrid.py

Grafica la distribución de carga del fluido (nodos qsi) junto con los
subpuntos generados por subdivisión subgrid, comparando ambos con la
curva teórica de Debye-Hückel en función de la distancia a la partícula.

Muestra:
  - Puntos de nodo (ρ_net por voxel) en gris
  - Subpuntos (q_sub / volumen_sub) en color por voxel o por shell
  - Curva teórica DH: ρ_net(r) = -κ²·q/(4π·r)·exp(-κr)

Uso:
    python plot_qsi_subgrid.py -n 1000 -L 32 -q 1.0 -eps 1e4 \\
        --rho-el 8e-3 --rcut 3.0 --nsub 3

    python plot_qsi_subgrid.py -n 1000 -L 32 -q 1.0 -eps 1e4 \\
        --rho-el 8e-3 --rcut 3.0 --nsub 5 --quadratic-full \\
        --rmin-sub 0.5 --output subgrid_plot
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import argparse
import os
import sys

font_size = 14


# ---------------------------------------------------------------------------
# Lector de archivos qsi
# ---------------------------------------------------------------------------

class QsiReader:
    def __init__(self, filename, grid_size, num_species=2):
        self.filename = filename
        self.nx, self.ny, self.nz = grid_size
        self.num_species = num_species
        self.species = []

    def read(self):
        data = np.loadtxt(self.filename)
        n_points = self.nx * self.ny * self.nz
        if data.ndim == 1:
            data = data.reshape(-1, 1)
            self.num_species = 1
        elif data.shape[0] != n_points:
            raise ValueError(f"El archivo tiene {data.shape[0]} filas, "
                             f"pero se esperaban {n_points}")
        if data.shape[1] != self.num_species:
            self.num_species = data.shape[1]
        self.species = []
        for i in range(self.num_species):
            self.species.append(
                data[:, i].reshape((self.nx, self.ny, self.nz)))
        return self.species

    def net_charge_3d(self):
        if len(self.species) < 2:
            return self.species[0]
        return self.species[0] - self.species[1]


class ColloidReader:
    @staticmethod
    def read_colloid_position_cds(filename):
        with open(filename, 'r') as f:
            lines = f.readlines()
        pos_line = lines[35].strip()
        parts = pos_line.split()
        if len(parts) < 3:
            raise ValueError(f"Línea 36 de {filename} no tiene 3 valores: '{pos_line}'")
        return np.array([float(parts[0]), float(parts[1]), float(parts[2])])


# ---------------------------------------------------------------------------
# Funciones de interpolación (copiadas de qsi_subgrid_force.py)
# ---------------------------------------------------------------------------

def _rho_net_at(species_3d_list, i, j, k, nx, ny, nz):
    ii, jj, kk = i % nx, j % ny, k % nz
    if len(species_3d_list) >= 2:
        return species_3d_list[0][ii, jj, kk] - species_3d_list[1][ii, jj, kk]
    return species_3d_list[0][ii, jj, kk]


def _quad_weight_1d(t, rm1, r0, rp1):
    return rm1 * t * (t - 1.0) / 2.0 - r0 * (t + 1.0) * (t - 1.0) + rp1 * t * (t + 1.0) / 2.0


def _quartic_weight_1d(t, rm2, rm1, r0, rp1, rp2):
    L_m2 = (t + 1.0) * t * (t - 1.0) * (t - 2.0) / 24.0
    L_m1 = (t + 2.0) * t * (t - 1.0) * (t - 2.0) / (-6.0)
    L_0  = (t + 2.0) * (t + 1.0) * (t - 1.0) * (t - 2.0) / 4.0
    L_p1 = (t + 2.0) * (t + 1.0) * t * (t - 2.0) / (-6.0)
    L_p2 = (t + 2.0) * (t + 1.0) * t * (t - 1.0) / 24.0
    return rm2 * L_m2 + rm1 * L_m1 + r0 * L_0 + rp1 * L_p1 + rp2 * L_p2


def _lagrange_bases_1d(t, nodes):
    n = len(nodes)
    L = np.ones(n)
    for k in range(n):
        for j in range(n):
            if j != k:
                L[k] *= (t - nodes[j]) / (nodes[k] - nodes[j])
    return L


def _tensor_lagrange_weights(offsets, nodes_1d, rho_cube):
    nodes = np.array(nodes_1d, dtype=float)
    n_sub = len(offsets)
    Lx = np.array([_lagrange_bases_1d(offsets[s, 0], nodes) for s in range(n_sub)])
    Ly = np.array([_lagrange_bases_1d(offsets[s, 1], nodes) for s in range(n_sub)])
    Lz = np.array([_lagrange_bases_1d(offsets[s, 2], nodes) for s in range(n_sub)])
    return np.einsum('sa,sb,sc,abc->s', Lx, Ly, Lz, rho_cube)


def subgrid_offsets(nsub):
    step = 1.0 / (nsub + 1)
    positions = np.arange(1, nsub + 1) * step - 0.5
    return np.array([[dx, dy, dz]
                     for dx in positions
                     for dy in positions
                     for dz in positions])


def compute_q_subs(species_3d_list, i, j, k, nx, ny, nz, rho_net, offsets,
                   gradient_dist, quadratic_dist, quartic_dist,
                   quadratic_full, quartic_full):
    """Devuelve array q_subs (sin redistribución por exclusión)."""
    n_sub = len(offsets)
    q_vox = rho_net

    if quartic_full:
        nodes = [-2, -1, 0, 1, 2]
        rho_cube = np.array([[[
            _rho_net_at(species_3d_list, i+a, j+b, k+c, nx, ny, nz)
            for c in nodes] for b in nodes] for a in nodes])
        w_raw = _tensor_lagrange_weights(offsets, nodes, rho_cube)
    elif quadratic_full:
        nodes = [-1, 0, 1]
        rho_cube = np.array([[[
            _rho_net_at(species_3d_list, i+a, j+b, k+c, nx, ny, nz)
            for c in nodes] for b in nodes] for a in nodes])
        w_raw = _tensor_lagrange_weights(offsets, nodes, rho_cube)
    elif quartic_dist:
        rxm2 = _rho_net_at(species_3d_list, i-2, j,   k,   nx, ny, nz)
        rxm1 = _rho_net_at(species_3d_list, i-1, j,   k,   nx, ny, nz)
        rxp1 = _rho_net_at(species_3d_list, i+1, j,   k,   nx, ny, nz)
        rxp2 = _rho_net_at(species_3d_list, i+2, j,   k,   nx, ny, nz)
        rym2 = _rho_net_at(species_3d_list, i,   j-2, k,   nx, ny, nz)
        rym1 = _rho_net_at(species_3d_list, i,   j-1, k,   nx, ny, nz)
        ryp1 = _rho_net_at(species_3d_list, i,   j+1, k,   nx, ny, nz)
        ryp2 = _rho_net_at(species_3d_list, i,   j+2, k,   nx, ny, nz)
        rzm2 = _rho_net_at(species_3d_list, i,   j,   k-2, nx, ny, nz)
        rzm1 = _rho_net_at(species_3d_list, i,   j,   k-1, nx, ny, nz)
        rzp1 = _rho_net_at(species_3d_list, i,   j,   k+1, nx, ny, nz)
        rzp2 = _rho_net_at(species_3d_list, i,   j,   k+2, nx, ny, nz)
        wx = _quartic_weight_1d(offsets[:, 0], rxm2, rxm1, rho_net, rxp1, rxp2)
        wy = _quartic_weight_1d(offsets[:, 1], rym2, rym1, rho_net, ryp1, ryp2)
        wz = _quartic_weight_1d(offsets[:, 2], rzm2, rzm1, rho_net, rzp1, rzp2)
        w_raw = wx * wy * wz
    elif quadratic_dist:
        rxm1 = _rho_net_at(species_3d_list, i-1, j,   k,   nx, ny, nz)
        rxp1 = _rho_net_at(species_3d_list, i+1, j,   k,   nx, ny, nz)
        rym1 = _rho_net_at(species_3d_list, i,   j-1, k,   nx, ny, nz)
        ryp1 = _rho_net_at(species_3d_list, i,   j+1, k,   nx, ny, nz)
        rzm1 = _rho_net_at(species_3d_list, i,   j,   k-1, nx, ny, nz)
        rzp1 = _rho_net_at(species_3d_list, i,   j,   k+1, nx, ny, nz)
        wx = _quad_weight_1d(offsets[:, 0], rxm1, rho_net, rxp1)
        wy = _quad_weight_1d(offsets[:, 1], rym1, rho_net, ryp1)
        wz = _quad_weight_1d(offsets[:, 2], rzm1, rho_net, rzp1)
        w_raw = wx * wy * wz
    elif gradient_dist:
        gx = (_rho_net_at(species_3d_list, i+1, j,   k,   nx, ny, nz)
            - _rho_net_at(species_3d_list, i-1, j,   k,   nx, ny, nz)) / 2.0
        gy = (_rho_net_at(species_3d_list, i,   j+1, k,   nx, ny, nz)
            - _rho_net_at(species_3d_list, i,   j-1, k,   nx, ny, nz)) / 2.0
        gz = (_rho_net_at(species_3d_list, i,   j,   k+1, nx, ny, nz)
            - _rho_net_at(species_3d_list, i,   j,   k-1, nx, ny, nz)) / 2.0
        q_base = q_vox / n_sub
        return q_base * (1.0 + gx * offsets[:, 0] + gy * offsets[:, 1] + gz * offsets[:, 2])
    else:
        return np.full(n_sub, q_vox / n_sub)

    # Para los métodos de interpolación: normalizar
    w_sum = w_raw.sum()
    if w_sum == 0.0:
        return np.full(n_sub, q_vox / n_sub)
    return q_vox * w_raw / w_sum


# ---------------------------------------------------------------------------
# Extracción de puntos: nodos y subpuntos
# ---------------------------------------------------------------------------

def extract_points(species_3d_list, particle_pos, grid_size,
                   rcut, rmin, nsub, rmin_sub, rcut_sub,
                   gradient_dist, quadratic_dist, quartic_dist,
                   quadratic_full, quartic_full):
    """
    Devuelve:
      node_dists, node_rho   : arrays 1D con distancia y ρ_net de cada nodo
      sub_dists, sub_rho_eff : arrays 1D con distancia y densidad efectiva de
                               cada subpunto activo (q_sub / vol_sub, donde
                               vol_sub = 1/n_sub lu³)
      sub_voxel_dist         : distancia del voxel padre para cada subpunto
                               (para colorear)
    """
    nx, ny, nz = grid_size
    px, py, pz = particle_pos
    offsets = subgrid_offsets(nsub)
    n_sub = len(offsets)
    vol_sub = 1.0 / n_sub  # volumen de cada subpartícula en lu³

    node_dists = []
    node_rho   = []
    sub_dists  = []
    sub_rho_eff = []
    sub_voxel_dist = []

    for i in range(nx):
        xi = i + 1
        for j in range(ny):
            yj = j + 1
            for k in range(nz):
                zk = k + 1
                r_vec = np.array([xi - px, yj - py, zk - pz])
                dist = np.linalg.norm(r_vec)

                if dist > rcut or dist < rmin:
                    continue

                if len(species_3d_list) >= 2:
                    rho_net = (species_3d_list[0][i, j, k]
                               - species_3d_list[1][i, j, k])
                else:
                    rho_net = species_3d_list[0][i, j, k]

                node_dists.append(dist)
                node_rho.append(rho_net)

                # Subpuntos
                q_subs = compute_q_subs(
                    species_3d_list, i, j, k, nx, ny, nz, rho_net, offsets,
                    gradient_dist, quadratic_dist, quartic_dist,
                    quadratic_full, quartic_full)

                sub_positions = np.array([xi, yj, zk]) + offsets
                r_to_part = np.array([px, py, pz]) - sub_positions
                dists_sub = np.linalg.norm(r_to_part, axis=1)

                active = (dists_sub >= rmin_sub) & (dists_sub > 1e-10)
                if rcut_sub:
                    active &= (dists_sub <= rcut)

                n_active = np.sum(active)
                if n_active == 0:
                    continue

                # Redistribuir carga de excluidas
                q_excluida = np.sum(q_subs[~active])
                q_subs_active = q_subs.copy()
                q_subs_active[~active] = 0.0
                if q_excluida != 0.0:
                    q_subs_active[active] += q_excluida / n_active

                for s_idx in np.where(active)[0]:
                    sub_dists.append(dists_sub[s_idx])
                    # q_sub: carga asignada al subpunto (mismas unidades que rho_net del nodo,
                    # ya que rho_net es densidad×1lu³ = carga)
                    sub_rho_eff.append(q_subs_active[s_idx])
                    sub_voxel_dist.append(dist)

    return (np.array(node_dists), np.array(node_rho),
            np.array(sub_dists), np.array(sub_rho_eff),
            np.array(sub_voxel_dist))


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

def plot_subgrid(node_dists, node_rho, sub_dists, sub_rho_eff,
                 sub_voxel_dist, kappa, charge, rmin, rmax,
                 nsub, dist_mode,
                 show_nodes=True, show_subs=True, show_dh=True,
                 log_y=False, output_file='qsi_subgrid.png'):

    fig, ax = plt.subplots(figsize=(10, 6))

    # Rango de la curva DH: desde la distancia mínima real de los datos
    all_dists = np.concatenate([node_dists, sub_dists]) if (len(node_dists) > 0 and len(sub_dists) > 0) \
                else (node_dists if len(node_dists) > 0 else sub_dists)
    r_start = max(np.min(all_dists) * 0.95, 0.1) if len(all_dists) > 0 else max(rmin, 0.1)
    r_theory = np.linspace(r_start, rmax, 500)

    # Curva DH: ρ_net(r) × 1lu³ para comparar con q_sub y rho_net de nodos
    if show_dh and kappa > 0:
        rho_dh = -kappa**2 * charge / (4.0 * np.pi * r_theory) * np.exp(-kappa * r_theory)
        ax.plot(r_theory, rho_dh, 'r-', linewidth=2,
                label=f'Teoría DH (κ={kappa:.4f}, λ_D={1/kappa:.2f})')

    # Nodos
    if show_nodes and len(node_dists) > 0:
        idx = np.argsort(node_dists)
        ax.scatter(node_dists[idx], node_rho[idx],
                   s=15, color='gray', alpha=0.5, zorder=2,
                   label=f'Nodos qsi ({len(node_dists)} pts)')

    # Subpuntos coloreados por distancia del voxel padre
    if show_subs and len(sub_dists) > 0:
        sc = ax.scatter(sub_dists, sub_rho_eff,
                        c=sub_voxel_dist, cmap='viridis',
                        s=6, alpha=0.4, zorder=3,
                        label=f'Subpuntos nsub={nsub} ({len(sub_dists)} pts)')
        cbar = plt.colorbar(sc, ax=ax)
        cbar.set_label('Distancia del voxel padre (lu)', fontsize=font_size - 2)

    if kappa > 0:
        ax.axvline(x=1.0 / kappa, color='purple', linestyle=':', linewidth=1.5,
                   alpha=0.7, label=f'λ_D = {1/kappa:.2f}')

    ax.set_xlabel('r (lu)', fontsize=font_size)
    ax.set_ylabel('q (lu⁻³·lu³ = carga)', fontsize=font_size)
    ax.set_title(
        f'Distribución de carga: nodos vs subpuntos vs DH\n'
        f'nsub={nsub} ({nsub**3} subpts/voxel), distribución={dist_mode}',
        fontsize=font_size)
    ax.legend(fontsize=font_size - 3, loc='upper right')
    ax.grid(True, alpha=0.3)

    if log_y:
        # Graficar |ρ| en escala log
        ax.set_yscale('log')
        ax.set_ylabel('|ρ_net| (lu⁻³)', fontsize=font_size)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nGráfico guardado en: {output_file}")
    plt.close()


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Grafica nodos qsi y subpuntos subgrid vs. teoría DH',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Ejemplos:
    python plot_qsi_subgrid.py -n 1000 -L 32 -q 1.0 -eps 1e4 \\
        --rho-el 8e-3 --rcut 3.0 --nsub 3

    python plot_qsi_subgrid.py -n 1000 -L 32 -q 1.0 -eps 1e4 \\
        --rho-el 8e-3 --rcut 3.0 --nsub 5 --quadratic-full --rmin-sub 0.5
""")

    # Parámetros de simulación
    parser.add_argument('-n', '--nstep', type=int, required=True)
    parser.add_argument('-L', '--grid-size', type=int, required=True)
    parser.add_argument('-q', '--charge', type=float, default=1.0)
    parser.add_argument('-eps', '--epsilon', type=float, default=1.0)
    parser.add_argument('--kappa', type=float, default=None,
                        help='Parámetro de Debye kappa (si no se da, se calcula desde --rho-el)')
    parser.add_argument('--rho-el', type=float, default=None,
                        help='Densidad de electrolito para calcular kappa')
    parser.add_argument('-kt', '--kt', type=float, default=1.0)
    parser.add_argument('--num-species', type=int, default=2)
    parser.add_argument('-d', '--directory', default='.')
    parser.add_argument('--qsi-file', default=None)
    parser.add_argument('--cds-file', default=None)

    # Radios
    parser.add_argument('--rcut', type=float, required=True,
                        help='Radio de corte para voxels')
    parser.add_argument('--rmin', type=float, default=0.5,
                        help='Radio mínimo para nodos (default: 0.5)')
    parser.add_argument('--rmax', type=float, default=None,
                        help='Radio máximo para el plot (default: rcut)')
    parser.add_argument('--rmin-sub', type=float, default=0.0,
                        help='Radio mínimo de exclusión para subpuntos (default: 0.0)')
    parser.add_argument('--rcut-sub', action='store_true',
                        help='Aplicar rcut también a subpuntos individuales')

    # Subdivisión
    parser.add_argument('--nsub', type=int, default=3)

    # Distribución de carga
    parser.add_argument('--gradient-dist', action='store_true')
    parser.add_argument('--quadratic-dist', action='store_true')
    parser.add_argument('--quartic-dist', action='store_true')
    parser.add_argument('--quadratic-full', action='store_true')
    parser.add_argument('--quartic-full', action='store_true')

    # Visualización
    parser.add_argument('--no-nodes', action='store_true',
                        help='No mostrar los nodos qsi')
    parser.add_argument('--no-subs', action='store_true',
                        help='No mostrar los subpuntos')
    parser.add_argument('--no-dh', action='store_true',
                        help='No mostrar la curva DH')
    parser.add_argument('--log-y', action='store_true',
                        help='Escala logarítmica en Y (usa |ρ|)')
    parser.add_argument('--output', default='qsi_subgrid',
                        help='Nombre base del archivo de salida (default: qsi_subgrid)')

    args = parser.parse_args()

    # Calcular kappa
    if args.kappa is not None:
        kappa = args.kappa
    elif args.rho_el is not None:
        lb = 1.0 / (4.0 * np.pi * args.epsilon * args.kt)
        V = float(args.grid_size ** 3)
        rho_counterions = abs(args.charge) / V
        rho_total = args.rho_el + rho_counterions
        kappa = np.sqrt(4.0 * np.pi * lb * rho_total * 2.0)
        print(f"Parámetros DH: λ_B={lb:.4e}, κ={kappa:.6e}, λ_D={1/kappa:.4f} lu")
    else:
        kappa = 0.0
        print("AVISO: kappa=0, sin curva DH.")

    grid_size = (args.grid_size, args.grid_size, args.grid_size)
    rmax = args.rmax if args.rmax is not None else args.rcut

    # Archivos de entrada
    base_dir = args.directory
    nstep_str9 = f"{args.nstep:09d}"
    nstep_str8 = f"{args.nstep:08d}"
    qsi_file = args.qsi_file or os.path.join(base_dir, f"qsi-{nstep_str9}.001-001")
    cds_file = args.cds_file or os.path.join(base_dir, f"config.cds{nstep_str8}.001-001")

    for fname in [qsi_file, cds_file]:
        if not os.path.exists(fname):
            print(f"ERROR: No se encontró: {fname}")
            sys.exit(1)

    print(f"Leyendo posición de la partícula desde: {cds_file}")
    particle_pos = ColloidReader.read_colloid_position_cds(cds_file)
    print(f"  Posición: {particle_pos}")

    print(f"Leyendo densidades de carga desde: {qsi_file}")
    reader = QsiReader(qsi_file, grid_size, num_species=args.num_species)
    species_list = reader.read()
    print(f"  Especies: {reader.num_species}")

    # Modo de distribución
    if args.quartic_full:
        dist_mode = 'cuártica tensorial (5x5x5)'
    elif args.quadratic_full:
        dist_mode = 'cuadrática tensorial (3x3x3)'
    elif args.quartic_dist:
        dist_mode = 'cuártica separable'
    elif args.quadratic_dist:
        dist_mode = 'cuadrática separable'
    elif args.gradient_dist:
        dist_mode = 'gradiente lineal'
    else:
        dist_mode = 'uniforme'

    print(f"\nExtrayendo puntos (rcut={args.rcut}, nsub={args.nsub}, "
          f"rmin_sub={args.rmin_sub}, dist={dist_mode})...")

    node_dists, node_rho, sub_dists, sub_rho_eff, sub_voxel_dist = extract_points(
        species_list, particle_pos, grid_size,
        rcut=args.rcut, rmin=args.rmin, nsub=args.nsub,
        rmin_sub=args.rmin_sub, rcut_sub=args.rcut_sub,
        gradient_dist=args.gradient_dist,
        quadratic_dist=args.quadratic_dist,
        quartic_dist=args.quartic_dist,
        quadratic_full=args.quadratic_full,
        quartic_full=args.quartic_full)

    print(f"  Nodos: {len(node_dists)}")
    print(f"  Subpuntos activos: {len(sub_dists)}")

    output_file = f"{args.output}-n_{args.nstep:09d}.png"
    plot_subgrid(
        node_dists, node_rho, sub_dists, sub_rho_eff, sub_voxel_dist,
        kappa=kappa, charge=args.charge,
        rmin=args.rmin, rmax=rmax,
        nsub=args.nsub, dist_mode=dist_mode,
        show_nodes=not args.no_nodes,
        show_subs=not args.no_subs,
        show_dh=not args.no_dh,
        log_y=args.log_y,
        output_file=output_file)


if __name__ == '__main__':
    main()
