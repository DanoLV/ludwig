#!/usr/bin/env python3
"""
qsi_subgrid_force.py

Calcula la fuerza electrostática que ejerce la distribución de carga del
fluido sobre la partícula coloidal, usando un refinamiento subgrid:

  - Cada voxel de la malla de Ludwig contiene una densidad de carga ρ_net.
  - La carga total del voxel (ρ_net × 1 lu³) se redistribuye uniformemente
    en n×n×n subpartículas equiespaciadas dentro del voxel.
  - Las subpartículas se ubican en posiciones interiores al voxel: la
    distancia desde el centro del voxel al borde es 0.5 lu, y los n puntos
    se colocan de modo que la distancia de cada subpartícula a la pared más
    cercana sea 1/(2*(n+1)) lu (división del arista en n+1 segmentos iguales).
  - Solo se procesan los voxels dentro del radio de corte (rcut) respecto a
    la partícula definida en config.cds.
  - La fuerza de Coulomb de cada subpartícula sobre la partícula es:
      F = q_colloid * q_sub / (4π·ε·r²) * r̂
    donde r̂ apunta desde la subpartícula hacia la partícula.
  - El resultado se guarda en un CSV y se imprime la fuerza total.

Uso:
    python qsi_subgrid_force.py -n 1000 -L 32 -q 1.0 -eps 1e4 \\
        --rcut 5.0 --nsub 3 --csv fuerza_subgrid

Formato del archivo qsi: qsi-NNNNNNNNN.001-001
Formato del archivo cds: config.cds<NNNNNNNNN>.001-001
"""

import numpy as np
import argparse
import os
import sys


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
            raise ValueError(
                f"El archivo tiene {data.shape[0]} filas, "
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


# ---------------------------------------------------------------------------
# Lector de posición del coloide
# ---------------------------------------------------------------------------

class ColloidReader:
    @staticmethod
    def read_colloid_position_cds(filename):
        with open(filename, 'r') as f:
            lines = f.readlines()
        pos_line = lines[35].strip()
        parts = pos_line.split()
        if len(parts) < 3:
            raise ValueError(
                f"Línea 36 de {filename} no tiene 3 valores: '{pos_line}'")
        return np.array([float(parts[0]), float(parts[1]), float(parts[2])])


# ---------------------------------------------------------------------------
# Cálculo subgrid
# ---------------------------------------------------------------------------

def subgrid_offsets(nsub):
    """
    Genera los nsub³ desplazamientos subgrid dentro de un voxel unitario.
    Los n puntos se ubican a distancias 1/(n+1), 2/(n+1), ..., n/(n+1)
    a lo largo de cada arista, centrado en 0 (rango [-0.5, +0.5]).
    """
    step = 1.0 / (nsub + 1)
    positions = np.arange(1, nsub + 1) * step - 0.5  # de -0.5+step/... a +0.5-step/...
    offsets = np.array([[dx, dy, dz]
                        for dx in positions
                        for dy in positions
                        for dz in positions])
    return offsets  # shape (nsub^3, 3)


def _rho_net_at(species_3d_list, i, j, k, nx, ny, nz):
    """ρ_net en nodo (i,j,k) con índices periódicos."""
    ii, jj, kk = i % nx, j % ny, k % nz
    if len(species_3d_list) >= 2:
        return species_3d_list[0][ii, jj, kk] - species_3d_list[1][ii, jj, kk]
    return species_3d_list[0][ii, jj, kk]


def _quad_weight_1d(t, rm1, r0, rp1):
    """
    Interpolación de Lagrange cuadrática en x = -1, 0, +1.
    Evalúa el polinomio en t ∈ [-0.5, 0.5].
    Puede devolver valores negativos si la parábola lo requiere.
    """
    return rm1 * t * (t - 1.0) / 2.0 - r0 * (t + 1.0) * (t - 1.0) + rp1 * t * (t + 1.0) / 2.0


def _quartic_weight_1d(t, rm2, rm1, r0, rp1, rp2):
    """
    Interpolación de Lagrange de grado 4 en x = -2, -1, 0, +1, +2.
    Evalúa el polinomio en t ∈ [-0.5, 0.5].

    Bases de Lagrange:
      L_{-2}(t) = t(t-1)(t+1)(t-2) / ((-2-(-1))(-2-0)(-2-1)(-2-2))
               = t(t-1)(t+1)(t-2) / ((-1)(-2)(-3)(-4)) = t(t-1)(t+1)(t-2)/24
      L_{-1}(t) = t(t-2)(t+1)(t+2) / ((-1-(-2))(-1-0)(-1-1)(-1-2))
               = t(t-2)(t+1)(t+2) / ((1)(-1)(-2)(-3)) = -t(t-2)(t+1)(t+2)/6
      L_0(t)  = (t-2)(t-1)(t+1)(t+2) / ((0-(-2))(0-(-1))(0-1)(0-2))
               = (t-2)(t-1)(t+1)(t+2) / ((2)(1)(-1)(-2)) = (t-2)(t-1)(t+1)(t+2)/4
      L_{+1}(t) = t(t-2)(t-1)(t+2) / ((1-(-2))(1-(-1))(1-0)(1-2))
               = t(t-2)(t-1)(t+2) / ((3)(2)(1)(-1)) = -t(t-2)(t-1)(t+2)/6
      L_{+2}(t) = t(t-1)(t+1)(t-2) corregido:
               = t(t+1)(t-1)(t+2) / ((2-(-2))(2-(-1))(2-0)(2-1))  ... ver abajo
      L_{+2}(t) = t(t-1)(t+1)(t+2) / ((4)(3)(2)(1)) = t(t-1)(t+1)(t+2)/24
                  pero evaluado en el polinomio correcto para x=+2:
               numerador: (t-(-2))(t-(-1))(t-0)(t-1) = (t+2)(t+1)t(t-1)
               denominador: (2+2)(2+1)(2-0)(2-1) = 4*3*2*1 = 24
    """
    L_m2 =  t * (t - 1.0) * (t + 1.0) * (t - 2.0) / 24.0   # (t+2 cancelado: nodo -2 → base sin (t-(-2)))
    # Rehaciendo correctamente con numeradores explícitos:
    # L_{-2}: producto de (t-xj) para j≠-2: (t+1)(t-0)(t-1)(t-2) / denom
    # denom: (-2+1)(-2-0)(-2-1)(-2-2) = (-1)(-2)(-3)(-4) = 24
    L_m2 = (t + 1.0) * t * (t - 1.0) * (t - 2.0) / 24.0
    # L_{-1}: (t+2)(t-0)(t-1)(t-2) / denom; denom: (-1+2)(-1-0)(-1-1)(-1-2)=(1)(-1)(-2)(-3)=-6
    L_m1 = (t + 2.0) * t * (t - 1.0) * (t - 2.0) / (-6.0)
    # L_0:  (t+2)(t+1)(t-1)(t-2) / denom; denom: (0+2)(0+1)(0-1)(0-2)=(2)(1)(-1)(-2)=4
    L_0  = (t + 2.0) * (t + 1.0) * (t - 1.0) * (t - 2.0) / 4.0
    # L_{+1}: (t+2)(t+1)(t-0)(t-2) / denom; denom: (1+2)(1+1)(1-0)(1-2)=(3)(2)(1)(-1)=-6
    L_p1 = (t + 2.0) * (t + 1.0) * t * (t - 2.0) / (-6.0)
    # L_{+2}: (t+2)(t+1)(t-0)(t-1) / denom; denom: (2+2)(2+1)(2-0)(2-1)=(4)(3)(2)(1)=24
    L_p2 = (t + 2.0) * (t + 1.0) * t * (t - 1.0) / 24.0
    return rm2 * L_m2 + rm1 * L_m1 + r0 * L_0 + rp1 * L_p1 + rp2 * L_p2


def _lagrange_bases_1d(t, nodes):
    """
    Calcula las bases de Lagrange L_k(t) para un conjunto de nodos dados.
    nodes: lista o array de posiciones de los nodos (e.g. [-1,0,1] o [-2,-1,0,1,2])
    Devuelve array de len(nodes) bases evaluadas en t.
    """
    n = len(nodes)
    L = np.ones(n)
    for k in range(n):
        for j in range(n):
            if j != k:
                L[k] *= (t - nodes[j]) / (nodes[k] - nodes[j])
    return L


def _tensor_lagrange_weights(offsets, nodes_1d, rho_cube):
    """
    Interpolación de Lagrange producto tensorial completo.

    offsets:   array (n_sub, 3) con posiciones de subpartículas en [-0.5, 0.5]
    nodes_1d:  posiciones de los nodos vecinos, e.g. [-1,0,1] o [-2,-1,0,1,2]
    rho_cube:  array (len(nodes_1d), len(nodes_1d), len(nodes_1d)) con ρ_net
               en cada nodo (a, b, c) del cubo de vecinos

    Devuelve array (n_sub,) con los pesos interpolados (sin normalizar).
    """
    nodes = np.array(nodes_1d, dtype=float)
    n_sub = len(offsets)
    # Precalcular bases 1D para cada subpartícula en cada eje
    Lx = np.array([_lagrange_bases_1d(offsets[s, 0], nodes) for s in range(n_sub)])  # (n_sub, n_nodes)
    Ly = np.array([_lagrange_bases_1d(offsets[s, 1], nodes) for s in range(n_sub)])
    Lz = np.array([_lagrange_bases_1d(offsets[s, 2], nodes) for s in range(n_sub)])
    # Suma tensorial: w(s) = Σ_{a,b,c} rho_cube[a,b,c] * Lx[s,a] * Ly[s,b] * Lz[s,c]
    w = np.einsum('sa,sb,sc,abc->s', Lx, Ly, Lz, rho_cube)
    return w


def compute_subgrid_forces(species_3d_list, particle_pos, grid_size,
                            charge, epsilon, rcut, nsub, rmin_sub=0.0,
                            gradient_dist=False, quadratic_dist=False,
                            quartic_dist=False, quadratic_full=False,
                            quartic_full=False, rcut_sub=False):
    """
    Para cada voxel dentro de rcut:
      1. Calcula ρ_net y la carga total del voxel (q_vox = ρ_net × 1 lu³).
      2. Distribuye q_vox entre nsub³ subpartículas:
         - uniforme (default): q_s = q_vox / n_sub para todas.
         - gradiente (--gradient-dist): q_s varía linealmente según el gradiente
           de ρ_net por diferencias centradas: q_s = q_base*(1 + gx*dx + gy*dy + gz*dz).
           La suma total se conserva igual a q_vox.
         - cuadrática (--quadratic-dist): peso de cada subpartícula proporcional al
           producto de las interpolaciones de Lagrange cuadráticas 1D en cada eje,
           usando ρ_net en los nodos vecinos. Se normaliza para conservar q_vox.
      3. Excluye subpartículas a distancia < rmin_sub de la partícula;
         su carga se redistribuye uniformemente entre las restantes.
      4. Calcula la fuerza de Coulomb de cada subpartícula activa sobre la partícula.

    Devuelve lista de dicts con info por voxel.
    """
    nx, ny, nz = grid_size
    px, py, pz = particle_pos
    num_species = len(species_3d_list)

    offsets = subgrid_offsets(nsub)   # (nsub^3, 3)
    n_sub = len(offsets)              # nsub^3

    factor = charge / (4.0 * np.pi * epsilon)  # constante de fuerza

    rows = []

    for i in range(nx):
        xi = i + 1  # posición del nodo (coordenada 1-based)
        for j in range(ny):
            yj = j + 1
            for k in range(nz):
                zk = k + 1

                # Vector desde la partícula al centro del voxel
                r_vec = np.array([xi - px, yj - py, zk - pz])
                dist = np.linalg.norm(r_vec)

                if dist > rcut:
                    continue

                # Carga neta del voxel (densidad × volumen = densidad × 1)
                if num_species >= 2:
                    rho_net = (species_3d_list[0][i, j, k]
                               - species_3d_list[1][i, j, k])
                else:
                    rho_net = species_3d_list[0][i, j, k]

                q_vox = rho_net  # volumen = 1 lu³

                # Distribución de carga entre subpartículas
                if quartic_full:
                    # Lagrange grado 4, producto tensorial completo con 5x5x5=125 nodos.
                    # Captura todas las correlaciones cruzadas diagonales.
                    nodes = [-2, -1, 0, 1, 2]
                    rho_cube = np.array([[[
                        _rho_net_at(species_3d_list, i+a, j+b, k+c, nx, ny, nz)
                        for c in nodes] for b in nodes] for a in nodes])
                    w_raw = _tensor_lagrange_weights(offsets, nodes, rho_cube)
                    w_sum = w_raw.sum()
                    if w_sum == 0.0:
                        q_subs = np.full(n_sub, q_vox / n_sub)
                    else:
                        q_subs = q_vox * w_raw / w_sum
                elif quadratic_full:
                    # Lagrange grado 2, producto tensorial completo con 3x3x3=27 nodos.
                    # Captura todas las correlaciones cruzadas diagonales.
                    nodes = [-1, 0, 1]
                    rho_cube = np.array([[[
                        _rho_net_at(species_3d_list, i+a, j+b, k+c, nx, ny, nz)
                        for c in nodes] for b in nodes] for a in nodes])
                    w_raw = _tensor_lagrange_weights(offsets, nodes, rho_cube)
                    w_sum = w_raw.sum()
                    if w_sum == 0.0:
                        q_subs = np.full(n_sub, q_vox / n_sub)
                    else:
                        q_subs = q_vox * w_raw / w_sum
                elif quartic_dist:
                    # Interpolación de Lagrange grado 4 separable por eje.
                    # Usa los 5 nodos: i-2, i-1, i, i+1, i+2 en cada dirección.
                    rxm2 = _rho_net_at(species_3d_list, i-2, j,   k,   nx, ny, nz)
                    rxm1 = _rho_net_at(species_3d_list, i-1, j,   k,   nx, ny, nz)
                    rx0  = rho_net
                    rxp1 = _rho_net_at(species_3d_list, i+1, j,   k,   nx, ny, nz)
                    rxp2 = _rho_net_at(species_3d_list, i+2, j,   k,   nx, ny, nz)
                    rym2 = _rho_net_at(species_3d_list, i,   j-2, k,   nx, ny, nz)
                    rym1 = _rho_net_at(species_3d_list, i,   j-1, k,   nx, ny, nz)
                    ry0  = rho_net
                    ryp1 = _rho_net_at(species_3d_list, i,   j+1, k,   nx, ny, nz)
                    ryp2 = _rho_net_at(species_3d_list, i,   j+2, k,   nx, ny, nz)
                    rzm2 = _rho_net_at(species_3d_list, i,   j,   k-2, nx, ny, nz)
                    rzm1 = _rho_net_at(species_3d_list, i,   j,   k-1, nx, ny, nz)
                    rz0  = rho_net
                    rzp1 = _rho_net_at(species_3d_list, i,   j,   k+1, nx, ny, nz)
                    rzp2 = _rho_net_at(species_3d_list, i,   j,   k+2, nx, ny, nz)
                    wx = _quartic_weight_1d(offsets[:, 0], rxm2, rxm1, rx0, rxp1, rxp2)
                    wy = _quartic_weight_1d(offsets[:, 1], rym2, rym1, ry0, ryp1, ryp2)
                    wz = _quartic_weight_1d(offsets[:, 2], rzm2, rzm1, rz0, rzp1, rzp2)
                    w_raw = wx * wy * wz
                    w_sum = w_raw.sum()
                    if w_sum == 0.0:
                        q_subs = np.full(n_sub, q_vox / n_sub)
                    else:
                        q_subs = q_vox * w_raw / w_sum
                elif quadratic_dist:
                    # Interpolación cuadrática de Lagrange separable por eje.
                    # Para cada eje se usan los 3 nodos vecinos: i-1, i, i+1.
                    # El peso de cada subpartícula es wx*wy*wz (evaluado en su offset),
                    # normalizado para que la suma total sea q_vox.
                    rxm1 = _rho_net_at(species_3d_list, i-1, j,   k,   nx, ny, nz)
                    rx0  = rho_net
                    rxp1 = _rho_net_at(species_3d_list, i+1, j,   k,   nx, ny, nz)
                    rym1 = _rho_net_at(species_3d_list, i,   j-1, k,   nx, ny, nz)
                    ry0  = rho_net
                    ryp1 = _rho_net_at(species_3d_list, i,   j+1, k,   nx, ny, nz)
                    rzm1 = _rho_net_at(species_3d_list, i,   j,   k-1, nx, ny, nz)
                    rz0  = rho_net
                    rzp1 = _rho_net_at(species_3d_list, i,   j,   k+1, nx, ny, nz)
                    wx = _quad_weight_1d(offsets[:, 0], rxm1, rx0, rxp1)
                    wy = _quad_weight_1d(offsets[:, 1], rym1, ry0, ryp1)
                    wz = _quad_weight_1d(offsets[:, 2], rzm1, rz0, rzp1)
                    w_raw = wx * wy * wz
                    w_sum = w_raw.sum()
                    if w_sum == 0.0:
                        q_subs = np.full(n_sub, q_vox / n_sub)
                    else:
                        q_subs = q_vox * w_raw / w_sum
                elif gradient_dist:
                    # Gradiente por diferencias centradas con condiciones periódicas
                    gx = (_rho_net_at(species_3d_list, i+1, j,   k,   nx, ny, nz)
                        - _rho_net_at(species_3d_list, i-1, j,   k,   nx, ny, nz)) / 2.0
                    gy = (_rho_net_at(species_3d_list, i,   j+1, k,   nx, ny, nz)
                        - _rho_net_at(species_3d_list, i,   j-1, k,   nx, ny, nz)) / 2.0
                    gz = (_rho_net_at(species_3d_list, i,   j,   k+1, nx, ny, nz)
                        - _rho_net_at(species_3d_list, i,   j,   k-1, nx, ny, nz)) / 2.0
                    q_base = q_vox / n_sub
                    # sum(offsets) = 0 por simetría → sum(q_subs) = q_vox
                    q_subs = q_base * (1.0
                                       + gx * offsets[:, 0]
                                       + gy * offsets[:, 1]
                                       + gz * offsets[:, 2])
                else:
                    # Distribución uniforme
                    q_subs = np.full(n_sub, q_vox / n_sub)

                # Posiciones de las subpartículas (absolutas)
                sub_positions = np.array([xi, yj, zk]) + offsets  # (n_sub, 3)

                # Vectores desde cada subpartícula hacia la partícula
                r_to_part = np.array([px, py, pz]) - sub_positions  # (n_sub, 3)
                dists_sub = np.linalg.norm(r_to_part, axis=1)       # (n_sub,)

                # Subpartículas activas: fuera del radio mínimo y con r > 0
                # Si rcut_sub=True, también se excluyen las que superan rcut
                active = (dists_sub >= rmin_sub) & (dists_sub > 1e-10)
                if rcut_sub:
                    active &= (dists_sub <= rcut)
                n_active = np.sum(active)
                if n_active == 0:
                    continue

                # Redistribuir carga de excluidas entre activas
                q_excluida = np.sum(q_subs[~active])
                q_subs_active = q_subs.copy()
                q_subs_active[~active] = 0.0
                if q_excluida != 0.0:
                    q_subs_active[active] += q_excluida / n_active

                # Fuerza de Coulomb: F = q_coll * q_s / (4π·ε·r²) * r̂
                # r̂ apunta desde subpartícula hacia la partícula
                f_mag = np.zeros(n_sub)
                f_mag[active] = factor * q_subs_active[active] / dists_sub[active]**2

                r_hat = np.zeros((n_sub, 3))
                r_hat[active] = (r_to_part[active]
                                 / dists_sub[active, np.newaxis])

                # Fuerza total del voxel sobre la partícula (suma de subpartículas activas)
                f_vec_total = np.sum(f_mag[:, np.newaxis] * r_hat, axis=0)
                f_mod_total = np.linalg.norm(f_vec_total)

                # Fuerza sin subdivisión: voxel como punto único en su centro
                r_hat_nodo = -r_vec / dist  # apunta del nodo hacia la partícula
                f_nodo_mag = factor * q_vox / dist**2
                f_nodo_vec = f_nodo_mag * r_hat_nodo

                # Guardar info de todas las subpartículas para el CSV
                subparts = []
                for s_idx in range(n_sub):
                    sp = sub_positions[s_idx]
                    is_active = bool(active[s_idx])
                    fv = f_mag[s_idx] * r_hat[s_idx] if is_active else np.zeros(3)
                    subparts.append({
                        'x': sp[0], 'y': sp[1], 'z': sp[2],
                        'dist': dists_sub[s_idx],
                        'q': q_subs_active[s_idx] if is_active else 0.0,
                        'Fx': fv[0], 'Fy': fv[1], 'Fz': fv[2],
                        'F': abs(f_mag[s_idx]),
                        'excluido': 0 if is_active else 1,
                    })

                rows.append({
                    'ix': i, 'iy': j, 'iz': k,
                    'x': xi, 'y': yj, 'z': zk,
                    'dist': dist,
                    'rho_net': rho_net,
                    'q_vox': q_vox,
                    'Fx': f_vec_total[0],
                    'Fy': f_vec_total[1],
                    'Fz': f_vec_total[2],
                    'F': f_mod_total,
                    'Fx_nodo': f_nodo_vec[0],
                    'Fy_nodo': f_nodo_vec[1],
                    'Fz_nodo': f_nodo_vec[2],
                    'F_nodo': np.linalg.norm(f_nodo_vec),
                    'subparts': subparts,
                })

    # Ordenar por distancia creciente
    rows.sort(key=lambda r: r['dist'])
    return rows


# ---------------------------------------------------------------------------
# Exportar CSV
# ---------------------------------------------------------------------------

def export_csv(rows, particle_pos, csv_file):
    """
    CSV con dos tipos de filas (columna 'tipo': VOXEL o SUB).
    Columnas compartidas: tipo, px, py, pz, ix, iy, iz
    Columnas solo VOXEL:  x, y, z, dist, rho_net, q_vox,
                          Fx, Fy, Fz, F, Fx_nodo, Fy_nodo, Fz_nodo, F_nodo
    Columnas solo SUB:    sub_x, sub_y, sub_z, sub_dist, sub_q,
                          sub_Fx, sub_Fy, sub_Fz, sub_F
    """
    px, py, pz = particle_pos
    header = ("tipo,px,py,pz,ix,iy,iz,"
              "x,y,z,dist,rho_net,q_vox,"
              "Fx,Fy,Fz,F,Fx_nodo,Fy_nodo,Fz_nodo,F_nodo,"
              "sub_x,sub_y,sub_z,sub_dist,sub_q,sub_Fx,sub_Fy,sub_Fz,sub_F,excluido")

    n_sub_total = sum(len(r['subparts']) for r in rows)

    with open(csv_file, 'w') as f:
        f.write(header + '\n')
        for r in rows:
            # Fila del voxel: columnas SUB y excluido vacías
            f.write(
                f"VOXEL,{px:.16e},{py:.16e},{pz:.16e},"
                f"{r['ix']:d},{r['iy']:d},{r['iz']:d},"
                f"{r['x']:.16e},{r['y']:.16e},{r['z']:.16e},"
                f"{r['dist']:.16e},{r['rho_net']:.16e},{r['q_vox']:.16e},"
                f"{r['Fx']:.16e},{r['Fy']:.16e},{r['Fz']:.16e},{r['F']:.16e},"
                f"{r['Fx_nodo']:.16e},{r['Fy_nodo']:.16e},{r['Fz_nodo']:.16e},{r['F_nodo']:.16e},"
                f",,,,,,,,,\n"
            )
            # Filas SUB: cols 7-12 vacías (x,y,z,dist,rho_net,q_vox = 6 campos)
            #            cols 13-20 vacías (Fx,Fy,Fz,F,Fx_nodo,Fy_nodo,Fz_nodo,F_nodo = 8 campos)
            for sp in r['subparts']:
                f.write(
                    f"SUB,{px:.16e},{py:.16e},{pz:.16e},"
                    f"{r['ix']:d},{r['iy']:d},{r['iz']:d},"
                    f",,,,,,"
                    f",,,,,,,,"
                    f"{sp['x']:.16e},{sp['y']:.16e},{sp['z']:.16e},"
                    f"{sp['dist']:.16e},{sp['q']:.16e},"
                    f"{sp['Fx']:.16e},{sp['Fy']:.16e},{sp['Fz']:.16e},{sp['F']:.16e},"
                    f"{sp['excluido']}\n"
                )

    print(f"CSV exportado: {csv_file} ({len(rows)} voxels, {n_sub_total} subpartículas)")


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Calcula fuerza subgrid sobre la partícula desde densidades qsi',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Ejemplos:
    python qsi_subgrid_force.py -n 1000 -L 32 -q 1.0 -eps 1e4 \\
        --rcut 5.0 --nsub 3 --csv fuerza_subgrid
""")

    parser.add_argument('-n', '--nstep', type=int, required=True,
                        help='Número de paso de tiempo')
    parser.add_argument('-L', '--grid-size', type=int, required=True,
                        help='Tamaño de la malla cúbica (NxNxN)')
    parser.add_argument('-q', '--charge', type=float, default=1.0,
                        help='Carga de la partícula coloidal (default: 1.0)')
    parser.add_argument('-eps', '--epsilon', type=float, default=1.0,
                        help='Permitividad dieléctrica (default: 1.0)')
    parser.add_argument('--rcut', type=float, required=True,
                        help='Radio de corte: solo voxels dentro de este radio')
    parser.add_argument('--nsub', type=int, default=3,
                        help='Subdivisión por arista del voxel (default: 3 → 27 subpartículas)')
    parser.add_argument('--rmin-sub', type=float, default=0.0,
                        help='Radio mínimo desde la partícula: subpartículas más cercanas '
                             'se excluyen y su carga se redistribuye entre las restantes '
                             '(default: 0.0, sin exclusión)')
    parser.add_argument('--rcut-sub', action='store_true',
                        help='Aplicar rcut también a las subpartículas individuales. '
                             'Por defecto solo se filtra el voxel por rcut; con esta opción '
                             'las subpartículas que superen rcut también se excluyen '
                             'y su carga se redistribuye entre las activas del mismo voxel.')
    parser.add_argument('--gradient-dist', action='store_true',
                        help='Distribuir carga según gradiente lineal de ρ_net en el voxel '
                             '(diferencias centradas). Por defecto: distribución uniforme.')
    parser.add_argument('--quadratic-dist', action='store_true',
                        help='Distribuir carga según interpolación cuadrática de Lagrange '
                             'usando los 3 nodos vecinos en cada eje (i-1, i, i+1). '
                             'Separable en x, y, z; normalizado para conservar q_vox.')
    parser.add_argument('--quartic-dist', action='store_true',
                        help='Distribuir carga según interpolación de Lagrange de grado 4 '
                             'usando 5 nodos en cada eje (i-2, i-1, i, i+1, i+2). '
                             'Separable en x, y, z; normalizado para conservar q_vox.')
    parser.add_argument('--quadratic-full', action='store_true',
                        help='Lagrange grado 2 con producto tensorial completo (3x3x3=27 nodos). '
                             'Captura correlaciones cruzadas diagonales. Más costoso que --quadratic-dist.')
    parser.add_argument('--quartic-full', action='store_true',
                        help='Lagrange grado 4 con producto tensorial completo (5x5x5=125 nodos). '
                             'Captura correlaciones cruzadas diagonales. Más costoso que --quartic-dist.')
    parser.add_argument('--num-species', type=int, default=2,
                        help='Número de especies en el archivo qsi (default: 2)')
    parser.add_argument('--csv', default=None,
                        help='Nombre base del archivo CSV de salida')
    parser.add_argument('--no-csv', action='store_true',
                        help='No exportar CSV aunque se especifique --csv')
    parser.add_argument('--qsi-file', default=None,
                        help='Ruta explícita al archivo qsi (si no se usa el nombre estándar)')
    parser.add_argument('--cds-file', default=None,
                        help='Ruta explícita al archivo config.cds (si no se usa el nombre estándar)')

    args = parser.parse_args()

    grid_size = (args.grid_size, args.grid_size, args.grid_size)

    # Archivos de entrada
    if args.qsi_file:
        qsi_file = args.qsi_file
    else:
        qsi_file = f"./qsi-{args.nstep:09d}.001-001"

    if args.cds_file:
        cds_file = args.cds_file
    else:
        cds_file = f"./config.cds{args.nstep:08d}.001-001"

    # Leer posición de la partícula
    print(f"Leyendo posición de la partícula desde: {cds_file}")
    particle_pos = ColloidReader.read_colloid_position_cds(cds_file)
    print(f"  Posición: {particle_pos}")

    # Leer densidades qsi
    print(f"Leyendo densidades de carga desde: {qsi_file}")
    reader = QsiReader(qsi_file, grid_size, num_species=args.num_species)
    species_list = reader.read()
    num_species = reader.num_species
    print(f"  Número de especies: {num_species}")
    for s, sp in enumerate(species_list):
        print(f"  Especie {s}: min={sp.min():.4e}, max={sp.max():.4e}, sum={sp.sum():.4e}")

    # Calcular fuerzas subgrid
    if args.quartic_full:
        dist_mode = 'cuártica grado 4 tensorial completo (5x5x5=125 nodos)'
    elif args.quadratic_full:
        dist_mode = 'cuadrática tensorial completo (3x3x3=27 nodos)'
    elif args.quartic_dist:
        dist_mode = 'cuártica grado 4 separable (5 nodos/eje)'
    elif args.quadratic_dist:
        dist_mode = 'cuadrática separable (3 nodos/eje)'
    elif args.gradient_dist:
        dist_mode = 'gradiente lineal'
    else:
        dist_mode = 'uniforme'
    print(f"\nCalculando fuerzas subgrid (rcut={args.rcut}, nsub={args.nsub}, "
          f"{args.nsub**3} subpartículas/voxel, rmin_sub={args.rmin_sub}, "
          f"distribución={dist_mode})...")
    rows = compute_subgrid_forces(
        species_list, particle_pos, grid_size,
        charge=args.charge, epsilon=args.epsilon,
        rcut=args.rcut, nsub=args.nsub, rmin_sub=args.rmin_sub,
        gradient_dist=args.gradient_dist,
        quadratic_dist=args.quadratic_dist,
        quartic_dist=args.quartic_dist,
        quadratic_full=args.quadratic_full,
        quartic_full=args.quartic_full,
        rcut_sub=args.rcut_sub)

    print(f"  Voxels procesados: {len(rows)}")

    # Fuerza total subgrid (suma vectorial de las fuerzas de cada voxel)
    Fx_tot = sum(r['Fx'] for r in rows)
    Fy_tot = sum(r['Fy'] for r in rows)
    Fz_tot = sum(r['Fz'] for r in rows)
    F_tot  = np.sqrt(Fx_tot**2 + Fy_tot**2 + Fz_tot**2)

    # Fuerza total sin subdivisión (para comparación)
    Fx_nodo_tot = sum(r['Fx_nodo'] for r in rows)
    Fy_nodo_tot = sum(r['Fy_nodo'] for r in rows)
    Fz_nodo_tot = sum(r['Fz_nodo'] for r in rows)
    F_nodo_tot  = np.sqrt(Fx_nodo_tot**2 + Fy_nodo_tot**2 + Fz_nodo_tot**2)

    print(f"\nFuerza total sobre la partícula (subgrid, {args.nsub}³ subpartículas/voxel):")
    print(f"  Fx = {Fx_tot:.16e}")
    print(f"  Fy = {Fy_tot:.16e}")
    print(f"  Fz = {Fz_tot:.16e}")
    print(f"  |F| = {F_tot:.16e}")

    print(f"\nFuerza total sobre la partícula (sin subdivisión, nodo como punto único):")
    print(f"  Fx = {Fx_nodo_tot:.16e}")
    print(f"  Fy = {Fy_nodo_tot:.16e}")
    print(f"  Fz = {Fz_nodo_tot:.16e}")
    print(f"  |F| = {F_nodo_tot:.16e}")

    # Exportar CSV
    if args.csv is not None and not args.no_csv:
        csv_base, csv_ext = os.path.splitext(args.csv)
        if not csv_ext:
            csv_ext = '.csv'
        csv_file = f"{csv_base}-n_{args.nstep:09d}{csv_ext}"
        print(f"\nExportando a CSV: {csv_file}")
        export_csv(rows, particle_pos, csv_file)

    print("\nAnálisis completado!")


if __name__ == '__main__':
    main()
