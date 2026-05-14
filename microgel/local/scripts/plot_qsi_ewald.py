#!/usr/bin/env python3
"""
Script para graficar la densidad de carga (qsi) calculada por Ludwig (Ewald),
comparándola con la solución teórica de Debye-Hückel.

Los archivos qsi contienen la densidad de carga por especie en cada nodo:
  - Columna 0: densidad de cationes ρ₊ (≥ 0)
  - Columna 1: densidad de aniones  ρ₋ (≥ 0)
  Carga neta: ρ_net = ρ₊ - ρ₋

La densidad de carga neta teórica de Debye-Hückel es:
  ρ_net(r) = -κ² ε kT / e · ψ_DH(r) / (4π)
equivalente a:
  ρ_net(r) = -κ² q exp(-κ r) / (4π r)    [en unidades de Ludwig]

Uso básico:
    python plot_qsi_ewald.py -n 2000 -L 32 -q 1.0 -kt 0.00001 -eps 1.0e4 \\
        --rho-el 8.0e-3 --alpha 0.7 --rc 3.0

Mostrar especie 0 y carga neta:
    python plot_qsi_ewald.py -n 2000 -L 32 -q 1.0 -kt 0.00001 -eps 1.0e4 \\
        --rho-el 8.0e-3 --show-species0 --show-species1 --show-net

Con dirección y CSV:
    python plot_qsi_ewald.py -n 2000 -L 32 -q 1.0 -kt 0.00001 -eps 1.0e4 \\
        --rho-el 8.0e-3 --direction x --angle-tol 5 \\
        --csv qsi_data
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
import os

font_size = 16


# ---------------------------------------------------------------------------
# Lector de archivos qsi (dos escalares por nodo: rho+, rho-)
# ---------------------------------------------------------------------------

class QsiReader:
    """Lee y procesa archivos qsi de Ludwig (densidades por especie)"""

    def __init__(self, filename, grid_size, num_species=2):
        self.filename = filename
        self.nx, self.ny, self.nz = grid_size
        self.num_species = num_species
        self.species = []

    def read(self):
        """Lee el archivo qsi y organiza cada especie en malla 3D"""
        data = np.loadtxt(self.filename)
        n_points = self.nx * self.ny * self.nz

        if data.ndim == 1:
            # Archivo con una sola columna (carga neta de coloides)
            print("Archivo con 1 columna detectado (probablemente qsi_colloid)")
            data = data.reshape(-1, 1)
            self.num_species = 1
        elif data.shape[0] != n_points:
            raise ValueError(
                f"El archivo tiene {data.shape[0]} filas, "
                f"pero se esperaban {n_points}")

        if data.shape[1] != self.num_species:
            print(f"Advertencia: Se encontraron {data.shape[1]} columnas, "
                  f"pero se esperaban {self.num_species}")
            self.num_species = data.shape[1]

        # Ludwig escribe en orden x-major: para cada x, para cada y, todos los z
        self.species = []
        for i in range(self.num_species):
            species_data = data[:, i].reshape((self.nx, self.ny, self.nz))
            if self.num_species > 1:
                min_val = np.min(species_data)
                if min_val < -1e-15:
                    print(f"ADVERTENCIA: Especie {i} tiene valores negativos "
                          f"(min={min_val:.4e}). Las densidades deben ser ≥ 0.")
            self.species.append(species_data)

        return self.species

    def net_charge_3d(self):
        """ρ_net = ρ₊ - ρ₋"""
        if len(self.species) < 2:
            return self.species[0]
        return self.species[0] - self.species[1]


# ---------------------------------------------------------------------------
# Lector de posición del coloide
# ---------------------------------------------------------------------------

class ColloidReader:

    @staticmethod
    def read_colloid_position_cds(filename):
        """Lee posición desde config.cds*.001-001 (posición en línea 36)"""
        with open(filename, 'r') as f:
            lines = f.readlines()
        pos_line = lines[35].strip()
        parts = pos_line.split()
        if len(parts) < 3:
            raise ValueError(
                f"Línea 36 de {filename} no tiene 3 valores: '{pos_line}'")
        return np.array([float(parts[0]), float(parts[1]), float(parts[2])])


# ---------------------------------------------------------------------------
# Teoría de Debye-Hückel para densidad de carga
# ---------------------------------------------------------------------------

class TheoreticalChargeDensity:
    """
    Densidad de carga neta de Debye-Hückel:
      ρ_net(r) = -κ² · q · exp(-κ·r) / (4π·r)    [sistema infinito]

    Esta es la derivada del potencial DH:
      ψ_DH(r) = (q / 4π·ε·kT) · exp(-κ·r) / r
    mediante la ecuación de Poisson:
      ρ_net(r) = -ε·kT · ∇²ψ / (4π)   (en unidades Gausianas de Ludwig)
    """

    @staticmethod
    def compute_dh(distances, q=1.0, epsilon=None, kt=None, kappa=0.0):
        """
        Densidad de carga neta DH para sistema infinito.
        ρ_net(r) = -κ² · q/(4π·r) · exp(-κ·r)   si kappa > 0
        ρ_net(r) = 0                               si kappa = 0 (Coulomb puro)
        """
        r = np.maximum(distances, 0.01)
        if kappa <= 0:
            return np.zeros_like(r)
        return -kappa**2 * q / (4.0 * np.pi * r) * np.exp(-kappa * r)

    @staticmethod
    def compute_dh_finite_radius(distances, q=1.0, kappa=0.0, a=0.0):
        """
        Densidad de carga DH para partícula esférica de radio 'a' (sistema infinito).

        Solución exacta linealizada de Poisson-Boltzmann para r > a:
            ρ_net(r) = -κ² q exp(-κ(r-a)) / (4π r (1 + κa))

        Para r ≤ a (interior de la partícula) retorna NaN.
        """
        r = np.asarray(distances, dtype=float)
        if kappa <= 0:
            return np.zeros_like(r)
        result = np.full_like(r, np.nan)
        outside = r > a
        r_out = r[outside]
        result[outside] = (
            -kappa**2 * q * np.exp(-kappa * (r_out - a))
            / (4.0 * np.pi * r_out * (1.0 + kappa * a))
        )
        return result

    @staticmethod
    def compute_dh_periodic(distances, q=1.0, epsilon=None, kt=None, kappa=0.0,
                             L=32.0, n_shells=1, direction=None):
        """
        ρ_net periódico a lo largo de una dirección dada.
        Suma imagen central + imágenes periódicas.
        """
        if direction is None:
            direction = np.array([1.0, 0.0, 0.0])
        d = np.asarray(direction, dtype=float)
        d = d / np.linalg.norm(d)

        r = np.maximum(distances, 0.01)
        r_vec = r[:, np.newaxis] * d[np.newaxis, :]  # (N, 3)

        def rho_scalar(dist):
            dist = np.maximum(dist, 1e-10)
            if kappa <= 0:
                return np.zeros_like(dist)
            return -kappa**2 * q / (4.0 * np.pi * dist) * np.exp(-kappa * dist)

        rho_total = rho_scalar(r)

        for inx in range(-n_shells, n_shells + 1):
            for iny in range(-n_shells, n_shells + 1):
                for inz in range(-n_shells, n_shells + 1):
                    if inx == 0 and iny == 0 and inz == 0:
                        continue
                    r_img = np.array([inx * L, iny * L, inz * L])
                    dr = r_vec - r_img[np.newaxis, :]
                    dist_img = np.linalg.norm(dr, axis=1)
                    rho_total = rho_total + rho_scalar(dist_img)

        return rho_total

    @staticmethod
    def compute_dh_periodic_at_points(r_vecs, q=1.0, epsilon=None, kt=None,
                                       kappa=0.0, L=32.0, n_shells=1):
        """ρ_net periódico en puntos 3D arbitrarios r_vecs (N,3) relativos a la partícula"""
        r_vecs = np.asarray(r_vecs, dtype=float)

        def rho_scalar(dist):
            dist = np.maximum(dist, 1e-10)
            if kappa <= 0:
                return np.zeros_like(dist)
            return -kappa**2 * q / (4.0 * np.pi * dist) * np.exp(-kappa * dist)

        dist0 = np.linalg.norm(r_vecs, axis=1)
        rho_total = rho_scalar(dist0)

        for inx in range(-n_shells, n_shells + 1):
            for iny in range(-n_shells, n_shells + 1):
                for inz in range(-n_shells, n_shells + 1):
                    if inx == 0 and iny == 0 and inz == 0:
                        continue
                    r_img = np.array([inx * L, iny * L, inz * L])
                    dr = r_vecs - r_img[np.newaxis, :]
                    dist_img = np.linalg.norm(dr, axis=1)
                    rho_total = rho_total + rho_scalar(dist_img)

        return rho_total

    # CHANGE INIT - DHPeriodicFiniteRadius - DH periódico con partícula de radio finito a
    @staticmethod
    def compute_dh_periodic_finite_radius(distances, q=1.0, kappa=0.0, a=0.0,
                                           L=32.0, n_shells=1, direction=None):
        """
        ρ_net DH periódico con partícula de radio finito 'a' a lo largo de
        una dirección dada.

        Imagen central: ρ_DH(r; a) = -κ² q exp(-κ(r-a)) / (4π r (1+κa))  para r > a
        Imágenes periódicas: usan la fórmula de punto (r=0 → imagen es punto).

        Para r ≤ a retorna NaN.
        """
        if direction is None:
            direction = np.array([1.0, 0.0, 0.0])
        d = np.asarray(direction, dtype=float)
        d = d / np.linalg.norm(d)

        r = np.asarray(distances, dtype=float)
        r_vec = r[:, np.newaxis] * d[np.newaxis, :]  # (N, 3)

        def rho_point(dist):
            dist = np.maximum(dist, 1e-10)
            if kappa <= 0:
                return np.zeros_like(dist)
            return -kappa**2 * q / (4.0 * np.pi * dist) * np.exp(-kappa * dist)

        # Imagen central con radio finito
        result = np.full_like(r, np.nan)
        outside = r > a
        r_out = r[outside]
        if kappa > 0:
            result[outside] = (
                -kappa**2 * q * np.exp(-kappa * (r_out - a))
                / (4.0 * np.pi * r_out * (1.0 + kappa * a))
            )
            # Para r <= a dentro del rango, dejar NaN (interior de la partícula)
        else:
            result[outside] = 0.0

        # Imágenes periódicas (tratadas como punto)
        for inx in range(-n_shells, n_shells + 1):
            for iny in range(-n_shells, n_shells + 1):
                for inz in range(-n_shells, n_shells + 1):
                    if inx == 0 and iny == 0 and inz == 0:
                        continue
                    r_img = np.array([inx * L, iny * L, inz * L])
                    dr = r_vec - r_img[np.newaxis, :]
                    dist_img = np.linalg.norm(dr, axis=1)
                    contrib = rho_point(dist_img)
                    # Solo sumar en nodos fuera de la partícula central
                    result[outside] += contrib[outside]

        return result
    # CHANGE END - DHPeriodicFiniteRadius

    # CHANGE INIT - GaussianDH - Densidad de carga iónica DH con fuente gaussiana
    @staticmethod
    def compute_gaussian_dh(distances, q=1.0, kappa=0.0, sigma=1.0):
        """
        Densidad de carga IÓNICA de Debye-Hückel con fuente gaussiana:

            ρ_iones(r) = -κ² ψ_gauss(r)

        donde ψ_gauss(r) es el potencial DH para fuente gaussiana, calculado vía
        la función de Green de Yukawa en simetría esférica (numéricamente estable):

            ψ_gauss(r) = (q/4π) · (1/κr) · { e^{-κr} ∫₀ʳ sinh(κr') r' ρ̃(r') dr'
                                             + sinh(κr) ∫ᵣ^∞ e^{-κr'} r' ρ̃(r') dr' }

        con ρ̃(r') = 1/(π^{3/2}σ³) exp(-r'²/σ²)  [distribución gaussiana normalizada].

        Esto se compara con ρ_net = ρ₊ - ρ₋ de Ludwig (iones en solución,
        la carga del coloide NO está en qsi).

        Args:
            distances: array de distancias
            q:        carga total Q
            kappa:    parámetro de Debye κ
            sigma:    anchura de la distribución gaussiana σ
        """
        r = np.asarray(distances, dtype=float)
        r_safe = np.maximum(r, 1e-10)

        if kappa <= 0:
            return np.zeros_like(r_safe)

        # Pre-computar integranda sobre malla fina de r'
        prefactor = q / (4.0 * np.pi)
        rmax = max(15.0 * sigma, float(np.max(r_safe)) * 1.5)
        n_rho = 20000
        rp = np.linspace(1e-8, rmax, n_rho)
        drp = rp[1] - rp[0]
        rho_rp = 1.0 / (np.pi**1.5 * sigma**3) * np.exp(-rp**2 / sigma**2) * rp

        integrand_inner = np.sinh(kappa * rp) * rho_rp
        integrand_outer = np.exp(-kappa * rp) * rho_rp

        cumul_inner = np.cumsum(integrand_inner) * drp
        total_outer = np.sum(integrand_outer) * drp
        cumul_outer = total_outer - np.cumsum(integrand_outer) * drp

        inner_at_r = np.interp(r_safe, rp, cumul_inner)
        outer_at_r = np.interp(r_safe, rp, cumul_outer)

        # Factor 4π de la reducción de la integral 3D a 1D en simetría esférica
        psi_gauss = (4.0 * np.pi * prefactor / (kappa * r_safe)) * (
            np.exp(-kappa * r_safe) * inner_at_r
            + np.sinh(kappa * r_safe) * outer_at_r
        )

        # ρ_iones = -κ² ψ_gauss  (en unidades de Ludwig)
        return -kappa**2 * psi_gauss

    @staticmethod
    def compute_gaussian_dh_periodic(distances, q=1.0, kappa=0.0, sigma=1.0,
                                     L=32.0, n_shells=1, direction=None):
        """
        ρ_ions gaussiana-DH periódico a lo largo de una dirección.
        Imagen central con fuente gaussiana; imágenes periódicas como Yukawa puntual.
        """
        if direction is None:
            direction = np.array([1.0, 0.0, 0.0])
        d = np.asarray(direction, dtype=float)
        d = d / np.linalg.norm(d)

        r = np.asarray(distances, dtype=float)
        r_safe = np.maximum(r, 1e-10)
        r_vec = r_safe[:, np.newaxis] * d[np.newaxis, :]

        result = TheoreticalChargeDensity.compute_gaussian_dh(
            r_safe, q=q, kappa=kappa, sigma=sigma)

        def rho_point(dist):
            dist = np.maximum(dist, 1e-10)
            if kappa <= 0:
                return np.zeros_like(dist)
            return -kappa**2 * q / (4.0 * np.pi * dist) * np.exp(-kappa * dist)

        for inx in range(-n_shells, n_shells + 1):
            for iny in range(-n_shells, n_shells + 1):
                for inz in range(-n_shells, n_shells + 1):
                    if inx == 0 and iny == 0 and inz == 0:
                        continue
                    r_img = np.array([inx * L, iny * L, inz * L])
                    dr = r_vec - r_img[np.newaxis, :]
                    dist_img = np.linalg.norm(dr, axis=1)
                    result = result + rho_point(dist_img)

        return result

    @staticmethod
    def compute_gaussian_dh_periodic_at_points(r_vecs, q=1.0, kappa=0.0, sigma=1.0,
                                               L=32.0, n_shells=1):
        """
        ρ_ions gaussiana-DH periódico en puntos 3D arbitrarios r_vecs (N,3).
        """
        r_vecs = np.asarray(r_vecs, dtype=float)
        dist0 = np.linalg.norm(r_vecs, axis=1)
        result = TheoreticalChargeDensity.compute_gaussian_dh(
            dist0, q=q, kappa=kappa, sigma=sigma)

        def rho_point(dist):
            dist = np.maximum(dist, 1e-10)
            if kappa <= 0:
                return np.zeros_like(dist)
            return -kappa**2 * q / (4.0 * np.pi * dist) * np.exp(-kappa * dist)

        for inx in range(-n_shells, n_shells + 1):
            for iny in range(-n_shells, n_shells + 1):
                for inz in range(-n_shells, n_shells + 1):
                    if inx == 0 and iny == 0 and inz == 0:
                        continue
                    r_img = np.array([inx * L, iny * L, inz * L])
                    dr = r_vecs - r_img[np.newaxis, :]
                    dist_img = np.linalg.norm(dr, axis=1)
                    result = result + rho_point(dist_img)

        return result
    # CHANGE END - GaussianDH


# ---------------------------------------------------------------------------
# Filtro de dirección
# ---------------------------------------------------------------------------

def parse_direction(direction_str):
    if direction_str is None:
        return None
    aliases = {
        'x':   np.array([1.0, 0.0, 0.0]),
        'y':   np.array([0.0, 1.0, 0.0]),
        'z':   np.array([0.0, 0.0, 1.0]),
        'xy':  np.array([1.0, 1.0, 0.0]),
        'xz':  np.array([1.0, 0.0, 1.0]),
        'yz':  np.array([0.0, 1.0, 1.0]),
        'xyz': np.array([1.0, 1.0, 1.0]),
    }
    s = direction_str.strip().lower()
    if s in aliases:
        v = aliases[s]
    else:
        try:
            parts = [float(x) for x in s.split(',')]
            if len(parts) != 3:
                raise ValueError
            v = np.array(parts)
        except ValueError:
            raise ValueError(
                f"--direction '{direction_str}' no reconocido. "
                f"Use 'x','y','z','xy','xz','yz','xyz' o 'dx,dy,dz'.")
    norm = np.linalg.norm(v)
    if norm < 1e-10:
        raise ValueError(f"--direction '{direction_str}' es el vector cero.")
    return v / norm


def _direction_label(d):
    if d is None:
        return 'todas'
    return f"({d[0]:.3f},{d[1]:.3f},{d[2]:.3f})"


# ---------------------------------------------------------------------------
# Extracción de puntos
# ---------------------------------------------------------------------------

def extract_points_qsi(species_3d_list, particle_pos, grid_size,
                        min_radius=0.5, max_radius=None,
                        direction=None, angle_tol_deg=15.0):
    """
    Extrae distancias, densidades por especie, y carga neta para cada nodo.

    Returns:
        distances  : (N,)
        species_vals: list of (N,) arrays, one per species
        net_vals   : (N,)
        r_vecs     : (N,3)
    """
    nx, ny, nz = grid_size
    Lx, Ly, Lz = float(nx), float(ny), float(nz)
    # CHANGE INIT - PBCExtract - Aplicar PBC en dr y corregir max_radius default
    # Sin PBC: nodos del lado opuesto del box tienen dr grande y quedan excluidos.
    # El default L/2 también cortaba puntos accesibles en la diagonal del box.
    if max_radius is None:
        max_radius = np.sqrt((Lx/2)**2 + (Ly/2)**2 + (Lz/2)**2)
    # CHANGE END - PBCExtract

    cos_tol = np.cos(np.radians(angle_tol_deg)) if direction is not None else None
    px, py, pz = particle_pos
    num_species = len(species_3d_list)

    distances = []
    r_vecs = []
    species_vals = [[] for _ in range(num_species)]

    for i in range(nx):
        xi = i + 1
        for j in range(ny):
            yj = j + 1
            for k in range(nz):
                zk = k + 1
                # CHANGE INIT - PBCExtract - wrap dr a imagen más cercana
                drx = xi - px; drx -= Lx * round(drx / Lx)
                dry = yj - py; dry -= Ly * round(dry / Ly)
                drz = zk - pz; drz -= Lz * round(drz / Lz)
                dr = np.array([drx, dry, drz])
                # CHANGE END - PBCExtract
                dist = np.linalg.norm(dr)

                if dist < min_radius or dist > max_radius:
                    continue
                if direction is not None:
                    cos_angle = abs(np.dot(dr / dist, direction))
                    if cos_angle < cos_tol:
                        continue

                distances.append(dist)
                r_vecs.append(dr)
                for s in range(num_species):
                    species_vals[s].append(species_3d_list[s][i, j, k])

    distances = np.array(distances)
    r_vecs = np.array(r_vecs) if len(r_vecs) > 0 else np.zeros((0, 3))
    species_vals = [np.array(sv) for sv in species_vals]

    if num_species >= 2:
        net_vals = species_vals[0] - species_vals[1]
    else:
        net_vals = species_vals[0]

    return distances, species_vals, net_vals, r_vecs


def _sort_subsample(dist, *arrays, max_pts=10000, rmin=0.5):
    """Ordena por distancia, aplica rmin, submuestrea. Devuelve tupla de arrays."""
    idx = np.argsort(dist)
    dist = dist[idx]
    sorted_arrays = [a[idx] for a in arrays]

    mask = dist >= rmin
    dist = dist[mask]
    sorted_arrays = [a[mask] for a in sorted_arrays]

    if len(dist) > max_pts:
        step = len(dist) // max_pts
        dist = dist[::step]
        sorted_arrays = [a[::step] for a in sorted_arrays]

    return (dist, *sorted_arrays)


# ---------------------------------------------------------------------------
# Promedio radial en cascaras esféricas
# ---------------------------------------------------------------------------

def compute_radial_average(distances, species_vals, net_vals,
                            shell_width, rmin, rmax):
    """
    Calcula el promedio de carga en cascaras esféricas de grosor shell_width.

    Para cada cascaras [r, r+dr):
      - r_mid: radio central de la cascara
      - n_nodes: número de nodos en la cascara
      - mean_rho_s: promedio de cada especie
      - mean_rho_net: promedio de carga neta
      - std_rho_net: desviación estándar de la carga neta

    Returns:
        r_mids      : (M,) centros de cada cascara
        mean_species: list of (M,) arrays, uno por especie
        mean_net    : (M,)
        std_net     : (M,)
        n_nodes     : (M,) número de nodos por cascara
    """
    r_edges = np.arange(0.0, rmax + shell_width, shell_width)
    r_mids = []
    mean_species = [[] for _ in range(len(species_vals))]
    mean_net = []
    std_net = []
    n_nodes = []

    for r_lo, r_hi in zip(r_edges[:-1], r_edges[1:]):
        mask = (distances >= r_lo) & (distances < r_hi)
        count = np.sum(mask)
        r_mids.append(0.5 * (r_lo + r_hi))
        n_nodes.append(count)
        if count > 0:
            for s, sv in enumerate(species_vals):
                mean_species[s].append(np.mean(sv[mask]))
            mean_net.append(np.mean(net_vals[mask]))
            std_net.append(np.std(net_vals[mask]))
        else:
            for s in range(len(species_vals)):
                mean_species[s].append(np.nan)
            mean_net.append(np.nan)
            std_net.append(np.nan)

    r_mids = np.array(r_mids)
    mean_species = [np.array(ms) for ms in mean_species]
    mean_net = np.array(mean_net)
    std_net = np.array(std_net)
    n_nodes = np.array(n_nodes)

    return r_mids, mean_species, mean_net, std_net, n_nodes


def plot_radial_average(r_mids, mean_species, mean_net, std_net,
                         grid_size, charge, epsilon, kt, kappa,
                         shell_width, rmin,
                         distances=None, species_vals=None, net_vals=None,
                         r_vecs=None,
                         scatter=False,
                         output_file='qsi_radial.png',
                         log_scale=False, log_log=False,
                         xmin=None, ymin=None,
                         show_species0=True, show_species1=True, show_net=True,
                         show_dh=True, show_dh_periodic=False,
                         dh_periodic_shells=1,
                         show_dh_periodic_pts=False,
                         show_error_inf=True,
                         error_ymax=None,
                         particle_radius=0.0,
                         gaussian_source=False, sigma=1.0,
                         gaussian_error=False, gaussian_error_periodic=False):

    L = float(grid_size[0])
    use_log = log_scale or log_log
    abs_suffix = ' |·|' if use_log else ''

    # Mascara de cascaras con datos
    valid = ~np.isnan(mean_net)

    r_max_data = float(np.nanmax(r_mids[valid])) * 1.1 if valid.any() else 16.0
    r_max_theory = max(r_max_data, 16.0)
    r_theory = np.linspace(rmin, r_max_theory, 500)

    rho_dh_inf_th = TheoreticalChargeDensity.compute_dh(
        r_theory, q=charge, kappa=kappa)

    rho_dh_rad_th = None
    if particle_radius > 0 and kappa > 0:
        rho_dh_rad_th = TheoreticalChargeDensity.compute_dh_finite_radius(
            r_theory, q=charge, kappa=kappa, a=particle_radius)

    rho_dh_per_th = None
    if show_dh_periodic:
        print(f"  Calculando ρ_net DH periódico para promedio radial...")
        rho_dh_per_th = TheoreticalChargeDensity.compute_dh_periodic(
            r_theory, q=charge, kappa=kappa,
            L=L, n_shells=dh_periodic_shells)

    # CHANGE INIT - GaussianDH - Cálculo de curvas gaussiana-DH para qsi
    rho_gauss_dh_th = None
    rho_gauss_dh_per_th = None
    if gaussian_source:
        print(f"  Calculando ρ_ions gaussiana-DH (σ={sigma:.3f})...")
        rho_gauss_dh_th = TheoreticalChargeDensity.compute_gaussian_dh(
            r_theory, q=charge, kappa=kappa, sigma=sigma)
        if show_dh_periodic and gaussian_error_periodic:
            n_img = (2 * dh_periodic_shells + 1) ** 3 - 1
            print(f"  Calculando ρ_ions gaussiana-DH periódico ({dh_periodic_shells} capa(s), "
                  f"{n_img} imgs)...")
            rho_gauss_dh_per_th = TheoreticalChargeDensity.compute_gaussian_dh_periodic(
                r_theory, q=charge, kappa=kappa, sigma=sigma,
                L=L, n_shells=dh_periodic_shells)
    # CHANGE END - GaussianDH

    # DH evaluado en los centros de cascara (para panel de error en modo media)
    rho_dh_inf_pts = TheoreticalChargeDensity.compute_dh(
        r_mids, q=charge, kappa=kappa)
    if particle_radius > 0 and kappa > 0:
        rho_dh_rad_pts = TheoreticalChargeDensity.compute_dh_finite_radius(
            r_mids, q=charge, kappa=kappa, a=particle_radius)
    else:
        rho_dh_rad_pts = None

    # DH periódico evaluado en posiciones 3D exactas de cada nodo
    rho_dh_per_pts = None
    if show_dh_periodic_pts and r_vecs is not None and kappa > 0:
        n_img = (2 * dh_periodic_shells + 1) ** 3 - 1
        print(f"  Calculando ρ_net DH periódico en posiciones exactas: "
              f"{dh_periodic_shells} capa(s), {n_img} imágenes, {len(distances)} nodos...")
        rho_dh_per_pts = TheoreticalChargeDensity.compute_dh_periodic_at_points(
            r_vecs, q=charge, kappa=kappa, L=L, n_shells=dh_periodic_shells)

    # CHANGE INIT - GaussianDH - Calcular curvas gaussiana en centros de cascara para error
    rho_gauss_dh_pts = None
    rho_gauss_dh_per_pts = None
    if gaussian_source and gaussian_error:
        rho_gauss_dh_pts = TheoreticalChargeDensity.compute_gaussian_dh(
            r_mids, q=charge, kappa=kappa, sigma=sigma)
    if gaussian_source and gaussian_error_periodic and show_dh_periodic and r_vecs is not None and kappa > 0:
        n_img = (2 * dh_periodic_shells + 1) ** 3 - 1
        print(f"  Calculando ρ_ions gaussiana-DH periódico en posiciones exactas: "
              f"{dh_periodic_shells} capa(s), {n_img} imágenes, {len(distances)} nodos...")
        rho_gauss_dh_per_pts = TheoreticalChargeDensity.compute_gaussian_dh_periodic_at_points(
            r_vecs, q=charge, kappa=kappa, sigma=sigma, L=L, n_shells=dh_periodic_shells)
    # CHANGE END - GaussianDH

    show_panel2 = show_net and kappa > 0 and (
        (show_error_inf and not scatter) or
        (show_dh_periodic_pts and r_vecs is not None) or
        (gaussian_source and gaussian_error) or
        (gaussian_source and gaussian_error_periodic and show_dh_periodic)
    )

    _, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 11),
                                  gridspec_kw={'height_ratios': [2, 1]})

    species_colors = ['tab:blue', 'tab:red', 'tab:green', 'tab:orange']
    species_labels = ['ρ₊ (cationes)', 'ρ₋ (aniones)', 'ρ₂', 'ρ₃']

    def _pv(y):
        return np.abs(y) if use_log else y

    if scatter and distances is not None and net_vals is not None:
        # Modo scatter: todos los puntos individuales, coloreados por cascara
        r_max_edges = r_max_data / 1.1 + shell_width
        r_edges = np.arange(0.0, r_max_edges + shell_width, shell_width)
        cmap = plt.get_cmap('viridis', max(len(r_edges) - 1, 1))

        def _weighted_r(dist_in_shell, weights):
            """Distancia media pesada por |peso|; si todos cero devuelve media aritmética."""
            w = np.abs(weights)
            w_sum = w.sum()
            if w_sum == 0:
                return dist_in_shell.mean()
            return np.dot(w, dist_in_shell) / w_sum

        if show_species0 and species_vals is not None and len(species_vals) >= 1:
            for s_idx, (r_lo, r_hi) in enumerate(zip(r_edges[:-1], r_edges[1:])):
                mask = (distances >= r_lo) & (distances < r_hi)
                if not np.any(mask):
                    continue
                ax1.plot(distances[mask], _pv(species_vals[0][mask]), '.',
                         color=species_colors[0], markersize=3, alpha=0.4,
                         label=f'Ludwig: {species_labels[0]}{abs_suffix}' if s_idx == 0 else '_')

        if show_species1 and species_vals is not None and len(species_vals) >= 2:
            for s_idx, (r_lo, r_hi) in enumerate(zip(r_edges[:-1], r_edges[1:])):
                mask = (distances >= r_lo) & (distances < r_hi)
                if not np.any(mask):
                    continue
                ax1.plot(distances[mask], _pv(species_vals[1][mask]), '.',
                         color=species_colors[1], markersize=3, alpha=0.4,
                         label=f'Ludwig: {species_labels[1]}{abs_suffix}' if s_idx == 0 else '_')

        if show_net:
            for s_idx, (r_lo, r_hi) in enumerate(zip(r_edges[:-1], r_edges[1:])):
                mask = (distances >= r_lo) & (distances < r_hi)
                if not np.any(mask):
                    continue
                ax1.plot(distances[mask], _pv(net_vals[mask]), '.',
                         color='green', markersize=5, alpha=0.6,
                         label=f'Ludwig: ρ_net{abs_suffix} (Δr={shell_width:.2f})' if s_idx == 0 else '_')

        ylabel_mode = 'puntos individuales'
    else:
        # Modo media ± std: x en distancia media pesada por |carga|
        def _weighted_r_shell(dist_all, vals_all, r_lo, r_hi):
            mask = (dist_all >= r_lo) & (dist_all < r_hi)
            if not np.any(mask):
                return np.nan
            d = dist_all[mask]
            w = np.abs(vals_all[mask])
            w_sum = w.sum()
            if w_sum == 0:
                return d.mean()
            return np.dot(w, d) / w_sum

        if show_species0 and len(mean_species) >= 1 and distances is not None:
            ms0 = mean_species[0]
            r_edges_avg = np.arange(0.0, r_mids[-1] + shell_width, shell_width)
            xw0 = np.array([_weighted_r_shell(distances, species_vals[0], r_lo, r_hi)
                            for r_lo, r_hi in zip(r_edges_avg[:-1], r_edges_avg[1:])], dtype=float)
            xw0_valid = xw0[valid]
            ax1.plot(xw0_valid, _pv(ms0[valid]), 'o-', color=species_colors[0],
                     markersize=5, linewidth=1.5,
                     label=f'<{species_labels[0]}>{abs_suffix} por cascara')
        elif show_species0 and len(mean_species) >= 1:
            ms0 = mean_species[0]
            ax1.plot(r_mids[valid], _pv(ms0[valid]), 'o-', color=species_colors[0],
                     markersize=5, linewidth=1.5,
                     label=f'<{species_labels[0]}>{abs_suffix} por cascara')

        if show_species1 and len(mean_species) >= 2 and distances is not None:
            ms1 = mean_species[1]
            r_edges_avg = np.arange(0.0, r_mids[-1] + shell_width, shell_width)
            xw1 = np.array([_weighted_r_shell(distances, species_vals[1], r_lo, r_hi)
                            for r_lo, r_hi in zip(r_edges_avg[:-1], r_edges_avg[1:])], dtype=float)
            xw1_valid = xw1[valid]
            ax1.plot(xw1_valid, _pv(ms1[valid]), 's-', color=species_colors[1],
                     markersize=5, linewidth=1.5,
                     label=f'<{species_labels[1]}>{abs_suffix} por cascara')
        elif show_species1 and len(mean_species) >= 2:
            ms1 = mean_species[1]
            ax1.plot(r_mids[valid], _pv(ms1[valid]), 's-', color=species_colors[1],
                     markersize=5, linewidth=1.5,
                     label=f'<{species_labels[1]}>{abs_suffix} por cascara')

        if show_net and distances is not None:
            mn = mean_net[valid]
            sn = std_net[valid]
            r_edges_avg = np.arange(0.0, r_mids[-1] + shell_width, shell_width)
            xwn = np.array([_weighted_r_shell(distances, net_vals, r_lo, r_hi)
                            for r_lo, r_hi in zip(r_edges_avg[:-1], r_edges_avg[1:])], dtype=float)
            xwn_valid = xwn[valid]
            ax1.plot(xwn_valid, _pv(mn), 'D-', color='green',
                     markersize=5, linewidth=1.5,
                     label=f'<ρ_net>{abs_suffix} por cascara (Δr={shell_width:.2f})')
            if not use_log:
                ax1.fill_between(xwn_valid, mn - sn, mn + sn,
                                 color='green', alpha=0.15, label='±σ ρ_net')
        elif show_net:
            mn = mean_net[valid]
            sn = std_net[valid]
            rm = r_mids[valid]
            ax1.plot(rm, _pv(mn), 'D-', color='green',
                     markersize=5, linewidth=1.5,
                     label=f'<ρ_net>{abs_suffix} por cascara (Δr={shell_width:.2f})')
            if not use_log:
                ax1.fill_between(rm, mn - sn, mn + sn,
                                 color='green', alpha=0.15, label='±σ ρ_net')

        ylabel_mode = 'promedio por cascara'

    # Curvas teóricas (comunes a ambos modos)
    if show_dh and kappa > 0:
        label_dh = f'Teoría DH: ρ_net{abs_suffix} (λ_D={1.0/kappa:.2f}, κ={kappa:.4f})'
        lw_dh = 1.0 if scatter else 2.0
        ax1.plot(r_theory, _pv(rho_dh_inf_th), 'k-', linewidth=lw_dh, alpha=0.9,
                 label=label_dh)
    if rho_dh_rad_th is not None:
        lw_dh = 1.0 if scatter else 2.0
        ax1.plot(r_theory, _pv(rho_dh_rad_th), 'b--', linewidth=lw_dh, alpha=0.9,
                 label=f'DH radio finito{abs_suffix} (a={particle_radius:.2f} lu)')
        ax1.axvline(x=particle_radius, color='blue', linestyle=':', linewidth=1.0,
                    alpha=0.6, label=f'superficie (a={particle_radius:.2f})')
    if show_dh_periodic and rho_dh_per_th is not None:
        n_img = (2 * dh_periodic_shells + 1) ** 3 - 1
        ax1.plot(r_theory, _pv(rho_dh_per_th), 'r--', linewidth=1.5, alpha=0.9,
                 label=f'DH periódico{abs_suffix} ({dh_periodic_shells} capa(s), {n_img} imgs)')

    # CHANGE INIT - GaussianDH - Curvas gaussiana-DH en panel 1
    if gaussian_source and rho_gauss_dh_th is not None:
        if kappa > 0:
            label_gauss = (f'Gauss-DH ρ_ions{abs_suffix} (σ={sigma:.2f}, '
                           f'λ_D={1.0/kappa:.2f}, κ={kappa:.3f})')
        else:
            label_gauss = f'Gauss-Coulomb ρ{abs_suffix} (σ={sigma:.2f})'
        ax1.plot(r_theory, _pv(rho_gauss_dh_th), 'c-', linewidth=2.0, alpha=0.9,
                 label=label_gauss)
    if gaussian_source and rho_gauss_dh_per_th is not None:
        n_img = (2 * dh_periodic_shells + 1) ** 3 - 1
        ax1.plot(r_theory, _pv(rho_gauss_dh_per_th), 'c--', linewidth=1.5, alpha=0.9,
                 label=(f'Gauss-DH periódico ρ_ions{abs_suffix} ({dh_periodic_shells} '
                        f'cap., {n_img} imgs)'))
    # CHANGE END - GaussianDH

    if kappa > 0:
        ax1.axvline(x=1.0 / kappa, color='purple', linestyle=':', linewidth=2, alpha=0.7,
                    label=f'λ_D={1.0/kappa:.2f}')

    ylabel_base = '|Densidad de carga|' if use_log else 'Densidad de carga'
    ax1.set_xlabel('r (lu)', fontsize=font_size)
    ax1.set_ylabel(f'{ylabel_base} — {ylabel_mode} (adim.)', fontsize=font_size)
    title = (f'Distribución radial de carga (qsi) — Δr={shell_width:.2f}  κ={kappa:.4f}\n'
             f'L={grid_size[0]}×{grid_size[1]}×{grid_size[2]}, q={charge:.2e}, '
             f'ε={epsilon:.2e}, kT={kt:.2e}')
    ax1.set_title(title, fontsize=font_size)
    ax1.legend(fontsize=font_size - 3, loc='best')
    ax1.grid(True, alpha=0.3)

    if log_log:
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        x_left = xmin if xmin is not None else max(rmin, 0.1)
        ax1.set_xlim(x_left, r_max_data)
        if ymin is not None:
            ax1.set_ylim(bottom=ymin)
    elif log_scale:
        ax1.set_yscale('log')
        ax1.set_xlim(xmin if xmin is not None else 0.0, r_max_data)
        if ymin is not None:
            ax1.set_ylim(bottom=ymin)
    else:
        ax1.set_xlim(xmin if xmin is not None else 0.0, r_max_data)
        if ymin is not None:
            ax1.set_ylim(bottom=ymin)

    # Panel 2: error relativo vs DH
    if show_panel2:
        if show_error_inf and not scatter:
            rel_gap = (np.abs(mean_net - rho_dh_inf_pts)
                       / (np.abs(rho_dh_inf_pts) + 1e-15) * 100)
            ax2.plot(r_mids[valid], rel_gap[valid], 'purple',
                     markersize=5, alpha=0.8, marker='D', linestyle='-',
                     label='|<ρ_net> - DH inf| / DH inf × 100')
            if rho_dh_rad_pts is not None:
                valid_rad = valid & ~np.isnan(rho_dh_rad_pts)
                rel_gap_rad = (np.abs(mean_net - rho_dh_rad_pts)
                               / (np.abs(rho_dh_rad_pts) + 1e-15) * 100)
                ax2.plot(r_mids[valid_rad], rel_gap_rad[valid_rad], 'b',
                         markersize=5, alpha=0.8, marker='s', linestyle='-',
                         label=f'|<ρ_net> - DH a={particle_radius:.2f}| / DH_a × 100')
        if rho_dh_per_pts is not None and net_vals is not None:
            n_img = (2 * dh_periodic_shells + 1) ** 3 - 1
            rel_gap_per = (np.abs(net_vals - rho_dh_per_pts)
                           / (np.abs(rho_dh_per_pts) + 1e-15) * 100)
            ax2.plot(distances, rel_gap_per, 'tomato',
                     markersize=3, alpha=0.6, marker='.', linestyle='none',
                     label=f'|ρ_net - DH per exacto| / DH_per × 100 ({n_img} imgs)')

        # CHANGE INIT - GaussianDH - Error vs gaussiana-DH en panel 2
        if gaussian_source and gaussian_error and rho_gauss_dh_pts is not None:
            rel_gap_gauss = (np.abs(mean_net - rho_gauss_dh_pts)
                             / (np.abs(rho_gauss_dh_pts) + 1e-15) * 100)
            if kappa > 0:
                label_ge = f'|<ρ_net> - Gauss-DH (σ={sigma:.2f})| / Gauss-DH × 100'
            else:
                label_ge = f'|<ρ_net> - Gauss-Coulomb (σ={sigma:.2f})| / Gauss-Coulomb × 100'
            ax2.plot(r_mids[valid], rel_gap_gauss[valid], 'teal',
                     markersize=5, alpha=0.8, marker='D', linestyle='-', label=label_ge)
        if gaussian_source and gaussian_error_periodic and rho_gauss_dh_per_pts is not None and net_vals is not None:
            n_img = (2 * dh_periodic_shells + 1) ** 3 - 1
            rel_gap_gauss_per = (np.abs(net_vals - rho_gauss_dh_per_pts)
                                 / (np.abs(rho_gauss_dh_per_pts) + 1e-15) * 100)
            ax2.plot(distances, rel_gap_gauss_per, 'darkcyan',
                     markersize=3, alpha=0.6, marker='.', linestyle='none',
                     label=(f'|ρ_net - Gauss-DH per (σ={sigma:.2f})| / Gauss-DH_per × 100'
                            f' ({n_img} imgs)'))
        # CHANGE END - GaussianDH

        ax2.axhline(y=10, color='r', linestyle=':', alpha=0.5, linewidth=2, label='10%')
        ax2.axhline(y=1, color='g', linestyle=':', alpha=0.5, linewidth=2, label='1%')
        ax2.set_xlabel('r (lu)', fontsize=font_size)
        ax2.set_ylabel('Error relativo (%)', fontsize=font_size)
        ax2.set_title('Error relativo ρ_net vs. teoría DH', fontsize=font_size)
        ax2.legend(fontsize=font_size - 4, loc='upper right')
        ax2.grid(True, alpha=0.3)
        if error_ymax is not None:
            ax2.set_ylim(top=error_ymax)
        ax2.set_xlim(ax1.get_xlim())
    else:
        ax2.set_visible(False)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nGráfico radial guardado en: {output_file}")
    plt.close()


def export_csv_radial(r_mids, mean_species, mean_net, std_net, n_nodes,
                       charge, kappa, csv_file):
    """Exporta el promedio radial a CSV con comparación DH."""

    def rho_dh_scalar(dist):
        dist = max(float(dist), 1e-10)
        if kappa <= 0:
            return 0.0
        return -kappa**2 * charge / (4.0 * np.pi * dist) * np.exp(-kappa * dist)

    num_species = len(mean_species)
    species_header = ','.join(f'mean_rho_{s}' for s in range(num_species))
    header = f"r_mid,n_nodes,{species_header},mean_rho_net,std_rho_net,rho_dh_inf,err_inf%"

    with open(csv_file, 'w') as f:
        f.write(header + '\n')
        for i, r in enumerate(r_mids):
            if np.isnan(mean_net[i]):
                continue
            rho_inf = rho_dh_scalar(r)
            err_inf = abs(mean_net[i] - rho_inf) / (abs(rho_inf) + 1e-15) * 100.0
            species_str = ','.join(f"{mean_species[s][i]:.8e}" for s in range(num_species))
            f.write(f"{r:.8e},{n_nodes[i]:d},{species_str},"
                    f"{mean_net[i]:.8e},{std_net[i]:.8e},"
                    f"{rho_inf:.8e},{err_inf:.4f}\n")

    print(f"CSV radial exportado: {csv_file} ({np.sum(~np.isnan(mean_net))} cascaras)")


# ---------------------------------------------------------------------------
# Gráfico principal
# ---------------------------------------------------------------------------

def plot_qsi(distances, species_vals, net_vals, r_vecs,
             grid_size, charge, epsilon, kt, kappa,
             rmin=0.5, rmax=None, output_file='qsi_ewald.png', max_points=10000,
             log_scale=False, log_log=False, log_log_error=False,
             xmin=None, ymin=None,
             show_species0=True, show_species1=True, show_net=True,
             show_dh=True, show_dh_periodic=False, dh_periodic_shells=1,
             show_dh_periodic_pts=False,
             direction=None, angle_tol_deg=15.0,
             pointwise_error=False,
             show_error_inf=True, show_error_dir=True,
             show_error=False, error_ymax=None,
             particle_radius=0.0,
             only_error_dh_per_rad=False,
             gaussian_source=False, sigma=1.0,
             gaussian_error=False, gaussian_error_periodic=False):

    L = float(grid_size[0])
    num_species = len(species_vals)

    # En modo log se usa valor absoluto (rho_net es negativa cerca de carga +)
    use_log = log_scale or log_log
    abs_suffix = ' |·|' if use_log else ''

    # Preparar arrays para plot
    arrays_to_sort = list(species_vals) + [net_vals, r_vecs]
    sorted_out = _sort_subsample(distances, *arrays_to_sort,
                                  max_pts=max_points, rmin=rmin)
    d = sorted_out[0]
    sorted_species = list(sorted_out[1:1 + num_species])
    net = sorted_out[1 + num_species]
    rv = sorted_out[2 + num_species]

    r_max_data = float(np.max(d)) * 1.1 if len(d) > 0 else 16.0
    # CHANGE INIT - RmaxTheory - usar rmax del usuario y evitar singularidades periódicas
    # np.linspace puede caer exactamente en n*L donde la imagen periódica está a dist~0,
    # generando picos de ~1e8 que dominan el autoscale del eje Y.
    if rmax is not None:
        r_max_theory = rmax
    else:
        r_max_theory = r_max_data
    r_theory = np.linspace(rmin, r_max_theory, 500)
    if show_dh_periodic:
        safe = np.ones(len(r_theory), dtype=bool)
        for n in range(1, int(r_max_theory / L) + 2):
            safe &= np.abs(r_theory - n * L) > 0.1
        r_theory = r_theory[safe]
    # CHANGE END - RmaxTheory

    # Curvas teóricas
    rho_dh_inf = TheoreticalChargeDensity.compute_dh(
        r_theory, q=charge, epsilon=epsilon, kt=kt, kappa=kappa)

    rho_dh_rad = None
    if particle_radius > 0 and kappa > 0:
        rho_dh_rad = TheoreticalChargeDensity.compute_dh_finite_radius(
            r_theory, q=charge, kappa=kappa, a=particle_radius)

    rho_dh_per = None
    if show_dh_periodic:
        print(f"  Calculando ρ_net DH periódico: {dh_periodic_shells} capa(s), "
              f"dir={_direction_label(direction)}...")
        rho_dh_per = TheoreticalChargeDensity.compute_dh_periodic(
            r_theory, q=charge, epsilon=epsilon, kt=kt, kappa=kappa,
            L=L, n_shells=dh_periodic_shells, direction=direction)

    # Periódico en posiciones 3D exactas de los nodos (subsampled)
    rho_dh_per_pts = None
    if show_dh_periodic_pts and kappa > 0:
        n_img = (2 * dh_periodic_shells + 1) ** 3 - 1
        print(f"  Calculando ρ_net DH periódico en posiciones exactas: "
              f"{dh_periodic_shells} capa(s), {n_img} imágenes, {len(d)} nodos...")
        rho_dh_per_pts = TheoreticalChargeDensity.compute_dh_periodic_at_points(
            rv, q=charge, kappa=kappa, L=L, n_shells=dh_periodic_shells)

    # CHANGE INIT - OnlyErrorDHPerRad - Pre-calcular DH periódico con radio finito
    rho_dh_per_rad = None
    if only_error_dh_per_rad and particle_radius > 0 and show_dh_periodic and kappa > 0:
        n_img = (2 * dh_periodic_shells + 1) ** 3 - 1
        print(f"  Calculando ρ_net DH periódico radio finito (a={particle_radius:.2f}): "
              f"{dh_periodic_shells} capa(s), {n_img} imágenes...")
        rho_dh_per_rad = TheoreticalChargeDensity.compute_dh_periodic_finite_radius(
            d, q=charge, kappa=kappa, a=particle_radius,
            L=L, n_shells=dh_periodic_shells, direction=direction)
    # CHANGE END - OnlyErrorDHPerRad

    # CHANGE INIT - GaussianDH - Cálculo de curvas gaussiana-DH en plot_qsi
    rho_gauss_dh = None
    rho_gauss_dh_per = None
    rho_gauss_dh_per_pts = None
    if gaussian_source:
        print(f"  Calculando ρ_ions gaussiana-DH (σ={sigma:.3f}) en r_theory...")
        rho_gauss_dh = TheoreticalChargeDensity.compute_gaussian_dh(
            r_theory, q=charge, kappa=kappa, sigma=sigma)
        if show_dh_periodic:
            n_img = (2 * dh_periodic_shells + 1) ** 3 - 1
            print(f"  Calculando ρ_ions gaussiana-DH periódico ({dh_periodic_shells} capa(s), "
                  f"{n_img} imgs)...")
            rho_gauss_dh_per = TheoreticalChargeDensity.compute_gaussian_dh_periodic(
                r_theory, q=charge, kappa=kappa, sigma=sigma,
                L=L, n_shells=dh_periodic_shells, direction=direction)
        if gaussian_error_periodic and show_dh_periodic and rv is not None and len(rv) > 0 and kappa > 0:
            n_img = (2 * dh_periodic_shells + 1) ** 3 - 1
            print(f"  Calculando ρ_ions gaussiana-DH periódico en posiciones exactas: "
                  f"{dh_periodic_shells} capa(s), {n_img} imgs, {len(rv)} nodos...")
            rho_gauss_dh_per_pts = TheoreticalChargeDensity.compute_gaussian_dh_periodic_at_points(
                rv, q=charge, kappa=kappa, sigma=sigma, L=L, n_shells=dh_periodic_shells)
    # CHANGE END - GaussianDH

    # Determinar si mostrar panel de error
    show_panel2 = show_net and (
        (show_dh and show_error_inf and kappa > 0) or
        (show_dh_periodic and show_error_dir and kappa > 0) or
        pointwise_error or
        (show_error and kappa > 0) or
        (show_dh_periodic_pts and kappa > 0) or
        (only_error_dh_per_rad and rho_dh_per_rad is not None) or
        (gaussian_source and gaussian_error) or
        (gaussian_source and gaussian_error_periodic and show_dh_periodic)
    )

    _, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 11),
                                  gridspec_kw={'height_ratios': [2, 1]})

    dir_label_str = (f" dir={_direction_label(direction)} ±{angle_tol_deg:.0f}°"
                     if direction is not None else "")

    # Colores de las especies
    species_colors = ['tab:blue', 'tab:red', 'tab:green', 'tab:orange']
    species_labels = ['ρ₊ (cationes)', 'ρ₋ (aniones)', 'ρ₂', 'ρ₃']

    def _plot_vals(ax, x, y, *args, **kwargs):
        """Grafica y o |y| según si estamos en modo log."""
        ax.plot(x, np.abs(y) if use_log else y, *args, **kwargs)

    # Panel 1: densidades vs r
    if show_species0 and num_species >= 1 and len(sorted_species[0]) > 0:
        _plot_vals(ax1, d, sorted_species[0], '.', color=species_colors[0],
                   markersize=4, alpha=0.6,
                   label=f'Ludwig: {species_labels[0]}{abs_suffix}')
    if show_species1 and num_species >= 2 and len(sorted_species[1]) > 0:
        _plot_vals(ax1, d, sorted_species[1], '.', color=species_colors[1],
                   markersize=4, alpha=0.6,
                   label=f'Ludwig: {species_labels[1]}{abs_suffix}')
    if show_net and len(net) > 0:
        _plot_vals(ax1, d, net, 'g.', markersize=6, alpha=0.7,
                   label=f'Ludwig: ρ_net = ρ₊ - ρ₋{abs_suffix}{dir_label_str}')

    if show_dh and kappa > 0:
        label_dh = f'Teoría DH: ρ_net{abs_suffix} (λ_D={1.0/kappa:.2f}, κ={kappa:.4f})'
        _plot_vals(ax1, r_theory, rho_dh_inf, 'k-', linewidth=1.5, alpha=0.8,
                   label=label_dh)
    if rho_dh_rad is not None:
        _plot_vals(ax1, r_theory, rho_dh_rad, 'b--', linewidth=1.5, alpha=0.9,
                   label=f'DH radio finito{abs_suffix} (a={particle_radius:.2f} lu)')
        ax1.axvline(x=particle_radius, color='blue', linestyle=':', linewidth=1.0,
                    alpha=0.6, label=f'superficie (a={particle_radius:.2f})')
    if show_dh_periodic and rho_dh_per is not None:
        n_img = (2 * dh_periodic_shells + 1) ** 3 - 1
        _plot_vals(ax1, r_theory, rho_dh_per, 'r--', linewidth=1.5, alpha=0.9,
                   label=f'DH periódico{abs_suffix} ({dh_periodic_shells} capa(s), {n_img} imgs)')
    if rho_dh_per_pts is not None:
        n_img = (2 * dh_periodic_shells + 1) ** 3 - 1
        _plot_vals(ax1, d, rho_dh_per_pts, 'r.', markersize=3, alpha=0.5,
                   label=f'DH per. exacto{abs_suffix} ({n_img} imgs, por nodo)')

    # CHANGE INIT - GaussianDH - Curvas gaussiana-DH en panel 1 de plot_qsi
    if gaussian_source and rho_gauss_dh is not None:
        if kappa > 0:
            label_gauss = (f'Gauss-DH ρ_ions{abs_suffix} (σ={sigma:.2f}, '
                           f'λ_D={1.0/kappa:.2f}, κ={kappa:.3f})')
        else:
            label_gauss = f'Gauss-Coulomb ρ{abs_suffix} (σ={sigma:.2f})'
        _plot_vals(ax1, r_theory, rho_gauss_dh, 'c-', linewidth=2.0, alpha=0.9,
                   label=label_gauss)
    if gaussian_source and rho_gauss_dh_per is not None:
        n_img = (2 * dh_periodic_shells + 1) ** 3 - 1
        _plot_vals(ax1, r_theory, rho_gauss_dh_per, 'c--', linewidth=1.5, alpha=0.9,
                   label=(f'Gauss-DH periódico ρ_ions{abs_suffix} ({dh_periodic_shells} '
                          f'cap., {n_img} imgs)'))
    # CHANGE END - GaussianDH

    if kappa > 0:
        ax1.axvline(x=1.0 / kappa, color='purple', linestyle=':', linewidth=2, alpha=0.7,
                    label=f'λ_D={1.0/kappa:.2f}')

    ylabel = ('|Densidad de carga| (adimensional)' if use_log
              else 'Densidad de carga (adimensional)')
    ax1.set_xlabel('r (lu)', fontsize=font_size)
    ax1.set_ylabel(ylabel, fontsize=font_size)
    title = (f'Densidad de carga eléctrica (qsi) Ewald  κ={kappa:.4f}\n'
             f'L={grid_size[0]}×{grid_size[1]}×{grid_size[2]}, q={charge:.2e}, '
             f'ε={epsilon:.2e}, kT={kt:.2e}')
    if direction is not None:
        title += f'\nDirección: {_direction_label(direction)} (±{angle_tol_deg:.0f}°)'
    ax1.set_title(title, fontsize=font_size)
    ax1.legend(fontsize=font_size - 3, loc='best')
    ax1.grid(True, alpha=0.3)

    # Escala
    if log_log:
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        x_left = xmin if xmin is not None else max(rmin, 0.1)
        ax1.set_xlim(x_left, r_max_data)
        if ymin is not None:
            ax1.set_ylim(bottom=ymin)
    elif log_scale:
        ax1.set_yscale('log')
        ax1.set_xlim(xmin if xmin is not None else 0.0, r_max_data)
        if ymin is not None:
            ax1.set_ylim(bottom=ymin)
    else:
        ax1.set_xlim(xmin if xmin is not None else 0.0, r_max_data)
        if ymin is not None:
            ax1.set_ylim(bottom=ymin)

    # Panel 2: error relativo de ρ_net vs DH
    if show_panel2:
        # CHANGE INIT - OnlyErrorDHPerRad - Suprimir otras curvas si only_error_dh_per_rad
        if not only_error_dh_per_rad:
            if show_dh and show_error_inf and kappa > 0:
                rho_inf_at_pts = TheoreticalChargeDensity.compute_dh(
                    d, q=charge, epsilon=epsilon, kt=kt, kappa=kappa)
                rel_gap_inf = (np.abs(net - rho_inf_at_pts)
                               / (np.abs(rho_inf_at_pts) + 1e-15) * 100)
                ax2.plot(d, rel_gap_inf, 'purple', markersize=3.5, alpha=0.6,
                         marker='.', linestyle='none',
                         label='|ρ_net - DH inf| / DH inf × 100')
            if particle_radius > 0 and kappa > 0:
                rho_rad_at_pts = TheoreticalChargeDensity.compute_dh_finite_radius(
                    d, q=charge, kappa=kappa, a=particle_radius)
                valid_rad = ~np.isnan(rho_rad_at_pts)
                rel_gap_rad = (np.abs(net[valid_rad] - rho_rad_at_pts[valid_rad])
                               / (np.abs(rho_rad_at_pts[valid_rad]) + 1e-15) * 100)
                ax2.plot(d[valid_rad], rel_gap_rad, 'b', markersize=3.5, alpha=0.7,
                         marker='.', linestyle='none',
                         label=f'|ρ_net - DH a={particle_radius:.2f}| / DH_a × 100')

            if show_dh_periodic and show_error_dir and rho_dh_per is not None and kappa > 0:
                rho_per_at_pts = TheoreticalChargeDensity.compute_dh_periodic(
                    d, q=charge, epsilon=epsilon, kt=kt, kappa=kappa,
                    L=L, n_shells=dh_periodic_shells, direction=direction)
                rel_gap_per = (np.abs(net - rho_per_at_pts)
                               / (np.abs(rho_per_at_pts) + 1e-15) * 100)
                ax2.plot(d, rel_gap_per, 'darkorange', markersize=3.5, alpha=0.7,
                         marker='.', linestyle='none',
                         label='|ρ_net - DH per (dir)| / DH × 100')

            if rho_dh_per_pts is not None:
                rel_gap_per_pts = (np.abs(net - rho_dh_per_pts)
                                   / (np.abs(rho_dh_per_pts) + 1e-15) * 100)
                n_img = (2 * dh_periodic_shells + 1) ** 3 - 1
                ax2.plot(d, rel_gap_per_pts, 'tomato', markersize=3.5, alpha=0.7,
                         marker='.', linestyle='none',
                         label=f'|ρ_net - DH per exacto| / DH_per × 100 ({n_img} imgs)')

            if pointwise_error and rv is not None and len(rv) > 0 and kappa > 0:
                print(f"  Calculando error puntual en {len(rv)} nodos...")
                rho_pw = TheoreticalChargeDensity.compute_dh_periodic_at_points(
                    rv, q=charge, epsilon=epsilon, kt=kt, kappa=kappa,
                    L=L, n_shells=dh_periodic_shells)
                rel_gap_pw = (np.abs(net - rho_pw)
                              / (np.abs(rho_pw) + 1e-15) * 100)
                ax2.plot(d, rel_gap_pw, 'green', markersize=5, alpha=0.8,
                         marker='.', linestyle='none',
                         label='|ρ_net - DH per (punto exacto)| / DH × 100')

        if rho_dh_per_rad is not None:
            n_img = (2 * dh_periodic_shells + 1) ** 3 - 1
            valid_per_rad = ~np.isnan(rho_dh_per_rad)
            rel_gap_per_rad = (np.abs(net[valid_per_rad] - rho_dh_per_rad[valid_per_rad])
                               / (np.abs(rho_dh_per_rad[valid_per_rad]) + 1e-15) * 100)
            ax2.plot(d[valid_per_rad], rel_gap_per_rad, 'darkgreen',
                     markersize=3.5, alpha=0.8,
                     marker='.', linestyle='none',
                     label=(f'|ρ_net - DH per a={particle_radius:.2f}| / DH_per_a × 100'
                            f' ({dh_periodic_shells} capa(s), {n_img} imgs)'))
        # CHANGE END - OnlyErrorDHPerRad

        # CHANGE INIT - GaussianDH - Error vs gaussiana-DH en panel 2 de plot_qsi
        if gaussian_source and gaussian_error:
            rho_gauss_at_pts = TheoreticalChargeDensity.compute_gaussian_dh(
                d, q=charge, kappa=kappa, sigma=sigma)
            rel_gap_gauss = (np.abs(net - rho_gauss_at_pts)
                             / (np.abs(rho_gauss_at_pts) + 1e-15) * 100)
            if kappa > 0:
                label_ge = f'|ρ_net - Gauss-DH (σ={sigma:.2f})| / Gauss-DH × 100'
            else:
                label_ge = f'|ρ_net - Gauss-Coulomb (σ={sigma:.2f})| / Gauss-Coulomb × 100'
            ax2.plot(d, rel_gap_gauss, 'teal', markersize=3.5, alpha=0.8,
                     marker='.', linestyle='none', label=label_ge)

        if gaussian_source and gaussian_error_periodic and rho_gauss_dh_per_pts is not None:
            rel_gap_gauss_per = (np.abs(net - rho_gauss_dh_per_pts)
                                 / (np.abs(rho_gauss_dh_per_pts) + 1e-15) * 100)
            n_img = (2 * dh_periodic_shells + 1) ** 3 - 1
            ax2.plot(d, rel_gap_gauss_per, 'darkcyan', markersize=4, alpha=0.85,
                     marker='.', linestyle='none',
                     label=(f'|ρ_net - Gauss-DH per (σ={sigma:.2f})| / Gauss-DH_per × 100'
                            f' ({dh_periodic_shells} cap., {n_img} imgs)'))
        # CHANGE END - GaussianDH

        ax2.axhline(y=10, color='r', linestyle=':', alpha=0.5, linewidth=2,
                    label='10%')
        ax2.axhline(y=1, color='g', linestyle=':', alpha=0.5, linewidth=2,
                    label='1%')
        ax2.set_xlabel('r (lu)', fontsize=font_size)
        ax2.set_ylabel('Error relativo (%)', fontsize=font_size)
        ax2.set_title('Error relativo ρ_net vs. teoría DH', fontsize=font_size)
        ax2.legend(fontsize=font_size - 4, loc='upper right')
        ax2.grid(True, alpha=0.3)
        if log_log_error:
            ax2.set_yscale('log')
        if error_ymax is not None:
            ax2.set_ylim(top=error_ymax)
        ax2.set_xlim(ax1.get_xlim())
    else:
        ax2.set_visible(False)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nGráfico guardado en: {output_file}")
    plt.close()


# ---------------------------------------------------------------------------
# Exportar CSV
# ---------------------------------------------------------------------------

def export_csv_qsi(species_3d_list, particle_pos, grid_size,
                   charge, epsilon, kt, kappa, L, n_shells,
                   rmin, csv_file):
    """
    Exporta todos los nodos (distancia >= rmin) con:
      - Posición de la partícula (px, py, pz)
      - Índice del nodo (ix, iy, iz) y posición absoluta (x, y, z)
      - Distancia a la partícula
      - Densidad de cada especie (rho_0, rho_1, ...)
      - Carga neta: rho_net = rho_0 - rho_1
      - Teoría DH infinita: rho_dh_inf
      - Teoría DH periódica: rho_dh_per
      - Error relativo vs DH inf y periódico (%)
    """
    nx, ny, nz = grid_size
    px, py, pz = particle_pos
    num_species = len(species_3d_list)

    def rho_dh_scalar(dist):
        dist = max(float(dist), 1e-10)
        if kappa <= 0:
            return 0.0
        return -kappa**2 * charge / (4.0 * np.pi * dist) * np.exp(-kappa * dist)

    def rho_dh_periodic(r_vec):
        dist0 = np.linalg.norm(r_vec)
        val = rho_dh_scalar(dist0)
        for inx in range(-n_shells, n_shells + 1):
            for iny in range(-n_shells, n_shells + 1):
                for inz in range(-n_shells, n_shells + 1):
                    if inx == 0 and iny == 0 and inz == 0:
                        continue
                    r_img = np.array([inx * L, iny * L, inz * L])
                    dr = r_vec - r_img
                    d = np.linalg.norm(dr)
                    if d > 1e-10:
                        val += rho_dh_scalar(d)
        return val

    def rel_err(val, ref):
        return abs(val - ref) / (abs(ref) + 1e-15) * 100.0

    rows = []
    for i in range(nx):
        xi = i + 1
        for j in range(ny):
            yj = j + 1
            for k in range(nz):
                zk = k + 1
                r_vec = np.array([xi - px, yj - py, zk - pz])
                dist = np.linalg.norm(r_vec)
                if dist < rmin:
                    continue

                rho_species = [species_3d_list[s][i, j, k]
                               for s in range(num_species)]

                if num_species >= 2:
                    rho_net = rho_species[0] - rho_species[1]
                else:
                    rho_net = rho_species[0]

                rho_inf = rho_dh_scalar(dist)
                rho_per = rho_dh_periodic(r_vec)

                rows.append((dist, i, j, k, xi, yj, zk,
                             *rho_species, rho_net,
                             rho_inf, rho_per,
                             rel_err(rho_net, rho_inf),
                             rel_err(rho_net, rho_per)))

    rows.sort(key=lambda r: r[0])

    # Cabecera dinámica según número de especies
    species_header = ','.join(f'rho_{s}' for s in range(num_species))
    header = (f"px,py,pz,dist,ix,iy,iz,x,y,z,"
              f"{species_header},rho_net,"
              f"rho_dh_inf,rho_dh_per,"
              f"err_inf%,err_per%")

    with open(csv_file, 'w') as f:
        f.write(header + '\n')
        for r in rows:
            # Primeros 10 valores fijos, luego num_species densidades, luego 5 más
            base = r[:7]       # (dist, i, j, k, xi, yj, zk)
            rho_s = r[7:7 + num_species]  # densidades por especie
            tail = r[7 + num_species:]    # rho_net, rho_inf, rho_per, err_inf, err_per

            line = (f"{px:.6f},{py:.6f},{pz:.6f},"
                    f"{base[0]:.8e},{base[1]:d},{base[2]:d},{base[3]:d},"
                    f"{base[4]:.6f},{base[5]:.6f},{base[6]:.6f},")
            line += ','.join(f"{v:.8e}" for v in rho_s) + ','
            line += ','.join(f"{v:.8e}" for v in tail[:3]) + ','
            line += f"{tail[3]:.4f},{tail[4]:.4f}\n"
            f.write(line)

    print(f"CSV exportado: {csv_file} ({len(rows)} puntos)")


def export_csv_scatter(distances, species_vals, net_vals, r_vecs,
                       particle_pos, charge, epsilon, kappa,
                       csv_file):
    """
    Exporta todos los nodos en modo scatter con:
      - Posición de la partícula (px, py, pz)
      - Posición del nodo (x, y, z)
      - Distancia a la partícula (r)
      - Densidad de cada especie (rho_0, rho_1, ...)
      - Carga neta calculada (rho_net)
      - Carga teórica DH en ese nodo (rho_dh)
      - Fuerza radial calculada sobre la partícula: F_calc = q*rho_net / (4pi*eps*r^2)
      - Fuerza radial teórica DH:                  F_dh   = q*rho_dh  / (4pi*eps*r^2)
    """
    ppx, ppy, ppz = particle_pos
    num_species = len(species_vals)

    def _rho_dh(r):
        if kappa <= 0 or r < 1e-10:
            return 0.0
        return -kappa**2 * charge / (4.0 * np.pi * r) * np.exp(-kappa * r)

    def _force(rho, r):
        if r < 1e-10:
            return 0.0
        return charge * rho / (4.0 * np.pi * epsilon * r**2)

    species_header = ','.join(f'rho_{s}' for s in range(num_species))
    header = (f"px,py,pz,x,y,z,r,{species_header},rho_net,rho_dh,"
              f"F_calc_x,F_calc_y,F_calc_z,F_calc,"
              f"F_dh_x,F_dh_y,F_dh_z,F_dh")

    n = len(distances)
    sorted_idx = np.argsort(distances)

    with open(csv_file, 'w') as f:
        f.write(header + '\n')
        for idx in sorted_idx:
            r = distances[idx]
            rv = r_vecs[idx]
            x = ppx + rv[0]
            y = ppy + rv[1]
            z = ppz + rv[2]
            rho_net = net_vals[idx]
            rho_dh = _rho_dh(r)
            f_calc_mag = _force(rho_net, r)
            f_dh_mag = _force(rho_dh, r)
            r_hat = rv / r if r > 1e-10 else np.zeros(3)
            fc_vec = f_calc_mag * r_hat
            fd_vec = f_dh_mag * r_hat
            rho_s = [species_vals[s][idx] for s in range(num_species)]
            line = (f"{ppx:.16e},{ppy:.16e},{ppz:.16e},"
                    f"{x:.16e},{y:.16e},{z:.16e},"
                    f"{r:.16e},")
            line += ','.join(f"{v:.16e}" for v in rho_s) + ','
            line += (f"{rho_net:.16e},{rho_dh:.16e},"
                     f"{fc_vec[0]:.16e},{fc_vec[1]:.16e},{fc_vec[2]:.16e},{f_calc_mag:.16e},"
                     f"{fd_vec[0]:.16e},{fd_vec[1]:.16e},{fd_vec[2]:.16e},{f_dh_mag:.16e}\n")
            f.write(line)

    print(f"CSV scatter exportado: {csv_file} ({n} puntos)")


# ---------------------------------------------------------------------------
# Mapa angular en cáscara esférica
# ---------------------------------------------------------------------------

def plot_angular_map(r_vecs, net_vals, r_shell, dr,
                     charge, kappa,
                     output_file='qsi_angular.png',
                     n_theta=90, n_phi=180,
                     quantity='net'):
    """
    Mapa de calor (θ, φ) de la densidad de carga en una cáscara esférica.

    Selecciona nodos con r_shell - dr/2 < |r| < r_shell + dr/2,
    convierte a coordenadas esféricas y hace un mapa de calor interpolado
    en una grilla regular (θ, φ).

    Si el sistema fuera perfectamente esférico el mapa sería uniforme.
    Patrones cúbicos indican rotura de simetría por la red.

    Parámetros
    ----------
    r_vecs    : (N, 3) vectores posición relativa a la partícula
    net_vals  : (N,)   densidad de carga neta en cada nodo
    r_shell   : radio central de la cáscara
    dr        : grosor de la cáscara
    n_theta   : resolución en θ (polar, 0..π)
    n_phi     : resolución en φ (azimutal, 0..2π)
    quantity  : 'net' → carga neta; 'abs' → valor absoluto
    """
    r_norms = np.linalg.norm(r_vecs, axis=1)
    mask = (r_norms >= r_shell - dr / 2.0) & (r_norms <= r_shell + dr / 2.0)

    if np.sum(mask) == 0:
        print(f"AVISO: No hay nodos en la cáscara r={r_shell:.2f} ± {dr/2:.2f}")
        return

    rv = r_vecs[mask]
    rn = r_norms[mask]
    vals = net_vals[mask]
    if quantity == 'abs':
        vals = np.abs(vals)

    print(f"  Nodos en cáscara r={r_shell:.2f} ± {dr/2:.2f}: {np.sum(mask)}")

    # Coordenadas esféricas: θ ∈ [0,π], φ ∈ [0,2π)
    x_hat = rv[:, 0] / rn
    y_hat = rv[:, 1] / rn
    z_hat = rv[:, 2] / rn
    theta = np.arccos(np.clip(z_hat, -1.0, 1.0))          # polar
    phi   = np.arctan2(y_hat, x_hat) % (2.0 * np.pi)      # azimutal

    # Grilla regular
    theta_edges = np.linspace(0, np.pi,      n_theta + 1)
    phi_edges   = np.linspace(0, 2 * np.pi,  n_phi   + 1)
    theta_mid   = 0.5 * (theta_edges[:-1] + theta_edges[1:])
    phi_mid     = 0.5 * (phi_edges[:-1]   + phi_edges[1:])

    # Acumular valores en cada celda (media de los nodos que caen en ella)
    grid_sum   = np.zeros((n_theta, n_phi))
    grid_count = np.zeros((n_theta, n_phi), dtype=int)

    itheta = np.searchsorted(theta_edges[1:], theta, side='left')
    iphi   = np.searchsorted(phi_edges[1:],   phi,   side='left')
    itheta = np.clip(itheta, 0, n_theta - 1)
    iphi   = np.clip(iphi,   0, n_phi   - 1)

    np.add.at(grid_sum,   (itheta, iphi), vals)
    np.add.at(grid_count, (itheta, iphi), 1)

    with np.errstate(invalid='ignore'):
        grid_mean = np.where(grid_count > 0, grid_sum / grid_count, np.nan)

    # Teoría DH (valor escalar esperado a distancia r_shell)
    rho_dh_expected = None
    if kappa > 0:
        rho_dh_expected = TheoreticalChargeDensity.compute_dh(
            np.array([r_shell]), q=charge, kappa=kappa)[0]
        if quantity == 'abs':
            rho_dh_expected = abs(rho_dh_expected)

    # --- Figura ---
    fig, axes = plt.subplots(1, 2, figsize=(16, 6),
                              gridspec_kw={'width_ratios': [3, 1]})
    ax_map, ax_hist = axes

    # Escala de color centrada en el valor teórico si está disponible
    valid_vals = grid_mean[~np.isnan(grid_mean)]
    vmin = np.nanmin(valid_vals) if len(valid_vals) > 0 else -1
    vmax = np.nanmax(valid_vals) if len(valid_vals) > 0 else  1

    # Escala acotada a los valores de los datos
    if vmin < 0 and vmax > 0:
        cmap = 'RdBu_r'
    else:
        cmap = 'viridis'

    # φ en eje X (columnas), θ en eje Y (filas)
    phi_deg   = np.degrees(phi_mid)
    theta_deg = np.degrees(theta_mid)
    label_q = '|ρ_net|' if quantity == 'abs' else 'ρ_net'

    n_nodes = int(np.sum(mask))
    USE_SCATTER = n_nodes <= 20  # pocos nodos: scatter directo

    if USE_SCATTER:
        phi_deg_pts   = np.degrees(phi)
        theta_deg_pts = np.degrees(theta)
        norm = plt.Normalize(vmin=vmin, vmax=vmax)
        sc = ax_map.scatter(phi_deg_pts, theta_deg_pts, c=vals,
                            cmap=cmap, norm=norm, s=120, zorder=5,
                            edgecolors='k', linewidths=0.5)
        cbar = plt.colorbar(sc, ax=ax_map)
        # Anotar valor de cada punto
        for pd, td, v in zip(phi_deg_pts, theta_deg_pts, vals):
            ax_map.annotate(f'{v:.3e}', (pd, td),
                            textcoords='offset points', xytext=(6, 4),
                            fontsize=font_size - 6, color='k')
    else:
        im = ax_map.pcolormesh(phi_deg, theta_deg, grid_mean,
                                cmap=cmap, vmin=vmin, vmax=vmax,
                                shading='auto')
        cbar = plt.colorbar(im, ax=ax_map)
    cbar.set_label(f'{label_q} (adim.)', fontsize=font_size - 2)

    if rho_dh_expected is not None:
        ax_map.set_title(
            f'Mapa angular {label_q}  r={r_shell:.2f}±{dr/2:.2f} lu\n'
            f'κ={kappa:.4f}  ρ_DH(r)={rho_dh_expected:.4e}  '
            f'n_nodos={np.sum(mask)}',
            fontsize=font_size - 1)
    else:
        ax_map.set_title(
            f'Mapa angular {label_q}  r={r_shell:.2f}±{dr/2:.2f} lu  '
            f'n_nodos={np.sum(mask)}',
            fontsize=font_size - 1)

    ax_map.set_xlabel('φ (grados)', fontsize=font_size)
    ax_map.set_ylabel('θ (grados)', fontsize=font_size)
    ax_map.set_xlim(0, 360)
    ax_map.set_ylim(0, 180)
    ax_map.set_xticks([0, 90, 180, 270, 360])
    ax_map.set_yticks([0, 45, 90, 135, 180])

    # Marcar ejes cartesianos (+x, -x, +y, -y, +z, -z) en el mapa
    # +z: θ=0  -z: θ=180  +x: θ=90,φ=0  -x: θ=90,φ=180  +y: θ=90,φ=90  -y: θ=90,φ=270
    cart_axes = {
        '+z': (0,   0  ), '-z': (180, 0  ),
        '+x': (90,  0  ), '-x': (90,  180),
        '+y': (90,  90 ), '-y': (90,  270),
    }
    for label, (th, ph) in cart_axes.items():
        ax_map.plot(ph, th, 'w+', markersize=12, markeredgewidth=2)
        ax_map.text(ph + 4, th + 4, label, color='white',
                    fontsize=font_size - 4, fontweight='bold')

    # Líneas en múltiplos de 45°: patrón distinto para múltiplos de 90°
    # (10, (5, 2, 1, 2)) = guión largo, espacio, guión largo, espacio, punto, espacio
    ls_90 = (0, (8, 3, 2, 3, 2, 3))   # larga-larga-punto
    ls_45 = '--'
    for ph in range(0, 361, 45):
        if ph % 90 == 0:
            ax_map.axvline(ph, color='gray', linestyle=ls_90, linewidth=1.2, alpha=0.8)
        else:
            ax_map.axvline(ph, color='#b0b0b0', linestyle=ls_45, linewidth=0.5, alpha=0.8)
    for th in range(0, 181, 45):
        if th % 90 == 0:
            ax_map.axhline(th, color='gray', linestyle=ls_90, linewidth=1.2, alpha=0.8)
        else:
            ax_map.axhline(th, color='#b0b0b0', linestyle=ls_45, linewidth=0.5, alpha=0.8)

    # Histograma lateral de los valores en la cáscara
    if len(vals) <= 20:
        # Puntos individuales: X = índice, Y = valor
        ax_hist.scatter(np.arange(len(vals)), vals, color='steelblue',
                        edgecolors='k', s=60, zorder=5)
        ax_hist.set_xlabel('índice nodo', fontsize=font_size - 2)
        ax_hist.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
    else:
        n_bins = min(40, max(1, len(vals) // 2))
        ax_hist.hist(vals, bins=n_bins, orientation='horizontal',
                     color='steelblue', edgecolor='k', alpha=0.7)
    if rho_dh_expected is not None:
        ax_hist.axhline(rho_dh_expected, color='red', linewidth=2,
                        linestyle='--', label=f'DH: {rho_dh_expected:.3e}')
        ax_hist.legend(fontsize=font_size - 4)
    if len(vals) > 20:
        ax_hist.set_xlabel('N nodos', fontsize=font_size - 2)
    ax_hist.set_ylabel(f'{label_q}', fontsize=font_size - 2)
    ax_hist.set_title('Distribución', fontsize=font_size - 2)
    ax_hist.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Mapa angular guardado en: {output_file}")
    plt.close()


# ---------------------------------------------------------------------------
# Mapa de plano
# ---------------------------------------------------------------------------

def _parse_plane_normal(plane_str):
    """
    Parsea --plane y devuelve (normal_unitario, label).

    Formatos aceptados:
      'xy'         → normal = (0,0,1)
      'xz'         → normal = (0,1,0)
      'yz'         → normal = (1,0,0)
      'nx,ny,nz'   → normal arbitraria (se normaliza automáticamente)
    El plano siempre pasa por la posición de la partícula.
    """
    p = plane_str.strip().lower()
    presets = {'xy': ([0., 0., 1.], 'xy'),
               'xz': ([0., 1., 0.], 'xz'),
               'yz': ([1., 0., 0.], 'yz')}
    if p in presets:
        n_raw, label = presets[p]
        return np.array(n_raw), label
    parts = p.split(',')
    if len(parts) != 3:
        raise ValueError(
            f"Formato de plano no reconocido: '{plane_str}'. "
            "Usa 'xy','xz','yz' o 'nx,ny,nz'.")
    n = np.array([float(x) for x in parts])
    norm = np.linalg.norm(n)
    if norm < 1e-12:
        raise ValueError("La normal del plano no puede ser el vector cero.")
    n = n / norm
    label = f"n=({n[0]:.2f},{n[1]:.2f},{n[2]:.2f})"
    return n, label


def _build_plane_axes(normal):
    """Construye dos vectores ortonormales (e1, e2) en el plano perpendicular a normal."""
    aux = np.array([1., 0., 0.])
    if abs(np.dot(normal, aux)) > 0.9:
        aux = np.array([0., 1., 0.])
    e1 = np.cross(normal, aux)
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(normal, e1)
    e2 /= np.linalg.norm(e2)
    return e1, e2


def _axis_label(vec):
    """Etiqueta legible para un eje en el plano."""
    names = ['x', 'y', 'z']
    nonzero = [i for i in range(3) if abs(vec[i]) > 0.05]
    if len(nonzero) == 1:
        # Eje cartesiano puro
        return f'{names[nonzero[0]]} (lu)'
    # Eje mixto: mostrar qué combinación de ejes cartesianos representa
    parts = [f'{vec[i]:+.2f}{names[i]}' for i in nonzero]
    return f'{"".join(parts)} (lu)'


def plot_plane_map(r_vecs, net_vals, plane, dz,
                   charge, kappa,
                   output_file='qsi_plane.png',
                   plot_type='pcolor',
                   quantity='net',
                   n_contour=20):
    """
    Mapa 2D de densidad de carga sobre un plano que pasa por la partícula.

    Parámetros
    ----------
    r_vecs    : (N, 3)  vectores posición relativa a la partícula
    net_vals  : (N,)    densidad de carga neta en cada nodo
    plane     : str     'xy', 'xz', 'yz' o 'nx,ny,nz' (normal al plano)
    dz        : float   semiancho de la rebanada en la dirección normal
    plot_type : str     'pcolor' (sin interpolar) o 'contour'
    quantity  : str     'net' o 'abs'
    n_contour : int     número de niveles para contourf
    """
    try:
        normal, plane_label = _parse_plane_normal(plane)
    except ValueError as e:
        print(f"ERROR: {e}")
        return

    e1, e2 = _build_plane_axes(normal)

    # Distancia de cada nodo al plano (el plano pasa por el origen = partícula)
    dist_n = r_vecs @ normal
    mask = np.abs(dist_n) <= dz

    if np.sum(mask) == 0:
        print(f"AVISO: No hay nodos en plano {plane_label} con dist <= {dz:.2f}")
        return

    # Coordenadas en el plano
    a_vals = r_vecs[mask] @ e1
    b_vals = r_vecs[mask] @ e2
    vals   = net_vals[mask]
    if quantity == 'abs':
        vals = np.abs(vals)

    n_nodes = int(np.sum(mask))
    print(f"  Nodos en plano {plane_label} (dist <= {dz:.2f} lu): {n_nodes}")

    r_norms = np.linalg.norm(r_vecs[mask], axis=1)
    rho_dh_expected = None
    if kappa > 0:
        rho_dh_expected = TheoreticalChargeDensity.compute_dh(
            np.array([np.mean(r_norms)]), q=charge, kappa=kappa)[0]
        if quantity == 'abs':
            rho_dh_expected = abs(rho_dh_expected)

    vmin = float(np.nanmin(vals))
    vmax = float(np.nanmax(vals))
    cmap = 'RdBu_r' if (vmin < 0 and vmax > 0) else 'viridis'
    label_q = '|ρ_net|' if quantity == 'abs' else 'ρ_net'

    fig, ax = plt.subplots(figsize=(8, 7))

    if plot_type == 'contour':
        from scipy.interpolate import griddata
        margin = 0.5
        a_grid = np.linspace(a_vals.min() - margin, a_vals.max() + margin, 300)
        b_grid = np.linspace(b_vals.min() - margin, b_vals.max() + margin, 300)
        AA, BB = np.meshgrid(a_grid, b_grid)
        ZZ = griddata((a_vals, b_vals), vals, (AA, BB), method='linear')
        cf = ax.contourf(AA, BB, ZZ, levels=n_contour,
                         cmap=cmap, vmin=vmin, vmax=vmax)
        ax.contour(AA, BB, ZZ, levels=n_contour,
                   colors='k', linewidths=0.3, alpha=0.4)
        cbar = plt.colorbar(cf, ax=ax, shrink=0.8, pad=0.02)
    else:
        sc = ax.scatter(a_vals, b_vals, c=vals, cmap=cmap,
                        vmin=vmin, vmax=vmax,
                        s=40, marker='s', linewidths=0)
        cbar = plt.colorbar(sc, ax=ax, shrink=0.8, pad=0.02)

    cbar.set_label(f'{label_q} (adim.)', fontsize=font_size - 2)
    ax.plot(0, 0, 'o', color='red', markersize=5, label='partícula')
    ax.grid(True, which='major', color='#d0d0d0', linewidth=0.4, alpha=0.6)
    ax.set_axisbelow(True)

    r_min_nodes = float(np.min(r_norms))
    r_max_nodes = float(np.max(r_norms))

    def _vec_str(v):
        return f'({v[0]:+.2f},{v[1]:+.2f},{v[2]:+.2f})'

    ax.set_title(
        f'Plano  n̂={plane_label}  {label_q}  dist_normal≤{dz:.2f} lu  n={n_nodes}\n'
        f'e₁={_vec_str(e1)}  e₂={_vec_str(e2)}  '
        f'r∈[{r_min_nodes:.2f},{r_max_nodes:.2f}] lu',
        fontsize=font_size - 1)
    ax.set_xlabel("e₁ (lu)", fontsize=font_size)
    ax.set_ylabel("e₂ (lu)", fontsize=font_size)
    ax.set_aspect('equal')
    ax.legend(fontsize=font_size - 4, loc='upper right')

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Mapa de plano guardado en: {output_file}")
    plt.close()


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Grafica densidades de carga (qsi) de Ewald vs. teoría Debye-Hückel',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Ejemplos:
    python plot_qsi_ewald.py -n 2000 -L 32 -q 1.0 -kt 0.00001 -eps 1e4 \\
        --rho-el 8e-3

    python plot_qsi_ewald.py -n 2000 -L 32 -q 1.0 -kt 0.00001 -eps 1e4 \\
        --rho-el 8e-3 --show-species0 --show-species1 --show-net

    python plot_qsi_ewald.py -n 2000 -L 32 -q 1.0 -kt 0.00001 -eps 1e4 \\
        --rho-el 8e-3 --direction x --angle-tol 5 \\
        --show-dh-periodic --dh-periodic-shells 1 \\
        --pointwise-error --log-log --csv qsi_data

    # Promedio radial con cascara de 0.5 lu y exportar CSV radial:
    python plot_qsi_ewald.py -n 2000 -L 32 -q 1.0 -kt 0.00001 -eps 1e4 \\
        --rho-el 8e-3 --show-net --radial --shell-width 0.5 --csv-radial qsi_radial
""")

    parser.add_argument('-n', '--nstep', type=int, required=True,
                        help='Número de paso de tiempo')
    parser.add_argument('-L', '--grid-size', type=int, required=True,
                        help='Tamaño de malla cúbica L×L×L')
    parser.add_argument('-q', '--charge', type=float, default=1.0,
                        help='Carga de la partícula')
    parser.add_argument('-kt', '--kt', type=float, default=1.0,
                        help='Energía térmica kT')
    parser.add_argument('-eps', '--epsilon', type=float, default=1.0,
                        help='Permitividad relativa')
    parser.add_argument('--kappa', type=float, default=None,
                        help='Parámetro de Debye kappa (si no se da, se calcula desde --rho-el)')
    parser.add_argument('--rho-el', type=float, default=None,
                        help='Densidad de electrolito para calcular kappa')
    parser.add_argument('-d', '--directory', default='.',
                        help='Directorio con los archivos de Ludwig (default: .)')

    parser.add_argument('--rmin', type=float, default=0.5,
                        help='Radio mínimo a graficar (default: 0.5)')
    parser.add_argument('--rmax', type=float, default=None,
                        help='Radio máximo (default: L/2)')
    parser.add_argument('--particle-radius', type=float, default=0.0,
                        help='Radio de la partícula (lu). Si > 0, agrega curva DH '
                             'con radio finito: ρ(r)=-κ²q exp(-κ(r-a))/(4πr(1+κa)) '
                             '(default: 0 = DH puntual)')

    parser.add_argument('--num-species', type=int, default=2,
                        help='Número de especies iónicas en el archivo qsi (default: 2)')

    # Qué mostrar
    parser.add_argument('--show-species0', action='store_true', default=False,
                        help='Mostrar densidad de especie 0 (cationes)')
    parser.add_argument('--show-species1', action='store_true', default=False,
                        help='Mostrar densidad de especie 1 (aniones)')
    parser.add_argument('--show-net', action='store_true', default=None,
                        help='Mostrar carga neta ρ_net = ρ₊ - ρ₋')
    parser.add_argument('--no-net', action='store_true', default=False,
                        help='Ocultar carga neta')

    parser.add_argument('--show-dh', action='store_true', default=None,
                        help='Mostrar curva ρ_net DH infinita')
    parser.add_argument('--show-dh-periodic', action='store_true', default=None,
                        help='Mostrar curva ρ_net DH periódica (evaluada a lo largo '
                             'de la dirección seleccionada)')
    parser.add_argument('--show-dh-periodic-pts', action='store_true', default=False,
                        help='Mostrar ρ_net DH periódico evaluado en la posición 3D '
                             'exacta de cada nodo (más preciso; usa --dh-periodic-shells)')
    parser.add_argument('--dh-periodic-shells', type=int, default=1,
                        help='Capas de imágenes periódicas (default: 1)')
    parser.add_argument('--pointwise-error', action='store_true', default=False,
                        help='Calcular error en la posición exacta de cada nodo')
    parser.add_argument('--no-error-inf', action='store_true', default=False,
                        help='Ocultar error vs DH infinito')
    parser.add_argument('--no-error-dir', action='store_true', default=False,
                        help='Ocultar error vs DH periódico direccional')
    # CHANGE INIT - OnlyErrorDHPerRad - Opción para mostrar solo error vs DH periódico radio finito
    parser.add_argument('--only-error-dh-per-rad', action='store_true', default=False,
                        help='Mostrar en el panel de error SOLO la curva de error relativo '
                             'vs DH periódico con partícula de radio a. Requiere '
                             '--particle-radius > 0 y --show-dh-periodic. '
                             'Suprime las demás curvas de error.')
    # CHANGE END - OnlyErrorDHPerRad

    # CHANGE INIT - GaussianDH - Argumentos para fuente de carga gaussiana
    parser.add_argument('--gaussian-source', action='store_true', default=False,
                        help='Usar distribución gaussiana ρ(r)=Q/(π^(3/2)σ³)exp(-r²/σ²) '
                             'como fuente en la teoría DH. Activa curva ρ_ions gaussiana-DH.')
    parser.add_argument('--sigma', type=float, default=1.0,
                        help='Anchura σ de la distribución gaussiana de carga (lu). '
                             'Solo se usa con --gaussian-source (default: 1.0)')
    parser.add_argument('--gaussian-error', action='store_true', default=False,
                        help='Mostrar error relativo vs teoría gaussiana-DH infinita '
                             'en el panel 2. Requiere --gaussian-source.')
    parser.add_argument('--gaussian-error-periodic', action='store_true', default=False,
                        help='Mostrar error relativo vs teoría gaussiana-DH periódica '
                             'en el panel 2. Requiere --gaussian-source y --show-dh-periodic.')
    # CHANGE END - GaussianDH

    parser.add_argument('--direction', default=None,
                        help="Filtrar por dirección: 'x','y','z','xy','xyz' o 'dx,dy,dz'")
    parser.add_argument('--angle-tol', type=float, default=15.0,
                        help='Tolerancia angular en grados (default: 15)')

    parser.add_argument('--semilog', action='store_true', default=False,
                        help='Escala semilog en Y (default: lineal; usa |valor| para log)')
    parser.add_argument('--log-log', action='store_true', default=False,
                        help='Escala log-log en ambos ejes (usa |valor| automáticamente)')
    parser.add_argument('--log-log-error', action='store_true', default=False,
                        help='Escala log en eje Y del panel de error')
    parser.add_argument('--show-error', action='store_true', default=False,
                        help='Mostrar panel de error relativo punto a punto vs DH inf '
                             '(útil cuando se grafican todos los nodos sin --direction)')
    parser.add_argument('--error-ymax', type=float, default=None,
                        help='Límite superior del eje Y del panel de error en %% '
                             '(útil para recortar outliers, ej: --error-ymax 200)')
    parser.add_argument('--xmin', type=float, default=None)
    parser.add_argument('--ymin', type=float, default=None)

    parser.add_argument('--output', default='qsi_ewald',
                        help='Nombre base del archivo de salida (default: qsi_ewald)')
    parser.add_argument('--max-points', type=int, default=10000,
                        help='Máximo número de puntos a graficar (default: 10000)')
    parser.add_argument('--csv', default=None,
                        help='Exportar datos a CSV (nombre base; se añade -n_NNNNNN.csv)')

    # Promedio radial
    parser.add_argument('--radial', action='store_true', default=False,
                        help='Graficar promedio radial (media ± std) en cascaras esféricas')
    parser.add_argument('--radial-scatter', action='store_true', default=False,
                        help='Graficar todos los puntos individuales coloreados por cascara '
                             '(en lugar de la media). Requiere --radial o se activa solo.')
    parser.add_argument('--shell-width', type=float, default=1.0,
                        help='Grosor de cada cascara esférica en lu (default: 1.0)')
    parser.add_argument('--csv-radial', default=None,
                        help='Exportar promedio radial a CSV (nombre base)')
    parser.add_argument('--csv-scatter', default=None,
                        help='Exportar nodos individuales scatter a CSV (nombre base)')

    # Mapa angular
    parser.add_argument('--angular-map', action='store_true', default=False,
                        help='Generar mapa de calor (θ,φ) en una cáscara esférica')
    parser.add_argument('--angular-r', type=float, default=None,
                        help='Radio central de la cáscara para el mapa angular (lu)')
    parser.add_argument('--angular-dr', type=float, default=1.0,
                        help='Grosor de la cáscara para el mapa angular (default: 1.0 lu)')
    parser.add_argument('--angular-n-theta', type=int, default=90,
                        help='Resolución en θ del mapa angular (default: 90)')
    parser.add_argument('--angular-n-phi', type=int, default=180,
                        help='Resolución en φ del mapa angular (default: 180)')
    parser.add_argument('--angular-abs', action='store_true', default=False,
                        help='Usar |ρ_net| en lugar de ρ_net en el mapa angular')

    # Mapa de plano
    parser.add_argument('--plane-map', action='store_true', default=False,
                        help='Generar mapa 2D de densidad de carga sobre un plano')
    parser.add_argument('--plane', default='xy',
                        help="Normal al plano (el plano pasa por la partícula): "
                             "'xy','xz','yz' o 'nx,ny,nz' (default: xy)")
    parser.add_argument('--plane-dz', type=float, default=0.5,
                        help='Semiancho de la rebanada en la dirección normal (default: 0.5 lu)')
    parser.add_argument('--plane-type', default='pcolor',
                        choices=['pcolor', 'contour'],
                        help="Tipo de gráfico: 'pcolor' (sin interpolar) o 'contour' (default: pcolor)")
    parser.add_argument('--plane-contour-levels', type=int, default=20,
                        help='Número de niveles para contourf (default: 20)')
    parser.add_argument('--plane-abs', action='store_true', default=False,
                        help='Usar |ρ_net| en lugar de ρ_net en el mapa de plano')

    args = parser.parse_args()

    # Parsear dirección
    try:
        direction = parse_direction(args.direction)
    except ValueError as e:
        print(f"ERROR: {e}")
        sys.exit(1)
    if direction is not None:
        print(f"Dirección seleccionada: {_direction_label(direction)}, "
              f"tolerancia ±{args.angle_tol:.1f}°")

    # Defaults de visualización
    show_net = not args.no_net
    if args.show_net:
        show_net = True

    # Si no se especifica nada para mostrar, mostrar carga neta y DH inf
    nothing_requested = (not args.show_species0 and not args.show_species1
                         and args.show_net is None and not args.no_net)
    if nothing_requested:
        show_net = True

    if args.show_dh is None and args.show_dh_periodic is None:
        args.show_dh = True
        args.show_dh_periodic = False
    if args.show_dh is None:
        args.show_dh = False
    if args.show_dh_periodic is None:
        args.show_dh_periodic = False

    # Calcular kappa
    if args.kappa is not None:
        kappa = args.kappa
    elif args.rho_el is not None:
        lb = 1.0 / (4.0 * np.pi * args.epsilon * args.kt)
        V = float(args.grid_size ** 3)
        rho_counterions = abs(args.charge) / V
        rho_total = args.rho_el + rho_counterions
        kappa = np.sqrt(4.0 * np.pi * lb * rho_total * 2.0)
        print(f"\nParámetros de Debye-Hückel:")
        print(f"  Longitud de Bjerrum λ_B = {lb:.6e}")
        print(f"  ρ_el = {args.rho_el:.6e}, ρ_contraiones = {rho_counterions:.6e}")
        print(f"  kappa = {kappa:.6e}  (λ_D = {1.0/kappa:.4f} lu)")
    else:
        kappa = 0.0
        print("AVISO: kappa=0 (sin electrolito). La curva DH no estará disponible.")

    grid_size = (args.grid_size, args.grid_size, args.grid_size)
    base_dir = args.directory

    nstep_str9 = f"{args.nstep:09d}"
    nstep_str8 = f"{args.nstep:08d}"
    qsi_file     = os.path.join(base_dir, f"qsi-{nstep_str9}.001-001")
    colloid_file = os.path.join(base_dir, f"config.cds{nstep_str8}.001-001")

    for fname in [qsi_file, colloid_file]:
        if not os.path.exists(fname):
            print(f"ERROR: No se encontró el archivo: {fname}")
            sys.exit(1)

    print(f"Leyendo posición de la partícula desde: {colloid_file}")
    particle_pos = ColloidReader.read_colloid_position_cds(colloid_file)
    print(f"  Posición: {particle_pos}")

    print(f"\nLeyendo densidades de carga desde: {qsi_file}")
    reader = QsiReader(qsi_file, grid_size, num_species=args.num_species)
    try:
        species_list = reader.read()
    except Exception as e:
        print(f"ERROR: {e}")
        sys.exit(1)

    print(f"  Número de especies: {len(species_list)}")
    for s, sp in enumerate(species_list):
        sign = '+' if s == 0 else '-'
        print(f"  Especie {s} (ρ{sign}): min={sp.min():.4e}, max={sp.max():.4e}, "
              f"sum={sp.sum():.4e}")
    if len(species_list) >= 2:
        net_3d = reader.net_charge_3d()
        print(f"  ρ_net total = {net_3d.sum():.4e}")

    print(f"\nExtrayendo puntos (rmin={args.rmin}, rmax={args.rmax})...")
    distances, species_vals, net_vals, r_vecs = extract_points_qsi(
        species_list, particle_pos, grid_size,
        min_radius=args.rmin, max_radius=args.rmax,
        direction=direction, angle_tol_deg=args.angle_tol)
    print(f"  Puntos extraídos: {len(distances)}")

    output_base, output_ext = os.path.splitext(args.output)
    if not output_ext:
        output_ext = '.png'
    output_file = f"{output_base}-n_{args.nstep:09d}{output_ext}"

    log_scale = args.semilog and not args.log_log
    print(f"\nGenerando gráfico...")
    plot_qsi(distances, species_vals, net_vals, r_vecs,
             grid_size, args.charge, args.epsilon, args.kt, kappa,
             rmin=args.rmin, rmax=args.rmax, output_file=output_file, max_points=args.max_points,
             log_scale=log_scale, log_log=args.log_log,
             log_log_error=args.log_log_error,
             xmin=args.xmin, ymin=args.ymin,
             show_species0=args.show_species0,
             show_species1=args.show_species1,
             show_net=show_net,
             show_dh=args.show_dh,
             show_dh_periodic=args.show_dh_periodic,
             show_dh_periodic_pts=args.show_dh_periodic_pts,
             dh_periodic_shells=args.dh_periodic_shells,
             direction=direction, angle_tol_deg=args.angle_tol,
             pointwise_error=args.pointwise_error,
             show_error_inf=not args.no_error_inf,
             show_error_dir=not args.no_error_dir,
             show_error=args.show_error,
             error_ymax=args.error_ymax,
             particle_radius=args.particle_radius,
             only_error_dh_per_rad=args.only_error_dh_per_rad,
             gaussian_source=args.gaussian_source,
             sigma=args.sigma,
             gaussian_error=args.gaussian_error,
             gaussian_error_periodic=args.gaussian_error_periodic)

    if args.radial or args.radial_scatter:
        do_scatter = args.radial_scatter
        rmax_radial = args.rmax if args.rmax is not None else min(grid_size) / 2.0
        print(f"\nCalculando promedio radial (Δr={args.shell_width}, "
              f"rmin={args.rmin}, rmax={rmax_radial:.1f})...")
        r_mids, mean_species_r, mean_net_r, std_net_r, n_nodes_r = compute_radial_average(
            distances, species_vals, net_vals,
            shell_width=args.shell_width, rmin=args.rmin, rmax=rmax_radial)
        print(f"  Cascaras con datos: {np.sum(~np.isnan(mean_net_r))}")

        radial_base, radial_ext = os.path.splitext(args.output)
        if not radial_ext:
            radial_ext = '.png'
        suffix = '-scatter' if do_scatter else '-radial'
        radial_file = f"{radial_base}{suffix}-n_{args.nstep:09d}{radial_ext}"
        print(f"Generando gráfico {'scatter' if do_scatter else 'radial'}...")
        plot_radial_average(
            r_mids, mean_species_r, mean_net_r, std_net_r,
            grid_size, args.charge, args.epsilon, args.kt, kappa,
            shell_width=args.shell_width, rmin=args.rmin,
            distances=distances,
            species_vals=species_vals,
            net_vals=net_vals,
            r_vecs=r_vecs,
            scatter=do_scatter,
            output_file=radial_file,
            log_scale=log_scale, log_log=args.log_log,
            xmin=args.xmin, ymin=args.ymin,
            show_species0=args.show_species0,
            show_species1=args.show_species1,
            show_net=show_net,
            show_dh=args.show_dh,
            show_dh_periodic=args.show_dh_periodic,
            dh_periodic_shells=args.dh_periodic_shells,
            show_dh_periodic_pts=args.show_dh_periodic_pts,
            show_error_inf=not args.no_error_inf,
            error_ymax=args.error_ymax,
            particle_radius=args.particle_radius,
            gaussian_source=args.gaussian_source,
            sigma=args.sigma,
            gaussian_error=args.gaussian_error,
            gaussian_error_periodic=args.gaussian_error_periodic)

        if args.csv_radial is not None:
            csv_r_base, csv_r_ext = os.path.splitext(args.csv_radial)
            if not csv_r_ext:
                csv_r_ext = '.csv'
            csv_r_file = f"{csv_r_base}-n_{args.nstep:09d}{csv_r_ext}"
            print(f"\nExportando promedio radial a CSV: {csv_r_file}")
            export_csv_radial(r_mids, mean_species_r, mean_net_r, std_net_r, n_nodes_r,
                              charge=args.charge, kappa=kappa,
                              csv_file=csv_r_file)

        if do_scatter and args.csv_scatter is not None:
            csv_s_base, csv_s_ext = os.path.splitext(args.csv_scatter)
            if not csv_s_ext:
                csv_s_ext = '.csv'
            csv_s_file = f"{csv_s_base}-n_{args.nstep:09d}{csv_s_ext}"
            print(f"\nExportando scatter a CSV: {csv_s_file}")
            export_csv_scatter(distances, species_vals, net_vals, r_vecs,
                               particle_pos=particle_pos,
                               charge=args.charge, epsilon=args.epsilon,
                               kappa=kappa, csv_file=csv_s_file)

    if args.angular_map:
        angular_r = args.angular_r
        if angular_r is None:
            angular_r = min(grid_size) / 4.0
            print(f"  --angular-r no especificado; usando r={angular_r:.2f} (L/4)")
        # Extraer todos los nodos sin filtro de dirección para el mapa angular
        print(f"\nGenerando mapa angular (r={angular_r:.2f} ± {args.angular_dr/2:.2f} lu)...")
        _, _, net_vals_full, r_vecs_full = extract_points_qsi(
            species_list, particle_pos, grid_size,
            min_radius=0.5, max_radius=None,
            direction=None, angle_tol_deg=180.0)
        quantity = 'abs' if args.angular_abs else 'net'
        ang_base, ang_ext = os.path.splitext(args.output)
        if not ang_ext:
            ang_ext = '.png'
        ang_file = f"{ang_base}-angular-r{angular_r:.1f}-n_{args.nstep:09d}{ang_ext}"
        plot_angular_map(
            r_vecs_full, net_vals_full,
            r_shell=angular_r, dr=args.angular_dr,
            charge=args.charge, kappa=kappa,
            output_file=ang_file,
            n_theta=args.angular_n_theta,
            n_phi=args.angular_n_phi,
            quantity=quantity)

    if args.plane_map:
        print(f"\nGenerando mapa de plano {args.plane.upper()} "
              f"(|normal| <= {args.plane_dz:.2f} lu, tipo={args.plane_type})...")
        # Reusar los nodos completos si ya fueron extraídos, si no extraer
        try:
            r_vecs_full
        except NameError:
            _, _, net_vals_full, r_vecs_full = extract_points_qsi(
                species_list, particle_pos, grid_size,
                min_radius=0.0, max_radius=None,
                direction=None, angle_tol_deg=180.0)
        quantity_pl = 'abs' if args.plane_abs else 'net'
        pl_base, pl_ext = os.path.splitext(args.output)
        if not pl_ext:
            pl_ext = '.png'
        pl_file = (f"{pl_base}-plane_{args.plane}-dz{args.plane_dz:.2f}"
                   f"-{args.plane_type}-n_{args.nstep:09d}{pl_ext}")
        plot_plane_map(
            r_vecs_full, net_vals_full,
            plane=args.plane,
            dz=args.plane_dz,
            charge=args.charge, kappa=kappa,
            output_file=pl_file,
            plot_type=args.plane_type,
            quantity=quantity_pl,
            n_contour=args.plane_contour_levels)

    if args.csv is not None:
        csv_base, csv_ext = os.path.splitext(args.csv)
        if not csv_ext:
            csv_ext = '.csv'
        csv_file = f"{csv_base}-n_{args.nstep:09d}{csv_ext}"
        print(f"\nExportando datos a CSV: {csv_file}")
        export_csv_qsi(species_list, particle_pos, grid_size,
                       charge=args.charge, epsilon=args.epsilon,
                       kt=args.kt, kappa=kappa,
                       L=float(args.grid_size),
                       n_shells=args.dh_periodic_shells,
                       rmin=args.rmin, csv_file=csv_file)

    print("\nAnálisis completado!")


if __name__ == '__main__':
    main()
