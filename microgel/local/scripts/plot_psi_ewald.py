#!/usr/bin/env python3
"""
Script para graficar el potencial electrostático ψ calculado por Ludwig (Ewald),
comparándolo con la solución teórica de Debye-Hückel.
Soporta las componentes psi_real y psi_fourier de la misma manera que
plot_componentes_ewald.py maneja efield_real y efield_fourier.

Datos de entrada: archivos psi-NNNNNNNN.001-001 generados por Ludwig.
Opcionalmente: psi_real-NNNNNNNN.001-001 y psi_fourier-NNNNNNNN.001-001.

Uso básico:
    python plot_psi_ewald.py -n 2000 -L 32 -q 1.0 -kt 0.00001 -eps 1.0e4 \\
        --rho-el 8.0e-3 --alpha 0.7 --rc 3.0

Con componentes:
    python plot_psi_ewald.py -n 2000 -L 32 -q 1.0 -kt 0.00001 -eps 1.0e4 \\
        --rho-el 8.0e-3 --alpha 0.7 --rc 3.0 --show-components

Con dirección y error puntual:
    python plot_psi_ewald.py -n 2000 -L 32 -q 1.0 -kt 0.00001 -eps 1.0e4 \\
        --rho-el 8.0e-3 --alpha 0.7 --rc 3.0 \\
        --direction 1,1,1 --angle-tol 5 --pointwise-error \\
        --show-dh-periodic --dh-periodic-shells 1

Con fuente gaussiana (σ=1.0):
    python plot_psi_ewald.py -n 2000 -L 32 -q 1.0 -kt 0.00001 -eps 1.0e4 \\
        --rho-el 8.0e-3 --alpha 0.7 --rc 3.0 \\
        --gaussian-source --sigma 1.0

Con fuente gaussiana y error vs gaussiana-DH periódica:
    python plot_psi_ewald.py -n 2000 -L 32 -q 1.0 -kt 0.00001 -eps 1.0e4 \\
        --rho-el 8.0e-3 --alpha 0.7 --rc 3.0 \\
        --gaussian-source --sigma 1.0 --show-dh-periodic --dh-periodic-shells 1 \\
        --direction 1,0,0 --angle-tol 5 --gaussian-error --gaussian-error-periodic
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
import os

font_size = 16


# ---------------------------------------------------------------------------
# Lector de archivos psi (escalar)
# ---------------------------------------------------------------------------

class PsiReader:
    """Lee y procesa archivos psi de Ludwig (un escalar por nodo)"""

    def __init__(self, filename, grid_size):
        self.filename = filename
        self.nx, self.ny, self.nz = grid_size
        self.psi = None

    def read(self):
        """Lee el archivo psi y lo organiza en malla 3D"""
        data = np.loadtxt(self.filename)
        n_points = self.nx * self.ny * self.nz
        if data.ndim != 1 or len(data) != n_points:
            raise ValueError(
                f"El archivo tiene {data.size} valores, "
                f"pero se esperaban {n_points} (1 escalar por nodo)")
        # Ludwig escribe en orden x-major (para cada x, para cada y, todos los z)
        self.psi = data.reshape((self.nx, self.ny, self.nz))
        return self.psi


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
            raise ValueError(f"Línea 36 de {filename} no tiene 3 valores: '{pos_line}'")
        return np.array([float(parts[0]), float(parts[1]), float(parts[2])])


# ---------------------------------------------------------------------------
# Curvas teóricas de Debye-Hückel (potencial)
# ---------------------------------------------------------------------------

class TheoreticalPotential:
    """Potencial eléctrico teórico de Debye-Hückel: psi(r) = (q/4*pi*eps*kT)*exp(-kappa*r)/r"""

    @staticmethod
    def compute_dh(distances, q=1.0, epsilon=1.0, kt=1.0, kappa=0.0):
        """psi_DH(r) para sistema infinito"""
        prefactor = q / (4.0 * np.pi * epsilon * kt)
        r = np.maximum(distances, 0.01)
        if kappa > 0:
            return prefactor * np.exp(-kappa * r) / r
        else:
            return prefactor / r

    @staticmethod
    def compute_dh_periodic(distances, q=1.0, epsilon=1.0, kt=1.0, kappa=0.0,
                            L=32.0, n_shells=1, direction=None):
        """
        psi periodico evaluado a lo largo de una direccion d.
        Suma imagen central + imagenes periodicas (escalar).
        """
        if direction is None:
            direction = np.array([1.0, 0.0, 0.0])
        d = np.asarray(direction, dtype=float)
        d = d / np.linalg.norm(d)

        prefactor = q / (4.0 * np.pi * epsilon * kt)
        r = np.maximum(distances, 0.01)

        def psi_scalar(dist):
            dist = np.maximum(dist, 1e-10)
            if kappa > 0:
                return prefactor * np.exp(-kappa * dist) / dist
            else:
                return prefactor / dist

        # Punto de evaluacion: r_vec = r * d
        r_vec = r[:, np.newaxis] * d[np.newaxis, :]   # (N, 3)

        psi_total = psi_scalar(r)  # contribucion central

        for inx in range(-n_shells, n_shells + 1):
            for iny in range(-n_shells, n_shells + 1):
                for inz in range(-n_shells, n_shells + 1):
                    if inx == 0 and iny == 0 and inz == 0:
                        continue
                    r_img = np.array([inx * L, iny * L, inz * L])
                    dr = r_vec - r_img[np.newaxis, :]
                    dist_img = np.linalg.norm(dr, axis=1)
                    psi_total = psi_total + psi_scalar(dist_img)

        return psi_total

    @staticmethod
    def compute_dh_periodic_at_points(r_vecs, q=1.0, epsilon=1.0, kt=1.0, kappa=0.0,
                                      L=32.0, n_shells=1):
        """psi periodico en puntos 3D arbitrarios r_vecs (N,3) relativos a la particula"""
        prefactor = q / (4.0 * np.pi * epsilon * kt)
        r_vecs = np.asarray(r_vecs, dtype=float)

        def psi_scalar(dist):
            dist = np.maximum(dist, 1e-10)
            if kappa > 0:
                return prefactor * np.exp(-kappa * dist) / dist
            else:
                return prefactor / dist

        dist0 = np.linalg.norm(r_vecs, axis=1)
        psi_total = psi_scalar(dist0)

        for inx in range(-n_shells, n_shells + 1):
            for iny in range(-n_shells, n_shells + 1):
                for inz in range(-n_shells, n_shells + 1):
                    if inx == 0 and iny == 0 and inz == 0:
                        continue
                    r_img = np.array([inx * L, iny * L, inz * L])
                    dr = r_vecs - r_img[np.newaxis, :]
                    dist_img = np.linalg.norm(dr, axis=1)
                    psi_total = psi_total + psi_scalar(dist_img)

        return psi_total


    # CHANGE INIT - DHFiniteRadius - Potencial DH con partícula de radio finito a
    @staticmethod
    def compute_dh_finite_radius(distances, q=1.0, epsilon=1.0, kt=1.0, kappa=0.0, a=0.0):
        """
        Potencial DH para partícula esférica de radio 'a' (sistema infinito).

        Solución exacta linealizada de Poisson-Boltzmann para r > a:
            psi(r) = q / (4π ε kT) * exp(-κ(r-a)) / (r (1+κa))

        Para r ≤ a retorna NaN.
        """
        prefactor = q / (4.0 * np.pi * epsilon * kt)
        r = np.asarray(distances, dtype=float)
        result = np.full_like(r, np.nan)
        outside = r > a
        r_out = r[outside]
        if kappa > 0:
            result[outside] = prefactor * np.exp(-kappa * (r_out - a)) / (r_out * (1.0 + kappa * a))
        else:
            result[outside] = prefactor / r_out
        return result

    @staticmethod
    def compute_dh_periodic_finite_radius(distances, q=1.0, epsilon=1.0, kt=1.0, kappa=0.0,
                                          a=0.0, L=32.0, n_shells=1, direction=None):
        """
        psi DH periódico con partícula de radio finito 'a' a lo largo de una dirección.

        Imagen central: psi_DH(r; a)  para r > a (NaN para r <= a).
        Imágenes periódicas: tratadas como cargas puntuales.
        """
        if direction is None:
            direction = np.array([1.0, 0.0, 0.0])
        d = np.asarray(direction, dtype=float)
        d = d / np.linalg.norm(d)

        prefactor = q / (4.0 * np.pi * epsilon * kt)
        r = np.asarray(distances, dtype=float)
        r_vec = r[:, np.newaxis] * d[np.newaxis, :]  # (N, 3)

        def psi_point(dist):
            dist = np.maximum(dist, 1e-10)
            if kappa > 0:
                return prefactor * np.exp(-kappa * dist) / dist
            else:
                return prefactor / dist

        # Imagen central con radio finito
        result = np.full_like(r, np.nan)
        outside = r > a
        r_out = r[outside]
        if kappa > 0:
            result[outside] = prefactor * np.exp(-kappa * (r_out - a)) / (r_out * (1.0 + kappa * a))
        else:
            result[outside] = prefactor / r_out

        # Imágenes periódicas (puntuales)
        for inx in range(-n_shells, n_shells + 1):
            for iny in range(-n_shells, n_shells + 1):
                for inz in range(-n_shells, n_shells + 1):
                    if inx == 0 and iny == 0 and inz == 0:
                        continue
                    r_img = np.array([inx * L, iny * L, inz * L])
                    dr = r_vec - r_img[np.newaxis, :]
                    dist_img = np.linalg.norm(dr, axis=1)
                    result[outside] += psi_point(dist_img)[outside]

        return result
    # CHANGE END - DHFiniteRadius

    # CHANGE INIT - GaussianDH - Potencial DH con fuente de carga gaussiana (sistema infinito y periódico)
    @staticmethod
    def compute_gaussian_dh(distances, q=1.0, epsilon=1.0, kt=1.0, kappa=0.0, sigma=1.0):
        """
        Potencial electrostático de Poisson-Boltzmann linealizado (Debye-Hückel)
        para una fuente de carga con distribución gaussiana:

            ρ(r) = Q / (π^(3/2) σ³) · exp(-r²/σ²)

        Calculado vía la función de Green de Yukawa en simetría esférica:

            ψ(r) = (1/κr) · { e^{-κr} · ∫₀ʳ sinh(κr') r' ρ(r') dr'
                             + sinh(κr) · ∫ᵣ^∞ e^{-κr'} r' ρ(r') dr' }

        donde el prefactor 1/(ε kT) está incluido al normalizar ρ.

        Este método es numéricamente estable para todos los r, σ y κ.

        Casos límite:
          κ→0: ψ(r) → Q/(4πε kT) · erf(r/σ)/r   (Coulomb suavizado por gaussiana)
          σ→0: ψ(r) → Q/(4πε kT) · e^{-κr}/r     (Yukawa puntual estándar)
          r→∞: ψ(r) → exp(κ²σ²/4) · Q/(4πε kT) · e^{-κr}/r  (tail Yukawa con factor gaussiano)

        Args:
            distances: array de distancias radiales
            q:        carga total Q
            epsilon:  permitividad relativa
            kt:       energía térmica kT
            kappa:    parámetro de Debye κ = 1/λ_D (0 = Coulomb puro)
            sigma:    anchura de la distribución gaussiana σ

        Returns:
            psi: array de potencial eléctrico
        """
        from scipy.special import erf
        prefactor = q / (4.0 * np.pi * epsilon * kt)
        r = np.asarray(distances, dtype=float)
        r_safe = np.maximum(r, 1e-10)

        if kappa <= 0:
            # κ=0: ψ(r) = Q/(4πε kT) · erf(r/σ)/r
            return prefactor * erf(r_safe / sigma) / r_safe

        # Pre-computar integranda sobre malla fina de r'
        rmax = max(15.0 * sigma, float(np.max(r_safe)) * 1.5)
        n_rho = 20000
        rp = np.linspace(1e-8, rmax, n_rho)
        drp = rp[1] - rp[0]
        # rho normalizada (sin prefactor q/(4pi*eps*kT) - se aplica al final)
        rho_rp = 1.0 / (np.pi**1.5 * sigma**3) * np.exp(-rp**2 / sigma**2) * rp

        integrand_inner = np.sinh(kappa * rp) * rho_rp
        integrand_outer = np.exp(-kappa * rp) * rho_rp

        cumul_inner = np.cumsum(integrand_inner) * drp
        total_outer = np.sum(integrand_outer) * drp
        cumul_outer = total_outer - np.cumsum(integrand_outer) * drp

        inner_at_r = np.interp(r_safe, rp, cumul_inner)
        outer_at_r = np.interp(r_safe, rp, cumul_outer)

        # Factor 4π de la reducción de la integral 3D a 1D en simetría esférica:
        # ∫ρ(r')G(|r-r'|)d³r' = (4π/κr)[e^{-κr}∫₀ʳsinh(κr')r'ρ̃dr' + sinh(κr)∫ᵣ^∞e^{-κr'}r'ρ̃dr']
        result = (4.0 * np.pi * prefactor / (kappa * r_safe)) * (
            np.exp(-kappa * r_safe) * inner_at_r
            + np.sinh(kappa * r_safe) * outer_at_r
        )
        return result

    @staticmethod
    def compute_gaussian_dh_periodic(distances, q=1.0, epsilon=1.0, kt=1.0, kappa=0.0,
                                     sigma=1.0, L=32.0, n_shells=1, direction=None):
        """
        ψ gaussiana-DH periódico a lo largo de una dirección dada.

        Imagen central: usa compute_gaussian_dh (fuente suavizada).
        Imágenes periódicas: tratadas como cargas puntuales de Yukawa (σ→0,
        ya que están muy lejos y la gaussiana no tiene efecto).

        Args:
            distances: array de distancias desde la imagen central
            direction: vector unitario de evaluación (None → eje x)
            n_shells:  capas de imágenes periódicas
        """
        from scipy.special import erfc, erf
        if direction is None:
            direction = np.array([1.0, 0.0, 0.0])
        d = np.asarray(direction, dtype=float)
        d = d / np.linalg.norm(d)

        prefactor = q / (4.0 * np.pi * epsilon * kt)
        r = np.asarray(distances, dtype=float)
        r_safe = np.maximum(r, 1e-10)
        r_vec = r_safe[:, np.newaxis] * d[np.newaxis, :]  # (N, 3)

        # Imagen central con fuente gaussiana
        result = TheoreticalPotential.compute_gaussian_dh(
            r_safe, q=q, epsilon=epsilon, kt=kt, kappa=kappa, sigma=sigma)

        # Imágenes periódicas como Yukawa puntual
        def psi_point(dist):
            dist = np.maximum(dist, 1e-10)
            if kappa > 0:
                return prefactor * np.exp(-kappa * dist) / dist
            else:
                return prefactor / dist

        for inx in range(-n_shells, n_shells + 1):
            for iny in range(-n_shells, n_shells + 1):
                for inz in range(-n_shells, n_shells + 1):
                    if inx == 0 and iny == 0 and inz == 0:
                        continue
                    r_img = np.array([inx * L, iny * L, inz * L])
                    dr = r_vec - r_img[np.newaxis, :]
                    dist_img = np.linalg.norm(dr, axis=1)
                    result = result + psi_point(dist_img)

        return result

    @staticmethod
    def compute_gaussian_dh_periodic_at_points(r_vecs, q=1.0, epsilon=1.0, kt=1.0,
                                               kappa=0.0, sigma=1.0, L=32.0, n_shells=1):
        """
        ψ gaussiana-DH periódico en puntos 3D arbitrarios r_vecs (N,3) relativos a la partícula.
        Imagen central con fuente gaussiana; imágenes periódicas como Yukawa puntual.
        """
        prefactor = q / (4.0 * np.pi * epsilon * kt)
        r_vecs = np.asarray(r_vecs, dtype=float)

        dist0 = np.linalg.norm(r_vecs, axis=1)
        result = TheoreticalPotential.compute_gaussian_dh(
            dist0, q=q, epsilon=epsilon, kt=kt, kappa=kappa, sigma=sigma)

        def psi_point(dist):
            dist = np.maximum(dist, 1e-10)
            if kappa > 0:
                return prefactor * np.exp(-kappa * dist) / dist
            else:
                return prefactor / dist

        for inx in range(-n_shells, n_shells + 1):
            for iny in range(-n_shells, n_shells + 1):
                for inz in range(-n_shells, n_shells + 1):
                    if inx == 0 and iny == 0 and inz == 0:
                        continue
                    r_img = np.array([inx * L, iny * L, inz * L])
                    dr = r_vecs - r_img[np.newaxis, :]
                    dist_img = np.linalg.norm(dr, axis=1)
                    result = result + psi_point(dist_img)

        return result
    # CHANGE END - GaussianDH

    @staticmethod
    def compute_ewald_components(distances, q=1.0, epsilon=1.0, kt=1.0, alpha=0.5, rc=5.0):
        """
        Componentes teóricas de Ewald para el potencial:
          psi_real(r)  = prefactor * erfc(alpha*r) / r   (espacio real, corte en rc)
          psi_recip(r) = prefactor * erf(alpha*r)  / r   (complemento, largo alcance)
          psi_real + psi_recip = prefactor / r  (Coulomb)
        prefactor = q / (4*pi*eps*kt)
        """
        from scipy.special import erfc, erf
        prefactor = q / (4.0 * np.pi * epsilon * kt)
        r = np.maximum(distances, 0.01)
        psi_real  = prefactor * erfc(alpha * r) / r
        psi_recip = prefactor * erf(alpha * r)  / r
        psi_real  = np.where(r > rc, 0.0, psi_real)
        return psi_real, psi_recip


# ---------------------------------------------------------------------------
# Filtro de direccion
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
            raise ValueError(f"--direction '{direction_str}' no reconocido. "
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
# Extraccion de puntos
# ---------------------------------------------------------------------------

def compute_background(psi_3d, particle_pos, grid_size, bg_rmin=12.0):
    """
    Calcula el potencial de fondo promediando nodos a distancia >= bg_rmin.
    El fondo es la contribucion negativa uniforme de los contraiones.
    """
    nx, ny, nz = grid_size
    px, py, pz = particle_pos
    vals = []
    for i in range(nx):
        xi = i + 1
        for j in range(ny):
            yj = j + 1
            for k in range(nz):
                zk = k + 1
                dist = np.sqrt((xi-px)**2 + (yj-py)**2 + (zk-pz)**2)
                if dist >= bg_rmin:
                    vals.append(psi_3d[i, j, k])
    return float(np.mean(vals)) if vals else 0.0


def compute_rise_to_dhp(distances, psi_vals, r_vecs,
                        charge, epsilon, kt, kappa, L, n_shells):
    """
    Encuentra el punto con el valor minimo de psi, evalua la teoria DH-periodica
    en ese punto exacto, y devuelve el offset = psi_DH_per(r_min) - psi_min.
    Sumando este offset a todos los datos, el punto mas bajo toca la curva teorica.
    """
    idx_min = np.argmin(psi_vals)
    psi_min = psi_vals[idx_min]
    r_vec_min = r_vecs[idx_min]

    prefactor = charge / (4.0 * np.pi * epsilon * kt)

    def psi_scalar(d):
        d = max(float(d), 1e-10)
        return prefactor * (np.exp(-kappa * d) / d if kappa > 0 else 1.0 / d)

    dist0 = np.linalg.norm(r_vec_min)
    psi_theory = psi_scalar(dist0)
    for inx in range(-n_shells, n_shells + 1):
        for iny in range(-n_shells, n_shells + 1):
            for inz in range(-n_shells, n_shells + 1):
                if inx == 0 and iny == 0 and inz == 0:
                    continue
                r_img = np.array([inx * L, iny * L, inz * L])
                dr = r_vec_min - r_img
                d = np.linalg.norm(dr)
                if d > 1e-10:
                    psi_theory += psi_scalar(d)

    offset = psi_theory - psi_min
    print(f"  --rise-to-DHp: psi_min={psi_min:.6e} en r={distances[idx_min]:.3f}")
    print(f"    psi_DH_per en ese punto={psi_theory:.6e}, offset={offset:.6e}")
    return float(offset)


def extract_points(psi_3d, particle_pos, grid_size,
                   min_radius=0.5, max_radius=None,
                   direction=None, angle_tol_deg=15.0, return_rvecs=False,
                   background=0.0):
    """
    Extrae puntos de la malla con sus distancias y valores de potencial.
    El archivo psi de Ludwig ya esta en unidades adimensionales e*psi/kT,
    equivalente al prefactor q/(4*pi*eps*kT) de la teoria DH. No se divide por kT.
    Si background != 0, se resta ese valor a cada punto (correccion de fondo de contraiones).
    """
    nx, ny, nz = grid_size
    if max_radius is None:
        max_radius = min(nx, ny, nz) / 2.0

    cos_tol = np.cos(np.radians(angle_tol_deg)) if direction is not None else None
    px, py, pz = particle_pos

    distances, psi_vals, r_vecs = [], [], []

    for i in range(nx):
        xi = i + 1
        for j in range(ny):
            yj = j + 1
            for k in range(nz):
                zk = k + 1
                dr = np.array([xi - px, yj - py, zk - pz])
                dist = np.linalg.norm(dr)

                if dist < min_radius or dist > max_radius:
                    continue
                if direction is not None:
                    cos_angle = abs(np.dot(dr / dist, direction))
                    if cos_angle < cos_tol:
                        continue

                distances.append(dist)
                psi_vals.append(psi_3d[i, j, k] - background)
                r_vecs.append(dr)

    if return_rvecs:
        return np.array(distances), np.array(psi_vals), np.array(r_vecs)
    return np.array(distances), np.array(psi_vals)


# ---------------------------------------------------------------------------
# Grafico principal
# ---------------------------------------------------------------------------

def _sort_subsample_mask(dist, vals, vecs, max_pts, rmin):
    """Ordena por distancia, submuestrea y aplica rmin. Devuelve (dist, vals, vecs)."""
    idx = np.argsort(dist)
    dist = dist[idx]; vals = vals[idx]
    vecs = vecs[idx] if vecs is not None else None
    if len(dist) > max_pts:
        step = len(dist) // max_pts
        dist = dist[::step]; vals = vals[::step]
        vecs = vecs[::step] if vecs is not None else None
    mask = dist >= rmin
    dist = dist[mask]; vals = vals[mask]
    vecs = vecs[mask] if vecs is not None else None
    return dist, vals, vecs


def plot_psi(distances_total, psi_total,
             distances_real, psi_real,
             distances_fourier, psi_fourier,
             grid_size, charge, epsilon, kt, kappa, alpha, rc,
             rmin=0.1, output_file='psi_ewald.png', max_points=10000,
             log_scale=True, log_log=False, log_log_error=False,
             xmin=None, ymin=None,
             show_components=False, show_total=True,
             show_dh=True, show_ewald_theory=True,
             show_dh_periodic=False, dh_periodic_shells=1,
             direction=None, angle_tol_deg=15.0,
             r_vecs_total=None, pointwise_error=False,
             show_error_inf=True, show_error_dir=True,
             particle_radius=0.0, only_error_dh_per_rad=False,
             gaussian_source=False, sigma=1.0,
             gaussian_error=False, gaussian_error_periodic=False,
             show_gaussian_dh_inf=True):

    L = float(grid_size[0])

    # Ordenar, submuestrar, filtrar por rmin
    d_tot, p_tot, rv_tot = _sort_subsample_mask(
        distances_total, psi_total, r_vecs_total, max_points, rmin)
    d_real, p_real, _ = _sort_subsample_mask(
        distances_real, psi_real, None, max_points, rmin)
    d_four, p_four, _ = _sort_subsample_mask(
        distances_fourier, psi_fourier, None, max_points, rmin)

    all_dists = np.concatenate([d for d in [d_tot, d_real, d_four] if len(d) > 0])
    r_max_data = float(np.max(all_dists)) * 1.1 if len(all_dists) > 0 else 16.0
    r_max_theory = max(r_max_data, 16.0)
    r_theory = np.linspace(rmin, r_max_theory, 500)

    # Curvas teoricas
    psi_dh_inf = TheoreticalPotential.compute_dh(
        r_theory, q=charge, epsilon=epsilon, kt=kt, kappa=kappa)
    psi_theory_real, psi_theory_recip = TheoreticalPotential.compute_ewald_components(
        r_theory, q=charge, epsilon=epsilon, kt=kt, alpha=alpha, rc=rc)
    psi_theory_total = psi_theory_real + psi_theory_recip

    psi_dh_per = None
    if show_dh_periodic:
        print(f"  Calculando psi DH periodico: {dh_periodic_shells} capa(s), "
              f"dir={_direction_label(direction)}...")
        psi_dh_per = TheoreticalPotential.compute_dh_periodic(
            r_theory, q=charge, epsilon=epsilon, kt=kt, kappa=kappa,
            L=L, n_shells=dh_periodic_shells, direction=direction)

    # CHANGE INIT - DHFiniteRadius - Curvas teóricas con radio finito
    psi_dh_rad = None
    if particle_radius > 0:
        psi_dh_rad = TheoreticalPotential.compute_dh_finite_radius(
            r_theory, q=charge, epsilon=epsilon, kt=kt, kappa=kappa, a=particle_radius)

    psi_dh_per_rad = None
    if only_error_dh_per_rad and particle_radius > 0 and show_dh_periodic:
        n_img = (2 * dh_periodic_shells + 1) ** 3 - 1
        print(f"  Calculando psi DH periodico radio finito (a={particle_radius:.2f}): "
              f"{dh_periodic_shells} capa(s), {n_img} imagenes...")
        psi_dh_per_rad = TheoreticalPotential.compute_dh_periodic_finite_radius(
            d_tot, q=charge, epsilon=epsilon, kt=kt, kappa=kappa,
            a=particle_radius, L=L, n_shells=dh_periodic_shells, direction=direction)
    # CHANGE END - DHFiniteRadius

    # CHANGE INIT - GaussianDH - Cálculo de curvas teóricas con fuente gaussiana
    psi_gauss_dh = None
    psi_gauss_dh_per = None
    if gaussian_source:
        print(f"  Calculando psi gaussiana-DH (σ={sigma:.3f}) en r_theory...")
        psi_gauss_dh = TheoreticalPotential.compute_gaussian_dh(
            r_theory, q=charge, epsilon=epsilon, kt=kt, kappa=kappa, sigma=sigma)
        if show_dh_periodic and gaussian_error_periodic:
            n_img = (2 * dh_periodic_shells + 1) ** 3 - 1
            print(f"  Calculando psi gaussiana-DH periódico ({dh_periodic_shells} capa(s), "
                  f"{n_img} imgs) para error...")
            psi_gauss_dh_per = TheoreticalPotential.compute_gaussian_dh_periodic(
                r_theory, q=charge, epsilon=epsilon, kt=kt, kappa=kappa,
                sigma=sigma, L=L, n_shells=dh_periodic_shells, direction=direction)
    # CHANGE END - GaussianDH

    # Determinar si mostrar panel de error
    show_panel2 = show_total and (
        (show_dh and show_error_inf) or
        (show_dh_periodic and show_error_dir) or
        pointwise_error or
        (only_error_dh_per_rad and psi_dh_per_rad is not None) or
        (gaussian_source and gaussian_error) or
        (gaussian_source and gaussian_error_periodic and show_dh_periodic)
    )

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 11),
                                    gridspec_kw={'height_ratios': [2, 1]})

    dir_label_str = (f" dir={_direction_label(direction)} ±{angle_tol_deg:.0f}°"
                     if direction is not None else "")

    # Panel 1: psi vs r
    if show_components and len(d_real) > 0:
        ax1.plot(d_real, p_real, 'r.', markersize=4, alpha=0.7,
                 label='Ludwig: psi real (erfc)')
    if show_components and len(d_four) > 0:
        ax1.plot(d_four, p_four, 'b.', markersize=4, alpha=0.7,
                 label='Ludwig: psi Fourier')
    if show_total and len(d_tot) > 0:
        ax1.plot(d_tot, p_tot, 'g.', markersize=6, alpha=0.6,
                 label=f'Ludwig: psi total{dir_label_str}')

    if show_components:
        ax1.plot(r_theory, psi_theory_real, 'r-', linewidth=1.0, alpha=0.7,
                 label=f'Teoría: erfc(αr)/r (α={alpha:.2f}, rc={rc:.1f})')
        ax1.plot(r_theory, psi_theory_recip, 'b-', linewidth=1.0, alpha=0.7,
                 label=f'Teoría: erf(αr)/r (largo alcance, α={alpha:.2f})')
    if show_ewald_theory:
        ax1.plot(r_theory, psi_theory_total, 'm-', linewidth=1.5, alpha=0.8,
                 label=f'Teoría: Ewald total (real+recíp, α={alpha:.2f})')
    if show_dh:
        if kappa > 0:
            label_dh = f'Teoría: DH (λ_D={1.0/kappa:.2f}, κ={kappa:.3f})'
        else:
            label_dh = 'Teoría: Coulomb 1/r'
        ax1.plot(r_theory, psi_dh_inf, 'k-', linewidth=1.0, alpha=0.8, label=label_dh)
    if show_dh_periodic and psi_dh_per is not None:
        n_img = (2*dh_periodic_shells + 1)**3 - 1
        ax1.plot(r_theory, psi_dh_per, 'r--', linewidth=1.5, alpha=0.9,
                 label=f'DH periódico ({dh_periodic_shells} capa(s), {n_img} imgs)')
    # CHANGE INIT - DHFiniteRadius - Curvas radio finito en panel 1
    if psi_dh_rad is not None:
        ax1.plot(r_theory, psi_dh_rad, 'b--', linewidth=1.5, alpha=0.9,
                 label=f'DH radio finito (a={particle_radius:.2f} lu)')
        ax1.axvline(x=particle_radius, color='blue', linestyle=':', linewidth=1.0,
                    alpha=0.6, label=f'superficie (a={particle_radius:.2f})')
    # CHANGE END - DHFiniteRadius

    # CHANGE INIT - GaussianDH - Curva teórica gaussiana-DH en panel 1
    if gaussian_source and psi_gauss_dh is not None and show_gaussian_dh_inf:
        if kappa > 0:
            label_gauss = (f'Gauss-DH (σ={sigma:.2f}, λ_D={1.0/kappa:.2f}, κ={kappa:.3f})')
        else:
            label_gauss = f'Gauss-Coulomb (σ={sigma:.2f})'
        ax1.plot(r_theory, psi_gauss_dh, 'c-', linewidth=2.0, alpha=0.9, label=label_gauss)
    if gaussian_source and psi_gauss_dh_per is not None:
        n_img = (2 * dh_periodic_shells + 1) ** 3 - 1
        ax1.plot(r_theory, psi_gauss_dh_per, 'c--', linewidth=1.5, alpha=0.9,
                 label=f'Gauss-DH periódico ({dh_periodic_shells} cap., {n_img} imgs)')
    # CHANGE END - GaussianDH

    ax1.axvline(x=rc, color='orange', linestyle='--', linewidth=2, alpha=0.7,
                label=f'rc={rc:.1f}')
    if kappa > 0:
        ax1.axvline(x=1.0/kappa, color='purple', linestyle=':', linewidth=2, alpha=0.7,
                    label=f'λ_D={1.0/kappa:.2f}')

    # Ajuste de ylim basado en datos de Ludwig
    ludwig_vals = np.concatenate([v for v in [p_tot, p_real, p_four] if len(v) > 0])
    if len(ludwig_vals) > 0:
        pos_lud = ludwig_vals[ludwig_vals > 0]
        y_top = float(np.max(ludwig_vals)) * 2.0
        y_bot = float(np.min(pos_lud)) * 0.5 if len(pos_lud) > 0 else None
    else:
        y_top = y_bot = None

    ax1.set_xlabel('r (lu)', fontsize=font_size)
    ax1.set_ylabel('e·ψ/kT  (adimensional)', fontsize=font_size)
    title = (f'Potencial electrostatico Ewald  α={alpha:.3f}  rc={rc:.1f}  κ={kappa:.4f}\n'
             f'L={grid_size[0]}×{grid_size[1]}×{grid_size[2]}, q={charge:.2e}, ε={epsilon:.2e}, kT={kt:.2e}')
    if direction is not None:
        title += f'\nDirección: {_direction_label(direction)} (±{angle_tol_deg:.0f}°)'
    ax1.set_title(title, fontsize=font_size)
    ax1.legend(fontsize=font_size - 3, loc='best')
    ax1.grid(True, alpha=0.3)

    if log_log:
        ax1.set_xscale('log'); ax1.set_yscale('log')
        x_left = xmin if xmin is not None else max(rmin, 0.1)
        ax1.set_xlim(x_left, r_max_data)
        ax1.set_ylim(bottom=ymin if ymin is not None else y_bot, top=y_top)
    elif log_scale:
        ax1.set_yscale('log')
        ax1.set_xlim(xmin if xmin is not None else 0.0, r_max_data)
        ax1.set_ylim(bottom=ymin if ymin is not None else y_bot, top=y_top)
    else:
        ax1.set_xlim(xmin if xmin is not None else 0.0, r_max_data)
        if ymin is not None:
            ax1.set_ylim(bottom=ymin, top=y_top)

    # Panel 2: error relativo del total
    if show_panel2:
        # CHANGE INIT - DHFiniteRadius - Condicionar otras curvas de error con only_error_dh_per_rad
        if not only_error_dh_per_rad:
            if show_dh and show_error_inf:
                psi_inf_at_pts = TheoreticalPotential.compute_dh(
                    d_tot, q=charge, epsilon=epsilon, kt=kt, kappa=kappa)
                rel_gap_inf = np.abs(p_tot - psi_inf_at_pts) / (np.abs(psi_inf_at_pts) + 1e-15) * 100
                ax2.plot(d_tot, rel_gap_inf, 'purple', markersize=3.5, alpha=0.6,
                         marker='.', linestyle='none', label='|psi - DH inf| / DH inf × 100')

            if particle_radius > 0:
                psi_rad_at_pts = TheoreticalPotential.compute_dh_finite_radius(
                    d_tot, q=charge, epsilon=epsilon, kt=kt, kappa=kappa, a=particle_radius)
                valid_rad = ~np.isnan(psi_rad_at_pts)
                rel_gap_rad = (np.abs(p_tot[valid_rad] - psi_rad_at_pts[valid_rad])
                               / (np.abs(psi_rad_at_pts[valid_rad]) + 1e-15) * 100)
                ax2.plot(d_tot[valid_rad], rel_gap_rad, 'b', markersize=3.5, alpha=0.7,
                         marker='.', linestyle='none',
                         label=f'|psi - DH a={particle_radius:.2f}| / DH_a × 100')

            if show_dh_periodic and show_error_dir and psi_dh_per is not None:
                psi_per_at_pts = TheoreticalPotential.compute_dh_periodic(
                    d_tot, q=charge, epsilon=epsilon, kt=kt, kappa=kappa,
                    L=L, n_shells=dh_periodic_shells, direction=direction)
                rel_gap_per = np.abs(p_tot - psi_per_at_pts) / (np.abs(psi_per_at_pts) + 1e-15) * 100
                ax2.plot(d_tot, rel_gap_per, 'darkorange', markersize=3.5, alpha=0.7,
                         marker='.', linestyle='none',
                         label='|psi - DH per (dir)| / DH × 100')

            if pointwise_error and rv_tot is not None and len(rv_tot) > 0:
                print(f"  Calculando error puntual en {len(rv_tot)} nodos...")
                psi_pw = TheoreticalPotential.compute_dh_periodic_at_points(
                    rv_tot, q=charge, epsilon=epsilon, kt=kt, kappa=kappa,
                    L=L, n_shells=dh_periodic_shells)
                rel_gap_pw = np.abs(p_tot - psi_pw) / (np.abs(psi_pw) + 1e-15) * 100
                ax2.plot(d_tot, rel_gap_pw, 'green', markersize=5, alpha=0.8,
                         marker='.', linestyle='none',
                         label='|psi - DH per (punto exacto)| / DH × 100')

        if psi_dh_per_rad is not None:
            n_img = (2 * dh_periodic_shells + 1) ** 3 - 1
            valid_pr = ~np.isnan(psi_dh_per_rad)
            rel_gap_per_rad = (np.abs(p_tot[valid_pr] - psi_dh_per_rad[valid_pr])
                               / (np.abs(psi_dh_per_rad[valid_pr]) + 1e-15) * 100)
            ax2.plot(d_tot[valid_pr], rel_gap_per_rad, 'darkgreen',
                     markersize=3.5, alpha=0.8,
                     marker='.', linestyle='none',
                     label=(f'|psi - DH per a={particle_radius:.2f}| / DH_per_a × 100'
                            f' ({dh_periodic_shells} capa(s), {n_img} imgs)'))
        # CHANGE END - DHFiniteRadius

        # CHANGE INIT - GaussianDH - Curvas de error vs gaussiana-DH en panel 2
        if gaussian_source and gaussian_error:
            psi_gauss_at_pts = TheoreticalPotential.compute_gaussian_dh(
                d_tot, q=charge, epsilon=epsilon, kt=kt, kappa=kappa, sigma=sigma)
            rel_gap_gauss = (np.abs(p_tot - psi_gauss_at_pts)
                             / (np.abs(psi_gauss_at_pts) + 1e-15) * 100)
            if kappa > 0:
                label_ge = f'|psi - Gauss-DH (σ={sigma:.2f})| / Gauss-DH × 100'
            else:
                label_ge = f'|psi - Gauss-Coulomb (σ={sigma:.2f})| / Gauss-Coulomb × 100'
            ax2.plot(d_tot, rel_gap_gauss, 'teal', markersize=3.5, alpha=0.8,
                     marker='.', linestyle='none', label=label_ge)

        if gaussian_source and gaussian_error_periodic and show_dh_periodic:
            print(f"  Calculando error vs gaussiana-DH periódica en {len(d_tot)} nodos...")
            if rv_tot is not None and len(rv_tot) > 0:
                psi_gauss_per_pts = TheoreticalPotential.compute_gaussian_dh_periodic_at_points(
                    rv_tot, q=charge, epsilon=epsilon, kt=kt, kappa=kappa,
                    sigma=sigma, L=L, n_shells=dh_periodic_shells)
            else:
                psi_gauss_per_pts = TheoreticalPotential.compute_gaussian_dh_periodic(
                    d_tot, q=charge, epsilon=epsilon, kt=kt, kappa=kappa,
                    sigma=sigma, L=L, n_shells=dh_periodic_shells, direction=direction)
            rel_gap_gauss_per = (np.abs(p_tot - psi_gauss_per_pts)
                                 / (np.abs(psi_gauss_per_pts) + 1e-15) * 100)
            n_img = (2 * dh_periodic_shells + 1) ** 3 - 1
            ax2.plot(d_tot, rel_gap_gauss_per, 'darkcyan', markersize=4, alpha=0.85,
                     marker='.', linestyle='none',
                     label=(f'|psi - Gauss-DH per (σ={sigma:.2f})| / Gauss-DH_per × 100'
                            f' ({dh_periodic_shells} cap., {n_img} imgs)'))
        # CHANGE END - GaussianDH

        ax2.axhline(y=10, color='r', linestyle=':', alpha=0.5, linewidth=2, label='10%')
        ax2.axhline(y=1,  color='g', linestyle=':', alpha=0.5, linewidth=2, label='1%')
        ax2.axvline(x=rc, color='orange', linestyle='--', linewidth=2, alpha=0.7,
                    label=f'rc={rc:.1f}')
        ax2.set_xlabel('r (lu)', fontsize=font_size)
        ax2.set_ylabel('Error relativo (%)', fontsize=font_size)
        ax2.set_title('Error relativo psi total vs. teoría', fontsize=font_size)
        ax2.legend(fontsize=font_size - 4, loc='upper right')
        ax2.grid(True, alpha=0.3)
        if log_log_error:
            ax2.set_yscale('log')
        ax2.set_xlim(ax1.get_xlim())
    else:
        ax2.set_visible(False)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nGrafico guardado en: {output_file}")
    plt.close()


# ---------------------------------------------------------------------------
# Exportar CSV
# ---------------------------------------------------------------------------

def export_csv_psi(psi_3d, particle_pos, grid_size, kt,
                   charge, epsilon, kappa, L, n_shells, rmin, csv_file):
    """Exporta todos los nodos con psi/kT y error vs DH-inf y DH-periodico"""
    nx, ny, nz = grid_size
    px, py, pz = particle_pos
    prefactor = charge / (4.0 * np.pi * epsilon * kt)

    def psi_scalar(d):
        d = max(float(d), 1e-10)
        if kappa > 0:
            return prefactor * np.exp(-kappa * d) / d
        else:
            return prefactor / d

    def psi_periodic(r_vec):
        dist0 = np.linalg.norm(r_vec)
        val = psi_scalar(dist0)
        for inx in range(-n_shells, n_shells + 1):
            for iny in range(-n_shells, n_shells + 1):
                for inz in range(-n_shells, n_shells + 1):
                    if inx == 0 and iny == 0 and inz == 0:
                        continue
                    r_img = np.array([inx * L, iny * L, inz * L])
                    dr = r_vec - r_img
                    d = np.linalg.norm(dr)
                    if d > 1e-10:
                        val += psi_scalar(d)
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

                psi_lud = psi_3d[i, j, k]
                psi_inf = psi_scalar(dist)
                psi_per = psi_periodic(r_vec)

                rows.append((dist, i, j, k, xi, yj, zk,
                             psi_lud, psi_inf, psi_per,
                             rel_err(psi_lud, psi_inf),
                             rel_err(psi_lud, psi_per)))

    rows.sort(key=lambda r: r[0])

    header = ("px,py,pz,dist,ix,iy,iz,x,y,z,"
              "psi_ludwig,psi_dh_inf,psi_dh_per,"
              "err_inf%,err_per%")
    with open(csv_file, 'w') as f:
        f.write(header + '\n')
        for r in rows:
            f.write(
                f"{px:.6f},{py:.6f},{pz:.6f},"
                f"{r[0]:.8e},{r[1]:d},{r[2]:d},{r[3]:d},"
                f"{r[4]:.6f},{r[5]:.6f},{r[6]:.6f},"
                f"{r[7]:.8e},{r[8]:.8e},{r[9]:.8e},"
                f"{r[10]:.4f},{r[11]:.4f}\n"
            )

    print(f"CSV exportado: {csv_file} ({len(rows)} puntos)")


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Grafica el potencial psi de Ewald (total, real, Fourier) vs. Debye-Huckel',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Ejemplos:
    python plot_psi_ewald.py -n 2000 -L 32 -q 1.0 -kt 0.00001 -eps 1e4 \\
        --rho-el 8e-3 --alpha 0.7 --rc 3.0

    python plot_psi_ewald.py -n 2000 -L 32 -q 1.0 -kt 0.00001 -eps 1e4 \\
        --rho-el 8e-3 --alpha 0.7 --rc 3.0 --show-components

    python plot_psi_ewald.py -n 2000 -L 32 -q 1.0 -kt 0.00001 -eps 1e4 \\
        --rho-el 8e-3 --alpha 0.7 --rc 3.0 \\
        --show-dh-periodic --dh-periodic-shells 1 \\
        --direction 1,5,5 --angle-tol 5 --pointwise-error --log-log
""")

    parser.add_argument('-n', '--nstep', type=int, required=True,
                        help='Numero de paso de tiempo')
    parser.add_argument('-L', '--grid-size', type=int, required=True,
                        help='Tamano de la malla cubica L x L x L')
    parser.add_argument('-q', '--charge', type=float, default=1.0,
                        help='Carga de la particula')
    parser.add_argument('-kt', '--kt', type=float, default=1.0,
                        help='Energia termica kT')
    parser.add_argument('-eps', '--epsilon', type=float, default=1.0,
                        help='Permitividad relativa')
    parser.add_argument('--alpha', type=float, default=0.5,
                        help='Parametro de splitting de Ewald alpha')
    parser.add_argument('--rc', type=float, default=5.0,
                        help='Radio de corte del espacio real')
    parser.add_argument('--kappa', type=float, default=None,
                        help='Parametro de Debye kappa (si no se da, se calcula desde --rho-el)')
    parser.add_argument('--rho-el', type=float, default=None,
                        help='Densidad de electrolito para calcular kappa')
    parser.add_argument('-d', '--directory', default='.',
                        help='Directorio con los archivos de Ludwig (default: .)')

    parser.add_argument('--rmin', type=float, default=0.5,
                        help='Radio minimo a graficar (default: 0.5)')
    parser.add_argument('--rmax', type=float, default=None,
                        help='Radio maximo (default: L/2)')

    parser.add_argument('--show-components', action='store_true', default=False,
                        help='Mostrar psi_real y psi_fourier como series separadas')
    parser.add_argument('--no-total', action='store_true', default=False,
                        help='Ocultar psi total de Ludwig')
    parser.add_argument('--show-ewald-theory', action='store_true', default=None,
                        help='Mostrar curva teorica Ewald total (erfc+erf)/r')
    parser.add_argument('--show-dh', action='store_true', default=None,
                        help='Mostrar curva DH infinita')
    parser.add_argument('--show-dh-periodic', action='store_true', default=None,
                        help='Mostrar curva DH periodica')
    parser.add_argument('--dh-periodic-shells', type=int, default=1,
                        help='Capas de imagenes periodicas (default: 1)')
    parser.add_argument('--pointwise-error', action='store_true', default=False,
                        help='Calcular error en la posicion exacta de cada nodo (verde)')
    parser.add_argument('--no-error-inf', action='store_true', default=False,
                        help='Ocultar serie de error vs DH infinito (purpura)')
    parser.add_argument('--no-error-dir', action='store_true', default=False,
                        help='Ocultar serie de error vs DH periodico direccional (naranja)')
    # CHANGE INIT - DHFiniteRadius - Argumento radio de partícula y solo-error periódico con radio
    parser.add_argument('--particle-radius', type=float, default=0.0,
                        help='Radio de la partícula a (lu). Activa curva DH radio finito '
                             'y error vs DH(a) en el panel de error (default: 0 = desactivado)')
    parser.add_argument('--only-error-dh-per-rad', action='store_true', default=False,
                        help='Mostrar en el panel de error SOLO la curva de error relativo '
                             'vs DH periódico con partícula de radio a. Requiere '
                             '--particle-radius > 0 y --show-dh-periodic. '
                             'Suprime las demás curvas de error.')
    # CHANGE END - DHFiniteRadius

    # CHANGE INIT - GaussianDH - Argumentos para fuente de carga gaussiana
    parser.add_argument('--gaussian-source', action='store_true', default=False,
                        help='Usar distribución de carga gaussiana ρ(r)=Q/(π^(3/2)σ³)exp(-r²/σ²) '
                             'como fuente en la teoría DH (en lugar de carga puntual). '
                             'Activa la curva teórica ψ_gauss-DH en el panel 1.')
    parser.add_argument('--sigma', type=float, default=1.0,
                        help='Anchura σ de la distribución gaussiana de carga (lu). '
                             'Solo se usa con --gaussian-source (default: 1.0)')
    parser.add_argument('--gaussian-error', action='store_true', default=False,
                        help='Mostrar el error relativo vs la teoría gaussiana-DH (sistema infinito) '
                             'en el panel 2. Requiere --gaussian-source.')
    parser.add_argument('--gaussian-error-periodic', action='store_true', default=False,
                        help='Mostrar el error relativo vs la teoría gaussiana-DH periódica '
                             'en el panel 2. Requiere --gaussian-source y --show-dh-periodic. '
                             'Usa --dh-periodic-shells para el número de capas.')
    parser.add_argument('--no-gaussian-dh-inf', action='store_true', default=False,
                        help='No graficar la curva gaussiana-DH sistema infinito (no periódica).')
    # CHANGE END - GaussianDH

    parser.add_argument('--linear', action='store_true', default=False,
                        help='Escala lineal en Y (default: semilog)')
    parser.add_argument('--log-log', action='store_true', default=False,
                        help='Escala log-log en ambos ejes')
    parser.add_argument('--log-log-error', action='store_true', default=False,
                        help='Escala log en eje Y del panel de error')
    parser.add_argument('--xmin', type=float, default=None)
    parser.add_argument('--ymin', type=float, default=None)

    parser.add_argument('--subtract-background', action='store_true', default=False,
                        help='Restar el potencial de fondo de contraiones (promedio en nodos lejanos)')
    parser.add_argument('--bg-rmin', type=float, default=12.0,
                        help='Radio minimo para calcular el fondo (default: 12.0)')
    parser.add_argument('--rise-to-DHp', action='store_true', default=False,
                        help='Subir todos los valores de psi para que el minimo de los puntos '
                             'seleccionados coincida con la curva DH-periodica en ese punto. '
                             'Requiere --show-dh-periodic y --dh-periodic-shells.')

    parser.add_argument('--direction', default=None,
                        help="Filtrar por direccion: 'x','y','z','xy','xyz' o 'dx,dy,dz'")
    parser.add_argument('--angle-tol', type=float, default=15.0,
                        help='Tolerancia angular en grados (default: 15)')

    parser.add_argument('--output', default='psi_ewald',
                        help='Nombre base del archivo de salida (default: psi_ewald)')
    parser.add_argument('--max-points', type=int, default=10000,
                        help='Maximo numero de puntos a graficar (default: 10000)')
    parser.add_argument('--csv', default=None,
                        help='Exportar datos a CSV (nombre base; se anade -n_NNNNNN.csv)')

    args = parser.parse_args()

    # Parsear direccion
    try:
        direction = parse_direction(args.direction)
    except ValueError as e:
        print(f"ERROR: {e}")
        sys.exit(1)
    if direction is not None:
        print(f"Direccion seleccionada: {_direction_label(direction)} = {direction}, "
              f"tolerancia +/-{args.angle_tol:.1f} deg")

    # Flags de visualizacion: si no se especifica ninguno, mostrar DH inf y teoria Ewald
    if args.show_dh is None and args.show_dh_periodic is None and args.show_ewald_theory is None:
        args.show_dh = True
        args.show_dh_periodic = False
        args.show_ewald_theory = True
    if args.show_dh is None:
        args.show_dh = False
    if args.show_dh_periodic is None:
        args.show_dh_periodic = False
    if args.show_ewald_theory is None:
        args.show_ewald_theory = False

    # Calcular kappa
    if args.kappa is not None:
        kappa = args.kappa
    elif args.rho_el is not None:
        lb = 1.0 / (4.0 * np.pi * args.epsilon * args.kt)
        V = float(args.grid_size ** 3)
        rho_counterions = abs(args.charge) / V
        rho_total = args.rho_el + rho_counterions
        kappa = np.sqrt(4.0 * np.pi * lb * rho_total * 2.0)
        print(f"\nParametros de Debye-Huckel:")
        print(f"  Longitud de Bjerrum lambda_B = {lb:.6e}")
        print(f"  kappa = {kappa:.6e}  (lambda_D = {1.0/kappa:.4f} lu)")
    else:
        kappa = 0.0

    grid_size = (args.grid_size, args.grid_size, args.grid_size)
    base_dir = args.directory

    nstep_str9 = f"{args.nstep:09d}"
    nstep_str8 = f"{args.nstep:08d}"
    psi_file         = os.path.join(base_dir, f"psi-{nstep_str9}.001-001")
    psi_real_file    = os.path.join(base_dir, f"psi_real-{nstep_str9}.001-001")
    psi_fourier_file = os.path.join(base_dir, f"psi_fourier-{nstep_str9}.001-001")
    colloid_file     = os.path.join(base_dir, f"config.cds{nstep_str8}.001-001")

    for fname in [psi_file, colloid_file]:
        if not os.path.exists(fname):
            print(f"ERROR: No se encontro el archivo: {fname}")
            sys.exit(1)

    print(f"Leyendo posicion de la particula desde: {colloid_file}")
    particle_pos = ColloidReader.read_colloid_position_cds(colloid_file)
    print(f"  Posicion: {particle_pos}")

    # Leer psi total (obligatorio)
    print(f"\nLeyendo potencial total desde: {psi_file}")
    try:
        psi_3d = PsiReader(psi_file, grid_size).read()
        print(f"  Rango psi total: [{psi_3d.min():.6e}, {psi_3d.max():.6e}]")
    except Exception as e:
        print(f"ERROR: {e}"); sys.exit(1)

    # Leer psi_real (opcional)
    psi_real_3d = None
    if os.path.exists(psi_real_file):
        print(f"Leyendo psi_real desde: {psi_real_file}")
        try:
            psi_real_3d = PsiReader(psi_real_file, grid_size).read()
            print(f"  Rango psi_real: [{psi_real_3d.min():.6e}, {psi_real_3d.max():.6e}]")
        except Exception as e:
            print(f"AVISO: no se pudo leer psi_real: {e}")
    else:
        print(f"AVISO: no se encontro {psi_real_file} (se omitira componente real)")

    # Leer psi_fourier (opcional)
    psi_fourier_3d = None
    if os.path.exists(psi_fourier_file):
        print(f"Leyendo psi_fourier desde: {psi_fourier_file}")
        try:
            psi_fourier_3d = PsiReader(psi_fourier_file, grid_size).read()
            print(f"  Rango psi_fourier: [{psi_fourier_3d.min():.6e}, {psi_fourier_3d.max():.6e}]")
        except Exception as e:
            print(f"AVISO: no se pudo leer psi_fourier: {e}")
    else:
        print(f"AVISO: no se encontro {psi_fourier_file} (se omitira componente Fourier)")

    # Si se piden componentes pero no existen los archivos, desactivar
    if args.show_components and (psi_real_3d is None or psi_fourier_3d is None):
        print("AVISO: --show-components requiere psi_real y psi_fourier; se desactiva")
        args.show_components = False

    # Calcular y restar fondo si se solicita
    background_tot = 0.0
    background_real = 0.0
    background_four = 0.0
    if args.subtract_background:
        background_tot = compute_background(psi_3d, particle_pos, grid_size,
                                            bg_rmin=args.bg_rmin)
        print(f"\nPotencial de fondo (promedio r > {args.bg_rmin}):")
        print(f"  psi_total:   {background_tot:.6e}")
        if psi_real_3d is not None:
            background_real = compute_background(psi_real_3d, particle_pos, grid_size,
                                                 bg_rmin=args.bg_rmin)
            print(f"  psi_real:    {background_real:.6e}")
        if psi_fourier_3d is not None:
            background_four = compute_background(psi_fourier_3d, particle_pos, grid_size,
                                                 bg_rmin=args.bg_rmin)
            print(f"  psi_fourier: {background_four:.6e}")
        print(f"  (se restara a cada componente antes de comparar con DH)")

    print(f"\nExtrayendo puntos (rmin={args.rmin}, rmax={args.rmax})...")
    need_rvecs = args.pointwise_error

    # Total
    if need_rvecs:
        distances_tot, psi_tot, r_vecs_arr = extract_points(
            psi_3d, particle_pos, grid_size,
            min_radius=args.rmin, max_radius=args.rmax,
            direction=direction, angle_tol_deg=args.angle_tol,
            return_rvecs=True, background=background_tot)
    else:
        distances_tot, psi_tot = extract_points(
            psi_3d, particle_pos, grid_size,
            min_radius=args.rmin, max_radius=args.rmax,
            direction=direction, angle_tol_deg=args.angle_tol,
            background=background_tot)
        r_vecs_arr = None
    print(f"  Puntos total: {len(distances_tot)}")

    # --rise-to-DHp: calcular offset para que el minimo de psi_tot toque DH-periodico
    rise_offset = 0.0
    if getattr(args, 'rise_to_DHp', False):
        if not args.show_dh_periodic:
            print("AVISO: --rise-to-DHp requiere --show-dh-periodic; se ignora")
        elif need_rvecs and r_vecs_arr is not None and len(r_vecs_arr) > 0:
            print("\nCalculando rise-to-DHp...")
            rise_offset = compute_rise_to_dhp(
                distances_tot, psi_tot, r_vecs_arr,
                charge=args.charge, epsilon=args.epsilon, kt=args.kt, kappa=kappa,
                L=float(args.grid_size), n_shells=args.dh_periodic_shells)
            psi_tot = psi_tot + rise_offset
        else:
            print("AVISO: --rise-to-DHp requiere --pointwise-error (para tener r_vecs); se ignora")

    # Real y Fourier
    if psi_real_3d is not None:
        distances_real, psi_real_vals = extract_points(
            psi_real_3d, particle_pos, grid_size,
            min_radius=args.rmin, max_radius=args.rmax,
            direction=direction, angle_tol_deg=args.angle_tol,
            background=background_real)
        print(f"  Puntos real: {len(distances_real)}")
    else:
        distances_real = np.array([]); psi_real_vals = np.array([])

    if psi_fourier_3d is not None:
        distances_four, psi_four_vals = extract_points(
            psi_fourier_3d, particle_pos, grid_size,
            min_radius=args.rmin, max_radius=args.rmax,
            direction=direction, angle_tol_deg=args.angle_tol,
            background=background_four)
        print(f"  Puntos Fourier: {len(distances_four)}")
        if rise_offset != 0.0:
            psi_four_vals = psi_four_vals + rise_offset
    else:
        distances_four = np.array([]); psi_four_vals = np.array([])

    output_base, output_ext = os.path.splitext(args.output)
    if not output_ext:
        output_ext = '.png'
    output_file = f"{output_base}-n_{args.nstep:09d}{output_ext}"

    log_scale = not args.linear and not args.log_log
    print(f"\nGenerando grafico...")
    plot_psi(distances_tot, psi_tot,
             distances_real, psi_real_vals,
             distances_four, psi_four_vals,
             grid_size, args.charge, args.epsilon, args.kt,
             kappa, args.alpha, args.rc,
             rmin=args.rmin, output_file=output_file, max_points=args.max_points,
             log_scale=log_scale, log_log=args.log_log, log_log_error=args.log_log_error,
             xmin=args.xmin, ymin=args.ymin,
             show_components=args.show_components,
             show_total=not args.no_total,
             show_dh=args.show_dh, show_ewald_theory=args.show_ewald_theory,
             show_dh_periodic=args.show_dh_periodic,
             dh_periodic_shells=args.dh_periodic_shells,
             direction=direction, angle_tol_deg=args.angle_tol,
             r_vecs_total=r_vecs_arr, pointwise_error=args.pointwise_error,
             show_error_inf=not args.no_error_inf,
             show_error_dir=not args.no_error_dir,
             particle_radius=args.particle_radius,
             only_error_dh_per_rad=args.only_error_dh_per_rad,
             gaussian_source=args.gaussian_source,
             sigma=args.sigma,
             gaussian_error=args.gaussian_error,
             gaussian_error_periodic=args.gaussian_error_periodic,
             show_gaussian_dh_inf=not args.no_gaussian_dh_inf)

    if args.csv is not None:
        csv_base, csv_ext = os.path.splitext(args.csv)
        if not csv_ext:
            csv_ext = '.csv'
        csv_file = f"{csv_base}-n_{args.nstep:09d}{csv_ext}"
        print(f"\nExportando datos a CSV: {csv_file}")
        export_csv_psi(psi_3d, particle_pos, grid_size, args.kt,
                       charge=args.charge, epsilon=args.epsilon, kappa=kappa,
                       L=float(args.grid_size), n_shells=args.dh_periodic_shells,
                       rmin=args.rmin, csv_file=csv_file)

    print("\nAnalisis completado!")


if __name__ == '__main__':
    main()
