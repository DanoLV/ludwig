#!/usr/bin/env python3
"""
Script para plotear por separado las componentes Real y Recíproca del campo eléctrico
calculadas por Ewald en Ludwig.

Este script ayuda a diagnosticar el problema de desconexión entre las dos partes.

Para usar este script, necesitamos modificar temporalmente ewald_charge.c para
que genere dos archivos separados:
  - efield_real: Componente de espacio real
  - efield_fourier: Componente de espacio recíproco

Autor: Script para análisis de Ewald
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
import os
from pathlib import Path

font_size = 16

class EFieldReader:
    """Lee y procesa archivos efield de Ludwig"""

    def __init__(self, filename, grid_size):
        """
        Inicializa el lector de archivos efield

        Args:
            filename: Ruta al archivo efield
            grid_size: Tupla (nx, ny, nz) con el tamaño de la malla
        """
        self.filename = filename
        self.nx, self.ny, self.nz = grid_size
        self.Ex = None
        self.Ey = None
        self.Ez = None

    def read(self):
        """Lee el archivo efield y lo organiza en mallas 3D para cada componente"""
        data = np.loadtxt(self.filename)

        # El archivo efield tiene 3 componentes por fila (Ex, Ey, Ez por línea)
        n_points = self.nx * self.ny * self.nz
        expected_rows = n_points

        # np.loadtxt devuelve shape (n_points, 3) si hay 3 columnas
        if data.ndim == 1:
            # Si es un array 1D, reshapear a (n_points, 3)
            if len(data) != n_points * 3:
                raise ValueError(f"El archivo tiene {len(data)} valores, "
                               f"pero se esperaban {n_points * 3} (3 componentes por punto)")
            data = data.reshape((n_points, 3))
        elif data.shape[0] != expected_rows or data.shape[1] != 3:
            raise ValueError(f"El archivo tiene forma {data.shape}, "
                           f"pero se esperaba ({expected_rows}, 3)")

        # Extraer componentes: cada columna es Ex, Ey, Ez
        Ex_flat = data[:, 0]
        Ey_flat = data[:, 1]
        Ez_flat = data[:, 2]

        # Reshape a malla 3D: Ludwig escribe en orden z-major
        # Para cada x, para cada y, todos los z
        self.Ex = Ex_flat.reshape((self.nx, self.ny, self.nz))
        self.Ey = Ey_flat.reshape((self.nx, self.ny, self.nz))
        self.Ez = Ez_flat.reshape((self.nx, self.ny, self.nz))

        return self.Ex, self.Ey, self.Ez

    def compute_magnitude(self):
        """Calcula el módulo del campo eléctrico |E| = sqrt(Ex² + Ey² + Ez²)"""
        if self.Ex is None or self.Ey is None or self.Ez is None:
            raise ValueError("Primero debe leer el archivo con read()")

        return np.sqrt(self.Ex**2 + self.Ey**2 + self.Ez**2)


class ColloidReader:
    """Lee archivos de coloides de Ludwig"""

    @staticmethod
    def read_colloid_position(filename):
        """
        Lee la posición del coloide desde un archivo colloids-*.csv

        Args:
            filename: Ruta al archivo colloids-*.csv

        Returns:
            numpy array con posición [x, y, z]
        """
        with open(filename, 'r') as f:
            # Saltar líneas de comentario y encabezados
            for line in f:
                if line.startswith('#'):
                    continue
                if 'id' in line.lower() or 'index' in line.lower():
                    continue

                # Primera línea de datos: index,x,y,z,...
                parts = line.strip().split(',')
                if len(parts) >= 4:
                    x = float(parts[1])
                    y = float(parts[2])
                    z = float(parts[3])
                    return np.array([x, y, z])

        raise ValueError(f"No se pudo leer posición del coloide en {filename}")

    @staticmethod
    def read_colloid_position_cds(filename):
        """
        Lee la posición del coloide desde un archivo config.cds*.001-001 de Ludwig.

        El archivo tiene un coloide por bloque; la posición x,y,z está en la
        línea 36 (3 valores separados por espacios).

        Args:
            filename: Ruta al archivo config.cds*.001-001

        Returns:
            numpy array con posición [x, y, z]
        """
        with open(filename, 'r') as f:
            lines = f.readlines()

        # Línea 36 es índice 35 (base 0); contiene "x  y  z"
        pos_line = lines[35].strip()
        parts = pos_line.split()
        if len(parts) < 3:
            raise ValueError(f"Línea 36 de {filename} no tiene 3 valores: '{pos_line}'")
        return np.array([float(parts[0]), float(parts[1]), float(parts[2])])


class TheoreticalField:
    """Calcula el campo eléctrico teórico según modelo de Debye-Hückel"""

    @staticmethod
    def compute_theory_for_distances(distances, q=1.0, epsilon=1.0, kt=1.0, kappa=0.0):
        """
        Calcula el campo eléctrico teórico para un array de distancias usando Debye-Hückel

        Args:
            distances: Array de distancias radiales
            q: Carga
            epsilon: Permitividad relativa
            kt: Energía térmica k_B*T (INCLUIDO en el prefactor para Debye-Hückel)
            kappa: Parámetro de Debye κ = 1/λ_D (default: 0.0 = Coulomb)

        Returns:
            E_mag: Módulo del campo eléctrico para cada distancia

        Nota: Para Debye-Hückel con iones, el campo eléctrico es:
              E = (q/4πε*kT) · exp(-κr) · (κ/r + 1/r²) · r̂
        """
        # Prefactor (CON kT, como se espera en Debye-Hückel)
        prefactor = q / (4.0 * np.pi * epsilon * kt)

        # Evitar división por cero
        r = np.maximum(distances, 0.01)

        # Campo eléctrico según el modelo
        if kappa > 0:
            # Debye-Hückel: medio con iones
            # |E| = (q/4πε*kT) · exp(-κr) · (κ/r + 1/r²)
            E_mag = prefactor * np.exp(-kappa * r) * (kappa / r + 1.0 / r**2)
        else:
            # Coulomb: vacío sin screening
            # |E| = (q/4πε*kT) · 1/r²
            E_mag = prefactor / r**2

        return E_mag

    @staticmethod
    def compute_dh_periodic(distances, q=1.0, epsilon=1.0, kt=1.0, kappa=0.0, L=32.0, n_shells=1,
                            direction=None):
        """
        Calcula el módulo del campo eléctrico de Debye-Hückel incluyendo imágenes periódicas,
        evaluado a lo largo de una dirección específica.

        El punto de evaluación se fija en r_vec = r · d̂, donde d̂ es el vector unitario
        de dirección. Las imágenes de la partícula están en r_img = n·L (n entero 3D).
        El campo vectorial total es:

            E_vec(r) = E_central(r) · d̂  +  Σ_{n≠0}  E(|r·d̂ - n·L|) · (r·d̂ - n·L)/|r·d̂ - n·L|

        Se devuelve |E_vec_total|.

        Args:
            distances: Array de distancias desde la partícula central
            q: Carga
            epsilon: Permitividad relativa
            kt: Energía térmica kT
            kappa: Parámetro de Debye κ = 1/λ_D (0 = Coulomb puro)
            L: Tamaño de la caja periódica (cubo L×L×L)
            n_shells: Capas de imágenes (1 = 26 vecinos, 2 = 124, ...)
            direction: Vector unitario [dx,dy,dz] de evaluación (default: eje x = [1,0,0])

        Returns:
            E_mag: |E_total| para cada distancia
        """
        if direction is None:
            direction = np.array([1.0, 0.0, 0.0])
        d = np.asarray(direction, dtype=float)
        d = d / np.linalg.norm(d)  # asegurar unitario

        prefactor = q / (4.0 * np.pi * epsilon * kt)
        r = np.maximum(distances, 0.01)

        def dh_field_scalar(dist):
            """Magnitud escalar del campo DH a distancia dist"""
            if kappa > 0:
                return prefactor * np.exp(-kappa * dist) * (kappa / dist + 1.0 / dist**2)
            else:
                return prefactor / dist**2

        # Punto de evaluación: r_vec = r * d̂  (shape: (len(r), 3))
        r_vec = r[:, np.newaxis] * d[np.newaxis, :]   # (N, 3)

        # Campo vectorial acumulado: contribución de la partícula central
        E_vec = dh_field_scalar(r)[:, np.newaxis] * d[np.newaxis, :]  # (N, 3), apunta en d̂

        # Sumar imágenes periódicas
        for inx in range(-n_shells, n_shells + 1):
            for iny in range(-n_shells, n_shells + 1):
                for inz in range(-n_shells, n_shells + 1):
                    if inx == 0 and iny == 0 and inz == 0:
                        continue

                    # Posición de la imagen: r_img = (inx·L, iny·L, inz·L)
                    r_img = np.array([inx * L, iny * L, inz * L])  # (3,)

                    # Vector del punto al origen de la imagen: dr = r_vec - r_img
                    dr = r_vec - r_img[np.newaxis, :]   # (N, 3)
                    dist_img = np.linalg.norm(dr, axis=1)  # (N,)

                    # Magnitud del campo de la imagen
                    E_img = dh_field_scalar(dist_img)  # (N,)

                    # Contribución vectorial: E_img * dr/dist_img
                    E_vec += (E_img / dist_img)[:, np.newaxis] * dr  # (N, 3)

        return np.linalg.norm(E_vec, axis=1)   # |E_total| para cada r

    @staticmethod
    def compute_dh_periodic_at_points(r_vecs, q=1.0, epsilon=1.0, kt=1.0, kappa=0.0,
                                      L=32.0, n_shells=1):
        """
        Calcula |E_DH periódico| en puntos 3D arbitrarios r_vecs (relativos a la partícula).

        Para cada punto r_vec, suma el campo de la partícula central más el de todas
        sus imágenes periódicas:

            E_vec(r) = E(|r|)·r̂  +  Σ_{n≠0}  E(|r - n·L|)·(r - n·L)/|r - n·L|

        Args:
            r_vecs: Array (N, 3) de vectores posición desde la partícula
            q, epsilon, kt, kappa: parámetros físicos
            L: tamaño de caja (cubo L×L×L)
            n_shells: capas de imágenes

        Returns:
            E_mag: Array (N,) con |E_total| en cada punto
        """
        prefactor = q / (4.0 * np.pi * epsilon * kt)
        r_vecs = np.asarray(r_vecs, dtype=float)  # (N, 3)

        def dh_scalar(dist):
            d = np.maximum(dist, 1e-10)
            if kappa > 0:
                return prefactor * np.exp(-kappa * d) * (kappa / d + 1.0 / d**2)
            else:
                return prefactor / d**2

        dist0 = np.linalg.norm(r_vecs, axis=1)  # (N,)
        # Contribución partícula central: E·r̂
        with np.errstate(invalid='ignore', divide='ignore'):
            r_hat = np.where(dist0[:, np.newaxis] > 1e-10,
                             r_vecs / dist0[:, np.newaxis], 0.0)
        E_vec = dh_scalar(dist0)[:, np.newaxis] * r_hat  # (N, 3)

        # Sumar imágenes
        for inx in range(-n_shells, n_shells + 1):
            for iny in range(-n_shells, n_shells + 1):
                for inz in range(-n_shells, n_shells + 1):
                    if inx == 0 and iny == 0 and inz == 0:
                        continue
                    r_img = np.array([inx * L, iny * L, inz * L])
                    dr = r_vecs - r_img[np.newaxis, :]       # (N, 3)
                    dist_img = np.linalg.norm(dr, axis=1)    # (N,)
                    E_img = dh_scalar(dist_img)              # (N,)
                    E_vec += (E_img / np.maximum(dist_img, 1e-10))[:, np.newaxis] * dr

        return np.linalg.norm(E_vec, axis=1)  # (N,)

    # CHANGE INIT - DHFiniteRadius - |E| DH con partícula de radio finito a
    @staticmethod
    def compute_dh_finite_radius(distances, q=1.0, epsilon=1.0, kt=1.0, kappa=0.0, a=0.0):
        """
        Módulo del campo eléctrico DH para partícula esférica de radio 'a' (sistema infinito).

        Para r > a:
            |E(r; a)| = prefactor * exp(-κ(r-a)) / (r(1+κa)) * (κ + 1/r)
        donde prefactor = q / (4π ε kT).

        Para r ≤ a retorna NaN.
        """
        prefactor = q / (4.0 * np.pi * epsilon * kt)
        r = np.asarray(distances, dtype=float)
        result = np.full_like(r, np.nan)
        outside = r > a
        r_out = r[outside]
        if kappa > 0:
            result[outside] = (prefactor * np.exp(-kappa * (r_out - a))
                               / (r_out * (1.0 + kappa * a))
                               * (kappa + 1.0 / r_out))
        else:
            result[outside] = prefactor / r_out**2
        return result

    @staticmethod
    def compute_dh_periodic_finite_radius(distances, q=1.0, epsilon=1.0, kt=1.0, kappa=0.0,
                                          a=0.0, L=32.0, n_shells=1, direction=None):
        """
        |E_total| DH periódico con partícula de radio finito 'a', evaluado a lo largo
        de una dirección.

        Imagen central: usa la fórmula de radio finito (r > a).
        Imágenes periódicas: tratadas como cargas puntuales.

        Para r ≤ a retorna NaN.
        """
        if direction is None:
            direction = np.array([1.0, 0.0, 0.0])
        d = np.asarray(direction, dtype=float)
        d = d / np.linalg.norm(d)

        prefactor = q / (4.0 * np.pi * epsilon * kt)
        r = np.asarray(distances, dtype=float)
        r_vec = r[:, np.newaxis] * d[np.newaxis, :]  # (N, 3)

        def dh_field_scalar(dist):
            dist = np.maximum(dist, 1e-10)
            if kappa > 0:
                return prefactor * np.exp(-kappa * dist) * (kappa / dist + 1.0 / dist**2)
            else:
                return prefactor / dist**2

        outside = r > a
        r_out = r[outside]

        # Campo vectorial imagen central con radio finito (apunta en d̂)
        E_vec = np.zeros((len(r), 3))
        if kappa > 0:
            e_central = (prefactor * np.exp(-kappa * (r_out - a))
                         / (r_out * (1.0 + kappa * a))
                         * (kappa + 1.0 / r_out))
        else:
            e_central = prefactor / r_out**2
        E_vec[outside] = e_central[:, np.newaxis] * d[np.newaxis, :]

        # Imágenes periódicas (puntuales, vectoriales)
        for inx in range(-n_shells, n_shells + 1):
            for iny in range(-n_shells, n_shells + 1):
                for inz in range(-n_shells, n_shells + 1):
                    if inx == 0 and iny == 0 and inz == 0:
                        continue
                    r_img = np.array([inx * L, iny * L, inz * L])
                    dr = r_vec - r_img[np.newaxis, :]     # (N, 3)
                    dist_img = np.linalg.norm(dr, axis=1)  # (N,)
                    E_img = dh_field_scalar(dist_img)
                    E_vec[outside] += ((E_img / np.maximum(dist_img, 1e-10))
                                       [:, np.newaxis] * dr)[outside]

        result = np.full(len(r), np.nan)
        result[outside] = np.linalg.norm(E_vec[outside], axis=1)
        return result
    # CHANGE END - DHFiniteRadius

    # CHANGE INIT - GaussianDH - Campo eléctrico DH con fuente de carga gaussiana
    @staticmethod
    def compute_gaussian_dh(distances, q=1.0, epsilon=1.0, kt=1.0, kappa=0.0, sigma=1.0):
        """
        Módulo del campo eléctrico |E(r)| = -dψ/dr para una fuente gaussiana en DH.

        El potencial gaussiana-DH es:
            ψ(r) = Q/(4πε kT r) · ½ · exp(κ²σ²/4)
                   · [e^{-κr} erfc(r/σ - κσ/2) - e^{κr} erfc(r/σ + κσ/2)]

        El campo eléctrico se obtiene diferenciando numéricamente con paso fino,
        ya que la expresión analítica de -dψ/dr es larga aunque bien definida y finita.

        Comportamientos límite:
          r→0: |E(0)| finito (sin singularidad)
          r≫σ: |E| → E_DH estándar = prefactor·e^{-κr}·(κ/r + 1/r²)
          κ→0: |E| → E_Coulomb suavizado por gaussiana

        Args:
            distances: array de distancias radiales
            q:        carga total Q
            epsilon:  permitividad relativa
            kt:       energía térmica kT
            kappa:    parámetro de Debye κ (0 = Coulomb puro)
            sigma:    anchura de la distribución gaussiana σ
        """
        from scipy.special import erf
        prefactor = q / (4.0 * np.pi * epsilon * kt)
        r = np.asarray(distances, dtype=float)
        r_safe = np.maximum(r, 1e-10)

        if kappa <= 0:
            # κ=0: E = -d/dr [prefactor*erf(r/σ)/r]
            #     = prefactor * [erf(r/σ)/r² - 2/(sqrt(π)*σ*r)*exp(-r²/σ²)]
            from scipy.special import erf
            e_mag = prefactor * (erf(r_safe / sigma) / r_safe**2
                                 - 2.0 / (np.sqrt(np.pi) * sigma * r_safe)
                                 * np.exp(-r_safe**2 / sigma**2))
            return np.abs(e_mag)

        # Calcular E analíticamente desde la función de Green de Yukawa:
        # ψ(r) = (pref/(κr)) * [e^{-κr}·I_<(r) + sinh(κr)·I_>(r)]
        # E(r) = -dψ/dr = (pref/(κr)) * [κe^{-κr}·I_< + κcosh(κr)·I_> - r·ρ̃(r)·(e^{-κr}·sinh(κr)/κ - sinh(κr)·e^{-κr}/κ)]
        #   simplificando: los términos con r·ρ̃(r) se cancelan
        # E(r) = -dψ/dr = ψ(r)/r + (pref/(κr)) * [κe^{-κr}·I_<(r) + κcosh(κr)·I_>(r) - κψ_r_term]
        # La derivada exacta es:
        # dψ/dr = -(pref/(κr²))*[e^{-κr}·I_< + sinh(κr)·I_>]
        #          + (pref/(κr))*[-κe^{-κr}·I_< + κcosh(κr)·I_> + r·ρ̃(r)*(sinh(κr)·(-1) + sinh(κr))]
        # El término con r·ρ̃(r) cancela exactamente. Queda:
        # dψ/dr = -ψ(r)/r + (pref/(κr)) * κ * [-e^{-κr}·I_< + cosh(κr)·I_>]

        rmax = max(15.0 * sigma, float(np.max(r_safe)) * 1.5)
        n_rho = 50000
        rp = np.linspace(1e-8, rmax, n_rho)
        drp = rp[1] - rp[0]
        rho_rp = 1.0 / (np.pi**1.5 * sigma**3) * np.exp(-rp**2 / sigma**2) * rp

        integrand_inner = np.sinh(kappa * rp) * rho_rp
        integrand_outer = np.exp(-kappa * rp) * rho_rp

        cumul_inner = np.cumsum(integrand_inner) * drp
        total_outer = np.sum(integrand_outer) * drp
        cumul_outer = total_outer - np.cumsum(integrand_outer) * drp

        I_lt = np.interp(r_safe, rp, cumul_inner)   # ∫₀ʳ sinh(κr')r'ρ̃dr'
        I_gt = np.interp(r_safe, rp, cumul_outer)   # ∫ᵣ^∞ e^{-κr'}r'ρ̃dr'

        # Factor 4π de la reducción de la integral 3D a 1D en simetría esférica
        fac = 4.0 * np.pi * prefactor
        psi = (fac / (kappa * r_safe)) * (
            np.exp(-kappa * r_safe) * I_lt + np.sinh(kappa * r_safe) * I_gt
        )
        # dψ/dr = -ψ/r + (4π·pref/r) * [-e^{-κr}·I_< + cosh(κr)·I_>]
        dpsi_dr = (-psi / r_safe
                   + (fac / r_safe) * (
                       -np.exp(-kappa * r_safe) * I_lt
                       + np.cosh(kappa * r_safe) * I_gt
                   ))
        return np.abs(-dpsi_dr)

    @staticmethod
    def compute_gaussian_dh_periodic(distances, q=1.0, epsilon=1.0, kt=1.0, kappa=0.0,
                                     sigma=1.0, L=32.0, n_shells=1, direction=None):
        """
        |E_total| gaussiana-DH periódico a lo largo de una dirección.
        Imagen central con fuente gaussiana (derivada numérica de ψ_gauss);
        imágenes periódicas como Yukawa puntual vectorial.
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

        # Campo imagen central: E_central = -dψ_gauss/dr, apunta en d̂
        E_central = TheoreticalField.compute_gaussian_dh(
            r_safe, q=q, epsilon=epsilon, kt=kt, kappa=kappa, sigma=sigma)
        E_vec = E_central[:, np.newaxis] * d[np.newaxis, :]  # (N, 3)

        # Imágenes periódicas (Yukawa puntual vectorial)
        def dh_field_scalar(dist):
            dist = np.maximum(dist, 1e-10)
            if kappa > 0:
                return prefactor * np.exp(-kappa * dist) * (kappa / dist + 1.0 / dist**2)
            else:
                return prefactor / dist**2

        for inx in range(-n_shells, n_shells + 1):
            for iny in range(-n_shells, n_shells + 1):
                for inz in range(-n_shells, n_shells + 1):
                    if inx == 0 and iny == 0 and inz == 0:
                        continue
                    r_img = np.array([inx * L, iny * L, inz * L])
                    dr = r_vec - r_img[np.newaxis, :]
                    dist_img = np.linalg.norm(dr, axis=1)
                    E_img = dh_field_scalar(dist_img)
                    E_vec += (E_img / np.maximum(dist_img, 1e-10))[:, np.newaxis] * dr

        return np.linalg.norm(E_vec, axis=1)

    @staticmethod
    def compute_gaussian_dh_periodic_at_points(r_vecs, q=1.0, epsilon=1.0, kt=1.0,
                                               kappa=0.0, sigma=1.0, L=32.0, n_shells=1):
        """
        |E_total| gaussiana-DH periódico en puntos 3D arbitrarios r_vecs (N,3).
        Imagen central con fuente gaussiana; imágenes periódicas como Yukawa puntual.
        """
        prefactor = q / (4.0 * np.pi * epsilon * kt)
        r_vecs = np.asarray(r_vecs, dtype=float)
        dist0 = np.linalg.norm(r_vecs, axis=1)

        # Campo imagen central: radial con magnitud gaussiana-DH
        E_central = TheoreticalField.compute_gaussian_dh(
            dist0, q=q, epsilon=epsilon, kt=kt, kappa=kappa, sigma=sigma)
        with np.errstate(invalid='ignore', divide='ignore'):
            r_hat = np.where(dist0[:, np.newaxis] > 1e-10,
                             r_vecs / dist0[:, np.newaxis], 0.0)
        E_vec = E_central[:, np.newaxis] * r_hat

        def dh_field_scalar(dist):
            dist = np.maximum(dist, 1e-10)
            if kappa > 0:
                return prefactor * np.exp(-kappa * dist) * (kappa / dist + 1.0 / dist**2)
            else:
                return prefactor / dist**2

        for inx in range(-n_shells, n_shells + 1):
            for iny in range(-n_shells, n_shells + 1):
                for inz in range(-n_shells, n_shells + 1):
                    if inx == 0 and iny == 0 and inz == 0:
                        continue
                    r_img = np.array([inx * L, iny * L, inz * L])
                    dr = r_vecs - r_img[np.newaxis, :]
                    dist_img = np.linalg.norm(dr, axis=1)
                    E_img = dh_field_scalar(dist_img)
                    E_vec += (E_img / np.maximum(dist_img, 1e-10))[:, np.newaxis] * dr

        return np.linalg.norm(E_vec, axis=1)
    # CHANGE END - GaussianDH

    @staticmethod
    def compute_ewald_components(distances, q=1.0, epsilon=1.0, kt=1.0, alpha=0.5, rc=5.0):
        """
        Calcula las componentes teóricas de Ewald para un array de distancias

        Args:
            distances: Array de distancias radiales
            q: Carga
            epsilon: Permitividad relativa
            kt: Energía térmica k_B*T
            alpha: Parámetro de splitting de Ewald
            rc: Radio de corte

        Returns:
            E_real: Componente de espacio real  = prefactor * [erfc(αr)/r² + (2α/√π)*exp(-α²r²)/r]
                    → tiende a Coulomb cuando r→0, se aplica corte en rc
            E_recip: Componente recíproca (largo alcance) en espacio real:
                     = prefactor * [erf(αr)/r² - (2α/√π)*exp(-α²r²)/r]
                     → tiende a 0 cuando r→0 (largo alcance puro)
                     → tiende a Coulomb cuando r→∞ (erf→1, exp→0)
                     E_real + E_recip = E_coulomb = prefactor/r² (sin corte)
        """
        from scipy.special import erfc, erf

        # Prefactor (CON kT, como en Debye-Hückel)
        prefactor = q / (4.0 * np.pi * epsilon * kt)

        # Evitar división por cero
        r = np.maximum(distances, 0.01)

        # Componente de espacio real (parte con erfc)
        # E_real = (q/4πε*kT) · erfc(α·r)/r² · [1 + 2αr/√π · exp(-(αr)²)]
        alpha_r = alpha * r
        exp_term = np.exp(-alpha_r**2)
        real_term = erfc(alpha_r) / r**2 + (2.0 * alpha / np.sqrt(np.pi)) * exp_term / r
        E_real = prefactor * real_term

        # Componente recíproca expresada en espacio real:
        # E_recip(r) = E_coulomb - E_real_sin_corte
        #            = prefactor * [erf(αr)/r² - (2α/√π)*exp(-α²r²)/r]
        # Esta expresión tiende a 0 cuando r→0 (la singularidad se cancela),
        # y a prefactor/r² (Coulomb) para r grande donde erf(αr)→1 y exp→0.
        E_recip = prefactor * (erf(alpha_r) / r**2 - (2.0 * alpha / np.sqrt(np.pi)) * exp_term / r)

        # Aplicar corte en la parte real
        E_real = np.where(r > rc, 0.0, E_real)

        return E_real, E_recip


def _direction_label(direction):
    """Devuelve una string legible para el vector dirección dado."""
    if direction is None:
        return 'todas'
    d = np.asarray(direction)
    # Intentar reconocer casos canónicos
    canonical = {
        'x':   np.array([1., 0., 0.]),
        'y':   np.array([0., 1., 0.]),
        'z':   np.array([0., 0., 1.]),
        'xy':  np.array([1., 1., 0.]) / np.sqrt(2),
        'xz':  np.array([1., 0., 1.]) / np.sqrt(2),
        'yz':  np.array([0., 1., 1.]) / np.sqrt(2),
        'xyz': np.array([1., 1., 1.]) / np.sqrt(3),
    }
    for name, v in canonical.items():
        if np.allclose(d, v, atol=1e-6):
            return name
    return f'({d[0]:.3f},{d[1]:.3f},{d[2]:.3f})'


def parse_direction(direction_str):
    """
    Parsea el argumento --direction y devuelve un vector unitario normalizado.

    Formatos aceptados:
        'x'          → (1, 0, 0)
        'y'          → (0, 1, 0)
        'z'          → (0, 0, 1)
        'xy'         → (1, 1, 0) normalizado
        'xz'         → (1, 0, 1) normalizado
        'yz'         → (0, 1, 1) normalizado
        'xyz'        → (1, 1, 1) normalizado
        '1,0,0'      → vector arbitrario (se normaliza)
        '1,1,0'      → vector arbitrario (se normaliza)

    Returns:
        dir_vec: numpy array unitario [dx, dy, dz], o None si direction_str es None
    """
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


def extract_all_points(E_mag_3d, particle_pos, grid_size, min_radius=0.5, max_radius=None,
                       direction=None, angle_tol_deg=15.0, return_rvecs=False):
    """
    Extrae puntos de la malla con sus distancias y valores de campo.

    Si `direction` es None se devuelven todos los puntos.
    Si `direction` es un vector unitario [dx,dy,dz], solo se devuelven los puntos
    cuyo vector desde la partícula forma un ángulo ≤ angle_tol_deg con esa dirección.

    Args:
        E_mag_3d: Campo |E| en malla 3D
        particle_pos: Posición de la partícula [x, y, z]
        grid_size: Tupla (nx, ny, nz)
        min_radius: Radio mínimo (default: 0.5)
        max_radius: Radio máximo (default: None = L/2)
        direction: Vector unitario numpy [dx,dy,dz] o None
        angle_tol_deg: Tolerancia angular en grados (default: 15)
        return_rvecs: Si True, devuelve también array (N,3) de vectores r desde la partícula

    Returns:
        distances: Array (N,) con distancias de cada punto
        E_values: Array (N,) con |E| de cada punto
        r_vecs: Array (N,3) de vectores posición desde la partícula [solo si return_rvecs=True]
    """
    nx, ny, nz = grid_size

    if max_radius is None:
        max_radius = min(nx, ny, nz) / 2.0

    cos_tol = np.cos(np.radians(angle_tol_deg)) if direction is not None else None

    distances = []
    E_values = []
    r_vecs = []

    px, py, pz = particle_pos

    # Coordenada del nodo en Ludwig: ic = i+1
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
                E_values.append(E_mag_3d[i, j, k])
                r_vecs.append(dr)

    if return_rvecs:
        return np.array(distances), np.array(E_values), np.array(r_vecs)
    return np.array(distances), np.array(E_values)


def extract_all_points_full(efield_real_reader, efield_fourier_reader, efield_total_reader,
                             particle_pos, grid_size, kt,
                             min_radius=0.5, max_radius=None,
                             direction=None, angle_tol_deg=15.0):
    """
    Extrae todos los puntos con las componentes vectoriales completas de cada campo.

    Returns:
        list of dicts con claves: dist, ix, iy, iz, x, y, z, r_vec(3),
        ex_r, ey_r, ez_r, mag_r, ex_f, ey_f, ez_f, mag_f, ex_t, ey_t, ez_t, mag_t
    """
    nx, ny, nz = grid_size
    if max_radius is None:
        max_radius = min(nx, ny, nz) / 2.0

    cos_tol = np.cos(np.radians(angle_tol_deg)) if direction is not None else None
    px, py, pz = particle_pos

    rows = []
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

                ex_r = efield_real_reader.Ex[i, j, k] / kt if efield_real_reader else float('nan')
                ey_r = efield_real_reader.Ey[i, j, k] / kt if efield_real_reader else float('nan')
                ez_r = efield_real_reader.Ez[i, j, k] / kt if efield_real_reader else float('nan')

                ex_f = efield_fourier_reader.Ex[i, j, k] / kt if efield_fourier_reader else float('nan')
                ey_f = efield_fourier_reader.Ey[i, j, k] / kt if efield_fourier_reader else float('nan')
                ez_f = efield_fourier_reader.Ez[i, j, k] / kt if efield_fourier_reader else float('nan')

                ex_t = efield_total_reader.Ex[i, j, k] / kt if efield_total_reader else float('nan')
                ey_t = efield_total_reader.Ey[i, j, k] / kt if efield_total_reader else float('nan')
                ez_t = efield_total_reader.Ez[i, j, k] / kt if efield_total_reader else float('nan')

                rows.append({
                    'dist': dist, 'r_vec': dr,
                    'ix': i, 'iy': j, 'iz': k, 'x': xi, 'y': yj, 'z': zk,
                    'ex_r': ex_r, 'ey_r': ey_r, 'ez_r': ez_r,
                    'ex_f': ex_f, 'ey_f': ey_f, 'ez_f': ez_f,
                    'ex_t': ex_t, 'ey_t': ey_t, 'ez_t': ez_t,
                })

    rows.sort(key=lambda r: r['dist'])
    return rows


def export_csv_from_points(rows, particle_pos, charge, epsilon, kt, kappa, L, n_shells, csv_file):
    """
    Graba CSV con errores respecto a DH∞ y DH periódico, a partir de la lista
    de puntos devuelta por extract_all_points_full().
    Misma estructura que export_csv() pero calculada sobre los puntos exactos
    que se muestran en el gráfico.
    """
    prefactor = charge / (4.0 * np.pi * epsilon * kt)

    def dh_scalar(d):
        d = max(float(d), 1e-10)
        if kappa > 0:
            return prefactor * np.exp(-kappa * d) * (kappa / d + 1.0 / d**2)
        else:
            return prefactor / d**2

    def dh_vec_inf(r_vec):
        dist = np.linalg.norm(r_vec)
        if dist < 1e-10:
            return np.zeros(3)
        return dh_scalar(dist) * r_vec / dist

    def dh_vec_periodic(r_vec):
        E = dh_vec_inf(r_vec)
        for inx in range(-n_shells, n_shells + 1):
            for iny in range(-n_shells, n_shells + 1):
                for inz in range(-n_shells, n_shells + 1):
                    if inx == 0 and iny == 0 and inz == 0:
                        continue
                    r_img = np.array([inx * L, iny * L, inz * L])
                    dr = r_vec - r_img
                    dist_img = np.linalg.norm(dr)
                    if dist_img < 1e-10:
                        continue
                    E += dh_scalar(dist_img) * dr / dist_img
        return E

    def rel_err(val, ref):
        return abs(val - ref) / (abs(ref) + 1e-15) * 100.0

    px, py, pz = particle_pos

    header = ("px,py,pz,"
              "dist,ix,iy,iz,x,y,z,"
              "Ex_real/kT,Ey_real/kT,Ez_real/kT,|E_real|/kT,"
              "Ex_fourier/kT,Ey_fourier/kT,Ez_fourier/kT,|E_fourier|/kT,"
              "Ex_total/kT,Ey_total/kT,Ez_total/kT,|E_total|/kT,"
              "Ex_dh_inf/kT,Ey_dh_inf/kT,Ez_dh_inf/kT,|E_dh_inf|/kT,"
              "Ex_dh_per/kT,Ey_dh_per/kT,Ez_dh_per/kT,|E_dh_per|/kT,"
              "err_Ex_inf%,err_Ey_inf%,err_Ez_inf%,err_mag_inf%,"
              "err_Ex_per%,err_Ey_per%,err_Ez_per%,err_mag_per%")

    with open(csv_file, 'w') as f:
        f.write(header + '\n')
        for r in rows:
            ex_r, ey_r, ez_r = r['ex_r'], r['ey_r'], r['ez_r']
            ex_f, ey_f, ez_f = r['ex_f'], r['ey_f'], r['ez_f']
            ex_t, ey_t, ez_t = r['ex_t'], r['ey_t'], r['ez_t']
            mag_r = np.sqrt(ex_r**2 + ey_r**2 + ez_r**2)
            mag_f = np.sqrt(ex_f**2 + ey_f**2 + ez_f**2)
            mag_t = np.sqrt(ex_t**2 + ey_t**2 + ez_t**2)

            E_inf = dh_vec_inf(r['r_vec'])
            E_per = dh_vec_periodic(r['r_vec'])
            mag_inf = np.linalg.norm(E_inf)
            mag_per = np.linalg.norm(E_per)

            f.write(
                f"{px:.6f},{py:.6f},{pz:.6f},"
                f"{r['dist']:.8e},{r['ix']:d},{r['iy']:d},{r['iz']:d},"
                f"{r['x']:.6f},{r['y']:.6f},{r['z']:.6f},"
                f"{ex_r:.8e},{ey_r:.8e},{ez_r:.8e},{mag_r:.8e},"
                f"{ex_f:.8e},{ey_f:.8e},{ez_f:.8e},{mag_f:.8e},"
                f"{ex_t:.8e},{ey_t:.8e},{ez_t:.8e},{mag_t:.8e},"
                f"{E_inf[0]:.8e},{E_inf[1]:.8e},{E_inf[2]:.8e},{mag_inf:.8e},"
                f"{E_per[0]:.8e},{E_per[1]:.8e},{E_per[2]:.8e},{mag_per:.8e},"
                f"{rel_err(ex_t,E_inf[0]):.4f},{rel_err(ey_t,E_inf[1]):.4f},"
                f"{rel_err(ez_t,E_inf[2]):.4f},{rel_err(mag_t,mag_inf):.4f},"
                f"{rel_err(ex_t,E_per[0]):.4f},{rel_err(ey_t,E_per[1]):.4f},"
                f"{rel_err(ez_t,E_per[2]):.4f},{rel_err(mag_t,mag_per):.4f}\n"
            )

    print(f"CSV exportado: {csv_file} ({len(rows)} puntos)")


def plot_components_comparison(distances_real, E_real, distances_fourier, E_fourier,
                               distances_total, E_total, grid_size, rho_el, charge, epsilon, kt, kappa, alpha, rc,
                               rmin=0.1, output_file='ewald_components.png', max_points=10000,
                               log_scale=True, log_log=False, log_log_error=False, xmin=None, ymin=None,
                               show_ewald_components=True,
                               show_ewald_total=True, show_dh=True, show_ewald_theory_total=True,
                               show_dh_periodic=False, dh_periodic_shells=1,
                               direction=None, angle_tol_deg=15.0,
                               r_vecs_total=None, pointwise_error=False,
                               show_error_inf=True, show_error_dir=True,
                               particle_radius=0.0, only_error_dh_per_rad=False,
                               show_dh_periodic_rad=False,
                               gaussian_source=False, sigma=1.0,
                               gaussian_error=False, gaussian_error_periodic=False):
    """
    Genera gráfico comparando las componentes Real, Fourier, Total y Teórico

    Args:
        distances_real: Distancias para componente real
        E_real: Campo eléctrico de espacio real
        distances_fourier: Distancias para componente Fourier
        E_fourier: Campo eléctrico de espacio recíproco
        distances_total: Distancias para campo total
        E_total: Campo eléctrico total (real + fourier)
        charge: Carga de la partícula
        epsilon: Permitividad
        alpha: Parámetro de splitting de Ewald
        rc: Radio de corte
        output_file: Nombre del archivo de salida
        max_points: Máximo número de puntos a graficar
        log_scale: Si True, usa escala logarítmica
        show_ewald_components: Si True, grafica las componentes real y Fourier de Ewald por separado
        show_ewald_total: Si True, grafica la suma total de Ewald (real + Fourier)
        show_dh: Si True, grafica la curva teórica de Debye-Hückel
        show_ewald_theory_total: Si True, grafica la suma teórica total de Ewald (real + recíproco teóricos)
    """
    # Ordenar por distancia
    sort_idx_real = np.argsort(distances_real)
    distances_real_sorted = distances_real[sort_idx_real]
    E_real_sorted = E_real[sort_idx_real]

    sort_idx_fourier = np.argsort(distances_fourier)
    distances_fourier_sorted = distances_fourier[sort_idx_fourier]
    E_fourier_sorted = E_fourier[sort_idx_fourier]

    sort_idx_total = np.argsort(distances_total)
    distances_total_sorted = distances_total[sort_idx_total]
    E_total_sorted = E_total[sort_idx_total]
    r_vecs_total_sorted = r_vecs_total[sort_idx_total] if r_vecs_total is not None else None

    # Submuestrear si hay demasiados puntos
    def subsample(dist, vals, max_pts):
        if len(dist) > max_pts:
            step = len(dist) // max_pts
            return dist[::step], vals[::step]
        return dist, vals

    def subsample3(dist, vals, vecs, max_pts):
        if len(dist) > max_pts:
            step = len(dist) // max_pts
            return dist[::step], vals[::step], vecs[::step] if vecs is not None else None
        return dist, vals, vecs

    distances_real_sorted, E_real_sorted = subsample(distances_real_sorted, E_real_sorted, max_points)
    distances_fourier_sorted, E_fourier_sorted = subsample(distances_fourier_sorted, E_fourier_sorted, max_points)
    distances_total_sorted, E_total_sorted, r_vecs_total_sorted = subsample3(
        distances_total_sorted, E_total_sorted, r_vecs_total_sorted, max_points)

    # Filtrar datos para cortar en rmin (solo plotear r >= rmin)
    mask_real = distances_real_sorted >= rmin
    distances_real_sorted = distances_real_sorted[mask_real]
    E_real_sorted = E_real_sorted[mask_real]

    mask_fourier = distances_fourier_sorted >= rmin
    distances_fourier_sorted = distances_fourier_sorted[mask_fourier]
    E_fourier_sorted = E_fourier_sorted[mask_fourier]

    mask_total = distances_total_sorted >= rmin
    distances_total_sorted = distances_total_sorted[mask_total]
    E_total_sorted = E_total_sorted[mask_total]
    r_vecs_total_sorted = r_vecs_total_sorted[mask_total] if r_vecs_total_sorted is not None else None

    # Determinar rango de r para los datos y para la teoría por separado
    non_empty = [d for d in [distances_real_sorted, distances_fourier_sorted,
                              distances_total_sorted] if len(d) > 0]
    all_dists = np.concatenate(non_empty) if non_empty else np.array([])
    r_max_data = float(np.max(all_dists)) * 1.1 if len(all_dists) > 0 else 16.0
    # La curva teórica siempre se extiende hasta al menos r_max_data,
    # con un mínimo de 16 unidades para que sea informativa
    r_max_theory = max(r_max_data, 16.0)
    r_theory = np.linspace(rmin, r_max_theory, 500)

    # Calcular curvas teóricas
    L = float(grid_size[0])  # tamaño de caja (asume cubo)
    E_theory_dh = TheoreticalField.compute_theory_for_distances(r_theory, q=charge, epsilon=epsilon,
                                                                kt=kt, kappa=kappa)
    E_theory_real, E_theory_recip = TheoreticalField.compute_ewald_components(r_theory, q=charge,
                                                                               epsilon=epsilon, kt=kt,
                                                                               alpha=alpha, rc=rc)
    # Suma teórica Ewald = componente real teórica + componente recíproca teórica
    E_theory_ewald_total = E_theory_real + E_theory_recip

    # DH periódico
    E_theory_dh_periodic = None
    if show_dh_periodic:
        n_img = (2*dh_periodic_shells + 1)**3 - 1
        dir_label = _direction_label(direction)
        print(f"  Calculando DH periódico: {dh_periodic_shells} capa(s), {n_img} imágenes, dir={dir_label}...")
        E_theory_dh_periodic = TheoreticalField.compute_dh_periodic(
            r_theory, q=charge, epsilon=epsilon, kt=kt, kappa=kappa,
            L=L, n_shells=dh_periodic_shells, direction=direction)

    # CHANGE INIT - DHFiniteRadius - Curvas teóricas con radio finito
    E_theory_dh_rad = None
    if particle_radius > 0:
        E_theory_dh_rad = TheoreticalField.compute_dh_finite_radius(
            r_theory, q=charge, epsilon=epsilon, kt=kt, kappa=kappa, a=particle_radius)

    E_theory_dh_per_rad = None
    if particle_radius > 0 and (show_dh_periodic or show_dh_periodic_rad):
        n_img = (2 * dh_periodic_shells + 1) ** 3 - 1
        print(f"  Calculando DH periódico radio finito para panel 1 (a={particle_radius:.2f}): "
              f"{dh_periodic_shells} capa(s), {n_img} imágenes...")
        E_theory_dh_per_rad = TheoreticalField.compute_dh_periodic_finite_radius(
            r_theory, q=charge, epsilon=epsilon, kt=kt, kappa=kappa,
            a=particle_radius, L=L, n_shells=dh_periodic_shells, direction=direction)

    E_dh_per_rad_pts = None
    if only_error_dh_per_rad and particle_radius > 0 and show_dh_periodic_rad:
        n_img = (2 * dh_periodic_shells + 1) ** 3 - 1
        print(f"  Calculando DH periódico radio finito para panel 2 (a={particle_radius:.2f}): "
              f"{dh_periodic_shells} capa(s), {n_img} imágenes...")
        E_dh_per_rad_pts = TheoreticalField.compute_dh_periodic_finite_radius(
            distances_total_sorted, q=charge, epsilon=epsilon, kt=kt, kappa=kappa,
            a=particle_radius, L=L, n_shells=dh_periodic_shells, direction=direction)
    # CHANGE END - DHFiniteRadius

    # CHANGE INIT - GaussianDH - Cálculo de curvas gaussiana-DH para efield
    E_gauss_dh = None
    E_gauss_dh_per = None
    if gaussian_source:
        print(f"  Calculando |E| gaussiana-DH (σ={sigma:.3f}) en r_theory...")
        E_gauss_dh = TheoreticalField.compute_gaussian_dh(
            r_theory, q=charge, epsilon=epsilon, kt=kt, kappa=kappa, sigma=sigma)
        if show_dh_periodic and gaussian_error_periodic:
            n_img = (2 * dh_periodic_shells + 1) ** 3 - 1
            print(f"  Calculando |E| gaussiana-DH periódico ({dh_periodic_shells} capa(s), "
                  f"{n_img} imgs)...")
            E_gauss_dh_per = TheoreticalField.compute_gaussian_dh_periodic(
                r_theory, q=charge, epsilon=epsilon, kt=kt, kappa=kappa,
                sigma=sigma, L=L, n_shells=dh_periodic_shells, direction=direction)
    # CHANGE END - GaussianDH

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 12))

    # Panel 1: Componentes separadas
    if show_ewald_components:
        ax1.plot(distances_real_sorted, E_real_sorted, 'r.', markersize=6, alpha=0.8,
                 label='Ludwig: Espacio Real (erfc)')
        ax1.plot(distances_fourier_sorted, E_fourier_sorted, 'b.', markersize=6, alpha=0.8,
                 label='Ludwig: Espacio Recíproco (Fourier)')

    if show_ewald_total:
        ax1.plot(distances_total_sorted, E_total_sorted, 'g.', markersize=6, alpha=0.8,
                 label='Ludwig: Total (Real + Fourier)')

    # Teoría de componentes Ewald (solo si se muestran componentes)
    if show_ewald_components:
        ax1.plot(r_theory, E_theory_real, 'r-', linewidth=1.0, alpha=0.7,
                 label=f'Teoría: Real erfc(αr)/r² (α={alpha:.2f}, rc={rc:.1f})')
        ax1.plot(r_theory, E_theory_recip, 'b-', linewidth=1.0, alpha=0.7,
                 label=f'Teoría: Recíproco largo alcance [erf(αr)/r²-(2α/√π)e^(-α²r²)/r] (α={alpha:.2f})')

    # Suma teórica de Ewald (real + recíproco teóricos)
    if show_ewald_theory_total:
        ax1.plot(r_theory, E_theory_ewald_total, 'm-', linewidth=1.5, alpha=0.8,
                 label=f'Teoría: Ewald total (real+recíproco, α={alpha:.2f})')

    # Etiqueta para teoría total (Debye-Hückel o Coulomb)
    if kappa > 0:
        lambda_d = 1.0 / kappa
        theory_label = f'Teoría: Debye-Hückel (λ_D={lambda_d:.2f}, κ={kappa:.3f})'
    else:
        theory_label = 'Teoría: Coulomb puro 1/r²'
    if show_dh:
        ax1.plot(r_theory, E_theory_dh, 'k-', linewidth=1.0, alpha=0.8, label=theory_label)

    if show_dh_periodic and E_theory_dh_periodic is not None:
        n_img = (2*dh_periodic_shells + 1)**3 - 1
        dir_label = _direction_label(direction)
        if kappa > 0:
            periodic_label = (f'Teoría: DH periódico ({dh_periodic_shells} capa(s), '
                              f'{n_img} imgs, dir={dir_label})')
        else:
            periodic_label = (f'Teoría: Coulomb periódico ({dh_periodic_shells} capa(s), '
                              f'{n_img} imgs, L={L:.0f}, dir={dir_label})')
        ax1.plot(r_theory, E_theory_dh_periodic, 'r--', linewidth=1.5, alpha=0.9,
                 label=periodic_label)
    # CHANGE INIT - DHFiniteRadius - Curva radio finito en panel 1
    if E_theory_dh_rad is not None:
        ax1.plot(r_theory, E_theory_dh_rad, 'b--', linewidth=1.5, alpha=0.9,
                 label=f'Teoría: DH radio finito (a={particle_radius:.2f} lu)')
        ax1.axvline(x=particle_radius, color='blue', linestyle=':', linewidth=1.0,
                    alpha=0.6, label=f'superficie (a={particle_radius:.2f})')
    if E_theory_dh_per_rad is not None:
        n_img = (2*dh_periodic_shells + 1)**3 - 1
        dir_label = _direction_label(direction)
        ax1.plot(r_theory, E_theory_dh_per_rad, 'c--', linewidth=1.5, alpha=0.9,
                 label=(f'Teoría: DH periódico radio finito (a={particle_radius:.2f}, '
                        f'{dh_periodic_shells} capa(s), {n_img} imgs, dir={dir_label})'))
    # CHANGE END - DHFiniteRadius

    # CHANGE INIT - GaussianDH - Curvas gaussiana-DH en panel 1
    if gaussian_source and E_gauss_dh is not None:
        if kappa > 0:
            label_gauss = (f'Gauss-DH (σ={sigma:.2f}, λ_D={1.0/kappa:.2f}, κ={kappa:.3f})')
        else:
            label_gauss = f'Gauss-Coulomb (σ={sigma:.2f})'
        ax1.plot(r_theory, E_gauss_dh, 'c-', linewidth=2.0, alpha=0.9, label=label_gauss)
    if gaussian_source and E_gauss_dh_per is not None:
        n_img = (2 * dh_periodic_shells + 1) ** 3 - 1
        dir_label = _direction_label(direction)
        ax1.plot(r_theory, E_gauss_dh_per, 'c--', linewidth=1.5, alpha=0.9,
                 label=(f'Gauss-DH periódico (σ={sigma:.2f}, {dh_periodic_shells} cap., '
                        f'{n_img} imgs, dir={dir_label})'))
    # CHANGE END - GaussianDH

    # Marcar rc y λ_D
    ax1.axvline(x=rc, color='orange', linestyle='--', linewidth=2, alpha=0.7,
                label=f'Radio de corte rc={rc:.1f}')
    if kappa > 0:
        ax1.axvline(x=lambda_d, color='purple', linestyle=':', linewidth=2, alpha=0.7,
                    label=f'Longitud Debye λ_D={lambda_d:.2f}')

    ax1.set_xlabel('Distancia r [lattice units]', fontsize=font_size)
    ax1.set_ylabel('|E(r)|/kT [lattice units]', fontsize=font_size)

    title_str = f'|E(r)|/kT vs. distancia\nL=({grid_size[0]:d}×{grid_size[1]:d}×{grid_size[2]:d}), rho_el={rho_el:.2e}, q={charge:.2e}, ε={epsilon:.2e}, kT={kt:.2e}, α={alpha:.3f}, rc={rc:.1f}'
    if kappa > 0:
        title_str += f', κ={kappa:.3f}, λ_D={lambda_d:.2f}'
    if direction is not None:
        title_str += f'\nDirección datos: {_direction_label(direction)} (tol ±{angle_tol_deg:.0f}°)'
    ax1.set_title(title_str, fontsize=font_size, fontweight='bold')
    ax1.legend(fontsize=font_size, loc='best')
    ax1.grid(True, alpha=0.3)

    # Calcular ylim basado en los datos de Ludwig para que no queden aplastados
    # por la curva teórica que diverge cerca de r=0
    ludwig_vals = np.concatenate([v for v in [E_real_sorted, E_fourier_sorted, E_total_sorted]
                                  if len(v) > 0])
    if len(ludwig_vals) > 0:
        pos_ludwig = ludwig_vals[ludwig_vals > 0]
        y_top = float(np.max(ludwig_vals)) * 2.0
        y_bottom_auto = float(np.min(pos_ludwig)) * 0.5 if len(pos_ludwig) > 0 else None
    else:
        y_top = None
        y_bottom_auto = None

    # Configurar escala de ejes: X limitado al rango de datos, Y a los valores de Ludwig
    x_max = r_max_data
    if log_log:
        x_left = xmin if xmin is not None else max(rmin, 0.1)
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.set_xlim(x_left, x_max)
        ax1.set_ylim(bottom=ymin if ymin is not None else y_bottom_auto,
                     top=y_top)
    elif log_scale:
        x_left = xmin if xmin is not None else 0.0
        ax1.set_yscale('log')
        ax1.set_xlim(x_left, x_max)
        ax1.set_ylim(bottom=ymin if ymin is not None else y_bottom_auto,
                     top=y_top)
    else:
        x_left = xmin if xmin is not None else 0.0
        ax1.set_xlim(x_left, x_max)
        if ymin is not None:
            ax1.set_ylim(bottom=ymin, top=y_top)
        elif y_top is not None:
            ax1.set_ylim(top=y_top)

    # Panel 2: Diferencia entre total y teoría (solo si hay ambos datos)
    show_panel2 = show_ewald_total and (
        (show_dh and show_error_inf) or
        (show_dh_periodic and show_error_dir) or
        (show_dh_periodic_rad and show_error_dir) or
        pointwise_error or
        (only_error_dh_per_rad and E_dh_per_rad_pts is not None) or
        (gaussian_source and gaussian_error) or
        (gaussian_source and gaussian_error_periodic and show_dh_periodic))
    if show_panel2:
        # CHANGE INIT - DHFiniteRadius - Condicionar otras curvas con only_error_dh_per_rad
        if not only_error_dh_per_rad:
            if show_dh and show_error_inf:
                E_theory_at_ludwig = TheoreticalField.compute_theory_for_distances(
                    distances_total_sorted, q=charge, epsilon=epsilon, kt=kt, kappa=kappa)
                gap = E_total_sorted - E_theory_at_ludwig
                rel_gap = np.abs(gap) / (E_theory_at_ludwig + 1e-15) * 100
                ax2.plot(distances_total_sorted, rel_gap, 'purple', markersize=3.5, alpha=0.6,
                         marker='.', linestyle='none', label='|E_total - DH∞| / DH∞ × 100')

            if particle_radius > 0:
                E_rad_at_pts = TheoreticalField.compute_dh_finite_radius(
                    distances_total_sorted, q=charge, epsilon=epsilon, kt=kt,
                    kappa=kappa, a=particle_radius)
                valid_rad = ~np.isnan(E_rad_at_pts)
                rel_gap_rad = (np.abs(E_total_sorted[valid_rad] - E_rad_at_pts[valid_rad])
                               / (np.abs(E_rad_at_pts[valid_rad]) + 1e-15) * 100)
                ax2.plot(distances_total_sorted[valid_rad], rel_gap_rad, 'b',
                         markersize=3.5, alpha=0.7, marker='.', linestyle='none',
                         label=f'|E_total - DH a={particle_radius:.2f}| / DH_a × 100')

            if show_dh_periodic and show_error_dir and E_theory_dh_periodic is not None:
                E_dh_per_at_ludwig = TheoreticalField.compute_dh_periodic(
                    distances_total_sorted, q=charge, epsilon=epsilon, kt=kt, kappa=kappa,
                    L=L, n_shells=dh_periodic_shells, direction=direction)
                gap_per = E_total_sorted - E_dh_per_at_ludwig
                rel_gap_per = np.abs(gap_per) / (np.abs(E_dh_per_at_ludwig) + 1e-15) * 100
                ax2.plot(distances_total_sorted, rel_gap_per, 'darkorange', markersize=3.5, alpha=0.7,
                         marker='.', linestyle='none',
                         label=f'|E_total - DH periódico (dir)| / DH periódico × 100')

            if pointwise_error and r_vecs_total_sorted is not None and len(r_vecs_total_sorted) > 0:
                print(f"  Calculando error puntual en {len(r_vecs_total_sorted)} nodos...")
                E_dh_pointwise = TheoreticalField.compute_dh_periodic_at_points(
                    r_vecs_total_sorted, q=charge, epsilon=epsilon, kt=kt, kappa=kappa,
                    L=L, n_shells=dh_periodic_shells)
                gap_pw = E_total_sorted - E_dh_pointwise
                rel_gap_pw = np.abs(gap_pw) / (np.abs(E_dh_pointwise) + 1e-15) * 100
                ax2.plot(distances_total_sorted, rel_gap_pw, 'green', markersize=5, alpha=0.8,
                         marker='.', linestyle='none',
                         label=f'|E_total - DH periódico (punto exacto)| / DH × 100')

        if E_dh_per_rad_pts is not None:
            n_img = (2 * dh_periodic_shells + 1) ** 3 - 1
            valid_pr = ~np.isnan(E_dh_per_rad_pts)
            rel_gap_pr = (np.abs(E_total_sorted[valid_pr] - E_dh_per_rad_pts[valid_pr])
                          / (np.abs(E_dh_per_rad_pts[valid_pr]) + 1e-15) * 100)
            ax2.plot(distances_total_sorted[valid_pr], rel_gap_pr, 'darkgreen',
                     markersize=3.5, alpha=0.8, marker='.', linestyle='none',
                     label=(f'|E_total - DH per a={particle_radius:.2f}| / DH_per_a × 100'
                            f' ({dh_periodic_shells} capa(s), {n_img} imgs)'))
        # CHANGE END - DHFiniteRadius

        # CHANGE INIT - GaussianDH - Error vs gaussiana-DH en panel 2
        if gaussian_source and gaussian_error:
            E_gauss_at_pts = TheoreticalField.compute_gaussian_dh(
                distances_total_sorted, q=charge, epsilon=epsilon, kt=kt,
                kappa=kappa, sigma=sigma)
            rel_gap_gauss = (np.abs(E_total_sorted - E_gauss_at_pts)
                             / (np.abs(E_gauss_at_pts) + 1e-15) * 100)
            if kappa > 0:
                label_ge = f'|E_total - Gauss-DH (σ={sigma:.2f})| / Gauss-DH × 100'
            else:
                label_ge = f'|E_total - Gauss-Coulomb (σ={sigma:.2f})| / Gauss-Coulomb × 100'
            ax2.plot(distances_total_sorted, rel_gap_gauss, 'teal', markersize=3.5, alpha=0.8,
                     marker='.', linestyle='none', label=label_ge)

        if gaussian_source and gaussian_error_periodic and show_dh_periodic:
            print(f"  Calculando error vs gaussiana-DH periódica en {len(distances_total_sorted)} nodos...")
            if r_vecs_total_sorted is not None and len(r_vecs_total_sorted) > 0:
                E_gauss_per_pts = TheoreticalField.compute_gaussian_dh_periodic_at_points(
                    r_vecs_total_sorted, q=charge, epsilon=epsilon, kt=kt,
                    kappa=kappa, sigma=sigma, L=L, n_shells=dh_periodic_shells)
            else:
                E_gauss_per_pts = TheoreticalField.compute_gaussian_dh_periodic(
                    distances_total_sorted, q=charge, epsilon=epsilon, kt=kt,
                    kappa=kappa, sigma=sigma, L=L, n_shells=dh_periodic_shells,
                    direction=direction)
            rel_gap_gauss_per = (np.abs(E_total_sorted - E_gauss_per_pts)
                                 / (np.abs(E_gauss_per_pts) + 1e-15) * 100)
            n_img = (2 * dh_periodic_shells + 1) ** 3 - 1
            ax2.plot(distances_total_sorted, rel_gap_gauss_per, 'darkcyan', markersize=4, alpha=0.85,
                     marker='.', linestyle='none',
                     label=(f'|E_total - Gauss-DH per (σ={sigma:.2f})| / Gauss-DH_per × 100'
                            f' ({dh_periodic_shells} cap., {n_img} imgs)'))
        # CHANGE END - GaussianDH

        ax2.axhline(y=10, color='r', linestyle=':', alpha=0.5, linewidth=2, label='10% error')
        ax2.axhline(y=1, color='g', linestyle=':', alpha=0.5, linewidth=2, label='1% error')
        ax2.axvline(x=rc, color='orange', linestyle='--', linewidth=2, alpha=0.7,
                    label=f'rc={rc:.1f}')

        # problem_zone = (distances_total_sorted > 1.0) & (distances_total_sorted < 3.0)
        # if np.any(problem_zone):
        #     r_problem = distances_total_sorted[problem_zone]
        #     gap_problem = rel_gap[problem_zone]
        #     ax2.fill_between(r_problem, 0, gap_problem, color='red', alpha=0.2,
        #                     label='Zona problemática (r ≈ 1-3)')

        ax2.set_xlabel('Distancia r [lattice units]', fontsize=font_size)
        ax2.set_ylabel('Error relativo [%]', fontsize=font_size)
        ax2.set_title('Error relativo vs. distancia',
                      fontsize=font_size, fontweight='bold')
        ax2.legend(fontsize=font_size)
        ax2.grid(True, alpha=0.3)
        if log_log or log_log_error:
            x_left2 = xmin if xmin is not None else max(rmin, 0.1)
            ax2.set_xscale('log')
            ax2.set_xlim(x_left2, x_max)
        else:
            x_left2 = xmin if xmin is not None else 0.0
            ax2.set_xlim(x_left2, x_max)
        if log_log_error:
            ax2.set_yscale('log')
            if ymin is not None:
                ax2.set_ylim(bottom=ymin)
    else:
        ax2.set_visible(False)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nGráfico de componentes guardado en: {output_file}")
    plt.close()


def export_csv(efield_real_reader, efield_fourier_reader, efield_total_reader,
               particle_pos, grid_size, kt, rmin, rmax, csv_file,
               charge=1.0, epsilon=1.0, kappa=0.0, L=32.0, n_shells=1):
    """
    Exporta todos los puntos de la malla a CSV ordenados por distancia a la partícula.

    Columnas: ix, iy, iz, x, y, z, dist,
              Ex_real, Ey_real, Ez_real, |E_real|/kT,
              Ex_fourier, Ey_fourier, Ez_fourier, |E_fourier|/kT,
              Ex_total, Ey_total, Ez_total, |E_total|/kT,
              Ex_dh_inf, Ey_dh_inf, Ez_dh_inf, |E_dh_inf|,
              Ex_dh_per, Ey_dh_per, Ez_dh_per, |E_dh_per|,
              err_Ex_inf%, err_Ey_inf%, err_Ez_inf%, err_mag_inf%,
              err_Ex_per%, err_Ey_per%, err_Ez_per%, err_mag_per%
    """
    nx, ny, nz = grid_size
    # Sin límite de rmax: exportar todos los puntos de la malla
    max_possible = np.sqrt(nx**2 + ny**2 + nz**2)
    if rmax is None:
        rmax = max_possible

    px, py, pz = particle_pos
    prefactor = charge / (4.0 * np.pi * epsilon * kt)

    def dh_scalar(d):
        d = max(d, 1e-10)
        if kappa > 0:
            return prefactor * np.exp(-kappa * d) * (kappa / d + 1.0 / d**2)
        else:
            return prefactor / d**2

    def dh_vec_inf(r_vec):
        """Campo DH infinito: E = dh_scalar(|r|) * r/|r|"""
        dist = np.linalg.norm(r_vec)
        if dist < 1e-10:
            return np.zeros(3)
        return dh_scalar(dist) * r_vec / dist

    def dh_vec_periodic(r_vec):
        """Campo DH periódico: suma imagen central + shells imágenes periódicas"""
        E = dh_vec_inf(r_vec)
        for inx in range(-n_shells, n_shells + 1):
            for iny in range(-n_shells, n_shells + 1):
                for inz in range(-n_shells, n_shells + 1):
                    if inx == 0 and iny == 0 and inz == 0:
                        continue
                    r_img = np.array([inx * L, iny * L, inz * L])
                    dr = r_vec - r_img
                    dist_img = np.linalg.norm(dr)
                    if dist_img < 1e-10:
                        continue
                    E += dh_scalar(dist_img) * dr / dist_img
        return E

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
                if dist < rmin or dist > rmax:
                    continue

                ex_r = efield_real_reader.Ex[i, j, k] / kt if efield_real_reader else float('nan')
                ey_r = efield_real_reader.Ey[i, j, k] / kt if efield_real_reader else float('nan')
                ez_r = efield_real_reader.Ez[i, j, k] / kt if efield_real_reader else float('nan')
                mag_r = np.sqrt(ex_r**2 + ey_r**2 + ez_r**2)

                ex_f = efield_fourier_reader.Ex[i, j, k] / kt if efield_fourier_reader else float('nan')
                ey_f = efield_fourier_reader.Ey[i, j, k] / kt if efield_fourier_reader else float('nan')
                ez_f = efield_fourier_reader.Ez[i, j, k] / kt if efield_fourier_reader else float('nan')
                mag_f = np.sqrt(ex_f**2 + ey_f**2 + ez_f**2)

                ex_t = efield_total_reader.Ex[i, j, k] / kt if efield_total_reader else float('nan')
                ey_t = efield_total_reader.Ey[i, j, k] / kt if efield_total_reader else float('nan')
                ez_t = efield_total_reader.Ez[i, j, k] / kt if efield_total_reader else float('nan')
                mag_t = np.sqrt(ex_t**2 + ey_t**2 + ez_t**2)

                # Teoría DH en este punto exacto
                E_inf = dh_vec_inf(r_vec)
                E_per = dh_vec_periodic(r_vec)
                mag_inf = np.linalg.norm(E_inf)
                mag_per = np.linalg.norm(E_per)

                # Errores relativos (%) componente a componente y módulo
                err_ex_inf = rel_err(ex_t, E_inf[0])
                err_ey_inf = rel_err(ey_t, E_inf[1])
                err_ez_inf = rel_err(ez_t, E_inf[2])
                err_mag_inf = rel_err(mag_t, mag_inf)

                err_ex_per = rel_err(ex_t, E_per[0])
                err_ey_per = rel_err(ey_t, E_per[1])
                err_ez_per = rel_err(ez_t, E_per[2])
                err_mag_per = rel_err(mag_t, mag_per)

                rows.append((dist, i, j, k, xi, yj, zk,
                             ex_r, ey_r, ez_r, mag_r,
                             ex_f, ey_f, ez_f, mag_f,
                             ex_t, ey_t, ez_t, mag_t,
                             E_inf[0], E_inf[1], E_inf[2], mag_inf,
                             E_per[0], E_per[1], E_per[2], mag_per,
                             err_ex_inf, err_ey_inf, err_ez_inf, err_mag_inf,
                             err_ex_per, err_ey_per, err_ez_per, err_mag_per))

    rows.sort(key=lambda r: r[0])

    header = ("px,py,pz,"
              "dist,ix,iy,iz,x,y,z,"
              "Ex_real/kT,Ey_real/kT,Ez_real/kT,|E_real|/kT,"
              "Ex_fourier/kT,Ey_fourier/kT,Ez_fourier/kT,|E_fourier|/kT,"
              "Ex_total/kT,Ey_total/kT,Ez_total/kT,|E_total|/kT,"
              "Ex_dh_inf/kT,Ey_dh_inf/kT,Ez_dh_inf/kT,|E_dh_inf|/kT,"
              "Ex_dh_per/kT,Ey_dh_per/kT,Ez_dh_per/kT,|E_dh_per|/kT,"
              "err_Ex_inf%,err_Ey_inf%,err_Ez_inf%,err_mag_inf%,"
              "err_Ex_per%,err_Ey_per%,err_Ez_per%,err_mag_per%")

    with open(csv_file, 'w') as f:
        f.write(header + '\n')
        for r in rows:
            f.write(f"{px:.6f},{py:.6f},{pz:.6f},"
                    f"{r[0]:.8e},{r[1]:d},{r[2]:d},{r[3]:d},"
                    f"{r[4]:.6f},{r[5]:.6f},{r[6]:.6f},"
                    f"{r[7]:.8e},{r[8]:.8e},{r[9]:.8e},{r[10]:.8e},"
                    f"{r[11]:.8e},{r[12]:.8e},{r[13]:.8e},{r[14]:.8e},"
                    f"{r[15]:.8e},{r[16]:.8e},{r[17]:.8e},{r[18]:.8e},"
                    f"{r[19]:.8e},{r[20]:.8e},{r[21]:.8e},{r[22]:.8e},"
                    f"{r[23]:.8e},{r[24]:.8e},{r[25]:.8e},{r[26]:.8e},"
                    f"{r[27]:.4f},{r[28]:.4f},{r[29]:.4f},{r[30]:.4f},"
                    f"{r[31]:.4f},{r[32]:.4f},{r[33]:.4f},{r[34]:.4f}\n")

    print(f"CSV exportado: {csv_file} ({len(rows)} puntos)")


def main():
    parser = argparse.ArgumentParser(
        description='Plotea por separado las componentes Real y Recíproca de Ewald',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Ejemplos:
    # Mostrar todo (comportamiento por defecto si no se especifica ningún --show-*)
    python plot_componentes_ewald.py -n 1 -L 32 -q 1.0 -eps 10000 --alpha 0.5 --rc 5.0

    # Solo la suma total de Ewald y la curva de Debye-Hückel
    python plot_componentes_ewald.py -n 1 -L 32 -q 1.0 -eps 10000 --show-ewald-total --show-dh

    # Solo las componentes real y Fourier de Ewald (sin DH)
    python plot_componentes_ewald.py -n 1 -L 32 -q 1.0 -eps 10000 --show-ewald-components

    # Solo la curva teórica de Debye-Hückel
    python plot_componentes_ewald.py -n 1 -L 32 -q 1.0 -eps 10000 --show-dh

    # Con directorio personalizado
    python plot_componentes_ewald.py -n 1000 -d results/ -L 32 -q 1.0 -eps 10000

Flags de visualización (si no se especifica ninguno, se muestran todos):
    --show-ewald-components      Componentes real y Fourier de Ewald por separado
    --show-ewald-total           Suma total numérica de Ewald (real + Fourier de Ludwig)
    --show-ewald-theory-total    Suma teórica total de Ewald (real + recíproco teóricos)
    --show-dh                    Curva teórica de Debye-Hückel (sistema infinito)
    --show-dh-periodic           DH con imágenes periódicas (primeros vecinos a distancia L)
    --dh-periodic-shells N       Número de capas de imágenes (default: 1 = 26 imágenes a ~L)

NOTA: Este script requiere que ewald_charge.c genere archivos separados
      efield_real y efield_fourier (además del efield total).
      Los archivos deben seguir el formato:
        efield-NNNNNNNN.001-001
        efield_real-NNNNNNNN.001-001
        efield_fourier-NNNNNNNN.001-001
        colloids-NNNNNNNN.csv
        """
    )

    parser.add_argument('-n', '--nstep', type=int, required=True,
                       help='Número de paso temporal (ej: 1, 1000, etc)')
    parser.add_argument('-d', '--directory', default='.',
                       help='Directorio donde están los archivos (default: directorio actual)')
    parser.add_argument('-L', '--grid-size', type=int, required=True,
                       help='Tamaño de la malla (asume cubo LxLxL)')
    parser.add_argument('-q', '--charge', type=float, required=True,
                       help='Carga de la partícula')
    parser.add_argument('-eps', '--epsilon', type=float, required=True,
                       help='Permitividad relativa ε_r')
    parser.add_argument('-kt', '--kt', type=float, default=1.0,
                       help='Energía térmica kT (default: 1.0)')
    parser.add_argument('--kappa', type=float, default=None,
                       help='Parámetro de Debye κ = 1/λ_D (si no se especifica, se calcula de rho_el)')
    parser.add_argument('--rho-el', type=float, default=None,
                       help='Densidad electrolítica ρ_el para calcular κ (si no se da kappa)')
    parser.add_argument('--alpha', type=float, default=0.5,
                       help='Parámetro de splitting de Ewald α (default: 0.5)')
    parser.add_argument('--rc', type=float, default=5.0,
                       help='Radio de corte rc (default: 5.0)')
    parser.add_argument('-o', '--output', default='ewald_components.png',
                       help='Archivo de salida (default: ewald_components.png)')
    parser.add_argument('--rmin', type=float, default=0.1,
                       help='Radio mínimo (default: 0.1)')
    parser.add_argument('--rmax', type=float, default=None,
                       help='Radio máximo (default: L/2)')
    parser.add_argument('--max-points', type=int, default=10000,
                       help='Máximo número de puntos (default: 10000)')
    parser.add_argument('--linear', action='store_true',
                       help='Usar escala lineal (default: logarítmica)')
    parser.add_argument('--log-log', action='store_true',
                       help='Usar doble eje logarítmico log-log (default: semi-log)')
    parser.add_argument('--log-log-error', action='store_true',
                       help='Usar doble eje logarítmico en el panel de error relativo')
    parser.add_argument('--xmin', type=float, default=None,
                       help='Valor mínimo del eje X en ambos paneles (default: 0 o rmin en log-log)')
    parser.add_argument('--ymin', type=float, default=None,
                       help='Valor mínimo del eje Y (aplica a ax1 en log/semilog, y a ax2 en log-log-error)')
    parser.add_argument('--show-ewald-components', action='store_true', default=None,
                       help='Mostrar componentes real y Fourier de Ewald por separado')
    parser.add_argument('--show-ewald-total', action='store_true', default=None,
                       help='Mostrar la suma total de Ewald numérica (real + Fourier)')
    parser.add_argument('--show-ewald-theory-total', action='store_true', default=None,
                       help='Mostrar la suma teórica total de Ewald (real + recíproco teóricos)')
    parser.add_argument('--show-dh', action='store_true', default=None,
                       help='Mostrar la curva teórica de Debye-Hückel')
    parser.add_argument('--show-dh-periodic', action='store_true', default=None,
                       help='Mostrar Debye-Hückel con imágenes periódicas')
    parser.add_argument('--dh-periodic-shells', type=int, default=1,
                       help='Número de capas de imágenes periódicas (default: 1 = 26 vecinos a distancia ~L)')
    parser.add_argument('--direction', default=None,
                       help=('Dirección para filtrar datos y evaluar DH periódico. '
                             "Ejemplos: 'x', 'y', 'z', 'xy', 'xz', 'yz', 'xyz', '1,0,0', '1,1,0'. "
                             'Sin este argumento se usan todos los puntos (isotrópico).'))
    parser.add_argument('--angle-tol', type=float, default=15.0,
                       help='Tolerancia angular en grados para filtrar puntos por dirección (default: 15)')
    parser.add_argument('--pointwise-error', action='store_true', default=False,
                       help=('Calcula el error de cada nodo contra el DH periódico evaluado '
                             'en el mismo punto 3D exacto (no en r·d̂). '
                             'Requiere --show-dh-periodic y --direction.'))
    parser.add_argument('--no-error-inf', action='store_true', default=False,
                       help='Ocultar la serie de error vs DH infinito (púrpura) en el panel 2')
    parser.add_argument('--no-error-dir', action='store_true', default=False,
                       help='Ocultar la serie de error vs DH periódico direccional (naranja) en el panel 2')
    # CHANGE INIT - DHFiniteRadius - Radio de partícula y solo-error periódico con radio
    parser.add_argument('--particle-radius', type=float, default=0.0,
                       help='Radio de la partícula a (lu). Activa curva DH radio finito '
                            'en panel 1 y error vs DH(a) en panel 2 (default: 0 = desactivado)')
    parser.add_argument('--only-error-dh-per-rad', action='store_true', default=False,
                       help='Mostrar en el panel de error SOLO la curva de error relativo '
                            'vs DH periódico con partícula de radio a. Requiere '
                            '--particle-radius > 0 y --show-dh-periodic-rad. '
                            'Suprime las demás curvas de error.')
    parser.add_argument('--show-dh-periodic-rad', action='store_true', default=None,
                       help='Mostrar DH periódico con radio finito en panel 1, independientemente '
                            'de --show-dh-periodic. Requiere --particle-radius > 0.')
    # CHANGE END - DHFiniteRadius

    # CHANGE INIT - GaussianDH - Argumentos para fuente de carga gaussiana
    parser.add_argument('--gaussian-source', action='store_true', default=False,
                       help='Usar distribución gaussiana ρ(r)=Q/(π^(3/2)σ³)exp(-r²/σ²) '
                            'como fuente en la teoría DH. Activa curva |E| gaussiana-DH.')
    parser.add_argument('--sigma', type=float, default=1.0,
                       help='Anchura σ de la distribución gaussiana de carga (lu). '
                            'Solo se usa con --gaussian-source (default: 1.0)')
    parser.add_argument('--gaussian-error', action='store_true', default=False,
                       help='Mostrar el error relativo vs la teoría gaussiana-DH (infinita) '
                            'en el panel 2. Requiere --gaussian-source.')
    parser.add_argument('--gaussian-error-periodic', action='store_true', default=False,
                       help='Mostrar el error relativo vs la teoría gaussiana-DH periódica '
                            'en el panel 2. Requiere --gaussian-source y --show-dh-periodic.')
    # CHANGE END - GaussianDH

    parser.add_argument('--csv', default=None,
                       help='Exportar datos a CSV (nombre base sin extensión; se añade -n_NNNNNN.csv)')
    parser.add_argument('--csv-points', default=None,
                       help='Exportar a CSV los mismos puntos que aparecen en el gráfico '
                            '(con filtro de dirección), incluyendo error vs DH∞ y DH periódico. '
                            'Nombre base sin extensión; se añade -n_NNNNNN.csv')

    args = parser.parse_args()

    # Parsear dirección
    try:
        direction = parse_direction(args.direction)
    except ValueError as e:
        print(f"ERROR: {e}")
        sys.exit(1)
    if direction is not None:
        print(f"Dirección seleccionada: {_direction_label(direction)} = {direction}, "
              f"tolerancia ±{args.angle_tol:.1f}°")

    # Si no se especifica ningún flag de visualización, mostrar todo por defecto
    if (args.show_ewald_components is None and args.show_ewald_total is None
            and args.show_ewald_theory_total is None and args.show_dh is None
            and args.show_dh_periodic is None):
        args.show_ewald_components = True
        args.show_ewald_total = True
        args.show_ewald_theory_total = True
        args.show_dh = True
        args.show_dh_periodic = False  # No por defecto: requiere L explícito
    else:
        args.show_ewald_components = bool(args.show_ewald_components)
        args.show_ewald_total = bool(args.show_ewald_total)
        args.show_ewald_theory_total = bool(args.show_ewald_theory_total)
        args.show_dh = bool(args.show_dh)
        args.show_dh_periodic = bool(args.show_dh_periodic)
    args.show_dh_periodic_rad = bool(args.show_dh_periodic_rad)

    grid_size = (args.grid_size, args.grid_size, args.grid_size)

    # Calcular κ (parámetro de Debye) si no se proporciona
    if args.kappa is None:
        if args.rho_el is not None:
            # Calcular κ a partir de la densidad electrolítica total
            # Para electroneutralidad: rho_total = rho_el *2 + |q|/V
            # donde rho_el es la densidad de iones del electrolito
            # y |q|/V es la densidad de contraiones para neutralizar la partícula
            volume = args.grid_size**3
            rho_counterions = abs(args.charge) / volume
            rho_total = args.rho_el * 2 + rho_counterions

            # κ² = 4πλ_B × ρ_total, donde λ_B = e²/(4πε kT)
            # λ_B es propiedad del solvente: se calcula con carga del ion = 1 (monovalente),
            # no con la carga del coloide (que puede ser arbitraria).
            lb = 1.0 / (4.0 * np.pi * args.epsilon * args.kt)  # Longitud de Bjerrum (ion unitario)
            kappa = np.sqrt(4.0 * np.pi * lb * rho_total)
            lambda_d = 1.0 / kappa if kappa > 0 else np.inf
            print(f"\nParámetros de Debye-Hückel calculados:")
            print(f"  Longitud de Bjerrum λ_B = {lb:.6e}")
            print(f"  Volumen del sistema V = {volume:.6e}")
            print(f"  Densidad electrolítica (input) ρ_el = {args.rho_el:.6e}")
            print(f"  Densidad de contraiones |q|/V = {rho_counterions:.6e}")
            print(f"  Densidad total (electroneutra) ρ_total = {rho_total:.6e}")
            print(f"  Parámetro de Debye κ = {kappa:.6e}")
            print(f"  Longitud de Debye λ_D = {lambda_d:.6e}")
        else:
            # Sin screening (Coulomb puro)
            kappa = 0.0
            print(f"\nUsando modelo de Coulomb puro (sin screening, κ=0)")
    else:
        kappa = args.kappa
        lambda_d = 1.0 / kappa if kappa > 0 else np.inf
        print(f"\nUsando κ = {kappa:.6e} (λ_D = {lambda_d:.6e})")

    # Construir nombres de archivos
    nstep_str = f"{args.nstep:09d}"
    base_dir = args.directory
    efield_real_file = os.path.join(base_dir, f"efield_real-{nstep_str}.001-001")
    efield_fourier_file = os.path.join(base_dir, f"efield_fourier-{nstep_str}.001-001")
    efield_total_file = os.path.join(base_dir, f"efield-{nstep_str}.001-001")
    nstep_str = f"{args.nstep:08d}"
    colloid_file = os.path.join(base_dir, f"config.cds{nstep_str}.001-001")

    # Determinar qué archivos son necesarios
    need_components = args.show_ewald_components or (args.csv is not None) or (args.csv_points is not None)
    need_total = (args.show_ewald_total or args.pointwise_error or args.only_error_dh_per_rad
                  or (not args.no_error_inf) or (not args.no_error_dir)
                  or (args.csv is not None) or (args.csv_points is not None))

    files_to_check = [colloid_file]
    if need_components:
        files_to_check += [efield_real_file, efield_fourier_file]
    if need_total:
        files_to_check.append(efield_total_file)

    for fname in files_to_check:
        if not os.path.exists(fname):
            print(f"ERROR: No se encontró el archivo: {fname}")
            sys.exit(1)

    # Leer posición de la partícula
    print(f"Leyendo posición de la partícula desde: {colloid_file}")
    particle_pos = ColloidReader.read_colloid_position_cds(colloid_file)
    print(f"  Posición: {particle_pos}")

    # Leer componente real
    efield_real_reader = None
    E_mag_real = None
    if args.show_ewald_components or args.csv is not None or args.csv_points is not None:
        print(f"\nLeyendo componente de espacio real desde: {efield_real_file}")
        efield_real_reader = EFieldReader(efield_real_file, grid_size)
        try:
            efield_real_reader.read()
            E_mag_real = efield_real_reader.compute_magnitude() / args.kt
            print(f"  Rango de |E_real|/kT: [{np.min(E_mag_real):.6e}, {np.max(E_mag_real):.6e}]")
        except Exception as e:
            print(f"ERROR: {e}")
            sys.exit(1)

    # Leer componente Fourier
    efield_fourier_reader = None
    E_mag_fourier = None
    if args.show_ewald_components or args.csv is not None or args.csv_points is not None:
        print(f"\nLeyendo componente de espacio recíproco desde: {efield_fourier_file}")
        efield_fourier_reader = EFieldReader(efield_fourier_file, grid_size)
        try:
            efield_fourier_reader.read()
            E_mag_fourier = efield_fourier_reader.compute_magnitude() / args.kt
            print(f"  Rango de |E_fourier|/kT: [{np.min(E_mag_fourier):.6e}, {np.max(E_mag_fourier):.6e}]")
        except Exception as e:
            print(f"ERROR: {e}")
            sys.exit(1)

    # Leer campo total
    efield_total_reader = None
    E_mag_total = None
    if need_total:
        print(f"\nLeyendo campo total desde: {efield_total_file}")
        efield_total_reader = EFieldReader(efield_total_file, grid_size)
        try:
            efield_total_reader.read()
            E_mag_total = efield_total_reader.compute_magnitude() / args.kt
            print(f"  Rango de |E_total|/kT: [{np.min(E_mag_total):.6e}, {np.max(E_mag_total):.6e}]")
        except Exception as e:
            print(f"ERROR: {e}")
            sys.exit(1)

    # Extraer puntos (arrays vacíos si no se necesitan)
    print(f"\nExtrayendo puntos (rmin={args.rmin}, rmax={args.rmax})...")
    empty = np.array([])
    distances_real, E_real = (extract_all_points(E_mag_real, particle_pos, grid_size,
                                                 min_radius=args.rmin, max_radius=args.rmax,
                                                 direction=direction, angle_tol_deg=args.angle_tol)
                              if E_mag_real is not None else (empty, empty))
    distances_fourier, E_fourier = (extract_all_points(E_mag_fourier, particle_pos, grid_size,
                                                       min_radius=args.rmin, max_radius=args.rmax,
                                                       direction=direction, angle_tol_deg=args.angle_tol)
                                    if E_mag_fourier is not None else (empty, empty))
    need_rvecs = (args.pointwise_error or (args.gaussian_source and args.gaussian_error_periodic)) and E_mag_total is not None
    if need_rvecs:
        distances_total, E_total, r_vecs_total = extract_all_points(
            E_mag_total, particle_pos, grid_size,
            min_radius=args.rmin, max_radius=args.rmax,
            direction=direction, angle_tol_deg=args.angle_tol,
            return_rvecs=True)
    else:
        distances_total, E_total = (extract_all_points(E_mag_total, particle_pos, grid_size,
                                                       min_radius=args.rmin, max_radius=args.rmax,
                                                       direction=direction, angle_tol_deg=args.angle_tol)
                                    if E_mag_total is not None else (empty, empty))
        r_vecs_total = None

    print(f"  Puntos extraídos - Real: {len(distances_real)}, Fourier: {len(distances_fourier)}, Total: {len(distances_total)}")

    # Construir nombre de archivo de salida con número de paso
    output_base, output_ext = os.path.splitext(args.output)
    if not output_ext:
        output_ext = '.png'
    output_file = f"{output_base}-n_{args.nstep:09d}{output_ext}"

    # Generar gráfico
    log_scale = not args.linear and not args.log_log
    print(f"\nGenerando gráfico de componentes...")
    plot_components_comparison(distances_real, E_real, distances_fourier, E_fourier,
                              distances_total, E_total, grid_size, args.rho_el, args.charge, args.epsilon,
                              args.kt, kappa, args.alpha, args.rc, args.rmin,
                              output_file, args.max_points, log_scale, args.log_log, args.log_log_error,
                              args.xmin, args.ymin, args.show_ewald_components, args.show_ewald_total,
                              args.show_dh, args.show_ewald_theory_total,
                              args.show_dh_periodic, args.dh_periodic_shells,
                              direction, args.angle_tol,
                              r_vecs_total, args.pointwise_error,
                              show_error_inf=not args.no_error_inf,
                              show_error_dir=not args.no_error_dir,
                              particle_radius=args.particle_radius,
                              only_error_dh_per_rad=args.only_error_dh_per_rad,
                              show_dh_periodic_rad=args.show_dh_periodic_rad,
                              gaussian_source=args.gaussian_source,
                              sigma=args.sigma,
                              gaussian_error=args.gaussian_error,
                              gaussian_error_periodic=args.gaussian_error_periodic)

    # Exportar CSV si se solicitó
    if args.csv is not None:
        csv_base, csv_ext = os.path.splitext(args.csv)
        if not csv_ext:
            csv_ext = '.csv'
        csv_file = f"{csv_base}-n_{args.nstep:09d}{csv_ext}"
        print(f"\nExportando datos a CSV...")
        export_csv(efield_real_reader, efield_fourier_reader, efield_total_reader,
                   particle_pos, grid_size, args.kt, args.rmin, args.rmax, csv_file,
                   charge=args.charge, epsilon=args.epsilon, kappa=kappa,
                   L=float(grid_size[0]), n_shells=args.dh_periodic_shells)

    # Exportar CSV con los mismos puntos del gráfico (filtro de dirección incluido)
    if args.csv_points is not None:
        csv_base, csv_ext = os.path.splitext(args.csv_points)
        if not csv_ext:
            csv_ext = '.csv'
        csv_file_pts = f"{csv_base}-n_{args.nstep:09d}{csv_ext}"
        print(f"\nExtrayendo puntos para CSV (con filtro de dirección)...")
        rows_full = extract_all_points_full(
            efield_real_reader, efield_fourier_reader, efield_total_reader,
            particle_pos, grid_size, args.kt,
            min_radius=args.rmin, max_radius=args.rmax,
            direction=direction, angle_tol_deg=args.angle_tol)
        print(f"  {len(rows_full)} puntos extraídos")
        print(f"  Calculando errores y exportando a {csv_file_pts}...")
        export_csv_from_points(rows_full, particle_pos,
                               charge=args.charge, epsilon=args.epsilon, kt=args.kt, kappa=kappa,
                               L=float(grid_size[0]), n_shells=args.dh_periodic_shells,
                               csv_file=csv_file_pts)

    print("\n¡Análisis de componentes completado!")


if __name__ == '__main__':
    main()
