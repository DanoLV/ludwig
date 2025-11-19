#!/usr/bin/env python3
"""
Script para verificar la simetría del campo eléctrico con respecto a la posición de la partícula.

Este script:
1. Lee archivos psi de Ludwig
2. Extrae la posición de la partícula del archivo colloids
3. Calcula el campo eléctrico E = -∇ψ
4. Verifica la simetría radial y en diferentes planos
5. Compara con modelo teórico de Poisson-Boltzmann (Debye-Hückel)
6. Genera gráficos comparativos y métricas cuantitativas

Autor: Script para análisis de simulaciones Ludwig
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
import os
from pathlib import Path
import re


class PsiReader:
    """Lee y procesa archivos psi de Ludwig"""

    def __init__(self, filename, grid_size):
        """
        Inicializa el lector de archivos psi

        Args:
            filename: Ruta al archivo psi
            grid_size: Tupla (nx, ny, nz) con el tamaño de la malla
        """
        self.filename = filename
        self.nx, self.ny, self.nz = grid_size
        self.psi = None

    def read(self):
        """Lee el archivo psi y lo organiza en una malla 3D"""
        data = np.loadtxt(self.filename)

        if len(data) != self.nx * self.ny * self.nz:
            raise ValueError(f"El archivo tiene {len(data)} valores, "
                           f"pero se esperaban {self.nx * self.ny * self.nz}")

        # Ludwig escribe en orden z-major: para cada x, para cada y, todos los z
        self.psi = data.reshape((self.nx, self.ny, self.nz))
        return self.psi

    def compute_electric_field(self):
        """
        Calcula el campo eléctrico E = -∇ψ usando diferencias finitas centradas

        Returns:
            Ex, Ey, Ez: Componentes del campo eléctrico
        """
        if self.psi is None:
            raise ValueError("Primero debe leer el archivo con read()")

        # Usar gradiente de numpy (diferencias finitas centradas en el interior)
        grad_psi = np.gradient(self.psi)

        # E = -∇ψ
        Ex = -grad_psi[0]
        Ey = -grad_psi[1]
        Ez = -grad_psi[2]

        return Ex, Ey, Ez


class ColloidReader:
    """Lee archivos de coloides de Ludwig"""

    @staticmethod
    def read_colloid_position(filename):
        """
        Lee la posición del coloide desde un archivo config.cds o colloids-*.csv

        Args:
            filename: Ruta al archivo (config.cds*.001-001 o colloids-*.csv)

        Returns:
            (x, y, z): Posición del coloide
        """
        # Determinar el formato del archivo
        if 'config.cds' in filename:
            # Formato config.cds de Ludwig (ASCII)
            # La posición está típicamente en la línea 40
            with open(filename, 'r') as f:
                lines = f.readlines()
                # Buscar línea con 3 valores flotantes (posición x, y, z)
                for i, line in enumerate(lines):
                    parts = line.strip().split()
                    if len(parts) == 3:
                        try:
                            x, y, z = float(parts[0]), float(parts[1]), float(parts[2])
                            # Verificar que sean valores razonables (no flags/enteros pequeños)
                            if abs(x) > 1.0 or abs(y) > 1.0 or abs(z) > 1.0:
                                return (x, y, z)
                        except ValueError:
                            continue
            raise ValueError(f"No se pudo leer posición del coloide en {filename}")

        else:
            # Formato colloids-*.csv
            with open(filename, 'r') as f:
                # Saltar líneas de comentario y encabezados
                for line in f:
                    if line.startswith('#'):
                        continue
                    if 'id' in line.lower() or 'index' in line.lower():
                        continue

                    # Primera línea de datos
                    parts = line.strip().split()
                    if len(parts) >= 3:
                        # Las columnas típicamente son: index x y z ...
                        x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                        return (x, y, z)
            raise ValueError(f"No se pudo leer posición del coloide en {filename}")


class TheoreticalField:
    """Calcula el campo eléctrico teórico según modelos analíticos"""

    @staticmethod
    def debye_huckel_field(charge_pos, grid_size, q=1.0, epsilon=1.0, kt=1.0,
                          kappa=0.0, scale_factor=None):
        """
        Calcula el campo eléctrico teórico usando el modelo de Debye-Hückel
        (solución de Poisson-Boltzmann linearizada)

        Args:
            charge_pos: Posición de la carga (x, y, z)
            grid_size: Tupla (nx, ny, nz)
            q: Carga (default: 1.0)
            epsilon: Permitividad relativa (default: 1.0)
            kt: Energía térmica k_B*T (default: 1.0)
            kappa: Parámetro de Debye κ = 1/λ_D (default: 0.0 = Coulomb)
            scale_factor: Factor de escala manual. Si None, usa q/(4πε*kt)

        Potenciales:
            - kappa = 0 (sin iones): φ = q/(4πε*kt·r) → E = q/(4πε*kt) · r/r³
            - kappa > 0 (con iones): φ = q/(4πε*kt·r) · exp(-κr) →
                                      E = q/(4πε*kt) · exp(-κr) · (1/r² + κ/r) · r̂

        Returns:
            Ex, Ey, Ez: Componentes del campo eléctrico teórico (arrays 3D)
        """
        nx, ny, nz = grid_size
        xc, yc, zc = charge_pos

        # Crear malla de coordenadas centradas en nodos
        x = np.arange(0.5, nx)
        y = np.arange(0.5, ny)
        z = np.arange(0.5, nz)

        X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

        # Vector distancia desde la carga
        dx = X - xc
        dy = Y - yc
        dz = Z - zc

        # Distancia radial
        r_squared = dx**2 + dy**2 + dz**2
        r = np.sqrt(r_squared)

        # Evitar división por cero en la posición de la carga
        r = np.where(r < 0.01, 0.01, r)
        r_squared = r**2

        # Prefactor
        if scale_factor is None:
            prefactor = q / (4.0 * np.pi * epsilon * kt)
        else:
            prefactor = scale_factor

        # Campo eléctrico según el modelo
        if kappa > 0:
            # Debye-Hückel: medio con iones
            # E = (q/4πε) · exp(-κr) · (κ/r + 1/r²) · r̂
            exp_factor = np.exp(-kappa * r)
            magnitude = prefactor * exp_factor * (kappa / r + 1.0 / r_squared)

            Ex = magnitude * dx / r
            Ey = magnitude * dy / r
            Ez = magnitude * dz / r
        else:
            # Coulomb: vacío sin screening
            # E = (q/4πε) · r/r³
            r_cubed = r_squared * r

            Ex = prefactor * dx / r_cubed
            Ey = prefactor * dy / r_cubed
            Ez = prefactor * dz / r_cubed

        return Ex, Ey, Ez

    @staticmethod
    def debye_huckel_potential(charge_pos, grid_size, q=1.0, epsilon=1.0, kt=1.0,
                              kappa=0.0, scale_factor=None):
        """
        Calcula el potencial eléctrico teórico usando el modelo de Debye-Hückel

        Args:
            charge_pos: Posición de la carga (x, y, z)
            grid_size: Tupla (nx, ny, nz)
            q: Carga (default: 1.0)
            epsilon: Permitividad relativa (default: 1.0)
            kt: Energía térmica k_B*T (default: 1.0)
            kappa: Parámetro de Debye κ = 1/λ_D (default: 0.0 = Coulomb)
            scale_factor: Factor de escala manual. Si None, usa q/(4πε*kt)

        Potenciales:
            - kappa = 0 (sin iones): φ = q/(4πε*kt·r)
            - kappa > 0 (con iones): φ = q/(4πε*kt·r) · exp(-κr)

        Returns:
            psi: Potencial eléctrico teórico (array 3D)
        """
        nx, ny, nz = grid_size
        xc, yc, zc = charge_pos

        # Crear malla de coordenadas centradas en nodos
        x = np.arange(0.5, nx)
        y = np.arange(0.5, ny)
        z = np.arange(0.5, nz)

        X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

        # Vector distancia desde la carga
        dx = X - xc
        dy = Y - yc
        dz = Z - zc

        # Distancia radial
        r = np.sqrt(dx**2 + dy**2 + dz**2)

        # Evitar división por cero en la posición de la carga
        r = np.where(r < 0.01, 0.01, r)

        # Prefactor
        if scale_factor is None:
            prefactor = q / (4.0 * np.pi * epsilon * kt)
        else:
            prefactor = scale_factor

        # Potencial según el modelo
        if kappa > 0:
            # Debye-Hückel: medio con iones
            # φ = (q/4πε*kt) · exp(-κr) / r
            psi = prefactor * np.exp(-kappa * r) / r
        else:
            # Coulomb: vacío sin screening
            # φ = (q/4πε*kt) / r
            psi = prefactor / r

        return psi

    @staticmethod
    def fit_scale_factor(Ex_sim, Ey_sim, Ez_sim, charge_pos, grid_size,
                        q=1.0, kappa=0.0, kt=1.0, exclude_radius=2.0):
        """
        Encuentra el factor de escala óptimo entre simulación y teoría

        Args:
            Ex_sim, Ey_sim, Ez_sim: Componentes del campo simulado
            charge_pos: Posición de la carga
            grid_size: Tupla (nx, ny, nz)
            q: Carga
            kappa: Parámetro de Debye
            kt: Energía térmica
            exclude_radius: Radio de exclusión cerca de la carga (lattice units)

        Returns:
            scale_factor: Factor de escala óptimo
        """
        # Calcular campo teórico con factor unitario
        Ex_unit, Ey_unit, Ez_unit = TheoreticalField.debye_huckel_field(
            charge_pos, grid_size, q=1.0, epsilon=1.0, kt=1.0,
            scale_factor=1.0, kappa=kappa
        )

        # Crear máscara para excluir región cercana a la carga
        nx, ny, nz = grid_size
        xc, yc, zc = charge_pos
        x = np.arange(0.5, nx)
        y = np.arange(0.5, ny)
        z = np.arange(0.5, nz)
        X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

        r_squared = (X - xc)**2 + (Y - yc)**2 + (Z - zc)**2
        mask = r_squared > exclude_radius**2

        # Ajuste de mínimos cuadrados
        # scale_factor = Σ(E_sim · E_unit) / Σ(E_unit²)
        numerator = (np.sum(Ex_sim[mask] * Ex_unit[mask]) +
                     np.sum(Ey_sim[mask] * Ey_unit[mask]) +
                     np.sum(Ez_sim[mask] * Ez_unit[mask]))

        denominator = (np.sum(Ex_unit[mask]**2) +
                       np.sum(Ey_unit[mask]**2) +
                       np.sum(Ez_unit[mask]**2))

        scale_factor = numerator / denominator

        return scale_factor


class SymmetryAnalyzer:
    """Analiza la simetría del campo eléctrico"""

    def __init__(self, psi, Ex, Ey, Ez, grid_size, particle_pos):
        """
        Inicializa el analizador de simetría

        Args:
            psi: Potencial eléctrico (array 3D)
            Ex, Ey, Ez: Componentes del campo eléctrico (arrays 3D)
            grid_size: Tupla (nx, ny, nz)
            particle_pos: Posición de la partícula (x, y, z)
        """
        self.psi = psi
        self.Ex = Ex
        self.Ey = Ey
        self.Ez = Ez
        self.nx, self.ny, self.nz = grid_size
        self.particle_pos = particle_pos

        # Coordenadas de la malla (centradas en nodos)
        self.x_coords = np.arange(0.5, self.nx)
        self.y_coords = np.arange(0.5, self.ny)
        self.z_coords = np.arange(0.5, self.nz)

        print(f"\nPosición de la partícula: ({particle_pos[0]:.2f}, {particle_pos[1]:.2f}, {particle_pos[2]:.2f})")

    def check_radial_symmetry(self, num_radii=10, max_radius=15.0,
                             Ex_theory=None, Ey_theory=None, Ez_theory=None,
                             psi_theory=None, shell_thickness=0.1):
        """
        Verifica la simetría radial del campo eléctrico y compara con teoría

        Args:
            num_radii: Número de radios a analizar
            max_radius: Radio máximo a analizar
            Ex_theory, Ey_theory, Ez_theory: Campos teóricos (opcional)
            psi_theory: Potencial teórico (opcional)
            shell_thickness: Grosor de la cáscara esférica para muestreo (default: 0.1)

        Returns:
            dict: Diccionario con resultados del análisis
        """
        print("\n=== ANÁLISIS DE SIMETRÍA RADIAL ===")
        print(f"Grosor de cáscaras esféricas: ±{shell_thickness/2:.3f} (grosor total: {shell_thickness:.3f})")

        radii = np.linspace(1.0, max_radius, num_radii)
        px, py, pz = self.particle_pos

        # Verificar si hay campo teórico para comparar
        has_theory = (Ex_theory is not None and Ey_theory is not None and
                     Ez_theory is not None)

        results = {
            'radii': radii,
            'E_magnitude_mean': [],
            'E_magnitude_std': [],
            'Ex_mean': [],
            'Ex_std': [],
            'Ey_mean': [],
            'Ey_std': [],
            'Ez_mean': [],
            'Ez_std': [],
            'psi_mean': [],
            'psi_std': [],
            'num_points': []
        }

        if has_theory:
            results['E_theory_mean'] = []
            results['E_theory_std'] = []
            results['Ex_theory_mean'] = []
            results['Ey_theory_mean'] = []
            results['Ez_theory_mean'] = []
            results['psi_theory_mean'] = []
            results['psi_theory_std'] = []
            results['relative_error'] = []

        for r in radii:
            # Buscar todos los puntos a distancia r de la partícula
            E_mag_at_r = []
            Ex_at_r = []
            Ey_at_r = []
            Ez_at_r = []
            psi_at_r = []
            E_theory_at_r = []
            Ex_theory_at_r = []
            Ey_theory_at_r = []
            Ez_theory_at_r = []
            psi_theory_at_r = []

            for i, x in enumerate(self.x_coords):
                for j, y in enumerate(self.y_coords):
                    for k, z in enumerate(self.z_coords):
                        # Calcular distancia al centro de la partícula
                        dist = np.sqrt((x - px)**2 + (y - py)**2 + (z - pz)**2)

                        # Si está en el rango r ± shell_thickness/2
                        if abs(dist - r) < shell_thickness / 2:
                            E_mag = np.sqrt(self.Ex[i, j, k]**2 +
                                          self.Ey[i, j, k]**2 +
                                          self.Ez[i, j, k]**2)
                            E_mag_at_r.append(E_mag)
                            Ex_at_r.append(self.Ex[i, j, k])
                            Ey_at_r.append(self.Ey[i, j, k])
                            Ez_at_r.append(self.Ez[i, j, k])
                            psi_at_r.append(self.psi[i, j, k])

                            if has_theory:
                                E_th_mag = np.sqrt(Ex_theory[i, j, k]**2 +
                                                  Ey_theory[i, j, k]**2 +
                                                  Ez_theory[i, j, k]**2)
                                E_theory_at_r.append(E_th_mag)
                                Ex_theory_at_r.append(Ex_theory[i, j, k])
                                Ey_theory_at_r.append(Ey_theory[i, j, k])
                                Ez_theory_at_r.append(Ez_theory[i, j, k])
                                if psi_theory is not None:
                                    psi_theory_at_r.append(psi_theory[i, j, k])

            if len(E_mag_at_r) > 0:
                E_mean = np.mean(E_mag_at_r)
                E_std = np.std(E_mag_at_r)

                results['E_magnitude_mean'].append(E_mean)
                results['E_magnitude_std'].append(E_std)
                results['Ex_mean'].append(np.mean(Ex_at_r))
                results['Ex_std'].append(np.std(Ex_at_r))
                results['Ey_mean'].append(np.mean(Ey_at_r))
                results['Ey_std'].append(np.std(Ey_at_r))
                results['Ez_mean'].append(np.mean(Ez_at_r))
                results['Ez_std'].append(np.std(Ez_at_r))
                results['psi_mean'].append(np.mean(psi_at_r))
                results['psi_std'].append(np.std(psi_at_r))
                results['num_points'].append(len(E_mag_at_r))

                # Calcular coeficiente de variación (CV = std/mean)
                cv_E = E_std / E_mean if E_mean != 0 else 0
                cv_psi = np.std(psi_at_r) / abs(np.mean(psi_at_r)) if np.mean(psi_at_r) != 0 else 0

                # Información de teoría si está disponible
                theory_str = ""
                if has_theory and len(E_theory_at_r) > 0:
                    E_th_mean = np.mean(E_theory_at_r)
                    E_th_std = np.std(E_theory_at_r)
                    rel_error = abs(E_mean - E_th_mean) / E_th_mean if E_th_mean != 0 else 0

                    results['E_theory_mean'].append(E_th_mean)
                    results['E_theory_std'].append(E_th_std)
                    results['Ex_theory_mean'].append(np.mean(Ex_theory_at_r))
                    results['Ey_theory_mean'].append(np.mean(Ey_theory_at_r))
                    results['Ez_theory_mean'].append(np.mean(Ez_theory_at_r))
                    results['relative_error'].append(rel_error)

                    if psi_theory is not None and len(psi_theory_at_r) > 0:
                        results['psi_theory_mean'].append(np.mean(psi_theory_at_r))
                        results['psi_theory_std'].append(np.std(psi_theory_at_r))
                    else:
                        results['psi_theory_mean'].append(np.nan)
                        results['psi_theory_std'].append(np.nan)

                    theory_str = f" | Teoría={E_th_mean:.6e} (error={rel_error*100:.2f}%)"

                print(f"Radio {r:.2f}: {len(E_mag_at_r)} puntos | "
                      f"|E| = {E_mean:.6e} ± {E_std:.6e} (CV={cv_E:.3f}) | "
                      f"ψ = {np.mean(psi_at_r):.6e} ± {np.std(psi_at_r):.6e} (CV={cv_psi:.3f}){theory_str}")
            else:
                results['E_magnitude_mean'].append(np.nan)
                results['E_magnitude_std'].append(np.nan)
                results['Ex_mean'].append(np.nan)
                results['Ex_std'].append(np.nan)
                results['Ey_mean'].append(np.nan)
                results['Ey_std'].append(np.nan)
                results['Ez_mean'].append(np.nan)
                results['Ez_std'].append(np.nan)
                results['psi_mean'].append(np.nan)
                results['psi_std'].append(np.nan)
                results['num_points'].append(0)

                if has_theory:
                    results['E_theory_mean'].append(np.nan)
                    results['E_theory_std'].append(np.nan)
                    results['Ex_theory_mean'].append(np.nan)
                    results['Ey_theory_mean'].append(np.nan)
                    results['Ez_theory_mean'].append(np.nan)
                    results['psi_theory_mean'].append(np.nan)
                    results['psi_theory_std'].append(np.nan)
                    results['relative_error'].append(np.nan)

                print(f"Radio {r:.2f}: Sin puntos en este radio")

        return results

    def check_plane_symmetry(self, plane='xy'):
        """
        Verifica la simetría en un plano específico

        Args:
            plane: 'xy', 'xz', o 'yz'

        Returns:
            dict: Diccionario con resultados del análisis
        """
        print(f"\n=== ANÁLISIS DE SIMETRÍA EN PLANO {plane.upper()} ===")

        px, py, pz = self.particle_pos

        # Determinar el índice del plano que pasa por la partícula
        if plane == 'xy':
            pos_idx = int(pz)
            E_plane = np.sqrt(self.Ex[:, :, pos_idx]**2 +
                            self.Ey[:, :, pos_idx]**2 +
                            self.Ez[:, :, pos_idx]**2)
            psi_plane = self.psi[:, :, pos_idx]
            coord1, coord2 = self.x_coords, self.y_coords
            p1, p2 = px, py
            label1, label2 = 'X', 'Y'

        elif plane == 'xz':
            pos_idx = int(py)
            E_plane = np.sqrt(self.Ex[:, pos_idx, :]**2 +
                            self.Ey[:, pos_idx, :]**2 +
                            self.Ez[:, pos_idx, :]**2)
            psi_plane = self.psi[:, pos_idx, :]
            coord1, coord2 = self.x_coords, self.z_coords
            p1, p2 = px, pz
            label1, label2 = 'X', 'Z'

        elif plane == 'yz':
            pos_idx = int(px)
            E_plane = np.sqrt(self.Ex[pos_idx, :, :]**2 +
                            self.Ey[pos_idx, :, :]**2 +
                            self.Ez[pos_idx, :, :]**2)
            psi_plane = self.psi[pos_idx, :, :]
            coord1, coord2 = self.y_coords, self.z_coords
            p1, p2 = py, pz
            label1, label2 = 'Y', 'Z'
        else:
            raise ValueError("plane debe ser 'xy', 'xz', o 'yz'")

        # Analizar simetría en diferentes direcciones radiales
        num_angles = 8
        angles = np.linspace(0, 2*np.pi, num_angles, endpoint=False)
        max_dist = min(self.nx, self.ny, self.nz) // 2

        results = {
            'plane': plane,
            'position': pos_idx,
            'angles': angles,
            'profiles': [],
            'particle_pos': (p1, p2)
        }

        for angle in angles:
            # Dirección radial
            dx = np.cos(angle)
            dy = np.sin(angle)

            # Muestrear a lo largo de esta dirección
            distances = np.linspace(0.5, max_dist, 50)
            E_profile = []
            psi_profile = []

            for dist in distances:
                # Posición en el plano
                x_pos = p1 + dist * dx
                y_pos = p2 + dist * dy

                # Interpolar valor en esta posición
                if (coord1[0] <= x_pos <= coord1[-1] and
                    coord2[0] <= y_pos <= coord2[-1]):
                    # Índices más cercanos
                    i = np.argmin(np.abs(coord1 - x_pos))
                    j = np.argmin(np.abs(coord2 - y_pos))

                    E_profile.append(E_plane[i, j])
                    psi_profile.append(psi_plane[i, j])
                else:
                    E_profile.append(np.nan)
                    psi_profile.append(np.nan)

            results['profiles'].append({
                'angle': angle,
                'distances': distances,
                'E_magnitude': np.array(E_profile),
                'psi': np.array(psi_profile)
            })

        print(f"Plano {plane.upper()} en índice {pos_idx}")
        print(f"Partícula en ({p1:.2f}, {p2:.2f})")
        print(f"Analizadas {num_angles} direcciones radiales")

        return results

    def plot_radial_symmetry(self, results, output=None, kappa=None, log_scale='none', psi_offset_factor=1.1, psi_theory_mean_offset=0.0):
        """Grafica los resultados del análisis de simetría radial con comparación teórica

        Args:
            results: Diccionario con resultados del análisis radial
            output: Ruta del archivo de salida (None para mostrar en pantalla)
            kappa: Parámetro de Debye (para mostrar en el título)
            log_scale: 'none' (lineal), 'y' (log solo en Y), 'xy' (log-log en ambos ejes)
            psi_offset_factor: Factor para calcular offset: offset = -psi_min * factor
            psi_theory_mean_offset: Offset de media de teoría para recuperar valores absolutos en log
        """
        # Verificar si hay datos teóricos
        has_theory = 'E_theory_mean' in results

        if has_theory:
            # Layout 2x2: Campo (arriba izq), Potencial (arriba der),
            #             Error E (abajo izq), Error psi (abajo der)
            fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
        else:
            # Sin teoría: solo campo y potencial lado a lado
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

        radii = results['radii']
        E_mean = results['E_magnitude_mean']
        E_std = results['E_magnitude_std']
        psi_mean = results['psi_mean']
        psi_std = results['psi_std']

        # Gráfico de magnitud del campo eléctrico
        ax1.errorbar(radii, E_mean, yerr=E_std, fmt='o-', linewidth=2,
                    markersize=8, capsize=5, label='Simulación |E| ± σ', color='blue')

        # Agregar campo teórico si está disponible
        if has_theory:
            E_theory = results['E_theory_mean']
            E_theory_std = results['E_theory_std']
            model_name = 'Debye-Hückel' if (kappa is not None and kappa > 0) else 'Coulomb'
            ax1.errorbar(radii, E_theory, yerr=E_theory_std, fmt='s--', linewidth=2,
                        markersize=6, capsize=5, label=f'Teoría ({model_name}) ± σ', color='red')

        ax1.set_xlabel('Distancia desde la partícula', fontsize=14)
        ax1.set_ylabel('Magnitud del campo eléctrico |E|', fontsize=14)
        title_str = 'Simetría Radial del Campo Eléctrico'
        if has_theory and kappa is not None and kappa > 0:
            title_str += f' (κ={kappa:.4f}, λ_D={1.0/kappa:.4f})'
        if log_scale == 'y':
            title_str += ' (escala log-Y)'
        elif log_scale == 'xy':
            title_str += ' (escala log-log)'
        ax1.set_title(title_str, fontsize=16)
        ax1.grid(True, alpha=0.3)
        ax1.legend(fontsize=12)

        # Aplicar escala logarítmica si se solicita
        if log_scale == 'y':
            ax1.set_yscale('log')
            ax1.set_ylabel('Magnitud del campo eléctrico |E| (escala log)', fontsize=14)
        elif log_scale == 'xy':
            ax1.set_xscale('log')
            ax1.set_yscale('log')
            ax1.set_xlabel('Distancia desde la partícula (escala log)', fontsize=14)
            ax1.set_ylabel('Magnitud del campo eléctrico |E| (escala log)', fontsize=14)

        # Añadir línea de ajuste 1/r² solo si NO hay teoría
        if not has_theory:
            valid_idx = ~np.isnan(E_mean)
            if np.any(valid_idx):
                r_valid = np.array(radii)[valid_idx]
                E_valid = np.array(E_mean)[valid_idx]
                # Ajustar a E = k/r²
                from scipy.optimize import curve_fit
                try:
                    def coulomb(r, k):
                        return k / r**2
                    popt, _ = curve_fit(coulomb, r_valid, E_valid)
                    r_fit = np.linspace(r_valid[0], r_valid[-1], 100)
                    ax1.plot(r_fit, coulomb(r_fit, *popt), '--',
                            label=f'Ajuste: k/r² (k={popt[0]:.4e})', linewidth=2, color='green')
                    ax1.legend(fontsize=12)
                except:
                    pass

        # Gráfico de potencial
        # Si se usa escala logarítmica, necesitamos desplazar el potencial para que sea positivo
        psi_offset = 0.0

        # Preparar datos para graficar
        if log_scale in ['y', 'xy']:
            # En escala log, queremos recuperar los valores absolutos originales
            # Primero sumamos el offset de la teoría a ambas curvas
            # Esto las pone en la misma referencia absoluta (la de la teoría original)
            psi_mean_with_theory_ref = np.array(psi_mean) + psi_theory_mean_offset
            if has_theory and 'psi_theory_mean' in results:
                psi_theory_mean_values = results['psi_theory_mean']
                psi_theory_with_ref = np.array(psi_theory_mean_values) + psi_theory_mean_offset
            else:
                psi_theory_with_ref = None

            # Ahora encontrar el mínimo para calcular un offset adicional si es necesario
            psi_min = np.nanmin(psi_mean_with_theory_ref)
            if psi_theory_with_ref is not None:
                psi_min = min(psi_min, np.nanmin(psi_theory_with_ref))

            # Si aún hay valores negativos, aplicar offset adicional
            if psi_min < 0:
                psi_offset = -psi_min * psi_offset_factor
            else:
                psi_offset = 0.0

            # Aplicar offset final
            psi_mean_plot = psi_mean_with_theory_ref + psi_offset
            if psi_theory_with_ref is not None:
                psi_theory_mean_plot = psi_theory_with_ref + psi_offset
            else:
                psi_theory_mean_plot = None

            # El offset total es la suma de ambos
            psi_offset = psi_theory_mean_offset + psi_offset
        else:
            # Escala lineal: usar valores originales
            psi_mean_plot = psi_mean
            if has_theory and 'psi_theory_mean' in results:
                psi_theory_mean_plot = results['psi_theory_mean']
            else:
                psi_theory_mean_plot = None

        ax2.errorbar(radii, psi_mean_plot, yerr=psi_std, fmt='o-', linewidth=2,
                    markersize=8, capsize=5, color='purple', label='Simulación ψ ± σ')

        # Agregar potencial teórico si está disponible
        if has_theory and psi_theory_mean_plot is not None:
            psi_theory_std = results['psi_theory_std']
            model_name = 'Debye-Hückel' if (kappa is not None and kappa > 0) else 'Coulomb'
            ax2.errorbar(radii, psi_theory_mean_plot, yerr=psi_theory_std, fmt='s--', linewidth=2,
                        markersize=6, capsize=5, label=f'Teoría ({model_name}) ± σ', color='orange')

        ax2.set_xlabel('Distancia desde la partícula', fontsize=14)
        ylabel_psi = 'Potencial eléctrico ψ'
        psi_title = 'Simetría Radial del Potencial'
        if log_scale == 'y':
            psi_title += ' (escala log-Y)'
            ylabel_psi += ' (escala log)'
            if psi_offset > 0:
                psi_title += f' [offset +{psi_offset:.4f}]'
        elif log_scale == 'xy':
            psi_title += ' (escala log-log)'
            ylabel_psi += ' (escala log)'
            if psi_offset > 0:
                psi_title += f' [offset +{psi_offset:.4f}]'
        ax2.set_ylabel(ylabel_psi, fontsize=14)
        ax2.set_title(psi_title, fontsize=16)
        ax2.grid(True, alpha=0.3)
        ax2.legend(fontsize=12)

        # Aplicar escala logarítmica si se solicita
        if log_scale == 'y':
            ax2.set_yscale('log')
        elif log_scale == 'xy':
            ax2.set_xscale('log')
            ax2.set_yscale('log')
            ax2.set_xlabel('Distancia desde la partícula (escala log)', fontsize=14)

        # Gráfico de error relativo del campo eléctrico (subplot inferior izquierdo)
        if has_theory:
            rel_error_E = results['relative_error']
            ax3.plot(radii, np.array(rel_error_E) * 100, 'o-', linewidth=2,
                    markersize=8, color='green')
            ax3.set_xlabel('Distancia desde la partícula', fontsize=14)
            ax3.set_ylabel('Error Relativo (%)', fontsize=14)
            ax3.set_title('Error Relativo del Campo: |E_sim - E_theory| / E_theory', fontsize=14)
            ax3.grid(True, alpha=0.3)
            ax3.axhline(y=0, color='k', linestyle='--', linewidth=1, alpha=0.5)

            # Líneas de referencia
            ax3.axhline(y=5, color='orange', linestyle=':', linewidth=1, alpha=0.5, label='±5%')
            ax3.axhline(y=10, color='red', linestyle=':', linewidth=1, alpha=0.5, label='±10%')
            ax3.legend(fontsize=10)

            # Gráfico de error relativo del potencial (subplot inferior derecho)
            # Calcular error relativo del potencial
            psi_sim = results['psi_mean']
            psi_theory = results['psi_theory_mean']
            rel_error_psi = []
            for ps, pt in zip(psi_sim, psi_theory):
                if not np.isnan(ps) and not np.isnan(pt) and pt != 0:
                    rel_error_psi.append(abs(ps - pt) / abs(pt))
                else:
                    rel_error_psi.append(np.nan)

            ax4.plot(radii, np.array(rel_error_psi) * 100, 'o-', linewidth=2,
                    markersize=8, color='purple')
            ax4.set_xlabel('Distancia desde la partícula', fontsize=14)
            ax4.set_ylabel('Error Relativo (%)', fontsize=14)
            ax4.set_title('Error Relativo del Potencial: |ψ_sim - ψ_theory| / |ψ_theory|', fontsize=14)
            ax4.grid(True, alpha=0.3)
            ax4.axhline(y=0, color='k', linestyle='--', linewidth=1, alpha=0.5)

            # Líneas de referencia
            ax4.axhline(y=5, color='orange', linestyle=':', linewidth=1, alpha=0.5, label='±5%')
            ax4.axhline(y=10, color='red', linestyle=':', linewidth=1, alpha=0.5, label='±10%')
            ax4.legend(fontsize=10)

        # Ajustar a ψ = k/r (Coulomb) solo si no hay teoría
        if not has_theory:
            valid_idx = ~np.isnan(psi_mean)
            if np.any(valid_idx):
                r_valid = np.array(radii)[valid_idx]
                psi_valid = np.array(psi_mean)[valid_idx]
                from scipy.optimize import curve_fit
                try:
                    def coulomb_potential(r, k):
                        return k / r
                    popt, _ = curve_fit(coulomb_potential, r_valid, psi_valid)
                    r_fit = np.linspace(r_valid[0], r_valid[-1], 100)
                    ax2.plot(r_fit, coulomb_potential(r_fit, *popt), '--',
                            label=f'Ajuste: k/r (k={popt[0]:.4e})', linewidth=2)
                    ax2.legend(fontsize=12)
                except:
                    pass

        plt.tight_layout()

        if output:
            plt.savefig(output, dpi=300, bbox_inches='tight')
            print(f"\nGráfico de simetría radial guardado en: {output}")
        else:
            plt.show()

        plt.close()

    def plot_components(self, results, output=None, kappa=None):
        """Grafica las componentes individuales del campo eléctrico con comparación teórica"""
        has_theory = 'Ex_theory_mean' in results

        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 14))

        radii = results['radii']
        Ex_mean = results['Ex_mean']
        Ex_std = results['Ex_std']
        Ey_mean = results['Ey_mean']
        Ey_std = results['Ey_std']
        Ez_mean = results['Ez_mean']
        Ez_std = results['Ez_std']

        # Componente Ex
        ax1.errorbar(radii, Ex_mean, yerr=Ex_std, fmt='o-', linewidth=2,
                    markersize=8, capsize=5, label='Simulación Ex ± σ', color='red')
        if has_theory:
            Ex_theory = results['Ex_theory_mean']
            model_name = 'Debye-Hückel' if (kappa is not None and kappa > 0) else 'Coulomb'
            ax1.plot(radii, Ex_theory, 's--', linewidth=2, markersize=6,
                    label=f'Teoría ({model_name})', color='darkred')
        ax1.set_xlabel('Distancia desde la partícula', fontsize=14)
        ax1.set_ylabel('Ex', fontsize=14)
        ax1.set_title('Componente Ex del Campo Eléctrico', fontsize=16)
        ax1.grid(True, alpha=0.3)
        ax1.legend(fontsize=12)
        ax1.axhline(y=0, color='k', linestyle='--', linewidth=0.5, alpha=0.5)

        # Componente Ey
        ax2.errorbar(radii, Ey_mean, yerr=Ey_std, fmt='o-', linewidth=2,
                    markersize=8, capsize=5, label='Simulación Ey ± σ', color='green')
        if has_theory:
            Ey_theory = results['Ey_theory_mean']
            ax2.plot(radii, Ey_theory, 's--', linewidth=2, markersize=6,
                    label=f'Teoría ({model_name})', color='darkgreen')
        ax2.set_xlabel('Distancia desde la partícula', fontsize=14)
        ax2.set_ylabel('Ey', fontsize=14)
        ax2.set_title('Componente Ey del Campo Eléctrico', fontsize=16)
        ax2.grid(True, alpha=0.3)
        ax2.legend(fontsize=12)
        ax2.axhline(y=0, color='k', linestyle='--', linewidth=0.5, alpha=0.5)

        # Componente Ez
        ax3.errorbar(radii, Ez_mean, yerr=Ez_std, fmt='o-', linewidth=2,
                    markersize=8, capsize=5, label='Simulación Ez ± σ', color='blue')
        if has_theory:
            Ez_theory = results['Ez_theory_mean']
            ax3.plot(radii, Ez_theory, 's--', linewidth=2, markersize=6,
                    label=f'Teoría ({model_name})', color='darkblue')
        ax3.set_xlabel('Distancia desde la partícula', fontsize=14)
        ax3.set_ylabel('Ez', fontsize=14)
        ax3.set_title('Componente Ez del Campo Eléctrico', fontsize=16)
        ax3.grid(True, alpha=0.3)
        ax3.legend(fontsize=12)
        ax3.axhline(y=0, color='k', linestyle='--', linewidth=0.5, alpha=0.5)

        title_str = 'Componentes del Campo Eléctrico vs Radio'
        if has_theory and kappa is not None and kappa > 0:
            title_str += f' (κ={kappa:.4f}, λ_D={1.0/kappa:.4f})'
        plt.suptitle(title_str, fontsize=16, y=0.995)

        plt.tight_layout()

        if output:
            plt.savefig(output, dpi=300, bbox_inches='tight')
            print(f"Gráfico de componentes guardado en: {output}")
        else:
            plt.show()

        plt.close()

    def save_radial_profiles(self, results, output_prefix):
        """
        Guarda los perfiles radiales en archivos CSV separados

        Args:
            results: Diccionario con resultados del análisis radial
            output_prefix: Prefijo para los archivos de salida
        """
        import pandas as pd

        has_theory = 'Ex_theory_mean' in results
        radii = results['radii']

        # Crear DataFrame para cada componente
        # Archivo para Ex
        df_Ex = pd.DataFrame({
            'radius': radii,
            'Ex_mean': results['Ex_mean'],
            'Ex_std': results['Ex_std'],
            'num_points': results['num_points']
        })
        if has_theory:
            df_Ex['Ex_theory'] = results['Ex_theory_mean']

        filename_Ex = f"{output_prefix}Ex_profile.csv"
        df_Ex.to_csv(filename_Ex, index=False)
        print(f"Perfil de Ex guardado en: {filename_Ex}")

        # Archivo para Ey
        df_Ey = pd.DataFrame({
            'radius': radii,
            'Ey_mean': results['Ey_mean'],
            'Ey_std': results['Ey_std'],
            'num_points': results['num_points']
        })
        if has_theory:
            df_Ey['Ey_theory'] = results['Ey_theory_mean']

        filename_Ey = f"{output_prefix}Ey_profile.csv"
        df_Ey.to_csv(filename_Ey, index=False)
        print(f"Perfil de Ey guardado en: {filename_Ey}")

        # Archivo para Ez
        df_Ez = pd.DataFrame({
            'radius': radii,
            'Ez_mean': results['Ez_mean'],
            'Ez_std': results['Ez_std'],
            'num_points': results['num_points']
        })
        if has_theory:
            df_Ez['Ez_theory'] = results['Ez_theory_mean']

        filename_Ez = f"{output_prefix}Ez_profile.csv"
        df_Ez.to_csv(filename_Ez, index=False)
        print(f"Perfil de Ez guardado en: {filename_Ez}")

        # Archivo para módulo
        df_mag = pd.DataFrame({
            'radius': radii,
            'E_magnitude_mean': results['E_magnitude_mean'],
            'E_magnitude_std': results['E_magnitude_std'],
            'num_points': results['num_points']
        })
        if has_theory:
            df_mag['E_theory'] = results['E_theory_mean']
            df_mag['relative_error'] = results['relative_error']

        filename_mag = f"{output_prefix}E_magnitude_profile.csv"
        df_mag.to_csv(filename_mag, index=False)
        print(f"Perfil de |E| guardado en: {filename_mag}")

        # Archivo para potencial
        df_psi = pd.DataFrame({
            'radius': radii,
            'psi_mean': results['psi_mean'],
            'psi_std': results['psi_std'],
            'num_points': results['num_points']
        })

        filename_psi = f"{output_prefix}psi_profile.csv"
        df_psi.to_csv(filename_psi, index=False)
        print(f"Perfil de ψ guardado en: {filename_psi}")

    def plot_plane_symmetry(self, results, output=None):
        """Grafica los resultados del análisis de simetría en plano"""
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

        plane = results['plane']
        profiles = results['profiles']
        p1, p2 = results['particle_pos']

        # Gráfico 1: Perfiles radiales de |E|
        for profile in profiles:
            angle = profile['angle']
            distances = profile['distances']
            E_mag = profile['E_magnitude']
            label = f'{np.degrees(angle):.0f}°'
            ax1.plot(distances, E_mag, '-', linewidth=2, label=label, alpha=0.7)

        ax1.set_xlabel('Distancia desde la partícula', fontsize=14)
        ax1.set_ylabel('Magnitud del campo |E|', fontsize=14)
        ax1.set_title(f'Perfiles Radiales de |E| en Plano {plane.upper()}', fontsize=16)
        ax1.grid(True, alpha=0.3)
        ax1.legend(fontsize=10, ncol=2)

        # Gráfico 2: Perfiles radiales de ψ
        for profile in profiles:
            angle = profile['angle']
            distances = profile['distances']
            psi = profile['psi']
            label = f'{np.degrees(angle):.0f}°'
            ax2.plot(distances, psi, '-', linewidth=2, label=label, alpha=0.7)

        ax2.set_xlabel('Distancia desde la partícula', fontsize=14)
        ax2.set_ylabel('Potencial ψ', fontsize=14)
        ax2.set_title(f'Perfiles Radiales de ψ en Plano {plane.upper()}', fontsize=16)
        ax2.grid(True, alpha=0.3)
        ax2.legend(fontsize=10, ncol=2)

        plt.tight_layout()

        if output:
            plt.savefig(output, dpi=300, bbox_inches='tight')
            print(f"Gráfico de simetría en plano guardado en: {output}")
        else:
            plt.show()

        plt.close()

    def compute_symmetry_metrics(self, radial_results, plane_results):
        """
        Calcula métricas cuantitativas de simetría

        Args:
            radial_results: Resultados del análisis radial
            plane_results: Lista de resultados de análisis en planos

        Returns:
            dict: Métricas de simetría
        """
        print("\n=== MÉTRICAS DE SIMETRÍA ===")

        metrics = {}

        # Métrica 1: Coeficiente de variación promedio en función del radio
        E_std = np.array(radial_results['E_magnitude_std'])
        E_mean = np.array(radial_results['E_magnitude_mean'])
        valid = ~np.isnan(E_std) & ~np.isnan(E_mean) & (E_mean != 0)

        if np.any(valid):
            cv_radial = E_std[valid] / E_mean[valid]
            metrics['cv_radial_mean'] = np.mean(cv_radial)
            metrics['cv_radial_max'] = np.max(cv_radial)
            print(f"Coeficiente de variación radial (|E|):")
            print(f"  Promedio: {metrics['cv_radial_mean']:.4f}")
            print(f"  Máximo: {metrics['cv_radial_max']:.4f}")

        # Métrica 2: Dispersión entre perfiles angulares
        metrics['plane_dispersion'] = {}
        for plane_result in plane_results:
            plane = plane_result['plane']
            profiles = plane_result['profiles']

            # Calcular la dispersión entre perfiles a cada distancia
            num_distances = len(profiles[0]['distances'])
            dispersions = []

            for i in range(num_distances):
                values = [p['E_magnitude'][i] for p in profiles]
                values = [v for v in values if not np.isnan(v)]
                if len(values) > 1:
                    dispersions.append(np.std(values) / np.mean(values) if np.mean(values) != 0 else 0)

            if dispersions:
                metrics['plane_dispersion'][plane] = {
                    'mean': np.mean(dispersions),
                    'max': np.max(dispersions)
                }
                print(f"\nDispersión angular en plano {plane.upper()}:")
                print(f"  Promedio: {metrics['plane_dispersion'][plane]['mean']:.4f}")
                print(f"  Máximo: {metrics['plane_dispersion'][plane]['max']:.4f}")

        return metrics


def main():
    parser = argparse.ArgumentParser(
        description='Verificar simetría del campo eléctrico respecto a la posición de la partícula',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Ejemplos de uso:

ANÁLISIS BÁSICO:
1. Análisis completo usando archivo config.cds:
   ./verificar_simetria_psi.py -f psi-000005000.001-001 -c config.cds00005000.001-001 -s 32 32 32

2. Análisis usando archivo colloids.csv:
   ./verificar_simetria_psi.py -f psi-000005000.001-001 -c colloids-000005000.csv -s 32 32 32

3. Especificar posición de la partícula manualmente:
   ./verificar_simetria_psi.py -f psi-000005000.001-001 -s 32 32 32 -p 16.1 16.5 16.5

COMPARACIÓN CON MODELO TEÓRICO (POISSON-BOLTZMANN):
4. Comparar con modelo de Coulomb (vacío, sin iones):
   ./verificar_simetria_psi.py -f psi-000005000.001-001 -c config.cds00005000.001-001 -s 32 32 32 \
       --compare-theory --charge 1.0 --epsilon 1.0 --kappa 0.0 -o simetria_coulomb_

5. Comparar con modelo de Debye-Hückel (electrolito con iones):
   ./verificar_simetria_psi.py -f psi-000005000.001-001 -c config.cds00005000.001-001 -s 32 32 32 \
       --compare-theory --charge 1.0 --epsilon 1.0 --kappa 0.1 -o simetria_DH_

6. Comparar con parámetros de Ludwig (epsilon=1e4, rho_el=2e-3):
   ./verificar_simetria_psi.py -f psi-000005000.001-001 -c config.cds00005000.001-001 -s 32 32 32 \
       --compare-theory --charge 1.0 --epsilon 1e4 --kt 0.0005 --kappa 0.0 -o simetria_

OPCIONES AVANZADAS:
7. Análisis con rango de radios personalizado:
   ./verificar_simetria_psi.py -f psi-000005000.001-001 -c config.cds00005000.001-001 -s 32 32 32 \
       --max-radius 10 --num-radii 15

8. Analizar solo planos específicos:
   ./verificar_simetria_psi.py -f psi-000005000.001-001 -c config.cds00005000.001-001 -s 32 32 32 \
       --planes xy xz

9. Usar escala logarítmica en el eje Y del gráfico:
   ./verificar_simetria_psi.py -f psi-000005000.001-001 -c config.cds00005000.001-001 -s 32 32 32 \
       --compare-theory --kappa 0.1 --log-scale y -o simetria_DH_log_

   Usar escala log-log (ambos ejes) - muestra leyes de potencia como líneas rectas:
   ./verificar_simetria_psi.py -f psi-000005000.001-001 -c config.cds00005000.001-001 -s 32 32 32 \
       --compare-theory --kappa 0.1 --log-scale xy -o simetria_DH_loglog_

10. Usar cáscaras esféricas más finas para mayor precisión:
   ./verificar_simetria_psi.py -f psi-000005000.001-001 -c config.cds00005000.001-001 -s 32 32 32 \
       --compare-theory --kappa 0.1 --shell-thickness 0.1 -o simetria_DH_

11. Solo realizar análisis radial (salta planos, componentes y métricas). Solo guarda el gráfico:
   ./verificar_simetria_psi.py -f psi-000005000.001-001 -c config.cds00005000.001-001 -s 32 32 32 \
       --compare-theory --kappa 0.1 --log-scale xy --radial-only -o simetria_DH_

12. Análisis radial con exportación de datos CSV además del gráfico:
   ./verificar_simetria_psi.py -f psi-000005000.001-001 -c config.cds00005000.001-001 -s 32 32 32 \
       --compare-theory --kappa 0.1 --log-scale xy --radial-only --save-csv -o simetria_DH_

13. Ajustar el offset del potencial en escala logarítmica:
   ./verificar_simetria_psi.py -f psi-000005000.001-001 -c config.cds00005000.001-001 -s 32 32 32 \
       --compare-theory --kappa 0.1 --log-scale xy --psi-offset-factor 1.5 -o simetria_DH_

NOTAS:
  - kappa = 0: Modelo de Coulomb (vacío, sin screening iónico)
  - kappa > 0: Modelo de Debye-Hückel / Poisson-Boltzmann (electrolito)
  - λ_D = 1/κ es la longitud de Debye
  - --compare-theory ajusta automáticamente el factor de escala
  - --shell-thickness: grosor de las cáscaras esféricas. Menor grosor = más precisión pero menos puntos por cáscara
        """
    )

    # Argumentos obligatorios
    parser.add_argument('-f', '--file', required=True,
                       help='Archivo psi a analizar')
    parser.add_argument('-s', '--size', nargs=3, type=int, required=True,
                       metavar=('NX', 'NY', 'NZ'),
                       help='Tamaño de la malla (nx ny nz)')

    # Posición de la partícula
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-c', '--colloid-file',
                       help='Archivo config.cds*.001-001 o colloids*.csv para leer posición de la partícula')
    group.add_argument('-p', '--position', nargs=3, type=float,
                       metavar=('X', 'Y', 'Z'),
                       help='Posición de la partícula (x y z)')

    # Parámetros de análisis
    parser.add_argument('--max-radius', type=float, default=15.0,
                       help='Radio máximo para análisis radial (default: 15.0)')
    parser.add_argument('--num-radii', type=int, default=10,
                       help='Número de radios a analizar (default: 10)')
    parser.add_argument('--shell-thickness', type=float, default=0.1,
                       help='Grosor de las cáscaras esféricas para muestreo radial (default: 0.1)')
    parser.add_argument('--planes', nargs='+', choices=['xy', 'xz', 'yz'],
                       default=['xy', 'xz', 'yz'],
                       help='Planos a analizar (default: xy xz yz)')
    parser.add_argument('--log-scale', choices=['y', 'xy', 'none'], default='none',
                       help='Escala logarítmica: "y" (solo eje Y), "xy" (ambos ejes log-log), "none" (lineal)')
    parser.add_argument('--radial-only', action='store_true',
                       help='Solo realizar análisis radial (salta planos, componentes y métricas). Por defecto solo guarda el gráfico radial.')
    parser.add_argument('--save-csv', action='store_true',
                       help='Guardar perfiles radiales en archivos CSV (además de los gráficos). Por defecto solo se guardan gráficos.')
    parser.add_argument('--psi-offset-factor', type=float, default=1.1,
                       help='Factor para calcular offset del potencial en escala log: offset = -psi_min * factor (default: 1.1)')

    # Parámetros del modelo teórico (Poisson-Boltzmann / Debye-Hückel)
    parser.add_argument('--compare-theory', action='store_true',
                       help='Comparar con modelo teórico de Poisson-Boltzmann (Debye-Hückel)')
    parser.add_argument('--charge', type=float, default=1.0,
                       help='Carga de la partícula (default: 1.0)')
    parser.add_argument('--epsilon', type=float, default=1.0,
                       help='Permitividad relativa del medio (default: 1.0)')
    parser.add_argument('--kt', type=float, default=1.0,
                       help='Energía térmica k_B*T (default: 1.0)')
    parser.add_argument('--kappa', type=float, default=0.0,
                       help='Parámetro de Debye κ = 1/λ_D. κ=0 usa Coulomb (vacío), κ>0 usa Debye-Hückel (electrolito). (default: 0.0)')

    # Salida
    parser.add_argument('-o', '--output-prefix',
                       help='Prefijo para archivos de salida (default: mostrar en pantalla)')

    args = parser.parse_args()

    # Verificar que el archivo existe
    if not os.path.exists(args.file):
        print(f"Error: No se encuentra el archivo {args.file}")
        sys.exit(1)

    # Obtener posición de la partícula
    if args.colloid_file:
        if not os.path.exists(args.colloid_file):
            print(f"Error: No se encuentra el archivo {args.colloid_file}")
            sys.exit(1)
        particle_pos = ColloidReader.read_colloid_position(args.colloid_file)
        print(f"Posición leída desde {args.colloid_file}: {particle_pos}")
    else:
        particle_pos = tuple(args.position)
        print(f"Posición especificada: {particle_pos}")

    # Leer archivo psi
    print(f"\nLeyendo archivo {args.file}...")
    reader = PsiReader(args.file, tuple(args.size))
    psi = reader.read()
    print(f"Archivo leído exitosamente. Malla: {args.size[0]}x{args.size[1]}x{args.size[2]}")

    # Calcular campo eléctrico
    print("Calculando campo eléctrico E = -∇ψ...")
    Ex, Ey, Ez = reader.compute_electric_field()
    print("Campo eléctrico calculado.")

    # Calcular campo teórico si se solicita
    Ex_theory, Ey_theory, Ez_theory = None, None, None
    psi_theory = None
    scale_factor = None

    if args.compare_theory:
        print("\n=== COMPARACIÓN CON MODELO TEÓRICO ===")
        model_name = 'Debye-Hückel (Poisson-Boltzmann linearizado)' if args.kappa > 0 else 'Coulomb'
        print(f"Modelo: {model_name}")
        print(f"Parámetros: q={args.charge}, ε={args.epsilon}, kT={args.kt}, κ={args.kappa}")
        if args.kappa > 0:
            print(f"Longitud de Debye: λ_D = {1.0/args.kappa:.4f}")

        # Ajustar factor de escala
        print("\nAjustando factor de escala entre simulación y teoría...")
        scale_factor = TheoreticalField.fit_scale_factor(
            Ex, Ey, Ez, particle_pos, tuple(args.size),
            q=args.charge, kappa=args.kappa, kt=args.kt
        )
        print(f"Factor de escala encontrado: {scale_factor:.6e}")
        print(f"Factor teórico q/(4πε*kt): {args.charge / (4.0 * np.pi * args.epsilon * args.kt):.6e}")
        print(f"Ratio simulado/teórico: {scale_factor / (args.charge / (4.0 * np.pi * args.epsilon * args.kt)):.4f}")

        # Calcular campo teórico
        print("\nCalculando campo teórico...")
        Ex_theory, Ey_theory, Ez_theory = TheoreticalField.debye_huckel_field(
            particle_pos, tuple(args.size),
            q=args.charge, epsilon=args.epsilon, kt=args.kt,
            kappa=args.kappa, scale_factor=scale_factor
        )
        print("Campo teórico calculado.")

        # Calcular potencial teórico
        print("Calculando potencial teórico...")
        psi_theory = TheoreticalField.debye_huckel_potential(
            particle_pos, tuple(args.size),
            q=args.charge, epsilon=args.epsilon, kt=args.kt,
            kappa=args.kappa, scale_factor=scale_factor
        )

        # Ajustar potencial teórico para tener media cero (igual que la simulación)
        psi_theory_mean_offset = np.mean(psi_theory)
        psi_theory = psi_theory - psi_theory_mean_offset
        print(f"Potencial teórico calculado (ajustado a media cero, offset={psi_theory_mean_offset:.6e}).")

    # Crear analizador de simetría
    analyzer = SymmetryAnalyzer(psi, Ex, Ey, Ez, tuple(args.size), particle_pos)

    # Análisis de simetría radial
    print(f"\nAnalizando simetría radial (radio máximo: {args.max_radius}, {args.num_radii} radios)...")
    radial_results = analyzer.check_radial_symmetry(
        num_radii=args.num_radii,
        max_radius=args.max_radius,
        Ex_theory=Ex_theory,
        Ey_theory=Ey_theory,
        Ez_theory=Ez_theory,
        psi_theory=psi_theory,
        shell_thickness=args.shell_thickness
    )

    # Graficar simetría radial
    if args.output_prefix:
        output_radial = f"{args.output_prefix}radial.png"
    else:
        output_radial = None
    analyzer.plot_radial_symmetry(radial_results, output=output_radial,
                                  kappa=args.kappa if args.compare_theory else None,
                                  log_scale=args.log_scale,
                                  psi_offset_factor=args.psi_offset_factor,
                                  psi_theory_mean_offset=psi_theory_mean_offset if args.compare_theory else 0.0)

    # Si solo se pidió análisis radial, terminar aquí
    if args.radial_only:
        # Guardar perfiles radiales en archivos CSV solo si se especificó --save-csv
        if args.output_prefix and args.save_csv:
            print("\nGuardando perfiles radiales en archivos CSV...")
            analyzer.save_radial_profiles(radial_results, args.output_prefix)
        print("\n¡Análisis radial completado!")
        return

    # Graficar componentes individuales
    if args.output_prefix:
        output_components = f"{args.output_prefix}components.png"
    else:
        output_components = None
    print("\nGenerando gráfico de componentes individuales...")
    analyzer.plot_components(radial_results, output=output_components, kappa=args.kappa if args.compare_theory else None)

    # Guardar perfiles radiales en archivos CSV solo si se especificó --save-csv
    if args.output_prefix and args.save_csv:
        print("\nGuardando perfiles radiales en archivos CSV...")
        analyzer.save_radial_profiles(radial_results, args.output_prefix)

    # Análisis de simetría en planos
    plane_results = []
    for plane in args.planes:
        print(f"\nAnalizando simetría en plano {plane.upper()}...")
        plane_result = analyzer.check_plane_symmetry(plane=plane)
        plane_results.append(plane_result)

        # Graficar simetría en plano
        if args.output_prefix:
            output_plane = f"{args.output_prefix}plane_{plane}.png"
        else:
            output_plane = None
        analyzer.plot_plane_symmetry(plane_result, output=output_plane)

    # Calcular y mostrar métricas de simetría
    metrics = analyzer.compute_symmetry_metrics(radial_results, plane_results)

    # Guardar métricas en archivo si se especificó prefijo
    if args.output_prefix:
        metrics_file = f"{args.output_prefix}metrics.txt"
        with open(metrics_file, 'w') as f:
            f.write("=== MÉTRICAS DE SIMETRÍA DEL CAMPO ELÉCTRICO ===\n\n")
            f.write(f"Archivo PSI: {args.file}\n")
            f.write(f"Posición de la partícula: {particle_pos}\n\n")

            if 'cv_radial_mean' in metrics:
                f.write("Simetría Radial:\n")
                f.write(f"  CV promedio: {metrics['cv_radial_mean']:.6f}\n")
                f.write(f"  CV máximo: {metrics['cv_radial_max']:.6f}\n\n")

            if 'plane_dispersion' in metrics:
                f.write("Simetría Angular por Plano:\n")
                for plane, values in metrics['plane_dispersion'].items():
                    f.write(f"  Plano {plane.upper()}:\n")
                    f.write(f"    Dispersión promedio: {values['mean']:.6f}\n")
                    f.write(f"    Dispersión máxima: {values['max']:.6f}\n")

        print(f"\nMétricas guardadas en: {metrics_file}")

    print("\n¡Análisis de simetría completado!")


if __name__ == '__main__':
    main()
