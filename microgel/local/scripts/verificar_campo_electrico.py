#!/usr/bin/env python3
"""
Script para verificar el campo eléctrico calculado por Ludwig contra el teórico.

Este script:
1. Lee archivos efield de Ludwig
2. Extrae la posición de la partícula del archivo colloids
3. Calcula el módulo del campo eléctrico |E| para TODOS los puntos
4. Ordena por distancia y compara con el campo teórico de Debye-Hückel
5. Genera gráficos scatter de todos los puntos vs teoría

Modelo teórico (campo NO reducido, como en Ludwig):
    - kappa = 0 (sin iones): E = (q/4πε) · 1/r²
    - kappa > 0 (con iones): E = (q/4πε) · exp(-κr) · (κ/r + 1/r²)

NOTA: Ludwig calcula el campo eléctrico E en unidades de fuerza/carga,
      NO el campo reducido E* = E·kT. Por lo tanto, NO se divide por kT.

Autor: Script para análisis de simulaciones Ludwig con Ewald
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
import os
from pathlib import Path


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


class TheoreticalField:
    """Calcula el campo eléctrico teórico según modelo de Debye-Hückel"""

    @staticmethod
    def compute_theory_for_distances(distances, q=1.0, epsilon=1.0, kt=1.0, kappa=0.0):
        """
        Calcula el campo eléctrico teórico para un array de distancias

        Args:
            distances: Array de distancias radiales
            q: Carga
            epsilon: Permitividad relativa
            kt: Energía térmica k_B*T (usado para calcular kappa, NO en el prefactor)
            kappa: Parámetro de Debye κ = 1/λ_D (default: 0.0 = Coulomb)

        Returns:
            E_mag: Módulo del campo eléctrico para cada distancia

        NOTA: Ludwig calcula E sin dividir por kT (unidades de fuerza/carga),
              no el campo reducido E* = E·kT. Por lo tanto, NO incluimos kT
              en el denominador del prefactor, como en ewald_charge.c línea 2824.
        """
        # Prefactor (SIN kT, como en Ludwig ewald_charge.c)
        prefactor = q / (4.0 * np.pi * epsilon)

        # Evitar división por cero
        r = np.maximum(distances, 0.01)

        # Campo eléctrico según el modelo
        if kappa > 0:
            # Debye-Hückel: medio con iones
            # |E| = (q/4πε) · exp(-κr) · (κ/r + 1/r²)
            E_mag = prefactor * np.exp(-kappa * r) * (kappa / r + 1.0 / r**2)
        else:
            # Coulomb: vacío sin screening
            # |E| = (q/4πε) · 1/r²
            E_mag = prefactor / r**2

        return E_mag


def extract_all_points(E_mag_3d, particle_pos, grid_size, min_radius=0.5, max_radius=None):
    """
    Extrae TODOS los puntos de la malla con sus distancias y valores de campo

    Args:
        E_mag_3d: Campo |E| en malla 3D
        particle_pos: Posición de la partícula [x, y, z]
        grid_size: Tupla (nx, ny, nz)
        min_radius: Radio mínimo (default: 0.5)
        max_radius: Radio máximo (default: None = L/2)

    Returns:
        distances: Array con distancias de cada punto
        E_values: Array con |E| de cada punto
    """
    nx, ny, nz = grid_size

    if max_radius is None:
        max_radius = min(nx, ny, nz) / 2.0

    # Crear coordenadas centradas en nodos (0.5, 1.5, 2.5, ...)
    x = np.arange(0.5, nx)
    y = np.arange(0.5, ny)
    z = np.arange(0.5, nz)

    distances = []
    E_values = []

    px, py, pz = particle_pos

    # Iterar sobre todos los puntos
    for i, xi in enumerate(x):
        for j, yj in enumerate(y):
            for k, zk in enumerate(z):
                # Calcular distancia a la partícula
                dist = np.sqrt((xi - px)**2 + (yj - py)**2 + (zk - pz)**2)

                # Filtrar por rango
                if min_radius <= dist <= max_radius:
                    distances.append(dist)
                    E_values.append(E_mag_3d[i, j, k])

    return np.array(distances), np.array(E_values)


def plot_comparison(distances, E_simul, E_theory, charge, epsilon, kt, kappa,
                   output_file='efield_comparison.png', max_points=10000, log_scale=True):
    """
    Genera gráfico comparando TODOS los puntos simulados con el teórico

    Args:
        distances: Distancias de todos los puntos
        E_simul: Campo eléctrico simulado para cada punto
        E_theory: Campo eléctrico teórico para cada punto
        charge: Carga de la partícula
        epsilon: Permitividad
        kt: Energía térmica
        kappa: Parámetro de Debye
        output_file: Nombre del archivo de salida
        max_points: Máximo número de puntos a graficar (para evitar sobrecarga)
        log_scale: Si True, usa escala logarítmica; si False, usa escala lineal
    """
    # Ordenar por distancia
    sort_idx = np.argsort(distances)
    distances_sorted = distances[sort_idx]
    E_simul_sorted = E_simul[sort_idx]
    E_theory_sorted = E_theory[sort_idx]

    # Si hay demasiados puntos, submuestrear uniformemente
    if len(distances_sorted) > max_points:
        step = len(distances_sorted) // max_points
        distances_sorted = distances_sorted[::step]
        E_simul_sorted = E_simul_sorted[::step]
        E_theory_sorted = E_theory_sorted[::step]
        print(f"  Submuestreo: graficando {len(distances_sorted)} de {len(distances)} puntos")

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))

    # Panel 1: Comparación |E|(r) - todos los puntos
    ax1.plot(distances_sorted, E_simul_sorted, 'b.', markersize=1, alpha=0.3,
             label='Ludwig (Ewald) - todos los puntos')
    ax1.plot(distances_sorted, E_theory_sorted, 'r-', linewidth=2, alpha=0.8,
             label='Teórico (Debye-Hückel)')

    ax1.set_xlabel('Distancia r [lattice units]', fontsize=12)
    ax1.set_ylabel('|E(r)| [lattice units]', fontsize=12)

    title = f'Campo Eléctrico: q={charge:.2e}, ε={epsilon:.2e}'
    if kappa > 0:
        title += f', κ={kappa:.2e}, λD={1.0/kappa:.2e}'
    ax1.set_title(title, fontsize=14)

    ax1.legend(fontsize=11, loc='upper right')
    ax1.grid(True, alpha=0.3)
    if log_scale:
        ax1.set_yscale('log')
        ax1.set_xscale('log')
    else:
        ax1.set_yscale('linear')
        ax1.set_xscale('linear')

    # Panel 2: Error relativo punto a punto
    rel_error = np.abs(E_simul_sorted - E_theory_sorted) / (E_theory_sorted + 1e-15) * 100

    ax2.plot(distances_sorted, rel_error, 'bo', markersize=1, alpha=0.5)
    ax2.axhline(y=10, color='r', linestyle=':', alpha=0.5, linewidth=2, label='10% error')
    ax2.axhline(y=1, color='g', linestyle=':', alpha=0.5, linewidth=2, label='1% error')

    ax2.set_xlabel('Distancia r [lattice units]', fontsize=12)
    ax2.set_ylabel('Error relativo [%]', fontsize=12)
    ax2.set_title('Error Relativo punto a punto: |E_simul - E_teórico| / E_teórico × 100', fontsize=14)
    ax2.legend(fontsize=11)
    ax2.grid(True, alpha=0.3)
    if log_scale:
        ax2.set_xscale('log')
        ax2.set_yscale('log')
    else:
        ax2.set_xscale('linear')
        ax2.set_yscale('linear')

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nGráfico guardado en: {output_file}")
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description='Verifica el campo eléctrico de Ludwig contra el modelo de Debye-Hückel',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Ejemplos:
    # Análisis básico (vacío, kappa=0, Coulomb)
    python verificar_campo_electrico.py -e efield-00000001.001-001 -c colloids-00000001.csv \\
        -L 32 -q 1.0 -eps 10000

    # Con screening de Debye (medio con iones, kappa > 0)
    python verificar_campo_electrico.py -e efield-00000001.001-001 -c colloids-00000001.csv \\
        -L 32 -q 1.0 -eps 10000 --kappa 0.1

    # Personalizar rango radial
    python verificar_campo_electrico.py -e efield-00000001.001-001 -c colloids-00000001.csv \\
        -L 32 -q 1.0 -eps 10000 --rmin 1.0 --rmax 15

    # Usar escala lineal en lugar de logarítmica
    python verificar_campo_electrico.py -e efield-00000001.001-001 -c colloids-00000001.csv \\
        -L 32 -q 1.0 -eps 10000 --linear

NOTA: Ludwig calcula el campo E en unidades de fuerza/carga, NO el campo reducido
      E* = E·kT. Por lo tanto, el parámetro -kt NO se usa en la fórmula del campo,
      solo está disponible para referencia.
        """
    )

    parser.add_argument('-e', '--efield', required=True,
                       help='Archivo efield de Ludwig (ej: efield-00000001.001-001)')
    parser.add_argument('-c', '--colloid', required=True,
                       help='Archivo de coloides (ej: colloids-00000001.csv)')
    parser.add_argument('-L', '--grid-size', type=int, required=True,
                       help='Tamaño de la malla (asume cubo LxLxL)')
    parser.add_argument('-q', '--charge', type=float, required=True,
                       help='Carga de la partícula (en unidades de Ludwig)')
    parser.add_argument('-eps', '--epsilon', type=float, required=True,
                       help='Permitividad relativa ε_r')
    parser.add_argument('-kt', '--kt', type=float, default=1.0,
                       help='Energía térmica k_B*T (solo para referencia, NO se usa en fórmula de E)')
    parser.add_argument('--kappa', type=float, default=0.0,
                       help='Parámetro de Debye κ = 1/λ_D (default: 0.0 = Coulomb)')
    parser.add_argument('-o', '--output', default='efield_comparison.png',
                       help='Archivo de salida para el gráfico (default: efield_comparison.png)')
    parser.add_argument('--rmin', type=float, default=0.5,
                       help='Radio mínimo de análisis (default: 0.5)')
    parser.add_argument('--rmax', type=float, default=None,
                       help='Radio máximo de análisis (default: L/2)')
    parser.add_argument('--max-points', type=int, default=10000,
                       help='Máximo número de puntos a graficar (default: 10000)')
    parser.add_argument('--linear', action='store_true',
                       help='Usar escala lineal en los gráficos (default: logarítmica)')

    args = parser.parse_args()

    # Tamaño de la malla
    grid_size = (args.grid_size, args.grid_size, args.grid_size)

    # Leer posición de la partícula
    print(f"Leyendo posición de la partícula desde: {args.colloid}")
    particle_pos = ColloidReader.read_colloid_position(args.colloid)
    print(f"  Posición: {particle_pos}")

    # Leer campo eléctrico de Ludwig
    print(f"\nLeyendo campo eléctrico desde: {args.efield}")
    efield_reader = EFieldReader(args.efield, grid_size)
    try:
        Ex, Ey, Ez = efield_reader.read()
    except ValueError as e:
        print(f"ERROR: {e}")
        sys.exit(1)

    E_mag_simul = efield_reader.compute_magnitude()
    print(f"  Rango de |E|: [{np.min(E_mag_simul):.6e}, {np.max(E_mag_simul):.6e}]")

    # Extraer TODOS los puntos con sus distancias
    print(f"\nExtrayendo todos los puntos (rmin={args.rmin}, rmax={args.rmax})...")
    distances, E_values = extract_all_points(E_mag_simul, particle_pos, grid_size,
                                            min_radius=args.rmin, max_radius=args.rmax)
    print(f"  Total de puntos extraídos: {len(distances)}")
    print(f"  Distancia mínima: {np.min(distances):.4f}")
    print(f"  Distancia máxima: {np.max(distances):.4f}")

    # Calcular campo eléctrico teórico para TODAS las distancias
    modelo = "Debye-Hückel (con screening)" if args.kappa > 0 else "Coulomb (vacío)"
    print(f"\nCalculando campo eléctrico teórico ({modelo}):")
    print(f"  Carga q = {args.charge:.6e}")
    print(f"  Permitividad ε = {args.epsilon:.6e}")
    if args.kappa > 0:
        print(f"  Parámetro Debye κ = {args.kappa:.6e}")
        print(f"  Longitud Debye λD = {1.0/args.kappa:.6e}")
    print(f"  NOTA: Campo NO reducido (E, no E*=E·kT), como en Ludwig")

    E_theory = TheoreticalField.compute_theory_for_distances(
        distances, q=args.charge, epsilon=args.epsilon, kt=args.kt, kappa=args.kappa
    )

    # Estadísticas
    print(f"\nEstadísticas de comparación:")
    rel_error = np.abs(E_values - E_theory) / (E_theory + 1e-15) * 100
    print(f"  Error relativo promedio: {np.mean(rel_error):.2f}%")
    print(f"  Error relativo mediano: {np.median(rel_error):.2f}%")
    print(f"  Error relativo máximo: {np.max(rel_error):.2f}%")
    print(f"  Error relativo mínimo: {np.min(rel_error):.2f}%")

    # Generar gráfico
    log_scale = not args.linear
    scale_str = "lineal" if args.linear else "logarítmica"
    print(f"\nGenerando gráfico con escala {scale_str}...")
    plot_comparison(distances, E_values, E_theory, args.charge, args.epsilon,
                   args.kt, args.kappa, args.output, args.max_points, log_scale)

    print("\n¡Análisis completado!")


if __name__ == '__main__':
    main()
