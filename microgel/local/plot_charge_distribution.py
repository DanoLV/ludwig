#!/usr/bin/env python3
"""
Script para graficar la distribución de carga desde archivos qsi de Ludwig.

Los archivos qsi contienen la densidad de carga en cada punto de la malla.
Formato típico (2 columnas):
  - Columna 0: Densidad de carga positiva ρ₊ (cationes, carga +1)
  - Columna 1: Densidad de carga negativa ρ₋ (aniones, carga -1)

IMPORTANTE: Las densidades son siempre ≥ 0 (representan concentraciones).
La carga neta en cada nodo es: ρ_net = ρ₊ - ρ₋

Modos de visualización:
1. Plano 2D: Muestra la distribución de carga en un plano específico (XY, XZ, o YZ)
2. Línea 1D: Muestra la distribución de carga a lo largo de una línea
3. Plano 3D: Muestra superficie 3D de la densidad de carga en un plano

Funcionalidades adicionales:
- Cálculo de carga total del sistema
- Verificación de electroneutralidad
- Visualización de carga neta (ρ₊ - ρ₋)
- Visualización de densidad total (ρ₊ + ρ₋)
- Visualización de especies individuales

Autor: Script generado para análisis de simulaciones Ludwig
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
import os
from matplotlib.colors import Normalize
from matplotlib import cm

class QsiReader:
    """Lee y procesa archivos qsi de Ludwig"""

    def __init__(self, filename, grid_size, num_species=2):
        """
        Inicializa el lector de archivos qsi

        Args:
            filename: Ruta al archivo qsi
            grid_size: Tupla (nx, ny, nz) con el tamaño de la malla
            num_species: Número de especies iónicas (por defecto: 2)
        """
        self.filename = filename
        self.nx, self.ny, self.nz = grid_size
        self.num_species = num_species
        self.charge_density = None
        self.species = []

    def read(self):
        """Lee el archivo qsi y lo organiza en mallas 3D"""
        data = np.loadtxt(self.filename)

        if data.shape[0] != self.nx * self.ny * self.nz:
            raise ValueError(f"El archivo tiene {data.shape[0]} filas, "
                           f"pero se esperaban {self.nx * self.ny * self.nz}")

        if data.shape[1] != self.num_species:
            print(f"Advertencia: Se encontraron {data.shape[1]} columnas, "
                  f"pero se esperaban {self.num_species}")
            self.num_species = data.shape[1]

        # Organizar cada especie en una malla 3D
        # Ludwig escribe en orden z-major: para cada x, para cada y, todos los z
        self.species = []
        for i in range(self.num_species):
            species_data = data[:, i].reshape((self.nx, self.ny, self.nz))

            # Validar que no haya valores negativos (densidades deben ser ≥ 0)
            min_value = np.min(species_data)
            if min_value < 0:
                print(f"ADVERTENCIA: Especie {i} tiene valores negativos (mín: {min_value:.6e})")
                print(f"Las densidades de carga deben ser no-negativas.")
                print(f"Esto puede indicar un error en la simulación o en el archivo.")

            self.species.append(species_data)

        return self.species

    def compute_net_charge(self):
        """
        Calcula la carga neta en cada punto

        La carga neta es: ρ_net = ρ₊ - ρ₋
        donde:
          ρ₊ = densidad de carga positiva (columna 0, cationes con carga +1)
          ρ₋ = densidad de carga negativa (columna 1, aniones con carga -1)

        Returns:
            Carga neta en cada punto de la malla
        """
        if len(self.species) < 2:
            print("Advertencia: Se necesitan al menos 2 especies para calcular carga neta")
            return self.species[0]

        # ρ_net = ρ₊ - ρ₋
        net_charge = self.species[0] - self.species[1]
        return net_charge

    def compute_total_charge(self):
        """
        Calcula la carga total del sistema para cada especie

        Returns:
            Lista con la carga total de cada especie
        """
        if self.species is None or len(self.species) == 0:
            raise ValueError("Primero debe leer el archivo con read()")

        total_charges = []
        for i, species in enumerate(self.species):
            total = np.sum(species)
            total_charges.append(total)

        return total_charges

    def compute_charge_stats(self):
        """
        Calcula estadísticas de la distribución de carga

        Returns:
            Diccionario con estadísticas por especie
        """
        if self.species is None or len(self.species) == 0:
            raise ValueError("Primero debe leer el archivo con read()")

        stats = {}
        for i, species in enumerate(self.species):
            stats[f'species_{i}'] = {
                'total': np.sum(species),
                'mean': np.mean(species),
                'std': np.std(species),
                'min': np.min(species),
                'max': np.max(species)
            }

        # Estadísticas de carga neta
        if len(self.species) >= 2:
            net = self.compute_net_charge()
            stats['net_charge'] = {
                'total': np.sum(net),
                'mean': np.mean(net),
                'std': np.std(net),
                'min': np.min(net),
                'max': np.max(net)
            }

        return stats


class ChargeDistributionPlotter:
    """Genera gráficos de la distribución de carga"""

    def __init__(self, species, grid_size):
        """
        Inicializa el graficador

        Args:
            species: Lista de arrays 3D con la densidad de carga de cada especie
            grid_size: Tupla (nx, ny, nz)
        """
        self.species = species
        self.nx, self.ny, self.nz = grid_size
        self.num_species = len(species)

    def plot_plane(self, plane='xy', position=None, component='net',
                   species_index=0, show_colorbar=True, output=None):
        """
        Grafica la distribución de carga en un plano 2D

        Args:
            plane: 'xy', 'xz', o 'yz'
            position: Posición del plano (índice). Si es None, usa el centro
            component: 'net' (carga neta), 'species0', 'species1', 'total' (suma absoluta)
            species_index: Índice de la especie si component='speciesN'
            show_colorbar: Mostrar barra de colores
            output: Nombre del archivo de salida (None = mostrar en pantalla)
        """
        plane = plane.lower()

        # Determinar la posición del plano si no se especifica
        if position is None:
            if plane == 'xy':
                position = self.nz // 2
            elif plane == 'xz':
                position = self.ny // 2
            elif plane == 'yz':
                position = self.nx // 2

        # Seleccionar qué mostrar
        if component == 'net':
            if len(self.species) < 2:
                raise ValueError("Se necesitan al menos 2 especies para mostrar carga neta")
            charge = self.species[0] - self.species[1]
            title = f'Carga Neta (plano {plane.upper()}, posición {position})'
            clabel = 'ρ_net'
            cmap = 'RdBu_r'
        elif component.startswith('species'):
            idx = int(component.replace('species', ''))
            if idx >= len(self.species):
                raise ValueError(f"Especie {idx} no existe. Hay {len(self.species)} especies.")
            charge = self.species[idx]
            title = f'Densidad de Carga Especie {idx} (plano {plane.upper()}, posición {position})'
            clabel = f'ρ_{idx}'
            cmap = 'viridis'
        elif component == 'total':
            # Densidad iónica total = ρ₊ + ρ₋ (ambas son ≥ 0)
            charge = sum(s for s in self.species)
            title = f'Densidad Iónica Total (plano {plane.upper()}, posición {position})'
            clabel = 'ρ_total'
            cmap = 'plasma'
        else:
            raise ValueError("component debe ser 'net', 'species0', 'species1', ..., o 'total'")

        # Extraer el slice del plano
        if plane == 'xy':
            field = charge[:, :, position]
            xlabel, ylabel = 'X', 'Y'
            x = np.arange(self.nx)
            y = np.arange(self.ny)
        elif plane == 'xz':
            field = charge[:, position, :]
            xlabel, ylabel = 'X', 'Z'
            x = np.arange(self.nx)
            y = np.arange(self.nz)
        elif plane == 'yz':
            field = charge[position, :, :]
            xlabel, ylabel = 'Y', 'Z'
            x = np.arange(self.ny)
            y = np.arange(self.nz)
        else:
            raise ValueError("plane debe ser 'xy', 'xz', o 'yz'")

        # Crear meshgrid
        X, Y = np.meshgrid(x, y, indexing='ij')

        # Crear figura
        fig, ax = plt.subplots(figsize=(10, 8))

        # Plot de contorno
        im = ax.contourf(X, Y, field, levels=20, cmap=cmap)

        if show_colorbar:
            cbar = plt.colorbar(im, ax=ax, label=clabel)

        ax.set_xlabel(xlabel, fontsize=14)
        ax.set_ylabel(ylabel, fontsize=14)
        ax.set_title(title, fontsize=16)
        ax.set_aspect('equal')
        ax.grid(True, alpha=0.3)

        plt.tight_layout()

        if output:
            plt.savefig(output, dpi=300, bbox_inches='tight')
            print(f"Gráfico guardado en: {output}")
        else:
            plt.show()

        plt.close()

    def plot_plane_3d(self, plane='xy', position=None, component='net',
                      species_index=0, colormap='viridis', elevation=30,
                      azimuth=-60, output=None):
        """
        Grafica la distribución de carga como superficie 3D

        Args:
            plane: 'xy', 'xz', o 'yz'
            position: Posición del plano (índice). Si es None, usa el centro
            component: 'net', 'species0', 'species1', 'total'
            species_index: Índice de la especie si component='speciesN'
            colormap: Mapa de colores
            elevation: Ángulo de elevación de la vista (grados)
            azimuth: Ángulo azimutal de la vista (grados)
            output: Nombre del archivo de salida (None = mostrar en pantalla)
        """
        from mpl_toolkits.mplot3d import Axes3D

        plane = plane.lower()

        # Determinar la posición del plano si no se especifica
        if position is None:
            if plane == 'xy':
                position = self.nz // 2
            elif plane == 'xz':
                position = self.ny // 2
            elif plane == 'yz':
                position = self.nx // 2

        # Seleccionar qué mostrar
        if component == 'net':
            if len(self.species) < 2:
                raise ValueError("Se necesitan al menos 2 especies para mostrar carga neta")
            charge = self.species[0] - self.species[1]
            title = f'Carga Neta 3D (plano {plane.upper()}, posición {position})'
            clabel = 'ρ_net'
        elif component.startswith('species'):
            idx = int(component.replace('species', ''))
            if idx >= len(self.species):
                raise ValueError(f"Especie {idx} no existe. Hay {len(self.species)} especies.")
            charge = self.species[idx]
            title = f'Densidad Especie {idx} 3D (plano {plane.upper()}, posición {position})'
            clabel = f'ρ_{idx}'
        elif component == 'total':
            # Densidad iónica total = ρ₊ + ρ₋ (ambas son ≥ 0)
            charge = sum(s for s in self.species)
            title = f'Densidad Iónica Total 3D (plano {plane.upper()}, posición {position})'
            clabel = 'ρ_total'
        else:
            raise ValueError("component debe ser 'net', 'species0', 'species1', ..., o 'total'")

        # Extraer el slice del plano
        if plane == 'xy':
            field = charge[:, :, position]
            xlabel, ylabel = 'X', 'Y'
            x = np.arange(self.nx)
            y = np.arange(self.ny)
        elif plane == 'xz':
            field = charge[:, position, :]
            xlabel, ylabel = 'X', 'Z'
            x = np.arange(self.nx)
            y = np.arange(self.nz)
        elif plane == 'yz':
            field = charge[position, :, :]
            xlabel, ylabel = 'Y', 'Z'
            x = np.arange(self.ny)
            y = np.arange(self.nz)
        else:
            raise ValueError("plane debe ser 'xy', 'xz', o 'yz'")

        # Crear meshgrid
        X, Y = np.meshgrid(x, y, indexing='ij')

        # Crear figura 3D
        fig = plt.figure(figsize=(12, 9))
        ax = fig.add_subplot(111, projection='3d')

        # Graficar superficie
        surf = ax.plot_surface(X, Y, field, cmap=colormap,
                              linewidth=0, antialiased=True, alpha=0.9)

        # Configurar vista
        ax.view_init(elev=elevation, azim=azimuth)

        # Etiquetas y título
        ax.set_xlabel(xlabel, fontsize=12, labelpad=10)
        ax.set_ylabel(ylabel, fontsize=12, labelpad=10)
        ax.set_zlabel(clabel, fontsize=12, labelpad=10)
        ax.set_title(title, fontsize=14, pad=20)

        # Barra de colores
        cbar = fig.colorbar(surf, ax=ax, shrink=0.5, aspect=5, pad=0.1)
        cbar.set_label(clabel, fontsize=12)

        plt.tight_layout()

        if output:
            plt.savefig(output, dpi=300, bbox_inches='tight')
            print(f"Gráfico guardado en: {output}")
        else:
            plt.show()

        plt.close()

    def plot_line(self, start, end, num_points=100, component='net', output=None):
        """
        Grafica la distribución de carga a lo largo de una línea

        Args:
            start: Punto inicial (x, y, z) en coordenadas de la malla
            end: Punto final (x, y, z) en coordenadas de la malla
            num_points: Número de puntos a interpolar
            component: 'net', 'species0', 'species1', 'total', 'all'
            output: Nombre del archivo de salida (None = mostrar en pantalla)
        """
        from scipy.interpolate import RegularGridInterpolator

        # Crear puntos a lo largo de la línea
        points = np.array([np.linspace(start[i], end[i], num_points)
                          for i in range(3)]).T

        # Coordenadas de la malla
        x = np.arange(self.nx)
        y = np.arange(self.ny)
        z = np.arange(self.nz)

        # Interpolar cada especie
        species_lines = []
        for species in self.species:
            interp = RegularGridInterpolator((x, y, z), species)
            species_lines.append(interp(points))

        # Calcular distancia a lo largo de la línea
        distance = np.sqrt(np.sum((points - start)**2, axis=1))

        # Crear gráfico
        if component == 'all':
            fig, axes = plt.subplots(2, 1, figsize=(10, 10))

            # Especies individuales
            for i, species_line in enumerate(species_lines):
                axes[0].plot(distance, species_line, linewidth=2, label=f'Especie {i}')

            axes[0].set_xlabel('Distancia a lo largo de la línea', fontsize=14)
            axes[0].set_ylabel('Densidad de carga', fontsize=14)
            axes[0].set_title(f'Densidades por Especie\ndesde ({start[0]:.1f}, {start[1]:.1f}, {start[2]:.1f}) '
                             f'hasta ({end[0]:.1f}, {end[1]:.1f}, {end[2]:.1f})', fontsize=16)
            axes[0].legend(fontsize=12)
            axes[0].grid(True, alpha=0.3)

            # Carga neta
            if len(species_lines) >= 2:
                net_charge = species_lines[0] - species_lines[1]
                axes[1].plot(distance, net_charge, 'purple', linewidth=2)
                axes[1].set_xlabel('Distancia a lo largo de la línea', fontsize=14)
                axes[1].set_ylabel('Carga neta', fontsize=14)
                axes[1].set_title('Carga Neta', fontsize=16)
                axes[1].grid(True, alpha=0.3)

        else:
            fig, ax = plt.subplots(1, 1, figsize=(10, 6))

            if component == 'net':
                if len(species_lines) < 2:
                    raise ValueError("Se necesitan al menos 2 especies para carga neta")
                data = species_lines[0] - species_lines[1]
                ylabel = 'Carga neta'
                title = 'Carga Neta'
                color = 'purple'
            elif component.startswith('species'):
                idx = int(component.replace('species', ''))
                if idx >= len(species_lines):
                    raise ValueError(f"Especie {idx} no existe")
                data = species_lines[idx]
                ylabel = f'Densidad especie {idx}'
                title = f'Densidad de Carga Especie {idx}'
                color = f'C{idx}'
            elif component == 'total':
                # Densidad iónica total = ρ₊ + ρ₋ (ambas son ≥ 0)
                data = sum(s for s in species_lines)
                ylabel = 'Densidad iónica total'
                title = 'Densidad Iónica Total'
                color = 'black'
            else:
                raise ValueError("component debe ser 'net', 'species0', 'species1', 'total', o 'all'")

            ax.plot(distance, data, color=color, linewidth=2)
            ax.set_xlabel('Distancia a lo largo de la línea', fontsize=14)
            ax.set_ylabel(ylabel, fontsize=14)
            ax.set_title(f'{title}\ndesde ({start[0]:.1f}, {start[1]:.1f}, {start[2]:.1f}) '
                        f'hasta ({end[0]:.1f}, {end[1]:.1f}, {end[2]:.1f})', fontsize=16)
            ax.grid(True, alpha=0.3)

        plt.tight_layout()

        if output:
            plt.savefig(output, dpi=300, bbox_inches='tight')
            print(f"Gráfico guardado en: {output}")
        else:
            plt.show()

        plt.close()


def print_charge_statistics(stats, grid_size):
    """
    Imprime estadísticas de la distribución de carga

    Args:
        stats: Diccionario con estadísticas
        grid_size: Tupla (nx, ny, nz)
    """
    nx, ny, nz = grid_size
    volume = nx * ny * nz

    print("\n" + "="*70)
    print("ESTADÍSTICAS DE DISTRIBUCIÓN DE CARGA")
    print("="*70)
    print(f"Tamaño de malla: {nx} × {ny} × {nz} = {volume} puntos\n")

    for key, values in stats.items():
        if key == 'species_0':
            print("Densidad de Carga Positiva ρ₊ (Cationes, carga +1):")
        elif key == 'species_1':
            print("Densidad de Carga Negativa ρ₋ (Aniones, carga -1):")
        elif key.startswith('species'):
            species_num = key.split('_')[1]
            print(f"Especie {species_num}:")
        else:
            print("Carga Neta (ρ₊ - ρ₋):")

        print(f"  Carga/Densidad total: {values['total']:15.6e}")
        print(f"  Promedio:             {values['mean']:15.6e}")
        print(f"  Desv. std:            {values['std']:15.6e}")
        print(f"  Mínimo:               {values['min']:15.6e}")

        # Advertencia si hay valores negativos en densidades
        if key.startswith('species') and values['min'] < 0:
            print(f"  Máximo:               {values['max']:15.6e}")
            print("  ⚠️  ADVERTENCIA: Valores negativos detectados!")
            print("      Las densidades deben ser ≥ 0")
        else:
            print(f"  Máximo:               {values['max']:15.6e}")
        print()

    # Verificar electroneutralidad
    if 'net_charge' in stats:
        net_total = stats['net_charge']['total']
        if abs(net_total) < 1e-10:
            print("✓ Sistema electroneutro (carga neta total ≈ 0)")
        else:
            print(f"⚠️  Carga neta total = {net_total:.6e}")
            print("   (Verificar si esto es esperado)")
        print()

    print("="*70 + "\n")


def main():
    parser = argparse.ArgumentParser(
        description='Graficar distribución de carga desde archivos qsi de Ludwig',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Ejemplos de uso:

PLANOS 2D:
1. Carga neta en plano XY central:
   ./plot_charge_distribution.py -f qsi-000010250.001-001 -s 32 32 32 -m plane -p xy -c net

2. Especie 0 en plano XZ:
   ./plot_charge_distribution.py -f qsi-000010250.001-001 -s 32 32 32 -m plane -p xz -c species0

3. Densidad total en plano YZ:
   ./plot_charge_distribution.py -f qsi-000010250.001-001 -s 32 32 32 -m plane -p yz -c total

SUPERFICIES 3D:
4. Superficie 3D de carga neta:
   ./plot_charge_distribution.py -f qsi-000010250.001-001 -s 32 32 32 -m plane3d -p xy -c net

5. Superficie 3D de especie 1 con colormap personalizado:
   ./plot_charge_distribution.py -f qsi-000010250.001-001 -s 32 32 32 -m plane3d -p xz -c species1 --colormap plasma

LÍNEAS 1D:
6. Todas las especies a lo largo del eje X:
   ./plot_charge_distribution.py -f qsi-000010250.001-001 -s 32 32 32 -m line --start 0 16 16 --end 31 16 16 -c all

7. Solo carga neta a lo largo de diagonal:
   ./plot_charge_distribution.py -f qsi-000010250.001-001 -s 32 32 32 -m line --start 0 0 0 --end 31 31 31 -c net

ESTADÍSTICAS:
8. Calcular carga total sin graficar:
   ./plot_charge_distribution.py -f qsi-000010250.001-001 -s 32 32 32 -m stats
        """
    )

    # Argumentos obligatorios
    parser.add_argument('-f', '--file', required=True,
                       help='Archivo qsi a leer')
    parser.add_argument('-s', '--size', nargs=3, type=int, required=True,
                       metavar=('NX', 'NY', 'NZ'),
                       help='Tamaño de la malla (nx ny nz)')
    parser.add_argument('-m', '--mode', choices=['plane', 'line', 'plane3d', 'stats'], required=True,
                       help='Modo: plane (2D), line (1D), plane3d (superficie 3D), stats (solo estadísticas)')

    # Argumentos para modos gráficos
    parser.add_argument('-p', '--plane', choices=['xy', 'xz', 'yz'],
                       help='Plano a graficar (obligatorio para mode=plane/plane3d)')
    parser.add_argument('-pos', '--position', type=int,
                       help='Posición del plano (índice). Por defecto: centro')

    # Argumentos para modo línea
    parser.add_argument('--start', nargs=3, type=float, metavar=('X', 'Y', 'Z'),
                       help='Punto inicial de la línea (obligatorio para mode=line)')
    parser.add_argument('--end', nargs=3, type=float, metavar=('X', 'Y', 'Z'),
                       help='Punto final de la línea (obligatorio para mode=line)')
    parser.add_argument('--num-points', type=int, default=100,
                       help='Número de puntos a lo largo de la línea (por defecto: 100)')

    # Argumentos compartidos
    parser.add_argument('-c', '--component',
                       choices=['net', 'species0', 'species1', 'total', 'all'],
                       default='net',
                       help='Componente a graficar: net (carga neta), species0/1 (especie individual), '
                            'total (densidad total), all (todas en modo line). Default: net')

    # Argumentos para modo 3D
    parser.add_argument('--colormap', default='RdBu_r',
                       help='Mapa de colores para modo plane3d (default: RdBu_r)')
    parser.add_argument('--elevation', type=float, default=30,
                       help='Ángulo de elevación de la vista 3D (default: 30)')
    parser.add_argument('--azimuth', type=float, default=-60,
                       help='Ángulo azimutal de la vista 3D (default: -60)')

    # Argumentos generales
    parser.add_argument('--num-species', type=int, default=2,
                       help='Número de especies iónicas (default: 2)')
    parser.add_argument('-o', '--output',
                       help='Archivo de salida (por defecto: mostrar en pantalla)')

    args = parser.parse_args()

    # Validar argumentos según el modo
    if args.mode in ['plane', 'plane3d'] and args.plane is None:
        parser.error(f"mode={args.mode} requiere especificar --plane")
    if args.mode == 'line' and (args.start is None or args.end is None):
        parser.error("mode=line requiere especificar --start y --end")

    # Verificar que el archivo existe
    if not os.path.exists(args.file):
        print(f"Error: No se encuentra el archivo {args.file}")
        sys.exit(1)

    # Leer archivo qsi
    print(f"Leyendo archivo {args.file}...")
    reader = QsiReader(args.file, tuple(args.size), args.num_species)
    species = reader.read()
    print(f"Archivo leído exitosamente. Malla: {args.size[0]}x{args.size[1]}x{args.size[2]}")
    print(f"Número de especies: {len(species)}")

    # Calcular y mostrar estadísticas
    stats = reader.compute_charge_stats()
    print_charge_statistics(stats, tuple(args.size))

    # Si solo se requieren estadísticas, terminar aquí
    if args.mode == 'stats':
        print("¡Listo!")
        return

    # Crear graficador
    plotter = ChargeDistributionPlotter(species, tuple(args.size))

    # Generar gráfico según el modo
    if args.mode == 'plane':
        print(f"Graficando plano {args.plane.upper()}...")
        plotter.plot_plane(
            plane=args.plane,
            position=args.position,
            component=args.component,
            output=args.output
        )
    elif args.mode == 'plane3d':
        print(f"Graficando superficie 3D en plano {args.plane.upper()}...")
        plotter.plot_plane_3d(
            plane=args.plane,
            position=args.position,
            component=args.component,
            colormap=args.colormap,
            elevation=args.elevation,
            azimuth=args.azimuth,
            output=args.output
        )
    else:  # mode == 'line'
        print(f"Graficando línea desde {args.start} hasta {args.end}...")
        plotter.plot_line(
            start=args.start,
            end=args.end,
            num_points=args.num_points,
            component=args.component,
            output=args.output
        )

    print("¡Listo!")


if __name__ == '__main__':
    main()
