#!/usr/bin/env python3
"""
Script para graficar el campo eléctrico desde archivos psi de Ludwig.

El campo eléctrico se calcula como E = -∇ψ, donde ψ es el potencial eléctrico.

Modos de visualización:
1. Plano 2D: Muestra el campo eléctrico en un plano específico (XY, XZ, o YZ)
2. Línea 1D: Muestra el campo eléctrico a lo largo de una línea en el espacio
3. Plano 3D: Muestra superficie 3D del potencial o campo eléctrico en un plano

Autor: Script generado para análisis de simulaciones Ludwig
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
import os
from matplotlib.colors import Normalize
from matplotlib import cm

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


class ElectricFieldPlotter:
    """Genera gráficos del campo eléctrico"""

    def __init__(self, psi, Ex, Ey, Ez, grid_size, external_field=None):
        """
        Inicializa el graficador

        Args:
            psi: Potencial eléctrico (array 3D)
            Ex, Ey, Ez: Componentes del campo eléctrico (arrays 3D)
            grid_size: Tupla (nx, ny, nz)
            external_field: Campo externo constante (Ex_ext, Ey_ext, Ez_ext) o None
        """
        self.psi = psi
        self.nx, self.ny, self.nz = grid_size

        # Si hay campo externo, restar del campo total
        if external_field is not None:
            Ex_ext, Ey_ext, Ez_ext = external_field
            self.Ex = Ex - Ex_ext
            self.Ey = Ey - Ey_ext
            self.Ez = Ez - Ez_ext
            self.external_field = external_field
            print(f"Campo externo restado: E_ext = ({Ex_ext:.6e}, {Ey_ext:.6e}, {Ez_ext:.6e})")
        else:
            self.Ex = Ex
            self.Ey = Ey
            self.Ez = Ez
            self.external_field = None

    def plot_plane(self, plane='xy', position=None, component='magnitude',
                   show_vectors=True, vector_stride=2, output=None):
        """
        Grafica el campo eléctrico en un plano 2D

        Args:
            plane: 'xy', 'xz', o 'yz'
            position: Posición del plano (índice). Si es None, usa el centro
            component: 'magnitude', 'x', 'y', 'z', o 'psi'
            show_vectors: Si True, superpone vectores del campo
            vector_stride: Espaciado entre vectores (1 = todos los puntos)
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

        # Extraer el slice del plano
        if plane == 'xy':
            E1 = self.Ex[:, :, position]
            E2 = self.Ey[:, :, position]
            E3 = self.Ez[:, :, position]
            psi_slice = self.psi[:, :, position]
            xlabel, ylabel = 'X', 'Y'
            x = np.arange(self.nx)
            y = np.arange(self.ny)
        elif plane == 'xz':
            E1 = self.Ex[:, position, :]
            E2 = self.Ez[:, position, :]
            E3 = self.Ey[:, position, :]
            psi_slice = self.psi[:, position, :]
            xlabel, ylabel = 'X', 'Z'
            x = np.arange(self.nx)
            y = np.arange(self.nz)
        elif plane == 'yz':
            E1 = self.Ey[position, :, :]
            E2 = self.Ez[position, :, :]
            E3 = self.Ex[position, :, :]
            psi_slice = self.psi[position, :, :]
            xlabel, ylabel = 'Y', 'Z'
            x = np.arange(self.ny)
            y = np.arange(self.nz)
        else:
            raise ValueError("plane debe ser 'xy', 'xz', o 'yz'")

        # Seleccionar qué mostrar
        if component == 'magnitude':
            if plane == 'xy':
                field = np.sqrt(self.Ex[:, :, position]**2 +
                               self.Ey[:, :, position]**2 +
                               self.Ez[:, :, position]**2)
            elif plane == 'xz':
                field = np.sqrt(self.Ex[:, position, :]**2 +
                               self.Ey[:, position, :]**2 +
                               self.Ez[:, position, :]**2)
            else:
                field = np.sqrt(self.Ex[position, :, :]**2 +
                               self.Ey[position, :, :]**2 +
                               self.Ez[position, :, :]**2)
            title = f'Magnitud del Campo Eléctrico (plano {plane.upper()}, posición {position})'
            clabel = '|E|'
        elif component == 'x':
            field = E1 if plane in ['xy', 'xz'] else E3
            title = f'Componente Ex del Campo Eléctrico (plano {plane.upper()}, posición {position})'
            clabel = 'Ex'
        elif component == 'y':
            field = E2 if plane == 'xy' else E1 if plane == 'yz' else E3
            title = f'Componente Ey del Campo Eléctrico (plano {plane.upper()}, posición {position})'
            clabel = 'Ey'
        elif component == 'z':
            field = E2 if plane in ['xz', 'yz'] else E3
            title = f'Componente Ez del Campo Eléctrico (plano {plane.upper()}, posición {position})'
            clabel = 'Ez'
        elif component == 'psi':
            field = psi_slice
            title = f'Potencial Eléctrico ψ (plano {plane.upper()}, posición {position})'
            clabel = 'ψ'
        else:
            raise ValueError("component debe ser 'magnitude', 'x', 'y', 'z', o 'psi'")

        # Crear figura
        fig, ax = plt.subplots(figsize=(10, 8))

        # Crear meshgrid para el contour plot
        X, Y = np.meshgrid(x, y, indexing='ij')

        # Plot de contorno del campo
        im = ax.contourf(X, Y, field, levels=20, cmap='RdBu_r')
        cbar = plt.colorbar(im, ax=ax, label=clabel)

        # Superponer vectores si se solicita
        if show_vectors and component != 'psi':
            # Submuestreo para no saturar el gráfico
            X_sub = X[::vector_stride, ::vector_stride]
            Y_sub = Y[::vector_stride, ::vector_stride]
            E1_sub = E1[::vector_stride, ::vector_stride]
            E2_sub = E2[::vector_stride, ::vector_stride]

            ax.quiver(X_sub, Y_sub, E1_sub, E2_sub,
                     color='black', alpha=0.6, scale=None, width=0.003)

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

    def plot_line(self, start, end, num_points=100, components='all', output=None):
        """
        Grafica el campo eléctrico a lo largo de una línea en el espacio

        Args:
            start: Punto inicial (x, y, z) en coordenadas de la malla
            end: Punto final (x, y, z) en coordenadas de la malla
            num_points: Número de puntos a interpolar
            components: 'all', 'x', 'y', 'z', 'magnitude', o 'psi'
            output: Nombre del archivo de salida (None = mostrar en pantalla)
        """
        # Crear puntos a lo largo de la línea
        t = np.linspace(0, 1, num_points)
        points = np.array([np.linspace(start[i], end[i], num_points)
                          for i in range(3)]).T

        # Interpolar los valores en estos puntos
        from scipy.interpolate import RegularGridInterpolator

        x = np.arange(self.nx)
        y = np.arange(self.ny)
        z = np.arange(self.nz)

        interp_Ex = RegularGridInterpolator((x, y, z), self.Ex)
        interp_Ey = RegularGridInterpolator((x, y, z), self.Ey)
        interp_Ez = RegularGridInterpolator((x, y, z), self.Ez)
        interp_psi = RegularGridInterpolator((x, y, z), self.psi)

        Ex_line = interp_Ex(points)
        Ey_line = interp_Ey(points)
        Ez_line = interp_Ez(points)
        psi_line = interp_psi(points)

        # Calcular magnitud
        E_mag = np.sqrt(Ex_line**2 + Ey_line**2 + Ez_line**2)

        # Calcular distancia a lo largo de la línea
        distance = np.sqrt(np.sum((points - start)**2, axis=1))

        # Determinar qué componentes graficar
        if components == 'all':
            # Gráfico original con todas las componentes
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))

            # Componentes del campo eléctrico
            ax1.plot(distance, Ex_line, 'r-', label='Ex', linewidth=2)
            ax1.plot(distance, Ey_line, 'g-', label='Ey', linewidth=2)
            ax1.plot(distance, Ez_line, 'b-', label='Ez', linewidth=2)
            ax1.plot(distance, E_mag, 'k--', label='|E|', linewidth=2)
            ax1.set_xlabel('Distancia a lo largo de la línea', fontsize=14)
            ax1.set_ylabel('Campo Eléctrico', fontsize=14)
            ax1.set_title(f'Campo Eléctrico desde ({start[0]:.1f}, {start[1]:.1f}, {start[2]:.1f}) '
                         f'hasta ({end[0]:.1f}, {end[1]:.1f}, {end[2]:.1f})', fontsize=16)
            ax1.legend(fontsize=12)
            ax1.grid(True, alpha=0.3)

            # Potencial eléctrico
            ax2.plot(distance, psi_line, 'purple', linewidth=2)
            ax2.set_xlabel('Distancia a lo largo de la línea', fontsize=14)
            ax2.set_ylabel('Potencial ψ', fontsize=14)
            ax2.set_title('Potencial Eléctrico', fontsize=16)
            ax2.grid(True, alpha=0.3)

        else:
            # Gráfico simple con componente seleccionada
            fig, ax = plt.subplots(1, 1, figsize=(10, 6))

            if components == 'x':
                ax.plot(distance, Ex_line, 'r-', linewidth=2)
                ylabel = 'Ex'
                title = f'Componente Ex del Campo Eléctrico'
            elif components == 'y':
                ax.plot(distance, Ey_line, 'g-', linewidth=2)
                ylabel = 'Ey'
                title = f'Componente Ey del Campo Eléctrico'
            elif components == 'z':
                ax.plot(distance, Ez_line, 'b-', linewidth=2)
                ylabel = 'Ez'
                title = f'Componente Ez del Campo Eléctrico'
            elif components == 'magnitude':
                ax.plot(distance, E_mag, 'k-', linewidth=2)
                ylabel = '|E|'
                title = f'Magnitud del Campo Eléctrico'
            elif components == 'psi':
                ax.plot(distance, psi_line, 'purple', linewidth=2)
                ylabel = 'ψ'
                title = f'Potencial Eléctrico'
            else:
                raise ValueError("components debe ser 'all', 'x', 'y', 'z', 'magnitude', o 'psi'")

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

    def plot_plane_3d(self, plane='xy', position=None, component='psi',
                      colormap='viridis', elevation=30, azimuth=-60, output=None):
        """
        Grafica el campo eléctrico o potencial en 3D sobre un plano

        Args:
            plane: 'xy', 'xz', o 'yz'
            position: Posición del plano (índice). Si es None, usa el centro
            component: 'magnitude', 'x', 'y', 'z', o 'psi'
            colormap: Mapa de colores ('viridis', 'plasma', 'RdBu_r', etc.)
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

        # Extraer el slice del plano y preparar datos
        if plane == 'xy':
            if component == 'magnitude':
                field = np.sqrt(self.Ex[:, :, position]**2 +
                               self.Ey[:, :, position]**2 +
                               self.Ez[:, :, position]**2)
                clabel = '|E|'
                title = f'Magnitud del Campo Eléctrico 3D (plano XY, z={position})'
            elif component == 'x':
                field = self.Ex[:, :, position]
                clabel = 'Ex'
                title = f'Componente Ex del Campo Eléctrico 3D (plano XY, z={position})'
            elif component == 'y':
                field = self.Ey[:, :, position]
                clabel = 'Ey'
                title = f'Componente Ey del Campo Eléctrico 3D (plano XY, z={position})'
            elif component == 'z':
                field = self.Ez[:, :, position]
                clabel = 'Ez'
                title = f'Componente Ez del Campo Eléctrico 3D (plano XY, z={position})'
            elif component == 'psi':
                field = self.psi[:, :, position]
                clabel = 'ψ'
                title = f'Potencial Eléctrico ψ 3D (plano XY, z={position})'
            else:
                raise ValueError("component debe ser 'magnitude', 'x', 'y', 'z', o 'psi'")

            xlabel, ylabel = 'X', 'Y'
            x = np.arange(self.nx)
            y = np.arange(self.ny)

        elif plane == 'xz':
            if component == 'magnitude':
                field = np.sqrt(self.Ex[:, position, :]**2 +
                               self.Ey[:, position, :]**2 +
                               self.Ez[:, position, :]**2)
                clabel = '|E|'
                title = f'Magnitud del Campo Eléctrico 3D (plano XZ, y={position})'
            elif component == 'x':
                field = self.Ex[:, position, :]
                clabel = 'Ex'
                title = f'Componente Ex del Campo Eléctrico 3D (plano XZ, y={position})'
            elif component == 'y':
                field = self.Ey[:, position, :]
                clabel = 'Ey'
                title = f'Componente Ey del Campo Eléctrico 3D (plano XZ, y={position})'
            elif component == 'z':
                field = self.Ez[:, position, :]
                clabel = 'Ez'
                title = f'Componente Ez del Campo Eléctrico 3D (plano XZ, y={position})'
            elif component == 'psi':
                field = self.psi[:, position, :]
                clabel = 'ψ'
                title = f'Potencial Eléctrico ψ 3D (plano XZ, y={position})'
            else:
                raise ValueError("component debe ser 'magnitude', 'x', 'y', 'z', o 'psi'")

            xlabel, ylabel = 'X', 'Z'
            x = np.arange(self.nx)
            y = np.arange(self.nz)

        elif plane == 'yz':
            if component == 'magnitude':
                field = np.sqrt(self.Ex[position, :, :]**2 +
                               self.Ey[position, :, :]**2 +
                               self.Ez[position, :, :]**2)
                clabel = '|E|'
                title = f'Magnitud del Campo Eléctrico 3D (plano YZ, x={position})'
            elif component == 'x':
                field = self.Ex[position, :, :]
                clabel = 'Ex'
                title = f'Componente Ex del Campo Eléctrico 3D (plano YZ, x={position})'
            elif component == 'y':
                field = self.Ey[position, :, :]
                clabel = 'Ey'
                title = f'Componente Ey del Campo Eléctrico 3D (plano YZ, x={position})'
            elif component == 'z':
                field = self.Ez[position, :, :]
                clabel = 'Ez'
                title = f'Componente Ez del Campo Eléctrico 3D (plano YZ, x={position})'
            elif component == 'psi':
                field = self.psi[position, :, :]
                clabel = 'ψ'
                title = f'Potencial Eléctrico ψ 3D (plano YZ, x={position})'
            else:
                raise ValueError("component debe ser 'magnitude', 'x', 'y', 'z', o 'psi'")

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
                              linewidth=0, antialiased=True,
                              alpha=0.9)

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

        # Ajustar diseño
        plt.tight_layout()

        if output:
            plt.savefig(output, dpi=300, bbox_inches='tight')
            print(f"Gráfico guardado en: {output}")
        else:
            plt.show()

        plt.close()


def main():
    parser = argparse.ArgumentParser(
        description='Graficar campo eléctrico desde archivos psi de Ludwig',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Ejemplos de uso:

PLANOS 2D:
1. Graficar campo eléctrico en plano XY central:
   ./plot_electric_field.py -f psi-000010050.001-001 -s 32 32 32 -m plane -p xy

2. Graficar componente Ez en plano XZ en z=10 con vectores:
   ./plot_electric_field.py -f psi-000010050.001-001 -s 32 32 32 -m plane -p xz -pos 10 -c z -v

3. Graficar potencial eléctrico en plano YZ:
   ./plot_electric_field.py -f psi-000010050.001-001 -s 32 32 32 -m plane -p yz -c psi

SUPERFICIES 3D:
4. Superficie 3D del potencial en plano XY:
   ./plot_electric_field.py -f psi-000010050.001-001 -s 32 32 32 -m plane3d -p xy -c psi

5. Superficie 3D de la magnitud del campo con mapa de colores personalizado:
   ./plot_electric_field.py -f psi-000010050.001-001 -s 32 32 32 -m plane3d -p xz -c magnitude --colormap plasma

6. Superficie 3D con vista personalizada:
   ./plot_electric_field.py -f psi-000010050.001-001 -s 32 32 32 -m plane3d -p xy -c psi --elevation 45 --azimuth -30

LÍNEAS 1D:
7. Todas las componentes a lo largo de una diagonal:
   ./plot_electric_field.py -f psi-000010050.001-001 -s 32 32 32 -m line --start 0 0 0 --end 31 31 31 -c all

8. Solo componente Ez a lo largo del eje X:
   ./plot_electric_field.py -f psi-000010050.001-001 -s 32 32 32 -m line --start 0 16 16 --end 31 16 16 -c z

9. Magnitud del campo a lo largo del eje Z:
   ./plot_electric_field.py -f psi-000010050.001-001 -s 32 32 32 -m line --start 16 16 0 --end 16 16 31 -c magnitude

RESTAR CAMPO EXTERNO:
10. Restar campo externo uniforme en dirección Z (ej: E_ext = 0.001 en Z):
   ./plot_electric_field.py -f psi-000010050.001-001 -s 32 32 32 -m plane -p xy --external-field 0 0 0.001

11. Restar campo externo en dirección X:
   ./plot_electric_field.py -f psi-000010050.001-001 -s 32 32 32 -m plane3d -p xz -c magnitude --external-field 0.005 0 0

12. Visualizar solo el campo generado por cargas (restar campo aplicado):
   ./plot_electric_field.py -f psi-000010050.001-001 -s 32 32 32 -m line --start 0 16 16 --end 31 16 16 -c all --external-field 0 0 0.001
        """
    )

    # Argumentos obligatorios
    parser.add_argument('-f', '--file', required=True,
                       help='Archivo psi a leer')
    parser.add_argument('-s', '--size', nargs=3, type=int, required=True,
                       metavar=('NX', 'NY', 'NZ'),
                       help='Tamaño de la malla (nx ny nz)')
    parser.add_argument('-m', '--mode', choices=['plane', 'line', 'plane3d'], required=True,
                       help='Modo de graficación: plane (plano 2D), line (línea 1D), o plane3d (superficie 3D)')

    # Argumentos para modo plano
    parser.add_argument('-p', '--plane', choices=['xy', 'xz', 'yz'],
                       help='Plano a graficar (obligatorio para mode=plane)')
    parser.add_argument('-pos', '--position', type=int,
                       help='Posición del plano (índice). Por defecto: centro')
    parser.add_argument('-v', '--vectors', action='store_true',
                       help='Mostrar vectores del campo (solo para mode=plane)')
    parser.add_argument('--vector-stride', type=int, default=2,
                       help='Espaciado entre vectores (por defecto: 2)')

    # Argumentos para modo línea
    parser.add_argument('--start', nargs=3, type=float, metavar=('X', 'Y', 'Z'),
                       help='Punto inicial de la línea (obligatorio para mode=line)')
    parser.add_argument('--end', nargs=3, type=float, metavar=('X', 'Y', 'Z'),
                       help='Punto final de la línea (obligatorio para mode=line)')
    parser.add_argument('--num-points', type=int, default=100,
                       help='Número de puntos a lo largo de la línea (por defecto: 100)')

    # Argumentos compartidos
    parser.add_argument('-c', '--component',
                       choices=['magnitude', 'x', 'y', 'z', 'psi', 'all'],
                       default='magnitude',
                       help='Componente a graficar. Plane/Plane3D: magnitude/x/y/z/psi (default: magnitude). '
                            'Line: all/magnitude/x/y/z/psi (default: magnitude). '
                            'Usar "all" en modo line muestra todas las componentes juntas')

    # Argumentos para modo 3D
    parser.add_argument('--colormap', default='viridis',
                       help='Mapa de colores para modo plane3d (default: viridis). '
                            'Opciones: viridis, plasma, inferno, magma, RdBu_r, coolwarm, seismic')
    parser.add_argument('--elevation', type=float, default=30,
                       help='Ángulo de elevación de la vista 3D en grados (default: 30)')
    parser.add_argument('--azimuth', type=float, default=-60,
                       help='Ángulo azimutal de la vista 3D en grados (default: -60)')

    # Argumentos para campo externo
    parser.add_argument('--external-field', nargs=3, type=float, metavar=('Ex', 'Ey', 'Ez'),
                       help='Campo eléctrico externo constante a restar (Ex Ey Ez). '
                            'Útil para visualizar solo el campo generado por cargas, '
                            'sin el campo aplicado externamente.')

    # Argumentos generales
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

    # Leer archivo psi
    print(f"Leyendo archivo {args.file}...")
    reader = PsiReader(args.file, tuple(args.size))
    psi = reader.read()
    print(f"Archivo leído exitosamente. Malla: {args.size[0]}x{args.size[1]}x{args.size[2]}")

    # Calcular campo eléctrico
    print("Calculando campo eléctrico E = -∇ψ...")
    Ex, Ey, Ez = reader.compute_electric_field()
    print("Campo eléctrico calculado.")

    # Preparar campo externo si se especificó
    external_field = None
    if args.external_field is not None:
        external_field = tuple(args.external_field)
        print(f"\nRestando campo externo aplicado:")
        print(f"  E_ext = ({external_field[0]:.6e}, {external_field[1]:.6e}, {external_field[2]:.6e})")
        print(f"  Se graficará: E_total - E_ext (campo generado por cargas)\n")

    # Crear graficador
    plotter = ElectricFieldPlotter(psi, Ex, Ey, Ez, tuple(args.size), external_field=external_field)

    # Generar gráfico según el modo
    if args.mode == 'plane':
        print(f"Graficando plano {args.plane.upper()}...")
        plotter.plot_plane(
            plane=args.plane,
            position=args.position,
            component=args.component,
            show_vectors=args.vectors,
            vector_stride=args.vector_stride,
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
            components=args.component,
            output=args.output
        )

    print("¡Listo!")


if __name__ == '__main__':
    main()
