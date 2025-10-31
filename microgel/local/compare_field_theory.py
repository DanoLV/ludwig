#!/usr/bin/env python3
"""
Script para comparar el campo eléctrico simulado con el campo teórico
de una carga puntual.

Para una carga puntual q en el vacío (o medio uniforme), el campo eléctrico es:
E = (q / (4πε)) * (r / |r|³)

donde r es el vector desde la carga al punto de observación.

En Ludwig, el sistema está discretizado en una malla y se usa el sistema de
unidades de lattice Boltzmann.
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
import os

class PsiReader:
    """Lee y procesa archivos psi de Ludwig"""

    def __init__(self, filename, grid_size):
        self.filename = filename
        self.nx, self.ny, self.nz = grid_size
        self.psi = None

    def read(self):
        """Lee el archivo psi y lo organiza en una malla 3D"""
        data = np.loadtxt(self.filename)

        if len(data) != self.nx * self.ny * self.nz:
            raise ValueError(f"El archivo tiene {len(data)} valores, "
                           f"pero se esperaban {self.nx * self.ny * self.nz}")

        self.psi = data.reshape((self.nx, self.ny, self.nz))
        return self.psi

    def compute_electric_field(self):
        """Calcula E = -∇ψ"""
        if self.psi is None:
            raise ValueError("Primero debe leer el archivo con read()")

        grad_psi = np.gradient(self.psi)
        Ex = -grad_psi[0]
        Ey = -grad_psi[1]
        Ez = -grad_psi[2]

        return Ex, Ey, Ez


def theoretical_field(charge_pos, grid_size, q=1.0, epsilon=1.0, scale_factor=None):
    """
    Calcula el campo eléctrico teórico de una carga puntual

    Args:
        charge_pos: Posición de la carga (x, y, z)
        grid_size: Tupla (nx, ny, nz)
        q: Carga (default: 1.0)
        epsilon: Permitividad (default: 1.0)
        scale_factor: Factor de escala manual. Si es None, usa q/(4πε)

    Returns:
        Ex_theory, Ey_theory, Ez_theory: Componentes del campo teórico
    """
    nx, ny, nz = grid_size
    xc, yc, zc = charge_pos

    # Crear malla de coordenadas centradas en nodos (0.5, 1.5, 2.5, ...)
    x = np.arange(0.5, nx)
    y = np.arange(0.5, ny)
    z = np.arange(0.5, nz)

    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

    # Vector distancia desde la carga
    dx = X - xc
    dy = Y - yc
    dz = Z - zc

    # Distancia al cuadrado
    r_squared = dx**2 + dy**2 + dz**2

    # Evitar división por cero en la posición de la carga
    r_squared = np.where(r_squared < 0.01, 0.01, r_squared)

    r_cubed = r_squared * np.sqrt(r_squared)

    # Campo eléctrico teórico: E = prefactor * (r / r³)
    if scale_factor is None:
        # Usar factor teórico estándar
        prefactor = q / (4.0 * np.pi * epsilon)
    else:
        # Usar factor de escala proporcionado
        prefactor = scale_factor

    Ex_theory = prefactor * dx / r_cubed
    Ey_theory = prefactor * dy / r_cubed
    Ez_theory = prefactor * dz / r_cubed

    return Ex_theory, Ey_theory, Ez_theory


def fit_scale_factor(Ex_sim, Ey_sim, Ez_sim, charge_pos, grid_size, q=1.0):
    """
    Encuentra el factor de escala óptimo comparando el campo simulado con el patrón teórico

    Args:
        Ex_sim, Ey_sim, Ez_sim: Componentes del campo simulado
        charge_pos: Posición de la carga
        grid_size: Tupla (nx, ny, nz)
        q: Carga

    Returns:
        scale_factor: Factor de escala óptimo
    """
    # Calcular campo teórico con factor unitario
    Ex_unit, Ey_unit, Ez_unit = theoretical_field(charge_pos, grid_size, q=1.0, epsilon=1.0, scale_factor=1.0)

    # Crear máscara para excluir región cercana a la carga (donde la discretización afecta)
    nx, ny, nz = grid_size
    xc, yc, zc = charge_pos
    x = np.arange(0.5, nx)
    y = np.arange(0.5, ny)
    z = np.arange(0.5, nz)
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

    r_squared = (X - xc)**2 + (Y - yc)**2 + (Z - zc)**2
    mask = r_squared > 4.0  # Excluir región dentro de 2 unidades de lattice

    # Calcular factor de escala mediante ajuste de mínimos cuadrados
    # scale_factor = sum(E_sim * E_unit) / sum(E_unit^2)
    numerator = (np.sum(Ex_sim[mask] * Ex_unit[mask]) +
                 np.sum(Ey_sim[mask] * Ey_unit[mask]) +
                 np.sum(Ez_sim[mask] * Ez_unit[mask]))

    denominator = (np.sum(Ex_unit[mask]**2) +
                   np.sum(Ey_unit[mask]**2) +
                   np.sum(Ez_unit[mask]**2))

    scale_factor = numerator / denominator

    return scale_factor


def compare_fields_line(Ex_sim, Ey_sim, Ez_sim, Ex_theory, Ey_theory, Ez_theory,
                        start, end, grid_size, num_points=100, output=None):
    """
    Compara los campos simulado y teórico a lo largo de una línea
    """
    from scipy.interpolate import RegularGridInterpolator

    nx, ny, nz = grid_size

    # Coordenadas de la malla
    x = np.arange(0.5, nx)
    y = np.arange(0.5, ny)
    z = np.arange(0.5, nz)

    # Crear interpoladores
    interp_Ex_sim = RegularGridInterpolator((x, y, z), Ex_sim, bounds_error=False, fill_value=0)
    interp_Ey_sim = RegularGridInterpolator((x, y, z), Ey_sim, bounds_error=False, fill_value=0)
    interp_Ez_sim = RegularGridInterpolator((x, y, z), Ez_sim, bounds_error=False, fill_value=0)

    interp_Ex_theory = RegularGridInterpolator((x, y, z), Ex_theory, bounds_error=False, fill_value=0)
    interp_Ey_theory = RegularGridInterpolator((x, y, z), Ey_theory, bounds_error=False, fill_value=0)
    interp_Ez_theory = RegularGridInterpolator((x, y, z), Ez_theory, bounds_error=False, fill_value=0)

    # Crear puntos a lo largo de la línea
    points = np.array([np.linspace(start[i], end[i], num_points) for i in range(3)]).T

    # Interpolar
    Ex_sim_line = interp_Ex_sim(points)
    Ey_sim_line = interp_Ey_sim(points)
    Ez_sim_line = interp_Ez_sim(points)

    Ex_theory_line = interp_Ex_theory(points)
    Ey_theory_line = interp_Ey_theory(points)
    Ez_theory_line = interp_Ez_theory(points)

    # Magnitudes
    E_sim_mag = np.sqrt(Ex_sim_line**2 + Ey_sim_line**2 + Ez_sim_line**2)
    E_theory_mag = np.sqrt(Ex_theory_line**2 + Ey_theory_line**2 + Ez_theory_line**2)

    # Distancia a lo largo de la línea
    distance = np.sqrt(np.sum((points - start)**2, axis=1))

    # Crear gráfico
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Componente X
    axes[0, 0].plot(distance, Ex_sim_line, 'b-', label='Simulado', linewidth=2)
    axes[0, 0].plot(distance, Ex_theory_line, 'r--', label='Teórico', linewidth=2)
    axes[0, 0].set_xlabel('Distancia', fontsize=12)
    axes[0, 0].set_ylabel('Ex', fontsize=12)
    axes[0, 0].set_title('Componente X del Campo Eléctrico', fontsize=14)
    axes[0, 0].legend(fontsize=10)
    axes[0, 0].grid(True, alpha=0.3)

    # Componente Y
    axes[0, 1].plot(distance, Ey_sim_line, 'b-', label='Simulado', linewidth=2)
    axes[0, 1].plot(distance, Ey_theory_line, 'r--', label='Teórico', linewidth=2)
    axes[0, 1].set_xlabel('Distancia', fontsize=12)
    axes[0, 1].set_ylabel('Ey', fontsize=12)
    axes[0, 1].set_title('Componente Y del Campo Eléctrico', fontsize=14)
    axes[0, 1].legend(fontsize=10)
    axes[0, 1].grid(True, alpha=0.3)

    # Componente Z
    axes[1, 0].plot(distance, Ez_sim_line, 'b-', label='Simulado', linewidth=2)
    axes[1, 0].plot(distance, Ez_theory_line, 'r--', label='Teórico', linewidth=2)
    axes[1, 0].set_xlabel('Distancia', fontsize=12)
    axes[1, 0].set_ylabel('Ez', fontsize=12)
    axes[1, 0].set_title('Componente Z del Campo Eléctrico', fontsize=14)
    axes[1, 0].legend(fontsize=10)
    axes[1, 0].grid(True, alpha=0.3)

    # Magnitud
    axes[1, 1].plot(distance, E_sim_mag, 'b-', label='Simulado', linewidth=2)
    axes[1, 1].plot(distance, E_theory_mag, 'r--', label='Teórico', linewidth=2)
    axes[1, 1].set_xlabel('Distancia', fontsize=12)
    axes[1, 1].set_ylabel('|E|', fontsize=12)
    axes[1, 1].set_title('Magnitud del Campo Eléctrico', fontsize=14)
    axes[1, 1].legend(fontsize=10)
    axes[1, 1].grid(True, alpha=0.3)

    plt.suptitle(f'Comparación Campo Simulado vs Teórico\n'
                 f'Línea desde ({start[0]:.1f}, {start[1]:.1f}, {start[2]:.1f}) '
                 f'hasta ({end[0]:.1f}, {end[1]:.1f}, {end[2]:.1f})',
                 fontsize=16, y=0.995)

    plt.tight_layout()

    if output:
        plt.savefig(output, dpi=300, bbox_inches='tight')
        print(f"Gráfico guardado en: {output}")
    else:
        plt.show()

    plt.close()

    # Calcular error relativo promedio
    error_x = np.mean(np.abs(Ex_sim_line - Ex_theory_line) / (np.abs(Ex_theory_line) + 1e-10))
    error_y = np.mean(np.abs(Ey_sim_line - Ey_theory_line) / (np.abs(Ey_theory_line) + 1e-10))
    error_z = np.mean(np.abs(Ez_sim_line - Ez_theory_line) / (np.abs(Ez_theory_line) + 1e-10))
    error_mag = np.mean(np.abs(E_sim_mag - E_theory_mag) / (E_theory_mag + 1e-10))

    print("\n" + "="*60)
    print("ANÁLISIS DE ERROR RELATIVO PROMEDIO")
    print("="*60)
    print(f"Error relativo Ex:  {error_x:.4f} ({error_x*100:.2f}%)")
    print(f"Error relativo Ey:  {error_y:.4f} ({error_y*100:.2f}%)")
    print(f"Error relativo Ez:  {error_z:.4f} ({error_z*100:.2f}%)")
    print(f"Error relativo |E|: {error_mag:.4f} ({error_mag*100:.2f}%)")
    print("="*60 + "\n")


def compare_fields_plane(Ex_sim, Ey_sim, Ez_sim, Ex_theory, Ey_theory, Ez_theory,
                         plane='xy', position=None, grid_size=(32, 32, 32), output=None):
    """
    Compara los campos simulado y teórico en un plano
    """
    nx, ny, nz = grid_size
    plane = plane.lower()

    if position is None:
        if plane == 'xy':
            position = nz // 2
        elif plane == 'xz':
            position = ny // 2
        elif plane == 'yz':
            position = nx // 2

    # Extraer slices
    if plane == 'xy':
        E_sim_mag = np.sqrt(Ex_sim[:, :, position]**2 +
                           Ey_sim[:, :, position]**2 +
                           Ez_sim[:, :, position]**2)
        E_theory_mag = np.sqrt(Ex_theory[:, :, position]**2 +
                              Ey_theory[:, :, position]**2 +
                              Ez_theory[:, :, position]**2)
        x = np.arange(0.5, nx)
        y = np.arange(0.5, ny)
        xlabel, ylabel = 'X', 'Y'
    elif plane == 'xz':
        E_sim_mag = np.sqrt(Ex_sim[:, position, :]**2 +
                           Ey_sim[:, position, :]**2 +
                           Ez_sim[:, position, :]**2)
        E_theory_mag = np.sqrt(Ex_theory[:, position, :]**2 +
                              Ey_theory[:, position, :]**2 +
                              Ez_theory[:, position, :]**2)
        x = np.arange(0.5, nx)
        y = np.arange(0.5, nz)
        xlabel, ylabel = 'X', 'Z'
    elif plane == 'yz':
        E_sim_mag = np.sqrt(Ex_sim[position, :, :]**2 +
                           Ey_sim[position, :, :]**2 +
                           Ez_sim[position, :, :]**2)
        E_theory_mag = np.sqrt(Ex_theory[position, :, :]**2 +
                              Ey_theory[position, :, :]**2 +
                              Ez_theory[position, :, :]**2)
        x = np.arange(0.5, ny)
        y = np.arange(0.5, nz)
        xlabel, ylabel = 'Y', 'Z'

    X, Y = np.meshgrid(x, y, indexing='ij')

    # Crear figura con 3 subplots
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # Campo simulado
    im1 = axes[0].contourf(X, Y, E_sim_mag, levels=20, cmap='viridis')
    axes[0].set_xlabel(xlabel, fontsize=12)
    axes[0].set_ylabel(ylabel, fontsize=12)
    axes[0].set_title('Campo Simulado |E|', fontsize=14)
    axes[0].set_aspect('equal')
    plt.colorbar(im1, ax=axes[0])

    # Campo teórico
    im2 = axes[1].contourf(X, Y, E_theory_mag, levels=20, cmap='viridis')
    axes[1].set_xlabel(xlabel, fontsize=12)
    axes[1].set_ylabel(ylabel, fontsize=12)
    axes[1].set_title('Campo Teórico |E|', fontsize=14)
    axes[1].set_aspect('equal')
    plt.colorbar(im2, ax=axes[1])

    # Diferencia
    diff = np.abs(E_sim_mag - E_theory_mag)
    im3 = axes[2].contourf(X, Y, diff, levels=20, cmap='RdYlBu_r')
    axes[2].set_xlabel(xlabel, fontsize=12)
    axes[2].set_ylabel(ylabel, fontsize=12)
    axes[2].set_title('Diferencia Absoluta |E_sim - E_theory|', fontsize=14)
    axes[2].set_aspect('equal')
    plt.colorbar(im3, ax=axes[2])

    plt.suptitle(f'Comparación en Plano {plane.upper()} (posición {position})', fontsize=16)
    plt.tight_layout()

    if output:
        plt.savefig(output, dpi=300, bbox_inches='tight')
        print(f"Gráfico guardado en: {output}")
    else:
        plt.show()

    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description='Comparar campo eléctrico simulado con campo teórico de carga puntual',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Ejemplos de uso:

1. Comparar a lo largo del eje X pasando por la carga:
   ./compare_field_theory.py -f psi-000030000.001-001 -s 32 32 32 \\
       --charge-pos 16 16 16 --charge 1.0 \\
       -m line --start 0.5 16 16 --end 31.5 16 16

2. Comparar a lo largo de una diagonal:
   ./compare_field_theory.py -f psi-000030000.001-001 -s 32 32 32 \\
       --charge-pos 16 16 16 --charge 1.0 \\
       -m line --start 0.5 0.5 0.5 --end 31.5 31.5 31.5

3. Comparar en un plano XY:
   ./compare_field_theory.py -f psi-000030000.001-001 -s 32 32 32 \\
       --charge-pos 16 16 16 --charge 1.0 \\
       -m plane -p xy --position 16

4. Especificar permitividad y guardar resultado:
   ./compare_field_theory.py -f psi-000030000.001-001 -s 32 32 32 \\
       --charge-pos 16 16 16 --charge 1.0 --epsilon 4.0 \\
       -m line --start 0.5 16 16 --end 31.5 16 16 -o comparison.png
        """
    )

    parser.add_argument('-f', '--file', required=True,
                       help='Archivo psi a leer')
    parser.add_argument('-s', '--size', nargs=3, type=int, required=True,
                       metavar=('NX', 'NY', 'NZ'),
                       help='Tamaño de la malla')
    parser.add_argument('--charge-pos', nargs=3, type=float, required=True,
                       metavar=('X', 'Y', 'Z'),
                       help='Posición de la carga puntual')
    parser.add_argument('--charge', type=float, default=1.0,
                       help='Magnitud de la carga (default: 1.0)')
    parser.add_argument('--epsilon', type=float, default=1.0,
                       help='Permitividad del medio (default: 1.0)')

    parser.add_argument('-m', '--mode', choices=['line', 'plane'], required=True,
                       help='Modo de comparación: line o plane')

    # Para modo línea
    parser.add_argument('--start', nargs=3, type=float, metavar=('X', 'Y', 'Z'),
                       help='Punto inicial de la línea (obligatorio para mode=line)')
    parser.add_argument('--end', nargs=3, type=float, metavar=('X', 'Y', 'Z'),
                       help='Punto final de la línea (obligatorio para mode=line)')
    parser.add_argument('--num-points', type=int, default=100,
                       help='Número de puntos a lo largo de la línea (default: 100)')

    # Para modo plano
    parser.add_argument('-p', '--plane', choices=['xy', 'xz', 'yz'],
                       help='Plano a comparar (obligatorio para mode=plane)')
    parser.add_argument('--position', type=int,
                       help='Posición del plano (índice). Default: centro')

    parser.add_argument('-o', '--output',
                       help='Archivo de salida (default: mostrar en pantalla)')

    args = parser.parse_args()

    # Validar argumentos
    if args.mode == 'line' and (args.start is None or args.end is None):
        parser.error("mode=line requiere --start y --end")
    if args.mode == 'plane' and args.plane is None:
        parser.error("mode=plane requiere --plane")

    if not os.path.exists(args.file):
        print(f"Error: No se encuentra el archivo {args.file}")
        sys.exit(1)

    # Leer archivo psi
    print(f"Leyendo archivo {args.file}...")
    reader = PsiReader(args.file, tuple(args.size))
    psi = reader.read()
    print(f"Archivo leído. Malla: {args.size[0]}x{args.size[1]}x{args.size[2]}")

    # Calcular campo simulado
    print("Calculando campo eléctrico simulado E = -∇ψ...")
    Ex_sim, Ey_sim, Ez_sim = reader.compute_electric_field()

    # Ajustar factor de escala automáticamente
    print("\nAjustando factor de escala entre simulación y teoría...")
    scale_factor = fit_scale_factor(Ex_sim, Ey_sim, Ez_sim, args.charge_pos, tuple(args.size), q=args.charge)
    print(f"Factor de escala encontrado: {scale_factor:.6e}")
    print(f"Factor teórico q/(4πε):      {args.charge / (4.0 * np.pi * args.epsilon):.6e}")
    print(f"Ratio simulado/teórico:      {scale_factor / (args.charge / (4.0 * np.pi * args.epsilon)):.4f}")

    # Calcular campo teórico con factor de escala ajustado
    print(f"\nCalculando campo teórico para carga q={args.charge} en "
          f"({args.charge_pos[0]}, {args.charge_pos[1]}, {args.charge_pos[2]}) "
          f"con factor de escala ajustado...")
    Ex_theory, Ey_theory, Ez_theory = theoretical_field(
        args.charge_pos, tuple(args.size), q=args.charge, epsilon=args.epsilon,
        scale_factor=scale_factor
    )

    # Comparar según el modo
    if args.mode == 'line':
        print(f"\nComparando a lo largo de línea desde {args.start} hasta {args.end}...")
        compare_fields_line(
            Ex_sim, Ey_sim, Ez_sim,
            Ex_theory, Ey_theory, Ez_theory,
            args.start, args.end, tuple(args.size),
            num_points=args.num_points,
            output=args.output
        )
    else:  # mode == 'plane'
        print(f"\nComparando en plano {args.plane.upper()}...")
        compare_fields_plane(
            Ex_sim, Ey_sim, Ez_sim,
            Ex_theory, Ey_theory, Ez_theory,
            plane=args.plane,
            position=args.position,
            grid_size=tuple(args.size),
            output=args.output
        )

    print("\n¡Comparación completada!")


if __name__ == '__main__':
    main()
