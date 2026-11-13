#!/usr/bin/env python3
"""
Script para comparar el potencial eléctrico simulado con el potencial teórico
de una carga puntual.

Para una carga puntual q, el potencial eléctrico es:
ψ = q / (4πε|r|)

donde |r| es la distancia desde la carga al punto de observación.
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


def theoretical_potential(charge_pos, grid_size, q=1.0, epsilon=1.0, scale_factor=None):
    """
    Calcula el potencial eléctrico teórico de una carga puntual

    Args:
        charge_pos: Posición de la carga (x, y, z)
        grid_size: Tupla (nx, ny, nz)
        q: Carga (default: 1.0)
        epsilon: Permitividad (default: 1.0)
        scale_factor: Factor de escala manual. Si es None, usa q/(4πε)

    Returns:
        psi_theory: Potencial teórico en la malla
    """
    nx, ny, nz = grid_size
    xc, yc, zc = charge_pos

    # Crear malla de coordenadas centradas en nodos (0.5, 1.5, 2.5, ...)
    x = np.arange(0.5, nx)
    y = np.arange(0.5, ny)
    z = np.arange(0.5, nz)

    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

    # Distancia desde la carga
    r = np.sqrt((X - xc)**2 + (Y - yc)**2 + (Z - zc)**2)

    # Evitar división por cero en la posición de la carga
    r = np.where(r < 0.1, 0.1, r)

    # Potencial teórico: ψ = (q / 4πε) / r
    if scale_factor is None:
        # Usar factor teórico estándar
        prefactor = q / (4.0 * np.pi * epsilon)
    else:
        # Usar factor de escala proporcionado
        prefactor = scale_factor

    psi_theory = prefactor / r

    return psi_theory


def fit_scale_factor(psi_sim, charge_pos, grid_size, q=1.0):
    """
    Encuentra el factor de escala óptimo comparando el potencial simulado con el patrón teórico

    Args:
        psi_sim: Potencial simulado
        charge_pos: Posición de la carga
        grid_size: Tupla (nx, ny, nz)
        q: Carga

    Returns:
        scale_factor: Factor de escala óptimo
    """
    # Calcular potencial teórico con factor unitario
    psi_unit = theoretical_potential(charge_pos, grid_size, q=1.0, epsilon=1.0, scale_factor=1.0)

    # Crear máscara para excluir región cercana a la carga
    nx, ny, nz = grid_size
    xc, yc, zc = charge_pos
    x = np.arange(0.5, nx)
    y = np.arange(0.5, ny)
    z = np.arange(0.5, nz)
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

    r = np.sqrt((X - xc)**2 + (Y - yc)**2 + (Z - zc)**2)
    mask = r > 2.0  # Excluir región dentro de 2 unidades de lattice

    # Calcular factor de escala mediante ajuste de mínimos cuadrados
    # scale_factor = sum(psi_sim * psi_unit) / sum(psi_unit^2)
    numerator = np.sum(psi_sim[mask] * psi_unit[mask])
    denominator = np.sum(psi_unit[mask]**2)

    scale_factor = numerator / denominator

    return scale_factor


def compare_potential_line(psi_sim, psi_theory, start, end, grid_size,
                           num_points=100, charge_pos=None, output=None):
    """
    Compara los potenciales simulado y teórico a lo largo de una línea
    """
    from scipy.interpolate import RegularGridInterpolator

    nx, ny, nz = grid_size

    # Coordenadas de la malla
    x = np.arange(0.5, nx)
    y = np.arange(0.5, ny)
    z = np.arange(0.5, nz)

    # Crear interpoladores
    interp_psi_sim = RegularGridInterpolator((x, y, z), psi_sim, bounds_error=False, fill_value=0)
    interp_psi_theory = RegularGridInterpolator((x, y, z), psi_theory, bounds_error=False, fill_value=0)

    # Crear puntos a lo largo de la línea
    points = np.array([np.linspace(start[i], end[i], num_points) for i in range(3)]).T

    # Interpolar
    psi_sim_line = interp_psi_sim(points)
    psi_theory_line = interp_psi_theory(points)

    # Distancia a lo largo de la línea
    distance = np.sqrt(np.sum((points - start)**2, axis=1))

    # Si se proporciona posición de carga, calcular distancia desde la carga
    if charge_pos is not None:
        r_from_charge = np.sqrt(np.sum((points - charge_pos)**2, axis=1))
    else:
        r_from_charge = distance

    # Crear gráfico
    fig, axes = plt.subplots(2, 1, figsize=(12, 10))

    # Potencial vs distancia
    axes[0].plot(distance, psi_sim_line, 'b-', label='Simulado', linewidth=2)
    axes[0].plot(distance, psi_theory_line, 'r--', label='Teórico', linewidth=2)
    axes[0].set_xlabel('Distancia a lo largo de la línea', fontsize=12)
    axes[0].set_ylabel('Potencial ψ', fontsize=12)
    axes[0].set_title('Potencial Eléctrico', fontsize=14)
    axes[0].legend(fontsize=11)
    axes[0].grid(True, alpha=0.3)

    # Diferencia absoluta
    diff = np.abs(psi_sim_line - psi_theory_line)
    axes[1].plot(distance, diff, 'g-', linewidth=2)
    axes[1].set_xlabel('Distancia a lo largo de la línea', fontsize=12)
    axes[1].set_ylabel('|ψ_sim - ψ_theory|', fontsize=12)
    axes[1].set_title('Diferencia Absoluta', fontsize=14)
    axes[1].grid(True, alpha=0.3)

    plt.suptitle(f'Comparación Potencial Simulado vs Teórico\n'
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

    # Calcular error relativo promedio (excluyendo región muy cercana a la carga)
    mask = r_from_charge > 2.0
    if np.any(mask):
        error_rel = np.mean(np.abs(psi_sim_line[mask] - psi_theory_line[mask]) /
                           (np.abs(psi_theory_line[mask]) + 1e-10))

        print("\n" + "="*60)
        print("ANÁLISIS DE ERROR RELATIVO PROMEDIO")
        print("="*60)
        print(f"Error relativo ψ (r > 2.0): {error_rel:.4f} ({error_rel*100:.2f}%)")
        print("="*60 + "\n")


def compare_potential_plane(psi_sim, psi_theory, plane='xy', position=None,
                            grid_size=(32, 32, 32), output=None):
    """
    Compara los potenciales simulado y teórico en un plano
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
        psi_sim_slice = psi_sim[:, :, position]
        psi_theory_slice = psi_theory[:, :, position]
        x = np.arange(0.5, nx)
        y = np.arange(0.5, ny)
        xlabel, ylabel = 'X', 'Y'
    elif plane == 'xz':
        psi_sim_slice = psi_sim[:, position, :]
        psi_theory_slice = psi_theory[:, position, :]
        x = np.arange(0.5, nx)
        y = np.arange(0.5, nz)
        xlabel, ylabel = 'X', 'Z'
    elif plane == 'yz':
        psi_sim_slice = psi_sim[position, :, :]
        psi_theory_slice = psi_theory[position, :, :]
        x = np.arange(0.5, ny)
        y = np.arange(0.5, nz)
        xlabel, ylabel = 'Y', 'Z'

    X, Y = np.meshgrid(x, y, indexing='ij')

    # Crear figura con 3 subplots
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # Potencial simulado
    im1 = axes[0].contourf(X, Y, psi_sim_slice, levels=20, cmap='viridis')
    axes[0].set_xlabel(xlabel, fontsize=12)
    axes[0].set_ylabel(ylabel, fontsize=12)
    axes[0].set_title('Potencial Simulado ψ', fontsize=14)
    axes[0].set_aspect('equal')
    plt.colorbar(im1, ax=axes[0])

    # Potencial teórico
    im2 = axes[1].contourf(X, Y, psi_theory_slice, levels=20, cmap='viridis')
    axes[1].set_xlabel(xlabel, fontsize=12)
    axes[1].set_ylabel(ylabel, fontsize=12)
    axes[1].set_title('Potencial Teórico ψ', fontsize=14)
    axes[1].set_aspect('equal')
    plt.colorbar(im2, ax=axes[1])

    # Diferencia relativa
    diff_rel = np.abs((psi_sim_slice - psi_theory_slice) / (psi_theory_slice + 1e-10))
    im3 = axes[2].contourf(X, Y, diff_rel, levels=20, cmap='RdYlBu_r')
    axes[2].set_xlabel(xlabel, fontsize=12)
    axes[2].set_ylabel(ylabel, fontsize=12)
    axes[2].set_title('Error Relativo |Δψ/ψ_theory|', fontsize=14)
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
        description='Comparar potencial eléctrico simulado con potencial teórico de carga puntual',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Ejemplos de uso:

1. Comparar a lo largo del eje X pasando por la carga:
   ./compare_potential_theory.py -f psi-000030000.001-001 -s 32 32 32 \\
       --charge-pos 16 16 16 --charge 1.0 --epsilon 1.0e4 \\
       -m line --start 0.5 16 16 --end 31.5 16 16

2. Comparar a lo largo de una diagonal:
   ./compare_potential_theory.py -f psi-000030000.001-001 -s 32 32 32 \\
       --charge-pos 16 16 16 --charge 1.0 --epsilon 1.0e4 \\
       -m line --start 0.5 0.5 0.5 --end 31.5 31.5 31.5

3. Comparar en un plano XY:
   ./compare_potential_theory.py -f psi-000030000.001-001 -s 32 32 32 \\
       --charge-pos 16 16 16 --charge 1.0 --epsilon 1.0e4 \\
       -m plane -p xy --position 16

4. Guardar resultado:
   ./compare_potential_theory.py -f psi-000030000.001-001 -s 32 32 32 \\
       --charge-pos 16 16 16 --charge 1.0 --epsilon 1.0e4 \\
       -m line --start 0.5 16 16 --end 31.5 16 16 -o comparison_psi.png
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
    psi_sim = reader.read()
    print(f"Archivo leído. Malla: {args.size[0]}x{args.size[1]}x{args.size[2]}")

    # Ajustar factor de escala automáticamente
    print("\nAjustando factor de escala entre simulación y teoría...")
    scale_factor = fit_scale_factor(psi_sim, args.charge_pos, tuple(args.size), q=args.charge)
    print(f"Factor de escala encontrado: {scale_factor:.6e}")
    print(f"Factor teórico q/(4πε):      {args.charge / (4.0 * np.pi * args.epsilon):.6e}")
    print(f"Ratio simulado/teórico:      {scale_factor / (args.charge / (4.0 * np.pi * args.epsilon)):.4f}")

    # Calcular potencial teórico con factor de escala ajustado
    print(f"\nCalculando potencial teórico para carga q={args.charge} en "
          f"({args.charge_pos[0]}, {args.charge_pos[1]}, {args.charge_pos[2]}) "
          f"con factor de escala ajustado...")
    psi_theory = theoretical_potential(
        args.charge_pos, tuple(args.size), q=args.charge, epsilon=args.epsilon,
        scale_factor=scale_factor
    )

    # Comparar según el modo
    if args.mode == 'line':
        print(f"\nComparando a lo largo de línea desde {args.start} hasta {args.end}...")
        compare_potential_line(
            psi_sim, psi_theory,
            args.start, args.end, tuple(args.size),
            num_points=args.num_points,
            charge_pos=args.charge_pos,
            output=args.output
        )
    else:  # mode == 'plane'
        print(f"\nComparando en plano {args.plane.upper()}...")
        compare_potential_plane(
            psi_sim, psi_theory,
            plane=args.plane,
            position=args.position,
            grid_size=tuple(args.size),
            output=args.output
        )

    print("\n¡Comparación completada!")


if __name__ == '__main__':
    main()
