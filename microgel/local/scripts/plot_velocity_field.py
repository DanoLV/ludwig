#!/usr/bin/env python3
"""
Script para visualizar campos de velocidad de archivos vel-* de Ludwig
Permite graficar en planos xy, xz, yz o cualquier plano definido por el usuario
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import argparse
import os
import glob

def read_velocity_file(filename, grid_size=None):
    """
    Lee un archivo de velocidades y devuelve un array 3D

    Args:
        filename: nombre del archivo vel-*
        grid_size: tupla (nx, ny, nz) opcional. Si None, asume malla cúbica

    Returns:
        velocity_field: array de forma (nx, ny, nz, 3) con componentes vx, vy, vz
        grid_size: tupla (nx, ny, nz)
    """
    # Leer datos
    data = np.loadtxt(filename)

    if grid_size is not None:
        # Usar dimensiones especificadas por el usuario
        nx, ny, nz = grid_size
        expected_points = nx * ny * nz
        actual_points = data.shape[0]

        if actual_points != expected_points:
            raise ValueError(
                f"El número de puntos en el archivo ({actual_points}) no coincide "
                f"con las dimensiones especificadas {nx}x{ny}x{nz} = {expected_points}"
            )

        print(f"Malla especificada: {nx}x{ny}x{nz} = {actual_points} puntos")
    else:
        # Calcular tamaño de la malla (asumiendo malla cúbica)
        total_points = data.shape[0]
        nx = int(round(total_points ** (1/3)))
        ny = nz = nx

        print(f"Malla detectada (cúbica): {nx}x{ny}x{nz} = {total_points} puntos")

    # Reorganizar datos en malla 3D
    # Asumiendo orden: z varía más rápido, luego y, luego x
    velocity_field = data.reshape(nx, ny, nz, 3)

    return velocity_field, (nx, ny, nz)

def plot_slice_xy(velocity_field, z_index, title="", save_path=None, skip=1, scale=None, colormap='viridis', magnitude_only=False, single_plot=False):
    """
    Grafica un plano XY en un índice z dado

    Args:
        velocity_field: array (nx, ny, nz, 3)
        z_index: índice z del plano a graficar
        title: título del gráfico
        save_path: ruta para guardar la figura (opcional)
        skip: factor de submuestreo para las flechas (default: 1, graficar todas)
        scale: escala para las flechas (opcional, auto si None)
        colormap: mapa de colores para la magnitud
        magnitude_only: si True, grafica solo la magnitud en plano sin vectores (default: False)
        single_plot: si True, grafica solo magnitud+vectores en un plot (default: False, 2 plots)
    """
    nx, ny, nz, _ = velocity_field.shape

    # Extraer plano
    vx = velocity_field[:, :, z_index, 0]
    vy = velocity_field[:, :, z_index, 1]
    vz = velocity_field[:, :, z_index, 2]

    # Crear malla
    x = np.arange(nx)
    y = np.arange(ny)
    X, Y = np.meshgrid(x, y, indexing='ij')

    # Calcular magnitud del campo de velocidad en el plano
    v_magnitude = np.sqrt(vx**2 + vy**2)
    v_total_magnitude = np.sqrt(vx**2 + vy**2 + vz**2)

    if magnitude_only:
        # Graficar solo la magnitud en plano (sin vectores)
        fig, ax = plt.subplots(1, 1, figsize=(10, 8))
        c = ax.contourf(X, Y, v_magnitude, levels=20, cmap=colormap)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_title(f'{title}\nPlano XY (z={z_index}) - Magnitud en plano')
        ax.set_aspect('equal')
        cbar = plt.colorbar(c, ax=ax, label='|v| en plano XY')
        plt.tight_layout()
    elif single_plot:
        # Graficar solo magnitud + vectores (panel izquierdo del caso default)
        fig, ax = plt.subplots(1, 1, figsize=(10, 8))
        c = ax.contourf(X, Y, v_magnitude, levels=20, cmap=colormap, alpha=0.6)
        q = ax.quiver(X[::skip, ::skip], Y[::skip, ::skip],
                      vx[::skip, ::skip], vy[::skip, ::skip],
                      scale=scale, scale_units='xy', alpha=0.8)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_title(f'{title}\nPlano XY (z={z_index}) - Magnitud en plano')
        ax.set_aspect('equal')
        plt.colorbar(c, ax=ax, label='|v| en plano XY')
        ax.quiverkey(q, 0.9, 1.05, np.max(v_magnitude), f'{np.max(v_magnitude):.2e}', labelpos='E')
        plt.tight_layout()
    else:
        # Crear figura con dos subplots
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))

        # Subplot 1: Vectores en el plano XY con magnitud en plano
        c1 = ax1.contourf(X, Y, v_magnitude, levels=20, cmap=colormap, alpha=0.6)
        q1 = ax1.quiver(X[::skip, ::skip], Y[::skip, ::skip],
                        vx[::skip, ::skip], vy[::skip, ::skip],
                        scale=scale, scale_units='xy', alpha=0.8)
        ax1.set_xlabel('X')
        ax1.set_ylabel('Y')
        ax1.set_title(f'{title}\nPlano XY (z={z_index}) - Magnitud en plano')
        ax1.set_aspect('equal')
        plt.colorbar(c1, ax=ax1, label='|v| en plano XY')
        ax1.quiverkey(q1, 0.9, 1.05, np.max(v_magnitude), f'{np.max(v_magnitude):.2e}', labelpos='E')

        # Subplot 2: Magnitud total con componente z
        c2 = ax2.contourf(X, Y, v_total_magnitude, levels=20, cmap=colormap, alpha=0.6)
        # Colores basados en vz
        colors_vz = vz[::skip, ::skip]
        q2 = ax2.quiver(X[::skip, ::skip], Y[::skip, ::skip],
                        vx[::skip, ::skip], vy[::skip, ::skip],
                        colors_vz, cmap='RdBu_r', scale=scale, scale_units='xy', alpha=0.8)
        ax2.set_xlabel('X')
        ax2.set_ylabel('Y')
        ax2.set_title(f'{title}\nPlano XY (z={z_index}) - Color: Vz')
        ax2.set_aspect('equal')
        plt.colorbar(c2, ax=ax2, label='|v| total')
        plt.colorbar(q2, ax=ax2, label='Vz')

        plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Figura guardada en: {save_path}")

    return fig

def plot_slice_xz(velocity_field, y_index, title="", save_path=None, skip=1, scale=None, colormap='viridis', magnitude_only=False, single_plot=False):
    """
    Grafica un plano XZ en un índice y dado
    """
    nx, ny, nz, _ = velocity_field.shape

    # Extraer plano
    vx = velocity_field[:, y_index, :, 0]
    vy = velocity_field[:, y_index, :, 1]
    vz = velocity_field[:, y_index, :, 2]

    # Crear malla
    x = np.arange(nx)
    z = np.arange(nz)
    X, Z = np.meshgrid(x, z, indexing='ij')

    # Calcular magnitud
    v_magnitude = np.sqrt(vx**2 + vz**2)
    v_total_magnitude = np.sqrt(vx**2 + vy**2 + vz**2)

    if magnitude_only:
        # Graficar solo la magnitud en plano (sin vectores)
        fig, ax = plt.subplots(1, 1, figsize=(10, 8))
        c = ax.contourf(X, Z, v_magnitude, levels=20, cmap=colormap)
        ax.set_xlabel('X')
        ax.set_ylabel('Z')
        ax.set_title(f'{title}\nPlano XZ (y={y_index}) - Magnitud en plano')
        ax.set_aspect('equal')
        cbar = plt.colorbar(c, ax=ax, label='|v| en plano XZ')
        plt.tight_layout()
    elif single_plot:
        # Graficar solo magnitud + vectores
        fig, ax = plt.subplots(1, 1, figsize=(10, 8))
        c = ax.contourf(X, Z, v_magnitude, levels=20, cmap=colormap, alpha=0.6)
        q = ax.quiver(X[::skip, ::skip], Z[::skip, ::skip],
                      vx[::skip, ::skip], vz[::skip, ::skip],
                      scale=scale, scale_units='xy', alpha=0.8)
        ax.set_xlabel('X')
        ax.set_ylabel('Z')
        ax.set_title(f'{title}\nPlano XZ (y={y_index}) - Magnitud en plano')
        ax.set_aspect('equal')
        plt.colorbar(c, ax=ax, label='|v| en plano XZ')
        ax.quiverkey(q, 0.9, 1.05, np.max(v_magnitude), f'{np.max(v_magnitude):.2e}', labelpos='E')
        plt.tight_layout()
    else:
        # Crear figura
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))

        # Subplot 1: Vectores en el plano XZ
        c1 = ax1.contourf(X, Z, v_magnitude, levels=20, cmap=colormap, alpha=0.6)
        q1 = ax1.quiver(X[::skip, ::skip], Z[::skip, ::skip],
                        vx[::skip, ::skip], vz[::skip, ::skip],
                        scale=scale, scale_units='xy', alpha=0.8)
        ax1.set_xlabel('X')
        ax1.set_ylabel('Z')
        ax1.set_title(f'{title}\nPlano XZ (y={y_index}) - Magnitud en plano')
        ax1.set_aspect('equal')
        plt.colorbar(c1, ax=ax1, label='|v| en plano XZ')
        ax1.quiverkey(q1, 0.9, 1.05, np.max(v_magnitude), f'{np.max(v_magnitude):.2e}', labelpos='E')

        # Subplot 2: Magnitud total con componente y
        c2 = ax2.contourf(X, Z, v_total_magnitude, levels=20, cmap=colormap, alpha=0.6)
        colors_vy = vy[::skip, ::skip]
        q2 = ax2.quiver(X[::skip, ::skip], Z[::skip, ::skip],
                        vx[::skip, ::skip], vz[::skip, ::skip],
                        colors_vy, cmap='RdBu_r', scale=scale, scale_units='xy', alpha=0.8)
        ax2.set_xlabel('X')
        ax2.set_ylabel('Z')
        ax2.set_title(f'{title}\nPlano XZ (y={y_index}) - Color: Vy')
        ax2.set_aspect('equal')
        plt.colorbar(c2, ax=ax2, label='|v| total')
        plt.colorbar(q2, ax=ax2, label='Vy')

        plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Figura guardada en: {save_path}")

    return fig

def plot_slice_yz(velocity_field, x_index, title="", save_path=None, skip=1, scale=None, colormap='viridis', magnitude_only=False, single_plot=False):
    """
    Grafica un plano YZ en un índice x dado
    """
    nx, ny, nz, _ = velocity_field.shape

    # Extraer plano
    vx = velocity_field[x_index, :, :, 0]
    vy = velocity_field[x_index, :, :, 1]
    vz = velocity_field[x_index, :, :, 2]

    # Crear malla
    y = np.arange(ny)
    z = np.arange(nz)
    Y, Z = np.meshgrid(y, z, indexing='ij')

    # Calcular magnitud
    v_magnitude = np.sqrt(vy**2 + vz**2)
    v_total_magnitude = np.sqrt(vx**2 + vy**2 + vz**2)

    if magnitude_only:
        # Graficar solo la magnitud en plano (sin vectores)
        fig, ax = plt.subplots(1, 1, figsize=(10, 8))
        c = ax.contourf(Y, Z, v_magnitude, levels=20, cmap=colormap)
        ax.set_xlabel('Y')
        ax.set_ylabel('Z')
        ax.set_title(f'{title}\nPlano YZ (x={x_index}) - Magnitud en plano')
        ax.set_aspect('equal')
        cbar = plt.colorbar(c, ax=ax, label='|v| en plano YZ')
        plt.tight_layout()
    elif single_plot:
        # Graficar solo magnitud + vectores
        fig, ax = plt.subplots(1, 1, figsize=(10, 8))
        c = ax.contourf(Y, Z, v_magnitude, levels=20, cmap=colormap, alpha=0.6)
        q = ax.quiver(Y[::skip, ::skip], Z[::skip, ::skip],
                      vy[::skip, ::skip], vz[::skip, ::skip],
                      scale=scale, scale_units='xy', alpha=0.8)
        ax.set_xlabel('Y')
        ax.set_ylabel('Z')
        ax.set_title(f'{title}\nPlano YZ (x={x_index}) - Magnitud en plano')
        ax.set_aspect('equal')
        plt.colorbar(c, ax=ax, label='|v| en plano YZ')
        ax.quiverkey(q, 0.9, 1.05, np.max(v_magnitude), f'{np.max(v_magnitude):.2e}', labelpos='E')
        plt.tight_layout()
    else:
        # Crear figura
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))

        # Subplot 1: Vectores en el plano YZ
        c1 = ax1.contourf(Y, Z, v_magnitude, levels=20, cmap=colormap, alpha=0.6)
        q1 = ax1.quiver(Y[::skip, ::skip], Z[::skip, ::skip],
                        vy[::skip, ::skip], vz[::skip, ::skip],
                        scale=scale, scale_units='xy', alpha=0.8)
        ax1.set_xlabel('Y')
        ax1.set_ylabel('Z')
        ax1.set_title(f'{title}\nPlano YZ (x={x_index}) - Magnitud en plano')
        ax1.set_aspect('equal')
        plt.colorbar(c1, ax=ax1, label='|v| en plano YZ')
        ax1.quiverkey(q1, 0.9, 1.05, np.max(v_magnitude), f'{np.max(v_magnitude):.2e}', labelpos='E')

        # Subplot 2: Magnitud total con componente x
        c2 = ax2.contourf(Y, Z, v_total_magnitude, levels=20, cmap=colormap, alpha=0.6)
        colors_vx = vx[::skip, ::skip]
        q2 = ax2.quiver(Y[::skip, ::skip], Z[::skip, ::skip],
                        vy[::skip, ::skip], vz[::skip, ::skip],
                        colors_vx, cmap='RdBu_r', scale=scale, scale_units='xy', alpha=0.8)
        ax2.set_xlabel('Y')
        ax2.set_ylabel('Z')
        ax2.set_title(f'{title}\nPlano YZ (x={x_index}) - Color: Vx')
        ax2.set_aspect('equal')
        plt.colorbar(c2, ax=ax2, label='|v| total')
        plt.colorbar(q2, ax=ax2, label='Vx')

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Figura guardada en: {save_path}")

    return fig

def plot_all_planes(velocity_field, grid_size, output_dir=".", prefix="vel", skip=1, scale=None, magnitude_only=False, single_plot=False):
    """
    Crea gráficos para los tres planos principales en el centro de la malla
    """
    nx, ny, nz = grid_size

    # Plano XY en el centro
    z_mid = nz // 2
    fig1 = plot_slice_xy(velocity_field, z_mid,
                         title=f"Campo de velocidad",
                         save_path=f"{output_dir}/{prefix}_xy_z{z_mid}.png",
                         skip=skip, scale=scale, magnitude_only=magnitude_only, single_plot=single_plot)

    # Plano XZ en el centro
    y_mid = ny // 2
    fig2 = plot_slice_xz(velocity_field, y_mid,
                         title=f"Campo de velocidad",
                         save_path=f"{output_dir}/{prefix}_xz_y{y_mid}.png",
                         skip=skip, scale=scale, magnitude_only=magnitude_only, single_plot=single_plot)

    # Plano YZ en el centro
    x_mid = nx // 2
    fig3 = plot_slice_yz(velocity_field, x_mid,
                         title=f"Campo de velocidad",
                         save_path=f"{output_dir}/{prefix}_yz_x{x_mid}.png",
                         skip=skip, scale=scale, magnitude_only=magnitude_only, single_plot=single_plot)

    return fig1, fig2, fig3

def find_velocity_file(directory, timestep):
    """
    Busca archivos de velocidad para un timestep dado

    Args:
        directory: directorio donde buscar
        timestep: número de paso temporal

    Returns:
        lista de archivos que coinciden
    """
    pattern = f"{directory}/vel-{timestep:09d}.*"
    files = glob.glob(pattern)
    return files

def main():
    parser = argparse.ArgumentParser(
        description='Graficar campos de velocidad de archivos vel-* de Ludwig',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Ejemplos de uso:
  # Graficar planos principales para el paso 50 (malla cúbica auto-detectada)
  %(prog)s -d ./directorio -t 50

  # Graficar con dimensiones específicas (por ejemplo, malla 2D: 64x64x4)
  %(prog)s -f vel-000000050.001-001 -s 64 64 4

  # Graficar solo plano XY en z=2 con dimensiones específicas
  %(prog)s -f vel-000000050.001-001 -s 64 64 4 --plane xy --index 2

  # Graficar con submuestreo de flechas (cada 2 puntos)
  %(prog)s -f vel-000000050.001-001 -s 64 64 4 --skip 2

  # Graficar plano XZ en y=32 con escala personalizada
  %(prog)s -f vel-000000050.001-001 -s 64 64 4 --plane xz --index 32 --scale 100

  # Graficar solo la magnitud en plano (sin vectores)
  %(prog)s -f vel-000000050.001-001 -s 64 64 4 --plane xy --magnitude-only

  # Graficar solo un panel con magnitud y vectores
  %(prog)s -f vel-000000050.001-001 -s 64 64 4 --plane xy --single-plot --skip 2
        """
    )

    # Argumentos principales
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-f', '--file', type=str,
                       help='Archivo de velocidad específico')
    group.add_argument('-d', '--directory', type=str,
                       help='Directorio con archivos de velocidad')

    parser.add_argument('-t', '--timestep', type=int,
                        help='Número de paso temporal (requerido si se usa -d)')

    # Dimensiones del grid
    parser.add_argument('-s', '--size', nargs=3, type=int,
                        metavar=('NX', 'NY', 'NZ'),
                        help='Tamaño de la malla (nx ny nz). Si no se especifica, asume malla cúbica')

    # Opciones de plano
    parser.add_argument('--plane', type=str, choices=['xy', 'xz', 'yz', 'all'],
                        default='all',
                        help='Plano a graficar (default: all)')
    parser.add_argument('--index', type=int,
                        help='Índice del plano (default: centro de la malla)')

    # Opciones de visualización
    parser.add_argument('--skip', type=int, default=1,
                        help='Factor de submuestreo para flechas (default: 1)')
    parser.add_argument('--scale', type=float,
                        help='Escala para las flechas (default: auto)')
    parser.add_argument('--colormap', type=str, default='viridis',
                        help='Mapa de colores (default: viridis)')
    parser.add_argument('--magnitude-only', action='store_true',
                        help='Graficar solo la magnitud en plano sin vectores')
    parser.add_argument('--single-plot', action='store_true',
                        help='Graficar solo un panel con magnitud y vectores (en lugar de 2 paneles)')

    # Opciones de salida
    parser.add_argument('-o', '--output', type=str, default='.',
                        help='Directorio de salida (default: .)')
    parser.add_argument('--prefix', type=str, default='vel',
                        help='Prefijo para archivos de salida (default: vel)')
    parser.add_argument('--show', action='store_true',
                        help='Mostrar gráficos interactivamente')

    args = parser.parse_args()

    # Validar argumentos
    if args.directory and not args.timestep:
        parser.error("-d/--directory requiere -t/--timestep")

    # Determinar archivo a leer
    if args.file:
        vel_file = args.file
    else:
        files = find_velocity_file(args.directory, args.timestep)
        if not files:
            print(f"Error: No se encontraron archivos para el timestep {args.timestep}")
            return 1
        if len(files) > 1:
            print(f"Se encontraron múltiples archivos: {files}")
            print(f"Usando: {files[0]}")
        vel_file = files[0]

    # Verificar que el archivo existe
    if not os.path.exists(vel_file):
        print(f"Error: El archivo {vel_file} no existe")
        return 1

    print(f"Leyendo archivo: {vel_file}")

    # Leer datos
    grid_size_input = tuple(args.size) if args.size else None
    velocity_field, grid_size = read_velocity_file(vel_file, grid_size=grid_size_input)
    nx, ny, nz = grid_size

    # Crear directorio de salida si no existe
    os.makedirs(args.output, exist_ok=True)

    # Determinar índices
    if args.plane == 'all':
        # Graficar todos los planos en el centro
        print("\nGraficando todos los planos...")
        plot_all_planes(velocity_field, grid_size,
                       output_dir=args.output, prefix=args.prefix,
                       skip=args.skip, scale=args.scale,
                       magnitude_only=args.magnitude_only,
                       single_plot=args.single_plot)
    else:
        # Graficar plano específico
        if args.index is None:
            # Usar centro si no se especifica
            if args.plane == 'xy':
                idx = nz // 2
            elif args.plane == 'xz':
                idx = ny // 2
            else:  # yz
                idx = nx // 2
        else:
            idx = args.index

        save_path = f"{args.output}/{args.prefix}_{args.plane}_{idx}.png"

        print(f"\nGraficando plano {args.plane.upper()} en índice {idx}...")

        if args.plane == 'xy':
            plot_slice_xy(velocity_field, idx,
                         title=f"Campo de velocidad",
                         save_path=save_path,
                         skip=args.skip, scale=args.scale,
                         colormap=args.colormap,
                         magnitude_only=args.magnitude_only,
                         single_plot=args.single_plot)
        elif args.plane == 'xz':
            plot_slice_xz(velocity_field, idx,
                         title=f"Campo de velocidad",
                         save_path=save_path,
                         skip=args.skip, scale=args.scale,
                         colormap=args.colormap,
                         magnitude_only=args.magnitude_only,
                         single_plot=args.single_plot)
        else:  # yz
            plot_slice_yz(velocity_field, idx,
                         title=f"Campo de velocidad",
                         save_path=save_path,
                         skip=args.skip, scale=args.scale,
                         colormap=args.colormap,
                         magnitude_only=args.magnitude_only,
                         single_plot=args.single_plot)

    print("\n¡Listo!")

    if args.show:
        plt.show()

    return 0

if __name__ == "__main__":
    exit(main())
