#!/usr/bin/env python3
"""
Script para graficar campo_electrico_medio.csv de múltiples simulaciones
variando el parámetro rhoel
"""

import os
import re
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def extract_rhoel(dirname):
    """Extrae el valor de rhoel del nombre del directorio"""
    match = re.search(r'rhoel_([\d.e+-]+)', dirname)
    if match:
        return float(match.group(1))
    return None

def plot_electric_field_rhoel(base_path, log_scale=True, symlog_threshold=1e-20):
    """
    Lee todos los archivos campo_electrico_medio.csv y los grafica
    según el parámetro rhoel

    Args:
        base_path: Ruta base donde buscar las subcarpetas
        log_scale: Si es True, usa escala logarítmica en el eje Y. Si es False, usa escala lineal
        symlog_threshold: Umbral para la escala symlog (usado en las componentes del campo)
    """
    # Buscar todos los archivos campo_electrico_medio.csv
    csv_files = []
    base_path = Path(base_path)

    for csv_file in base_path.rglob('campo_electrico_medio.csv'):
        parent_dir = csv_file.parent.name
        rhoel = extract_rhoel(parent_dir)
        if rhoel is not None:
            csv_files.append((rhoel, csv_file))

    if not csv_files:
        print(f"No se encontraron archivos campo_electrico_medio.csv en {base_path}")
        return

    # Ordenar por rhoel
    csv_files.sort(key=lambda x: x[0])

    print(f"Se encontraron {len(csv_files)} archivos:")
    for rhoel, filepath in csv_files:
        print(f"  rhoel = {rhoel}: {filepath.parent.name}")

    # Crear figura
    fig, ax = plt.subplots(figsize=(12, 8))

    # Colormap para diferentes valores de rhoel
    colors = plt.cm.viridis(np.linspace(0, 1, len(csv_files)))

    # Leer y graficar cada archivo
    for (rhoel, csv_file), color in zip(csv_files, colors):
        try:
            # Leer CSV
            data = np.loadtxt(csv_file, delimiter=',', skiprows=1)

            # Extraer posición x y campo eléctrico modular
            pos_x = data[:, 0]
            mean_Emod = data[:, 3]

            # Graficar
            ax.plot(pos_x, mean_Emod, 'o-', color=color,
                   label=f'rhoel = {rhoel}', linewidth=2, markersize=4, alpha=0.7)

        except Exception as e:
            print(f"Error al leer {csv_file}: {e}")

    # Configurar gráfico
    ax.set_xlabel('Posición X', fontsize=14)
    ax.set_ylabel('Campo Eléctrico Medio |E|', fontsize=14)
    title = 'Campo Eléctrico Medio vs Posición para diferentes valores de rhoel'
    if log_scale:
        title += ' (escala log)'
    ax.set_title(title, fontsize=16)
    ax.legend(loc='best', fontsize=10, ncol=2)
    ax.grid(True, alpha=0.3)

    if log_scale:
        ax.set_yscale('log')  # Escala logarítmica

    plt.tight_layout()

    # Guardar figura del módulo
    scale_suffix = '_log' if log_scale else '_linear'
    output_file = base_path / f'campo_electrico_rhoel_comparison{scale_suffix}.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nGráfico del módulo guardado en: {output_file}")
    plt.close()

    # Crear figura para las componentes del campo eléctrico
    fig, axes = plt.subplots(3, 1, figsize=(12, 14), sharex=True)

    # Leer y graficar cada archivo para las componentes
    for (rhoel, csv_file), color in zip(csv_files, colors):
        try:
            # Leer CSV
            data = np.loadtxt(csv_file, delimiter=',', skiprows=1)

            # Extraer posición x y componentes del campo eléctrico
            pos_x = data[:, 0]
            mean_Ex = data[:, 4]  # mean_Esub_X
            mean_Ey = data[:, 5]  # mean_Esub_Y
            mean_Ez = data[:, 6]  # mean_Esub_Z

            # Graficar cada componente
            axes[0].plot(pos_x, mean_Ex, 'o-', color=color,
                        label=f'rhoel = {rhoel}', linewidth=2, markersize=4, alpha=0.7)
            axes[1].plot(pos_x, mean_Ey, 'o-', color=color,
                        label=f'rhoel = {rhoel}', linewidth=2, markersize=4, alpha=0.7)
            axes[2].plot(pos_x, mean_Ez, 'o-', color=color,
                        label=f'rhoel = {rhoel}', linewidth=2, markersize=4, alpha=0.7)

        except Exception as e:
            print(f"Error al leer {csv_file} para componentes: {e}")

    # Configurar cada subplot
    comp_title = 'Componentes del Campo Eléctrico vs Posición para diferentes valores de rhoel'
    if log_scale:
        comp_title += ' (escala log)'

    axes[0].set_ylabel('Ex', fontsize=14)
    axes[0].set_title(comp_title, fontsize=16)
    axes[0].legend(loc='best', fontsize=10, ncol=2)
    axes[0].grid(True, alpha=0.3)
    axes[0].axhline(y=0, color='k', linestyle='--', alpha=0.3)
    if log_scale:
        axes[0].set_yscale('symlog', linthresh=symlog_threshold)

    axes[1].set_ylabel('Ey', fontsize=14)
    axes[1].legend(loc='best', fontsize=10, ncol=2)
    axes[1].grid(True, alpha=0.3)
    axes[1].axhline(y=0, color='k', linestyle='--', alpha=0.3)
    if log_scale:
        axes[1].set_yscale('symlog', linthresh=symlog_threshold)

    axes[2].set_ylabel('Ez', fontsize=14)
    axes[2].set_xlabel('Posición X', fontsize=14)
    axes[2].legend(loc='best', fontsize=10, ncol=2)
    axes[2].grid(True, alpha=0.3)
    axes[2].axhline(y=0, color='k', linestyle='--', alpha=0.3)
    if log_scale:
        axes[2].set_yscale('symlog', linthresh=symlog_threshold)

    plt.tight_layout()

    # Guardar figura de las componentes
    output_file_comp = base_path / f'campo_electrico_componentes_rhoel_comparison{scale_suffix}.png'
    plt.savefig(output_file_comp, dpi=300, bbox_inches='tight')
    print(f"Gráfico de componentes guardado en: {output_file_comp}")
    plt.close()

if __name__ == "__main__":
    import sys

    # Ruta base donde están las subcarpetas
    base_path = '../'

    # Determinar si usar escala logarítmica o lineal
    # Por defecto usa logarítmica, se puede pasar 'linear' como argumento
    use_log_scale = True
    symlog_thresh = 1e-20  # Valor por defecto

    if len(sys.argv) > 1:
        if sys.argv[1].lower() in ['linear', 'lin', 'l', 'false']:
            use_log_scale = False
            print("Usando escala lineal")
        else:
            print("Usando escala logarítmica")
    else:
        print("Usando escala logarítmica (por defecto)")
        print("Para usar escala lineal, ejecuta: python3 plot_campo_electrico_rhoel.py linear")

    # Parámetro opcional para el threshold del symlog
    if len(sys.argv) > 2:
        try:
            symlog_thresh = float(sys.argv[2])
            print(f"Usando symlog threshold: {symlog_thresh}")
        except ValueError:
            print(f"Advertencia: '{sys.argv[2]}' no es un número válido. Usando threshold por defecto: {symlog_thresh}")

    plot_electric_field_rhoel(base_path, log_scale=use_log_scale, symlog_threshold=symlog_thresh)
