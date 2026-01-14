#!/usr/bin/env python3
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
import sys
from pathlib import Path

paso=1000

# Parsear argumentos
show_plot = True  # Por defecto mostrar el gráfico
base_dir = None

for arg in sys.argv[1:]:
    if arg == '--no-show':
        show_plot = False
    else:
        base_dir = arg

# Directorio base
if base_dir is None:
    # Si no se pasa directorio, usar el comportamiento original
    os.chdir("..")
    base_dir = os.getcwd()

# Lista para almacenar los resultados
results = []

# Recorrer todas las subcarpetas
subdirs = sorted([d for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d))])

print(f"Procesando {len(subdirs)} carpetas...")

for subdir in subdirs:
    csv_file = os.path.join(base_dir, subdir, 'proceced_data', 'particle_Esub.csv')
    if not os.path.isfile(csv_file):
        continue
    
    # Verificar que el archivo existe
    if not os.path.exists(csv_file):
        print(f"Archivo no encontrado: {csv_file}")
        continue

    # Extraer posiciones x, y, z del nombre de la carpeta
    # Formato: single-r_0.8-v_0.02-kT_0.0005-ah_1.605000000000000e-02-lb_fluctuation_0-fe_electro-q_1.0-e_0.0_0.0_0.0-p200-N_50k_16.00_16.50_16.50
    match = re.search(r'pos_([\d\.]+)_([\d\.]+)_([\d\.]+)', subdir)
    if match:
        pos_x = float(match.group(1))
        pos_y = float(match.group(2))
        pos_z = float(match.group(3))
    else:
        print(f"No se pudo extraer posición de: {subdir}")
        continue

    print(f"Procesando {subdir} (x={pos_x}, y={pos_y}, z={pos_z})...")

    # Leer el archivo CSV (usando coma como separador decimal y punto y coma como separador de columna)
    try:
        # Leer archivo, saltando la primera línea (comentario con #)
        df = pd.read_csv(csv_file, sep=';', decimal=',', skiprows=1, header=None,
                         names=['Step', 'Index', 'Emod', 'Esub_X', 'Esub_Y', 'Esub_Z', 'EmodPB', 'EPB_X', 'EPB_Y', 'EPB_Z'] )

        # Limpiar espacios en blanco de todas las columnas
        for col in df.columns:
            if df[col].dtype == 'object':
                df[col] = df[col].str.strip()

        # Filtrar datos desde el paso 1000
        df_filtered = df[df['Step'] >= paso]

        if len(df_filtered) == 0:
            print(f"  No hay datos después del paso {paso} en {subdir}")
            continue

        # Calcular medias (excluyendo Step e Index)
        mean_Emod = df_filtered['Emod'].mean()
        mean_Esub_X = df_filtered['Esub_X'].mean()
        mean_Esub_Y = df_filtered['Esub_Y'].mean()
        mean_Esub_Z = df_filtered['Esub_Z'].mean()
        mean_EPB_X = df_filtered['EPB_X'].mean()
        mean_EPB_Y = df_filtered['EPB_Y'].mean()
        mean_EPB_Z = df_filtered['EPB_Z'].mean()

        # Agregar resultados
        results.append({
            'pos_x': pos_x,
            'pos_y': pos_y,
            'pos_z': pos_z,
            'mean_Emod': mean_Emod,
            'mean_Esub_X': mean_Esub_X,
            'mean_Esub_Y': mean_Esub_Y,
            'mean_Esub_Z': mean_Esub_Z,
            'mean_EPB_X': mean_EPB_X,
            'mean_EPB_Y': mean_EPB_Y,
            'mean_EPB_Z': mean_EPB_Z
        })

        print(f"  Procesados {len(df_filtered)} pasos (de {len(df)} totales)")

    except Exception as e:
        print(f"Error procesando {csv_file}: {e}")
        continue

# Crear DataFrame con resultados
df_results = pd.DataFrame(results)

# Ordenar por posición x
df_results = df_results.sort_values('pos_x')

# Guardar en archivo CSV
output_file = os.path.join(base_dir, 'campo_electrico_medio.csv')
df_results.to_csv(output_file, index=False)
print(f"\nResultados guardados en: {output_file}")
print(f"Total de carpetas procesadas: {len(results)}")

# Calcular módulo del campo eléctrico a partir de las componentes medias
df_results['E_modulo_calculado'] = np.sqrt(
    df_results['mean_Esub_X']**2 +
    df_results['mean_Esub_Y']**2 +
    df_results['mean_Esub_Z']**2
)

# Crear gráficos
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))

# Subplot 1: Módulo del campo eléctrico vs posición x
ax1.plot(df_results['pos_x'], df_results['mean_Emod'], 'o-', linewidth=2, markersize=8, label='Emod (mean)')
ax1.plot(df_results['pos_x'], df_results['E_modulo_calculado'], 's--', linewidth=2, markersize=6, label='|E| calculated')
ax1.set_xlabel('Position x', fontsize=12)
ax1.set_ylabel('Eletric field module', fontsize=12)
ax1.set_title('Eletric field module vs position x', fontsize=14)
ax1.grid(True, alpha=0.3)
ax1.legend()

# Subplot 2: Componentes del campo eléctrico vs posición x
ax2.plot(df_results['pos_x'], df_results['mean_Esub_X'], 'o-', linewidth=2, markersize=8, label='E_x')
ax2.plot(df_results['pos_x'], df_results['mean_Esub_Y'], 's-', linewidth=2, markersize=8, label='E_y')
ax2.plot(df_results['pos_x'], df_results['mean_Esub_Z'], '^-', linewidth=2, markersize=8, label='E_z')
ax2.set_xlabel('Position x', fontsize=12)
ax2.set_ylabel('Eletric field components', fontsize=12)
ax2.set_title('Eletric field components vs position x', fontsize=14)
ax2.grid(True, alpha=0.3)
ax2.legend()
ax2.axhline(y=0, color='k', linestyle='-', linewidth=0.5, alpha=0.5)

plt.tight_layout()

# Guardar figura
plot_file = os.path.join(base_dir, 'campo_electrico_vs_x.png')
plt.savefig(plot_file, dpi=300, bbox_inches='tight')
print(f"Gráfico guardado en: {plot_file}")

if show_plot:
    plt.show()

print("\n¡Análisis completado!")
