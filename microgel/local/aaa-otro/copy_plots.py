#!/usr/bin/env python3
#-------------------------------------------------------------------------------
import os
import shutil
import argparse

# Configurar argumentos de línea de comandos
parser = argparse.ArgumentParser(description='Copiar carpetas plots de origen a destino')
parser.add_argument('-origen', type=str, help='Carpeta origen')
parser.add_argument('-destino', type=str, help='Carpeta destino')
args = parser.parse_args()

origen = args.origen
destino = args.destino

os.makedirs(destino, exist_ok=True)

for hijo in os.listdir(origen):
    ruta_hijo_origen = os.path.join(origen, hijo)
    ruta_hijo_destino = os.path.join(destino, hijo)
    if os.path.isdir(ruta_hijo_origen):
        os.makedirs(ruta_hijo_destino, exist_ok=True)
        ruta_plots_origen = os.path.join(ruta_hijo_origen, "plots")
        ruta_plots_destino = os.path.join(ruta_hijo_destino, "plots")
        if os.path.isdir(ruta_plots_origen):
            # Copia todo el contenido de 'plots' al nuevo destino
            shutil.copytree(ruta_plots_origen, ruta_plots_destino, dirs_exist_ok=True, 
                            ignore=shutil.ignore_patterns('*.py', '*.csv'))