#!/usr/bin/env python3
"""
Corrida larga H-microgel-ruido (2026-11-13, pedido del usuario): neutro y superficie, pared
Steele eps_w = 1e-5 (la que resultó estable con ruido en la prueba de pared), 50000 pasos,
posiciones cada 100 pasos, volcados completos cada 10000 (permiten reiniciar / extender).
Copia de <gel>-steele1e-5 (misma configuración inicial, misma semilla) cambiando solo la salida
y la duración.

Uso:  python3 capacitor/scripts/H-microgel-ruido/crear_largo.py
Crea H-microgel-ruido/<gel>-largo/ (no sobrescribe).
"""
import os
import re
import shutil
import sys

H = "/home/danolv/Sim/ludwig/capacitor/local/WALL-DIRICHLET/H-microgel-ruido"
REP = [(r"^N_cycles .*$", "N_cycles 50000     # corrida larga (equilibrio por la serie de Rg)"),
       (r"^colloid_io_freq .*$", "colloid_io_freq 100        # posiciones (config.cds*)"),
       (r"^vel_io_freq .*$", "vel_io_freq     10000      # campo de velocidad: <u_a^2> = kT/rho"),
       (r"^freq_config .*$", "freq_config     10000      # volcados completos (reinicio)"),
       (r"^config_at_end .*$", "config_at_end   yes")]

for gel in ("neutro", "superficie"):
    src, d = f"{H}/{gel}-steele1e-5", f"{H}/{gel}-largo"
    if os.path.exists(d):
        sys.exit(f"{d} ya existe")
    os.makedirs(f"{d}/logs")
    os.makedirs(f"{d}/proceced_data")
    txt = open(f"{src}/input").read()
    # cabecera en dos líneas: Ludwig parte las líneas de más de 127 bytes (NKEY_LENGTH en src/runtime.c)
    txt = re.sub(r"^# H-microgel-ruido/.*$", f"# H-microgel-ruido/{gel}-largo (crear_largo.py, 2026-11-13): corrida larga.\n"
                 "# Pared Steele eps_w = 1e-5 (estable con ruido, ver README).", txt, count=1, flags=re.M)
    txt = re.sub(r"^# Prueba corta: 5000 pasos.*$", "# 50000 pasos, posiciones cada 100. Indicadores como en la prueba corta:", txt,
                 count=1, flags=re.M)
    for pat, new in REP:
        txt, n = re.subn(pat, new, txt, count=1, flags=re.M)
        assert n == 1, pat
    open(f"{d}/input", "w").write(txt)
    for f in ("config.cds.init.001-001", ".petscrc", "Ludwig.exe"):
        shutil.copy2(f"{src}/{f}", f"{d}/{f}")
    print("creada", d)
