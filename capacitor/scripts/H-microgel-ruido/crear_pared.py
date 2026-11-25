#!/usr/bin/env python3
"""
Prueba de pared con ruido térmico (H-microgel-ruido, 2026-11-13, pedido del usuario).

Un monómero libre, empezando a d = 4.0 del plano sólido/fluido inferior (z = 5.5; x = y = 16.5),
mismas unidades y parámetros que las corridas H (escala de Ingo: kT 1e-4, eta 0.01, rho 1,
a_h 0.308, al 0.474; epsilon 1e3), sin sal, placas a psi = 0. 2000 pasos, posición en cada paso.
Cuatro corridas: carga q = 0 y q = +1 (la imagen lo atrae a la pared) x las dos paredes Steele
(eps_w = 6.2572e-19 "steeleG" y 1e-5 "steele1e-5", sigma_w = 1).
Base: el input de H-microgel-ruido/superficie-<pared> y el bloque del monómero 1 de su
config.cds.init (q = 1, a_h ya cambiado), sin enlaces.

Uso:  python3 capacitor/scripts/H-microgel-ruido/crear_pared.py [z0]   (z0 = 5.5 por defecto)
Nota: con eps_w = 1e-5, z0 = 5.5 (h = 0.5) cae dentro de la zona repulsiva y aborta en el paso 1:
usar z0 = 6.5 (h = 1.5, fuera del alcance z* = 0.987).
Crea H-microgel-ruido/pared-q{0,1}-<pared>/ (no sobrescribe).
"""
import os
import re
import shutil
import sys

H = "/home/danolv/Sim/ludwig/capacitor/local/WALL-DIRICHLET/H-microgel-ruido"
Z0 = float(sys.argv[1]) if len(sys.argv) > 1 else 5.5   # 5.5 -> d = 4.0; 6.5 -> d = 5.0
SUF = "" if Z0 == 5.5 else f"-d{Z0 - 1.5:g}"


def main():
    for wall in ("steeleG", "steele1e-5"):
        src = f"{H}/superficie-{wall}"
        L = open(f"{src}/config.cds.init.001-001").read().split("\n")
        block = L[1:69]
        assert abs(float(block[45].split()[0]) - 1.0) < 1e-12, "el monómero 1 debe tener q = 1"
        for q in (0, 1):
            name = f"pared-q{q}-{wall}{SUF}"
            d = f"{H}/{name}"
            if os.path.exists(d):
                sys.exit(f"{d} ya existe")
            os.makedirs(f"{d}/logs")
            os.makedirs(f"{d}/proceced_data")
            b = list(block)
            b[0] = f"{1:24d}"                               # index
            b[2] = f"{0:24d}"                               # nbonds
            for j in range(9, 15):                          # bond[]
                b[j] = f"{0:24d}"
            b[34] = f"{16.5:24.15e} {16.5:24.15e} {Z0:24.15e} "
            b[45] = f"{float(q):24.15e}"                    # q0
            open(f"{d}/config.cds.init.001-001", "w").write("\n".join([f"{1:24d}"] + b) + "\n")
            txt = open(f"{src}/input").read()
            txt = re.sub(r"^# H-microgel-ruido/.*$", f"# H-microgel-ruido/{name}  (crear_pared.py, 2026-11-13): UN monómero "
                         f"libre q = {q} empezando a d = {Z0 - 1.5:g} del plano inferior (z = {Z0:g}); prueba de la pared con ruido.",
                         txt, count=1, flags=re.M)
            txt = re.sub(r"^N_cycles .*$", "N_cycles 2000      # prueba de pared", txt, count=1, flags=re.M)
            open(f"{d}/input", "w").write(txt)
            shutil.copy2(f"{src}/.petscrc", f"{d}/.petscrc")
            shutil.copy2(f"{src}/Ludwig.exe", f"{d}/Ludwig.exe")
            print("creada", d)


if __name__ == "__main__":
    main()
