#!/usr/bin/env python3
"""
LONG-SALT-LD2LB-CEN-3E4 (2026-11-16, pedido del usuario): LONG-MASS-3E5 (parámetros de Ingo, una carga q = +1
libre, LB térmico) pero CON SAL 1:1 tal que lambda_D = 2 l_B (elección del usuario entre l_B, 2 l_B y 3 l_B),
advección centred (elección del usuario) y 3e4 pasos (prueba corta, elección del usuario).

Sal: lambda_D = (8 pi l_B n_s)^(-1/2) para sal 1:1 (fórmula que usa Ludwig para el "Debye length" del log)
     -> n_s = 1/(8 pi l_B lambda_D^2) = 1/(8 pi 2.84 5.68^2) = 4.3426e-4 por especie.
Los contraiones del +1 (1/27^3 = 5.08e-5) los agrega Ludwig por electroneutralidad: con ellos, lambda_D efectivo
(4 pi l_B sum z^2 n)^(-1/2) = 5.52 (no 5.68); Ludwig informa solo el de la sal.
Todo lo demás idéntico a LONG-MASS-3E5 (config.cds.init, semilla, kernel Hann-8, D de Ingo, salida cada 1000).
Ejecutable: el actual de src/ (269ec1b8: pesos advectivos + no-flujo en GPU).

Uso:  python3 microgel/scripts/INGO-1PART-LBTERM-DIAG/crear_sal.py
Crea microgel/local/INGO-1PART-LBTERM-DIAG/LONG-SALT-LD2LB-CEN-3E4/ (no sobrescribe).
"""
import math
import os
import re
import shutil
import sys

ROOT = "/home/danolv/Sim/ludwig"
BASE = f"{ROOT}/microgel/local/INGO-1PART-LBTERM-DIAG"
SRC = f"{BASE}/LONG-MASS-3E5"
DST = f"{BASE}/LONG-SALT-LD2LB-CEN-3E4"
EXE = f"{ROOT}/src/Ludwig.exe"

LB, LD = 2.84, 2 * 2.84
NS = 1.0 / (8 * math.pi * LB * LD ** 2)

HEADER = f"""##############################################################################
#
#  LONG-SALT-LD2LB-CEN-3E4  (crear_sal.py, 2026-11-16, pedido del usuario)
#  LONG-MASS-3E5 con SAL 1:1: lambda_D = 2 l_B = {LD:.2f} (l_B = {LB}),
#  n_s = 1/(8 pi l_B lambda_D^2) = {NS:.4e} por especie (+ contraiones del +1).
#  Adveccion centred (no flla_mass), 3e4 pasos (prueba corta). Resto igual.
#
##############################################################################
"""


def sub(t, pat, new):
    t2, n = re.subn(pat, new, t, flags=re.M)
    assert n == 1, pat
    return t2


def main():
    if os.path.exists(DST):
        sys.exit(DST + " ya existe: no se sobrescribe")
    t = open(f"{SRC}/input").read()
    t = sub(t, r"^N_cycles .*$", "N_cycles 30000            # prueba corta (eleccion del usuario)")
    t = sub(t, r"^electrokinetics_init_rho_el .*$",
            f"electrokinetics_init_rho_el {NS:.4e}  # sal 1:1 por especie: lambda_D = 2 l_B = {LD:.2f}")
    t = sub(t, r"^electrokinetics_advection .*$",
            "electrokinetics_advection centred   # pesos W_p = 2 wgradients (20261112) + no-flujo (20261116)")
    t = HEADER + "#  --- input original de LONG-MASS-3E5 (modificado arriba) ---\n" + t
    for l in t.split("\n"):
        assert len(l.encode()) <= 120, l        # NKEY_LENGTH 128 (src/runtime.c)
    os.makedirs(f"{DST}/proceced_data")
    open(f"{DST}/input", "w").write(t)
    for f in ("config.cds.init.001-001", "posicion_inicial.txt"):
        shutil.copy2(f"{SRC}/{f}", DST)
    shutil.copy2(EXE, f"{DST}/Ludwig.exe")
    print("creada", DST, "n_s = %.4e" % NS)


if __name__ == "__main__":
    main()
