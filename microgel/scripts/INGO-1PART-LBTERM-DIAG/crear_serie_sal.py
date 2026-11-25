#!/usr/bin/env python3
"""
SAL-SERIE (2026-11-16, pedido del usuario): cómo cambia la nube de iones alrededor de la carga q = +1 con sal
(lambda_D = 2 l_B) al variar D_ion, la temperatura y la resolución, y si la nube sigue a la partícula.
Base: LONG-SALT-LD2LB-CEN-3E4 (parámetros de Ingo, 27^3, kT 1e-4, eta 0.01, l_B 2.84, n_s 4.3426e-4,
lambda_D 5.68, Hann-8, advección centred, LB térmico, partícula libre, 3e4 pasos). Esa corrida ES el caso
BASE-LIBRE y no se repite aquí.

Decisiones del usuario (2026-11-16):
  * D_ion x10 y x30; la partícula conserva su movilidad desnuda mu_b = 1/gamma (a_h, a_L sin cambio).
  * kT/4 y kT/16 con epsilon ~ 1/kT: l_B = 1/(4 pi eps kT) = 2.84 fijo, misma sal y lambda_D; solo bajan el ruido
    del LB, el empujón de la partícula y D_p = kT mu_b. D_ion no cambia (Ludwig usa D y el potencial en kT).
  * dx/2, 1.5e4 pasos físicos: mismo sistema físico con escalado difusivo del tiempo (dt/4):
        caja 54^3; Hann-16 (mismo ancho físico); D_LU y nu_LU iguales (d, viscosity sin cambio, tau_LB 0.53);
        kT_LU x2 (unidad de energía M dx^2/dt^2 con M = rho dx^3 -> /2); l_B = 5.68 nodos -> eps = eps/4 = 70.0506;
        n_s/8 = 5.4283e-5 (lambda_D 11.36 nodos); mu_b en LU /2 = 3.2575 con a_L = 0.4950778 fijo ->
        a_h = 1/(6 pi eta mu_b + 1/a_L) = 0.379664 (el a_h físico no se puede conservar: daría mu_b < 0);
        posición física igual: X' = 0.5 + 2 (X - 0.5); 6e4 pasos, salida cada 4000.
  * Gemelas: partícula FIJA con LB térmico (separa "la nube no sigue a la partícula" de "el ruido deforma la
    nube") en todos los casos, y FIJA SIN RUIDO (lb_fluctuations off) solo en la base (valida la PB con sal).
Todo lo demás idéntico a la base (misma semilla del ruido, mismo ejecutable que la base: src/Ludwig.exe 269ec1b8).
Duración 3e4 pasos (como la base). Ojo: con D_ion x1 el tiempo de Debye lambda_D^2/D ~ 4.7e4 pasos es mayor.

Uso:  python3 microgel/scripts/INGO-1PART-LBTERM-DIAG/crear_serie_sal.py      (crea; NO corre)
      bash microgel/scripts/INGO-1PART-LBTERM-DIAG/correr_serie_sal.sh [caso ...]   (para correr, cuando se decida)
Crea microgel/local/INGO-1PART-LBTERM-DIAG/SAL-SERIE/<caso>/ (no sobrescribe).
"""
import math
import os
import re
import shutil
import sys

ROOT = "/home/danolv/Sim/ludwig"
DIAG = f"{ROOT}/microgel/local/INGO-1PART-LBTERM-DIAG"
SRC = f"{DIAG}/LONG-SALT-LD2LB-CEN-3E4"
DST = f"{DIAG}/SAL-SERIE"
EXE = f"{SRC}/Ludwig.exe"           # mismo ejecutable que la base

D0, EPS0, KT0, NS0 = 1.0614504e-4, 280.2024, 1.0e-4, 4.3426e-4
A0 = 6.137828
AL = 0.4950778
ETA = 0.01


def ah_from_mub(mub):
    return 1.0 / (6 * math.pi * ETA * mub + 1.0 / AL)


# caso: (descripción, dict de cambios)
CASES = {
    "BASE-FIJA":         ("base, particula fija, LB termico", dict(fixed=True)),
    "BASE-FIJA-NOFLUCT": ("base, particula fija, sin ruido (valida PB con sal)", dict(fixed=True, nofluct=True)),
    "D10-LIBRE":  ("D_ion x10", dict(dfac=10)),
    "D10-FIJA":   ("D_ion x10, particula fija", dict(dfac=10, fixed=True)),
    "D30-LIBRE":  ("D_ion x30", dict(dfac=30)),
    "D30-FIJA":   ("D_ion x30, particula fija", dict(dfac=30, fixed=True)),
    "T4-LIBRE":   ("kT/4, eps x4 (l_B fijo)", dict(tfac=4)),
    "T4-FIJA":    ("kT/4, eps x4, particula fija", dict(tfac=4, fixed=True)),
    "T16-LIBRE":  ("kT/16, eps x16 (l_B fijo)", dict(tfac=16)),
    "T16-FIJA":   ("kT/16, eps x16, particula fija", dict(tfac=16, fixed=True)),
    "DXH-LIBRE":  ("dx/2 (mismo sistema fisico), 1.5e4 pasos fisicos", dict(dxh=True)),
    "DXH-FIJA":   ("dx/2, particula fija", dict(dxh=True, fixed=True)),
}

HEADER = """##############################################################################
#
#  SAL-SERIE/{case}  (crear_serie_sal.py, 2026-11-16, pedido del usuario)
#  {desc}.
#  Base: LONG-SALT-LD2LB-CEN-3E4 (Ingo + sal, lambda_D = 2 l_B, centred).
#  Cambios respecto de la base: {chg}
#
##############################################################################
"""


def sub(t, pat, new):
    t2, n = re.subn(pat, new, t, flags=re.M)
    assert n == 1, pat
    return t2


def make(case, desc, o):
    t = open(f"{SRC}/input").read()
    t = t[t.index("#  --- input original"):]          # quitar la cabecera de la base
    chg = []
    if o.get("dfac"):
        f = o["dfac"]
        d = D0 * f
        t = sub(t, r"^electrokinetics_d0 .*$", f"electrokinetics_d0          {d:.7e}  # D_ion x{f}: D_fisico = d*A0 = {d*A0:.4e}")
        t = sub(t, r"^electrokinetics_d1 .*$", f"electrokinetics_d1          {d:.7e}  # A0(D3Q27) = {A0}")
        chg.append(f"d0 = d1 = {d:.7e} (x{f})")
    if o.get("tfac"):
        f = o["tfac"]
        kt, eps = KT0 / f, EPS0 * f
        t = sub(t, r"^temperature .*$", f"temperature     {kt:.6e}   # kT/{f} (ruido del LB y de la particula menor)")
        t = sub(t, r"^electrokinetics_epsilon .*$",
                f"electrokinetics_epsilon     {eps:.4f} # x{f}: l_B = 1/(4 pi eps kT) = 2.84 fijo")
        chg.append(f"temperature {kt:.6e}, epsilon {eps:.4f}")
    if o.get("dxh"):
        mub = 6.515 / 2      # mu_b(LU) con dx/2, dt/4, M/8: x (M/8)/(dt/4) = /2
        ah = ah_from_mub(mub)
        t = sub(t, r"^N_cycles .*$", "N_cycles 60000            # 1.5e4 pasos fisicos x 4 (dt/4)")
        t = sub(t, r"^size .*$", "size 54_54_54            # 27 x 2 (dx/2)")
        t = sub(t, r"^temperature .*$", "temperature     2.0e-4   # kT_LU x2 con dx/2, dt/4, M = rho dx^3")
        t = sub(t, r"^electrokinetics_epsilon .*$",
                "electrokinetics_epsilon     70.0506  # eps/4: l_B = 1/(4 pi eps kT) = 5.68 nodos = 2.84 fisico")
        t = sub(t, r"^electrokinetics_init_rho_el .*$",
                f"electrokinetics_init_rho_el {NS0/8:.4e}  # n_s/8: lambda_D = 11.36 nodos = 5.68 fisico")
        t = sub(t, r"^subgrid_hann_order .*$", "subgrid_hann_order 16     # mismo ancho fisico que Hann-8 con dx = 1")
        for key in ("freq_statistics", "colloid_io_freq", "psi_io_freq", "rho_io_freq"):
            t = sub(t, rf"^{key}( +)1000\b", rf"{key}\g<1>4000")
        chg.append(f"dx/2: size 54, Hann-16, kT 2e-4, eps 70.0506, rho_el {NS0/8:.4e}, a_h {ah:.6f}, 6e4 pasos")
        o["ah"] = ah
    if o.get("nofluct"):
        t = sub(t, r"^lb_fluctuations .*$", "lb_fluctuations off      # sin ruido termico (ni LB ni empujon de la particula)")
        chg.append("lb_fluctuations off")
    if o.get("fixed"):
        chg.append("particula fija (isfixedr = 1)")
    t = HEADER.format(case=case, desc=desc, chg="\n#    ".join([""] + chg)) + t
    for line in t.split("\n"):
        assert len(line.encode()) <= 120, line       # NKEY_LENGTH 128 (src/runtime.c)

    L = open(f"{SRC}/config.cds.init.001-001").read().split("\n")
    assert int(L[0].split()[0]) == 1
    assert int(float(L[5])) == 0 and abs(float(L[34]) - 0.3078881242253297) < 1e-12
    if o.get("fixed"):
        L[5] = f"{1:24d}"
    if o.get("dxh"):
        L[34] = f"{o['ah']:24.15e}"
        x = [float(v) for v in L[35].split()]
        L[35] = "".join(f"{0.5 + 2 * (c - 0.5):24.15e}" for c in x)
    return t, "\n".join(L)


def main():
    for case, (desc, o) in CASES.items():
        d = f"{DST}/{case}"
        if os.path.exists(d):
            print(f"{d}: ya existe, no se modifica")
            continue
        t, cfg = make(case, desc, dict(o))
        os.makedirs(f"{d}/proceced_data")
        open(f"{d}/input", "w").write(t)
        open(f"{d}/config.cds.init.001-001", "w").write(cfg)
        shutil.copy2(EXE, f"{d}/Ludwig.exe")
        print("creada", d)


if __name__ == "__main__":
    main()
