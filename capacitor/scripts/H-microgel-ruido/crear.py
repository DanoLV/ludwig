#!/usr/bin/env python3
"""
H-microgel-ruido (2026-11-13, pedido del usuario): microgel de G-microgel-relax entre placas
a psi = 0, ahora CON ruido térmico (LB fluctuante + empujón aleatorio de las partículas
subgrid), en la escala de unidades de Ingo (escala "C"):

  kT = 1e-4, eta = 0.01, rho = 1 (como INGO-1PART-LBTERM), a_h = 0.3078881 (el de Ingo;
  al = 0.4743549 se deja como en G).  Para conservar el problema físico de G se re-escalan
  con kT:  epsilon = 1/(4 pi lB kT) con lB = 0.795775 (igual que G) -> 1000;
           bond_harmonic_k = 10 kT / LU^2 = 1e-3 (G: 1e-4 con kT 1e-5);
           soft_sphere_epsilon = kT = 1e-4 (G: 1e-5 = kT).
  Cocientes que se conservan: lB/r0 = 0.199, sqrt(kT/k)/r0 = 7.9 %, eps_rep/kT = 1, caja.

Prueba corta: 5000 pasos, posiciones de los coloides en CADA paso (para el desplazamiento
por paso). Cuatro corridas = {neutro, superficie} x {pared Steele de G re-escalada x10
(eps_w = 6.2572e-19), pared Steele estable (eps_w = 1e-5)}, ambas con sigma_w = 1.
Configuración inicial: config.cds.init.001-001 de G (neutro / superficie), con a_h cambiado.

Uso:  python3 capacitor/scripts/H-microgel-ruido/crear.py
Crea capacitor/local/WALL-DIRICHLET/H-microgel-ruido/<corrida>/ (no sobrescribe).
"""
import os
import re
import shutil
import sys

ROOT = "/home/danolv/Sim/ludwig"
G = f"{ROOT}/capacitor/local/WALL-DIRICHLET/G-microgel-relax"
H = f"{ROOT}/capacitor/local/WALL-DIRICHLET/H-microgel-ruido"
EXE = f"{ROOT}/src/Ludwig.exe"
AH = "3.078881242253297e-01"          # a_h de Ingo (INGO-1PART-LBTERM config.cds.init, línea 35)
WALLS = {"steeleG": ("6.2572e-19", "Steele de G re-escalada x10 (eps_w ~ kT-escalado); según el criterio "
                                   "mu*V'' < 1 es INESTABLE con ruido (mu*V'' ~ 40 donde V = kT): se prueba igual"),
         "steele1e-5": ("1.0e-5", "Steele estable: mu*V'' ~ 0.9 donde V = 5 kT; los monómeros quedan a "
                                  "d >~ 4.2 del plano sólido/fluido (con la de G llegaban a ~3.58)")}

HEADER = """# ===========================================================================
# H-microgel-ruido/{name}  (crear.py, 2026-11-13, pedido del usuario)
# Microgel de 20 monómeros de G-microgel-relax ({gel}) entre placas a psi = 0,
# CON RUIDO TÉRMICO, escala de unidades de Ingo (kT 1e-4, eta 0.01, rho 1, a_h 0.308).
# Re-escalado con kT para conservar el problema de G: epsilon 1e3 (lB = 0.796),
# k = 1e-3 (sqrt(kT/k)/r0 = 7.9 %), soft_sphere_epsilon = kT.
# Pared: {wall}.
# Prueba corta: 5000 pasos, posiciones en cada paso. Indicadores a evaluar:
# desplazamiento por paso vs drmax 0.8, equipartición de enlaces <U> ~ 35/2 kT,
# <u_a^2> = kT/rho del fluido, distancia mínima a las paredes.
# ===========================================================================
"""


def make_input(src, name, gel, wall_key):
    eps_w, wall_txt = WALLS[wall_key]
    rep = [
        (r"^N_start .*", "N_start  0"),
        (r"^N_cycles .*", "N_cycles 5000      # prueba corta"),
        (r"^fluid_rho0 .*", "fluid_rho0 1.0              # escala de Ingo"),
        (r"^viscosity .*", "viscosity      0.01        # escala de Ingo (tau_LB = 3 eta/rho + 1/2 = 0.53)"),
        (r"^viscosity_bulk .*", "viscosity_bulk 0.01"),
        (r"^lb_fluctuations .*", "lb_fluctuations on         # RUIDO TERMICO (LB + empujon de las particulas subgrid)"),
        (r"^temperature .*", "temperature 1.0e-4         # kT (escala de Ingo)"),
        (r"^electrokinetics_epsilon .*", "electrokinetics_epsilon     1.0e3          # lB = 1/(4 pi eps kT) = 0.796, igual que G"),
        (r"^electrokinetics_d0 .*", "electrokinetics_d0 0.05    # sin iones (rho_el 0): inerte. A0(D3Q27) = 6.138 -> D_fisico = 0.307"),
        (r"^electrokinetics_d1 .*", "electrokinetics_d1 0.05    # idem"),
        (r"^bond_harmonic_k .*", "bond_harmonic_k  1.0e-3     # 10 kT/LU^2, igual cociente que G"),
        (r"^soft_sphere_epsilon .*", "soft_sphere_epsilon 1.0e-4  # = kT, igual que G"),
        (r"^freq_statistics .*", "freq_statistics 100"),
        (r"^freq_psi_resid .*", "freq_psi_resid  500"),
        (r"^vel_io_freq .*", "vel_io_freq     1000       # campo de velocidad: <u_a^2> = kT/rho"),
        (r"^freq_config .*", "freq_config     5000"),
        (r"^colloid_io_freq .*", "colloid_io_freq 1          # posiciones en cada paso (desplazamiento por paso)"),
        (r"^pforce_io_freq .*", "pforce_io_freq  100"),
        (r"^wall_ss_cut_epsilon .*", f"wall_ss_cut_epsilon   {eps_w}"),
    ]
    out, seen = [HEADER.format(name=name, gel=gel, wall=wall_txt)], set()
    lines = open(src).read().splitlines()
    k = 0
    while k < len(lines) and (lines[k].startswith("#") or not lines[k].strip()):
        k += 1                                   # quitar la cabecera de G (no corresponde)
    for line in lines[k:]:
        for pat, new in rep:
            if re.match(pat, line):
                line = new
                seen.add(pat)
                break
        out.append(line)
    # Ludwig lee líneas de a 127 bytes (NKEY_LENGTH 128 en src/runtime.c): una línea más larga
    # se parte y el pedazo se lee como clave (aborta si tiene 1 sola palabra). Partir comentarios.
    wrapped = []
    for line in out:
        for sub in line.split("\n"):
            if len(sub.encode()) > 120 and sub.lstrip().startswith("#"):
                import textwrap
                wrapped += ["# " + x for x in textwrap.wrap(sub.lstrip()[1:].strip(), width=90)]
            else:
                wrapped.append(sub)
    out = wrapped
    missing = [p for p, _ in rep if p not in seen]
    assert not missing, missing
    return "\n".join(out) + "\n"


def make_config(src):
    L = open(src).read().split("\n")
    n = int(L[0].split()[0])
    for c in range(n):
        b = 1 + 68 * c
        assert abs(float(L[b + 33].split()[0]) - 0.05) < 1e-12, "a_h inesperado"
        L[b + 33] = f"{float(AH):24.15e}"
        assert int(float(L[b + 4].split()[0])) == 0, "isfixedr debe ser 0"
    return "\n".join(L)


def main():
    src_input = f"{G}/superficie-pared-steele/input"
    for gel in ("neutro", "superficie"):
        for wk in WALLS:
            name = f"{gel}-{wk}"
            d = f"{H}/{name}"
            if os.path.exists(d):
                sys.exit(f"{d} ya existe: no se sobrescribe")
            os.makedirs(f"{d}/logs")
            os.makedirs(f"{d}/proceced_data")
            open(f"{d}/input", "w").write(make_input(src_input, name, gel, wk))
            open(f"{d}/config.cds.init.001-001", "w").write(make_config(f"{G}/{gel}/config.cds.init.001-001"))
            shutil.copy2(f"{ROOT}/.petscrc", f"{d}/.petscrc")
            shutil.copy2(EXE, f"{d}/Ludwig.exe")
            print("creada", d)


if __name__ == "__main__":
    main()
