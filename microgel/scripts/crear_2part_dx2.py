#!/usr/bin/env python3
"""
Two fixed charges in the periodic box with dx = 2 (2026-11-09, user's request): same physical system as
SLIDES-LD8-hann08-corr_no-N10000 (dx = 1, box 48, Hann-8, lambda_D = 8, q1 = q2 = +1), refined by r = 1/2 with the
recipe of capacitor/local/WALL-DIRICHLET/crear_refinado.py:
  - lengths x r: box 48 -> 24, Hann-8 -> Hann-4 (physical width 8)
  - electrokinetics_epsilon / r (= 2e4): q/(eps kT) in lattice units scales with r, psi(physical x) unchanged
  - electrokinetics_init_rho_el / r^3 (= 6.17504e-3): with eps/r gives lambda_D = 8 r = 4 nodes (8 physical)
  - D per step unchanged -> diffusion time in steps scales as r^2: N = 10000 r^2 = 2500, rounded up to 3000
  - freq_config 1000 r^2 = 250 (steady state: psi(3000) - psi(2750))
Positions: the same PHYSICAL positions as dx = 1 (user's choice). Ludwig node i is at coordinate i, the box spans
[0.5, N + 0.5], so the physical coordinate is X - 0.5 and X_r = 0.5 + r (X_1 - 0.5): the fixed charge at X_1 = 8 ->
4.25, the other at X_1 = 8 + d -> 4.25 + d/2, y = z = 24 -> 12.25. No charge sits on a node of the coarse lattice
(with dx = 1 both are on nodes); psi_electroneutral therefore removes no "near" nodes here (it removes the 2 occupied
nodes with dx = 1).
Executable: the one of SLIDES-LD8-hann08-corr_no (md5 1c30df4b..., 2026-09-15), as the dx = 1 runs.
Creates SLIDES-LD8-hann04-corr_no-dx2/r_DD (d = 1 ... 16); existing directories are not touched.
Usage: python3 crear_2part_dx2.py   and then   cd SLIDES-LD8-hann04-corr_no-dx2 && ./correr.sh 4 > bateria.log
"""
import os, re, shutil

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "SLIDES-LD8-hann08-corr_no-N10000")
BAT = os.path.join(HERE, "SLIDES-LD8-hann04-corr_no-dx2")
R = 0.5; NSTEPS = 3000; FREQ = 250

HEADER = """\
##############################################################################
#
#  SLIDES-LD8-hann04-corr_no-dx2/r_{d:02d} (crear_2part_dx2.py, 2026-11-09)
#  Dos cargas q = +1 fijas, caja periodica, dx = 2: mismo sistema fisico que
#  SLIDES-LD8-hann08-corr_no-N10000/pos_{x1:.2f}_24.00_24.00 (dx = 1), separacion d = {d} fisica.
#  Caja 24^3, Hann-4 (ancho fisico 8), epsilon 2e4, rho_el 6.17504e-3 (lambda_D = 4 nodos = 8 fisico).
#  Cargas en x = {xa:.2f} y {xb:.2f}, y = z = 12.25 (misma posicion fisica que dx = 1; ninguna sobre un nodo).
#  N = 3000 (10000 r^2 = 2500 redondeado), psi cada 250 pasos: estacionario con psi(3000) - psi(2750).
#
##############################################################################
"""

def sub(t, pat, new):
    t2, n = re.subn(pat, new, t, flags=re.M)
    assert n == 1, pat
    return t2

if __name__ == "__main__":
    os.makedirs(BAT, exist_ok=True)
    for d in range(1, 17):
        x1 = 8.0 + d; src = os.path.join(SRC, "pos_%.2f_24.00_24.00" % x1); dst = os.path.join(BAT, "r_%02d" % d)
        if os.path.exists(dst):
            print("%s: ya existe, no se modifica" % dst); continue
        xa = 0.5 + R * (8.0 - 0.5); xb = 0.5 + R * (x1 - 0.5); yc = 0.5 + R * (24.0 - 0.5)
        t = open(os.path.join(src, "input")).read()
        t = HEADER.format(d=d, x1=x1, xa=xa, xb=xb) + t[t.index("\nN_start"):]
        t = sub(t, r"^N_cycles .*$", "N_cycles %d             # 10000 r^2 = 2500 redondeado a 1000 hacia arriba (D por paso igual)" % NSTEPS)
        t = sub(t, r"^size .*$", "size 24_24_24          # 48 r")
        t = sub(t, r"^electrokinetics_epsilon .*$", "electrokinetics_epsilon 2.0e4   # 1e4/r: l_B = 0.3978874 nodos = 0.7957747 fisico")
        t = sub(t, r"^electrokinetics_init_rho_el .*$", "electrokinetics_init_rho_el 6.17504e-03   # 7.7188e-4/r^3: lambda_D = 4 nodos = 8 fisico")
        t = sub(t, r"^subgrid_hann_order .*$", "subgrid_hann_order 4      # orden Hann: ancho fisico 4/r = 8")
        t = sub(t, r"^pm_sr_rcut .*$", "pm_sr_rcut        8.5      # no se usa (pm_sr_correction no); 17 r")
        t = sub(t, r"^vel_io_freq .*$", "vel_io_freq %d" % NSTEPS)
        t = sub(t, r"^freq_config .*$", "freq_config %d        # 1000 r^2: estacionario con psi(N) - psi(N - %d)" % (FREQ, FREQ))
        lines = open(os.path.join(src, "config.cds.init.001-001")).read().split("\n")
        assert [float(v) for v in lines[35].split()] == [x1, 24.0, 24.0], lines[35]
        assert [float(v) for v in lines[103].split()] == [8.0, 24.0, 24.0], lines[103]
        lines[35] = "".join("{:>25.15e}".format(v) for v in (xb, yc, yc))
        lines[103] = "".join("{:>25.15e}".format(v) for v in (xa, yc, yc))
        for sd in ("", "logs", "colloid_data", "proceced_data"):
            os.makedirs(os.path.join(dst, sd), exist_ok=True)
        open(os.path.join(dst, "input"), "w").write(t)
        open(os.path.join(dst, "config.cds.init.001-001"), "w").write("\n".join(lines))
        shutil.copy(os.path.join(src, "Ludwig.exe"), os.path.join(dst, "Ludwig.exe"))
        print("%s: cargas en x = %.2f y %.2f, y = z = %.2f" % (dst, xa, xb, yc))
    t = open(os.path.join(SRC, "correr.sh")).read().replace("pos_*", "r_*")
    open(os.path.join(BAT, "correr.sh"), "w").write(t); os.chmod(os.path.join(BAT, "correr.sh"), 0o755)
