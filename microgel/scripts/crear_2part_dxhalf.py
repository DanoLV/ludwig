#!/usr/bin/env python3
"""
Two fixed charges in the periodic box with dx = 1/2 (2026-11-09, user's request): ONLY separation d = 4, 20000 steps,
psi/efield every 2000 steps to see when the steady state is reached (user's choice). Same physical system as
SLIDES-LD8-hann08-corr_no-N10000 (dx = 1) and SLIDES-LD8-hann04-corr_no-dx2 (dx = 2), refined by r = 2 with the recipe of
crear_2part_dx2.py:
  - lengths x r: box 48 -> 96, Hann-8 -> Hann-16 (physical width 8)
  - electrokinetics_epsilon / r (= 5e3), electrokinetics_init_rho_el / r^3 (= 9.6485e-5): lambda_D = 16 nodes = 8 physical
  - D per step unchanged: the dx = 1 runs needed ~9000-10000 steps -> ~36000-40000 here; 20000 = t 5000 (user's choice)
Positions: same PHYSICAL positions, X_r = 0.5 + r (X_1 - 0.5): x = 8 -> 15.5, x = 8 + d -> 15.5 + 2d, y = z = 24 -> 47.5
(between nodes in the three coordinates; with dx = 1 both charges are on nodes).
Executable (user's choice): the CURRENT src/Ludwig.exe (md5 4f0c20e7..., 2026-11-04, with the halo double-count fix of
2026-09-26), NOT the one of the dx = 1 and dx = 2 runs (1c30df4b..., 2026-09-15).
Creates SLIDES-LD8-hann16-corr_no-dxhalf/r_04; an existing directory is not touched.
2026-11-10: optional arguments d and N (python3 crear_2part_dxhalf.py 5 12000 -> r_05 with 12000 steps); without
arguments it does exactly what created r_04 (d = 4, 20000 steps).
Usage: python3 crear_2part_dxhalf.py [d N]   and then   cd SLIDES-LD8-hann16-corr_no-dxhalf && ./correr.sh 1 > bateria.log
"""
import os, re, shutil, sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "SLIDES-LD8-hann08-corr_no-N10000")
BAT = os.path.join(HERE, "SLIDES-LD8-hann16-corr_no-dxhalf")
R = 2.0; NSTEPS = 20000; FREQ = 2000; DS = (4,)
if len(sys.argv) == 3:
    DS = (int(sys.argv[1]),); NSTEPS = int(sys.argv[2])

HEADER = """\
##############################################################################
#
#  SLIDES-LD8-hann16-corr_no-dxhalf/r_{d:02d} (crear_2part_dxhalf.py, 2026-11-09)
#  Dos cargas q = +1 fijas, caja periodica, dx = 1/2: mismo sistema fisico que
#  SLIDES-LD8-hann08-corr_no-N10000/pos_{x1:.2f}_24.00_24.00 (dx = 1), separacion d = {d} fisica.
#  Caja 96^3, Hann-16 (ancho fisico 8), epsilon 5e3, rho_el 9.6485e-5 (lambda_D = 16 nodos = 8 fisico).
#  Cargas en x = {xa:.2f} y {xb:.2f}, y = z = 47.50 (misma posicion fisica que dx = 1; entre nodos).
#  N = {n} (t = {t:.0f} en pasos de dx = 1), psi/efield cada 2000 pasos para ver cuando llega al estacionario.
#  Ejecutable actual (src/Ludwig.exe, 2026-11-04), no el de las corridas dx = 1 y 2.
#
##############################################################################
"""

def sub(t, pat, new):
    t2, n = re.subn(pat, new, t, flags=re.M)
    assert n == 1, pat
    return t2

if __name__ == "__main__":
    os.makedirs(BAT, exist_ok=True)
    for d in DS:
        x1 = 8.0 + d; src = os.path.join(SRC, "pos_%.2f_24.00_24.00" % x1); dst = os.path.join(BAT, "r_%02d" % d)
        if os.path.exists(dst):
            print("%s: ya existe, no se modifica" % dst); continue
        xa = 0.5 + R * (8.0 - 0.5); xb = 0.5 + R * (x1 - 0.5); yc = 0.5 + R * (24.0 - 0.5)
        t = open(os.path.join(src, "input")).read()
        t = HEADER.format(d=d, x1=x1, xa=xa, xb=xb, n=NSTEPS, t=NSTEPS / R ** 2) + t[t.index("\nN_start"):]
        t = sub(t, r"^N_cycles .*$", "N_cycles %d             # t = %.0f en pasos de dx = 1 (D por paso igual: pasos x r^2)" % (NSTEPS, NSTEPS / R ** 2))
        t = sub(t, r"^size .*$", "size 96_96_96          # 48 r")
        t = sub(t, r"^electrokinetics_epsilon .*$", "electrokinetics_epsilon 5.0e3   # 1e4/r: l_B = 1.5915494 nodos = 0.7957747 fisico")
        t = sub(t, r"^electrokinetics_init_rho_el .*$", "electrokinetics_init_rho_el 9.6485e-05   # 7.7188e-4/r^3: lambda_D = 16 nodos = 8 fisico")
        t = sub(t, r"^subgrid_hann_order .*$", "subgrid_hann_order 16      # orden Hann: ancho fisico 16/r = 8")
        t = sub(t, r"^pm_sr_sigma .*$", "pm_sr_sigma       2.8920      # no se usa (pm_sr_correction no); 1.4460 r")
        t = sub(t, r"^pm_sr_rcut .*$", "pm_sr_rcut        34.0      # no se usa (pm_sr_correction no); 17 r")
        t = sub(t, r"^vel_io_freq .*$", "vel_io_freq %d" % NSTEPS)
        t = sub(t, r"^freq_config .*$", "freq_config %d        # psi/efield cada 2000 pasos (t = 500): seguir la convergencia al estacionario" % FREQ)
        lines = open(os.path.join(src, "config.cds.init.001-001")).read().split("\n")
        assert [float(v) for v in lines[35].split()] == [x1, 24.0, 24.0], lines[35]
        assert [float(v) for v in lines[103].split()] == [8.0, 24.0, 24.0], lines[103]
        lines[35] = "".join("{:>25.15e}".format(v) for v in (xb, yc, yc))
        lines[103] = "".join("{:>25.15e}".format(v) for v in (xa, yc, yc))
        for sd in ("", "logs", "colloid_data", "proceced_data", "plots", "scripts"):
            os.makedirs(os.path.join(dst, sd), exist_ok=True)
        open(os.path.join(dst, "input"), "w").write(t)
        open(os.path.join(dst, "config.cds.init.001-001"), "w").write("\n".join(lines))
        shutil.copy(os.path.join(HERE, "..", "..", "src", "Ludwig.exe"), os.path.join(dst, "Ludwig.exe"))   # current exe
        print("%s: cargas en x = %.2f y %.2f, y = z = %.2f" % (dst, xa, xb, yc))
    t = open(os.path.join(SRC, "correr.sh")).read().replace("pos_*", "r_*")
    open(os.path.join(BAT, "correr.sh"), "w").write(t); os.chmod(os.path.join(BAT, "correr.sh"), 0o755)
