#!/usr/bin/env python3
"""
Two fixed charges, dx = 2, FIXED separation r = 4 (physical), both charges shifted together inside the cell
(2026-11-11, user's request: force vs intra-cell position of the pair).
  Shift delta (physical) = 0, 0.2, ..., 2.0 (one full cell of the dx = 2 lattice = 1 node; delta = 0 and 2.0 must be
  equivalent by translation invariance: control), along x only (user's choice). In lattice units both charges move by
  delta/2: charge 1 6.25 -> 6.25 + delta/2, charge 2 4.25 -> 4.25 + delta/2; y = z = 12.25 unchanged.
Copies of SLIDES-LD8-hann04-corr_no-dx2/r_04 (box 24^3, Hann-4, lambda_D = 8 physical, 3000 steps); only the two charge
positions change. Executable: the CURRENT src/Ludwig.exe (45dc94bc); delta = 0 is redone with it (r_04 used 1c30df4b).
Creates SLIDES-LD8-hann04-corr_no-dx2-r4-shift/d_X.X; existing directories are not touched.
Usage: python3 crear_2part_dx2_r4_shift.py   and then
       cd SLIDES-LD8-hann04-corr_no-dx2-r4-shift && ./correr.sh 4 > bateria.log
"""
import os, shutil

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "SLIDES-LD8-hann04-corr_no-dx2", "r_04")
BAT = os.path.join(HERE, "SLIDES-LD8-hann04-corr_no-dx2-r4-shift")
EXE = os.path.join(HERE, "..", "..", "src", "Ludwig.exe")
DELTAS = [round(0.2 * i, 1) for i in range(11)]          # 0.0 ... 2.0 physical
X1, X2, YZ = 6.25, 4.25, 12.25

HEADER = """\
##############################################################################
# SLIDES-LD8-hann04-corr_no-dx2-r4-shift/{name} (crear_2part_dx2_r4_shift.py, 2026-11-11, pedido del usuario)
# Copia de SLIDES-LD8-hann04-corr_no-dx2/r_04 (dx = 2, separacion fisica 4): las DOS cargas corridas en x
# delta = {d:.1f} fisico = {h:.2f} de red: x = {x2:.2f} y {x1:.2f} (antes 4.25 y 6.25), y = z = 12.25.
# Ejecutable ACTUAL (src/Ludwig.exe, md5 45dc94bc).
##############################################################################
"""

if __name__ == "__main__":
    os.makedirs(BAT, exist_ok=True)
    for d in DELTAS:
        name = "d_%.1f" % d; dst = os.path.join(BAT, name); h = d / 2.0
        if os.path.exists(dst):
            print("%s: ya existe, no se modifica" % dst); continue
        lines = open(os.path.join(SRC, "config.cds.init.001-001")).read().split("\n")
        assert [float(v) for v in lines[35].split()] == [X1, YZ, YZ], lines[35]
        assert [float(v) for v in lines[103].split()] == [X2, YZ, YZ], lines[103]
        lines[35] = "".join("{:>25.15e}".format(v) for v in (X1 + h, YZ, YZ))
        lines[103] = "".join("{:>25.15e}".format(v) for v in (X2 + h, YZ, YZ))
        for sd in ("", "logs", "colloid_data", "proceced_data", "plots", "scripts"):
            os.makedirs(os.path.join(dst, sd), exist_ok=True)
        t = HEADER.format(name=name, d=d, h=h, x1=X1 + h, x2=X2 + h) + open(os.path.join(SRC, "input")).read()
        open(os.path.join(dst, "input"), "w").write(t)
        open(os.path.join(dst, "config.cds.init.001-001"), "w").write("\n".join(lines))
        shutil.copy(EXE, os.path.join(dst, "Ludwig.exe"))
        print("%s: cargas en x = %.2f y %.2f" % (dst, X2 + h, X1 + h))
    t = open(os.path.join(HERE, "SLIDES-LD8-hann04-corr_no-dx2", "correr.sh")).read().replace("r_*", "d_*")
    open(os.path.join(BAT, "correr.sh"), "w").write(t); os.chmod(os.path.join(BAT, "correr.sh"), 0o755)
