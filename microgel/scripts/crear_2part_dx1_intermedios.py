#!/usr/bin/env python3
"""
Two fixed charges, dx = 1, intermediate separations d = 4.1, 4.2, ..., 4.9 (2026-11-11, user's request: points every 0.1
between d = 4 and 5). Copies of SLIDES-LD8-hann08-corr_no-N10000/pos_12.00_24.00_24.00 (box 48^3, lambda_D = 8, Hann-8,
10000 steps, freq_config 1000) where ONLY the x of the moving charge changes: 12 -> 8 + d (the other stays at x = 8;
y = z = 24). With d non-integer the moving charge is between nodes in x.
Executable (user's choice): the CURRENT src/Ludwig.exe (md5 4f0c20e7, 2026-11-04, with the halo fix of 2026-09-26), NOT
the 1c30df4b of the d = 4 and d = 5 runs of SLIDES-LD8-hann08-corr_no-N10000 (those are not redone).
Creates SLIDES-LD8-hann08-corr_no-N10000-d4a5/pos_XX.XX_24.00_24.00; an existing directory is not touched.
Usage: python3 crear_2part_dx1_intermedios.py   and then   cd SLIDES-LD8-hann08-corr_no-N10000-d4a5 && ./correr.sh 3 > bateria.log
"""
import os, shutil

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "SLIDES-LD8-hann08-corr_no-N10000")
BAT = SRC + "-d4a5"
DS = [round(4.0 + 0.1 * i, 1) for i in range(1, 10)]

HEADER = """\
##############################################################################
# SLIDES-LD8-hann08-corr_no-N10000-d4a5/{name} (crear_2part_dx1_intermedios.py, 2026-11-11, pedido del usuario)
# Copia de SLIDES-LD8-hann08-corr_no-N10000/pos_12.00_24.00_24.00: solo cambia la x de la carga movil,
# 12 -> {x1:.2f} (separacion d = {d:.1f}; la otra carga sigue en x = 8, y = z = 24).
# Ejecutable ACTUAL (src/Ludwig.exe, md5 4f0c20e7, 2026-11-04), no el 1c30df4b de las corridas originales.
##############################################################################
"""

if __name__ == "__main__":
    src = os.path.join(SRC, "pos_12.00_24.00_24.00")
    os.makedirs(BAT, exist_ok=True)
    for d in DS:
        x1 = 8.0 + d; name = "pos_%.2f_24.00_24.00" % x1; dst = os.path.join(BAT, name)
        if os.path.exists(dst):
            print("%s: ya existe, no se modifica" % dst); continue
        lines = open(os.path.join(src, "config.cds.init.001-001")).read().split("\n")
        assert [float(v) for v in lines[35].split()] == [12.0, 24.0, 24.0], lines[35]
        assert [float(v) for v in lines[103].split()] == [8.0, 24.0, 24.0], lines[103]
        lines[35] = "".join("{:>25.15e}".format(v) for v in (x1, 24.0, 24.0))
        for sd in ("", "logs", "colloid_data", "proceced_data", "plots", "scripts"):
            os.makedirs(os.path.join(dst, sd), exist_ok=True)
        t = HEADER.format(name=name, x1=x1, d=d) + open(os.path.join(src, "input")).read()
        open(os.path.join(dst, "input"), "w").write(t)
        open(os.path.join(dst, "config.cds.init.001-001"), "w").write("\n".join(lines))
        shutil.copy(os.path.join(HERE, "..", "..", "src", "Ludwig.exe"), os.path.join(dst, "Ludwig.exe"))   # current exe
        print("%s: cargas en x = 8.00 y %.2f" % (dst, x1))
    shutil.copy(os.path.join(SRC, "correr.sh"), os.path.join(BAT, "correr.sh"))
