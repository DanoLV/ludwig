#!/usr/bin/env python3
"""
Two fixed charges, intermediate separations d = 8.2, 8.4, ..., 9.8 for dx = 2 and dx = 1 (2026-11-11, user's request:
look for a fluctuation of the force with the position of the charge inside the cell; the dx = 2 error jumps between
integer d and dx = 1 does not).
  dx = 1: copies of SLIDES-LD8-hann08-corr_no-N10000/pos_16.00_24.00_24.00 (d = 8; box 48, Hann-8, 10000 steps);
          only the x of the moving charge changes: 16 -> 8 + d (the other stays at x = 8; y = z = 24).
          -> SLIDES-LD8-hann08-corr_no-N10000-d8a10/pos_XX.XX_24.00_24.00
  dx = 2: copies of SLIDES-LD8-hann04-corr_no-dx2/r_08 (d = 8; box 24, Hann-4, 3000 steps); the moving charge
          8.25 -> 4.25 + d/2 (same physical position as dx = 1; the other stays at 4.25; y = z = 12.25).
          -> SLIDES-LD8-hann04-corr_no-dx2-d8a10/r_XX.X
  Inside one coarse cell: dx = 2 moves the charge by 0.1 lattice units per point, dx = 1 by 0.2.
Executable: the CURRENT src/Ludwig.exe (md5 45dc94bc; with hydrodynamics on it gives results byte-identical to 4f0c20e7,
checked 2026-11-11 with TEST-HIDRO-dx2/hidro_si), as the user chose for d = 4.1 ... 4.9. The d = 8 and d = 10 runs
(old exe 1c30df4b) are not redone. Existing directories are not touched.
Usage: python3 crear_2part_d8a10.py   and then   ./correr_d8a10.sh > bateria_d8a10.log
"""
import os, shutil

HERE = os.path.dirname(os.path.abspath(__file__))
DS = [round(8.0 + 0.2 * i, 1) for i in range(1, 10)]   # 8.2 ... 9.8 (9.0 included: also redone with the current exe)
EXE = os.path.join(HERE, "..", "..", "src", "Ludwig.exe")
CASES = (  # dx, source run, destination dir, name format, moving-charge x(d), its x in the source, other charge, yz
    (1.0, "SLIDES-LD8-hann08-corr_no-N10000/pos_16.00_24.00_24.00", "SLIDES-LD8-hann08-corr_no-N10000-d8a10",
     lambda d: "pos_%.2f_24.00_24.00" % (8 + d), lambda d: 8.0 + d, 16.0, 8.0, 24.0),
    (2.0, "SLIDES-LD8-hann04-corr_no-dx2/r_08", "SLIDES-LD8-hann04-corr_no-dx2-d8a10",
     lambda d: "r_%04.1f" % d, lambda d: 4.25 + d / 2, 8.25, 4.25, 12.25))

HEADER = """\
##############################################################################
# {dst}/{name} (crear_2part_d8a10.py, 2026-11-11, pedido del usuario)
# Copia de {src}: solo cambia la x de la carga movil, {x0:.2f} -> {x1:.2f}
# (separacion fisica d = {d:.1f}, dx = {dx:g}; la otra carga sigue en x = {xo:.2f}, y = z = {yz:.2f}).
# Ejecutable ACTUAL (src/Ludwig.exe, md5 45dc94bc), no el 1c30df4b de la corrida original.
##############################################################################
"""

if __name__ == "__main__":
    for dx, src, dstdir, fname, fx, x0, xo, yz in CASES:
        src = os.path.join(HERE, src); B = os.path.join(HERE, dstdir); os.makedirs(B, exist_ok=True)
        for d in DS:
            name = fname(d); dst = os.path.join(B, name); x1 = fx(d)
            if os.path.exists(dst):
                print("%s: ya existe, no se modifica" % dst); continue
            lines = open(os.path.join(src, "config.cds.init.001-001")).read().split("\n")
            assert [float(v) for v in lines[35].split()] == [x0, yz, yz], lines[35]
            assert [float(v) for v in lines[103].split()] == [xo, yz, yz], lines[103]
            lines[35] = "".join("{:>25.15e}".format(v) for v in (x1, yz, yz))
            for sd in ("", "logs", "colloid_data", "proceced_data", "plots", "scripts"):
                os.makedirs(os.path.join(dst, sd), exist_ok=True)
            t = HEADER.format(dst=dstdir, name=name, src=os.path.relpath(src, HERE), x0=x0, x1=x1, d=d, dx=dx, xo=xo, yz=yz)
            open(os.path.join(dst, "input"), "w").write(t + open(os.path.join(src, "input")).read())
            open(os.path.join(dst, "config.cds.init.001-001"), "w").write("\n".join(lines))
            shutil.copy(EXE, os.path.join(dst, "Ludwig.exe"))
            print("%s: carga movil en x = %.2f" % (dst, x1))
