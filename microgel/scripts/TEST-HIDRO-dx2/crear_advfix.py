#!/usr/bin/env python3
"""
Rehace el barrido de viscosidad de TEST-HIDRO-dx2 (7 corridas, 2026-11-10) con el ejecutable ACTUAL
(src/Ludwig.exe, md5 221f3229, compilado 2026-11-12 19:35), que incluye el arreglo del peso de la
advección de NP (CHANGE 20261112: W_p = 2 wgradients[p]; antes los iones se advectaban 9x con stencil 27).
Pedido del usuario 2026-11-16: comparar el potencial químico con hidrodinámica en función de la
viscosidad antes y después del arreglo.

Cada corrida nueva es una copia del input y de config.cds.init de la original, sin otro cambio físico
(misma advección "centred" por defecto, mismos 3000 pasos, mismo solver fft); solo se agregan líneas de
comentario en la cabecera.

Uso:  python3 scripts/crear_advfix.py   (desde TEST-HIDRO-dx2). Crea advfix/<corrida>/ (no sobrescribe).
"""
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = "/home/danolv/Sim/ludwig/microgel/local/TEST-HIDRO-dx2"   # scripts/ es symlink: no usar HERE/..
EXE = "/home/danolv/Sim/ludwig/src/Ludwig.exe"
RUNS = ("hidro_si_eta0.02", "hidro_si", "hidro_si_eta0.1", "hidro_si_eta0.2", "hidro_si_eta0.5",
        "hidro_si_eta1.0", "hidro_si_eta2.0")
HEAD = ("#  advfix (2026-11-16, crear_advfix.py): copia de ../../{run} con el exe actual (md5 221f3229),\n"
        "#  que incluye el peso W_p = 2 wgradients[p] de la advección de NP (CHANGE 20261112).\n"
        "#  Único cambio respecto de la original: el ejecutable. Comparar Delta mu antes/después.\n")

for run in RUNS:
    src, dst = os.path.join(BASE, run), os.path.join(BASE, "advfix", run)
    if os.path.exists(dst):
        sys.exit(dst + " ya existe")
    os.makedirs(dst)
    for sub in ("logs", "proceced_data", "colloid_data", "plots"):
        os.makedirs(os.path.join(dst, sub))
    txt = open(os.path.join(src, "input")).read()
    lines = txt.split("\n")
    assert lines[0].startswith("####")
    txt = "\n".join([lines[0], HEAD.format(run=run).rstrip("\n")] + lines[1:])
    for l in txt.split("\n"):
        assert len(l.encode()) <= 120, l          # NKEY_LENGTH 128 en src/runtime.c
    open(os.path.join(dst, "input"), "w").write(txt)
    shutil.copy2(os.path.join(src, "config.cds.init.001-001"), dst)
    shutil.copy2(EXE, os.path.join(dst, "Ludwig.exe"))
    print("creada", dst)
