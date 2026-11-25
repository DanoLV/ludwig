#!/usr/bin/env python3
"""
Barrido de viscosidad de TEST-HIDRO-dx2 (7 corridas, 2026-11-10) con el ejecutable actual (md5 221f3229)
y advección de NP "flla_mass" (variante de flla de DeVincentis & Thomas 2023, limitador superbee; CHANGE 20261112 en
src/nernst_planck.c) en lugar de "centred". Pedido del usuario 2026-11-16: comparar Delta mu vs viscosidad
con la advección vieja, la centrada corregida (advfix/), flla (advfix_flla/) y flla_mass.

Cada corrida es una copia del input y de config.cds.init de la original; cambian el ejecutable y la
línea agregada "electrokinetics_advection flla_mass" (requisitos: nhalo >= 2 (aquí 3), un solo proceso MPI,
enlaces D3Q19 del LB contenidos en el stencil 27 de NP).

Uso:  python3 scripts/crear_advfix_flla_mass.py   (desde TEST-HIDRO-dx2). Crea advfix_flla_mass/<corrida>/ (no sobrescribe).
"""
import os
import shutil
import sys

BASE = "/home/danolv/Sim/ludwig/microgel/local/TEST-HIDRO-dx2"   # scripts/ es symlink: rutas absolutas
EXE = "/home/danolv/Sim/ludwig/src/Ludwig.exe"
RUNS = ("hidro_si_eta0.02", "hidro_si", "hidro_si_eta0.1", "hidro_si_eta0.2", "hidro_si_eta0.5",
        "hidro_si_eta1.0", "hidro_si_eta2.0")
HEAD = ("#  advfix_flla_mass (2026-11-16, crear_advfix_flla_mass.py): copia de ../../{run} con el exe actual\n"
        "#  (md5 221f3229) y advección de NP flla_mass (superbee) en vez de centred. Comparar Delta mu.\n")
ADV = ("\n# Advección de NP (agregado por crear_advfix_flla_mass.py): flujo de masa LB sin normalizar + superbee\n"
       "electrokinetics_advection flla_mass\n")

for run in RUNS:
    src, dst = os.path.join(BASE, run), os.path.join(BASE, "advfix_flla_mass", run)
    if os.path.exists(dst):
        sys.exit(dst + " ya existe")
    txt = open(os.path.join(src, "input")).read()
    assert "electrokinetics_advection" not in txt
    lines = txt.split("\n")
    assert lines[0].startswith("####")
    txt = "\n".join([lines[0], HEAD.format(run=run).rstrip("\n")] + lines[1:]).rstrip("\n") + "\n" + ADV
    for l in txt.split("\n"):
        assert len(l.encode()) <= 120, l          # NKEY_LENGTH 128 en src/runtime.c
    os.makedirs(dst)
    for sub in ("logs", "proceced_data", "colloid_data", "plots"):
        os.makedirs(os.path.join(dst, sub))
    open(os.path.join(dst, "input"), "w").write(txt)
    shutil.copy2(os.path.join(src, "config.cds.init.001-001"), dst)
    shutil.copy2(EXE, os.path.join(dst, "Ludwig.exe"))
    print("creada", dst)
