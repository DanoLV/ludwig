#!/usr/bin/env python3
"""
Continuacion de C-electrolito-hann8-LD6-desplazada (crear_desplazada.py) desde el paso 5000 hasta el 10000, para
verificar el estacionario. Las corridas originales no se tocan: cada z_* se copia entero (con las salidas del paso 5000:
dist, rho, vel, psi, qsi, config.cds) a C-electrolito-hann8-LD6-desplazada-cont/, y en la copia el input pasa a
  N_start 5000, N_cycles 5000 (pasos adicionales), freq_config 1000 (psi/efield cada 1000 pasos).
Ludwig reanuda leyendo esos archivos (ludwig.c: "Re-starting simulation at step 5000"; colloids_rt.c lee
config.cds00005000). Probado el 2026-11-05 con 1 paso: max|psi(5001) - psi(5000)| = 2e-15, y los archivos del paso 5000
no se modifican.
correr_cont.sh NO borra salidas previas (correr.sh si: borraria los archivos del paso 5000 que hacen falta para
reanudar). Escribe run_cont.log.
Las posiciones ya continuadas se saltean (se pueden agregar posiciones nuevas despues).
Uso: python3 crear_continuacion.py   y despues   cd C-electrolito-hann8-LD6-desplazada-cont && bash correr_cont.sh 3 > bateria.log
"""
import os, shutil

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "C-electrolito-hann8-LD6-desplazada")
BAT = SRC + "-cont"
N0, NADD = 5000, 5000

RUNNER = """\
#!/bin/bash
# Reanuda todas las posiciones z_* desde el paso {n0} (N en paralelo). No borra nada; log en run_cont.log.
export LD_LIBRARY_PATH=/usr/local/petsc-cuda-hypre/lib:/usr/local/ompi/lib:$LD_LIBRARY_PATH
cd "$(dirname "$(readlink -f "$0")")"
N=${{1:-3}}
uno() {{
  cd "$1" || return 1
  if grep -aqs "finished normally" run_cont.log; then echo "[SKIP] $1"; return 0; fi
  stdbuf -o0 ./Ludwig.exe > run_cont.log 2>&1
  grep -aq "finished normally" run_cont.log && grep -aq "Re-starting simulation at step {n0}" run_cont.log \\
    && echo "[OK] $1" || echo "[FALLO] $1"
}}
export -f uno
ls -d z_* | xargs -P "$N" -I{{}} bash -c 'uno {{}}'
echo BATERIA_LISTA
"""

for z in sorted(d for d in os.listdir(SRC) if d.startswith("z_")):
    src, dst = os.path.join(SRC, z), os.path.join(BAT, z)
    if os.path.exists(dst):                           # continuaciones ya creadas: no se tocan
        print("%s: ya existe, no se modifica" % dst); continue
    assert os.path.exists(os.path.join(src, "psi-%09d.001-001" % N0)), src
    shutil.copytree(src, dst, symlinks=True)
    t = open(os.path.join(dst, "input")).read()
    for old, new in (("\nN_start 0\n", "\nN_start %d            # continuacion desde el paso %d (crear_continuacion.py)\n" % (N0, N0)),
                     ("\nN_cycles 5000\n", "\nN_cycles %d            # pasos adicionales: termina en el paso %d\n" % (NADD, N0 + NADD)),
                     ("\nfreq_config 5000\n", "\nfreq_config 1000        # psi/efield cada 1000 pasos: deriva psi(N) - psi(N - 1000)\n")):
        assert t.count(old) == 1, old
        t = t.replace(old, new)
    open(os.path.join(dst, "input"), "w").write(t)
    print("%s -> %s" % (src, dst))
open(os.path.join(BAT, "correr_cont.sh"), "w").write(RUNNER.format(n0=N0))
