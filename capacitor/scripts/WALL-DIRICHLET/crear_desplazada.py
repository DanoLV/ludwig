#!/usr/bin/env python3
"""
Bateria C-electrolito-hann8-LD6-desplazada: el mismo sistema que C-electrolito-hann8-LD6/z_17.50 (Hann-8, lambda_D = 6,
dx = 1, placas a tierra, q = +1 fija en x = y = 16), con la carga en el centro del canal y corrida en z hacia la placa
superior por fracciones de celda:
  z = 17.50 + delta,  delta = 0, 0.10, 0.25, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90
(z = 18.50, delta = 1, seria de nuevo un punto medio entre nodos, simetrico; no se corre). Las corridas que ya
existen no se recrean ni se modifican. Las de delta > 0.5 (agregadas despues) guardan psi/efield cada 1000 pasos
(freq_config 1000, solo cambia la salida) para verificar el estacionario en el paso 5000 con psi(5000) - psi(4000).
Acercamiento a la placa (2026-11-05): 10 posiciones de delta = 1 a 12, equiespaciadas y redondeadas a 0.1 (= dx/10):
  delta = 1.0, 2.2, 3.4, 4.7, 5.9, 7.1, 8.3, 9.6, 10.8, 12.0
delta = 12 es la mas cercana sin que el soporte del kernel, [z0 - 4, z0 + 4], cruce la placa superior (z = +16,
a mitad del enlace): el borde del soporte queda justo en la placa, donde W = 0.
Puntos medios entre esas posiciones (2026-11-05), multiplos de 0.1 (4.05 -> 4.0 y 8.95 -> 9.0, eleccion del usuario):
  delta = 1.6, 2.8, 4.0, 5.3, 6.5, 7.7, 9.0, 10.2, 11.4
z = 17.50 es el plano medio (mitad de camino entre los nodos 17 y 18) y z = 18.00 es un nodo: las dos son posiciones de
simetria de la red en z; 17.60, 17.75 y 17.90 no lo son (primer momento discreto del Hann no nulo).
delta = 0 se repite con el ejecutable actual (la corrida de C-electrolito-hann8-LD6/z_17.50 usa un enlace a src/ y se
hizo con una version anterior del codigo); esa corrida vieja no se toca.

Solo cambian, respecto de C-electrolito-hann8-LD6/z_17.50: la coordenada z de la particula (config.cds.init) y la
cabecera del input. Ludwig.exe se copia (no enlace) para que recompilar no afecte estas corridas.
Uso: python3 crear_desplazada.py   y despues   cd C-electrolito-hann8-LD6-desplazada && bash correr.sh 3 > bateria.log
     python3 crear_desplazada.py N=10000   las mismas corridas desde cero con N_cycles = 10000 y psi cada 1000 pasos
     (freq_config 1000), para verificar el estacionario: C-electrolito-hann8-LD6-desplazada-N10000.
"""
import os, shutil, sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "C-electrolito-hann8-LD6")
BASE = os.path.join(SRC, "z_17.50")
EXE = os.path.join(HERE, "../../../src/Ludwig.exe")
NCYC = int(([a[2:] for a in sys.argv[1:] if a.startswith("N=")] or ["5000"])[0])
BAT = os.path.join(HERE, "C-electrolito-hann8-LD6-desplazada" + ("-N%d" % NCYC if NCYC != 5000 else ""))
DELTAS = (0.0, 0.10, 0.25, 0.40, 0.50,
          0.60, 0.70, 0.80, 0.90,    # 2026-11-05: completan la celda, cada 0.1 (posiciones de no simetria)
          1.0, 2.2, 3.4, 4.7, 5.9, 7.1, 8.3, 9.6, 10.8, 12.0,   # acercamiento a la placa (ver arriba)
          1.6, 2.8, 4.0, 5.3, 6.5, 7.7, 9.0, 10.2, 11.4,        # puntos medios del acercamiento
          # 2026-11-06 (pedido del usuario): delta = 1 ... 4 cada 0.1 (completa 3 celdas de dx = 1)
          1.1, 1.2, 1.3, 1.4, 1.5, 1.7, 1.8, 1.9, 2.0, 2.1, 2.3, 2.4, 2.5, 2.6, 2.7, 2.9, 3.0, 3.1, 3.2, 3.3, 3.5, 3.6, 3.7, 3.8, 3.9,
          # 2026-11-06 (pedido del usuario): resto de la grilla 0 ... 12 cada 0.1, igual que Hann-4 (crear_desplazada_ref.py r=0.5)
          0.2, 0.3, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.8, 4.9, 5.0, 5.1, 5.2, 5.4, 5.5, 5.6, 5.7, 5.8,
          6.0, 6.1, 6.2, 6.3, 6.4, 6.6, 6.7, 6.8, 6.9, 7.0, 7.2, 7.3, 7.4, 7.5, 7.6, 7.8, 7.9, 8.0, 8.1, 8.2,
          8.4, 8.5, 8.6, 8.7, 8.8, 8.9, 9.1, 9.2, 9.3, 9.4, 9.5, 9.7, 9.8, 9.9, 10.0, 10.1, 10.3, 10.4, 10.5,
          10.6, 10.7, 10.9, 11.0, 11.1, 11.2, 11.3, 11.5, 11.6, 11.7, 11.8, 11.9)

HEADER = """\
##############################################################################
#
#  C-electrolito-hann8-LD6-desplazada (crear_desplazada.py)
#  Igual que C-electrolito-hann8-LD6/z_17.50 salvo la posicion de la carga:
#  z = {z:.2f} = centro del canal (17.50) + {dl:.2f} hacia la placa superior.
#  Hann-8, dx = 1, L = 32, lambda_D = 6, q = +1 fija en x = y = 16, placas a tierra.
#  Fuera de z = 17.50 (plano medio) y z = 18.00 (nodo) la posicion no es de
#  simetria de la red en z.
#
##############################################################################
"""

for dl in DELTAS:
    z = 17.5 + dl
    d = os.path.join(BAT, "z_%05.2f" % z)
    if os.path.exists(d):                             # no tocar corridas ya creadas
        print("%s: ya existe, no se modifica" % d); continue
    for sd in ("", "logs", "colloid_data", "proceced_data"):
        os.makedirs(os.path.join(d, sd), exist_ok=True)
    t = open(os.path.join(BASE, "input")).read()
    t = HEADER.format(z=z, dl=dl) + t[t.index("\nN_start"):]
    if NCYC != 5000:
        assert t.count("\nN_cycles 5000\n") == 1 and t.count("\nfreq_config 5000\n") == 1
        t = t.replace("\nN_cycles 5000\n", "\nN_cycles %d            # verificacion del estacionario (base: 5000)\n" % NCYC)
        t = t.replace("\nfreq_config 5000\n", "\nfreq_config 1000        # psi cada 1000 pasos: deriva psi(N) - psi(N - 1000)\n")
    if (dl > 0.5 or dl in (0.2, 0.3)) and NCYC == 5000:   # 0.2, 0.3 (2026-11-06): tambien, para verificar el estacionario
        assert t.count("\nfreq_config 5000\n") == 1
        t = t.replace("\nfreq_config 5000\n", "\nfreq_config 1000        # psi/efield cada 1000 pasos: estacionario con psi(5000) - psi(4000)\n")
    open(os.path.join(d, "input"), "w").write(t)
    lines = open(os.path.join(BASE, "config.cds.init.001-001")).read().split("\n")
    assert [float(v) for v in lines[35].split()] == [16.0, 16.0, 17.5]
    lines[35] = "    %.15e    %.15e    %.15e" % (16.0, 16.0, z)
    open(os.path.join(d, "config.cds.init.001-001"), "w").write("\n".join(lines))
    shutil.copy(EXE, os.path.join(d, "Ludwig.exe"))
    shutil.copy(os.path.join(BASE, ".petscrc"), os.path.join(d, ".petscrc"))
    print("%s: carga en (16, 16, %.2f)" % (d, z))
shutil.copy(os.path.join(SRC, "correr.sh"), os.path.join(BAT, "correr.sh"))
