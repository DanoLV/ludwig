#!/usr/bin/env python3
"""
Carga desplazada (crear_desplazada.py) en la red refinada dx = 1/2: mismo sistema fisico con Hann-16 (ancho fisico 8),
con la misma receta de escalado que crear_refinado.py para r = 2:
  - caja 64 x 64 x 66 (32r x 32r x (32r + 2)), carga en x = y = 16r = 32
  - electrokinetics_epsilon 1e4/r, electrokinetics_init_rho_el 1.388889e-3/r^3 (lambda_D = 6r = 12 nodos)
  - D por paso igual -> N = 5000 r^2 = 20000
Posiciones (eleccion del usuario, 2026-11-05): extremos de la serie de acercamiento a la placa y dos intermedias lo mas
equiespaciadas posible entre las ya corridas con dx = 1 (ideales 4.67 y 8.33):
  delta = 1.0, 4.7, 8.3, 12.0   (z0 = delta fisico desde el plano medio; Ludwig z = 1.5 + r (16 + delta))
Para verificar el estacionario psi/efield se guardan cada 4000 pasos (freq_config 4000): mismo tiempo fisico que 1000
pasos con dx = 1 (el paso de tiempo fisico escala como 1/r^2); criterio sobre psi(20000) - psi(16000).
Se parte del input de C-electrolito-hann8-LD6/z_17.50 (el mismo de las corridas desplazadas con dx = 1).
Uso: python3 crear_desplazada_ref.py   y despues   cd C-electrolito-hann16-LD12-ref2-desplazada && bash correr.sh 2 > bateria.log
     python3 crear_desplazada_ref.py r=0.5   red dx = 2 (2026-11-06): Hann-4, caja 16 x 16 x 18, lambda_D = 3 nodos,
     N = 2000 (5000 r^2 = 1250 redondeado a 1000 hacia arriba, como crear_refinado.py), psi cada 250 pasos (1000 r^2):
     estacionario con psi(2000) - psi(1750). Carga en Ludwig z = 10.00, 11.85, 13.65, 15.50.
     -> C-electrolito-hann4-LD3-ref0.5-desplazada. Con r=0.5 corre ademas todas las delta de crear_desplazada.py (28 en total).
"""
import math, os, re, shutil, sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "C-electrolito-hann8-LD6")
BASE = os.path.join(SRC, "z_17.50")
EXE = os.path.join(HERE, "../../../src/Ludwig.exe")
R = float(([a[2:] for a in sys.argv[1:] if a.startswith("r=")] or ["2"])[0])
N = int(round(8 * R)); L = int(round(32 * R))
EPS0, RHO0, SIG8 = 1.0e4, 1.38888889e-03, 1.44612
NSTEPS = int(math.ceil(5000 * R * R / 1000.0) * 1000); FREQ = int(1000 * R * R)
BAT = os.path.join(HERE, "C-electrolito-hann%d-LD%g-ref%g-desplazada" % (N, 6 * R, R))
DELTAS = (1.0, 4.7, 8.3, 12.0)
if R == 2:                                     # 2026-11-07 (pedido del usuario): celda delta = 1 ... 3 cada 0.1 (solo Hann-16)
    DELTAS += (1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0)
if R == 0.5:                                   # 2026-11-06 (pedido del usuario): todas las delta corridas con Hann-8 (crear_desplazada.py)
    DELTAS += (0.0, 0.10, 0.25, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90,
               2.2, 3.4, 5.9, 7.1, 9.6, 10.8, 1.6, 2.8, 4.0, 5.3, 6.5, 7.7, 9.0, 10.2, 11.4,
               # 2026-11-06 (pedido del usuario): delta = 1 ... 4 cada 0.1, igual que con Hann-8
               1.1, 1.2, 1.3, 1.4, 1.5, 1.7, 1.8, 1.9, 2.0, 2.1, 2.3, 2.4, 2.5, 2.6, 2.7, 2.9, 3.0, 3.1, 3.2, 3.3, 3.5, 3.6, 3.7, 3.8, 3.9,
               # 2026-11-06 (pedido del usuario): resto de la grilla 0 ... 12 cada 0.1 (solo Hann-4)
               0.2, 0.3, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.8, 4.9, 5.0, 5.1, 5.2, 5.4, 5.5, 5.6, 5.7, 5.8,
               6.0, 6.1, 6.2, 6.3, 6.4, 6.6, 6.7, 6.8, 6.9, 7.0, 7.2, 7.3, 7.4, 7.5, 7.6, 7.8, 7.9, 8.0, 8.1, 8.2,
               8.4, 8.5, 8.6, 8.7, 8.8, 8.9, 9.1, 9.2, 9.3, 9.4, 9.5, 9.7, 9.8, 9.9, 10.0, 10.1, 10.3, 10.4, 10.5,
               10.6, 10.7, 10.9, 11.0, 11.1, 11.2, 11.3, 11.5, 11.6, 11.7, 11.8, 11.9)

HEADER = """\
##############################################################################
#
#  {bat} (crear_desplazada_ref.py)
#  Carga desplazada con dx = {dx:g}: Hann-{n} (ancho fisico 8), L = {L} nodos (= 32 fisicos),
#  lambda_D = {ld:g} nodos (= 6). Carga q = +1 fija en x = y = {xc:g} nodos, z = {zc:.3f} nodos
#  = delta {dl:.2f} fisico desde el plano medio hacia la placa superior. Placas a tierra.
#
##############################################################################
"""

def sub(txt, pat, new):
    out, k = re.subn(pat, new, txt, flags=re.M)
    assert k == 1, pat
    return out

for dl in DELTAS:
    zc = 1.5 + R * (16 + dl)
    d = os.path.join(BAT, "z_%05.2f" % zc)
    if os.path.exists(d):
        print("%s: ya existe, no se modifica" % d); continue
    for sd in ("", "logs", "colloid_data", "proceced_data"):
        os.makedirs(os.path.join(d, sd), exist_ok=True)
    t = open(os.path.join(BASE, "input")).read()
    t = HEADER.format(bat=os.path.basename(BAT), dx=1 / R, n=N, L=L, ld=6 * R, xc=16 * R, zc=zc, dl=dl) + t[t.index("\nN_start"):]
    t = sub(t, r"^N_cycles .*$", "N_cycles %d             # 5000 r^2 redondeado a 1000 hacia arriba (D por paso igual que dx = 1)" % NSTEPS)
    t = sub(t, r"^size .*$", "size %d_%d_%d          # 32r x 32r x (32r + 2 nodos de pared)" % (L, L, L + 2))
    t = sub(t, r"^electrokinetics_epsilon .*$",
            "electrokinetics_epsilon %.8e   # 1e4/r: q/(eps kT) = 10 r en unidades de red, l_B = %.6f" % (EPS0 / R, 0.7957747 * R))
    t = sub(t, r"^electrokinetics_init_rho_el .*$",
            "electrokinetics_init_rho_el %.8e   # 1.388889e-3/r^3 -> lambda_D = 6 r = %g" % (RHO0 / R ** 3, 6 * R))
    t = sub(t, r"^pm_sr_sigma .*$", "pm_sr_sigma       %.5f      # = sigma_kernel (no se usa: pm_sr_correction no)" % (SIG8 * R))
    t = sub(t, r"^pm_sr_rcut .*$", "pm_sr_rcut        %.1f      # 2 + krange (no se usa)" % (2 + N // 2))
    t = sub(t, r"^subgrid_hann_order .*$", "subgrid_hann_order %d      # orden Hann: ancho fisico n/r = 8" % N)
    t = sub(t, r"^vel_io_freq .*$", "vel_io_freq %d" % NSTEPS)
    t = sub(t, r"^freq_config .*$", "freq_config %d        # psi/efield cada 1000 r^2 pasos: estacionario con psi(N) - psi(N - %d)" % (FREQ, FREQ))
    t = sub(t, r"^# lambda_D = 6 -> .*$",
            "# lambda_D = 6r = %g -> rho_s = 1/(8 pi l_B lambda_D^2) = %.6e por especie (l_B = %.6f)"
            % (6 * R, RHO0 / R ** 3, 0.7957747 * R))
    t = sub(t, r"^# Solo importa el ESTACIONARIO.*$",
            "# Solo importa el ESTACIONARIO (no depende de D). D por paso igual que dx = 1: en unidades fisicas el paso\n"
            "# de tiempo es 1/r^2 del de dx = 1, por eso N = 5000 r^2.")
    open(os.path.join(d, "input"), "w").write(t)
    lines = open(os.path.join(BASE, "config.cds.init.001-001")).read().split("\n")
    assert [float(v) for v in lines[35].split()] == [16.0, 16.0, 17.5]
    lines[35] = "".join("{:>25.15e}".format(v) for v in (16 * R, 16 * R, zc))
    open(os.path.join(d, "config.cds.init.001-001"), "w").write("\n".join(lines))
    shutil.copy(EXE, os.path.join(d, "Ludwig.exe")); shutil.copy(os.path.join(BASE, ".petscrc"), os.path.join(d, ".petscrc"))
    print("%s: size %d_%d_%d, N = %d, carga en (%g, %g, %.2f)" % (d, L, L, L + 2, NSTEPS, 16 * R, 16 * R, zc))
shutil.copy(os.path.join(SRC, "correr.sh"), os.path.join(BAT, "correr.sh"))
