#!/usr/bin/env python3
"""
Carga desplazada con Hann cos^4 de sigma (casi) igual al Hann-8 cos^2 (2026-11-09, pedido del usuario): prueba en
delta = 1 con dx = 1 y dx = 2, para la convergencia y para comparar con las corridas cos^2 de la misma posicion.

Ancho: sigma^2 = n^2 (1/12 - 1/(2 pi^2)) para cos^2 y n^2 (1/12 - 5/(8 pi^2)) para cos^4 (comentario de d_hann, subgrid.c).
  cos^2, n = 8:  sigma = 1.44607 (fisico)
  cos^4: n = 10.223 daria el mismo sigma, pero n debe ser entero (suma 1 exacta). n = 10 (eleccion del usuario):
         sigma = 1.41449 fisico (2.2 % menor); n = 11 daria 1.55594 (7.6 % mayor).
Escalado con el refinamiento igual que cos^2 (ancho fisico fijo 10): dx = 1 -> n = 10, dx = 2 -> n = 5 (impar: permitido,
la suma 1 sigue exacta para n entero >= p/2 + 1 = 3; solo cambia entre 5 y 4 el numero de nodos activos con el offset).

Cada corrida es una copia de la corrida cos^2 de la misma delta y la misma red; solo cambian subgrid_hann_order,
subgrid_hann_power (nuevo, 4), pm_sr_sigma (no se usa: pm_sr_correction no) y la cabecera. Mismo Ludwig.exe (se copia
el de la corrida base, md5 4f0c20e7..., = src/Ludwig.exe del 2026-11-04).
  dx = 1: base C-electrolito-hann8-LD6-desplazada/z_18.50          -> C-electrolito-hann10p4-LD6-desplazada/z_18.50
  dx = 2: base C-electrolito-hann4-LD3-ref0.5-desplazada/z_10.00   -> C-electrolito-hann5p4-LD3-ref0.5-desplazada/z_10.00
Las corridas que ya existen no se tocan.
Uso: python3 crear_desplazada_p4.py   y despues, en cada bateria,   bash correr.sh 1 > bateria.log
"""
import os, re, shutil

HERE = os.path.dirname(os.path.abspath(__file__))
SIG4 = 0.1414490                     # sigma cos^4 / n = sqrt(1/12 - 5/(8 pi^2))
RUNS = (  # (base, bateria nueva, z, n, dx)
    ("C-electrolito-hann8-LD6-desplazada", "C-electrolito-hann10p4-LD6-desplazada", "z_18.50", 10, 1.0),
    ("C-electrolito-hann4-LD3-ref0.5-desplazada", "C-electrolito-hann5p4-LD3-ref0.5-desplazada", "z_10.00", 5, 2.0))

HEADER = """\
##############################################################################
#
#  {bat}/{z} (crear_desplazada_p4.py, 2026-11-09)
#  Copia de {base}/{z} (delta = 1, dx = {dx:g}) con kernel Hann cos^4:
#  subgrid_hann_order {n} (ancho fisico 10), subgrid_hann_power 4,
#  sigma = {sl:.5f} nodos = {sf:.5f} fisico (Hann-8 cos^2: 1.44607 fisico).
#  Todo lo demas igual a la corrida base.
#
##############################################################################
"""

def sub(t, pat, new):
    t2, k = re.subn(pat, new, t, flags=re.M)
    assert k == 1, pat
    return t2

for base, bat, z, n, dx in RUNS:
    src = os.path.join(HERE, base, z); dst = os.path.join(HERE, bat, z)
    if os.path.exists(dst):
        print("%s: ya existe, no se modifica" % dst); continue
    sl = SIG4 * n
    t = open(os.path.join(src, "input")).read()
    t = HEADER.format(bat=bat, base=base, z=z, dx=dx, n=n, sl=sl, sf=sl * dx) + t[t.index("\nN_start"):]
    t = sub(t, r"^subgrid_hann_order .*$",
            "subgrid_hann_order %d      # ancho fisico n dx = 10 (cos^4 con sigma ~ Hann-8 cos^2)\n"
            "subgrid_hann_power 4       # w ~ cos^4(pi u/n): C^3, suma 1 exacta (n >= 3)" % n)
    t = sub(t, r"^pm_sr_sigma .*$", "pm_sr_sigma       %.5f      # = sigma_kernel cos^4 (no se usa: pm_sr_correction no)" % sl)
    for sd in ("", "logs", "colloid_data", "proceced_data"):
        os.makedirs(os.path.join(dst, sd), exist_ok=True)
    open(os.path.join(dst, "input"), "w").write(t)
    for f in ("config.cds.init.001-001", "Ludwig.exe", ".petscrc"):
        shutil.copy2(os.path.join(src, f), os.path.join(dst, f))
    shutil.copy(os.path.join(HERE, "C-electrolito-hann8-LD6", "correr.sh"), os.path.join(HERE, bat, "correr.sh"))
    print("%s: n = %d, sigma = %.5f fisico" % (dst, n, sl * dx))
