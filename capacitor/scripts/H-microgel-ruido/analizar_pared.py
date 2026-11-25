#!/usr/bin/env python3
"""
Análisis de la prueba de pared con ruido (H-microgel-ruido/pared-*, 2026-11-13).

Para cada corrida (un monómero libre que empieza a d = 4.0 del plano z = 1.5):
  * si terminó o abortó (mensaje de Ludwig);
  * h = d - 3.5 (d = distancia del centro al plano sólido/fluido más cercano; 3.5 = radio del
    kernel Hann-8: h < 0 significa que el kernel toca nodos de pared): mínimo, fracción h < 0;
  * desplazamiento por paso: rms por componente y máximo frente a drmax = 0.8;
  * V(h)/kT de la Steele (Hoang Ngoc Minh et al. 2023, ecs. 7-8, cortada en z* = 0.98698 sigma_w):
    máximo alcanzado y, dentro del alcance de la pared (0 < h < z*), comparación del promedio
    <V>/kT de la corrida con el de Boltzmann, int V e^{-V/kT} dh / int e^{-V/kT} dh en (0, z*).
    Esa comparación supone equilibrio local junto a la pared (densidad uniforme sin pared);
    con 2000 pasos el muestreo es corto: es un indicador, no un test estadístico.
    Para q = 1 la imagen agrega atracción: no se compara con Boltzmann de la pared sola.
Uso:  python3 capacitor/scripts/H-microgel-ruido/analizar_pared.py [carpetas...]
"""
import glob
import math
import os
import re
import sys

import numpy as np
from scipy.integrate import quad

H = "/home/danolv/Sim/ludwig/capacitor/local/WALL-DIRICHLET/H-microgel-ruido"
KT, ZLO, ZHI, RK, SIG, DRMAX = 1e-4, 1.5, 33.5, 3.5, 1.0, 0.8


def vads(h, e):
    return 2 * math.pi * e * ((2 / 5) * (SIG / h) ** 10 - (SIG / h) ** 4
                              - math.sqrt(2) / (3 * (h / SIG + 0.61 / math.sqrt(2)) ** 3))


def V(h, e):
    zs = 0.98698 * SIG
    if h <= 0:
        return float("inf")
    return vads(h, e) - vads(zs, e) if h < zs else 0.0


def analyse(d):
    name = os.path.basename(d)
    eps = float(re.search(r"^wall_ss_cut_epsilon\s+(\S+)", open(f"{d}/input").read(), re.M).group(1))
    q = re.search(r"pared-q(\d)", name).group(1)
    files = sorted(glob.glob(f"{d}/config.cds0*.001-001"))
    log = open(f"{d}/logs/output.txt", errors="replace").read()
    fin = "terminó normalmente" if "finished normally" in log else "ABORTÓ / no terminó"
    err = re.findall(r"dr= *[0-9.e+-]+.{0,40}exceeded max *[0-9.e+-]+", log)[:1]
    z = np.array([float(open(f).read().split("\n")[35].split()[2]) for f in files])
    xyz = np.array([[float(v) for v in open(f).read().split("\n")[35].split()[:3]] for f in files])
    if len(files) < 2:
        print(f"\n== {name}: {len(files)} paso(s), {fin}" + (f"; {err}" if err else ""))
        return
    dr = np.diff(xyz, axis=0)
    dist = np.minimum(z - ZLO, ZHI - z)
    h = dist - RK
    Vs = np.array([V(x, eps) for x in h]) / KT
    zs = 0.98698 * SIG
    print(f"\n== {name} (q = {q}, eps_w = {eps:g}): {len(files)} pasos, {fin}" + (f"; {err}" if err else ""))
    print(f"   h: inicial {h[0]:.3f}, mínimo {h.min():.4f}, fracción con h < 0: {(h < 0).mean():.2e}; "
          f"fracción dentro del alcance (h < {zs:.3f}): {(h < zs).mean():.2f}")
    print(f"   paso: rms por componente {np.sqrt((dr**2).mean()):.4f}, máximo |dr| {np.linalg.norm(dr,axis=1).max():.4f} "
          f"({np.linalg.norm(dr,axis=1).max()/DRMAX*100:.0f} % de drmax)")
    fin_V = Vs[np.isfinite(Vs)]
    print(f"   V/kT máximo alcanzado {fin_V.max() if fin_V.size else float('nan'):.2f}")
    inside = (h > 0) & (h < zs)
    if q == "0" and inside.sum() > 20:
        num = quad(lambda x: V(x, eps) / KT * math.exp(-V(x, eps) / KT), 1e-4, zs, limit=200)[0]
        den = quad(lambda x: math.exp(-V(x, eps) / KT), 1e-4, zs, limit=200)[0]
        print(f"   dentro del alcance: <V>/kT corrida {Vs[inside].mean():.3f} vs Boltzmann {num/den:.3f} "
              f"({inside.sum()} muestras; indicador, muestreo corto)")


def main():
    for d in (sys.argv[1:] or sorted(glob.glob(f"{H}/pared-*/"))):
        analyse(d.rstrip("/"))


if __name__ == "__main__":
    main()
