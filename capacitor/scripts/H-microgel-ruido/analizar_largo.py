#!/usr/bin/env python3
"""
Análisis de las corridas largas H-microgel-ruido/{neutro,superficie}-largo (2026-11-13).

Posiciones cada 100 pasos (no hay pasos consecutivos: el desplazamiento por paso de
analizar_ruido.py no se puede medir aquí). Para cada corrida imprime:
  * serie por bloques de 5000 pasos (heurístico: 50 volcados por bloque) de
      Rg (radio de giro, sin pesos), <U_enlaces>/kT (ref. equipartición armónica 35/2 = 17.5),
      z del centro de masa, d mínimo a los planos z = 1.5 / 33.5 (h = d - 3.5);
    sirve para decidir a ojo desde cuándo la corrida está estacionaria: NO elige la ventana.
  * fluido: <u_a^2>/(kT/rho) en cada vel-* (equipartición LB, Dünweg & Ladd 2009).
  * máximo desplazamiento de un monómero entre volcados (100 pasos), solo como control de saltos.
Reutiliza read_cds y las constantes de analizar_ruido.py. Solo lee archivos.

Uso:  python3 capacitor/scripts/H-microgel-ruido/analizar_largo.py [carpetas...]
"""
import glob
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from analizar_ruido import H, L, KT, RHO, K, R0, ZLO, ZHI, RK, read_cds  # noqa: E402

BLOQUE = 5000


def analyse(d):
    files = sorted(glob.glob(f"{d}/config.cds0*.001-001"))
    steps = np.array([int(os.path.basename(f)[10:18]) for f in files])
    idx, _, bonds = read_cds(files[0])
    bset = sorted({tuple(sorted(b)) for b in bonds})
    i_of = {v: k for k, v in enumerate(idx)}
    pos = np.array([read_cds(f)[1] for f in files])
    # enlaces
    a = np.array([i_of[p] for p, _ in bset]); b = np.array([i_of[q] for _, q in bset])
    v = pos[:, a] - pos[:, b]
    v[..., :2] -= L * np.round(v[..., :2] / L)
    U = 0.5 * K * ((np.linalg.norm(v, axis=-1) - R0) ** 2).sum(-1) / KT
    # Rg, centro de masa, paredes (imagen mínima en x,y respecto del monómero 0 para Rg)
    rel = pos - pos[:, :1]
    rel[..., :2] -= L * np.round(rel[..., :2] / L)
    rg = np.sqrt(((rel - rel.mean(1, keepdims=True)) ** 2).sum(-1).mean(-1))
    zcm = pos[..., 2].mean(1)
    dmin = np.minimum(pos[..., 2] - ZLO, ZHI - pos[..., 2]).min(1)
    dr = np.diff(pos, axis=0)
    dr[..., :2] -= L * np.round(dr[..., :2] / L)
    jump = np.linalg.norm(dr, axis=-1).max()
    log = open(f"{d}/logs/output.txt", errors="replace").read()
    fin = "terminó normalmente" if "finished normally" in log else "NO terminó"
    print(f"\n== {os.path.basename(d)}: pasos {steps[0]}..{steps[-1]} ({len(files)} volcados), {fin}")
    print(f"   {len(bset)} enlaces; máximo |dr| entre volcados (100 pasos) = {jump:.3f}")
    print("   bloque [pasos)      Rg     <U>/kT   z_cm    d_min (h)")
    print(f"   paso {steps[0]:<13d} {rg[0]:6.3f} {U[0]:8.2f} {zcm[0]:7.2f}  {dmin[0]:6.2f}  (inicial)")
    for t0 in range(0, steps[-1], BLOQUE):
        m = (steps > t0) & (steps <= t0 + BLOQUE)
        if m.any():
            print(f"   ({t0:5d},{t0+BLOQUE:5d}]   {rg[m].mean():6.3f} {U[m].mean():8.2f} {zcm[m].mean():7.2f}  "
                  f"{dmin[m].min():6.2f} ({dmin[m].min()-RK:5.2f})")
    for vf in sorted(glob.glob(f"{d}/vel-0*.001-001")):
        u = np.loadtxt(vf)
        if u.ndim != 2 or u.shape[1] < 3:
            continue
        u = u[:, -3:].reshape(32, 32, 34, 3)[:, :, 1:33, :]
        print(f"   fluido {os.path.basename(vf)[:13]}: <u_a^2>/(kT/rho) = "
              + ", ".join(f"{(u[..., c]**2).mean()/(KT/RHO):.3f}" for c in range(3)))


def main():
    for d in (sys.argv[1:] or sorted(glob.glob(f"{H}/*-largo/"))):
        analyse(d.rstrip("/"))


if __name__ == "__main__":
    main()
