#!/usr/bin/env python3
"""
Figura comparativa de las baterias C: error relativo del campo de imagen en
funcion de la distancia a la pared, con/sin fondo de contraiones, con/sin
correccion pm_sr, Hann-8 y Hann-6. Escribe graficos/C_comparacion.png.
"""
import os, sys, glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.argv = sys.argv[:1]
import analizar_C as ac

ROOT = ac.ROOT
casos = [  # carpeta, orden, etiqueta, color, marcador, relleno
    ("C-carga-imagen-hann08",   8, "Hann-8, background, no correction", "#2a78d6", "o", True),
    ("C-corr-hann08",           8, "Hann-8, background, correction",    "#2a78d6", "x", False),
    ("C-sin-fondo-hann08",      8, "Hann-8, net charge, no correction", "#eb6834", "s", True),
    ("C-sin-fondo-corr-hann08", 8, "Hann-8, net charge, correction",    "#eb6834", "+", False),
    ("C-carga-imagen-hann06",   6, "Hann-6, background, no correction", "#1baf7a", "^", True),
]
plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.color": "#e4e2dc",
                     "axes.spines.top": False, "axes.spines.right": False,
                     "legend.frameon": False, "savefig.dpi": 160, "savefig.bbox": "tight"})
fig, ax = plt.subplots(figsize=(8.5, 5))
for carpeta, n, lab, c, mk, fill in casos:
    dmin = n / 2.0 - 0.5
    D, R = [], []
    for d in sorted(glob.glob(os.path.join(ROOT, carpeta, "z_*")),
                    key=lambda p: float(p.rsplit("_", 1)[1])):
        z0 = float(d.rsplit("_", 1)[1]); dd = z0 - 1.5
        if dd < dmin or dd > 15.9: continue
        r = ac.read_run(d)
        th = ac.image_field(dd, k0=not r["fondo"])
        D.append(dd); R.append(abs(r["ez"] - th) / abs(th))
    kw = dict(mfc=c if fill else "none", mec=c, color=c)
    if mk in "x+": kw = dict(color=c)
    ax.semilogy(D, R, mk, ms=6 if mk not in "x+" else 8, mew=1.4, label=lab, **kw)
for n, c in ((8, "#2a78d6"), (6, "#1baf7a")):
    ax.axvline(n / 2.0 - 0.5, color=c, lw=0.9, ls=":")
ax.text(3.55, 3e-2, "Hann-8\nd_min", color="#2a78d6", fontsize=8.5)
ax.text(2.55, 3e-2, "Hann-6\nd_min", color="#1baf7a", fontsize=8.5)
ax.set(xlabel="distance to the wall plane d (lattice units)",
       ylabel=r"$|E_z - E_{theory}| / |E_{theory}|$", ylim=(1e-6, 0.1),
       title="Image field on one charged subgrid particle between grounded walls:\n"
             "relative error vs theory (kernel clear of the wall)")
ax.legend(loc="upper right", fontsize=9)
out = os.path.join(ROOT, "graficos", "C_comparacion.png")
fig.savefig(out); print("  ->", os.path.relpath(out, ROOT))
