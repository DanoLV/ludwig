#!/usr/bin/env python3
"""Two-panel plot vs viscosity: chemical-potential residual (top) and fluid speed (bottom) (2026-11-10, user's request).

Runs: the 7-viscosity sweep of analizar_visc.py (dx = 2, r = 4, 24^3), last dump (step 3000).
Top: Delta mu_a = max_n - min_n of mu_a/kT = ln rho_a + Z_a psi (as grafico_dmu_visc.py).
Bottom: max_n |u_n| and RMS_n |u_n| (vel dump, lattice units), as analizar_visc.py.
Shared log eta axis; the points are joined by straight segments only as a guide (no fit).
Usage: python3 grafico_dmu_u_visc.py   (from this directory) -> ../graficos/dmu_u_vs_viscosidad.png
"""
import os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from analizar_visc import RUNS, ETA, Z

ST = 3000
dmu, uu = [], []
for run in RUNS:
    p = os.path.join(HERE, "..", run)
    psi = np.loadtxt(os.path.join(p, "psi-%09d.001-001" % ST))
    q = np.loadtxt(os.path.join(p, "qsi-%09d.001-001" % ST))
    dmu.append([np.ptp(np.log(q[:, a]) + Z[a] * psi) for a in (0, 1)])
    u = np.sqrt((np.loadtxt(os.path.join(p, "vel-%09d.001-001" % ST)) ** 2).sum(1))
    uu.append([u.max(), np.sqrt(np.mean(u ** 2))])
dmu, uu = np.array(dmu), np.array(uu)
for e, (a, b), (c, d) in zip(ETA, dmu, uu):
    print("eta %.2f  dmu_+ %.3e  dmu_- %.3e  max|u| %.3e  RMS|u| %.3e" % (e, a, b, c, d))

BLUE, ORANGE = "#2a78d6", "#eb6834"
fig, (a1, a2) = plt.subplots(2, 1, figsize=(6, 7), sharex=True)
for ax, data, labs, ylab in ((a1, dmu, (r"$\Delta\mu_+$", r"$\Delta\mu_-$"), r"$\Delta\mu = \max - \min$  ($k_BT$)"),
                             (a2, uu, (r"max$|u|$", r"RMS$|u|$"), r"$|u|$ (lattice units)")):
    for i, (col, mk) in enumerate(((BLUE, "o"), (ORANGE, "s"))):
        ax.plot(ETA, data[:, i], "-", marker=mk, color=col, lw=2, ms=7, label=labs[i])
        ax.annotate(labs[i], (ETA[-1], data[-1, i]), xytext=(8, 0), textcoords="offset points", va="center", color="#333333")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_ylabel(ylab)
    ax.grid(True, which="major", color="#dddddd", lw=0.6); ax.grid(True, which="minor", color="#eeeeee", lw=0.4)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    ax.legend(frameon=False)
a2.set_xlim(0.015, 3.5); a2.set_xlabel(r"viscosity $\eta$ (lattice units)")
a1.set_title(r"$\Delta\mu$ and $|u|$ vs viscosity ($\Delta x=2$, $r=4$, step 3000)", fontsize=10)
fig.tight_layout()
out = os.path.join(HERE, "..", "graficos", "dmu_u_vs_viscosidad.png")
fig.savefig(out, dpi=150); print(out)
