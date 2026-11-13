#!/usr/bin/env python3
"""Plot of the chemical-potential residual vs viscosity (2026-11-10, user's request).

Runs: the 7-viscosity sweep of analizar_visc.py (dx = 2, r = 4, 24^3, 3000 steps), last dump (step 3000).
Delta mu_a = max_n - min_n of mu_a/kT = ln rho_a + Z_a psi (same definition as energia_2part.py / analizar_visc.py).
Log-log axes; the points are joined by straight segments only as a guide (no fit).
Usage: python3 grafico_dmu_visc.py   (from this directory) -> ../graficos/dmu_vs_viscosidad.png
"""
import os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from analizar_visc import RUNS, ETA, Z

ST = 3000
dmu = []
for run in RUNS:
    p = os.path.join(HERE, "..", run)
    psi = np.loadtxt(os.path.join(p, "psi-%09d.001-001" % ST))
    q = np.loadtxt(os.path.join(p, "qsi-%09d.001-001" % ST))
    dmu.append([np.ptp(np.log(q[:, a]) + Z[a] * psi) for a in (0, 1)])
dmu = np.array(dmu)
for e, (a, b) in zip(ETA, dmu):
    print("eta %.2f  dmu_+ %.3e  dmu_- %.3e" % (e, a, b))

fig, ax = plt.subplots(figsize=(6, 4.2))
for i, (lab, col, mk) in enumerate(((r"$\Delta\mu_+$", "#2a78d6", "o"), (r"$\Delta\mu_-$", "#eb6834", "s"))):
    ax.plot(ETA, dmu[:, i], "-", marker=mk, color=col, lw=2, ms=7, label=lab)
    ax.annotate(lab, (ETA[-1], dmu[-1, i]), xytext=(8, 0), textcoords="offset points", va="center", color="#333333")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel(r"viscosity $\eta$ (lattice units)")
ax.set_ylabel(r"$\Delta\mu = \max - \min$  ($k_BT$)")
ax.set_title(r"Floor of $\Delta\mu$ vs viscosity ($\Delta x=2$, $r=4$, step 3000)", fontsize=10)
ax.grid(True, which="major", color="#dddddd", lw=0.6); ax.grid(True, which="minor", color="#eeeeee", lw=0.4)
for s in ("top", "right"): ax.spines[s].set_visible(False)
ax.set_xlim(0.015, 3.5); ax.legend(frameon=False)
fig.tight_layout()
out = os.path.join(HERE, "..", "graficos", "dmu_vs_viscosidad.png")
fig.savefig(out, dpi=150); print(out)
