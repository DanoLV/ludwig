#!/usr/bin/env python3
"""Delta mu and |u| vs viscosity, BEFORE vs AFTER the NP advection-weight fix (2026-11-16, user's request).

Before: the 7-viscosity sweep of 2026-11-10 (../hidro_si*, exe 4f0c20e7: unweighted centred advection, ions advected
9x too fast with stencil 27). After: the same 7 inputs rerun with the current exe (../advfix/hidro_si*, md5 221f3229,
CHANGE 20261112 W_p = 2 wgradients[p]); see crear_advfix.py. Same quantities and dump as grafico_dmu_u_visc.py:
Delta mu_a = max_n - min_n of mu_a/kT = ln rho_a + Z_a psi (top), max_n|u_n| and RMS_n|u_n| (bottom), step 3000.
Straight segments only guide the eye (no fit). Prints a table with the before/after ratios.
Usage: python3 scripts/grafico_dmu_u_visc_advfix.py  -> graficos/dmu_u_vs_viscosidad_advfix.png
"""
import os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "/home/danolv/Sim/ludwig/microgel/local/TEST-HIDRO-dx2"   # scripts/ is a symlink: absolute paths
sys.path.insert(0, os.path.join(BASE, "scripts"))
from analizar_visc import RUNS, ETA, Z

ST = 3000


def collect(root):
    dmu, uu = [], []
    for run in RUNS:
        p = os.path.join(root, run)
        psi = np.loadtxt(os.path.join(p, "psi-%09d.001-001" % ST))
        q = np.loadtxt(os.path.join(p, "qsi-%09d.001-001" % ST))
        dmu.append([np.ptp(np.log(q[:, a]) + Z[a] * psi) for a in (0, 1)])
        u = np.sqrt((np.loadtxt(os.path.join(p, "vel-%09d.001-001" % ST)) ** 2).sum(1))
        uu.append([u.max(), np.sqrt(np.mean(u ** 2))])
    return np.array(dmu), np.array(uu)


old = collect(BASE)
new = collect(os.path.join(BASE, "advfix"))
print(" eta  | dmu_+ old  new    ratio | dmu_- old  new    ratio | max|u| old new  ratio | RMS|u| old new  ratio")
for i, e in enumerate(ETA):
    cells = []
    for blk, j in ((0, 0), (0, 1), (1, 0), (1, 1)):
        a, b = old[blk][i, j], new[blk][i, j]
        cells.append("%.3e %.3e %6.2f" % (a, b, a / b))
    print("%5.2f | " % e + " | ".join(cells))

BLUE, ORANGE = "#2a78d6", "#eb6834"
fig, (a1, a2) = plt.subplots(2, 1, figsize=(6.4, 7.6), sharex=True)
panels = ((a1, 0, (r"$\Delta\mu_+$", r"$\Delta\mu_-$"), r"$\Delta\mu = \max - \min$  ($k_BT$)"),
          (a2, 1, (r"max$|u|$", r"RMS$|u|$"), r"$|u|$ (lattice units)"))
for ax, blk, labs, ylab in panels:
    for i, (col, mk) in enumerate(((BLUE, "o"), (ORANGE, "s"))):
        ax.plot(ETA, old[blk][:, i], "--", marker=mk, color=col, lw=1.5, ms=7, mfc="white",
                label=labs[i] + " before fix")
        ax.plot(ETA, new[blk][:, i], "-", marker=mk, color=col, lw=2, ms=7, label=labs[i] + " after fix")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_ylabel(ylab)
    ax.grid(True, which="major", color="#dddddd", lw=0.6); ax.grid(True, which="minor", color="#eeeeee", lw=0.4)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    ax.legend(frameon=False, fontsize=8, ncol=2)
a2.text(0.03, 0.04, "before and after coincide (same to 4 digits):\nthe fix does not change the flow",
        transform=a2.transAxes, fontsize=8, color="#555555")
a2.set_xlim(0.015, 3.5); a2.set_xlabel(r"viscosity $\eta$ (lattice units)")
a1.set_title("$\\Delta\\mu$ and $|u|$ vs viscosity, before/after the NP advection-weight fix\n"
             r"($\Delta x=2$, $r=4$, step 3000; open = old exe, filled = current exe)", fontsize=9.5)
fig.tight_layout()
out = os.path.join(BASE, "graficos", "dmu_u_vs_viscosidad_advfix.png")
if os.path.exists(out):
    sys.exit(out + " already exists (not overwritten)")
fig.savefig(out, dpi=150); print(out)
