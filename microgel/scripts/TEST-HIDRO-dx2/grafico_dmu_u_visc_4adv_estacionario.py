#!/usr/bin/env python3
"""Delta mu and |u| vs viscosity for four NP advection variants, STEADY STATE (step 6000) (2026-11-16, user's request).

Same as grafico_dmu_u_visc_4adv.py (step 3000) but read from the continuations cont6000/<variant>/<run>
(crear_cont_visc.py; steady by estacionario_visc.py: relative change of Delta mu over steps 5000-6000
<= 1.7e-6 in all 28 runs). Variants as listed below (the cont6000 copies use the same exe as the source).

  old     ../hidro_si*              exe 4f0c20e7, centred advection WITHOUT link weights (ions 9x too fast)
  centred ../advfix/hidro_si*       exe 221f3229, centred advection with W_p = 2 wgradients[p] (crear_advfix.py)
  flla    ../advfix_flla/hidro_si*  exe 221f3229, electrokinetics_advection flla (crear_advfix_flla.py)
  flla_mass ../advfix_flla_mass/hidro_si*  exe 221f3229, electrokinetics_advection flla_mass
          (crear_advfix_flla_mass.py)
Same quantities and dump as grafico_dmu_u_visc.py: Delta mu_a = max_n - min_n of mu_a/kT = ln rho_a + Z_a psi
(top), max_n|u_n| and RMS_n|u_n| (bottom), step 3000. Straight segments only guide the eye (no fit).
Prints the table of values and of the ratios centred/variant.
Usage: python3 scripts/grafico_dmu_u_visc_4adv_estacionario.py  -> graficos/dmu_u_vs_viscosidad_4adv_estacionario.png
"""
import os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "/home/danolv/Sim/ludwig/microgel/local/TEST-HIDRO-dx2"   # scripts/ is a symlink: absolute paths
sys.path.insert(0, os.path.join(BASE, "scripts"))
from analizar_visc import RUNS, ETA, Z

ST = 6000
SETS = (("old exe (unweighted centred)", os.path.join(BASE, "cont6000", "old"), "--", "white"),
        ("centred (weighted)", os.path.join(BASE, "cont6000", "centred"), "-", None),
        ("flla", os.path.join(BASE, "cont6000", "flla"), ":", None),
        ("flla_mass", os.path.join(BASE, "cont6000", "flla_mass"), "-.", "white"))


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


data = [collect(root) for _, root, _, _ in SETS]
names = ("dmu_+", "dmu_-", "max|u|", "RMS|u|")
for j, nm in enumerate(names):
    blk, c = divmod(j, 2)
    print("\n%s:  eta | %s | centred/old  centred/flla  centred/flla_mass" % (nm, " | ".join("%-12s" % s[0][:12] for s in SETS)))
    for i, e in enumerate(ETA):
        v = [d[blk][i, c] for d in data]
        print("      %5.2f | %s | %8.3f  %8.3f  %8.3f" % (e, " | ".join("%.4e  " % x for x in v),
                                                         v[1] / v[0], v[1] / v[2], v[1] / v[3]))

BLUE, ORANGE = "#2a78d6", "#eb6834"
MK = ("o", "s")
fig, (a1, a2) = plt.subplots(2, 1, figsize=(7.0, 8.6), sharex=True)
panels = ((a1, 0, (r"$\Delta\mu_+$", r"$\Delta\mu_-$"), r"$\Delta\mu = \max - \min$  ($k_BT$)"),
          (a2, 1, (r"max$|u|$", r"RMS$|u|$"), r"$|u|$ (lattice units)"))
for ax, blk, labs, ylab in panels:
    for i, col in enumerate((BLUE, ORANGE)):
        for (lab, _, ls, mfc), d in zip(SETS, data):
            ax.plot(ETA, d[blk][:, i], ls, marker=MK[i], color=col, lw=1.8, ms=6.5,
                    mfc=(mfc or col), label="%s, %s" % (labs[i], lab))
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_ylabel(ylab)
    ax.grid(True, which="major", color="#dddddd", lw=0.6); ax.grid(True, which="minor", color="#eeeeee", lw=0.4)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    ax.legend(frameon=False, fontsize=7, ncol=2)
uu = np.array([d[1] for d in data])
if np.allclose(uu, uu[0], rtol=1e-3, atol=0):
    a2.text(0.03, 0.04, "the four variants coincide (rel. diff. < 1e-3):\nthe advection scheme does not change the flow",
            transform=a2.transAxes, fontsize=8, color="#555555")
a2.set_xlim(0.015, 3.5); a2.set_xlabel(r"viscosity $\eta$ (lattice units)")
a1.set_title("$\\Delta\\mu$ and $|u|$ vs viscosity for four NP advection variants\n"
             r"($\Delta x=2$, $r=4$, steady state: step 6000)" "\n" "dashed = old exe, solid = centred, dotted = flla, dash-dot = flla_mass", fontsize=9.5)
fig.tight_layout()
out = os.path.join(BASE, "graficos", "dmu_u_vs_viscosidad_4adv_estacionario.png")
if os.path.exists(out):
    sys.exit(out + " already exists (not overwritten)")
fig.savefig(out, dpi=150); print(out)
