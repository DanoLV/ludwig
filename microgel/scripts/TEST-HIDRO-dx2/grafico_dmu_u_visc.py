#!/usr/bin/env python3
"""Two-panel plot vs viscosity: chemical-potential residual (top) and fluid speed (bottom) (2026-11-10, user's request).

Runs: the 7-viscosity sweep of analizar_visc.py (dx = 2, r = 4, 24^3), last dump (step 3000).
Top: Delta mu_a = max_n - min_n of mu_a/kT = ln rho_a + Z_a psi (as grafico_dmu_visc.py).
Bottom: max_n |u_n| and RMS_n |u_n| (vel dump, lattice units), as analizar_visc.py.
Shared log eta axis; the points are joined by straight segments only as a guide (no fit).
Usage: python3 grafico_dmu_u_visc.py   (from this directory) -> ../graficos/dmu_u_vs_viscosidad.png

2026-11-16 (user's request): added to the top panel
  * Delta mu with the weighted centred advection (advfix/, exe 221f3229, W_p = 2 wgradients[p]), same step 3000;
    open markers + dashed = old exe (unweighted advection, the original curves), filled + solid = weighted;
  * Delta mu WITHOUT hydrodynamics (hidro_no_exe45dc94bc, exe 45dc94bc, hydrodynamics no, eta irrelevant) at step
    3000 as horizontal lines; its value at step 9000 (hidro_no_exe45dc94bc-cont2) is written as text, because at
    step 3000 it is still decreasing and 4e-14 would compress the axis by 6 decades.
  |u| (bottom) is unchanged: identical in the old and weighted runs (rel. diff < 1e-3, grafico_dmu_u_visc_3adv.py).
  The previous version of the figure is kept as ../graficos/dmu_u_vs_viscosidad_v1_20261110.png.

2026-11-16, later (user's request: points in the steady state): the hydrodynamic curves are now taken at step 6000
from the continuations cont6000/old/<run> (old exe) and cont6000/centred/<run> (weighted), steady by
estacionario_visc.py (|dmu(6000) - dmu(5000)|/dmu(6000) < 1e-4); |u| from cont6000/centred at step 6000 (the
relative difference max/RMS|u| old vs weighted is printed). The no-hydrodynamics value is the one at step 9000
(hidro_no_exe45dc94bc-cont2), constant from step 6250 to 9000 (dumps every 250) at ~4.4e-14, the round-off level.
The step-3000 version of this figure is kept as ../graficos/dmu_u_vs_viscosidad_v2_paso3000.png.
"""
import os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = "/home/danolv/Sim/ludwig/microgel/local/TEST-HIDRO-dx2/scripts"   # scripts/ is a symlink: absolute path
sys.path.insert(0, HERE)
from analizar_visc import RUNS, ETA, Z

ST = 3000
BASE = "/home/danolv/Sim/ludwig/microgel/local/TEST-HIDRO-dx2"   # not HERE/..: the OS resolves the symlink


def dmu_of(p, st):
    psi = np.loadtxt(os.path.join(p, "psi-%09d.001-001" % st))
    q = np.loadtxt(os.path.join(p, "qsi-%09d.001-001" % st))
    return [np.ptp(np.log(q[:, a]) + Z[a] * psi) for a in (0, 1)]


ST = 6000                       # steady state (cont6000), see the docstring
dmu, dmuw, uu, uo = [], [], [], []
for run in RUNS:
    po, pw = os.path.join(BASE, "cont6000", "old", run), os.path.join(BASE, "cont6000", "centred", run)
    dmu.append(dmu_of(po, ST))
    dmuw.append(dmu_of(pw, ST))
    for p, lst in ((pw, uu), (po, uo)):
        u = np.sqrt((np.loadtxt(os.path.join(p, "vel-%09d.001-001" % ST)) ** 2).sum(1))
        lst.append([u.max(), np.sqrt(np.mean(u ** 2))])
dmu, dmuw, uu, uo = np.array(dmu), np.array(dmuw), np.array(uu), np.array(uo)
print("max rel. diff of (max|u|, RMS|u|) old vs weighted at step %d: %.2e" % (ST, np.abs(uo / uu - 1).max()))
dno = dmu_of(os.path.join(BASE, "hidro_no_exe45dc94bc-cont2"), 9000)
dno6 = dmu_of(os.path.join(BASE, "hidro_no_exe45dc94bc-cont2"), 6250)
for e, (a, b), (aw, bw), (c, d) in zip(ETA, dmu, dmuw, uu):
    print("eta %.2f  dmu_+ %.3e  dmu_- %.3e | weighted dmu_+ %.3e  dmu_- %.3e | max|u| %.3e  RMS|u| %.3e"
          % (e, a, b, aw, bw, c, d))
print("hydrodynamics no: step 6250 dmu_+ %.3e dmu_- %.3e; step 9000 dmu_+ %.3e dmu_- %.3e" % (*dno6, *dno))

BLUE, ORANGE = "#2a78d6", "#eb6834"
fig, (a1, a2) = plt.subplots(2, 1, figsize=(9.8, 8), sharex=True)
# Simulation parameters (from hidro_si/input and run.log; identical in all runs of the figure except viscosity,
# executable and hydrodynamics on/off, checked with diff of the inputs and cmp of config.cds.init, 2026-11-16).
PARAMS = "\n".join([
    r"$\bf{Simulation\ parameters}$",
    r"(lattice units; 1 lattice spacing = 2 physical units)",
    "",
    r"box $24^3$ (48$^3$ physical), periodic",
    r"2 fixed charges $q=+1$, Hann-4 kernel",
    r"  (width 8 physical), at $x = 4.25,\ 6.25$",
    r"  $y=z=12.25$: separation 2 (4 physical)",
    r"1:1 salt, $n_\pm = 6.175\times10^{-3}$",
    r"  per node ($7.72\times10^{-4}$ physical)",
    r"$\ell_B = 0.398$ (0.796 physical)",
    r"$\lambda_D = 4$ (8 physical), $\kappa\ell_B = 0.10$",
    r"$k_BT = 10^{-5}$, $\varepsilon = 2\times10^{4}$",
    r"ion $D$: input 0.01, physical",
    r"  $0.01\,A_0 = 0.0614$ ($A_0$ = 6.138)",
    r"$\rho_{LB} = 0.8$, $\eta = \eta_{bulk}$ (x axis)",
    r"LB M10, no thermal noise",
    r"Poisson: FFT; NP stencil 27, centred",
    r"Hydrodynamics off: same, $\eta$ irrelevant",
])
fig.subplots_adjust(right=0.66)
fig.text(0.685, 0.5, PARAMS, fontsize=8.5, va="center", ha="left", color="#333333", linespacing=1.45)
for i, (col, mk) in enumerate(((BLUE, "o"), (ORANGE, "s"))):
    sp = (r"$\Delta\mu_+$", r"$\Delta\mu_-$")[i]
    a1.plot(ETA, dmu[:, i], "--", marker=mk, color=col, lw=1.6, ms=7, mfc="white", label=sp + ", old exe (unweighted)")
    a1.plot(ETA, dmuw[:, i], "-", marker=mk, color=col, lw=2, ms=7, label=sp + r", weighted $W_p$")
    a1.axhline(dno[i], color=col, lw=1.4, ls=":", label=sp + ", hydrodynamics off (step 9000)")
a1.set_ylim(1e-14, 3e-5)
a1.text(0.02, 0.10, "hydrodynamics off: constant from step 6250 to 9000\n(%.1e, %.1e $k_BT$; probably round-off, not verified)"
        % (dno[0], dno[1]), transform=a1.transAxes, fontsize=7.5, color="#555555")
for i, (col, mk) in enumerate(((BLUE, "o"), (ORANGE, "s"))):
    lab = (r"max$|u|$", r"RMS$|u|$")[i]
    a2.plot(ETA, uu[:, i], "-", marker=mk, color=col, lw=2, ms=7, label=lab)
    a2.annotate(lab, (ETA[-1], uu[-1, i]), xytext=(8, 0), textcoords="offset points", va="center", color="#333333")
a2.text(0.02, 0.03, "max / RMS over the lattice nodes of $|u|$ (step 6000, weighted exe);\nold exe: same within %.0e"
        % np.abs(uo / uu - 1).max(),
        transform=a2.transAxes, fontsize=7.5, color="#555555")
for ax, ylab in ((a1, r"$\Delta\mu = \max - \min$  ($k_BT$)"), (a2, r"$|u|$ (lattice units)")):
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_ylabel(ylab)
    ax.grid(True, which="major", color="#dddddd", lw=0.6); ax.grid(True, which="minor", color="#eeeeee", lw=0.4)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    ax.legend(frameon=False, fontsize=7.5 if ax is a1 else 9, ncol=2 if ax is a1 else 1)
a2.set_xlim(0.015, 3.5); a2.set_xlabel(r"viscosity $\eta$ (lattice units)")
a1.set_title(r"$\Delta\mu$ and $|u|$ vs viscosity ($\Delta x=2$, charge separation $r=4$ physical), steady state:" "\n"
             r"step 6000 with hydrodynamics, step 9000 without", fontsize=10)
fig.tight_layout(rect=(0, 0, 0.67, 1))
out = os.path.join(BASE, "graficos", "dmu_u_vs_viscosidad.png")
fig.savefig(out, dpi=150); print(out)
