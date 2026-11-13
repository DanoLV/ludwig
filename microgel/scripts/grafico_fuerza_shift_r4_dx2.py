#!/usr/bin/env python3
"""Force between two fixed charges at FIXED separation r = 4 vs the shift delta of the pair inside the cell, dx = 2
(2026-11-11, user's request; runs of crear_2part_dx2_r4_shift.py: SLIDES-LD8-hann04-corr_no-dx2-r4-shift/d_X.X,
both charges moved along x by delta (physical) = 0 ... 2.0, one full cell; exe 45dc94bc; 3000 steps).
Force on each charge, physical units: F = kT Esub_x / dx (dx = 2), last step (proceced_data/particle_Esub.csv;
index 1 = charge at larger x, index 2 = the other). Plotted: F_x on charge 1 and -F_x on charge 2 (equal for an
isolated pair), with the reference F_ref(4) of error_fuerza_2part.py (periodic DH, continuous Hann of width 8; it does not
depend on delta). Panels: force, relative error (F - F_ref)/F_ref in %, and (2026-11-11) variation about the mean of one period
(F - <F>)/<F> in %, <F> over delta = 0 ... 1.8.
Printed: table, and the control delta = 0 vs delta = 2.0 (same configuration translated by one node).
Usage: python3 grafico_fuerza_shift_r4_dx2.py   -> ../FUERZA-2PART-MEAS/FUERZA_shift_r4_dx2.png
"""
import os
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from error_fuerza_2part import F_ref, KT, HERE, OUT

DX = 2.0; R = 4.0
B = os.path.join(HERE, "..", "SLIDES-LD8-hann04-corr_no-dx2-r4-shift")
DELTAS = [round(0.2 * i, 1) for i in range(11)]

def esub_last(path):
    rows = [l.strip().replace(",", ".").split(";") for l in open(os.path.join(path, "proceced_data", "particle_Esub.csv"))
            if not l.startswith("#")]
    st = max(int(v[0]) for v in rows)
    return {int(v[1]): float(v[3]) for v in rows if int(v[0]) == st}, st

if __name__ == "__main__":
    fr = F_ref(np.array([R]))[0]; F1, F2 = [], []
    print("delta |   F1 (charge 1)  |  -F2 (charge 2)  | e1 (%)   e2 (%)  | step")
    for d in DELTAS:
        e, st = esub_last(os.path.join(B, "d_%.1f" % d)); assert st == 3000, (d, st)
        F1.append(KT * e[1] / DX); F2.append(-KT * e[2] / DX)
        print("%5.1f | %.9e | %.9e | %+7.3f %+7.3f | %d" % (d, F1[-1], F2[-1], 100 * (F1[-1] - fr) / fr, 100 * (F2[-1] - fr) / fr, st))
    F1, F2 = np.array(F1), np.array(F2)
    print("F_ref(4) = %.9e; control delta 0 vs 2.0: rel. diff F1 %.1e, F2 %.1e"
          % (fr, abs(F1[-1] - F1[0]) / fr, abs(F2[-1] - F2[0]) / fr))
    # 2026-11-11 (user's request): third panel, variation about the mean force of the interval,
    # (F - <F>)/<F> in %, <F> = plain mean over ONE period, delta = 0 ... 1.8 (10 points; delta = 2.0 is the same
    # configuration as delta = 0 and is left out of the mean so it is not counted twice; choice of this script).
    m1, m2 = F1[:-1].mean(), F2[:-1].mean()
    print("<F> over delta = 0 ... 1.8: charge 1 %.9e (%+.3f %% vs F_ref), charge 2 %.9e; max |F - <F>|/<F>: %.3f %%, %.3f %%"
          % (m1, 100 * (m1 - fr) / fr, m2, 100 * np.abs(F1 - m1).max() / m1, 100 * np.abs(F2 - m2).max() / m2))
    fig, ax = plt.subplots(3, 1, figsize=(7.5, 10.5), sharex=True, gridspec_kw=dict(hspace=0.08))
    ax[0].axhline(fr, color="#8a8880", lw=1.0, label="reference $F_{ref}(4)$: periodic DH, continuous Hann (width 8)")
    # 2026-11-11 (user's request): only the force on charge 1 is drawn (charge 2 is still read and printed: -F2 = F1 to
    # <= 5e-6), points joined by lines.
    kw = dict(color="#eb6834", ms=7, lw=1.0)
    ax[0].plot(DELTAS, F1, "^-", label="Ludwig, $\\Delta x = 2$ (Hann-4): $F_x$ on charge 1", **kw)
    ax[1].plot(DELTAS, 100 * (F1 - fr) / fr, "^-", label="charge 1", **kw)
    ax[1].axhline(0, color="#8a8880", lw=0.8)
    ax[2].plot(DELTAS, 100 * (F1 - m1) / m1, "^-", label="charge 1", **kw)
    ax[2].axhline(0, color="#8a8880", lw=0.8)
    ax[2].set_ylabel("$(F - \\langle F \\rangle) / \\langle F \\rangle$ (%)\n$\\langle F \\rangle$: mean over $\\delta = 0 \\ldots 1.8$")
    ax[0].set_ylabel("force (lattice units)"); ax[1].set_ylabel("$(F - F_{ref}) / F_{ref}$ (%)")
    ax[2].set_xlabel("shift $\\delta$ of the pair along $x$ (physical units; one cell = 2)")
    ax[2].set_xticks(DELTAS)
    for a in ax:
        a.grid(True, axis="x", color="#e5e5e5", lw=0.6)
    ax[0].legend(fontsize=8)   # one series: no legend needed in the error panels
    ax[0].set_title("Two charges $q = +1$ at fixed separation $r = 4$, both shifted by $\\delta$\n$\\Delta x = 2$, "
                    "box 48, $\\lambda_D = 8$, $k_BT = 10^{-5}$ (lattice units)", fontsize=10, loc="left")
    out = os.path.join(OUT, "FUERZA_shift_r4_dx2.png"); fig.savefig(out, dpi=300, bbox_inches="tight"); print(out)
