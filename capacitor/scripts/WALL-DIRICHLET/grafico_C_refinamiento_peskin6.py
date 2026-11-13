#!/usr/bin/env python3
"""Nested refinement dx = 2, 1, 1/2 with Ludwig's stretched Peskin 6-point kernel (C^3; subgrid_peskin6_stretch
s = 2r, fixed physical support |u| < 6, sigma = 1.690), runs built by "crear_refinado.py peskin6 [q=0.25] 4 8 16",
against the analytic Peskin-6 solution (grafico_C_hoja_peskin6.py); side by side with the Hann-8 runs
(C^1, support |u| <= 4) of grafico_C_refinamiento_q.py, and with the generalised Hann window cos^4 (subgrid_hann_power 4,
C^3, the SAME support |u| <= 4 and lattices as Hann-8, sigma = 1.132; "crear_refinado.py p=4 [q=0.25] 4 8 16", analytic
solution grafico_C_hoja_hannp.py).

Question: is the loss of order near the support edge (order ~1 instead of 2) caused by the jump of W'' of Hann?
Peskin-6 is C^3 (W'''' is the first derivative that jumps), so if the edge order is a kernel-smoothness effect it
must go back to 2 here.

Per lattice, band RMS of the nodal residual on bands of width 2 with ends on faces of dx = 2 (whole cells in all
three lattices, exact midpoint rule; see grafico_C_refinamiento_anidado.py), and the split of the residual with the two
charges (grafico_C_refinamiento_q.py): d_q/q = L + q^2 N, so at each node
    L = [16 (d/q)_{0.25} - (d/q)_1]/15     (part linear in the source: discretisation error)
    N = 16 [(d/q)_1 - (d/q)_{0.25}]/15      (part ~ q^2: nonlinear Poisson-Boltzmann response at q = 1)
Usage: python3 grafico_C_refinamiento_peskin6.py
"""
import os
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import grafico_C_formulas as F
import grafico_C_hoja_hann as HH
import grafico_C_hoja_peskin6 as HP
import grafico_C_hoja_hannp as H4
from grafico_C_refinamiento import Z0, L, psi_files, load, line
from grafico_C_refinamiento_anidado import band_rms
from grafico_C_refinamiento_q import QRUNS

PRUNS = {q: [(0.5, "C-electrolito-peskin6s1-LD3-ref0.5%s/z_05.50" % sfx), (1.0, "C-electrolito-peskin6s2-LD6-ref1%s/z_09.50" % sfx),
             (2.0, "C-electrolito-peskin6s4-LD12-ref2%s/z_17.50" % sfx)] for q, sfx in ((1.0, ""), (0.25, "-q0.25"))}
H4RUNS = {q: [(0.5, "C-electrolito-hann4p4-LD3-ref0.5%s/z_05.50" % sfx), (1.0, "C-electrolito-hann8p4-LD6-ref1%s/z_09.50" % sfx),
              (2.0, "C-electrolito-hann16p4-LD12-ref2%s/z_17.50" % sfx)] for q, sfx in ((1.0, ""), (0.25, "-q0.25"))}
KERNELS = (("Hann-8 ($C^1$, support $\\pm4$)", HH, QRUNS, 4.0, "#eb6834"),
           ("Peskin-6 $s=2r$ ($C^3$, support $\\pm6$)", HP, PRUNS, 6.0, "#2a78d6"),
           ("Hann-8 $\\cos^4$ ($C^3$, support $\\pm4$)", H4, H4RUNS, 4.0, "#1baf7a"))
W = 2.0

def residuals(model, runs, q, mu):
    """{r: (z, (psi - q phi_lin)/q)} for the finished runs."""
    out = {}
    for r, path in runs[q]:
        fs = psi_files(path)
        if not fs: print("  (no data yet) %s" % path); continue
        a = load(fs[-1], r); z, ps, _ = line(a, r, 0)
        if len(fs) > 2:
            print("  %-45s drift over the last N/4: %.1e" % (path, np.abs(a - load(fs[-2], r)).max()))
        out[r] = (z, (ps - q * model.phi(0.0, z, 0, 0, mu)) / q)
    return out

def bands(S):
    """{r: band RMS} and the band centres."""
    R = {r: band_rms(z, d, r, W, W) for r, (z, d) in S.items()}
    return (next(iter(R.values()))[0] if R else None), {r: v for r, (_, v) in R.items()}

def orders(R):
    return {(r0, r1): np.log2(R[r0] / R[r1]) for r0, r1 in ((0.5, 1.0), (1.0, 2.0)) if r0 in R and r1 in R}

if __name__ == "__main__":
    F.style()
    fig, ax = plt.subplots(3, 1, figsize=(10.0, 11.0), sharex=True, gridspec_kw=dict(hspace=0.08))
    for tag, model, runs, sup, col in KERNELS:
        mu = model.mean_U(); print(tag)
        S1 = residuals(model, runs, 1.0, mu); S4 = residuals(model, runs, 0.25, mu)
        zc, R1 = bands(S1)
        if zc is None: continue
        both = sorted(set(S1) & set(S4))
        Ls = {r: (S1[r][0], (16 * S4[r][1] - S1[r][1]) / 15) for r in both}
        Ns = {r: (S1[r][0], 16 * (S1[r][1] - S4[r][1]) / 15) for r in both}
        _, RL = bands(Ls); _, RN = bands(Ns)
        print("  band centre        " + " ".join("%7.1f" % c for c in zc))
        for name, R in (("d (q=1)", R1), ("L      ", RL)):
            for (r0, r1), p in orders(R).items():
                print("  order %-7s %g->%-3s" % (name, 1 / r0, "1/2" if r1 == 2 else "%g" % (1 / r1)) + " ".join("%7.2f" % v for v in p))
        for r in sorted(RN):
            print("  RMS N   dx=%-4g     " % (1 / r) + " ".join("%7.1e" % v for v in RN[r]))
            print("  RMS L   dx=%-4g     " % (1 / r) + " ".join("%7.1e" % v for v in RL[r]))
        for k, (R, a) in enumerate(((R1, ax[0]), (RL, ax[1]))):
            for (r0, r1), p in orders(R).items():
                fine = r1 == 2.0
                a.plot(zc, p, "o-" if fine else "s--", ms=5 if fine else 3.5, lw=1.8 if fine else 0.9, color=col,
                       alpha=1.0 if fine else 0.5, label="%s, $\\Delta x = %g \\to %s$" % (tag, 1 / r0, "1/2" if fine else "%g" % (1 / r1)))
        for r in sorted(RN):
            ax[2].semilogy(zc, RN[r], "o-", ms=3 + 1.2 * np.log2(2 * r), lw=1.2, color=col, alpha=0.45 + 0.27 * np.log2(2 * r),
                           label="%s, N, $\\Delta x = %s$" % (tag, "1/2" if r == 2 else "%g" % (1 / r)))
        if 2.0 in RL:
            ax[2].semilogy(zc, RL[2.0], ":", lw=1.6, color=col, label="%s, L, $\\Delta x = 1/2$" % tag)
        for a in ax:
            a.axvline(Z0 - sup, color=col, lw=0.9, ls="--", alpha=0.6); a.axvline(Z0 + sup, color=col, lw=0.9, ls="--", alpha=0.6)
    for a in ax:
        a.axvline(Z0, color="#b9b7ae", lw=0.8, ls=":")
        for x in np.arange(-L / 2, L / 2 + 1e-9, 2.0):
            a.axvline(x, color="#e4e2dc", lw=0.5, zorder=0)
    for a in ax[:2]:
        a.axhline(2, color="#55544f", lw=0.9, ls=":"); a.axhline(1, color="#b9b7ae", lw=0.9, ls="-."); a.set_ylim(-1, 4)
    ax[0].set_ylabel("local order of RMS$(d)$, $q=1$\n(residual against linear DH)")
    ax[1].set_ylabel("local order of RMS$(L)$\n(part linear in the source)")
    ax[2].set_ylabel("band RMS of $N$ (nonlinear part, $q=1$)\nand of $L$ at $\\Delta x=1/2$")
    ax[2].set_xlabel("band centre z (physical units, from the mid-plane)"); ax[0].set_xlim(-16, 16)
    ax[0].legend(fontsize=8, loc="lower right", ncol=2); ax[1].legend(fontsize=8, loc="lower right", ncol=2)
    ax[2].legend(fontsize=7.5, loc="upper right", ncol=2)
    ax[0].set_title("Nested lattices $\\Delta x = 2, 1, 1/2$, bands of width 2 (dashed verticals: kernel support edges)", fontsize=10.5)
    fig.savefig(os.path.join(F.ROOT, "graficos", "C_refinamiento_peskin6.png"), dpi=170, bbox_inches="tight")
