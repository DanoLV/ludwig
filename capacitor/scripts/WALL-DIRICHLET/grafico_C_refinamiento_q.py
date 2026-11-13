#!/usr/bin/env python3
"""Linearity test of the dx-independent floor C(z) in the refinement residual (Hann, nested lattices dx = 2, 1, 1/2):
same runs with charge q = 1 (grafico_C_refinamiento_anidado.py) and q = 0.25 ("crear_refinado.py q=0.25 4 8 16").

The analytic solution is linear Debye-Hueckel, so it scales with q: residual d_q = psi_Ludwig(q) - q phi_lin.
Ludwig solves Nernst-Planck, whose equilibrium is nonlinear Poisson-Boltzmann (ion densities ~ exp(-+psi)):
psi = q phi_lin + q^3 phi_3 + ..., so the nonlinear part of d_q / q scales as q^2. If the floor is nonlinearity,
    C_{0.25}/0.25 = C_1/1 / 16.
Discretisation errors are linear in the source: (d_q/q) at a fixed dx is the same for both q except for the floor.

Band means (signed): M_[a,b] = (1/(b-a)) sum_{z_k in [a,b]} d_k dx, bands of width 2 with ends on faces of dx = 2
(nested: whole cells in all three lattices, exact midpoint rule, no interpolation). M is linear in d, so the model
M = C + a dx^2 holds band by band. C from Richardson with dx = 1 and 1/2: C = (4 M_{1/2} - M_1)/3; the 3-lattice
least-squares fit (with dx = 2) is printed as a check.
Usage: python3 grafico_C_refinamiento_q.py
"""
import os
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import grafico_C_formulas as F
import grafico_C_hoja_hann as HH
from grafico_C_refinamiento import Z0, L, psi_files, load, line
from grafico_C_refinamiento_anidado import NRUNS, band_rms

QRUNS = {1.0: NRUNS,
         0.25: [(0.5, "C-electrolito-hann4-LD3-ref0.5-q0.25/z_05.50"), (1.0, "C-electrolito-hann8-LD6-ref1-q0.25/z_09.50"),
                (2.0, "C-electrolito-hann16-LD12-ref2-q0.25/z_17.50")]}
QCOL = {1.0: "#2b2a27", 0.25: "#2a78d6"}
W = 2.0

def band_mean(z, d, r):
    S = np.arange(-L / 2, L / 2 - W + 1e-9, W); out = []
    for s in S:
        m = (z > s) & (z < s + W)
        assert m.sum() == int(round(W * r))
        out.append(d[m].mean())
    return S + W / 2, np.array(out)

def residuals(q, mu):
    """[(r, z, (psi - q phi_lin)/q)] for the finished runs of one q."""
    out = []
    for r, path in QRUNS[q]:
        fs = psi_files(path)
        if not fs: print("  (no data yet) %s" % path); continue
        a = load(fs[-1], r); z, ps, _ = line(a, r, 0)
        if len(fs) > 2:                                   # needs an intermediate dump (not the initial one)
            print("  q = %g  dx = %g  drift over the last N/4: %.1e" % (q, 1 / r, np.abs(a - load(fs[-2], r)).max()))
        out.append((r, z, (ps - q * HH.phi(0.0, z, 0, 0, mu)) / q))
    return out

if __name__ == "__main__":
    F.style(); mu = HH.mean_U()
    fig, ax = plt.subplots(2, 1, figsize=(10.0, 8.2), sharex=True, gridspec_kw=dict(hspace=0.08))
    C = {}
    for q in (1.0, 0.25):
        Sx = residuals(q, mu)
        if len(Sx) < 3: continue
        dx = np.array([1 / r for r, _, _ in Sx])
        M = np.array([band_mean(z, d, r)[1] for r, z, d in Sx]); zc = band_mean(Sx[0][1], Sx[0][2], Sx[0][0])[0]
        R = np.array([band_rms(z, d, r, W, W)[1] for r, z, d in Sx])
        i1, i2 = list(dx).index(1.0), list(dx).index(0.5)
        C[q] = (4 * M[i2] - M[i1]) / 3
        A = np.vstack([np.ones_like(dx), dx ** 2]).T; co = np.linalg.lstsq(A, M, rcond=None)[0]
        p12 = np.log(R[i1] / R[i2]) / np.log(2.0)
        print("q = %g   band centre       " % q + " ".join("%9.1f" % c for c in zc))
        print("   order RMS 1 -> 1/2      " + " ".join("%9.2f" % v for v in p12))
        print("   C/q (Richardson 1, 1/2) " + " ".join("%+.2e" % v for v in C[q]))
        print("   C/q (fit dx = 2, 1, 1/2)" + " ".join("%+.2e" % v for v in co[0]))
        print("   max |M - fit| / |M|_max " + " ".join("%9.1e" % v for v in np.abs(M - A @ co).max(0) / np.abs(M).max(0)))
        ax[0].plot(zc, p12, "o-", ms=5, lw=1.6, color=QCOL[q], label="$q=%g$" % q)
        ax[1].plot(zc, C[q], "o-", ms=5, lw=1.6, color=QCOL[q], label="$C/q$, $q=%g$" % q)
    if 0.25 in C:
        ax[1].plot(zc, 16 * C[0.25], "s--", ms=4, lw=1.0, color=QCOL[0.25], alpha=0.6, label="$16\\times C/q$, $q=0.25$")
        print("   ratio (C/q)_1 / (C/q)_0.25 " + " ".join("%9.1f" % v for v in C[1.0] / C[0.25]))
    for a in ax:
        a.axvspan(Z0 - 4, Z0 + 4, color="#f1efe8", zorder=0); a.axvline(Z0, color="#b9b7ae", lw=0.8, ls=":")
    ax[0].axhline(2, color="#55544f", lw=0.9, ls=":"); ax[0].axhline(1, color="#b9b7ae", lw=0.9, ls="-.")
    ax[0].set_ylim(-1, 4); ax[0].set_ylabel("local order of RMS$(d/q)$, $\\Delta x = 1 \\to 1/2$")
    ax[1].axhline(0, color="#8a8880", lw=0.8); ax[1].set_ylabel("floor $C/q$ (band mean, Richardson)")
    ax[1].set_xlabel("band centre z (physical units, from the mid-plane)"); ax[0].set_xlim(-16, 16)
    ax[0].legend(fontsize=9, loc="lower right"); ax[1].legend(fontsize=9, loc="lower right")
    ax[0].set_title("Hann, nested lattices $\\Delta x=2,1,1/2$, bands of width 2: linearity test of the $\\Delta x$-independent floor",
                    fontsize=10.5)
    fig.savefig(os.path.join(F.ROOT, "graficos", "C_refinamiento_q_hann.png"), dpi=170, bbox_inches="tight")
