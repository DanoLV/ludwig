#!/usr/bin/env python3
"""Local convergence order from the RMS of the nodal residual on NARROW overlapping bands (rho = 0, z axis through
the charge), for the Hann-n and Gaussian refinement series (grafico_C_refinamiento.py, grafico_C_refinamiento_gauss.py).

Bands [s, s + W] with W = 2 physical, s every STEP = 1 (consecutive bands overlap by 1).
Each node stands for its cell [z_k - dx/2, z_k + dx/2] (cell-centred lattice, wall at half a link). The band integral
is the sum over cells of d_k^2 times the length of the cell that lies inside the band (exact for whole cells, i.e. the
midpoint rule, and clipped at the two ends):
    RMS_[s, s+W] = ( (1/W) sum_k d_k^2 |cell_k ∩ [s, s+W]| )^(1/2)
Unlike the 4-unit bands, the ends of a 2-unit band are not cell faces of every lattice (common faces of dx = 1, 0.8,
2/3, 1/2 are 4 apart), so the clipped end cells add a quadrature error O(dx^2) that depends on where the end falls in
the cell. Its size is measured with a smooth test function (QUAD check below) and printed: up to 10 % for dx = 0.8 and
2/3, enough to fake orders of 3 or -4 where the residual changes by 10x inside a band (support edges).
dx = 1 and dx = 1/2 have cell faces at every integer, so for that pair every band end is a face of both lattices and
the quadrature is exact (midpoint rule): the order 1 -> 1/2 (thick black) is free of this artefact.
Usage: python3 grafico_C_refinamiento_bandas.py
"""
import os
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import grafico_C_formulas as F
import grafico_C_hoja_hann as HH
import grafico_C_hoja_gauss as HG
from grafico_C_refinamiento import RUNS, Z0, L
from grafico_C_refinamiento_gauss import GRUNS, RCOL, series

W, STEP = 2.0, 1.0
S = np.arange(-L / 2, L / 2 - W + 1e-9, STEP)                  # band starts
ZC = S + W / 2                                                  # band centres

def band_rms(z, d, r):
    h = 1 / r; lo, hi = z - h / 2, z + h / 2
    out = []
    for s in S:
        w = np.clip(np.minimum(hi, s + W) - np.maximum(lo, s), 0, None)
        out.append(np.sqrt(np.sum(w * d ** 2) / W))
    return np.array(out)

def quad_check(runs):
    """Relative error of the clipped-cell quadrature on f = psi_analytic'' (smooth, scales like the residual)."""
    mu = HH.mean_U(); zf = np.linspace(-L / 2, L / 2, 4001)
    f = lambda z: (HH.phi(0.0, z + 1e-2, 0, 0, mu) - 2 * HH.phi(0.0, z, 0, 0, mu) + HH.phi(0.0, z - 1e-2, 0, 0, mu)) / 1e-4
    ff = f(zf) ** 2
    exact = np.array([np.sqrt(np.trapezoid(ff[(zf >= s) & (zf <= s + W)], zf[(zf >= s) & (zf <= s + W)]) / W) for s in S])
    print("quadrature check (f = psi''): max relative error of the band RMS over all bands")
    for r, _ in runs:
        z = (np.arange(1, int(round(32 * r)) + 1) - 0.5) / r - L / 2
        print("  dx = %.3f   %.2e" % (1 / r, np.abs(band_rms(z, f(z), r) / exact - 1).max()))

def orders(Sx):
    R = np.array([band_rms(z, d, r) for r, z, d in Sx]); dx = np.array([1 / r for r, _, _ in Sx])
    return dx, R, np.log(R[:-1] / R[1:]) / np.log(dx[:-1] / dx[1:])[:, None]

if __name__ == "__main__":
    F.style()
    quad_check([(r, p) for _, r, p in RUNS])
    SH = series([(r, p) for _, r, p in RUNS], HH, HH.mean_U()); SG = series(GRUNS, HG, HG.mean_U())
    fig, ax = plt.subplots(2, 2, figsize=(13.0, 7.6), sharex=True, gridspec_kw=dict(height_ratios=[1, 1], hspace=0.08, wspace=0.22))
    for j, (tag, Sx, ksup) in enumerate((("Hann-$n$ (support $|u|\\leq 4$)", SH, 4.0), ("Gaussian $\\sigma=1.446$ (cut $|u|=7.95$)", SG, 7.95))):
        dx, R, P = orders(Sx)
        for k, (r, _, _) in enumerate(Sx):
            ax[0, j].semilogy(ZC, R[k], "o-", ms=3.5, lw=1.1, color=RCOL[r], label="$\\Delta x=%s$" % ("1" if r == 1 else "1/%g" % r))
        for k in range(len(dx) - 1):
            ax[1, j].plot(ZC, P[k], "o-", ms=3.0, lw=0.9, color=RCOL[Sx[k + 1][0]], alpha=0.45,
                          label="$%s \\to %s$ (clipped cells)" % tuple("1" if x == 1 else "1/%g" % (1 / x) for x in dx[k:k + 2]))
        i1, i2 = list(dx).index(1.0), list(dx).index(0.5)
        P12 = np.log(R[i1] / R[i2]) / np.log(2.0)
        ax[1, j].plot(ZC, P12, "o-", ms=4.5, lw=2.0, color="#2b2a27", label="$1 \\to 1/2$ (exact quadrature)", zorder=5)
        print("  order 1 -> 1/2 (exact)  " + " ".join("%6.2f" % v for v in P12))
        print(tag); print("  band centre " + " ".join("%6.1f" % c for c in ZC))
        for k in range(len(dx) - 1):
            print("  order %.3f->%.3f " % (dx[k], dx[k + 1]) + " ".join("%6.2f" % v for v in P[k]))
        for a in ax[:, j]:
            a.axvspan(Z0 - 4, Z0 + 4, color="#f1efe8", zorder=0); a.axvline(Z0, color="#b9b7ae", lw=0.8, ls=":")
            if ksup > 4: a.axvline(Z0 - ksup, color="#b9b7ae", lw=0.8, ls="--"); a.axvline(Z0 + ksup, color="#b9b7ae", lw=0.8, ls="--")
        ax[1, j].axhline(2, color="#55544f", lw=0.9, ls=":"); ax[1, j].axhline(1, color="#b9b7ae", lw=0.9, ls="-.")
        ax[1, j].set_ylim(-1, 4); ax[1, j].set_xlabel("band centre z (physical units, from the mid-plane)")
        ax[0, j].set_title(tag, fontsize=11); ax[0, j].legend(fontsize=8.5, loc="upper right")
        ax[1, j].legend(fontsize=8.5, loc="lower right", title="local order", title_fontsize=8.5)
    ax[0, 0].set_ylabel("RMS of the nodal residual\non 2-unit bands, $\\rho=0$"); ax[1, 0].set_ylabel("local order (log ratio)")
    ax[0, 0].set_xlim(-16, 16)
    fig.suptitle("Bands of width 2, one every 1 unit (overlap 1); shaded: Hann support, dashed: Gaussian cut", fontsize=11, y=0.95)
    fig.savefig(os.path.join(F.ROOT, "graficos", "C_refinamiento_bandas_finas.png"), dpi=170, bbox_inches="tight")
