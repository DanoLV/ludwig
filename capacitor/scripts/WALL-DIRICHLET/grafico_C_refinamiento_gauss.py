#!/usr/bin/env python3
"""Smoothness control for the order ~1 near the edge of the Hann support (grafico_C_refinamiento.py): the same lattice
refinement with Ludwig's truncated Gaussian kernel (subgrid_kernel gauss; runs built by "crear_refinado.py gauss"),
sigma = 1.44612 physical (same variance as Hann-8), cut at 5.5 sigma, vs the analytic Gaussian solution
(grafico_C_hoja_gauss.py, full Gaussian; the cut changes the kernel by exp(-5.5^2/2) = 2.7e-7 relative).
Hann is only C^1 (W'' jumps at |u| = 4); the Gaussian is smooth. If the order-1 bands of Hann become order 2 with the
Gaussian, the loss of order comes from the non-smooth edge of the Hann kernel.

Metric: RMS of the nodal residual on the same bands as grafico_C_refinamiento.py (ends on cell faces common to every
lattice, so no interpolation), rho = 0, z axis through the charge.
Usage: python3 grafico_C_refinamiento_gauss.py ;  LAYOUT=left python3 grafico_C_refinamiento_gauss.py
"""
import os
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import grafico_C_formulas as F
import grafico_C_hoja_hann as HH
import grafico_C_hoja_gauss as HG
from grafico_C_refinamiento import RUNS, BANDS, Z0, psi_files, load, line, band_rms

GRUNS = [(1.0, "C-electrolito-gauss-LD6-ref1/z_09.50"), (1.25, "C-electrolito-gauss-LD7.5-ref1.25/z_11.50"),
         (1.5, "C-electrolito-gauss-LD9-ref1.5/z_13.50"), (2.0, "C-electrolito-gauss-LD12-ref2/z_17.50")]
RCOL = {1.0: "#eb6834", 1.25: "#2a78d6", 1.5: "#1baf7a", 2.0: "#8e44ad"}

def series(runs, model, mu):
    """[(r, z, residual)] for the finished runs of one kernel."""
    out = []
    for r, path in runs:
        fs = psi_files(path)
        if not fs: print("  (no data yet) %s" % path); continue
        z, ps, _ = line(load(fs[-1], r), r, 0)
        out.append((r, z, ps - model.phi(0.0, z, 0, 0, mu)))
    return out

def table(tag, S):
    R = np.array([band_rms(z, d, r) for r, z, d in S]); dxs = np.array([1 / r for r, _, _ in S])
    print("%s: RMS of the nodal residual on bands " % tag + " | ".join("[%g, %g]" % (a, b) for a, b, _, _ in BANDS))
    for k, (r, _, _) in enumerate(S):
        print("  dx = %.3f   " % (1 / r) + "  ".join("%.3e" % v for v in R[k]))
    for k in range(1, len(S)):
        print("  order %.3f -> %.3f  " % (dxs[k - 1], dxs[k])
              + "  ".join("%7.2f" % (np.log(R[k - 1, j] / R[k, j]) / np.log(dxs[k - 1] / dxs[k])) for j in range(len(BANDS))))
    return dxs, R

def conv(a, dxs, R, title):
    for j, (_, _, lab, c) in enumerate(BANDS):
        a.loglog(dxs, R[:, j], "D-", color=c, ms=5.0, lw=1.3, label=lab)
    x = np.array([0.5, 1.0])
    for y0 in (2e-2, 6e-5):
        a.loglog(x, y0 * x ** 2, ":", color="#55544f", lw=1.0)
        a.loglog(x, y0 * 0.45 * x, "-.", color="#b9b7ae", lw=0.9)
    a.text(0.99, 2e-2, r"$\propto\Delta x^{2}$", ha="right", va="bottom", fontsize=9, color="#55544f")
    a.text(0.99, 6e-5 * 0.45, r"$\propto\Delta x$", ha="right", va="top", fontsize=9, color="#8a8880")
    a.set_xticks([0.5, 1 / 1.5, 0.8, 1.0]); a.set_xticklabels(["1/2", "1/1.5", "1/1.25", "1"]); a.minorticks_off()
    a.set_xlabel(r"lattice spacing $\Delta x$ (physical units)"); a.set_title(title, fontsize=11)

if __name__ == "__main__":
    F.style()
    SH = series([(r, p) for _, r, p in RUNS], HH, HH.mean_U())
    SG = series(GRUNS, HG, HG.mean_U())
    dH, RH = table("Hann (C^1)", SH); dG, RG = table("Gaussian (smooth)", SG)
    fig = plt.figure(figsize=(15.5, 5.6) if F.LAYOUT != "left" else (15.5, 8.2))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.35, 1, 1], wspace=0.30)
    a0 = fig.add_subplot(gs[0]); a1 = fig.add_subplot(gs[1]); a2 = fig.add_subplot(gs[2], sharey=a1, sharex=a1)
    for r, z, d in SG:
        a0.plot(z, d, "o-", ms=3.0, lw=1.0, color=RCOL[r], mec="white", mew=0.3,
                label="Ludwig gauss, $\\Delta x=%s$" % ("1" if r == 1 else "1/%g" % r))
    for a, b, _, c in BANDS:
        a0.axvspan(a, b, color=c, alpha=0.07, lw=0)
    a0.axhline(0, color="#8a8880", lw=0.8); a0.set_xlim(-16, 16)
    a0.set_xlabel("z (physical units, from the mid-plane)"); a0.set_ylabel("Ludwig $-$ Gaussian analytic, $\\rho=0$")
    a0.set_title("Gaussian kernel: residual on the axis (bands shaded)", fontsize=11); a0.legend(fontsize=8.5, loc="upper right")
    conv(a1, dH, RH, "Hann-$n$ ($C^1$: $W''$ jumps at $|u|=4$)")
    conv(a2, dG, RG, "Gaussian $\\sigma=1.446$ (smooth)")
    a1.set_ylabel(r"RMS$_{\rm band}$ of the nodal residual, $\rho=0$"); a2.legend(fontsize=8.5, loc="center right", bbox_to_anchor=(1.0, 0.62), title="band")
    if F.LAYOUT == "left":
        fig.text(0.5, 0.0, r"Same physical system and lattices for both kernels ($\Delta x = 1/r$, $\lambda_D=6r$,"
              r" $\varepsilon\to\varepsilon/r$, $\rho_{el}\to\rho_{el}/r^3$). RMS over bands whose ends are cell faces of every lattice.",
              ha="center", va="top", fontsize=10)
    fig.savefig(os.path.join(F.ROOT, "graficos", F.out_name("C_refinamiento_hann_vs_gauss.png")), dpi=170, bbox_inches="tight")
