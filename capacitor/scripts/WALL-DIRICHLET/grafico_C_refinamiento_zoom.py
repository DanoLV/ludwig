#!/usr/bin/env python3
"""Zoom of grafico_C_refinamiento.py on |z - z0| <= 2 (physical units) along the z axis through the charge (rho = 0):
potential psi and field E_z = -dpsi/dz for Ludwig Hann-8, -10, -12, -16 vs the analytic continuous-Hann solution
(grafico_C_hoja_hann.py), with the residuals below. With the argument "full" the whole channel -16 <= z <= 16.

On the axis E_x = E_y = 0 by symmetry, so E_z is the whole field.
Ludwig E_z: the field Ludwig itself writes (efield-*, psi_electric_field in psi_gradients.c): 27-point gradient stencil
(electrokinetics_solver_stencil 27), difference form with the wall links doubled (Asta walls). The file is in lattice
units of the potential (kT/e) per lattice spacing: E_reduced, physical = efield / kT * r, kT = 1e-5 (input "temperature").
Third row: relative error |Ludwig - analytic| / |analytic| (%, log scale). It grows where psi -> 0 at the grounded
plates; at z0 it stays finite because E_z and its error are both odd about the charge (both vanish linearly).
Fourth row: the same relative error with its sign, (Ludwig - analytic) / analytic in %, on a linear scale.
For reference the residual panel also shows the plain central difference -(psi_{k+1} - psi_{k-1})/(2 dx) (open, faint).
Analytic E_z: central difference of the analytic phi with h = 1e-4 (error ~1e-8, negligible).
Usage: python3 grafico_C_refinamiento_zoom.py [full]   (-> C_refinamiento_hann_zoom.png / C_refinamiento_hann_full.png)
"""
import os, sys, glob
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import grafico_C_formulas as F
import grafico_C_hoja_hann as HH
from grafico_C_refinamiento import RUNS, COL, Z0, L, psi_files, load, line

FULL = "full" in sys.argv[1:]
KT = 1.0e-5                                           # "temperature" in the inputs (efield file is in units kT/e)

def ez_file(path, fname_psi, r):
    """E_z (reduced, physical units) on the axis from Ludwig's efield dump of the same step as fname_psi."""
    n = int(round(32 * r)); ic = int(round(16 * r)) - 1
    fe = os.path.join(F.ROOT, path, "efield-" + os.path.basename(fname_psi)[4:])
    e = np.loadtxt(fe).reshape(n, n, n + 2, 3)
    return e[ic, ic, 1:n + 1, 2] / KT * r
W = 2.0                                               # half width of the zoom (physical units)
ZLO, ZHI = (-L / 2, L / 2) if FULL else (Z0 - W, Z0 + W)

def ez_an(z, mu, h=1e-4):
    return -(HH.phi(0.0, z + h, 0, 0, mu) - HH.phi(0.0, z - h, 0, 0, mu)) / (2 * h)

if __name__ == "__main__":
    F.style()
    muH = HH.mean_U()
    fig, ax = plt.subplots(4, 2, figsize=(11.0, 11.8), sharex=True, gridspec_kw=dict(height_ratios=[1.5, 1, 1, 1], hspace=0.08, wspace=0.28))
    zz = np.linspace(ZLO, ZHI, 1281 if FULL else 401)
    ax[0, 0].plot(zz, HH.phi(0.0, zz, 0, 0, muH), "-", color="#2b2a27", lw=1.5, label="analytic, continuous Hann (width 8)")
    ax[0, 1].plot(zz, ez_an(zz, muH), "-", color="#2b2a27", lw=1.5)
    print(" n   dx     | max|dpsi| in range | max|dE_z| file (27-pt) | central diff")
    for n, r, path in RUNS:
        fs = psi_files(path)
        if not fs: continue
        z, ps, _ = line(load(fs[-1], r), r, 0)
        ez = ez_file(path, fs[-1], r)                                       # Ludwig's own 27-point field
        ecd = -(ps[2:] - ps[:-2]) * r / 2; zc = z[1:-1]                     # plain central difference (reference)
        m = (z >= ZLO - 1e-9) & (z <= ZHI + 1e-9); mc = (zc >= ZLO - 1e-9) & (zc <= ZHI + 1e-9)
        pa = HH.phi(0.0, z[m], 0, 0, muH); ea = ez_an(z[m], muH)
        dp = ps[m] - pa; de = ez[m] - ea; dcd = ecd[mc] - ez_an(zc[mc], muH)
        lab = "Ludwig Hann-%d, $\\Delta x=1/%g$" % (n, r) if r != 1 else "Ludwig Hann-8, $\\Delta x=1$"
        kw = dict(ms=5.0, color=COL[n], mec="white", mew=0.5)
        ax[0, 0].plot(z[m], ps[m], "o", label=lab, **kw); ax[0, 1].plot(z[m], ez[m], "o", **kw)
        ax[1, 0].plot(z[m], dp, "o-", lw=1.0, **kw); ax[1, 1].plot(z[m], de, "o-", lw=1.0, **kw)
        ax[1, 1].plot(zc[mc], dcd, "o", ms=3.5, mfc="none", mec=COL[n], mew=0.7, alpha=0.45, zorder=1)
        ax[2, 0].semilogy(z[m], 100 * np.abs(dp / pa), "o-", lw=1.0, **kw)
        ax[2, 1].semilogy(z[m], 100 * np.abs(de / ea), "o-", lw=1.0, **kw)
        ax[3, 0].plot(z[m], 100 * dp / pa, "o-", lw=1.0, **kw)
        ax[3, 1].plot(z[m], 100 * de / ea, "o-", lw=1.0, **kw)
        print("%2d  %.3f  |     %.2e      |    %.2e     |     %.2e" % (n, 1 / r, np.abs(dp).max(), np.abs(de).max(), np.abs(dcd).max()))
    for a in ax[1]:
        a.axhline(0, color="#8a8880", lw=0.8)
    for a in ax[2]:
        a.set_ylim(1e-3 if FULL else 1e-2, 1e2)
    for a in ax[3]:
        a.axhline(0, color="#8a8880", lw=0.8); a.set_xlabel("z (physical units, from the mid-plane)")
    ax[2, 0].set_ylabel(r"$\psi$: |relative error| (%)"); ax[2, 1].set_ylabel(r"$E_z$: |relative error| (%)")
    ax[3, 0].set_ylabel(r"$\psi$: relative error (%)"); ax[3, 1].set_ylabel(r"$E_z$: relative error (%)")
    for a in ax.ravel():
        a.axvline(Z0, color="#b9b7ae", lw=0.8, ls=":"); a.set_xlim(ZLO - 0.05 * (1 + 3 * FULL), ZHI + 0.05 * (1 + 3 * FULL))
        if FULL: a.axvspan(Z0 - 4, Z0 + 4, color="#f1efe8", zorder=0)
    ax[0, 0].set_ylabel(r"$\psi = e\phi/k_BT$ at $\rho=0$"); ax[0, 1].set_ylabel(r"$E_z=-\partial_z\psi$ at $\rho=0$")
    ax[1, 0].set_ylabel(r"$\psi$: Ludwig $-$ analytic"); ax[1, 1].set_ylabel(r"$E_z$: Ludwig $-$ analytic")
    ax[1, 1].plot([], [], "o", ms=3.5, mfc="none", mec="#8a8880", label="central difference of $\\psi$ (reference)")
    ax[1, 1].plot([], [], "o-", ms=4, color="#8a8880", label="Ludwig efield (27-point stencil)")
    ax[1, 1].legend(fontsize=7.5, loc="lower right")
    ax[0, 0].legend(loc="upper right" if FULL else "lower center", bbox_to_anchor=(1.04, 1.02) if FULL else None, fontsize=8 if FULL else 8.5, ncol=1)
    if FULL: ax[1, 0].text(Z0, 0.96, "kernel support", transform=ax[1, 0].get_xaxis_transform(), ha="center", va="top", fontsize=8, color="#6b6a64")
    fig.suptitle((r"Whole channel $-16\leq z\leq 16$" if FULL else r"Zoom $|z-z_0|\leq 2$")
                 + r" on the axis through the charge ($z_0=-8$, $\lambda_D=6$, grounded plates)", fontsize=11, y=0.95)
    F.put(fig, 0.5, 0.0, r"Ludwig $E_z$: its own output (efield, 27-point gradient stencil, wall links doubled), $E/k_BT$ in"
          r" physical units; analytic $E_z$: exact derivative. On the axis $E_x=E_y=0$ by symmetry.", ha="center", va="top", fontsize=9.5)
    fig.savefig(os.path.join(F.ROOT, "graficos", "C_refinamiento_hann_full.png" if FULL else "C_refinamiento_hann_zoom.png"), dpi=170, bbox_inches="tight")
