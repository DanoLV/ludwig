#!/usr/bin/env python3
"""Lattice refinement of the electrolyte capacitor with one Gaussian-like charge: Ludwig with Hann-8, -10, -12, -16 on
lattices dx = 1, 1/1.25, 1/1.5, 1/2 (same physical system, built by crear_refinado.py) vs the analytic solution with the
continuous Hann charge (grafico_C_hoja_hann.py), whose physical width n dx = 8 is the same for all runs.

Convergence panel: residual Ludwig - analytic on the z axis through the charge (rho = 0) at FIXED physical distances
z - z0 (cubic spline through the nodes: inside the support with the nodes |z - z0| < 4, outside with the nodes
z - z0 > 4 on the channel side), so every run is compared at the same points.

Physical units: L = 32, lambda_D = 6, z0 = -8 (d = 8 from the lower plate), q/(eps kT) = 10, plates grounded.
Ludwig node k (1-based, z) is at physical z = (k - 1.5)/r - L/2; the charge sits at x = y = 16 r.
Field line y = 0 (through the charge), x = rho: rho_lat = round(r rho), evaluated at the exact rho_lat/r.
Usage: python3 grafico_C_refinamiento.py ;  LAYOUT=left python3 grafico_C_refinamiento.py  (slide version, *_ppt.png)
"""
import os, glob, math
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import grafico_C_formulas as F
import grafico_C_hoja_hann as HH
import grafico_C_hoja_gauss as HG

L = F.L; Z0 = F.D
RUNS = [(8, 1.0, "C-electrolito-hann8-LD6/z_09.50"),
        (10, 1.25, "C-electrolito-hann10-LD7.5-ref1.25/z_11.50"),
        (12, 1.5, "C-electrolito-hann12-LD9-ref1.5/z_13.50"),
        (16, 2.0, "C-electrolito-hann16-LD12-ref2/z_17.50")]
COL = {8: "#eb6834", 10: "#2a78d6", 12: "#1baf7a", 16: "#8e44ad"}
PTS = ((0.0, "inside"), (2.0, "inside"), (4.5, "outside"), (6.0, "outside"), (12.0, "outside"), (20.0, "outside"))
PCOL = {0.0: "#c0392b", 2.0: "#e67e22", 4.5: "#2a78d6", 6.0: "#5dade2", 12.0: "#1baf7a", 20.0: "#7f8c8d"}
# interpolation-free metric: RMS of the nodal residual over bands [a, b] of z (physical). The ends are multiples of 4
# from the wall (z + 16 = 4j), which are cell faces of every lattice (4 = 4*1 = 5*0.8 = 6*(2/3) = 8*0.5): each band is
# tiled by whole cells, so the mean over its nodes is the midpoint rule of (1/(b-a)) int d^2 dz, O(dx^2) accurate.
BANDS = ((-12.0, -4.0, "inside $|z-z_0|\\leq 4$", "#c0392b"), (-4.0, 0.0, "edge $4\\leq z-z_0\\leq 8$", "#2a78d6"),
         (-16.0, -12.0, "plate gap $-8\\leq z-z_0\\leq -4$", "#e67e22"), (0.0, 16.0, "far $8\\leq z-z_0\\leq 24$", "#1baf7a"))

def band_rms(z, d, r):
    """RMS of the nodal residual on each band (no interpolation); checks the band holds a whole number of cells."""
    out = []
    for a, b, _, _ in BANDS:
        m = (z > a) & (z < b)
        assert m.sum() == int(round((b - a) * r)), (a, b, r, m.sum())
        out.append(math.sqrt(np.mean(d[m] ** 2)))
    return np.array(out)

def at_points(z, d):
    """Residual d(z) at the fixed distances PTS (spline on the matching side of the support edge)."""
    from scipy.interpolate import CubicSpline
    u = z - Z0; out = []
    for p, side in PTS:
        m = (np.abs(u) < 4) if side == "inside" else (u > 4)
        out.append(float(CubicSpline(u[m], d[m])(p)))
    return np.array(out)

def psi_files(path):
    """psi dumps of a run, or [] if the final step (N_cycles of its input) has not been written yet."""
    fs = sorted(glob.glob(os.path.join(F.ROOT, path, "psi-0*")))
    ncyc = [int(l.split()[1]) for l in open(os.path.join(F.ROOT, path, "input")) if l.startswith("N_cycles")][0]
    return fs if fs and int(os.path.basename(fs[-1])[4:13]) == ncyc else []

def load(fname, r):
    n = int(round(32 * r)); a = np.loadtxt(fname).reshape(n, n, n + 2)
    return a

def line(a, r, rho):
    """Ludwig psi on fluid nodes along y = y_charge, x = x_charge + round(r rho); physical z and exact rho."""
    n = int(round(32 * r)); ic = int(round(16 * r)) - 1; dx = int(round(r * rho))
    z = (np.arange(1, n + 3) - 1.5) / r - L / 2
    return z[1:n + 1], a[(ic + dx) % n, ic, 1:n + 1], dx / r

if __name__ == "__main__":
    F.style()
    muH = HH.mean_U(); muG = HG.mean_U()
    res = {}
    print(" n    r   step(last) | max|psi(N) - psi(3N/4)| | rho  | max|Ludwig - Hann| out / in kernel | max|Ludwig - Gauss| out / in")
    for n, r, path in RUNS:
        fs = psi_files(path)
        if not fs:
            print("%2d  (no data yet)" % n); continue
        a = load(fs[-1], r); a_prev = load(fs[-2], r) if len(fs) > 2 else None
        drift = np.abs(a - a_prev).max() if a_prev is not None else float("nan")
        for rho in (0, 6, 10):
            z, ps, rr = line(a, r, rho)
            dh = ps - HH.phi(rr, z, 0, 0, muH); dg = ps - HG.phi(rr, z, 0, 0, muG)
            out = np.abs(z - Z0) > 4.0
            res[(n, rho)] = (z, ps, rr, dh, dg)
            print("%2d %5.2f %s |       %.2e        | %5.2f |  %.2e / %.2e |  %.2e / %.2e"
                  % (n, r, os.path.basename(fs[-1])[4:13], drift, rr, np.abs(dh[out]).max(), np.abs(dh[~out]).max(),
                     np.abs(dg[out]).max(), np.abs(dg[~out]).max()))
    done = [x for x in RUNS if (x[0], 0) in res]

    if F.LAYOUT == "left":
        fig = plt.figure(figsize=(24.0, 9.0)); gs = fig.add_gridspec(1, 2, width_ratios=[0.95, 1.6], wspace=0.08)
        F._AXT[0] = fig.add_subplot(gs[0]); F._AXT[0].axis("off"); F._YMAX[0] = 1.0
        g2 = gs[1].subgridspec(2, 3, height_ratios=[1.5, 1], width_ratios=[1.5, 1, 1], hspace=0.30, wspace=0.36)
    else:
        fig = plt.figure(figsize=(15.0, 7.4)); g2 = fig.add_gridspec(2, 3, height_ratios=[1.5, 1], width_ratios=[1.5, 1, 1],
                                                                     hspace=0.30, wspace=0.32)
    a0 = fig.add_subplot(g2[0, 0]); a1 = fig.add_subplot(g2[1, 0], sharex=a0); a2 = fig.add_subplot(g2[:, 1])
    a3 = fig.add_subplot(g2[:, 2], sharex=a2)
    if F.LAYOUT != "left": fig.subplots_adjust(left=0.08, right=0.98, top=0.95, bottom=0.08)   # before reading positions
    zz = np.linspace(-L / 2, L / 2, 641)
    a0.plot(zz, HH.phi(0.0, zz, 0, 0, muH), "-", color="#2b2a27", lw=1.5, label="analytic, continuous Hann (width 8)", zorder=2)
    a0.plot(zz, HG.phi(0.0, zz, 0, 0, muG), "--", color="#8a8880", lw=1.1, dashes=(4, 3), label="analytic, Gaussian ($\\sigma=1.446$)", zorder=2)
    for n, r, _ in done:
        z, ps, rr, dh, dg = res[(n, 0)]
        a0.plot(z, ps, "o", ms=4.0, color=COL[n], mec="white", mew=0.4, zorder=4,
                label="Ludwig Hann-%d, $\\Delta x=1/%g$" % (n, r) if r != 1 else "Ludwig Hann-8, $\\Delta x=1$")
        a1.plot(z, dh, "o-", ms=3.0, lw=1.0, color=COL[n], mec="white", mew=0.3)
    a0.set_xlim(-16, 0); a0.set_ylabel(r"$\psi = e\phi/k_BT$ at $\rho=0$")
    a0.set_title("Same physical system, %d lattice spacings (grounded plates, $\\lambda_D=6$)" % len(done), fontsize=11, pad=58)
    a0.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=2, fontsize=8.5); plt.setp(a0.get_xticklabels(), visible=False)
    a1.axvspan(Z0 - 4, Z0 + 4, color="#f1efe8", zorder=0); a1.axhline(0, color="#8a8880", lw=0.8)
    a1.set_ylabel("Ludwig $-$ Hann analytic"); a1.set_xlabel("z (physical units, from the mid-plane)")
    a1.text(Z0, 0.96, "kernel support", transform=a1.get_xaxis_transform(), ha="center", va="top", fontsize=8, color="#6b6a64")
    # convergence on the z axis through the charge (rho = 0), at fixed physical distances z - z0
    dxs = np.array([1 / r for _, r, _ in done])
    P = np.array([at_points(res[(n, 0)][0], res[(n, 0)][3]) for n, _, _ in done])        # (runs, points)
    print("rho = 0, residual at fixed z - z0:   " + "  ".join("%9.1f" % p for p, _ in PTS))
    for k, (n, r, _) in enumerate(done):
        print("  Hann-%-2d dx = %.3f            " % (n, 1 / r) + "  ".join("%+.2e" % v for v in P[k]))
    for k in range(1, len(done)):
        print("  order %.3f -> %.3f            " % (dxs[k - 1], dxs[k])
              + "  ".join("%9.2f" % (np.log(P[k - 1, j] / P[k, j]) / np.log(dxs[k - 1] / dxs[k])) for j in range(len(PTS))))
    for j, (p, side) in enumerate(PTS):
        a2.loglog(dxs, np.abs(P[:, j]), "o-" if side == "inside" else "s--", color=PCOL[p], ms=5.5, lw=1.3,
                  mfc=PCOL[p] if side == "inside" else "white", label="$z-z_0=%g$%s" % (p, " (inside)" if side == "inside" else ""))
    if len(done) > 1:
        x = np.array([dxs.min(), dxs.max()])
        for y0 in (np.abs(P[0, 0]) * 1.5, np.abs(P[0, 5]) * 0.6):
            a2.loglog(x, y0 * (x / dxs[0]) ** 2, ":", color="#55544f", lw=1.0)
        a2.text(dxs[0] * 0.97, np.abs(P[0, 0]) * 1.5, r"$\propto\Delta x^{2}$", ha="right", va="bottom", fontsize=9, color="#55544f")
    a2.set_xlabel(r"lattice spacing $\Delta x$ (physical units)"); a2.set_ylabel(r"$|$Ludwig $-$ Hann analytic$|$ on the axis $\rho=0$")
    a2.set_xticks(dxs); a2.set_xticklabels(["1" if r == 1 else "1/%g" % r for _, r, _ in done]); a2.minorticks_off()
    a2.legend(loc="center right", fontsize=8.5, title="distance from the charge", title_fontsize=8.5)
    a2.set_title("At fixed points (spline of the residual)", fontsize=11)
    # interpolation-free alternative: RMS of the nodal residual over bands bounded by common cell faces
    R = np.array([band_rms(res[(n, 0)][0], res[(n, 0)][3], r) for n, r, _ in done])        # (runs, bands)
    print("rho = 0, RMS of the nodal residual on bands: " + " | ".join("[%g, %g]" % (a, b) for a, b, _, _ in BANDS))
    for k, (n, r, _) in enumerate(done):
        print("  Hann-%-2d dx = %.3f            " % (n, 1 / r) + "  ".join("%.3e" % v for v in R[k]))
    for k in range(1, len(done)):
        print("  order %.3f -> %.3f            " % (dxs[k - 1], dxs[k])
              + "  ".join("%9.2f" % (np.log(R[k - 1, j] / R[k, j]) / np.log(dxs[k - 1] / dxs[k])) for j in range(len(BANDS))))
    for j, (a, b, lab, c) in enumerate(BANDS):
        a3.loglog(dxs, R[:, j], "D-", color=c, ms=5.0, lw=1.3, label=lab)
    if len(done) > 1:
        for y0 in (R[0, 0] * 1.5, R[0, 3] * 0.6):
            a3.loglog(x, y0 * (x / dxs[0]) ** 2, ":", color="#55544f", lw=1.0)
        a3.text(dxs[0] * 0.97, R[0, 0] * 1.5, r"$\propto\Delta x^{2}$", ha="right", va="bottom", fontsize=9, color="#55544f")
    a3.set_xlabel(r"lattice spacing $\Delta x$ (physical units)")
    a3.set_ylabel(r"RMS$_{\rm band}$ of the nodal residual, $\rho=0$")
    a3.set_xticks(dxs); a3.set_xticklabels(["1" if r == 1 else "1/%g" % r for _, r, _ in done]); a3.minorticks_off()
    a3.legend(loc="lower right", fontsize=8.5, title="band", title_fontsize=8.5)
    a3.set_title("Over bands (nodes only, no interpolation)", fontsize=11)
    fig.text(0.5 * (a2.get_position().x0 + a3.get_position().x1), a2.get_position().y1 + 0.045,
             "Convergence on the axis through the charge", ha="center", va="bottom", fontsize=12)
    fsz = 10.5; put = lambda y, s, **kw: F.put(fig, 0.5, y, s, ha="center", va="top", **kw)
    if F.LAYOUT != "left": put = lambda y, s, **kw: F.put(fig, 0.5, y * 0.75 - 0.02, s, ha="center", va="top", **kw)
    put(-0.02, "Lattice refinement at fixed physical system (r = n/8)", fontsize=fsz, fontweight="bold")
    put(-0.08, r"$\Delta x=1/r,\quad \ell,\,L\to 32r,\quad z_0\to r z_0,\quad \lambda_D\to 6r,\quad"
        r"\varepsilon\to\varepsilon/r\ (q/\varepsilon k_BT\to 10r),\quad \rho_{el}\to\rho_{el}/r^{3},\quad N_{steps}\to 5000\,r^{2}$", fontsize=fsz)
    put(-0.15, r"Hann-$n$ with $n=8r$: physical support $n\,\Delta x=8$ for all runs $\Rightarrow$ one analytic curve", fontsize=fsz)
    put(-0.22, r"$\phi=\phi_C(z)+U+\dfrac{\kappa L}{2\tanh(\kappa L/2)}\langle U\rangle\left[1-\dfrac{\cosh\kappa z}{\cosh(\kappa L/2)}\right]$",
        fontsize=fsz)
    put(-0.32, r"$\rho_{\mathrm{part}}=q\,W(x)W(y)W(z-z_0),\qquad W(u)=\dfrac{1}{8}\left[1+\cos\dfrac{2\pi u}{8}\right],\ |u|\leq 4$",
        fontsize=fsz)
    put(-0.42, r"(full formulas: grafico_C_hoja_hann.py). Ludwig places the charge on nodes and uses the 27-point Laplacian;"
        r" both errors vanish as $\Delta x\to 0$.", fontsize=fsz - 0.5)
    name = "C_refinamiento_hann.png"
    fig.savefig(os.path.join(F.ROOT, "graficos", F.out_name(name)), dpi=170, bbox_inches="tight")
