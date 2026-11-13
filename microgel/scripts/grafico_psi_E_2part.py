#!/usr/bin/env python3
"""Potential and electric field of two fixed charges (r = 4 physical), same figure as
TEST-HIDRO-dx1/scripts/grafico_psi_E_hidro_no.py but for any lattice (2026-11-11, user's request: Hann-16 and Hann-8 WITH
hydrodynamics, to compare with the Hann-8 run without hydrodynamics).

Cases (python3 grafico_psi_E_2part.py <case>):
  hann8_hidro : SLIDES-LD8-hann08-corr_no-N10000/pos_12.00_24.00_24.00, dx = 1, box 48^3, Hann-8, charges x = 8, 12,
                y = z = 24 (on nodes), step 10000, exe 1c30df4b (old), hydrodynamics on.
  hann16_hidro: SLIDES-LD8-hann16-corr_no-dxhalf/r_04, dx = 1/2, box 96^3, Hann-16, charges x = 15.5, 23.5,
                y = z = 47.5 (BETWEEN nodes), step 20000, exe 4f0c20e7, hydrodynamics on.
Everything is plotted in PHYSICAL units (dx = 1 lattice): x_phys = 0.5 + (X - 0.5) dx (X = Ludwig coordinate),
E per physical length = (efield / kT) / dx; psi in e psi/kT (lattice independent).
Line: the node line y = z = YL closest to the charge axis. For dx = 1 it is the axis itself; for dx = 1/2 the axis
y = z = 47.5 has no nodes, the line y = z = 47 is used (offset -0.25 physical in y and z; choice of this script).
Reference (as error_capas_2part.py): linear DH, periodic box, continuous Hann sources, lattice units of each lattice
(l_B / dx, kappa dx, Hann width 8 / dx), |n| <= M = 3N per direction:
  on the nodes: dh_nodes (zero node mean);
  on a fine x grid along the line (offsets dy = YL - Y_c, dz = YL - Z_c from the charges):
    phi(x) = (4 pi l_B/V) sum_kx What(kx) S1(kx) sum_p cos(kx (x - X_p)),  E_x = (4 pi l_B/V) sum_kx kx What S1 sum_p sin(.)
    S1(kx) = sum_{ky,kz} What(ky) What(kz) cos(ky dy) cos(kz dz) / (k^2 + kappa^2),
  its constant fixed to match dh_nodes on the nodes (printed: spread of the constant and E_x mismatch, ~ round-off).
Rows: profiles on the line (Ludwig points, reference line); pointwise error Ludwig - reference on the nodes (signed,
absolute); maps on the node plane z = YL: psi, and |E_xy| with (E_x, E_y) arrows every 3 physical units.
Usage: python3 grafico_psi_E_2part.py hann16_hidro   -> ../FUERZA-2PART-MEAS/PSI_E_<case>.png
"""
import os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
from error_capas_2part import dh_nodes, what, LB, KAPPA, KT, WHANN, OUT

CASES = {
    "hann8_hidro": dict(path="SLIDES-LD8-hann08-corr_no-N10000/pos_12.00_24.00_24.00", N=48, dx=1.0, st=10000,
                        pos=((12.0, 24.0, 24.0), (8.0, 24.0, 24.0)), YL=24.0,
                        lab="$\\Delta x = 1$, Hann-8, WITH hydrodynamics (step 10000, old exe 1c30df4b)"),
    "hann16_hidro": dict(path="SLIDES-LD8-hann16-corr_no-dxhalf/r_04", N=96, dx=0.5, st=20000,
                         pos=((23.5, 47.5, 47.5), (15.5, 47.5, 47.5)), YL=47.0,
                         lab="$\\Delta x = 1/2$, Hann-16, WITH hydrodynamics (step 20000)")}

case = sys.argv[1]; C = CASES[case]; N, dx, ST, POS, YL = C["N"], C["dx"], C["st"], C["pos"], C["YL"]
J = int(YL) - 1; RUN = os.path.join(HERE, "..", C["path"])
nh, lb, kap = WHANN / dx, LB / dx, KAPPA * dx
phys = lambda X: 0.5 + (np.asarray(X, dtype=float) - 0.5) * dx

psi = np.loadtxt(os.path.join(RUN, "psi-%09d.001-001" % ST)).reshape(N, N, N)
E = (np.loadtxt(os.path.join(RUN, "efield-%09d.001-001" % ST)) / KT / dx).reshape(N, N, N, 3)
phi, Er = dh_nodes(N, kap, lb, nh, POS, 3 * N); Er = Er / dx
X = np.arange(1, N + 1, dtype=float); x = phys(X)

M = 3 * N; n1 = np.arange(-M, M + 1); k1 = 2 * np.pi * n1 / N; w1 = what(k1, nh)
dy, dz = YL - POS[0][1], YL - POS[0][2]; assert all(p[1] == POS[0][1] and p[2] == POS[0][2] for p in POS)
ky, kz = np.meshgrid(k1, k1, indexing="ij"); wyz = np.outer(w1 * np.cos(k1 * dy), w1 * np.cos(k1 * dz))
S1 = np.array([np.sum(wyz / (kx * kx + ky * ky + kz * kz + kap ** 2)) for kx in k1]); A = 4 * np.pi * lb / N ** 3 * w1 * S1
def line(Xs):
    c = sum(np.cos(np.outer(Xs - p[0], k1)) for p in POS); s = sum(np.sin(np.outer(Xs - p[0], k1)) for p in POS)
    return c @ A, (s @ (k1 * A)) / dx
Xf = np.linspace(0.5, N + 0.5, 50 * N + 1); phif, Exf = line(Xf); phin, Exn = line(X)
c0 = phi[:, J, J] - phin; phif = phif + c0.mean(); xf = phys(Xf)
print("%s: line y = z = %g (lattice) = %g (physical)" % (case, YL, phys(YL)))
print("fine reference: spread of the constant on the nodes %.1e, max|E_x(line) - E_x(dh_nodes)| on the nodes %.1e"
      % (np.ptp(c0), np.abs(Exn - Er[:, J, J, 0]).max()))
for name, a, b in (("psi", psi[:, J, J], phi[:, J, J]), ("E_x", E[:, J, J, 0], Er[:, J, J, 0])):
    print("%s on the line: max|Ludwig - ref| = %.3e, max|ref| = %.3e" % (name, np.abs(a - b).max(), np.abs(b).max()))

INK = "#2b2a27"; REF = "#8a8880"; L = N * dx; xc = [phys(p[0]) for p in POS]
fig = plt.figure(figsize=(11, 13))
gs = fig.add_gridspec(3, 2, height_ratios=(1, 0.6, 1.15), hspace=0.32, wspace=0.28)
for c, (lab, a, b, bn) in enumerate(((r"$e\psi/k_BT$", psi[:, J, J], phif, phi[:, J, J]),
                                     (r"$eE_x/k_BT$ (per physical length)", E[:, J, J, 0], Exf, Er[:, J, J, 0]))):
    ax = fig.add_subplot(gs[0, c])
    ax.plot(xf, b, "-", color=REF, lw=1.0, label="reference: periodic DH, continuous Hann (width 8)")
    ax.plot(x, a, "o", color=INK, ms=3.5 if dx == 1 else 2.5, label="Ludwig, Hann-%d (%d nodes = width 8, $\\Delta x = %g$)" % (nh, nh, dx))
    ax.set_ylabel(lab); ax.legend(fontsize=7.5, frameon=False, loc="lower left", bbox_to_anchor=(0, 1.01), borderaxespad=0)   # outside the data
    ae = fig.add_subplot(gs[1, c], sharex=ax)
    ae.plot(x, a - bn, "o-", color=INK, ms=3.5 if dx == 1 else 2.5, lw=0.8)
    ae.set_ylabel("Ludwig $-$ reference"); ae.set_xlabel("physical x (line $y = z = %g$)" % phys(YL))
    plt.setp(ax.get_xticklabels(), visible=False)
    for q in (ax, ae):
        for p in xc:
            q.axvline(p, color=REF, lw=0.8, ls="--", zorder=0)
        q.axhline(0, color=REF, lw=0.6); q.set_xlim(0, L)
        for s in ("top", "right"): q.spines[s].set_visible(False)
ext = (phys(0.5), phys(N + 0.5), phys(0.5), phys(N + 0.5))   # pixel edges: node X spans X +- 1/2 (lattice) -> physical
# (2026-11-11 fix: was (0, L, 0, L), which shifted the maps by half a node with respect to the "+" marks and the arrows)
m0 = ax = fig.add_subplot(gs[2, 0]); v = np.abs(psi[:, :, J]).max()
im = ax.imshow(psi[:, :, J].T, origin="lower", extent=ext, cmap="RdBu_r", vmin=-v, vmax=v)
fig.colorbar(im, ax=ax, shrink=0.85, label=r"$e\psi/k_BT$")
ax.set_title(r"$\psi$, plane $z = %g$" % phys(YL), fontsize=10)
m1 = ax = fig.add_subplot(gs[2, 1]); Em = np.linalg.norm(E[:, :, J, :2], axis=-1)
im = ax.imshow(Em.T, origin="lower", extent=ext, cmap="Blues")
fig.colorbar(im, ax=ax, shrink=0.85, label=r"$e|E_{xy}|/k_BT$")
step = int(round(3 / dx)); s = slice(int(round(1 / dx)), N, step); XX, YY = np.meshgrid(x[s], x[s], indexing="ij")
ax.quiver(XX, YY, E[s, s, J, 0], E[s, s, J, 1], color=INK, scale_units="xy", angles="xy", width=0.003,
          scale=Em.max() / 2.5)     # arrow length proportional to |E|: the largest is 2.5 physical units long
ax.set_title(r"$|E|$ in the plane and $(E_x, E_y)$ arrows, $z = %g$" % phys(YL), fontsize=10)
for a in (m0, m1):
    a.set_xlabel("physical x"); a.set_ylabel("physical y")
    for p in POS:
        a.plot(phys(p[0]), phys(p[1]), "+", color="k", ms=8, mew=1.5)
fig.suptitle(r"Two charges $q = +1$ at physical $x = 8, 12$, $y = z = 24$; " + C["lab"] + r", $\lambda_D = 8$, box 48",
             fontsize=10.5, x=0.06, y=0.925, ha="left")   # just above the legends (default y = 0.98 left a large gap)
out = os.path.join(OUT, "PSI_E_%s.png" % case)
fig.savefig(out, dpi=200, bbox_inches="tight"); print(out)
