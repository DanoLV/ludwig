#!/usr/bin/env python3
"""Potential and electric field of the dx = 1, Hann-8, r = 4 run WITHOUT hydrodynamics (2026-11-11, user's request).

Run: ../hidro_no_exe45dc94bc (box 48^3, charges q = +1 at x = 8 and 12, y = z = 24; exe 45dc94bc), last dump (step 30000).
Units: psi dump in e psi/kT; efield dump / kT -> e E/kT per lattice unit (dx = 1: per physical length).
Reference (user's choice): linear Debye-Hueckel, periodic box, continuous Hann-8 sources = dh_nodes of
error_capas_2part.py (zero node mean, E_DH = -grad phi_DH), evaluated on the nodes, M = 3N modes per direction.
Top: profiles on the line y = z = 24 (through both charges), Ludwig (points) and reference (line).
Middle: pointwise error on that line, Ludwig - reference on the nodes (signed, absolute).
Bottom: maps on the plane z = 24 (contains both charges): psi, and |E| with arrows of (E_x, E_y) every 3 nodes.
Node at array index j <-> Ludwig coordinate j + 1.
Usage: python3 grafico_psi_E_hidro_no.py   -> ../graficos/psi_E_hidro_no_paso30000.png
"""
import os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "scripts"))
from error_capas_2part import dh_nodes, LB, KAPPA, KT, WHANN

N = 48; ST = 30000; POS = ((12.0, 24.0, 24.0), (8.0, 24.0, 24.0)); J = 23     # y = z = 24 -> index 23
RUN = os.path.join(HERE, "..", "hidro_no_exe45dc94bc")

psi = np.loadtxt(os.path.join(RUN, "psi-%09d.001-001" % ST)).reshape(N, N, N)
E = (np.loadtxt(os.path.join(RUN, "efield-%09d.001-001" % ST)) / KT).reshape(N, N, N, 3)
phi, Er = dh_nodes(N, KAPPA, LB, WHANN, POS, 3 * N)
x = np.arange(1, N + 1)
for name, a, b in (("psi", psi[:, J, J], phi[:, J, J]), ("E_x", E[:, J, J, 0], Er[:, J, J, 0])):
    print("%s on the axis: max|Ludwig - ref| = %.3e, max|ref| = %.3e" % (name, np.abs(a - b).max(), np.abs(b).max()))

# 2026-11-11 (user's request): reference on a FINE x grid along the axis y = z = 24 (it passes through both charges, so the
# ky, kz phases vanish): same modes as dh_nodes (|n| <= M = 3N per direction),
#   phi(x) = (4 pi l_B / V) sum_kx What(kx) S1(kx) sum_p cos(kx (x - X_p)),   S1(kx) = sum_{ky,kz} What(ky) What(kz)/(k^2 + kappa^2)
#   E_x(x) = -dphi/dx = (4 pi l_B / V) sum_kx kx What(kx) S1(kx) sum_p sin(kx (x - X_p)).
# The additive constant is fixed so that phi matches dh_nodes on the nodes (zero node mean); printed: the spread of that
# constant over the nodes and the E_x mismatch on the nodes (both must be ~ round-off).
from error_capas_2part import what
M = 3 * N; n1 = np.arange(-M, M + 1); k1 = 2 * np.pi * n1 / N; w1 = what(k1, WHANN)
ky, kz = np.meshgrid(k1, k1, indexing="ij"); wyz = np.outer(w1, w1)
S1 = np.array([np.sum(wyz / (kx * kx + ky * ky + kz * kz + KAPPA ** 2)) for kx in k1]); A = 4 * np.pi * LB / N ** 3 * w1 * S1
def line(xs):
    xs = np.atleast_1d(xs)
    c = sum(np.cos(np.outer(xs - p[0], k1)) for p in POS); s = sum(np.sin(np.outer(xs - p[0], k1)) for p in POS)
    return c @ A, s @ (k1 * A)
def line_one(xs, p):
    """Same as line() for ONE charge at p (2026-11-11): its own periodic DH field on the axis."""
    xs = np.atleast_1d(xs); c = np.cos(np.outer(xs - p[0], k1)); s = np.sin(np.outer(xs - p[0], k1))
    return c @ A, s @ (k1 * A)
xf = np.linspace(0.5, N + 0.5, 2401); phif, Exf = line(xf); phin, Exn = line(x.astype(float))
c0 = phi[:, J, J] - phin; phif = phif + c0.mean()
print("fine reference: spread of the constant on the nodes %.1e, max|E_x(line) - E_x(dh_nodes)| on the nodes %.1e"
      % (np.ptp(c0), np.abs(Exn - Er[:, J, J, 0]).max()))

INK = "#2b2a27"; REF = "#8a8880"
fig = plt.figure(figsize=(11, 17))
gs = fig.add_gridspec(4, 2, height_ratios=(1, 0.6, 1, 1.15), hspace=0.32, wspace=0.28)
# 2026-11-11 (user's request): row 2 = pointwise error on the axis nodes, Ludwig - reference (signed, absolute, same units
# as the profile above; reference = dh_nodes on the nodes).
for c, (lab, a, b, bn) in enumerate(((r"$e\psi/k_BT$", psi[:, J, J], phif, phi[:, J, J]),
                                     (r"$eE_x/k_BT$ (per lattice unit)", E[:, J, J, 0], Exf, Er[:, J, J, 0]))):
    ax = fig.add_subplot(gs[0, c])
    ax.plot(xf, b, "-", color=REF, lw=1.0, label="reference: periodic DH, continuous Hann (width 8)")
    ax.plot(x, a, "o", color=INK, ms=3.5, label="Ludwig, Hann-8 (8 nodes = width 8, $\\Delta x = 1$), no hydrodynamics (step %d)" % ST)
    ax.set_ylabel(lab); ax.legend(fontsize=7.5, frameon=False, loc="lower left", bbox_to_anchor=(0, 1.01), borderaxespad=0)   # outside the data
    ae = fig.add_subplot(gs[1, c], sharex=ax)
    ae.plot(x, a - bn, "o-", color=INK, ms=3.5, lw=0.8)
    ae.set_ylabel("Ludwig $-$ reference"); ae.set_xlabel("x (lattice = physical, $y = z = 24$)")
    plt.setp(ax.get_xticklabels(), visible=False)
    for q in (ax, ae):
        for p in POS:
            q.axvline(p[0], color=REF, lw=0.8, ls="--", zorder=0)
        q.axhline(0, color=REF, lw=0.6); q.set_xlim(0.5, N + 0.5)
        for s in ("top", "right"): q.spines[s].set_visible(False)
# 2026-11-11 (user's request): row 3 = analytic solution of EACH charge alone and their sum (linear DH: superposition,
# the sum is exactly the reference of row 1), on the axis y = z = 24, fine grid. Constant of psi: each charge gets half
# of the constant c0 of the total (both charges are on nodes and identical, so their node means are equal; choice of
# this script); E_x needs no constant.
one = [line_one(xf, p) for p in POS]
for c, (lab, k) in enumerate(((r"$e\psi/k_BT$", 0), (r"$eE_x/k_BT$ (per lattice unit)", 1))):
    ax = fig.add_subplot(gs[2, c])
    for (p, col, nm) in zip(POS, ("#2a78d6", "#eb6834"), ("charge 1", "charge 2")):
        f = one[POS.index(p)][k] + (0.5 * c0.mean() if k == 0 else 0.0)
        ax.plot(xf, f, "-", color=col, lw=1.2, label="%s alone ($x = %g$)" % (nm, p[0]))
    tot = (phif, Exf)[k]
    ax.plot(xf, tot, "-", color=INK, lw=1.2, label="sum (= reference above)")
    for p in POS:
        ax.axvline(p[0], color=REF, lw=0.8, ls="--", zorder=0)
    ax.axhline(0, color=REF, lw=0.6); ax.set_xlim(0.5, N + 0.5); ax.set_ylabel(lab)
    ax.set_xlabel("x (lattice = physical, $y = z = 24$)")
    ax.set_title("analytic: periodic DH, continuous Hann (width 8)", fontsize=10, loc="left")
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    ax.legend(fontsize=7.5, frameon=False, loc="upper right")
chk = max(np.abs(one[0][0] + one[1][0] + c0.mean() - phif).max(), np.abs(one[0][1] + one[1][1] - Exf).max())
print("superposition check: max|charge 1 + charge 2 - total| = %.1e" % chk)
ext = (0.5, N + 0.5, 0.5, N + 0.5)
m0 = ax = fig.add_subplot(gs[3, 0]); v = np.abs(psi[:, :, J]).max()
im = ax.imshow(psi[:, :, J].T, origin="lower", extent=ext, cmap="RdBu_r", vmin=-v, vmax=v)
fig.colorbar(im, ax=ax, shrink=0.85, label=r"$e\psi/k_BT$")
ax.set_title(r"$\psi$, plane $z = 24$", fontsize=10)
m1 = ax = fig.add_subplot(gs[3, 1]); Em = np.linalg.norm(E[:, :, J, :2], axis=-1)
im = ax.imshow(Em.T, origin="lower", extent=ext, cmap="Blues")
fig.colorbar(im, ax=ax, shrink=0.85, label=r"$e|E_{xy}|/k_BT$")
s = slice(1, N, 3); X, Y = np.meshgrid(x[s], x[s], indexing="ij")
ax.quiver(X, Y, E[s, s, J, 0], E[s, s, J, 1], color=INK, scale_units="xy", angles="xy", width=0.003,
          scale=Em.max() / 2.5)     # arrow length proportional to |E|: the largest is 2.5 nodes long
ax.set_title(r"$|E|$ in the plane and $(E_x, E_y)$ arrows, $z = 24$", fontsize=10)
for a in (m0, m1):
    a.set_xlabel("x"); a.set_ylabel("y")
    for p in POS:
        a.plot(p[0], p[1], "+", color="k", ms=8, mew=1.5)
fig.suptitle(r"Two charges $q = +1$ at $x = 8, 12$; $\Delta x = 1$, Hann-8, $\lambda_D = 8$, box 48, no hydrodynamics",
             fontsize=11, x=0.06, y=0.915, ha="left")   # just above the legends (default y = 0.98 left a large gap)
out = os.path.join(HERE, "..", "graficos", "psi_E_hidro_no_paso%d.png" % ST)
os.makedirs(os.path.dirname(out), exist_ok=True); fig.savefig(out, dpi=200, bbox_inches="tight"); print(out)
