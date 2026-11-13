#!/usr/bin/env python3
"""Force between two fixed charges in the periodic box vs separation, dx = 2 and 1 (2026-11-09, user's request).

Runs (as error_capas_2part.py): dx = 1 SLIDES-LD8-hann08-corr_no-N10000 (box 48, Hann-8), dx = 2 SLIDES-LD8-hann04-corr_no-dx2
(box 24, Hann-4, same physical positions); q1 = q2 = +1, lambda_D = 8; charge 2 fixed, charge 1 at separation r along +x.
Ludwig force on charge 1 (physical units, kT per physical length times kT): F_L = kT Esub_x / dx, Esub_x of particle 1 at
the last step (proceced_data/particle_Esub.csv). Esub, not particle_force.csv: in the dx = 2 runs (old exe 1c30df4b) the
force of particle 2 (x = 4.25, edge colloid cell) is doubled by the halo bug fixed 2026-09-26; particle 1 and Esub are right.

Reference (user's choice: option 2, the same as for psi and E): linear Debye-Hueckel, periodic box L = 48, both charges spread
with the continuous Hann window of width 8; force = kernel-weighted field at charge 1, F_ref = kT q int W(r - r1) E_DH,x(r) d^3r.
With phi(r) = sum_k A(k) sum_p e^{i k.(r - r_p)}, A = (4 pi l_B/V) What3(k)/(k^2 + kappa^2), What3 = What(kx) What(ky) What(kz),
and int W(r - r1) e^{i k.r} d^3r = What3(k) e^{i k.r1}:
  F_ref/kT = (4 pi l_B/V) sum_k kx What3(k)^2 / (k^2 + kappa^2) sin(kx r)     (r = x1 - x2; the self term p = 1 is odd in k: 0)
           = (4 pi l_B/V) sum_kx kx What(kx)^2 sin(kx r) S(kx),  S(kx) = sum_{ky,kz} What(ky)^2 What(kz)^2/(k^2 + kappa^2)
l_B = 0.7957747, kappa = 0.124975 (as error_capas_2part.py). Modes |n| <= M per direction; checks printed: M -> 2M, and
F_ref against -dU/dr by central differences of the interaction energy U(r) = (4 pi l_B/V) sum_k What3^2/(k^2+kappa^2) cos(kx r).
Error e = (F_L - F_ref)/F_ref; observed order p = log2 |e(dx = 2)/e(dx = 1)| (two lattices: one number per r, it assumes
e ~ C dx^p already holds at dx = 1).
Not separated: nonlinear Poisson-Boltzmann part (Ludwig nonlinear, reference linear).
Usage: python3 error_fuerza_2part.py   -> ../FUERZA-2PART-MEAS/FUERZA_LD08.0_Hann-8_corr-no_dx2-1.png (+ table printed)
"""
import os, glob, math
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from error_capas_2part import what, LB, KAPPA, KT, WHANN, OUT, HERE

L = 48.0; M = 400
# 2026-11-11 (user's request): kernel scales marked on every force-vs-r plot. sigma = std. dev. of the continuous Hann window
# of width n = WHANN: sigma^2 = n^2 (1/12 - 1/(2 pi^2)) -> 1.44605 (= pm_sr_sigma of the Hann-8 runs). Lines: r = sigma/2 and
# r = sigma (user's choice: sigma = minimum separation taken as physically possible without real overlap), r = 8 (supports
# of width 8 stop overlapping).
SIG = WHANN * math.sqrt(1.0 / 12.0 - 1.0 / (2.0 * math.pi ** 2))
def kernel_lines(a):
    a.axvline(SIG / 2, color="#8a8880", lw=0.8, ls=":", zorder=0, label="$r = \\sigma/2 = %.3f$" % (SIG / 2))
    a.axvline(SIG, color="#8a8880", lw=0.8, ls="-.", zorder=0,
              label="$r = \\sigma = %.3f$: kernel std. dev." % SIG)
    # 2026-11-11 (user's request): multiples k sigma inside the kernel width (k sigma < 8: k = 2 ... 5), one legend entry
    ks = [k for k in range(2, 100) if k * SIG < WHANN]
    for i, k in enumerate(ks):
        a.axvline(k * SIG, color="#b4b2a9", lw=0.7, ls="-.", zorder=0,
                  label=("$r = k\\sigma$, $k = %d \\ldots %d$ ($k\\sigma < 8$)" % (ks[0], ks[-1])) if i == 0 else None)
    a.axvline(WHANN, color="#8a8880", lw=0.8, ls="--", zorder=0, label="$r = 8$: kernel supports (width 8) stop overlapping")

def S_of_kx(M):
    n = np.arange(-M, M + 1); k = 2 * np.pi * n / L; w2 = what(k, WHANN) ** 2
    ky, kz = np.meshgrid(k, k, indexing="ij"); wyz = np.outer(w2, w2)
    return k, w2, np.array([np.sum(wyz / (kx * kx + ky * ky + kz * kz + KAPPA ** 2)) for kx in k])

def F_ref(r, M=M):
    k, w2, S = S_of_kx(M)
    return KT * 4 * np.pi * LB / L ** 3 * np.sum((k * w2 * S)[:, None] * np.sin(np.outer(k, np.atleast_1d(r))), axis=0)

def U_ref(r, M=M):
    k, w2, S = S_of_kx(M)
    return KT * 4 * np.pi * LB / L ** 3 * np.sum((w2 * S)[:, None] * np.cos(np.outer(k, np.atleast_1d(r))), axis=0)

def esub_x(path, idx):
    rows = [l.strip().replace(",", ".").split(";") for l in open(os.path.join(path, "proceced_data", "particle_Esub.csv"))
            if not l.startswith("#")]
    st = max(int(v[0]) for v in rows)
    return [float(v[3]) for v in rows if int(v[0]) == st and int(v[1]) == idx][0], st

if __name__ == "__main__":
    runs = {1.0: {float(os.path.basename(p).split("_")[1]) - 8: p for p in glob.glob(os.path.join(HERE, "..", "SLIDES-LD8-hann08-corr_no-N10000", "pos_*"))},
            2.0: {float(os.path.basename(p)[2:]): p for p in glob.glob(os.path.join(HERE, "..", "SLIDES-LD8-hann04-corr_no-dx2", "r_*"))}}
    rs = np.array(sorted(runs[1.0])); assert list(rs) == sorted(runs[2.0])
    fr = F_ref(rs); fr2 = F_ref(rs, M // 2)
    h = 1e-4; du = -(U_ref(rs + h) - U_ref(rs - h)) / (2 * h)
    print("checks: max rel |F_ref(M) - F_ref(M/2)| %.1e | max rel |F_ref + dU/dr| %.1e (M = %d)"
          % (np.abs((fr - fr2) / fr).max(), np.abs((fr - du) / fr).max(), M))
    FL = {}
    for dx in (2.0, 1.0):
        FL[dx] = np.array([KT * esub_x(runs[dx][r], 1)[0] / dx for r in rs])
    e = {dx: (FL[dx] - fr) / fr for dx in FL}; p = np.log2(np.abs(e[2.0] / e[1.0]))
    print("  r  |     F_ref     | F_L dx = 2    e (%%)    | F_L dx = 1    e (%%)    |  p(2->1)")
    for i, r in enumerate(rs):
        print("%4.0f | %.6e | %.6e %+8.3f | %.6e %+8.3f | %6.2f" % (r, fr[i], FL[2.0][i], 100 * e[2.0][i], FL[1.0][i], 100 * e[1.0][i], p[i]))
    rr = np.linspace(0.05, 24, 480); fc = F_ref(rr)
    fig, ax = plt.subplots(3, 1, figsize=(8.5, 11.0), sharex=True, gridspec_kw=dict(hspace=0.08))
    ax[0].plot(rr, fc, "-", color="#8a8880", lw=1.0, label="reference: periodic DH, continuous Hann (width 8)")
    LAB = {2.0: ("$\\Delta x = 2$ (Hann-4)", "#eb6834", "^"), 1.0: ("$\\Delta x = 1$ (Hann-8)", "#2b2a27", "o")}
    for dx in (2.0, 1.0):
        lab, c, mk = LAB[dx]
        ax[0].plot(rs, FL[dx], mk, color=c, ms=5, label="Ludwig, " + lab)
        ax[1].plot(rs, 100 * e[dx], mk + "-", color=c, lw=1.0, ms=5, label=lab)
    ax[2].plot(rs, p, "s-", color="#2a78d6", lw=1.0, ms=5, label="$p(\\Delta x = 2 \\to 1)$")
    for q in (1, 2):
        ax[2].axhline(q, color="#8a8880", lw=0.8, ls=":")
    ax[1].axhline(0, color="#8a8880", lw=0.8)
    for a in ax:
        kernel_lines(a)
        if a is ax[2]:   # order panel: only the p entry (the vertical lines are explained in the panels above)
            a.legend(handles=a.lines[:1], fontsize=8, loc="upper left")
        else:
            a.legend(fontsize=8)
    ax[0].set_ylabel("$F_x$ on charge 1 (lattice units)")
    ax[1].set_ylabel("$(F_L - F_{ref}) / F_{ref}$ (%)")
    ax[2].set_ylabel("observed order $p$")
    ax[2].set_xlabel("separation $r$ between the charges")
    ax[0].set_title("Two charges $q = +1$, periodic box 48, $\\lambda_D = 8$, Hann of width 8, $k_BT = 10^{-5}$ (lattice units)", fontsize=10, loc="left")
    ax[0].set_xlim(0, 24)
    fig.savefig(os.path.join(OUT, "FUERZA_LD08.0_Hann-8_corr-no_dx2-1.png"), dpi=300, bbox_inches="tight")
    # 2026-11-11 (user's request): same figure + the dx = 1/2 force points (box 96, Hann-16, current exe 4f0c20e7; r = 1..8),
    # last step of each run (20000; r = 5, 6, 7 from the continuations -cont). Saved to a SEPARATE file (the plot above is kept).
    # Panel 3 (order p) unchanged: only p(2 -> 1).
    H = os.path.join(HERE, "..", "SLIDES-LD8-hann16-corr_no-dxhalf")
    rh = np.arange(1.0, 9.0)
    paths = [H + ("-cont" if r in (5, 6, 7) else "") + "/r_%02d" % r for r in rh]
    st = [esub_x(p, 1) for p in paths]; assert all(s == 20000 for _, s in st), st
    FLh = np.array([KT * v / 0.5 for v, _ in st]); eh = (FLh - F_ref(rh)) / F_ref(rh)
    print("dx = 1/2 (step 20000):"); [print("%4.0f | %.6e %+8.3f" % (r, f, 100 * x)) for r, f, x in zip(rh, FLh, eh)]
    lab = "$\\Delta x = 1/2$ (Hann-16, current exe)"
    ax[0].plot(rh, FLh, "D", color="#2a78d6", ms=5, label="Ludwig, " + lab)
    ax[1].plot(rh, 100 * eh, "D-", color="#2a78d6", lw=1.0, ms=5, label=lab)
    for a in ax[:2]:
        a.legend(fontsize=8)
    # fig.savefig(os.path.join(OUT, "FUERZA_LD08.0_Hann-8_corr-no_dx2-1-half.png"), dpi=300, bbox_inches="tight")
    # 2026-11-11 (user's request): the "-half" figure WITHOUT the order panel: drawn on its own 2-panel figure (force,
    # relative error) with the same series; the 3-panel figure above continues for "-half-int".
    fg, axh = plt.subplots(2, 1, figsize=(8.5, 7.6), sharex=True, gridspec_kw=dict(hspace=0.08))
    # 2026-11-11 (user's request): thinner reference line; vertical line at sigma/2 (user's choice of position), sigma = second
    # moment width of the continuous Hann window of width n = 8: sigma^2 = n^2 (1/12 - 1/(2 pi^2)) -> sigma = 1.44605
    # (= pm_sr_sigma of the Hann-8 runs); legend entries for both vertical lines.
    axh[0].plot(rr, fc, "-", color="#8a8880", lw=1.0, label="reference: periodic DH, continuous Hann (width 8)")
    for r_, F_, (lab_, c, mk) in ((rs, FL[2.0], LAB[2.0]), (rs, FL[1.0], LAB[1.0]), (rh, FLh, (lab, "#2a78d6", "D"))):
        axh[0].plot(r_, F_, mk, color=c, ms=5, label="Ludwig, " + lab_)
        axh[1].plot(r_, 100 * (F_ - F_ref(r_)) / F_ref(r_), mk + "-", color=c, lw=1.0, ms=5, label=lab_)
    axh[1].axhline(0, color="#8a8880", lw=0.8)
    for a in axh:
        kernel_lines(a)
        a.legend(fontsize=8)
    axh[0].set_ylabel("$F_x$ on charge 1 (lattice units)"); axh[1].set_ylabel("$(F_L - F_{ref}) / F_{ref}$ (%)")
    axh[1].set_xlabel("separation $r$ between the charges")
    axh[0].set_title("Two charges $q = +1$, periodic box 48, $\\lambda_D = 8$, Hann of width 8, $k_BT = 10^{-5}$ (lattice units)", fontsize=10, loc="left")
    axh[0].set_xlim(0, 24)
    fg.savefig(os.path.join(OUT, "FUERZA_LD08.0_Hann-8_corr-no_dx2-1-half.png"), dpi=300, bbox_inches="tight")
    # 2026-11-11 (user's request): + dx = 1 intermediate separations d = 4.1 ... 4.9 (SLIDES-LD8-hann08-corr_no-N10000-d4a5,
    # exe 4f0c20e7, step 10000), open circles. Saved to a SEPARATE file (the two plots above are kept).
    I = os.path.join(HERE, "..", "SLIDES-LD8-hann08-corr_no-N10000-d4a5")
    ri = np.round(4.0 + 0.1 * np.arange(1, 10), 1)
    st = [esub_x(os.path.join(I, "pos_%.2f_24.00_24.00" % (8 + r)), 1) for r in ri]; assert all(s == 10000 for _, s in st), st
    FLi = np.array([KT * v for v, _ in st]); ei = (FLi - F_ref(ri)) / F_ref(ri)
    print("dx = 1 intermediate (step 10000):"); [print("%4.1f | %.6e %+8.3f" % (r, f, 100 * x)) for r, f, x in zip(ri, FLi, ei)]
    # 2026-11-11 (user's request, same figure): + intermediate d = 8.2 ... 9.8 every 0.2 for dx = 1
    # (SLIDES-LD8-hann08-corr_no-N10000-d8a10, step 10000) and dx = 2 (SLIDES-LD8-hann04-corr_no-dx2-d8a10, step 3000,
    # F = kT Esub_x / 2), exe 45dc94bc (crear_2part_d8a10.py). Open markers = intermediate separations.
    r8 = np.round(8.0 + 0.2 * np.arange(1, 10), 1)
    st = [esub_x(os.path.join(HERE, "..", "SLIDES-LD8-hann08-corr_no-N10000-d8a10", "pos_%.2f_24.00_24.00" % (8 + r)), 1)
          for r in r8]; assert all(s == 10000 for _, s in st), st
    ri = np.concatenate([ri, r8]); FLi = np.concatenate([FLi, [KT * v for v, _ in st]]); ei = (FLi - F_ref(ri)) / F_ref(ri)
    st = [esub_x(os.path.join(HERE, "..", "SLIDES-LD8-hann04-corr_no-dx2-d8a10", "r_%04.1f" % r), 1) for r in r8]
    assert all(s == 3000 for _, s in st), st
    FL2 = np.array([KT * v / 2.0 for v, _ in st]); e2 = (FL2 - F_ref(r8)) / F_ref(r8)
    print("dx = 1 intermediate, all:"); [print("%4.1f | %.6e %+8.3f" % (r, f, 100 * x)) for r, f, x in zip(ri, FLi, ei)]
    print("dx = 2 intermediate (step 3000):"); [print("%4.1f | %.6e %+8.3f" % (r, f, 100 * x)) for r, f, x in zip(r8, FL2, e2)]
    lab = "$\\Delta x = 1$, intermediate $r$ (4.1 ... 4.9, 8.2 ... 9.8; current exe)"
    ax[0].plot(ri, FLi, "o", mfc="none", mec="#2b2a27", ms=5, label="Ludwig, " + lab)
    ax[1].plot(ri, 100 * ei, "o", mfc="none", mec="#2b2a27", ms=5, label=lab)
    lab = "$\\Delta x = 2$, intermediate $r$ (8.2 ... 9.8; current exe)"
    ax[0].plot(r8, FL2, "^", mfc="none", mec="#eb6834", ms=5, label="Ludwig, " + lab)
    ax[1].plot(r8, 100 * e2, "^", mfc="none", mec="#eb6834", ms=5, label=lab)
    for a in ax[:2]:
        a.legend(fontsize=8)
    fig.savefig(os.path.join(OUT, "FUERZA_LD08.0_Hann-8_corr-no_dx2-1-half-int.png"), dpi=300, bbox_inches="tight")
    # 2026-11-11 (user's request): NEW figure, same series as the "-half-int" one, without the order panel and with the
    # ABSOLUTE error F_L - F_ref (signed, force units of the top panel) above the relative error.
    fig, ax = plt.subplots(3, 1, figsize=(8.5, 11.0), sharex=True, gridspec_kw=dict(hspace=0.08))
    ax[0].plot(rr, fc, "-", color="#8a8880", lw=1.0, label="reference: periodic DH, continuous Hann (width 8)")
    SER = [(rs, FL[2.0], "^", "#eb6834", True, "$\\Delta x = 2$ (Hann-4)"),
           (rs, FL[1.0], "o", "#2b2a27", True, "$\\Delta x = 1$ (Hann-8)"),
           (rh, FLh, "D", "#2a78d6", True, "$\\Delta x = 1/2$ (Hann-16, current exe)"),
           (ri, FLi, "o", "#2b2a27", False, "$\\Delta x = 1$, intermediate $r$ (4.1 ... 4.9, 8.2 ... 9.8; current exe)"),
           (r8, FL2, "^", "#eb6834", False, "$\\Delta x = 2$, intermediate $r$ (8.2 ... 9.8; current exe)")]
    for r_, F_, mk, c, full, lab in SER:
        fref = F_ref(r_); kw = dict(color=c, ms=5) if full else dict(mfc="none", mec=c, ms=5, ls="none")
        ls = "-" if full else ""
        ax[0].plot(r_, F_, mk, label="Ludwig, " + lab, **kw)
        ax[1].plot(r_, F_ - fref, mk + ls, lw=1.0, label=lab, **kw)
        ax[2].plot(r_, 100 * (F_ - fref) / fref, mk + ls, lw=1.0, label=lab, **kw)
    for a in ax:
        kernel_lines(a)
        a.legend(fontsize=8)
    for a in ax[1:]:
        a.axhline(0, color="#8a8880", lw=0.8)
    ax[0].set_ylabel("$F_x$ on charge 1 (lattice units)")
    ax[1].set_ylabel("$F_L - F_{ref}$ (same units)")
    ax[2].set_ylabel("$(F_L - F_{ref}) / F_{ref}$ (%)")
    ax[2].set_xlabel("separation $r$ between the charges")
    ax[0].set_title("Two charges $q = +1$, periodic box 48, $\\lambda_D = 8$, Hann of width 8, $k_BT = 10^{-5}$ (lattice units)", fontsize=10, loc="left")
    ax[0].set_xlim(0, 24)
    fig.savefig(os.path.join(OUT, "FUERZA_LD08.0_Hann-8_corr-no_abs-rel.png"), dpi=300, bbox_inches="tight")
    # 2026-11-11 (user's request): NEW figure, zoom of the "abs-rel" figure on 8 <= r <= 10 (same series and panels; only
    # the points inside the window are drawn, so each panel's y range fits them; reference on a fine grid).
    X0, X1 = 7.8, 10.2
    fig, ax = plt.subplots(3, 1, figsize=(8.5, 11.0), sharex=True, gridspec_kw=dict(hspace=0.08))
    rz = np.linspace(X0, X1, 481)
    ax[0].plot(rz, F_ref(rz), "-", color="#8a8880", lw=1.0, label="reference: periodic DH, continuous Hann (width 8)")
    # 2026-11-11 (user's request): one series per lattice (integer and intermediate r together, same marker), no dx = 1/2.
    ZSER = [(np.concatenate([rs, r8]), np.concatenate([FL[2.0], FL2]), "^", "#eb6834", "$\\Delta x = 2$ (Hann-4)"),
            (np.concatenate([rs, ri]), np.concatenate([FL[1.0], FLi]), "o", "#2b2a27", "$\\Delta x = 1$ (Hann-8)")]
    for r_, F_, mk, c, lab in ZSER:
        m = (r_ >= X0) & (r_ <= X1); o = np.argsort(r_[m]); r_, F_ = r_[m][o], F_[m][o]
        fref = F_ref(r_); kw = dict(color=c, ms=6, ls="-", lw=1.0)   # 2026-11-11 (user): points joined by lines (sorted in r)
        ax[0].plot(r_, F_, mk, label="Ludwig, " + lab, **kw)
        ax[1].plot(r_, F_ - fref, mk, label=lab, **kw)
        ax[2].plot(r_, 100 * (F_ - fref) / fref, mk, label=lab, **kw)
    for a in ax[1:]:
        a.axhline(0, color="#8a8880", lw=0.8)
    for a in ax:
        a.legend(fontsize=8, loc="upper right") if a is ax[0] else None   # inside the top panel (free corner)
    ax[0].set_ylabel("$F_x$ on charge 1 (lattice units)")
    ax[1].set_ylabel("$F_L - F_{ref}$ (same units)")
    ax[2].set_ylabel("$(F_L - F_{ref}) / F_{ref}$ (%)")
    ax[2].set_xlabel("separation $r$ between the charges")
    ax[0].set_title("Zoom $8 \\leq r \\leq 10$\nTwo charges $q = +1$, periodic box 48, $\\lambda_D = 8$, Hann of width 8, $k_BT = 10^{-5}$ (lattice units)",
                    fontsize=10, loc="left")
    ax[0].set_xlim(X0, X1); ax[2].set_xticks(np.arange(8.0, 10.01, 0.2))
    for a in ax:
        a.grid(True, axis="x", color="#e5e5e5", lw=0.6)
    fig.savefig(os.path.join(OUT, "FUERZA_LD08.0_Hann-8_corr-no_abs-rel_zoom8-10.png"), dpi=300, bbox_inches="tight")
    # 2026-11-11 (user's request): NEW figure, same as the 8-10 zoom but for 4 <= r <= 5, Hann-8 (dx = 1) only:
    # r = 4.0 and 5.0 (old exe 1c30df4b) and 4.1 ... 4.9 (exe 4f0c20e7), one series joined by lines.
    X0, X1 = 3.95, 5.05
    fig, ax = plt.subplots(3, 1, figsize=(8.5, 11.0), sharex=True, gridspec_kw=dict(hspace=0.08))
    rz = np.linspace(X0, X1, 441)
    ax[0].plot(rz, F_ref(rz), "-", color="#8a8880", lw=1.0, label="reference: periodic DH, continuous Hann (width 8)")
    r_, F_ = np.concatenate([rs, ri]), np.concatenate([FL[1.0], FLi])
    m = (r_ >= X0) & (r_ <= X1); o = np.argsort(r_[m]); r_, F_ = r_[m][o], F_[m][o]
    fref = F_ref(r_); kw = dict(color="#2a78d6", ms=6, ls="-", lw=1.0); lab = "$\\Delta x = 1$ (Hann-8)"   # 2026-11-11 (user): coloured points
    ax[0].plot(r_, F_, "o", label="Ludwig, " + lab, **kw)
    ax[1].plot(r_, F_ - fref, "o", label=lab, **kw)
    ax[2].plot(r_, 100 * (F_ - fref) / fref, "o", label=lab, **kw)
    for a in ax[1:]:
        a.axhline(0, color="#8a8880", lw=0.8)
    ax[0].legend(fontsize=8, loc="upper right")
    ax[0].set_ylabel("$F_x$ on charge 1 (lattice units)")
    ax[1].set_ylabel("$F_L - F_{ref}$ (same units)")
    ax[2].set_ylabel("$(F_L - F_{ref}) / F_{ref}$ (%)")
    ax[2].set_xlabel("separation $r$ between the charges")
    ax[0].set_title("Zoom $4 \\leq r \\leq 5$\nTwo charges $q = +1$, periodic box 48, $\\lambda_D = 8$, Hann of width 8, $k_BT = 10^{-5}$ (lattice units)",
                    fontsize=10, loc="left")
    ax[0].set_xlim(X0, X1); ax[2].set_xticks(np.round(np.arange(4.0, 5.001, 0.1), 1))
    for a in ax:
        a.grid(True, axis="x", color="#e5e5e5", lw=0.6)
    fig.savefig(os.path.join(OUT, "FUERZA_LD08.0_Hann-8_corr-no_abs-rel_zoom4-5.png"), dpi=300, bbox_inches="tight")
