#!/usr/bin/env python3
"""Convergence of the displaced charge between dx = 1 (Hann-8, C-electrolito-hann8-LD6-desplazada, crear_desplazada.py)
and dx = 1/2 (Hann-16, C-electrolito-hann16-LD12-ref2-desplazada, crear_desplazada_ref.py): same physical system and
the same physical kernel (continuous Hann of width 8), delta = 1.0, 4.7, 8.3, 12.0 (z0 = delta from the mid-plane).

Errors (all against the analytic linear Debye-Hueckel solution with the continuous Hann charge, as in
grafico_C_desplazada.py / grafico_C_fuerza_desplazada.py):
  psi, E_z  on the z axis through the charge: d = Ludwig - analytic at each node. The z nodes of the two lattices do
            not coincide (dx = 1: z = +-0.5, +-1.5, ...; dx = 1/2: z = +-0.25, +-0.75, ...), so the comparison is by
            bands (user's choice, 2026-11-06): RMS of d on bands of width 2 with ends at z + 16 = 2j (2 nodes with
            dx = 1, 4 nodes with dx = 1/2; grafico_C_refinamiento_anidado.band_rms).
  F         F_L - F_kernel, F_L = Esub_z * R (Ludwig), F_kernel = int W E_DH (continuous limit, the target force).
Observed order between the two lattices: p = log2(e(dx = 1) / e(dx = 1/2)); with two lattices it is a single number
per band (or per delta), not a fit. It assumes e ~ C dx^p already holds at dx = 1; it does not test that.
Guide lines of slope 1 and 2 in the log-log panel are visual references only (through the mean of the dx = 1 points).
Not separated: the nonlinear Poisson-Boltzmann part (Ludwig is nonlinear, the analytic solution linear).
Steady state of the dx = 1 and 1/2 batches was checked by grafico_C_desplazada.py (placa, ref2); not repeated here.
dx = 2 (2026-11-06, user's request; force figure only): Hann-4, C-electrolito-hann4-LD3-ref0.5-desplazada
(crear_desplazada_ref.py r=0.5), N = 2000, psi every 250 steps. Its steady state is checked here with the same criterion
(max|psi(2000) - psi(1750)| and max|E(2000) - E(1750)| < 1 % of the smallest band RMS of d_psi, d_E; 1 node per band).
In the force figure p is given for each pair of lattices, 2 -> 1 and 1 -> 1/2.
All deltas (2026-11-06, user's request: understand why the observed order fluctuates): every delta run with Hann-8
(cell 0 ... 0.9 and approach to the plate 1.0 ... 12.0) also run with Hann-4 (dx = 2). Third figure,
graficos/C_convergencia_desplazada_fuerza_todos.png: (F_L - F_kernel)/F_kernel vs delta for dx = 2, 1 (and 1/2 at its
4 deltas), and p(2 -> 1) [p(1 -> 1/2)] vs delta. Lattice nodes are marked on the delta axis (dx = 1: delta = k + 1/2;
dx = 2: delta = 2k + 1). dx = 1 cell runs with delta <= 0.5 use the -cont run (10000 steps), as mode celda of
grafico_C_desplazada.py. Steady state of every dx = 2 run checked here as above. Usage: python3
grafico_C_convergencia_desplazada.py [todos]   (todos: only the third figure)
Output: graficos/C_convergencia_desplazada_campo.png, graficos/C_convergencia_desplazada_fuerza.png (dpi 300).
Usage: python3 grafico_C_convergencia_desplazada.py
"""
import os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import grafico_C_formulas as F
import grafico_C_hoja_hann as HH
from grafico_C_refinamiento import load, line
from grafico_C_refinamiento_anidado import band_rms
from grafico_C_refinamiento_zoom import ez_file, KT
from grafico_C_desplazada import set_charge, ez_an, psi_files
from grafico_C_fuerza_desplazada import F_kernel

DELTAS = (1.0, 4.7, 8.3, 12.0)
BATS = {0.5: "C-electrolito-hann4-LD3-ref0.5-desplazada", 1: "C-electrolito-hann8-LD6-desplazada",
        2: "C-electrolito-hann16-LD12-ref2-desplazada"}
COL = {d: matplotlib.colormaps["viridis"](0.9 * (d - min(DELTAS)) / (max(DELTAS) - min(DELTAS))) for d in DELTAS}
DPI = 300

ALL = tuple(sorted({round(0.1 * k, 1) for k in range(121)} | {0.25}))   # 0 ... 12 every 0.1, plus 0.25 (2026-11-06)

def run_path(dl, r):
    z = "z_%05.2f" % (1.5 + r * (16 + dl))
    cont = os.path.join(BATS[1] + "-cont", z)
    return cont if r == 1 and os.path.isdir(os.path.join(F.ROOT, cont)) else os.path.join(BATS[r], z)

def steady_dx2(path, fs, mu):
    """Steady state of a dx = 2 run: (D_psi, threshold, D_E, threshold), 1 % of the smallest band RMS of d."""
    r = 0.5; n = 16
    z, ps, _ = line(load(fs[-1], r), r, 0); ez = ez_file(path, fs[-1], r)
    _, Rp = band_rms(z, ps - HH.phi(0.0, z, 0, 0, mu), r, 2.0, 2.0); _, Re = band_rms(z, ez - ez_an(z, mu), r, 2.0, 2.0)
    fe = lambda f: np.loadtxt(os.path.join(F.ROOT, path, "efield-" + os.path.basename(f)[4:])).reshape(n, n, n + 2, 3) / KT * r
    return (np.abs(load(fs[-1], r) - load(fs[-2], r)).max(), 0.01 * Rp.min(),
            np.abs(fe(fs[-1]) - fe(fs[-2])).max(), 0.01 * Re.min())

def figure_todos():
    rows = []; absd = []                               # absd: (delta, F_kernel, F_L - F_kernel for dx = 2, 1, 1/2)
    for dl in ALL:
        set_charge(dl); mu = HH.mean_U(); fk = F_kernel(dl, mu); e = {}; d = {}
        for r in (0.5, 1, 2):
            path = run_path(dl, r)
            if not os.path.isdir(os.path.join(F.ROOT, path)): continue
            fs = psi_files(path)
            if not fs and r == 2:                       # 2026-11-07: partial analysis while the dx = 1/2 batch runs
                print("delta %5.2f dx = 1/2: not finished, left out (%s)" % (dl, path)); continue
            assert fs, "not finished: " + path
            if r == 0.5:
                s = steady_dx2(path, fs, mu)
                print("delta %5.2f dx = 2: D_psi %.1e vs %.1e %s | D_E %.1e vs %.1e %s" % (dl, s[0], s[1],
                      "OK" if s[0] < s[1] else "NOT steady", s[2], s[3], "OK" if s[2] < s[3] else "NOT steady"))
            d[r] = f_ludwig(path, r)[1] - fk
            e[r] = d[r] / fk if fk != 0 else np.nan
        rows.append((dl, e.get(0.5, np.nan), e.get(1, np.nan), e.get(2, np.nan)))
        absd.append((dl, fk, d.get(0.5, np.nan), d.get(1, np.nan), d.get(2, np.nan)))
        print("delta %5.2f  rel. F_L - F_kernel (%%): dx = 2 %+.4f | dx = 1 %+.4f | dx = 1/2 %s" % (dl, 100 * rows[-1][1],
              100 * rows[-1][2], "%+.4f" % (100 * rows[-1][3]) if not np.isnan(rows[-1][3]) else "-"))
    r = np.array(rows)[1:]                              # delta = 0: F = 0, relative error undefined
    p0 = np.log2(np.abs(r[:, 1] / r[:, 2])); p1 = np.log2(np.abs(r[:, 2] / r[:, 3]))
    F.style()
    fig, ax = plt.subplots(4, 1, figsize=(12.0, 17.0), sharex=True, gridspec_kw=dict(hspace=0.08))
    ad = np.array(absd)                                # 2026-11-06 (user's request): F_L over F_kernel, and F_L - F_kernel
    ax[0].plot(ad[:, 0], ad[:, 1], "-", color="#8a8880", lw=2.5, alpha=0.6, label="$F_{kernel} = \\int W E_{DH}$")
    ax[1].axhline(0, color="#8a8880", lw=0.8)
    for k, lab, c, mk in ((1, "$\\Delta x = 2$ (Hann-4)", "#eb6834", "^"), (2, "$\\Delta x = 1$ (Hann-8)", "#2b2a27", "o"),
                          (3, "$\\Delta x = 1/2$ (Hann-16)", "#2a78d6", "s")):
        m = ~np.isnan(ad[:, k + 1])
        ax[0].plot(ad[m, 0], ad[m, 1] + ad[m, k + 1], mk, color=c, ms=5, label="$F_L$, " + lab)
        m = ~np.isnan(r[:, k])
        ax[1].plot(r[m, 0], 100 * r[m, k], mk + "-", color=c, lw=1.0, ms=5, label=lab)   # signed, delta = 0 left out
        ax[2].semilogy(r[m, 0], 100 * np.abs(r[m, k]), mk + "-", color=c, lw=1.0, ms=5, label=lab)
    ax[3].plot(r[:, 0], p0, "^-", color="#eb6834", lw=1.0, ms=5, label="$p(\\Delta x = 2 \\to 1)$")
    m = ~np.isnan(p1)
    ax[3].plot(r[m, 0], p1[m], "s", color="#2a78d6", ms=6, label="$p(\\Delta x = 1 \\to 1/2)$")
    for q in (1, 2):
        ax[3].axhline(q, color="#8a8880", lw=0.8, ls=":")
    for a in ax:
        for zn in np.arange(0.5, 12.6, 1.0):            # dx = 1 nodes (delta = k + 1/2)
            a.axvline(zn, color="#c9c6bd", lw=0.6, ls="--", zorder=0)
        for zn in np.arange(1.0, 12.1, 2.0):            # dx = 2 nodes (delta = 2k + 1)
            a.axvline(zn, color="#eb6834", lw=0.9, ls=":", zorder=0)
    ax[0].set_ylabel("$F_z/q$ (reduced, $e E/k_BT$)"); ax[1].set_ylabel("$(F_L - F_{kernel}) / F_{kernel}$ (%)")
    ax[2].set_ylabel("$|F_L - F_{kernel}| / F_{kernel}$ (%)"); ax[3].set_ylabel("observed order $p$")
    ax[3].set_xlabel("$\\delta$ (charge at $z_0 = \\delta$ from the mid-plane); vertical lines: nodes of $\\Delta x = 1$ "
                     "(grey dashed) and $\\Delta x = 2$ (orange dotted)")
    ax[0].legend(fontsize=9, loc="upper left"); ax[1].legend(fontsize=9, loc="upper left")
    ax[2].legend(fontsize=9); ax[3].legend(fontsize=9)
    ax[0].set_title("Force on the displaced charge, all positions: error vs $F_{kernel} = \\int W E_{DH}$ and observed order",
                    fontsize=10, loc="left")
    fig.savefig(os.path.join(F.ROOT, "graficos", "C_convergencia_desplazada_fuerza_todos.png"), dpi=DPI, bbox_inches="tight")

def f_ludwig(path, r):
    """(step, Esub_z * r) of the last line of particle_Esub.csv (as grafico_C_fuerza_desplazada.F_ludwig)."""
    last = [l for l in open(os.path.join(F.ROOT, path, "proceced_data", "particle_Esub.csv")) if not l.startswith("#")][-1]
    v = last.strip().replace(",", ".").split(";")
    return int(v[0]), float(v[5]) * r

if __name__ == "__main__":
    if "todos" in sys.argv[1:]:
        figure_todos(); sys.exit()
    band = {}; force = {}
    for dl in DELTAS:
        set_charge(dl); mu = HH.mean_U()
        fk = F_kernel(dl, mu)
        for r in (0.5, 1, 2):
            path = run_path(dl, r); fs = psi_files(path)
            assert fs, "not finished: " + path
            z, ps, _ = line(load(fs[-1], r), r, 0)
            ez = ez_file(path, fs[-1], r)
            dp = ps - HH.phi(0.0, z, 0, 0, mu); de = ez - ez_an(z, mu)
            zc, Rp = band_rms(z, dp, r, 2.0, 2.0); _, Re = band_rms(z, de, r, 2.0, 2.0)
            band[dl, r] = (Rp, Re)
            if r == 0.5:                                   # steady state of the dx = 2 runs (not checked elsewhere)
                n = int(round(32 * r)); fe = lambda f: np.loadtxt(os.path.join(F.ROOT, path, "efield-" + os.path.basename(f)[4:])).reshape(n, n, n + 2, 3) / KT * r
                dps = np.abs(load(fs[-1], r) - load(fs[-2], r)).max(); des = np.abs(fe(fs[-1]) - fe(fs[-2])).max()
                print("delta %5.2f  dx = 2  dumps %s -> %s: D_psi %.1e vs %.1e %s | D_E %.1e vs %.1e %s"
                      % (dl, os.path.basename(fs[-2])[4:13], os.path.basename(fs[-1])[4:13], dps, 0.01 * Rp.min(),
                         "OK" if dps < 0.01 * Rp.min() else "NOT steady", des, 0.01 * Re.min(), "OK" if des < 0.01 * Re.min() else "NOT steady"))
            step, fl = f_ludwig(path, r)
            force[dl, r] = (fl, fk)
            print("delta %5.2f  dx = %-3g  step %5d  F_L %+.6e  F_kernel %+.6e  F_L - F_kernel %+.3e (%+.3f %%)"
                  % (dl, 1 / r, step, fl, fk, fl - fk, 100 * (fl - fk) / fk))
    p = {dl: [np.log2(band[dl, 1][k] / band[dl, 2][k]) for k in (0, 1)] for dl in DELTAS}
    ef = {(dl, r): abs(force[dl, r][0] - force[dl, r][1]) / abs(force[dl, r][1]) for dl in DELTAS for r in (0.5, 1, 2)}
    pf = {dl: np.log2(ef[dl, 1] / ef[dl, 2]) for dl in DELTAS}         # 1 -> 1/2
    pf0 = {dl: np.log2(ef[dl, 0.5] / ef[dl, 1]) for dl in DELTAS}      # 2 -> 1
    print("band centre      " + " ".join("%7.1f" % c for c in zc))
    for k, name in ((0, "p psi"), (1, "p E  ")):
        for dl in DELTAS:
            print("%s  %5.2f  " % (name, dl) + " ".join("%7.2f" % v for v in p[dl][k]))
    print("p F (2 -> 1):   " + "  ".join("delta %.1f: %.2f" % (dl, pf0[dl]) for dl in DELTAS))
    print("p F (1 -> 1/2): " + "  ".join("delta %.1f: %.2f" % (dl, pf[dl]) for dl in DELTAS))

    # ---- fields: band RMS at both lattices and band order, psi (left) and E_z (right)
    F.style()
    fig, ax = plt.subplots(len(DELTAS) + 1, 2, figsize=(12.0, 16.0), sharex=True,
                           gridspec_kw=dict(hspace=0.1, wspace=0.25))
    for i, dl in enumerate(DELTAS):
        for k in (0, 1):
            a = ax[i, k]
            a.axvspan(dl - 4, dl + 4, color="#f1efe8", zorder=0)
            a.semilogy(zc, band[dl, 1][k], "o--", color=COL[dl], mfc="white", lw=1.0, label="$\\Delta x = 1$ (Hann-8)")
            a.semilogy(zc, band[dl, 2][k], "o-", color=COL[dl], lw=1.2, label="$\\Delta x = 1/2$ (Hann-16)")
            a.text(0.01, 0.95, "$\\delta = %.1f$" % dl, transform=a.transAxes, va="top", fontsize=10)
        ax[i, 0].set_ylabel("band RMS of $d_\\psi$"); ax[i, 1].set_ylabel("band RMS of $d_E$")
    ax[0, 0].legend(fontsize=9, loc="upper right")
    for k in (0, 1):
        a = ax[-1, k]
        for dl in DELTAS:
            a.plot(zc, p[dl][k], "o-", color=COL[dl], lw=1.0, label="$\\delta = %.1f$" % dl)
        for q in (1, 2):
            a.axhline(q, color="#8a8880", lw=0.8, ls=":")
        a.set_xlabel("z (physical units, from the mid-plane; band centres)")
    ax[-1, 0].set_ylabel("order $p$ = log$_2$(RMS$_{\\Delta x=1}$ / RMS$_{\\Delta x=1/2}$), $\\psi$")
    ax[-1, 1].set_ylabel("order $p$, $E_z$")
    ax[-1, 0].legend(fontsize=9, ncol=2)
    for a in ax.ravel():
        a.set_xlim(-16.5, 16.5)
    ax[0, 0].set_title("Displaced charge: convergence $\\Delta x = 1 \\to 1/2$ by bands of width 2 (shaded: kernel support; "
                       "dotted: p = 1, 2)", fontsize=10, loc="left")
    fig.savefig(os.path.join(F.ROOT, "graficos", "C_convergencia_desplazada_campo.png"), dpi=DPI, bbox_inches="tight")

    # ---- force: F_L vs dx with the continuous-limit target, log-log error, order vs delta
    fig = plt.figure(figsize=(12.0, 11.0))
    gs = fig.add_gridspec(3, len(DELTAS), height_ratios=[1, 1.3, 1], hspace=0.35, wspace=0.45)
    dxs = np.array([2.0, 1.0, 0.5]); RS = (0.5, 1, 2)
    for i, dl in enumerate(DELTAS):
        a = fig.add_subplot(gs[0, i]); fk = force[dl, 1][1]
        a.plot(dxs, [force[dl, r][0] for r in RS], "o-", color=COL[dl], label="Ludwig $F_L$")
        a.axhline(fk, color="#2a78d6", ls="--", lw=1.0, label="$F_{kernel}$ (continuous limit)")
        a.plot([0], [fk], "s", color="#2a78d6")
        a.set_xlim(-0.12, 2.15); a.set_xticks([0, 0.5, 1, 2]); a.set_xlabel("$\\Delta x$")
        a.set_title("$\\delta = %.1f$" % dl, fontsize=10)
        a.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        if i == 0:
            a.set_ylabel("$F_z/q$ (reduced)"); a.legend(fontsize=7, loc="upper left")
    a = fig.add_subplot(gs[1, :])
    for dl in DELTAS:
        a.loglog(dxs, [ef[dl, r] for r in RS], "o-", color=COL[dl],
                 label="$\\delta = %.1f$ (p: %.2f, %.2f)" % (dl, pf0[dl], pf[dl]))
    e1 = np.exp(np.mean([np.log(ef[dl, 1]) for dl in DELTAS]))
    xg = np.array([0.45, 2.2])
    a.loglog(xg, e1 * xg, ":", color="#8a8880", lw=1.0, label="slope 1 (guide)")
    a.loglog(xg, e1 * xg ** 2, "--", color="#8a8880", lw=1.0, label="slope 2 (guide)")
    a.set_xlabel("$\\Delta x$"); a.set_ylabel("$|F_L - F_{kernel}| / |F_{kernel}|$")
    a.set_xticks([0.5, 1.0, 2.0]); a.set_xticklabels(["1/2", "1", "2"]); a.legend(fontsize=8, ncol=2)
    a = fig.add_subplot(gs[2, :])
    a.plot(DELTAS, [pf0[dl] for dl in DELTAS], "s--", color="#8a8880", label="$\\Delta x = 2 \\to 1$")
    a.plot(DELTAS, [pf[dl] for dl in DELTAS], "o-", color="#2b2a27", label="$\\Delta x = 1 \\to 1/2$")
    for q in (1, 2):
        a.axhline(q, color="#8a8880", lw=0.8, ls=":")
    a.set_xlabel("$\\delta$ (charge at $z_0 = \\delta$ from the mid-plane)")
    a.set_ylabel("order $p$ of $F_L - F_{kernel}$"); a.legend(fontsize=8)
    fig.suptitle("Force on the displaced charge: $\\Delta x = 2$ (Hann-4) $\\to$ 1 (Hann-8) $\\to$ 1/2 (Hann-16), "
                 "target $F_{kernel} = \\int W E_{DH}$", fontsize=11, x=0.08, ha="left", y=0.93)
    fig.savefig(os.path.join(F.ROOT, "graficos", "C_convergencia_desplazada_fuerza.png"), dpi=DPI, bbox_inches="tight")
