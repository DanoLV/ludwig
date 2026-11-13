#!/usr/bin/env python3
"""Potential and field at the centre of the displaced charge, all positions and the three lattices (2026-11-08, user's
request; same runs as grafico_C_convergencia_desplazada.py todos): dx = 2 (Hann-4), 1 (Hann-8), 1/2 (Hann-16).

Ludwig: psi and E_z (efield dump) of the last step on the vertical line x = y = x_charge (the charge is on a node in x, y
in the three lattices), linearly interpolated in z at z0 = delta. Trilinear interpolation at r0 reduces to this. The
interpolation adds its own error (user's choice, 2026-11-08: pointwise value at r0, interpolated).
Reference: linear Debye-Hueckel solution with the continuous Hann source (grafico_C_hoja_hann.py) at r0:
  phi_DH(r0) = HH.phi(0, delta, 0, 0, mu),   E_DH,z(r0) = -d phi_DH/dz (closed form, HH.dphi via grafico_C_desplazada.ez_an).
Not separated: the nonlinear Poisson-Boltzmann part (Ludwig is nonlinear, phi_DH linear).
Figure (one column per quantity): value of Ludwig over the analytic curve; (Ludwig - DH)/DH in % with sign;
observed order p = log2 |e(coarse)/e(fine)| at the same delta. delta = 0: E = 0, relative error of E undefined.
dx = 1/2 runs that are not finished are left out.
Usage: python3 grafico_C_carga_desplazada_puntual.py   -> graficos/C_carga_desplazada_puntual.png
"""
import os
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import grafico_C_formulas as F
import grafico_C_hoja_hann as HH
from grafico_C_refinamiento import load, line
from grafico_C_refinamiento_zoom import ez_file
from grafico_C_desplazada import set_charge, ez_an, psi_files
from grafico_C_convergencia_desplazada import ALL, run_path, DPI

LATS = ((0.5, "$\\Delta x = 2$ (Hann-4)", "#eb6834", "^"), (1, "$\\Delta x = 1$ (Hann-8)", "#2b2a27", "o"),
        (2, "$\\Delta x = 1/2$ (Hann-16)", "#2a78d6", "s"))

def at_charge(path, r, dl):
    """(psi, E_z) of Ludwig at z0 = delta, linear interpolation in z along the line through the charge."""
    fs = psi_files(path)
    z, ps, _ = line(load(fs[-1], r), r, 0); ez = ez_file(path, fs[-1], r)
    return np.interp(dl, z, ps), np.interp(dl, z, ez)

if __name__ == "__main__":
    rows = []                                          # (delta, phi_DH, E_DH, psi_L and E_L for r = 0.5, 1, 2)
    for dl in ALL:
        set_charge(dl); mu = HH.mean_U()
        pa = float(HH.phi(0.0, np.array([dl]), 0, 0, mu)[0]); ea = float(ez_an(np.array([dl]), mu)[0])
        v = [dl, pa, ea]
        for r, *_ in LATS:
            path = run_path(dl, r)
            if not os.path.isdir(os.path.join(F.ROOT, path)) or not psi_files(path):
                v += [np.nan, np.nan]; continue
            v += list(at_charge(path, r, dl))
        rows.append(v)
        print("delta %5.2f | phi_DH %+.6e  psi_L %s | E_DH %+.6e  E_L %s" % (dl, pa, " ".join("%+.6e" % x for x in v[3::2]),
              ea, " ".join("%+.6e" % x for x in v[4::2])))
    a = np.array(rows)
    F.style()
    fig, ax = plt.subplots(3, 2, figsize=(16.0, 12.0), sharex=True, gridspec_kw=dict(hspace=0.08, wspace=0.22))
    for j, (ka, name, unit) in enumerate(((1, "\\psi", "$e\\psi/k_BT$"), (2, "E", "$e E_z/k_BT$, $z$ component"))):
        ref = a[:, ka]; ok = ref != 0                  # E at delta = 0: relative error undefined
        ax[0, j].plot(a[:, 0], ref, "-", color="#8a8880", lw=2.5, alpha=0.6, label="DH analytic at $r_0$")
        e = {}
        for i, (r, lab, c, mk) in enumerate(LATS):
            val = a[:, 3 + 2 * i + (ka - 1)]; m = ~np.isnan(val)
            ax[0, j].plot(a[m, 0], val[m], mk, color=c, ms=4, label="Ludwig, " + lab)
            e[r] = np.where(ok, (val - ref) / np.where(ok, ref, 1.0), np.nan)
            m = ~np.isnan(e[r])
            ax[1, j].plot(a[m, 0], 100 * e[r][m], mk + "-", color=c, lw=1.0, ms=4, label=lab)
        for (rc, rf), lab, c, mk in (((0.5, 1), "$p(\\Delta x = 2 \\to 1)$", "#eb6834", "^-"),
                                     ((1, 2), "$p(\\Delta x = 1 \\to 1/2)$", "#2a78d6", "s")):
            p = np.log2(np.abs(e[rc] / e[rf])); m = ~np.isnan(p)
            ax[2, j].plot(a[m, 0], p[m], mk, color=c, lw=1.0, ms=4, label=lab)
        for q in (1, 2):
            ax[2, j].axhline(q, color="#8a8880", lw=0.8, ls=":")
        ax[1, j].axhline(0, color="#8a8880", lw=0.8)
        ax[0, j].set_ylabel("$%s(r_0)$ (reduced, %s)" % (name, unit))
        ax[1, j].set_ylabel("$(%s_L - %s_{DH}) / %s_{DH}$ at $r_0$ (%%)" % (name, name, name))
        ax[2, j].set_ylabel("observed order $p$")
        ax[2, j].set_xlabel("$\\delta$ (charge at $z_0 = \\delta$ from the mid-plane)")
        ax[0, j].set_title("$%s$ at the centre of the charge (Ludwig: linear interpolation in $z$)" % name, fontsize=10, loc="left")
        for k in range(3):
            ax[k, j].legend(fontsize=8)
            for zn in np.arange(0.5, 12.6, 1.0):       # dx = 1 nodes (delta = k + 1/2)
                ax[k, j].axvline(zn, color="#c9c6bd", lw=0.6, ls="--", zorder=0)
            for zn in np.arange(1.0, 12.1, 2.0):       # dx = 2 nodes (delta = 2k + 1)
                ax[k, j].axvline(zn, color="#eb6834", lw=0.9, ls=":", zorder=0)
    fig.text(0.5, 0.035, "vertical lines: nodes of $\\Delta x = 1$ (grey dashed) and $\\Delta x = 2$ (orange dotted)",
             ha="center", fontsize=9)
    fig.savefig(os.path.join(F.ROOT, "graficos", "C_carga_desplazada_puntual.png"), dpi=DPI, bbox_inches="tight")
