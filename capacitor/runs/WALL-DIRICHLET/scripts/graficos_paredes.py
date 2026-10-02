#!/usr/bin/env python3
"""
Figuras de validacion de las paredes a potencial constante (PETSc, Asta et al.
2019). Lee las salidas de Ludwig de A-sin-iones/ y B-debye-huckel/ (con sus
subcarpetas previo/, que guardan el campo ANTES de corregir psi_electric_field)
y recalcula las pruebas numericas con verificar_carga_pared.py y
convergencia_campo.py. Escribe PNG en ../graficos/.

Textos de las figuras en ingles (para presentaciones).
Uso:  python3 graficos_paredes.py
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "graficos")
sys.path.insert(0, HERE)
import verificar_carga_pared as vcp          # noqa: E402
import convergencia_campo as conv            # noqa: E402

# Paleta de referencia (ranuras 1-3 en orden fijo) + neutros para referencias
C1, C2, C3 = "#2a78d6", "#eb6834", "#1baf7a"   # blue, orange, aqua
REF = "#55544f"                                # reference curves (ink)
GRID = "#e4e2dc"
plt.rcParams.update({
    "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
    "axes.edgecolor": "#8a887f", "axes.linewidth": 0.8,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False,
    "legend.frameon": True, "legend.facecolor": "white", "legend.framealpha": 1.0,
    "legend.edgecolor": "#8a887f", "legend.fancybox": False,
    "lines.linewidth": 2.0,
    "savefig.dpi": 160, "savefig.bbox": "tight",
})

NX, NY, NZ = 4, 4, 34
BETA = 1.0e5            # efield sale como -grad(psi)/beta
DPSI = 1.0e-6
L = 32.0                # distancia entre planos medios (z = 1.5 y 33.5)
Z = np.arange(1, NZ + 1)
FL = slice(1, NZ - 1)   # nodos fluidos z = 2..33


def load(path, ncol=1):
    a = np.loadtxt(path)
    return a.reshape(NX, NY, NZ, ncol) if ncol > 1 else a.reshape(NX, NY, NZ)


def wall_bands(ax):
    """Sombrea los nodos de pared y marca los planos medios."""
    for z0 in (1, NZ):
        ax.axvspan(z0 - 0.5, z0 + 0.5, color="#d9d6cc", lw=0, zorder=0)
    for zm in (1.5, NZ - 0.5):
        ax.axvline(zm, color="#8a887f", lw=0.8, ls=(0, (2, 2)), zorder=1)


def save(fig, name):
    os.makedirs(OUT, exist_ok=True)
    fig.savefig(os.path.join(OUT, name))
    plt.close(fig)
    print("  ->", os.path.join("graficos", name))


# --------------------------------------------------------------------------
def fig_test_A():
    d = os.path.join(ROOT, "A-sin-iones")
    psi = load(os.path.join(d, "psi-000000001.001-001"))[0, 0, :]
    psi_naive = load(os.path.join(d, "naive", "psi-000000001.001-001"))[0, 0, :]
    E = -load(os.path.join(d, "efield-000000001.001-001"), 3)[0, 0, :, 2] * BETA
    Eo = -load(os.path.join(d, "previo", "efield-000000001.001-001"), 3)[0, 0, :, 2] * BETA
    Eex = DPSI / L
    zz = np.linspace(1.5, NZ - 0.5, 200)
    exact = (zz - 1.5) / L                              # psi / dpsi
    zzn = np.linspace(1, NZ, 200)
    naive = (zzn - 1) / (NZ - 1)                        # wall at the solid node

    fig, ax = plt.subplot_mosaic([["a", "a"], ["b", "c"]], figsize=(9.5, 8.2),
                                  height_ratios=[1, 1])
    a = ax["a"]
    wall_bands(a)
    a.plot(zz, exact, color=REF, lw=1.2, ls="--", label="exact (wall at mid-plane)")
    a.plot(zzn, naive, color=C2, lw=1.4, ls=":", label="wall at solid node (theory)")
    a.plot(Z, psi_naive / DPSI, "^", ms=5, color=C2, mfc="white", mew=1.4,
           label="Ludwig, PETSc: matrix without the wall correction")
    a.plot(Z, psi / DPSI, "o", ms=5, color=C1,
           label="Ludwig, PETSc: matrix with the wall correction")
    a.set(xlabel="z (lattice units)", ylabel=r"$\psi / \Delta\psi$",
          title="(a) Potential, no ions")
    a.legend(loc="upper left")

    # (b) Zoom at the lower wall, up to z = 5. No legend: same series/colors as (a).
    a = ax["b"]
    wall_bands(a)
    zn = np.linspace(1, 5, 50)
    a.plot(zn, (zn - 1.5) / L, color=REF, lw=1.2, ls="--")
    a.plot(zn, (zn - 1) / (NZ - 1), color=C2, lw=1.4, ls=":")
    a.plot(Z[:5], psi_naive[:5] / DPSI, "^", ms=6, color=C2, mfc="white", mew=1.6)
    a.plot(Z[:5], psi[:5] / DPSI, "o", ms=6, color=C1)
    a.set(xlim=(0.5, 5), ylim=(-0.01, 0.14), xlabel="z (lattice units)",
          ylabel=r"$\psi / \Delta\psi$", title="(b) Zoom at the lower wall")

    a = ax["c"]
    wall_bands(a)
    a.axhline(1.0, color=REF, lw=1.2, ls="--", label="exact")
    a.plot(Z[FL], Eo[FL] / Eex, "s", ms=5, color=C2, mfc="white", mew=1.6,
           label="same $\\psi$, plain central difference")
    a.plot(Z[FL], E[FL] / Eex, "o", ms=4, color=C1,
           label="same $\\psi$, wall-consistent gradient")
    a.set(ylim=(0.7, 1.05), xlabel="z (lattice units)", ylabel=r"$E_z / E_{exact}$",
          title="(c) Electric field")
    a.legend(loc="center")
    fig.suptitle("Test A - plate capacitor without ions (walls at z = 1 and z = 34, "
                 r"$\psi_{low} = 0$, $\psi_{up} = 10^{-6}$)", y=1.0, fontsize=11)
    fig.tight_layout()
    save(fig, "A_capacitor_sin_iones.png")


# --------------------------------------------------------------------------
def fig_test_B():
    d = os.path.join(ROOT, "B-debye-huckel")
    step = 60000
    psi = load(os.path.join(d, "psi-%09d.001-001" % step))[0, 0, :]
    q = load(os.path.join(d, "qsi-%09d.001-001" % step), 2)[0, 0, :, :]
    E = -load(os.path.join(d, "efield-%09d.001-001" % step), 3)[0, 0, :, 2] * BETA
    Eo = -load(os.path.join(d, "previo", "efield-%09d.001-001" % step), 3)[0, 0, :, 2] * BETA
    kap, rhos = 0.25, 3.125e-3
    zz = np.linspace(1.5, NZ - 0.5, 400)
    zc = lambda z: z - 17.5
    cont = lambda z: 0.5 * DPSI + 0.5 * DPSI * np.sinh(kap * zc(z)) / np.sinh(kap * L / 2)
    dcont = lambda z: 0.5 * DPSI * kap * np.cosh(kap * zc(z)) / np.sinh(kap * L / 2)
    # DH discreta con la misma regla de pared
    n = NZ - 2
    disc = conv.disc_rule(n, kap, "pared") * DPSI      # misma ecuacion (L=n, h=1)

    fig, ax = plt.subplots(2, 2, figsize=(11.5, 7.6))
    a = ax[0, 0]
    wall_bands(a)
    a.plot(zz, cont(zz) / DPSI, color=REF, lw=1.2, ls="--", label="Debye-Hückel, continuum (eq. 17)")
    a.plot(Z[FL], disc / DPSI, color=C3, lw=1.6, label="Debye-Hückel, discrete (same lattice)")
    a.plot(Z, psi / DPSI, "o", ms=4, color=C1, label="Ludwig, PETSc")
    a.set(xlabel="z (lattice units)", ylabel=r"$\psi / \Delta\psi$",
          title=r"(a) Potential, $\lambda_D = 4$, steady state")
    a.legend(loc="upper left")

    a = ax[0, 1]
    wall_bands(a)
    e_disc = np.abs(psi[FL] - disc) / DPSI
    e_cont = np.abs(psi[FL] - cont(Z[FL])) / DPSI
    a.semilogy(Z[FL], e_cont, "s", ms=4, color=C2, mfc="white", mew=1.4,
               label="vs continuum (lattice error)")
    a.semilogy(Z[FL], np.maximum(e_disc, 1e-16), "o", ms=4, color=C1,
               label="vs discrete (implementation)")
    a.legend(loc="lower center")
    a.set(ylim=(1e-12, 1e-1), xlabel="z (lattice units)",
          ylabel=r"$|\psi - \psi_{DH}| / \Delta\psi$", title="(b) Deviation from Debye-Hückel")

    a = ax[1, 0]
    wall_bands(a)
    pb = 0.5 * DPSI            # potencial de bulk por simetria
    a.plot(zz, -(cont(zz) - pb) * 1e7, color=REF, lw=1.2, ls="--",
           label=r"Boltzmann: $\mp(\psi - \psi_b)$")
    a.plot(zz, (cont(zz) - pb) * 1e7, color=REF, lw=1.2, ls="--")
    a.plot(Z[FL], (q[FL, 0] / rhos - 1) * 1e7, "o", ms=4, color=C1, label=r"cation $\rho_+$")
    a.plot(Z[FL], (q[FL, 1] / rhos - 1) * 1e7, "^", ms=4, color=C2, label=r"anion $\rho_-$")
    a.set(xlabel="z (lattice units)", ylabel=r"$(\rho_\pm / \rho_s - 1) \times 10^7$",
          title="(c) Ion densities")
    a.legend(loc="upper center", ncol=1)

    a = ax[1, 1]
    wall_bands(a)
    a.axhline(1.0, color=REF, lw=1.2, ls="--", label="exact")
    a.plot(Z[FL], Eo[FL] / dcont(Z[FL]), "s", ms=5, color=C2, mfc="white", mew=1.6,
           label="same $\\psi$, plain central difference")
    a.plot(Z[FL], E[FL] / dcont(Z[FL]), "o", ms=4, color=C1,
           label="same $\\psi$, wall-consistent gradient")
    a.set(ylim=(0.68, 1.06), xlabel="z (lattice units)",
          ylabel=r"$E_z / E_{DH,\,continuum}$", title="(d) Electric field")
    a.annotate("%.3f" % (Eo[1] / dcont(2.0)), (2, Eo[1] / dcont(2.0)), xytext=(6, 0.72),
               color="#3b3a36", arrowprops=dict(arrowstyle="-", color="#8a887f", lw=0.8))
    a.annotate("%.5f" % (E[1] / dcont(2.0)), (2, E[1] / dcont(2.0)), xytext=(6, 0.93),
               color="#3b3a36", arrowprops=dict(arrowstyle="-", color="#8a887f", lw=0.8))
    a.legend(loc="lower center")
    fig.suptitle(r"Test B - capacitor with 1:1 electrolyte ($\lambda_D = 4$, $L = 32$, "
                 r"$\Delta\psi = 10^{-6}$ kT/e)", y=1.0, fontsize=11)
    fig.tight_layout()
    save(fig, "B_debye_huckel.png")


# --------------------------------------------------------------------------
def fig_charge_rules():
    nf, ns = 20, 3
    phi, chi, nz, _ = vcp.solve(nf, ns, 0.0, 1.0, lambda z: 0.0, "modified")
    exact = 1.0 / nf
    rules = [("eq15", "Eq. 15 read literally"),
             ("eq12", "Eq. 12, unmodified"),
             ("sym", "factor 2 on every wall link\n(symmetric rule, implemented)")]
    vals, neut = [], []
    rng = np.random.default_rng(1)
    rz = rng.normal(size=nf) * 1e-2
    phi2, chi2, nz2, _ = vcp.solve(nf, ns, 0.1, 0.2, lambda z: rz[int(z - 0.5)], "modified")
    for key, _ in rules:
        lo, hi = vcp.charges(phi, chi, nz, 2, ns, nf, key)
        vals.append(hi / exact)
        lo2, hi2 = vcp.charges(phi2, chi2, nz2, 2, ns, nf, key)
        neut.append(abs(lo2 + hi2 + rz.sum()))

    fig, ax = plt.subplots(1, 2, figsize=(11, 3.2))
    y = np.arange(len(rules))[::-1]
    cols = ["#b9b6ab", "#b9b6ab", C1]
    a = ax[0]
    a.barh(y, vals, height=0.5, color=cols)
    a.axvline(1.0, color=REF, lw=1.2, ls="--")
    for yi, v in zip(y, vals):
        a.text(max(v, 0) + 0.02, yi, "%.4f" % (abs(v) if abs(v) < 5e-5 else v), va="center",
               color="#3b3a36")
    a.set(yticks=y, yticklabels=[r for _, r in rules], xlim=(0, 1.25),
          xlabel=r"$Q_{electrode} / Q_{Gauss}$", title="(a) Induced charge, no ions")
    a.grid(axis="y", visible=False)
    a = ax[1]
    a.barh(y, np.maximum(neut, 1e-16), height=0.5, color=cols)
    a.set_xscale("log")
    for yi, v in zip(y, neut):
        a.text(max(v, 1e-16) * 1.6, yi, "%.1e" % v, va="center", color="#3b3a36")
    a.set(yticks=y, yticklabels=["", "", ""], xlim=(1e-16, 1),
          xlabel=r"$|Q_{low} + Q_{up} + Q_{ions}|$", title="(b) Charge neutrality, random ionic charge")
    a.grid(axis="y", visible=False)
    fig.suptitle("Electrode charge from Eq. 16 of Asta et al. - which Laplacian at the electrode nodes "
                 "(D3Q27, plate capacitor)", y=1.04, fontsize=11)
    save(fig, "carga_inducida_reglas.png")


# --------------------------------------------------------------------------
def fig_convergence():
    hs, field, prof = conv.convergence()
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    a = ax[0]
    a.loglog(hs, prof["pared"], "o-", color=C1, ms=6, label="wall at mid-plane (implemented)")
    a.loglog(hs, prof["ingenuo"], "s--", color=C2, ms=6, mfc="white", mew=1.6,
             label="wall at solid node (naive)")
    a.loglog(hs, prof["pared"][0] * (hs / hs[0]) ** 2, color=REF, lw=0.9, ls=":")
    a.loglog(hs, prof["ingenuo"][0] * (hs / hs[0]), color=REF, lw=0.9, ls=":")
    a.text(hs[2], prof["pared"][0] * (hs[2] / hs[0]) ** 2 * 0.4, r"$\propto h^2$", color="#3b3a36")
    a.text(hs[2], prof["ingenuo"][0] * (hs[2] / hs[0]) * 1.6, r"$\propto h$", color="#3b3a36")
    a.set(xlabel=r"lattice spacing $h$ ($\kappa L = 8$ fixed)", ylabel="max profile error / Δψ",
          title="(a) Potential profile")
    a.legend(loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=1)
    a = ax[1]
    lab = {"central (antes)": ("central difference (before)", C2, "s--"),
           "factor 2 (ahora)": ("factor 2, difference form (implemented)", C1, "o-"),
           "cuadratica": ("quadratic interpolation (rejected)", C3, "^-.")}
    for k, (name, c, st) in lab.items():
        kw = dict(mfc="white", mew=1.6) if "s" in st else {}
        a.loglog(hs, field[k], st, color=c, ms=6, label=name, **kw)
    a.set(xlabel=r"lattice spacing $h$", ylabel="relative error of E at the first fluid node",
          title="(b) Electric field next to the wall", ylim=(5e-8, 1))
    a.legend(loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=1)
    fig.suptitle("Convergence on the discrete Debye-Hückel solution (Ludwig reproduces it to 8e-9)",
                 y=1.02, fontsize=11)
    fig.tight_layout()
    save(fig, "convergencia.png")


if __name__ == "__main__":
    print("Figuras de validacion de paredes:")
    fig_test_A()
    fig_test_B()
    fig_charge_rules()
    fig_convergence()
