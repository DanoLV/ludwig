#!/usr/bin/env python3
"""Effect of the lateral images as the screening decreases: two-potential electrolyte capacitor (plates at psi = 0.1 and 0.3),
point charge at d = 8 from the lower plate, lambda_D = 6, 8, 10, 12, 14.  Ludwig potential vs the closed-system formula

  phi = phi0 + [kL/(2 tanh(kL/2))] <phi0> [1 - cosh(kz)/cosh(kL/2)],   phi0 = phi_CAP + phi_q^p,   k = 1/lambda_D
Dashed: phi0 of the single cell (infinite plates). Solid: phi0 summed over lateral images (10 per side).
Lines rho = 6 and 10 (outside the kernel support).   Usage: python3 grafico_C_dospot_vs_LD.py
"""
import os, glob
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import grafico_C_formulas as F

L, A, FOUR, D = F.L, F.A, F.FOUR, F.D
P1, P2 = 0.1, 0.3
LDS = (6, 8, 10, 12, 14)
BAT = {6: "C-electrolito-dospotenciales-hann8-LD6"}
BAT.update({ld: "C-electrolito-dospotenciales-hann8-LD%d" % ld for ld in LDS if ld != 6})

def mean_phi0(kap):
    n = np.arange(0, 20000); m = 2 * n + 1
    mq = 4 * FOUR / (np.pi * A * L) * np.sum((-1.0) ** n / m * np.cos(m * np.pi * D / L) / ((m * np.pi / L) ** 2 + kap ** 2))
    return mq + (P1 + P2) * np.tanh(kap * L / 2) / (kap * L)

def closed(z, kap):
    return (kap * L) / (2 * np.tanh(kap * L / 2)) * mean_phi0(kap) * (1 - np.cosh(kap * z) / np.cosh(kap * L / 2))

def ready(bat):
    p = glob.glob(os.path.join(F.ROOT, bat, "z_09.50", "psi-0*"))
    return len(p) > 1 and "finished normally" in open(os.path.join(F.ROOT, bat, "z_09.50", "run.log"), errors="ignore").read()

if __name__ == "__main__":
    plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.color": "#e4e2dc",
                         "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False})
    fig, axs = plt.subplots(2, 3, figsize=(15, 8.4)); axs = axs.ravel()
    cols = {6: "#2a78d6", 10: "#eb6834"}; zz = np.linspace(-L / 2 + 0.2, L / 2 - 0.2, 160)
    res = {}
    print(" lD   phibar   rho | max|Ludwig - formula| no images | with images  | max(Ludwig) | error ratio")
    for k, ld in enumerate(LDS):
        a = axs[k]; kap = 1.0 / ld
        if not ready(BAT[ld]):
            a.set_title("$\\lambda_D=%d$: simulation not finished" % ld); continue
        for dx, c in cols.items():
            z, ps = F.ludwig(os.path.join(F.ROOT, BAT[ld]), dx)
            f0 = lambda zv: F.phi_cap_elec(zv, P1, P2, kap) + F.phi_q_AB(float(dx), zv, D, kap)
            rep = lambda zv: F.replicas_AB(float(dx), zv, kap)
            a.plot(zz, f0(zz) + closed(zz, kap), "--", color=c, lw=1.3, dashes=(5, 3), zorder=2)
            a.plot(zz, f0(zz) + rep(zz) + closed(zz, kap), "-", color=c, lw=1.4, alpha=0.55, zorder=2)
            a.plot(z, ps, "o", ms=4.2, color=c, mec="white", mew=0.5, zorder=4, label="Ludwig, $\\rho=%d$" % dx)
            e0 = np.abs(ps - f0(z) - closed(z, kap)).max(); e1 = np.abs(ps - f0(z) - rep(z) - closed(z, kap)).max()
            res[(ld, dx)] = (e0, e1)
            print("%3d  %.4f   %3d |     %.4f              |   %.4f     |   %.4f    |  %.1f"
                  % (ld, (kap * L) / (2 * np.tanh(kap * L / 2)) * mean_phi0(kap), dx, e0, e1, ps.max(), e0 / e1))
        a.set_title("$\\lambda_D=%d$" % ld); a.set_xlabel("z (lattice units)")
        if k % 3 == 0: a.set_ylabel(r"reduced potential $\psi$")
    axs[0].plot([], [], "--", color="#55544f", lw=1.3, dashes=(5, 3), label="formula, infinite plates")
    axs[0].plot([], [], "-", color="#55544f", lw=1.4, alpha=0.55, label="formula + lateral images")
    axs[0].legend(loc="upper left", fontsize=8.5)
    b = axs[5]
    if res:
        lds = sorted({k[0] for k in res})
        for dx, c in cols.items():
            xs = [l for l in lds if (l, dx) in res]
            b.semilogy(xs, [res[(l, dx)][0] for l in xs], "--o", color=c, lw=1.3, dashes=(5, 3), ms=5, label="$\\rho=%d$, infinite plates" % dx)
            b.semilogy(xs, [res[(l, dx)][1] for l in xs], "-s", color=c, lw=1.6, alpha=0.7, ms=5, label="$\\rho=%d$, with images" % dx)
        b.set_xlabel("$\\lambda_D$ (lattice units)"); b.set_ylabel(r"max $|\psi_{\mathrm{Ludwig}}-\psi_{\mathrm{formula}}|$")
        b.set_title("Error of the formula vs the simulation"); b.legend(fontsize=8.5)
    fig.suptitle("1:1 electrolyte, plates at $\\psi=0.1$ and $0.3$, point charge at d = 8: Ludwig vs closed-system formula, "
                 "without and with lateral images", y=0.995)
    fig.tight_layout(); fig.savefig(os.path.join(F.ROOT, "graficos", "C_dospot_electrolyte_vs_lambdaD.png"), dpi=160, bbox_inches="tight")
