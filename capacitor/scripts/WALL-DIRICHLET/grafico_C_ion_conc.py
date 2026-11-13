#!/usr/bin/env python3
"""Ion concentrations n+/n0 and n-/n0 of the closed-system electrolyte capacitor (plates at psi=0.1 and 0.3, lambda_D=6,
point charge q/eps=10 at d=8 from the lower plate) -- Ludwig (qsi file, step 5000) vs the linearised closed-system theory

   n_pm = n0 [ 1 -/+ (phi - <phi>) ],   phi = phi0 + sum of lateral images + closed-system correction  (same phi as the potential figure)

Reduced units psi = e phi / kBT.  n0 = mean ion density (= initial rho_el, ion number conserved).
Usage: python3 grafico_C_ion_conc.py -> ../graficos/C_ion_concentration_closed_2pot.png
"""
import os, sys
BAT = "C-electrolito-dospotenciales-hann8-LD6"
sys.argv = [sys.argv[0], BAT, "0.1", "0.3", "x.png"]
import glob
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import grafico_C_formulas as F
import grafico_C_formulas_closed as C

L, D, KAP, P1, P2 = F.L, F.D, C.KAP, C.P1, C.P2
N0 = 1.38888889e-03

def qsi(bat, dx):
    p = sorted(glob.glob(os.path.join(bat, "z_%05.2f" % (F.D_LOWER + 1.5), "qsi-0*")))[-1]
    a = np.loadtxt(p).reshape(32, 32, 34, 2)
    zc = np.arange(1, 35) - 1.5
    return zc[1:33] - L / 2, a[(15 + dx) % 32, 15, 1:33, :] / N0

if __name__ == "__main__":
    plt.rcParams.update({"font.size": 11, "axes.grid": True, "grid.color": "#e4e2dc",
                         "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False})
    bat = os.path.join(F.ROOT, BAT)
    phibar = (KAP * L) / (2 * np.tanh(KAP * L / 2)) * C.mean_phi0()
    cols = {0: "#1baf7a", 6: "#2a78d6", 10: "#eb6834"}
    zz = np.linspace(-L / 2 + 0.2, L / 2 - 0.2, 300)
    fig, axs = plt.subplots(1, 2, figsize=(12.5, 5.0), sharey=True)
    for k, (a, sgn, lab) in enumerate(zip(axs, (-1, +1), ("$n_+/n_0$ (positive ions)", "$n_-/n_0$ (negative ions)"))):
        for dx, c in cols.items():
            z, n = qsi(bat, dx)
            phi0 = lambda zv: F.phi_cap_elec(zv, P1, P2, KAP) + F.phi_q_AB(float(dx), zv, D, KAP)
            rep = lambda zv: F.replicas_AB(float(dx), zv, KAP)
            a.plot(zz, 1 + sgn * (phi0(zz) + C.closed(zz) - phibar), "--", color=c, lw=1.3, dashes=(5, 3))
            a.plot(zz, 1 + sgn * (phi0(zz) + rep(zz) + C.closed(zz) - phibar), "-", color=c, lw=1.4, alpha=0.55)
            a.plot(z, n[:, k], "o", color=c, ms=4.2, mec="white", mew=0.5, label=r"Ludwig, $\rho=%d$" % dx)
        a.axhline(1, color="#888", lw=0.8)
        a.set_xlabel("z (distance from the mid-plane, lattice units)"); a.set_title(lab)
    axs[0].set_ylim(0.7, 1.3)    # theory is singular at the point charge (rho=0, z=d): clipped
    axs[0].set_ylabel(r"ion concentration $n_\pm/n_0$")
    axs[0].plot([], [], "--", color="#55544f", lw=1.3, dashes=(5, 3), label="linearised closed-system theory, infinite plates")
    axs[0].plot([], [], "-", color="#55544f", lw=1.4, alpha=0.55, label="theory + lateral images (10 per side)")
    axs[0].legend(loc="upper right", fontsize=9)
    fig.suptitle(r"1:1 electrolyte, $\lambda_D=6$, plates at $\psi=0.1$ and $0.3$, charge at d = 8 from the lower plate", fontsize=12)
    fig.tight_layout()
    out = os.path.join(F.ROOT, "graficos", "C_ion_concentration_closed_2pot.png")
    fig.savefig(out, dpi=170); print(out)
    # residuals outside the kernel support
    for dx in cols:
        z, n = qsi(bat, dx); sel = np.abs(z - D) > 4.5
        phi = F.phi_cap_elec(z, P1, P2, KAP) + F.phi_q_AB(float(dx), z, D, KAP) + F.replicas_AB(float(dx), z, KAP) + C.closed(z) - phibar
        print(" rho=%2d  max|n+ - theory| = %.4f   max|n- - theory| = %.4f   range n+ = %.3f..%.3f" %
              (dx, np.abs(n[:, 0] - (1 - phi))[sel].max(), np.abs(n[:, 1] - (1 + phi))[sel].max(), n[:, 0].min(), n[:, 0].max()))
