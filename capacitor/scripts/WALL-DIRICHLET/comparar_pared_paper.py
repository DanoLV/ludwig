#!/usr/bin/env python3
"""Wall potential of Ludwig (wall_ss_cut, soft sphere with energy+force shift) vs the ion-wall potential of
Hoang Ngoc Minh et al., J. Chem. Phys. 158, 104103 (2023), Eqs. (6)-(8):
   V_ads(z) = 2 pi eps_w [ (2/5)(s/z)^10 - (s/z)^4 - sqrt(2) / (3 (z/s + 0.61/sqrt(2))^3) ]     (10-4-3 Steele, fcc (100) face, lattice parameter s = sigma_w)
   V_rep(z) = V_ads(z) - V_ads(z*)  for z <= z* = 0.987 s (minimum of V_ads);  0 for z > z*       (purely repulsive wall)
Ludwig: v(h) = eps (s/h)^nu - eps (s/hc)^nu + (eps nu/s)(s/hc)^(nu+1) (h - hc), h < hc;  eps/kT = 1.3105e-16, s = 1, nu = 12, hc = 3.
Usage: python3 comparar_pared_paper.py -> ../graficos/comparacion_potencial_pared_paper.png"""
import os
import numpy as np
from scipy.optimize import minimize_scalar, brentq
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

def vads(z, s=1.0):                         # in units of eps_w
    return 2 * np.pi * (0.4 * (s / z) ** 10 - (s / z) ** 4 - np.sqrt(2) / (3 * (z / s + 0.61 / np.sqrt(2)) ** 3))
zmin = minimize_scalar(vads, bounds=(0.5, 1.5), method="bounded", options={"xatol": 1e-12}).x
def vrep(z, s=1.0):                         # units of eps_w
    z = np.asarray(z, float); return np.where(z <= zmin * s, vads(z, s) - vads(zmin * s, s), 0.0)
def frep(z, s=1.0, dz=1e-6):                # -dV/dz  (units eps_w / sigma_w)
    return np.where(np.asarray(z) <= zmin * s, -(vads(z + dz, s) - vads(z - dz, s)) / (2 * dz), 0.0)

EPS, SIG, NU, HC = 1.3105e-16, 1.0, 12.0, 3.0      # eps in kT
def vss(h):
    h = np.asarray(h, float); vc = EPS * (SIG / HC) ** NU; dvc = -EPS * NU / SIG * (SIG / HC) ** (NU + 1)
    return np.where(h < HC, EPS * (SIG / h) ** NU - vc - (h - HC) * dvc, 0.0)
def fss(h):
    h = np.asarray(h, float); dvc = -EPS * NU / SIG * (SIG / HC) ** (NU + 1)
    return np.where(h < HC, EPS * NU / SIG * (SIG / h) ** (NU + 1) + dvc, 0.0)

if __name__ == "__main__":
    print("z* (minimum of V_ads) = %.4f sigma_w ;  V_ads(z*) = %.4f eps_w   (paper: z* ~ 0.987 sigma_w)" % (zmin, vads(zmin)))
    for be in (0.17, 0.50):
        z1 = brentq(lambda z: be * vrep(z) - 1, 0.5, zmin - 1e-9)
        print("paper beta eps_w = %.2f : V_rep = 1 kT at z = %.4f sigma_w ; V_rep(0.8 s) = %.3f kT" % (be, z1, be * float(vrep(0.8))))
    h1 = brentq(lambda h: vss(h) - 1, 1e-3, 1.0)
    print("Ludwig: v = 1 kT at h = %.4f ; v(0.1) = %.3e kT ; F(0.1) = %.3e kT/length ; range hc = %.1f" % (h1, vss(0.1), fss(0.1), HC))
    plt.rcParams.update({"font.size": 10.5, "axes.grid": True, "grid.color": "#e4e2dc",
                         "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False})
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))
    zs = np.linspace(0.45, 1.6, 600); hs = np.logspace(-3, np.log10(3.4), 600)
    cs = {0.17: "#2a78d6", 0.50: "#1baf7a"}
    for be, c in cs.items():
        ax[0].plot(zs, be * vads(zs), color=c, lw=1.2, ls=":", label=r"paper $V_w^{ads}$, $\beta\varepsilon_w=%.2f$" % be)
        ax[0].plot(zs[zs <= zmin], be * vrep(zs[zs <= zmin]), color=c, lw=2.0, label=r"paper $V_w^{rep}$ (truncated, shifted), $\beta\varepsilon_w=%.2f$" % be)
    ax[0].set(ylim=(-1.5, 12), xlim=(0.45, 1.6), xlabel=r"distance to the wall, $z/\sigma_w$ (paper)", ylabel=r"$\beta V$", title="(a) Paper: ion-wall potential (linear)")
    ax[0].legend(fontsize=8, loc="upper right")
    for be, c in cs.items():
        zz = np.linspace(0.3, zmin, 400)
        ax[1].semilogy(zz, be * vrep(zz) + 1e-30, color=c, lw=2.0, label=r"paper $V_w^{rep}$, $\beta\varepsilon_w=%.2f$ (range $0.987\,\sigma_w$)" % be)
    ax[1].semilogy(hs, vss(hs), color="#eb6834", lw=2.0, label=r"Ludwig soft sphere, $\varepsilon/k_BT=1.3\times10^{-16}$, $\nu=12$, $h_c=3$")
    ax[1].axhline(1, color="#888", lw=0.8); ax[1].text(0.011, 3, r"$1\,k_BT$", color="#666")
    ax[1].axvline(0.1, color="#aaa", ls="--", lw=0.9); ax[1].text(0.105, 1e-30 * 3e15, "Ludwig stop\n(h = 0.1)", fontsize=8.5, color="#666")
    ax[1].set(xscale="log", xlim=(0.01, 3.4), ylim=(1e-30, 1e8), xlabel=r"distance: $z/\sigma_w$ (paper)   |   $h=d-R$, in lattice units (Ludwig)",
              ylabel=r"$\beta V$ (log)", title="(b) Both potentials (log-log)")
    ax[1].legend(fontsize=8, loc="upper right")
    for be, c in cs.items():
        zz = np.linspace(0.3, zmin, 400)
        ax[2].semilogy(zz, be * frep(zz) + 1e-30, color=c, lw=2.0)
    ax[2].semilogy(hs, fss(hs), color="#eb6834", lw=2.0)
    ax[2].axhline(1.5e-7 / 1e-5, color="#aaa", ls="--", lw=0.9); ax[2].text(0.011, 2.6e-2 * 1.4, r"net electrostatic force on the particle $\approx 0.015\,k_BT$/length", fontsize=8.5, color="#666")
    ax[2].set(xscale="log", xlim=(0.01, 3.4), ylim=(1e-12, 1e10), xlabel=ax[1].get_xlabel(), ylabel=r"repulsive force $\beta F$ (log), units $1$/length", title="(c) Repulsive force")
    fig.suptitle("Wall repulsion: Ludwig soft sphere (wall_ss_cut) vs ion-wall potential of Hoang Ngoc Minh et al., J. Chem. Phys. 158, 104103 (2023)", y=1.02)
    fig.tight_layout(); out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "graficos", "comparacion_potencial_pared_paper.png")
    fig.savefig(out, dpi=160, bbox_inches="tight"); print(out)
