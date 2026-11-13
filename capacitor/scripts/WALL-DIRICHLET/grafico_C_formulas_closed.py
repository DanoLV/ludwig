#!/usr/bin/env python3
"""Electrolyte capacitor with a point charge: Ludwig potential vs the CLOSED-system (fixed ion number) formula

  phi(rho,z) = phi0(rho,z) + [kL / (2 tanh(kL/2))] <phi0> [1 - cosh(k z)/cosh(k L/2)],     k = kappa
  phi0 = phi_CAP(z) + sum_n [phi^A_n + phi^B_n]       (open-system solution, point charge)
  <phi0> = <phi_q> + <phi_CAP>,  <phi_CAP> = (p1+p2) tanh(kL/2)/(kL),
  <phi_q> = (4q/(pi eps V)) sum_{n>=0} (-1)^n/(2n+1) cos((2n+1) pi d/L)/(((2n+1)pi/L)^2 + k^2),  V = l^2 L
Dashed: phi0 of the single cell (infinite plates); solid: phi0 summed over lateral images (10 per side).
Reduced units, q/eps = 10.  Data: C-electrolito-hann8-LD6 (Hann-8, d = 8 from the lower plate).
Usage: python3 grafico_C_formulas_closed.py -> ../graficos/C_formula_electrolyte_closed_vs_ludwig.png
"""
import os, sys, math
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import grafico_C_formulas as F

L, A, FOUR, D = F.L, F.A, F.FOUR, F.D
KAP = 1 / 6.0
# optional arguments: battery folder, psi_1, psi_2, output name   (defaults: grounded plates, original figure)
BAT_NAME = sys.argv[1] if len(sys.argv) > 1 else "C-electrolito-hann8-LD6"
P1 = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0
P2 = float(sys.argv[3]) if len(sys.argv) > 3 else 0.0
OUT_NAME = sys.argv[4] if len(sys.argv) > 4 else "C_formula_electrolyte_closed_vs_ludwig.png"
GROUNDED = (P1 == 0.0 and P2 == 0.0)

def _sy(y):
    if F.LAYOUT != "left": return y
    if not GROUNDED:        # PPT two-potential figure: compact (no "change vs open" heading)
        if y > -0.42: return y
        return -0.46 + 0.80 * (y + 0.42) - (0.04 if y <= -1.37 else 0) - (0.08 if y <= -1.43 else 0)
    return y - (0.08 if y <= -0.33 else 0) - (0.1 if y <= -1.37 else 0) - (0.08 if y <= -1.43 else 0)

def mean_phi_q():
    n = np.arange(0, 20000); m = 2 * n + 1; V = A * L
    return 4 * FOUR / (np.pi * V) * np.sum((-1.0) ** n / m * np.cos(m * np.pi * D / L) / ((m * np.pi / L) ** 2 + KAP ** 2))

def mean_phi_cap():
    return (P1 + P2) * np.tanh(KAP * L / 2) / (KAP * L)

def mean_phi0():
    return mean_phi_q() + mean_phi_cap()

def closed(z):
    fac = (KAP * L) / (2 * np.tanh(KAP * L / 2)) * mean_phi0()
    return fac * (1 - np.cosh(KAP * z) / np.cosh(KAP * L / 2))

if __name__ == "__main__":
    plt.rcParams.update({"font.size": 11, "axes.grid": True, "grid.color": "#e4e2dc",
                         "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False})
    fig, a = F.new_fig(*((1.75, 12.2) if GROUNDED else (1.42, 8.2)))
    cols = {0: "#1baf7a", 6: "#2a78d6", 10: "#eb6834"}
    zz = np.linspace(-L / 2 + 0.2, L / 2 - 0.2, 200)
    bat = os.path.join(F.ROOT, BAT_NAME)
    print("<phi0> = %.6f ; phibar = %.6f" % (mean_phi0(), (KAP * L) / (2 * np.tanh(KAP * L / 2)) * mean_phi0()))
    print(" rho   max|Ludwig - formula|: infinite plates | with images   (outside kernel support)   max(Ludwig)")
    for dx, c in cols.items():
        z, ps = F.ludwig(bat, dx)
        f0 = lambda zv: F.phi_cap_elec(zv, P1, P2, KAP) + F.phi_q_AB(float(dx), zv, D, KAP)
        rep = lambda zv: F.replicas_AB(float(dx), zv, KAP)
        a.plot(zz, f0(zz) + closed(zz), "--", color=c, lw=1.3, dashes=(5, 3), zorder=2)
        a.plot(zz, f0(zz) + rep(zz) + closed(zz), "-", color=c, lw=1.4, alpha=0.55, zorder=2)
        a.plot(z, ps, "o", ms=4.6, color=c, mec="white", mew=0.5, zorder=4, label="Ludwig, $\\rho=%d$" % dx)
        sel = np.abs(z - D) > 4.5
        print("  %2d     %.4f | %.4f   (%.4f)" % (dx, np.abs(ps - f0(z) - closed(z))[sel].max(),
              np.abs(ps - f0(z) - rep(z) - closed(z))[sel].max(), ps.max()))
    a.plot([], [], "--", color="#55544f", lw=1.3, dashes=(5, 3), label="closed-system formula, infinite plates")
    a.plot([], [], "-", color="#55544f", lw=1.4, alpha=0.55, label="closed-system formula + lateral images (10 per side)")
    a.set_xlabel("z (distance from the mid-plane, lattice units)"); a.set_ylabel(r"reduced potential $\psi = e\phi/k_BT$")
    a.set_title("1:1 electrolyte, $\\lambda_D=6$, " + ("grounded plates" if GROUNDED else "plates at $\\psi=%g$ and $%g$" % (P1, P2)) + " (d = 8 from the lower plate)")
    a.set_ylim(*((-0.01, 0.45) if GROUNDED else (0.05, 0.62))); a.legend(loc="upper right", fontsize=(8.5 if F.LAYOUT == "left" and not GROUNDED else 9), bbox_to_anchor=(1.0, 1.0 if F.LAYOUT == "left" and not GROUNDED else 1.0))
    fig.tight_layout()
    if F.LAYOUT == "left" and not GROUNDED:      # shorter plot, centred vertically
        p = a.get_position(); h = p.height * 0.72
        a.set_position([p.x0, p.y0 + (p.height - h) / 2, p.width, h])
    fs = 10.5
    F.put(fig, 0.5, _sy(-0.02), r"Poisson equation with ion-number conservation (closed system: $N_{\pm}$ fixed)",
             ha="center", va="top", fontsize=fs, fontweight="bold")
    F.put(fig, 0.5, _sy(-0.07), r"$\left(\nabla^{2}-\kappa^{2}\right)\phi=\;-\kappa^{2}\langle\phi\rangle\;-\dfrac{\rho_{\mathrm{part}}}{\varepsilon}"
             r"\;\Longleftrightarrow\;\nabla^{2}\phi-\kappa^{2}\left(\phi-\langle\phi\rangle\right)=-\dfrac{\rho_{\mathrm{part}}}{\varepsilon}$",
             ha="center", va="top", fontsize=fs)
    F.put(fig, 0.5, _sy(-0.20), r"$\rho_{\mathrm{part}}=q\,\delta(\mathbf{r}-\mathbf{r}_0)$, $\phi(z=\mp L/2)=\phi_{1,2}$, periodic in $x,y$; "
             r"$n_{\pm}=n_0\left[1\mp\dfrac{e}{k_BT}\left(\phi-\langle\phi\rangle\right)\right]$, $\int n_{\pm}\,dV=N$, "
             r"$\langle\phi\rangle=\dfrac{1}{V}\int_V\phi\,dV$", ha="center", va="top", fontsize=fs - 0.5)
    if not (F.LAYOUT == "left" and not GROUNDED):   # PPT two-potential figure: formulas only
      F.put(fig, 0.5, _sy(-0.33), r"Change vs the open case: extra source term $-\kappa^{2}\langle\phi\rangle$ in the differential equation",
               ha="center", va="top", fontsize=fs, color="#c0392b", fontweight="bold")
    F.put(fig, 0.5, _sy(-0.42), r"$\phi(x,y,z)=\phi_0(x,y,z)+\dfrac{\kappa L}{2\tanh(\kappa L/2)}\,\langle\phi_0\rangle"
             r"\left[1-\dfrac{\cosh(\kappa z)}{\cosh(\kappa L/2)}\right]$", ha="center", va="top", fontsize=fs)
    F.put(fig, 0.5, _sy(-0.55), r"$\phi_0(x,y,z)=\phi_{\mathrm{CAP}}(z)+\sum_{i,j}\phi_q^{p}(\rho_{ij},z),\qquad "
             r"\phi_{\mathrm{CAP}}(z)=\dfrac{\phi_1\sinh\left[\kappa(L/2-z)\right]+\phi_2\sinh\left[\kappa(L/2+z)\right]}{\sinh(\kappa L)}$",
             ha="center", va="top", fontsize=fs)
    F.put(fig, 0.5, _sy(-0.69), r"$\phi_q^{p}=\sum_{n=0}^{\infty}\left[\phi^{A}_{q_n}+\phi^{B}_{q_n}\right],\qquad "
             r"\kappa^{2}=\dfrac{2e^{2}n_0}{\varepsilon k_B T}$", ha="center", va="top", fontsize=fs)
    F.put(fig, 0.5, _sy(-0.81), r"$\phi^{A}_{q_n}=\dfrac{q}{\pi L\varepsilon}\,\sin\left(\dfrac{2\pi n}{L}d\right)\sin\left(\dfrac{2\pi n}{L}z\right)"
             r"K_0\left(\sqrt{\left(\dfrac{2\pi n}{L}\right)^{2}+\kappa^{2}}\,\rho\right)$", ha="center", va="top", fontsize=fs)
    F.put(fig, 0.5, _sy(-0.94), r"$\phi^{B}_{q_n}=\dfrac{q}{\pi L\varepsilon}\,\cos\left(\dfrac{(2n+1)\pi}{L}d\right)\cos\left(\dfrac{(2n+1)\pi}{L}z\right)"
             r"K_0\left(\sqrt{\left(\dfrac{(2n+1)\pi}{L}\right)^{2}+\kappa^{2}}\,\rho\right)$", ha="center", va="top", fontsize=fs)
    F.put(fig, 0.5, _sy(-1.08), r"$\langle\phi_0\rangle=\langle\phi_q\rangle+\langle\phi_{\mathrm{CAP}}\rangle,\qquad "
             r"\langle\phi_{\mathrm{CAP}}\rangle=\dfrac{(\phi_1+\phi_2)\tanh(\kappa L/2)}{\kappa L}$",
             ha="center", va="top", fontsize=fs)
    F.put(fig, 0.5, _sy(-1.21), r"$\langle\phi_q\rangle=\dfrac{4q}{\pi\varepsilon V}\sum_{n=0}^{\infty}\dfrac{(-1)^{n}}{2n+1}\,"
             r"\dfrac{\cos\left(\frac{(2n+1)\pi d}{L}\right)}{\left(\frac{(2n+1)\pi}{L}\right)^{2}+\kappa^{2}},\qquad V=\ell^{2}L$",
             ha="center", va="top", fontsize=fs)
    F.put(fig, 0.5, _sy(-1.37), r"$\rho_{ij}=\sqrt{(x-i\ell)^{2}+(y-j\ell)^{2}}$, $\ell=32$. Solid: $|i|,|j|\leq 10$ (lateral images); "
             r"dashed: only $(i,j)=(0,0)$, i.e. $\rho_{00}=\sqrt{x^{2}+y^{2}}$ (infinite plates)",
             ha="center", va="top", fontsize=fs - 0.5)
    F.put(fig, 0.5, _sy(-1.43), (r"Plates at $z=\mp L/2$ ($L=32$), $\psi_1=%g$, $\psi_2=%g$; " % (P1, P2)) + r"$\lambda_D=1/\kappa=6$; $d=-8$ (charge position); "
             r"point charge, $q/\varepsilon=10$. Field points on the line $y=0$, $x=\rho$", ha="center", va="top", fontsize=fs - 0.5)
    fig.savefig(os.path.join(F.ROOT, "graficos", F.out_name(OUT_NAME)), dpi=170, bbox_inches="tight")
