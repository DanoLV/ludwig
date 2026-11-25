#!/usr/bin/env python3
"""Ludwig potential vs the analytical formulas (plates at z = -/+ L/2, charge at z = d measured from the centre).

(a) No electrolyte, plates at psi_1 = 0.1 (lower) and psi_2 = 0.3 (upper), Gaussian charge (sigma of Hann-8):
      phi_TOT = phi_CAP(z) + phi_G(rho, z),   phi_CAP = (p1+p2)/2 + (p2-p1) z/L
      phi_G   = l_B q sum_n [ erf(r_n^+/(sqrt2 sigma))/r_n^+ - erf(r_n^-/(sqrt2 sigma))/r_n^- ]
      r_n^+ = sqrt(rho^2 + (z-(2nL+d))^2),  r_n^- = sqrt(rho^2 + (z-((2n+1)L-d))^2)
(b) 1:1 electrolyte, grounded plates (as written on the sheet, NO mean-potential term):
      phi_CAP = [p1 sinh k(L/2-z) + p2 sinh k(L/2+z)]/sinh kL ,  phi_TOT = phi_CAP + sum_n [phi^A_n + phi^B_n]
Data: C-dospotenciales-hann08 and C-electrolito-hann8-LD6 (Hann-8, d = 8 from the lower plate, L = 32, box 32x32).
Reduced units (psi = e phi / k_B T, 4 pi l_B = 10).   Usage: python3 grafico_C_formulas.py
"""
import os, glob, math
import numpy as np
from scipy.special import erf as verf, k0 as _k0       # erf and modified Bessel K0 (scipy.special)
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
# 2026-11-16: absolute path. Since the 2026-11-13 move to capacitor/scripts/ Python resolves the scripts/ symlink,
# so dirname(HERE) pointed to capacitor/scripts instead of the runs. Original: ROOT = os.path.dirname(HERE)
ROOT = "/home/danolv/Sim/ludwig/capacitor/local/WALL-DIRICHLET"
L = 32.0; A = 32.0 ** 2; FOUR = 10.0; LB = FOUR / (4 * np.pi)
SIG = 8 * math.sqrt(1 / 12 - 1 / (2 * math.pi ** 2))
D_LOWER = 8.0; D = D_LOWER - L / 2            # charge z-coordinate measured from the centre

def phi_cap_lin(z, p1, p2):
    return (p1 + p2) / 2 + (p2 - p1) * z / L

def phi_G(rho, z, d=D, nmax=3000):
    n = np.arange(-nmax, nmax + 1); out = np.empty(len(z))
    for i, zz in enumerate(z):
        rp = np.hypot(rho, zz - (2 * n * L + d)); rm = np.hypot(rho, zz - ((2 * n + 1) * L - d))
        fp = np.where(rp < 8 * SIG, verf(np.minimum(rp, 8 * SIG) / (math.sqrt(2) * SIG)) / np.maximum(rp, 1e-9), 1 / rp)
        fm = np.where(rm < 8 * SIG, verf(np.minimum(rm, 8 * SIG) / (math.sqrt(2) * SIG)) / np.maximum(rm, 1e-9), 1 / rm)
        out[i] = LB * (fp.sum() - fm.sum())
    return out

def phi_q_AB(rho, z, d, kap, N=400):
    rho = max(rho, 0.02); N = 2500 if rho < 0.5 else N      # rho = 0: series evaluated at rho = 0.02
    out = np.zeros(len(z))
    for n in range(1, N + 1):                  # family A (n = 0 term vanishes)
        g = np.sqrt((2 * np.pi * n / L) ** 2 + kap ** 2); k0 = float(_k0(g * rho))
        out += np.sin(2 * np.pi * n * d / L) * np.sin(2 * np.pi * n * z / L) * k0
    for n in range(0, N + 1):                  # family B
        a = (2 * n + 1) * np.pi / L; g = np.sqrt(a ** 2 + kap ** 2); k0 = float(_k0(g * rho))
        out += np.cos(a * d) * np.cos(a * z) * k0
    return FOUR / (np.pi * L) * out

NREP = 10                                          # lateral replicas on each side: |i|,|j| <= 10

def replicas_G(dx, z, nmax=1500):
    n = np.arange(-nmax, nmax + 1)[None, :]; zc = z[:, None]; out = np.zeros(len(z))
    for i in range(-NREP, NREP + 1):
        for j in range(-NREP, NREP + 1):
            if i == 0 and j == 0: continue
            rho = math.hypot(dx - i * 32.0, j * 32.0)
            rp = np.hypot(rho, zc - (2 * n * L + D)); rm = np.hypot(rho, zc - ((2 * n + 1) * L - D))
            out += LB * (1 / rp - 1 / rm).sum(axis=1)
    return out

def replicas_AB(dx, z, kap, N=60):
    out = np.zeros(len(z))
    for i in range(-NREP, NREP + 1):
        for j in range(-NREP, NREP + 1):
            if i == 0 and j == 0: continue
            rho = math.hypot(dx - i * 32.0, j * 32.0)
            out += phi_q_AB(rho, z, D, kap, N=N)
    return out

def phi_cap_elec(z, p1, p2, kap):
    return (p1 * np.sinh(kap * (L / 2 - z)) + p2 * np.sinh(kap * (L / 2 + z))) / np.sinh(kap * L)

def closed_term(z, d, kap):
    mo = np.arange(1, 4001, 2); g2 = (mo * np.pi / L) ** 2 + kap ** 2
    mean = FOUR / (np.pi * L) * (2 * np.pi / (A * L)) * np.sum(np.sin(mo * np.pi * (d + L / 2) / L) * (2 * L / (mo * np.pi)) / g2)
    pbar = (kap * L / 2) / np.tanh(kap * L / 2) * mean
    return pbar * (1 - np.cosh(kap * z) / np.cosh(kap * L / 2))

def ludwig(bat, dx):
    p = sorted(glob.glob(os.path.join(bat, "z_%05.2f" % (D_LOWER + 1.5), "psi-0*")))[-1]
    a = np.loadtxt(p).reshape(32, 32, 34)
    zc = np.arange(1, 35) - 1.5
    return zc[1:33] - L / 2, a[(15 + dx) % 32, 15, 1:33]


# --- layout: formulas under the plot (default) or at the LEFT of the plot (slides): LAYOUT=left python3 ... -> *_ppt.png ---
LAYOUT = os.environ.get("LAYOUT", "bottom")
_AXT = [None]; _YMAX = [1.0]

def new_fig(ymax, height=7.2):
    _YMAX[0] = ymax
    if LAYOUT == "left":
        fig = plt.figure(figsize=(16.5, height)); gs = fig.add_gridspec(1, 2, width_ratios=[1.3, 1], wspace=0.04)
        _AXT[0] = fig.add_subplot(gs[0]); _AXT[0].axis("off"); a = fig.add_subplot(gs[1])
        return fig, a
    return plt.subplots(figsize=(7.2, 5.0))

def _wrap(s, n=92):
    """Break a long note at spaces outside $...$ so it fits in the left panel."""
    toks, cur, inm = [], "", False
    for ch in s:
        if ch == "$": inm = not inm
        if ch == " " and not inm: toks.append(cur); cur = ""
        else: cur += ch
    toks.append(cur); lines, line = [], ""
    for t in toks:
        cand = (line + " " + t).strip()
        if len(cand.replace("$", "")) > n and line: lines.append(line); line = t
        else: line = cand
    lines.append(line)
    return "\n".join(lines)

def put(fig, x, y, s, **kw):
    if LAYOUT == "left":
        kw = dict(kw)
        if kw.get("fontsize", 10.5) < 10.4: s = _wrap(s)            # notes (smaller font) are wrapped
        kw["fontsize"] = kw.get("fontsize", 10.5) * 0.9
        _AXT[0].text(0.5, 1.0 + y / _YMAX[0], s, transform=_AXT[0].transAxes, linespacing=1.5, **kw)
    else:
        fig.text(x, y, s, **kw)

def out_name(name):
    return name.replace(".png", "_ppt.png") if LAYOUT == "left" else name

def style():
    plt.rcParams.update({"font.size": 11, "axes.grid": True, "grid.color": "#e4e2dc",
                         "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False})

COLS = {0: "#1baf7a", 6: "#2a78d6", 10: "#eb6834"}

def make(which):
    style(); fig, a = new_fig(0.90 if which == "a" else 1.12)
    zz = np.linspace(-L / 2 + 0.2, L / 2 - 0.2, 200)
    kap = 1 / 6.0; P1, P2 = 0.1, 0.3
    bat = os.path.join(ROOT, "C-dospotenciales-hann08" if which == "a" else "C-electrolito-hann8-LD6")
    print("figure", which, " rho   max|Ludwig - formula|   max(Ludwig)")
    for dx, c in COLS.items():
        z, ps = ludwig(bat, dx)
        if which == "a":
            f = lambda zv: phi_cap_lin(zv, P1, P2) + phi_G(float(dx), zv)
            rep = lambda zv: replicas_G(float(dx), zv)
        else:
            f = lambda zv: phi_cap_elec(zv, 0.0, 0.0, kap) + phi_q_AB(float(dx), zv, D, kap)
            rep = lambda zv: replicas_AB(float(dx), zv, kap)
        a.plot(zz, f(zz), "--", color=c, lw=1.3, dashes=(5, 3), zorder=2)            # infinite plates: thin dashed
        a.plot(zz, f(zz) + rep(zz), "-", color=c, lw=1.4, alpha=0.55, zorder=2)      # with lateral images: thick solid
        a.plot(z, ps, "o", ms=4.6, color=c, mec="white", mew=0.5, zorder=4, label="Ludwig, $\\rho=%d$" % dx)
        sel = np.abs(z - D) > 4.5
        print("   %3d   no images %.4f | with images %.4f (outside kernel support)   max(Ludwig) %.4f"
              % (dx, np.abs(ps - f(z))[sel].max(), np.abs(ps - f(z) - rep(z))[sel].max(), ps.max()))
    a.plot([], [], "--", color="#55544f", lw=1.3, dashes=(5, 3), label="analytical formula, infinite plates")
    a.plot([], [], "-", color="#55544f", lw=1.4, alpha=0.55, label="formula + lateral images (10 per side)")
    a.set_xlabel("z (distance from the mid-plane, lattice units)")
    a.set_ylabel(r"reduced potential $\psi = e\phi/k_BT$")
    if which == "a":
        a.set_title("No electrolyte, plates at $\\psi=0.1$ and $0.3$ (d = 8 from the lower plate)")
        a.legend(loc="upper right" if LAYOUT == "left" else "upper left", fontsize=9.5); name = "C_formula_noelectrolyte_vs_ludwig.png"
    else:
        a.set_title("1:1 electrolyte, $\\lambda_D=6$, grounded plates (d = 8 from the lower plate)")
        a.set_ylim(-0.01, 0.45); a.legend(loc="upper right", fontsize=9.5); name = "C_formula_electrolyte_vs_ludwig.png"
    fig.tight_layout()
    if which == "a":                           # formula used, written below the plot
        fs = 10.5
        put(fig, 0.5, -0.02, r"Poisson equation, no electrolyte (Gaussian charge between two plates at fixed potential)",
                 ha="center", va="top", fontsize=fs, fontweight="bold")
        put(fig, 0.5, -0.07, r"$\nabla^{2}\phi=-\dfrac{\rho_{\mathrm{part}}}{\varepsilon},\qquad "
                 r"\rho_{\mathrm{part}}=q\,\dfrac{e^{-|\mathbf{r}-\mathbf{r}_0|^{2}/(2\sigma^{2})}}{(2\pi\sigma^{2})^{3/2}}$",
                 ha="center", va="top", fontsize=fs)
        put(fig, 0.5, -0.20, r"$\phi(z=\mp L/2)=\phi_{1,2}$; laterally unbounded (dashed) or periodic with period $\ell$ (solid)",
                 ha="center", va="top", fontsize=fs - 0.5)
        put(fig, 0.5, -0.32, r"$\phi_{\mathrm{TOT}}(x,y,z)=\phi_{\mathrm{CAP}}(z)+\sum_{i,j}\phi_G(\rho_{ij},z),\qquad"
                 r"\phi_{\mathrm{CAP}}(z)=\dfrac{\phi_1+\phi_2}{2}+\dfrac{\phi_2-\phi_1}{L}\,z$",
                 ha="center", va="top", fontsize=fs)
        put(fig, 0.5, -0.43, r"$\phi_G(\rho,z)=\dfrac{q}{4\pi\varepsilon}\sum_{n=-\infty}^{\infty}"
                 r"\left[\dfrac{\mathrm{erf}\left(r_n^{+}/\sqrt{2}\,\sigma\right)}{r_n^{+}}"
                 r"-\dfrac{\mathrm{erf}\left(r_n^{-}/\sqrt{2}\,\sigma\right)}{r_n^{-}}\right]$",
                 ha="center", va="top", fontsize=fs)
        put(fig, 0.5, -0.57, r"$r_n^{+}=\sqrt{\rho^{2}+\left[z-(2nL+d)\right]^{2}},\qquad"
                 r"r_n^{-}=\sqrt{\rho^{2}+\left[z-\left((2n+1)L-d\right)\right]^{2}}$",
                 ha="center", va="top", fontsize=fs)
        put(fig, 0.5, -0.65, r"$\rho_{ij}=\sqrt{(x-i\ell)^{2}+(y-j\ell)^{2}}$, $\ell=32$. Solid: $|i|,|j|\leq 10$ (lateral images); "
                 r"dashed: only $(i,j)=(0,0)$, i.e. $\rho_{00}=\sqrt{x^{2}+y^{2}}$ (infinite plates)",
                 ha="center", va="top", fontsize=fs - 0.5)
        put(fig, 0.5, -0.71, r"Plates at $z=\mp L/2$ ($L=32$), $\psi_1=0.1$, $\psi_2=0.3$; $d=-8$ (charge position); "
                 r"$\sigma=1.446$ (Hann-8 width); $q/\varepsilon=10$. Field points on the line $y=0$, $x=\rho$",
                 ha="center", va="top", fontsize=fs - 0.5)
    if which == "b":                           # formula used, written below the plot
        fs = 10.5
        put(fig, 0.5, -0.02, r"Poisson equation, no ion-number conservation (open system: reservoir at fixed $\mu$)",
                 ha="center", va="top", fontsize=fs, fontweight="bold")
        put(fig, 0.5, -0.07, r"$\left(\nabla^{2}-\kappa^{2}\right)\phi=-\dfrac{\rho_{\mathrm{part}}}{\varepsilon},\qquad "
                 r"\rho_{\mathrm{part}}=q\,\delta(\mathbf{r}-\mathbf{r}_0),\quad \phi(z=\mp L/2)=\phi_{1,2},\quad \mathrm{periodic\ in}\ x,y$",
                 ha="center", va="top", fontsize=fs)
        put(fig, 0.5, -0.20, r"from $n_{\pm}=n_0\left[1\mp\dfrac{e\,\phi}{k_BT}\right]$,  $\kappa^{2}=\dfrac{2e^{2}n_0}{\varepsilon k_BT}$  "
                 r"(ion numbers are not conserved)", ha="center", va="top", fontsize=fs - 0.5)
        put(fig, 0.5, -0.34, r"$\phi_{\mathrm{TOT}}(x,y,z)=\phi_{\mathrm{CAP}}(z)+\sum_{i,j}\phi_q^{p}(\rho_{ij},z),\qquad"
                 r"\phi_{\mathrm{CAP}}(z)=\dfrac{\phi_1\sinh\left[\kappa(L/2-z)\right]+\phi_2\sinh\left[\kappa(L/2+z)\right]}{\sinh(\kappa L)}$",
                 ha="center", va="top", fontsize=fs)
        put(fig, 0.5, -0.47, r"$\phi_q^{p}(\rho,z)=\sum_{n=0}^{\infty}\left[\phi^{A}_{q_n}+\phi^{B}_{q_n}\right],\qquad"
                 r"\kappa^{2}=\dfrac{2e^{2}n_0}{\varepsilon k_B T}$", ha="center", va="top", fontsize=fs)
        put(fig, 0.5, -0.60, r"$\phi^{A}_{q_n}=\dfrac{q}{\pi L\varepsilon}\,\sin\left(\dfrac{2\pi n}{L}d\right)"
                 r"\sin\left(\dfrac{2\pi n}{L}z\right)K_0\left(\sqrt{\left(\dfrac{2\pi n}{L}\right)^{2}+\kappa^{2}}\,\rho\right)$",
                 ha="center", va="top", fontsize=fs)
        put(fig, 0.5, -0.73, r"$\phi^{B}_{q_n}=\dfrac{q}{\pi L\varepsilon}\,\cos\left(\dfrac{(2n+1)\pi}{L}d\right)"
                 r"\cos\left(\dfrac{(2n+1)\pi}{L}z\right)K_0\left(\sqrt{\left(\dfrac{(2n+1)\pi}{L}\right)^{2}+\kappa^{2}}\,\rho\right)$",
                 ha="center", va="top", fontsize=fs)
        put(fig, 0.5, -0.87, r"$\rho_{ij}=\sqrt{(x-i\ell)^{2}+(y-j\ell)^{2}}$, $\ell=32$. Solid: $|i|,|j|\leq 10$ (lateral images); "
                 r"dashed: only $(i,j)=(0,0)$, i.e. $\rho_{00}=\sqrt{x^{2}+y^{2}}$ (infinite plates)",
                 ha="center", va="top", fontsize=fs - 0.5)
        put(fig, 0.5, -0.93, r"Plates at $z=\mp L/2$ ($L=32$), $\psi_1=\psi_2=0$; $\lambda_D=1/\kappa=6$; $d=-8$ (charge position); "
                 r"point charge, $q/\varepsilon=10$. Field points on the line $y=0$, $x=\rho$",
                 ha="center", va="top", fontsize=fs - 0.5)
    fig.savefig(os.path.join(ROOT, "graficos", out_name(name)), dpi=170, bbox_inches="tight")

if __name__ == "__main__":
    make("a"); make("b")
