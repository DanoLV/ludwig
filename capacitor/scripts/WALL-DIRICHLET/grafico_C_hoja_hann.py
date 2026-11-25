#!/usr/bin/env python3
"""Electrolyte capacitor (closed system): the solution of the sheet (grafico_C_hoja_gauss.py) with the charge spread by
Ludwig's Hann-8 kernel instead of a Gaussian, vs Ludwig and vs the Gaussian solution.

Charge density (continuous version of Ludwig's spread, subgrid.c d_hann):
  rho_part = q W(x) W(y) W(z - z0),   W(u) = (1/n)[1 + cos(2 pi u/n)] for |u| <= n/2, 0 otherwise,  n = 8
Same derivation as the sheet; only the source transform changes:
  phi   = phi_C + U + [kL/(2 tanh(kL/2))] <U> [1 - cosh(k z)/cosh(k L/2)]    (phi_C: closed-system capacitor, as there)
  U     = (2 pi/A) sum_nm psit(kx, ky, z) exp(i k.r)
  psit  = q What(kx) What(ky) / (2 pi eps a sinh(aL)) [sinh(a(L/2-z)) I_< + sinh(a(L/2+z)) I_>],  a^2 = kx^2+ky^2+kappa^2
  What(k) = [sin(x)/x] / [1 - (x/pi)^2],  x = k n/2       (2D transform is NOT radial: depends on kx, ky separately)
  I_<   = int_{-L/2}^{z} W(s-z0) sinh(a(s+L/2)) ds,   I_> = int_{z}^{L/2} W(s-z0) sinh(a(L/2-s)) ds
  int e^{b s} [1 + cos(w(s-z0))] ds = e^{b s} [1/b + (b cos(w(s-z0)) + w sin(w(s-z0)))/(b^2 + w^2)],   w = 2 pi/n
The z support [z0-4, z0+4] lies inside the channel, so C1 = 1 (no truncation).
What(k): transform of the raised cosine, sinc(x) + [sinc(x-pi) + sinc(x+pi)]/2 (elementary; value 1/2 at x = pi).
The sinh products are evaluated as sums of exponentials with non-positive exponents (no overflow for large a).
Checks printed: psit against direct quadrature, and <U> against the closed form
  <U> = (q/eps)/(kappa^2 A L) [1 - int W(s-z0) cosh(kappa s) ds / cosh(kappa L/2)].

Reduced units, q/eps = 10, L = l = 32, kappa = 1/6, z0 = -8 (particle at z = 9.5, x = y = 16 in Ludwig).
Usage: python3 grafico_C_hoja_hann.py ;  LAYOUT=left python3 grafico_C_hoja_hann.py  (slide version, *_ppt.png)
"""
import os, math
import numpy as np
from scipy.integrate import quad
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import grafico_C_formulas as F
import grafico_C_hoja_gauss as H

L, A, FOUR = F.L, F.A, F.FOUR
KAP = 1 / 6.0; Z0 = F.D; NH = 8.0; OM = 2 * math.pi / NH
NK = 96                                                 # lateral modes |n|,|m| <= NK (What ~ 1/k^3: slower than the Gaussian)
ZA, ZB = max(-L / 2, Z0 - NH / 2), min(L / 2, Z0 + NH / 2)

def W(u):
    return np.where(np.abs(u) < NH / 2, (1 + np.cos(OM * u)) / NH, 0.0)

def What(k):
    x = np.abs(k) * NH / 2
    with np.errstate(invalid="ignore", divide="ignore"):
        v = np.sinc(x / np.pi) / (1 - (x / np.pi) ** 2)
    return np.where(np.abs(x - np.pi) < 1e-9, 0.5, v)

def _J(b, c, s1, s2):
    """int_{s1}^{s2} e^{b(s-c)} W(s-z0) ds (s1 <= s2 inside the support), closed form; arrays broadcast."""
    def P(s):
        u = s - Z0
        return np.exp(b * (s - c)) / NH * (1 / b + (b * np.cos(OM * u) + OM * np.sin(OM * u)) / (b * b + OM * OM))
    return np.where(s2 > s1, P(s2) - P(np.minimum(s1, s2)), 0.0)

def psit(kx, ky, z):
    """psit(kx, ky, z): kx, ky arrays (modes), z array -> shape (modes, z)."""
    a = np.sqrt(kx ** 2 + ky ** 2 + KAP ** 2)[:, None]; zz = z[None, :]
    lo, hi = ZA, np.clip(zz, ZA, ZB)                    # s < z part
    S = (_J(a, zz, lo, hi) - _J(-a, -(L + zz), lo, hi) - _J(a, L - zz, lo, hi) + _J(-a, -(2 * L - zz), lo, hi))
    lo, hi = np.clip(zz, ZA, ZB), ZB                    # s > z part
    S += (_J(-a, zz, lo, hi) - _J(-a, -(L + zz), lo, hi) - _J(a, L - zz, lo, hi) + _J(a, 2 * L + zz, lo, hi))
    # sinh(a(L/2-z_>)) sinh(a(z_<+L/2)) / sinh(aL) = e^{-aL} (sum of the 4 exponentials) / (2 (1 - e^{-2aL}))
    G = S / (2 * a * (1 - np.exp(-2 * a * L)))
    return FOUR / (2 * math.pi) * (What(kx) * What(ky))[:, None] * G

def dpsit(kx, ky, z):
    """d psit/dz, closed form (2026-11-06): derivative of each exponential of psit under the s integral; the terms from
    the moving limit s = z cancel (the Green's function is continuous at s = z). Shape (modes, z)."""
    a = np.sqrt(kx ** 2 + ky ** 2 + KAP ** 2)[:, None]; zz = z[None, :]
    lo, hi = ZA, np.clip(zz, ZA, ZB)                    # s < z part
    S = a * (-_J(a, zz, lo, hi) + _J(-a, -(L + zz), lo, hi) - _J(a, L - zz, lo, hi) + _J(-a, -(2 * L - zz), lo, hi))
    lo, hi = np.clip(zz, ZA, ZB), ZB                    # s > z part
    S += a * (_J(-a, zz, lo, hi) + _J(-a, -(L + zz), lo, hi) - _J(a, L - zz, lo, hi) - _J(a, 2 * L + zz, lo, hi))
    G = S / (2 * a * (1 - np.exp(-2 * a * L)))
    return FOUR / (2 * math.pi) * (What(kx) * What(ky))[:, None] * G

def U(rho, z):
    """(2 pi/A) sum_nm psit e^{i k.r} on the line y = 0, x = rho (charge at x = y = 0): sum of cos(kx rho)."""
    n, m = np.meshgrid(np.arange(-NK, NK + 1), np.arange(-NK, NK + 1), indexing="ij")
    kx = 2 * np.pi * n.ravel() / L; ky = 2 * np.pi * m.ravel() / L
    out = np.zeros(len(z))
    for i0 in range(0, len(kx), 4000):                  # in blocks (memory)
        sl = slice(i0, i0 + 4000)
        out += (np.cos(kx[sl] * rho)[:, None] * psit(kx[sl], ky[sl], z)).sum(axis=0)
    return 2 * np.pi / A * out

def dU(rho, z):
    """dU/dz on the line y = 0, x = rho: as U with dpsit."""
    n, m = np.meshgrid(np.arange(-NK, NK + 1), np.arange(-NK, NK + 1), indexing="ij")
    kx = 2 * np.pi * n.ravel() / L; ky = 2 * np.pi * m.ravel() / L
    out = np.zeros(len(z))
    for i0 in range(0, len(kx), 4000):
        sl = slice(i0, i0 + 4000)
        out += (np.cos(kx[sl] * rho)[:, None] * dpsit(kx[sl], ky[sl], z)).sum(axis=0)
    return 2 * np.pi / A * out

def mean_U():
    x, w = np.polynomial.legendre.leggauss(400); tot = 0.0
    for a, b in ((-L / 2, ZA), (ZA, Z0), (Z0, ZB), (ZB, L / 2)):
        z = 0.5 * (b - a) * x + 0.5 * (b + a); tot += 0.5 * (b - a) * np.sum(w * psit(np.zeros(1), np.zeros(1), z)[0])
    return 2 * np.pi / (A * L) * tot

def phi(rho, z, p1, p2, mu):
    fac = (KAP * L / 2) / math.tanh(KAP * L / 2)
    return H.cap_closed(z, p1, p2) + U(rho, z) + fac * mu * H.h(z)

def dphi(rho, z, p1, p2, mu):
    """d phi/dz, closed form: phi_C' + dU/dz + fac <U> h'(z)."""
    fac = (KAP * L / 2) / math.tanh(KAP * L / 2)
    dcap = 0.5 * (p2 - p1) * KAP * np.cosh(KAP * z) / np.sinh(KAP * L / 2)
    dh = -KAP * np.sinh(KAP * z) / np.cosh(KAP * L / 2)
    return dcap + dU(rho, z) + fac * mu * dh

def checks(mu):
    for kx, ky, z in ((0.0, 0.0, -8.0), (0.6, 0.2, -6.3), (2.4, 1.1, -9.0), (0.4, 0.0, 5.0)):
        a = math.sqrt(kx * kx + ky * ky + KAP ** 2)
        G = lambda s: (np.sinh(a * (min(z, s) + L / 2)) * np.sinh(a * (L / 2 - max(z, s))) / (a * np.sinh(a * L)))
        qd = quad(lambda s: float(W(s - Z0)) * G(s), ZA, ZB, points=[z, Z0], limit=200)[0]
        qd *= FOUR / (2 * math.pi) * float(What(kx) * What(ky))
        print("psit(%.1f,%.1f,z=%5.1f): closed form %.12e | quadrature %.12e"
              % (kx, ky, z, psit(np.array([kx]), np.array([ky]), np.array([z]))[0, 0], qd))
    ic = quad(lambda s: float(W(s - Z0)) * math.cosh(KAP * s), ZA, ZB)[0]
    print("<U>: quadrature of psit(0,z) %.10e | closed form %.10e"
          % (mu, FOUR / (KAP ** 2 * A * L) * (1 - ic / math.cosh(KAP * L / 2))))

def figure(bat, p1, p2, name, mu, mu_g):
    F.style(); grounded = (p1 == 0 and p2 == 0)
    if F.LAYOUT == "left":
        fig = plt.figure(figsize=(19.5, 9.6)); gs = fig.add_gridspec(1, 2, width_ratios=[1.65, 1], wspace=0.06)
        F._AXT[0] = fig.add_subplot(gs[0]); F._AXT[0].axis("off"); F._YMAX[0] = 1.62
        a = fig.add_subplot(gs[1])
    else:
        # 2026-11-16 (user's request): only the potential panel, the formulas and the parameters (the residual panel
        # "Ludwig - formula" and the note on the derivation were removed; the residuals are still printed)
        fig = plt.figure(figsize=(7.2, 4.6)); a = fig.add_subplot(111)
    zz = np.linspace(-L / 2, L / 2, 321)
    print("%s (psi_1 = %g, psi_2 = %g)" % (bat, p1, p2))
    print("  rho | max|Ludwig - Hann| out / in kernel | max|Ludwig - Gauss| out / in | max|Hann - Gauss|")
    for dx, c in F.COLS.items():
        z, ps = F.ludwig(os.path.join(F.ROOT, bat), dx)
        fh = lambda zv: phi(float(dx), zv, p1, p2, mu); fg = lambda zv: H.phi(float(dx), zv, p1, p2, mu_g)
        a.plot(zz, fh(zz), "-", color=c, lw=1.6, alpha=0.75, zorder=2)
        a.plot(zz, fg(zz), "--", color=c, lw=1.1, dashes=(4, 3), zorder=3)
        a.plot(z, ps, "o", ms=4.6, color=c, mec="white", mew=0.5, zorder=4, label="Ludwig, $\\rho=%d$" % dx)
        dh, dg = ps - fh(z), ps - fg(z); out = np.abs(z - Z0) > NH / 2
        print("   %2d |  %.2e / %.2e |  %.2e / %.2e |  %.2e"
              % (dx, np.abs(dh[out]).max(), np.abs(dh[~out]).max(), np.abs(dg[out]).max(), np.abs(dg[~out]).max(),
                 np.abs(fh(zz) - fg(zz)).max()))
    a.plot([], [], "-", color="#55544f", lw=1.6, alpha=0.75, label="derived solution, Hann-8 charge")
    a.plot([], [], "--", color="#55544f", lw=1.1, dashes=(4, 3), label="derived solution, Gaussian charge")
    a.set_ylabel(r"reduced potential $\psi = e\phi/k_BT$"); a.legend(loc="upper right", fontsize=9)
    a.set_title("1:1 electrolyte, $\\lambda_D=6$, " + ("grounded plates" if grounded else
                "plates at $\\psi=%g$ and $%g$" % (p1, p2)) + " (d = 8 from the lower plate)",
                fontsize=(10 if F.LAYOUT != "left" else None))
    a.set_xlabel("z (distance from the mid-plane, lattice units)")
    if F.LAYOUT != "left": fig.subplots_adjust(left=0.13, right=0.97, top=0.92, bottom=0.11)
    if F.LAYOUT != "left":
        # 2026-11-16 (user's request): parameters ABOVE the formulas and no line wider than the plot. One equation per
        # line, centred on the axes; vertical positions in inches below the figure, converted to figure fraction.
        H_IN = fig.get_size_inches()[1]; xc = 0.55; y_in = [0.42]
        def line(s, h, **kw):
            fig.text(xc, -y_in[0] / H_IN, s, ha="center", va="top", **kw); y_in[0] += h
        sm, eq = dict(fontsize=9.3, color="#333333"), dict(fontsize=10.5)
        line(r"$\bf{Parameters}$", 0.24, fontsize=9.8)
        line(r"plates at $z=\mp L/2$, $L=\ell=32$;  $\psi_1=%g$, $\psi_2=%g$;  $\lambda_D=1/\kappa=6$" % (p1, p2), 0.22, **sm)
        line(r"charge at $z_0=-8$ (d = 8 from the lower plate);  Hann kernel $n=8$;  $q/\varepsilon=10$", 0.22, **sm)
        line(r"field points on the line $y=0$, $x=\rho$;  dashed: Gaussian, same variance ($\sigma=1.446$)", 0.40, **sm)
        line(r"$\bf{Linearised\ Poisson\text{-}Boltzmann,\ closed\ system,\ Hann\text{-}8\ charge}$", 0.34, fontsize=9.8)
        line(r"$\left(\nabla^{2}-\kappa^{2}\right)\phi=-\kappa^{2}\langle\phi\rangle-\dfrac{\rho_{\mathrm{part}}}{\varepsilon},\qquad"
             r"\rho_{\mathrm{part}}=q\,W(x)\,W(y)\,W(z-z_0)$", 0.50, **eq)
        line(r"$W(u)=\dfrac{1}{n}\left[1+\cos\left(\dfrac{2\pi u}{n}\right)\right],\ |u|\leq\dfrac{n}{2},\qquad"
             r"\hat W(k)=\dfrac{\sin(kn/2)}{kn/2}\,\dfrac{1}{1-\left(kn/2\pi\right)^{2}}$", 0.52, **eq)
        line(r"(Gaussian: $\hat W(k)=e^{-\sigma^{2}k^{2}/2}$)", 0.32, **sm)
        line(r"$\phi=\phi_{C}(z)+U(x,y,z)+\dfrac{\kappa L}{2\tanh(\kappa L/2)}\,\langle U\rangle"
             r"\left[1-\dfrac{\cosh(\kappa z)}{\cosh(\kappa L/2)}\right]$", 0.55, **eq)
        line(r"$U=\dfrac{2\pi}{A}\sum_{n,m}\tilde\psi_\infty(\mathbf{k}_{nm},z)\,e^{i\mathbf{k}_{nm}\cdot\mathbf{r}},\qquad"
             r"\langle U\rangle=\dfrac{2\pi}{AL}\int_{-L/2}^{L/2}\tilde\psi_\infty(0,z)\,dz$", 0.58, **eq)
        line(r"$\tilde\psi_\infty(\mathbf{k},z)=\dfrac{q\,\hat W(k_x)\,\hat W(k_y)}{2\pi\varepsilon\,\alpha\sinh(\alpha L)}"
             r"\left[\sinh\left(\alpha(L/2-z)\right)I_<+\sinh\left(\alpha(L/2+z)\right)I_>\right]$", 0.55, **eq)
        line(r"$\alpha^{2}=k_x^{2}+k_y^{2}+\kappa^{2},\qquad I_<=\int_{-L/2}^{z}W(s-z_0)\sinh\left(\alpha(s+L/2)\right)ds$", 0.48, **eq)
        line(r"$I_>=\int_{z}^{L/2}W(s-z_0)\sinh\left(\alpha(L/2-s)\right)ds$", 0.48, **eq)
        line(r"$\int e^{\beta s}\left[1+\cos\omega(s-z_0)\right]ds=e^{\beta s}\left[\dfrac{1}{\beta}"
             r"+\dfrac{\beta\cos\omega(s-z_0)+\omega\sin\omega(s-z_0)}{\beta^{2}+\omega^{2}}\right]$", 0.55, **eq)
        line(r"$\beta=\pm\alpha,\quad\omega=\dfrac{2\pi}{n}$", 0.48, **eq)
        line(r"$\phi_{C}(z)=V_m+\dfrac{\Delta V}{2\sinh(\kappa L/2)}\,\sinh(\kappa z),\qquad"
             r"V_m=\dfrac{\phi_1+\phi_2}{2},\quad\Delta V=\phi_2-\phi_1$", 0.55, **eq)
        line(r"$\phi_C''=\kappa^{2}\left(\phi_C-\langle\phi_C\rangle\right),\quad\langle\phi_C\rangle=V_m$", 0.3, **eq)
        fig.savefig(os.path.join(F.ROOT, "graficos", F.out_name(name)), dpi=170, bbox_inches="tight")
        return
    fs = 10.5; put = lambda y, s, **kw: F.put(fig, 0.5, y, s, ha="center", va="top", **kw)
    if F.LAYOUT != "left": put = lambda y, s, **kw: F.put(fig, 0.5, y * 0.78 * 6.4 / 4.6 - 0.03, s, ha="center", va="top", **kw)
    put(-0.02, r"Linearised Poisson-Boltzmann, closed system, charge spread with Ludwig's Hann-8 kernel",
        fontsize=fs, fontweight="bold")
    put(-0.07, r"$\left(\nabla^{2}-\kappa^{2}\right)\phi=-\kappa^{2}\langle\phi\rangle-\dfrac{\rho_{\mathrm{part}}}{\varepsilon},\qquad"
        r"\rho_{\mathrm{part}}=q\,W(x)\,W(y)\,W(z-z_0),\qquad W(u)=\dfrac{1}{n}\left[1+\cos\left(\dfrac{2\pi u}{n}\right)\right],\ |u|\leq\dfrac{n}{2}$",
        fontsize=fs)
    put(-0.19, r"$\phi=\phi_{C}(z)+U(x,y,z)+\dfrac{\kappa L}{2\tanh(\kappa L/2)}\,\langle U\rangle"
        r"\left[1-\dfrac{\cosh(\kappa z)}{\cosh(\kappa L/2)}\right],\qquad"
        r"U=\dfrac{2\pi}{A}\sum_{n,m}\tilde\psi_\infty(\mathbf{k}_{nm},z)\,e^{i\mathbf{k}_{nm}\cdot\mathbf{r}}$", fontsize=fs)
    put(-0.31, r"$\tilde\psi_\infty(\mathbf{k},z)=\dfrac{q\,\hat W(k_x)\,\hat W(k_y)}{2\pi\varepsilon\,\alpha\sinh(\alpha L)}"
        r"\left[\sinh\left(\alpha(L/2-z)\right)I_<+\sinh\left(\alpha(L/2+z)\right)I_>\right],\qquad"
        r"\alpha^{2}=k_x^{2}+k_y^{2}+\kappa^{2}$", fontsize=fs)
    put(-0.43, r"$\hat W(k)=\dfrac{\sin(kn/2)}{kn/2}\,\dfrac{1}{1-\left(kn/2\pi\right)^{2}}\quad"
        r"\left(\mathrm{Gaussian:}\ e^{-\sigma^{2}k^{2}/2}\right),\qquad"
        r"I_<=\int_{-L/2}^{z}W(s-z_0)\sinh\left(\alpha(s+L/2)\right)ds,\qquad"
        r"I_>=\int_{z}^{L/2}W(s-z_0)\sinh\left(\alpha(L/2-s)\right)ds$", fontsize=fs)
    put(-0.55, r"$\int e^{\beta s}\left[1+\cos\omega(s-z_0)\right]ds=e^{\beta s}\left[\dfrac{1}{\beta}"
        r"+\dfrac{\beta\cos\omega(s-z_0)+\omega\sin\omega(s-z_0)}{\beta^{2}+\omega^{2}}\right],\qquad"
        r"\beta=\pm\alpha,\quad\omega=\dfrac{2\pi}{n}$", fontsize=fs)
    put(-0.67, H.CAP_TEX, fontsize=fs)
    put(-0.79, r"$\langle U\rangle=\dfrac{2\pi}{AL}\int_{-L/2}^{L/2}\tilde\psi_\infty(0,z)\,dz$", fontsize=fs)
    put(-0.91 if F.LAYOUT != "left" else -1.07,
        (r"Plates at $z=\mp L/2$ ($L=\ell=32$), $\psi_1=%g$, $\psi_2=%g$; $\lambda_D=1/\kappa=6$; $z_0=-8$; " % (p1, p2))
        + r"$n=8$; $q/\varepsilon=10$. Field points on the line $y=0$, $x=\rho$. Dashed: Gaussian with the same variance ($\sigma=1.446$)",
        fontsize=fs - 0.5)
    fig.savefig(os.path.join(F.ROOT, "graficos", F.out_name(name)), dpi=170, bbox_inches="tight")

if __name__ == "__main__":
    mu = mean_U(); mu_g = H.mean_U()
    checks(mu)
    print("<U>: Hann %.10e | Gaussian %.10e" % (mu, mu_g))
    figure("C-electrolito-hann8-LD6", 0.0, 0.0, "C_hoja_hann_closed_vs_ludwig.png", mu, mu_g)
    figure("C-electrolito-dospotenciales-hann8-LD6", 0.1, 0.3, "C_hoja_hann_closed_2pot_vs_ludwig.png", mu, mu_g)
