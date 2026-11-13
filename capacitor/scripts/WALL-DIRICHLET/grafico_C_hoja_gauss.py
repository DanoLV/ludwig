#!/usr/bin/env python3
"""Electrolyte capacitor (closed system) with a GAUSSIAN charge: the hand-derived solution of the sheet (2026-10-30),
evaluated term by term, vs Ludwig.

  phi   = phi_C(z) + U(x,y,z) + [kL/(2 tanh(kL/2))] <U> [1 - cosh(k z)/cosh(k L/2)],   k = kappa
  phi_C = V_m + (DV/2) sinh(k z)/sinh(k L/2),  V_m = (p1+p2)/2, DV = p2 - p1   (capacitor of the CLOSED system, see below)
  U     = (2 pi/A) sum_nm  psit(k_nm, z) exp(i k_nm . r),     k_nm = (2 pi/l)(n, m)      (exact periodic images)
  <U>   = (2 pi/(A L)) int_{-L/2}^{L/2} psit(0, z) dz                                    (k = 0 mode, by quadrature)
  psit  = q exp(-s^2 k^2/2) / (2 pi eps a sinh(aL)) [sinh(a(L/2 - z)) I_< + sinh(a(L/2 + z)) I_>],   a^2 = k^2 + kappa^2
  I_<   = (1/C1) int_{-L/2}^{z}  exp(-(s'-z0)^2/(2 s^2)) sinh(a(s'+L/2)) ds'
  I_>   = (1/C1) int_{z}^{L/2}   exp(-(s'-z0)^2/(2 s^2)) sinh(a(L/2-s')) ds'
  C1    = int_{-L/2}^{L/2} exp(-(s'-z0)^2/(2 s^2)) ds'        (Gaussian truncated to the channel, renormalised)
  phi_CAP = [p1 sinh k(L/2-z) + p2 sinh k(L/2+z)]/sinh kL   (open system: phi'' - k^2 phi = 0)
phi_C solves phi'' - k^2 (phi - <phi>) = 0 with phi(-/+L/2) = p1, p2 and <phi_C> = V_m. It equals phi_CAP plus the
closed-system term of the sheet applied to phi_CAP, [kL/(2 tanh(kL/2))] <phi_CAP> [1 - cosh/cosh] = V_m [1 - cosh/cosh]:
splitting phi_CAP = V_m cosh(kz)/cosh(kL/2) + (DV/2) sinh(kz)/sinh(kL/2), the cosh parts cancel.
I_< and I_> are evaluated with the erf closed forms of the sheet, written in log/erfcx form to avoid overflow and
cancellation for large a (algebraically identical).  Green function: Jackson, Classical Electrodynamics 3rd ed.,
Problem 3.17 (with k^2 -> k^2 + kappa^2).  Two corrections with respect to the sheet as written:
  (1) sign: the sheet has sinh(a(z - L/2)) I_< + sinh(a(z + L/2)) I_> with its I_> = -(the I_> above), i.e. -psit;
  (2) phi_CAP of the sheet solves the open-system equation; the closed-system capacitor is phi_C
      (with p1 = p2 = V it gives phi_C = V, no double layer: the ions only shift their reference).
Checks printed: the k = 0 mean against the closed form of grafico_C_electrolito_gauss.py, and the grounded
potential against analytic() of that script (sine series + Hankel + K0 point images).

Reduced units, q/eps = 10, L = l = 32, kappa = 1/6, z0 = -8, sigma = 1.446 (Hann-8 width).
Usage: python3 grafico_C_hoja_gauss.py            (grounded + two-potential figures)
       LAYOUT=left python3 grafico_C_hoja_gauss.py  (slide version, formulas at the left -> *_ppt.png)
"""
import os, math
import numpy as np
from scipy.special import erf, erfcx                   # erf and scaled erfc (scipy.special, Cephes)
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import grafico_C_formulas as F
import grafico_C_electrolito_gauss as G

L, A, FOUR = F.L, F.A, F.FOUR                          # FOUR = q/eps
KAP = 1 / 6.0; SIG = F.SIG; Z0 = F.D; NK = 24          # lateral modes |n|,|m| <= NK (exp(-s^2 k^2/2) < 1e-20)
SQ = math.sqrt(2) * SIG

def c1():
    return SIG * math.sqrt(math.pi / 2) * (erf((L / 2 - Z0) / SQ) - erf((-L / 2 - Z0) / SQ))

def _eerf(lnp, v, u):
    """exp(lnp) [erf(v) - erf(u)], stable when both arguments have the same sign (erfc = exp(-x^2) erfcx)."""
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        pos = np.exp(lnp - u * u) * erfcx(u) - np.exp(lnp - v * v) * erfcx(v)
        neg = np.exp(lnp - v * v) * erfcx(-v) - np.exp(lnp - u * u) * erfcx(-u)
        mix = np.exp(lnp) * (erf(v) - erf(u))
    return np.where(u >= 0, pos, np.where(v <= 0, neg, mix))

def _lnsinh(x):
    with np.errstate(divide="ignore"):
        return np.where(x > 0, x + np.log1p(-np.exp(-2 * np.maximum(x, 1e-300))) - math.log(2), -np.inf)

def psit(k, z):
    """psit(k, z) of the sheet (corrected sign), for one |k| and an array of z."""
    a = math.sqrt(k * k + KAP * KAP); s2 = a * SIG ** 2
    lr = _lnsinh(a * (L / 2 - z)) - _lnsinh(np.array(a * L))      # ln[sinh(a(L/2-z))/sinh(aL)]
    ll = _lnsinh(a * (L / 2 + z)) - _lnsinh(np.array(a * L))      # ln[sinh(a(L/2+z))/sinh(aL)]
    # I_< : e^{a(z0+L/2)}{erf((z-z0-a s^2)/sq) - erf((-L/2-z0-a s^2)/sq)} - e^{-a(z0+L/2)}{same with +a s^2}
    ilt = (_eerf(lr + a * (Z0 + L / 2), (z - Z0 - s2) / SQ, (-L / 2 - Z0 - s2) / SQ)
           - _eerf(lr - a * (Z0 + L / 2), (z - Z0 + s2) / SQ, (-L / 2 - Z0 + s2) / SQ))
    # I_> : e^{-a(z0-L/2)}{erf((L/2-z0+a s^2)/sq) - erf((z-z0+a s^2)/sq)} - e^{a(z0-L/2)}{same with -a s^2}
    igt = (_eerf(ll - a * (Z0 - L / 2), (L / 2 - Z0 + s2) / SQ, (z - Z0 + s2) / SQ)
           - _eerf(ll + a * (Z0 - L / 2), (L / 2 - Z0 - s2) / SQ, (z - Z0 - s2) / SQ))
    # q e^{-s^2k^2/2}/(2 pi eps a) * sqrt(pi/2) s e^{(a s)^2/2}/(2 C1);  e^{-s^2k^2/2} e^{a^2 s^2/2} = e^{kappa^2 s^2/2}
    pre = FOUR / (2 * math.pi * a) * math.sqrt(math.pi / 2) * SIG * math.exp(KAP ** 2 * SIG ** 2 / 2) / (2 * c1())
    return pre * (ilt + igt)

def U(rho, z):
    """(2 pi/A) sum_nm psit(k_nm, z) cos(k_x rho): line y = 0, x = rho; charge at x = y = 0."""
    n, m = np.meshgrid(np.arange(-NK, NK + 1), np.arange(-NK, NK + 1), indexing="ij")
    kx = 2 * np.pi * n.ravel() / L; k = 2 * np.pi * np.hypot(n, m).ravel() / L
    uk, inv = np.unique(np.round(k, 12), return_inverse=True)
    out = np.zeros(len(z))
    for i, kk in enumerate(uk):
        out += np.cos(kx[inv == i] * rho).sum() * psit(kk, z)
    return 2 * np.pi / A * out

def mean_U():
    """<U> = (2 pi/(A L)) int psit(0, z) dz  (Gauss-Legendre, split at z0)."""
    x, w = np.polynomial.legendre.leggauss(400); tot = 0.0
    for a, b in ((-L / 2, Z0), (Z0, L / 2)):
        z = 0.5 * (b - a) * x + 0.5 * (b + a); tot += 0.5 * (b - a) * np.sum(w * psit(0.0, z))
    return 2 * np.pi / (A * L) * tot

def mean_cap(p1, p2):
    return (p1 + p2) * math.tanh(KAP * L / 2) / (KAP * L)

def h(z):
    return 1 - np.cosh(KAP * z) / np.cosh(KAP * L / 2)

def cap_closed(z, p1, p2):
    """Capacitor of the closed system: V_m + (DV/2) sinh(kz)/sinh(kL/2), V_m = (p1+p2)/2, DV = p2 - p1."""
    return 0.5 * (p1 + p2) + 0.5 * (p2 - p1) * np.sinh(KAP * z) / np.sinh(KAP * L / 2)

def phi(rho, z, p1, p2, mu):
    fac = (KAP * L / 2) / math.tanh(KAP * L / 2)
    return cap_closed(z, p1, p2) + U(rho, z) + fac * mu * h(z)

def check_cap(p1=0.1, p2=0.3):
    """phi_C: residual of phi'' - k^2 (phi - <phi>) = 0, mean and wall values (and the same for the open phi_CAP)."""
    x, w = np.polynomial.legendre.leggauss(200); zq = L / 2 * x
    zz = np.linspace(-L / 2 + 0.5, L / 2 - 0.5, 200); e = 1e-3
    for name, f in (("open   phi_CAP", lambda z: F.phi_cap_elec(z, p1, p2, KAP)), ("closed phi_C  ", lambda z: cap_closed(z, p1, p2))):
        m = 0.5 * np.sum(w * f(zq))
        res = (f(zz + e) - 2 * f(zz) + f(zz - e)) / e ** 2 - KAP ** 2 * (f(zz) - m)
        print("%s (p1 = %g, p2 = %g): <phi> = %.6f, walls = %.4f %.4f, max|phi'' - k^2(phi - <phi>)| = %.1e"
              % (name, p1, p2, m, f(np.array([-L / 2]))[0], f(np.array([L / 2]))[0], np.abs(res).max()))
    d = np.abs(cap_closed(zz, p1, p2) - F.phi_cap_elec(zz, p1, p2, KAP) - 0.5 * (p1 + p2) * h(zz)).max()
    print("phi_C vs phi_CAP + V_m [1 - cosh/cosh]: max diff %.1e" % d)

def checks(mu):
    G.TRUNC = True
    print("C1/(sqrt(2 pi) sigma) = %.12f" % (c1() / (SIG * math.sqrt(2 * math.pi))))
    print("<U>: quadrature of psit(0,z) %.10e | closed form (grafico_C_electrolito_gauss) %.10e" % (mu, G.mean_q()))
    zz = np.linspace(-L / 2 + 0.2, L / 2 - 0.2, 61)
    for rho in (0, 6, 10):
        d = np.abs(phi(float(rho), zz, 0, 0, mu) - G.analytic(float(rho), zz)[0]).max()
        print("grounded, rho = %2d: max|sheet - analytic() of grafico_C_electrolito_gauss| = %.2e" % (rho, d))

CAP_TEX = (r"$\phi_{C}(z)=V_m+\dfrac{\Delta V}{2\sinh(\kappa L/2)}\,\sinh(\kappa z),\qquad"
           r"V_m=\dfrac{\phi_1+\phi_2}{2},\quad\Delta V=\phi_2-\phi_1,\qquad"
           r"\phi_C''=\kappa^{2}\left(\phi_C-\langle\phi_C\rangle\right),\quad\langle\phi_C\rangle=V_m$")

def figure(bat, p1, p2, name, mu):
    F.style()
    grounded = (p1 == 0 and p2 == 0)
    if F.LAYOUT == "left":                    # wider formula panel than F.new_fig (long I_< and I_> lines)
        fig = plt.figure(figsize=(19.5, 9.6)); gs = fig.add_gridspec(1, 2, width_ratios=[1.65, 1], wspace=0.04)
        F._AXT[0] = fig.add_subplot(gs[0]); F._AXT[0].axis("off"); a = fig.add_subplot(gs[1]); F._YMAX[0] = 1.62
    else:
        fig, a = F.new_fig(1.62, 9.6)
    zz = np.linspace(-L / 2, L / 2, 241)
    print("%s (psi_1 = %g, psi_2 = %g):  rho | max|Ludwig - sheet| outside kernel (|z-z0|>4.5) | inside | max(Ludwig)"
          % (bat, p1, p2))
    for dx, c in F.COLS.items():
        z, ps = F.ludwig(os.path.join(F.ROOT, bat), dx)
        a.plot(zz, phi(float(dx), zz, p1, p2, mu), "-", color=c, lw=1.6, alpha=0.75, zorder=2)
        a.plot(z, ps, "o", ms=4.6, color=c, mec="white", mew=0.5, zorder=4, label="Ludwig, $\\rho=%d$" % dx)
        d = np.abs(ps - phi(float(dx), z, p1, p2, mu)); out = np.abs(z - Z0) > 4.5
        print("                 %2d | %.2e | %.2e | %.4f" % (dx, d[out].max(), d[~out].max(), ps.max()))
    a.plot([], [], "-", color="#55544f", lw=1.6, alpha=0.75, label="derived solution (Gaussian charge, closed system)")
    a.set_xlabel("z (distance from the mid-plane, lattice units)"); a.set_ylabel(r"reduced potential $\psi = e\phi/k_BT$")
    a.set_title("1:1 electrolyte, $\\lambda_D=6$, " + ("grounded plates" if grounded else
                "plates at $\\psi=%g$ and $%g$" % (p1, p2)) + " (d = 8 from the lower plate)")
    a.legend(loc="upper right", fontsize=9)
    fig.tight_layout()
    fs = 10.5; put = lambda y, s, **kw: F.put(fig, 0.5, y, s, ha="center", va="top", **kw)
    put(-0.02, r"Linearised Poisson-Boltzmann with ion-number conservation (closed system), Gaussian charge",
        fontsize=fs, fontweight="bold")
    put(-0.07, r"$\left(\nabla^{2}-\kappa^{2}\right)\phi=-\kappa^{2}\langle\phi\rangle-\dfrac{\rho_{\mathrm{part}}}{\varepsilon},\qquad"
        r"\rho_{\mathrm{part}}=\dfrac{q}{2\pi\sigma^{2}}\,e^{-(x^{2}+y^{2})/(2\sigma^{2})}\;\dfrac{e^{-(z-z_0)^{2}/(2\sigma^{2})}}{C_1},\qquad"
        r"\phi(z=\mp L/2)=\phi_{1,2}$", fontsize=fs)
    put(-0.19, r"$\phi=\phi_{C}(z)+U(x,y,z)+\dfrac{\kappa L}{2\tanh(\kappa L/2)}\,\langle U\rangle"
        r"\left[1-\dfrac{\cosh(\kappa z)}{\cosh(\kappa L/2)}\right]$", fontsize=fs)
    put(-0.31, r"$U=\dfrac{2\pi}{A}\sum_{n,m}\tilde\psi_\infty(k_{nm},z)\,e^{i\mathbf{k}_{nm}\cdot\mathbf{r}},\qquad"
        r"\langle U\rangle=\dfrac{2\pi}{A L}\int_{-L/2}^{L/2}\tilde\psi_\infty(0,z)\,dz,\qquad"
        r"\mathbf{k}_{nm}=\dfrac{2\pi}{\ell}(n,m),\quad A=\ell^{2}$", fontsize=fs)
    put(-0.43, r"$\tilde\psi_\infty(k,z)=\dfrac{q\,e^{-\sigma^{2}k^{2}/2}}{2\pi\varepsilon\,\alpha\sinh(\alpha L)}"
        r"\left[\sinh\left(\alpha\left(L/2-z\right)\right)I_<+\sinh\left(\alpha\left(L/2+z\right)\right)I_>\right],\qquad"
        r"\alpha^{2}=k^{2}+\kappa^{2}$", fontsize=fs)
    put(-0.555, r"$I_<=\sqrt{\dfrac{\pi}{2}}\,\dfrac{\sigma e^{(\alpha\sigma)^{2}/2}}{2C_1}\left[e^{\alpha(z_0+L/2)}"
        r"\left\{\mathrm{erf}\left(\dfrac{z-z_0-\alpha\sigma^{2}}{\sqrt{2}\,\sigma}\right)-\mathrm{erf}\left(\dfrac{-L/2-z_0-\alpha\sigma^{2}}{\sqrt{2}\,\sigma}\right)\right\}"
        r"-e^{-\alpha(z_0+L/2)}\left\{\alpha\sigma^{2}\rightarrow-\alpha\sigma^{2}\right\}\right]$", fontsize=fs)
    put(-0.68, r"$I_>=\sqrt{\dfrac{\pi}{2}}\,\dfrac{\sigma e^{(\alpha\sigma)^{2}/2}}{2C_1}\left[e^{-\alpha(z_0-L/2)}"
        r"\left\{\mathrm{erf}\left(\dfrac{L/2-z_0+\alpha\sigma^{2}}{\sqrt{2}\,\sigma}\right)-\mathrm{erf}\left(\dfrac{z-z_0+\alpha\sigma^{2}}{\sqrt{2}\,\sigma}\right)\right\}"
        r"-e^{\alpha(z_0-L/2)}\left\{\alpha\sigma^{2}\rightarrow-\alpha\sigma^{2}\right\}\right]$", fontsize=fs)
    put(-0.81, r"$C_1=\int_{-L/2}^{L/2}e^{-(s-z_0)^{2}/(2\sigma^{2})}\,ds,\qquad"
        r"\kappa^{2}=\dfrac{2e^{2}n_0}{\varepsilon k_BT}$", fontsize=fs)
    put(-0.92, CAP_TEX, fontsize=fs)
    put(-1.05, r"Green function in $z$: Jackson, Classical Electrodynamics 3rd ed., Problem 3.17 ($k^{2}\to k^{2}+\kappa^{2}$). "
        r"Sum over $|n|,|m|\leq %d$: exact periodic (Gaussian) images in $x,y$." % NK, fontsize=fs - 0.5)
    put(-1.22 if F.LAYOUT == "left" else -1.15, (r"Plates at $z=\mp L/2$ ($L=\ell=32$), $\psi_1=%g$, $\psi_2=%g$; $\lambda_D=1/\kappa=6$; $z_0=-8$; " % (p1, p2))
        + r"$\sigma=1.446$ (Hann-8 width); $q/\varepsilon=10$. Field points on the line $y=0$, $x=\rho$", fontsize=fs - 0.5)
    fig.savefig(os.path.join(F.ROOT, "graficos", F.out_name(name)), dpi=170, bbox_inches="tight")

if __name__ == "__main__":
    mu = mean_U()
    checks(mu); check_cap()
    figure("C-electrolito-hann8-LD6", 0.0, 0.0, "C_hoja_gauss_closed_vs_ludwig.png", mu)
    figure("C-electrolito-dospotenciales-hann8-LD6", 0.1, 0.3, "C_hoja_gauss_closed_2pot_vs_ludwig.png", mu)
