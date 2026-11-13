#!/usr/bin/env python3
"""Check of the derived solution (electrolyte, closed box) with a GAUSSIAN charge, against
 (1) an independent numerical solver (lateral Fourier modes + finite differences in z, <phi - phibar> = 0 imposed numerically)
 (2) the Ludwig potential (C-electrolito-hann8-LD6, Hann-8 kernel; same variance as the gaussian).

Derived solution (plates at z = -/+ L/2 grounded, charge at z = d from the centre, box A = 32^2, Gaussian of std sigma):
  phi = phi_q(rho,z) + phibar [1 - cosh(kappa z)/cosh(kappa L/2)]
  phi_q = (q/(pi L eps)) sum_m c_m xi_m(z) I_m(rho) + lateral images,
  I_m(rho) = int_0^inf k J0(k rho) exp(-k^2 sigma^2/2)/(k^2 + g_m^2) dk,  alpha_m = m pi/L,  g_m^2 = alpha_m^2 + kappa^2
  c_m = int_{-L/2}^{L/2} g(s) xi_m(s) ds   (projection of the z profile g on the plate eigenfunctions)
  phibar  = (kappa L/2) coth(kappa L/2) <phi_q>

z profile g (CHANGE 20261030): with TRUNC = True the Gaussian is truncated to [-L/2, L/2] and renormalised,
  g(s) = exp(-(s-d)^2/(2 sigma^2)) / C1,  C1 = int_{-L/2}^{L/2} exp(-(s-d)^2/(2 sigma^2)) ds,
so the charge between the plates is exactly q (Ludwig loses the part of the kernel that falls on wall nodes instead,
see README "Pendiente"). c_m by Gauss-Legendre quadrature. With TRUNC = False, c_m = exp(-alpha_m^2 sigma^2/2) xi_m(d)
(full Gaussian, the tails beyond the plates are dropped without renormalising) as before.
The mean uses the closed form (derived by hand, no reference: int G(z,s) dz = h(s) solves (d2/ds2 - kappa^2) h = -1 with
h(+-L/2) = 0, so h = [1 - cosh(kappa s)/cosh(kappa L/2)]/kappa^2):
  <phi_q> = (q/eps)/(kappa^2 A L) [1 - int g(s) cosh(kappa s) ds / cosh(kappa L/2)]
Images stay point charges laterally (L >= 6 sigma) with the same z profile g.
Reduced units q/eps = 10.   Usage: python3 grafico_C_electrolito_gauss.py
"""
import os, glob, math
import numpy as np
from scipy.special import j0 as _j0, k0 as _k0       # Bessel J0 and modified Bessel K0 (scipy.special, Cephes)
from scipy.special import erf
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
L = 32.0; A = L * L; FOUR = 10.0; KAP = 1 / 6.0; D = 8.0 - L / 2
SIG = 8 * math.sqrt(1 / 12 - 1 / (2 * math.pi ** 2)); NREP = 10; M = 120
TRUNC = True                                          # Gaussian truncated to the channel and renormalised

def xi(m, z): return np.sin(m * np.pi * (z + L / 2) / L)

def _gint(a, b, c):
    """int_a^b exp(-(s-c)^2/(2 SIG^2)) ds"""
    sq = math.sqrt(2) * SIG
    return SIG * math.sqrt(math.pi / 2) * (erf((b - c) / sq) - erf((a - c) / sq))

def c1():
    """Normalisation of the z profile: truncated integral (TRUNC) or the full Gaussian."""
    return _gint(-L / 2, L / 2, D) if TRUNC else SIG * math.sqrt(2 * math.pi)

def zcoef(ms):
    """c_m = int g(s) xi_m(s) ds over the channel."""
    if not TRUNC:
        return np.exp(-(ms * np.pi / L) ** 2 * SIG ** 2 / 2) * xi(ms, D)
    a, b = max(-L / 2, D - 12 * SIG), min(L / 2, D + 12 * SIG)
    x, w = np.polynomial.legendre.leggauss(1000)
    s = 0.5 * (b - a) * x + 0.5 * (b + a); w = 0.5 * (b - a) * w
    g = np.exp(-(s - D) ** 2 / (2 * SIG ** 2)) / c1()
    return (xi(ms[:, None], s[None, :]) * (g * w)[None, :]).sum(axis=1)

def mean_q():
    """<phi_q> = (q/eps)/(kappa^2 A L) [1 - int g cosh(kappa s) ds / cosh(kappa L/2)]  (closed form, see docstring)."""
    k, h = KAP, L / 2
    if TRUNC:
        ic = 0.5 * math.exp(k * k * SIG ** 2 / 2) * (math.exp(k * D) * _gint(-h, h, D + k * SIG ** 2)
                                                    + math.exp(-k * D) * _gint(-h, h, D - k * SIG ** 2)) / c1()
    else:
        ic = math.exp(k * k * SIG ** 2 / 2) * math.cosh(k * D)
    return FOUR / (k * k * A * L) * (1 - ic / math.cosh(k * h))

def I_m(rho, ms):
    kk, ww = np.polynomial.legendre.leggauss(500); kmax = 9.0 / SIG
    k = 0.5 * (kk + 1) * kmax; w = 0.5 * kmax * ww
    j0 = _j0(k * rho)
    g2 = (ms * np.pi / L) ** 2 + KAP ** 2
    return np.array([np.sum(w * k * j0 * np.exp(-k ** 2 * SIG ** 2 / 2) / (k ** 2 + gg)) for gg in g2])

def analytic(rho, z, with_images=True):
    ms = np.arange(1, M + 1); al = ms * np.pi / L
    c = zcoef(ms) * I_m(rho, ms)
    phi = FOUR / (np.pi * L) * (xi(ms[None, :], z[:, None]) * c[None, :]).sum(axis=1)
    if with_images:                                   # replicas are >= 22 away: point-charge K0 series, 60 modes
        ms2 = np.arange(1, 61); al2 = ms2 * np.pi / L; g2 = np.sqrt(al2 ** 2 + KAP ** 2)
        cz2 = zcoef(ms2)                              # same z profile as the central charge
        for i in range(-NREP, NREP + 1):
            for j in range(-NREP, NREP + 1):
                if i == 0 and j == 0: continue
                r = math.hypot(rho - i * L, j * L)
                k0 = _k0(g2 * r)
                cc = cz2 * k0
                phi += FOUR / (np.pi * L) * (xi(ms2[None, :], z[:, None]) * cc[None, :]).sum(axis=1)
    pbar = (KAP * L / 2) / np.tanh(KAP * L / 2) * mean_q()
    return phi + pbar * (1 - np.cosh(KAP * z) / np.cosh(KAP * L / 2)), pbar

def numeric(rho_lines, z_out, Nz=800, Mx=20):
    """Independent solver: periodic box, modes (mx,my), finite differences in z, <phi - phibar> = 0 imposed."""
    dz = L / (Nz + 1); zg = -L / 2 + dz * np.arange(1, Nz + 1)
    gz = np.exp(-(zg - D) ** 2 / (2 * SIG ** 2)) / c1()          # same z profile g as analytic()
    mx, my = np.meshgrid(np.arange(-Mx, Mx + 1), np.arange(-Mx, Mx + 1), indexing="ij")
    k2 = ((2 * np.pi / L) ** 2 * (mx ** 2 + my ** 2)).ravel()
    uk, inv = np.unique(np.round(k2, 12), return_inverse=True)
    # (d2/dz2 - k2 - kap2) chi = -S ,  S = (FOUR/A) exp(-k2 sig2/2) gz ; Thomas algorithm, vectorised over modes
    diag = -(2 / dz ** 2 + uk + KAP ** 2)[:, None] * np.ones((1, Nz)); off = 1 / dz ** 2
    rhs = -(FOUR / A) * np.exp(-uk * SIG ** 2 / 2)[:, None] * gz[None, :]
    cp = np.zeros_like(diag); dp = np.zeros_like(rhs)
    cp[:, 0] = off / diag[:, 0]; dp[:, 0] = rhs[:, 0] / diag[:, 0]
    for j in range(1, Nz):
        den = diag[:, j] - off * cp[:, j - 1]
        cp[:, j] = off / den; dp[:, j] = (rhs[:, j] - off * dp[:, j - 1]) / den
    chi = np.zeros_like(rhs); chi[:, -1] = dp[:, -1]
    for j in range(Nz - 2, -1, -1):
        chi[:, j] = dp[:, j] - cp[:, j] * chi[:, j + 1]
    k0idx = int(np.argmin(uk)); mean_chi = np.trapezoid(chi[k0idx], zg) / L           # <chi_q> (k = 0 mode, includes 1/A)
    hmean = 2 / (KAP * L) * math.tanh(KAP * L / 2)
    pbar = mean_chi / hmean
    out = {}
    kx = (2 * np.pi / L) * mx.ravel()
    for rho in rho_lines:
        field = np.zeros(Nz)
        wgt = np.cos(kx * rho)                                  # line y = 0, x = rho
        for ui in range(len(uk)):
            field += wgt[inv == ui].sum() * chi[ui]
        field = field + pbar * (1 - np.cosh(KAP * zg) / np.cosh(KAP * L / 2))
        out[rho] = np.interp(z_out, zg, field)
    return out, pbar

def ludwig(dx):
    p = sorted(glob.glob(os.path.join(ROOT, "C-electrolito-hann8-LD6", "z_09.50", "psi-0*")))[-1]
    a = np.loadtxt(p).reshape(32, 32, 34); zc = np.arange(1, 35) - 1.5
    return zc[1:33] - L / 2, a[(15 + dx) % 32, 15, 1:33]

if __name__ == "__main__":
    rhos = (0, 6, 10); zz = np.linspace(-L / 2 + 0.2, L / 2 - 0.2, 120)
    z32 = ludwig(0)[0]
    num, pbar_n = numeric(rhos, np.concatenate([zz, z32]))
    print("phibar: numeric %.6f" % pbar_n)
    print("TRUNC = %s; fraction of the full Gaussian inside the plates = %.3e"
          % (TRUNC, _gint(-L / 2, L / 2, D) / (SIG * math.sqrt(2 * math.pi))))
    plt.rcParams.update({"font.size": 11, "axes.grid": True, "grid.color": "#e4e2dc",
                         "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False})
    fig, a = plt.subplots(figsize=(7.2, 5.0)); cols = {0: "#1baf7a", 6: "#2a78d6", 10: "#eb6834"}
    print(" rho | max|analytic - numeric solver| | max|Ludwig - analytic| (outside kernel) | (inside)  | max(Ludwig)")
    for rho in rhos:
        an_z, pbar = analytic(float(rho), np.concatenate([zz, z32]))
        n_full = np.concatenate([zz, z32]); nz = len(zz)
        z, ps = ludwig(rho); outside = np.abs(z - D) > 4.5
        a.plot(zz, an_z[:nz], "-", color=cols[rho], lw=1.6)
        a.plot(z, ps, "o", ms=3.8, color=cols[rho], label="Ludwig, $\\rho=%d$" % rho)
        d_an_num = np.abs(an_z - num[rho]).max()
        d_l = np.abs(ps - an_z[nz:]); print("  %2d |  %.5f (phibar analytic %.6f)  |  %.4f  | %.4f | %.4f"
              % (rho, d_an_num, pbar, d_l[outside].max(), d_l[~outside].max(), ps.max()))
    a.plot([], [], "-", color="#55544f", label="derived solution (Gaussian charge, mean potential, lateral images)")
    a.set_xlabel("z (distance from the mid-plane, lattice units)"); a.set_ylabel(r"reduced potential $\psi = e\phi/k_BT$")
    a.set_title("1:1 electrolyte, $\\lambda_D=6$, grounded plates (d = 8 from the lower plate)")
    a.legend(loc="center right", fontsize=9, bbox_to_anchor=(1.0, 0.62))
    fig.tight_layout(); fig.savefig(os.path.join(ROOT, "graficos", "C_electrolyte_gaussian_derived_vs_ludwig.png"), dpi=170, bbox_inches="tight")
