#!/usr/bin/env python3
"""Electrolyte capacitor (closed system): the solution of grafico_C_hoja_hann.py with Ludwig's generalised Hann window
w ~ cos^p (subgrid_hann_power p; subgrid.c d_hann), support n = 8 physical, p even (here p = 4: C^3).

Charge density (continuous version of Ludwig's spread):
  rho_part = q W(x) W(y) W(z - z0),   W(u) = cos^p(pi u/n) / (n C_p) for |u| < n/2,   C_p = (p-1)!!/p!!
For even p, cos^p(theta) = 2^-p [binom(p, p/2) + 2 sum_{m=1}^{p/2} binom(p, p/2 - m) cos(2 m theta)], so
  W(u) = sum_{m=0}^{p/2} c_m cos(m w u),   w = 2 pi/n
(p = 2: c = (1, 1)/n, the Hann sheet; p = 4: c = (3, 4, 1)/(3n)). Every term has elementary integrals, so the sheet's
closed forms carry over term by term:
  What(k) = sum_m c_m (n/2) [sinc((k - m w) n/2) + sinc((k + m w) n/2)],   sinc(x) = sin(x)/x
  int e^{b s} cos(m w (s - z0)) ds = e^{b s} (b cos(m w (s - z0)) + m w sin(m w (s - z0)))/(b^2 + m^2 w^2)
and psit, U, <U>, phi exactly as in grafico_C_hoja_hann.py (Green function in z: Jackson 3rd ed., Problem 3.17).
Checks printed: p = 2 reproduces grafico_C_hoja_hann.py, psit against scipy quad (also at large k), <U> against its
closed form, mode truncation NK vs 3NK/4.

Reduced units, q/eps = 10, L = l = 32, kappa = 1/6, z0 = -8.
Usage: python3 grafico_C_hoja_hannp.py   (checks; used as a module, P = 4 by default: set_power(p) to change)
"""
import math
import numpy as np
from scipy.integrate import quad
from scipy.special import comb
import grafico_C_formulas as F
import grafico_C_hoja_gauss as H

L, A, FOUR = F.L, F.A, F.FOUR
KAP = 1 / 6.0; Z0 = F.D; NH = 8.0; OM = 2 * math.pi / NH
NK = 64                                                 # lateral modes |n|,|m| <= NK (What ~ k^-(p+1))
ZA, ZB = Z0 - NH / 2, Z0 + NH / 2
P = 4; CM = None; SIG = None

def set_power(p):
    """Harmonic coefficients c_m of W = cos^p(pi u/n)/(n C_p) and the width sigma."""
    global P, CM, SIG
    assert p >= 2 and p % 2 == 0
    cp = np.prod([(k - 1) / k for k in range(2, p + 1, 2)])
    cm = np.array([comb(p, p // 2)] + [2 * comb(p, p // 2 - m) for m in range(1, p // 2 + 1)]) / 2.0 ** p
    P, CM = p, cm / (NH * cp)
    SIG = math.sqrt(quad(lambda u: u * u * float(W(u)), -NH / 2, NH / 2)[0])
def W(u):
    u = np.asarray(u, float)
    return np.where(np.abs(u) < NH / 2, sum(c * np.cos(m * OM * u) for m, c in enumerate(CM)), 0.0)

set_power(P)

def What(k):
    k = np.abs(np.asarray(k, float)); a = NH / 2
    return sum(c * a * (np.sinc((k - m * OM) * a / np.pi) + np.sinc((k + m * OM) * a / np.pi)) for m, c in enumerate(CM))

def _J(b, c, s1, s2):
    """int_{s1}^{s2} e^{b(s-c)} W(s-z0) ds (s1 <= s2 inside the support), closed form term by term."""
    def Pr(s):
        u = s - Z0; e = np.exp(b * (s - c)); out = CM[0] * e / b
        for m in range(1, len(CM)):
            wm = m * OM
            out = out + CM[m] * e * (b * np.cos(wm * u) + wm * np.sin(wm * u)) / (b * b + wm * wm)
        return out
    return np.where(s2 > s1, Pr(s2) - Pr(np.minimum(s1, s2)), 0.0)

def Gt(a, z):
    """int W(s - z0) g_a(s, z) ds, (modes, z): the four exponentials of grafico_C_hoja_hann.psit (non-positive exponents)."""
    a = a[:, None]; zz = z[None, :]
    lo, hi = ZA, np.clip(zz, ZA, ZB)
    S = (_J(a, zz, lo, hi) - _J(-a, -(L + zz), lo, hi) - _J(a, L - zz, lo, hi) + _J(-a, -(2 * L - zz), lo, hi))
    lo, hi = np.clip(zz, ZA, ZB), ZB
    S += (_J(-a, zz, lo, hi) - _J(-a, -(L + zz), lo, hi) - _J(a, L - zz, lo, hi) + _J(a, 2 * L + zz, lo, hi))
    return S / (2 * a * (1 - np.exp(-2 * a * L)))

def psit(kx, ky, z):
    a = np.sqrt(kx ** 2 + ky ** 2 + KAP ** 2)
    return FOUR / (2 * math.pi) * (What(kx) * What(ky))[:, None] * Gt(a, z)

def U(rho, z, nk=None):
    """(2 pi/A) sum_nm psit e^{i k.r} on the line y = 0, x = rho; Gt depends on |k| only."""
    nk = NK if nk is None else nk
    z = np.atleast_1d(np.asarray(z, float))
    n = np.arange(-nk, nk + 1); kv = 2 * np.pi * n / L
    wh = What(kv); cx = np.cos(kv * rho) * wh
    nn, mm = np.meshgrid(n, n, indexing="ij"); q2 = (nn ** 2 + mm ** 2).ravel()
    uq, inv = np.unique(q2, return_inverse=True)
    coef = np.bincount(inv, weights=(cx[:, None] * wh[None, :]).ravel())
    a = np.sqrt((2 * np.pi / L) ** 2 * uq + KAP ** 2); out = np.zeros(len(z))
    for i0 in range(0, len(uq), 2000):
        sl = slice(i0, i0 + 2000); out += (coef[sl, None] * Gt(a[sl], z)).sum(0)
    return 2 * np.pi / A * FOUR / (2 * math.pi) * out

def mean_U():
    """<U> = (q/eps)/(kappa^2 A L) [1 - int W(s - z0) cosh(kappa s) ds / cosh(kappa L/2)]."""
    ic = quad(lambda s: float(W(s - Z0)) * math.cosh(KAP * s), ZA, ZB, epsabs=1e-14, epsrel=1e-13)[0]
    return FOUR / (KAP ** 2 * A * L) * (1 - ic / math.cosh(KAP * L / 2))

def phi(rho, z, p1, p2, mu):
    fac = (KAP * L / 2) / math.tanh(KAP * L / 2)
    return H.cap_closed(z, p1, p2) + U(rho, z) + fac * mu * H.h(z)

def checks():
    import grafico_C_hoja_hann as HH
    p0 = P; set_power(2); zz = np.linspace(-L / 2 + 0.25, L / 2 - 0.25, 64)
    d = np.abs(U(0.0, zz, HH.NK) - HH.U(0.0, zz)).max()
    print("p = 2 vs grafico_C_hoja_hann.py (same NK): max|U| diff %.1e, <U> %.12e vs %.12e, sigma %.5f"
          % (d, mean_U(), HH.mean_U(), SIG))
    set_power(p0)
    g = lambda a, s, z: (np.sinh(a * (min(z, s) + L / 2)) * np.sinh(a * (L / 2 - max(z, s))) / (a * np.sinh(a * L)))
    print("p = %d: c_m = %s, sigma = %.5f (8 sqrt(1/12 - 5/(8 pi^2)) = %.5f), What(0) = %.15f"
          % (P, np.array2string(CM * NH, precision=6), SIG, 8 * math.sqrt(1 / 12 - 5 / (8 * math.pi ** 2)), What(0.0)))
    for kx, ky, z in ((0.0, 0.0, -8.0), (0.6, 0.2, -6.3), (2.4, 1.1, -9.0), (0.4, 0.0, 5.0), (9.0, 6.0, -7.3)):
        a = math.sqrt(kx * kx + ky * ky + KAP ** 2)
        qd = quad(lambda s: float(W(s - Z0)) * g(a, s, z), ZA, ZB, points=[z, Z0], limit=400, epsabs=0, epsrel=1e-12)[0]
        wk = [quad(lambda u: float(W(u)) * math.cos(k * u), -NH / 2, NH / 2, limit=400, epsabs=0, epsrel=1e-12)[0] for k in (kx, ky)]
        print("  psit(%.1f,%.1f,z=%5.1f): closed form %.12e | quad %.12e"
              % (kx, ky, z, psit(np.array([kx]), np.array([ky]), np.array([z]))[0, 0], FOUR / (2 * math.pi) * wk[0] * wk[1] * qd))
    x, w = np.polynomial.legendre.leggauss(400); tot = 0.0
    for lo, hi in ((-L / 2, ZA), (ZA, Z0), (Z0, ZB), (ZB, L / 2)):
        z = 0.5 * (hi - lo) * x + 0.5 * (hi + lo); tot += 0.5 * (hi - lo) * np.sum(w * psit(np.zeros(1), np.zeros(1), z)[0])
    print("  <U>: quadrature of psit(0,z) %.10e | closed form %.10e" % (2 * np.pi / (A * L) * tot, mean_U()))
    print("  mode truncation: max|U(NK=%d) - U(NK=%d)| on the axis = %.1e (max U = %.3f)"
          % (NK, 3 * NK // 4, np.abs(U(0.0, zz) - U(0.0, zz, 3 * NK // 4)).max(), U(0.0, zz).max()))

if __name__ == "__main__":
    checks()
