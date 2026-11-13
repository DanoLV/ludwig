#!/usr/bin/env python3
"""Analytic linear Debye-Hueckel reference of grafico_C_hoja_hann.py generalised to the Hann window cos^p (2026-11-09,
user's request: reference for the cos^4 kernel, subgrid.c d_hann with subgrid_hann_power 4).

Kernel (subgrid.c, CHANGE 20261104): w(u) = cos^p(pi u/n) / (n C_p), |u| < n/2, C_p = (p-1)!!/p!!.
For even p, cos^p(theta) = 2^-p [binom(p, p/2) + 2 sum_{m=1}^{p/2} binom(p, p/2 - m) cos(2 m theta)], so with
omega = 2 pi/n
  W(u) = (1/n) sum_{m=0}^{p/2} c_m cos(m omega u),   c_0 = 1,  c_m = 2 binom(p, p/2 - m)/binom(p, p/2)
  p = 2: c = (1, 1)   [(1/n)(1 + cos(omega u)), the original Hann]
  p = 4: c = (1, 4/3, 1/3)
Only the source changes in the derivation of grafico_C_hoja_hann.py; the three kernel-dependent pieces become sums over
the harmonics m omega:
  What(k) = sum_m (c_m / (1 + [m > 0])) [sinc(x - m pi) + sinc(x + m pi)] with x = k n/2, sinc(y) = sin(y)/y
            (m = 0: sinc(x); integral of cos(m omega u) cos(k u) over |u| < n/2)
  int e^{b s} W(s - z0) ds = e^{b s} (1/n) sum_m c_m (b cos(m omega u) + m omega sin(m omega u)) / (b^2 + (m omega)^2),
            u = s - z0   (m = 0 term: 1/b)
  W'(u)   = -(1/n) sum_m c_m m omega sin(m omega u)
set_kernel(n, p) replaces these functions (and NH, OM) in the module grafico_C_hoja_hann (HH) and W' in
grafico_C_fuerza_desplazada, so psit, dpsit, U, dU, mean_U, phi, dphi, checks and F_kernel use the new kernel without
modifying those files. After set_kernel, call grafico_C_desplazada.set_charge(delta) (support limits ZA, ZB use NH).
With p = 2, n = 8 the functions coincide with the originals (checked by kernel_checks()).
"""
import math
import numpy as np
from scipy.integrate import quad
import grafico_C_hoja_hann as HH
import grafico_C_fuerza_desplazada as FD

_ORIG = dict(W=HH.W, What=HH.What, _J=HH._J, NH=HH.NH, OM=HH.OM, Wp=FD.Wp)
COEF = None

def coefs(p):
    b0 = math.comb(p, p // 2)
    return np.array([1.0] + [2.0 * math.comb(p, p // 2 - m) / b0 for m in range(1, p // 2 + 1)])

def set_kernel(n, p):
    global COEF
    assert p >= 2 and p % 2 == 0
    COEF = coefs(p); HH.NH = float(n); HH.OM = 2 * math.pi / n
    HH.W, HH.What, HH._J, FD.Wp = W, What, _J, Wp

def W(u):
    n, om = HH.NH, HH.OM
    v = sum(c * np.cos(m * om * u) for m, c in enumerate(COEF)) / n
    return np.where(np.abs(u) < n / 2, v, 0.0)

def Wp(u):
    n, om = HH.NH, HH.OM
    v = -sum(c * m * om * np.sin(m * om * u) for m, c in enumerate(COEF)) / n
    return np.where(np.abs(u) < n / 2, v, 0.0)

def What(k):
    x = np.abs(k) * HH.NH / 2
    s = lambda y: np.sinc(y / np.pi)                    # sin(y)/y
    v = s(x)
    for m, c in enumerate(COEF[1:], start=1):
        v = v + 0.5 * c * (s(x - m * np.pi) + s(x + m * np.pi))
    return v

def _J(b, c, s1, s2):
    """int_{s1}^{s2} e^{b(s-c)} W(s-z0) ds (s1 <= s2 inside the support), closed form; arrays broadcast."""
    n, om = HH.NH, HH.OM
    def P(s):
        u = s - HH.Z0; tot = 1 / b
        for m, cm in enumerate(COEF[1:], start=1):
            w = m * om
            tot = tot + cm * (b * np.cos(w * u) + w * np.sin(w * u)) / (b * b + w * w)
        return np.exp(b * (s - c)) / n * tot
    return np.where(s2 > s1, P(s2) - P(np.minimum(s1, s2)), 0.0)

def kernel_checks(p):
    """Prints: W against cos^p/(n C_p); int W = 1; What against quadrature of the transform; _J against quadrature."""
    n = HH.NH; cp = np.prod([(k - 1) / k for k in range(2, p + 1, 2)])
    u = np.linspace(-n / 2 + 1e-9, n / 2 - 1e-9, 2001)
    d0 = np.abs(W(u) - np.cos(np.pi * u / n) ** p / (n * cp)).max()
    i0 = quad(lambda t: float(W(t)), -n / 2, n / 2)[0]
    dk = max(abs(float(What(np.array([k]))[0]) - quad(lambda t: float(W(t)) * math.cos(k * t), -n / 2, n / 2, limit=200)[0])
             for k in (0.0, 0.3, 2 * math.pi / n, 1.7, 4.0, 9.1))
    dw = max(abs(float(Wp(np.array([t]))[0]) - (float(W(np.array([t + 1e-6]))[0]) - float(W(np.array([t - 1e-6]))[0])) / 2e-6)
             for t in (-0.4 * n, -0.1 * n, 0.05 * n, 0.33 * n))
    z0 = HH.Z0; dj = 0.0
    for b, c, s1, s2 in ((0.7, 0.0, z0 - n / 3, z0 + n / 4), (-1.3, 2.0, z0 - n / 2, z0), (0.2, -5.0, z0, z0 + n / 2)):
        qd = quad(lambda s: math.exp(b * (s - c)) * float(W(np.array([s - z0]))[0]), s1, s2, limit=200)[0]
        dj = max(dj, abs(float(_J(b, c, np.array(s1), np.array(s2))) - qd) / abs(qd))
    print("kernel cos^%d n = %g: max|W - cos^p/(n C_p)| %.1e | int W - 1 = %.1e | max|What - quadrature| %.1e | "
          "max|W' - finite diff.| %.1e | max rel |J - quadrature| %.1e" % (p, n, d0, i0 - 1, dk, dw, dj))
