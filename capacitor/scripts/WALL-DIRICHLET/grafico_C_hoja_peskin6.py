#!/usr/bin/env python3
"""Electrolyte capacitor (closed system): the solution of grafico_C_hoja_hann.py with the charge spread by Ludwig's
STRETCHED Peskin 6-point kernel (subgrid_kernel peskin6, subgrid_peskin6_stretch s; subgrid.c d_peskin6).

Charge density (continuous version of Ludwig's spread, in physical units = lattice units of the dx = 1 run):
  rho_part = q W(x) W(y) W(z - z0),   W(u) = phi6(u/S)/S,   S = 2 (s = 2r lattice units = 2 physical for every r)
  phi6 = Bao, Kaplan & Peskin (2016) 6-point kernel, support |v| < 3, C^3 (J. Comput. Phys. 347 (2017) 183;
  arXiv:1505.07529, Eqs. 2.16-2.23 as coded in subgrid.c peskin6_phi_segment). Support |u| < 3S = 6, sigma = S sqrt(K)
  = 1.69006 with K = 59/60 - sqrt(29)/20 (its second-moment condition).
Same derivation as the Hann sheet; only the source transform changes:
  phi   = phi_C + U + [kL/(2 tanh(kL/2))] <U> [1 - cosh(k z)/cosh(k L/2)]
  U     = (2 pi/A) sum_nm psit(kx, ky, z) exp(i k.r)
  psit  = q What(kx) What(ky)/(2 pi eps) int W(s - z0) g_a(s, z) ds,   a^2 = kx^2 + ky^2 + kappa^2
  g_a   = sinh(a(z_< + L/2)) sinh(a(L/2 - z_>)) / (a sinh(aL))      (Jackson 3rd ed., Problem 3.17, k^2 -> k^2 + kappa^2)
  What(k) = int phi6(v) cos(k S v) dv
phi6 has no elementary antiderivative (it contains a square root), so What and the z integral are done by
Gauss-Legendre on each smooth piece: the unit segments of phi6 (breakpoints z0 + S m) and, for the z integral, split
at s = z (kink of g). On each piece the integrand is analytic, so 24 nodes give ~machine precision (checked below
against scipy quad). The support [z0 - 6, z0 + 6] = [-14, -2] lies inside the channel (C1 = 1).
Checks printed: What(0) = 1, psit against quad, <U> against its closed form, mode-truncation NK vs 3NK/4.

Reduced units, q/eps = 10, L = l = 32, kappa = 1/6, z0 = -8.
Usage: python3 grafico_C_hoja_peskin6.py   (checks only; used as a module by grafico_C_refinamiento_peskin6.py)
"""
import math
import numpy as np
from scipy.integrate import quad
import grafico_C_formulas as F
import grafico_C_hoja_gauss as H

L, A, FOUR = F.L, F.A, F.FOUR
KAP = 1 / 6.0; Z0 = F.D; S = 2.0
K6 = 59.0 / 60 - math.sqrt(29.0) / 20
SIG = S * math.sqrt(K6)
NK = 64                                                 # lateral modes |n|,|m| <= NK (checked against 3NK/4)
NG = 24                                                 # Gauss-Legendre nodes per smooth piece
ZA, ZB = Z0 - 3 * S, Z0 + 3 * S
BRK = Z0 + S * np.arange(-3, 4)                         # breakpoints of W(s - z0)
XG, WG = np.polynomial.legendre.leggauss(NG)

def phi6(v):
    """Bao et al. 6-point kernel (unit support |v| < 3), vectorised port of subgrid.c d_peskin6_unit."""
    v = np.asarray(v, float); r = v - np.floor(v); seg = 3 + np.floor(v).astype(int)
    beta = 9 / 4 - 1.5 * (K6 + r * r) + (22 / 3 - 7 * K6) * r - (7 / 3) * r ** 3
    t1 = (3 * K6 - 1) * r + r ** 3; t2 = (4 - 3 * K6) * r - r ** 3
    gam = -11 / 32 * r * r + 3 / 32 * (2 * K6 + r * r) * r * r + t1 * t1 / 72 + t2 * t2 / 18
    p = (-beta + np.sqrt(np.maximum(beta * beta - 112 * gam, 0))) / 56
    out = np.select([seg == 0, seg == 1, seg == 2, seg == 3, seg == 4, seg == 5],
                    [p, -3 * p - 1 / 16 + (K6 + r * r) / 8 + (3 * K6 - 1) * r / 12 + r ** 3 / 12,
                     2 * p + 1 / 4 + (4 - 3 * K6) * r / 6 - r ** 3 / 6, 2 * p + 5 / 8 - (K6 + r * r) / 4,
                     -3 * p + 1 / 4 - (4 - 3 * K6) * r / 6 + r ** 3 / 6,
                     p - 1 / 16 + (K6 + r * r) / 8 - (3 * K6 - 1) * r / 12 - r ** 3 / 12], 0.0)
    return np.where(np.abs(v) < 3, out, 0.0)

def W(u):
    return phi6(np.asarray(u, float) / S) / S

def _nodes(a, b):
    return 0.5 * (b - a) * XG + 0.5 * (b + a), 0.5 * (b - a) * WG

_VN = np.concatenate([_nodes(m, m + 1)[0] for m in range(-3, 3)])
_VW = np.concatenate([_nodes(m, m + 1)[1] for m in range(-3, 3)])
_PV = phi6(_VN)

def What(k):
    k = np.asarray(k, float)
    return (np.cos(np.multiply.outer(k * S, _VN)) * (_PV * _VW)).sum(-1)

def _g(a, s, z):
    """g_a(s, z), stable for large a: a (modes, 1, 1), s and z broadcast."""
    zl, zg = np.minimum(s, z), np.maximum(s, z)
    x1, x2 = a * (zl + L / 2), a * (L / 2 - zg)
    return (np.exp(x1 + x2 - a * L) * (-np.expm1(-2 * x1)) * (-np.expm1(-2 * x2))
            / (2 * a * (-np.expm1(-2 * a * L))))

def Gt(a, z):
    """int W(s - z0) g_a(s, z) ds for a (modes,) and z (nz,) -> (modes, nz); pieces split at the breakpoints and z."""
    out = np.zeros((len(a), len(z)))
    for j, zj in enumerate(z):
        b = np.unique(np.concatenate([BRK, [zj]] if ZA < zj < ZB else [BRK]))
        s, w = map(np.concatenate, zip(*[_nodes(b[i], b[i + 1]) for i in range(len(b) - 1)]))
        out[:, j] = (_g(a[:, None], s[None, :], zj) * (W(s - Z0) * w)[None, :]).sum(1)
    return out

def psit(kx, ky, z):
    a = np.sqrt(kx ** 2 + ky ** 2 + KAP ** 2)
    return FOUR / (2 * math.pi) * (What(kx) * What(ky))[:, None] * Gt(a, z)

def U(rho, z, nk=None):
    """(2 pi/A) sum_nm psit e^{i k.r} on the line y = 0, x = rho: Gt depends on |k| only, computed once per |k|."""
    nk = NK if nk is None else nk
    z = np.atleast_1d(np.asarray(z, float))
    n = np.arange(-nk, nk + 1); kv = 2 * np.pi * n / L
    wh = What(kv); cx = np.cos(kv * rho) * wh
    nn, mm = np.meshgrid(n, n, indexing="ij"); q2 = (nn ** 2 + mm ** 2).ravel()
    uq, inv = np.unique(q2, return_inverse=True)
    coef = np.bincount(inv, weights=(cx[:, None] * wh[None, :]).ravel())    # sum over modes with the same |k|
    a = np.sqrt((2 * np.pi / L) ** 2 * uq + KAP ** 2); out = np.zeros(len(z))
    for i0 in range(0, len(uq), 400):
        sl = slice(i0, i0 + 400); out += (coef[sl, None] * Gt(a[sl], z)).sum(0)
    return 2 * np.pi / A * FOUR / (2 * math.pi) * out

def mean_U():
    """<U> = (q/eps)/(kappa^2 A L) [1 - int W(s - z0) cosh(kappa s) ds / cosh(kappa L/2)] (closed form of the Hann sheet)."""
    s = np.concatenate([_nodes(BRK[i], BRK[i + 1])[0] for i in range(6)])
    w = np.concatenate([_nodes(BRK[i], BRK[i + 1])[1] for i in range(6)])
    ic = np.sum(w * W(s - Z0) * np.cosh(KAP * s))
    return FOUR / (KAP ** 2 * A * L) * (1 - ic / math.cosh(KAP * L / 2))

def phi(rho, z, p1, p2, mu):
    fac = (KAP * L / 2) / math.tanh(KAP * L / 2)
    return H.cap_closed(z, p1, p2) + U(rho, z) + fac * mu * H.h(z)

def checks(mu):
    print("phi6: sum over a unit grid %.15f, 2nd moment %.6f (K = %.6f), What(0) = %.15f, sigma = %.5f"
          % (phi6(np.arange(-3, 4) + 0.3).sum(), quad(lambda v: v * v * float(phi6(v)), -3, 3, points=range(-2, 3))[0],
             K6, What(np.zeros(1))[0], SIG))
    for kx, ky, z in ((0.0, 0.0, -8.0), (0.6, 0.2, -6.3), (2.4, 1.1, -9.0), (0.4, 0.0, 5.0)):
        a = math.sqrt(kx * kx + ky * ky + KAP ** 2)
        qd = quad(lambda s: float(W(s - Z0)) * float(_g(a, s, z)), ZA, ZB, points=list(BRK[1:-1]) + [z], limit=400)[0]
        wk = quad(lambda v: float(phi6(v)) * math.cos(kx * S * v), -3, 3, points=range(-2, 3))[0] * \
             quad(lambda v: float(phi6(v)) * math.cos(ky * S * v), -3, 3, points=range(-2, 3))[0]
        print("psit(%.1f,%.1f,z=%5.1f): Gauss pieces %.12e | quad %.12e"
              % (kx, ky, z, psit(np.array([kx]), np.array([ky]), np.array([z]))[0, 0], FOUR / (2 * math.pi) * wk * qd))
    x, w = np.polynomial.legendre.leggauss(200); tot = 0.0
    for lo, hi in ((-L / 2, ZA),) + tuple((BRK[i], BRK[i + 1]) for i in range(6)) + ((ZB, L / 2),):
        z = 0.5 * (hi - lo) * x + 0.5 * (hi + lo); tot += 0.5 * (hi - lo) * np.sum(w * psit(np.zeros(1), np.zeros(1), z)[0])
    print("<U>: quadrature of psit(0,z) %.10e | closed form %.10e" % (2 * np.pi / (A * L) * tot, mu))
    zz = np.linspace(-L / 2 + 0.25, L / 2 - 0.25, 64)
    print("mode truncation: max|U(NK=%d) - U(NK=%d)| on the axis = %.1e (max U = %.3f)"
          % (NK, 3 * NK // 4, np.abs(U(0.0, zz) - U(0.0, zz, 3 * NK // 4)).max(), U(0.0, zz).max()))

if __name__ == "__main__":
    checks(mean_U())
