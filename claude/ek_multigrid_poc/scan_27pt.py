#!/usr/bin/env python3
"""
scan_27pt.py -- pair force on the refined grid WITHOUT short-range correction,
using Ludwig's D3Q27 stencils (stencil_d3q27.c):

    laplacian(psi)(x) = 6 sum_p w_p [psi(x+c_p) - psi(x)]
    grad(psi)(x)      = 3 sum_p w_p c_p psi(x+c_p)

both applied in Fourier space (exactly the same operators on a periodic box).

Two spread/gather kernels:
  hann6  : Hann order 6 defined in FINE nodes (physical sigma = 1.085/r)
  gauss  : Ludwig's truncated Gaussian d_gauss (support 6 sigma, renormalised
           over the nodes), with a FIXED PHYSICAL sigma_c (coarse units),
           i.e. sigma_f = r sigma_c and support 6 r sigma_c fine nodes.

References (periodic box, +1/-1 pair):
  Coulomb  : point charges (Ewald)
  Gauss-G  : two Gaussian clouds of width sigma_c (reciprocal sum with
             form factor exp(-k^2 sigma^2))

Usage:  python3 scan_27pt.py [--sigma 1.0] [--nplace 12]
"""
import argparse
import math
import numpy as np
import ek_multigrid_poc as poc

# D3Q27 velocities and weights
CV = np.array([(i, j, k) for i in (-1, 0, 1) for j in (-1, 0, 1)
               for k in (-1, 0, 1)])
WV = np.array([{0: 8 / 27, 1: 2 / 27, 2: 1 / 54, 3: 1 / 216}[int(np.sum(c * c))]
               for c in CV])


def symbols(n):
    """Fourier symbols of the D3Q27 laplacian (as -lap >= 0) and gradient."""
    k = 2 * np.pi * np.fft.fftfreq(n)
    KX, KY, KZ = np.meshgrid(k, k, k, indexing="ij")
    neg_lap = np.zeros_like(KX)
    grad = [np.zeros_like(KX) for _ in range(3)]
    for c, w in zip(CV, WV):
        ph = KX * c[0] + KY * c[1] + KZ * c[2]
        neg_lap += -6 * w * (np.cos(ph) - 1)
        for d in range(3):
            grad[d] += 3 * w * c[d] * np.sin(ph)    # grad_k = i * grad[d]
    return neg_lap, grad


_SYM = {}


def hann_1d(n_order):
    def k1d(xf, n):
        half = n_order / 2
        m = np.arange(math.floor(xf - half) + 1, math.floor(xf + half) + 1)
        t = xf - m
        w = np.where(np.abs(t) < half, (1 + np.cos(2 * np.pi * t / n_order)) / n_order, 0.0)
        return (m - 1) % n, w
    return k1d


def gauss_1d(sigma_f, support_f):
    """Exactly Ludwig's d_gauss: truncated at support/2, renormalised."""
    rc = 0.5 * support_f
    a = 0.5 / sigma_f**2
    gcut = math.exp(-a * rc * rc)

    def k1d(xf, n):
        m = np.arange(math.floor(xf - rc) + 1, math.floor(xf + rc) + 1)
        t = xf - m
        g = np.where(np.abs(t) < rc, np.exp(-a * t * t) - gcut, 0.0)
        return (m - 1) % n, g / g.sum()
    return k1d


def pm_forces_27(pos_c, q, L, r, eps_c, k1d):
    n = r * L
    eps_f = eps_c / r
    if n not in _SYM:
        _SYM[n] = symbols(n)
    neg_lap, grad = _SYM[n]
    rho = np.zeros((n, n, n))
    st_all = []
    for p, qi in zip(pos_c, q):
        xf = poc.c2f_pos(p, r)
        st = [k1d(xf[d], n) for d in range(3)]
        st_all.append(st)
        (ix, wx), (iy, wy), (iz, wz) = st
        rho[np.ix_(ix, iy, iz)] += qi * np.einsum("i,j,k->ijk", wx, wy, wz)
    rk = np.fft.fftn(rho)
    nl = neg_lap.copy()
    nl[0, 0, 0] = 1.0
    pk = rk / (eps_f * nl)
    pk[0, 0, 0] = 0.0
    E = [np.real(np.fft.ifftn(-1j * grad[d] * pk)) for d in range(3)]
    F = np.zeros((len(q), 3))
    for a, (st, qi) in enumerate(zip(st_all, q)):
        (ix, wx), (iy, wy), (iz, wz) = st
        W = np.einsum("i,j,k->ijk", wx, wy, wz)
        for d in range(3):
            F[a, d] = r * qi * np.sum(W * E[d][np.ix_(ix, iy, iz)])
    return F


def gauss_pair_reference(pos, q, L, eps, sigma):
    """Periodic Gaussian clouds (width sigma each): reciprocal sum is exact."""
    C = 1 / (4 * np.pi * eps)
    kmax = 5.3 / sigma
    nmax = int(math.ceil(kmax * L / (2 * np.pi)))
    nv = np.arange(-nmax, nmax + 1)
    KX, KY, KZ = np.meshgrid(nv, nv, nv, indexing="ij")
    kv = 2 * np.pi / L * np.stack([KX.ravel(), KY.ravel(), KZ.ravel()], 1)
    k2 = np.sum(kv**2, 1)
    keep = (k2 > 0) & (k2 <= kmax**2)
    kv, k2 = kv[keep], k2[keep]
    pref = 4 * np.pi / L**3 * np.exp(-k2 * sigma**2) / k2
    F = np.zeros((len(q), 3))
    for i in range(len(q)):
        s = np.zeros(len(kv))
        for j in range(len(q)):
            s += q[j] * np.sin(kv @ (pos[i] - pos[j]))
        F[i] = C * q[i] * (pref * s) @ kv
    return F


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sigma", type=float, default=1.0,
                    help="physical Gaussian width (coarse units)")
    ap.add_argument("--nplace", type=int, default=12)
    args = ap.parse_args()

    L, eps, q = 16, 1 / (4 * np.pi), np.array([1.0, -1.0])
    ds = [1.0, 1.5, 2.0, 3.0, 4.0]
    rng = np.random.default_rng(7)
    plac = {}
    for d in ds:
        plac[d] = []
        for _ in range(args.nplace):
            p1 = 0.5 + L * rng.random(3)
            nh = rng.normal(size=3)
            nh /= np.linalg.norm(nh)
            pos = np.array([p1, p1 + d * nh])
            plac[d].append((pos, poc.ewald_forces(pos, q, L, eps),
                            gauss_pair_reference(pos, q, L, eps, args.sigma)))

    print("D3Q27 stencils, no short-range correction. Relative pair-force "
          f"error, mean (max) over {args.nplace} placements; d in coarse units")
    print(f"gauss kernel: physical sigma = {args.sigma} (sigma_f = r*sigma, "
          "support 6 sigma_f)\n")
    print("Gauss-G reference vs Coulomb at each d: " + "  ".join(
        f"d={d}: {np.mean([np.linalg.norm(G[0])/np.linalg.norm(C[0]) for _, C, G in plac[d]]):.3f}"
        for d in ds))
    print()
    hdr = f"{'kernel':6s} {'ref':8s} {'r':>2s} " + "".join(f"  d={d:<5}          " for d in ds)
    print(hdr)
    for name in ("hann6", "gauss"):
        for r in (1, 2, 4):
            k1d = hann_1d(6) if name == "hann6" else gauss_1d(r * args.sigma, 6 * r * args.sigma)
            res = {d: [] for d in ds}
            for d in ds:
                for pos, Fc, Fg in plac[d]:
                    F = pm_forces_27(pos, q, L, r, eps, k1d)
                    res[d].append((np.linalg.norm(F[0] - Fc[0]) / np.linalg.norm(Fc[0]),
                                   np.linalg.norm(F[0] - Fg[0]) / np.linalg.norm(Fg[0]),
                                   np.linalg.norm(F[0] + F[1]) / np.linalg.norm(Fc[0])))
            for iref, ref in ((0, "Coulomb"), (1, "Gauss-G")):
                row = " ".join(
                    f"{100*np.mean([e[iref] for e in res[d]]):7.2f}% ({100*np.max([e[iref] for e in res[d]]):6.2f}%)"
                    for d in ds)
                print(f"{name:6s} {ref:8s} {r:2d} {row}")
            mom = max(e[2] for d in ds for e in res[d])
            print(f"{'':6s} {'|F1+F2|':8s} {r:2d}  max {mom:.1e}")
        print()


if __name__ == "__main__":
    main()
