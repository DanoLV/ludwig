#!/usr/bin/env python3
"""Radial profile of the ionic (fluid) charge around one subgrid charge, with
and without the pm_sr short-range correction, against theory.

The qsi dump holds the two ion densities (rho0, rho1); the net node charge is
rho0 - rho1, in the same units as the particle charge (verified: the sum over
the box is exactly -q). In linearised Poisson-Boltzmann the induced charge is
rho_ion = -eps * kappa^2 * phi, so the same screened potential used for the
field gives the charge with no extra calibration.
"""
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import erfc


def phi_gauss(d, sigma, kappa, q, eps):
    B = sigma * kappa / np.sqrt(2.0)
    C = 1.0 / (sigma * np.sqrt(2.0))
    t1 = np.exp(-kappa * d) * erfc(B - C * d)
    t2 = np.exp(kappa * d) * erfc(B + C * d)
    return q / (8.0 * np.pi * eps) * np.exp(kappa ** 2 * sigma ** 2 / 2.0) * (t1 - t2) / d


def phi_point(d, kappa, q, eps):
    return q * np.exp(-kappa * d) / (4.0 * np.pi * eps * d)


def geom(n, pos, rmax):
    g = np.arange(1, n + 1, dtype=float)
    X, Y, Z = np.meshgrid(g, g, g, indexing='ij')
    dx, dy, dz = X - pos[0], Y - pos[1], Z - pos[2]
    for v in (dx, dy, dz):
        v -= n * np.round(v / n)
    r = np.sqrt(dx * dx + dy * dy + dz * dz)
    m = (r > 0.0) & (r <= rmax)
    return r[m], (dx[m], dy[m], dz[m]), m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--no-corr', default=None, help='qsi dump of the uncorrected run')
    ap.add_argument('--corr', default=None, help='qsi dump of the corrected run')
    ap.add_argument('-n', '--ntotal', type=int, default=48)
    ap.add_argument('--pos', type=float, nargs=3, default=[24.3, 24.0, 24.0])
    ap.add_argument('-q', '--charge', type=float, default=1.0)
    ap.add_argument('--epsilon', type=float, default=1.0e4)
    ap.add_argument('--ld', type=float, required=True)
    ap.add_argument('--sigma', type=float, required=True)
    ap.add_argument('--rcut', type=float, default=None)
    ap.add_argument('--rmax', type=float, default=18.0)
    ap.add_argument('--images', type=int, default=2)
    ap.add_argument('--kernel', default=None)
    ap.add_argument('-o', '--output', default='charge_vs_r.png')
    a = ap.parse_args()

    n, q, eps, kappa = a.ntotal, a.charge, a.epsilon, 1.0 / a.ld
    rr, dvec, m = geom(n, a.pos, a.rmax)

    def load(p):
        d = np.loadtxt(p)
        return (d[:, 0] - d[:, 1]).reshape(n, n, n)[m]

    if not a.no_corr and not a.corr:
        raise SystemExit('give --no-corr, --corr, or both')
    # each run is optional, so one file alone produces its own figure
    runs = []
    if a.no_corr:
        runs.append((load(a.no_corr), 'tab:red', 'no correction'))
    if a.corr:
        runs.append((load(a.corr), 'tab:blue', 'with correction'))

    off = np.arange(-a.images, a.images + 1) * n
    th = {}
    for kind in ('point', 'gauss'):
        acc = np.zeros_like(rr)
        for ox in off:
            for oy in off:
                for oz in off:
                    dd = np.sqrt((dvec[0] + ox) ** 2 + (dvec[1] + oy) ** 2 + (dvec[2] + oz) ** 2)
                    phi = (phi_point(dd, kappa, q, eps) if kind == 'point'
                           else phi_gauss(dd, a.sigma, kappa, q, eps))
                    acc += -eps * kappa ** 2 * phi
        th[kind] = acc

    fig, (ax, axe) = plt.subplots(2, 1, figsize=(9.5, 8.6), sharex=True,
                                  gridspec_kw={'height_ratios': [3, 1.15]})
    o = np.argsort(rr)
    rs = rr[o]
    for qv, c, lab in runs:
        ax.plot(rr, -qv, '.', color=c, ms=2.5, alpha=0.30, label=f'Lattice, {lab}')
    ax.plot(rs, -th['point'][o], '-.', color='tab:green', lw=1.8,
            label=r'Point Debye-Hückel ($\lambda_D$=' + f'{a.ld:.4g})')
    ax.plot(rs, -th['gauss'][o], ':', color='k', lw=2.4,
            label=r'Gaussian cloud, Debye-Hückel ($\sigma$=' + f'{a.sigma:g})  — target')
    if a.rcut:
        for b in (ax, axe):
            b.axvline(a.rcut, color='tab:orange', ls='--', lw=1.2)
        ax.text(a.rcut * 1.02, ax.get_ylim()[1] * 0.2, r'$r_{cut}$', color='tab:orange', fontsize=10)
    ax.set_yscale('log')
    ax.set_ylabel(r'Counter-charge density  $-(\rho_0-\rho_1)$', fontsize=12)
    sub = [f'single charge q={q:g}', f'box {n}$^3$']
    if a.kernel:
        sub.append(f'{a.kernel} kernel')
    what = ('with and without the short-range correction' if len(runs) == 2
            else runs[0][2])
    ax.set_title(f'Ionic charge around one particle, {what}\n'
                 + ',   '.join(sub), fontsize=13)
    ax.grid(True, alpha=0.3, which='both')
    ax.legend(fontsize=9.5)

    for qv, c, lab in runs:
        axe.plot(rr, 100.0 * (qv - th['gauss']) / th['gauss'], '.', color=c, ms=2.5,
                 alpha=0.30, label=lab)
    axe.axhline(0.0, color='k', ls='--', lw=1.0)
    for y in (1.0, -1.0):
        axe.axhline(y, color='tab:gray', ls=':', lw=0.9)
    axe.set_ylim(-60, 60)
    axe.set_xlabel('Distance $r$ from the particle (lattice units)', fontsize=12)
    axe.set_ylabel('relative error vs\nGaussian cloud (%)', fontsize=10)
    axe.grid(True, alpha=0.3)
    axe.legend(fontsize=9, markerscale=3)

    fig.tight_layout()
    fig.savefig(a.output, dpi=150)
    print(f'saved: {a.output}\n')
    head = ''.join(f"{lab:>18}" for _, _, lab in runs)
    print(f"{'r range':>12}{head}")
    for lo, hi in ((0, 2), (2, 4), (4, 6), (6, 8), (8, 12), (12, 18)):
        s = (rr >= lo) & (rr < hi)
        if s.any():
            cells = ''.join(
                f"{np.median(100.0*(qv[s]-th['gauss'][s])/th['gauss'][s]):+17.2f}%"
                for qv, _, _ in runs)
            print(f"  [{lo:2d},{hi:2d}){cells}")


if __name__ == '__main__':
    main()
