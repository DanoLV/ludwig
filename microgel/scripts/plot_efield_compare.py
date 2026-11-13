#!/usr/bin/env python3
"""Radial field around one subgrid charge, with and without the pm_sr
short-range correction, on the same axes.

Note the correction enters the stored efield only through
pm_sr_apply_force_correction (the particle<->node term). The particle-particle
pair correction never touches the field, and with a single particle it would
do nothing anyway.
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


def field_of(kind, d, sigma, kappa, q, eps, h=1e-4):
    if kind == 'coulomb':
        return q / (4.0 * np.pi * eps * d ** 2)
    if kind == 'dh':
        return q / (4.0 * np.pi * eps) * np.exp(-kappa * d) * (kappa / d + 1.0 / d ** 2)
    return -(phi_gauss(d + h, sigma, kappa, q, eps)
             - phi_gauss(d - h, sigma, kappa, q, eps)) / (2.0 * h)


def radial(path, n, pos, rmax):
    d = np.loadtxt(path).reshape(n, n, n, 3)
    g = np.arange(1, n + 1, dtype=float)
    X, Y, Z = np.meshgrid(g, g, g, indexing='ij')
    dx, dy, dz = X - pos[0], Y - pos[1], Z - pos[2]
    for v in (dx, dy, dz):
        v -= n * np.round(v / n)
    r = np.sqrt(dx * dx + dy * dy + dz * dz)
    m = (r > 0.0) & (r <= rmax)
    rr = r[m]
    u = (dx[m] / rr, dy[m] / rr, dz[m] / rr)
    Er = d[..., 0][m] * u[0] + d[..., 1][m] * u[1] + d[..., 2][m] * u[2]
    return rr, Er, (dx[m], dy[m], dz[m]), u


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--no-corr', required=True)
    ap.add_argument('--corr', required=True)
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
    ap.add_argument('-o', '--output', default='efield_compare.png')
    a = ap.parse_args()

    n, q, eps, kappa = a.ntotal, a.charge, a.epsilon, 1.0 / a.ld
    rr, Eno, dvec, u = radial(a.no_corr, n, a.pos, a.rmax)
    _, Eyes, _, _ = radial(a.corr, n, a.pos, a.rmax)

    off = np.arange(-a.images, a.images + 1) * n
    th = {}
    for kind in ('coulomb', 'dh', 'gauss'):
        acc = np.zeros_like(rr)
        for ox in off:
            for oy in off:
                for oz in off:
                    px, py, pz = dvec[0] + ox, dvec[1] + oy, dvec[2] + oz
                    dd = np.sqrt(px * px + py * py + pz * pz)
                    acc += field_of(kind, dd, a.sigma, kappa, q, eps) * \
                        (px * u[0] + py * u[1] + pz * u[2]) / dd
        th[kind] = acc

    fig, (ax, axe) = plt.subplots(2, 1, figsize=(9.5, 8.6), sharex=True,
                                  gridspec_kw={'height_ratios': [3, 1.15]})
    o = np.argsort(rr)
    rs = rr[o]
    ax.plot(rr, Eno, '.', color='tab:red', ms=2.5, alpha=0.30, label='Lattice, no correction')
    ax.plot(rr, Eyes, '.', color='tab:blue', ms=2.5, alpha=0.30, label='Lattice, with correction')
    ax.plot(rs, th['coulomb'][o], '--', color='dimgray', lw=1.6, label='Pure Coulomb')
    ax.plot(rs, th['dh'][o], '-.', color='tab:green', lw=1.8,
            label=r'Debye-Hückel ($\lambda_D$=' + f'{a.ld:.4g})')
    ax.plot(rs, th['gauss'][o], ':', color='k', lw=2.4,
            label=r'Gaussian cloud, Debye-Hückel ($\sigma$=' + f'{a.sigma:g})  — target')
    if a.rcut:
        for b in (ax, axe):
            b.axvline(a.rcut, color='tab:orange', ls='--', lw=1.2)
        ax.text(a.rcut * 1.02, ax.get_ylim()[1] * 0.25, r'$r_{cut}$', color='tab:orange', fontsize=10)
    ax.set_yscale('log')
    ax.set_ylabel(r'Radial field $kT\,E_r$  (lattice units)', fontsize=12)
    sub = [f'single charge q={q:g}', f'box {n}$^3$']
    if a.kernel:
        sub.append(f'{a.kernel} kernel')
    ax.set_title('Field around one particle, with and without the short-range correction\n'
                 + ',   '.join(sub), fontsize=13)
    ax.grid(True, alpha=0.3, which='both')
    ax.legend(fontsize=9.5)

    for E, c, lab in ((Eno, 'tab:red', 'no correction'), (Eyes, 'tab:blue', 'with correction')):
        axe.plot(rr, 100.0 * (E - th['gauss']) / th['gauss'], '.', color=c, ms=2.5,
                 alpha=0.30, label=lab)
    axe.axhline(0.0, color='k', ls='--', lw=1.0)
    for y in (1.0, -1.0):
        axe.axhline(y, color='tab:gray', ls=':', lw=0.9)
    axe.set_ylim(-35, 35)
    axe.set_xlabel('Distance $r$ from the particle (lattice units)', fontsize=12)
    axe.set_ylabel('relative error vs\nGaussian cloud (%)', fontsize=10)
    axe.grid(True, alpha=0.3)
    axe.legend(fontsize=9, markerscale=3)

    fig.tight_layout()
    fig.savefig(a.output, dpi=150)
    print(f'saved: {a.output}\n')
    print(f"{'r range':>12} {'no corr':>18} {'with corr':>18}")
    for lo, hi in ((0, 2), (2, 4), (4, 6), (6, 8), (8, 12), (12, 18)):
        s = (rr >= lo) & (rr < hi)
        if s.any():
            a1 = np.median(100.0 * (Eno[s] - th['gauss'][s]) / th['gauss'][s])
            a2 = np.median(100.0 * (Eyes[s] - th['gauss'][s]) / th['gauss'][s])
            print(f"  [{lo:2d},{hi:2d})   {a1:+17.2f}% {a2:+17.2f}%")


if __name__ == '__main__':
    main()
