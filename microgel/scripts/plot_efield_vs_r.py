#!/usr/bin/env python3
"""Radial profile of the lattice electric field around a single subgrid charge,
against theory.

The dump stores kT*E (psi_force.c), which is also the force per unit charge:
F = kt * reunit * q * Esub. So the stored field is compared directly against
the force-per-unit-charge expressions, in the same units used for the
two-particle plots.

Three references:
  - point Coulomb      q/(4 pi eps r^2)
  - point Debye-Huckel q/(4 pi eps) exp(-kr) (k/r + 1/r^2)
  - Gaussian cloud in a Debye-Huckel medium, width sigma (NOT sigma*sqrt(2):
    that factor belongs to the interaction between TWO clouds, while here the
    field of a SINGLE cloud is wanted).
Periodic images are summed as vectors and projected on the radial direction.
"""
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import erfc


def phi_gauss(d, sigma, kappa, q, eps):
    """Screened potential of a Gaussian charge cloud of width sigma."""
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('efield')
    ap.add_argument('-n', '--ntotal', type=int, default=48)
    ap.add_argument('--pos', type=float, nargs=3, default=[24.3, 24.0, 24.0])
    ap.add_argument('-q', '--charge', type=float, default=1.0)
    ap.add_argument('--epsilon', type=float, default=1.0e4)
    ap.add_argument('--ld', type=float, required=True, help='Debye length')
    ap.add_argument('--sigma', type=float, required=True, help='kernel sigma')
    ap.add_argument('--rmax', type=float, default=18.0)
    ap.add_argument('--images', type=int, default=2)
    ap.add_argument('--kernel', default=None)
    ap.add_argument('-o', '--output', default='efield_vs_r.png')
    a = ap.parse_args()

    n, q, eps = a.ntotal, a.charge, a.epsilon
    kappa = 1.0 / a.ld
    d = np.loadtxt(a.efield).reshape(n, n, n, 3)

    # node coordinates: sites sit at 1..n
    g = np.arange(1, n + 1, dtype=float)
    X, Y, Z = np.meshgrid(g, g, g, indexing='ij')
    dx = X - a.pos[0]
    dy = Y - a.pos[1]
    dz = Z - a.pos[2]
    for v in (dx, dy, dz):                      # minimum image
        v -= n * np.round(v / n)
    r = np.sqrt(dx * dx + dy * dy + dz * dz)

    m = (r > 0.0) & (r <= a.rmax)
    rr = r[m]
    ux, uy, uz = dx[m] / rr, dy[m] / rr, dz[m] / rr
    Er = d[..., 0][m] * ux + d[..., 1][m] * uy + d[..., 2][m] * uz

    # theory, images summed as vectors then projected on the primary direction
    off = np.arange(-a.images, a.images + 1) * n
    th = {}
    for kind in ('coulomb', 'dh', 'gauss'):
        acc = np.zeros_like(rr)
        for ox in off:
            for oy in off:
                for oz in off:
                    px, py, pz = dx[m] + ox, dy[m] + oy, dz[m] + oz
                    dd = np.sqrt(px * px + py * py + pz * pz)
                    f = field_of(kind, dd, a.sigma, kappa, q, eps)
                    acc += f * (px * ux + py * uy + pz * uz) / dd
        th[kind] = acc

    fig, (ax, axe) = plt.subplots(2, 1, figsize=(9, 8), sharex=True,
                                  gridspec_kw={'height_ratios': [3, 1]})
    ax.plot(rr, Er, '.', color='tab:blue', ms=2.5, alpha=0.35,
            label='Lattice (all nodes, radial component)')
    o = np.argsort(rr)
    rs = rr[o]
    ax.plot(rs, th['coulomb'][o], 'r--', lw=2, label='Pure Coulomb')
    ax.plot(rs, th['dh'][o], '-.', color='tab:green', lw=2,
            label=r'Debye-Hückel ($\lambda_D$=' + f'{a.ld:g}, {a.images} periodic images/side)')
    ax.plot(rs, th['gauss'][o], ':', color='tab:purple', lw=2.5,
            label=r'Gaussian cloud, Debye-Hückel ($\sigma$=' + f'{a.sigma:g})')
    ax.set_yscale('log')
    ax.set_ylabel(r'Radial field $kT\,E_r$  (lattice units)', fontsize=12)
    sub = [f'single charge q={q:g}', f'box {n}$^3$']
    if a.kernel:
        sub.append(f'{a.kernel} kernel')
    ax.set_title('Field around one particle vs theory\n' + ',   '.join(sub), fontsize=13)
    ax.grid(True, alpha=0.3, which='both')
    ax.legend(fontsize=10)

    err = 100.0 * (Er - th['gauss']) / th['gauss']
    axe.plot(rr, err, '.', color='tab:purple', ms=2.5, alpha=0.35)
    axe.axhline(0.0, color='k', ls='--', lw=1.0)
    for y in (1.0, -1.0):
        axe.axhline(y, color='tab:gray', ls=':', lw=0.9)
    axe.set_ylim(-25, 25)
    axe.set_xlabel('Distance $r$ from the particle (lattice units)', fontsize=12)
    axe.set_ylabel('relative error vs\nGaussian cloud (%)', fontsize=10)
    axe.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(a.output, dpi=150)
    print(f'saved: {a.output}')
    for lo, hi in ((0, 2), (2, 4), (4, 6), (6, 10), (10, 18)):
        s = (rr >= lo) & (rr < hi)
        if s.any():
            print(f'  r in [{lo:2d},{hi:2d}):  median error = {np.median(err[s]):+7.2f} %'
                  f'   spread (p5..p95) = {np.percentile(err[s],5):+7.2f} .. {np.percentile(err[s],95):+7.2f}')


if __name__ == '__main__':
    main()
