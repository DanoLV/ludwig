#!/usr/bin/env python3
"""Two-particle force: corrected and uncorrected runs on the same axes.

Both runs must be a set of pos_* directories at matching separations. The
theory curves and the reference are shared, so the only difference between
the two point sets is pm_sr_correction.

Reuses the readers and the screened-Gaussian potential of
compare_force_coulomb.py, so the numbers match those single-run figures.
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import erfc

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compare_force_coulomb import (find_sim_dirs, parse_physics,      # noqa: E402
                                   read_force_csv, read_colloid_positions)


def series(base, step, box_l):
    """(r, Fx) for the moving particle in every pos_* directory of `base`."""
    rows = []
    for d in find_sim_dirs(base):
        fcsv = os.path.join(d, 'proceced_data', 'particle_force.csv')
        if not os.path.isfile(fcsv):
            continue
        _, forces = read_force_csv(fcsv, step)
        _, pos = read_colloid_positions(d, step)
        if len(pos) < 2 or not forces:
            continue
        # The fixed particle always sits at lower x and the moving one at
        # higher x, as in compare_force_coulomb.py. Picking by index instead
        # of by x silently returns the wrong particle and flips the sign.
        ids = list(pos)
        fixed = min(ids, key=lambda i: pos[i][0])
        moving = max(ids, key=lambda i: pos[i][0])
        if moving not in forces:
            continue
        rows.append((pos[moving][0] - pos[fixed][0], forces[moving][0]))
    rows.sort()
    return np.array([r for r, _ in rows]), np.array([f for _, f in rows])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--no-corr', required=True, help='base dir of the uncorrected run')
    ap.add_argument('--corr', required=True, help='base dir of the corrected run')
    ap.add_argument('-n', '--step', type=int, default=3000)
    ap.add_argument('-q', '--charge', type=float, default=1.0)
    ap.add_argument('--gaussian-sigma', type=float, required=True)
    ap.add_argument('--rcut', type=float, default=None)
    ap.add_argument('--images', type=int, default=2)
    ap.add_argument('--kernel', default=None)
    ap.add_argument('-o', '--output', required=True)
    a = ap.parse_args()

    ph = parse_physics(os.path.join(find_sim_dirs(a.corr)[0], 'logs', 'output.txt'))
    eps, kappa, box_l = ph['epsilon'], ph['kappa'], ph['box_l']
    q = a.charge
    sp = a.gaussian_sigma * np.sqrt(2.0)          # two clouds -> sigma*sqrt(2)

    r_no, f_no = series(a.no_corr, a.step, box_l)
    r_ye, f_ye = series(a.corr, a.step, box_l)

    def phi_g(d):
        B = sp * kappa / np.sqrt(2.0)
        C = 1.0 / (sp * np.sqrt(2.0))
        t1 = np.exp(-kappa * d) * erfc(B - C * d)
        t2 = np.exp(kappa * d) * erfc(B + C * d)
        return q * q / (8 * np.pi * eps) * np.exp(kappa ** 2 * sp ** 2 / 2.0) * (t1 - t2) / d

    def iso(kind, d, h=1e-4):
        if kind == 'coulomb':
            return q * q / (4 * np.pi * eps * d ** 2)
        if kind == 'dh':
            return q * q / (4 * np.pi * eps) * np.exp(-kappa * d) * (kappa / d + 1 / d ** 2)
        return -(phi_g(d + h) - phi_g(d - h)) / (2 * h)

    def periodic(kind, r):
        acc = np.zeros_like(r)
        for i in range(-a.images, a.images + 1):
            for j in range(-a.images, a.images + 1):
                for k in range(-a.images, a.images + 1):
                    dxs = r + i * box_l
                    dd = np.sqrt(dxs ** 2 + (j * box_l) ** 2 + (k * box_l) ** 2)
                    acc += iso(kind, dd) * dxs / dd
        return acc

    rf = np.linspace(min(r_no.min(), r_ye.min()), max(r_no.max(), r_ye.max()), 400)

    # At r = L/2 the force vanishes by symmetry (the moving particle is
    # equidistant from the fixed one and its periodic image), so both the
    # simulated and the theoretical force are numerical zeros there and their
    # ratio is noise. Drop any point whose target is below 1e-6 of the peak,
    # as compare_force_coulomb.py already does, and blank the theory curves
    # below the same floor so the log axis is not stretched to ~1e-23.
    g_fine = periodic('gauss', rf)
    floor = 1e-6 * np.max(np.abs(g_fine))
    def keep(r):
        return np.abs(periodic('gauss', r)) > floor
    k_no, k_ye = keep(r_no), keep(r_ye)
    dropped = sorted(set(np.r_[r_no[~k_no], r_ye[~k_ye]].tolist()))
    r_no, f_no = r_no[k_no], f_no[k_no]
    r_ye, f_ye = r_ye[k_ye], f_ye[k_ye]
    def blank(y):
        y = np.array(y, dtype=float)
        y[np.abs(g_fine) <= floor] = np.nan
        return y
    fig, (ax, axe) = plt.subplots(2, 1, figsize=(9, 8.6), sharex=True,
                                  gridspec_kw={'height_ratios': [3, 1.15]})
    ax.plot(rf, blank(periodic('coulomb', rf)), '--', color='dimgray', lw=1.6, label='Pure Coulomb')
    ax.plot(rf, blank(periodic('dh', rf)), '-.', color='tab:green', lw=1.8,
            label=r'Debye-Hückel ($\lambda_D$=' + f'{1/kappa:.4g}, {a.images} images/side)')
    ax.plot(rf, blank(periodic('gauss', rf)), ':', color='k', lw=2.4,
            label=r'Gaussian Debye-Hückel ($\sigma$=' + f'{a.gaussian_sigma:g}, '
                  + r'$\sigma_{pair}$=' + f'{sp:.3g})  — target')
    ax.plot(r_no, f_no, 'o', color='tab:red', ms=7, mfc='none', mew=1.6,
            label='Simulation, no correction')
    ax.plot(r_ye, f_ye, 'o', color='tab:blue', ms=7, label='Simulation, with correction')
    if a.rcut:
        for b in (ax, axe):
            b.axvline(a.rcut, color='tab:orange', ls='--', lw=1.2)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylabel(r'Force $F_x$ on the moving particle', fontsize=12)
    sub = [r'$\lambda_D$ = ' + f'{1/kappa:.4g}']
    if a.kernel:
        sub.append(f'{a.kernel} kernel')
    sub.append(r'$\sigma$ = ' + f'{a.gaussian_sigma:g}')
    ax.set_title(f'Two-particle force, with and without the short-range correction '
                 f'(step={a.step})\n' + ',   '.join(sub) + f'\nq={q:g}, epsilon={eps:.4g}',
                 fontsize=13)
    ax.grid(True, alpha=0.3, which='both')
    ax.legend(fontsize=9.5)

    for r, f, c, lab in ((r_no, f_no, 'tab:red', 'no correction'),
                         (r_ye, f_ye, 'tab:blue', 'with correction')):
        g = periodic('gauss', r)
        axe.plot(r, 100 * (f - g) / g, 'o-', color=c, ms=5, lw=1.2, label=lab)
    axe.axhline(0.0, color='k', ls='--', lw=1.0)
    for y in (1.0, -1.0):
        axe.axhline(y, color='tab:gray', ls=':', lw=0.9)
    axe.set_xscale('log')
    axe.set_xlabel('Separation r between particles (lattice units)', fontsize=12)
    axe.set_ylabel('relative error vs\nGaussian-Gaussian (%)', fontsize=10)
    axe.grid(True, alpha=0.3, which='both')
    axe.legend(fontsize=9)

    fig.tight_layout()
    fig.savefig(a.output, dpi=150, bbox_inches='tight')
    print(f'saved: {a.output}')
    if dropped:
        print(f'  omitted (force ~0 by symmetry, ratio is noise): r = {dropped}')
    for r, f, lab in ((r_no, f_no, 'no corr  '), (r_ye, f_ye, 'with corr')):
        g = periodic('gauss', r)
        e = 100 * (f - g) / g
        print(f'  {lab}: err(r=1) = {e[0]:+7.2f} %   |err|max = {np.abs(e).max():6.2f} %')


if __name__ == '__main__':
    main()
