#!/usr/bin/env python3
"""Electric field around a single subgrid charge.

Reads Ludwig's ascii efield dump (3 components per site, z fastest, x slowest)
and draws the field in the plane through the particle, plus the profile along
the axis through it.

The dump stores kT*E in lattice units (see psi_force.c), so the colour scale is
that stored quantity, not E itself.
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('efield')
    ap.add_argument('-n', '--ntotal', type=int, default=48)
    ap.add_argument('--pos', type=float, nargs=3, default=[24.3, 24.0, 24.0],
                    help='particle position x y z')
    ap.add_argument('--ld', type=float, default=None, help='Debye length, for the label')
    ap.add_argument('--kernel', default=None)
    ap.add_argument('-o', '--output', default='efield_single.png')
    a = ap.parse_args()

    n = a.ntotal
    d = np.loadtxt(a.efield).reshape(n, n, n, 3)   # [x, y, z, component]
    px, py, pz = a.pos
    kz = int(round(pz)) - 1                        # nodes are at 1..n

    ex = d[:, :, kz, 0].T                          # transpose -> rows are y
    ey = d[:, :, kz, 1].T
    mag = np.hypot(ex, ey)
    ax_ = np.arange(1, n + 1)

    fig, (axm, axp) = plt.subplots(
        1, 2, figsize=(13.5, 5.4), gridspec_kw={'width_ratios': [1.15, 1]})

    im = axm.pcolormesh(ax_, ax_, mag, norm=LogNorm(vmin=max(mag.min(), mag.max() * 1e-5),
                                                    vmax=mag.max()),
                        cmap='inferno', shading='nearest')
    axm.streamplot(ax_, ax_, ex, ey, color='white', linewidth=0.6,
                   density=1.1, arrowsize=0.7)
    axm.plot(px, py, 'o', mfc='none', mec='cyan', ms=9, mew=1.6)
    axm.set_xlim(px - 14, px + 14)
    axm.set_ylim(py - 14, py + 14)
    axm.set_aspect('equal')
    axm.set_xlabel('$x$ (lattice units)')
    axm.set_ylabel('$y$ (lattice units)')
    axm.set_title(f'In-plane field magnitude, slice $z$ = {pz:g}')
    fig.colorbar(im, ax=axm, label=r'$|kT\,\mathbf{E}_{xy}|$  (lattice units)')

    # profile along x through the particle
    jy = int(round(py)) - 1
    line = d[:, jy, kz, 0]
    axp.plot(ax_, line, 'o-', color='tab:blue', ms=4, label=r'$kT\,E_x$ on the lattice')
    axp.axvline(px, color='tab:red', linestyle='--', linewidth=1.2,
                label=f'particle at $x$ = {px:g}')
    axp.axhline(0.0, color='k', linewidth=1.0)
    axp.set_xlim(px - 16, px + 16)
    axp.set_xlabel('$x$ (lattice units)')
    axp.set_ylabel(r'$kT\,E_x$  (lattice units)')
    axp.set_title(f'Field along the axis  ($y$ = {py:g}, $z$ = {pz:g})')
    axp.grid(True, alpha=0.3)
    axp.legend(fontsize=9)

    sub = [f'single subgrid charge $q$ = +1', f'box {n}$^3$']
    if a.ld:
        sub.append(r'$\lambda_D$ = ' + f'{a.ld:g}')
    if a.kernel:
        sub.append(f'{a.kernel} kernel')
    fig.suptitle('Electric field around one particle  —  ' + ',   '.join(sub), fontsize=13)
    fig.tight_layout()
    fig.savefig(a.output, dpi=150)
    print(f'saved: {a.output}')
    print(f'|kT E| max in plane = {mag.max():.4e}')
    print(f'E_x zero crossing near x = '
          f'{np.interp(0.0, [line[int(round(px))-1], line[int(round(px))]], [px, px+1]):.3f}'
          if False else '')


if __name__ == '__main__':
    main()
