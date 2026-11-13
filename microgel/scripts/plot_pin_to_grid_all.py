#!/usr/bin/env python3
"""Lattice pinning (pin-to-grid) for all spread/gather kernels, lambda_D = 8, L = 48.

Left: F_x on the isolated particle, linear axis. Right: |F_x|, log axis (zeros
by symmetry at the nodes and at the cell midpoint are round-off and are left out).
"""
import os
import re
import glob

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
STEP = 3000
RUNS = [('peskin4', 'Peskin-4 (support 4)', 'o-', 'tab:red'),
        ('bspline4', 'B-spline-4 (support 4)', 's-', 'tab:orange'),
        ('hann04', 'Hann-4', '^-', 'tab:pink'),
        ('bspline6', 'B-spline-6 (support 6)', 'D-', 'tab:cyan'),
        ('peskin6', 'Peskin-6 (support 6)', '<-', 'tab:olive'),
        ('hann06', 'Hann-6', 'v-', 'tab:green'),
        ('hann08', 'Hann-8 (support 8)', 'o-', 'tab:blue'),
        ('gauss08', r'Gaussian ($\sigma$=1.446, support 8)', 'D-', 'tab:purple'),
        ('peskin8', 'Peskin-8 (stretched, support 8)', 'v-', 'tab:brown')]


def read_fx(path, step=STEP):
    for line in open(path):
        if line.startswith('#'):
            continue
        p = line.strip().replace(',', '.').split(';')
        if len(p) >= 6 and int(p[0]) == step:
            return float(p[3])
    return None


def series(tag):
    pts = []
    for d in glob.glob(f'{BASE}/PIN2GRID-{tag}/x*/'):
        m = re.search(r'x([0-9.]+)/$', d)
        f = read_fx(d + 'proceced_data/particle_force.csv')
        if m and f is not None:
            pts.append((float(m.group(1)), f))
    pts.sort()
    return np.array([p[0] for p in pts]), np.array([p[1] for p in pts])


def main():
    data = {tag: series(tag) for tag, *_ in RUNS}
    fig, (axl, axr) = plt.subplots(1, 2, figsize=(14, 5.4))
    for tag, lab, sty, col in RUNS:
        x, f = data[tag]
        axl.plot(x, f, sty, color=col, mfc='none', ms=5, label=lab)
        y = np.where(np.abs(f) > 1e-18, np.abs(f), np.nan)
        axr.semilogy(x, y, sty, color=col, mfc='none', ms=5, label=lab)
    axl.axhline(0.0, color='k', lw=1.0)
    for ax in (axl, axr):
        ax.axvline(24.5, color='gray', ls='--', lw=1.0)
        ax.set_xlabel('Particle position $x$ (lattice units)')
        ax.grid(True, alpha=0.3, which='both')
    axl.set_ylabel('$F_x$ on the isolated particle')
    axl.set_title('Linear axis')
    axr.set_ylabel(r'$|F_x|$ on the isolated particle')
    axr.set_title('Log axis')
    axr.legend(fontsize=8, ncol=1)
    fig.suptitle(r'Lattice pinning, all kernels: spurious force on an isolated particle '
                 r'($L$=48, $\lambda_D$=8, no correction; zero by symmetry)', fontsize=12)
    fig.tight_layout()
    out = f'{BASE}/graficos/pin_to_grid_all_kernels_LD8.png'
    fig.savefig(out, dpi=150)
    print(f'saved: {out}')

    # Linear zoom without the B-splines (they are 1e2-1e3 times larger and flatten
    # everything else in the left panel above).
    fig2, ax = plt.subplots(figsize=(8, 5.4))
    for tag, lab, sty, col in RUNS:
        if tag.startswith('bspline'):
            continue
        x, f = data[tag]
        ax.plot(x, f, sty, color=col, mfc='none', ms=5, label=lab)
    ax.axhline(0.0, color='k', lw=1.0)
    ax.axvline(24.5, color='gray', ls='--', lw=1.0)
    ax.set_xlabel('Particle position $x$ (lattice units)')
    ax.set_ylabel('$F_x$ on the isolated particle')
    ax.set_title(r'Linear zoom, B-splines excluded ($L$=48, $\lambda_D$=8, no correction)')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    fig2.tight_layout()
    out2 = f'{BASE}/graficos/pin_to_grid_zoom_sin_bspline_LD8.png'
    fig2.savefig(out2, dpi=150)
    print(f'saved: {out2}')

    # Zoom to the Hann-6 scale: only the kernels whose pinning is <= that of Hann-6.
    fig3, ax = plt.subplots(figsize=(8, 5.4))
    for tag, lab, sty, col in RUNS:
        if tag not in ('hann06', 'hann08', 'gauss08', 'peskin8'):
            continue
        x, f = data[tag]
        ax.plot(x, f, sty, color=col, mfc='none', ms=5, label=lab)
    ax.axhline(0.0, color='k', lw=1.0)
    ax.axvline(24.5, color='gray', ls='--', lw=1.0)
    ax.set_xlabel('Particle position $x$ (lattice units)')
    ax.set_ylabel('$F_x$ on the isolated particle')
    ax.set_title(r'Linear zoom at the Hann-6 scale ($L$=48, $\lambda_D$=8, no correction)')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    fig3.tight_layout()
    out3 = f'{BASE}/graficos/pin_to_grid_zoom_hann06_LD8.png'
    fig3.savefig(out3, dpi=150)
    print(f'saved: {out3}')
    print(f"{'kernel':34} {'n':>3} {'max |F_x|':>11}")
    for tag, lab, *_ in RUNS:
        x, f = data[tag]
        print(f"{tag:34} {len(x):3d} {np.abs(f).max() if len(f) else float('nan'):11.2e}")


if __name__ == '__main__':
    main()
