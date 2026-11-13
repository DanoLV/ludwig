#!/usr/bin/env python3
"""Lattice pinning for every spread/gather kernel available.

Left: linear axis, Gaussian vs Hann-8 (same support, same width, same cost
class). Right: |F_x| on a log axis for all four kernels. The points on nodes
and at the cell midpoint are zero by symmetry (~1e-22, round-off) and are left
out of the log panel.
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
RUNS = [('peskin8', 'Peskin-8 (stretched, support 8)', 'v-', 'tab:brown'),
        ('gauss08', r'Gaussian ($\sigma_{eff}$=1.446, support 8)', 'D-', 'tab:purple'),
        ('hann08', 'Hann-8', 'o-', 'tab:blue'),
        ('peskin6', 'Peskin-6 (support 6)', '<-', 'tab:olive'),
        ('hann06', 'Hann-6', 's-', 'tab:green'),
        ('hann04', 'Hann-4', '^-', 'tab:red')]


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
    fig, (axl, axr) = plt.subplots(1, 2, figsize=(14, 5.2))

    for tag, lab, sty, col in RUNS[:3]:
        x, f = data[tag]
        axl.plot(x, f, sty, color=col, mfc='none', label=lab)
    axl.axhline(0.0, color='k', lw=1.0)
    axl.axvline(24.5, color='gray', ls='--', lw=1.0)
    axl.set_xlabel('Particle position $x$ (lattice units)')
    axl.set_ylabel('$F_x$ on the isolated particle')
    axl.set_title('Support 8: Peskin-8 vs Gaussian vs Hann-8 (linear axis)')
    axl.grid(True, alpha=0.3)
    axl.legend(fontsize=9)

    for tag, lab, sty, col in RUNS:
        x, f = data[tag]
        # Drop the symmetry zeros (round-off) and break the line there, so it
        # does not suggest a nonzero force through the node or the midpoint.
        y = np.where(np.abs(f) > 1e-18, np.abs(f), np.nan)
        axr.semilogy(x, y, sty, color=col, mfc='none', label=lab)
    axr.axvline(24.5, color='gray', ls='--', lw=1.0)
    axr.set_xlabel('Particle position $x$ (lattice units)')
    axr.set_ylabel(r'$|F_x|$ on the isolated particle')
    axr.set_title('All six kernels (log axis)')
    axr.grid(True, alpha=0.3, which='both')
    axr.legend(fontsize=9)

    fig.suptitle(r'Lattice pinning: spurious force on an isolated particle  '
                 r'($\lambda_D$=8, no correction — zero by symmetry)', fontsize=13)
    fig.tight_layout()
    out = f'{BASE}/PIN2GRID-gauss08/pin_to_grid_kernels.png'
    fig.savefig(out, dpi=150)
    print(f'saved: {out}')
    print(f"{'kernel':34} {'max |F_x|':>11}")
    for tag, lab, *_ in RUNS:
        print(f"{lab.replace('$','').replace(chr(92)+'sigma','sigma'):34} {np.abs(data[tag][1]).max():11.2e}")


if __name__ == '__main__':
    main()
