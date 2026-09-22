#!/usr/bin/env python3
"""Lattice pinning: spurious force on an isolated particle as it is moved
across one lattice cell. By symmetry the force on a single charge in an
otherwise uniform periodic box must vanish; whatever is left measures how
strongly the spread/gather kernel couples to the mesh.

Reads PIN2GRID-hann{04,06,08}/x*/proceced_data/particle_force.csv.
"""
import glob
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
STEP = 3000
STYLE = {'08': ('o-', 'tab:blue'), '06': ('s-', 'tab:green'), '04': ('^-', 'tab:red')}


def read_fx(path, step=STEP):
    """Last Fx recorded at `step` (the file uses ';' and a decimal comma)."""
    for line in open(path):
        if line.startswith('#'):
            continue
        p = line.strip().replace(',', '.').split(';')
        if len(p) >= 6 and int(p[0]) == step:
            return float(p[3])
    return None


def main():
    fig, ax = plt.subplots(figsize=(9, 5))
    for order in ('08', '06', '04'):
        xs, fs = [], []
        for d in sorted(glob.glob(f'{BASE}/PIN2GRID-hann{order}/x*/')):
            m = re.search(r'x([0-9.]+)/$', d)
            f = read_fx(d + 'proceced_data/particle_force.csv')
            if m and f is not None:
                xs.append(float(m.group(1)))
                fs.append(f)
        if not xs:
            continue
        o = sorted(zip(xs, fs))
        style, colour = STYLE[order]
        ax.plot([p[0] for p in o], [p[1] for p in o], style, color=colour,
                markerfacecolor='none', label=f'hann{int(order)}')

    ax.axhline(0.0, color='k', linewidth=1.2)
    for xv in (24.0, 25.0):
        ax.axvline(xv, color='gray', linestyle=':', linewidth=1.0)
    ax.axvline(24.5, color='gray', linestyle='--', linewidth=1.0)
    ax.text(24.02, ax.get_ylim()[1] * 0.93, 'lattice node', fontsize=9, color='gray')
    ax.text(24.52, ax.get_ylim()[1] * 0.93, 'cell midpoint', fontsize=9, color='gray')

    ax.set_xlabel('Particle position $x$ (lattice units)', fontsize=12)
    ax.set_ylabel('$F_x$ on the isolated particle', fontsize=12)
    ax.set_title('Lattice pinning: spurious force on an isolated particle\n'
                 r'($\lambda_D$=8, no correction — zero by symmetry)', fontsize=13)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=11)
    fig.tight_layout()
    out = f'{BASE}/PIN2GRID-hann08/pin_to_grid.png'
    fig.savefig(out, dpi=150)
    print(f'saved: {out}')


if __name__ == '__main__':
    main()
