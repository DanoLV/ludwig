#!/usr/bin/env python3
"""Two-particle force error for two kernels, with and without the pm_sr
correction, on the same axes. Same reference as compare_force_two_runs.py
(Gaussian-Gaussian Debye-Huckel, 3D periodic image sum); points where the
target is ~0 by symmetry are dropped."""
import argparse, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.special import erfc
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compare_force_two_runs import series
from compare_force_coulomb import parse_physics, find_sim_dirs

ap = argparse.ArgumentParser()
ap.add_argument('--a-no', required=True); ap.add_argument('--a-yes', required=True)
ap.add_argument('--b-no', required=True); ap.add_argument('--b-yes', required=True)
ap.add_argument('--a-label', default='kernel A'); ap.add_argument('--b-label', default='kernel B')
ap.add_argument('--sigma', type=float, required=True); ap.add_argument('--rcut', type=float, default=None)
ap.add_argument('-o', '--output', required=True)
a = ap.parse_args()

ph = parse_physics(os.path.join(find_sim_dirs(a.a_yes)[0], 'logs', 'output.txt'))
eps, kappa, L = ph['epsilon'], ph['kappa'], ph['box_l']
sp = a.sigma*np.sqrt(2.0)
def phi(d):
    B = sp*kappa/np.sqrt(2); C = 1/(sp*np.sqrt(2))
    return 1/(8*np.pi*eps)*np.exp(kappa**2*sp**2/2)*(np.exp(-kappa*d)*erfc(B-C*d)-np.exp(kappa*d)*erfc(B+C*d))/d
F = lambda d, h=1e-4: -(phi(d+h)-phi(d-h))/(2*h)
def target(r):
    acc = np.zeros_like(r)
    for i in range(-2, 3):
        for j in range(-2, 3):
            for k in range(-2, 3):
                dx = r+i*L; dd = np.sqrt(dx*dx+(j*L)**2+(k*L)**2); acc += F(dd)*dx/dd
    return acc

sets = [(a.a_no, f'{a.a_label}, no correction', 'o', 'tab:blue', 'none'),
        (a.a_yes, f'{a.a_label}, with correction', 'o', 'tab:blue', 'tab:blue'),
        (a.b_no, f'{a.b_label}, no correction', 'D', 'tab:purple', 'none'),
        (a.b_yes, f'{a.b_label}, with correction', 'D', 'tab:purple', 'tab:purple')]

def coulomb(d): return 1/(4*np.pi*eps*d**2)
def dh(d): return 1/(4*np.pi*eps)*np.exp(-kappa*d)*(kappa/d + 1/d**2)
def image_sum(fn, r):
    acc = np.zeros_like(r)
    for i in range(-2, 3):
        for j in range(-2, 3):
            for k in range(-2, 3):
                dx = r+i*L; dd = np.sqrt(dx*dx+(j*L)**2+(k*L)**2); acc += fn(dd)*dx/dd
    return acc

data = {}
for d, lab, *_ in sets:
    r, f = series(d, 3000, L)
    t = target(r); keep = np.abs(t) > 1e-6*np.abs(t).max()
    data[lab] = (r[keep], f[keep], 100*(f[keep]/t[keep]-1))

rall = np.concatenate([v[0] for v in data.values()])
rf = np.linspace(rall.min(), rall.max(), 400)
tf = target(rf); floor = 1e-6*np.abs(tf).max()
blank = lambda y: np.where(np.abs(tf) > floor, y, np.nan)

fig = plt.figure(figsize=(13, 10))
gsp = fig.add_gridspec(2, 2, height_ratios=[1.6, 1], hspace=0.28, wspace=0.18)
ax = fig.add_subplot(gsp[0, :]); ae = fig.add_subplot(gsp[1, 0]); az = fig.add_subplot(gsp[1, 1])

ax.plot(rf, blank(image_sum(coulomb, rf)), '--', color='dimgray', lw=1.5, label='Pure Coulomb')
ax.plot(rf, blank(image_sum(dh, rf)), '-.', color='tab:green', lw=1.7,
        label=r'Debye-Hückel ($\lambda_D$=' + f'{1/kappa:.4g})')
ax.plot(rf, blank(tf), ':', color='k', lw=2.4,
        label=r'Gaussian Debye-Hückel ($\sigma$=' + f'{a.sigma:g})  — target')
for d, lab, mk, col, mfc in sets:
    r, f, e = data[lab]
    ax.plot(r, f, mk, color=col, mfc=mfc, ms=7, mew=1.5, ls='none', label=f'Simulation: {lab}')
    for axis in (ae, az):
        axis.plot(r, e, mk + ('-' if mfc != 'none' else '--'), color=col, mfc=mfc, ms=5, lw=1.3, label=lab)

ax.set_xscale('log'); ax.set_yscale('log')
ax.set_ylabel(r'Force $F_x$ on the moving particle', fontsize=12)
ax.grid(True, alpha=0.3, which='both'); ax.legend(fontsize=9, ncol=2, loc='lower left')
if a.rcut: ax.axvline(a.rcut, color='tab:orange', ls='--', lw=1.2)

for axis in (ae, az):
    axis.axhline(0, color='k', ls='--', lw=1)
    if a.rcut: axis.axvline(a.rcut, color='tab:orange', ls='--', lw=1.2)
    axis.set_xscale('log')
    axis.set_xlabel('Separation r between particles (lattice units)')
    axis.grid(True, alpha=0.3, which='both')
ae.set_ylabel('relative error vs\nGaussian-Gaussian (%)'); ae.set_title('Relative error, full scale')
az.set_ylim(-1.5, 1.5); az.set_title('Zoom on the corrected runs')
for y in (1, -1): az.axhline(y, color='tab:gray', ls=':', lw=0.9)
ae.legend(fontsize=8)

fig.suptitle(f'Two-particle force: {a.a_label} vs {a.b_label}\n'
             + r'$\lambda_D$ = ' + f'{1/kappa:.4g},   target ' + r'$\sigma$ = ' + f'{a.sigma:g},   step 3000',
             fontsize=13)
fig.savefig(a.output, dpi=150, bbox_inches='tight')
print(f'saved: {a.output}')
