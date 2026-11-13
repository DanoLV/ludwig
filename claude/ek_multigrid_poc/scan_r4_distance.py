#!/usr/bin/env python3
"""Pair-force error vs separation for r = 1, 2, 4 and several kernels
(kernel width defined in FINE nodes). Reference: Ewald, point charges.
Question: with r = 4, can pm_sr be dropped for d >= 2 coarse?"""
import math, numpy as np
import ek_multigrid_poc as poc

def hann(n):
    def f(x):
        return np.where(np.abs(x) < n / 2, (1 + np.cos(2 * np.pi * x / n)) / n, 0.0)
    return f, n / 2

KERNELS = {"peskin4": (poc.peskin4, 2.0), "hann4": hann(4), "hann6": hann(6)}

def make_kernel_1d(phi, half):
    def k1d(xf, s, n):
        m = np.arange(math.floor(xf - half * s) + 1, math.floor(xf + half * s) + 1)
        return (m - 1) % n, phi((xf - m) / s) / s
    return k1d

L, eps, q = 16, 1 / (4 * np.pi), np.array([1.0, -1.0])
rng = np.random.default_rng(7)
nplace = 12
ds = [1.0, 1.5, 2.0, 3.0, 4.0]
plac = {}
for d in ds:
    plac[d] = []
    for _ in range(nplace):
        p1 = 0.5 + L * rng.random(3)
        nh = rng.normal(size=3); nh /= np.linalg.norm(nh)
        pos = np.array([p1, p1 + d * nh])
        plac[d].append((pos, poc.ewald_forces(pos, q, L, eps)))

print("relative pair-force error |F_PM-F_Ewald|/|F_Ewald|: mean (max) over",
      nplace, "placements; separations in coarse units")
print(f"{'kernel':8s} {'r':>2s} " + "".join(f"   d={d:<4}        " for d in ds))
for name, (phi, half) in KERNELS.items():
    poc.kernel_1d = make_kernel_1d(phi, half)
    for r in (1, 2, 4):
        row = []
        for d in ds:
            e = []
            for pos, Fe in plac[d]:
                F = poc.pm_forces(pos, q, L, r, eps, 1)
                e.append(np.linalg.norm(F[0] - Fe[0]) / np.linalg.norm(Fe[0]))
            row.append(f"{100*np.mean(e):6.2f}% ({100*np.max(e):5.2f}%)")
        print(f"{name:8s} {r:2d} " + " ".join(row))
