#!/usr/bin/env python3
"""Hann-6 vs Gaussian kernel with the SAME width in fine nodes
(sigma_f = 1.085 = rms width of Hann-6), D3Q27 stencils, no short-range
correction. References: point Coulomb (Ewald) and Gaussian clouds of
physical width sigma_f/r."""
import numpy as np
import scan_27pt as s
import ek_multigrid_poc as poc

SIG_F = 1.085
L, eps, q = 16, 1 / (4 * np.pi), np.array([1.0, -1.0])
ds = [1.0, 1.5, 2.0, 3.0, 4.0]
nplace = 12
rng = np.random.default_rng(7)
geo = {d: [] for d in ds}
for d in ds:
    for _ in range(nplace):
        p1 = 0.5 + L * rng.random(3)
        nh = rng.normal(size=3); nh /= np.linalg.norm(nh)
        geo[d].append(np.array([p1, p1 + d * nh]))
fc = {d: [poc.ewald_forces(p, q, L, eps) for p in geo[d]] for d in ds}

kernels = [("hann6", lambda r: s.hann_1d(6)),
           ("gauss sup6", lambda r: s.gauss_1d(SIG_F, 6.0)),
           ("gauss sup8s", lambda r: s.gauss_1d(SIG_F, 8 * SIG_F))]

print(f"sigma_f = {SIG_F} fine nodes for both kernels; mean (max) relative "
      f"pair-force error over {nplace} placements; d in coarse units\n")
print(f"{'kernel':11s} {'ref':7s} {'r':>2s} " + "".join(f"  d={d:<5}          " for d in ds))
for r in (1, 2, 4):
    sig_c = SIG_F / r
    fg = {d: [s.gauss_pair_reference(p, q, L, eps, sig_c) for p in geo[d]] for d in ds}
    for name, mk in kernels:
        k1d = mk(r)
        F = {d: [s.pm_forces_27(p, q, L, r, eps, k1d) for p in geo[d]] for d in ds}
        for ref, R in (("Coulomb", fc), ("Gauss", fg)):
            row = []
            for d in ds:
                e = [np.linalg.norm(Fm[0] - Fr[0]) / np.linalg.norm(Fr[0])
                     for Fm, Fr in zip(F[d], R[d])]
                row.append(f"{100*np.mean(e):7.2f}% ({100*np.max(e):6.2f}%)")
            print(f"{name:11s} {ref:7s} {r:2d} " + " ".join(row))
    print(f"  (r={r}: physical sigma = {sig_c:.3f}; Gauss/Coulomb force ratio: " +
          ", ".join(f"d={d}: {np.mean([np.linalg.norm(a[0])/np.linalg.norm(b[0]) for a,b in zip(fg[d],fc[d])]):.3f}" for d in ds) + ")\n")
