#!/usr/bin/env python3
"""Viscosity test (2026-11-10, user's request): does the chemical-potential floor come from the steady flow?

Runs (dx = 2, r = 4, current exe 4f0c20e7): ../hidro_si (viscosity 0.05, as SLIDES-LD8-hann04-corr_no-dx2/r_04) and
../hidro_si_eta0.5 (viscosity and bulk viscosity 0.5). If the floor of Delta mu comes from the flow, u ~ 1/eta should drop
~10x and Delta mu with it.
Per dump (definitions of microgel/local/scripts/energia_2part.py): F/kT and dF between dumps, Delta mu_a = max - min of
mu_a/kT = ln rho_a + Z_a psi, max|u| and RMS|u| (vel dump, lattice units).
2026-11-10: extended to the 7-viscosity sweep (0.02 ... 2) + summary table of the last dump.
Last dump: RMS of delta mu = mu - <mu> per shell of physical distance to the nearest charge (width 1, as mu_capas_2part.py).
Usage: python3 analizar_visc.py   (from this directory)
"""
import os, sys, glob
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "scripts"))
from energia_2part import rho_fix, Z
from error_capas_2part import nearest_distance

N = 24; DX = 2.0; NH = 4.0; POS = ((6.25, 12.25, 12.25), (4.25, 12.25, 12.25))
RUNS = ("hidro_si_eta0.02", "hidro_si", "hidro_si_eta0.1", "hidro_si_eta0.2", "hidro_si_eta0.5", "hidro_si_eta1.0",
        "hidro_si_eta2.0")   # 2026-11-10: barrido de 7 viscosidades (hidro_si = 0.05)
ETA = (0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0)

if __name__ == "__main__":
    rf = rho_fix(N, NH, POS); k = np.floor((nearest_distance(N, POS) * DX).ravel()).astype(int); last = {}; resumen = {}
    for run in RUNS:
        path = os.path.join(HERE, "..", run); prev = None
        print("\n== %s" % run)
        print(" step |   F/kT              dF (kT)   | dmu_+ (kT)  dmu_- (kT) | max|u|     RMS|u|")
        for f in sorted(glob.glob(os.path.join(path, "psi-0*"))):
            st = int(os.path.basename(f)[4:13]); psi = np.loadtxt(f)
            q = np.loadtxt(os.path.join(path, "qsi-%09d.001-001" % st)); rp, rm = q[:, 0], q[:, 1]
            F = np.sum(rp * (np.log(rp) - 1)) + np.sum(rm * (np.log(rm) - 1)) + 0.5 * np.sum((rp - rm + rf.ravel()) * psi)
            mu = [np.log(rp) + Z[0] * psi, np.log(rm) + Z[1] * psi]
            fv = os.path.join(path, "vel-%09d.001-001" % st)
            if os.path.exists(fv):
                u = np.sqrt((np.loadtxt(fv) ** 2).sum(1)); us = "%.3e  %.3e" % (u.max(), np.sqrt(np.mean(u ** 2)))
            else:
                us = "(no vel)"
            print("%5d | %+.10e %+.3e | %.3e   %.3e  | %s" % (st, F, np.nan if prev is None else F - prev, np.ptp(mu[0]), np.ptp(mu[1]), us))
            prev = F; last[run] = [m - m.mean() for m in mu]
            resumen[run] = (st, F, np.ptp(mu[0]), np.ptp(mu[1]), us)
    print("\nLast dump: RMS delta mu per shell (kT)")
    print("   d   | " + " | ".join("%-25s" % r for r in RUNS))
    print("       | " + " | ".join("%-12s %-12s" % ("mu_+", "mu_-") for r in RUNS))
    for s in range(k.max() + 1):
        m = k == s
        if not m.any(): continue
        print("%5.1f  | " % (s + 0.5) + " | ".join("%.3e    %.3e   " % tuple(np.sqrt(np.mean(x[m] ** 2)) for x in last[r]) for r in RUNS))
    print("\nResumen, ultimo dump (tau = 3 eta/0.8 + 1/2)")
    print("  eta    tau   | step |   F/kT              | dmu_+ (kT)  dmu_- (kT) | max|u|     RMS|u|")
    for e, r in zip(ETA, RUNS):
        st, F, a, b, us = resumen[r]
        print("%5.2f  %6.4f | %4d | %+.10e | %.3e   %.3e  | %s" % (e, 3 * e / 0.8 + 0.5, st, F, a, b, us))
