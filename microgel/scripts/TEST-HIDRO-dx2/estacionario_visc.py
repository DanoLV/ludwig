#!/usr/bin/env python3
"""Steady-state check of the 4 x 7 viscosity runs (2026-11-16, user's request).

For each variant (old, centred, flla, flla_mass) and viscosity, the time series of Delta mu_a = max - min of
mu_a/kT = ln rho_a + Z_a psi (as analizar_visc.py) is built from the dumps every 250 steps of the original run
(0..3000) and of its continuations cont6000/, cont9000/ (if present).

Steady-state criterion (HEURISTIC, proposed by Claude, not a standard; stated to the user):
    r_a = |Delta mu_a(N) - Delta mu_a(N - 1000)| / Delta mu_a(N) < TOL = 1 %   for both species,
with N the last available dump. Printed: Delta mu_a at 3000, at N-1000 and at N, r_a, and STEADY yes/no.
Writes nothing. Usage: python3 scripts/estacionario_visc.py
"""
import glob, os, sys
import numpy as np

BASE = "/home/danolv/Sim/ludwig/microgel/local/TEST-HIDRO-dx2"
sys.path.insert(0, os.path.join(BASE, "scripts"))
from analizar_visc import RUNS, ETA, Z

VARIANTS = {"old": BASE, "centred": BASE + "/advfix", "flla": BASE + "/advfix_flla",
            "flla_mass": BASE + "/advfix_flla_mass"}
TOL, LAG = 0.01, 1000


def series(v, r):
    dirs = [os.path.join(VARIANTS[v], r)] + sorted(glob.glob(os.path.join(BASE, "cont*", v, r)),
                                                    key=lambda d: int(d.split("/cont")[1].split("/")[0]))
    out = {}
    for d in dirs:
        for f in glob.glob(os.path.join(d, "psi-0*")):
            st = int(os.path.basename(f)[4:13])
            if st < 250 or st in out:
                continue
            q = np.loadtxt(os.path.join(d, "qsi-%09d.001-001" % st)); psi = np.loadtxt(f)
            out[st] = [np.ptp(np.log(q[:, a]) + Z[a] * psi) for a in (0, 1)]
    st = np.array(sorted(out)); return st, np.array([out[s] for s in st])


def check(v, r):
    st, d = series(v, r)
    n = st[-1]; i_n = -1; i_l = int(np.where(st == n - LAG)[0][0])
    rel = np.abs(d[i_n] - d[i_l]) / d[i_n]
    return st, d, n, rel, bool((rel < TOL).all())


if __name__ == "__main__":
    print("criterion: |dmu(N) - dmu(N-%d)|/dmu(N) < %.0f %% for both species (heuristic)" % (LAG, 100 * TOL))
    for v in VARIANTS:
        print("\n== %s" % v)
        print("  eta  |  N   | dmu_+ 3000  N-1000     N         r_+     | dmu_- 3000  N-1000     N         r_-     | steady")
        for e, r in zip(ETA, RUNS):
            st, d, n, rel, ok = check(v, r)
            i3, il = int(np.where(st == 3000)[0][0]), int(np.where(st == n - LAG)[0][0])
            print(" %5.2f | %4d | %.3e  %.3e  %.3e  %6.2f%% | %.3e  %.3e  %.3e  %6.2f%% | %s" % (
                e, n, d[i3, 0], d[il, 0], d[-1, 0], 100 * rel[0], d[i3, 1], d[il, 1], d[-1, 1], 100 * rel[1],
                "yes" if ok else "NO"))
