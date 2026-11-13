#!/usr/bin/env python3
"""Chemical-potential residual around the charges, two fixed charges r = 4, three lattices (2026-11-10, user's request).

Same runs and definitions as energia_2part.py: mu_a/kT = ln rho_a + Z_a psi (qsi = fluid only, psi in kT/e).
delta mu_a,n = mu_a,n - <mu_a>, <.> = plain mean over all nodes (choice of this script; the residual max - min of
energia_2part.py does not depend on it).
Per shell of physical distance d to the NEAREST charge (minimum image, width 1 physical, as error_capas_2part.py):
signed mean of delta mu, RMS and max |delta mu| (kT), for each species.
Dumps: physical time t = 5000 (steps 1250, 5000, 20000 for dx = 2, 1, 1/2: same time in the three lattices) and the last
dump of each run (dx = 2: t = 12000, dx = 1: t = 10000; dx = 1/2: t = 5000, the same).
Usage: python3 mu_capas_2part.py   (from microgel/local/scripts) -> tables printed
"""
import os
import numpy as np
from energia_2part import RUNS, HERE, Z
from error_capas_2part import nearest_distance

DUMPS = {2.0: (1250, 3000), 1.0: (5000, 10000), 0.5: (20000,)}

if __name__ == "__main__":
    for dx, rel, N, n, pos in RUNS:
        path = os.path.join(HERE, "..", rel); d = (nearest_distance(N, pos) * dx).ravel(); k = np.floor(d).astype(int)
        for st in DUMPS[dx]:
            psi = np.loadtxt(os.path.join(path, "psi-%09d.001-001" % st))
            q = np.loadtxt(os.path.join(path, "qsi-%09d.001-001" % st))
            dm = [np.log(q[:, a]) + Z[a] * psi for a in (0, 1)]; dm = [m - m.mean() for m in dm]
            print("\n== dx = %g, step %d (t = %g): max - min  mu_+ %.3e  mu_- %.3e (kT)"
                  % (dx, st, st * dx * dx, np.ptp(dm[0]), np.ptp(dm[1])))
            print("   d    nodes |  mean dmu_+   RMS dmu_+   max|dmu_+| |  mean dmu_-   RMS dmu_-   max|dmu_-|")
            for s in range(k.max() + 1):
                m = k == s
                if not m.any(): continue
                v = [(x[m].mean(), np.sqrt(np.mean(x[m] ** 2)), np.abs(x[m]).max()) for x in dm]
                print("%5.1f %7d | %+.3e  %.3e  %.3e | %+.3e  %.3e  %.3e" % ((s + 0.5), m.sum(), *v[0], *v[1]))
