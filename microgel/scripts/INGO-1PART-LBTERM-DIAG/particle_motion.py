"""Mean velocity / position of the free subgrid particle under thermal LB (2026-11-12).
Reads config.cds* (line 36: r, line 37: v = displacement per step). Unwraps periodic jumps.
"""
import glob, sys
import numpy as np
L = 27

def traj(d):
    fs = sorted(glob.glob(f"{d}/config.cds0*.001-001"))
    st = np.array([int(f.split("config.cds")[1][:8]) for f in fs])
    R, V = [], []
    for f in fs:
        t = open(f).read().split("\n")
        R.append([float(x) for x in t[35].split()]); V.append([float(x) for x in t[36].split()])
    R, V = np.array(R), np.array(V)
    dR = np.diff(R, axis=0); dR -= L*np.round(dR/L)
    Ru = np.vstack([R[:1], R[0] + np.cumsum(dR, axis=0)])
    return st, Ru, V

# per-step velocities (BASE-MASS)
st, R, V = traj("BASE-MASS")
n = len(V)
print(f"BASE-MASS: {n} steps; mean v = {V.mean(0)} ; std v = {V.std(0)}")
# blocking estimate of the standard error of the mean (accounts for correlations)
for b in [1, 10, 50, 100]:
    m = n//b; blk = V[:m*b].reshape(m, b, 3).mean(1)
    print(f"  block {b:4d}: s.e.(mean v) = {blk.std(0, ddof=1)/np.sqrt(m)}")
# velocity autocorrelation (x,y,z averaged), normalised
vc = V - V.mean(0)
acf = [np.mean(np.sum(vc[:n-k]*vc[k:], 1))/np.mean(np.sum(vc*vc, 1)) for k in range(0, 6)]
print("  normalised velocity autocorrelation, lags 0..5:", np.round(acf, 3))

# long trajectory every 100 steps
for d in ["LONG-MASS", "LONG-MASS-D10", "LONG-MASS-D30"]:
    st, R, V = traj(d)
    print(f"{d}: steps {st[0]}..{st[-1]}, net displacement {R[-1]-R[0]} (|.|={np.linalg.norm(R[-1]-R[0]):.2f})")
st, R, V = traj("LONG-MASS-D10")
dt = st[1] - st[0]
print("MSD (time-averaged over the single trajectory) -> D = MSD/(6 t):")
for lag in [1, 2, 5, 10, 20, 50]:
    dd = R[lag:] - R[:-lag]
    msd = np.mean(np.sum(dd**2, 1))
    print(f"  t = {lag*dt:5d}: MSD = {msd:.3f}, D = {msd/(6*lag*dt):.2e}, independent windows ~ {len(R)//lag}")
