"""Radial counterion profile around the subgrid particle (particle frame).

For each dump step s (psi/rho every 500): read counterion density n (qsi column 1),
fluid density rho, particle position (config.cds at step s, colloid_io_freq 100).
Distances: periodic minimum image in a cubic box L=27, nodes at integer coords 1..L.
Outputs (text): profile n(r)/nbar in bins of 1 LU, averaged over a time window;
time series of the mean excess in r<3 and the excess counterion number within r<R;
far-field (r>9) ion vs fluid statistics.
"""
import sys, glob, os
import numpy as np

L = 27
NBAR = 5.0805263425290857e-05   # psi_electroneutral rhoi (run.log DIAG line)
EDGES = np.arange(0.0, 14.0, 1.0)

i = np.arange(L**3)
X = np.array([i // (L*L) + 1, (i // L) % L + 1, i % L + 1]).T.astype(float)


def load(d, s):
    q = np.loadtxt(f"{d}/qsi-{s:09d}.001-001")[:, 1]
    rho = np.loadtxt(f"{d}/rho-{s:09d}.001-001")
    r0 = np.array([float(v) for v in open(f"{d}/config.cds{s:08d}.001-001").read().split("\n")[35].split()])
    dv = X - r0
    dv -= L*np.round(dv/L)
    return q, rho, np.linalg.norm(dv, axis=1), r0


def profile(n, dist):
    idx = np.digitize(dist, EDGES) - 1
    return np.array([n[idx == k].mean() for k in range(len(EDGES) - 1)])


for d in sys.argv[1:]:
    steps = sorted(int(f.split("/qsi-")[-1].split(".")[0]) for f in glob.glob(f"{d}/qsi-0*.001-001"))
    steps = [s for s in steps if s > 1]
    print(f"\n=== {d}: dumps {steps[0]}..{steps[-1]}")
    profs = {}
    print(" step  particle(x,y,z)        excess<3   Nex(r<3)  Nex(r<6)  far: ion_std fluid_std corr   total/nbar")
    for s in steps:
        q, rho, dist, r0 = load(d, s)
        n = q/NBAR
        dn = n - 1.0
        dr = rho/rho.mean() - 1.0
        far = dist > 9
        profs[s] = profile(n, dist)
        print(f"{s:6d}  ({r0[0]:6.2f},{r0[1]:6.2f},{r0[2]:6.2f})  {dn[dist<3].mean():+.4f}  "
              f"{(q[dist<3]-NBAR).sum():+.4f}  {(q[dist<6]-NBAR).sum():+.4f}   "
              f"{dn[far].std():.4f}   {dr[far].std():.4f}  {np.corrcoef(dn[far], dr[far])[0,1]:+.3f}  {n.mean():.6f}")
    win = [s for s in steps if s >= int(os.environ.get("WIN0", "5000"))]  # averaging window start (env WIN0)
    avg = np.mean([profs[s] for s in win], axis=0)
    sd = np.std([profs[s] for s in win], axis=0)
    print(f" profile n(r)/nbar, mean (std) over steps {win[0]}-{win[-1]} ({len(win)} dumps):")
    for k in range(len(EDGES) - 1):
        print(f"   r in [{EDGES[k]:4.1f},{EDGES[k+1]:4.1f}):  {avg[k]:.4f}  ({sd[k]:.4f})")
