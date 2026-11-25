"""Radial profiles of BOTH ion species around the subgrid particle, with salt, vs the discrete PB reference
(LONG-SALT-LD2LB-CEN-3E4, 2026-11-16). Two-species copy of cloud_profile.py.

For each dump s > 1 (qsi/rho/config.cds every 1000):
  * n_+(r)/nbar_+ and n_-(r)/nbar_-, nbar_k = N_k/L^3 (totals from step 1; NP conserves them), bins of 1 LU,
    periodic minimum image around the particle position of that dump (particle frame);
  * PB reference (pb_reference_sal.py) solved at THAT dump's particle position, binned the same way
    (so the reference sees the same particle-vs-grid offsets as the data);
  * time series: excess of n_- and n_+ in r < 3, net ion charge within r < 3 and r < 6
    (Q_ion(R) = sum_{r<R} (n_+ - n_-), screening of the +1 particle), same for PB;
  * far field r > 9: std of dn_k/nbar_k, std of drho/rho and their correlation.
Window average of the profiles over dumps with step >= WIN0 (env, default 15000 = second half of the run;
HEURISTIC choice, stated to the user), also for the per-dump PB profiles. Also the PB profile at the INITIAL
position (as in pb_reference.py).
Writes <run>/cloud_profile_sal.txt (series + profiles) and <run>/cloud_profile_sal.csv (window profiles).
Usage: [WIN0=15000] python3 cloud_profile_sal.py <run dir>
"""
import glob
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pb_reference_sal import L, EDGES, solve_pb, distances, radial, totals_from_run, position  # noqa: E402

d = sys.argv[1]
WIN0 = int(os.environ.get("WIN0", "15000"))
Np, Nm = totals_from_run(d)
NBP, NBM = Np / L ** 3, Nm / L ** 3
nb = len(EDGES) - 1
out = []


def say(s=""):
    print(s)
    out.append(s)


steps = sorted(int(f.split("/qsi-")[-1].split(".")[0]) for f in glob.glob(f"{d}/qsi-0*.001-001"))
steps = [s for s in steps if s > 1]
say(f"=== {d}: dumps {steps[0]}..{steps[-1]}; N+ = {Np:.8f} (nbar {NBP:.6e}), N- = {Nm:.8f} (nbar {NBM:.6e})")
say(" step  particle(x,y,z)          ex-(r<3)  ex+(r<3)  Qion(r<3) Qion(r<6) | PB: ex-(r<3) Qion(r<3) Qion(r<6)"
    " | far: std- std+ std_rho corr-")
prof, pbprof = {}, {}
for s in steps:
    q = np.loadtxt(f"{d}/qsi-{s:09d}.001-001").reshape(L, L, L, 2)
    rho = np.loadtxt(f"{d}/rho-{s:09d}.001-001").reshape(L, L, L)
    r0 = position(d, s)
    dist = distances(r0)
    npl, nmi = q[..., 0], q[..., 1]
    ppb, mpb, it, err = solve_pb(r0, Np, Nm)
    assert err < 1e-12, (s, it, err)
    prof[s] = (radial(npl / NBP, dist), radial(nmi / NBM, dist))
    pbprof[s] = (radial(ppb / NBP, dist), radial(mpb / NBM, dist))
    c3, c6 = dist < 3, dist < 6
    far = dist > 9
    dm, dp, dr = nmi[far] / NBM - 1, npl[far] / NBP - 1, rho[far] / rho.mean() - 1
    say(f"{s:6d}  ({r0[0]:6.2f},{r0[1]:6.2f},{r0[2]:6.2f})  {(nmi[c3] / NBM).mean() - 1:+.4f}  "
        f"{(npl[c3] / NBP).mean() - 1:+.4f}  {(npl - nmi)[c3].sum():+.4f}   {(npl - nmi)[c6].sum():+.4f}  | "
        f"{(mpb[c3] / NBM).mean() - 1:+.4f}    {(ppb - mpb)[c3].sum():+.4f}   {(ppb - mpb)[c6].sum():+.4f}  | "
        f"{dm.std():.4f} {dp.std():.4f} {dr.std():.4f} {np.corrcoef(dm, dr)[0, 1]:+.3f}")

win = [s for s in steps if s >= WIN0]
A = np.array([prof[s] for s in win])          # [dump, species(+,-), bin]
P = np.array([pbprof[s] for s in win])
r0i = position(d, 1)
pi_p, pi_m, _, _ = solve_pb(r0i, Np, Nm)
di = distances(r0i)
P0 = (radial(pi_p / NBP, di), radial(pi_m / NBM, di))
say(f"\n profiles n_k(r)/nbar_k, mean (std) over dumps {win[0]}-{win[-1]} ({len(win)} dumps, WIN0 = {WIN0}, heuristic);"
    " PB at each dump's position (same window) and PB at the initial position")
say("   r bin       |  n- Ludwig        n- PB(dumps) n- PB(r0) |  n+ Ludwig        n+ PB(dumps) n+ PB(r0)")
rows = []
for b in range(nb):
    lm, sm, lp, sp = A[:, 1, b].mean(), A[:, 1, b].std(), A[:, 0, b].mean(), A[:, 0, b].std()
    pm, pp = P[:, 1, b].mean(), P[:, 0, b].mean()
    say(f"   [{EDGES[b]:4.1f},{EDGES[b+1]:4.1f}) |  {lm:.4f} ({sm:.4f})  {pm:.4f}       {P0[1][b]:.4f}    |  "
        f"{lp:.4f} ({sp:.4f})  {pp:.4f}       {P0[0][b]:.4f}")
    rows.append([EDGES[b], EDGES[b + 1], lm, sm, pm, P0[1][b], lp, sp, pp, P0[0][b]])

open(f"{d}/cloud_profile_sal.txt", "w").write("\n".join(out) + "\n")
np.savetxt(f"{d}/cloud_profile_sal.csv", np.array(rows), delimiter=",", fmt="%.6g",
           header="r_lo,r_hi,nminus_ludwig,nminus_ludwig_std,nminus_pb_dumps,nminus_pb_r0,"
                  "nplus_ludwig,nplus_ludwig_std,nplus_pb_dumps,nplus_pb_r0  (n_k/nbar_k; window "
                  f"{win[0]}-{win[-1]})", comments="# ")
print(f"-> {d}/cloud_profile_sal.txt, cloud_profile_sal.csv")
