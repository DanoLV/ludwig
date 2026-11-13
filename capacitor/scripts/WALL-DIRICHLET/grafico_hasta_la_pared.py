#!/usr/bin/env python3
"""Free particle (no electrolyte, eta = 0.0005) falling onto the lower grounded wall: full trajectory (steps 0-30000 from
C-particula-libre-hann08-eta0.0005 and the continuation C-particula-libre-hann08-eta0.0005-hastapared until d < 0.1).
Panels: (a) d(t); (b) reduced field on the particle vs d, with the fixed-particle values of C-sin-fondo-hann08 (static);
(c) induced charge on the walls vs d (sum = -q_eff: charge lost when the kernel reaches wall nodes).
Usage: python3 grafico_hasta_la_pared.py -> ../graficos/particula_libre_hasta_la_pared.png
"""
import os, re, glob
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
R1 = ("C-particula-libre-hann08-eta0.0005", 0); R2 = ("C-particula-libre-hann08-eta0.0005-hastapared", 30000)
pat = re.compile(r"Wall charge lower\s+(\S+)\s+upper\s+(\S+)")

def series(bat, t0):
    d = os.path.join(ROOT, bat, "z_09.50"); z = {}
    for f in glob.glob(os.path.join(d, "config.cds0*")):
        st = int(re.search(r"cds(\d+)", f).group(1)); L = [l.strip() for l in open(f)]
        z[st] = (float(L[36 - 1].split()[2]), float(L[37 - 1].split()[2]))
    fz = {}
    for l in open(os.path.join(d, "proceced_data", "particle_force.csv")):
        if l.startswith("#"): continue
        p = l.strip().split(";"); fz[int(p[0])] = float(p[5].replace(",", "."))
    wq = pat.findall(open(os.path.join(d, "run.log"), errors="ignore").read())
    wc = {t0 + 100 * (k + 1): (float(a), float(b)) for k, (a, b) in enumerate(wq)}
    return z, fz, wc

if __name__ == "__main__":
    S = {}
    for bat, t0 in (R1, R2):
        z, fz, wc = series(bat, t0)
        for k, v in z.items(): S.setdefault(k, {})["z"] = v
        for k, v in fz.items(): S.setdefault(k, {})["f"] = v
        for k, v in wc.items(): S.setdefault(k, {})["w"] = v
    st = sorted(k for k in S if "z" in S[k] and k > 0 and "f" in S[k] and "w" in S[k])
    t = np.array(st); d = np.array([S[k]["z"][0] - 1.5 for k in st]); E = np.array([S[k]["f"] for k in st]) * 1e5
    qlo = np.array([S[k]["w"][0] for k in st]); qhi = np.array([S[k]["w"][1] for k in st])
    # static fixed-particle values (C-sin-fondo-hann08): E_z(d) and wall charges
    sd, sE, sQ = [], [], []
    for p in sorted(glob.glob(os.path.join(ROOT, "C-sin-fondo-hann08", "z_*"))):
        dd = float(p.rsplit("_", 1)[1]) - 1.5
        rows = [l for l in open(os.path.join(p, "proceced_data", "particle_force.csv")) if not l.startswith("#")]
        sd.append(dd); sE.append(float(rows[0].strip().split(";")[5].replace(",", ".")) * 1e5)
        w = pat.findall(open(os.path.join(p, "run.log"), errors="ignore").read())[0]; sQ.append(float(w[0]) + float(w[1]))
    plt.rcParams.update({"font.size": 10.5, "axes.grid": True, "grid.color": "#e4e2dc",
                         "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False})
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
    ax[0].plot(t, d, color="#1baf7a", lw=1.8)
    ax[0].axhline(3.5, color="#999", ls="--", lw=1); ax[0].text(t[0], 3.6, "kernel reaches wall nodes (d = 3.5)", fontsize=8.5, color="#666")
    ax[0].axhline(0.0, color="#444", lw=1); ax[0].text(t[0], 0.1, "wall plane (d = 0)", fontsize=8.5, color="#444")
    ax[0].set(xlabel="step", ylabel="distance to the lower wall, d (lattice units)", title="(a) Trajectory")
    ax[1].plot(d, E, color="#1baf7a", lw=1.8, label="moving particle")
    ax[1].plot(sd, sE, "o", ms=4.5, color="#2a78d6", label="fixed particle (static)")
    ax[1].axvline(3.5, color="#999", ls="--", lw=1)
    ax[1].set(xlabel="d (lattice units)", ylabel=r"reduced field on the particle $E_z=F_z\beta/q$", title="(b) Force vs distance"); ax[1].legend(fontsize=9)
    ax[2].plot(d, qlo, color="#2a78d6", lw=1.6, label="lower wall"); ax[2].plot(d, qhi, color="#eb6834", lw=1.6, label="upper wall")
    ax[2].plot(d, qlo + qhi, color="#1baf7a", lw=1.8, label="sum")
    ax[2].plot(sd, sQ, "o", ms=4, mfc="none", color="#1baf7a", label="sum, fixed particle")
    ax[2].axvline(3.5, color="#999", ls="--", lw=1)
    ax[2].set(xlabel="d (lattice units)", ylabel="induced wall charge (units of e)", title="(c) Wall charge (sum = -q_eff)"); ax[2].legend(fontsize=9)
    fig.suptitle("Free particle falling onto the lower grounded wall (no electrolyte, viscosity 0.0005)", y=1.02)
    fig.tight_layout(); fig.savefig(os.path.join(ROOT, "graficos", "particula_libre_hasta_la_pared.png"), dpi=160, bbox_inches="tight")
    for dd in (5.0, 4.0, 3.5, 3.0, 2.5, 2.0, 1.5, 1.0, 0.5, 0.2):
        i = np.argmin(np.abs(d - dd)); print("d=%.2f (paso %d)  E_z=%+.3e  Q_lo=%+.4f Q_hi=%+.4f suma=%+.4f" % (d[i], t[i], E[i], qlo[i], qhi[i], qlo[i] + qhi[i]))
    print("ultimo: paso %d d=%.4f" % (t[-1], d[-1]))
