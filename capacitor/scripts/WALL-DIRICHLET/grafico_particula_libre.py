#!/usr/bin/env python3
"""Free charged particle between grounded plates (L = 32, Hann-8, q = +1, eta = 0.0005, initial d = 8 from the lower plate).
Three runs: no electrolyte; electrolyte (lambda_D = 6) started from scratch (uniform ions, particle free from step 0);
electrolyte (lambda_D = 6) released from the equilibrated state (restart at step 5000 of the fixed-particle run).
Reads config.cds* (position/velocity every 100 steps) and proceced_data/particle_force.csv.
Usage: python3 grafico_particula_libre.py -> ../graficos/particula_libre_eta0.0005.png
"""
import os, glob, re
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
RUNS = [("No electrolyte", "C-particula-libre-hann08-eta0.0005", 0, "#1baf7a"),
        (r"Electrolyte $\lambda_D=6$, started from scratch", "C-particula-libre-electrolito-hann8-LD6-eta0.0005", 0, "#eb6834"),
        (r"Electrolyte $\lambda_D=6$, released from equilibrium", "C-particula-libre-electrolito-hann8-LD6-eta0.0005-equil", 5000, "#2a78d6")]

def read_cfg(f):
    L = [l.strip() for l in open(f)]
    return float(L[1 + 34].split()[2]), float(L[1 + 35].split()[2])

def load(bat, t0):
    d = os.path.join(ROOT, bat, "z_09.50")
    rows = []
    for f in glob.glob(os.path.join(d, "config.cds0*")):
        st = int(re.search(r"cds(\d+)", f).group(1)); z, v = read_cfg(f); rows.append((st - t0, z, v))
    rows.sort(); t, z, v = (np.array(x) for x in zip(*rows))
    fz = {}
    for l in open(os.path.join(d, "proceced_data", "particle_force.csv")):
        if l.startswith("#"): continue
        p = l.strip().split(";"); fz[int(p[0]) - t0] = float(p[5].replace(",", "."))
    tf = np.array(sorted(fz)); ff = np.array([fz[k] for k in tf])
    return t, z - 1.5, v, tf, ff            # d = z - 1.5 (distance to the lower wall plane)

if __name__ == "__main__":
    plt.rcParams.update({"font.size": 10.5, "axes.grid": True, "grid.color": "#e4e2dc",
                         "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False})
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
    print("run | d(0) | d(final) | shift | v_z(final) | E_z(start) | E_z(final)   (E_z reduced = F*1e5/q)")
    for name, bat, t0, c in RUNS:
        t, d, v, tf, ff = load(bat, t0)
        m = t >= 0
        ax[0].plot(t[m], d[m], color=c, lw=1.8, label=name)
        ax[1].plot(tf[tf >= 0], ff[tf >= 0] * 1e5, color=c, lw=1.8)
        ax[2].plot(d[t > 0], -v[t > 0], color=c, lw=1.8)
        j0 = np.argmin(np.abs(tf - 100)); print("%-52s | %.4f | %.4f | %+.4f | %.3e | %.3e | %.3e" % (
            name, d[m][0], d[m][-1], d[m][-1] - d[m][0], v[m][-1], ff[j0] * 1e5, ff[-1] * 1e5))
    ax[0].set(xlabel="steps since release", ylabel="distance to the lower plate, d (lattice units)", title="(a) Position of the particle")
    ax[1].set(xlabel="steps since release", ylabel=r"reduced field on the particle $E_z=F_z\beta/q$", title="(b) Force on the particle")
    ax[2].set(xlabel="d (lattice units)", ylabel=r"approach speed $-v_z$ (lattice units per step)", title="(c) Speed vs distance to the wall")
    ax[0].legend(fontsize=9)
    fig.suptitle("Free charged particle between grounded plates, viscosity 0.0005, initial d = 8 (no wall repulsion in the model)", y=1.02)
    fig.tight_layout(); fig.savefig(os.path.join(ROOT, "graficos", "particula_libre_eta0.0005.png"), dpi=160, bbox_inches="tight")
