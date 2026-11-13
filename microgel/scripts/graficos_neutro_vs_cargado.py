#!/usr/bin/env python3
"""
Compara el microgel neutro en equilibrio con el mismo microgel cargado en la
superficie, sin sal, tambien en equilibrio (misma estructura de red de
partida: MGEL-RED-neutro/config.cds00060000 == MGEL-RED-SURF-NOSALT
config.cds.init, antes de aplicar la carga).

(a) Trayectoria de Rg: relajacion neutra (0->60000) y relajacion bajo carga
    (0->500000), cada una con su propia asintota extrapolada.
(b) Radio de cada monomero al centro de masa, antes (=neutro) y despues
    (=cargado en equilibrio), separado en capa cargada (65 monomeros
    externos) y nucleo sin carga (93 monomeros), para ver DONDE hincha.

Escribe graficos/neutro_vs_cargado.png
"""
import json
import math
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def read_config(path):
    lines = open(path).readlines()
    n = int(lines[0]); rest = lines[1:]; b = len(rest) // n
    out = {}
    for i in range(n):
        blk = rest[i * b:(i + 1) * b]
        out[int(blk[0])] = ([float(x) for x in blk[34].split()], float(blk[45]))
    return out


# --- paleta de referencia (orden fijo) ---------------------------------
C_CORE, C_SHELL, REF = "#2a78d6", "#eb6834", "#55544f"
GRID = "#e4e2dc"
plt.rcParams.update({
    "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
    "axes.edgecolor": "#8a887f", "axes.linewidth": 0.8,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False,
    "legend.frameon": False, "lines.linewidth": 2.0,
    "savefig.dpi": 160, "savefig.bbox": "tight",
})

fig, ax = plt.subplots(1, 2, figsize=(12, 4.6), gridspec_kw={"width_ratios": [1.15, 1]})

# ---- (a) trayectorias --------------------------------------------------
neu = json.load(open(os.path.join(ROOT, "MGEL-RED-neutro", "serie.json")))
car = json.load(open(os.path.join(ROOT, "MGEL-RED-SURF-NOSALT", "serie.json")))
neu_steps = sorted(int(k) for k in neu)
car_steps = sorted(int(k) for k in car)
neu_rg = [neu[str(s)]["Rg"] for s in neu_steps]
car_rg = [car[str(s)]["Rg"] for s in car_steps]
rg_inf_neu, rg_inf_car = 6.581, 7.303   # ajuste exponencial (mgel_equilibrio.py)

a = ax[0]
a.semilogx(neu_steps, neu_rg, color=C_CORE, label="neutral relaxation")
a.semilogx(car_steps, car_rg, color=C_SHELL, label="charged, no salt (from the neutral structure)")
a.axhline(rg_inf_neu, color=C_CORE, lw=1.0, ls=":", zorder=0)
a.axhline(rg_inf_car, color=C_SHELL, lw=1.0, ls=":", zorder=0)
a.annotate("neutral: %.3f" % rg_inf_neu, (1.3, rg_inf_neu), xytext=(1.3, rg_inf_neu - 0.09),
           color=C_CORE, fontsize=9)
a.annotate("charged: %.3f (+%.1f%%)" % (rg_inf_car, 100 * (rg_inf_car / rg_inf_neu - 1)),
           (1.3, rg_inf_car), xytext=(1.3, rg_inf_car + 0.04), color=C_SHELL, fontsize=9)
a.set(xlabel="step (each run's own clock, from its starting structure)",
      ylabel=r"$R_g$ (lattice units)", xlim=(1, 5.5e5), ylim=(6.4, 7.4),
      title="(a) Radius of gyration: neutral vs. charged relaxation")
a.legend(loc="lower right")

# ---- (b) radios por grupo, antes/despues -------------------------------
before = read_config(os.path.join(ROOT, "MGEL-RED-SURF-NOSALT", "config.cds00000001.001-001"))
after = read_config(os.path.join(ROOT, "MGEL-RED-SURF-NOSALT", "config.cds00500000.001-001"))
ids = sorted(before)
n = len(ids)
com_b = [sum(before[i][0][k] for i in ids) / n for k in range(3)]
com_a = [sum(after[i][0][k] for i in ids) / n for k in range(3)]
shell_ids = [i for i in ids if before[i][1] != 0]
core_ids = [i for i in ids if before[i][1] == 0]


def radii(cfg, com, sel):
    return np.array([math.dist(cfg[i][0], com) for i in sel])


groups = [
    ("shell, charged\n(65 outer monomers)", radii(before, com_b, shell_ids), radii(after, com_a, shell_ids), C_SHELL),
    ("core, uncharged\n(93 inner monomers)", radii(before, com_b, core_ids), radii(after, com_a, core_ids), C_CORE),
]

a = ax[1]
rng = np.random.default_rng(0)
for row, (label, rb, ra, col) in enumerate(groups):
    y0, y1 = row * 1.0 + 0.72, row * 1.0 + 0.28
    jb = y0 + rng.uniform(-0.07, 0.07, size=rb.size)
    ja = y1 + rng.uniform(-0.07, 0.07, size=ra.size)
    a.plot(rb, jb, "o", ms=5, mfc="none", mec=col, mew=1.3, alpha=0.75, zorder=2)
    a.plot(ra, ja, "o", ms=5, color=col, alpha=0.85, zorder=3)
    a.plot([rb.mean()] * 2, [y0 - 0.1, y0 + 0.1], color=col, lw=2.5, zorder=4)
    a.plot([ra.mean()] * 2, [y1 - 0.1, y1 + 0.1], color=col, lw=2.5, zorder=4)
    a.annotate("%.2f" % rb.mean(), (rb.mean(), y0 + 0.14), ha="center", fontsize=8, color=col)
    a.annotate("%.2f" % ra.mean(), (ra.mean(), y1 - 0.20), ha="center", fontsize=8, color=col)

a.set_yticks([0.72, 0.28, 1.72, 1.28])
a.set_yticklabels(["before\n(neutral)", "after\n(charged)", "before\n(neutral)", "after\n(charged)"], fontsize=8.5)
for row, (label, *_ ) in enumerate(groups):
    a.text(-0.5, row * 1.0 + 0.5, label, fontsize=9.5, color="#3b3a36", va="center", ha="left",
           fontweight="medium")
a.set(xlabel="distance to center of mass (lattice units)", ylim=(-0.15, 2.15), xlim=(0, 13.5),
      title="(b) Where the swelling happens (per-monomer radius)")
a.set_xlim(left=-4.3)
a.axvline(0, color="#c9c6bb", lw=0.7)

fig.suptitle("Neutral microgel vs. surface-charged microgel, no added salt (same crosslinked network, z = 4.72)",
             y=1.03, fontsize=11.5)
fig.tight_layout()
out = os.path.join(ROOT, "graficos", "neutro_vs_cargado.png")
os.makedirs(os.path.dirname(out), exist_ok=True)
fig.savefig(out)
print("->", os.path.relpath(out, ROOT))
