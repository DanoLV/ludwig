#!/usr/bin/env python3
"""
Version de grafico_paper_asta_carga.py con SOLO el panel (a) (Q(t)/Q_inf vs
t/tau_norm) y SIN el ajuste exponencial superpuesto -- solo los datos crudos
de Ludwig. Reutiliza la lectura de datos de grafico_paper_asta_carga.py.

Uso: python3 grafico_paper_asta_carga_sinfit.py [carpeta] [nombre_salida] [titulo_extra]
  (por defecto: E-paper-asta-carga-A0 / E_paper_asta_carga_A0_sinfit / "")
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import grafico_paper_asta_carga as base   # noqa: E402

ROOT = base.ROOT
OUT = base.OUT
DIRNAME = sys.argv[1] if len(sys.argv) > 1 else "E-paper-asta-carga-A0"
OUTNAME = sys.argv[2] if len(sys.argv) > 2 else "E_paper_asta_carga_A0_sinfit"
SUBTITLE = sys.argv[3] if len(sys.argv) > 3 else ""
E = os.path.join(ROOT, DIRNAME)

C1 = base.C1
TAU_NORM = base.TAU_NORM
FREQ = base.FREQ

plt.rcParams.update({
    "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
    "axes.edgecolor": "#8a887f", "axes.linewidth": 0.8,
    "axes.grid": True, "grid.color": "#e4e2dc", "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False,
    "legend.frameon": True, "legend.facecolor": "white", "legend.framealpha": 1.0,
    "legend.edgecolor": "#8a887f", "legend.fancybox": False,
    "lines.linewidth": 2.0,
    "savefig.dpi": 160, "savefig.bbox": "tight",
})


def main():
    lo, hi = base.read_charge_series(os.path.join(E, "run.log"))
    n = len(hi)
    t = (np.arange(n) + 1) * FREQ
    Q = hi
    Qinf = Q[-1]
    tn = t / TAU_NORM

    fig, a = plt.subplots(1, 1, figsize=(6.5, 4.6))
    a.plot(tn, Q / Qinf, "o", ms=1.8, mew=0, color=C1, label="Ludwig, PETSc (this work)")
    a.set(xlabel=r"$t / (L\lambda_D/2D)$", ylabel=r"$Q(t)/Q_\infty$",
          title="Charging a parallel plate capacitor")
    a.legend(loc="lower right")

    extra = (" " + SUBTITLE) if SUBTITLE else ""
    fig.suptitle("Asta et al. 2019, Fig. 4a "
                  r"(charging dynamics, $L=76\,\Delta x$, $\lambda_D=6\,\Delta x$)" + extra,
                  y=1.0, fontsize=11)
    fig.tight_layout()
    os.makedirs(OUT, exist_ok=True)
    out = os.path.join(OUT, OUTNAME + ".png")
    fig.savefig(out)
    plt.close(fig)
    print("  ->", os.path.join("graficos", OUTNAME + ".png"))


if __name__ == "__main__":
    main()
