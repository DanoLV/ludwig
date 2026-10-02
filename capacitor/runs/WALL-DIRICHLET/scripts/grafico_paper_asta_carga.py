#!/usr/bin/env python3
"""
Comparacion DIRECTA con la Fig. 4a de Asta et al. 2019 (dinamica de carga del
capacitor de placas paralelas, seccion III.A.3): Q(t)/Q_inf vs t/tau_norm,
con tau_norm = L*lambda_D/(2D), y el ajuste exponencial
Q(t) = Q_inf + (Q0 - Q_inf) * exp(-t/tau).

Q(t) se saca de las lineas "Wall charge lower/upper" que Ludwig imprime cada
freq_psi_resid pasos (ver E-paper-asta-carga/input). Se usa la carga de la
pared superior (mismo signo que Dpsi>0); Q0 es el primer valor (capacitor
neutro, antes de que los iones se muevan) y Q_inf el ultimo (estacionario).

Lee <carpeta>/run.log. Escribe graficos/<nombre>.png.
Uso: python3 grafico_paper_asta_carga.py [carpeta] [nombre_salida] [titulo_extra]
  (por defecto: E-paper-asta-carga / E_paper_asta_carga / "")
"""
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "graficos")
DIRNAME = sys.argv[1] if len(sys.argv) > 1 else "E-paper-asta-carga"
OUTNAME = sys.argv[2] if len(sys.argv) > 2 else "E_paper_asta_carga"
SUBTITLE = sys.argv[3] if len(sys.argv) > 3 else ""
E = os.path.join(ROOT, DIRNAME)

C1 = "#2a78d6"
REF = "#55544f"
GRID = "#e4e2dc"
plt.rcParams.update({
    "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
    "axes.edgecolor": "#8a887f", "axes.linewidth": 0.8,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False,
    "legend.frameon": True, "legend.facecolor": "white", "legend.framealpha": 1.0,
    "legend.edgecolor": "#8a887f", "legend.fancybox": False,
    "lines.linewidth": 2.0,
    "savefig.dpi": 160, "savefig.bbox": "tight",
})

# --- Parametros, identicos a D-paper-asta/ (mismo capacitor) ----------------
L = 76.0
LAMBDA_D = 6.0
D0 = 0.05
FREQ = 50          # freq_psi_resid del input: paso entre lineas sucesivas
TAU_NORM = L * LAMBDA_D / (2 * D0)


def read_charge_series(path):
    pat = re.compile(r"Wall charge lower\s+([\-0-9.eE+]+)\s+upper\s+([\-0-9.eE+]+)")
    lo, hi = [], []
    with open(path) as f:
        for line in f:
            m = pat.search(line)
            if m:
                lo.append(float(m.group(1)))
                hi.append(float(m.group(2)))
    return np.array(lo), np.array(hi)


def main():
    lo, hi = read_charge_series(os.path.join(E, "run.log"))
    n = len(hi)
    # El primer "Wall charge" impreso en el log NO es t=0: ntimestep empieza
    # en 1 (no en 0) para esta cuenta, asi que la condicion
    # "ntimestep % freq_psi_resid == 0" recien se cumple en el paso FREQ, no
    # en el paso 0 (verificado contando las lineas "Current step =0" antes
    # del primer "Wall charge": salen 49, con FREQ=50). O sea la muestra i
    # (i=0,1,2,...) corresponde al paso (i+1)*FREQ, no a i*FREQ.
    t = (np.arange(n) + 1) * FREQ

    Q = hi                            # carga de la pared superior (signo de Dpsi>0)
    Q0, Qinf = Q[0], Q[-1]
    print("Dinamica de carga (Asta et al. 2019, Fig. 4a):")
    print("  tau_norm = L*lambda_D/(2D) = %.1f pasos" % TAU_NORM)
    print("  Q0  (capacitor neutro, t=0)   = %.6e" % Q0)
    print("  Qinf (estacionario, t=%d)     = %.6e" % (t[-1], Qinf))
    print("  N muestras = %d, hasta t/tau_norm = %.2f" % (n, t[-1] / TAU_NORM))

    # Ajuste exponencial: ln[(Q-Qinf)/(Q0-Qinf)] = -t/tau (lineal, sin el primer
    # punto t=0 que es exacto por construccion, y recortando cuando Q-Qinf
    # se vuelve puro ruido de iteracion del solver).
    resid = (Q - Qinf) / (Q0 - Qinf)
    good = resid > 1e-3                      # antes de tocar el piso de ruido
    slope, intercept = np.polyfit(t[good], np.log(resid[good]), 1)
    tau_fit = -1.0 / slope
    print("  tau (ajuste exponencial)      = %.1f pasos" % tau_fit)
    print("  tau / tau_norm                = %.4f" % (tau_fit / TAU_NORM))

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))

    a = ax[0]
    tn = t / TAU_NORM
    tt = np.linspace(0, tn.max(), 300)
    a.plot(tn, Q / Qinf, "o", ms=4, color=C1, label="Ludwig, PETSc (this work)")
    a.plot(tt, (Qinf + (Q0 - Qinf) * np.exp(-tt * TAU_NORM / tau_fit)) / Qinf,
           color=REF, lw=1.6, ls="--",
           label=r"$Q_\infty + (Q_0-Q_\infty)e^{-t/\tau}$, $\tau=%.0f$" % tau_fit)
    a.set(xlabel=r"$t / (L\lambda_D/2D)$", ylabel=r"$Q(t)/Q_\infty$",
          title="(a) Charging a parallel plate capacitor")
    a.legend(loc="lower right")

    a = ax[1]
    a.semilogy(tn[good], resid[good], "o", ms=4, color=C1, label="Ludwig, PETSc")
    a.semilogy(tt, np.exp(-tt * TAU_NORM / tau_fit), color=REF, lw=1.6, ls="--",
               label=r"$e^{-t/\tau}$")
    a.set(xlabel=r"$t / (L\lambda_D/2D)$", ylabel=r"$(Q(t)-Q_\infty)/(Q_0-Q_\infty)$",
          title=r"(b) Exponential relaxation, $\tau/\tau_{norm}=%.3f$" % (tau_fit / TAU_NORM))
    a.legend(loc="upper right")

    extra = (" " + SUBTITLE) if SUBTITLE else ""
    fig.suptitle("Direct comparison with Asta et al. 2019, Fig. 4a "
                  r"(charging dynamics, $L=76\,\Delta x$, $\lambda_D=6\,\Delta x$)" + extra,
                  y=1.03, fontsize=11)
    fig.tight_layout()
    os.makedirs(OUT, exist_ok=True)
    out = os.path.join(OUT, OUTNAME + ".png")
    fig.savefig(out)
    plt.close(fig)
    print("  ->", os.path.join("graficos", OUTNAME + ".png"))


if __name__ == "__main__":
    main()
