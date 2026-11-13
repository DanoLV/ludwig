#!/usr/bin/env python3
"""
Comparacion DIRECTA con la Fig. 3 de Asta et al. 2019 (capacitor de placas
paralelas, seccion III.A.1): mismos parametros reducidos del paper
(L = 76 Dx, lambda_D = 6 Dx, beta*e*Dpsi = 0.1, D = 0.05), corridos con nuestra
implementacion de paredes a potencial constante. Mismo formato de la figura
del paper: panel (a) perfil de potencial, panel (b) concentracion de iones.

A diferencia de graficos_paredes.py (pruebas A y B), aca NO se calcula la
solucion discreta de referencia ni un panel de desviacion: se compara solo
contra la formula analitica continua (ec. 17 del paper) y la distribucion de
Boltzmann que de ella se deriva, que es lo unico que el paper mismo grafica
en su Fig. 3.

Lee la salida de D-paper-asta/. Escribe graficos/D_paper_asta.png.
Uso: python3 grafico_paper_asta.py
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "graficos")
D = os.path.join(ROOT, "D-paper-asta")

C1, C2 = "#2a78d6", "#eb6834"
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

# --- Parametros, identicos a los del input (ver comentarios alli) ----------
NX, NY, NZ = 4, 4, 78
STEP = 10000
BETA, EUNIT, EPS = 1.0e5, 1.0, 1.0e4
RHOS = 1.388889e-3
D0 = 0.05
L = 76.0                       # planos medios en z=1.5 y z=77.5 (NZ=78)
ZC = 39.5                      # centro del capacitor
# psi1, psi2 del paper (0.1, 0.2 kT/e); Ludwig corrio con 0_0.1 (misma Dpsi,
# ver nota en el input) y se le suma 0.1 solo para graficar en los mismos
# ejes que la Fig. 3a -- un corrimiento aditivo no cambia la fisica (el paper
# aclara esto explicitamente).
PSI1, PSI2 = 0.1, 0.2
OFFSET = PSI1

KAPPA = np.sqrt(8 * np.pi * (BETA * EUNIT**2 / (4 * np.pi * EPS)) * RHOS)
LAMBDA_D = 1.0 / KAPPA


def load(path, ncol=1):
    a = np.loadtxt(path)
    return a.reshape(NX, NY, NZ, ncol) if ncol > 1 else a.reshape(NX, NY, NZ)


def dh_continuum(z):
    """Ec. 17 del paper, con z medido desde el centro del capacitor."""
    return 0.5 * (PSI1 + PSI2) + 0.5 * (PSI2 - PSI1) * np.sinh(KAPPA * z) / np.sinh(KAPPA * L / 2)


def main():
    psi = load(os.path.join(D, "psi-%09d.001-001" % STEP))[0, 0, :] + OFFSET
    q = load(os.path.join(D, "qsi-%09d.001-001" % STEP), 2)[0, 0, :, :]
    Z = np.arange(1, NZ + 1)
    FL = slice(1, NZ - 1)           # nodos fluidos z = 2..77
    zL = (Z - ZC) / L               # eje horizontal del paper: z/L

    zz = np.linspace(-L / 2, L / 2, 400)
    dh = dh_continuum(zz)
    pb = 0.5 * (PSI1 + PSI2)
    # Boltzmann COMPLETO (no linealizado) con el psi de DH: con beta*e*Dpsi=0.1
    # la no linealidad ya no es despreciable como en B-debye-huckel/ (Dpsi=1e-6).
    boltz_cat = np.exp(-(dh - pb))
    boltz_anu = np.exp(+(dh - pb))

    err = np.abs(psi[FL] - dh_continuum(Z[FL] - ZC)) / (PSI2 - PSI1)
    print("Replica del capacitor de Asta et al. 2019 (Fig. 3):")
    print("  L = %.1f Dx (paper: 76)" % L)
    print("  lambda_D = %.4f Dx (paper: 6)" % LAMBDA_D)
    print("  L/lambda_D = %.3f (paper: 76/6 = %.3f)" % (L / LAMBDA_D, 76 / 6))
    print("  beta*e*Dpsi = %.3f (paper: 0.1)" % (PSI2 - PSI1))
    print("  max |psi - psi_DH| / Dpsi = %.3e" % err.max())

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))

    a = ax[0]
    a.plot(zz / L, dh, color=REF, lw=1.6, ls="--", label="Debye-Hückel (eq. 17, Asta et al. 2019)")
    a.plot(zL[FL], psi[FL], "o", ms=4.5, color=C1, label="Ludwig, PETSc (this work)")
    a.set(xlabel=r"$z/L$", ylabel=r"$\beta e\,\psi(z)$", xlim=(-0.5, 0.5),
          title=r"(a) Potential profile, $L/\lambda_D = %.2f$, $\beta e\Delta\psi = 0.1$" % (L / LAMBDA_D))
    a.legend(loc="upper left")

    a = ax[1]
    a.plot(zz / L, boltz_cat, color=REF, lw=1.6, ls="--", label="Debye-Hückel (Boltzmann)")
    a.plot(zz / L, boltz_anu, color=REF, lw=1.6, ls="--")
    a.plot(zL[FL], q[FL, 0] / RHOS, "o", ms=4.5, color=C1, label=r"Ludwig, $\rho_+/\rho_s$")
    a.plot(zL[FL], q[FL, 1] / RHOS, "^", ms=4.5, color=C2, label=r"Ludwig, $\rho_-/\rho_s$")
    a.set(xlabel=r"$z/L$", ylabel=r"$\rho_\pm(z)/\rho_s$", xlim=(-0.5, 0.5),
          title="(b) Ionic concentration profiles")
    a.legend(loc="upper center")

    fig.suptitle("Direct comparison with Asta et al. 2019, Fig. 3 "
                  r"(parallel plate capacitor, $L=76\,\Delta x$, $\lambda_D=6\,\Delta x$)",
                  y=1.03, fontsize=11)
    fig.tight_layout()
    os.makedirs(OUT, exist_ok=True)
    out = os.path.join(OUT, "D_paper_asta.png")
    fig.savefig(out)
    plt.close(fig)
    print("  ->", os.path.join("graficos", "D_paper_asta.png"))


if __name__ == "__main__":
    main()
