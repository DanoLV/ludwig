#!/usr/bin/env python3
"""
Prueba C: una particula subgrid cargada (q = +1, Hann-8) entre dos paredes a
tierra, sin sal, en funcion de su distancia a la pared inferior.

Geometria: caja 32 x 32 x 34, paredes MAP_BOUNDARY en z = 1 y z = 34, planos
medios en z = 1.5 y 33.5 -> L = 32; d = z0 - 1.5. Periodica en x, y (A = 32^2).

Teoria (unidades de Ludwig: psi reducido, 4 pi l_B = beta/epsilon = 10):

 * Carga inducida (Shockley-Ramo): una carga Q a distancia d de la placa
   inferior induce -Q (L-d)/L en la inferior y -Q d/L en la superior. Se suma
   la particula y el fondo de contraiones que agrega psi_electroneutral
   (uniforme en los nodos fluidos, salvo n_near nodos pegados a la particula,
   que quedan en cero -- se incluye exactamente).

 * Campo de imagen sobre la particula (serie de Fourier en x, y con Dirichlet
   en z; promedio de las derivadas laterales en z0, que elimina el termino
   propio). Se toma de particle_force.csv (E = F beta / q), que incluye la
   correccion pm_sr cuando esta activa:
     E_z = -(4 pi l_B q / A) [ (L-2d)/(2L) + sum_{k != 0} sinh(k(L-2d)) / (2 sinh kL) ]
   El fondo uniforme -q aporta exactamente +(4 pi l_B q/A)(L-2d)/(2L) y cancela
   el termino k = 0. Con electrokinetics_electroneutral no (sin fondo) el
   termino k = 0 queda, y la carga inducida total es -q. Para una carga con simetria esferica que no toca la pared
   la imagen actua como sobre una puntual, asi que vale para el kernel mientras
   su soporte no alcance los nodos de pared (Hann-8: z0 >= 5, d >= 3.5).

Uso: analizar_C.py [directorio de la bateria] [orden Hann n]
     -> tabla + graficos/C_carga_imagen_hannNN.png
"""
import glob
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
BAT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "C-carga-imagen-hann08")
NHANN = int(sys.argv[2]) if len(sys.argv) > 2 else 8

NXY, NZ = 32, 34
L = 32.0
A = float(NXY * NXY)
FOURPI_LB = 10.0          # beta / epsilon = 1e5 / 1e4
BETA = 1.0e5
Q = 1.0
# Hann-n: peso nulo para |x| >= n/2 -> no toca el nodo de pared z = 1 si
# z0 - 1 >= n/2, es decir d = z0 - 1.5 >= n/2 - 1/2
D_MIN = NHANN / 2.0 - 0.5


def image_field(d, M=80, k0=False):
    """E_z de imagen; k0=True agrega el termino k = 0 (corrida sin fondo)."""
    m = np.arange(-M, M + 1)
    mx, my = np.meshgrid(m, m, indexing="ij")
    k = 2 * np.pi * np.hypot(mx, my) / NXY
    k = k[k > 0]
    a, b = k * (L - 2 * d), k * L
    # sinh(a)/sinh(b) estable: e^{a-b} (1 - e^{-2a}) / (1 - e^{-2b})
    ratio = np.exp(a - b) * (-np.expm1(-2 * a)) / (-np.expm1(-2 * b))
    s = np.sum(ratio / 2.0) + ((L - 2 * d) / (2 * L) if k0 else 0.0)
    return -(FOURPI_LB * Q / A) * s


def read_run(d):
    log = open(os.path.join(d, "run.log"), errors="ignore").read()
    wq = re.findall(r"Wall charge lower\s+(\S+) upper\s+(\S+)", log)
    fq = re.findall(r"Fluid charge\s+(\S+) total\s+(\S+)", log)
    nn = re.search(r"n_near_global=(\d+) vf=(\d+) vf_eff=(\d+)", log)
    rows = [l for l in open(os.path.join(d, "proceced_data", "particle_force.csv"))
            if not l.startswith("#")]
    f = rows[0].strip().split(";")
    # E reducido = F beta / q  (F incluye la correccion pm_sr si esta activa)
    ez = float(f[5].replace(",", ".")) * BETA / Q
    ex = float(f[3].replace(",", ".")) * BETA / Q
    ey = float(f[4].replace(",", ".")) * BETA / Q
    fondo = nn is not None
    corr = bool(re.search(r"^PM short-range correction tables built", log, re.M))
    return dict(qlo=float(wq[0][0]), qhi=float(wq[0][1]), qfl=float(fq[0][0]),
                ez=ez, exy=max(abs(ex), abs(ey)), fondo=fondo, corr=corr,
                nnear=int(nn.group(1)) if fondo else 0,
                vf=int(nn.group(2)) if fondo else 0,
                vfeff=int(nn.group(3)) if fondo else 1)


def charge_theory(z0, r):
    d = z0 - 1.5
    q_lo = -Q * (L - d) / L
    q_hi = -Q * d / L
    if not r["fondo"]:
        return q_lo, q_hi
    c = -Q / r["vfeff"]                         # carga de fondo por nodo
    # suma sobre todos los nodos fluidos: <(L-d)/L> = 1/2
    lo_bg = -c * (0.5 * r["vf"] - r["nnear"] * (L - d) / L)
    hi_bg = -c * (0.5 * r["vf"] - r["nnear"] * d / L)
    return q_lo + lo_bg, q_hi + hi_bg


def main():
    dirs = sorted(glob.glob(os.path.join(BAT, "z_*")),
                  key=lambda p: float(p.rsplit("_", 1)[1]))
    Z, D, EZ, TH, QL, QH, TL, TH2, QF, EXY = ([] for _ in range(10))
    print("%6s %6s %13s %13s %9s %12s %12s %10s %9s" %
          ("z0", "d", "Ez Ludwig", "Ez teoria", "err rel", "Qinf Ludwig",
           "Qinf teoria", "err Q", "Q perdida"))
    for dd in dirs:
        z0 = float(dd.rsplit("_", 1)[1])
        r = read_run(dd)
        d = z0 - 1.5
        th = image_field(d, k0=not r["fondo"])
        tl, th_hi = charge_theory(z0, r)
        FONDO, CORR = r["fondo"], r["corr"]
        Z.append(z0); D.append(d); EZ.append(r["ez"]); TH.append(th)
        QL.append(r["qlo"]); QH.append(r["qhi"]); TL.append(tl); TH2.append(th_hi)
        # carga neta esperada en el fluido: 0 con fondo, q sin fondo; lo que
        # falta respecto de eso es carga depositada en nodos de pared (perdida)
        QF.append(r["qfl"] - (0.0 if r["fondo"] else Q)); EXY.append(r["exy"])
        rel = (r["ez"] - th) / abs(th) if abs(th) > 1e-12 else float("nan")
        print("%6.2f %6.2f %13.6e %13.6e %+9.2e %12.8f %12.8f %+10.2e %+9.2e%s" %
              (z0, d, r["ez"], th, rel, r["qlo"], tl, r["qlo"] - tl, QF[-1],
               "   <- kernel toca la pared" if d < D_MIN else ""))
    D, EZ, TH = map(np.array, (D, EZ, TH))
    QL, QH, TL, TH2, QF = map(np.array, (QL, QH, TL, TH2, QF))
    ok = D >= D_MIN

    # ---- grafico -------------------------------------------------------
    C1, C2, C3, REF = "#2a78d6", "#eb6834", "#1baf7a", "#55544f"
    plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.color": "#e4e2dc",
                         "axes.spines.top": False, "axes.spines.right": False,
                         "legend.frameon": False, "lines.linewidth": 2.0,
                         "savefig.dpi": 160, "savefig.bbox": "tight"})
    fig, ax = plt.subplots(2, 2, figsize=(11.5, 7.8))

    def overlap(a):
        a.axvspan(D.min() - 0.3, D_MIN, color="#efe9dd", lw=0, zorder=0)

    dd = np.linspace(0.8, 16, 400)
    a = ax[0, 0]; overlap(a)
    a.plot(dd, [image_field(x, k0=not FONDO) for x in dd], color=REF, lw=1.3, ls="--",
           label="theory: images, periodic in x,y")
    a.plot(dd, -FOURPI_LB * Q / (4 * np.pi) / (2 * dd) ** 2, color=C3, lw=1.2, ls=":",
           label="single wall, isolated charge")
    a.plot(D[ok], EZ[ok], "o", ms=5, color=C1, label="Ludwig, PETSc (Hann-%d)" % NHANN)
    a.plot(D[~ok], EZ[~ok], "s", ms=5, color=C2, mfc="white", mew=1.5,
           label="Ludwig, kernel reaches the wall")
    a.set(xlabel="distance to the wall plane d (lattice units)",
          ylabel=r"$E_z$ on the particle (reduced units)", ylim=(-0.04, 0.003),
          title="(a) Image field on the particle")
    a.text(0.9, -0.037, "kernel\noverlaps\nwall", fontsize=8.5, color="#6b6960")
    a.legend(loc="lower right")

    a = ax[0, 1]; overlap(a)
    m = D < 15.9                     # en el centro la fuerza es 0
    rel = np.full_like(EZ, np.nan)
    rel[m] = np.abs(EZ[m] - TH[m]) / np.abs(TH[m])                     # en el centro la fuerza es 0
    a.semilogy(D[m & ok], rel[m & ok], "o", ms=5, color=C1, label="kernel clear of the wall")
    a.semilogy(D[m & ~ok], rel[m & ~ok], "s", ms=5, color=C2, mfc="white", mew=1.5,
               label="kernel reaches the wall")
    a.set(xlabel="distance to the wall plane d (lattice units)",
          ylabel=r"$|E_z - E_{theory}| / |E_{theory}|$", title="(b) Relative error of the image field")
    a.legend(loc="upper right")

    a = ax[1, 0]; overlap(a)
    bg = 0.5 if FONDO else 0.0
    a.plot(dd, [-(L - x) / L + bg for x in dd], color=REF, lw=1.3, ls="--",
           label="theory (Shockley-Ramo)")
    a.plot(dd, [-x / L + bg for x in dd], color=REF, lw=1.3, ls="--")
    a.plot(D, QL, "o", ms=5, color=C1, label="lower wall, Ludwig")
    a.plot(D, QH, "^", ms=5, color=C2, label="upper wall, Ludwig")
    a.set(xlabel="distance to the wall plane d (lattice units)",
          ylabel="induced charge (units of e)", title="(c) Charge induced on each wall")
    a.legend(loc="center left", bbox_to_anchor=(0.14, 0.44))

    a = ax[1, 1]; overlap(a)
    a.semilogy(D, np.maximum(np.abs(QL - TL), 1e-17), "o", ms=5, color=C1,
               label="lower wall vs theory")
    a.semilogy(D, np.maximum(np.abs(QH - TH2), 1e-17), "^", ms=5, color=C2,
               label="upper wall vs theory")
    a.semilogy(D, np.maximum(np.abs(QF), 1e-17), "x", ms=6, color=REF, mew=1.5,
               label="charge lost into the wall")
    a.set(xlabel="distance to the wall plane d (lattice units)", ylabel="absolute error (units of e)",
          ylim=(1e-17, 1), title="(d) Charge error and charge lost")
    a.legend(loc="upper right")

    fig.suptitle("Test C - one charged subgrid particle (q = +1, Hann-%d) between two grounded walls "
                 "(L = 32, box 32 x 32, no salt)\ncounterion background: %s   |   pm_sr correction: %s"
                 % (NHANN, "yes (neutral)" if FONDO else "no (net charge +1)",
                    "yes (measured reference)" if CORR else "no"), y=1.03, fontsize=11)
    fig.tight_layout()
    out = os.path.join(ROOT, "graficos",
                       "C_" + os.path.basename(os.path.normpath(BAT))[2:].replace("-", "_") + ".png")
    fig.savefig(out); print("  ->", os.path.relpath(out, ROOT))

    print("\nresumen (kernel sin tocar la pared, d >= %.1f, sin el centro):" % D_MIN)
    mm = ok & m
    print("  campo: error relativo max %.2e, mediana %.2e" % (rel[mm].max(), np.median(rel[mm])))
    print("  carga: error absoluto max %.2e" % max(np.abs(QL - TL)[ok].max(), np.abs(QH - TH2)[ok].max()))
    print("  carga perdida max %.2e" % np.abs(QF[ok]).max())


if __name__ == "__main__":
    main()
