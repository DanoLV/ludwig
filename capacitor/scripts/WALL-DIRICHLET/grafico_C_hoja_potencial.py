#!/usr/bin/env python3
"""Potencial electrostatico de Ludwig (carga q = +1 subgrid, Hann-8, entre placas a tierra, electrolito 1:1)
contra la solucion analitica de la hoja (separacion de variables, Fourier en z + K0 en rho), CORREGIDA:

  phi_q(rho, z) = (q/(pi L eps)) sum_{m>=1} sen[m pi (d+L/2)/L] sen[m pi (z+L/2)/L] K0( sqrt((m pi/L)^2 + kappa^2) rho )
  phi_TOT(abierta)  = phi_1 sinh k(L/2-z)/sinh kL + phi_2 sinh k(L/2+z)/sinh kL + phi_q     (placas en +-L/2)
  phi_TOT(cerrada)  = phi_TOT(abierta) + psibar * [1 - cosh(kappa z)/cosh(kappa L/2)],
        psibar = (kappa L/2) coth(kappa L/2) <phi_q>,  <phi_q> = (1/A L) int phi_q dV
               = (q/(pi L eps)) (2 pi/(A L)) sum_{m impar} sen(m pi (d+L/2)/L) (2L/(m pi)) / g_m^2,  g_m^2 = (m pi/L)^2 + kappa^2
La "abierta" supone sal en contacto con un reservorio (n0 a potencial cero); Ludwig conserva el numero de iones en la caja
de area A = 32^2, y eso agrega el desplazamiento psibar (orden 1/A) y h(z). Mas replicas laterales (|i|,|j| <= 2).
Unidades reducidas: q/eps = 4 pi l_B q = 10.   Uso: python3 grafico_C_hoja_potencial.py
"""
import os, glob
import numpy as np
from scipy.special import k0 as _k0                   # modified Bessel K0 (scipy.special, Cephes)
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
L = 32.0; ELL = 32.0; A = ELL ** 2; FOUR = 10.0; M = 400
m_ = np.arange(1, M + 1)

def phi_q(rho, z, d, kap):
    g = np.sqrt((m_ * np.pi / L) ** 2 + kap ** 2)
    k0 = _k0(g * rho)
    sd = np.sin(m_ * np.pi * (d + L / 2) / L)
    sz = np.sin(np.outer(np.atleast_1d(z) + L / 2, m_) * np.pi / L)
    return FOUR / (np.pi * L) * (sz * (sd * k0)).sum(axis=1)

def phi_q_box(x, y, z, d, kap, NR=2):
    tot = 0.0
    for i in range(-NR, NR + 1):
        for j in range(-NR, NR + 1):
            tot = tot + phi_q(max(np.hypot(x - i * ELL, y - j * ELL), 1e-6), z, d, kap)
    return tot

def h_closed(z, d, kap):
    mo = np.arange(1, 4001, 2); g2 = (mo * np.pi / L) ** 2 + kap ** 2
    mean = FOUR / (np.pi * L) * (2 * np.pi / (A * L)) * np.sum(np.sin(mo * np.pi * (d + L / 2) / L) * (2 * L / (mo * np.pi)) / g2)
    pbar = (kap * L / 2) / np.tanh(kap * L / 2) * mean
    return pbar * (1 - np.cosh(kap * z) / np.cosh(kap * L / 2)), pbar

def ludwig(bat, dl, dx):
    p = sorted(glob.glob(os.path.join(bat, "z_%05.2f" % (dl + 1.5), "psi-0*")))[-1]
    a = np.loadtxt(p).reshape(32, 32, 34)
    zc = np.arange(1, 35) - 1.5
    return zc[1:33] - L / 2, a[(15 + dx) % 32, 15, 1:33]

if __name__ == "__main__":
    plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.color": "#e4e2dc",
                         "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False})
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6), sharey=True)
    cfg = [(8, 6, 8.0, "Hann-8, $\\lambda_D=6$, $d=8$"), (8, 4, 8.0, "Hann-8, $\\lambda_D=4$, $d=8$"),
           (8, 8, 8.0, "Hann-8, $\\lambda_D=8$, $d=8$")]
    cols = {6: "#2a78d6", 10: "#eb6834"}
    print("hann LD  rho   max|Ludwig-teoria(abierta)|   max|Ludwig-teoria(cerrada)|   (max ψ en la línea)   psibar")
    for a, (hann, ld, dl, tit) in zip(ax, cfg):
        bat = os.path.join(ROOT, "C-electrolito-hann%d-LD%d" % (hann, ld)); kap = 1.0 / ld; d = dl - L / 2
        zz = np.linspace(-L / 2 + 0.3, L / 2 - 0.3, 60)
        for dx, c in cols.items():
            z, ps = ludwig(bat, dl, dx)
            op = phi_q_box(dx, 0.0, z, d, kap); sh, pbar = h_closed(z, d, kap)
            opf = phi_q_box(dx, 0.0, zz, d, kap); shf, _ = h_closed(zz, d, kap)
            a.plot(zz, opf, "--", color=c, lw=1.2)
            a.plot(zz, opf + shf, "-", color=c, lw=1.6)
            a.plot(z, ps, "o", ms=3.5, color=c, label="Ludwig, $\\rho=%d$" % dx)
            print("%4d %2d %4d   %.4f      %.4f      (%.4f)   %.4f" % (hann, ld, dx, np.abs(ps - op).max(),
                  np.abs(ps - op - sh).max(), ps.max(), pbar))
        a.axvline(d, color="#999", lw=0.8, ls=":")
        a.set(xlabel="z (centro entre placas = 0)", title=tit)
    ax[0].set_ylabel(r"potencial reducido $\psi$")
    ax[0].plot([], [], "--", color="#55544f", label="hoja corregida (abierta)")
    ax[0].plot([], [], "-", color="#55544f", label="+ conservación de iones (cerrada)")
    ax[0].legend(loc="upper right", fontsize=8.5)
    fig.suptitle("Potencial de una carga entre placas a tierra con electrolito 1:1 (L = 32, caja lateral 32): Ludwig contra la solución analítica", y=1.02)
    fig.tight_layout(); fig.savefig(os.path.join(ROOT, "graficos", "C_hoja_potencial.png"), dpi=160, bbox_inches="tight")
