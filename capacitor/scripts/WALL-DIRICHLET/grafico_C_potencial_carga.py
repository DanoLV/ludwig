#!/usr/bin/env python3
"""Prueba C (sin fondo): potencial sobre la linea de la particula y carga
inducida en las paredes, medidos, contra el sistema ideal de una carga puntual
entre dos placas a tierra INFINITAS (imagenes solo en z). Se agrega la teoria
periodica en x,y de la caja (Fourier) como referencia.
psi reducido = l_B * sum q_i/|r - r_i|, l_B = 10/(4 pi). Particula en x=y=16 (nodo).
Uso: python3 grafico_C_potencial_carga.py -> ../graficos/C_potencial_carga_vs_placas_infinitas.png
"""
import os, sys, glob
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
BAT = os.path.join(ROOT, "C-sin-fondo-hann08")
sys.argv = ["x", BAT, "8"]; sys.path.insert(0, HERE)
import analizar_C as C
LB = C.FOURPI_LB / (4 * np.pi); L = C.L; NZ = 34; NXY = 32; A = float(NXY * NXY)

def psi_ideal(z, d, nmax=4000):
    """z, d medidos desde el plano medio inferior."""
    n = np.arange(-nmax, nmax + 1)
    zp = d + 2 * L * n; zn = -d + 2 * L * n
    z = np.atleast_1d(z)[:, None]
    return LB * (np.sum(1 / np.abs(z - zp), axis=1) - np.sum(1 / np.abs(z - zn), axis=1))

def psi_periodic(z, d, M=60):
    m = np.arange(-M, M + 1); mx, my = np.meshgrid(m, m, indexing="ij")
    k = 2 * np.pi * np.hypot(mx, my) / NXY; k = k[k > 0]
    z = np.atleast_1d(z); out = np.zeros_like(z, dtype=float)
    for i, zz in enumerate(z):
        zl, zg = min(zz, d), max(zz, d)
        s = np.sum(np.sinh(k * zl) * np.sinh(k * (L - zg)) / (k * np.sinh(k * L)))
        s += zl * (L - zg) / L
        out[i] = C.FOURPI_LB * s / A
    return out

def read_psi(path):
    return np.loadtxt(path).reshape(NXY, NXY, NZ)

plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.color": "#e4e2dc",
                     "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False})
fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.8))
cols = {4.0: "#eb6834", 8.0: "#2a78d6", 12.0: "#1baf7a"}
a = ax[0]
zc = np.arange(1, NZ + 1) - 1.5                      # distancia al plano medio inferior
zz = np.linspace(0, L, 600)
for d, c in cols.items():
    p = read_psi(os.path.join(BAT, "z_%05.2f" % (d + 1.5), "psi-000000001.001-001"))[15, 15, :]
    fl = (zc > 0) & (zc < L)
    near = np.abs(zc - d) < 4.5                      # soporte del kernel Hann-8: no comparable
    a.plot(zz, np.clip(psi_ideal(zz, d), None, 0.5), "--", color=c, lw=1.2)
    a.plot(zz, np.clip(psi_periodic(zz, d), None, 0.5), ":", color=c, lw=1.4)
    a.plot(zc[fl & ~near], p[fl & ~near], "o", ms=4, color=c, label="d = %g" % d)
    a.plot(zc[fl & near], p[fl & near], "o", ms=4, mfc="none", color=c)
a.set(xlabel="z medido desde el plano medio inferior (unidades de red)", ylim=(0, 0.45),
      ylabel=r"potencial reducido $\psi$ sobre la línea de la partícula",
      title="(a) Potencial")
a.plot([], [], "--", color="#55544f", label="ideal: placas infinitas")
a.plot([], [], ":", color="#55544f", label="teoría con periodicidad lateral")
a.plot([], [], "o", mfc="none", color="#55544f", label="dentro del soporte del kernel")
a.legend(loc="upper right", fontsize=8.5)

b = ax[1]
dirs = sorted(glob.glob(os.path.join(BAT, "z_*")), key=lambda p: float(p.rsplit("_", 1)[1]))
d = np.array([float(p.rsplit("_", 1)[1]) - 1.5 for p in dirs])
r = [C.read_run(p) for p in dirs]
qlo = np.array([x["qlo"] for x in r]); qhi = np.array([x["qhi"] for x in r])
dd = np.linspace(0.5, 16.5, 200)
b.axvspan(0, C.D_MIN, color="#eee8dc", zorder=0)
b.plot(dd, -(L - dd) / L, "--", color="#55544f", label="ideal (placas infinitas): $-q(L-d)/L$, $-qd/L$")
b.plot(dd, -dd / L, "--", color="#55544f")
ok = d >= C.D_MIN
b.plot(d[ok], qlo[ok], "o", ms=5, color="#2a78d6", label="pared inferior, Ludwig")
b.plot(d[ok], qhi[ok], "^", ms=5, color="#eb6834", label="pared superior, Ludwig")
b.plot(d[~ok], qlo[~ok], "o", ms=5, mfc="none", color="#2a78d6")
b.plot(d[~ok], qhi[~ok], "^", ms=5, mfc="none", color="#eb6834")
b.set(xlabel="distancia d a la pared inferior (unidades de red)", ylabel="carga inducida (unidades de e)",
      title="(b) Carga inducida en cada pared")
b.legend(loc="upper right", fontsize=8.5, bbox_to_anchor=(1.0, 0.88))
fig.suptitle("Una carga q = +1 entre dos paredes a tierra, L = 32, sin fondo de contraiones (Hann-8, caja 32×32)", y=1.0)
fig.tight_layout()
out = os.path.join(ROOT, "graficos", "C_potencial_carga_vs_placas_infinitas.png"); fig.savefig(out, dpi=160)
print("max |Q medida - ideal| (kernel libre): inf %.2e sup %.2e" %
      (np.max(np.abs(qlo[ok] + (L - d[ok]) / L)), np.max(np.abs(qhi[ok] + d[ok] / L))))
