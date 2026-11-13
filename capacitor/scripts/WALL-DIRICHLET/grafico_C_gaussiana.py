#!/usr/bin/env python3
"""Prueba C (sin fondo): potencial sobre la linea de la particula y carga
inducida, medidos, contra una carga GAUSSIANA entre placas a tierra infinitas.
Gaussiana isotropa con la misma varianza que el kernel Hann-n:
  sigma^2 = n^2 (1/12 - 1/(2 pi^2)) -> sigma = 1.446 para n = 8.
Potencial de una gaussiana en el vacio: psi = l_B q erf(r/(sqrt2 sigma))/r
(valor central l_B q sqrt(2/pi)/sigma). Imagenes: la misma gaussiana, +q en
d+2nL y -q en -d+2nL. Referencia adicional: caja periodica (Fourier, cada modo
lateral multiplicado por exp(-k^2 sigma^2/2), convolucion con la gaussiana en z).
Uso: python3 grafico_C_gaussiana.py -> ../graficos/C_gaussiana_vs_placas_infinitas.png
"""
import os, sys, glob, math
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
BAT = os.path.join(ROOT, "C-sin-fondo-hann08")
sys.argv = ["x", BAT, "8"]; sys.path.insert(0, HERE)
import analizar_C as C
NH = 8
SIG = NH * math.sqrt(1 / 12 - 1 / (2 * math.pi ** 2))
LB = C.FOURPI_LB / (4 * np.pi); L = C.L; NZ = 34; NXY = 32; A = float(NXY * NXY)
from scipy.special import erf as verf                 # vectorised erf (scipy.special)

def kern(r):                      # erf(r/(sqrt2 sig))/r, finito en 0
    r = np.abs(np.asarray(r, float)); out = np.empty_like(r)
    s = r < 1e-9
    out[s] = math.sqrt(2 / math.pi) / SIG
    out[~s] = verf(r[~s] / (math.sqrt(2) * SIG)) / r[~s]
    return out

def psi_ideal(z, d, nmax=3000):
    z = np.atleast_1d(np.asarray(z, float)); n = np.arange(-nmax, nmax + 1)
    out = np.zeros_like(z)
    for i, zz in enumerate(z):
        zp = d + 2 * L * n; zn = -d + 2 * L * n
        rp, rn = np.abs(zz - zp), np.abs(zz - zn)
        pos = np.where(rp < 10 * SIG, kern(np.minimum(rp, 1e3)), 1 / np.maximum(rp, 1e-9))
        neg = np.where(rn < 10 * SIG, kern(np.minimum(rn, 1e3)), 1 / np.maximum(rn, 1e-9))
        out[i] = LB * (pos.sum() - neg.sum())
    return out

def psi_periodic(z, d, M=16):
    m = np.arange(-M, M + 1); mx, my = np.meshgrid(m, m, indexing="ij")
    k = 2 * np.pi * np.hypot(mx, my) / NXY
    w = np.exp(-k ** 2 * SIG ** 2 / 2)
    zs = d + np.linspace(-5 * SIG, 5 * SIG, 41); zs = zs[(zs > 0) & (zs < L)]
    g = np.exp(-(zs - d) ** 2 / (2 * SIG ** 2)); g /= g.sum()
    z = np.atleast_1d(z); out = np.zeros_like(z, dtype=float)
    kz = k[k > 0]; wz = w[k > 0]
    for i, zz in enumerate(z):
        zl = np.minimum(zz, zs); zg = np.maximum(zz, zs)
        # k = 0
        s0 = np.sum(g * zl * (L - zg) / L)
        sk = np.sum(g[None, :] * np.sinh(kz[:, None] * zl[None, :]) *
                    np.sinh(kz[:, None] * (L - zg[None, :])) / (kz[:, None] * np.sinh(kz[:, None] * L)) * wz[:, None])
        out[i] = C.FOURPI_LB * (s0 + sk) / A
    return out

plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.color": "#e4e2dc",
                     "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False})
fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.8))
cols = {4.0: "#eb6834", 8.0: "#2a78d6", 12.0: "#1baf7a"}
zc = np.arange(1, NZ + 1) - 1.5; zz = np.linspace(0, L, 400)
a = ax[0]; rms = []
for d, c in cols.items():
    p = np.loadtxt(os.path.join(BAT, "z_%05.2f" % (d + 1.5), "psi-000000001.001-001")).reshape(NXY, NXY, NZ)[15, 15, :]
    fl = (zc > 0) & (zc < L)
    a.plot(zz, psi_ideal(zz, d), "--", color=c, lw=1.2)
    a.plot(zz, psi_periodic(zz, d), ":", color=c, lw=1.4)
    a.plot(zc[fl], p[fl], "o", ms=4, color=c, label="d = %g" % d)
    ref = psi_periodic(zc[fl], d)
    rms.append((d, np.max(np.abs(p[fl] - ref)) / np.max(ref)))
a.set(xlabel="z medido desde el plano medio inferior (unidades de red)", ylim=(0, 0.5),
      ylabel=r"potencial reducido $\psi$ sobre la línea de la partícula", title="(a) Potencial")
a.plot([], [], "--", color="#55544f", label="gaussiana σ = %.3f, placas infinitas" % SIG)
a.plot([], [], ":", color="#55544f", label="gaussiana, con periodicidad lateral")
a.legend(loc="upper right", fontsize=8.5)

b = ax[1]
dirs = sorted(glob.glob(os.path.join(BAT, "z_*")), key=lambda p: float(p.rsplit("_", 1)[1]))
d = np.array([float(p.rsplit("_", 1)[1]) - 1.5 for p in dirs]); r = [C.read_run(p) for p in dirs]
qlo = np.array([x["qlo"] for x in r]); qhi = np.array([x["qhi"] for x in r]); ok = d >= C.D_MIN
dd = np.linspace(0.5, 16.5, 200); b.axvspan(0, C.D_MIN, color="#eee8dc", zorder=0)
b.plot(dd, -(L - dd) / L, "--", color="#55544f", label="gaussiana entre placas infinitas: $-q(L-d)/L$, $-qd/L$")
b.plot(dd, -dd / L, "--", color="#55544f")
b.plot(d[ok], qlo[ok], "o", ms=5, color="#2a78d6", label="pared inferior, Ludwig")
b.plot(d[ok], qhi[ok], "^", ms=5, color="#eb6834", label="pared superior, Ludwig")
b.plot(d[~ok], qlo[~ok], "o", ms=5, mfc="none", color="#2a78d6")
b.plot(d[~ok], qhi[~ok], "^", ms=5, mfc="none", color="#eb6834")
b.set(xlabel="distancia d a la pared inferior (unidades de red)", ylabel="carga inducida (unidades de e)",
      title="(b) Carga inducida en cada pared")
b.legend(loc="center", fontsize=8.5, bbox_to_anchor=(0.55, 0.5))
fig.suptitle("Carga q = +1 gaussiana (σ del kernel Hann-8) entre dos paredes a tierra, L = 32, sin fondo", y=1.0)
fig.tight_layout()
fig.savefig(os.path.join(ROOT, "graficos", "C_gaussiana_vs_placas_infinitas.png"), dpi=160)
print("sigma = %.4f ; centro gaussiana = %.4f" % (SIG, LB * math.sqrt(2 / math.pi) / SIG))
for dd_, e in rms: print("d=%g: max|psi_Ludwig - gaussiana periodica| / max(psi) = %.3f" % (dd_, e))
