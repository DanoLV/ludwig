#!/usr/bin/env python3
"""Prueba C (sin fondo): potencial sobre la linea de la particula y carga
inducida, medidos, contra una carga GAUSSIANA entre placas a tierra infinitas.
Gaussiana isotropa con la misma varianza que el kernel Hann-n:
  sigma^2 = n^2 (1/12 - 1/(2 pi^2)) -> sigma = 1.446 para n = 8.
Potencial de una gaussiana en el vacio: psi = l_B q erf(r/(sqrt2 sigma))/r
(valor central l_B q sqrt(2/pi)/sigma). Imagenes: la misma gaussiana, +q en
d+2nL y -q en -d+2nL. Referencia adicional: caja periodica (Fourier, cada modo
lateral multiplicado por exp(-k^2 sigma^2/2), convolucion con la gaussiana en z).
Paredes a DOS potenciales (V1 abajo, V2 arriba): por linealidad, psi = V1 + (V2-V1) z/L +
psi_carga(paredes a tierra). Corridas reales en C-dospotenciales-hann08/.
Uso: python3 grafico_C_campo_dospotenciales.py -> ../graficos/C_campo_dospotenciales.png
"""
import os, sys, glob, math
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
BAT = os.path.join(ROOT, "C-dospotenciales-hann08")
V1, V2 = 0.1, 0.3
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
    zs = d + np.linspace(-5 * SIG, 5 * SIG, 201); zs = zs[(zs > 0) & (zs < L)]
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



lin = lambda z: V1 + (V2 - V1) * np.asarray(z) / L
E0 = -(V2 - V1) / L                       # campo uniforme de las paredes
plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.color": "#e4e2dc",
                     "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False})
fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.8))
cols = {4.0: "#eb6834", 8.0: "#2a78d6", 12.0: "#1baf7a"}
zc = np.arange(1, NZ + 1) - 1.5; zz = np.linspace(0.3, L - 0.3, 500)
a = ax[0]; a.axhline(E0, color="#999", lw=1, label="sin carga: $-(V_2-V_1)/L$")
h = 0.05
for d, c in cols.items():
    p = np.loadtxt(os.path.join(BAT, "z_%05.2f" % (d + 1.5), "psi-000000001.001-001")).reshape(NXY, NXY, NZ)[15, 15, :]
    Em = -(p[2:] - p[:-2]) / 2.0; zm = zc[1:-1]; fl = (zm > 0.5) & (zm < L - 0.5)
    ps = lambda z, f: f(z, d)
    Ei = -(ps(zz + h, psi_ideal) - ps(zz - h, psi_ideal)) / (2 * h) + E0
    Ep = -(ps(zz + h, psi_periodic) - ps(zz - h, psi_periodic)) / (2 * h) + E0
    a.plot(zz, Ei, "--", color=c, lw=1.2); a.plot(zz, Ep, ":", color=c, lw=1.4)
    a.plot(zm[fl], Em[fl], "o", ms=4, color=c, label="Ludwig, d = %g" % d)
a.plot([], [], "--", color="#55544f", label="gaussiana σ = %.3f, placas infinitas" % SIG)
a.plot([], [], ":", color="#55544f", label="gaussiana, con periodicidad lateral")
a.set(xlabel="z medido desde el plano medio inferior (unidades de red)", ylim=(-0.2, 0.2),
      ylabel=r"$E_z=-\partial_z\psi$ sobre la línea de la partícula", title="(a) Perfil del campo")
a.legend(loc="upper right", fontsize=8.5)

b = ax[1]
dirs = sorted(glob.glob(os.path.join(BAT, "z_*")), key=lambda p: float(p.rsplit("_", 1)[1]))
d = np.array([float(p.rsplit("_", 1)[1]) - 1.5 for p in dirs]); ez = np.array([C.read_run(p)["ez"] for p in dirs])
ok = d >= C.D_MIN
def E_ideal(dd_):
    n = np.arange(1, 5001); pos = dd_ + 2 * L * np.concatenate([-n, n]); neg = -dd_ + 2 * L * np.arange(-5000, 5001)
    f = lambda z: (dd_ - z) / np.abs(dd_ - z) ** 3
    return LB * (f(pos).sum() - f(neg).sum())
dd = np.linspace(0.8, 16, 300); b.axvspan(0, C.D_MIN, color="#eee8dc", zorder=0)
b.axhline(E0, color="#999", lw=1)
b.plot(dd, [E_ideal(x) + E0 for x in dd], "--", color="#55544f", lw=1.5, label="ideal: carga entre placas infinitas + campo uniforme")
b.plot(dd, [C.image_field(x, k0=True) + E0 for x in dd], ":", color="#1baf7a", lw=1.6, label="con periodicidad lateral (Fourier) + campo uniforme")
b.plot(d[ok], ez[ok], "o", ms=5, color="#2a78d6", label="Ludwig, fuerza sobre la partícula / q")
b.plot(d[~ok], ez[~ok], "s", ms=6, mfc="none", color="#eb6834", label="Ludwig, kernel alcanza la pared")
b.set(xlabel="distancia d a la pared inferior (unidades de red)", ylim=(-0.05, 0.003),
      ylabel=r"$E_z$ sobre la partícula", title="(b) Campo sobre la partícula, paredes a ψ = %g y %g" % (V1, V2))
b.legend(loc="lower right", fontsize=8.5)
fig.tight_layout()
fig.savefig(os.path.join(ROOT, "graficos", "C_campo_dospotenciales.png"), dpi=160)
r = np.abs(ez[ok] - np.array([C.image_field(x, k0=True) + E0 for x in d[ok]])) / np.abs(ez[ok])
print("error relativo vs teoria periodica (kernel libre): mediana %.1e, max %.1e" % (np.median(r), r.max()))
