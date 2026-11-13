#!/usr/bin/env python3
"""Prueba C, particula con el soporte del kernel cruzando la pared (d < n/2 - 1/2 = 3.5, Hann-8):
Ludwig (C-sin-fondo-hann08) contra DOS tratamientos teoricos del kernel que cruza el plano medio:
  T0 "truncado": solo los nodos fluidos reciben carga y fuerza; la carga que cae en nodos de pared se pierde.
  T1 "imagen espejada" (DISCOS, Ladiges et al., PRE 106, 035104 (2022), Sec. II.B.1): se agrega una particula
     imagen de carga opuesta espejada en el plano medio, cuyos pesos en el dominio se restan a los de la real;
     el campo se interpola con celdas fantasma pares, es decir F = q sum_n [w(z_n-d) + w(z_n+d)] E(z_n).
Los dos usan el kernel Hann-8 MUESTREADO en los nodos (x, y, z), Fourier lateral con ell = 32 y G_k con Dirichlet en
z = 0 y L = 32. Carga inducida: reciprocidad de Green (Jackson 1.13) con las cargas nodales de cada tratamiento.
Uso: python3 grafico_C_subkernel.py -> ../graficos/C_subkernel_pared.png (y tabla por stdout)
"""
import os, sys, glob
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
BAT = os.path.join(ROOT, "C-sin-fondo-hann08")
sys.argv = ["x", BAT, "8"]; sys.path.insert(0, HERE)
import analizar_C as C

L = 32.0; N = 8; LB4 = C.FOURPI_LB; ELL = 32; A = ELL * ELL
w = lambda r: np.where(np.abs(np.asarray(r, float)) < N / 2, (1 + np.cos(2 * np.pi * np.asarray(r, float) / N)) / N, 0.0)
nn = np.arange(-3, 4); wn = w(nn)                         # kernel muestreado, particula sobre un nodo en x, y
m = np.arange(0, ELL // 2 + 1)
whd = np.array([(wn * np.cos(2 * np.pi * mi * nn / ELL)).sum() for mi in m])
mult = np.where((m > 0) & (m < ELL // 2), 2.0, 1.0)
mx, my = np.meshgrid(m, m, indexing="ij")
kk = (2 * np.pi * np.hypot(mx, my) / ELL).ravel()
w2 = ((mult[:, None] * mult[None, :]) * (whd[:, None] * whd[None, :]) ** 2).ravel()
zn = np.arange(32) + 0.5                                   # nodos fluidos, medidos desde el plano medio inferior

def dG(k, z, zp):
    kq = np.where(k > 0, k, 1.0)
    below = np.cosh(kq * z) * np.sinh(kq * (L - zp)) / np.sinh(kq * L)
    above = -np.sinh(kq * zp) * np.cosh(kq * (L - z)) / np.sinh(kq * L)
    eq = np.abs(z - zp) < 1e-9                             # diagonal: promedio de lados (el salto es el autocampo)
    d = np.where(eq, 0.5 * (below + above), np.where(z < zp, below, above))
    d0 = np.where(eq, 0.5 * ((L - zp) / L - zp / L), np.where(z < zp, (L - zp) / L, -zp / L))
    return np.where(k > 0, d, d0)
S = dG(kk[:, None, None], zn[None, :, None], zn[None, None, :])
def force(src, fw):
    return -LB4 / A * np.einsum("k,ka,a->", w2, np.einsum("kab,b->ka", S, src), fw)
def treat(d):
    wr, wi = w(zn - d), w(zn + d)
    return {"T0": (wr, wr), "T1": (wr - wi, wr + wi)}      # (carga nodal, pesos de fuerza)

dirs = sorted(glob.glob(os.path.join(BAT, "z_*")), key=lambda p: float(p.rsplit("_", 1)[1]))
rows = []
for p in dirs:
    d = float(p.rsplit("_", 1)[1]) - 1.5
    if d > 8.01: break
    r = C.read_run(p); t = treat(d)
    row = dict(d=d, ez=r["ez"], qlo=r["qlo"], qhi=r["qhi"])
    for k, (src, fw) in t.items():
        row[k] = force(src, fw); row[k + "_qlo"] = -(src * (L - zn) / L).sum(); row[k + "_qhi"] = -(src * zn / L).sum()
        row[k + "_qeff"] = src.sum()
    rows.append(row)
R = {k: np.array([r[k] for r in rows]) for k in rows[0]}
print("   d   Ez_Ludwig     T0/Ludwig  T1/Ludwig | q_eff(T0) q_eff(T1) | Q_inf: Ludwig   T0        T1")
for i, d in enumerate(R["d"]):
    print("%5.1f %+.4e   %.4f     %.4f   |  %.4f    %.4f   |  %+.4f %+.4f %+.4f" % (
        d, R["ez"][i], R["T0"][i] / R["ez"][i], R["T1"][i] / R["ez"][i], R["T0_qeff"][i], R["T1_qeff"][i],
        R["qlo"][i], R["T0_qlo"][i], R["T1_qlo"][i]))

plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.color": "#e4e2dc",
                     "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False})
fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.5))
a = ax[0]; a.axvspan(0, C.D_MIN, color="#eee8dc", zorder=0)
df = np.linspace(R["d"].min(), R["d"].max(), 160)       # teoria en grilla fina
tt = [treat(x) for x in df]
a.plot(df, [force(*t["T0"]) for t in tt], "-", color="#55544f", lw=1.6, label="T0: truncado (carga en pared perdida)")
a.plot(df, [force(*t["T1"]) for t in tt], "--", color="#eb6834", lw=1.6, label="T1: imagen espejada (DISCOS 2022)")
a.plot(R["d"], R["ez"], "o", ms=5, color="#2a78d6", label="Ludwig")
a.set(xlabel="distancia d al plano medio de la pared (unidades de red)", ylabel=r"$E_z$ sobre la partícula",
      title="(a) Campo sobre la partícula, Hann-8"); a.legend(loc="lower right")
b = ax[1]; b.axvspan(0, C.D_MIN, color="#eee8dc", zorder=0); b.axhline(1, color="#999", lw=1)
b.plot(R["d"], R["T0"] / R["ez"], "o-", color="#55544f", label="T0 / Ludwig")
b.plot(R["d"], R["T1"] / R["ez"], "s--", color="#eb6834", label="T1 / Ludwig")
b.set(xlabel="distancia d al plano medio de la pared", ylabel="teoría / Ludwig", title="(b) Cociente (sombreado: el soporte llega a nodos de pared)")
b.legend(loc="upper right")
fig.tight_layout(); fig.savefig(os.path.join(ROOT, "graficos", "C_subkernel_pared.png"), dpi=160)
