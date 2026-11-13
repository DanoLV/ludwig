#!/usr/bin/env python3
"""Prueba C (sin fondo): campo medido sobre la particula vs el del sistema ideal
de una carga puntual entre dos placas a tierra INFINITAS (sin periodicidad
lateral): serie de imagenes a lo largo de z, +q en d+2nL (n!=0), -q en -d+2nL.
E_z = l_B * sum q_i (d - z_i)/|d - z_i|^3, unidades de Ludwig (4 pi l_B = 10).
Uso: python3 grafico_C_placas_infinitas.py   -> ../graficos/C_vs_placas_infinitas.png
"""
import os, sys, glob
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
BAT = os.path.join(ROOT, "C-sin-fondo-hann08")
sys.argv = ["x", BAT, "8"]
sys.path.insert(0, HERE)
import analizar_C as C
LB = C.FOURPI_LB / (4 * np.pi); L = C.L; NMAX = 5000

def ideal(d):
    n = np.arange(1, NMAX + 1)
    # imagenes +q en d+2nL (n != 0) y -q en -d+2nL (todo n); campo sobre la carga en z=d
    pos = d + 2 * L * np.concatenate([-n, n]); neg = -d + 2 * L * np.arange(-NMAX, NMAX + 1)
    f = lambda z: (d - z) / np.abs(d - z) ** 3
    return LB * (np.sum(f(pos)) - np.sum(f(neg)))

dirs = sorted(glob.glob(os.path.join(BAT, "z_*")), key=lambda p: float(p.rsplit("_", 1)[1]))
d = np.array([float(p.rsplit("_", 1)[1]) - 1.5 for p in dirs])
ez = np.array([C.read_run(p)["ez"] for p in dirs])
th = np.array([ideal(x) for x in d]); ok = d >= C.D_MIN

plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.color": "#e4e2dc",
                     "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False})
fig, ax = plt.subplots(2, 1, figsize=(7.5, 7), sharex=True, gridspec_kw={"height_ratios": [2.2, 1]})
dd = np.linspace(0.8, 16, 400)
a = ax[0]; a.axvspan(0, C.D_MIN, color="#eee8dc", zorder=0)
a.plot(dd, [ideal(x) for x in dd], "--", color="#55544f", lw=1.5,
       label="ideal: carga puntual entre placas infinitas (imágenes en z)")
a.plot(dd, [C.image_field(x, k0=True) for x in dd], ":", color="#1baf7a", lw=1.6,
       label="teoría con la periodicidad lateral de la caja (serie de Fourier)")
a.plot(d[ok], ez[ok], "o", ms=5, color="#2a78d6", label="Ludwig (Hann-8, caja periódica 32×32)")
a.plot(d[~ok], ez[~ok], "s", ms=6, mfc="none", color="#eb6834", label="Ludwig, kernel alcanza la pared")
a.set(ylabel=r"$E_z$ sobre la partícula (unidades reducidas)", ylim=(-0.045, 0.003),
      title="Una carga q = +1 entre dos paredes a tierra, L = 32 (sin fondo de contraiones)")
a.legend(loc="lower right")
b = ax[1]; b.axvspan(0, C.D_MIN, color="#eee8dc", zorder=0)
rel = np.abs(ez - th) / np.abs(th)
b.semilogy(d[ok], rel[ok], "o", ms=5, color="#2a78d6")
b.semilogy(d[~ok], rel[~ok], "s", ms=6, mfc="none", color="#eb6834")
b.set(xlabel="distancia d a la pared inferior (unidades de red)", ylabel="|medido − ideal| / |ideal|")
fig.tight_layout()
out = os.path.join(ROOT, "graficos", "C_vs_placas_infinitas.png"); fig.savefig(out, dpi=160)
for x, e, t in zip(d, ez, th): print("d=%5.1f  Ludwig=%+.5e  ideal=%+.5e  rel=%.1e" % (x, e, t, abs(e-t)/abs(t) if t else 0))
