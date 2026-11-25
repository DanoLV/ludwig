#!/usr/bin/env python3
"""Prueba C (sin fondo): potencial sobre la linea de la particula y carga
inducida en las paredes, medidos, contra el sistema ideal de una carga puntual
entre dos placas a tierra INFINITAS (imagenes solo en z). Se agrega la teoria
periodica en x,y de la caja (Fourier) como referencia.
psi reducido = l_B * sum q_i/|r - r_i|, l_B = 10/(4 pi). Particula en x=y=16 (nodo).
Uso: python3 grafico_C_potencial_carga.py -> ../graficos/C_potencial_carga_vs_placas_infinitas.png

2026-11-16 (pedido del usuario): en el panel (a) la referencia teórica es ahora la carga repartida con el kernel
Hann-8 CONTINUO (psi_hann), no la carga puntual: misma serie de Fourier en x, y que psi_periodic (placas a tierra,
caja periódica 32x32, sin fondo), con la fuente transformada  W^(kx) W^(ky),  W^(k) = [sin(kn/2)/(kn/2)] / [1-(kn/2pi)^2]
(la misma de grafico_C_hoja_hann.py), y convolución en z con W(s - d):
   psi(z) = (4 pi l_B / A) sum_k W^(kx) W^(ky) int W(s-d) G_k(z, s) ds,
   G_k = [e^{-k|z-s|} - e^{-k(z+s)} - e^{-k(2L-z-s)} + e^{-k(2L-|z-s|)}] / [2k(1 - e^{-2kL})],  G_0 = z_<(L - z_>)/L.
Comprobado: lejos del soporte coincide con psi_periodic (puntual) a ~1e-6; M = 60 vs 80 modos: igual a 1e-7.
Aproximación: kernel continuo (Ludwig lo evalúa en los nodos). Textos de la figura en inglés (regla del usuario).
"""
import os, sys, glob
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = "/home/danolv/Sim/ludwig/capacitor/local/WALL-DIRICHLET"   # absolute: scripts/ is a symlink since 2026-11-13
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

NH = 8                                               # Hann-8 kernel


def What(k, n=NH):
    x = k * n / 2; s = np.where(np.abs(x) < 1e-12, 1.0, np.sin(x) / np.where(x == 0, 1, x))
    den = 1 - (k * n / (2 * np.pi)) ** 2
    return np.where(np.abs(den) < 1e-12, 0.5, s / np.where(np.abs(den) < 1e-12, 1, den))


def W(u, n=NH):
    return np.where(np.abs(u) <= n / 2, (1 + np.cos(2 * np.pi * u / n)) / n, 0.0)


def psi_hann(z, d, M=60, Ns=801):
    """Potential on the particle axis of the continuous Hann-8 charge between grounded plates, periodic in x, y."""
    m = np.arange(-M, M + 1); mx, my = np.meshgrid(m, m, indexing="ij")
    kx = 2 * np.pi * mx / NXY; ky = 2 * np.pi * my / NXY
    k = np.hypot(kx, ky).ravel(); wt = (What(kx) * What(ky)).ravel(); nz = k > 0
    ku, inv = np.unique(np.round(k[nz], 12), return_inverse=True); wu = np.bincount(inv, weights=wt[nz])
    s = np.linspace(d - NH / 2, d + NH / 2, Ns); ws = W(s - d) * (s[1] - s[0]); ws[[0, -1]] *= 0.5
    z = np.atleast_1d(z).astype(float)
    zl = np.minimum(z[:, None], s[None, :]); zg = np.maximum(z[:, None], s[None, :])
    out = (zl * (L - zg) / L) @ ws
    dz = np.abs(z[:, None] - s[None, :]); sp = z[:, None] + s[None, :]
    for i in range(0, len(ku), 40):
        kk = ku[i:i + 40][:, None, None]
        G = (np.exp(-kk * dz) - np.exp(-kk * sp) - np.exp(-kk * (2 * L - sp)) + np.exp(-kk * (2 * L - dz))) \
            / (2 * kk * (1 - np.exp(-2 * kk * L)))
        out += np.einsum("k,kzs,s->z", wu[i:i + 40], G, ws)
    return C.FOURPI_LB * out / A


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
    th = psi_hann(zz, d)
    a.plot(zz, th, "-", color=c, lw=1.2, alpha=0.8)
    a.plot(zc[fl & ~near], p[fl & ~near], "o", ms=4, color=c, label="d = %g" % d)
    a.plot(zc[fl & near], p[fl & near], "o", ms=4, mfc="none", color=c)
    pt = psi_hann(zc[fl], d); res = p[fl] - pt
    print("d = %g: max|Ludwig - Hann theory| outside / inside the kernel support: %.2e / %.2e"
          % (d, np.abs(res[~near[fl]]).max(), np.abs(res[near[fl]]).max()))
a.set(xlabel="z from the lower mid-plane (lattice units)", ylim=(0, 0.45),
      ylabel=r"reduced potential $\psi$ on the particle axis",
      title="(a) Potential")
a.plot([], [], "-", color="#55544f", label="theory: Hann-8 charge, periodic in x, y")
a.plot([], [], "o", mfc="none", color="#55544f", label="inside the kernel support")
a.legend(loc="upper right", fontsize=8.5)

b = ax[1]
dirs = sorted(glob.glob(os.path.join(BAT, "z_*")), key=lambda p: float(p.rsplit("_", 1)[1]))
d = np.array([float(p.rsplit("_", 1)[1]) - 1.5 for p in dirs])
r = [C.read_run(p) for p in dirs]
qlo = np.array([x["qlo"] for x in r]); qhi = np.array([x["qhi"] for x in r])
dd = np.linspace(0.5, 16.5, 200)
b.axvspan(0, C.D_MIN, color="#eee8dc", zorder=0)
b.plot(dd, -(L - dd) / L, "--", color="#55544f", label="theory: $-q(L-d)/L$, $-qd/L$")
b.plot(dd, -dd / L, "--", color="#55544f")
ok = d >= C.D_MIN
b.plot(d[ok], qlo[ok], "o", ms=5, color="#2a78d6", label="lower wall, Ludwig")
b.plot(d[ok], qhi[ok], "^", ms=5, color="#eb6834", label="upper wall, Ludwig")
b.plot(d[~ok], qlo[~ok], "o", ms=5, mfc="none", color="#2a78d6")
b.plot(d[~ok], qhi[~ok], "^", ms=5, mfc="none", color="#eb6834")
b.set(xlabel="distance d to the lower wall (lattice units)", ylabel="induced charge (units of e)",
      title="(b) Charge induced on each wall")
b.legend(loc="upper right", fontsize=8.5, bbox_to_anchor=(1.0, 0.88))
fig.suptitle("One charge q = +1 between two grounded walls, L = 32, no counterion background (Hann-8, box 32×32)", y=1.0)
fig.tight_layout()
out = os.path.join(ROOT, "graficos", "C_potencial_carga_vs_placas_infinitas.png"); fig.savefig(out, dpi=160)
print("max |Q medida - ideal| (kernel libre): inf %.2e sup %.2e" %
      (np.max(np.abs(qlo[ok] + (L - d[ok]) / L)), np.max(np.abs(qhi[ok] + d[ok] / L))))
