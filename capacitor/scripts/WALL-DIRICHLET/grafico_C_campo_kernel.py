#!/usr/bin/env python3
"""Prueba C, paredes a DOS potenciales (C-dospotenciales-hann08): campo medido
contra la teoria cuyo objetivo es la FORMA REAL DEL KERNEL, no una puntual ni
una gaussiana: carga q w(x) w(y) w(z), w(r) = (1/n)(1 + cos(2 pi r/n)), |r| <= n/2
(Hann-n, src/subgrid.c: d_hann(x) d_hann(y) d_hann(z)), n = 8.

Teoria (Fourier lateral, Dirichlet en z, constante de Green de la sec. C del README):
  psi(r) = (4 pi l_B q / A) sum_k e^{ik.r} wh(kx) wh(ky) int G_k(z,z') w(z') dz'
  wh(k) = int w(x) cos(kx) dx = sinc(u)/(1-u^2), u = k n/(2 pi)
  G_k(z,z') = sinh(k z<) sinh(k(L-z>)) / (k sinh kL),  G_0 = z<(L-z>)/L
Campo sobre la distribucion (fuerza/q; el autocampo se cancela por simetria):
  E = E0 - (4 pi l_B/A) sum_k wh(kx)^2 wh(ky)^2 II w(z) w(z') d_z G_k(z,z') dz dz'
Caja de lado ell = 32 (la de la simulacion) o ell = 200 (placas infinitas en la
practica: la imagen lateral decae exponencialmente con ell/L).
Uso: python3 grafico_C_campo_kernel.py -> ../graficos/C_campo_dospotenciales_kernel.png
"""
import os, sys, glob, math
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
BAT = os.path.join(ROOT, "C-dospotenciales-hann08")
sys.argv = ["x", BAT, "8"]; sys.path.insert(0, HERE)
import analizar_C as C
NH = 8; V1, V2 = 0.1, 0.3
L = C.L; NZ = 34; NXY = 32
LB4 = C.FOURPI_LB                       # 4 pi l_B
E0 = -(V2 - V1) / L
lin = lambda z: V1 + (V2 - V1) * np.asarray(z) / L
s_ = np.linspace(-NH / 2, NH / 2, 41); ds = s_[1] - s_[0]
wz = (1 + np.cos(2 * np.pi * s_ / NH)) / NH; wz = wz / wz.sum()      # pesos de cuadratura (suma 1)

def what(k):
    u = k * NH / (2 * np.pi)
    den = 1 - u ** 2
    safe = np.abs(den) > 1e-9
    return np.where(safe, np.sinc(u) / np.where(safe, den, 1.0), 0.5)

def modes(ell):
    M = int(3.2 * ell / (2 * np.pi)) + 1
    m = np.arange(0, M + 1)
    mult = np.where(m > 0, 2.0, 1.0)
    mx, my = np.meshgrid(m, m, indexing="ij")
    k = 2 * np.pi * np.hypot(mx, my) / ell
    mm = mult[:, None] * mult[None, :]
    wxy = what(2 * np.pi * mx / ell) * what(2 * np.pi * my / ell)
    # w1: potencial en el eje (un factor wh por direccion); w2: campo sobre la carga (dos)
    return k.ravel(), (mm * wxy ** 2).ravel(), ell ** 2, (mm * wxy).ravel()

def Gk(k, zl, zg):
    kk = np.where(k > 0, k, 1.0)
    g = np.sinh(kk * zl) * np.sinh(kk * (L - zg)) / (kk * np.sinh(kk * L))
    return np.where(k > 0, g, zl * (L - zg) / L)

def dGk(k, z, zp):                       # d/dz de G_k(z,z')
    kk = np.where(k > 0, k, 1.0)
    below = np.cosh(kk * z) * np.sinh(kk * (L - zp)) / np.sinh(kk * L)
    above = -np.sinh(kk * zp) * np.cosh(kk * (L - z)) / np.sinh(kk * L)
    eq = np.abs(z - zp) < 1e-9          # diagonal: promedio de los dos lados (el salto es el autocampo)
    d = np.where(eq, 0.5 * (below + above), np.where(z < zp, below, above))
    d0 = np.where(eq, 0.5 * ((L - zp) / L - zp / L), np.where(z < zp, (L - zp) / L, -zp / L))
    return np.where(k > 0, d, d0)

def psi_q(zs, d, mod):
    k, w2, A, w1 = mod; zp = d + s_
    out = []
    for z in zs:
        G = Gk(k[:, None], np.minimum(z, zp)[None, :], np.maximum(z, zp)[None, :])
        out.append(LB4 / A * np.sum(w1 * (G @ wz)))
    return np.array(out)

def E_part(d, mod):
    k, w2, A, w1 = mod; zp = d + s_
    tot = 0.0
    for i, z in enumerate(zp):                       # int w(z) d_z psi
        tot += wz[i] * np.sum(w2 * (dGk(k[:, None], z, zp[None, :]) @ wz))
    return E0 - LB4 / A * tot

if __name__ == "__main__":
    P_box, P_inf = modes(32.0), modes(200.0)
    plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.color": "#e4e2dc",
                         "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False})
    fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.8))
    cols = {4.0: "#eb6834", 8.0: "#2a78d6", 12.0: "#1baf7a"}
    zc = np.arange(1, NZ + 1) - 1.5; zz = np.linspace(0.4, L - 0.4, 160); h = 0.05
    a = ax[0]; a.axhline(E0, color="#999", lw=1, label="sin carga: $-(V_2-V_1)/L$")
    for d, c in cols.items():
        p = np.loadtxt(os.path.join(BAT, "z_%05.2f" % (d + 1.5), "psi-000000001.001-001")).reshape(NXY, NXY, NZ)[15, 15, :]
        Em = -(p[2:] - p[:-2]) / 2.0; zm = zc[1:-1]; fl = (zm > 0.5) & (zm < L - 0.5)
        for P, st, lw in ((P_inf, "--", 1.2), (P_box, ":", 1.5)):
            Ef = -(psi_q(zz + h, d, P) - psi_q(zz - h, d, P)) / (2 * h) + E0
            a.plot(zz, Ef, st, color=c, lw=lw)
        a.plot(zm[fl], Em[fl], "o", ms=4, color=c, label="Ludwig, d = %g" % d)
    a.plot([], [], "--", color="#55544f", label="kernel Hann-8, placas infinitas (ℓ = 200)")
    a.plot([], [], ":", color="#55544f", label="kernel Hann-8, caja periódica ℓ = 32")
    a.set(xlabel="z medido desde el plano medio inferior (unidades de red)", ylim=(-0.2, 0.2),
          ylabel=r"$E_z=-\partial_z\psi$ sobre la línea de la partícula", title="(a) Perfil del campo")
    a.legend(loc="upper right", fontsize=8.5)

    b = ax[1]
    dirs = sorted(glob.glob(os.path.join(BAT, "z_*")), key=lambda p: float(p.rsplit("_", 1)[1]))
    d = np.array([float(p.rsplit("_", 1)[1]) - 1.5 for p in dirs]); ez = np.array([C.read_run(p)["ez"] for p in dirs])
    ok = d >= C.D_MIN
    dd = np.linspace(4.0, 16, 25)
    b.axvspan(0, C.D_MIN, color="#eee8dc", zorder=0); b.axhline(E0, color="#999", lw=1)
    b.plot(dd, [E_part(x, P_inf) for x in dd], "--", color="#55544f", lw=1.5, label="kernel Hann-8 entre placas infinitas + campo uniforme")
    Eb = np.array([E_part(x, P_box) for x in d[ok]])
    b.plot(d[ok], Eb, ":", color="#1baf7a", lw=1.8, label="kernel Hann-8, caja periódica ℓ = 32")
    b.plot(d[ok], ez[ok], "o", ms=5, color="#2a78d6", label="Ludwig, fuerza sobre la partícula / q")
    b.plot(d[~ok], ez[~ok], "s", ms=6, mfc="none", color="#eb6834", label="Ludwig, kernel alcanza la pared")
    b.set(xlabel="distancia d a la pared inferior (unidades de red)", ylim=(-0.05, 0.003),
          ylabel=r"$E_z$ sobre la partícula", title="(b) Campo sobre la partícula, paredes a ψ = %g y %g" % (V1, V2))
    b.legend(loc="lower right", fontsize=8.5)
    fig.tight_layout()
    fig.savefig(os.path.join(ROOT, "graficos", "C_campo_dospotenciales_kernel.png"), dpi=160)
    r = np.abs(ez[ok] - Eb) / np.abs(ez[ok])
    print("Ludwig vs kernel periodico: mediana %.1e, max %.1e (d=%.1f)" % (np.median(r), r.max(), d[ok][r.argmax()]))
    print("convergencia en ell: d=8  E(ell=200)=%.6e  E(ell=400)=%.6e" % (E_part(8.0, P_inf), E_part(8.0, modes(400.0))))
