#!/usr/bin/env python3
"""Prueba C con electrolito (C-electrolito-hann{6,8}-LD{4,6,8}): particula subgrid q = +1 entre paredes a tierra,
sal 1:1, N iones por especie conservado, estacionario (5000 pasos). Teoria de Debye-Huckel con numero fijo de iones.

Linealizando rho_+- = rho_s [1 -+ (psi - <psi>)] (el promedio <psi> sobre el volumen fluido sale de conservar N_+-):
   nabla^2 psi - kappa^2 (psi - <psi>) = -4 pi l_B rho_p,   kappa = 1/lambda_D = sqrt(8 pi l_B rho_s)
 * modos laterales k != 0:   (k^2 + kappa^2) en lugar de k^2 en G_k(z,z') (misma G, q = sqrt(k^2+kappa^2)).
   Es la G de Zhou, Bruch, Wang (arXiv 2405.13261) Eq. (19) con kappa uniforme y G = 0 en la placa metalica.
 * modo k = 0: psi0 = phi_a + <psi> [1 - cosh(kappa(z-L/2))/cosh(kappa L/2)], phi_a con G_kappa y
   <psi> = <phi_a> (kappa l)/tanh(kappa l), l = L/2 (autoconsistencia del promedio).
Cargas inducidas: Q_inf = -(A/4 pi l_B) psi0'(0), Q_sup = +(A/4 pi l_B) psi0'(L).
Fuente y pesos de fuerza: kernel Hann-n muestreado en nodos fluidos (tratamiento T0, ver README).
Uso: python3 grafico_C_electrolito.py -> ../graficos/C_electrolito.png + tabla
"""
import os, re, glob, sys
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
L = 32.0; ELL = 32; A = float(ELL * ELL); LB4 = 10.0; BETA = 1.0e5
zn = np.arange(32) + 0.5

def wfun(r, N):
    r = np.asarray(r, float)
    return np.where(np.abs(r) < N / 2, (1 + np.cos(2 * np.pi * r / N)) / N, 0.0)

def setup(N):
    nn = np.arange(-(N // 2) + 1, N // 2); wn = wfun(nn, N)
    m = np.arange(0, ELL // 2 + 1)
    wh = np.array([(wn * np.cos(2 * np.pi * mi * nn / ELL)).sum() for mi in m])
    mult = np.where((m > 0) & (m < ELL // 2), 2.0, 1.0)
    mx, my = np.meshgrid(m, m, indexing="ij")
    k = (2 * np.pi * np.hypot(mx, my) / ELL).ravel()
    w2 = ((mult[:, None] * mult[None, :]) * (wh[:, None] * wh[None, :]) ** 2).ravel()
    return k, w2

def dG(q, z, zp):
    below = np.cosh(q * z) * np.sinh(q * (L - zp)) / np.sinh(q * L)
    above = -np.sinh(q * zp) * np.cosh(q * (L - z)) / np.sinh(q * L)
    eq = np.abs(z - zp) < 1e-9
    return np.where(eq, 0.5 * (below + above), np.where(z < zp, below, above))

def theory(d, N, kap, kmodes):
    k, w2 = kmodes
    wr = wfun(zn - d, N)                              # carga nodal (T0) = pesos de fuerza
    q = np.sqrt(k ** 2 + kap ** 2)
    iszero = k == 0
    S = dG(q[:, None, None], zn[None, :, None], zn[None, None, :])
    inner = np.einsum("kab,b->ka", S, wr)
    F = -LB4 / A * np.einsum("k,ka,a->", w2, inner, wr)          # incluye el modo k = 0 con G_kappa
    # modo k = 0: <psi> y la parte homogenea
    if kap > 0:
        intG = (1 - (np.sinh(kap * zn) + np.sinh(kap * (L - zn))) / np.sinh(kap * L)) / kap ** 2
        phia_mean = LB4 / A * (wr * intG).sum() / L
        l = L / 2
        psi_m = phia_mean * (kap * l) / np.tanh(kap * l)
        F += psi_m * kap * (wr * np.sinh(kap * (zn - l)) / np.cosh(kap * l)).sum()
        dphi0 = LB4 / A * (wr * np.sinh(kap * (L - zn)) / np.sinh(kap * L)).sum()
        dphiL = -LB4 / A * (wr * np.sinh(kap * zn) / np.sinh(kap * L)).sum()
        Qlo = -(A / LB4) * (dphi0 + psi_m * kap * np.tanh(kap * l))
        Qhi = +(A / LB4) * (dphiL - psi_m * kap * np.tanh(kap * l))
    else:
        Qlo = -(wr * (L - zn) / L).sum(); Qhi = -(wr * zn / L).sum()
    return F, Qlo, Qhi

pat = re.compile(r"Wall charge lower\s+(\S+)\s+upper\s+(\S+)")
def medido(p):
    wq = pat.findall(open(os.path.join(p, "run.log"), errors="ignore").read())[-1]
    rows = [l for l in open(os.path.join(p, "proceced_data", "particle_force.csv")) if not l.startswith("#")]
    last = rows[-1].strip().split(";")
    return float(last[5].replace(",", ".")) * BETA, float(wq[0]), float(wq[1]), int(last[0])

if __name__ == "__main__":
    plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.color": "#e4e2dc",
                         "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False})
    cols = {4: "#eb6834", 6: "#2a78d6", 8: "#1baf7a"}
    fig, ax = plt.subplots(2, 2, figsize=(12, 8.2))
    print("hann LD    d    Ez_Ludwig    Ez_teoria  teo/Lud | Qinf_Lud  Qinf_teo  dQ  | paso")
    for ci, N in enumerate((6, 8)):
        kmodes = setup(N)
        for ld in (4, 6, 8):
            kap = 1.0 / ld
            bat = os.path.join(ROOT, "C-electrolito-hann%d-LD%d" % (N, ld))
            ds, em, et, qm, qt = [], [], [], [], []
            for p in sorted(glob.glob(bat + "/z_*"), key=lambda s: float(s.rsplit("_", 1)[1])):
                d = float(p.rsplit("_", 1)[1]) - 1.5
                e, ql, qh, step = medido(p); F, Ql, Qh = theory(d, N, kap, kmodes)
                ds.append(d); em.append(e); et.append(F); qm.append(ql); qt.append(Ql)
                print("%4d %2d %5.1f %+.5e %+.5e %.4f | %+.4f %+.4f %+.1e | %d" % (N, ld, d, e, F, F / e, ql, Ql, ql - Ql, step))
            ds = np.array(ds)
            df = np.linspace(ds.min(), ds.max(), 140)               # teoria en grilla fina (no solo en los d simulados)
            tf = np.array([theory(x, N, kap, kmodes) for x in df])
            ax[0, ci].plot(ds, em, "o", ms=5, color=cols[ld]); ax[0, ci].plot(df, tf[:, 0], "-", color=cols[ld], lw=1.4, label="λD = %d" % ld)
            ax[1, ci].plot(ds, qm, "o", ms=5, color=cols[ld]); ax[1, ci].plot(df, tf[:, 1], "-", color=cols[ld], lw=1.4)
        ax[0, ci].set(title="Hann-%d: campo sobre la partícula (puntos: Ludwig)" % N, ylabel=r"$E_z$", xlabel="d"); ax[0, ci].legend()
        ax[1, ci].set(title="Hann-%d: carga inducida en la pared inferior" % N, ylabel=r"$Q_{inf}$", xlabel="d")
    fig.tight_layout(); fig.savefig(os.path.join(ROOT, "graficos", "C_electrolito.png"), dpi=160)
