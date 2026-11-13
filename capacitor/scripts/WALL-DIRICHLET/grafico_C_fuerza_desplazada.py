#!/usr/bin/env python3
"""Spurious force on the displaced charge (C-electrolito-hann8-LD6-desplazada[-cont], crear_desplazada.py /
crear_continuacion.py): Hann-8, dx = 1, q = +1 at x = y = 16 (on nodes), z0 = delta from the mid-plane,
delta = 0, 0.10, 0.25, 0.40, 0.50. Units: reduced field e E/kT per physical length; with q = 1 it is also the force
in kT per length.

Two quantities, z component (x, y vanish by symmetry):
  F_L      Ludwig: Esub_z of proceced_data/particle_Esub.csv at the last step. subgrid.c subgrid_update_Esub:
           Esub = sum_i W(r_i - r0) e(r_i), e = psi_electric_field (27-point stencil, wall links doubled), and the force
           is q kT/e Esub (subgrid_update_forces_electrokinetics).
           Check printed: the same sum recomputed here from the efield dump (efield/kT) at that step.
  F_kernel analytic force distributed by the kernel, the continuum limit of F_L (since 2026-11-06 computed directly,
           Fbar = -int W(z - z0) dPhibar/dz dz with the closed-form derivative HH.dpsit; the form below, integrated
           by parts, is kept as F_kernel_partes and printed as a check):
           Fbar = int W(r - r0) E_DH,z(r) d^3r,   E_DH = -grad phi_DH (linear Debye-Hueckel, grafico_C_hoja_hann.py)
           Integrating by parts in z (W vanishes at the support edges) and doing the x, y integrals in Fourier space,
             Fbar = int W'(z - z0) Phibar(z) dz,
             Phibar(z) = int int W(x) W(y) phi_DH(x, y, z) dx dy
                       = cap_closed(z) + (2 pi/A) sum_k What(kx) What(ky) psit(kx, ky, z) + fac <U> h(z)
           (int W(x) e^{i kx x} dx = What(kx), real and even; cap_closed and h depend on z only and int int W W = 1).
           z integral: Gauss-Legendre on [z0 - 4, z0] and [z0, z0 + 4].
Spurious force: F_L - F_kernel. Figure: the two forces, their difference, and (F_L - F_kernel)/F_kernel in %
(delta = 0 left out: F = 0). F_point = E_DH,z(r0) (analytic field at the centre) was removed on 2026-11-06 (user's
decision): it is not the force on the distributed charge.
Not separated: the nonlinear Poisson-Boltzmann part (Ludwig is nonlinear, phi_DH linear).
Usage: python3 grafico_C_fuerza_desplazada.py [cont | celda | placa | ref2]   (same run selection as grafico_C_desplazada.py;
       -> graficos/C_fuerza_desplazada[_cont | _celda | _placa | _ref2].png)
"""
import os, sys, glob, math
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import grafico_C_formulas as F
import grafico_C_hoja_hann as HH
import grafico_C_hoja_gauss as H
from grafico_C_desplazada import DELTAS, COL, set_charge, psi_files, run_path, SUF, R, KLAB
from grafico_C_refinamiento_zoom import KT


def Wp(u):
    """W'(u) of the Hann-8 window (1/n)(1 + cos(w u))."""
    return np.where(np.abs(u) < HH.NH / 2, -HH.OM * np.sin(HH.OM * u) / HH.NH, 0.0)

def Ubar(z, nk=HH.NK):
    """(2 pi/A) sum_k What(kx) What(ky) psit(kx, ky, z): U averaged over x, y with W(x) W(y)."""
    n, m = np.meshgrid(np.arange(-nk, nk + 1), np.arange(-nk, nk + 1), indexing="ij")
    kx = 2 * np.pi * n.ravel() / HH.L; ky = 2 * np.pi * m.ravel() / HH.L; out = np.zeros(len(z))
    for i0 in range(0, len(kx), 4000):
        sl = slice(i0, i0 + 4000)
        out += ((HH.What(kx[sl]) * HH.What(ky[sl]))[:, None] * HH.psit(kx[sl], ky[sl], z)).sum(axis=0)
    return 2 * np.pi / HH.A * out

def Phibar(z, mu, nk=HH.NK):
    fac = (HH.KAP * HH.L / 2) / math.tanh(HH.KAP * HH.L / 2)
    return H.cap_closed(z, 0, 0) + Ubar(z, nk) + fac * mu * H.h(z)

def dUbar(z, nk=HH.NK):
    """d Ubar/dz: as Ubar with the closed-form HH.dpsit."""
    n, m = np.meshgrid(np.arange(-nk, nk + 1), np.arange(-nk, nk + 1), indexing="ij")
    kx = 2 * np.pi * n.ravel() / HH.L; ky = 2 * np.pi * m.ravel() / HH.L; out = np.zeros(len(z))
    for i0 in range(0, len(kx), 4000):
        sl = slice(i0, i0 + 4000)
        out += ((HH.What(kx[sl]) * HH.What(ky[sl]))[:, None] * HH.dpsit(kx[sl], ky[sl], z)).sum(axis=0)
    return 2 * np.pi / HH.A * out

def dPhibar(z, mu, nk=HH.NK):
    """d Phibar/dz, closed form (grounded plates: phi_C' = 0)."""
    fac = (HH.KAP * HH.L / 2) / math.tanh(HH.KAP * HH.L / 2)
    return dUbar(z, nk) + fac * mu * (-HH.KAP * np.sinh(HH.KAP * z) / np.cosh(HH.KAP * HH.L / 2))

def F_kernel(z0, mu, ng=48, nk=HH.NK):
    """int W(z - z0) Ebar_z(z) dz with Ebar_z = -d Phibar/dz in closed form (2026-11-06, user's request)."""
    x, w = np.polynomial.legendre.leggauss(ng); tot = 0.0
    for a, b in ((z0 - HH.NH / 2, z0), (z0, z0 + HH.NH / 2)):
        z = 0.5 * (b - a) * x + 0.5 * (b + a)
        tot += 0.5 * (b - a) * np.sum(w * HH.W(z - z0) * -dPhibar(z, mu, nk))
    return tot

def F_kernel_partes(z0, mu, ng=48, nk=HH.NK):
    """Previous form of F_kernel (kept as a check): integrated by parts, int W'(z - z0) Phibar(z) dz."""
    x, w = np.polynomial.legendre.leggauss(ng); tot = 0.0
    for a, b in ((z0 - HH.NH / 2, z0), (z0, z0 + HH.NH / 2)):
        z = 0.5 * (b - a) * x + 0.5 * (b + a)
        tot += 0.5 * (b - a) * np.sum(w * Wp(z - z0) * Phibar(z, mu, nk))
    return tot

def F_ludwig(path):
    """(step, Esub_z * R) of the last line of particle_Esub.csv (decimal commas): Esub is per lattice spacing, the
    physical field is R times it (R = 1/dx)."""
    last = [l for l in open(os.path.join(F.ROOT, path, "proceced_data", "particle_Esub.csv")) if not l.startswith("#")][-1]
    v = last.strip().replace(",", ".").split(";")
    return int(v[0]), float(v[5]) * R

def F_ludwig_from_dump(path, step, z0lat):
    """sum_i W(r_i - r0) efield_z(r_i)/kT recomputed from the efield dump of that step (r0 = (16, 16, z0lat))."""
    n = int(round(32 * R)); nh = 8 * R                 # Hann order in nodes (physical width 8)
    e = np.loadtxt(os.path.join(F.ROOT, path, "efield-%09d.001-001" % step)).reshape(n, n, n + 2, 3)[..., 2] / KT * R
    k = np.arange(1, n + 3); w = lambda u: np.where(np.abs(u) < nh / 2, (1 + np.cos(2 * np.pi * u / nh)) / nh, 0.0)
    wx = w(16.0 * R - k[:n]); wz = w(z0lat - k)
    return np.einsum("i,j,k,ijk->", wx, wx, wz, e)

if __name__ == "__main__":
    rows = []
    print("delta  step  | F_L (Ludwig Esub)  [same from efield dump] | F_kernel (analytic)  [by parts - direct] | F_L - F_kernel | rel. (%)")
    for dl in DELTAS:
        path = run_path(dl)
        if not psi_files(path): print("(not finished) %s" % path); continue
        set_charge(dl); mu = HH.mean_U()
        step, fl = F_ludwig(path); fd = F_ludwig_from_dump(path, step, 1.5 + R * (16 + dl))
        fk = F_kernel(dl, mu); fkp = F_kernel_partes(dl, mu)
        rows.append((dl, fl, fk))
        print("%.2f  %5d  | %+.6e  [%+.6e]  | %+.6e  [%+.1e]  | %+.3e   | %s" % (dl, step, fl, fd, fk, fkp - fk, fl - fk,
              "%+.3f" % (100 * (fl - fk) / fk) if fk != 0 else "-"))
    # numerical checks of F_kernel at delta = 0.25: Gauss order, lateral modes, and the form integrated by parts
    set_charge(0.25); mu = HH.mean_U()
    f0 = F_kernel(0.25, mu); f1 = F_kernel(0.25, mu, ng=96); f2 = F_kernel(0.25, mu, nk=3 * HH.NK // 4)
    f3 = F_kernel_partes(0.25, mu)
    print("checks F_kernel(0.25): ng 48 -> 96: %.1e | NK -> 3NK/4: %.1e | analytic derivative vs int W' Phibar: %.1e  (F = %.3e)"
          % (abs(f1 - f0), abs(f2 - f0), abs(f3 - f0), f0))
    if rows:
        F.style(); r = np.array(rows)
        fig, ax = plt.subplots(3, 1, figsize=(7.5, 10.0), sharex=True, gridspec_kw=dict(hspace=0.08))
        m = r[:, 2] != 0                                  # delta = 0: F = 0, relative error undefined
        ax[0].plot(r[:, 0], r[:, 1], "o-", color="#2b2a27", label="Ludwig $\\sum_i W(r_i-r_0)E(r_i)$")
        ax[0].plot(r[:, 0], r[:, 2], "s--", color="#2a78d6", label="analytic, distributed by the kernel $\\int W E_{DH}$")
        ax[1].plot(r[:, 0], r[:, 1] - r[:, 2], "s-", color="#2a78d6")
        ax[2].plot(r[m, 0], 100 * (r[m, 1] - r[m, 2]) / r[m, 2], "s-", color="#2a78d6")
        for a in ax[1:]:
            a.axhline(0, color="#8a8880", lw=0.8)
        ax[0].set_ylabel("$F_z/q$ (reduced, $e E/k_BT$)"); ax[1].set_ylabel("$F_L - F_{kernel}$")
        ax[2].set_ylabel("$(F_L - F_{kernel})/F_{kernel}$ (%)")
        ax[2].set_xlabel("$\\delta$ (charge at $z_0 = \\delta$ from the mid-plane; 0 = mid-cell, 0.5 = node)")
        ax[0].legend(fontsize=8)
        ax[0].set_title(KLAB + ", charge displaced towards the upper plate", fontsize=10, loc="left")
        fig.savefig(os.path.join(F.ROOT, "graficos", "C_fuerza_desplazada%s.png" % SUF.replace("-", "_")), dpi=300, bbox_inches="tight")
