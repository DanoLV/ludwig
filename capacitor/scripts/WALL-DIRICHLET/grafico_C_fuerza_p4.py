#!/usr/bin/env python3
"""Force error of the displaced charge with the Hann cos^4 kernel vs the original Hann cos^2, delta = 1, dx = 2 and 1
(2026-11-09, user's request).

Runs (crear_desplazada_p4.py; the cos^4 runs are copies of the cos^2 ones with only the kernel changed):
  cos^2, physical width 8:  dx = 2 C-electrolito-hann4-LD3-ref0.5-desplazada/z_10.00 (Hann-4)
                            dx = 1 C-electrolito-hann8-LD6-desplazada/z_18.50          (Hann-8)
  cos^4, physical width 10: dx = 2 C-electrolito-hann5p4-LD3-ref0.5-desplazada/z_10.00 (n = 5)
                            dx = 1 C-electrolito-hann10p4-LD6-desplazada/z_18.50        (n = 10)
Each kernel against ITS OWN continuum target (same definitions as grafico_C_fuerza_desplazada.py):
  F_L      = Esub_z / dx (last line of proceced_data/particle_Esub.csv, physical units)
  F_kernel = int W(r - r0) E_DH,z(r) d^3r, E_DH the linear Debye-Hueckel field of the charge spread with the same W
             (grafico_C_hoja_hann_pot.py: grafico_C_hoja_hann.py generalised to cos^p)
  e = (F_L - F_kernel)/F_kernel,   observed order p = log2 |e(dx = 2)/e(dx = 1)| (two lattices: one number; it assumes
  e ~ C dx^p already holds at dx = 1, it does not test it).
Checks printed:
  - kernel functions vs direct formulas / quadrature (grafico_C_hoja_hann_pot.kernel_checks), psit and <U> (HH.checks);
  - F_kernel against its by-parts form F_kernel_partes;
  - cos^2: F_kernel with the generalised functions equals the original module's value;
  - F_L recomputed from the efield dump with the lattice weights of the kernel (sum_i W W W E_z), against Esub;
  - steady state: max|psi(N) - psi(N - freq)| and max|E(N) - E(N - freq)| on the whole lattice against 1 % of the smallest
    RMS of (Ludwig - DH) on the z axis through the charge in bands of width 2 (criterion of grafico_C_convergencia_desplazada).
Not separated: nonlinear Poisson-Boltzmann part (Ludwig nonlinear, phi_DH linear).
Usage: python3 grafico_C_fuerza_p4.py   (text output only)
"""
import os, math
import numpy as np
import grafico_C_formulas as F
import grafico_C_hoja_hann as HH
import grafico_C_fuerza_desplazada as FD
import grafico_C_hoja_hann_pot as HP
from grafico_C_desplazada import set_charge, ez_an, psi_files
from grafico_C_refinamiento import load, line
from grafico_C_refinamiento_anidado import band_rms
from grafico_C_refinamiento_zoom import ez_file, KT

DL = 1.0
KERNELS = (  # (label, physical width, p, {r = 1/dx: run})
    ("cos^2", 8, 2, {0.5: "C-electrolito-hann4-LD3-ref0.5-desplazada/z_10.00", 1: "C-electrolito-hann8-LD6-desplazada/z_18.50"}),
    ("cos^4", 10, 4, {0.5: "C-electrolito-hann5p4-LD3-ref0.5-desplazada/z_10.00", 1: "C-electrolito-hann10p4-LD6-desplazada/z_18.50"}))

def f_ludwig(path, r):
    last = [l for l in open(os.path.join(F.ROOT, path, "proceced_data", "particle_Esub.csv")) if not l.startswith("#")][-1]
    v = last.strip().replace(",", ".").split(";")
    return int(v[0]), float(v[5]) * r

def f_from_dump(path, step, r, nphys, p):
    """sum_i w(x_i) w(y_i) w(z_i - z0) E_z(r_i) from the efield dump, w = Ludwig's cos^p window of order nphys r (nodes)."""
    n = int(round(32 * r)); nh = nphys * r; cp = np.prod([(k - 1) / k for k in range(2, p + 1, 2)])
    e = np.loadtxt(os.path.join(F.ROOT, path, "efield-%09d.001-001" % step)).reshape(n, n, n + 2, 3)[..., 2] / KT * r
    w = lambda u: np.where(np.abs(u) < nh / 2, np.cos(np.pi * u / nh) ** p / (nh * cp), 0.0)
    k = np.arange(1, n + 3); z0lat = 1.5 + r * (16 + DL)
    wx = w(16.0 * r - k[:n]); wz = w(z0lat - k)
    return np.einsum("i,j,k,ijk->", wx, wx, wz, e)

def steady(path, r, mu):
    fs = psi_files(path); assert fs, "not finished: " + path
    n = int(round(32 * r))
    z, ps, _ = line(load(fs[-1], r), r, 0); ez = ez_file(path, fs[-1], r)
    _, Rp = band_rms(z, ps - HH.phi(0.0, z, 0, 0, mu), r, 2.0, 2.0); _, Re = band_rms(z, ez - ez_an(z, mu), r, 2.0, 2.0)
    fe = lambda f: np.loadtxt(os.path.join(F.ROOT, path, "efield-" + os.path.basename(f)[4:])).reshape(n, n, n + 2, 3) / KT * r
    return (os.path.basename(fs[-2])[4:13], os.path.basename(fs[-1])[4:13], np.abs(load(fs[-1], r) - load(fs[-2], r)).max(),
            0.01 * Rp.min(), np.abs(fe(fs[-1]) - fe(fs[-2])).max(), 0.01 * Re.min())

if __name__ == "__main__":
    set_charge(DL); fk_orig = FD.F_kernel(DL, HH.mean_U())          # original cos^2 module, before any replacement
    res = {}
    for lab, nphys, p, runs in KERNELS:
        HP.set_kernel(nphys, p); set_charge(DL); mu = HH.mean_U()
        print("==== %s, physical width %d, delta = %g" % (lab, nphys, DL))
        HP.kernel_checks(p); HH.checks(mu)
        fk = FD.F_kernel(DL, mu); fkp = FD.F_kernel_partes(DL, mu)
        print("F_kernel %.10e | by parts %.10e (rel diff %.1e)" % (fk, fkp, abs(fk - fkp) / fk))
        if p == 2:
            print("F_kernel generalised vs original module: %.10e vs %.10e (rel diff %.1e)" % (fk, fk_orig, abs(fk - fk_orig) / fk_orig))
        for r in (0.5, 1):
            st, fl = f_ludwig(runs[r], r); fd = f_from_dump(runs[r], st, r, nphys, p)
            s = steady(runs[r], r, mu); e = (fl - fk) / fk; res[lab, r] = e
            print("dx = %g: F_L %.10e (step %d; from the efield dump %.10e, rel diff %.1e) | (F_L - F_kernel)/F_kernel = %+.4f %%"
                  % (1 / r, fl, st, fd, abs(fd - fl) / abs(fl), 100 * e))
            print("        steady %s->%s: max|dpsi| %.1e vs %.1e %s, max|dE| %.1e vs %.1e %s"
                  % (s[0], s[1], s[2], s[3], "OK" if s[2] < s[3] else "NOT", s[4], s[5], "OK" if s[4] < s[5] else "NOT"))
        print("p(dx 2 -> 1) = %.3f" % math.log2(abs(res[lab, 0.5] / res[lab, 1])))
