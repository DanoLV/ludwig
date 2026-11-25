#!/usr/bin/env python3
"""Charge at the centre of the channel, displaced towards the upper plate by fractions of a cell (Hann-8, dx = 1,
lambda_D = 6, plates grounded; runs C-electrolito-hann8-LD6-desplazada built by crear_desplazada.py):
  z0 = delta (physical, from the mid-plane),  delta = 0, 0.10, 0.25, 0.40, 0.50, 0.60, 0.75, 0.90
delta = 0 is the mid-plane (half way between two nodes) and delta = 0.5 is a node; 0.10, 0.25, 0.40 are not lattice
symmetry positions in z.

On the z axis through the charge (x = y = 16, on nodes) and for each delta:
  psi residual   d_psi(z) = psi_Ludwig(z) - phi_DH(z)
  field residual d_E(z)   = E_z,Ludwig(z) - E_z,DH(z),   E_z,DH = -d phi_DH/dz
phi_DH: analytic linear Debye-Hueckel solution with the continuous Hann-8 charge centred at z0 = delta
(grafico_C_hoja_hann.py, its globals Z0, ZA, ZB set to the displaced charge). E_z,DH: closed-form derivative of phi_DH
(grafico_C_hoja_hann.dphi, since 2026-11-06; before, central difference with h = 1e-4, agreement ~5e-11). E_z,Ludwig: Ludwig's own efield dump (27-point stencil), read with
grafico_C_refinamiento_zoom.ez_file. The residuals contain both the lattice error and the nonlinear Poisson-Boltzmann
part (Ludwig is nonlinear, phi_DH is linear); this script does not separate them.
Band RMS: bands of width 2 with ends at z + 16 = 2j (whole cells), fixed in the lab frame (the charge moves inside
them). Pointwise errors (2026-11-05): absolute |d_i| and relative |d_i / analytic_i| at each node (log scale; the
relative one diverges where the analytic value crosses zero: E_z at the charge, psi at the grounded plates).
Relative band error (user's choice, 2026-11-05): e_b = RMS_b(d) / RMS_b(analytic), the band RMS of the residual
over the band RMS of the analytic phi_DH (or E_z,DH) in the same band.
Steady state (user's criterion, 2026-11-05): with the last two dumps N and N - 1000 (whole 3D lattice)
  D_psi = max|psi(N) - psi(N-1000)| < 1 % of the smallest band RMS of d_psi,   D_E = max|E(N) - E(N-1000)| (all three
  components, same units as d_E) < 1 % of the smallest band RMS of d_E.
Usage: python3 grafico_C_desplazada.py [cont | celda]
  (none)  runs C-electrolito-hann8-LD6-desplazada (step 5000)          -> graficos/C_desplazada.png
  cont    runs C-electrolito-hann8-LD6-desplazada-cont (restart from 5000 to 10000, crear_continuacion.py)
                                                                        -> graficos/C_desplazada_cont.png
  celda   the whole cell delta = 0 ... 0.9: for each delta the -cont run if it exists (delta <= 0.5, steady state
          checked at 9000 -> 10000), otherwise the base run (delta > 0.5, psi every 1000 steps, checked at
          4000 -> 5000)                                                 -> graficos/C_desplazada_celda.png
  placa   approach to the upper plate, delta = 1.0, 1.6, 2.2, ..., 12.0 (base runs, psi every 1000 steps; delta = 12 puts
          the support edge on the plate)                                -> graficos/C_desplazada_placa.png
  ref2    lattice dx = 1/2 (Hann-16, C-electrolito-hann16-LD12-ref2-desplazada, crear_desplazada_ref.py),
          delta = 1.0, 4.7, 8.3, 12.0; psi every 4000 steps (same physical time as 1000 steps with dx = 1),
          steady state checked at 16000 -> 20000                       -> graficos/C_desplazada_ref2.png
"""
import os, sys, glob
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import grafico_C_formulas as F
import grafico_C_hoja_hann as HH
from grafico_C_refinamiento import load, line
from grafico_C_refinamiento_anidado import band_rms
from grafico_C_refinamiento_zoom import ez_file, KT

MODE = ([a for a in sys.argv[1:] if a in ("cont", "celda", "placa", "ref2")] or [""])[0]
R = 2.0 if MODE == "ref2" else 1.0                   # nodes per physical length unit (dx = 1/R)
SUF = "-" + MODE if MODE else ""
KLAB = "Hann-16, $\\Delta x = 1/2$" if MODE == "ref2" else "Hann-8, $\\Delta x = 1$"   # figure titles
BAT = "C-electrolito-hann8-LD6-desplazada"
DELTAS = ((1.0, 4.7, 8.3, 12.0) if MODE == "ref2" else
          (1.0, 1.6, 2.2, 2.8, 3.4, 4.0, 4.7, 5.3, 5.9, 6.5, 7.1, 7.7, 8.3, 9.0, 9.6, 10.2, 10.8, 11.4, 12.0)
          if MODE == "placa" else
          (0.0, 0.10, 0.25, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90))

def run_path(dl):
    """Run directory of position delta for the chosen MODE."""
    z = "z_%05.2f" % (1.5 + R * (16 + dl))
    if MODE == "ref2":
        return os.path.join("C-electrolito-hann16-LD12-ref2-desplazada", z)
    if MODE == "cont" or (MODE == "celda" and os.path.isdir(os.path.join(F.ROOT, BAT + "-cont", z))):
        return os.path.join(BAT + "-cont", z)
    return os.path.join(BAT, z)
COL = {d: matplotlib.colormaps["viridis"](0.9 * (d - min(DELTAS)) / (max(DELTAS) - min(DELTAS))) for d in DELTAS}   # ordered: dark (0) -> yellow (0.9)

def psi_files(path):
    """psi dumps of a run, or [] if the last step N_start + N_cycles of its input has not been written yet."""
    fs = sorted(glob.glob(os.path.join(F.ROOT, path, "psi-0*")))
    key = {l.split()[0]: int(l.split()[1]) for l in open(os.path.join(F.ROOT, path, "input")) if l.startswith(("N_start", "N_cycles"))}
    return fs if fs and int(os.path.basename(fs[-1])[4:13]) == key["N_start"] + key["N_cycles"] else []

def set_charge(z0):
    HH.Z0 = z0; HH.ZA = max(-HH.L / 2, z0 - HH.NH / 2); HH.ZB = min(HH.L / 2, z0 + HH.NH / 2)

def ez_an(z, mu):
    """E_z,DH = -d phi_DH/dz on the axis, closed-form derivative (HH.dphi, 2026-11-06; before: ez_fd)."""
    return -HH.dphi(0.0, z, 0, 0, mu)

def ez_fd(z, mu, h=1e-4):
    """Previous version of ez_an: central difference of phi_DH (agrees with ez_an to ~5e-11)."""
    return -(HH.phi(0.0, z + h, 0, 0, mu) - HH.phi(0.0, z - h, 0, 0, mu)) / (2 * h)

def efield_3d(path, fname_psi):
    """Ludwig's whole efield dump (3 components) of the same step as fname_psi, in the units of ez_file."""
    fe = os.path.join(F.ROOT, path, "efield-" + os.path.basename(fname_psi)[4:])
    n = int(round(32 * R))
    return np.loadtxt(fe).reshape(n, n, n + 2, 3) / KT * R

if __name__ == "__main__":
    F.style()
    # 2026-11-16 (user's request): in the ref2 figure the two band-RMS rows are not drawn (still computed and printed,
    # they are used by the steady-state check)
    BANDS = MODE != "ref2"; NR = 6 if BANDS else 4
    fig, ax = plt.subplots(NR, 2, figsize=(11.0, 17.0 * NR / 6 + (0 if BANDS else 0.6)), sharex=True,
                           gridspec_kw=dict(height_ratios=[1.3] + [1] * (NR - 1), hspace=0.08, wspace=0.3))
    zz = np.linspace(-16, 16, 1281); tab = {}; steady = {}
    for dl in DELTAS:
        path = run_path(dl)
        fs = psi_files(path)
        if not fs: print("(no data yet) %s" % path); continue
        set_charge(dl); mu = HH.mean_U()
        a = load(fs[-1], R); z, ps, _ = line(a, R, 0)
        drift = np.abs(a - load(fs[-2], R)).max() if len(fs) > 2 else float("nan")
        drift_e = np.abs(efield_3d(path, fs[-1]) - efield_3d(path, fs[-2])).max() if len(fs) > 2 else float("nan")
        ez = ez_file(path, fs[-1], R)
        pa = HH.phi(0.0, z, 0, 0, mu); ea = ez_an(z, mu)
        dp = ps - pa; de = ez - ea
        zc, Rp = band_rms(z, dp, R, 2.0, 2.0); _, Re = band_rms(z, de, R, 2.0, 2.0)
        _, Ap = band_rms(z, pa, R, 2.0, 2.0); _, Ae = band_rms(z, ea, R, 2.0, 2.0)
        tab[dl] = (Rp, Re, Rp / Ap, Re / Ae)
        steady[dl] = (drift < 0.01 * Rp.min(), drift_e < 0.01 * Re.min())
        i0 = np.argmin(np.abs(z - dl))
        print("delta = %.2f: last dumps %s, %s" % (dl, os.path.basename(fs[-2])[4:13], os.path.basename(fs[-1])[4:13]))
        print("  D_psi = %.1e vs 1%% of min band RMS d_psi = %.1e -> %s | D_E = %.1e vs 1%% of min band RMS d_E = %.1e -> %s"
              % (drift, 0.01 * Rp.min(), "OK" if steady[dl][0] else "NOT steady",
                 drift_e, 0.01 * Re.min(), "OK" if steady[dl][1] else "NOT steady"))
        print("  max|d_psi| %.2e | max|d_E| %.2e | d_psi, d_E at the nearest node z = %+.1f: %+.2e, %+.2e"
              % (np.abs(dp).max(), np.abs(de).max(), z[i0], dp[i0], de[i0]))
        kw = dict(ms=4.0, color=COL[dl], mec="white", mew=0.4)
        lab = "$\\delta = %.2f$" % dl
        ax[0, 0].plot(zz, HH.phi(0.0, zz, 0, 0, mu), "-", lw=0.9, color=COL[dl], alpha=0.6)
        ax[0, 0].plot(z, ps, "o", label=lab, **kw)
        ax[0, 1].plot(zz, ez_an(zz, mu), "-", lw=0.9, color=COL[dl], alpha=0.6)
        ax[0, 1].plot(z, ez, "o", **kw)
        ax[1, 0].plot(z, dp, "o-", lw=0.9, **kw); ax[1, 1].plot(z, de, "o-", lw=0.9, **kw)
        ax[2, 0].semilogy(z, np.abs(dp), "o-", lw=0.9, **kw); ax[2, 1].semilogy(z, np.abs(de), "o-", lw=0.9, **kw)
        ax[3, 0].semilogy(z, 100 * np.abs(dp / pa), "o-", lw=0.9, **kw); ax[3, 1].semilogy(z, 100 * np.abs(de / ea), "o-", lw=0.9, **kw)
        if BANDS:
            ax[4, 0].semilogy(zc, Rp, "o-", lw=0.9, **kw); ax[4, 1].semilogy(zc, Re, "o-", lw=0.9, **kw)
            ax[5, 0].semilogy(zc, 100 * Rp / Ap, "o-", lw=0.9, **kw); ax[5, 1].semilogy(zc, 100 * Re / Ae, "o-", lw=0.9, **kw)
    if tab:
        print("band centre    " + " ".join("%8.1f" % c for c in zc))
        for k, name in ((0, "RMS d_psi   "), (1, "RMS d_E     "), (2, "rel. psi (%)"), (3, "rel. E   (%)")):
            for dl, R in tab.items():
                print("%s %.2f " % (name, dl) + " ".join("%8.1e" % (v * (100 if k >= 2 else 1)) for v in R[k]))
    for a in ax[1]:
        a.axhline(0, color="#8a8880", lw=0.8)
    for a in ax.ravel():
        a.axvspan(min(DELTAS) - 4, max(DELTAS) + 4, color="#f1efe8", zorder=0); a.set_xlim(-16.5, 16.5)
    for a in ax[NR - 1]:
        a.set_xlabel("z (physical units, from the mid-plane)")
    if BANDS:
        ax[5, 0].set_ylabel("RMS$_b(d_\\psi)$ / RMS$_b(\\phi_{DH})$ (%)"); ax[5, 1].set_ylabel("RMS$_b(d_E)$ / RMS$_b(E_{z,DH})$ (%)")
    ax[2, 0].set_ylabel("$|d_\\psi|$ at each node"); ax[2, 1].set_ylabel("$|d_E|$ at each node")
    ax[3, 0].set_ylabel("$|d_\\psi / \\phi_{DH}|$ at each node (%)"); ax[3, 1].set_ylabel("$|d_E / E_{z,DH}|$ at each node (%)")
    ax[0, 0].set_ylabel(r"$\psi$ (Ludwig: points; DH: lines)"); ax[0, 1].set_ylabel(r"$E_z$ (Ludwig efield: points; DH: lines)")
    ax[1, 0].set_ylabel(r"$d_\psi = \psi_{Ludwig} - \phi_{DH}$"); ax[1, 1].set_ylabel(r"$d_E = E_{z,Ludwig} - E_{z,DH}$")
    if BANDS:
        ax[4, 0].set_ylabel(r"band RMS of $d_\psi$ (width 2)"); ax[4, 1].set_ylabel(r"band RMS of $d_E$ (width 2)")
    ax[0, 0].legend(fontsize=8, loc="upper left")
    ax[0, 0].set_title(KLAB + ", charge at $z_0 = \\delta$ (shaded: union of the kernel supports)", fontsize=10, loc="left")
    fig.savefig(os.path.join(F.ROOT, "graficos", "C_desplazada%s.png" % SUF.replace("-", "_")), dpi=300, bbox_inches="tight")
