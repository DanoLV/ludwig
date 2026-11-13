#!/usr/bin/env python3
"""Free energy of the two-charge system vs time and chemical-potential residual (2026-11-10, user's request: steps 0-2).

Runs: separation r = 4 (physical) in the three lattices of the same physical system (box 48, lambda_D = 8, Hann of
physical width 8, cos^2, q1 = q2 = +1):
  dx = 2  : SLIDES-LD8-hann04-corr_no-dx2/r_04                 (24^3, Hann-4,  charges at x = 4.25, 6.25;  y = z = 12.25)
  dx = 1  : SLIDES-LD8-hann08-corr_no-N10000/pos_12.00_24.00_24.00 (48^3, Hann-8, x = 8, 12;  y = z = 24)
  dx = 1/2: SLIDES-LD8-hann16-corr_no-dxhalf/r_04              (96^3, Hann-16, x = 15.5, 23.5;  y = z = 47.5)
Every psi/qsi dump (psi in e psi/kT; qsi = rho_+, rho_- of the fluid only, per node; valencies Z = +1, -1).

Step 0 (checks): N_a = sum_n rho_a,n per dump (must be constant: NP is in flux form); ion charge Q_ion = N_+ - N_-;
  fixed charge Q_fix = sum_n rho_fix,n (2 if the lattice Hann weights sum to 1); total charge Q_ion + Q_fix.
Step 1: free energy in kT, the functional of fe_electro.c (Rotenberg et al. 2010), psi already in units of kT/e:
  F/kT = sum_n [ sum_a rho_a (ln rho_a - 1) + (1/2) rho_e psi ],  rho_e = rho_+ - rho_- + rho_fix,
  rho_fix,n = sum_i q_i w(x - X_i) w(y - Y_i) w(z - Z_i),  w = d_hann of subgrid.c (cos^2: (1 + cos(2 pi r/n))/n, |r| < n/2),
  node at array index j = Ludwig coordinate j + 1.
  Excess (reference stated, choice of this script): F_exc = F - sum_a N_a (ln(N_a/V) - 1), i.e. minus the entropy of the
  same ions spread uniformly; only a constant per run, it does not change the time differences.
  Printed: Delta F between consecutive dumps (kT) and Delta F / |F_exc|. Not included: LB kinetic energy.
Step 2: mu_a,n / kT = ln rho_a,n + Z_a psi_n (as calc_mui.py); residual Delta mu_a = max_n - min_n (kT). At the minimum of F
  with N_a fixed, mu_a is uniform (Delta mu_a = 0).
  Caveat: psi and qsi are written at the same step but psi comes from the Poisson solve with the densities of the start of
  the step; the mismatch is of the order of one NP step (zero at steady state).
Usage: python3 energia_2part.py   (from microgel/local/scripts) -> table printed and ../FUERZA-2PART-MEAS/ENERGIA_r04.csv
"""
import os, glob
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "FUERZA-2PART-MEAS")
RUNS = ((2.0, "SLIDES-LD8-hann04-corr_no-dx2/r_04", 24, 4.0, ((6.25, 12.25, 12.25), (4.25, 12.25, 12.25))),
        (1.0, "SLIDES-LD8-hann08-corr_no-N10000/pos_12.00_24.00_24.00", 48, 8.0, ((12.0, 24.0, 24.0), (8.0, 24.0, 24.0))),
        (0.5, "SLIDES-LD8-hann16-corr_no-dxhalf/r_04", 96, 16.0, ((23.5, 47.5, 47.5), (15.5, 47.5, 47.5))))
Q = (1.0, 1.0); Z = (1.0, -1.0)

def w_hann(r, n):
    return np.where(np.abs(r) < 0.5 * n, (1.0 + np.cos(2 * np.pi * r / n)) / n, 0.0)

def rho_fix(N, n, pos):
    j = np.arange(1, N + 1, dtype=float); out = np.zeros((N, N, N))
    for q, p in zip(Q, pos):
        w = [w_hann((j - c) - N * np.round((j - c) / N), n) for c in p]     # minimum image
        out += q * w[0][:, None, None] * w[1][None, :, None] * w[2][None, None, :]
    return out

if __name__ == "__main__":
    rows = []
    for dx, rel, N, n, pos in RUNS:
        path = os.path.join(HERE, "..", rel); rf = rho_fix(N, n, pos); V = float(N ** 3)
        print("\n== dx = %g: %s  (Q_fix = %.15f)" % (dx, rel, rf.sum()))
        print(" step  |  t(dx=1) |      N_+          N_-       | Q_ion + Q_fix |   F/kT              F_exc/kT    "
              "|  dF (kT)     dF/|F_exc| | dmu_+ (kT)  dmu_- (kT)")
        prev = None
        for f in sorted(glob.glob(os.path.join(path, "psi-0*"))):
            st = int(os.path.basename(f)[4:13])
            psi = np.loadtxt(f).reshape(N, N, N)
            q = np.loadtxt(os.path.join(path, "qsi-%09d.001-001" % st)); rp = q[:, 0].reshape(N, N, N); rm = q[:, 1].reshape(N, N, N)
            Np, Nm = rp.sum(), rm.sum()
            S = np.sum(rp * (np.log(rp) - 1)) + np.sum(rm * (np.log(rm) - 1))
            E = 0.5 * np.sum((rp - rm + rf) * psi)
            F = S + E; Fexc = F - sum(Na * (np.log(Na / V) - 1) for Na in (Np, Nm))
            mup = np.log(rp) + Z[0] * psi; mum = np.log(rm) + Z[1] * psi
            dmp, dmm = mup.max() - mup.min(), mum.max() - mum.min()
            dF = np.nan if prev is None else F - prev
            print("%6d | %8.0f | %.10e %.10e | %+.3e | %+.10e %+.8e | %+.3e %+.2e | %.3e  %.3e"
                  % (st, st * dx * dx, Np, Nm, Np - Nm + rf.sum(), F, Fexc, dF, dF / abs(Fexc), dmp, dmm))
            rows.append((dx, st, st * dx * dx, Np, Nm, Np - Nm + rf.sum(), S, E, F, Fexc, dF, dmp, dmm))
            prev = F
    os.makedirs(OUT, exist_ok=True)
    np.savetxt(os.path.join(OUT, "ENERGIA_r04.csv"), np.array(rows), delimiter=",", fmt="%.15e",
               header="dx,step,t_dx1,N_plus,N_minus,Q_total,S_entropy_kT,E_elec_kT,F_kT,F_exc_kT,dF_kT,dmu_plus_kT,dmu_minus_kT")
