#!/usr/bin/env python3
"""Error of psi and E at the lattice nodes, two fixed charges in the periodic box (2026-11-09, user's request; same metrics as
capacitor/local/WALL-DIRICHLET/scripts/grafico_C_error_capas_desplazada.py), two lattices and the observed order.

Runs (same physical system, box 48, lambda_D = 8, Hann of physical width 8, q1 = q2 = +1, separation r = 1 ... 16):
  dx = 1: SLIDES-LD8-hann08-corr_no-N10000/pos_X.00_24.00_24.00 (box 48^3, Hann-8, charges at x = 8 and X = 8 + r,
          y = z = 24), 10000 steps; stationarity psi(10000) - psi(9000).
  dx = 2: SLIDES-LD8-hann04-corr_no-dx2/r_DD (box 24^3, Hann-4, charges at x = 4.25 and 4.25 + r/2, y = z = 12.25: the
          same physical positions, crear_2part_dx2.py), 3000 steps; stationarity psi(3000) - psi(2750).
  (The first version used SLIDES-LD8-hann08-corr_no, 3000 steps: not stationary; and kappa of its run.log, 0 at r = 16.)
Units: psi dump in e psi/kT; efield dump / kT in lattice units, divided by dx -> per physical length. Distances, shells and
norms in physical units (dx = 1 lattice).

Reference: linear Debye-Hueckel, periodic box, continuous Hann sources (user's choice), in lattice units of each lattice:
  phi(r) = (4 pi l_B / V) sum_k What(kx) What(ky) What(kz) / (k^2 + kappa^2) sum_p q_p e^{i k.(r - r_p)}
  What(k) = sinc(k n/2)/(1 - (k n/(2 pi))^2) (Fourier transform of the Hann-n window), l_B = 0.7957747/dx nodes,
  kappa = 0.124975 dx (physical kappa = sqrt(4 pi l_B (rho0 + rho1)) from the mean equilibrium densities; Ludwig's own
  "Kappa calculado" is wrong when a charge sits on the central node, see the 2026-11-09 diary).
  Modes |n_i| <= M = 3N (per direction) folded mod N and summed with a 3D FFT: exact on the nodes up to the truncation at M
  (printed: change from 2N to 3N). E_DH = -grad phi (modes times -i k).
  The constant of psi is arbitrary in the periodic box: Ludwig's psi has zero node mean, so the node mean of the reference
  is removed too (E is not affected).
Error per node: e_psi = psi_L - phi_DH, e_E = |E_L - E_DH|.
1) Shells of d = distance to the NEAREST charge (minimum image), width W_SHELL (physical, w=): RMS(e)/RMS(ref) and
   max|e|/RMS(ref) per shell, and order p = log2 RMS(e, dx = 2)/RMS(e, dx = 1) in the same shell.
2) Global norms over all nodes: ||e||_2 = sqrt(sum dx^3 e^2), ||e||_inf, absolute and divided by the same norm of the
   reference, versus the separation r, and their orders p(2 -> 1).
Usage: python3 error_capas_2part.py [w=1]   (from microgel/local/scripts)
       -> ../FUERZA-2PART-MEAS/ERR_CAPAS_LD08.0_Hann-8_corr-no_dx2-1.png (shells), ..._dx2-1_orden.png (orders per shell),
          ERR_NORMAS_LD08.0_Hann-8_corr-no_dx2-1.png (global norms and orders)
"""
import os, sys, math, glob
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "FUERZA-2PART-MEAS")
KT = 1.0e-5; LB = 0.7957747; KAPPA = 0.124975; WHANN = 8.0     # physical l_B, kappa, Hann width
W_SHELL = float(([a[2:] for a in sys.argv[1:] if a.startswith("w=")] or ["1"])[0])
R_SHELLS = (1, 4, 8, 16)                                         # separations shown in the shell figures
LATS = ((2.0, "$\\Delta x = 2$ (Hann-4)", "#eb6834", "^"), (1.0, "$\\Delta x = 1$ (Hann-8)", "#2b2a27", "o"))

def runs(dx):
    """[(r, path, charge positions in lattice coordinates)] of the lattice dx, N."""
    if dx == 1.0:
        out = []
        for p in glob.glob(os.path.join(HERE, "..", "SLIDES-LD8-hann08-corr_no-N10000", "pos_*")):
            X = float(os.path.basename(p).split("_")[1]); out.append((X - 8.0, p, ((X, 24.0, 24.0), (8.0, 24.0, 24.0))))
        return sorted(out), 48
    out = []
    for p in glob.glob(os.path.join(HERE, "..", "SLIDES-LD8-hann04-corr_no-dx2", "r_*")):
        d = float(os.path.basename(p)[2:]); out.append((d, p, ((4.25 + d / 2, 12.25, 12.25), (4.25, 12.25, 12.25))))
    return sorted(out), 24

def what(k, n):
    x = np.abs(k) * n / 2
    with np.errstate(invalid="ignore", divide="ignore"):
        v = np.sinc(x / np.pi) / (1 - (x / np.pi) ** 2)
    return np.where(np.abs(x - np.pi) < 1e-9, 0.5, v)

def dh_nodes(N, kappa, lb, nh, pos, M):
    """phi_DH (zero node mean) and E_DH on the N^3 nodes, lattice units (array index j at Ludwig coordinate j + 1)."""
    L = float(N); n1 = np.arange(-M, M + 1); k1 = 2 * np.pi * n1 / L; w1 = what(k1, nh)
    C = np.zeros((4, N, N, N), dtype=complex)
    ky, kz = np.meshgrid(k1, k1, indexing="ij"); wyz = np.outer(w1, w1)
    iy, iz = np.meshgrid(n1 % N, n1 % N, indexing="ij")
    for a in range(len(n1)):                              # loop over kx
        kx = k1[a]; k2 = kx * kx + ky * ky + kz * kz
        amp = 4 * np.pi * lb / L ** 3 * w1[a] * wyz / (k2 + kappa ** 2)
        ph = sum(np.exp(1j * (kx * (1 - p[0]) + ky * (1 - p[1]) + kz * (1 - p[2]))) for p in pos)
        c = amp * ph
        for m, v in enumerate((c, -1j * kx * c, -1j * ky * c, -1j * kz * c)):
            np.add.at(C[m], (n1[a] % N, iy, iz), v)
    G = np.real(np.fft.ifftn(C, axes=(1, 2, 3))) * N ** 3
    return G[0] - G[0].mean(), np.stack(G[1:], axis=-1)

def nearest_distance(N, pos):
    """Distance (lattice units) of every node to the nearest charge, minimum image."""
    j = np.arange(1, N + 1, dtype=float); d2 = None
    for p in pos:
        dx = [(j - c) - N * np.round((j - c) / N) for c in p]
        dd = dx[0][:, None, None] ** 2 + dx[1][None, :, None] ** 2 + dx[2][None, None, :] ** 2
        d2 = dd if d2 is None else np.minimum(d2, dd)
    return np.sqrt(d2)

def shells(d, e, ref):
    """(centre, RMS(e)/RMS(ref), max|e|/RMS(ref), RMS(e)) per shell of width W_SHELL (d physical)."""
    k = np.floor(d / W_SHELL).astype(int).ravel(); e = np.abs(e).ravel(); ref = np.abs(ref).ravel(); out = []
    for s in range(k.max() + 1):
        m = k == s
        if not m.any(): continue
        rr = np.sqrt(np.mean(ref[m] ** 2)); re_ = np.sqrt(np.mean(e[m] ** 2))
        out.append(((s + 0.5) * W_SHELL, re_ / rr, e[m].max() / rr, re_))
    return np.array(out)

if __name__ == "__main__":
    res = {}
    print("shell width %g (physical); kappa %.6f (physical); reference modes |n| <= 3N" % (W_SHELL, KAPPA))
    for dx, *_ in LATS:
        rl, N = runs(dx); M = 3 * N; nh = WHANN / dx; lb = LB / dx; kap = KAPPA * dx
        for r, path, pos in rl:
            fs = sorted(glob.glob(os.path.join(path, "psi-0*")))
            psiL = np.loadtxt(fs[-1]).reshape(N, N, N); psiP = np.loadtxt(fs[-2]).reshape(N, N, N)
            EL = np.loadtxt(os.path.join(path, "efield-" + os.path.basename(fs[-1])[4:])).reshape(N, N, N, 3) / KT / dx
            phi, E = dh_nodes(N, kap, lb, nh, pos, M); E = E / dx
            ep = psiL - phi; eE = np.linalg.norm(EL - E, axis=-1); Em = np.linalg.norm(E, axis=-1)
            d = nearest_distance(N, pos) * dx; sp = shells(d, ep, phi); sE = shells(d, eE, Em); dv = dx ** 3
            g = dict(p2=math.sqrt(dv * np.sum(ep ** 2)), p2r=math.sqrt(np.sum(ep ** 2) / np.sum(phi ** 2)), pi=np.abs(ep).max(),
                     pir=np.abs(ep).max() / np.abs(phi).max(), E2=math.sqrt(dv * np.sum(eE ** 2)),
                     E2r=math.sqrt(np.sum(eE ** 2) / np.sum(Em ** 2)), Ei=eE.max(), Eir=eE.max() / Em.max())
            drift = np.abs(psiL - psiP).max(); thr = 0.01 * sp[:, 3].min()
            res[dx, r] = dict(sp=sp, sE=sE, g=g)
            msg = ""
            if r in (1, 16):                              # truncation check of the reference
                phi2, E2 = dh_nodes(N, kap, lb, nh, pos, 2 * N)
                msg = " | modes 2N -> 3N: phi %.1e, E %.1e" % (np.abs(phi - phi2).max(), np.abs(E - E2 / dx).max())
            print("dx %g r %4.1f | %s->%s drift %.1e vs %.1e %s | psi: L2 %.3e (rel %.3e) Linf %.3e (rel %.3e) | "
                  "E: L2 %.3e (rel %.3e) Linf %.3e (rel %.3e)%s"
                  % (dx, r, os.path.basename(fs[-2])[4:13], os.path.basename(fs[-1])[4:13], drift, thr,
                     "OK" if drift < thr else "NOT steady", g["p2"], g["p2r"], g["pi"], g["pir"], g["E2"], g["E2r"], g["Ei"], g["Eir"], msg))
    rs = np.array(sorted(r for (dx, r) in res if dx == 1.0))
    for r in rs:
        gc, gf = res[2.0, r]["g"], res[1.0, r]["g"]
        print("r %4.1f p(dx 2 -> 1): psi L2 %.2f Linf %.2f | E L2 %.2f Linf %.2f"
              % (r, *(math.log2(gc[k] / gf[k]) for k in ("p2", "pi", "E2", "Ei"))))
    # global norms and their orders
    fig, ax = plt.subplots(3, 2, figsize=(14.0, 11.0), sharex=True, gridspec_kw=dict(hspace=0.1, wspace=0.22))
    for j, (q, name) in enumerate((("p", "\\psi"), ("E", "E"))):
        for k, (s, lab) in enumerate((("2r", "$\\|e\\|_2 / \\|%s_{DH}\\|_2$"), ("ir", "$\\|e\\|_\\infty / \\|%s_{DH}\\|_\\infty$"))):
            for dx, ll, c, mk in LATS:
                ax[k, j].plot(rs, [100 * res[dx, r]["g"][q + s] for r in rs], mk + "-", color=c, lw=1.0, ms=4, label=ll)
            ax[k, j].set_ylabel((lab % name) + " (%)"); ax[k, j].legend(fontsize=8)
        for s, lab, mk in (("2", "$L_2$", "o-"), ("i", "$L_\\infty$", "s--")):
            ax[2, j].plot(rs, [math.log2(res[2.0, r]["g"][q + s] / res[1.0, r]["g"][q + s]) for r in rs], mk, color="#2a78d6",
                          lw=1.0, ms=4, label="$p(\\Delta x = 2 \\to 1)$, " + lab)
        for p in (1, 2):
            ax[2, j].axhline(p, color="#8a8880", lw=0.8, ls=":")
        ax[2, j].set_ylabel("observed order $p$"); ax[2, j].legend(fontsize=8)
        ax[0, j].set_title("$%s$: global norms over all nodes, Hann width 8, $\\lambda_D = 8$" % name, fontsize=10, loc="left")
        ax[2, j].set_xlabel("separation $r$ between the charges")
    fig.savefig(os.path.join(OUT, "ERR_NORMAS_LD08.0_Hann-8_corr-no_dx2-1.png"), dpi=300, bbox_inches="tight")
    # shells: error and order per shell
    for kind in ("err", "orden"):
        fig, ax = plt.subplots(len(R_SHELLS), 2, figsize=(14.0, 4.0 * len(R_SHELLS)), sharex=True,
                               gridspec_kw=dict(hspace=0.12, wspace=0.22))
        for i, r in enumerate(R_SHELLS):
            for j, (key, name) in enumerate((("sp", "\\psi"), ("sE", "E"))):
                a = ax[i, j]
                if kind == "err":
                    for dx, ll, c, mk in LATS:
                        s = res[dx, float(r)][key]
                        a.semilogy(s[:, 0], s[:, 1], mk + "-", color=c, lw=1.0, ms=4, label="RMS, " + ll)
                        a.semilogy(s[:, 0], s[:, 2], mk + ":", color=c, lw=0.8, ms=3, alpha=0.6, label="max, " + ll)
                    a.set_ylabel("$|e_{%s}|$ / RMS($%s_{DH}$) in the shell" % (name, name))
                else:
                    sc, sf = res[2.0, float(r)][key], res[1.0, float(r)][key]
                    _, ic, jf = np.intersect1d(np.round(sc[:, 0], 9), np.round(sf[:, 0], 9), return_indices=True)  # by shell centre
                    a.plot(sc[ic, 0], np.log2(sc[ic, 3] / sf[jf, 3]), "^-", color="#eb6834", lw=1.0, ms=4,
                           label="$p(\\Delta x = 2 \\to 1)$")
                    for p in (1, 2):
                        a.axhline(p, color="#8a8880", lw=0.8, ls=":")
                    a.set_ylabel("order $p$ of RMS($e_{%s}$) per shell" % name)
                a.axvline(WHANN / 2, color="#8a8880", lw=0.8, ls="--", zorder=0)
                a.set_title("$%s$, separation $r = %d$ (grey: kernel half width 4)" % (name, r), fontsize=9, loc="left")
                if i == 0: a.legend(fontsize=7, ncol=2 if kind == "err" else 1)
                if i == len(R_SHELLS) - 1: a.set_xlabel("distance $d$ to the nearest charge (shell centre, width %g)" % W_SHELL)
        fig.savefig(os.path.join(OUT, "ERR_CAPAS_LD08.0_Hann-8_corr-no_dx2-1%s.png" % ("" if kind == "err" else "_orden")),
                    dpi=300, bbox_inches="tight")
