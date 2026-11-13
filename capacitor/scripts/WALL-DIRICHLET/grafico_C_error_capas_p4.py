#!/usr/bin/env python3
"""Error of psi and E at the lattice nodes near the kernel edge: Hann cos^4 (physical width 10) vs Hann cos^2 (width 8),
delta = 1, dx = 2 and 1 (2026-11-09, user's request: convergence at the edge of the kernel).

Same metric as grafico_C_error_capas_desplazada.py (whose functions dh_on_nodes, distances, shells are used):
  e_psi(r_i) = psi_L(r_i) - phi_DH(r_i),  e_E(r_i) = |E_L(r_i) - E_DH(r_i)| at every fluid node, last dump, no interpolation;
  phi_DH: linear Debye-Hueckel solution of the charge spread with the SAME kernel as the run (cos^4: grafico_C_hoja_hann_pot);
  shells of d = |r_i - r0| (minimum image in x, y) of width W_SHELL (physical, parameter w=):
  RMS(e)/RMS(ref) per shell and the order p = log2 RMS(e, dx = 2)/RMS(e, dx = 1) in the same shell.
Kernel edge: the support is the cube |x|, |y|, |z - z0| < n/2 (tensor product), so in d it starts at n/2 (faces) and ends
at n sqrt(3)/2 (corners): 4 ... 6.93 for cos^2, 5 ... 8.66 for cos^4 (lines in the figure).
Smoothness at the edge: cos^2 window C^1 (W'' jumps), cos^4 C^3 (W'''' jumps) (subgrid.c d_hann comment).
Global norms (all fluid nodes, ||e||_2 = sqrt(sum dx^3 e^2), ||e||_inf) printed with their orders.
Runs: crear_desplazada_p4.py (cos^4) and the cos^2 runs of the same delta and lattices.
Usage: python3 grafico_C_error_capas_p4.py [w=1]   -> graficos/C_error_capas_p4.png
"""
import os, sys, math
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import grafico_C_formulas as F
import grafico_C_hoja_hann as HH
import grafico_C_hoja_hann_pot as HP
import grafico_C_error_capas_desplazada as EC
from grafico_C_refinamiento import load
from grafico_C_refinamiento_zoom import KT
from grafico_C_desplazada import set_charge, psi_files
from grafico_C_fuerza_p4 import KERNELS, DL
from grafico_C_convergencia_desplazada import DPI

EC.W_SHELL = float(([a[2:] for a in sys.argv[1:] if a.startswith("w=")] or ["1"])[0])
COLS = {"cos^2": "#2b2a27", "cos^4": "#2a78d6"}

if __name__ == "__main__":
    res = {}
    print("shell width %g (physical), delta = %g" % (EC.W_SHELL, DL))
    for lab, nphys, p, runs in KERNELS:
        HP.set_kernel(nphys, p); set_charge(DL); mu = HH.mean_U()
        for r in (0.5, 1):
            path = runs[r]; fs = psi_files(path); assert fs, "not finished: " + path
            n = int(round(32 * r)); ic = int(round(16 * r)) - 1
            psiL = load(fs[-1], r)[:, :, 1:n + 1]
            EL = np.loadtxt(os.path.join(F.ROOT, path, "efield-" + os.path.basename(fs[-1])[4:])).reshape(n, n, n + 2, 3)[:, :, 1:n + 1] / KT * r
            phi, E, z = EC.dh_on_nodes(r, mu)
            if r == 1:                                    # check of the FFT evaluation on the axis
                print("%s check on the axis: |phi_FFT - HH.phi| %.1e, |dphi_FFT - HH.dphi| %.1e" % (lab,
                      np.abs(phi[ic, ic] - HH.phi(0.0, z, 0, 0, mu)).max(), np.abs(-E[ic, ic, :, 2] - HH.dphi(0.0, z, 0, 0, mu)).max()))
            ep = psiL - phi; eE = np.linalg.norm(EL - E, axis=-1); Em = np.linalg.norm(E, axis=-1)
            d = EC.distances(r, DL); dv = (1 / r) ** 3
            g = dict(p2=math.sqrt(dv * np.sum(ep ** 2)), pi=np.abs(ep).max(), E2=math.sqrt(dv * np.sum(eE ** 2)), Ei=eE.max(),
                     p2r=math.sqrt(np.sum(ep ** 2) / np.sum(phi ** 2)), E2r=math.sqrt(np.sum(eE ** 2) / np.sum(Em ** 2)))
            res[lab, r] = dict(sp=EC.shells(d, ep, phi), sE=EC.shells(d, eE, Em), g=g)
            print("%s dx = %g | psi: L2 rel %.3e, Linf %.3e | E: L2 rel %.3e, Linf %.3e" % (lab, 1 / r, g["p2r"], g["pi"], g["E2r"], g["Ei"]))
        gc, gf = res[lab, 0.5]["g"], res[lab, 1]["g"]
        print("%s p(dx 2 -> 1) global: psi L2 %.2f Linf %.2f | E L2 %.2f Linf %.2f"
              % (lab, *(math.log2(gc[k] / gf[k]) for k in ("p2", "pi", "E2", "Ei"))))
    print("\norder p(dx 2 -> 1) of RMS(e) per shell (d = shell centre)")
    print("   d  | psi cos^2  psi cos^4 |  E cos^2   E cos^4")
    P = {}
    for lab, *_ in KERNELS:
        for key in ("sp", "sE"):
            P[lab, key] = EC.shell_order(res[lab, 0.5][key], res[lab, 1][key])
    for dc in P["cos^2", "sp"][0]:
        if dc > 12: break
        row = [P[k][1][np.isclose(P[k][0], dc)] for k in (("cos^2", "sp"), ("cos^4", "sp"), ("cos^2", "sE"), ("cos^4", "sE"))]
        print("%5.2f | %8.2f  %8.2f  | %8.2f  %8.2f" % (dc, *(v[0] if len(v) else np.nan for v in row)))
    F.style()
    fig, ax = plt.subplots(2, 2, figsize=(14.0, 9.0), sharex=True, gridspec_kw=dict(hspace=0.1, wspace=0.22))
    for j, (key, name) in enumerate((("sp", "\\psi"), ("sE", "E"))):
        for lab, nphys, p, _ in KERNELS:
            c = COLS[lab]
            for r, mk, ls in ((0.5, "^", "--"), (1, "o", "-")):
                s = res[lab, r][key]
                ax[0, j].semilogy(s[:, 0], s[:, 1], mk + ls, color=c, lw=1.0, ms=4,
                                  label="%s (width %d), $\\Delta x = %g$" % (lab.replace("^", "$^").replace("2", "2$").replace("4", "4$"), nphys, 1 / r))
            dd, pp = P[lab, key]
            ax[1, j].plot(dd, pp, "o-", color=c, lw=1.0, ms=4, label="%s (width %d)" % (lab.replace("^", "$^").replace("2", "2$").replace("4", "4$"), nphys))
            for a in ax[:, j]:
                a.axvspan(nphys / 2, nphys * math.sqrt(3) / 2, color=c, alpha=0.08, zorder=0)
                a.axvline(nphys / 2, color=c, lw=0.8, ls=":", zorder=0)
        for q in (1, 2):
            ax[1, j].axhline(q, color="#8a8880", lw=0.8, ls=":")
        ax[0, j].set_ylabel("RMS($e_{%s}$) / RMS($%s_{DH}$) in the shell" % (name, name))
        ax[1, j].set_ylabel("order $p(\\Delta x = 2 \\to 1)$ of RMS($e_{%s}$)" % name)
        ax[1, j].set_xlabel("distance $d = |r_i - r_0|$ (shell centre, width %g)" % EC.W_SHELL)
        ax[0, j].set_title("$%s$, $\\delta = %g$ (shaded: kernel edge, from the cube faces $n/2$ to the corners $n\\sqrt{3}/2$)"
                           % (name, DL), fontsize=9, loc="left")
        ax[0, j].legend(fontsize=7); ax[1, j].legend(fontsize=8)
        ax[1, j].set_xlim(0, 20)
    fig.savefig(os.path.join(F.ROOT, "graficos", "C_error_capas_p4.png"), dpi=DPI, bbox_inches="tight")
