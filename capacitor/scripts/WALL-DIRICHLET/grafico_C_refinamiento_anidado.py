#!/usr/bin/env python3
"""Nested lattice refinement with the Hann kernel: dx = 2, 1, 1/2 (Hann-4, -8, -16; same physical system, built by
crear_refinado.py 4 / the Hann-8 base / crear_refinado.py 16), local convergence order from the RMS of the nodal
residual on bands (rho = 0, z axis through the charge), vs the analytic continuous-Hann solution (grafico_C_hoja_hann.py).

Why bands and not nodes: the wall sits at half a link, so nodes are cell centres z + 16 = (m + 1/2) dx and never
coincide between lattices (halving dx splits each cell in two, with new centres at 1/4 and 3/4 of the old cell).
What does nest are the CELL FACES z + 16 = m dx: every face of dx = 2 is a face of dx = 1 and of dx = 1/2.
A band whose ends are faces of the coarse lattice (even z + 16) is tiled by whole cells in all three lattices, so
    RMS_[a,b] = ( (1/(b-a)) sum_{z_k in [a,b]} d_k^2 dx )^(1/2)
is the midpoint rule of (1/(b-a)) int d^2 dz in every lattice, with no clipped cells and no interpolation; the same
band covers exactly 1 + ... coarse cells, twice as many medium and four times as many fine ones.
Two families: width 2 every 2 (one coarse cell per band: contiguous) and width 4 every 2 (overlap 2).
Usage: python3 grafico_C_refinamiento_anidado.py
"""
import os
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import grafico_C_formulas as F
import grafico_C_hoja_hann as HH
from grafico_C_refinamiento import Z0, L
from grafico_C_refinamiento_gauss import series

NRUNS = [(0.5, "C-electrolito-hann4-LD3-ref0.5/z_05.50"), (1.0, "C-electrolito-hann8-LD6/z_09.50"),
         (2.0, "C-electrolito-hann16-LD12-ref2/z_17.50")]
NCOL = {0.5: "#c0392b", 1.0: "#eb6834", 2.0: "#8e44ad"}
FAMS = ((2.0, 2.0, "o-", "width 2, every 2"), (4.0, 2.0, "s--", "width 4, every 2 (overlap 2)"))

def band_rms(z, d, r, W, step):
    """RMS on bands [s, s+W], s = -16, -16+step, ...: ends on faces of dx = 2, so every band holds whole cells."""
    S = np.arange(-L / 2, L / 2 - W + 1e-9, step); out = []
    for s in S:
        m = (z > s) & (z < s + W)
        assert m.sum() == int(round(W * r)), (s, W, r, m.sum())
        out.append(np.sqrt(np.mean(d[m] ** 2)))
    return S + W / 2, np.array(out)

if __name__ == "__main__":
    F.style()
    SX = series(NRUNS, HH, HH.mean_U())
    fig, ax = plt.subplots(2, 1, figsize=(10.0, 8.0), sharex=True, gridspec_kw=dict(hspace=0.08))
    for W, step, sty, lab in FAMS:
        R = []
        for r, z, d in SX:
            zc, v = band_rms(z, d, r, W, step); R.append(v)
            ax[0].semilogy(zc, v, sty, ms=4.5 if W == 2 else 3.5, lw=1.3 if W == 2 else 0.9, color=NCOL[r],
                           alpha=1.0 if W == 2 else 0.6, label="$\\Delta x=%g$, %s" % (1 / r, lab))
        R = np.array(R)
        print("bands %s:  centre " % lab + " ".join("%6.1f" % c for c in zc))
        for k in range(len(SX) - 1):
            p = np.log(R[k] / R[k + 1]) / np.log(2.0)
            dx0, dx1 = 1 / SX[k][0], 1 / SX[k + 1][0]
            print("   order %g -> %g          " % (dx0, dx1) + " ".join("%6.2f" % v for v in p))
            ax[1].plot(zc, p, sty, ms=5 if W == 2 else 4, lw=1.6 if W == 2 else 1.0, color=NCOL[SX[k + 1][0]],
                       alpha=1.0 if W == 2 else 0.6, label="$%g \\to %s$, %s" % (dx0, "1/2" if dx1 == 0.5 else "%g" % dx1, lab))
    for a in ax:
        a.axvspan(Z0 - 4, Z0 + 4, color="#f1efe8", zorder=0); a.axvline(Z0, color="#b9b7ae", lw=0.8, ls=":")
        for x in np.arange(-L / 2, L / 2 + 1e-9, 2.0):
            a.axvline(x, color="#e4e2dc", lw=0.5, zorder=0)
    ax[1].axhline(2, color="#55544f", lw=0.9, ls=":"); ax[1].axhline(1, color="#b9b7ae", lw=0.9, ls="-.")
    ax[1].set_ylim(-1, 4); ax[0].set_xlim(-16, 16)
    ax[0].set_ylabel("band RMS of $\\psi_{\\rm Ludwig}-\\psi_{\\rm analytic}$\n"
                     "(potential $\\psi = e\\phi/k_BT$, axis $\\rho=0$)")
    ax[1].set_ylabel("local order $\\log_2({\\rm RMS}_{\\Delta x}/{\\rm RMS}_{\\Delta x/2})$")
    ax[1].set_xlabel("band centre z (physical units, from the mid-plane)")
    ax[0].legend(fontsize=8, loc="upper right", ncol=2); ax[1].legend(fontsize=8, loc="lower right", ncol=2)
    ax[0].set_title("Electrostatic potential: error of Ludwig vs the analytic solution (continuous Hann charge), "
                    "electrolyte capacitor\nnested lattices $\\Delta x = 2, 1, 1/2$ (Hann-4, -8, -16): band ends on "
                    "faces of all three (thin lines: faces of $\\Delta x=2$); shaded: charge support", fontsize=10.5)
    fig.savefig(os.path.join(F.ROOT, "graficos", "C_refinamiento_anidado_hann.png"), dpi=170, bbox_inches="tight")
