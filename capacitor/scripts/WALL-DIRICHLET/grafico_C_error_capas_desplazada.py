#!/usr/bin/env python3
"""Error of psi and E at the lattice nodes, displaced charge, three lattices (2026-11-09, user's request): dx = 2 (Hann-4),
1 (Hann-8), 1/2 (Hann-16), delta = 1.0, 4.7, 8.3, 12.0 (the positions with the three lattices).

Pointwise error at every fluid node r_i (no interpolation): e_psi(r_i) = psi_L(r_i) - phi_DH(r_i),
e_E(r_i) = |E_L(r_i) - E_DH(r_i)| (modulus of the vector difference). psi_L, E_L: last dump (efield/kT, as elsewhere).
Reference: linear Debye-Hueckel solution with the continuous Hann source (grafico_C_hoja_hann.py) at the nodes:
  phi_DH = cap(z) + U(x, y, z) + fac <U> h(z),  U = (2 pi/A) sum_{|n|,|m| <= NK} psit(kx, ky, z) e^{i(kx x + ky y)}.
  On the nodes x = j dx the modes n and n + N (N = 32/dx) coincide, so the sum is folded mod N and done with a 2D FFT
  (exact on the nodes). E_DH,x/y: modes times -i kx, -i ky; E_DH,z: closed-form dpsit (HH.dU), cap' and h'.
  Check printed: U and dU/dz on the axis through the charge vs HH.U, HH.dU.
1) Shells of distance d = |r_i - r0| (minimum image in x, y), width W_SHELL (physical units, parameter w=):
   RMS(e)/RMS(ref) in each shell (absolute error divided by the RMS of the reference in the same shell; user's choice),
   max|e|/RMS(ref), and the order p = log2 of the ratio of RMS(e) of two lattices in the same shell.
2) Global norms over all fluid nodes: ||e||_2 = sqrt(sum dx^3 e^2), ||e||_inf = max|e|; printed absolute and divided by
   the same norm of the reference, with their orders.
Not separated: nonlinear Poisson-Boltzmann part (Ludwig nonlinear, phi_DH linear).
Usage: python3 grafico_C_error_capas_desplazada.py [w=1]
       -> graficos/C_error_capas_desplazada.png (shells), graficos/C_error_capas_desplazada_orden.png (orders)
"""
import os, sys, math
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import grafico_C_formulas as F
import grafico_C_hoja_hann as HH
import grafico_C_hoja_gauss as H
from grafico_C_refinamiento import load
from grafico_C_refinamiento_zoom import KT
from grafico_C_desplazada import set_charge, psi_files
from grafico_C_convergencia_desplazada import run_path, DPI

W_SHELL = float(([a[2:] for a in sys.argv[1:] if a.startswith("w=")] or ["1"])[0])
DELTAS = (1.0, 4.7, 8.3, 12.0)
LATS = ((0.5, "$\\Delta x = 2$ (Hann-4)", "#eb6834", "^"), (1, "$\\Delta x = 1$ (Hann-8)", "#2b2a27", "o"),
        (2, "$\\Delta x = 1/2$ (Hann-16)", "#2a78d6", "s"))

def dh_on_nodes(r, mu):
    """phi_DH and E_DH (3 comps) on the fluid nodes, arrays (n, n, n) in Ludwig's index order, charge at x = y = 16r."""
    n = int(round(32 * r)); ic = int(round(16 * r)) - 1
    z = (np.arange(2, n + 2) - 1.5) / r - HH.L / 2          # fluid nodes k = 1..n (array index), physical z
    nn, mm = np.meshgrid(np.arange(-HH.NK, HH.NK + 1), np.arange(-HH.NK, HH.NK + 1), indexing="ij")
    nn = nn.ravel(); mm = mm.ravel(); kx = 2 * np.pi * nn / HH.L; ky = 2 * np.pi * mm / HH.L
    C = np.zeros((4, n, n, n), dtype=complex)               # U, dU/dx, dU/dy, dU/dz folded mod n, per z
    for i0 in range(0, len(kx), 4000):
        sl = slice(i0, i0 + 4000)
        p = HH.psit(kx[sl], ky[sl], z); dp = HH.dpsit(kx[sl], ky[sl], z)
        for c, v in enumerate((p, 1j * kx[sl][:, None] * p, 1j * ky[sl][:, None] * p, dp)):
            np.add.at(C[c], (nn[sl] % n, mm[sl] % n), v)
    # U(x_j, y_l, z) = (2 pi/A) sum_nm C_nm e^{2 pi i (n j + m l)/n}: numpy ifft2 has the 1/n^2 factor
    G = np.real(np.fft.ifft2(C, axes=(1, 2))) * n * n * 2 * np.pi / HH.A
    G = np.roll(G, (ic, ic), axis=(1, 2))                    # j = 0 (x = charge) -> array index ic
    fac = (HH.KAP * HH.L / 2) / math.tanh(HH.KAP * HH.L / 2)
    phi = H.cap_closed(z, 0, 0)[None, None, :] + G[0] + (fac * mu * H.h(z))[None, None, :]
    dzc = fac * mu * (-HH.KAP * np.sinh(HH.KAP * z) / np.cosh(HH.KAP * HH.L / 2))   # grounded plates: cap' = 0
    E = -np.stack((G[1], G[2], G[3] + dzc[None, None, :]), axis=-1)
    return phi, E, z

def distances(r, dl):
    n = int(round(32 * r)); ic = int(round(16 * r)) - 1
    d1 = (np.arange(n) - ic) / r; d1 = d1 - HH.L * np.round(d1 / HH.L)   # minimum image in x, y
    z = (np.arange(2, n + 2) - 1.5) / r - HH.L / 2
    return np.sqrt(d1[:, None, None] ** 2 + d1[None, :, None] ** 2 + (z - dl)[None, None, :] ** 2)

def shells(d, e, ref):
    """(centre, RMS(e)/RMS(ref), max|e|/RMS(ref), RMS(e)) per shell of width W_SHELL."""
    k = np.floor(d / W_SHELL).astype(int).ravel(); e = np.abs(e).ravel(); ref = np.abs(ref).ravel(); out = []
    for s in range(k.max() + 1):
        m = k == s
        if not m.any(): continue
        re, rr = np.sqrt(np.mean(e[m] ** 2)), np.sqrt(np.mean(ref[m] ** 2))
        out.append(((s + 0.5) * W_SHELL, re / rr, e[m].max() / rr, re))
    return np.array(out)

def shell_order(sc, sf):
    """(centres, p = log2 RMS(e) coarse / fine) on the shells present in BOTH lattices, matched by shell centre.
    (2026-11-09 fix: the shells were matched by index, wrong when a shell has no node in one lattice, e.g. the
    coarse lattice has no node at 1 <= d < 2 for delta = 1 and none at d < 1 for delta = 12.)"""
    cc, ic, jf = np.intersect1d(np.round(sc[:, 0], 9), np.round(sf[:, 0], 9), return_indices=True)
    return sc[ic, 0], np.log2(sc[ic, 3] / sf[jf, 3])

if __name__ == "__main__":
    res = {}
    print("shell width %g (physical units)" % W_SHELL)
    for dl in DELTAS:
        set_charge(dl); mu = HH.mean_U()
        for r, *_ in LATS:
            path = run_path(dl, r); fs = psi_files(path); assert fs, "not finished: " + path
            n = int(round(32 * r)); ic = int(round(16 * r)) - 1
            psiL = load(fs[-1], r)[:, :, 1:n + 1]
            EL = np.loadtxt(os.path.join(F.ROOT, path, "efield-" + os.path.basename(fs[-1])[4:])).reshape(n, n, n + 2, 3)[:, :, 1:n + 1] / KT * r
            phi, E, z = dh_on_nodes(r, mu)
            if r == 1:                                       # check of the FFT evaluation on the axis
                cu = np.abs(phi[ic, ic] - HH.phi(0.0, z, 0, 0, mu)).max(); cd = np.abs(-E[ic, ic, :, 2] - HH.dphi(0.0, z, 0, 0, mu)).max()
                print("delta %.1f check on the axis: |phi_FFT - HH.phi| max %.1e, |dphi_FFT - HH.dphi| max %.1e" % (dl, cu, cd))
            ep = psiL - phi; eE = np.linalg.norm(EL - E, axis=-1); Em = np.linalg.norm(E, axis=-1)
            d = distances(r, dl); dv = (1 / r) ** 3
            g = dict(p2=math.sqrt(dv * np.sum(ep ** 2)), p2r=math.sqrt(np.sum(ep ** 2) / np.sum(phi ** 2)),
                     pi=np.abs(ep).max(), pir=np.abs(ep).max() / np.abs(phi).max(),
                     E2=math.sqrt(dv * np.sum(eE ** 2)), E2r=math.sqrt(np.sum(eE ** 2) / np.sum(Em ** 2)),
                     Ei=eE.max(), Eir=eE.max() / Em.max())
            res[dl, r] = dict(sp=shells(d, ep, phi), sE=shells(d, eE, Em), g=g)
            print("delta %4.1f dx = %-3g | psi: L2 %.3e (rel %.3e)  Linf %.3e (rel %.3e) | E: L2 %.3e (rel %.3e)  Linf %.3e (rel %.3e)"
                  % (dl, 1 / r, g["p2"], g["p2r"], g["pi"], g["pir"], g["E2"], g["E2r"], g["Ei"], g["Eir"]))
        for (rc, rf) in ((0.5, 1), (1, 2)):
            gc, gf = res[dl, rc]["g"], res[dl, rf]["g"]
            print("delta %4.1f p(dx %g -> %g): psi L2 %.2f Linf %.2f | E L2 %.2f Linf %.2f" % (dl, 1 / rc, 1 / rf,
                  *(math.log2(gc[k] / gf[k]) for k in ("p2", "pi", "E2", "Ei"))))
    F.style()
    for kind in ("err", "orden"):
        fig, ax = plt.subplots(len(DELTAS), 2, figsize=(14.0, 4.2 * len(DELTAS)), sharex=True, gridspec_kw=dict(hspace=0.12, wspace=0.2))
        for i, dl in enumerate(DELTAS):
            for j, (key, name) in enumerate((("sp", "\\psi"), ("sE", "E"))):
                a = ax[i, j]
                if kind == "err":
                    for r, lab, c, mk in LATS:
                        s = res[dl, r][key]
                        a.semilogy(s[:, 0], s[:, 1], mk + "-", color=c, lw=1.0, ms=4, label="RMS, " + lab)
                        a.semilogy(s[:, 0], s[:, 2], mk + ":", color=c, lw=0.8, ms=3, alpha=0.6, label="max, " + lab)
                    a.set_ylabel("$|e_{%s}|$ / RMS($%s_{DH}$) in the shell" % (name, name))
                else:
                    for (rc, rf), lab, c, mk in (((0.5, 1), "$p(\\Delta x = 2 \\to 1)$", "#eb6834", "^-"),
                                                 ((1, 2), "$p(\\Delta x = 1 \\to 1/2)$", "#2a78d6", "s-")):
                        a.plot(*shell_order(res[dl, rc][key], res[dl, rf][key]), mk, color=c, lw=1.0, ms=4, label=lab)
                    for q in (1, 2):
                        a.axhline(q, color="#8a8880", lw=0.8, ls=":")
                    a.set_ylabel("order $p$ of RMS($e_{%s}$) per shell" % name)
                a.axvline(HH.NH / 2, color="#8a8880", lw=0.8, ls="--", zorder=0)    # kernel half width
                a.axvline(16 - dl, color="#c0392b", lw=0.8, ls="--", zorder=0)     # distance to the upper plate
                a.set_title("$%s$, $\\delta = %g$  (grey: kernel half width 4; red: distance to the upper plate)" % (name, dl),
                            fontsize=9, loc="left")
                if i == 0: a.legend(fontsize=7, ncol=2 if kind == "err" else 1)
                if i == len(DELTAS) - 1: a.set_xlabel("distance $d = |r_i - r_0|$ (shell centre, width %g)" % W_SHELL)
        fig.savefig(os.path.join(F.ROOT, "graficos", "C_error_capas_desplazada%s.png" % ("" if kind == "err" else "_orden")),
                    dpi=DPI, bbox_inches="tight")
