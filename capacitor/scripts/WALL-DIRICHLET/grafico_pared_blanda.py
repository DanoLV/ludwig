#!/usr/bin/env python3
"""Free particle (no electrolyte, eta = 0.0005) with and without the wall soft-sphere repulsion (wall_ss_cut, R = 4, nu = 12,
eps = 1.044e-12, hc = 3), both restarted from step 30000 (d = 5.04).  (a) d(t); (b) net force on the particle vs d;
(c) induced charge on the walls (sum) vs d.   Usage: python3 grafico_pared_blanda.py"""
import os, re, glob
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
# the run without wall potential and the soft-sphere run that gets closest to the wall (the R = 4 run is not shown)
STEELE = os.environ.get("STEELE", "0") == "1"
RUNS = [("without wall repulsion", "C-particula-libre-hann08-eta0.0005-hastapared", "#eb6834"),
        ("wall soft sphere, stops at d = 3.6 (R = 3.5)", "C-particula-libre-hann08-eta0.0005-pared-d3.5", "#2a78d6")]
if STEELE:
    RUNS.append(("paper wall (Steele 10-4-3, R = 3.5)", "C-particula-libre-hann08-eta0.0005-pared-steele-d3.5", "#8e44ad"))
pat = re.compile(r"Wall charge lower\s+(\S+)\s+upper\s+(\S+)")

def load(bat):
    d = os.path.join(ROOT, bat, "z_09.50"); z = {}
    for f in glob.glob(os.path.join(d, "config.cds0*")):
        st = int(re.search(r"cds(\d+)", f).group(1)); L = [l.strip() for l in open(f)]
        z[st] = float(L[35].split()[2]) - 1.5
    F = {}
    for l in open(os.path.join(d, "proceced_data", "particle_force.csv")):
        if l.startswith("#"): continue
        p = l.strip().split(";"); F[int(p[0])] = float(p[5].replace(",", "."))
    wq = pat.findall(open(os.path.join(d, "run.log"), errors="ignore").read())
    W = {30000 + 100 * (k + 1): float(a) + float(b) for k, (a, b) in enumerate(wq)}
    return z, F, W

if __name__ == "__main__":
    plt.rcParams.update({"font.size": 10.5, "axes.grid": True, "grid.color": "#e4e2dc",
                         "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False})
    fig = plt.figure(figsize=(16, 7.4)); gs = fig.add_gridspec(2, 3, height_ratios=[1, 0.62], hspace=0.38, wspace=0.32)
    ax = [fig.add_subplot(gs[0, i]) for i in range(3)]
    for name, bat, c in RUNS:
        z, F, W = load(bat); st = sorted(k for k in z if k >= 30000 and k in F and k in W)
        t = np.array(st); d = np.array([z[k] for k in st]); f = np.array([F[k] for k in st]) * 1e5; w = np.array([W[k] for k in st])
        ax[0].plot(t, d, color=c, lw=1.8, label=name); ax[1].plot(d, f, color=c, lw=1.8); ax[2].plot(d, w, color=c, lw=1.8)
    ax[0].axhline(3.5, color="#999", ls="--", lw=1); ax[0].text(30100, 3.58, "kernel reaches wall nodes (3.5)", fontsize=8.5, color="#666")
    ax[0].set(xlabel="step", ylabel="distance to the lower wall, d", title="(a) Trajectory"); ax[0].legend(fontsize=8.5, loc="upper right")
    ax[1].set(xlabel="d", ylabel=r"reduced force on the particle $F_z\beta/q$ (total, incl. wall)", title="(b) Net force", xlim=(0, 5.2))
    ax[2].set(xlabel="d", ylabel="total induced wall charge", title="(c) Wall charge (-1 = no charge lost)", xlim=(0, 5.2), ylim=(-1.0, -0.58))
    fig.suptitle("Free particle approaching the lower grounded wall, viscosity 0.0005, no electrolyte: " + ("no wall repulsion, soft sphere and paper wall" if STEELE else "with and without the wall soft-sphere repulsion"), y=0.98)
    axt = fig.add_subplot(gs[1, :]); axt.axis("off")
    if STEELE:
        axt.text(0.5, 1.0, "Paper wall repulsion (Hoang Ngoc Minh et al., J. Chem. Phys. 158, 104103, Eqs. 7-8): Steele 10-4-3, truncated and shifted", ha="center", va="top", fontsize=11.5, fontweight="bold")
        axt.text(0.5, 0.80, r"$z=h=d-R$,   $V_{rep}(z)=V_{ads}(z)-V_{ads}(z^*)$ for $z\leq z^*=0.98698\,\sigma_w$,   $V_{rep}=0$ beyond;   "
                 r"$V_{ads}(z)=2\pi\varepsilon_w\left[\dfrac{2}{5}\left(\dfrac{\sigma_w}{z}\right)^{10}-\left(\dfrac{\sigma_w}{z}\right)^{4}-\dfrac{\sqrt{2}}{3\left(z/\sigma_w+0.61/\sqrt{2}\right)^{3}}\right]$",
                 ha="center", va="top", fontsize=13)
        axt.text(0.5, 0.52, r"$F_z=-\dfrac{dV}{dz}=2\pi\varepsilon_w\left[\dfrac{4\sigma_w^{10}}{z^{11}}-\dfrac{4\sigma_w^{4}}{z^{5}}-\dfrac{\sqrt{2}}{\sigma_w\left(z/\sigma_w+0.61/\sqrt{2}\right)^{4}}\right]$ (pushes the particle away from the wall; $F_z=0$ at $z^*$)",
                 ha="center", va="top", fontsize=13)
        axt.text(0.5, 0.24, r"Parameters: $\varepsilon_w=6.2572\times10^{-20}$ ($\beta\varepsilon_w=6.3\times10^{-15}$, $k_BT=10^{-5}$),  $\sigma_w=1.0$,  $R=3.5$ (Hann-8 kernel radius)   |   planes at $z=1.5$ and $33.5$   |   "
                 r"$F(h=0.1)=1.5726\times10^{-7}$ = electrostatic attraction at $d=3.6$ (same calibration as the soft sphere)" "\n"
                 r"Soft sphere (blue): $\varepsilon=1.3105\times10^{-21}$, $\sigma=1$, $\nu=12$, $h_c=3$.   Range of the Steele wall: $0.987$ (force active from $d=4.487$).", ha="center", va="top", fontsize=11, linespacing=1.6)
    else:
        axt.text(0.5, 1.0, "Wall repulsion (wall_ss_cut, soft sphere shifted to zero force at the cut-off)", ha="center", va="top", fontsize=11.5, fontweight="bold")
        axt.text(0.5, 0.80, r"$h=d-R$ (gap between the kernel edge and the wall plane),  $v(h)=\varepsilon\left(\dfrac{\sigma}{h}\right)^{\nu}-\varepsilon\left(\dfrac{\sigma}{h_c}\right)^{\nu}"
                 r"+\dfrac{\varepsilon\nu}{\sigma}\left(\dfrac{\sigma}{h_c}\right)^{\nu+1}(h-h_c),\qquad h<h_c;\quad v=0\;\;(h\geq h_c)$",
                 ha="center", va="top", fontsize=13)
        axt.text(0.5, 0.52, r"$F_z=-\dfrac{dv}{dh}=\dfrac{\varepsilon\nu}{\sigma}\left[\left(\dfrac{\sigma}{h}\right)^{\nu+1}-\left(\dfrac{\sigma}{h_c}\right)^{\nu+1}\right]$ (pushes the particle away from the wall; "
                 r"$F_z=0$ at $h_c$, continuous)", ha="center", va="top", fontsize=13)
        axt.text(0.5, 0.24, r"Parameters: $\varepsilon=1.3105\times10^{-21}$,  $\sigma=1.0$,  $\nu=12$,  $h_c=3.0$,  $R=3.5$ (Hann-8 kernel radius at which charge is lost)   "
                 r"|   planes at $z=1.5$ and $33.5$ (lattice units)   |   $h=0.1\Rightarrow d_{eq}\approx3.6$ (observed stop: $d=3.6055$)" "\n"
                 r"The force is already active at the start ($d=5.04\Rightarrow h=1.54<h_c$).", ha="center", va="top", fontsize=11, linespacing=1.6)
    fig.savefig(os.path.join(ROOT, "graficos", "particula_libre_pared_steele.png" if STEELE else "particula_libre_pared_blanda.png"), dpi=160, bbox_inches="tight")
    for _, bat, _c in RUNS:
        z, F, W = load(bat); print("%s: carga total de las paredes entre %.6f y %.6f" % (bat, min(W.values()), max(W.values())))
