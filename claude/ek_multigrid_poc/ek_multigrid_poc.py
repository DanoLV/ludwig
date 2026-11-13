#!/usr/bin/env python3
"""
ek_multigrid_poc.py -- proof of concept for running the electrokinetics
(Poisson + Nernst-Planck + subgrid charges) on a grid refined by an integer
factor r over the LB grid (see memory ek-multigrid-plan, 2026-10-02).

The fine grid is treated as an ordinary lattice with its OWN lattice units
(dx_f = 1), exactly as the C plan proposes, so the solvers are written once
and know nothing about r.  Everything that depends on r lives in the
conversion layer:

    length      x_f = r (x_c - 0.5) + 0.5        (cell-centred alignment)
    epsilon     eps_f = eps_c / r
    diffusion   D_f   = D_c r^2
    density     c_f   = c_c / r^3                (per site)
    velocity    u_f   = u_c r                    (prolongation)
    field       E_f   = E_c / r
    force       F_c   = r * F_f                  (restriction: sum children)

Coarse (LB) units: node ic at x = ic, cell [ic-0.5, ic+0.5], beta = e = 1.

Tests
  1  1D capacitor (Dirichlet walls, Asta factor-2 rule): PNP transient and
     equilibrium vs r.  Checks eps/r, D r^2, c/r^3 and the wall alignment.
  2  Electro-osmotic flow: Stokes on the COARSE grid driven by the force
     restricted from the fine grid.  Checks the r*sum(children) restriction.
  3  2D advection of an ion blob with a prolonged Taylor-Green velocity:
     trilinear vs piecewise-constant prolongation; mass conservation.
  4  3D two-particle force, FFT Poisson on the fine grid, Peskin-4 kernel
     defined in fine nodes (A) or keeping the coarse physical width (B),
     against an Ewald reference.

Usage:  python3 ek_multigrid_poc.py [--quick] [--outdir DIR]
"""

import argparse
import math
import os
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# Conversion layer (the only place where r appears)
# ---------------------------------------------------------------------------

def c2f_pos(x_c, r):
    """Coarse coordinate -> fine coordinate (cell-centred alignment)."""
    return r * (np.asarray(x_c) - 0.5) + 0.5


def f2c_pos(x_f, r):
    return (np.asarray(x_f) - 0.5) / r + 0.5


def fine_centres_in_coarse(n_coarse, r):
    """Coarse coordinates of the fine nodes m = 1..r*n_coarse."""
    m = np.arange(1, r * n_coarse + 1)
    return f2c_pos(m, r)


# ---------------------------------------------------------------------------
# Test 1: 1D PNP capacitor
# ---------------------------------------------------------------------------

def pnp1d(r, N, eps_c, D_c, c0_c, psi_lo, psi_hi, nsteps, record_every=1):
    """
    Symmetric 1:1 electrolyte between two Dirichlet walls located at the
    faces z = 0.5 and z = N + 0.5 (coarse units), i.e. where LB bounce-back
    puts them.  Fine grid of n = r N cells in fine units.

    Poisson: finite volume, Asta rule (factor 2 on the fluid-wall link).
    NP: Ludwig-like link flux j = -D [dc + z c_avg dpsi], no flux at walls,
    explicit Euler with multisteps chosen from D_f (so they grow as r^2).

    Returns dict with time series (coarse units) and final fields.
    """
    n = r * N
    eps_f = eps_c / r
    D_f = D_c * r * r
    c_f = c0_c / r**3

    # Poisson matrix (fine units): sum_links f (psi_j - psi_i) = -rho_i/eps_f
    A = np.zeros((n, n))
    for i in range(n):
        A[i, i] = -2.0
        if i > 0:
            A[i, i - 1] = 1.0
        if i < n - 1:
            A[i, i + 1] = 1.0
    A[0, 0] = -3.0        # link to wall counted twice (factor 2): -1 -2
    A[-1, -1] = -3.0
    Ainv = np.linalg.inv(A)
    bc = np.zeros(n)
    bc[0] = 2.0 * psi_lo
    bc[-1] = 2.0 * psi_hi

    def poisson(cp, cm):
        return Ainv @ (-(cp - cm) / eps_f - bc)

    cp = np.full(n, c_f)
    cm = np.full(n, c_f)

    nsub = max(1, math.ceil(D_f / 0.25))
    dt = 1.0 / nsub

    t_rec, sig_rec, qion_rec, neut_rec = [], [], [], []

    def record(t, psi):
        # Wall charge from Gauss (coarse units): sigma = eps_c * E_c,
        # E_c = r * E_f, E_f at the wall face = (psi_lo - psi_0)/0.5
        sigma_lo = eps_c * r * (psi_lo - psi[0]) / 0.5
        # Ionic charge per unit COARSE area in the lower half:
        # r^2 fine columns per coarse area.
        qion_lo = r * r * np.sum((cp - cm)[: n // 2])
        # Global neutrality: sigma_lo + sigma_hi + Q_ion = 0
        sigma_hi = eps_c * r * (psi_hi - psi[-1]) / 0.5
        neut_rec.append(sigma_lo + sigma_hi + r * r * np.sum(cp - cm))
        t_rec.append(t)
        sig_rec.append(sigma_lo)
        qion_rec.append(qion_lo)

    psi = poisson(cp, cm)
    record(0, psi)
    for step in range(1, nsteps + 1):
        for _ in range(nsub):
            psi = poisson(cp, cm)
            dpsi = psi[1:] - psi[:-1]
            for c, z in ((cp, +1.0), (cm, -1.0)):
                j = -D_f * ((c[1:] - c[:-1]) + z * 0.5 * (c[1:] + c[:-1]) * dpsi)
                div = np.zeros(n)
                div[:-1] += j
                div[1:] -= j
                c -= dt * div
        if step % record_every == 0:
            record(step, psi)

    psi = poisson(cp, cm)
    return dict(r=r, n=n, nsub=nsub, t=np.array(t_rec),
                sigma=np.array(sig_rec), qion=np.array(qion_rec), neut=np.array(neut_rec),
                psi=psi, cp=cp, cm=cm, eps_f=eps_f,
                mass=r * r * (cp.sum() + cm.sum()),
                zc=fine_centres_in_coarse(N, r))


def dh_capacitor(z, N, lam, psi0):
    """Debye-Hueckel potential for walls at 0.5 (psi0) and N+0.5 (-psi0)."""
    zc = 0.5 * (N + 1)
    L = float(N)
    return psi0 * np.sinh((zc - z) / lam) / np.sinh(L / (2 * lam))


def dh_capacitor_dpsi(z, N, lam, psi0):
    zc = 0.5 * (N + 1)
    L = float(N)
    return -psi0 * np.cosh((zc - z) / lam) / (lam * np.sinh(L / (2 * lam)))


def test1(args, out):
    print("\n=== Test 1: 1D PNP capacitor (units: eps/r, D r^2, c/r^3) ===")
    N = 24
    lam = 1.5                      # Debye length, coarse units (under-resolved)
    c0 = 0.01
    eps = 2.0 * c0 * lam**2        # lambda^2 = eps/(2 c0) with beta=e=1
    D = 0.05
    psi0 = 0.02                    # small: Debye-Hueckel regime
    nsteps = 1500 if args.quick else 4000
    rs = [1, 2, 4, 8] if args.quick else [1, 2, 4, 8, 16]

    print(f"N={N} lambda_D={lam} c0={c0} eps={eps:.4g} D={D} psi_w=+-{psi0}")
    print(f"tau_Debye = {lam**2/D:.0f}  tau_RC = {lam*N/(2*D):.0f} LB steps,"
          f" run {nsteps}")

    res = {}
    for r in rs:
        t0 = time.time()
        res[r] = pnp1d(r, N, eps, D, c0, psi0, -psi0, nsteps)
        print(f"  r={r:2d}: n_fine={res[r]['n']:4d}  NP substeps/LB step="
              f"{res[r]['nsub']:3d}  ({time.time()-t0:.1f} s)")

    sigma_dh = eps * psi0 / (lam * math.tanh(N / (2 * lam)))
    rref = rs[-1]
    print(f"\n  sigma_wall Debye-Hueckel (continuum) = {sigma_dh:.6e}")
    print("   r   sigma_eq        err_vs_DH   err_vs_rmax  |psi-DH|_inf  "
          "neutrality   mass drift")
    m0 = 2 * N * c0
    for r in rs:
        R = res[r]
        s = R["sigma"][-1]
        psi_dh = dh_capacitor(R["zc"], N, lam, psi0)
        print(f"  {r:2d}  {s:.6e}  {(s-sigma_dh)/sigma_dh:+.3e}  "
              f"{(s-res[rref]['sigma'][-1])/res[rref]['sigma'][-1]:+.3e}  "
              f"{np.max(np.abs(R['psi']-psi_dh))/psi0:.3e}     "
              f"{np.max(np.abs(R['neut']))/s:.1e}  {(R['mass']-m0)/m0:+.1e}")

    # transient: compare sigma(t) with the finest grid at the same physical t
    print("\n  transient sigma(t)/sigma_DH  (same LB times for every r)")
    tt = [10, 50, 200, 800, nsteps]
    print("   r  " + "".join(f"  t={t:<6d}" for t in tt))
    for r in rs:
        print(f"  {r:2d}  " + "".join(
            f"  {res[r]['sigma'][t]/sigma_dh:.5f}" for t in tt))

    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    for r in rs:
        R = res[r]
        ax[0].plot(R["t"], R["sigma"] / sigma_dh, label=f"r={r}")
        ax[1].plot(R["zc"], R["psi"] / psi0, ".-", ms=3, label=f"r={r}")
    zz = np.linspace(0.5, N + 0.5, 400)
    ax[1].plot(zz, dh_capacitor(zz, N, lam, psi0) / psi0, "k--", lw=1,
               label="DH continuum")
    ax[0].set_xscale("log")
    ax[0].set_xlabel("t (LB steps)")
    ax[0].set_ylabel("sigma_wall / sigma_DH")
    ax[0].legend()
    ax[1].set_xlim(0.5, 6)
    ax[1].set_xlabel("z (coarse units)")
    ax[1].set_ylabel("psi / psi_wall")
    ax[1].legend()
    fig.suptitle("Test 1: capacitor PNP on fine grid")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "test1_capacitor.png"), dpi=130)
    plt.close(fig)
    return dict(res=res, N=N, lam=lam, eps=eps, psi0=psi0)


# ---------------------------------------------------------------------------
# Test 2: electro-osmotic flow, force restriction
# ---------------------------------------------------------------------------

def stokes_coarse(f, eta):
    """eta u'' = -f on N coarse cells, no slip at the wall faces."""
    N = len(f)
    A = np.zeros((N, N))
    for i in range(N):
        A[i, i] = -2.0
        if i > 0:
            A[i, i - 1] = 1.0
        if i < N - 1:
            A[i, i + 1] = 1.0
    A[0, 0] = A[-1, -1] = -3.0     # ghost u = -u_1 (mid-link wall)
    return np.linalg.solve(eta * A, -f)


def restrict_force_1d(fz_fine_per_site, r):
    """
    Fine force per fine site (fine units) -> coarse force per coarse site
    (coarse units).  A coarse cell has r^3 children: r along z times r^2
    identical columns in the 1D slab.  F_c = r * sum(children).
    """
    N = len(fz_fine_per_site) // r
    return r * (r * r) * fz_fine_per_site.reshape(N, r).sum(axis=1)


def test2(args, out, t1):
    print("\n=== Test 2: electro-osmosis, force restriction F_c = r*sum ===")
    N, lam, eps, psi0 = t1["N"], t1["lam"], t1["eps"], t1["psi0"]
    E_c = 1.0e-2                     # applied field along x, coarse units
    eta = 1.0 / 6.0
    z_nodes = np.arange(1, N + 1, dtype=float)

    # analytic: u = (eps E/eta)(psi - psi_linear)
    psi_lin = psi0 * (1.0 - 2.0 * (z_nodes - 0.5) / N)
    u_exact = (eps * E_c / eta) * (dh_capacitor(z_nodes, N, lam, psi0) - psi_lin)

    # floor: coarse Stokes with the EXACT cell-integrated force
    #   int_cell rho E dz = -eps E [psi'(z+1/2) - psi'(z-1/2)]
    f_exact = -eps * E_c * (dh_capacitor_dpsi(z_nodes + 0.5, N, lam, psi0)
                            - dh_capacitor_dpsi(z_nodes - 0.5, N, lam, psi0))
    u_floor = stokes_coarse(f_exact, eta)
    umax = np.max(np.abs(u_exact))
    print(f"  max|u| analytic = {umax:.4e}; coarse-Stokes floor error "
          f"(exact cell force) = {np.max(np.abs(u_floor-u_exact))/umax:.3e}")
    print("   r   |u-u_exact|/max   |u-u_floor|/max   |f_c-f_exact|/max|f|"
          "   (no x r factor: |u-u_floor|/max)")

    fig, ax = plt.subplots(figsize=(6, 4))
    for r, R in t1["res"].items():
        rho_f = R["cp"] - R["cm"]                 # charge per fine site
        E_f = E_c / r                             # field in fine units
        f_c = restrict_force_1d(rho_f * E_f, r)
        f_wrong = f_c / r                         # forgetting the r factor
        u = stokes_coarse(f_c, eta)
        u_w = stokes_coarse(f_wrong, eta)
        print(f"  {r:2d}   {np.max(np.abs(u-u_exact))/umax:.3e}         "
              f"{np.max(np.abs(u-u_floor))/umax:.3e}         "
              f"{np.max(np.abs(f_c-f_exact))/np.max(np.abs(f_exact)):.3e}"
              f"             {np.max(np.abs(u_w-u_floor))/umax:.3e}")
        ax.plot(z_nodes, u / umax, "o-", ms=3, label=f"r={r}")
    zz = np.linspace(0.5, N + 0.5, 400)
    pl = psi0 * (1.0 - 2.0 * (zz - 0.5) / N)
    ax.plot(zz, (eps * E_c / eta) * (dh_capacitor(zz, N, lam, psi0) - pl) / umax,
            "k--", lw=1, label="analytic")
    ax.set_xlabel("z (coarse)")
    ax.set_ylabel("u_x / max")
    ax.legend()
    ax.set_title("Test 2: EOF, LB-grid Stokes + restricted force")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "test2_eof.png"), dpi=130)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Test 3: 2D advection with prolonged velocity
# ---------------------------------------------------------------------------

def tg_velocity(x, y, U, N):
    k = 2 * np.pi / N
    return (U * np.sin(k * x) * np.cos(k * y),
            -U * np.cos(k * x) * np.sin(k * y))


def prolong_bilinear(uc, r):
    """Coarse node values (node ic at x=ic, periodic) -> fine centres."""
    N = uc.shape[0]
    X = fine_centres_in_coarse(N, r) - 1.0      # 0-based node units
    i0 = np.floor(X).astype(int)
    w = X - i0
    i0 %= N
    i1 = (i0 + 1) % N
    # separable bilinear
    a = (1 - w)[:, None] * uc[i0, :] + w[:, None] * uc[i1, :]
    return (1 - w)[None, :] * a[:, i0] + w[None, :] * a[:, i1]


def prolong_constant(uc, r):
    return np.repeat(np.repeat(uc, r, axis=0), r, axis=1)


def advect2d(c, ux, uy, nsteps, nsub):
    """Conservative upwind advection, link velocity = mean of the two nodes."""
    dt = 1.0 / nsub
    for _ in range(nsteps * nsub):
        for u, ax in ((ux, 0), (uy, 1)):
            ub = 0.5 * (u + np.roll(u, -1, axis=ax))
            cu = np.where(ub > 0, c, np.roll(c, -1, axis=ax))
            flux = ub * cu
            c = c - dt * (flux - np.roll(flux, 1, axis=ax))
    return c


def div2d(ux, uy):
    fx = 0.5 * (ux + np.roll(ux, -1, 0))
    fy = 0.5 * (uy + np.roll(uy, -1, 1))
    return (fx - np.roll(fx, 1, 0)) + (fy - np.roll(fy, 1, 1))


def test3(args, out):
    print("\n=== Test 3: 2D advection, velocity prolongation (u_f = r u_c) ===")
    N = 16
    U = 0.02
    nsteps = 100 if args.quick else 300
    xc = np.arange(1, N + 1, dtype=float)
    Xc, Yc = np.meshgrid(xc, xc, indexing="ij")
    uxc, uyc = tg_velocity(Xc, Yc, U, N)
    print(f"  Taylor-Green N={N} U={U}, {nsteps} LB steps")
    print("   r  prolong    rms div/(rU)   L2 err vs exact-u   mass drift")
    for r in (2, 4):
        xf = fine_centres_in_coarse(N, r)
        Xf, Yf = np.meshgrid(xf, xf, indexing="ij")
        uxe, uye = tg_velocity(Xf, Yf, U, N)
        c0 = np.exp(-((Xf - 0.6 * N) ** 2 + (Yf - 0.35 * N) ** 2) / (2 * 2.0**2))
        c0 /= r**3                                 # per fine site
        nsub = max(1, math.ceil(2 * r * U / 0.2))
        cref = advect2d(c0.copy(), r * uxe, r * uye, nsteps, nsub)
        for name, P in (("bilinear", prolong_bilinear),
                        ("constant", prolong_constant)):
            ux, uy = r * P(uxc, r), r * P(uyc, r)
            c = advect2d(c0.copy(), ux, uy, nsteps, nsub)
            dv = np.sqrt(np.mean(div2d(ux, uy) ** 2)) / (r * U)
            err = np.linalg.norm(c - cref) / np.linalg.norm(cref - c0)
            print(f"  {r:2d}  {name:9s}  {dv:.3e}      {err:.3e}"
                  f"          {(c.sum()-c0.sum())/c0.sum():+.1e}")
        dv = np.sqrt(np.mean(div2d(r * uxe, r * uye) ** 2)) / (r * U)
        print(f"  {r:2d}  exact u    {dv:.3e}      (reference)")
    print("  (L2 err normalised by how far the blob moved)")


# ---------------------------------------------------------------------------
# Test 4: 3D particle forces on the fine grid
# ---------------------------------------------------------------------------

def peskin4(x):
    a = np.abs(x)
    w = np.zeros_like(a)
    m1 = a < 1
    m2 = (a >= 1) & (a < 2)
    w[m1] = (3 - 2 * a[m1] + np.sqrt(1 + 4 * a[m1] - 4 * a[m1] ** 2)) / 8
    w[m2] = (5 - 2 * a[m2] - np.sqrt(-7 + 12 * a[m2] - 4 * a[m2] ** 2)) / 8
    return w


def kernel_1d(xf, s, n):
    """Nodes and weights for a Peskin-4 of width s fine nodes at fine coord xf."""
    m = np.arange(math.floor(xf - 2 * s) + 1, math.floor(xf + 2 * s) + 1)
    w = peskin4((xf - m) / s) / s
    return (m - 1) % n, w                # node m (1-based) -> array index


def pm_forces(pos_c, q, L, r, eps_c, width_fine):
    """Particle-mesh forces (coarse units) with all mesh work on the fine grid."""
    n = r * L
    eps_f = eps_c / r
    rho = np.zeros((n, n, n))
    stencils = []
    for p, qi in zip(pos_c, q):
        xf = c2f_pos(p, r)
        st = [kernel_1d(xf[d], width_fine, n) for d in range(3)]
        stencils.append(st)
        (ix, wx), (iy, wy), (iz, wz) = st
        rho[np.ix_(ix, iy, iz)] += qi * np.einsum("i,j,k->ijk", wx, wy, wz)

    k = 2 * np.pi * np.fft.fftfreq(n)
    kz = 2 * np.pi * np.fft.rfftfreq(n)
    K2 = (2 - 2 * np.cos(k))[:, None, None] + (2 - 2 * np.cos(k))[None, :, None] \
        + (2 - 2 * np.cos(kz))[None, None, :]
    rk = np.fft.rfftn(rho)
    K2[0, 0, 0] = 1.0
    pk = rk / (eps_f * K2)
    pk[0, 0, 0] = 0.0
    psi = np.fft.irfftn(pk, s=rho.shape, axes=(0, 1, 2))
    E = [-(np.roll(psi, -1, d) - np.roll(psi, 1, d)) / 2.0 for d in range(3)]

    F = np.zeros((len(q), 3))
    for a, (st, qi) in enumerate(zip(stencils, q)):
        (ix, wx), (iy, wy), (iz, wz) = st
        W = np.einsum("i,j,k->ijk", wx, wy, wz)
        for d in range(3):
            F[a, d] = r * qi * np.sum(W * E[d][np.ix_(ix, iy, iz)])  # F_c = r F_f
    return F


def ewald_forces(pos, q, L, eps, alpha=0.35, kmax=10):
    C = 1.0 / (4 * np.pi * eps)
    npart = len(q)
    F = np.zeros((npart, 3))
    rng = range(-2, 3)
    for i in range(npart):
        for j in range(npart):
            for nx in rng:
                for ny in rng:
                    for nz in rng:
                        if i == j and nx == ny == nz == 0:
                            continue
                        rv = pos[i] - pos[j] + L * np.array([nx, ny, nz])
                        d = np.linalg.norm(rv)
                        g = (math.erfc(alpha * d) / d**2
                             + 2 * alpha / math.sqrt(math.pi)
                             * math.exp(-(alpha * d) ** 2) / d)
                        F[i] += C * q[i] * q[j] * g * rv / d
    nv = np.arange(-kmax, kmax + 1)
    KX, KY, KZ = np.meshgrid(nv, nv, nv, indexing="ij")
    kv = 2 * np.pi / L * np.stack([KX.ravel(), KY.ravel(), KZ.ravel()], 1)
    kv = kv[np.any(kv != 0, axis=1)]
    k2 = np.sum(kv**2, 1)
    pref = 4 * np.pi / L**3 * np.exp(-k2 / (4 * alpha**2)) / k2
    for i in range(npart):
        s = np.zeros(len(kv))
        for j in range(npart):
            s += q[j] * np.sin(kv @ (pos[i] - pos[j]))
        F[i] += C * q[i] * (pref * s) @ kv
    return F


def test4(args, out):
    print("\n=== Test 4: 3D particle force, FFT Poisson on fine grid ===")
    L = 16
    eps = 1.0 / (4 * np.pi)          # Bjerrum length 1, Coulomb F = q1 q2 / r^2
    q = np.array([1.0, -1.0])
    nplace = 6 if args.quick else 16
    rsA = [1, 2, 4] if args.quick else [1, 2, 4, 8]
    rsB = [1, 2, 4]
    rng = np.random.default_rng(1)
    print(f"  L={L}, +1/-1 pair, Peskin-4, {nplace} random placements per d")
    print("  A: kernel width = 1 fine node   B: width = r fine nodes (coarse width)")
    print("   d  opt  r   mean|dF|/|F|   std|dF|/|F|   max|F1+F2|/|F|")
    table = {}
    for dsep in (2.0, 4.0):
        placements = []
        for _ in range(nplace):
            p1 = 0.5 + L * rng.random(3)
            nh = rng.normal(size=3)
            nh /= np.linalg.norm(nh)
            p2 = p1 + dsep * nh
            pos = np.array([p1, p2])
            placements.append((pos, ewald_forces(pos, q, L, eps)))
        for opt, rs in (("A", rsA), ("B", rsB)):
            for r in rs:
                width = 1 if opt == "A" else r
                errs, mom = [], []
                for pos, Fe in placements:
                    F = pm_forces(pos, q, L, r, eps, width)
                    nF = np.linalg.norm(Fe[0])
                    errs.append(np.linalg.norm(F[0] - Fe[0]) / nF)
                    mom.append(np.linalg.norm(F[0] + F[1]) / nF)
                table[(dsep, opt, r)] = (np.mean(errs), np.std(errs))
                print(f"  {dsep:.0f}   {opt}  {r:2d}   {np.mean(errs):.3e}     "
                      f"{np.std(errs):.3e}     {np.max(mom):.1e}")

    fig, ax = plt.subplots(figsize=(6, 4))
    for dsep, ls in ((2.0, "-"), (4.0, "--")):
        for opt, rs in (("A", rsA), ("B", rsB)):
            ax.errorbar(rs, [table[(dsep, opt, r)][0] for r in rs],
                        yerr=[table[(dsep, opt, r)][1] for r in rs],
                        ls=ls, marker="o", capsize=3,
                        label=f"d={dsep:.0f}, {opt}")
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xlabel("r")
    ax.set_ylabel("|F_PM - F_Ewald| / |F_Ewald|")
    ax.legend()
    ax.set_title("Test 4: pair force vs refinement")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "test4_pair_force.png"), dpi=130)
    plt.close(fig)


# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quick", action="store_true", help="smaller runs")
    ap.add_argument("--outdir", default=os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "out"))
    ap.add_argument("--tests", default="1234")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    t1 = None
    if "1" in args.tests or "2" in args.tests:
        t1 = test1(args, args.outdir)
    if "2" in args.tests:
        test2(args, args.outdir, t1)
    if "3" in args.tests:
        test3(args, args.outdir)
    if "4" in args.tests:
        test4(args, args.outdir)
    print(f"\nFigures in {args.outdir}")


if __name__ == "__main__":
    main()
