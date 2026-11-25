#!/usr/bin/env python3
"""Check of the hypothesis "flla advects the ions with u + F/(2 rho), not with u" (2026-11-16, user's request).

Background: with electrokinetics_advection flla the Delta mu floor does not decrease with viscosity
(grafico_dmu_u_visc_3adv.py). flla takes the link mass flux from the POST-collision distributions f*
(nernst_planck_lattice_flux_compute, src/nernst_planck.c), whose momentum is
    j* = sum_a f*_a c_a = j_pre + F = rho u + F/2      (collision.c: mode j += F; u = (j_pre + F/2)/rho),
while centred advection uses u. F/2 does not depend on the viscosity.

Exact reconstruction, no steady-state assumption: dist-<step> is written after propagation, so
f*_a(x) = f_dump_a(x + c_a) (propagation moves f*_a(x) to x + c_a). vel-<step> is u of the collision of the
same step. Then, per node:
    F/2 = j* - rho u                       (half the body force of that step)
    m_a(x) = f*_abar(x+c_a)/rho*(x+c_a) - f*_a(x)/rho*(x)   (flla mode 1, np_lattice_mflux_kernel)
    S(x) = sum_a m_a(x)   = net flla mass-fraction flux into x; flla changes a UNIFORM ion field n0 by
                            n0 * S per step (mode 1 conserves a uniform field only up to LB density differences)
Printed per viscosity (flla runs; the flow is identical in all variants): RMS and max of |u|, |F/(2 rho)|,
|j*/rho|, the RMS of the discrete divergence (central differences) of rho u and of F/2, RMS and max |S|,
and (cross-check) max|u| at step 1, where j_pre = 0 so u = F/(2 rho) exactly.
Reads files only. Usage: python3 scripts/comprobar_flla_fuerza.py [step]   (default 3000)
"""
import os, sys
import numpy as np

BASE = "/home/danolv/Sim/ludwig/microgel/local/TEST-HIDRO-dx2"
RUNS = ("hidro_si_eta0.02", "hidro_si", "hidro_si_eta0.1", "hidro_si_eta0.2", "hidro_si_eta0.5",
        "hidro_si_eta1.0", "hidro_si_eta2.0")
ETA = (0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0)
N = 24
CV = np.array([(0, 0, 0), (1, 1, 0), (1, 0, 1), (1, 0, 0), (1, 0, -1), (1, -1, 0), (0, 1, 1), (0, 1, 0), (0, 1, -1),
               (0, 0, 1), (0, 0, -1), (0, -1, 1), (0, -1, 0), (0, -1, -1), (-1, 1, 0), (-1, 0, 1), (-1, 0, 0),
               (-1, 0, -1), (-1, -1, 0)])          # LB_CV_D3Q19, src/lb_d3q19.h
OPP = np.array([int(np.where((CV == -c).all(1))[0][0]) for c in CV])
ST = int(sys.argv[1]) if len(sys.argv) > 1 else 3000


def load_grid(path, ncol):
    a = np.loadtxt(path).reshape(-1)          # dist: one value per line, 19 consecutive per site
    a = a.reshape(N ** 3, ncol)
    assert a.shape == (N ** 3, ncol), (path, a.shape)
    return a.reshape(N, N, N, ncol)                # x slow, z fast (as the other scripts of this series)


def shift(a, c):
    """b[x] = a[x + c] (periodic)."""
    return np.roll(a, shift=(-c[0], -c[1], -c[2]), axis=(0, 1, 2))


def div(v):
    return sum((np.roll(v[..., k], -1, axis=k) - np.roll(v[..., k], 1, axis=k)) / 2 for k in range(3))


rms = lambda x: np.sqrt(np.mean(x ** 2))
print("step %d (flla runs; flow identical in old/centred/flla)" % ST)
print(" eta  | RMS|u|     RMS|F/2rho| RMS|j*/rho| max|u|     max|F/2rho| | RMS div(rho u) RMS div(F/2) |"
      " RMS|S|     max|S|    | max|u| step 1")
for e, run in zip(ETA, RUNS):
    p = os.path.join(BASE, "advfix_flla", run)
    fd = load_grid(os.path.join(p, "dist-%09d.001-001" % ST), 19)
    u = load_grid(os.path.join(p, "vel-%09d.001-001" % ST), 3)
    fs = np.stack([shift(fd[..., a], CV[a]) for a in range(19)], axis=-1)       # post-collision f*
    rho = fs.sum(-1)
    js = np.einsum("xyza,ab->xyzb", fs, CV.astype(float))
    half_f = js - rho[..., None] * u
    S = np.zeros_like(rho)
    for a in range(1, 19):
        S += shift(fs[..., OPP[a]], CV[a]) / shift(rho, CV[a]) - fs[..., a] / rho
    u1 = load_grid(os.path.join(p, "vel-%09d.001-001" % 1), 3)
    mag = lambda v: np.sqrt((v ** 2).sum(-1))
    print("%5.2f | %.3e  %.3e  %.3e  %.3e  %.3e  | %.3e      %.3e     | %.3e  %.3e | %.3e" % (
        e, rms(mag(u)), rms(mag(half_f / rho[..., None])), rms(mag(js / rho[..., None])), mag(u).max(),
        mag(half_f / rho[..., None]).max(), rms(div(rho[..., None] * u)), rms(div(half_f)), rms(S), np.abs(S).max(),
        mag(u1).max()))
