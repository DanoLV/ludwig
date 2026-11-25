#!/usr/bin/env python3
"""Exact decomposition of the flla (mode 1) uniform-field source S, and its correlations (2026-11-16, user's request).

Follows comprobar_flla_fuerza.py (same reconstruction of the post-collision f* from dist-<step>, reused from it).
    S_i = D_i/rho_i + R_i,
    D_i = sum_a [f*_abar(x+c_a) - f*_a(x)]                 net LB mass flux into x (= rho(t+1) - rho*)
    R_i = sum_a f*_abar(x+c_a) (1/rho(x+c_a) - 1/rho(x))   part due to LB density differences between neighbours
Prints per viscosity: RMS of S, D/rho, R, corr(S, R), RMS and max of rho - <rho>; and for eta = 0.02 and 2 the
node-by-node correlation of S with -div(F/2), -div(rho u), -div(j*) (central differences; F/2 = j* - rho u).
Outcome 2026-11-16: D ~ 0, S = R exactly -> the hypothesis "flla advects with u + F/2rho" was REJECTED; the source
comes from the normalisation by the local LB density (documented limitation of flla in src/nernst_planck.c).
Reads files only. Usage: python3 scripts/comprobar_flla_descomposicion.py
"""
import os
import numpy as np

ns = {}
exec(open("/home/danolv/Sim/ludwig/microgel/local/TEST-HIDRO-dx2/scripts/comprobar_flla_fuerza.py").read()
     .split("rms = lambda")[0], ns)
rms = lambda x: np.sqrt(np.mean(x ** 2))


def fields(run):
    p = os.path.join(ns["BASE"], "advfix_flla", run)
    fd = ns["load_grid"](p + "/dist-000003000.001-001", 19)
    u = ns["load_grid"](p + "/vel-000003000.001-001", 3)
    fs = np.stack([ns["shift"](fd[..., a], ns["CV"][a]) for a in range(19)], -1)
    rho = fs.sum(-1)
    S, D, R = (np.zeros_like(rho) for _ in range(3))
    for a in range(1, 19):
        fin = ns["shift"](fs[..., ns["OPP"][a]], ns["CV"][a]); rj = ns["shift"](rho, ns["CV"][a])
        S += fin / rj - fs[..., a] / rho; D += fin - fs[..., a]; R += fin * (1 / rj - 1 / rho)
    assert np.allclose(S, D / rho + R, rtol=0, atol=1e-15)
    js = np.einsum("xyza,ab->xyzb", fs, ns["CV"].astype(float))
    return rho, u, js, S, D, R


print(" eta | RMS S      RMS D/rho  RMS R      | corr(S,R)  | RMS(rho-<rho>)  max|rho-<rho>|")
keep = {}
for e, run in zip(ns["ETA"], ns["RUNS"]):
    rho, u, js, S, D, R = fields(run)
    keep[e] = (rho, u, js, S)
    dr = rho - rho.mean()
    print("%4.2f | %.3e  %.3e  %.3e  | %+.4f    | %.3e       %.3e" % (e, rms(S), rms(D / rho), rms(R),
          np.corrcoef(S.ravel(), R.ravel())[0, 1], rms(dr), np.abs(dr).max()))
print()
for e in (0.02, 2.0):
    rho, u, js, S = keep[e]
    hf = js - rho[..., None] * u
    for nm, v in (("-div(F/2)", -ns["div"](hf)), ("-div(rho u)", -ns["div"](rho[..., None] * u)), ("-div(j*)", -ns["div"](js))):
        print("eta %4.2f: corr(S, %-11s) = %+.4f   slope S/%s = %+.3f" % (e, nm, np.corrcoef(S.ravel(), v.ravel())[0, 1],
              nm, (S * v).sum() / (v * v).sum()))
