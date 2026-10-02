"""
Capacitor de placas: electrodos solidos (ns capas) en ambos extremos de z,
fluido en el medio (nf capas), periodico en x,y. Stencil D3Q27, laplaciano
nabla^2 phi = (2/cs^2) sum_i w_i [phi(r+c_i) - phi(r)] * factor_i.
Se resuelve SOLO en nodos fluidos (Dirichlet en solidos, columnas eliminadas),
y se comparan formulas de carga inducida por electrodo.
Unidades: nabla^2 phi = -4 pi lB rho_el ; tomamos 4 pi lB = 1.
"""
import itertools, numpy as np
import scipy.sparse as sp, scipy.sparse.linalg as spla

cs2 = 1.0/3.0
cv = [c for c in itertools.product((-1,0,1), repeat=3)]
w = {0: 8/27, 1: 2/27, 2: 1/54, 3: 1/216}
wv = [w[sum(abs(x) for x in c)] for c in cv]

def solve(nf, ns, psi_lo, psi_hi, rho_fn, mode, nx=2):
    nz = nf + 2*ns
    chi = np.zeros(nz, int); chi[:ns] = 1; chi[ns+nf:] = 1
    psis = np.zeros(nz); psis[:ns] = psi_lo; psis[ns+nf:] = psi_hi
    idx = lambda i,j,k: (i*nx + j)*nz + k
    N = nx*nx*nz
    fluid = [n for n in range(N) if chi[n % nz] == 0]
    fmap = {n: m for m, n in enumerate(fluid)}
    rows, cols, vals = [], [], []
    b = np.zeros(len(fluid))
    for n in fluid:
        k = n % nz; j = (n//nz) % nx; i = n//(nz*nx)
        m = fmap[n]
        z = k - ns + 0.5            # z medido desde el plano medio izquierdo
        b[m] = -rho_fn(z)           # sum ... = -rho  (4 pi lB = 1)
        diag = 0.0
        for c, wi in zip(cv, wv):
            if c == (0,0,0): continue
            kk = k + c[2]
            nb = idx((i+c[0]) % nx, (j+c[1]) % nx, kk)
            if mode == "modified":
                fac = 1 + chi[kk] - 0
            else:                   # "naive": ec. 12 con el solido fijo
                fac = 1
            coef = (2/cs2)*wi*fac
            diag -= coef
            if chi[kk]:
                b[m] -= coef*psis[kk]           # columna Dirichlet al RHS
            else:
                rows.append(m); cols.append(fmap[nb]); vals.append(coef)
        rows.append(m); cols.append(m); vals.append(diag)
    A = sp.csr_matrix((vals, (rows, cols)), shape=(len(fluid),)*2)
    x = spla.spsolve(A.tocsc(), b)
    phi = psis.copy().repeat(1)
    full = np.array([psis[n % nz] for n in range(N)], float)
    for n in fluid: full[n] = x[fmap[n]]
    asym = abs(A - A.T).max()
    return full.reshape(nx, nx, nz), chi, nz, asym

def lap_at(phi, chi, nx, nz, i, j, k, rule):
    s = 0.0
    for c, wi in zip(cv, wv):
        if c == (0,0,0): continue
        kk = k + c[2]
        if kk < 0 or kk >= nz: continue
        d = phi[(i+c[0]) % nx, (j+c[1]) % nx, kk] - phi[i,j,k]
        if rule == "eq15":     fac = 1 + chi[kk] - chi[k]
        elif rule == "eq12":   fac = 1
        elif rule == "sym":    fac = 2 if chi[kk] != chi[k] else 1
        s += (2/cs2)*wi*d*fac
    return s

def charges(phi, chi, nz, nx, ns, nf, rule):
    """Q por unidad de area, cada electrodo: Q = -(1/4pi lB) sum_elec lap."""
    lo = sum(lap_at(phi, chi, nx, nz, i, j, k, rule)
             for i in range(nx) for j in range(nx) for k in range(ns))
    hi = sum(lap_at(phi, chi, nx, nz, i, j, k, rule)
             for i in range(nx) for j in range(nx) for k in range(ns+nf, nz))
    return -lo/(nx*nx), -hi/(nx*nx)

def main():
    print("=== Prueba 1: sin iones, psi_lo=0, psi_hi=1 ===")
    nf, ns = 20, 3
    phi, chi, nz, asym = solve(nf, ns, 0.0, 1.0, lambda z: 0.0, "modified")
    prof = phi[0,0,:]
    print("asimetria de la matriz fluido-fluido (con columnas eliminadas): %.1e" % asym)
    print("perfil en los fluidos extremos: phi[1er fluido]=%.6f  phi[ultimo]=%.6f" % (prof[ns], prof[ns+nf-1]))
    print("  lineal con planos medios a distancia nf=%d  => esperado %.6f y %.6f" % (nf, 0.5/nf, 1-0.5/nf))
    sigma = 1.0/nf       # E/(4 pi lB) con 4 pi lB = 1: sigma_hi = +E, sigma_lo = -E
    print("carga exacta por area (Gauss, planos medios): lo=%+.6f  hi=%+.6f" % (-sigma, sigma))
    for rule, name in (("eq15","ec.15 literal en nodos del electrodo"),
                       ("eq12","ec.12 sin corregir (lo que propuse antes)"),
                       ("sym", "factor 2 en enlaces de frontera, simetrico")):
        lo, hi = charges(phi, chi, nz, 2, ns, nf, rule)
        print("  %-45s lo=%+.6f  hi=%+.6f   hi/exacto=%.4f" % (name, lo, hi, hi/sigma))

    print("\n=== Prueba 2: carga ionica arbitraria en el fluido, neutralidad ===")
    rng = np.random.default_rng(1)
    rz = rng.normal(size=nf)*1e-2
    rho_fn = lambda z: rz[int(z-0.5)]
    phi, chi, nz, _ = solve(nf, ns, 0.1, 0.2, rho_fn, "modified")
    Qion = rz.sum()      # por unidad de area (una celda por capa)
    for rule in ("eq15","eq12","sym"):
        lo, hi = charges(phi, chi, nz, 2, ns, nf, rule)
        print("  %-5s  Q_lo+Q_hi+Q_ion = %+.3e   (Q_ion=%+.4e)" % (rule, lo+hi+Qion, Qion))

    print("\n=== Prueba 3: orden de precision (rho uniforme, solucion parabolica) ===")
    # continuo: phi'' = -rho0 en (0,L), phi(0)=phi(L)=0 con L = distancia entre
    # planos medios => phi = rho0/2 z (L - z)
    rho0 = 1e-3
    for mode in ("naive","modified"):
        errs = []
        for nf in (10, 20, 40):
            L = nf if mode == "modified" else nf + 1   # naive: la pared "esta" en el nodo solido
            phi, chi, nz, _ = solve(nf, 3, 0.0, 0.0, lambda z: rho0, mode)
            z = np.arange(nf) + (0.5 if mode == "modified" else 1.0)
            # comparar contra el continuo con la pared FISICA en el plano medio (L=nf)
            zf = np.arange(nf) + 0.5
            exact = rho0/2*zf*(nf - zf)
            errs.append(np.abs(phi[0,0,3:3+nf] - exact).max()/exact.max())
        print("  %-9s error rel. max (nf=10,20,40): %s   razones: %.2f %.2f" %
              (mode, " ".join("%.2e" % e for e in errs), errs[0]/errs[1], errs[1]/errs[2]))


if __name__ == "__main__":
    main()
