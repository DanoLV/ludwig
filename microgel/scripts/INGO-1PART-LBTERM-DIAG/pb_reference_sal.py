"""Discrete Poisson-Boltzmann reference WITH SALT (two species) for LONG-SALT-LD2LB-CEN-3E4 (2026-11-16).

Copy of pb_reference.py (counterions only) extended to a 1:1 electrolyte; same discretisation as the Ludwig run,
so it is the fixed point of Ludwig's NP + Poisson for a particle held at r0:
  * particle charge q = +1, separable Hann kernel of order 8 (subgrid.c), nodes at integer coords 1..L (periodic);
  * Poisson with the FFT solver's 27-point Laplacian, lambda(k) = 6 sum_p w_p (1 - cos k.c_p), phi = beta e psi:
        lambda(k) phi_hat(k) = 4 pi lB rho_hat(k),  k != 0,  rho = q_p + n_+ - n_-;
  * zero Capuani flux <=> n_k,i proportional to exp(-z_k phi_i): n_+ ~ exp(-phi), n_- ~ exp(+phi);
  * CANONICAL, one constraint per species (NP conserves each total): N_+ and N_- read from the run's step-1 qsi
        n_+,i = N_+ exp(-phi_i)/sum_j exp(-phi_j),   n_-,i = N_- exp(phi_i)/sum_j exp(phi_j).
Damped Picard iteration (n <- 0.5 n + 0.5 n_new), stop when max relative change < 1e-12 (same as pb_reference.py).

Usage (functions are imported by cloud_profile_sal.py):
    python3 pb_reference_sal.py <run dir>      -> solves at the INITIAL particle position (step 1) and prints
                                                 the radial profile n_k(r)/nbar_k.
"""
import sys
import numpy as np

L, LB, NORD = 27, 2.84, 8.0
EDGES = np.arange(0.0, 14.0, 1.0)

# discrete Laplacian eigenvalues (D3Q27 weights)
_k = 2 * np.pi * np.fft.fftfreq(L)
_KX, _KY, _KZ = np.meshgrid(_k, _k, _k, indexing="ij")
LAM = np.zeros_like(_KX)
for cx in (-1, 0, 1):
    for cy in (-1, 0, 1):
        for cz in (-1, 0, 1):
            s = abs(cx) + abs(cy) + abs(cz)
            if s:
                LAM += 6 * {1: 2 / 27, 2: 1 / 54, 3: 1 / 216}[s] * (1 - np.cos(_KX * cx + _KY * cy + _KZ * cz))
LAM[0, 0, 0] = 1.0


def hann1d(center):
    x = np.arange(1, L + 1, dtype=float)
    d = x - center
    d -= L * np.round(d / L)
    return np.where(np.abs(d) <= NORD / 2, (1.0 + np.cos(2 * np.pi * d / NORD)) / NORD, 0.0)


def particle_charge(r0, q=1.0):
    wx, wy, wz = (hann1d(c) for c in r0)
    return q * np.einsum("i,j,k->ijk", wx, wy, wz)


def solve_phi(rho):
    rh = np.fft.fftn(4 * np.pi * LB * rho) / LAM
    rh[0, 0, 0] = 0.0
    return np.real(np.fft.ifftn(rh))


def solve_pb(r0, Np, Nm, tol=1e-12, itmax=5000):
    qp = particle_charge(r0)
    npl = np.full((L, L, L), Np / L ** 3)
    nmi = np.full((L, L, L), Nm / L ** 3)
    for it in range(itmax):
        phi = solve_phi(qp + npl - nmi)
        ep = np.exp(-(phi - phi.min()))
        em = np.exp(phi - phi.max())
        np_new, nm_new = Np * ep / ep.sum(), Nm * em / em.sum()
        err = max(np.abs(np_new - npl).max() / npl.mean(), np.abs(nm_new - nmi).max() / nmi.mean())
        npl, nmi = 0.5 * npl + 0.5 * np_new, 0.5 * nmi + 0.5 * nm_new
        if err < tol:
            break
    return npl, nmi, it + 1, err


def distances(r0):
    g = np.arange(1, L + 1, dtype=float)
    X, Y, Z = np.meshgrid(g, g, g, indexing="ij")
    dv = np.stack([X - r0[0], Y - r0[1], Z - r0[2]], -1)
    dv -= L * np.round(dv / L)
    return np.linalg.norm(dv, axis=-1)


def radial(n, dist):
    idx = np.digitize(dist, EDGES) - 1
    return np.array([n[idx == b].mean() for b in range(len(EDGES) - 1)])


def totals_from_run(d):
    q = np.loadtxt(f"{d}/qsi-000000001.001-001")
    return q[:, 0].sum(), q[:, 1].sum()


def position(d, s):
    return np.array([float(v) for v in open(f"{d}/config.cds{s:08d}.001-001").read().split("\n")[35].split()])


if __name__ == "__main__":
    d = sys.argv[1]
    Np, Nm = totals_from_run(d)
    r0 = position(d, 1)
    npl, nmi, it, err = solve_pb(r0, Np, Nm)
    print(f"N+ = {Np:.8f}, N- = {Nm:.8f}, r0 = {r0}; Picard iterations {it}, last rel change {err:.1e}")
    dist = distances(r0)
    pp, pm = radial(npl / (Np / L ** 3), dist), radial(nmi / (Nm / L ** 3), dist)
    print("PB profile n_k(r)/nbar_k at r0:")
    for b in range(len(EDGES) - 1):
        print(f"   r in [{EDGES[b]:4.1f},{EDGES[b+1]:4.1f}):  n+ {pp[b]:.4f}   n- {pm[b]:.4f}")
