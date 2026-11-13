"""Discrete Poisson-Boltzmann equilibrium reference for INGO-1PART-LBTERM-DIAG
(2026-11-12). Counterions only, canonical (fixed total), same discretisation as
the Ludwig run, so it is the exact fixed point of Ludwig's NP + Poisson:

  * Particle charge q = +1 spread with the separable Hann kernel (subgrid.c, p = 2):
        w(r) = (1/n) [1 + cos(2 pi r/n)],  |r| <= n/2,  n = 8,
    on nodes at integer coordinates 1..L (periodic), particle at r0.
  * Poisson with the FFT solver's discrete 27-point Laplacian (psi_fft.c,
    stencil_d3q27.c after the 20251203 fix): lambda(k) = 6 sum_p w_p (1 - cos k.c_p),
    and phi = beta e psi (psi_fft.c multiplies rho by beta):
        lambda(k) phi_hat(k) = 4 pi lB rho_hat(k),   k != 0;  phi_hat(0) = 0.
  * NP zero-flux condition of Ludwig's link flux (Capuani form) is exactly
    n_i proportional to exp(-z phi_i); counterions z = -1, total N = 1:
        n_i = N exp(phi_i) / sum_j exp(phi_j).
Self-consistent solution by damped Picard iteration.
Output: radial profile (same bins as cloud_profile.py) and node-by-node comparison
with the Ludwig dumps given on the command line (NOFLUCT runs, particle fixed).
"""
import sys
import numpy as np

L, LB, NORD, NTOT = 27, 2.84, 8.0, 1.0
R0 = np.array([14.085456533601072, 14.302999340732786, 13.56327660187159])  # posicion_inicial.txt
EDGES = np.arange(0.0, 14.0, 1.0)

# --- particle charge on the grid (array index a <-> coordinate a+1) ---
def hann1d(center):
    x = np.arange(1, L + 1, dtype=float)
    d = x - center
    d -= L*np.round(d/L)
    w = np.where(np.abs(d) <= NORD/2, (1.0 + np.cos(2*np.pi*d/NORD))/NORD, 0.0)
    return w
wx, wy, wz = (hann1d(c) for c in R0)
qp = np.einsum("i,j,k->ijk", wx, wy, wz)          # [x, y, z]
print("particle charge on grid:", qp.sum())

# --- discrete Laplacian eigenvalues (D3Q27 weights) ---
k = 2*np.pi*np.fft.fftfreq(L)
KX, KY, KZ = np.meshgrid(k, k, k, indexing="ij")
lam = np.zeros_like(KX)
for cx in (-1, 0, 1):
    for cy in (-1, 0, 1):
        for cz in (-1, 0, 1):
            s = abs(cx) + abs(cy) + abs(cz)
            if s == 0:
                continue
            wp = {1: 2/27, 2: 1/54, 3: 1/216}[s]
            lam += 6*wp*(1 - np.cos(KX*cx + KY*cy + KZ*cz))
lam[0, 0, 0] = 1.0

def solve_phi(rho):
    rh = np.fft.fftn(4*np.pi*LB*rho)/lam
    rh[0, 0, 0] = 0.0
    return np.real(np.fft.ifftn(rh))

n = np.full((L, L, L), NTOT/L**3)
for it in range(2000):
    phi = solve_phi(qp - n)
    e = np.exp(phi - phi.max())
    nnew = NTOT*e/e.sum()
    err = np.abs(nnew - n).max()/n.mean()
    n = 0.5*n + 0.5*nnew
    if err < 1e-12:
        break
print(f"Picard iterations {it+1}, last rel change {err:.1e}")
nbar = NTOT/L**3
np.save("pb_reference_n.npy", n)

# --- radial profile ---
g = np.arange(1, L + 1, dtype=float)
X, Y, Z = np.meshgrid(g, g, g, indexing="ij")
dv = np.stack([X - R0[0], Y - R0[1], Z - R0[2]], -1)
dv -= L*np.round(dv/L)
dist = np.linalg.norm(dv, axis=-1)
idx = np.digitize(dist, EDGES) - 1
prof = np.array([(n/nbar)[idx == b].mean() for b in range(len(EDGES) - 1)])
print(f"PB excess r<3: {(n/nbar)[dist<3].mean()-1:+.4f};  Nex(r<6) = {(n-nbar)[dist<6].sum():+.4f} e")
print("profile n(r)/nbar (PB):")
for b in range(len(EDGES) - 1):
    print(f"   r in [{EDGES[b]:4.1f},{EDGES[b+1]:4.1f}):  {prof[b]:.4f}")

# --- node-by-node comparison with Ludwig dumps (ascii, x slowest, z fastest) ---
for f in sys.argv[1:]:
    q = np.loadtxt(f)[:, 1].reshape(L, L, L)
    rel = q/n - 1
    print(f"{f}: node-by-node rel. diff vs PB: max |.| {np.abs(rel).max():.4f}, rms {np.sqrt((rel**2).mean()):.4f}; "
          f"r<3: mean {rel[dist<3].mean():+.4f}; total ratio {q.sum()/n.sum():.6f}")
