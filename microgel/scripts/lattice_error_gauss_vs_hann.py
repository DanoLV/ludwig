import numpy as np
N = 48; L = 48.0; eps = 1.0e4; q = 1.0
sig = 0.18075603 * 8          # mismo ancho que hann8

# --- stencil D3Q27 de Ludwig: wlap = 6*wv (vecinos), wgrad = 3*wv ---
cv, wv = [], []
for a in (-1, 0, 1):
    for b in (-1, 0, 1):
        for c in (-1, 0, 1):
            n2 = a*a + b*b + c*c
            cv.append((a, b, c)); wv.append({0: 8/27, 1: 2/27, 2: 1/54, 3: 1/216}[n2])
cv = np.array(cv, float); wv = np.array(wv)
k1 = 2*np.pi*np.fft.fftfreq(N, d=L/N)
KX, KY, KZ = np.meshgrid(k1, k1, k1, indexing='ij')
lap = np.zeros_like(KX); gx = np.zeros_like(KX)
for (a, b, c), w in zip(cv, wv):
    ph = KX*a + KY*b + KZ*c
    if (a, b, c) != (0, 0, 0):
        lap += 6*w*(np.cos(ph) - 1)          # autovalor del Laplaciano discreto (<=0)
    gx += 3*w*a*np.sin(ph)                   # autovalor del gradiente discreto (x)
lap[0, 0, 0] = 1.0

def w_hann(x, n=8):
    return np.where(np.abs(x) < n/2, (1/n)*(1 + np.cos(2*np.pi*x/n)), 0.0)

def w_gauss(x, T=6):
    return np.where(np.abs(x) <= T, np.exp(-x*x/(2*sig*sig)), 0.0)

def axis_w(kern, c):
    w = kern(np.arange(N) - c)
    return w / w.sum()                       # renormaliza (el Hann ya da 1 exacto)

def lattice_force(kern, x1, x2):
    """F_x sobre la particula 2 por la 1, pipeline de Ludwig en la red."""
    wx, wy, wz = axis_w(kern, x1), axis_w(kern, 24.0), axis_w(kern, 24.0)
    rho = q * wx[:, None, None] * wy[None, :, None] * wz[None, None, :]
    psih = np.fft.fftn(rho) / (eps * (-lap)); psih[0, 0, 0] = 0
    Ex = np.real(np.fft.ifftn(-1j * gx * psih))
    g2 = axis_w(kern, x2)
    return q * np.einsum('i,j,k,ijk->', g2, wy, wz, Ex)

# --- referencia exacta en el continuo, periodica, para cada forma de nube ---
M = 40
km = 2*np.pi*np.arange(-M, M + 1)/L
QX, QY, QZ = np.meshgrid(km, km, km, indexing='ij')
Q2 = QX**2 + QY**2 + QZ**2; Q2[M, M, M] = 1.0

def ft_gauss(k):  return np.exp(-k*k*sig*sig/2)
def ft_hann(k, n=8):
    u = k*n/2
    out = np.ones_like(k)
    m = np.abs(u) > 1e-12
    s = np.sin(u[m])/u[m]
    d = 1 - (k[m]*n/(2*np.pi))**2
    out[m] = np.where(np.abs(d) > 1e-12, s/np.where(np.abs(d) > 1e-12, d, 1), 0.5)
    return out

def continuum_force(ft, r):
    S = ft(QX)*ft(QY)*ft(QZ)
    T = (S**2) * QX * np.sin(QX*r) / Q2
    T[M, M, M] = 0
    return q*q/(eps*L**3) * T.sum()

print("Error de la red contra la solucion EXACTA del continuo de la MISMA forma de nube")
print("(0 % = no hace falta correccion de corto alcance)\n")
print(f"{'r':>4} {'Gauss red vs Gauss cont':>24} {'Hann red vs Hann cont':>22} {'Hann red vs Gauss cont':>23}")
x1 = 24.0
for r in (1, 2, 3, 4, 5, 6, 8):
    fg  = lattice_force(w_gauss, x1, x1 + r)
    fh  = lattice_force(w_hann,  x1, x1 + r)
    cg  = continuum_force(ft_gauss, r)
    ch  = continuum_force(ft_hann,  r)
    print(f"{r:4d} {100*(fg/cg-1):+23.2f}% {100*(fh/ch-1):+21.2f}% {100*(fh/cg-1):+22.2f}%")
