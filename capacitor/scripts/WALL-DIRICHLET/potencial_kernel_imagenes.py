#!/usr/bin/env python3
"""Potencial sobre la linea de la particula (rho = 0) de la carga q con kernel Hann-n entre placas a tierra
(placas infinitas, sin imagenes laterales): potencial LIBRE del kernel (representacion de una integral 1D) +
imagenes con la correccion de hexadecapolo.  Placas en z = 0 y z = L, carga en z = d.

 phi_libre(r) = l_B q (2/sqrt(pi)) int_0^inf g(x;t) g(y;t) g(z;t) dt,   [1/r = (2/sqrt(pi)) int_0^inf e^{-r^2 t^2} dt]
 g(x;t) = int w(x') e^{-t^2 (x-x')^2} dx'  (cerrada con erf compleja, kappa = 2 pi/n):
   g = (1/n)(sqrt(pi)/(2t)) { erf(t(n/2-x)) + erf(t(n/2+x))
         + Re[ e^{i kappa x - kappa^2/(4t^2)} ( erf(t(n/2-x) - i kappa/(2t)) + erf(t(n/2+x) + i kappa/(2t)) ) ] }
 phi_img(z) = l_B q sum_j s_j [ 1/|z-z_j| + (7/4) kappa4 / |z-z_j|^5 ]     (imagenes, sobre el eje)
Uso: python3 potencial_kernel_imagenes.py -> tabla contra la teoria de Fourier con el kernel real (ell = 200)
"""
import sys, os
import numpy as np
from scipy.special import erf as _erf                  # complex-argument erf (scipy.special)
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
L = 32.0; LB = 10.0 / (4 * np.pi); NH = 8
K4 = NH ** 4 * (3 / (4 * np.pi ** 4) - 1 / 120)

def w(x, n=NH):
    x = np.asarray(x, float); return np.where(np.abs(x) < n / 2, (1 + np.cos(2 * np.pi * x / n)) / n, 0.0)

xq, wq = np.polynomial.legendre.leggauss(160)
def g_quad(x, t, n=NH):                      # integral directa (comprobacion de la forma cerrada)
    xp = 0.5 * n * xq; wp = 0.5 * n * wq
    return np.sum(wp * w(xp, n) * np.exp(-t ** 2 * (x - xp) ** 2))

def g_closed(x, t, n=NH):
    k = 2 * np.pi / n
    a = _erf(t * (n / 2 - x)) + _erf(t * (n / 2 + x))
    b = np.exp(1j * k * x - k ** 2 / (4 * t ** 2)) * (_erf(t * (n / 2 - x) - 1j * k / (2 * t)) +
                                                        _erf(t * (n / 2 + x) + 1j * k / (2 * t)))
    return float(np.sqrt(np.pi) / (2 * t * n) * (a + b.real))

def phi_libre(rho, z, n=NH):
    # integral en t por cuadratura (log) de g(x)g(y)g(z), con x = rho, y = 0 (eje) y z
    tt, tw = np.polynomial.legendre.leggauss(400)
    u = 0.5 * (tt + 1) * 12.0                # t en (0, 12)
    wu = 0.5 * 12.0 * tw
    s = 0.0
    for ti, wi in zip(u, wu):
        s += wi * g_quad(rho, ti, n) * g_quad(0.0, ti, n) * g_quad(z, ti, n)
    # cola t > 12: g ~ sqrt(pi) w(x)/t -> producto ~ pi^{3/2} w w w / t^3 ; contribuye pi^{3/2} wwww/(2*144)
    s += np.pi ** 1.5 * w(rho, n) * w(0.0, n) * w(z, n) / (2 * 12.0 ** 2)
    return LB * 2 / np.sqrt(np.pi) * s

def phi_img(z, d, NN=4000):
    n = np.arange(-NN, NN + 1)
    zp = d + 2 * L * n; zm = -d + 2 * L * n
    rp = np.abs(z - zp); rm = np.abs(z - zm)
    mp_ = (rp > 1e-9) & (n != 0)          # n = 0 positiva es la carga real: su potencial es phi_libre
    pos = np.where(mp_, 1 / np.where(mp_, rp, 1) + 1.75 * K4 / np.where(mp_, rp, 1) ** 5, 0.0)
    neg = 1 / rm + 1.75 * K4 / rm ** 5
    return LB * (pos.sum() - neg.sum())

if __name__ == "__main__":
    sys.argv = ["x", os.path.join(ROOT, "C-sin-fondo-hann08"), "8"]; sys.path.insert(0, HERE)
    import grafico_C_campo_kernel as K
    P = K.modes(200.0)
    print("comprobacion de la forma cerrada de g (x=1.3, t=0.7):", g_closed(1.3, 0.7), g_quad(1.3, 0.7))
    d = 8.0
    print("   z    phi_libre   phi_imagenes   suma      Fourier(kernel, ell=200)   dif")
    for z in (0.5, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 16.0, 24.0):
        pl = phi_libre(0.0, z - d); pi_ = phi_img(z, d); ref = K.psi_q(np.array([z]), d, P)[0]
        print("%5.1f  %.6f   %+.6f    %.6f   %.6f              %+.1e" % (z, pl, pi_, pl + pi_, ref, pl + pi_ - ref))
