#!/usr/bin/env python3
"""Prueba C: campo sobre la particula entre placas a tierra, con el kernel Hann-n, por SUMA DE IMAGENES
puntuales + correccion de hexadecapolo (cumulante kappa4 del kernel). Formula analitica para comparar con
la simulacion (ver README, 'Formula analitica para el kernel Hann').

  phi(r) = l_B q [ 1/r + (kappa4/24) sum_i d_i^4 (1/r) + O(r^-7) ]        (fuera del soporte)
  kappa4 = m4 - 3 m2^2 = n^4 (3/(4 pi^4) - 1/120),  m2 = n^2 (1/12 - 1/(2 pi^2))
  E_z(d) = -l_B q sum_{(i,j,n)} s d_Z [ 1/R + kappa4 (35 (X^4+Y^4+Z^4) - 21 R^4)/(4 R^9) ]
con imagenes +q en (i l, j l, d+2nL) (salvo la propia) y -q en (i l, j l, -d+2nL); Z = d - z_imagen.
Uso: python3 formula_imagenes_kernel.py [n=8]
"""
import sys
import numpy as np
L = 32.0; ELL = 32.0; LB = 10.0 / (4 * np.pi)

def kappa4(n):
    return n ** 4 * (3 / (4 * np.pi ** 4) - 1 / 120)

def Ez(d, n=8, corr=True, NL=6, NN=400):
    k4 = kappa4(n); idx = np.arange(-NL, NL + 1); nn = np.arange(-NN, NN + 1)
    X, Y, N = np.meshgrid(idx * ELL, idx * ELL, nn, indexing="ij")
    tot = 0.0
    for s, z0 in ((+1, d + 2 * L * N), (-1, -d + 2 * L * N)):
        Z = d - z0; R = np.sqrt(X ** 2 + Y ** 2 + Z ** 2); m = R > 1e-9
        R = np.where(m, R, 1.0); Zm = np.where(m, Z, 0.0)
        dz = -Zm / R ** 3
        if corr:
            Nn = 35 * (X ** 4 + Y ** 4 + Zm ** 4) - 21 * R ** 4
            dz = dz + k4 / 4 * ((140 * Zm ** 3 - 84 * R ** 2 * Zm) / R ** 9 - 9 * Zm * Nn / R ** 11)
        tot += s * np.where(m, dz, 0.0).sum()
    return -LB * tot

if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    print("kappa4(n=%d) = %.4f" % (n, kappa4(n)))
    print("  d    E_z puntual    E_z + hexadecapolo")
    for d in (n / 2 - 0.5, n / 2, 5, 6, 8, 10):
        print("%5.1f  %+.6e  %+.6e" % (d, Ez(d, n, False), Ez(d, n, True)))
