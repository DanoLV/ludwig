#!/usr/bin/env python3
"""
Analiza la prueba B (capacitor con electrolito, lambda_D = 4, estacionario):
  1) residuo de Poisson nodo a nodo con la regla de pared (ec. 15)
  2) equilibrio de Boltzmann de los iones: ln(rho+/rho-) = -2 psi
  3) perfil vs Debye-Hückel discreta (misma red) y continua (ec. 17)
  4) campo en el nodo junto a la pared vs derivada de DH continua
Uso: analizar_B.py [directorio de B] [paso]
"""
import sys, numpy as np
d = sys.argv[1] if len(sys.argv) > 1 else "../B-debye-huckel"
step = int(sys.argv[2]) if len(sys.argv) > 2 else 60000
nx, ny, nz = 4, 4, 34
eps, beta, e, rhos = 1.0e4, 1.0e5, 1.0, 3.125e-3
plo, phi_hi = 0.0, 1.0e-6
kap = np.sqrt(8*np.pi*(beta*e**2/(4*np.pi*eps))*rhos); L = 32.0
psi = np.loadtxt("%s/psi-%09d.001-001" % (d, step)).reshape(nx, ny, nz)
q = np.loadtxt("%s/qsi-%09d.001-001" % (d, step)).reshape(nx, ny, nz, 2)
E = np.loadtxt("%s/efield-%09d.001-001" % (d, step)).reshape(nx, ny, nz, 3)
p = psi[0,0,:]; rp = q[0,0,:,0]; rm = q[0,0,:,1]; f = slice(1, nz-1)
lap = np.zeros(nz)
for k in range(1, nz-1):
    lap[k] = (p[k+1]-p[k])*(2 if k+1 == nz-1 else 1) + (p[k-1]-p[k])*(2 if k-1 == 0 else 1)
res = eps*lap[f] + e*beta*e*(rp-rm)[f]
print("1) residuo Poisson max %.2e (escala %.2e)" % (abs(res).max(), abs(eps*lap[f]).max()))
y = np.log(rp[f]/rm[f]); s, c = np.polyfit(p[f], y, 1)
print("2) ln(rho+/rho-) = %.6f psi + %.1e  (Boltzmann exacto: -2)" % (s, c))
zc = np.arange(1, nz+1) - 17.5
cont = 0.5*(plo+phi_hi) + 0.5*(phi_hi-plo)*np.sinh(kap*zc)/np.sinh(kap*L/2)
n = nz-2; M = np.zeros((n, n)); r = np.zeros(n); pb = 0.5*(plo+phi_hi)
for i in range(n):
    M[i, i] = -2-kap**2; r[i] = -kap**2*pb
    if i == 0: M[i, i] -= 1; r[i] -= 2*plo
    else: M[i, i-1] = 1
    if i == n-1: M[i, i] -= 1; r[i] -= 2*phi_hi
    else: M[i, i+1] = 1
disc = np.linalg.solve(M, r); dps = phi_hi-plo
print("3) perfil: vs DH discreta %.2e | vs DH continua %.2e (relativo a dpsi)"
      % (abs(p[f]-disc).max()/dps, abs(p[f]-cont[f]).max()/dps))
dcont = lambda z: 0.5*dps*kap*np.cosh(kap*(z-17.5))/np.sinh(kap*L/2)
Ez = -E[0,0,:,2]/1e-5                      # efield sale como -grad(psi)/beta
print("4) campo junto a la pared (z=2): %.5f del continuo | z=33: %.5f | en pared: %.1e"
      % (Ez[1]/dcont(2.0), Ez[-2]/dcont(33.0), abs(E[...,[0,-1],:]).max()))
