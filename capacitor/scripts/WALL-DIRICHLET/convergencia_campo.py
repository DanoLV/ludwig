#!/usr/bin/env python3
"""
Orden de convergencia del campo en el nodo junto a la pared, sobre la solucion
de Debye-Hückel discreta (la misma que Ludwig reproduce a 1e-8), con kappa*L=8.
Compara: diferencia central actual, factor 2 en diferencias (implementado) e
interpolacion cuadratica (descartada).
"""
import numpy as np
def grads(p0, p1, ps, h):
    return {"central (antes)": (p1-ps)/(2*h),
            "factor 2 (ahora)": ((p1-p0)+2*(p0-ps))/(2*h),
            "cuadratica": (p1+3*p0-4*ps)/(3*h)}
def disc(n, kh):
    M = np.zeros((n, n)); r = np.zeros(n)
    for i in range(n):
        M[i, i] = -2-kh**2; r[i] = -kh**2*0.5
        if i == 0: M[i, i] -= 1
        else: M[i, i-1] = 1
        if i == n-1: M[i, i] -= 1; r[i] -= 2.0
        else: M[i, i+1] = 1
    return np.linalg.solve(M, r)
def disc_rule(n, kh, rule):
    """DH discreta con la pared en el plano medio ("pared", factor 2) o en el
    nodo solido ("ingenuo", factor 1)."""
    M = np.zeros((n, n)); r = np.zeros(n); w = 2 if rule == "pared" else 1
    for i in range(n):
        M[i, i] = -2-kh**2; r[i] = -kh**2*0.5
        if i == 0: M[i, i] -= (w-1)
        else: M[i, i-1] = 1
        if i == n-1: M[i, i] -= (w-1); r[i] -= w*1.0
        else: M[i, i+1] = 1
    return np.linalg.solve(M, r)


def convergence(ns=(32, 64, 128, 256, 512)):
    """Errores relativos del campo en el nodo junto a la pared y del perfil,
    kappa = 0.25 y L = 32 fijos (unidades fisicas), h = 32/n."""
    field = {}; prof = {"pared": [], "ingenuo": []}; hs = []
    for n in ns:
        h = 32.0/n; hs.append(h); phi = disc(n, 0.25*h)
        ex = 0.5*0.25*np.cosh(0.25*(0.5*h-16))/np.sinh(0.25*16)
        for k, v in grads(phi[0], phi[1], 0.0, h).items():
            field.setdefault(k, []).append(abs(v/ex-1))
        for rule in prof:
            ph = disc_rule(n, 0.25*h, rule)
            z = (np.arange(n)+0.5)*h - 16
            cont = 0.5 + 0.5*np.sinh(0.25*z)/np.sinh(0.25*16)
            prof[rule].append(np.abs(ph-cont).max())
    return np.array(hs), field, prof


def main():
    hs, rows, prof = convergence()
    for k, e in rows.items():
        print("%-17s error: %s  razones: %s" % (k, " ".join("%.1e" % x for x in e),
              " ".join("%.2f" % (e[i]/e[i+1]) for i in range(len(e)-1))))
    for k, e in prof.items():
        print("perfil %-10s error: %s  razones: %s" % (k, " ".join("%.1e" % x for x in e),
              " ".join("%.2f" % (e[i]/e[i+1]) for i in range(len(e)-1))))


if __name__ == "__main__":
    main()
