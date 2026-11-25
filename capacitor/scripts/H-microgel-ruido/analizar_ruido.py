#!/usr/bin/env python3
"""
Análisis de las pruebas cortas H-microgel-ruido (2026-11-13): ¿el ruido térmico es demasiado fuerte?

Para cada corrida (argumentos = carpetas; por defecto las 4 de H-microgel-ruido):
 1. Desplazamiento por paso de cada monómero (posiciones de config.cds en cada paso; imagen
    mínima en x, y; z no periódico). rms por componente y máximo |dr| frente a drmax = 0.8
    (límite de Ludwig). Referencias: sqrt(2 D_desnudo) con D_desnudo = kT/(6 pi eta)(1/a_h - 1/al)
    (empujón aleatorio propio de la partícula) y sqrt(2 kT/(6 pi eta a_h)) (difusión total de Stokes).
 2. Equipartición en los enlaces: U = sum_enlaces k/2 (r - r0)^2. Para fluctuaciones chicas la red
    es armónica y <U> = (rango de la matriz de rigidez)/2 kT = 35/2 kT (rango 35 según el README de
    G-microgel-relax). Teorema de equipartición; la aproximación armónica requiere sqrt(kT/k) << r0.
 3. Fluido: <u_a^2> = kT/rho por componente (equipartición de las fluctuaciones del LB; Dünweg &
    Ladd, Adv. Polym. Sci. 221, 89 (2009)), sobre nodos de fluido z = 2..33 de los vel-*.
 4. Distancia mínima de cada monómero a los planos sólido/fluido (z = 1.5 y 33.5): d y h = d - 3.5
    (radio del kernel Hann-8: con h < 0 el kernel toca nodos de pared).
Solo lee archivos; imprime un resumen.

Uso:  python3 capacitor/scripts/H-microgel-ruido/analizar_ruido.py [carpetas...]
"""
import glob
import math
import os
import sys

import numpy as np

H = "/home/danolv/Sim/ludwig/capacitor/local/WALL-DIRICHLET/H-microgel-ruido"
L = 32.0
KT, ETA, RHO, AH, AL, K, R0 = 1e-4, 0.01, 1.0, 0.3078881242253297, 0.4743549210637387, 1e-3, 4.0
ZLO, ZHI, RK = 1.5, 33.5, 3.5
DRMAX = 0.8


def read_cds(path):
    lines = open(path).read().split("\n")
    n = int(lines[0].split()[0])
    pos = np.empty((n, 3))
    idx = np.empty(n, int)
    bonds = []
    for c in range(n):
        b = 1 + 68 * c
        idx[c] = int(lines[b].split()[0])
        pos[c] = [float(x) for x in lines[b + 34].split()[:3]]
        for j in range(9, 15):
            p = int(lines[b + j].split()[0])
            if p:
                bonds.append((idx[c], p))
    order = np.argsort(idx)
    return idx[order], pos[order], bonds


def analyse(d):
    files = sorted(glob.glob(f"{d}/config.cds0*.001-001"))
    if not files:
        print(f"== {d}: sin config.cds"); return
    steps = np.array([int(os.path.basename(f)[10:18]) for f in files])
    idx, p0, bonds = read_cds(files[0])
    bset = sorted({tuple(sorted(b)) for b in bonds})
    pos = np.array([read_cds(f)[1] for f in files])           # [t, n, 3]
    # 1. desplazamiento por paso (solo pares de volcados consecutivos)
    cons = np.diff(steps) == 1
    dr = np.diff(pos, axis=0)[cons]
    dr[..., :2] -= L * np.round(dr[..., :2] / L)
    mag = np.linalg.norm(dr, axis=-1)
    Dbare = KT / (6 * math.pi * ETA) * (1 / AH - 1 / AL)
    Dst = KT / (6 * math.pi * ETA * AH)
    t_max = np.unravel_index(np.argmax(mag), mag.shape)
    # 2. energía de enlaces
    i_of = {v: k for k, v in enumerate(idx)}
    U = []
    for P in pos:
        u = 0.0
        for a, b in bset:
            v = P[i_of[a]] - P[i_of[b]]
            v[:2] -= L * np.round(v[:2] / L)
            u += 0.5 * K * (np.linalg.norm(v) - R0) ** 2
        U.append(u / KT)
    U = np.array(U)
    half = len(U) // 2
    # 4. paredes
    dlo, dhi = pos[..., 2] - ZLO, ZHI - pos[..., 2]
    dmin = np.minimum(dlo, dhi)
    cm = pos.mean(axis=1)
    rg = np.sqrt(((pos - cm[:, None, :]) ** 2).sum(-1).mean(-1))
    log = open(f"{d}/logs/output.txt", errors="replace").read()
    fin = "terminó normalmente" if "finished normally" in log else "NO terminó (ver logs/output.txt)"
    print(f"\n== {os.path.basename(d)}: pasos {steps[0]}..{steps[-1]} ({len(files)} volcados), {fin}")
    print(f" 1. desplazamiento por paso: rms por componente {np.sqrt((dr**2).mean()):.4f} "
          f"(x {np.sqrt((dr[...,0]**2).mean()):.4f}, y {np.sqrt((dr[...,1]**2).mean()):.4f}, z {np.sqrt((dr[...,2]**2).mean()):.4f}); "
          f"ref. sqrt(2 D_desnudo) = {math.sqrt(2*Dbare):.4f}, sqrt(2 D_Stokes) = {math.sqrt(2*Dst):.4f}")
    print(f"    máximo |dr| = {mag.max():.4f} = {mag.max()/DRMAX*100:.1f} % de drmax (paso {steps[1:][cons][t_max[0]]}, "
          f"monómero {idx[t_max[1]]}); percentil 99.99 = {np.percentile(mag, 99.99):.4f}")
    print(f" 2. enlaces ({len(bset)}): <U>/kT = {U[half:].mean():.2f} ± {U[half:].std()/math.sqrt(max(1,half/50)):.2f} "
          f"(segunda mitad; error con bloques de ~50 pasos, heurístico) vs 35/2 = 17.5; primera mitad {U[:half].mean():.2f}")
    print(f" 4. paredes: d mínimo {dmin.min():.3f} (h = d - 3.5 = {dmin.min()-RK:.3f}); "
          f"fracción de muestras con h < 0: {(dmin < RK).mean():.2e}")
    print(f"    Rg: {rg[0]:.3f} -> {rg[-1]:.3f}; z_cm: {cm[0,2]:.2f} -> {cm[-1,2]:.2f}")
    # 3. fluido
    for vf in sorted(glob.glob(f"{d}/vel-0*.001-001")):
        u = np.loadtxt(vf)
        if u.ndim != 2 or u.shape[1] < 3:
            continue
        u = u[:, -3:].reshape(32, 32, 34, 3)[:, :, 1:33, :]       # x lento, z rápido; fluido z = 2..33
        print(f" 3. fluido {os.path.basename(vf)[:13]}: <u_a^2>/(kT/rho) = "
              + ", ".join(f"{(u[..., a]**2).mean()/(KT/RHO):.3f}" for a in range(3)))


def main():
    dirs = sys.argv[1:] or sorted(glob.glob(f"{H}/*/"))
    for d in dirs:
        analyse(d.rstrip("/"))


if __name__ == "__main__":
    main()
