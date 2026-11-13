#!/usr/bin/env python3
"""Test de linealidad (Hann-8, L = 32), q = 0.25, 1, 2 en la misma configuracion.
La fuerza medida es F = q E, y el campo de imagen es proporcional a q, asi que F/q^2 y Q_inf/q son
independientes de q si Debye-Huckel linealizado es exacto; la teoria lineal da un unico valor de F/q^2."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import numpy as np
import grafico_C_electrolito as G
kmodes = G.setup(8)
print("LD    d     q     F/q^2         Qinf/q     |  teoria F/q^2  teoria/medido")
for ld, ds in ((4, (3.5, 6.0, 10.0)), (8, (6.0,))):
    for d in ds:
        z = d + 1.5
        F, Ql, Qh = G.theory(d, 8, 1.0 / ld, kmodes)
        for q, bat in ((0.25, "C-lineal-hann8-LD%d-q0.25"), (1.0, "C-electrolito-hann8-LD%d"), (2.0, "C-lineal-hann8-LD%d-q2")):
            p = os.path.join(ROOT, bat % ld, "z_%05.2f" % z)
            if not os.path.exists(os.path.join(p, "proceced_data", "particle_force.csv")) or "finished normally" not in open(os.path.join(p, "run.log"), errors="ignore").read():
                print("%d  %4.1f  %4.2f  (sin terminar)" % (ld, d, q)); continue
            e, ql, qh, step = G.medido(p)
            print("%d  %4.1f  %4.2f  %+.6e  %+.5f  |  %+.6e  %.5f" % (ld, d, q, e / q**2, ql / q, F, F / (e / q**2)))
