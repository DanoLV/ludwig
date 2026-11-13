#!/usr/bin/env python3
"""Aditividad: F_i(N cuerpos) vs suma de F_i medidas en los pares aislados.
Fuerzas del paso 1 (nube ionica aun uniforme => electrostatica pura)."""
import os, itertools, numpy as np
D=os.path.dirname(os.path.abspath(__file__))

def forces(tag, step=1):
    """dict {id_particula: F} leyendo particle_force.csv del run 'p<ids>'"""
    ids=[int(c) for c in tag[1:]]
    out={}
    for l in open(f'{D}/{tag}/proceced_data/particle_force.csv'):
        if not l[0].isdigit(): continue
        p=l.strip().split(';')
        if int(p[0])!=step: continue
        k=int(p[1])                       # indice en el archivo = orden en el config
        out[ids[k-1]]=np.array([float(x.replace(',','.')) for x in p[3:6]])
    return out

print(f"{'config':>8} {'part':>5}  {'|F| N cuerpos':>14} {'|F| suma pares':>15} {'|dif|':>11} {'rel':>10}")
for n in (2,3,4):
    for ids in itertools.combinations([1,2,3,4], n):
        tag='p'+''.join(map(str,ids))
        if not os.path.exists(f'{D}/{tag}/proceced_data/particle_force.csv'): continue
        F=forces(tag)
        for i in ids:
            s=np.zeros(3)
            for j in ids:
                if j==i: continue
                pt='p'+''.join(map(str,sorted((i,j))))
                s+=forces(pt)[i]
            d=np.linalg.norm(F[i]-s); rel=d/np.linalg.norm(F[i])
            print(f"{tag:>8} {i:>5}  {np.linalg.norm(F[i]):14.6e} {np.linalg.norm(s):15.6e} {d:11.3e} {rel:10.2e}")
