#!/usr/bin/env python3
"""Exporta una configuracion de coloides de Ludwig (config.cds*.001-001) a ParaView.
Uso: python3 exportar_csv.py [config] [prefijo_salida] [L]
Genera <prefijo>_puntos.csv (id,x,y,z,nbonds,isfixed), <prefijo>_enlaces.csv
(id1,id2,longitud con imagen minima,periodico) y <prefijo>.vtk (puntos + enlaces
no periodicos como lineas; los que cruzan la caja se omiten para no dibujar
lineas de lado a lado)."""
import sys
import numpy as np

cfg = sys.argv[1] if len(sys.argv) > 1 else "config.cds00003000.001-001"
out = sys.argv[2] if len(sys.argv) > 2 else "solido_final"
L = float(sys.argv[3]) if len(sys.argv) > 3 else 48.0

tok = open(cfg).read().split()
n = int(tok[0]); a = tok[1:]
assert len(a) == 80 * n, "formato inesperado (se esperan 32 enteros + 48 reales por coloide)"
I = np.array([[int(x) for x in a[80*i:80*i+32]] for i in range(n)])
D = np.array([[float(x) for x in a[80*i+32:80*i+80]] for i in range(n)])
ids, nb, bonds, fixed = I[:, 0], I[:, 2], I[:, 9:15], I[:, 4]
pos = D[:, 2:5]
order = np.argsort(ids); ids, nb, bonds, fixed, pos = ids[order], nb[order], bonds[order], fixed[order], pos[order]
row = {int(i): m for m, i in enumerate(ids)}

with open(out + "_puntos.csv", "w") as f:
    f.write("id,x,y,z,nbonds,isfixed\n")
    for m in range(n):
        f.write(f"{ids[m]},{pos[m,0]:.8f},{pos[m,1]:.8f},{pos[m,2]:.8f},{nb[m]},{fixed[m]}\n")

pairs = sorted({(min(m, row[int(bonds[m, k])]), max(m, row[int(bonds[m, k])]))
                for m in range(n) for k in range(nb[m])})
lines = []
with open(out + "_enlaces.csv", "w") as f:
    f.write("id1,id2,longitud,periodico\n")
    for p, q in pairs:
        d = pos[q] - pos[p]; dm = d - L * np.round(d / L)
        per = int(np.abs(d - dm).max() > 1e-9)
        f.write(f"{ids[p]},{ids[q]},{np.linalg.norm(dm):.8f},{per}\n")
        if not per: lines.append((p, q))

with open(out + ".vtk", "w") as f:
    f.write("# vtk DataFile Version 3.0\nsolido poroso\nASCII\nDATASET POLYDATA\n")
    f.write(f"POINTS {n} double\n")
    for m in range(n): f.write(f"{pos[m,0]:.8f} {pos[m,1]:.8f} {pos[m,2]:.8f}\n")
    f.write(f"VERTICES {n} {2*n}\n")
    for m in range(n): f.write(f"1 {m}\n")
    f.write(f"LINES {len(lines)} {3*len(lines)}\n")
    for p, q in lines: f.write(f"2 {p} {q}\n")
    f.write(f"POINT_DATA {n}\nSCALARS nbonds int 1\nLOOKUP_TABLE default\n")
    for m in range(n): f.write(f"{nb[m]}\n")
print(f"{n} puntos, {len(pairs)} enlaces ({len(pairs)-len(lines)} periodicos, omitidos en el .vtk)")
