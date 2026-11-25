#!/usr/bin/env python3
"""
TEST-NOFLUX-ADV (2026-11-16): conservación de los iones con advección centred + LB térmico + paredes.

Para cada corrida <ver>/<caso> y cada volcado qsi (cada 100 pasos):
  N_a(t) = sum de n_a sobre los nodos FLUIDOS (con pared: z = 2..L-1; sin pared: todos), a = +, -
  dN_a(t)/N_a(0) = cambio relativo respecto del paso 1 (la masa de iones de NP debería conservarse a redondeo)
También: n en los nodos de pared (deben quedar en el valor inicial) y RMS de u_z en la primera capa de fluido
(velocidad normal a la pared, la que según la lectura del código produce la pérdida).
Al final, si existen antes/ y despues/ del caso sinpared, compara qsi byte a byte (control de regresión).
Orden de los datos: x lento, z rápido (como los demás scripts del proyecto).
Solo lee archivos. Uso: python3 capacitor/scripts/TEST-NOFLUX-ADV/analizar.py
"""
import filecmp
import glob
import os
import numpy as np

BASE = "/home/danolv/Sim/ludwig/capacitor/local/WALL-DIRICHLET/TEST-NOFLUX-ADV"


def size(d):
    for l in open(os.path.join(d, "input")):
        if l.startswith("size"):
            return [int(x) for x in l.split()[1].split("_")]


for ver in ("antes", "despues"):
    for caso in ("pared", "sinpared"):
        d = os.path.join(BASE, ver, caso)
        if not glob.glob(os.path.join(d, "qsi-0*")):
            continue
        L = size(d)
        wall = caso == "pared"
        zf = slice(1, L[2] - 1) if wall else slice(0, L[2])
        print("\n== %s/%s  (size %s)" % (ver, caso, L))
        print("  step | N_+ fluid        dN_+/N_+(1)  | N_- fluid        dN_-/N_-(1)  |"
              " n_wall (+,-)          | RMS u_z capa 1")
        N1 = None
        for q in sorted(glob.glob(os.path.join(d, "qsi-0*"))):
            st = int(os.path.basename(q)[4:13])
            n = np.loadtxt(q).reshape(L[0], L[1], L[2], 2)
            N = n[:, :, zf].sum(axis=(0, 1, 2))
            N1 = N if N1 is None else N1
            rel = (N - N1) / N1
            nw = (n[:, :, 0].mean(axis=(0, 1)) if wall else [np.nan, np.nan])
            v = os.path.join(d, "vel-%09d.001-001" % st)
            uz = np.nan
            if os.path.exists(v):
                u = np.loadtxt(v).reshape(L[0], L[1], L[2], 3)
                uz = np.sqrt(np.mean(u[:, :, 1 if wall else 0, 2] ** 2))
            print("  %4d | %.10e  %+.3e   | %.10e  %+.3e   | %.6e %.6e | %.3e" % (
                st, N[0], rel[0], N[1], rel[1], nw[0], nw[1], uz))

a, b = os.path.join(BASE, "antes", "sinpared"), os.path.join(BASE, "despues", "sinpared")
qa = sorted(glob.glob(os.path.join(a, "qsi-0*")))
if qa and os.path.isdir(b):
    same = [filecmp.cmp(q, os.path.join(b, os.path.basename(q)), shallow=False)
            for q in qa if os.path.exists(os.path.join(b, os.path.basename(q)))]
    print("\nsinpared antes vs despues: %d/%d archivos qsi identicos byte a byte" % (sum(same), len(same)))
