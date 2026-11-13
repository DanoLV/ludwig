#!/usr/bin/env python3
"""
Sigue la evolucion estructural de un microgel a lo largo de una corrida.

Lee los config.cds volcados por Ludwig y, para cada paso, mide:
  Rg     radio de giro (RMS de la distancia al centro de masa)
  diam   2 * distancia maxima al centro de masa
  bond   longitud media de enlace (y su minimo/maximo)
  hmin   separacion minima entre monomeros no ligados
  vmax   modulo de velocidad maxima
  z      enlaces por monomero (constante: la topologia no cambia)

El equilibrio se juzga por dRg/dpaso: mientras crezca de forma sostenida el
microgel se sigue hinchando. Con la estructura de cadena ramificada (z=1.92)
esa pendiente nunca satura porque no hay elasticidad que la frene; con una
red entrecruzada debe caer a cero.

Uso:
  mgel_serie.py DIRECTORIO [-o serie.json] [--every N]
"""

import argparse
import itertools
import json
import math
import os
import re
import sys

I_ID, I_NBOND, I_BOND0, I_POS, I_VEL = 0, 2, 9, 34, 35


def read_config(path):
    """Devuelve {id: (pos, vel, [ids ligados])}."""
    with open(path) as fh:
        lines = fh.readlines()
    n = int(lines[0])
    rest = lines[1:]
    b = len(rest) // n
    out = {}
    for i in range(n):
        blk = rest[i * b:(i + 1) * b]
        nb = int(blk[I_NBOND])
        out[int(blk[I_ID])] = (
            [float(x) for x in blk[I_POS].split()],
            [float(x) for x in blk[I_VEL].split()],
            [int(blk[I_BOND0 + k]) for k in range(nb)],
        )
    return out


def measure(cfg):
    ids = sorted(cfg)
    n = len(ids)
    pos = {i: cfg[i][0] for i in ids}
    c = [sum(pos[i][k] for i in ids) / n for k in range(3)]
    rr = [math.dist(pos[i], c) for i in ids]

    bonds = set()
    for i in ids:
        for j in cfg[i][2]:
            bonds.add((min(i, j), max(i, j)))
    bl = [math.dist(pos[i], pos[j]) for i, j in bonds]

    # separacion minima entre pares NO ligados: es la que vigila el potencial
    # repulsivo, y la que avisa si la red se esta colapsando sobre si misma
    hmin = math.inf
    for i, j in itertools.combinations(ids, 2):
        if (i, j) in bonds:
            continue
        d = math.dist(pos[i], pos[j])
        if d < hmin:
            hmin = d

    vmax = max(math.dist(cfg[i][1], [0, 0, 0]) for i in ids)

    return {
        "n": n,
        "Rg": math.sqrt(sum(x * x for x in rr) / n),
        "diam": 2.0 * max(rr),
        "bond": sum(bl) / len(bl) if bl else 0.0,
        "bond_min": min(bl) if bl else 0.0,
        "bond_max": max(bl) if bl else 0.0,
        "z": 2.0 * len(bonds) / n,
        "hmin": hmin,
        "vmax": vmax,
    }


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("dirs", nargs="+", help="directorio(s) de corrida")
    p.add_argument("-o", "--out", default="serie.json",
                   help="nombre del json dentro de cada directorio")
    p.add_argument("--every", type=int, default=1,
                   help="procesar 1 de cada N volcados")
    a = p.parse_args()

    for d in a.dirs:
        files = []
        for f in os.listdir(d):
            m = re.fullmatch(r"config\.cds(\d{8})\.001-001", f)
            if m:
                files.append((int(m.group(1)), os.path.join(d, f)))
        files.sort()
        files = files[::a.every]
        if not files:
            print("%s: sin volcados" % d, file=sys.stderr)
            continue

        serie = {}
        for step, path in files:
            serie[str(step)] = measure(read_config(path))

        with open(os.path.join(d, a.out), "w") as fh:
            json.dump(serie, fh)

        steps = [s for s, _ in files]
        print("\n%s  (%d volcados, %d monomeros, z=%.3f)"
              % (d, len(steps), serie[str(steps[0])]["n"], serie[str(steps[0])]["z"]))
        print("%9s %8s %8s %8s %8s %8s %11s"
              % ("paso", "Rg", "diam", "<enlace>", "hmin", "vmax", "dRg/1000"))
        prev = None
        for s in steps:
            v = serie[str(s)]
            slope = ""
            if prev is not None and s != prev[0]:
                slope = "%+11.5f" % (1000.0 * (v["Rg"] - prev[1]) / (s - prev[0]))
            print("%9d %8.4f %8.3f %8.4f %8.4f %8.2e %11s"
                  % (s, v["Rg"], v["diam"], v["bond"], v["hmin"], v["vmax"], slope))
            prev = (s, v["Rg"])

    return 0


if __name__ == "__main__":
    sys.exit(main())
