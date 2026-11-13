#!/usr/bin/env python3
"""
Escribe un config.cds a partir de la estructura base equilibrada, aplicando un
perfil de cargas. La geometria (posiciones, enlaces) nunca se toca: todas las
corridas comparten exactamente la misma estructura, de modo que cualquier
diferencia que se mida viene del perfil de carga o de la sal, no del punto de
partida.

Indices dentro del bloque de cada coloide en el formato config.cds:
  0  = indice (id) del coloide      34 = posicion
  2  = numero de enlaces            35 = velocidad lineal
  9+ = ids de los coloides ligados  36 = velocidad angular
                                    45 = q0
                                    46 = q1

Uso:
  set_charges.py --mode all        -o config.cds.init.001-001
  set_charges.py --mode surface    -o ...
  set_charges.py --mode random --frac 0.4 --seed 7  -o ...
  set_charges.py --mode shell --rmin 5.5            -o ...
  set_charges.py --mode none       -o ...
"""

import argparse
import math
import os
import random
import sys

BASE_DEFAULT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "..", "MICROGEL-nuevo", "base_equilibrada.cds")
SURFACE_REF = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "..", "MICROGEL-nuevo", "superficie",
                           "config.cds.init.001-001")

I_ID, I_NBOND, I_BOND0, I_POS, I_Q0, I_Q1 = 0, 2, 9, 34, 45, 46


def read_config(path):
    with open(path) as fh:
        lines = fh.readlines()
    n = int(lines[0])
    rest = lines[1:]
    b = len(rest) // n
    blocks = {}
    for i in range(n):
        blk = rest[i * b:(i + 1) * b]
        blocks[int(blk[I_ID])] = list(blk)
    return lines[0], blocks, b


def positions(blocks):
    return {i: [float(x) for x in blk[I_POS].split()]
            for i, blk in blocks.items()}


def surface_ids():
    """Los monomeros que el generador marco como superficie, tomados del
    config que produjo en modo surface. Se usa esa lista y no un criterio
    geometrico nuevo para que la definicion de superficie sea la misma que
    en las corridas anteriores."""
    _, blocks, _ = read_config(SURFACE_REF)
    return {i for i, blk in blocks.items() if float(blk[I_Q0]) != 0.0}


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--mode", required=True,
                   choices=("all", "surface", "random", "shell", "none", "outer"))
    p.add_argument("--q", type=float, default=1.0, help="carga por monomero cargado")
    p.add_argument("--frac", type=float, default=0.4,
                   help="fraccion cargada en modo random (solo interior)")
    p.add_argument("--seed", type=int, default=12345)
    p.add_argument("--rmin", type=float, default=5.5,
                   help="radio interno en modo shell: se carga todo r >= rmin")
    p.add_argument("--nq", type=int, default=65,
                   help="modo outer: cuantos monomeros externos se cargan. "
                        "Fija la carga total (y por lo tanto la concentracion "
                        "de contraiones) sin depender de cuantos monomeros "
                        "tenga la estructura.")
    p.add_argument("--base", default=BASE_DEFAULT)
    p.add_argument("-o", "--out", required=True)
    a = p.parse_args()

    hdr, blocks, _ = read_config(a.base)
    ids = sorted(blocks)
    pos = positions(blocks)

    if a.mode == "none":
        chosen = set()
    elif a.mode == "all":
        chosen = set(ids)
    elif a.mode == "surface":
        chosen = surface_ids()
    elif a.mode == "random":
        interior = sorted(set(ids) - surface_ids())
        rng = random.Random(a.seed)
        chosen = set(i for i in interior if rng.random() < a.frac)
    elif a.mode == "shell":
        n = len(ids)
        c = [sum(pos[i][k] for i in ids) / n for k in range(3)]
        chosen = {i for i in ids if math.dist(pos[i], c) >= a.rmin}
    elif a.mode == "outer":
        # Los nq monomeros mas alejados del centro de masa. El modo "surface"
        # depende de una lista que solo existe para el microgel de cadena; este
        # es geometrico, sirve para cualquier estructura, y al fijar el NUMERO
        # de cargas deja la carga total (y lambda_D de los contraiones) igual
        # entre corridas con estructuras distintas.
        n = len(ids)
        c = [sum(pos[i][k] for i in ids) / n for k in range(3)]
        chosen = set(sorted(ids, key=lambda i: -math.dist(pos[i], c))[:a.nq])

    for i in ids:
        q = a.q if i in chosen else 0.0
        blocks[i][I_Q0] = "   %.15e\n" % q
        blocks[i][I_Q1] = "   %.15e\n" % 0.0

    out = [hdr]
    for i in ids:
        out.extend(blocks[i])
    with open(a.out, "w") as fh:
        fh.writelines(out)

    n = len(ids)
    c = [sum(pos[i][k] for i in ids) / n for k in range(3)]
    rr = [math.dist(pos[i], c) for i in ids]
    rg = math.sqrt(sum(x * x for x in rr) / n)
    rc = [math.dist(pos[i], c) for i in sorted(chosen)]
    print("%s  modo=%s  cargados=%d/%d  q=%g  Q_total=%+g  Rg=%.4f%s"
          % (a.out, a.mode, len(chosen), n, a.q, len(chosen) * a.q, rg,
             ("  <r>_cargados=%.2f" % (sum(rc) / len(rc))) if rc else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
