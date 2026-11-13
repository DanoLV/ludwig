#!/usr/bin/env python3
"""
Estado de equilibrio de una corrida de microgel a partir de su serie de Rg.

Ajusta la pendiente dRg/dpaso de los ultimos volcados a una exponencial
  s(t) = s0 exp(-(t - t0)/tau)
y reporta:
  - Rg actual y pendiente (por 1000 pasos, promediada en 10000 pasos)
  - tau y el Rg asintotico proyectado  Rg_inf = Rg + s * tau
  - si se cumple el criterio de equilibrio y, si no, en que paso se cumpliria

Criterio (por defecto): pendiente < 2e-4 por 1000 pasos, lo que con tau ~1.3e5
deja menos de ~0.4% de hinchamiento por delante.

Uso:
  mgel_equilibrio.py DIRECTORIO [--umbral 2e-4] [--ventana 50000] [-o archivo]
"""
import argparse
import math
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mgel_serie import measure, read_config   # noqa: E402


def serie_rg(d, cada):
    files = []
    for f in os.listdir(d):
        m = re.fullmatch(r"config\.cds(\d{8})\.001-001", f)
        if m and int(m.group(1)) % cada == 0:
            files.append((int(m.group(1)), os.path.join(d, f)))
    files.sort()
    return [(s, measure(read_config(p))["Rg"]) for s, p in files]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dir")
    ap.add_argument("--umbral", type=float, default=2e-4,
                    help="pendiente dRg/1000 pasos que se considera equilibrio")
    ap.add_argument("--ventana", type=int, default=50000,
                    help="pasos finales usados en el ajuste")
    ap.add_argument("--cada", type=int, default=5000,
                    help="usar solo volcados multiplos de este paso")
    ap.add_argument("-o", "--out", default=None)
    a = ap.parse_args()

    rg = serie_rg(a.dir, a.cada)
    lines = []
    if len(rg) < 4:
        lines.append("todavia no hay suficientes volcados")
    else:
        # pendientes centradas cada 10000 pasos (diferencias entre volcados
        # separados por 10000) dentro de la ventana final
        t_end = rg[-1][0]
        d = dict(rg)
        pts = []
        for s, r in rg:
            if s - 10000 in d and s >= t_end - a.ventana:
                pts.append((s - 5000, 1000.0 * (r - d[s - 10000]) / 10000.0))
        lines.append("paso actual         : %d" % t_end)
        lines.append("Rg actual           : %.4f" % rg[-1][1])
        if len(pts) >= 3 and all(p[1] > 0 for p in pts):
            n = len(pts)
            xs = [p[0] for p in pts]
            ys = [math.log(p[1]) for p in pts]
            mx, my = sum(xs) / n, sum(ys) / n
            b = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / \
                sum((x - mx) ** 2 for x in xs)
            tau = -1.0 / b if b < 0 else float("inf")
            s_now = math.exp(my + b * (t_end - mx))      # pendiente ajustada hoy
            rg_inf = rg[-1][1] + s_now * tau / 1000.0
            lines.append("pendiente dRg/1000  : %.5f (ajuste), %.5f (ultimos 10000)"
                         % (s_now, pts[-1][1]))
            lines.append("tau                 : %.0f pasos" % tau)
            lines.append("Rg asintotico       : %.4f (falta %.2f%%)"
                         % (rg_inf, 100 * (rg_inf - rg[-1][1]) / rg_inf))
            if s_now < a.umbral:
                lines.append("EQUILIBRIO          : SI (pendiente < %.1e)" % a.umbral)
            else:
                t_eq = t_end + tau * math.log(s_now / a.umbral)
                lines.append("EQUILIBRIO          : no; criterio (%.1e) estimado en el paso %.0f"
                             % (a.umbral, t_eq))
        else:
            lines.append("pendientes no monotonas o negativas; revisar la serie a mano")
            for p in pts[-5:]:
                lines.append("  t=%d  dRg/1000=%.5f" % p)

    # Aviso: el binario de esta corrida es anterior al arreglo de la fuerza
    # duplicada en celdas de borde (subgrid.c, 2026-09-26). Con celdas de 12 en
    # una caja de 48, una particula con alguna coordenada < 12.5 o > 36.5
    # recibiria 2x (o mas) la fuerza electrostatica.
    try:
        fs = sorted(f for f in os.listdir(a.dir) if re.fullmatch(r"config\.cds\d{8}\.001-001", f))
        cfg = read_config(os.path.join(a.dir, fs[-1]))
        P = [v[0] for v in cfg.values()]
        lo = min(min(p) for p in P); hi = max(max(p) for p in P)
        nb = sum(1 for p in P if any(x < 12.5 or x > 36.5 for x in p))
        lines.append("extension (x,y,z)   : %.2f .. %.2f   (celdas de borde: <12.5 o >36.5)" % (lo, hi))
        if nb:
            lines.append("AVISO               : %d particulas en celdas de borde -> fuerza duplicada con este binario" % nb)
    except Exception as exc:
        lines.append("extension           : no se pudo calcular (%s)" % exc)

    txt = "\n".join(lines) + "\n"
    if a.out:
        with open(a.out, "w") as fh:
            fh.write(txt)
    print(txt, end="")
    return 0


if __name__ == "__main__":
    sys.exit(main())
