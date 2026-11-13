#!/usr/bin/env python3
"""
Guarda en <proyecto>/saved_data/ lo necesario para rehacer cada corrida (2026-11-13,
pedido del usuario), sin guardar un input ni un archivo de coloides por corrida:

  saved_data/<carpeta-serie>/input.base          input real (con comentarios) de la corrida
                                                 que menos difiere del resto de la carpeta
  saved_data/<carpeta-serie>/coloides_<k>.cds    archivo(s) de coloides base de la carpeta
  saved_data/<carpeta-serie>/deltas.json         por corrida: claves del input que cambian /
                                                 se agregan / se quitan respecto de input.base,
                                                 base de coloides que usa y líneas que cambian

Coloides: dentro de una carpeta-serie el contenido más frecuente es la primera base; otro
contenido con el mismo número de líneas se guarda como cambios de línea si difiere en
<= MAX_LINEAS líneas (posición, carga, ...); si difiere en más, pasa a ser otra base.
(MAX_LINEAS es una elección práctica, no un criterio físico.)

Reconstrucción: util/rehacer_corrida.py. Verificación: al final este script reconstruye
TODAS las corridas y compara (coloides byte a byte; input clave a clave).
No modifica local/ ni runs/; reescribe saved_data/<carpeta-serie>/ de las carpetas que procesa.

Uso:  python3 util/guardar_bases.py [capacitor] [microgel]
"""
import collections
import hashlib
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import registro_simulaciones as R  # noqa: E402
import rehacer_corrida as H        # noqa: E402

MAX_LINEAS = 40
COLL = R.COLLOID_FILE


def key_lines(text):
    """{clave: valor} con la misma regla que Ludwig/registro (la última aparición gana)."""
    d = {}
    for line in text.splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        s = s.split("#", 1)[0].strip()
        if not s:
            continue
        tok = s.split(None, 1)
        d[tok[0]] = tok[1].strip() if len(tok) > 1 else ""
    return d


def group_dir(project, group):
    return os.path.join(R.ROOT, project, "saved_data", group.replace(" + ", "__"))


def process(project):
    runs = R.find_runs(project)
    groups = collections.defaultdict(list)
    for rel, loc in runs.items():
        d = loc["local"] or loc["runs"]
        groups[R.SERIES[project](rel)[1]].append((rel, d))
    tot_bytes = 0
    report = []
    for group, lst in sorted(groups.items()):
        out = group_dir(project, group)
        os.makedirs(out, exist_ok=True)
        for f in os.listdir(out):                      # rehacer la carpeta desde cero
            if f == "input.base" or f.startswith("coloides_") or f == "deltas.json":
                os.remove(os.path.join(out, f))
        texts, params, colls, skipped = {}, {}, {}, []
        for rel, d in lst:
            raw = open(os.path.join(d, "input"), "rb").read()
            if b"\0" in raw[:4096]:
                skipped.append(rel)
                continue
            texts[rel] = raw.decode(errors="replace")
            params[rel] = key_lines(texts[rel])
            cf = os.path.join(d, COLL)
            if os.path.isfile(cf):
                colls[rel] = open(cf, "rb").read()
        if not texts:
            continue
        # --- input base: la corrida con menos diferencias respecto del resto
        def ndiff(a, b):
            ka, kb = params[a], params[b]
            return sum(1 for k in set(ka) | set(kb) if ka.get(k) != kb.get(k))
        rels = sorted(texts)
        if len(rels) > 60:                     # grupos grandes: comparar contra una muestra
            sample = rels[:: max(1, len(rels) // 60)]
        else:
            sample = rels
        medoid = min(rels, key=lambda a: sum(ndiff(a, b) for b in sample))
        base_text = texts[medoid]
        open(os.path.join(out, "input.base"), "w").write(base_text)
        base_p = params[medoid]
        # --- coloides: bases y deltas
        freq = collections.Counter(colls.values())
        bases = []                                  # [(nombre, bytes, lineas)]
        assign = {}
        for content, _ in freq.most_common():
            lines = content.decode("latin-1").splitlines(keepends=True)
            best = None
            for name, bb, bl in bases:
                if len(bl) == len(lines):
                    nd = sum(1 for x, y in zip(bl, lines) if x != y)
                    if nd <= MAX_LINEAS and (best is None or nd < best[1]):
                        best = (name, nd, bl)
            if best is None:
                name = f"coloides_{len(bases) + 1}.cds"
                bases.append((name, content, lines))
                open(os.path.join(out, name), "wb").write(content)
                assign[content] = (name, {})
            else:
                name, _, bl = best
                assign[content] = (name, {str(i + 1): y for i, (x, y) in enumerate(zip(bl, lines)) if x != y})
        deltas = {"_info": {"input_base_de": medoid, "max_lineas_delta_coloides": MAX_LINEAS,
                            "sin_reconstruccion (input ilegible)": skipped}}
        for rel in rels:
            p = params[rel]
            dd = {"input_set": {k: v for k, v in p.items() if base_p.get(k) != v},
                  "input_remove": sorted(k for k in base_p if k not in p)}
            if rel in colls:
                name, ch = assign[colls[rel]]
                dd["coloides_base"] = name
                dd["coloides_lineas"] = ch
            deltas[rel] = dd
        json.dump(deltas, open(os.path.join(out, "deltas.json"), "w"), indent=1, ensure_ascii=False)
        size = sum(os.path.getsize(os.path.join(out, f)) for f in os.listdir(out))
        tot_bytes += size
        report.append((group, len(rels), len(bases), size))
    # --- verificación: reconstruir todas y comparar
    bad = []
    for rel, loc in runs.items():
        d = loc["local"] or loc["runs"]
        if b"\0" in open(os.path.join(d, "input"), "rb").read(4096):
            continue
        inp, col = H.rehacer(project, rel)
        if key_lines(inp) != key_lines(open(os.path.join(d, "input"), errors="replace").read()):
            bad.append((rel, "input"))
        cf = os.path.join(d, COLL)
        if os.path.isfile(cf):
            if col is None or hashlib.md5(col).hexdigest() != hashlib.md5(open(cf, "rb").read()).hexdigest():
                bad.append((rel, "coloides"))
    print(f"{project}: {len(report)} carpetas-serie, {sum(r[1] for r in report)} corridas, "
          f"{tot_bytes / 1e6:.2f} MB en saved_data")
    for g, n, nb, sz in sorted(report, key=lambda r: -r[3])[:8]:
        print(f"   {g[:60]:60s} {n:4d} corridas  {nb:3d} base(s) de coloides  {sz / 1e3:8.1f} kB")
    print(f"   verificación: {len(bad)} corridas NO reconstruidas exactamente" + (f": {bad[:10]}" if bad else ""))


if __name__ == "__main__":
    for p in (sys.argv[1:] or ["capacitor", "microgel"]):
        process(p)
