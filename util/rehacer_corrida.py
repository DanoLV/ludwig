#!/usr/bin/env python3
"""
Reconstruye el input y el config.cds.init.001-001 de una corrida a partir de
<proyecto>/saved_data/<carpeta-serie>/ (input.base, coloides_*.cds, deltas.json),
generados por util/guardar_bases.py (2026-11-13).

  python3 util/rehacer_corrida.py <proyecto> <ruta de la corrida> [carpeta destino]

<ruta de la corrida> es la columna "Ruta" del Registro_simulaciones.xlsx
(p. ej. SLIDES-LD8-hann08-corr_yes/pos_12.00_24.00_24.00). Sin carpeta destino solo
muestra los cambios respecto de la base. Con destino escribe ahí "input" y
"config.cds.init.001-001" (no sobrescribe archivos existentes).

El input reconstruido tiene las mismas claves y valores que el original; los
comentarios son los de input.base (las líneas cambiadas llevan "# [esta corrida]").
El archivo de coloides se reconstruye byte a byte.
Archivos auxiliares que el input pueda nombrar (.petscrc, tablas, anillo.cds, ...) no
se guardan aquí.
"""
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _group_dir(project, rel):
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import registro_simulaciones as R
    group = R.SERIES[project](rel)[1]
    return os.path.join(ROOT, project, "saved_data", group.replace(" + ", "__"))


def rehacer(project, rel):
    """Devuelve (texto del input, bytes del archivo de coloides o None)."""
    gd = _group_dir(project, rel)
    deltas = json.load(open(os.path.join(gd, "deltas.json")))
    if rel not in deltas:
        raise KeyError(f"{rel} no está en {gd}/deltas.json")
    dd = deltas[rel]
    lines = open(os.path.join(gd, "input.base")).read().splitlines()
    sets, remove = dict(dd["input_set"]), set(dd["input_remove"])
    out, done = [], set()
    for line in lines:
        s = line.strip()
        body = s.split("#", 1)[0].strip()
        key = body.split(None, 1)[0] if body and not s.startswith("#") else None
        if key in remove:
            out.append("# [quitada en esta corrida] " + line)
        elif key in sets:
            out.append(f"{key} {sets[key]}".rstrip() + "   # [esta corrida]")
            done.add(key)
        else:
            out.append(line)
    extra = [k for k in sets if k not in done]
    if extra:
        out.append("")
        out.append("# --- claves propias de esta corrida (no están en input.base) ---")
        out += [f"{k} {sets[k]}".rstrip() for k in extra]
    text = "\n".join(out) + "\n"
    col = None
    if "coloides_base" in dd:
        blines = open(os.path.join(gd, dd["coloides_base"]), "rb").read().decode("latin-1").splitlines(keepends=True)
        for i, v in dd["coloides_lineas"].items():
            blines[int(i) - 1] = v
        col = "".join(blines).encode("latin-1")
    return text, col


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)
    project, rel = sys.argv[1], sys.argv[2]
    text, col = rehacer(project, rel)
    if len(sys.argv) < 4:
        dd = json.load(open(os.path.join(_group_dir(project, rel), "deltas.json")))[rel]
        print(json.dumps(dd, indent=1, ensure_ascii=False))
        return
    dest = sys.argv[3]
    os.makedirs(dest, exist_ok=True)
    for name, data, mode in (("input", text, "w"), ("config.cds.init.001-001", col, "wb")):
        if data is None:
            continue
        p = os.path.join(dest, name)
        if os.path.exists(p):
            sys.exit(f"{p} ya existe: no se sobrescribe")
        open(p, mode).write(data)
        print("escrito", p)


if __name__ == "__main__":
    main()
