#!/usr/bin/env python3
"""
Genera code_changes/CAMBIOS_CODIGO.md a partir de las marcas de cambio local
del código de Ludwig (convención del proyecto, ver CLAUDE.md):

    /*CHANGE INIT - <AAAAMMDD> <descripción> */  ...  /*CHANGE END - ... */

Para cada marca: fecha, descripción, archivo y rango de líneas (INIT..END).
Agrupa por fecha (más reciente primero) y, dentro de cada fecha, por descripción.
Las entradas escritas a mano en code_changes/NOTAS_MANUALES.md (si existe) se
copian al principio del documento: sirven para explicar el porqué de un cambio
con más detalle del que cabe en la marca.

Uso (desde la raíz de Ludwig o desde cualquier lugar):
    python3 code_changes/extraer_cambios.py
Solo lee el código; escribe únicamente code_changes/CAMBIOS_CODIGO.md.
"""
import collections
import datetime as dt
import os
import re
import subprocess

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "code_changes", "CAMBIOS_CODIGO.md")
MANUAL = os.path.join(ROOT, "code_changes", "NOTAS_MANUALES.md")
DIRS = ["src", "target", "tests", "util", "mpi_s", "config"]
EXTS = (".c", ".h", ".cu", ".cpp", ".hpp", ".mk")
INIT_RE = re.compile(r"CHANGE\s+INIT\s*-?\s*(.*)")
END_RE = re.compile(r"CHANGE\s+END")
DATE_RE = re.compile(r"^\s*(20\d{6})\b[\s:,-]*(.*)$")


def clean(text):
    text = text.split("*/")[0]
    return " ".join(text.replace("/*", " ").replace("*", " ").split()).strip(" -:")


def scan():
    marks = []
    for d in DIRS:
        base = os.path.join(ROOT, d)
        for dirpath, dirnames, filenames in os.walk(base):
            dirnames[:] = [x for x in dirnames if not x.startswith(".")]
            for fn in filenames:
                if not (fn.endswith(EXTS) or fn.startswith("Makefile")):
                    continue
                path = os.path.join(dirpath, fn)
                rel = os.path.relpath(path, ROOT)
                with open(path, errors="replace") as fh:
                    lines = fh.readlines()
                for i, line in enumerate(lines):
                    m = INIT_RE.search(line)
                    if not m:
                        continue
                    text = m.group(1)
                    # descripción en varias líneas: seguir hasta cerrar el comentario (máx. 4)
                    j = i
                    while "*/" not in text and j + 1 < len(lines) and j < i + 4:
                        j += 1
                        text += " " + lines[j].strip()
                    text = clean(text)
                    dm = DATE_RE.match(text)
                    date, desc = (dm.group(1), dm.group(2).strip()) if dm else ("sin fecha", text)
                    end = next((k for k in range(i + 1, len(lines)) if END_RE.search(lines[k])), None)
                    marks.append(dict(date=date, desc=desc or "(sin descripción)", file=rel, raw=line.strip(),
                                      git=False,
                                      l0=i + 1, l1=(end + 1) if end is not None else None))
    return marks


def git_date(raw, rel):
    """Fecha (AAAAMMDD) del primer commit que introdujo esa línea de marca en ese archivo, o None."""
    try:
        r = subprocess.run(["git", "-C", ROOT, "log", "--reverse", "--format=%ad", "--date=format:%Y%m%d",
                            "-S" + raw, "--", rel], capture_output=True, text=True, timeout=120)
        out = r.stdout.split()
        return out[0] if out else None
    except Exception:
        return None


def date_undated(marks):
    """Marcas sin fecha: se fechan con git (fecha del commit, no de la escritura del código)."""
    cache = {}
    for m in marks:
        if m["date"] != "sin fecha":
            continue
        key = (m["raw"], m["file"])
        if key not in cache:
            cache[key] = git_date(*key)
        if cache[key]:
            m["date"], m["git"] = cache[key], True
        else:
            m["date"] = "sin commitear"


def fmt_date(d):
    return f"{d[:4]}-{d[4:6]}-{d[6:]}" if d[:1].isdigit() else d


def main():
    marks = scan()
    n_undated = sum(1 for m in marks if m["date"] == "sin fecha")
    date_undated(marks)
    by_date = collections.defaultdict(list)
    for m in marks:
        by_date[m["date"]].append(m)
    dates = sorted((d for d in by_date if d[:1].isdigit()), reverse=True)
    dates += [d for d in ("sin commitear", "sin fecha") if d in by_date]

    out = []
    out.append("# Cambios de código en Ludwig (registro)\n")
    out.append(f"Generado por `code_changes/extraer_cambios.py` el {dt.datetime.now():%Y-%m-%d %H:%M} "
               "a partir de las marcas `/*CHANGE INIT - <fecha> <descripción> */ … /*CHANGE END */` del código "
               f"(carpetas: {', '.join(DIRS)}). No editar a mano: volver a correr el script. "
               "Las explicaciones largas van en `code_changes/NOTAS_MANUALES.md`, que se copia abajo.\n")
    out.append(f"**{len(marks)} bloques de cambio** en **{len({m['file'] for m in marks})} archivos**, "
               f"{len(dates)} fechas.\n")
    n_git = sum(1 for m in marks if m["git"])
    out.append(f"Fechas: {len(marks) - n_undated} bloques tienen la fecha en la marca; {n_git} marcas sin fecha se "
               "fecharon con git (primer commit que introdujo esa línea; señaladas *(fecha de git)*: es la fecha "
               "del commit, el código pudo escribirse antes)"
               + (f"; {sum(1 for m in marks if m['date'] == 'sin commitear')} sin fecha y todavía sin commitear" if "sin commitear" in by_date else "")
               + ".\n")
    out.append("Limitación: solo registra cambios marcados con `CHANGE INIT`. Un cambio sin marca no aparece "
               "(para eso, `git log` / `git diff`).\n")

    if os.path.exists(MANUAL):
        out.append("\n---\n")
        out.append(open(MANUAL, encoding="utf-8").read().rstrip() + "\n")

    out.append("\n---\n\n## Índice por fecha\n")
    out.append("| Fecha | Bloques | Archivos | Descripciones (resumen) |")
    out.append("|---|---|---|---|")
    for d in dates:
        ms = by_date[d]
        files = sorted({m["file"] for m in ms})
        descs = list(dict.fromkeys(m["desc"] for m in ms))
        short = [x if len(x) <= 90 else x[:87] + "…" for x in descs]
        summ = "; ".join(short[:3]) + (f" (+{len(descs) - 3} más)" if len(descs) > 3 else "")
        out.append(f"| [{fmt_date(d)}](#{fmt_date(d)}) | {len(ms)} | {len(files)} | {summ.replace('|', '/')} |")

    out.append("\n## Detalle\n")
    for d in dates:
        ms = by_date[d]
        out.append(f"### {fmt_date(d)}\n")
        out.append("Archivos: " + ", ".join(f"`{f}`" for f in sorted({m['file'] for m in ms})) + "\n")
        by_desc = collections.OrderedDict()
        for m in sorted(ms, key=lambda m: (m["file"], m["l0"])):
            by_desc.setdefault(m["desc"], []).append(m)
        for desc, lst in by_desc.items():
            locs = ", ".join(f"`{m['file']}:{m['l0']}" + (f"-{m['l1']}`" if m["l1"] else "` (sin CHANGE END)")
                             for m in lst)
            tag = " *(fecha de git)*" if any(m["git"] for m in lst) else ""
            out.append(f"- **{desc}**{tag} — {locs}")
        out.append("")

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", encoding="utf-8") as fh:
        fh.write("\n".join(out) + "\n")
    print(f"{len(marks)} bloques, {len(dates)} fechas -> {OUT}")


if __name__ == "__main__":
    main()
