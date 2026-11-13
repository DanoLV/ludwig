#!/usr/bin/env python3
"""
Genera SCRIPTS.md en la raíz de Ludwig: qué hace y cómo se usa cada script
(2026-11-13, pedido del usuario).

Fuente: la documentación que trae cada archivo (no se ejecuta nada):
  * Python: docstring del módulo (o comentarios iniciales); "Uso"/"Usage" y argumentos
    de argparse si los hay; sys.argv[...] usados.
  * Shell: bloque de comentarios inicial (tras el shebang).
  * C (util/*.c): bloque de comentario inicial.
Alcance:
  * herramientas generales: util/, code_changes/, scripts de la raíz;
  * scripts de cada proyecto: <proyecto>/scripts/ (desde 2026-11-13; en local/ quedan
    enlaces a estos archivos).
Un script sin documentación queda marcado "SIN DOCUMENTACIÓN EN EL ARCHIVO".

Uso:  python3 util/documentar_scripts.py
"""
import ast
import datetime as dt
import glob
import os
import re

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "SCRIPTS.md")
PROJECTS = ["capacitor", "microgel"]
MAXLINES = 12


def py_doc(path):
    src = open(path, errors="replace").read()
    doc, args, argv = None, [], sorted(set(re.findall(r"sys\.argv\[(\d+)\]", src)))
    try:
        tree = ast.parse(src)
        doc = ast.get_docstring(tree)
        for node in ast.walk(tree):
            if isinstance(node, ast.Call) and getattr(node.func, "attr", "") == "add_argument":
                names = [a.value for a in node.args if isinstance(a, ast.Constant) and isinstance(a.value, str)]
                hlp = next((k.value.value for k in node.keywords if k.arg == "help"
                            and isinstance(k.value, ast.Constant)), "")
                if names:
                    args.append(f"`{' / '.join(names)}`" + (f": {hlp}" if hlp else ""))
    except SyntaxError:
        pass
    if not doc:
        doc = lead_comments(src.splitlines(), "#")
    return doc, args, argv


def lead_comments(lines, mark):
    out = []
    for line in lines:
        s = line.strip()
        if s.startswith("#!"):
            continue
        if s.startswith(mark):
            t = s.lstrip(mark).strip()
            if t and not set(t) <= set("=-*#_ "):
                out.append(t)
            continue
        if not s and not out:
            continue
        break
    return "\n".join(out) or None


def c_doc(path):
    txt = open(path, errors="replace").read()
    m = re.match(r"\s*/\*(.*?)\*/", txt, re.S)
    if not m:
        return None
    lines = [l.strip().lstrip("*").strip() for l in m.group(1).splitlines()]
    lines = [l for l in lines if l and not set(l) <= set("*-=_ ")]
    return "\n".join(lines) or None


def describe(path):
    ext = os.path.splitext(path)[1]
    args, argv = [], []
    if ext == ".py":
        doc, args, argv = py_doc(path)
    elif ext == ".sh" or open(path, "rb").read(2) == b"#!":
        doc = lead_comments(open(path, errors="replace").read().splitlines(), "#")
    elif ext == ".c":
        doc = c_doc(path)
    else:
        doc = None
    return doc, args, argv


def fmt(path, doc, args, argv):
    rel = os.path.relpath(path, ROOT)
    out = [f"#### `{rel}`"]
    if doc:
        lines = doc.strip().splitlines()
        usage = [l for l in lines if re.search(r"\b(uso|usage|python3?\s|\./|bash\s)", l, re.I)]
        body = lines[:MAXLINES]
        out.append("\n".join("> " + l if l.strip() else ">" for l in body))
        if len(lines) > MAXLINES:
            out.append(f"> … ({len(lines) - MAXLINES} líneas más en el archivo)")
        extra = [u for u in usage if u not in body][:4]
        if extra:
            out.append("\nUso (del archivo):\n" + "\n".join(f"    {u.strip()}" for u in extra))
    else:
        out.append("**SIN DOCUMENTACIÓN EN EL ARCHIVO.**")
    if args:
        out.append("\nArgumentos (argparse): " + "; ".join(args))
    elif argv:
        out.append(f"\nLee argumentos de línea de comandos: sys.argv[{', '.join(argv)}]")
    return "\n".join(out) + "\n"


def collect():
    sections = []
    gen = sorted(glob.glob(os.path.join(ROOT, "util", "*.py")) + glob.glob(os.path.join(ROOT, "util", "*.sh"))
                 + glob.glob(os.path.join(ROOT, "code_changes", "*.py")) + glob.glob(os.path.join(ROOT, "*.sh")))
    sections.append(("Herramientas generales (Python / shell)", gen))
    sections.append(("Programas en C de util/ (se compilan con util/Makefile)",
                     sorted(glob.glob(os.path.join(ROOT, "util", "*.c")))))
    for p in PROJECTS:
        sdir = os.path.join(ROOT, p, "scripts")
        files = []
        for r, ds, fs in os.walk(sdir):
            ds[:] = sorted(x for x in ds if x != "__pycache__")
            files += sorted(os.path.join(r, f) for f in fs if f.endswith((".py", ".sh")))
        sections.append((f"Proyecto {p} ({p}/scripts/)", files))
    return sections


def main():
    sections = collect()
    out = ["# Scripts de Ludwig: qué hace cada uno y cómo se usa\n",
           f"Generado por `util/documentar_scripts.py` el {dt.datetime.now():%Y-%m-%d %H:%M} a partir de la "
           "documentación que trae cada archivo (docstring o comentarios iniciales). No editar a mano: "
           "documentar el script en su propio archivo y volver a correr el generador.\n",
           "Alcance: `util/`, `code_changes/`, scripts de la raíz y `<proyecto>/scripts/` (scripts de "
           "cada proyecto; en `local/` quedan enlaces a ellos, ver `<proyecto>/scripts/README_ubicacion.md`). "
           "No incluye las copias de scripts dentro de cada corrida.\n"]
    summary, body = [], []
    for title, files in sections:
        docs = [(f, *describe(f)) for f in files]
        n_no = sum(1 for d in docs if not d[1])
        summary.append(f"| {title} | {len(files)} | {n_no} |")
        body.append(f"\n## {title}\n")
        cur = None
        for f, doc, args, argv in docs:
            folder = os.path.relpath(os.path.dirname(f), ROOT)
            if folder != cur:
                body.append(f"\n### `{folder}/`\n")
                cur = folder
            body.append(fmt(f, doc, args, argv))
    out.append("| Sección | Scripts | Sin documentación |\n|---|---|---|")
    out += summary
    out += body
    open(OUT, "w", encoding="utf-8").write("\n".join(out) + "\n")
    tot = sum(len(f) for _, f in sections)
    print(f"{tot} scripts -> {OUT}")


if __name__ == "__main__":
    main()
