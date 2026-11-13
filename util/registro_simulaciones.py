#!/usr/bin/env python3
"""
Registro de simulaciones por proyecto (2026-11-13, pedido del usuario).

Genera <proyecto>/Registro_simulaciones.xlsx para capacitor/ y microgel/:
  * una pestaña por SERIE (agrupación científica de carpetas, ver SERIES_*),
  * dentro de cada pestaña, una fila-título por carpeta-serie y una fila por
    simulación (= cada carpeta que contiene un archivo "input"),
  * columnas: ruta, descripción, fecha, ubicación (local / runs / borrada),
    ID del archivo de coloides y los parámetros del input que varían en la
    pestaña; los que no varían se listan en una fila "Comunes" arriba,
  * hoja "Coloides": un ID por contenido distinto de config.cds.init.001-001,
  * hoja "Leyenda".

Se puede volver a correr cuando haya corridas nuevas:
    python3 util/registro_simulaciones.py            # ambos proyectos
    python3 util/registro_simulaciones.py microgel   # uno solo
Las descripciones editadas a mano se conservan: si en el Excel existente la
columna "Descripción" difiere de la columna oculta "Descripción (auto)", se
mantiene la versión editada (clave = pestaña no importa, se busca por Ruta).
Los IDs de coloides se reasignan en cada corrida del script (orden por fecha).
No modifica nada en local/ ni en runs/; solo escribe el .xlsx.
"""
import collections
import datetime as dt
import glob
import hashlib
import json
import os
import re
import sys

import openpyxl
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
from openpyxl.utils import get_column_letter
from openpyxl.cell.cell import ILLEGAL_CHARACTERS_RE

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from registro_propositos import PROPOSITOS, CONCLUSIONES  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
XLSX = "Registro_simulaciones.xlsx"
COLLOID_FILE = "config.cds.init.001-001"
ALWAYS = ["size", "N_cycles", "temperature", "viscosity", "subgrid_kernel",
          "subgrid_hann_order", "electrokinetics_solver_type"]
GENERIC_COMMENT = ("Microgel polymer test", "Monomers are subgrid particles")
TEMPLATE_TITLES = {"System and MPI", "Fluid parameters", "Sistema y MPI", "Fluido", "Corrida"}


# --------------------------------------------------------------------------
# Series (pestañas) y carpetas-serie (filas-título). Acordado con el usuario
# el 2026-11-13. Devuelven (pestaña, grupo) a partir de la ruta relativa.
# --------------------------------------------------------------------------
def series_microgel(rel):
    parts = rel.split("/") if rel else [""]
    top = parts[0]
    if top == "old":
        sub = parts[1] if len(parts) > 1 else ""
        if sub.startswith(("PM-FFT", "PETSC", "test_size")):
            return "Historico (old)", "old/PM-FFT + PETSC + test_size"
        return "Historico (old)", "old/" + sub
    rules = [
        (("FFT-Force-", "SLIDES-", "PETSC-LD8-", "TEST-HIDRO-"), "Fuerza 2 particulas"),
        (("SS-SCAN-",), "Repulsion corto alcance"),
        (("NBODY-",), "Aditividad N cuerpos"),
        (("PIN2GRID-",), "1 part - autofuerza vs pos"),
        (("EFIELD-1PART-",), "1 part - campo electrico"),
        (("RING-", "CROSS-"), "Anillo y borde"),
        (("MGEL-",), "Microgel"),
        (("INGO-",), "Ingo 1 part LB termico"),
    ]
    for prefixes, tab in rules:
        if top.startswith(prefixes):
            return tab, top
    return "Varios", top or "(raiz de local)"


def series_capacitor(rel):
    parts = rel.split("/")
    if parts[0] != "WALL-DIRICHLET":
        if parts[0].startswith("SOLIDO"):
            return "Relajacion microgel-solido", parts[0]
        return "Varios", parts[0]
    sub = parts[1] if len(parts) > 1 else ""
    group = "WALL-DIRICHLET/" + sub
    if sub.startswith(("A-", "B-", "D-", "E-", "regresion")):
        return "Capacitor sin particula", group
    if sub.startswith("C-electrolito"):
        return "C particula + electrolito", group
    if sub.startswith("C-particula-libre"):
        return "C particula libre", group
    if sub.startswith("C-"):
        return "C particula fija", group
    if sub.startswith("F-"):
        return "F campo tangencial", group
    if sub.startswith("G-"):
        return "Relajacion microgel-solido", group
    return "Varios", group


SERIES = {"microgel": series_microgel, "capacitor": series_capacitor}
TAB_ORDER = {
    "microgel": ["Fuerza 2 particulas", "Repulsion corto alcance", "Aditividad N cuerpos",
                 "1 part - autofuerza vs pos", "1 part - campo electrico", "Anillo y borde",
                 "Microgel", "Ingo 1 part LB termico", "Historico (old)", "Varios"],
    "capacitor": ["Capacitor sin particula", "C particula fija", "C particula + electrolito",
                  "C particula libre", "F campo tangencial", "Relajacion microgel-solido", "Varios"],
}


# --------------------------------------------------------------------------
def find_runs(project):
    """{rel_dir: {"local": abs or None, "runs": abs or None}}"""
    runs = {}
    for where in ("local", "runs"):
        base = os.path.join(ROOT, project, where)
        for dirpath, dirnames, filenames in os.walk(base):
            if "input" in filenames and os.path.isfile(os.path.join(dirpath, "input")):
                rel = os.path.relpath(dirpath, base)
                rel = "" if rel == "." else rel
                runs.setdefault(rel, {"local": None, "runs": None})[where] = dirpath
    return runs


def parse_input(path):
    params, header, hdr_lines, past_header = {}, None, [], False
    if b"\0" in open(path, "rb").read(4096):
        return {}, "INPUT ILEGIBLE (archivo binario/corrupto, no es texto)"
    with open(path, errors="replace") as fh:
        for i, line in enumerate(fh):
            s = line.strip()
            if s.startswith("#"):
                # cabecera: lineas de comentario con texto real (no '=====', no la generica)
                txt = s.strip("#=-*_ \t")
                if (not past_header and len(txt) > 8 and re.search(r"[A-Za-z]", txt)
                        and not any(g in txt for g in GENERIC_COMMENT)
                        and txt not in TEMPLATE_TITLES):
                    hdr_lines.append(txt)
                continue
            s = s.split("#", 1)[0].strip()
            if not s:
                continue
            past_header = True      # la cabecera es solo el bloque de comentarios inicial
            tok = s.split(None, 1)
            params[tok[0]] = tok[1].strip() if len(tok) > 1 else ""
    header = " ".join(hdr_lines[:4]) or None
    return params, header


DATE_RE = re.compile(r"Start time:\s+(.*\d{4})")


def run_date(dirs):
    """Fecha de inicio desde run.log / logs/*.txt / *.log; si no, mtime del input."""
    for d in dirs:
        cands = [os.path.join(d, "run.log")] + sorted(glob.glob(os.path.join(d, "*.log"))) \
            + sorted(glob.glob(os.path.join(d, "logs", "*")))
        for c in cands:
            if not os.path.isfile(c) or os.path.getsize(c) > 50e6:
                continue
            with open(c, errors="replace") as fh:
                for k, line in enumerate(fh):
                    m = DATE_RE.search(line)
                    if m:
                        try:
                            t = dt.datetime.strptime(m.group(1).strip(), "%a %b %d %H:%M:%S %Y")
                            return t, "Start time (" + os.path.basename(c) + ")"
                        except ValueError:
                            pass
                    if k > 400:
                        break
    t = dt.datetime.fromtimestamp(os.path.getmtime(os.path.join(dirs[0], "input")))
    return t, "fecha de modificacion del input"


def colloid_info(path):
    raw = open(path, "rb").read()
    info = {"md5": hashlib.md5(raw).hexdigest(), "N": None, "q0": None, "q1": None,
            "a0": "", "ah": ""}
    try:
        lines = raw.decode().splitlines()
        n = int(lines[0].split()[0])
        q0 = q1 = 0.0
        a0s, ahs = set(), set()
        pos, fixed, qs = [], [], []
        for c in range(n):
            b = 1 + 68 * c          # first line of colloid c (0-based)
            if c < 4:
                pos.append([float(x) for x in lines[b + 34].split()[:3]])
                qs.append(float(lines[b + 45].split()[0]))
            fixed.append(int(float(lines[b + 4].split()[0])))
            a0s.add(round(float(lines[b + 32].split()[0]), 6))
            ahs.add(round(float(lines[b + 33].split()[0]), 6))
            q0 += float(lines[b + 45].split()[0])
            q1 += float(lines[b + 46].split()[0])
        info.update(N=n, q0=q0, q1=q1, pos=pos, fixed=fixed, qs=qs,
                    a0=", ".join(f"{v:g}" for v in sorted(a0s)),
                    ah=", ".join(f"{v:g}" for v in sorted(ahs)))
    except Exception:
        info["N"] = "no legible (binario?)"
    return info


def group_title(project, group):
    for where in ("local", "runs"):
        d = os.path.join(ROOT, project, where, group.split(" + ")[0])
        for r in sorted(glob.glob(os.path.join(d, "README*"))):
            with open(r, errors="replace") as fh:
                for line in fh:
                    t = line.strip().strip("#").strip()
                    if t:
                        return t[:250]
    return ""


def group_readmes(project, group):
    out = []
    for where in ("local", "runs"):
        d = os.path.join(ROOT, project, where, group.split(" + ")[0])
        out += sorted(glob.glob(os.path.join(d, "README*")))
        if out:
            break
    return out


def colloid_case(inf):
    """Texto del caso a partir del archivo de coloides (posiciones, separaciones, fijas/libres)."""
    if not inf or not isinstance(inf.get("N"), int):
        return ""
    n, pos, fx = inf["N"], inf.get("pos", []), inf.get("fixed", [])
    nf = sum(1 for f in fx if f)
    estado = "fija" if (n == 1 and nf == 1) else "libre" if (n == 1) else \
        f"{nf} fijas, {n - nf} libres" if 0 < nf < n else ("todas fijas" if nf else "todas libres")
    if n == 0:
        return "sin partículas"
    if n <= 4:
        ps = "; ".join(f"q={q:g} en ({x:.2f}, {y:.2f}, {z:.2f})" for (x, y, z), q in zip(pos, inf["qs"]))
        txt = f"{n} partícula{'s' if n > 1 else ''} ({estado}): {ps}"
        if n >= 2:
            ds = [f"r{i+1}{j+1} = {sum((a - b) ** 2 for a, b in zip(pos[i], pos[j])) ** 0.5:.3f}"
                  for i in range(n) for j in range(i + 1, n)]
            txt += "; separación " + ", ".join(ds)
        return txt
    return f"{n} partículas ({estado}), Σq0 = {inf['q0']:g}"


def readme_case(readmes, rel):
    """Linea del README de la carpeta-serie que empieza con el nombre de la corrida."""
    name = os.path.basename(rel)
    if not name:
        return None
    pat = re.compile(r"^\s*" + re.escape(name) + r"(\s{2,}|\s+[-=]|\s*[:(])\s*(.+)$")
    for r in readmes:
        with open(r, errors="replace") as fh:
            for line in fh:
                m = pat.match(line.rstrip())
                if m:
                    return m.group(2).strip()
    return None


def purpose(project, group):
    for p, rx, idea, src in PROPOSITOS:
        if p == project and re.search(rx, group):
            return idea, src
    return "Sin documentar: completar a mano.", "—"


def conclusion(project, group):
    for p, rx, txt, src in CONCLUSIONES:
        if p == project and re.search(rx, group):
            return f"{txt} [fuente: {src}]"
    return ""


EDITABLE = {"Descripción": "Descripción (auto)", "Conclusión": "Conclusión (auto)"}


def old_descriptions(path):
    """{(columna, ruta): texto editado a mano} para las columnas editables."""
    if not os.path.exists(path):
        return {}
    keep = {}
    wb = openpyxl.load_workbook(path)
    for ws in wb.worksheets:
        hdr = [c.value for c in ws[1]]
        if "Ruta" not in hdr:
            continue
        ir = hdr.index("Ruta")
        for col, auto in EDITABLE.items():
            if col not in hdr or auto not in hdr:
                continue
            idd, ia = hdr.index(col), hdr.index(auto)
            for row in ws.iter_rows(min_row=2, values_only=True):
                if row[ir] and row[ia] is not None and row[idd] is not None and row[idd] != row[ia]:
                    keep[(col, row[ir])] = row[idd]
    return keep


# --------------------------------------------------------------------------
NON_PARAM = {"Ruta", "Descripción", "Conclusión", "Fuente de la idea", "Fecha", "Fuente de la fecha",
             "Ubicación", "Coloides (ID)", "N coloides", "Rehacer desde saved_data", "Descripción (auto)", "Conclusión (auto)"}
GONE = "BORRADA (no está en local ni en runs)"


def old_rows(path):
    """Filas de simulación del Excel anterior: {ruta: {columna: valor}} (para no perder las borradas)."""
    if not os.path.exists(path):
        return {}
    rows = {}
    wb = openpyxl.load_workbook(path)
    for ws in wb.worksheets:
        hdr = [c.value for c in ws[1]]
        if "Ruta" not in hdr or "Ubicación" not in hdr:
            continue
        for vals in ws.iter_rows(min_row=2, values_only=True):
            d = dict(zip(hdr, vals))
            if d.get("Ruta") and d.get("Ubicación"):
                rel = "" if d["Ruta"] == "(raiz)" else d["Ruta"]
                rows[rel] = d
    return rows


_DELTAS = {}


def rehacer_info(project, group, rel):
    """Texto: de qué base de saved_data/ (util/guardar_bases.py) se rehace la corrida."""
    gname = group.replace(" + ", "__")
    path = os.path.join(ROOT, project, "saved_data", gname, "deltas.json")
    if path not in _DELTAS:
        _DELTAS[path] = json.load(open(path)) if os.path.exists(path) else {}
    dd = _DELTAS[path].get(rel)
    if dd is None:
        return "no guardada (correr util/guardar_bases.py)"
    n_in = len(dd["input_set"]) + len(dd["input_remove"])
    txt = f"saved_data/{gname}: input.base + {n_in} cambio(s)"
    if "coloides_base" in dd:
        txt += f"; {dd['coloides_base']} + {len(dd['coloides_lineas'])} línea(s)"
    return txt


def build(project):
    out = os.path.join(ROOT, project, XLSX)
    edited = old_descriptions(out)
    previous = old_rows(out)
    runs = find_runs(project)
    recs = []
    for rel, loc in runs.items():
        dirs = [d for d in (loc["local"], loc["runs"]) if d]
        params, header = parse_input(os.path.join(dirs[0], "input"))
        date, dsrc = run_date(dirs)
        if loc["local"] and loc["runs"]:
            where = "local + runs"
        elif loc["runs"]:
            where = "solo runs (borrada de local)"
        else:
            where = "solo local (falta espejar en runs)"
        cf = next((os.path.join(d, COLLOID_FILE) for d in dirs
                   if os.path.isfile(os.path.join(d, COLLOID_FILE))), None)
        tab, group = SERIES[project](rel)
        recs.append(dict(rel=rel, params=params, header=header, date=date, dsrc=dsrc,
                         where=where, cf=cf, tab=tab, group=group))
    # corridas registradas antes cuya carpeta ya no existe: se conservan, marcadas como borradas
    for rel, d in previous.items():
        if rel in runs:
            continue
        try:
            date = dt.datetime.strptime(str(d.get("Fecha")), "%Y-%m-%d %H:%M")
        except ValueError:
            date = dt.datetime(1970, 1, 1)
        tab, group = SERIES[project](rel)
        params = {k: str(v) for k, v in d.items() if k not in NON_PARAM and v not in (None, "—")}
        recs.append(dict(rel=rel, params=params, header=None, date=date,
                         dsrc=d.get("Fuente de la fecha") or "", cf=None, tab=tab, group=group,
                         where=GONE if not str(d.get("Ubicación", "")).startswith("BORRADA") else d["Ubicación"],
                         ghost=d))
    recs.sort(key=lambda r: r["date"])

    # colloid IDs by content, in order of first use
    cinfo, cid = {}, {}
    for r in recs:
        if r.get("ghost"):
            r["cid"] = "archivo borrado"
            r["cN"] = r["ghost"].get("N coloides") or ""
            continue
        if not r["cf"]:
            r["cid"] = "sin archivo"
            continue
        inf = colloid_info(r["cf"])
        if inf["md5"] not in cid:
            cid[inf["md5"]] = f"C{len(cid) + 1:03d}"
            inf.update(id=cid[inf["md5"]], first=r["rel"], nuse=0)
            cinfo[inf["md5"]] = inf
        cinfo[inf["md5"]]["nuse"] += 1
        r["cid"] = cid[inf["md5"]]
        r["cinf"] = inf
        r["cN"] = cinfo[inf["md5"]]["N"]

    wb = openpyxl.Workbook()
    wb.remove(wb.active)
    F = Font(name="Arial", size=10)
    FB = Font(name="Arial", size=10, bold=True)
    FH = Font(name="Arial", size=10, bold=True, color="FFFFFF")
    thin = Side(style="thin", color="BFBFBF")
    B = Border(left=thin, right=thin, top=thin, bottom=thin)
    fill_h = PatternFill("solid", fgColor="2F5496")
    fill_g = PatternFill("solid", fgColor="D9E1F2")
    fill_c = PatternFill("solid", fgColor="FFF2CC")
    fill_del = PatternFill("solid", fgColor="F8CBAD")

    tabs = collections.OrderedDict((t, []) for t in TAB_ORDER[project])
    for r in recs:
        tabs.setdefault(r["tab"], []).append(r)
    summary = []
    for tab, rs in tabs.items():
        if not rs:
            continue
        ws = wb.create_sheet(tab[:31])
        keys = sorted({k for r in rs for k in r["params"]})
        vals = {k: {r["params"].get(k, "—") for r in rs} for k in keys}
        varying = [k for k in ALWAYS if k in keys] + \
                  [k for k in keys if len(vals[k]) > 1 and k not in ALWAYS]
        common = [k for k in keys if len(vals[k]) == 1 and k not in ALWAYS]
        base = ["Ruta", "Descripción", "Conclusión", "Fuente de la idea", "Fecha", "Fuente de la fecha", "Ubicación",
                "Coloides (ID)", "N coloides", "Rehacer desde saved_data"]
        hdr = base + varying + ["Descripción (auto)", "Conclusión (auto)"]
        for j, h in enumerate(hdr, 1):
            c = ws.cell(1, j, h)
            c.font, c.fill, c.border = FH, fill_h, B
            c.alignment = Alignment(wrap_text=True, vertical="center", horizontal="center")
        row = 2
        txt = "; ".join(f"{k} = {next(iter(vals[k]))}" for k in common)
        ws.cell(row, 1, "Comunes a toda la pestaña").font = FB
        c = ws.cell(row, 2, txt or "(ninguno)")
        c.font, c.alignment = F, Alignment(wrap_text=True, vertical="top")
        ws.merge_cells(start_row=row, start_column=2, end_row=row, end_column=min(len(hdr), 7))
        for j in range(1, len(hdr) + 1):
            ws.cell(row, j).fill = fill_c
        ws.row_dimensions[row].height = min(400, 15 * (1 + len(txt) // 180))
        row += 1
        groups = collections.OrderedDict()
        for r in sorted(rs, key=lambda r: (min(x["date"] for x in rs if x["group"] == r["group"]),
                                           r["group"], r["date"])):
            groups.setdefault(r["group"], []).append(r)
        for g, grs in groups.items():
            idea, isrc = purpose(project, g)
            readmes = group_readmes(project, g)
            concl = conclusion(project, g)
            d0, d1 = min(r["date"] for r in grs), max(r["date"] for r in grs)
            title = f"{g}   ({len(grs)} sim., {d0:%Y-%m-%d}" + \
                    (f" a {d1:%Y-%m-%d})" if d1.date() != d0.date() else ")") + \
                    f" — {idea}"
            c = ws.cell(row, 1, title)
            c.font = FB
            for j in range(1, len(hdr) + 1):
                ws.cell(row, j).fill = fill_g
            row += 1
            # description: input header line + keys differing from the group's mode
            mode = {}
            for k in varying:
                cnt = collections.Counter(r["params"].get(k, "—") for r in grs)
                mode[k] = cnt.most_common(1)[0][0]
            for r in grs:
                diff = [f"{k}={r['params'].get(k, '—')}" for k in varying
                        if r["params"].get(k, "—") != mode[k] and len(grs) > 1]
                rcase = readme_case(readmes, r["rel"])
                ccase = colloid_case(r.get("cinf"))
                if rcase:
                    case = rcase
                elif r["header"]:
                    case = r["header"][:400]
                else:
                    case = (os.path.basename(r["rel"]) or r["rel"]) + (": " + ccase if ccase else "")
                if diff:
                    case += "  [parámetros distintivos: " + "; ".join(diff[:6]) + (" …" if len(diff) > 6 else "") + "]"
                auto = f"IDEA: {idea}\nCASO: {case}"
                src = isrc + (" / caso: README" if rcase else " / caso: cabecera input" if r["header"] else " / caso: archivo de coloides y parámetros")
                if r.get("ghost"):      # fila conservada de una corrida borrada: textos tal como estaban
                    g = r["ghost"]
                    auto = g.get("Descripción (auto)") or g.get("Descripción") or auto
                    src = g.get("Fuente de la idea") or src
                desc = edited.get(("Descripción", r["rel"]), auto)
                conc = edited.get(("Conclusión", r["rel"]), concl)
                line = [r["rel"] or "(raiz)", desc, conc, src, r["date"].strftime("%Y-%m-%d %H:%M"),
                        r["dsrc"], r["where"], r["cid"], r.get("cN", ""),
                        (r["ghost"].get("Rehacer desde saved_data") or "") if r.get("ghost")
                        else rehacer_info(project, r["group"], r["rel"])] + \
                       [r["params"].get(k, "—") for k in varying] + [auto, concl]
                for j, v in enumerate(line, 1):
                    if isinstance(v, str):
                        v = ILLEGAL_CHARACTERS_RE.sub("?", v)
                    c = ws.cell(row, j, v)
                    c.font, c.border = F, B
                    c.alignment = Alignment(wrap_text=(j in (2, 3)), vertical="top")
                if r["where"].startswith(("solo runs", "BORRADA")):
                    ws.cell(row, 7).fill = fill_del
                row += 1
        widths = [45, 90, 70, 22, 16, 22, 22, 11, 9, 40] + [14] * len(varying) + [30]
        for j, w in enumerate(widths, 1):
            ws.column_dimensions[get_column_letter(j)].width = w
        ws.column_dimensions[get_column_letter(len(hdr))].hidden = True
        ws.column_dimensions[get_column_letter(len(hdr) - 1)].hidden = True
        ws.freeze_panes = "C2"
        ws.auto_filter.ref = f"A1:{get_column_letter(len(hdr))}{row - 1}"
        summary.append((tab, len(rs), len(groups), len(varying)))

    wc = wb.create_sheet("Coloides")
    ch = ["ID", "md5 (8)", "N coloides", "Σ q0", "Σ q1", "a0 (valores)", "ah (valores)",
          "N simulaciones", "Primera simulación que lo usa"]
    for j, h in enumerate(ch, 1):
        c = wc.cell(1, j, h)
        c.font, c.fill = FH, fill_h
    for i, inf in enumerate(sorted(cinfo.values(), key=lambda x: x["id"]), 2):
        for j, v in enumerate([inf["id"], inf["md5"][:8], inf["N"], inf["q0"], inf["q1"],
                               inf["a0"], inf["ah"], inf["nuse"], inf["first"]], 1):
            wc.cell(i, j, v).font = F
    for j, w in enumerate([8, 11, 11, 10, 10, 18, 18, 14, 80], 1):
        wc.column_dimensions[get_column_letter(j)].width = w
    wc.freeze_panes = "B2"
    wc.auto_filter.ref = f"A1:I{len(cinfo) + 1}"

    wl = wb.create_sheet("Leyenda")
    legend = [
        ("Qué es", f"Registro de todas las simulaciones de {project}/ (carpetas con un archivo 'input' en "
                   f"{project}/local/ o {project}/runs/). Generado por util/registro_simulaciones.py el "
                   f"{dt.datetime.now():%Y-%m-%d %H:%M}. Volver a correrlo tras corridas nuevas."),
        ("Pestañas", "Una por serie científica (agrupación de carpetas acordada el 2026-11-13; ver SERIES_* "
                     "en el script). Fila amarilla: parámetros comunes a toda la pestaña. Filas azules: "
                     "carpeta-serie (n.º de simulaciones, fechas y primera línea de su README si existe)."),
        ("Ruta", f"Carpeta de la simulación relativa a {project}/local/ (o a {project}/runs/)."),
        ("Descripción", "IDEA: qué se intenta simular en la carpeta-serie (texto de util/registro_propositos.py). "
                        "CASO: qué distingue a esta corrida (línea del README de la carpeta con su nombre, o la cabecera "
                        "del input, o el nombre de la carpeta) + parámetros que difieren del valor más frecuente en la "
                        "carpeta (hasta 6). Editable: al regenerar se conserva lo editado (se compara con la columna "
                        "oculta 'Descripción (auto)')."),
        ("Conclusión", "Qué se puede concluir de la carpeta-serie, solo donde hay un resultado registrado (README, diario, "
                       "memoria o sesiones previas; la fuente va entre corchetes). Vacía = sin conclusión registrada: completar a mano. "
                       "Texto en util/registro_propositos.py (CONCLUSIONES); lo editado en el Excel se conserva al regenerar."),
        ("Fuente de la idea", "De dónde sale la IDEA: README, cabecera input, sesiones (conversaciones previas con Claude) "
                              "o INFERIDO (deducido del nombre y los parámetros: verificar)."),
        ("Fecha", "Start time del run.log (o de *.log / logs/*). Si no hay log (p.ej. corridas solo en runs/, "
                  "donde no se copian los .log), fecha de modificación del input: ver 'Fuente de la fecha'."),
        ("Ubicación", "local + runs | solo runs (borrada de local; fondo rojo) | BORRADA (no está en local ni en "
                      "runs: la fila se conserva del registro anterior con sus valores; fondo rojo) | solo local (falta correr "
                      "util/archive_runs.py)."),
        ("Rehacer desde saved_data", "Base de <proyecto>/saved_data/<carpeta>/ (input.base, coloides_k.cds) y número de "
                                     "cambios guardados en deltas.json. Reconstruir con: python3 util/rehacer_corrida.py "
                                     "<proyecto> <Ruta> <carpeta destino>. Bases y deltas: python3 util/guardar_bases.py."),
        ("Coloides (ID)", "Un ID por contenido distinto de config.cds.init.001-001 (hash md5); mismas ID = "
                          "mismo archivo exacto. Detalle en la hoja 'Coloides'. Los ID se reasignan al regenerar."),
        ("Parámetros", "Columnas = claves del input que varían dentro de la pestaña, más siempre: "
                       + ", ".join(ALWAYS) + ". '—' = clave ausente en ese input (vale el default de Ludwig). "
                       "Si una clave aparece dos veces en un input se toma la última."),
        ("Coloides: Σ q0, Σ q1", "Suma sobre coloides de las líneas q0 y q1 del archivo (líneas 47 y 48 del "
                                 "primer coloide; 68 líneas por coloide). No se interpreta el signo."),
    ]
    for i, (a, b) in enumerate(legend, 1):
        wl.cell(i, 1, a).font = FB
        c = wl.cell(i, 2, b)
        c.font, c.alignment = F, Alignment(wrap_text=True, vertical="top")
    wl.column_dimensions["A"].width = 20
    wl.column_dimensions["B"].width = 120

    wb.save(out)
    n_where = collections.Counter(r["where"] for r in recs)
    print(f"{project}: {len(recs)} simulaciones -> {out}")
    for t in summary:
        print(f"   {t[0]:30s} {t[1]:5d} sim  {t[2]:3d} carpetas  {t[3]:3d} col. de parámetros")
    print("   ubicación:", dict(n_where), "| coloides distintos:", len(cinfo),
          "| ediciones conservadas:", len(edited))


if __name__ == "__main__":
    for p in (sys.argv[1:] or ["capacitor", "microgel"]):
        build(p)
