#!/usr/bin/env python3
"""
Para cada proyecto de la forma <proyecto>/local/ en la raiz de Ludwig (p.ej.
capacitor/, microgel/), crea/actualiza <proyecto>/runs/ como un espejo
liviano de <proyecto>/local/: copia la estructura completa de carpetas, pero
excluyendo los volcados de datos pesados de Ludwig (dist-/psi-/qsi-/vel-/
rho-/efield-/config.cds-<paso>, sus "-metadata", binarios, logs) y cualquier
archivo >2MB que no sea .png. Al espejar el arbol completo (no solo las
carpetas con "input"), las carpetas compartidas de graficos/scripts que
viven al lado de varias corridas (como en capacitor/local/WALL-DIRICHLET/
graficos/) se copian solas, en su lugar natural.

Pensado para que <proyecto>/runs/ sea lo unico que se sube a git de cada
proyecto: lo minimo necesario para volver a correr cada simulacion
(input, .petscrc, condiciones iniciales chicas, scripts) mas los graficos
de resultado ya generados.

No modifica nada bajo <proyecto>/local/ -- solo crea/actualiza
<proyecto>/runs/. Se puede correr las veces que haga falta (p.ej. despues de
cada corrida nueva); los archivos ya copiados se vuelven a copiar si
cambiaron (sobrescribe, no hace merge ni borra lo que sobra en runs/).

Uso (desde cualquier directorio):
  python3 util/archive_runs.py --dry-run   # solo reporta, no copia nada
  python3 util/archive_runs.py             # copia de verdad
  python3 util/archive_runs.py microgel    # solo ese proyecto (atajo)
"""
import os
import re
import sys
import shutil

# Carpetas en la raiz de Ludwig que NO son proyectos de simulacion con una
# subcarpeta local/ de corridas (codigo fuente, tests, utilidades, etc.).
NOT_PROJECTS = {
    "src", "target", "tests", "util", "config", "docs", "claude", "mpi_s",
    "__pycache__", ".git",
}

SIZE_CAP = 2 * 1024 * 1024  # 2 MB

# Nombre de archivo = volcado de paso de Ludwig: prefijo + corrida de >=4
# digitos (el numero de paso) + ".RANK-NRANKS" (p.ej. dist-000003500.001-001,
# config.cds00001000.001-001, config.cds-000012000.001-001)
STEP_OUTPUT_RE = re.compile(r"^.+\d{4,}\.\d+-\d+$")
METADATA_RE = re.compile(r".*-metadata\.\d+-\d+$")

EXCLUDE_EXACT_NAMES = {"Ludwig.exe", "__pycache__"}
EXCLUDE_SUFFIXES = (".exe", ".pyc", ".o", ".so", ".log")


def should_exclude(name, full_path):
    if name in EXCLUDE_EXACT_NAMES:
        return True
    if name.endswith(EXCLUDE_SUFFIXES):
        return True
    if STEP_OUTPUT_RE.match(name) or METADATA_RE.match(name):
        return True
    if os.path.islink(full_path) and not os.path.exists(full_path):
        return True
    if name.lower().endswith(".png"):
        return False
    try:
        if os.path.getsize(full_path) > SIZE_CAP:
            return True
    except OSError:
        return True
    return False


def find_projects(ludwig_root):
    projects = []
    for name in sorted(os.listdir(ludwig_root)):
        if name in NOT_PROJECTS or name.startswith("."):
            continue
        local_dir = os.path.join(ludwig_root, name, "local")
        if os.path.isdir(local_dir):
            projects.append(name)
    return projects


def archive_project(ludwig_root, project, dry_run):
    local_root = os.path.join(ludwig_root, project, "local")
    runs_root = os.path.join(ludwig_root, project, "runs")

    n_files_copied = 0
    n_files_skipped = 0
    bytes_copied = 0

    for dirpath, dirnames, filenames in os.walk(local_root):
        rel_dir = os.path.relpath(dirpath, local_root)
        dest_dir = runs_root if rel_dir == "." else os.path.join(runs_root, rel_dir)
        for fn in filenames:
            full = os.path.join(dirpath, fn)
            if should_exclude(fn, full):
                n_files_skipped += 1
                continue
            dest_file = os.path.join(dest_dir, fn)
            try:
                sz = os.path.getsize(full)
            except OSError:
                sz = 0
            bytes_copied += sz
            n_files_copied += 1
            if not dry_run:
                os.makedirs(dest_dir, exist_ok=True)
                shutil.copy2(full, dest_file)

    return n_files_copied, n_files_skipped, bytes_copied


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    dry_run = "--dry-run" in sys.argv

    ludwig_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    projects = args if args else find_projects(ludwig_root)

    total_copied = total_skipped = total_bytes = 0
    for project in projects:
        if not os.path.isdir(os.path.join(ludwig_root, project, "local")):
            print(f"{project}: no existe {project}/local/, salteado")
            continue
        n_copied, n_skipped, n_bytes = archive_project(ludwig_root, project, dry_run)
        print(f"{project}: {n_copied} archivos "
              f"{'a copiar' if dry_run else 'copiados'} "
              f"({n_bytes/1e6:.1f} MB), {n_skipped} excluidos")
        total_copied += n_copied
        total_skipped += n_skipped
        total_bytes += n_bytes
    print(f"\nTOTAL: {total_copied} archivos ({total_bytes/1e6:.1f} MB), "
          f"{total_skipped} excluidos")


if __name__ == "__main__":
    main()
