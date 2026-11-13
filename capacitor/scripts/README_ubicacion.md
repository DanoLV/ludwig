# capacitor/scripts — ubicación de los scripts (2026-11-13)

Los scripts de capacitor se movieron acá desde `capacitor/local/` para que queden en git
(`local/` no se sube). Para no romper rutas relativas, en el lugar original quedó:

- **carpetas `scripts/`** (`local/scripts`, `local/<serie>/scripts`): un **enlace simbólico**
  a la carpeta de acá. Los scripts que usan `__file__` o `../` funcionan si se llaman por
  la ruta vieja (p. ej. `python3 local/WALL-DIRICHLET/scripts/analizar_C.py ...`); llamados
  desde acá, `../` apunta a otro lugar.
- **scripts sueltos de cada serie** (`local/<serie>/run.sh`, `correr*.sh`, `crear_*.py`, …):
  un **enlace duro** (el mismo archivo con dos nombres), porque varios lanzadores hacen
  `cd "$(dirname "$(readlink -f "$0")")"` y con un enlace simbólico entrarían acá.
  Editar cualquiera de los dos nombres cambia el mismo archivo. Ojo: un editor que guarda
  escribiendo un archivo nuevo y renombrándolo rompe el enlace duro.

Lista completa de lo movido: `code_changes/scripts_movidos_20261113.json`.
Documentación de cada script: `SCRIPTS.md` en la raíz (`python3 util/documentar_scripts.py`).
