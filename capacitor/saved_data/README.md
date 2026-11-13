# capacitor/saved_data

Datos de capacitor que se guardan en git a futuro (creada 2026-11-13).

- `capacitor/local/` (corridas completas) y `capacitor/runs/` (espejo liviano) **no** van a git.
- El registro de todas las corridas es `capacitor/Registro_simulaciones.xlsx` (sí va a git;
  se regenera con `python3 util/registro_simulaciones.py`).
- Acá va solo lo que deba conservarse realmente: resultados finales, figuras de
  referencia, condiciones iniciales necesarias para rehacer una corrida, etc.
- Dentro de esta carpeta no se aplican las reglas globales del `.gitignore`
  (`*.csv`, `*.vtk`, `*001-001`, …): todo lo que se ponga acá se sube.
  Evitar archivos > 50 MB (GitHub advierte en 50 MB y rechaza desde 100 MB).

## Bases para rehacer corridas (util/guardar_bases.py, 2026-11-13)

Una subcarpeta por carpeta-serie (mismo nombre que en local/), con:
- `input.base`: input real (con comentarios) de la corrida que menos difiere del resto.
- `coloides_<k>.cds`: archivo(s) de coloides base de la carpeta.
- `deltas.json`: por corrida, las claves del input que cambian/se agregan/se quitan y
  las líneas del archivo de coloides que cambian (posición, carga, ...).

Rehacer una corrida (ruta = columna "Ruta" del Registro_simulaciones.xlsx):

    python3 util/rehacer_corrida.py <proyecto> <ruta> <carpeta destino>

Regenerar bases y deltas tras corridas nuevas (verifica todas al final):

    python3 util/guardar_bases.py

No se guardan: archivos auxiliares que nombre el input (.petscrc, tablas, anillo.cds…),
ni el ejecutable (ver code_changes/ y git log para la versión del código).
