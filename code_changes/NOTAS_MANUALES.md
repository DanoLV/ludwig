## Notas manuales (por qué se hizo cada cambio importante)

Este archivo se copia al principio de `CAMBIOS_CODIGO.md`. Agregar una entrada
por cada cambio de código, la más reciente arriba, con este formato:

```
### AAAA-MM-DD — título corto
- **Qué:** qué hace el cambio.
- **Por qué:** qué problema resuelve y cómo se encontró.
- **Archivos:** lista de archivos tocados.
- **Cómo se verificó:** prueba realizada (y qué NO se verificó).
- **Corridas afectadas:** qué resultados anteriores pueden cambiar.
```

Las entradas siguientes se escribieron el 2026-11-13 a partir de las notas de
trabajo (diario y memoria de Claude). Los cambios anteriores a 2026-09 solo
figuran en el listado automático de abajo.

### 2026-11-16 — Condición de no-flujo en el driver de GPU de Nernst-Planck
- **Qué:** nuevo kernel `np_no_flux_boundary_kernel`, versión de GPU de `np_no_flux_boundary`
  del driver de CPU. Pone a cero todo flujo de enlace (advectivo + difusivo) cuyo nodo o vecino
  no es fluido. Se ejecuta entre los kernels difusivo y de antisimetrización.
- **Por qué:** el driver de GPU (el que se usa desde 2026-03-26) no tenía este paso. El kernel
  difusivo ya omitía los enlaces no fluidos, pero el advectivo `centred` no: por un enlace
  fluido–sólido pasaba W_p (u_i·c_p)/2 · (n_i + n_sólido)/2, que salía del nodo fluido y no
  entraba en ningún nodo actualizado. Encontrado al revisar la §7.2 del documento del modelo
  (pregunta del usuario).
- **Archivos:** `src/nernst_planck.c`.
- **Cómo se verificó:** `capacitor/local/WALL-DIRICHLET/TEST-NOFLUX-ADV`, caja 16×16×18 con
  paredes en z = 1 y 18, electrolito 1:1, LB térmico (kT = 1e-5), advección `centred`, 2000 pasos.
  - Total de iones en el fluido respecto del paso 1: antes (exe 221f3229) entre +8e-6 y −4.1e-5,
    con deriva; después (exe 269ec1b8) |Δ| ≤ 7e-15.
  - Control sin paredes: los 21 archivos qsi antes/después son idénticos byte a byte.
  - No verificado: coloides BBL, `flla`/`flla_mass` con paredes (ya omitían esos enlaces),
    varios procesos MPI, efecto sobre perfiles cerca de la pared.
- **Corridas afectadas:** advección `centred` + nodos sólidos (paredes o coloides BBL) + iones +
  velocidad normal a la superficie no nula en la primera capa de fluido (LB térmico, flujo hacia
  la pared). Con fluido en reposo o flujo paralelo a la pared el efecto es despreciable (p. ej.
  C-electrolito: |u| ≤ 2e-9, iones conservados a 1e-10).

### 2026-11-13 — Interpretación física de la rampa de pared con e0 tangencial (solo comentario)
- **Qué:** se corrigió el comentario del bloque `CHANGE 20261028` en `petsc_wall_value`
  (`src/psi_petsc.c`); sin cambios de código. También `claude/PAREDES_POTENCIAL_CONSTANTE.md` §8
  y el README de `capacitor/local/WALL-DIRICHLET/F-campo-tangencial`.
- **Por qué:** decía "un conductor ideal no admite campo tangencial, así que la pared sigue la
  rampa (el canal electroosmótico estándar)". Es al revés: un conductor ideal es equipotencial. La
  condición ψ_pared = ψ_w − e0·(r − r_c) representa un electrodo resistivo con corriente a lo largo
  de la placa (o un arreglo fino de electrodos a potenciales escalonados), con imagen de conductor
  para las cargas cercanas; no modela paredes aislantes (imagen (ε − ε_p)/(ε + ε_p), carga fija).
- **Archivos:** `src/psi_petsc.c` (comentario), `claude/PAREDES_POTENCIAL_CONSTANTE.md`.
- **Cómo se verificó:** revisión del razonamiento; no se recompiló (solo cambia un comentario).
- **Corridas afectadas:** ninguna numéricamente; cambia el uso físico de la prueba F.

### 2026-11-12 — Advección de iones en Nernst-Planck: normalización y esquema FLLA
- **Qué:** (1) el flujo advectivo por enlace ½(u0+u1)·c_p·½(ρ0+ρ1) ahora se pesa con
  W_p = 2·wgradients[p]; (2) nueva opción de entrada
  `electrokinetics_advection centred | flla | flla_mass` (default `centred`): advección
  con los flujos de masa de la red LB + reconstrucción lineal con limitador superbee
  (DeVincentis & Thomas 2023, ecs. 11–21); `flla_mass` usa el flujo de masa sin
  normalizar (Dubois & Lallemand 2008, ec. 17).
- **Por qué:** sin peso, la suma sobre todos los enlaces daba (M/2)∇·(ρu): iones
  advectados 9× (D3Q27) o 5× (D3Q19) demasiado rápido; con LB térmico la densidad de
  contraiones se volvía negativa (parámetros de Ingo, paso 72). La advección con la
  velocidad nodal además no es consistente con el transporte de masa del LB
  (∇·u_LB ≠ ∂tρ_LB); `flla_mass` hace que los iones sigan la densidad del fluido.
- **Archivos:** `src/nernst_planck.c`, `src/nernst_planck.h`, `src/ludwig.c`.
- **Cómo se verificó:** D3Q7 idéntico bit a bit (peso 1); con LB térmico sin densidades
  negativas; `flla_mass` correlación 1.000 con la densidad del fluido; test del escalón
  (u = 0.01): sin sobrepaso, difusión numérica ≈ 9e-5. NO verificado: MPI (FLLA solo
  con un proceso; `flla_mass` se rechaza con más de uno) ni el driver de CPU (FLLA solo GPU).
- **Corridas afectadas:** toda corrida con hidrodinámica activa y stencil EK de 19 o 27
  puntos advectó iones 5–9× de más (importa solo con flujo: electroósmosis,
  sedimentación, LB térmico). El equilibrio con fluido en reposo no cambia.

### 2026-10-28 — Fuerza de enlaces duplicada en celdas de borde (subgrid)
- **Qué:** se ponen a cero fex/tex de las copias de halo antes de la segunda suma en
  `subgrid_update_forces_electrokinetics`.
- **Por qué:** enlaces/pares de partículas subgrid en celdas de borde recibían 2× la
  fuerza. Encontrado con una cadena que cruza la frontera periódica
  (`capacitor/local/WALL-DIRICHLET/F-campo-tangencial/T4`).
- **Archivos:** `src/subgrid.c`.
- **Cómo se verificó:** suma de fuerzas de enlace 1e-21; borde vs referencia 1e-13.
- **Corridas afectadas:** cadenas/redes con monómeros en celdas de borde anteriores al arreglo.

### 2026-10-14 — σ por defecto de la corrección pm_sr = ancho del kernel
- **Qué:** si `pm_sr_sigma` no está en el input, se usa el segundo momento del kernel
  activo (`pm_sr_kernel_sigma`), en lugar del 0.7 fijo; el log informa de dónde salió σ.
- **Archivos:** `src/subgrid.c`, `src/ludwig.c`.

### 2026-09-26 — Paredes a potencial constante (Dirichlet) con PETSc + fuerza eléctrica duplicada en celdas de borde
- **Qué:** (1) paredes a potencial fijo en el solver PETSc (receta de Asta et al. 2019);
  documentación completa en `claude/PAREDES_POTENCIAL_CONSTANTE.md`. (2) en
  `subgrid_update_forces_electrokinetics` solo la copia propietaria (celdas 1..ncell)
  acumula la fuerza eléctrica.
- **Por qué (2):** las copias de halo sumaban la fuerza completa: ×2 en cara, ×4 en
  arista, ×8 en esquina. Encontrado en la prueba C de las paredes.
- **Archivos:** `src/ludwig.c`, `src/psi_petsc.c`, `src/psi_petsc.h`, `src/psi_gradients.c`, `src/psi.h`, `src/subgrid.c`.
- **Cómo se verificó:** pruebas A–E de `capacitor/local/WALL-DIRICHLET` (ver su README);
  después del arreglo (2) F·β/Esub = 1.000000000000 en todas las posiciones.
- **Corridas afectadas (2):** partículas cargadas en celdas de borde: `microgel/local/CROSS-borde`
  (×2 uniforme), `old/`, y `SLIDES-LD8-hann04-corr_no-dx2` (ejecutable viejo, partícula 2 ×2).

### 2025-10-29 a 2026-05 — Suma compensada (Kahan / Klein) en sumas de punto flotante
- **Qué:** sumas compensadas en lugar de sumas directas de doubles, usando los tipos
  `kahan_t` y `klein_t` de `src/util_sum.h`/`util_sum.c` (esa infraestructura es del
  Ludwig original, K. Stratford 2021; lo propio es su uso en los lugares siguientes):
  - **Transformación de modos del LB** (la más importante): `f2mode` y `mode2f` en
    `src/collision.c` (sumas sobre velocidades / modos con Kahan). Versión D3Q19
    (`d3q19_mode2f_chunk_kahan`): commit `b5bf67b0`, 2025-10-30 ("f2mode and mode2f kahan
    inplementation"). Versión D3Q27 (`d3q27_*_kahan`, bloques "D3Q27 Kahan summation for
    f2mode/mode2f") y reescritura de las marcas actuales: commit `354f9e7b`, 2026-01-14.
  - **Momentos y equilibrio del LB** (commit `354f9e7b`, 2026-01-14): momento 0 (densidad)
    y 1 en `src/lb_data.h`, inicialización del equilibrio en `src/lb_data.c`,
    normalizadores de `src/lb_d3q27.c` (como ya hacía D3Q19).
  - **Fuerzas subgrid** con `klein_t`: `src/psi_force.c` y `src/subgrid.c`
    (2025-10-29 / 2025-11-19); fuga de memoria de los `klein_t` de `entry->force[]`
    corregida el 2026-09-12 (`src/subgrid.c`).
  - **Nernst-Planck**: `src/nernst_planck.c` (2026-02-10; actualización de cargas con
    Kahan, fecha de git 2026-05-14).
  - **Otros**: fuerza total en `src/stats_velocity.c` (Klein), corrección de campo en
    `src/psi_fft_pn.c` (Klein), sumas de Ewald en `src/ewald_charge.c`.
- **Por qué (según el usuario, 2026-11-13):** en un sistema con una partícula cargada,
  una cierta concentración de sal y SIN campo externo, el momento lineal total del fluido
  debía conservarse, pero el error de redondeo de las sumas en la transformación de modos
  (f2mode / mode2f) se acumulaba paso a paso y el fluido se aceleraba de forma
  descontrolada. Con la suma compensada el problema desapareció.
- **Evidencia guardada** (las corridas originales ya no existen):
  - `microgel/runs/old/sim/particula cargada/particula.png` y `fluido.png`: partícula cargada
    con sal, sin campo: Vx/Vy/Vz crecen sin parar hasta ~3.5e-13 en 10^6 pasos, en forma de
    dientes de sierra (en vez de quedar en 0 por conservación del momento).
  - `microgel/runs/old/sim/Liquido solo/`: fluido solo (sin partícula, sin fuerzas, sin
    fluctuaciones; caja 32³, ρ = 0.8, η = 0.05, M10, 20 pasos). Los inputs son idénticos; solo
    cambia el código. `Sin_correcion/`: la velocidad media del fluido salta a
    ~1.7e-17 (x) y ~0.4e-17 (y) en los primeros 3 pasos y queda ahí (debería ser 0).
    `Correcion_f2m_m2f/` (Kahan en f2mode/mode2f; ρ = 0.5…1.0): velocidad media
    exactamente 0 en todos los casos.
  - `microgel/Graficos/` (en git, commit `e93f2fb4` del 2025-11-10, posterior al arreglo):
    `Fluid/fluido-none-rho_*-velocidades.png` (sin corrección: la velocidad media del fluido
    solo queda en ~±1e-17 en vez de 0) vs `Fluid/*-momento_kahan*.png` (con corrección: 0
    exacto). `Sin campo/single-…-q_1.0-e_0.0…-N_{250k,270k,500k}-*.png`: partícula cargada
    sin campo; tras un transitorio la velocidad de la partícula queda en ~1e-16 y la media
    del fluido deambula en ~1e-17 sin crecer sistemáticamente (vs ~3.5e-13 creciendo de
    `particula cargada/`). Que estas corridas sean con la corrección se INFIERE de la fecha
    del commit; el nombre del archivo no lo dice.
- **Historia:** el problema se planteó en el issue de Ludwig
  [ludwig-cf/ludwig#355](https://github.com/ludwig-cf/ludwig/issues/355) ("Charged
  polymer", DanoLV). Comentarios del 2025-08-08 (la velocidad de la partícula crece aun sin
  campo) y 2025-10-10 (el fluido mismo se acelera; con SOR llega a un estado estacionario,
  con PETSc no). La solución (commit `b5bf67b0`, 2025-10-30) no se publicó en el issue,
  que se cerró el 2026-01-29 sin ella.
- **Por completar:** medida cuantitativa de la mejora con partícula cargada (momento total
  vs paso con y sin Kahan).
- **Archivos:** `src/collision.c`, `src/lb_data.h`, `src/lb_data.c`, `src/lb_d3q27.c`,
  `src/psi_force.c`, `src/subgrid.c`, `src/nernst_planck.c`, `src/stats_velocity.c`,
  `src/psi_fft_pn.c`, `src/ewald_charge.c`.
- **Nota sobre las fechas:** muchas marcas de estos cambios no tienen fecha; en el
  listado de abajo aparecen con la fecha del commit que las introdujo *(fecha de git)*.

### 2025-12-03 — Normalización del laplaciano de los stencils D3Q19 y D3Q27
- **Qué:** se corrigió el peso del laplaciano en los stencils de 19 y 27 puntos.
- **Archivos:** `src/stencil_d3q19.c`, `src/stencil_d3q27.c`.
- **Por qué / cómo se verificó:** Poisson con los stencils de 7, 19 y 27 puntos daba
  soluciones distintas; el factor de normalización correcto es 6 en todos los casos (no 36 y
  216): `s->wlaplacian[p] = -6.0*wv[p];`. Verificado aplicando el stencil a
  f = x² + y² + z² (laplaciano analítico 6; el stencil da −6), con el test unitario
  `test_stencil_d3q27_analytical_laplacian_quadratic`.
- **Publicado en el Ludwig original:** issue
  [ludwig-cf/ludwig#389](https://github.com/ludwig-cf/ludwig/issues/389) (2026-01-12). La
  corrección entró en la rama `develop` del Ludwig original por el PR #411 y el issue se
  cerró el 2026-05-07 ("aparecerá en la próxima versión").
- **Nota relacionada (no es un cambio de código):** el flujo difusivo de Nernst-Planck
  suma sobre todos los enlaces sin dividir por A0, así que `electrokinetics_d0/d1` no es
  el D físico: D_físico = d·A0 (A0 = 3.828 en D3Q19, 6.138 en D3Q27).
