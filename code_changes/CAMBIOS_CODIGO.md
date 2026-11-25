# Cambios de código en Ludwig (registro)

Generado por `code_changes/extraer_cambios.py` el 2026-11-16 20:19 a partir de las marcas `/*CHANGE INIT - <fecha> <descripción> */ … /*CHANGE END */` del código (carpetas: src, target, tests, util, mpi_s, config). No editar a mano: volver a correr el script. Las explicaciones largas van en `code_changes/NOTAS_MANUALES.md`, que se copia abajo.

**529 bloques de cambio** en **60 archivos**, 62 fechas.

Fechas: 343 bloques tienen la fecha en la marca; 186 marcas sin fecha se fecharon con git (primer commit que introdujo esa línea; señaladas *(fecha de git)*: es la fecha del commit, el código pudo escribirse antes).

Limitación: solo registra cambios marcados con `CHANGE INIT`. Un cambio sin marca no aparece (para eso, `git log` / `git diff`).


---

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


---

## Índice por fecha

| Fecha | Bloques | Archivos | Descripciones (resumen) |
|---|---|---|---|
| [2026-11-16](#2026-11-16) | 2 | 1 | No-flux boundary in the GPU NP driver. GPU port of np_no_flux_boundary (CPU driver, abo…; No-flux boundary (GPU port of np_no_flux_boundary). Original: (no such step in the GPU … |
| [2026-11-12](#2026-11-12) | 11 | 3 | Lattice (FLLA) advection: state shared by nernst_planck_lattice_flux_compute() (called …; FLLA advection is implemented in the GPU driver only; Advective link-flux normalisation. The unweighted sum over all stencil links gives sum_… (+7 más) |
| [2026-11-11](#2026-11-11) | 5 | 2 | hydrodynamics no: without fluid there is nothing to spread the particle force onto. The…; hydrodynamics no: hydro may be NULL; only the fluid-side copies and the spreading of th…; hydrodynamics no (+1 más) |
| [2026-11-04](#2026-11-04) | 12 | 4 | stretched Peskin6: support /r/ < 3s, range 3s - 1 (same convention as the original s = …; Peskin6 stretched by an integer factor. The unstretched kernel is now d_peskin6_unit();…; Peskin6 stretched by an integer factor (+8 más) |
| [2026-10-28](#2026-10-28) | 7 | 3 | halo copies must not re-add the earlier fex. subgrid_force_from_particles() already ran…; tangential e0 with walls. Original: return solver->block->wall_psi[petsc_wall_side(...)…; tangential e0 with walls. Potential jump picked up going from DMDA global site (i,j,k) … (+4 más) |
| [2026-10-21](#2026-10-21) | 9 | 3 | alternative (Steele, paper) wall model; Steele 10-4-3 repulsive wall of Hoang Ngoc Minh et al., J. Chem. Phys. 158, 104103 (202…; Steele model: the loops use hc as the cut-off (+5 más) |
| [2026-10-18](#2026-10-18) | 6 | 3 | wall planes given explicitly (see wall_ss_cut.h); wall planes given explicitly; map-based walls: the planes are given explicitly and the separation is h = (distance ce… (+2 más) |
| [2026-10-14](#2026-10-14) | 5 | 3 | default pm_sr_sigma from the kernel's 2nd moment; default sigma = width of the active kernel, from its 2nd moment (variance-matched Gauss…; remember whether sigma came from the input. Original: if (rt_double_parameter(ludwig->r… (+1 más) |
| [2026-09-29](#2026-09-29) | 2 | 1 | stale check, unrelated to the OS update. NPAD_INT went from 7 to 3 in colloid.h on 2026…; same stale-constant issue as NPAD_INT above: NBOND_MAX went from 2 to 6 on 20260114. Or… |
| [2026-09-26](#2026-09-26) | 25 | 6 | constant-potential walls. Set by psi_solver_petsc_wall_set(); NULL (calloc) = no walls.…; add the force to the OWNER copy only. The cell loop includes the halo cells (0 and ncel…; constant-potential walls (psi_petsc.c). The particle-node correction reaches r_cut (~6)… (+20 más) |
| [2026-09-24](#2026-09-24) | 19 | 1 | spherical (microgel) mode; spherical mode; spherical (microgel) mode: a free-floating sphere has no walls and no faces, so both ar… (+8 más) |
| [2026-09-23](#2026-09-23) | 3 | 2 | Peskin stretched to support 8; Peskin 4-point function stretched to support 8 |
| [2026-09-22](#2026-09-22) | 23 | 5 | truncated Gaussian kernel parameters; Gaussian kernel; Gaussian kernel range = ceil(n/2) (+8 más) |
| [2026-09-21](#2026-09-21) | 2 | 2 | no DeviceToHost after the solve. Every solver leaves its result on the HOST: the FFT so…; DMDA stencil width = reach of the Poisson stencil, not the lattice halo. The width only… |
| [2026-09-20](#2026-09-20) | 8 | 2 | output frequency for the per-particle csv files particle_force.csv and particle_Esub.cs…; output frequency for the per-particle csv files See pforce_io_freq_ at file scope for t…; (sin descripción) (+2 más) |
| [2026-09-19](#2026-09-19) | 1 | 1 | minimum image convention. The raw coordinate difference is wrong for pairs straddling a… |
| [2026-09-18](#2026-09-18) | 4 | 3 | measured mesh reference for the pair table pm_sr_mesh_pair_esub_x Reference term comput…; measured-reference pair table; measured mesh reference for the pair table. "pm_sr_pair_ref measured" builds the refere… (+1 más) |
| [2026-09-17](#2026-09-17) | 1 | 1 | linear interpolation in the pair table. Previously the nearest bin was used as-is (ir =… |
| [2026-09-16](#2026-09-16) | 3 | 2 | continuum kernel-shaped target Gauss-Legendre nodes/weights on [-1,1] (Newton iteration…; continuum kernel-shaped target setup; selectable pair-correction target shape. target_shape: 0 = Gaussian (analytic erf form,… |
| [2026-09-12](#2026-09-12) | 1 | 1 | fuga de memoria: los tres klein_t de entry->force[] se asignan en subgrid.c (3 sizeof(k… |
| [2026-09-10](#2026-09-10) | 9 | 1 | minimum image for the dual real-space kernels. These kernels computed node<->particle a…; minimum image (see MIN_IMAGE_DUAL above). Correct for a single rank; an MPI build would…; use the halo instead of discarding it. nhalo is sized to the Ewald cutoff (ludwig.c:309… (+4 más) |
| [2026-09-02](#2026-09-02) | 5 | 3 | direct colloid-colloid PM short-range correction pm_sr_apply_force_correction (above) c…; direct colloid-colloid PM short-range correction; direct colloid-colloid PM short-range correction Corrects the discrete PM (kernel-on-me… |
| [2026-08-23](#2026-08-23) | 8 | 3 | per-particle force output for theory comparison; per-particle force output for theory comparison particle_Esub.csv records the field fel…; per-particle force output for theory comparison fex now holds the final averaged F = 1/… |
| [2026-07-29](#2026-07-29) | 17 | 8 | Gaussian_Ewald_Dual - per-species widths and precomputed σ_eff for pair interactions; Gaussian_Ewald_Dual — Dual-sigma kernels (Paso 3) Each kernel below is the analogue of …; Gaussian_Ewald_Dual (Paso 3b) — Fourier kernels with receiver form factor The structure… (+8 más) |
| [2026-07-15](#2026-07-15) | 1 | 1 | diagnostic: L2 norm of the Nyquist-x plane of a field. Returns sqrt(sum over k with mx=… |
| [2026-07-14](#2026-07-14) | 3 | 1 | for the interlacing random-offset hash; interlacing offset from input + random per-step mode. Input keys (all optional, default…; interlacing configuration from input |
| [2026-07-10](#2026-07-10) | 8 | 6 | trilinear: exact relabelling at integer offsets; trilinear (interpolating) kernel; spectral (FFT phase) shift of a scalar lattice field psi_fft_shift_field Translates a p… (+5 más) |
| [2026-07-06](#2026-07-06) | 3 | 1 | Hann kernel in PM self-force utility |
| [2026-06-30](#2026-06-30) | 22 | 7 | Hann spread/gather kernel order (configurable); Hann kernel; Hann kernel range = ceil(n/2) (+12 más) |
| [2026-06-29](#2026-06-29) | 1 | 1 | select subgrid kernel from input file |
| [2026-06-25](#2026-06-25) | 7 | 3 | PM short-range correction table; pm_sr_correct_phi / pm_sr_correct_efield; PM short-range correction tables (+2 más) |
| [2026-06-11](#2026-06-11) | 5 | 1 | Macro for dphi_dr with sigma->0 limit guard Used by all _dual real-space kernels. When …; sigma->0 guard for _pf potential; Handle sigma_ff=0 limit: erf(r/0)/r -> 1/r (+2 más) |
| [2026-06-10](#2026-06-10) | 1 | 1 | dual solver writes to ewald_charge_forces.csv, not particle_Esub.csv |
| [2026-05-18](#2026-05-18) | 2 | 2 | fix sprintf: remove dangling PB format specifiers; PM self-force utility using Peskin4 spread/gather |
| [2026-05-14](#2026-05-14) | 73 | 16 | psi_rho_zero; subgrid_charge_from_grid; Use Poisson-Boltzmann weight instead of Peskin (+63 más) |
| [2026-05-13](#2026-05-13) | 6 | 4 | setter for colloid info pointer; (sin descripción) |
| [2026-05-12](#2026-05-12) | 10 | 6 | new file; add FFT and EWALD solver ids; Ewald solver wrapper (+4 más) |
| [2026-04-27](#2026-04-27) | 2 | 2 | subgrid_scatter_fluid_offset: scatter fluid rho into auxiliary buffer; fluid scatter with offset into auxiliary buffer |
| [2026-04-26](#2026-04-26) | 2 | 1 | psi_force_gradmu_offset for interlacing; psi_force_gradmu_e_ewald_offset |
| [2026-04-25](#2026-04-25) | 6 | 2 | forward declare for interlacing offset variant; subgrid_charge_from_grid_offset for interlacing; Interlacing mode selector: 0 = No interlacing, use only Grid A (offset=0.0) 1 = Grid B … (+3 más) |
| [2026-04-24](#2026-04-24) | 3 | 1 | Kaiser-Bessel globals; Peskin 6-point kernel (Bao et al. 2016); Kaiser-Bessel order-4 kernel |
| [2026-04-22](#2026-04-22) | 23 | 7 | kernel parameter for subgrid_charge_from_particles; kernel parameter for subgrid_charge_from_grid; kernel parameter for subgrid_update_forces_electrokinetics (+15 más) |
| [2026-04-18](#2026-04-18) | 2 | 1 | Scale Peskin kernel with drange_ to support variable stencil width. drange_=1 (original…; Scale Peskin kernel with drange_ |
| [2026-03-26](#2026-03-26) | 9 | 4 | GPU implementation of nernst_planck_driver_d3qx nernst_planck_driver_d3qx_gpu GPU-paral…; Fix SEGV: target is a GPU struct, cannot deref ->data from host. Use tdpMemcpy(DeviceTo…; Use GPU vec/mat types when -vec_type cuda is set via .petscrc. DMSetVecType/DMSetMatTyp… (+6 más) |
| [2026-03-10](#2026-03-10) | 2 | 1 | Self-field correction |
| [2026-03-09](#2026-03-09) | 2 | 2 | GPU kernels for self-field correction table; Self-field correction for lattice-discretisation artefact |
| [2026-03-02](#2026-03-02) | 1 | 1 | Peskin-spread FFT Ewald version |
| [2026-02-12](#2026-02-12) | 4 | 2 | New ewald_charge_sum_FFT_full_gpu using cuFFT for Fourier-space |
| [2026-02-10](#2026-02-10) | 2 | 1 | Kahan compensated sum from util_sum; Use kahan_t from util_sum.h |
| [2026-02-06](#2026-02-06) | 1 | 1 | Antisymmetric flux assignment for exact conservation |
| [2026-01-21](#2026-01-21) | 2 | 2 | Refined grid Poisson solver |
| [2026-01-19](#2026-01-19) | 6 | 3 | New FFT-based Poisson solver using cuFFT; New FFT-based Poisson solver using cuFFT Edinburgh Soft Matter and Statistical Physics …; FFT Poisson solver |
| [2026-01-17](#2026-01-17) | 4 | 2 | cambio de soporte de peskin // static const double drange_ = 1.0; Max. range of interpo…; Short-range corrections infrastructure (Level 0) See docs/SHORT_RANGE_CORRECTIONS_ANALY…; Short-range corrections Level 1: Self-field See docs/SHORT_RANGE_CORRECTIONS_ANALYSIS.m… (+1 más) |
| [2026-01-14](#2026-01-14) | 39 | 17 | Subgrid charge // int psi_electroneutral(psi_t obj, map_t map); int psi_electroneutral(…; Extra bonds for microgel polymer networks (6 bonds, intpad reduced from 7 to 3) // #def…; Momentum conservation correction (+27 más) |
| [2025-12-03](#2025-12-03) | 26 | 9 | Electric field output; Fix Laplacian normalization for consistency; Electric field output --- IGNORE --- test_psi_initialise(pe); test_psi_create(pe); test… |
| [2025-11-23](#2025-11-23) | 6 | 2 | Poisson-Boltzmann kappa (inverse Debye length); Poisson-Boltzmann force calculation; Write electric field data to file (+1 más) |
| [2025-11-19](#2025-11-19) | 33 | 14 | Subgrid charge // for (int ia = 0; ia < model->nvel; ia++) { // wip += model->wv[ia] mo…; CUDA C++ compatibility fix; Write particle data after colloid_sums_halo, only for local (non-halo) particles (+10 más) |
| [2025-11-10](#2025-11-10) | 10 | 4 | Subgrid charge static int colloid_sums_m8(colloid_sum_t sum, int, int, int, int); //CHA…; Subgrid charge // static const int msize_[COLLOID_SUM_MAX] = { 0, 10, 35, 7, 4, 6, 7, 1…; Subgrid charge if (sum->mtype == COLLOID_SUM_ELECTRIC_FIELD) mloader_forw = colloid_sum… (+4 más) |
| [2025-10-29](#2025-10-29) | 3 | 1 | Subgrid charge klein_t flocal_k[3]; flocal_k[X] = klein_zero(); flocal_k[Y] = klein_zer…; Subgrid charge psi_electric_field(psi, index, e); // CHANGE psi->efield->data[addr_rank…; Subgrid charge psi_force_gradmu_e_subgrid |
| [2025-07-23](#2025-07-23) | 3 | 1 | Subgrid charge; (sin descripción) |
| [2025-07-01](#2025-07-01) | 7 | 5 | Subgrid charge |
| [2025-01-12](#2025-01-12) | 1 | 1 | Fix NULL pointer handling for test compatibility |

## Detalle

### 2026-11-16

Archivos: `src/nernst_planck.c`

- **No-flux boundary in the GPU NP driver. GPU port of np_no_flux_boundary (CPU driver, above): zero every link flux (advective + diffusive) whose site OR neighbour is not fluid. Before this, the GPU driver had no such step: the diffusive kernel already skipped non-fluid links, but the centred advective kernel did not, so a** — `src/nernst_planck.c:1418-1458`
- **No-flux boundary (GPU port of np_no_flux_boundary). Original: (no such step in the GPU driver)** — `src/nernst_planck.c:1928-1934`

### 2026-11-12

Archivos: `src/ludwig.c`, `src/nernst_planck.c`, `src/nernst_planck.h`

- **NP advection scheme: electrokinetics_advection centred | flla | flla_mass (default centred) flla = lattice advection with superbee limiter (DeVincentis & Thomas 2023), see nernst_planck.c.** — `src/ludwig.c:591-600`
- **Store the LB link fluxes for FLLA ion advection (read-only on the distributions; no-op unless electrokinetics_advection flla). Post-collision, post-bounce-back, pre-propagation.** — `src/ludwig.c:2268-2274`
- **Lattice (FLLA) advection: state shared by nernst_planck_lattice_flux_compute() (called in ludwig.c after bounce-back, just before propagation) and the GPU NP driver (next time step).** — `src/nernst_planck.c:74-86`
- **FLLA advection is implemented in the GPU driver only** — `src/nernst_planck.c:389-393`
- **Advective link-flux normalisation. The unweighted sum over all stencil links gives sum_p flux_p = (M/2) div(rho u), M = sum_p c_x^2 (D3Q7: 1, D3Q19: 5, D3Q27: 9), so ions were advected M/2 times faster than the fluid on D3Q19/D3Q27. Weight each link with W_p = 2 wgradients[p], so that** — `src/nernst_planck.c:1123-1132`
- **Advective link-flux normalisation (see np_advective_fluxes)** — `src/nernst_planck.c:1257-1259`, `src/nernst_planck.c:1288-1291`
- **Lattice (FLLA) advection of the ions. Reference: B. DeVincentis, J. Thomas, "Application of flux limiters to passive scalar advection for the lattice Boltzmann method", Comput. Math. Appl. 144, 1-11 (2023), doi:10.1016/j.camwa.2023.05.031.** — `src/nernst_planck.c:1505-1795`
- **Advective link-flux weights W_p = 2 wgradients[p] (see np_advective_fluxes).** — `src/nernst_planck.c:1843-1854`
- **Lattice (FLLA) advection, if selected. Before the first lattice flux exists (first call) no advection is applied.** — `src/nernst_planck.c:1896-1913`
- **Lattice (FLLA) advection of the ions. DeVincentis & Thomas, Comput. Math. Appl. 144, 1 (2023): advective fluxes from the LB distribution (their eq. 21) + linear reconstruction with the superbee flux limiter (their eqs. 13-19). Selected with the input key electrokinetics_advection centred | flla | flla_mass (default centred)** — `src/nernst_planck.h:32-41`

### 2026-11-11

Archivos: `src/ludwig.c`, `src/subgrid.c`

- **hydrodynamics no: total force needs the fluid** — `src/ludwig.c:2303-2310`
- **hydrodynamics no: without fluid there is nothing to spread the particle force onto. The total external force on each particle (lubrication + halo sum above) is kept.** — `src/subgrid.c:249-253`
- **hydrodynamics no: hydro may be NULL; only the fluid-side copies and the spreading of the force to the fluid are skipped (the force on the particle is computed as before).** — `src/subgrid.c:1443-1449`
- **hydrodynamics no** — `src/subgrid.c:1677-1680`, `src/subgrid.c:1689-1692`

### 2026-11-04

Archivos: `src/ludwig.c`, `src/psi_fft.c`, `src/subgrid.c`, `src/subgrid.h`

- **Hann exponent p (w ~ cos^p) from input** — `src/ludwig.c:1095-1113`
- **Peskin6 integer stretch s (support 6s) from input** — `src/ludwig.c:1114-1124`
- **Peskin6 stretch must be known before sizing the halo** — `src/ludwig.c:3451-3459`
- **stretched Peskin6: support |m| < 3s. Original: for (int m = -3; m <= 3; m++) {** — `src/psi_fft.c:1158-1162`
- **stretched Peskin6: support |r| < 3s, range 3s - 1 (same convention as the original s = 1 value 2, see the floor/ceil in subgrid_get_lattice_index_range). Original: else if (kernel == SUBGRID_KERNEL_PESKIN6) return 2;** — `src/subgrid.c:3166-3171`
- **Peskin6 stretched by an integer factor. The unstretched kernel is now d_peskin6_unit(); d_peskin6() below maps w(r) = phi6(r/s)/s, s = subgrid_peskin6_stretch_ (default 1: unchanged). Original: double d_peskin6(double s) { ... } (body unchanged below).** — `src/subgrid.c:4028-4033`
- **Peskin6 stretched by an integer factor** — `src/subgrid.c:4045-4081`, `src/subgrid.h:59-62`
- **Hann window with an even exponent p. w(r) = cos^p(pi r/n) / (n C_p), |r| < n/2, C_p = (p-1)!!/p!! (C_p n = int cos^p(pi r/n) dr over the support, Wallis integral). p = 2 is** — `src/subgrid.c:4180-4225`
- **Peskin6 stretch in the cache key (s = 1 keeps the old key). Original: return (kernel == SUBGRID_KERNEL_GAUSS) ? subgrid_gauss_support_ : subgrid_hann_order_;** — `src/subgrid.c:6460-6466`
- **Hann exponent in the cache key (p = 2 keeps the old key)** — `src/subgrid.c:6467-6470`
- **Hann window exponent p (cos^p)** — `src/subgrid.h:66-69`

### 2026-10-28

Archivos: `src/psi_gradients.c`, `src/psi_petsc.c`, `src/subgrid.c`

- **tangential e0 with walls. The doubled link reads the wall at the link MIDPOINT (half a spacing away), but the wall node stores its value at r + c. With the ramp psi_wall = psi_w - e0.r (psi_petsc.c) they differ by e0.c/2 on a diagonal link; without it the tangential field in the first fluid** — `src/psi_gradients.c:86-97`
- **tangential e0 with walls. Original: return solver->block->wall_psi[petsc_wall_side(...)]; A tangential external field enters as Ludwig's usual potential jump psi(x + L) = psi(x) - e0 L (psi_halo_psijump). The wall follows the applied ramp** — `src/psi_petsc.c:260-291`
- **tangential e0 with walls. Potential jump picked up going from DMDA global site (i,j,k) along displacement c: the stored (periodic) neighbour value plus this is the physical value, with the convention of psi_halo_psijump.** — `src/psi_petsc.c:294-311`
- **tangential e0 with walls: the jump term, link by link and for fluid neighbours only (a wall neighbour already carries the ramp in petsc_wall_value). Without walls the per-face sum below gives the same, since the wlaplacian weights with c_a = -1 add up to -1.** — `src/psi_petsc.c:709-718`
- **tangential e0 with walls: skip when walls are set (jump already in the loop above). Original: if (e0[X] || e0[Y] || e0[Z])** — `src/psi_petsc.c:750-753`
- **tangential e0 with walls: only the component normal to the walls is refused (that is a voltage between the plates, set it with wall_psi). Original: any non-zero e0 returned -2.** — `src/psi_petsc.c:1167-1171`
- **halo copies must not re-add the earlier fex. subgrid_force_from_particles() already ran this same (symmetric) halo sum on fex = bonds + pairs + angles + lubrication, which leaves the owner AND each periodic copy holding the total T. The owner then gets the electrostatic force E above, and the sum below added the copy's T** — `src/subgrid.c:1635-1662`

### 2026-10-21

Archivos: `src/colloids_rt.c`, `src/wall_ss_cut.c`, `src/wall_ss_cut.h`

- **wall_ss_cut_model soft_sphere (default) | steele steele: repulsive 10-4-3 wall of Hoang Ngoc Minh et al., J. Chem. Phys. 158, 104103 (2023), Eqs. (7)-(8); epsilon = eps_w, sigma = sigma_w, and wall_ss_cut_nu / wall_ss_cut_hc are not needed (cut-off 0.987 sigma_w). Original lines (all four keys required) kept below for soft_sphere.** — `src/colloids_rt.c:1446-1463`
- **nu and hc only checked for the soft sphere** — `src/colloids_rt.c:1470-1474`
- **alternative (Steele, paper) wall model** — `src/wall_ss_cut.c:58-60`
- **Steele 10-4-3 repulsive wall of Hoang Ngoc Minh et al., J. Chem. Phys. 158, 104103 (2023): Vads(z) = 2 pi eps_w [ (2/5)(s/z)^10 - (s/z)^4 - sqrt(2)/(3 (z/s + a)^3) ], a = 0.61/sqrt(2), s = sigma_w (Lennard-Jones integrated over a fcc (100) face) Vrep(z) = Vads(z) - Vads(z ) for z <= z , 0 otherwise,** — `src/wall_ss_cut.c:63-80`
- **Steele model: the loops use hc as the cut-off** — `src/wall_ss_cut.c:119-122`
- **(sin descripción)** — `src/wall_ss_cut.c:159-168`, `src/wall_ss_cut.c:175-177`
- **Steele repulsive wall (paper): f = -dV/dz > 0** — `src/wall_ss_cut.c:374-390`
- **alternative wall model: purely repulsive 10-4-3 (Steele) ion-wall potential of Hoang Ngoc Minh et al., J. Chem. Phys. 158, 104103 (2023), Eqs. (7)-(8). model 0 = original soft sphere (epsilon, sigma, nu, hc); model 1 = Steele, using epsilon = eps_w and sigma = sigma_w (nu and hc are ignored; the cut-off is z sigma_w).** — `src/wall_ss_cut.h:41-47`

### 2026-10-18

Archivos: `src/colloids_rt.c`, `src/wall_ss_cut.c`, `src/wall_ss_cut.h`

- **map-based walls (porous_media_init wall_x|y|z): wall_ss_cut_planes yes|no use the planes below (default no) wall_ss_cut_axis x|y|z axis normal to the planes (default z) wall_ss_cut_plane_lo <double> coordinate of the lower wall plane wall_ss_cut_plane_hi <double> coordinate of the upper wall plane** — `src/colloids_rt.c:1476-1500`
- **wall planes given explicitly (see wall_ss_cut.h)** — `src/wall_ss_cut.c:51-57`
- **wall planes given explicitly** — `src/wall_ss_cut.c:112-118`, `src/wall_ss_cut.c:181-188`
- **map-based walls: the planes are given explicitly and the separation is h = (distance centre-plane) - radius. h is kept above a tiny floor so that wall_ss_cut_single (assert h > 0) never aborts; the repulsion is then huge but finite. The original loop over wall->param->isboundary below is untouched for use_planes = 0.** — `src/wall_ss_cut.c:239-287`
- **wall planes given explicitly (map-based walls). The original potential only acts on the boundary walls of the wall_t object (isboundary); constant-potential walls built with porous_media_init wall_x|y|z are map sites and wall_t knows nothing of them. Original struct ended at 'hc'.** — `src/wall_ss_cut.h:30-40`

### 2026-10-14

Archivos: `src/ludwig.c`, `src/subgrid.c`, `src/subgrid.h`

- **default sigma = width of the active kernel, from its 2nd moment (variance-matched Gaussian), which is the value the production runs set by hand. pm_sr_sigma in the input still overrides it. Original: double pm_sr_sigma = 0.7;** — `src/ludwig.c:1396-1403`
- **remember whether sigma came from the input. Original: if (rt_double_parameter(ludwig->rt, "pm_sr_sigma", &v) == 1) pm_sr_sigma = v;** — `src/ludwig.c:1412-1419`
- **report where sigma came from; warn when r_cut is too short for the Gaussian reference to reach point Coulomb (needs r_cut >~ 2.9 sigma, see the 20260724 note above). Original: pe_info(ludwig->pe, " sigma = %.4g\n", pm_sr_sigma);** — `src/ludwig.c:1428-1438`
- **default pm_sr_sigma from the kernel's 2nd moment** — `src/subgrid.c:5390-5425`, `src/subgrid.h:189-191`

### 2026-09-29

Archivos: `tests/unit/test_colloid.c`

- **stale check, unrelated to the OS update. NPAD_INT went from 7 to 3 in colloid.h on 20260114 ("extra bonds for microgel polymer networks"), which also updated sref above (intpad {30,31,32} has 3 entries), but left this hardcoded assert unchanged -- it aborted the whole unit test binary before any suite could run.** — `tests/unit/test_colloid.c:82-89`
- **same stale-constant issue as NPAD_INT above: NBOND_MAX went from 2 to 6 on 20260114. Original: assert(NBOND_MAX == 2);** — `tests/unit/test_colloid.c:214-217`

### 2026-09-26

Archivos: `src/ludwig.c`, `src/psi.h`, `src/psi_gradients.c`, `src/psi_petsc.c`, `src/psi_petsc.h`, `src/subgrid.c`

- **constant-potential walls: parsed in ludwig_rt(), attached just before the time loop (see ludwig_walls_attach).** — `src/ludwig.c:201-207`
- **constant-potential (Dirichlet) walls for PETSc. electrokinetics_wall_potential psi_lo_psi_hi (units of psi) electrokinetics_wall_axis x|y|z (lower/upper half) Absent keys: the original fully periodic solver. Only read and checked here. The walls are attached to the solver just** — `src/ludwig.c:383-411`
- **optional with constant-potential walls. electrokinetics_electroneutral yes|no (default yes) With Dirichlet walls the Poisson problem is well posed for any net charge (the walls carry the compensating charge), so the counterions that make the fluid neutral can be left out. Without walls a net charge** — `src/ludwig.c:563-586`
- **constant-potential walls** — `src/ludwig.c:946-982`, `src/psi_petsc.c:61-63`, `src/psi_petsc.c:1140-1171`
- **with constant-potential walls the measured reference is still calibrated with the PERIODIC solver: the walls are attached after this block (ludwig_walls_attach). So the table is the bulk one and the pm_sr_meshref cache stays valid for periodic runs. It does NOT describe the particle-wall (image) short range: that is** — `src/ludwig.c:1488-1494`
- **constant-potential walls (after the pm_sr tables)** — `src/ludwig.c:1513-1515`
- **constant-potential walls. Set by psi_solver_petsc_wall_set(); NULL (calloc) = no walls. Lets psi_electric_field() apply the same wall rule as the Poisson solver.** — `src/psi.h:77-81`
- **constant-potential walls (psi_petsc.c). With walls (psi->wall_map set), a fluid node uses the DIFFERENCE form E_a = - sum_p wgradients[p] c_pa [psi(r+c_p) - psi(r)] (1 + chi(r+c_p)) with chi = 1 at MAP_BOUNDARY sites: the same factor 2 on wall links as the wall Laplacian (Asta et al. 2019, eq. 15), because the physical wall** — `src/psi_gradients.c:32-49`
- **constant-potential walls (see comment above)** — `src/psi_gradients.c:69-97`
- **constant-potential walls: stubs without PETSc** — `src/psi_petsc.c:191-201`
- **constant-potential walls (all zero = periodic)** — `src/psi_petsc.c:216-223`
- **constant-potential walls: helpers. Arguments are LOCAL lattice coordinates; halo values are allowed (the map halo is refreshed in psi_solver_petsc_wall_set()).** — `src/psi_petsc.c:226-291`
- **constant-potential walls. Rebuilding (wall_set runs after initialise) must start from zero values: the nonzero pattern is frozen (MAT_NEW_NONZERO_LOCATIONS false) and a row that now writes fewer entries would otherwise keep the old ones.** — `src/psi_petsc.c:468-475`
- **constant-potential walls (Asta et al. 2019, eq. 15). Wall site: decoupled row D x = D psi_wall (rhs_set). Fluid site: each link to a wall site is doubled, and its column is ELIMINATED (moved to the right-hand side in rhs_set) so the matrix stays symmetric -- CG + boomeramg need that. The diagonal can no** — `src/psi_petsc.c:490-530`
- **constant-potential walls: a Dirichlet wall makes the operator non-singular, so the constant-mode nullspace must go. Original (unconditional): MatNullSpaceCreate(...); MatSetNullSpace(a, nullsp); ...** — `src/psi_petsc.c:554-562`
- **... but not with walls: there is no nullspace then, and projecting out the mean would shift the prescribed wall potential.** — `src/psi_petsc.c:595-598`
- **constant-potential walls: induced charge** — `src/psi_petsc.c:626-637`
- **constant-potential walls. Wall site: the row is D x = b, so b = D psi_wall. Fluid site: move each eliminated wall column to the RHS, b -= A_fw psi_wall with A_fw = 2 wlaplacian[p] epsilon.** — `src/psi_petsc.c:684-718`
- **constant-potential walls: total fluid charge, for the neutrality check wall_q[0] + wall_q[1] + fluid_q = 0. With walls the e0 jump was already added link by link above (CHANGE 20261028); the per-face block below is for the periodic system only.** — `src/psi_petsc.c:731-741`
- **constant-potential (Dirichlet) walls, following Asta, Palaia, Trizac, Levesque & Rotenberg, arXiv:1907.04732, sec. II-C. MAP_BOUNDARY sites are held at a prescribed potential; the Laplacian at a fluid site doubles every link that reaches a wall site (their eq. 15), which puts the wall at the mid-plane between the last fluid and first solid node** — `src/psi_petsc.h:49-61`
- **add the force to the OWNER copy only. The cell loop includes the halo cells (0 and ncell+1), which hold periodic copies of particles in the outermost cell. That is right for the fluid reaction above: each copy scatters only to its LOCAL nodes (the node range is clipped** — `src/subgrid.c:1608-1628`
- **constant-potential walls (psi_petsc.c). The particle-node correction reaches r_cut (~6) from the particle, past the kernel support. With walls it must not use: (1) wall sites, and (2) fluid sites reached THROUGH the periodic boundary when that boundary is a wall -- their halo copy is the fluid on the other side of the electrode.** — `src/subgrid.c:5925-5958`
- **walls: see pm_sr_node_behind_wall** — `src/subgrid.c:6040-6042`

### 2026-09-24

Archivos: `util/porous_solid_cube.c`

- **spherical (microgel) mode** — `util/porous_solid_cube.c:141-144`, `util/porous_solid_cube.c:227-230`
- **spherical mode** — `util/porous_solid_cube.c:156-158`, `util/porous_solid_cube.c:191-193`, `util/porous_solid_cube.c:249-251`, `util/porous_solid_cube.c:408-416`, `util/porous_solid_cube.c:643-666`, `util/porous_solid_cube.c:1218-1220`, `util/porous_solid_cube.c:1247-1250`, `util/porous_solid_cube.c:1285-1288`
- **spherical (microgel) mode: a free-floating sphere has no walls and no faces, so both are refused rather than silently ignored. The box itself stays periodic (that is how Ludwig runs it); only the BONDING is made non-periodic, below.** — `util/porous_solid_cube.c:297-313`
- **the RSA minimum separation is now independent of the bonding window. Coupling them (min = lbond - delta) ties the monomer density to the bond tolerance, which caps the coordination number: a narrow bond window then forces a sparse packing. Default unchanged.** — `util/porous_solid_cube.c:350-356`
- **spherical mode: bonds never wrap for an isolated object, and the sphere is guaranteed to fit inside the box above.** — `util/porous_solid_cube.c:361-368`
- **spherical mode: sample in the sphere's bounding box and reject outside it, so the density applies to the sphere's volume.** — `util/porous_solid_cube.c:386-395`
- **the 3x3x3 cell scan is only exhaustive if the cell is at least as wide as every range searched: bonding (lbond + delta) and rejection (min_distance), which are no longer the same thing.** — `util/porous_solid_cube.c:436-442`
- **spherical mode: reject the bounding-box corners** — `util/porous_solid_cube.c:844-846`
- **the lower bond bound is lbond - delta, which used to be the same number as min_distance. Now that the RSA separation is set independently, it must be spelled out, or a smaller min_distance would silently widen the bond window (and pre-stress every short bond).** — `util/porous_solid_cube.c:922-928`
- **dmin, not min_distance -- see above** — `util/porous_solid_cube.c:959-961`
- **spherical mode: added R: and n** — `util/porous_solid_cube.c:1256-1258`

### 2026-09-23

Archivos: `src/ludwig.c`, `src/subgrid.h`

- **Peskin stretched to support 8** — `src/ludwig.c:1065-1067`, `src/ludwig.c:3431-3433`
- **Peskin 4-point function stretched to support 8** — `src/subgrid.h:96-98`

### 2026-09-22

Archivos: `src/ludwig.c`, `src/psi_force.c`, `src/subgrid.c`, `src/subgrid.h`, `util/PM_self_force.c`

- **truncated Gaussian kernel parameters from input. subgrid_gauss_support n full support width (lattice units, > 1, need not be integer). Default 6. subgrid_gauss_sigma sigma Gaussian width. Default n/6 (cut at 3 sigma). Called from both kernel parses (halo sizing and run setup) so that** — `src/ludwig.c:161-190`
- **truncated Gaussian kernel** — `src/ludwig.c:1063-1067`, `src/ludwig.c:3429-3433`
- **truncated Gaussian kernel parameters** — `src/ludwig.c:1125-1127`, `src/subgrid.c:102-107`
- **Gaussian support must be known before sizing the halo** — `src/ludwig.c:3460-3462`
- **Gaussian kernel** — `src/psi_force.c:511-515`, `src/psi_force.c:663-667`, `src/subgrid.c:534-541`, `src/subgrid.c:765-768`, `src/subgrid.c:981-988`, `src/subgrid.c:1590-1597`, `src/subgrid.c:1834-1841`, `src/subgrid.c:2556-2563`, `src/subgrid.c:5380-5383`, `util/PM_self_force.c:257-262`
- **Gaussian kernel range = ceil(n/2)** — `src/subgrid.c:3178-3181`
- **truncated Gaussian spread/gather kernel** — `src/subgrid.c:4228-4290`, `src/subgrid.h:71-76`
- **kernel shape parameter in the cache key. The "hann" slot of the header holds the kernel's shape parameter: the Hann order, or the support width for the Gaussian (whose sigma goes into the file name). Other kernels keep writing the Hann order, so existing cache files remain valid.** — `src/subgrid.c:6450-6466`
- **truncated Gaussian kernel id** — `src/subgrid.h:93-95`
- **Gaussian kernel in PM self-force utility** — `util/PM_self_force.c:105-107`
- **Gaussian kernel: argv[2] = support width n (default 6), argv[3] = sigma (default n/6).** — `util/PM_self_force.c:160-174`

### 2026-09-21

Archivos: `src/psi_petsc.c`, `src/subgrid.c`

- **DMDA stencil width = reach of the Poisson stencil, not the lattice halo. The width only sizes the matrix preallocation, (2w+1)^3 entries per row for DMDA_STENCIL_BOX. The Poisson matrix is assembled from cv[p] with components in {-1,0,1} for D3Q7/19/27, so w=1 is exact. Passing nhalo** — `src/psi_petsc.c:363-381`
- **no DeviceToHost after the solve. Every solver leaves its result on the HOST: the FFT solver copies DeviceToHost itself at the end of its solve (and re-uploads), and PETSc (psi_solver_petsc_da_to_psi) writes the host array directly. The main loop accordingly does no copy after solve(). The copy that** — `src/subgrid.c:6393-6404`

### 2026-09-20

Archivos: `src/ludwig.c`, `util/microgel_poly_cross_density.c`

- **output frequency for the per-particle csv files particle_force.csv and particle_Esub.csv were written every step, which is one line per colloid per step: ~370 MB over a 30k-step, 127-monomer run, for data normally sampled far more coarsely. "pforce_io_freq n" writes them only when (step % n) == 0, plus step 1. File scope rather than a local of** — `src/ludwig.c:192-207`
- **output frequency for the per-particle csv files See pforce_io_freq_ at file scope for the rationale.** — `src/ludwig.c:1199-1204`
- **(sin descripción)** — `util/microgel_poly_cross_density.c:230-234`, `util/microgel_poly_cross_density.c:910-915`, `util/microgel_poly_cross_density.c:983-1019`, `util/microgel_poly_cross_density.c:1063-1066`
- **modos de carga y semilla** — `util/microgel_poly_cross_density.c:885-889`
- **modos de carga. 0 = todos los monomeros; 1 = solo los de superficie; 2 = una fraccion charge_frac elegida al azar entre los interiores.** — `util/microgel_poly_cross_density.c:1075-1090`

### 2026-09-19

Archivos: `src/subgrid.c`

- **minimum image convention. The raw coordinate difference is wrong for pairs straddling a periodic boundary: two colloids genuinely 1.0 apart across the edge of a box of side 48 gave dx = -47, which falls outside r_cut, so the correction switched off silently exactly where it is largest (~26%** — `src/subgrid.c:6793-6806`

### 2026-09-18

Archivos: `src/ludwig.c`, `src/subgrid.c`, `src/subgrid.h`

- **measured mesh reference for the pair table. "pm_sr_pair_ref measured" builds the reference by running the mesh (deposit + Poisson solve + gather) instead of the analytic model.** — `src/ludwig.c:1473-1494`
- **measured mesh reference for the pair table pm_sr_mesh_pair_esub_x Reference term computed by RUNNING THE MESH, not by modelling it: deposits** — `src/subgrid.c:6319-6404`
- **measured-reference pair table** — `src/subgrid.c:6542-6632`
- **measured mesh reference (pm_sr_pair_ref measured). Same table, but the reference term is obtained by running the mesh — deposit two test charges, solve Poisson with the simulation's solver, gather the field back with the kernel — instead of the analytic phi_PM_pair model, which underestimates the real mesh pair force by ~9% at r=1. Serial only.** — `src/subgrid.h:234-243`

### 2026-09-17

Archivos: `src/subgrid.c`

- **linear interpolation in the pair table. Previously the nearest bin was used as-is (ir = (int)(r/dr)), which makes the correction piecewise constant in r. Measured with non-integer separations (hann8, r=0.5..3.0 step 0.1): the bin value departs from F(r) by up to ~4% below r=1, where the correction** — `src/subgrid.c:6810-6832`

### 2026-09-16

Archivos: `src/subgrid.c`, `src/subgrid.h`

- **continuum kernel-shaped target Gauss-Legendre nodes/weights on [-1,1] (Newton iteration on P_n).** — `src/subgrid.c:6231-6317`
- **continuum kernel-shaped target setup** — `src/subgrid.c:6680-6716`
- **selectable pair-correction target shape. target_shape: 0 = Gaussian (analytic erf form, original behaviour) 1 = continuum self-convolution of the real kernel, scaled so that its per-particle width equals sigma. Same shape family the mesh actually produces, so the near-field** — `src/subgrid.h:221-232`

### 2026-09-12

Archivos: `src/subgrid.c`

- **fuga de memoria: los tres klein_t de entry->force[] se asignan en subgrid.c (3 sizeof(klein_t)) pero solo se liberaba el array de punteros, perdiendo 3 x 32 B por nodo cubierto y por paso de tiempo. Original: if (c->entries[i]->force != NULL) free(c->entries[i]->force);** — `src/subgrid.c:1370-1379`

### 2026-09-10

Archivos: `src/ewald_charge.c`

- **minimum image for the dual real-space kernels. These kernels computed node<->particle and particle<->particle displacements as raw differences, with no periodic wrapping, so a pair separated across a box face was seen at distance ~L and dropped by the cutoff. Measured as a violation of translation invariance: shifting the whole system by an integer** — `src/ewald_charge.c:9543-9556`
- **minimum image (see MIN_IMAGE_DUAL above). Correct for a single rank; an MPI build would also need the neighbour-rank particle list, which this path does not build.** — `src/ewald_charge.c:9595-9601`, `src/ewald_charge.c:9727-9733`, `src/ewald_charge.c:9790-9796`
- **use the halo instead of discarding it. nhalo is sized to the Ewald cutoff (ludwig.c:3093, nhalo = ceil(ewald_rc) = irc) and is refreshed every step by psi_halo_rho() (ludwig.c:1732), so the halo already holds exactly the periodic images this loop needs; the addressing below is** — `src/ewald_charge.c:10000-10016`
- **minimum image for the particle-particle pair. Only bites when the pair separation exceeds L/2, but kept for consistency with the node loop above.** — `src/ewald_charge.c:10056-10062`
- **sign of reciprocal field on particles. Was "+=", which made the measured pp force equal A (F_real - F_recip) instead of A (F_real + F_recip): confirmed to 0.1-0.9% with no fitted parameters, and at large alpha the force converged to MINUS bare Coulomb. E = -grad(phi) = Gk k (Sk_cos sin - Sk_sin cos) = -Gk k im_part, which is** — `src/ewald_charge.c:10313-10326`
- **d_ltot was never uploaded by this driver, so the minimum-image macro in the real-space kernels would read zeros.** — `src/ewald_charge.c:11044-11047`
- **Ewald split convention: drop reciprocal form factors. The form factors were applied at BOTH source and receiver, so the reciprocal delivered erf(r/sqrt(sigma_eff^2 + 1/eta^2))/r, while the real-space kernels subtract erf(eta r)/r — the complement of a reciprocal WITHOUT form factors. The two halves therefore did not sum** — `src/ewald_charge.c:11079-11096`

### 2026-09-02

Archivos: `src/ludwig.c`, `src/subgrid.c`, `src/subgrid.h`

- **direct colloid-colloid PM short-range correction** — `src/ludwig.c:1002-1004`, `src/ludwig.c:1455-1494`, `src/ludwig.c:2092-2095`
- **direct colloid-colloid PM short-range correction pm_sr_apply_force_correction (above) can only ever see the ionic cloud: by the time it runs, the caller has already restored psi->rho to fluid-only (see ludwig.c), so it never corrects the discretisation error** — `src/subgrid.c:6135-6317`
- **direct colloid-colloid PM short-range correction Corrects the discrete PM (kernel-on-mesh) interaction directly BETWEEN nearby colloids (not colloid-vs-ion-cloud, which is what pm_sr_apply_force_correction above does). No screening (kappa) anywhere: both the target (two Gaussian charges of width sigma, pure Coulomb) and** — `src/subgrid.h:209-232`

### 2026-08-23

Archivos: `src/ludwig.c`, `src/subgrid.c`, `src/subgrid.h`

- **per-particle force output for theory comparison** — `src/ludwig.c:1195-1197`, `src/ludwig.c:2060-2064`, `src/ludwig.c:2097-2101`, `src/ludwig.c:2442-2446`, `src/subgrid.c:2172-2232`, `src/subgrid.h:134-136`
- **per-particle force output for theory comparison particle_Esub.csv records the field felt at each subgrid particle, but comparing the simulated force against the theoretical (screened) Coulomb force between colloids requires the force actually applied to each particle (pc->fex), not just the field. This mirrors the** — `src/ludwig.c:1240-1263`
- **per-particle force output for theory comparison fex now holds the final averaged F = 1/2(F_A + F_B) for this step.** — `src/ludwig.c:1769-1774`

### 2026-07-29

Archivos: `src/ewald_charge.c`, `src/ewald_charge.h`, `src/fe_force_method.c`, `src/fe_force_method.h`, `src/ludwig.c`, `src/psi_options.h`, `src/psi_solver_ewald.c`, `src/psi_solver_ewald.h`

- **Gaussian_Ewald_Dual - per-species widths and precomputed σ_eff for pair interactions** *(fecha de git)* — `src/ewald_charge.c:8954-8960`
- **Gaussian_Ewald_Dual — Dual-sigma kernels (Paso 3) Each kernel below is the analogue of its "_gaussian" counterpart but uses the appropriate σ_eff for the pair type: _pp suffix → particle-particle pair → d_gauss_sigma_eff_pp** *(fecha de git)* — `src/ewald_charge.c:9499-9541`
- **Gaussian_Ewald_Dual (Paso 3b) — Fourier kernels with receiver form factor The structure factor S(k) is built with f_p(k) for particle charges and f_f(k) for lattice charges (Paso 3 already does this). But when computing the field/potential AT a node or AT a particle, we must also multiply by** *(fecha de git)* — `src/ewald_charge.c:10091-10326`
- **Gaussian_Ewald_Dual Same as ewald_charge_sum_full_gaussian_gpu but with independent Gaussian widths for particles (σ_p) and fluid nodes (σ_f).** *(fecha de git)* — `src/ewald_charge.c:10943-11047`
- **Gaussian_Ewald_Dual - dual-sigma Gaussian Ewald sum Each particle is a Gaussian of width sigma_p; each fluid node is a Gaussian of width sigma_f. Pair interactions use sigma_eff = sqrt(sigma_i^2 + sigma_j^2): particle-particle : sqrt(2) sigma_p fluid-fluid : sqrt(2) sigma_f** *(fecha de git)* — `src/ewald_charge.h:156-167`
- **Gaussian_Ewald_Dual** *(fecha de git)* — `src/fe_force_method.c:63-67`, `src/fe_force_method.h:35-37`, `src/ludwig.c:3393-3403`, `src/ludwig.c:3858-3862`, `src/psi_options.h:36-38`, `src/psi_solver_ewald.c:172-240`, `src/psi_solver_ewald.h:67-91`
- **EWALD** *(fecha de git)* — `src/ludwig.c:995-997`
- **20260715d DETERMINISTIC 2-grid interlacing. Grid A at u=0 (= exact mode 0: no FFT shift; its psi feeds NP, so the ion cloud is perfectly stable). Grid B at u=0.5 (spectral). Forces are AVERAGED F = 1/2(F_A + F_B) on BOTH fluid and particle BEFORE the fluid integrates them: the u=0.5 grid cancels the first** *(fecha de git)* — `src/ludwig.c:1687-1774`
- **Update particle field before correctin phi** *(fecha de git)* — `src/ludwig.c:1871-1882`
- **Gaussian_Ewald: read sigma from input if present. Single-sigma mode key: ewald_gaussian_sigma** *(fecha de git)* — `src/ludwig.c:3551-3561`
- **Gaussian_Ewald_Dual: per-species sigmas Input keys: ewald_gaussian_sigma_p (particle) and ewald_gaussian_sigma_f (fluid). Defaults to ewald_gauss_sigma.** *(fecha de git)* — `src/ludwig.c:3563-3579`

### 2026-07-15

Archivos: `src/psi_fft.c`

- **diagnostic: L2 norm of the Nyquist-x plane of a field. Returns sqrt(sum over k with mx=nx/2 of |field_hat(k)|^2) / N, i.e. the real-space RMS amplitude carried by the Nyquist-x plane. Used to trace where u-dependent Nyquist content enters/leaves rho and psi. Standalone FFT (own cached plan) so it does not disturb the shift buffers.** — `src/psi_fft.c:1440-1499`

### 2026-07-14

Archivos: `src/ludwig.c`

- **for the interlacing random-offset hash** — `src/ludwig.c:22-24`
- **interlacing offset from input + random per-step mode. Input keys (all optional, defaults below): interlacing_mode 0|1 (overrides the static above) interlacing_offset <value> | random interlacing_seed <int> (used only for "random")** — `src/ludwig.c:212-233`
- **interlacing configuration from input** — `src/ludwig.c:1131-1157`

### 2026-07-10

Archivos: `src/ludwig.c`, `src/psi_fft.c`, `src/psi_fft.h`, `src/psi_force.c`, `src/subgrid.c`, `src/subgrid.h`

- **un-shift psi back to the original frame. The offset cycle solves Poisson with all charge placed at x + u, so the resulting potential lives in the shifted frame: Psi(y) = psi_true(y - u). Downstream consumers (Nernst-Planck, diagnostics) work in the original frame, so we recover psi_true(i) = Psi(i + u) by linear interpolation in X** — `src/ludwig.c:612-655`
- **spectral (FFT phase) shift of a scalar lattice field psi_fft_shift_field Translates a periodic scalar lattice field by a fractional displacement** — `src/psi_fft.c:1227-1499`
- **spectral shift stub (no CUDA)** — `src/psi_fft.c:1567-1594`
- **spectral (FFT phase) shift of a scalar lattice field** — `src/psi_fft.h:143-151`
- **trilinear: pointwise E at integer offsets** — `src/psi_force.c:668-672`
- **trilinear: exact relabelling at integer offsets** — `src/subgrid.c:759-761`
- **trilinear (interpolating) kernel** — `src/subgrid.c:3175-3177`
- **interpolating kernel for shifted-grid fluid ops. Trilinear is a delta at integer offsets, so fluid scatter/gather with an integer mesh shift reduces to an exact node relabelling (required for the integer-shift identity test of the interlacing machinery).** — `src/subgrid.h:87-92`

### 2026-07-06

Archivos: `util/PM_self_force.c`

- **Hann kernel in PM self-force utility** — `util/PM_self_force.c:102-104`, `util/PM_self_force.c:147-174`, `util/PM_self_force.c:252-262`

### 2026-06-30

Archivos: `src/ludwig.c`, `src/psi_rt.c`, `src/psi_solver.c`, `src/psi_solver_options.c`, `src/psi_solver_options.h`, `src/subgrid.c`, `src/subgrid.h`

- **Hann spread/gather kernel** — `src/ludwig.c:1060-1062`, `src/subgrid.h:63-69`
- **Hann kernel order (n = support width) from input. n must be a positive integer: partition of unity holds for any integer n (the FT of the truncated raised cosine has zeros at all non-zero integer frequencies). Odd n is allowed; it just makes the number of active nodes vary between n and n-1 with the sub-grid offset (even n gives a more** — `src/ludwig.c:1075-1094`
- **grow nhalo to cover the subgrid kernel support. The particle<->mesh scatter/gather reads/writes halo nodes up to krange = subgrid_get_range(kernel) deep. nhalo must be >= krange or the coupling accesses stale/incorrect halo data (silently wrong results). The kernel (and Hann order) are read here, before cs_nhalo_set, and set** — `src/ludwig.c:3413-3433`
- **select FFT variant by solver_type solver_type "fft" (not ewald) => solve the full Poisson equation with the discrete 27-point stencil, consistent with the FD stencil used elsewhere and directly comparable to the PETSc solver. The Ewald reciprocal-only path (create_ewald) is kept commented for** — `src/ludwig.c:3634-3652`
- **FFT Laplacian variant selectable from input. electrokinetics_fft_laplacian = discrete (default, matches PETSc) | analytic** — `src/psi_rt.c:263-282`
- **FFT Laplacian variant selectable from input via electrokinetics_fft_laplacian (default discrete = matches PETSc). ANALYTIC uses continuum k^2 and differs from PETSc by O(h^2).** — `src/psi_solver.c:65-73`
- **default FFT Laplacian = discrete (matches PETSc)** — `src/psi_solver_options.c:128-130`
- **FFT Laplacian variant selectable from input** — `src/psi_solver_options.h:35-44`
- **FFT Laplacian variant (only used when psolver==FFT)** — `src/psi_solver_options.h:54-56`
- **Hann spread/gather kernel order (configurable)** — `src/subgrid.c:94-100`
- **Hann kernel** — `src/subgrid.c:529-533`, `src/subgrid.c:762-764`, `src/subgrid.c:976-980`, `src/subgrid.c:1585-1589`, `src/subgrid.c:1829-1833`, `src/subgrid.c:2551-2555`, `src/subgrid.c:5377-5379`
- **Hann kernel range = ceil(n/2)** — `src/subgrid.c:3172-3174`
- **Fix d_peskin6 node-segment mapping. Previous version derived (r, segment) per integer band of s with a piecewise r definition, which broke partition of unity (Sum != 1) and symmetry (d6(s) != d6(-s)) at integer offsets (particle exactly on a node). That injected ~2.7% charge error wherever the particle / fluid** — `src/subgrid.c:4018-4033`
- **Hann (raised-cosine) spread/gather kernel** — `src/subgrid.c:4138-4225`
- **Hann kernel id** — `src/subgrid.h:84-86`

### 2026-06-29

Archivos: `src/ludwig.c`

- **select subgrid kernel from input file** — `src/ludwig.c:1051-1062`

### 2026-06-25

Archivos: `src/ludwig.c`, `src/subgrid.c`, `src/subgrid.h`

- **PM short-range correction tables** — `src/ludwig.c:998-1001`, `src/ludwig.c:1362-1403`
- **PM short-range phi correction (option C step 2)** — `src/ludwig.c:1883-1905`
- **PM short-range force/field correction** — `src/ludwig.c:2080-2090`
- **PM short-range correction table** — `src/subgrid.c:5322-5379`, `src/subgrid.h:179-191`
- **pm_sr_correct_phi / pm_sr_correct_efield** — `src/subgrid.c:5816-5958`

### 2026-06-11

Archivos: `src/ewald_charge.c`

- **Macro for dphi_dr with sigma->0 limit guard Used by all _dual real-space kernels. When sigma<1e-14 (point charge limit): dg_sig = d/dr[1/r] = -1/r^2 (Coulomb) otherwise use the standard Gaussian derivative formula. DPHI_DR_DUAL(dphi_dr_, sigma_, dist_, r2_, r_inv_) sets dphi_dr_ given** — `src/ewald_charge.c:9510-9541`
- **sigma->0 guard for _pf potential** — `src/ewald_charge.c:9606-9611`
- **Handle sigma_ff=0 limit: erf(r/0)/r -> 1/r** — `src/ewald_charge.c:9667-9671`
- **Force on fluid nodes from lattice nodes — fluid-fluid pair → σ_eff_ff** — `src/ewald_charge.c:9881-9955`
- **Real-space force from lattice nodes (fluid-fluid pair, σ_eff_ff)** — `src/ewald_charge.c:11301-11305`

### 2026-06-10

Archivos: `src/ludwig.c`

- **dual solver writes to ewald_charge_forces.csv, not particle_Esub.csv** — `src/ludwig.c:1275-1288`

### 2026-05-18

Archivos: `src/subgrid.c`, `util/PM_self_force.c`

- **fix sprintf: remove dangling PB format specifiers** — `src/subgrid.c:2634-2645`
- **PM self-force utility using Peskin4 spread/gather** — `util/PM_self_force.c:40-104`

### 2026-05-14

Archivos: `src/colloids.h`, `src/ewald_charge.c`, `src/ewald_charge.h`, `src/fe_force_method.c`, `src/fe_force_method.h`, `src/ludwig.c`, `src/nernst_planck.c`, `src/psi.c`, `src/psi.h`, `src/psi_fft.c`, `src/psi_fft_pn.c`, `src/psi_force.c`, `src/psi_options.h`, `src/psi_solver_ewald.c`, `src/psi_solver_ewald.h`, `src/subgrid.c`

- **Subgrid particle electric field interpolated fron mesh** *(fecha de git)* — `src/colloids.h:80` (sin CHANGE END)
- **Particle-mesh force** *(fecha de git)* — `src/colloids.h:81` (sin CHANGE END)
- **Subgrid particle force correction** *(fecha de git)* — `src/colloids.h:82` (sin CHANGE END)
- **PoissonVerification - Arrays for face sampling, Poisson check, Maxwell stress force** *(fecha de git)* — `src/ewald_charge.c:82-89`
- **PoissonVerification - Forward declarations for face-sampling helper functions** *(fecha de git)* — `src/ewald_charge.c:98-107`
- **PoissonVerification - Allocate face-sampling arrays** *(fecha de git)* — `src/ewald_charge.c:212-231`
- **PoissonVerification - Free face-sampling arrays** *(fecha de git)* — `src/ewald_charge.c:343-349`
- **SyncSinCos_GPU - Copy Sk to global sinx_/cosx_ so ewald_efield_at_r (CPU) can use the Fourier contribution. In the GPU path these are never populated otherwise.** *(fecha de git)* — `src/ewald_charge.c:4183-4188`
- **ForceFluidOutput - Store Ewald force density in psi->force_fluid for output** *(fecha de git)* — `src/ewald_charge.c:4513-4517`
- **ForceFluidOutput - Sync force_fluid back to device** *(fecha de git)* — `src/ewald_charge.c:4578-4582`
- **FaceDiv_GPU - Deduplicated face-sampling for div(E) Each unique face is evaluated exactly once by one thread. Unique faces for a window of size ni×nj×nk_win: ⊥X: (ni+1) × nj × nk_win faces, flat id = fi nj nk_win + fj nk_win + fk** *(fecha de git)* — `src/ewald_charge.c:4761-4864`
- **EvalE_Halo - include halo nodes so border face points get full real-space sum** *(fecha de git)* — `src/ewald_charge.c:4856-4864`
- **PoissonForceFunction - Poisson verification via GPU face-sampling (replaces CPU path)** *(fecha de git)* — `src/ewald_charge.c:5165-5403`
- **ForceFluidOutput_PoissonForce - Write force_fluid to file at same frequency as efield** *(fecha de git)* — `src/ewald_charge.c:5401-5403`
- **Gaussian_Ewald - Gaussian Ewald sum implementation All kernels and the main function are copies of the standard Ewald versions with suffix _gaussian. The physics changes are** *(fecha de git)* — `src/ewald_charge.c:8816-8960`
- **EfieldLatticGaussian - Gaussian real-space efield nodo->nodo** *(fecha de git)* — `src/ewald_charge.c:9258-9343`
- **EfieldLatticGaussian - Gaussian efield nodo->nodo** *(fecha de git)* — `src/ewald_charge.c:10751-10756`
- **PoissonVerification_FlagNodes - Diagnose nodes near particle faces** *(fecha de git)* — `src/ewald_charge.c:11486-11536`
- **PoissonVerification_EfieldAtPoint - E field at arbitrary continuous position** *(fecha de git)* — `src/ewald_charge.c:11553-11635`
- **RealSpace_ExactCenter - Use floor(r) as loop origin instead of round(r). round(r) was correct only when r coincides with a node center (integer coords). For face-sampling points r lies on a cell face (non-integer), so round(r) shifts the search box and misses nodes within rc on the opposite side. floor(r) + irc+1 guarantees all integer nodes within rc from the exact point r are visited.** *(fecha de git)* — `src/ewald_charge.c:11576-11635`
- **RhoUnits_EfieldAtR - psi_rho returns NUMBER density; multiply by eunit_ to get charge density [e/Δx³]. Colloid charges (pc->s.q0 etc.) are already in charge units, so this makes ions consistent.** *(fecha de git)* — `src/ewald_charge.c:11631-11635`
- **EfieldAtR_Sign - correct sign of Fourier E field** *(fecha de git)* — `src/ewald_charge.c:11721-11727`
- **PoissonVerification_FaceSampling - Sample E field on cube faces (CPU)** *(fecha de git)* — `src/ewald_charge.c:11763-11780`
- **FaceSample_Window - Restrict sampling to ±WINDOW nodes around domain centre** *(fecha de git)* — `src/ewald_charge.c:11772-11780`
- **PoissonVerification_CalcDivAndForce - div(E), Poisson error, Maxwell stress force** *(fecha de git)* — `src/ewald_charge.c:11921-11938`
- **CalcDiv_Window - Same window as face_sample so we only read valid E_face data** *(fecha de git)* — `src/ewald_charge.c:11933-11938`
- **RhoUnits_CalcDiv - psi_rho returns NUMBER density; multiply by eunit_ to get charge density. rho_eff must be in charge units [e/Δx³] so that rho_scaled = beta eunit rho_eff/eps matches div(E ) in LB units where E = beta eunit E_phys.** *(fecha de git)* — `src/ewald_charge.c:11977-11982`
- **ColoidAssignment - Distribute charge over all voxels containing the colloid A colloid exactly on a face/edge/corner is shared by 2/4/8 voxels equally. Count how many voxels claim this particle, then divide q by that count.** *(fecha de git)* — `src/ewald_charge.c:11992-12002`
- **DivEUnits_CalcDiv - corrected unit comparison for face-sampled div_E** *(fecha de git)* — `src/ewald_charge.c:12009-12012`
- **PoissonOnly - Skip Maxwell stress when EWALD_POISSON_ONLY is defined** *(fecha de git)* — `src/ewald_charge.c:12014-12084`
- **PoissonVerification_Output - Write Poisson and force data to files** *(fecha de git)* — `src/ewald_charge.c:12108-12126`
- **WritePoissonWindow - Only write nodes that were sampled** *(fecha de git)* — `src/ewald_charge.c:12121-12126`
- **RhoUnits_WritePoissonVer - psi_rho returns NUMBER density; multiply by eunit_ to get charge density** *(fecha de git)* — `src/ewald_charge.c:12144-12147`
- **ColoidAssignment - Distribute charge over all voxels containing the colloid** *(fecha de git)* — `src/ewald_charge.c:12157-12165`
- **WritePoissonVer_RelErr - add relative error** *(fecha de git)* — `src/ewald_charge.c:12169-12180`
- **WriteForceWindow - Only write nodes that were sampled** *(fecha de git)* — `src/ewald_charge.c:12201-12206`
- **PoissonStencil - Poisson check using discrete stencil Two independent verifications written to poisson_stencil_NNNNNN.dat: (A) div(E) via wgradients applied to psi->efield** *(fecha de git)* — `src/ewald_charge.c:12228-12302`
- **DivE_LBunits_Stencil - scale div_E to LB units (beta eunit div(E_phys))** *(fecha de git)* — `src/ewald_charge.c:12288-12302`
- **RhoUnits_PoissonStencil - psi_rho returns NUMBER density; multiply by eunit_ to get charge density. rho_scaled = beta eunit rho_eff/eps must match lap(psi ) and div(E ) in LB units.** *(fecha de git)* — `src/ewald_charge.c:12329-12333`
- **PoissonStencil_RelErr - add relative errors** *(fecha de git)* — `src/ewald_charge.c:12355-12370`
- **PoissonVerification - Compile-time parameters for face sampling** *(fecha de git)* — `src/ewald_charge.h:121-145`
- **Gaussian_Ewald - Gaussian Ewald sum: particles and nodes as Gaussian charge distributions** *(fecha de git)* — `src/ewald_charge.h:150-154`
- **PoissonVerification - New function: Poisson verification + Maxwell stress force** *(fecha de git)* — `src/ewald_charge.h:177-183`
- **Gaussian_Ewald** *(fecha de git)* — `src/fe_force_method.c:58-62`, `src/fe_force_method.h:32-34`, `src/ludwig.c:3382-3392`, `src/ludwig.c:3853-3857`, `src/psi_options.h:33-35`, `src/psi_solver_ewald.c:101-159`, `src/psi_solver_ewald.h:40-63`
- **Ewald sumation** *(fecha de git)* — `src/fe_force_method.h:29-31`
- **Ewald** *(fecha de git)* — `src/ludwig.c:1864-1869`
- **Subgrid charge debug data // // Charge conservation diagnostic: before Nernst-Planck** *(fecha de git)* — `src/ludwig.c:1939-1963`
- **Subgrid charge debug data TIMER_start(TIMER_ELECTRO_NPEQ); CHANGE INIT - 20260326 Use GPU-parallel NP driver** *(fecha de git)* — `src/ludwig.c:1966-1975`
- **Subgrid charge debug data // Charge conservation diagnostic: before Nernst-Planck** *(fecha de git)* — `src/ludwig.c:1978-2007`
- **Subgrid charge debug data // Necesario si no uso GPU pero esta activo el device if (!ewald)field_memcpy(ludwig->psi->efield, tdpMemcpyHostToDevice); // CHANGE END - Subgrid charge debug data psi_io_write(ludwig->psi, step);** *(fecha de git)* — `src/ludwig.c:2349-2352`
- **EWALD case FE_FORCE_METHOD_EWALD: char str_ewald_rc[BUFSIZ]; double ewald_rc; n = rt_string_parameter(rt, "ewald_rc", str_ewald_rc, BUFSIZ);** *(fecha de git)* — `src/ludwig.c:3371-3392`
- **END CHANGE INIT - Gaussian_Ewald** *(fecha de git)* — `src/ludwig.c:3381-3392`
- **Gaussian_Ewald: select standard or Gaussian Ewald wrapper** *(fecha de git)* — `src/ludwig.c:3596-3626`
- **Update charges (with Kahan summation for improved accuracy)** *(fecha de git)* — `src/nernst_planck.c:532-794`
- **psi_rho_zero** *(fecha de git)* — `src/psi.c:330-339`, `src/psi.h:106-108`
- **P3M Peskin deconvolution** *(fecha de git)* — `src/psi_fft.c:39-41`
- **eunit fix: was missing eunit factor (harmless when eunit=1)** *(fecha de git)* — `src/psi_fft.c:588-591`
- **Peskin continuous FT 1D continuous Fourier transform of the Peskin 4-point kernel: P̂_1d(k) = ∫_{-2}^{2} φ_peskin(r) cos(k·r) dr (real+symmetric) Computed analytically by integrating the two-piece formula over [-2,-1], [-1,0], [0,1], [1,2]. The 3D transform is separable: P̂(kx,ky,kz) = P̂_1d(kx)·P̂_1d(ky)·P̂_1d(kz).** *(fecha de git)* — `src/psi_fft_pn.c:198-226`
- **use continuous FT of Peskin kernel (not sinc^p approximation)** *(fecha de git)* — `src/psi_fft_pn.c:230-249`
- **Subgrid charge // klein_t flocal_k[3]; // flocal_k[X] = klein_zero(); // flocal_k[Y] = klein_zero(); // flocal_k[Z] = klein_zero();** *(fecha de git)* — `src/psi_force.c:316-321`
- **psi_force_gradmu_e fix: preserve node index** *(fecha de git)* — `src/psi_force.c:350-409`
- **subgrid_charge_from_grid** *(fecha de git)* — `src/subgrid.c:674-761`
- **Use Poisson-Boltzmann weight instead of Peskin** *(fecha de git)* — `src/subgrid.c:2355-2361`
- **Moved fprintf to after colloid_sums_halo to avoid writing halo particles** *(fecha de git)* — `src/subgrid.c:2449-2457`
- **Peskin scatter kernel on GPU** *(fecha de git)* — `src/subgrid.c:4851-5145`
- **subgrid_peskin_scatter_rho_buf** *(fecha de git)* — `src/subgrid.c:4983-5145`

### 2026-05-13

Archivos: `src/ewald_charge.c`, `src/ewald_charge.h`, `src/psi_solver_ewald.c`, `src/psi_solver_ewald.h`

- **setter for colloid info pointer** — `src/ewald_charge.c:282-310`
- **(sin descripción)** — `src/ewald_charge.h:189-193`, `src/psi_solver_ewald.c:79-89`, `src/psi_solver_ewald.c:149-159`, `src/psi_solver_ewald.h:36-38`, `src/psi_solver_ewald.h:60-63`

### 2026-05-12

Archivos: `src/ludwig.c`, `src/psi_solver.c`, `src/psi_solver_ewald.c`, `src/psi_solver_ewald.h`, `src/psi_solver_options.c`, `src/psi_solver_options.h`

- **Ewald solver wrapper** — `src/ludwig.c:121-123`
- **create Ewald psi_solver_t wrapper after ewald/self_table ready** — `src/ludwig.c:1267-1288`
- **select Poisson solver from electrokinetics_solver_type** — `src/ludwig.c:3507-3561`
- **add FFT solver to psi_solver_create** — `src/psi_solver.c:22-24`
- **FFT solver selectable from input for fe_electro** — `src/psi_solver.c:61-73`
- **new file** — `src/psi_solver_ewald.c:8-89`, `src/psi_solver_ewald.h:7-38`
- **add FFT and EWALD solver ids** — `src/psi_solver_options.c:64-71`, `src/psi_solver_options.c:96-99`, `src/psi_solver_options.h:29-32`

### 2026-04-27

Archivos: `src/subgrid.c`, `src/subgrid.h`

- **subgrid_scatter_fluid_offset: scatter fluid rho into auxiliary buffer** — `src/subgrid.c:683-761`
- **fluid scatter with offset into auxiliary buffer** — `src/subgrid.h:113-117`

### 2026-04-26

Archivos: `src/psi_force.c`

- **psi_force_gradmu_offset for interlacing** — `src/psi_force.c:79-101`
- **psi_force_gradmu_e_ewald_offset** — `src/psi_force.c:548-667`

### 2026-04-25

Archivos: `src/ludwig.c`, `src/subgrid.c`

- **Interlacing mode selector: 0 = No interlacing, use only Grid A (offset=0.0) 1 = Grid B only (offset=0.5, no interlacing) 2 = Full interlacing: average of Grid A and Grid B** — `src/ludwig.c:153-159`
- **Interlacing: one full scatter/Poisson/force cycle** — `src/ludwig.c:657-939`
- **Interlacing // DIAG: print rho_saved charge per X column (step 1 only)** — `src/ludwig.c:1596-1774`
- **Interlacing: psi_force_gradmu runs here for all modes (0, 1, 2). Uses psi->psi from interlacing_one_grid (or original Poisson for mode=0).** — `src/ludwig.c:1847-1869`
- **forward declare for interlacing offset variant** — `src/subgrid.c:109-113`
- **subgrid_charge_from_grid_offset for interlacing** — `src/subgrid.c:815-980`

### 2026-04-24

Archivos: `src/subgrid.c`

- **Kaiser-Bessel globals** — `src/subgrid.c:88-92`
- **Peskin 6-point kernel (Bao et al. 2016)** — `src/subgrid.c:3979-4033`
- **Kaiser-Bessel order-4 kernel** — `src/subgrid.c:4084-4136`

### 2026-04-22

Archivos: `src/ludwig.c`, `src/psi_fft.c`, `src/psi_fft.h`, `src/psi_force.c`, `src/psi_force.h`, `src/subgrid.c`, `src/subgrid.h`

- **assign FE_ELECTRO_EWALD id for ewald variant** — `src/ludwig.c:3509-3561`
- **B-spline deconvolution kernels** — `src/psi_fft.c:42-47`
- **deconv parameter in create_ewald** — `src/psi_fft.c:1013-1030`
- **deconv parameter in set_influence_ewald** — `src/psi_fft.c:1057-1162`, `src/psi_fft.h:126-130`
- **deconvolution kernel selection Ĝ(k) = exp(-k²/4α²) / (ε · k² · |P̂(k)|²) P̂(k) = Πᵢ Σ_m w(m) cos(kᵢ·m) separable in 3D** — `src/psi_fft.c:1120-1162`
- **Ewald deconvolution kernel selection** — `src/psi_fft.h:45-56`
- **deconv parameter in ewald** — `src/psi_fft.h:108-113`
- **psi_force_gradmu_e_ewald with kernel parameter** — `src/psi_force.c:36-41`
- **kernel parameter in psi_force_gradmu** — `src/psi_force.c:52-77`, `src/psi_force.h:22-30`
- **psi_force_gradmu_e_ewald** — `src/psi_force.c:427-515`
- **kernel parameter for subgrid_charge_from_particles** — `src/subgrid.c:446-533`
- **kernel parameter for subgrid_charge_from_grid** — `src/subgrid.c:675-761`
- **kernel parameter for subgrid_update_forces_electrokinetics** — `src/subgrid.c:1398-1449`, `src/subgrid.h:124-127`
- **kernel parameter for subgrid_update_Esub** — `src/subgrid.c:2240-2361`, `src/subgrid.h:129-133`
- **B-spline order-4 kernel** — `src/subgrid.c:3906-3935`
- **B-spline order-6 kernel** — `src/subgrid.c:3937-3977`
- **subgrid_get_drange** — `src/subgrid.c:4831-4835`, `src/subgrid.h:159-161`
- **B-spline order-4 kernel and kernel selection type** — `src/subgrid.h:52-62`

### 2026-04-18

Archivos: `src/subgrid.c`

- **Scale Peskin kernel with drange_ to support variable stencil width. drange_=1 (original, support=2): s=1, no scaling, identical to original. drange_=2 (support=4): r mapped to r/2, result divided by 2 to preserve partition-of-unity.** — `src/subgrid.c:3853-3859`
- **Scale Peskin kernel with drange_** — `src/subgrid.c:3868-3871`

### 2026-03-26

Archivos: `src/ludwig.c`, `src/nernst_planck.c`, `src/nernst_planck.h`, `src/psi_petsc.c`

- **Use GPU-parallel NP driver** — `src/ludwig.c:1969-1975`
- **GPU implementation of nernst_planck_driver_d3qx nernst_planck_driver_d3qx_gpu GPU-parallel version of nernst_planck_driver_d3qx. Uses a flat** — `src/nernst_planck.c:1209-1259`
- **Fix SEGV: target is a GPU struct, cannot deref ->data from host. Use tdpMemcpy(DeviceToHost) to read the device pointer value stored inside the GPU struct.** — `src/nernst_planck.c:1876-1886`
- **GPU driver for Nernst-Planck d3qx** — `src/nernst_planck.h:26-29`
- **Use GPU vec/mat types when -vec_type cuda is set via .petscrc. DMSetVecType/DMSetMatType must match the requested type BEFORE DMSetUp and DMCreateGlobalVector, otherwise -vec_type cuda in .petscrc is ignored and hypre PC fails with "HYPRE_MEMORY_DEVICE expects a device vector". Original (CPU only)** — `src/psi_petsc.c:383-402`
- **Allow .petscrc / command-line options for KSP solver** — `src/psi_petsc.c:422-424`
- **Fix DIVERGED_INDEFINITE_MAT with aijcusparse: MatSetNullSpace does not project the null space from the RHS/solution when using GPU matrices. Explicit projection via MatNullSpaceRemove on b and x ensures CG sees a consistent SPD system (no constant-mode contamination).** — `src/psi_petsc.c:591-598`
- **If the solution vector is on GPU (VECCUDA), bind it to CPU so that DMDAVecGetArray can access it from the host.** — `src/psi_petsc.c:1108-1111`
- **Restore GPU binding after host read.** — `src/psi_petsc.c:1133-1135`

### 2026-03-10

Archivos: `src/ewald_charge.c`

- **Self-field correction** — `src/ewald_charge.c:4632-4634`, `src/ewald_charge.c:4638-4669`

### 2026-03-09

Archivos: `src/ewald_charge.c`, `src/ewald_charge.h`

- **GPU kernels for self-field correction table** — `src/ewald_charge.c:3478-3562`
- **Self-field correction for lattice-discretisation artefact** — `src/ewald_charge.h:73-119`

### 2026-03-02

Archivos: `src/ewald_charge.c`

- **Peskin-spread FFT Ewald version** — `src/ewald_charge.c:6795-7511`

### 2026-02-12

Archivos: `src/ewald_charge.c`, `src/ewald_charge.h`

- **New ewald_charge_sum_FFT_full_gpu using cuFFT for Fourier-space** — `src/ewald_charge.c:3564-3969`, `src/ewald_charge.c:5925-6793`, `src/ewald_charge.c:7533-7542`, `src/ewald_charge.h:169-175`

### 2026-02-10

Archivos: `src/nernst_planck.c`

- **Kahan compensated sum from util_sum** — `src/nernst_planck.c:94-109`
- **Use kahan_t from util_sum.h** — `src/nernst_planck.c:954-972`

### 2026-02-06

Archivos: `src/nernst_planck.c`

- **Antisymmetric flux assignment for exact conservation** — `src/nernst_planck.c:597-794`

### 2026-01-21

Archivos: `src/psi_refined.c`, `src/psi_refined.h`

- **Refined grid Poisson solver** — `src/psi_refined.c:19` (sin CHANGE END), `src/psi_refined.h:15` (sin CHANGE END)

### 2026-01-19

Archivos: `src/ludwig.c`, `src/psi_fft.c`, `src/psi_fft.h`

- **FFT Poisson solver** — `src/ludwig.c:114-120`, `src/ludwig.c:143-159`, `src/ludwig.c:2456-2462`, `src/ludwig.c:3631-3652`
- **New FFT-based Poisson solver using cuFFT** — `src/psi_fft.c:21-41`
- **New FFT-based Poisson solver using cuFFT Edinburgh Soft Matter and Statistical Physics Group and Edinburgh Parallel Computing Centre** — `src/psi_fft.h:17-56`

### 2026-01-17

Archivos: `src/subgrid.c`, `src/subgrid.h`

- **cambio de soporte de peskin // static const double drange_ = 1.0; Max. range of interpolation - 1** — `src/subgrid.c:78-81`
- **Short-range corrections infrastructure (Level 0) See docs/SHORT_RANGE_CORRECTIONS_ANALYSIS.md for full documentation. This implements P3M-style short-range corrections for electric interactions of subgrid particles to eliminate spurious velocities and improve accuracy.** — `src/subgrid.c:115-208`
- **Short-range corrections Level 1: Self-field See docs/SHORT_RANGE_CORRECTIONS_ANALYSIS.md Section 7, Level 1. This corrects the self-field: the electric field that a particle generates on itself through its distributed charge on the mesh.** — `src/subgrid.c:3260-3542`
- **Short-range corrections declarations See docs/SHORT_RANGE_CORRECTIONS_ANALYSIS.md for full documentation. CHANGE 20260118: level1 disabled (no momentum conservation) level3_pn implements PN correction with Newton III** — `src/subgrid.h:166-177`

### 2026-01-14

Archivos: `src/collision.c`, `src/colloid.h`, `src/colloid_state_io.c`, `src/hydro.c`, `src/hydro.h`, `src/lb_d3q27.c`, `src/lb_data.c`, `src/lb_data.h`, `src/ludwig.c`, `src/psi.c`, `src/psi.h`, `src/psi_petsc.c`, `src/psi_petsc.h`, `src/subgrid.c`, `tests/unit/test_colloid.c`, `tests/unit/test_colloid_state_io.c`, `tests/unit/tests.c`

- **D3Q27 Kahan functions** *(fecha de git)* — `src/collision.c:71-74`
- **Use D3Q27 specific Kahan function** *(fecha de git)* — `src/collision.c:341-343`, `src/collision.c:589-591`, `src/collision.c:841-849`, `src/collision.c:1051-1059`
- **Kahan compensated summation for f2mode (D3Q27 and others)** *(fecha de git)* — `src/collision.c:345-382`, `src/collision.c:851-892`
- **Kahan compensated summation for mode2f (D3Q27 and others)** *(fecha de git)* — `src/collision.c:593-630`, `src/collision.c:1061-1101`
- **D3Q27 Kahan summation implementation** *(fecha de git)* — `src/collision.c:5532-5539`
- **D3Q27 constant definitions** *(fecha de git)* — `src/collision.c:5535-5539`
- **D3Q27 Kahan summation for f2mode** *(fecha de git)* — `src/collision.c:5546-9322`
- **D3Q27 Kahan summation for mode2f** *(fecha de git)* — `src/collision.c:9325-9964`
- **Extra bonds for microgel polymer networks (6 bonds, intpad reduced from 7 to 3) // #define NPAD_INT 7 // #define NPAD_DBL 4 // #define NBOND_MAX 2 #define NPAD_INT 3** *(fecha de git)* — `src/colloid.h:50-57`
- **Extra bonds for microgel polymer networks (6 bonds, intpad reduced from 7 to 3)** *(fecha de git)* — `src/colloid_state_io.c:74-122`, `src/colloid_state_io.c:236-283`
- **Resta velocidad media del sistema** *(fecha de git)* — `src/hydro.c:39-44`, `src/hydro.h:61-63`
- **Resta velocidad media del sistema Calcula la velocidad media del sistema y la resta de todos los puntos. Esto elimina cualquier movimiento neto del sistema.** *(fecha de git)* — `src/hydro.c:757-830`
- **Kernel para acumular velocidades Acumula la suma total de velocidades en el sistema.** *(fecha de git)* — `src/hydro.c:1021-1113`
- **Predictor-corrector: copy and extrapolate velocity** *(fecha de git)* — `src/hydro.c:1152-1233`, `src/hydro.h:64-67`
- **For Kahan summation** *(fecha de git)* — `src/lb_d3q27.c:22-84`
- **Use Kahan summation for normalisers like D3Q19** *(fecha de git)* — `src/lb_d3q27.c:75-84`
- **Kahan summation for equilibrium initialization** *(fecha de git)* — `src/lb_data.c:820-873`
- **Kahan summation for 0th moment (density)** *(fecha de git)* — `src/lb_data.h:254-276`
- **Kahan summation for 1st moment calculation** *(fecha de git)* — `src/lb_data.h:296-350`
- **PETSc solver with subgrid** *(fecha de git)* — `src/ludwig.c:111-113`
- **Subgrid charge debug data** *(fecha de git)* — `src/ludwig.c:2435-2446`
- **Stencil 27 subgrid // nhalo = 1; nhalo = 3; // CHANGE END - Stencil 27 subgrid psi_method = PSI_FORCE_GRADMU;** *(fecha de git)* — `src/ludwig.c:3365-3368`
- **Subgrid charge // int psi_electroneutral(psi_t psi, map_t map) { int psi_electroneutral(psi_t psi, map_t map, colloids_info_t cinfo, double rho_el) { // CHANGE END - Subgrid charge int ic, jc, kc, index;** *(fecha de git)* — `src/psi.c:935-938`
- **Subgrid charge int ncell[3]; double qcoll = 0.0; colloid_t p_colloid = NULL; colloids_info_ncell(cinfo, ncell);** *(fecha de git)* — `src/psi.c:990-1060`
- **Subgrid charge // int psi_electroneutral(psi_t obj, map_t map); int psi_electroneutral(psi_t psi, map_t map,colloids_info_t cinfo, double rho_el); // CHANGE END - Subgrid charge** *(fecha de git)* — `src/psi.h:131-134`
- **Add subgrid particle charges to PETSc RHS** *(fecha de git)* — `src/psi_petsc.c:821-1035`, `src/psi_petsc.h:40-47`
- **Momentum conservation correction** *(fecha de git)* — `src/subgrid.c:2952-3010`
- **Extra bonds for microgel polimer networks // {10, 11}, 12, {13, 14, 15}, {16, 17, 18}, 19, // 20, 21, 22, 23, 24, 25, // {26, 27, 28, 29, 30, 31, 32}, {10, 11, 12, 13, 14, 15}, 16, {17, 18, 19},** *(fecha de git)* — `tests/unit/test_colloid.c:52-59`
- **Extra bonds for microgel polimer networks // .bond = {10, 11}, // .rng = 12, // .isfixedrxyz = {13, 14, 15}, // .isfixedvxyz = {16, 17, 18},** *(fecha de git)* — `tests/unit/test_colloid_state_io.c:182-207`
- **PETSc solver tests - Test fails when PETSc is installed // test_psi_solver_petsc_suite() // CHANGE END - PETSc solver tests - Test fails when PETSc is installed test_psi_sor_suite(); test_nernst_planck_suite();** *(fecha de git)* — `tests/unit/tests.c:140-142`

### 2025-12-03

Archivos: `src/ewald_charge.c`, `src/psi.c`, `src/psi.h`, `src/psi_options.c`, `src/psi_options.h`, `src/psi_rt.c`, `src/stencil_d3q19.c`, `src/stencil_d3q27.c`, `tests/unit/test_psi.c`

- **Electric field output** — `src/ewald_charge.c:2909-3151`, `src/ewald_charge.c:4331-4392`, `src/ewald_charge.c:4407-4458`, `src/ewald_charge.c:4478-4498`, `src/ewald_charge.c:4500-4504`, `src/ewald_charge.c:4519-4540`, `src/ewald_charge.c:4549-4555`, `src/ewald_charge.c:4557-4577`, `src/psi.c:26-28`, `src/psi.c:136-142`, `src/psi.c:164-170`, `src/psi.c:1118-1158`, `src/psi.c:1173-1179`, `src/psi.c:1187-1193`, `src/psi.c:1196-1202`, `src/psi.c:1207-1227`, `src/psi.h:48-54`, `src/psi.h:100-102`, `src/psi_options.c:53-59`, `src/psi_options.c:63-69`, `src/psi_options.h:72-78`, `src/psi_rt.c:318-324`, `src/psi_rt.c:328-334`
- **Fix Laplacian normalization for consistency** — `src/stencil_d3q19.c:59-69`, `src/stencil_d3q27.c:63-72`
- **Electric field output --- IGNORE --- test_psi_initialise(pe); test_psi_create(pe); test_psi_psi_set(pe);** — `tests/unit/test_psi.c:43` (sin CHANGE END)

### 2025-11-23

Archivos: `src/subgrid.c`, `src/subgrid.h`

- **Poisson-Boltzmann kappa (inverse Debye length)** *(fecha de git)* — `src/subgrid.c:83-86`
- **Poisson-Boltzmann force calculation** *(fecha de git)* — `src/subgrid.c:2656-2931`, `src/subgrid.h:137-141`
- **Write electric field data to file** *(fecha de git)* — `src/subgrid.c:2901-2931`
- **Poisson-Boltzmann weight** *(fecha de git)* — `src/subgrid.c:4470-4835`, `src/subgrid.h:152-161`

### 2025-11-19

Archivos: `src/collision.c`, `src/distribution_rt.h`, `src/lb_d3q19.c`, `src/lb_data.c`, `src/lb_data.h`, `src/ludwig.c`, `src/psi_force.c`, `src/psi_petsc.c`, `src/stats_colloid.c`, `src/stats_colloid.h`, `src/stats_velocity.c`, `src/stats_velocity.h`, `src/subgrid.c`, `util/microgel_poly_cross_density.c`

- **Subgrid charge MODIFIED: Using Kahan version to reduce roundoff error accumulation** *(fecha de git)* — `src/collision.c:583-587`
- **Correct initial momentum to exactly zero** *(fecha de git)* — `src/distribution_rt.h:21-25`, `src/lb_data.c:34-37`, `src/lb_data.c:878-999`, `src/lb_data.h:172-178`
- **Subgrid charge // for (int ia = 0; ia < model->nvel; ia++) { // wip += model->wv[ia] model->ma[p][ia] model->ma[p][ia]; // }** *(fecha de git)* — `src/lb_d3q19.c:78-90`
- **Total force calculation** *(fecha de git)* — `src/ludwig.c:2301-2310`, `src/ludwig.c:2544-2546`, `src/ludwig.c:2613-2631`, `src/stats_velocity.c:25-28`, `src/stats_velocity.c:131-304`, `src/stats_velocity.h:18-20`, `src/stats_velocity.h:32-38`
- **Subgrid momentum** *(fecha de git)* — `src/ludwig.c:2532-2534`, `src/ludwig.c:2539-2541`, `src/ludwig.c:2560-2562`, `src/ludwig.c:2569-2574`, `src/ludwig.c:2590-2608`, `src/stats_colloid.c:63-70`, `src/stats_colloid.c:191-243`, `src/stats_colloid.h:19-21`
- **Subgrid charge // // flocal[X] += force[X]; // // flocal[Y] += force[Y]; // // flocal[Z] += force[Z]; // klein_add_double(&flocal_k[X], force[X]);** *(fecha de git)* — `src/psi_force.c:251-258`
- **Subgrid charge // flocal[X] = klein_sum(&flocal_k[X]); // flocal[Y] = klein_sum(&flocal_k[Y]); // flocal[Z] = klein_sum(&flocal_k[Z]); // // CHANGE END - Subgrid charge** *(fecha de git)* — `src/psi_force.c:264-268`
- **CUDA C++ compatibility fix** — `src/psi_petsc.c:481-488`, `src/psi_petsc.c:1417-1424`, `src/subgrid.c:1336-1340`, `util/microgel_poly_cross_density.c:1717-1721`, `util/microgel_poly_cross_density.c:1741-1745`
- **Write particle data after colloid_sums_halo, only for local (non-halo) particles** *(fecha de git)* — `src/subgrid.c:2622-2645`
- **Subgrid charge Loop through all cells (including the halo cells) and set the velocity at each particle to zero for this step.** *(fecha de git)* — `src/subgrid.c:3668-3689`
- **Subgrid charge klein_t force_k[3]; force_k[X] = klein_zero(); force_k[Y] = klein_zero(); force_k[Z] = klein_zero();** *(fecha de git)* — `src/subgrid.c:3702-3707`
- **Subgrid charge // p_colloid->fsub[X] += u[X] dr; // p_colloid->fsub[Y] += u[Y] dr; // p_colloid->fsub[Z] += u[Z] dr;** *(fecha de git)* — `src/subgrid.c:3743-3752`
- **Subgrid charge p_colloid->fsub[X] = klein_sum(&force_k[X]); p_colloid->fsub[Y] = klein_sum(&force_k[Y]); p_colloid->fsub[Z] = klein_sum(&force_k[Z]); //CHANGE END - Subgrid charge** *(fecha de git)* — `src/subgrid.c:3758-3762`

### 2025-11-10

Archivos: `src/colloid_sums.c`, `src/colloid_sums.h`, `src/ludwig.c`, `src/psi.c`

- **Subgrid charge static int colloid_sums_m8(colloid_sum_t sum, int, int, int, int); //CHANGE END - Subgrid charge Message sizes (doubles)** *(fecha de git)* — `src/colloid_sums.c:87-89`
- **Subgrid charge // static const int msize_[COLLOID_SUM_MAX] = { 0, 10, 35, 7, 4, 6, 7, 13}; static const int msize_[COLLOID_SUM_MAX] = { 0, 10, 35, 7, 4, 6, 7, 13, 4 }; //CHANGE END - Subgrid charge** *(fecha de git)* — `src/colloid_sums.c:93-96`
- **Subgrid charge if (sum->mtype == COLLOID_SUM_ELECTRIC_FIELD) mloader_forw = colloid_sums_m8; //CHANGE END - Subgrid charge assert(mloader_forw); mloader_back = mloader_forw;** *(fecha de git)* — `src/colloid_sums.c:401-403`
- **Subgrid charge colloid_sums_m8** *(fecha de git)* — `src/colloid_sums.c:877` (sin CHANGE END)
- **Subgrid charge COLLOID_SUM_MAX = 9} colloid_sum_enum_t; int colloid_sums_create(colloids_info_t cinfo, colloid_sum_t psum); void colloid_sums_free(colloid_sum_t sum);** *(fecha de git)* — `src/colloid_sums.h:36` (sin CHANGE END)
- **Subgrid charge debug data FILE fp, fpforce; double kappa_aux; CHANGE INIT - 20260823 per-particle force output for theory comparison** *(fecha de git)* — `src/ludwig.c:1192-1197`
- **Subgrid charge** *(fecha de git)* — `src/psi.c:29-31`, `src/psi.c:623-626`, `src/psi.c:646-649`, `src/psi.c:654-656`

### 2025-10-29

Archivos: `src/psi_force.c`

- **Subgrid charge klein_t flocal_k[3]; flocal_k[X] = klein_zero(); flocal_k[Y] = klein_zero(); flocal_k[Z] = klein_zero();** *(fecha de git)* — `src/psi_force.c:129-134`
- **Subgrid charge psi_electric_field(psi, index, e); // CHANGE psi->efield->data[addr_rank1(psi->efield->nsites, psi->efield->nf, index, X)] = kt e[X]; //kt e[X]; psi->efield->data[addr_rank1(psi->efield->nsites, psi->efield->nf, index, Y)] = kt e[Y];** *(fecha de git)* — `src/psi_force.c:169-229`
- **Subgrid charge psi_force_gradmu_e_subgrid** *(fecha de git)* — `src/psi_force.c:961-1024`

### 2025-07-23

Archivos: `src/ludwig.c`

- **Subgrid charge** *(fecha de git)* — `src/ludwig.c:4263-4265`, `src/ludwig.c:4272-4274`
- **(sin descripción)** *(fecha de git)* — `src/ludwig.c:4278` (sin CHANGE END)

### 2025-07-01

Archivos: `src/build.c`, `src/interaction.c`, `src/ludwig.c`, `src/subgrid.c`, `src/subgrid.h`

- **Subgrid charge** *(fecha de git)* — `src/build.c:1543-1545`, `src/interaction.c:31-35`, `src/ludwig.c:559-586`, `src/ludwig.c:2043-2064`, `src/subgrid.c:68-74`, `src/subgrid.c:346-533`, `src/subgrid.h:19-62`

### 2025-01-12

Archivos: `src/psi_petsc.c`

- **Fix NULL pointer handling for test compatibility** — `src/psi_petsc.c:581-586`

