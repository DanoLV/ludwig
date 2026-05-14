# Electrostatic FFT Solver Integration — Ludwig (develop-LaValle branch)

This document describes all functions and logic in `src/ludwig.c` (lines 53–1198) that are
relevant to solving the electrostatic problem using FFT + deconvolution (i.e., the reciprocal-
space Particle Mesh Ewald / P3M approach).  It is intended to allow a Claude Code instance to
transplant this functionality into a different Ludwig repository that has a slightly different
codebase.

---

## 1. Overview

The electrostatic solver pipeline computes:

1. **Scatter** — spread particle and fluid charge density onto the lattice mesh using a
   regularized delta-function kernel (spread / *scatter*).
2. **Poisson solve** — solve ∇²ψ = −ρ/ε in reciprocal space via FFT (cuFFT or PETSc FFT).
   Optionally apply **deconvolution** to correct for the smearing introduced by the kernel
   (making it equivalent to P3M optimal influence function).
3. **Gather forces** — compute the electric field E = −∇ψ at each particle location by
   gathering from the lattice using the same kernel (adjoint of scatter), and compute the
   fluid force via ∇μ_e.
4. **Interlacing** (optional) — average two shifted-grid solutions to cancel aliasing errors,
   analogous to the interlaced P3M approach.

---

## 2. New static globals and includes (lines 108–152)

These are additions relative to the upstream Ludwig codebase.

### 2.1 New includes (lines 108–117)

```c
/*CHANGE INIT - PETSc solver with subgrid */
#include "psi_petsc.h"
/*CHANGE END - PETSc solver with subgrid */
/*CHANGE INIT - 20260119 FFT Poisson solver */
#include "psi_fft.h"       /* FFT-based Poisson solver (psi_solver_fft_t) */
#include "psi_fft_pn.h"    /* PN short-range correction (psi_fft_pn_t) */
#include "psi_refined.h"   /* Refined-grid FFT solver */
#include "ewald_charge.h"  /* Ewald summation for charges (ewald_charge_t) */
#include "psi_exclusion.h" /* Charge exclusion / electroneutrality helpers */
/*CHANGE END - 20260119 FFT Poisson solver */
```

### 2.2 Static state variables (lines 137–152)

```c
// static psi_solver_fft_t* fft_solver_ = NULL;   /* (unused, kept for reference) */
static psi_fft_pn_t* fft_pn_ = NULL;              /* PN short-range correction object */
static ewald_charge_t* ewald = NULL;               /* Ewald summation object */
static psi_refined_t* psi_refined_ = NULL;         /* Refined-grid FFT object */
static int use_refined_grid_ = 0;   /* 0 = normal FFT, 1 = refined grid FFT */
static int refined_grid_factor_ = 2; /* Refinement factor (2, 4, 8, …) */

/* Active subgrid kernel — choose ONE of: */
subgrid_kernel_t kernel_g = SUBGRID_KERNEL_PESKIN4;
/* Other options: SUBGRID_KERNEL_BSPLINE4, SUBGRID_KERNEL_BSPLINE6,
                  SUBGRID_KERNEL_KB4, SUBGRID_KERNEL_PESKIN6 */

/* Deconvolution mode — choose ONE of: */
psi_fft_deconv_t deconv_g = PSI_FFT_DECONV_NONE;
/* Other options: PSI_FFT_DECONV_PESKIN, PSI_FFT_DECONV_BSPLINE4,
                  PSI_FFT_DECONV_BSPLINE6 */

/* Interlacing mode:
 *   0 = Grid A only (offset=0.0, no interlacing)
 *   1 = Grid B only (offset=0.5, equivalent to one shifted grid)
 *   2 = Full interlacing: average of Grid A + Grid B (commented out) */
static int interlacing_mode_ = 0;
```

These variables control the global electrostatic pipeline.  They should be declared at file
scope in `ludwig.c` in the target repository.

---

## 3. `interlacing_one_grid` function (lines 487–638)

```c
static int interlacing_one_grid(
    ludwig_t* lud,
    const double* rho_saved,   /* fluid-only charge density (no particle charge) */
    int rho_ndata,             /* = psi->nsites * psi->nk */
    double mesh_offset,        /* fractional offset in lattice units (0.0 or 0.5) */
    subgrid_kernel_t kernel,
    int step,
    FILE* fp,                  /* output file for Esub diagnostics */
    int flag)                  /* psi_force_method flag */
```

**Purpose:** One complete electrostatic cycle for a single mesh grid.  Called once per
timestep in mode 0/1, or twice (offset=0 and offset=0.5) in mode 2.

**Algorithm:**
1. Restore `psi->rho->data` from `rho_saved` (fluid-only baseline).
2. Shift particle positions by `mesh_offset` in X.
3. Scatter fluid charge onto auxiliary buffer `rho_buf` via
   `subgrid_scatter_fluid_offset(cinfo, psi, kernel, mesh_offset, rho_buf, n)`.
   Save a copy in `rho_fluid_offset`.
4. Copy `rho_buf` into `psi->rho->data`.
5. Scatter particle charges on top via
   `subgrid_charge_from_particles(cinfo, psi, &charge_sg, kernel)`.
6. `psi_halo_rho(psi)` — synchronize charge density halos.
7. Restore particle positions (undo shift).
8. Call `lud->poisson->impl->solve(lud->poisson, step)` — FFT Poisson solve.
9. Restore `psi->rho->data` to `rho_fluid_offset` (fluid-only-with-offset)
   so that the gather is the exact adjoint of the scatter.
10. Free `rho_fluid_offset` and `charge_sg`.
11. `hydro_u_halo`, `psi_halo_psi`, `psi_halo_psijump`, `psi_halo_rho`,
    `field_memcpy(psi->psi, HostToDevice)`.
12. `psi_force_gradmu_offset(psi, fe, phi, hydro, map, cinfo, kernel, mesh_offset)`
    — compute force on fluid using adjoint gather with the same offset.
13. Halo updates again, `hydro_memcpy(HostToDevice)`, `psi_zero_mean`.
14. Shift particle positions by `mesh_offset` again for gather.
15. `subgrid_update_Esub(cinfo, psi, step, fp, pe, kernel)` — gather E field at
    particle positions.
16. `subgrid_update_forces_electrokinetics(cinfo, map, phys, psi, hydro, kernel)`
    — convert E field to force on particles and fluid.
17. Restore particle positions.

**Key invariant:** The scatter (step 3) and gather (step 12 / step 15) use the *same*
`mesh_offset` so they are exact adjoints.  This is the core of the PME approach.

---

## 4. Main electrokinetic block inside `ludwig_run` (lines 853–1197)

This block executes when `ludwig->psi != NULL`.

### 4.1 Save baseline charge density (lines 860–863)

```c
int rho_ndata = ludwig->psi->nsites * ludwig->psi->nk;
double* rho_saved = (double*)malloc(rho_ndata * sizeof(double));
memcpy(rho_saved, ludwig->psi->rho->data, rho_ndata * sizeof(double));
```

`rho_saved` holds the fluid-only charge density at the start of each timestep, before any
particle contributions are scattered.  It is the baseline for both interlacing grids.

### 4.2 Mode 0 — Standard path (lines 882–900)

```c
if (interlacing_mode_ == 0) {
    distributed_charge_klein_t* charge = NULL;
    distributed_charge_klein_t* charge_sg = NULL;

    /* Scatter fluid charge from colloid grid to psi (used for Nernst-Planck) */
    subgrid_charge_from_grid(cinfo, psi, &charge, kernel_g);

    /* Scatter particle charges onto psi->rho */
    subgrid_charge_from_particles(cinfo, psi, &charge_sg, kernel_g);
    psi_halo_rho(psi);

    TIMER_start(TIMER_ELECTRO_POISSON);
    ludwig->poisson->impl->solve(ludwig->poisson, step);
    TIMER_stop(TIMER_ELECTRO_POISSON);

    subgrid_free_distributed_charge_t(&charge);
    subgrid_free_distributed_charge_t(&charge_sg);

    /* Restore rho to fluid-only (forces computed later in im==0 block) */
    memcpy(psi->rho->data, rho_saved, rho_ndata * sizeof(double));
    free(rho_saved);
}
```

In mode 0 the force calculation (`psi_force_gradmu`, `subgrid_update_Esub`,
`subgrid_update_forces_electrokinetics`) happens later at lines 1033–1196.

### 4.3 Mode 1/2 — Interlacing path (lines 901–996)

```c
else {
    int il_flag = 0;
    psi_force_method(psi, &il_flag);

    double mesh_offset_A = 1.0;  /* Grid A offset (in lattice units) */
    double mesh_offset_B = 0.0;  /* Grid B offset (TEST: both same → must equal mode 0) */

    /* Grid A: one full scatter/Poisson/force cycle */
    TIMER_start(TIMER_ELECTRO_POISSON);
    interlacing_one_grid(ludwig, rho_saved, rho_ndata,
                         mesh_offset_A, kernel_g, step, fp, il_flag);
    TIMER_stop(TIMER_ELECTRO_POISSON);

    /* NOTE: mode 2 (Grid B + average) is commented out — development in progress */

    /* Restore rho and sync forces to device */
    /* memcpy(psi->rho->data, rho_saved, rho_ndata * sizeof(double)); */
    /* if (hydro) field_memcpy(hydro->force, HostToDevice); */
}
```

**Note:** The commented-out mode 2 block would save Grid A forces, compute Grid B forces,
and average them.  It is not yet active.

### 4.4 Post-Poisson hydrodynamics (lines 999–1008)

```c
if (hydro) {
    hydro_u_halo(hydro);
    hydro_memcpy(hydro, tdpMemcpyDeviceToHost);
}
```

### 4.5 Nernst-Planck multistep loop (lines 1020–1149)

```c
psi_multisteps(psi, &multisteps);
for (im = 0; im < multisteps; im++) {
    psi_halo_psi(psi);
    psi_halo_psijump(psi);
    psi_halo_rho(psi);

    if (im == 0) {
        /* Mode 0: compute forces here */
        if (interlacing_mode_ == 0) {
            psi_force_method(psi, &flag);
            if (flag == PSI_FORCE_GRADMU) {
                psi_force_gradmu(psi, fe, phi, hydro, map, cinfo, kernel_g);
            }
            if (flag == PSI_FORCE_DIVERGENCE) {
                psi_force_divstress(psi, fe, hydro, cinfo);
            }
        }
        /* Modes 1/2: forces already computed inside interlacing_one_grid */
    }

    /* GPU-parallel Nernst-Planck driver */
    nernst_planck_driver_d3qx_gpu(psi, fe, hydro, map, cinfo);
}

psi_halo_psi(psi);
psi_halo_psijump(psi);
psi_halo_rho(psi);
if (hydro) hydro_memcpy(hydro, tdpMemcpyHostToDevice);
nernst_planck_adjust_multistep(psi);
psi_zero_mean(psi);
```

### 4.6 Mode 0 particle forces (lines 1177–1197)

```c
if (ludwig->psi) {
    if (interlacing_mode_ == 0) {
        subgrid_update_Esub(cinfo, psi, step, fp, pe, kernel_g);
        subgrid_update_forces_electrokinetics(cinfo, map, phys, psi, hydro, kernel_g);
    }
    /* Modes 1/2: already done inside interlacing_one_grid */
}
```

---

## 5. Function signatures required from other modules

The following functions are called from the electrostatic pipeline.  Each must exist
(with compatible signature) in the target repository.

### 5.1 `src/subgrid.h` / `src/subgrid.c`

```c
/* Kernel type enum */
typedef enum {
  SUBGRID_KERNEL_PESKIN4 = 0,
  SUBGRID_KERNEL_BSPLINE4 = 1,
  SUBGRID_KERNEL_BSPLINE6 = 2,
  SUBGRID_KERNEL_KB4 = 3,
  SUBGRID_KERNEL_PESKIN6 = 4
} subgrid_kernel_t;

/* Klein-compensated charge accumulation structures */
typedef struct distributed_charge_klein_entry_s { ... } distributed_charge_klein_entry_t;
typedef struct distributed_charge_klein_s { ... } distributed_charge_klein_t;

/* Scatter fluid charge from colloid lattice region into psi->rho */
int subgrid_charge_from_grid(colloids_info_t* cinfo, psi_t* obj,
                              distributed_charge_klein_t** charge,
                              subgrid_kernel_t kernel);

/* Scatter particle charge onto psi->rho using kernel */
int subgrid_charge_from_particles(colloids_info_t* cinfo, psi_t* obj,
                                   distributed_charge_klein_t** charge,
                                   subgrid_kernel_t kernel);

/* Free a distributed_charge_klein_t */
void subgrid_free_distributed_charge_t(distributed_charge_klein_t** charge);

/* Scatter fluid charge onto auxiliary buffer with mesh offset */
int subgrid_scatter_fluid_offset(colloids_info_t* cinfo, psi_t* obj,
                                  subgrid_kernel_t kernel, double mesh_offset,
                                  double* rho_buf, int ndata);

/* Gather E field at each particle's position (fills pc->Esub) */
int subgrid_update_Esub(colloids_info_t* cinfo, psi_t* psi, int step, FILE* fp,
                         pe_t* pe, subgrid_kernel_t kernel);

/* Convert pc->Esub to force on particles (pc->fex) and reaction on fluid */
int subgrid_update_forces_electrokinetics(colloids_info_t* cinfo, map_t* map,
                                           physics_t* phys, psi_t* psi,
                                           hydro_t* hydro,
                                           subgrid_kernel_t kernel);

/* Compute inverse Debye length from psi */
int subgrid_compute_kappa(psi_t* psi, double* kappa);
```

### 5.2 `src/psi_force.h` / `src/psi_force.c`

```c
/* Force on fluid via gradient of electrochemical potential (no offset) */
int psi_force_gradmu(psi_t* psi, fe_t* fe, field_t* phi, hydro_t* hydro,
                     map_t* map, colloids_info_t* cinfo, subgrid_kernel_t kernel);

/* Same as psi_force_gradmu but uses mesh_offset for the gather
 * (adjoint of subgrid_scatter_fluid_offset) */
int psi_force_gradmu_offset(psi_t* psi, fe_t* fe, field_t* phi, hydro_t* hydro,
                             map_t* map, colloids_info_t* cinfo,
                             subgrid_kernel_t kernel, double mesh_offset);

/* Alternative: force via divergence of Maxwell stress tensor */
int psi_force_divstress(psi_t* psi, fe_t* fe, hydro_t* hydro,
                         colloids_info_t* cinfo);
```

### 5.3 `src/psi_fft.h` and `src/psi_fft_pn.h`

```c
/* FFT Poisson solver object (concrete type) */
typedef struct psi_solver_fft_s psi_solver_fft_t;

/* Deconvolution mode */
typedef enum { PSI_FFT_DECONV_NONE, PSI_FFT_DECONV_PESKIN,
               PSI_FFT_DECONV_BSPLINE4, PSI_FFT_DECONV_BSPLINE6 } psi_fft_deconv_t;

/* PN short-range correction structure */
typedef struct psi_fft_pn_s psi_fft_pn_t;

int psi_fft_pn_create(psi_solver_fft_t* fft_solver, psi_t* psi, psi_fft_pn_t** pn);
int psi_fft_pn_free(psi_fft_pn_t** pn);
int psi_fft_pn_compute_green(psi_fft_pn_t* pn);
int psi_fft_pn_correction(psi_fft_pn_t* pn, colloids_info_t* cinfo, hydro_t* hydro);
int psi_fft_pn_ewald_correction_full(psi_fft_pn_t* pn, colloids_info_t* cinfo, hydro_t* hydro);
```

### 5.4 `src/nernst_planck.h`

```c
/* GPU-parallel Nernst-Planck solver (added in this branch) */
int nernst_planck_driver_d3qx_gpu(psi_t* psi, fe_t* fe, hydro_t* hydro,
                                   map_t* map, colloids_info_t* cinfo);

/* Original CPU solver (kept for reference) */
int nernst_planck_driver_d3qx(psi_t* psi, fe_t* fe, hydro_t* hydro,
                               map_t* map, colloids_info_t* cinfo);

int nernst_planck_adjust_multistep(psi_t* psi);
```

### 5.5 `src/psi.h`

```c
int psi_halo_psi(psi_t* psi);
int psi_halo_rho(psi_t* psi);
int psi_halo_psijump(psi_t* psi);
int psi_rho(psi_t* psi, int index, int n, double* rho);
int psi_zero_mean(psi_t* psi);
int psi_multisteps(psi_t* psi, int* multisteps);
int psi_force_method(psi_t* psi, int* flag);
/* flag values: PSI_FORCE_GRADMU, PSI_FORCE_DIVERGENCE */
```

### 5.6 `src/psi_solver.h`

```c
/* Polymorphic solver interface */
typedef struct psi_solver_s {
    struct psi_solver_vt_s* impl;
} psi_solver_t;

typedef struct psi_solver_vt_s {
    int (*solve)(psi_solver_t* solver, int timestep);
    int (*free)(psi_solver_t** solver);
} psi_solver_vt_t;
```

### 5.7 `src/ewald_charge.h`

```c
typedef struct ewald_charge_s ewald_charge_t;
/* Used as a static global ewald = NULL in this branch */
```

### 5.8 `src/colloids.h` (fields accessed)

```c
typedef struct colloid_s {
    colloid_state_t s;     /* includes s.r[3] = position */
    double fex[3];         /* external force (accumulated) */
    double Esub[3];        /* gathered electric field at particle */
    struct colloid_s* nextlocal;
    /* ... */
} colloid_t;

int colloids_info_local_head(colloids_info_t* cinfo, colloid_t** pc);
int colloids_info_nlocal(colloids_info_t* cinfo, int* n);
```

---

## 6. Data flow diagram

```
timestep start
    │
    ├─ memcpy rho_saved ← psi->rho->data   (fluid-only charge)
    │
    ├─ [mode 0] ─────────────────────────────────────────────────────
    │   subgrid_charge_from_grid(cinfo, psi, ...)      ← fluid scatter
    │   subgrid_charge_from_particles(cinfo, psi, ...) ← particle scatter
    │   psi_halo_rho
    │   poisson->solve(step)                           ← FFT Poisson
    │   restore psi->rho ← rho_saved
    │   [in im==0 of NP loop]:
    │     psi_force_gradmu(...)   ← fluid force via ∇μ_e
    │   [after NP loop]:
    │     subgrid_update_Esub     ← gather E at particles
    │     subgrid_update_forces_electrokinetics ← particle+fluid force
    │
    ├─ [mode 1/2] ──────────────────────────────────────────────────
    │   interlacing_one_grid(offset_A):
    │     restore rho ← rho_saved
    │     shift particles by offset_A
    │     subgrid_scatter_fluid_offset(offset_A) → rho_buf
    │     copy rho_buf → psi->rho
    │     subgrid_charge_from_particles → add particle charge
    │     psi_halo_rho
    │     poisson->solve(step)
    │     restore rho ← rho_fluid_offset
    │     psi_force_gradmu_offset(offset_A) ← fluid force
    │     shift particles by offset_A
    │     subgrid_update_Esub   ← gather E at shifted particles
    │     subgrid_update_forces_electrokinetics ← particle+fluid force
    │     restore particles
    │
    ├─ Nernst-Planck multistep loop
    │   nernst_planck_driver_d3qx_gpu × multisteps
    │
    └─ timestep end
```

---

## 7. Modified initialization (lines 446–456)

The original call to `psi_electroneutral(psi, map)` was replaced:

```c
/* CHANGE INIT - Subgrid charge */
double rho_el;
rt_double_parameter(rt, "electrokinetics_init_rho_el", &rho_el);
psi_electroneutral(ludwig->psi, ludwig->map, ludwig->collinfo, rho_el);
psi_halo_rho(ludwig->psi);
/* CHANGE END */
```

The new signature requires `colloids_info_t*` and an explicit charge density `rho_el`.

---

## 8. Modified initial conditions call (lines 329–332)

```c
/* CHANGE INIT - Correct initial momentum to exactly zero */
/* Original: lb_rt_initial_conditions(pe, rt, ludwig->lb, ludwig->phys); */
lb_rt_initial_conditions(pe, rt, ludwig->lb, ludwig->phys,
                          ludwig->map, ludwig->collinfo);
/* CHANGE END */
```

---

## 9. Diagnostic output files

Two CSV files are opened at the start of `ludwig_run`:

- `./proceced_data/ewald_charge_forces.csv` — total/fluid/particle force components
  per timestep.  Header: `|F_total|,F_total_X,...,F_particle_Z`
- `./proceced_data/particle_Esub.csv` — electric sub-field at each particle per step.
  Header: `# Step;Index;Emod;Esub_X;Esub_Y;Esub_Z;EmodPB;EPB_X;EPB_Y;EPB_Z`

Both files append if they exist.  The `fp` handle is passed to `subgrid_update_Esub` and
`interlacing_one_grid`.

---

## 10. Key constants and compile-time choices

| Variable | Location | Default | Notes |
|----------|----------|---------|-------|
| `kernel_g` | `ludwig.c` static | `SUBGRID_KERNEL_PESKIN4` | Active spread/gather kernel |
| `deconv_g` | `ludwig.c` static | `PSI_FFT_DECONV_NONE` | FFT deconvolution mode |
| `interlacing_mode_` | `ludwig.c` static | `0` | 0=standard, 1=shifted, 2=average |
| `use_refined_grid_` | `ludwig.c` static | `0` | 1=refined-grid FFT |
| `refined_grid_factor_` | `ludwig.c` static | `2` | Refinement factor |
| `PSI_FFT_PN_RC` | `psi_fft_pn.h` | `2.5` | PN correction cutoff radius |

---

## 11. Modo 0 — Flujo completo sin interlacing (`interlacing_mode_ == 0`)

Esta sección describe el camino de ejecución completo cuando `interlacing_mode_ = 0`.
Es el camino original, sin desplazamiento de grilla, y es la referencia para verificar
que los modos 1 y 2 producen resultados equivalentes.

### 11.1 Separación temporal del pipeline en dos fases

A diferencia de los modos con interlacing (donde scatter, Poisson y fuerzas ocurren
dentro de una misma función), en modo 0 el pipeline está **dividido en dos bloques
separados** del bucle principal:

- **Fase 1** (líneas 882–900): scatter de cargas + solve de Poisson.
- **Fase 2** (líneas 1034–1046 y 1183–1193): cálculo de fuerzas sobre fluido y partículas.

La separación existe porque el cálculo de fuerzas ocurre dentro del bucle `im == 0`
de Nernst-Planck, que también maneja el `psi_halo_*` necesario antes de llamar a
`psi_force_gradmu`.

### 11.2 Fase 1 — Scatter + Poisson (líneas 882–900)

```c
/* 1. Guardar rho antes de la dispersión */
int rho_ndata = psi->nsites * psi->nk;
double* rho_saved = malloc(rho_ndata * sizeof(double));
memcpy(rho_saved, psi->rho->data, rho_ndata * sizeof(double));
/* rho_saved contiene sólo la carga del fluido (sin contribución de partículas) */

/* 2. Dispersar carga del fluido desde la región de coloides hacia psi->rho
 *    (modifica psi->rho en los nodos ocupados por la partícula) */
distributed_charge_klein_t* charge = NULL;
subgrid_charge_from_grid(cinfo, psi, &charge, kernel_g);

/* 3. Dispersar carga de las partículas sobre psi->rho usando el kernel activo */
distributed_charge_klein_t* charge_sg = NULL;
subgrid_charge_from_particles(cinfo, psi, &charge_sg, kernel_g);

/* 4. Sincronizar halos de carga entre procesos MPI */
psi_halo_rho(psi);

/* 5. Resolver la ecuación de Poisson en espacio recíproco: ∇²ψ = -ρ/ε
 *    El resultado queda en psi->psi */
poisson->impl->solve(poisson, step);

/* 6. Liberar estructuras de carga dispersada */
subgrid_free_distributed_charge_t(&charge);
subgrid_free_distributed_charge_t(&charge_sg);

/* 7. IMPORTANTE: restaurar psi->rho al estado fluido-only (rho_saved)
 *    El solve de Poisson usó rho = fluido + partícula, pero las fuerzas de fluido
 *    en psi_force_gradmu deben ver sólo la carga del fluido para que el gather
 *    sea el adjunto exacto del scatter. */
memcpy(psi->rho->data, rho_saved, rho_ndata * sizeof(double));
free(rho_saved);
```

**Por qué se restaura `rho` después del solve:**
`psi_force_gradmu` calcula la fuerza como −ρ∇(δF/δρ), donde la densidad de carga ρ
que se usa en el *gather* del fluido debe corresponder exactamente a la misma distribución
que se usó para el *scatter* original.  Si se dejara la carga de las partículas incluida,
la fuerza resultante no sería el adjunto exacto del operador de scatter y se violaría la
conservación de momento.

### 11.3 Sincronización hidrodinámica (líneas 999–1006)

Tras el solve de Poisson y antes del bucle NP:

```c
hydro_u_halo(hydro);          /* halo de velocidades */
hydro_memcpy(hydro, tdpMemcpyDeviceToHost);  /* sincronizar device→host */
```

### 11.4 Fase 2a — Fuerza sobre el fluido, primera iteración NP (`im == 0`, líneas 1034–1046)

Dentro del bucle `for (im = 0; im < multisteps; im++)`, sólo en la primera iteración:

```c
psi_halo_psi(psi);       /* halo del potencial ψ */
psi_halo_psijump(psi);   /* halo de saltos de ψ en paredes */
psi_halo_rho(psi);       /* halo de densidad de carga */

/* Seleccionar método de fuerza configurado en el archivo de entrada */
psi_force_method(psi, &flag);

if (flag == PSI_FORCE_GRADMU) {
    /* Fuerza sobre fluido via gradiente del potencial electroquímico:
     *   f_fluido = -ρ ∇μ_e = -ρ ∇(eψ + kT ln ρ)
     * Usa kernel_g para el gather de ψ en cada nodo del fluido */
    psi_force_gradmu(psi, fe, phi, hydro, map, cinfo, kernel_g);
}
if (flag == PSI_FORCE_DIVERGENCE) {
    /* Alternativa: fuerza via divergencia del tensor de Maxwell */
    psi_force_divstress(psi, fe, hydro, cinfo);
}
```

En modo 0, `psi_force_gradmu` es la versión **sin offset** (offset implícito = 0.0).
Internamente realiza un gather del gradiente de ψ sobre los nodos del fluido con el mismo
kernel `kernel_g` y acumula la fuerza en `hydro->force`.

### 11.5 Bucle Nernst-Planck (líneas 1020–1149)

Todas las iteraciones `im = 0 … multisteps-1`:

```c
/* Transporte de iones en campo electrostático */
nernst_planck_driver_d3qx_gpu(psi, fe, hydro, map, cinfo);
```

Al finalizar el bucle:

```c
psi_halo_psi(psi);
psi_halo_psijump(psi);
psi_halo_rho(psi);
hydro_memcpy(hydro, tdpMemcpyHostToDevice);
nernst_planck_adjust_multistep(psi);   /* ajusta dt del multistep */
psi_zero_mean(psi);                    /* fija media de ψ a cero */
```

### 11.6 Fase 2b — Fuerza sobre partículas (líneas 1183–1193)

Fuera del bloque electrocinético, en la sección de dinámica del parámetro de orden:

```c
/* Gather del campo E en la posición de cada partícula:
 *   pc->Esub[3] ← −∇ψ evaluado en r_p via kernel_g */
subgrid_update_Esub(cinfo, psi, step, fp, pe, kernel_g);

/* Convertir pc->Esub en fuerza sobre partícula y reacción sobre fluido:
 *   pc->fex += q_p * pc->Esub
 *   hydro->force en nodos vecinos recibe -q_p * pc->Esub (Newton III) */
subgrid_update_forces_electrokinetics(cinfo, map, phys, psi, hydro, kernel_g);
```

### 11.7 Diagrama de secuencia completo — modo 0

```
inicio del paso temporal
│
├─ [1] memcpy rho_saved ← psi->rho->data      (fluido puro, sin partículas)
│
├─ [2] subgrid_charge_from_grid(...)           (scatter fluido → psi->rho)
├─ [3] subgrid_charge_from_particles(...)      (scatter partículas → psi->rho)
├─ [4] psi_halo_rho
├─ [5] poisson->solve(step)                    (FFT: ∇²ψ = -ρ/ε → psi->psi)
├─ [6] free charge structs
├─ [7] memcpy psi->rho->data ← rho_saved       (restaurar fluido puro)
│
├─ [8] hydro_u_halo
├─ [9] hydro_memcpy DeviceToHost
│
├─ bucle NP (im = 0 … multisteps-1)
│   ├─ psi_halo_psi / psi_halo_psijump / psi_halo_rho
│   ├─ [im==0] psi_force_gradmu(...)           (fuerza fluido: f += -ρ ∇μ_e)
│   └─ nernst_planck_driver_d3qx_gpu(...)      (transporte de iones)
│
├─ [10] psi_halo_* / hydro_memcpy / nernst_planck_adjust / psi_zero_mean
│
├─ [11] subgrid_update_Esub(...)               (gather E en partículas)
└─ [12] subgrid_update_forces_electrokinetics  (f_part += q*E; reacción fluido)
```

### 11.8 Diferencias clave respecto a los modos con interlacing

| Aspecto | Modo 0 (sin interlacing) | Modos 1/2 (con interlacing) |
|---------|--------------------------|------------------------------|
| Offset de partículas | Ninguno (posición real) | Se desplaza ±`mesh_offset` en X |
| Scatter de fluido | `subgrid_charge_from_grid` directo | `subgrid_scatter_fluid_offset` con offset |
| Función de fuerza fluido | `psi_force_gradmu` (sin offset) | `psi_force_gradmu_offset` (con offset) |
| Cuándo se calculan fuerzas | En `im==0` del bucle NP (fase 2a) y después del bucle NP (fase 2b) | Dentro de `interlacing_one_grid`, antes del bucle NP |
| Restauración de `rho` | Inmediata tras el solve (paso 7) | Usa `rho_fluid_offset` separado |
| Aliasing del kernel | Presente (un solo grid) | Reducido (promedio de dos grids desplazados) |

---

## 12. Files to copy / adapt

To transplant the electrostatic FFT pipeline into another Ludwig repository the following
source files are required (beyond the unmodified upstream files):

| File | Description |
|------|-------------|
| `src/psi_fft.c` / `.h` | FFT Poisson solver (cuFFT-based) |
| `src/psi_fft_pn.c` / `.h` | PN short-range correction, Ewald real-space |
| `src/psi_refined.c` / `.h` | Refined-grid FFT solver |
| `src/ewald_charge.c` / `.h` | Ewald summation for charges |
| `src/psi_exclusion.c` / `.h` | Charge exclusion helpers |
| `src/subgrid.c` / `.h` | Subgrid kernels, scatter/gather (heavily modified) |
| `src/psi_force.c` / `.h` | Force with mesh_offset variant |
| `src/nernst_planck.c` / `.h` | GPU NP driver |
| `src/psi.c` / `.h` | Modified `psi_electroneutral` signature |
| `src/util_sum.c` / `.h` | Klein compensated summation |

The changes to `src/ludwig.c` itself are the static globals (section 2) and the
electrokinetic block in `ludwig_run` (section 4), wrapped in `/*CHANGE INIT*/` /
`/*CHANGE END*/` markers.

---

## 13. Interpolation kernel functions (`src/subgrid.c`)

These functions implement the regularized delta-function kernels φ(r) used for both
scatter (charge assignment) and gather (field interpolation).  **All kernels must satisfy
the partition-of-unity property**: Σ_j φ(r_p − j) = 1 for any particle position r_p.
The 3D kernel is the tensor product: φ(r) = φ(r_x) φ(r_y) φ(r_z).

### 13.1 `d_peskin` — 4-point Peskin kernel (support [-2, 2])

```c
double d_peskin(double r) {
    double rmod;
    double delta = 0.0;

    rmod = fabs(r);

    if (rmod <= 1.0) {
        delta = 0.125 * (3.0 - 2.0 * rmod + sqrt(1.0 + 4.0 * rmod - 4.0 * rmod * rmod));
    }
    else if (rmod <= 2.0) {
        delta = 0.125 * (5.0 - 2.0 * rmod - sqrt(-7.0 + 12.0 * rmod - 4.0 * rmod * rmod));
    }

    return delta;
}
```

### 13.2 `d_peskin_derivative` — derivative of Peskin kernel

```c
double d_peskin_derivative(double x) {
    double r = fabs(x / drange_);
    double sign;
    double val = 0.0;

    if (x > 0.0) sign = 1.0;
    else if (x < 0.0) sign = -1.0;
    else sign = 0.0;

    if (r <= 1.0) {
        double tmp = sqrt(1.0 + 4.0 * r - 4.0 * r * r);
        val = 0.125 * (-2.0 + (4.0 - 8.0 * r) / (2.0 * tmp));
    }
    else if (r <= 2.0) {
        double tmp = sqrt(-7.0 + 12.0 * r - 4.0 * r * r);
        val = 0.125 * (-2.0 - (12.0 - 8.0 * r) / (2.0 * tmp));
    }

    return val * sign;
}
```

### 13.3 `d_bspline4` — cubic B-spline, support [-2, 2]

```c
double d_bspline4(double r) {
    double t = fabs(r);
    double val = 0.0;

    if (t < 1.0) {
        val = 2.0 / 3.0 - t * t + 0.5 * t * t * t;
    }
    else if (t < 2.0) {
        double u = 2.0 - t;
        val = u * u * u / 6.0;
    }

    return val;
}
```

### 13.4 `d_bspline6` — quintic B-spline, support [-3, 3]

```c
double d_bspline6(double r) {
    double t = fabs(r);
    double val = 0.0;

    if (t < 1.0) {
        double t2 = t * t;
        val = 11.0 / 20.0 - t2 / 2.0 + t2 * t2 / 4.0 - t2 * t2 * t / 12.0;
    }
    else if (t < 2.0) {
        double u = t - 1.0;
        val = 13.0 / 60.0 + u * (-5.0 / 12.0 + u * (1.0 / 6.0 + u * (1.0 / 6.0 + u * (-1.0 / 6.0 + u * (1.0 / 24.0)))));
    }
    else if (t < 3.0) {
        double u = 3.0 - t;
        double u2 = u * u;
        val = u2 * u2 * u / 120.0;
    }

    return val;
}
```

### 13.5 `d_peskin6` — Bao et al. 2016 6-point C³ kernel, support [-3, 3]

K = 59/60 − √29/20 ≈ 0.7141.  Particle at position x_p; node offset s = x_node − x_p.

```c
double d_peskin6(double s) {
    if (fabs(s) >= 3.0) return 0.0;

    static const double K = 0.71407520893979593; /* 59/60 - sqrt(29)/20 */

    double r;
    int segment;
    if (s > -3.0 && s <= -2.0) { r = -2.0 - s; segment = 0; }
    else if (s > -2.0 && s <= -1.0) { r = -1.0 - s; segment = 1; }
    else if (s > -1.0 && s <= 0.0)  { r = -s;        segment = 2; }
    else if (s >  0.0 && s <= 1.0)  { r = 1.0 - s;   segment = 3; }
    else if (s >  1.0 && s <= 2.0)  { r = 2.0 - s;   segment = 4; }
    else                             { r = 3.0 - s;   segment = 5; }

    double beta = 9.0/4.0 - 1.5*(K + r*r) + (22.0/3.0 - 7.0*K)*r - (7.0/3.0)*r*r*r;
    double t1 = (3.0*K - 1.0)*r + r*r*r;
    double t2 = (4.0 - 3.0*K)*r - r*r*r;
    double gamma_r = -11.0/32.0*r*r + 3.0/32.0*(2.0*K + r*r)*r*r
                   + t1*t1/72.0 + t2*t2/18.0;
    double phi_m3 = (-beta + sqrt(beta*beta - 112.0*gamma_r)) / 56.0;

    switch (segment) {
    case 0: return phi_m3;
    case 1: return -3.0*phi_m3 - 1.0/16.0 + (K + r*r)/8.0
                   + (3.0*K - 1.0)*r/12.0 + r*r*r/12.0;
    case 2: return  2.0*phi_m3 + 1.0/4.0 + (4.0 - 3.0*K)*r/6.0 - r*r*r/6.0;
    case 3: return  2.0*phi_m3 + 5.0/8.0 - (K + r*r)/4.0;
    case 4: return -3.0*phi_m3 + 1.0/4.0 - (4.0 - 3.0*K)*r/6.0 + r*r*r/6.0;
    case 5: return  phi_m3     - 1.0/16.0 + (K + r*r)/8.0
                   - (3.0*K - 1.0)*r/12.0 - r*r*r/12.0;
    default: return 0.0;
    }
}
```

### 13.6 `d_kb4` — Kaiser-Bessel kernel W=4, support [-2, 2]

```c
static double i0_series(double x) {
    double ax = fabs(x);
    double y, ans;
    if (ax < 3.75) {
        y = x / 3.75;
        y *= y;
        ans = 1.0 + y*(3.5156229 + y*(3.0899424 + y*(1.2067492
            + y*(0.2659732 + y*(0.0360768 + y*0.0045813)))));
    } else {
        y = 3.75 / ax;
        ans = (exp(ax) / sqrt(ax)) * (0.39894228 + y*(0.01328592
            + y*(0.00225319 + y*(-0.00157565 + y*(0.00916281
            + y*(-0.02057706 + y*(0.02635537 + y*(-0.01647633
            + y*0.00392377))))))));
    }
    return ans;
}

/* subgrid_kb4_beta_ is a file-scope static double, set via subgrid_set_kb4_beta() */
/* subgrid_kb4_norm_fluid_ is precomputed = (sum_m phi(m))^3 for normalisation */

void subgrid_set_kb4_beta(double beta) {
    subgrid_kb4_beta_ = beta;
    double i0b = i0_series(beta);
    double s = 0.0;
    for (int m = -2; m <= 2; m++) {
        double t = fabs((double)m);
        if (t < 2.0) s += i0_series(beta * sqrt(1.0 - (m/2.0)*(m/2.0))) / (2.0 * i0b);
    }
    subgrid_kb4_norm_fluid_ = s * s * s;
}

double d_kb4(double r) {
    double t = fabs(r);
    if (t >= 2.0) return 0.0;
    double arg = sqrt(1.0 - (r / 2.0) * (r / 2.0));
    return i0_series(subgrid_kb4_beta_ * arg) / (2.0 * i0_series(subgrid_kb4_beta_));
}
```

**Note on KB4 normalisation for fluid nodes:**  `d_kb4` does NOT satisfy partition of unity
on integer grids, so a normalisation factor `subgrid_kb4_norm_fluid_` is applied in scatter
and gather: `dr = d_kb4(dx)*d_kb4(dy)*d_kb4(dz) / subgrid_kb4_norm_fluid_`.

### 13.7 `d_trilinear` — trilinear, support [-1, 1]

```c
double d_trilinear(double r) {
    double rmod = fabs(r);
    double weight = 0.0;
    if (rmod <= 1.0) { weight = 1.0 - rmod; }
    return weight;
}
```

---

## 14. Range and index utilities (`src/subgrid.c`)

### 14.1 `subgrid_get_range`

Returns the integer half-support radius for the chosen kernel.

```c
int subgrid_get_range(subgrid_kernel_t kernel) {
    if      (kernel == SUBGRID_KERNEL_BSPLINE6) return 2;
    else if (kernel == SUBGRID_KERNEL_BSPLINE4) return 1;
    else if (kernel == SUBGRID_KERNEL_PESKIN4)  return 1;
    else if (kernel == SUBGRID_KERNEL_KB4)      return 1;
    else if (kernel == SUBGRID_KERNEL_PESKIN6)  return 2;
    else                                        return drange_;
}
```

### 14.2 `subgrid_get_lattice_index_range`

Clips to interior `[1, nlocal]`; used for particle scatter/gather.

```c
void subgrid_get_lattice_index_range(double r0[3], int range, int nlocal[3],
    int* i_min, int* i_max, int* j_min, int* j_max, int* k_min, int* k_max) {
    *i_min = imax(1, (int)floor(r0[X] - range));
    *i_max = imin(nlocal[X], (int)ceil(r0[X] + range));
    *j_min = imax(1, (int)floor(r0[Y] - range));
    *j_max = imin(nlocal[Y], (int)ceil(r0[Y] + range));
    *k_min = imax(1, (int)floor(r0[Z] - range));
    *k_max = imin(nlocal[Z], (int)ceil(r0[Z] + range));
}
```

### 14.3 `subgrid_get_lattice_index_range_halo`

Extends range into halos (`± 2*range`); used for fluid scatter/gather where the kernel
support of border nodes crosses the periodic boundary.

```c
void subgrid_get_lattice_index_range_halo(double r0[3], int range, int nlocal[3],
    int* i_min, int* i_max, int* j_min, int* j_max, int* k_min, int* k_max) {
    *i_min = imax(1 - 2*range, (int)floor(r0[X] - range));
    *i_max = imin(nlocal[X] + 2*range, (int)ceil(r0[X] + range));
    *j_min = imax(1 - 2*range, (int)floor(r0[Y] - range));
    *j_max = imin(nlocal[Y] + 2*range, (int)ceil(r0[Y] + range));
    *k_min = imax(1 - 2*range, (int)floor(r0[Z] - range));
    *k_max = imin(nlocal[Z] + 2*range, (int)ceil(r0[Z] + range));
}
```

---

## 15. Charge scatter — particle to mesh (`src/subgrid.c`)

### 15.1 Klein accumulator helpers

`add_charge_to_array` maintains a sorted array of `(cs_index, rho0_sum, rho1_sum)` entries
with Klein compensated summation and binary-search insertion:

```c
static int binary_search_charge_index(distributed_charge_klein_t* charge, int cs_index) {
    int left = 0, right = charge->count - 1;
    while (left <= right) {
        int mid = left + (right - left) / 2;
        if (charge->entries[mid]->cs_index == cs_index) return mid;
        if (charge->entries[mid]->cs_index < cs_index) left = mid + 1;
        else right = mid - 1;
    }
    return left;
}

void add_charge_to_array(distributed_charge_klein_t** charge_ptr, int cs_index,
                         double q0_dr, double q1_dr, psi_t* obj) {
    distributed_charge_klein_t* charge = *charge_ptr;

    if (charge == NULL) {
        charge = (distributed_charge_klein_t*)malloc(sizeof(distributed_charge_klein_t));
        charge->entries = (distributed_charge_klein_entry_t**)malloc(16 * sizeof(distributed_charge_klein_entry_t*));
        charge->count = 0;
        charge->capacity = 16;
        *charge_ptr = charge;
    }

    int pos = binary_search_charge_index(charge, cs_index);

    if (pos < charge->count && charge->entries[pos]->cs_index == cs_index) {
        klein_add_double(charge->entries[pos]->rho0_sum, q0_dr);
        klein_add_double(charge->entries[pos]->rho1_sum, q1_dr);
    }
    else {
        if (charge->count >= charge->capacity) {
            charge->capacity *= 2;
            charge->entries = (distributed_charge_klein_entry_t**)realloc(charge->entries,
                               charge->capacity * sizeof(distributed_charge_klein_entry_t*));
        }
        for (int i = charge->count; i > pos; i--) charge->entries[i] = charge->entries[i - 1];

        distributed_charge_klein_entry_t* entry = (distributed_charge_klein_entry_t*)malloc(sizeof(distributed_charge_klein_entry_t));
        entry->cs_index = cs_index;
        psi_rho(obj, cs_index, 0, &entry->rho0_original);
        psi_rho(obj, cs_index, 1, &entry->rho1_original);
        entry->rho0_sum = (klein_t*)malloc(sizeof(klein_t));
        entry->rho1_sum = (klein_t*)malloc(sizeof(klein_t));
        *entry->rho0_sum = klein_zero();
        *entry->rho1_sum = klein_zero();
        klein_add_double(entry->rho0_sum, q0_dr);
        klein_add_double(entry->rho1_sum, q1_dr);
        charge->entries[pos] = entry;
        charge->count++;
    }
}
```

### 15.2 `subgrid_charge_from_particles` — scatter particle charge to mesh

For each subgrid particle, distributes charge `q0`, `q1` onto surrounding lattice nodes
using kernel weights.  Stores original node values so they can be restored afterward.

```c
int subgrid_charge_from_particles(colloids_info_t* cinfo, psi_t* obj,
                                   distributed_charge_klein_t** charge,
                                   subgrid_kernel_t kernel) {
    int ic, jc, kc, i, j, k, i_min, i_max, j_min, j_max, k_min, k_max, index;
    int nlocal[3], offset[3], ncell[3];
    double r[3], r0[3], dr, q0_dr, q1_dr;
    colloid_t* p_colloid = NULL;

    if (cinfo->nsubgrid == 0) return 0;

    int krange = subgrid_get_range(kernel);
    cs_nlocal(cinfo->cs, nlocal);
    cs_nlocal_offset(cinfo->cs, offset);
    colloids_info_ncell(cinfo, ncell);

    for (ic = 0; ic <= ncell[X] + 1; ic++) {
    for (jc = 0; jc <= ncell[Y] + 1; jc++) {
    for (kc = 0; kc <= ncell[Z] + 1; kc++) {
        colloids_info_cell_list_head(cinfo, ic, jc, kc, &p_colloid);
        for (; p_colloid; p_colloid = p_colloid->next) {
            if (p_colloid->s.bc != COLLOID_BC_SUBGRID) continue;

            r0[X] = p_colloid->s.r[X] - 1.0 * offset[X];
            r0[Y] = p_colloid->s.r[Y] - 1.0 * offset[Y];
            r0[Z] = p_colloid->s.r[Z] - 1.0 * offset[Z];

            subgrid_get_lattice_index_range(r0, krange, nlocal,
                &i_min, &i_max, &j_min, &j_max, &k_min, &k_max);

            /* Precompute KB4 weights and normalisation if needed */
            double kb_w[4][4][4] = {{{0}}};
            double kb_w_sum = 1.0;
            if (kernel == SUBGRID_KERNEL_KB4) {
                kb_w_sum = 0.0;
                for (i = i_min; i <= i_max; i++)
                for (j = j_min; j <= j_max; j++)
                for (k = k_min; k <= k_max; k++) {
                    double wx = d_kb4(r0[X] - i);
                    double wy = d_kb4(r0[Y] - j);
                    double wz = d_kb4(r0[Z] - k);
                    kb_w[i-i_min][j-j_min][k-k_min] = wx * wy * wz;
                    kb_w_sum += wx * wy * wz;
                }
            }

            for (i = i_min; i <= i_max; i++) {
            for (j = j_min; j <= j_max; j++) {
            for (k = k_min; k <= k_max; k++) {
                index = cs_index(cinfo->cs, i, j, k);

                r[X] = r0[X] - (double)i;
                r[Y] = r0[Y] - (double)j;
                r[Z] = r0[Z] - (double)k;

                if      (kernel == SUBGRID_KERNEL_BSPLINE6)
                    dr = d_bspline6(r[X]) * d_bspline6(r[Y]) * d_bspline6(r[Z]);
                else if (kernel == SUBGRID_KERNEL_BSPLINE4)
                    dr = d_bspline4(r[X]) * d_bspline4(r[Y]) * d_bspline4(r[Z]);
                else if (kernel == SUBGRID_KERNEL_PESKIN4)
                    dr = d_peskin(r[X]) * d_peskin(r[Y]) * d_peskin(r[Z]);
                else if (kernel == SUBGRID_KERNEL_KB4)
                    dr = kb_w[i-i_min][j-j_min][k-k_min] / kb_w_sum;
                else if (kernel == SUBGRID_KERNEL_PESKIN6)
                    dr = d_peskin6(r[X]) * d_peskin6(r[Y]) * d_peskin6(r[Z]);

                q0_dr = p_colloid->s.q0 * dr;
                q1_dr = p_colloid->s.q1 * dr;
                add_charge_to_array(charge, index, q0_dr, q1_dr, obj);
            }}}
        }
    }}}

    /* Apply accumulated charges to psi */
    if (*charge != NULL) {
        for (int i = 0; i < (*charge)->count; i++) {
            distributed_charge_klein_entry_t* entry = (*charge)->entries[i];
            double new_rho0 = entry->rho0_original + klein_sum(entry->rho0_sum);
            double new_rho1 = entry->rho1_original + klein_sum(entry->rho1_sum);
            psi_rho_set(obj, entry->cs_index, 0, new_rho0);
            psi_rho_set(obj, entry->cs_index, 1, new_rho1);
        }
    }
    return 0;
}
```

### 15.3 `subgrid_charge_from_particles_restore`

Restores `psi->rho` to the values stored in the `rho*_original` fields when the charge
structure was built (i.e., before particle charge was added).

```c
int subgrid_charge_from_particles_restore(colloids_info_t* cinfo, psi_t* obj,
                                           distributed_charge_klein_t** charge) {
    if (*charge == NULL) return 0;
    for (int i = 0; i < (*charge)->count; i++) {
        distributed_charge_klein_entry_t* entry = (*charge)->entries[i];
        psi_rho_set(obj, entry->cs_index, 0, entry->rho0_original);
        psi_rho_set(obj, entry->cs_index, 1, entry->rho1_original);
    }
    return 0;
}
```

### 15.4 `subgrid_free_distributed_charge_t`

```c
void subgrid_free_distributed_charge_t(distributed_charge_klein_t** charge) {
    if (charge == NULL || *charge == NULL) return;
    distributed_charge_klein_t* c = *charge;
    for (int i = 0; i < c->count; i++) {
        if (c->entries[i]) {
            if (c->entries[i]->rho0_sum) free(c->entries[i]->rho0_sum);
            if (c->entries[i]->rho1_sum) free(c->entries[i]->rho1_sum);
            free(c->entries[i]);
        }
    }
    if (c->entries) free(c->entries);
    free(c);
    *charge = NULL;
}
```

---

## 16. Fluid charge scatter with mesh offset (`src/subgrid.c`)

### 16.1 `subgrid_scatter_fluid_offset`

Scatters fluid charge from `psi->rho` into an auxiliary buffer `rho_buf` using
`kernel` with an X-axis shift of `mesh_offset`.  Used for interlacing: mode 2 calls this
twice (offset=0 and offset=0.5) and averages the two Poisson solutions.  `psi->rho` is
NOT modified.

**Key design points:**
- For non-integer `mesh_offset`, border interior nodes whose kernel support crosses the
  periodic boundary are replaced with the corresponding halo nodes (`psi_halo_rho` must
  have been called beforehand).
- Destination indices outside `[1, nlocal]` are wrapped to interior via periodic BC.
- `rho_buf` is zeroed by the caller, size `= nsites * nk`.

```c
int subgrid_scatter_fluid_offset(colloids_info_t* cinfo, psi_t* obj,
                                  subgrid_kernel_t kernel, double mesh_offset,
                                  double* rho_buf, int ndata) {
    int i, j, k, i2, j2, k2;
    int i_min, i_max, j_min, j_max, k_min, k_max;
    int nlocal[3];

    cs_nlocal(cinfo->cs, nlocal);
    int krange = subgrid_get_range(kernel);

    double frac_x = mesh_offset - floor(mesh_offset);

    /* Determine which border interior nodes to replace with halo nodes */
    int ix_excl_lo = 0, ix_excl_hi = -1;
    int ix_halo_lo = 1, ix_halo_hi = 0;
    if (frac_x > 0.0) {
        ix_excl_lo = nlocal[X] - krange + 1; ix_excl_hi = nlocal[X];
        ix_halo_lo = 1 - krange;             ix_halo_hi = 0;
    } else if (frac_x < 0.0) {
        ix_excl_lo = 1;               ix_excl_hi = krange;
        ix_halo_lo = nlocal[X] + 1;  ix_halo_hi = nlocal[X] + krange;
    }

    int list_cap = nlocal[X] + 2 * krange + 4;
    int* i_list = (int*)malloc(list_cap * sizeof(int)); int i_list_n = 0;
    for (i = 1; i <= nlocal[X]; i++) {
        if (i < ix_excl_lo || i > ix_excl_hi) i_list[i_list_n++] = i;
    }
    for (i = ix_halo_lo; i <= ix_halo_hi; i++) i_list[i_list_n++] = i;

    for (int ii = 0; ii < i_list_n; ii++) {
        i = i_list[ii];
        /* Wrap halo source index to interior to read rho */
        int iw = i;
        if (iw < 1) iw += nlocal[X]; else if (iw > nlocal[X]) iw -= nlocal[X];

    for (j = 1; j <= nlocal[Y]; j++) {
    for (k = 1; k <= nlocal[Z]; k++) {
        int index_src = cs_index(cinfo->cs, iw, j, k);
        double rho0, rho1;
        psi_rho(obj, index_src, 0, &rho0);
        psi_rho(obj, index_src, 1, &rho1);
        if (rho0 == 0.0 && rho1 == 0.0) continue;

        /* Source position on shifted grid */
        double r0[3] = { (double)i + mesh_offset, (double)j, (double)k };
        subgrid_get_lattice_index_range_halo(r0, krange, nlocal,
            &i_min, &i_max, &j_min, &j_max, &k_min, &k_max);

        for (i2 = i_min; i2 <= i_max; i2++) {
        for (j2 = j_min; j2 <= j_max; j2++) {
        for (k2 = k_min; k2 <= k_max; k2++) {
            double dx = r0[X] - (double)i2;
            double dy = r0[Y] - (double)j2;
            double dz = r0[Z] - (double)k2;
            double dr;
            if      (kernel == SUBGRID_KERNEL_BSPLINE6)
                dr = d_bspline6(dx)*d_bspline6(dy)*d_bspline6(dz);
            else if (kernel == SUBGRID_KERNEL_BSPLINE4)
                dr = d_bspline4(dx)*d_bspline4(dy)*d_bspline4(dz);
            else if (kernel == SUBGRID_KERNEL_PESKIN4)
                dr = d_peskin(dx)*d_peskin(dy)*d_peskin(dz);
            else if (kernel == SUBGRID_KERNEL_KB4)
                dr = d_kb4(dx)*d_kb4(dy)*d_kb4(dz)/subgrid_kb4_norm_fluid();
            else if (kernel == SUBGRID_KERNEL_PESKIN6)
                dr = d_peskin6(dx)*d_peskin6(dy)*d_peskin6(dz);
            else dr = 0.0;
            if (dr == 0.0) continue;

            /* Periodic wrap for destination */
            int idw = i2, jdw = j2, kdw = k2;
            if (idw < 1) idw += nlocal[X]; else if (idw > nlocal[X]) idw -= nlocal[X];
            if (jdw < 1) jdw += nlocal[Y]; else if (jdw > nlocal[Y]) jdw -= nlocal[Y];
            if (kdw < 1) kdw += nlocal[Z]; else if (kdw > nlocal[Z]) kdw -= nlocal[Z];
            int index_dst = cs_index(cinfo->cs, idw, jdw, kdw);

            rho_buf[addr_rank1(obj->nsites, obj->nk, index_dst, 0)] += rho0 * dr;
            rho_buf[addr_rank1(obj->nsites, obj->nk, index_dst, 1)] += rho1 * dr;
        }}}
    }}}
    free(i_list);
    return 0;
}
```

### 16.2 `subgrid_charge_from_grid`

Thin wrapper — scatter fluid charge with zero offset (no interlacing):

```c
int subgrid_charge_from_grid(colloids_info_t* cinfo, psi_t* obj,
                              distributed_charge_klein_t** charge,
                              subgrid_kernel_t kernel) {
    return subgrid_charge_from_grid_offset(cinfo, obj, charge, kernel, 0.0);
}
```

---

## 17. E-field gather from mesh to particles (`src/subgrid.c`)

### 17.1 `subgrid_update_Esub`

Gather the electric field E = −∇ψ from surrounding lattice nodes onto each subgrid
particle.  The result is stored in `pc->Esub[3]`.  Uses the same kernel as scatter
(adjoint operator — required for Newton III / momentum conservation).

```c
int subgrid_update_Esub(colloids_info_t* cinfo, psi_t* psi, int step, FILE* fp,
                         pe_t* pe, subgrid_kernel_t kernel) {
    int i, j, k, ic, jc, kc, i_min, i_max, j_min, j_max, k_min, k_max;
    int ncell[3], index;
    int nlocal[3], offset[3];
    double dr, r[3], r0[3], e[3];
    double E_field[3] = { 0.0, 0.0, 0.0 };
    colloid_t* pc;
    MPI_Comm comm;

    if (cinfo->nsubgrid == 0) return 0;

    cs_nlocal(cinfo->cs, nlocal);
    cs_nlocal_offset(cinfo->cs, offset);
    cs_cart_comm(cinfo->cs, &comm);
    colloids_info_ncell(cinfo, ncell);

    int krange = subgrid_get_range(kernel);

    for (ic = 0; ic <= ncell[X] + 1; ic++) {
    for (jc = 0; jc <= ncell[Y] + 1; jc++) {
    for (kc = 0; kc <= ncell[Z] + 1; kc++) {
        colloids_info_cell_list_head(cinfo, ic, jc, kc, &pc);
        for (; pc; pc = pc->next) {
            if (pc->s.bc != COLLOID_BC_SUBGRID) continue;

            klein_t E_field_k[3];
            E_field_k[X] = klein_zero();
            E_field_k[Y] = klein_zero();
            E_field_k[Z] = klein_zero();

            r0[X] = pc->s.r[X] - 1.0 * offset[X];
            r0[Y] = pc->s.r[Y] - 1.0 * offset[Y];
            r0[Z] = pc->s.r[Z] - 1.0 * offset[Z];

            subgrid_get_lattice_index_range(r0, krange, nlocal,
                &i_min, &i_max, &j_min, &j_max, &k_min, &k_max);

            /* Precompute KB4 weights */
            double kb_w[4][4][4] = {{{0}}};
            double kb_w_sum = 1.0;
            if (kernel == SUBGRID_KERNEL_KB4) {
                kb_w_sum = 0.0;
                for (i = i_min; i <= i_max; i++)
                for (j = j_min; j <= j_max; j++)
                for (k = k_min; k <= k_max; k++) {
                    double wx = d_kb4(r0[X] - i);
                    double wy = d_kb4(r0[Y] - j);
                    double wz = d_kb4(r0[Z] - k);
                    kb_w[i-i_min][j-j_min][k-k_min] = wx * wy * wz;
                    kb_w_sum += wx * wy * wz;
                }
            }

            for (i = i_min; i <= i_max; i++) {
            for (j = j_min; j <= j_max; j++) {
            for (k = k_min; k <= k_max; k++) {
                index = cs_index(cinfo->cs, i, j, k);
                r[X] = r0[X] - (double)i;
                r[Y] = r0[Y] - (double)j;
                r[Z] = r0[Z] - (double)k;

                if      (kernel == SUBGRID_KERNEL_BSPLINE6)
                    dr = d_bspline6(r[X])*d_bspline6(r[Y])*d_bspline6(r[Z]);
                else if (kernel == SUBGRID_KERNEL_BSPLINE4)
                    dr = d_bspline4(r[X])*d_bspline4(r[Y])*d_bspline4(r[Z]);
                else if (kernel == SUBGRID_KERNEL_PESKIN4)
                    dr = d_peskin(r[X])*d_peskin(r[Y])*d_peskin(r[Z]);
                else if (kernel == SUBGRID_KERNEL_KB4)
                    dr = kb_w[i-i_min][j-j_min][k-k_min] / kb_w_sum;
                else if (kernel == SUBGRID_KERNEL_PESKIN6)
                    dr = d_peskin6(r[X])*d_peskin6(r[Y])*d_peskin6(r[Z]);

                psi_electric_field(psi, index, e);

                E_field[X] = e[X] * dr;
                E_field[Y] = e[Y] * dr;
                E_field[Z] = e[Z] * dr;

                klein_add_double(&E_field_k[X], E_field[X]);
                klein_add_double(&E_field_k[Y], E_field[Y]);
                klein_add_double(&E_field_k[Z], E_field[Z]);
            }}}

            pc->Esub[X] = klein_sum(&E_field_k[X]);
            pc->Esub[Y] = klein_sum(&E_field_k[Y]);
            pc->Esub[Z] = klein_sum(&E_field_k[Z]);
        }
    }}}

    colloid_sums_halo(cinfo, COLLOID_SUM_ELECTRIC_FIELD);
    /* ... diagnostic output to fp omitted for clarity ... */
    return 0;
}
```

After `colloid_sums_halo(COLLOID_SUM_ELECTRIC_FIELD)` the Esub value is consistent
across all MPI ranks.

---

## 18. Particle force from Esub (`src/subgrid.c`)

### 18.1 `subgrid_update_forces_electrokinetics`

Takes the gathered `pc->Esub`, computes `force = kt·reunit·Esub·(q0−q1)`, scatters it
back onto fluid nodes using the same kernel (Newton III), and accumulates `pc->fex`.

```c
int subgrid_update_forces_electrokinetics(colloids_info_t* cinfo, map_t* map,
                                           physics_t* phys, psi_t* psi,
                                           hydro_t* hydro, subgrid_kernel_t kernel) {
    int i, j, k, ic, jc, kc, i_min, i_max, j_min, j_max, k_min, k_max;
    int ncell[3], index;
    int nlocal[3], offset[3];
    double kt, eunit, reunit, dr;
    double r[3], r0[3];
    double force[3] = { 0.0, 0.0, 0.0 };
    distributed_force_klein_t* force_k_indexed = NULL;
    colloid_t* pc;
    MPI_Comm comm;

    if (cinfo->nsubgrid == 0) return 0;

    cs_nlocal(cinfo->cs, nlocal);
    cs_nlocal_offset(cinfo->cs, offset);
    cs_cart_comm(cinfo->cs, &comm);
    colloids_info_ncell(cinfo, ncell);
    physics_kt(phys, &kt);
    psi_unit_charge(psi, &eunit);
    reunit = 1.0 / eunit;

    hydro_memcpy(hydro, tdpMemcpyDeviceToHost);

    int krange = subgrid_get_range(kernel);

    for (ic = 0; ic <= ncell[X] + 1; ic++) {
    for (jc = 0; jc <= ncell[Y] + 1; jc++) {
    for (kc = 0; kc <= ncell[Z] + 1; kc++) {
        colloids_info_cell_list_head(cinfo, ic, jc, kc, &pc);
        for (; pc; pc = pc->next) {
            if (pc->s.bc != COLLOID_BC_SUBGRID) continue;

            /* Electrostatic force on particle */
            force[X] = kt * reunit * pc->Esub[X] * (pc->s.q0 - pc->s.q1);
            force[Y] = kt * reunit * pc->Esub[Y] * (pc->s.q0 - pc->s.q1);
            force[Z] = kt * reunit * pc->Esub[Z] * (pc->s.q0 - pc->s.q1);

            r0[X] = pc->s.r[X] - 1.0 * offset[X];
            r0[Y] = pc->s.r[Y] - 1.0 * offset[Y];
            r0[Z] = pc->s.r[Z] - 1.0 * offset[Z];

            subgrid_get_lattice_index_range(r0, krange, nlocal,
                &i_min, &i_max, &j_min, &j_max, &k_min, &k_max);

            /* Precompute KB4 weights */
            double kb_w[4][4][4] = {{{0}}};
            double kb_w_sum = 1.0;
            if (kernel == SUBGRID_KERNEL_KB4) {
                kb_w_sum = 0.0;
                for (i = i_min; i <= i_max; i++)
                for (j = j_min; j <= j_max; j++)
                for (k = k_min; k <= k_max; k++) {
                    double wx = d_kb4(r0[X] - i);
                    double wy = d_kb4(r0[Y] - j);
                    double wz = d_kb4(r0[Z] - k);
                    kb_w[i-i_min][j-j_min][k-k_min] = wx * wy * wz;
                    kb_w_sum += wx * wy * wz;
                }
            }

            for (i = i_min; i <= i_max; i++) {
            for (j = j_min; j <= j_max; j++) {
            for (k = k_min; k <= k_max; k++) {
                double force_aux[3] = { 0.0, 0.0, 0.0 };
                index = cs_index(cinfo->cs, i, j, k);
                r[X] = r0[X] - (double)i;
                r[Y] = r0[Y] - (double)j;
                r[Z] = r0[Z] - (double)k;

                if      (kernel == SUBGRID_KERNEL_BSPLINE6)
                    dr = d_bspline6(r[X])*d_bspline6(r[Y])*d_bspline6(r[Z]);
                else if (kernel == SUBGRID_KERNEL_BSPLINE4)
                    dr = d_bspline4(r[X])*d_bspline4(r[Y])*d_bspline4(r[Z]);
                else if (kernel == SUBGRID_KERNEL_PESKIN4)
                    dr = d_peskin(r[X])*d_peskin(r[Y])*d_peskin(r[Z]);
                else if (kernel == SUBGRID_KERNEL_KB4)
                    dr = kb_w[i-i_min][j-j_min][k-k_min] / kb_w_sum;
                else if (kernel == SUBGRID_KERNEL_PESKIN6)
                    dr = d_peskin6(r[X])*d_peskin6(r[Y])*d_peskin6(r[Z]);

                force_aux[X] = force[X] * dr;
                force_aux[Y] = force[Y] * dr;
                force_aux[Z] = force[Z] * dr;
                add_force_to_array(&force_k_indexed, index, force_aux);
            }}}

            pc->fex[X] += force[X];
            pc->fex[Y] += force[Y];
            pc->fex[Z] += force[Z];
        }
    }}}

    colloid_sums_halo(cinfo, COLLOID_SUM_FORCE_EXT_ONLY);

    /* Apply accumulated fluid forces */
    if (force_k_indexed != NULL) {
        for (int i = 0; i < force_k_indexed->count; i++) {
            distributed_force_klein_entry_t* entry = force_k_indexed->entries[i];
            double nf[3];
            nf[X] = klein_sum(entry->force[X]);
            nf[Y] = klein_sum(entry->force[Y]);
            nf[Z] = klein_sum(entry->force[Z]);
            hydro_f_local_add(hydro, entry->cs_index, nf);
        }
        subgrid_free_distributed_force_t(&force_k_indexed);
    }

    hydro_memcpy(hydro, tdpMemcpyHostToDevice);
    return 0;
}
```

---

## 19. Fluid force from E-field gradient (`src/psi_force.c`)

### 19.1 `psi_force_gradmu_e_ewald`

For each fluid node, computes the smoothed electric field by gathering from kernel
neighbours, and applies the body force `F = ρ_elec · kt·reunit · E_smooth`.

```c
int psi_force_gradmu_e_ewald(psi_t* psi, fe_t* fe, hydro_t* hydro,
                              colloids_info_t* cinfo, subgrid_kernel_t kernel) {
    int ic, jc, kc, ia, nlocal[3], index;
    double rho_elec, e[3], kt, eunit, reunit, force[3];

    physics_t* phys = NULL;
    cs_nlocal(psi->cs, nlocal);
    physics_ref(&phys);
    physics_kt(phys, &kt);
    psi_unit_charge(psi, &eunit);
    reunit = 1.0 / eunit;

    int krange = subgrid_get_range(kernel);

    for (ic = 1; ic <= nlocal[X]; ic++) {
    for (jc = 1; jc <= nlocal[Y]; jc++) {
    for (kc = 1; kc <= nlocal[Z]; kc++) {
        int index0 = cs_index(psi->cs, ic, jc, kc);
        psi_rho_elec(psi, index0, &rho_elec);

        double r[3], r0[3] = { ic, jc, kc };
        int i, j, k, i_min, i_max, j_min, j_max, k_min, k_max;
        klein_t E_field_k[3];
        E_field_k[X] = klein_zero();
        E_field_k[Y] = klein_zero();
        E_field_k[Z] = klein_zero();

        subgrid_get_lattice_index_range_halo(r0, krange, nlocal,
            &i_min, &i_max, &j_min, &j_max, &k_min, &k_max);

        for (i = i_min; i <= i_max; i++) {
        for (j = j_min; j <= j_max; j++) {
        for (k = k_min; k <= k_max; k++) {
            index = cs_index(psi->cs, i, j, k);
            r[X] = r0[X] - (double)i;
            r[Y] = r0[Y] - (double)j;
            r[Z] = r0[Z] - (double)k;

            double dr;
            if      (kernel == SUBGRID_KERNEL_BSPLINE6)
                dr = d_bspline6(r[X])*d_bspline6(r[Y])*d_bspline6(r[Z]);
            else if (kernel == SUBGRID_KERNEL_BSPLINE4)
                dr = d_bspline4(r[X])*d_bspline4(r[Y])*d_bspline4(r[Z]);
            else if (kernel == SUBGRID_KERNEL_PESKIN4)
                dr = d_peskin(r[X])*d_peskin(r[Y])*d_peskin(r[Z]);
            else if (kernel == SUBGRID_KERNEL_KB4)
                dr = d_kb4(r[X])*d_kb4(r[Y])*d_kb4(r[Z]) / subgrid_kb4_norm_fluid();
            else if (kernel == SUBGRID_KERNEL_PESKIN6)
                dr = d_peskin6(r[X])*d_peskin6(r[Y])*d_peskin6(r[Z]);

            psi_electric_field(psi, index, e);
            klein_add_double(&E_field_k[X], e[X] * dr);
            klein_add_double(&E_field_k[Y], e[Y] * dr);
            klein_add_double(&E_field_k[Z], e[Z] * dr);
        }}}

        e[X] = klein_sum(&E_field_k[X]);
        e[Y] = klein_sum(&E_field_k[Y]);
        e[Z] = klein_sum(&E_field_k[Z]);

        /* Store smoothed E-field */
        psi->efield->data[addr_rank1(psi->efield->nsites, psi->efield->nf, index0, X)] = kt * e[X];
        psi->efield->data[addr_rank1(psi->efield->nsites, psi->efield->nf, index0, Y)] = kt * e[Y];
        psi->efield->data[addr_rank1(psi->efield->nsites, psi->efield->nf, index0, Z)] = kt * e[Z];

        for (ia = 0; ia < 3; ia++) {
            e[ia] *= kt * reunit;
            force[ia] = rho_elec * e[ia];
        }
        if (hydro) hydro_f_local_add(hydro, index0, force);
    }}}
    return 0;
}
```

### 19.2 `psi_force_gradmu_e_ewald_offset`

Adjoint-of-scatter version for interlacing (mode 2).  Each fluid source node `(si,sj,sk)`
acts as if it sits at `r_src = {si + mesh_offset, sj, sk}` and gathers the smoothed
E-field from that virtual position.  The force is deposited on the real node `(si_dest, sj, sk)`.

```c
int psi_force_gradmu_e_ewald_offset(psi_t* psi, fe_t* fe, hydro_t* hydro,
                              colloids_info_t* cinfo, subgrid_kernel_t kernel,
                              double mesh_offset) {
    int si, sj, sk, ia, nlocal[3];
    double e[3], kt, eunit, reunit, force[3];
    physics_t* phys = NULL;

    cs_nlocal(psi->cs, nlocal);
    physics_ref(&phys);
    physics_kt(phys, &kt);
    psi_unit_charge(psi, &eunit);
    reunit = 1.0 / eunit;

    int krange = subgrid_get_range(kernel);

    /* Mirror border exclusion logic from scatter */
    double frac_x = mesh_offset - floor(mesh_offset);
    int ix_excl_lo = 0, ix_excl_hi = -1;
    int ix_halo_lo = 1, ix_halo_hi = 0;
    if (frac_x > 0.0) {
        ix_excl_lo = nlocal[X] - krange + 1; ix_excl_hi = nlocal[X];
        ix_halo_lo = 1 - krange;             ix_halo_hi = 0;
    } else if (frac_x < 0.0) {
        ix_excl_lo = 1;               ix_excl_hi = krange;
        ix_halo_lo = nlocal[X] + 1;  ix_halo_hi = nlocal[X] + krange;
    }

    int list_cap = nlocal[X] + 2 * krange + 4;
    int *i_list = (int*)malloc(list_cap*sizeof(int)); int i_list_n = 0;
    for (si = 1; si <= nlocal[X]; si++) {
        if (si < ix_excl_lo || si > ix_excl_hi) i_list[i_list_n++] = si;
    }
    for (si = ix_halo_lo; si <= ix_halo_hi; si++) i_list[i_list_n++] = si;

    for (int ii = 0; ii < i_list_n; ii++) {
        si = i_list[ii];
        int si_dest = si;
        if (si_dest < 1) si_dest += nlocal[X];
        else if (si_dest > nlocal[X]) si_dest -= nlocal[X];

    for (sj = 1; sj <= nlocal[Y]; sj++) {
    for (sk = 1; sk <= nlocal[Z]; sk++) {
        int index0 = cs_index(psi->cs, si_dest, sj, sk);

        double rho_elec;
        psi_rho_elec(psi, index0, &rho_elec);
        if (rho_elec == 0.0) continue;

        /* Shifted source position */
        double r_src[3] = { (double)si + mesh_offset, (double)sj, (double)sk };

        int i, j, k, i_min, i_max, j_min, j_max, k_min, k_max;
        klein_t E_field_k[3];
        E_field_k[X] = klein_zero();
        E_field_k[Y] = klein_zero();
        E_field_k[Z] = klein_zero();

        subgrid_get_lattice_index_range_halo(r_src, krange, nlocal,
            &i_min, &i_max, &j_min, &j_max, &k_min, &k_max);

        for (i = i_min; i <= i_max; i++) {
        for (j = j_min; j <= j_max; j++) {
        for (k = k_min; k <= k_max; k++) {
            /* Wrap halo to interior */
            int iw = i, jw = j, kw = k;
            if (iw < 1) iw += nlocal[X]; else if (iw > nlocal[X]) iw -= nlocal[X];
            if (jw < 1) jw += nlocal[Y]; else if (jw > nlocal[Y]) jw -= nlocal[Y];
            if (kw < 1) kw += nlocal[Z]; else if (kw > nlocal[Z]) kw -= nlocal[Z];
            int index = cs_index(psi->cs, iw, jw, kw);

            double rx = r_src[X] - (double)i;
            double ry = r_src[Y] - (double)j;
            double rz = r_src[Z] - (double)k;

            double dr;
            if      (kernel == SUBGRID_KERNEL_BSPLINE6)
                dr = d_bspline6(rx)*d_bspline6(ry)*d_bspline6(rz);
            else if (kernel == SUBGRID_KERNEL_BSPLINE4)
                dr = d_bspline4(rx)*d_bspline4(ry)*d_bspline4(rz);
            else if (kernel == SUBGRID_KERNEL_PESKIN4)
                dr = d_peskin(rx)*d_peskin(ry)*d_peskin(rz);
            else if (kernel == SUBGRID_KERNEL_KB4)
                dr = d_kb4(rx)*d_kb4(ry)*d_kb4(rz) / subgrid_kb4_norm_fluid();
            else if (kernel == SUBGRID_KERNEL_PESKIN6)
                dr = d_peskin6(rx)*d_peskin6(ry)*d_peskin6(rz);
            else { dr = 0.0; }
            if (dr == 0.0) continue;

            psi_electric_field(psi, index, e);
            klein_add_double(&E_field_k[X], e[X] * dr);
            klein_add_double(&E_field_k[Y], e[Y] * dr);
            klein_add_double(&E_field_k[Z], e[Z] * dr);
        }}}

        e[X] = klein_sum(&E_field_k[X]);
        e[Y] = klein_sum(&E_field_k[Y]);
        e[Z] = klein_sum(&E_field_k[Z]);

        /* Store smoothed E-field at source node */
        psi->efield->data[addr_rank1(psi->efield->nsites, psi->efield->nf, index0, X)] = kt * e[X];
        psi->efield->data[addr_rank1(psi->efield->nsites, psi->efield->nf, index0, Y)] = kt * e[Y];
        psi->efield->data[addr_rank1(psi->efield->nsites, psi->efield->nf, index0, Z)] = kt * e[Z];

        for (ia = 0; ia < 3; ia++) {
            e[ia] *= kt * reunit;
            force[ia] = rho_elec * e[ia];
        }
        if (hydro) hydro_f_local_add(hydro, index0, force);
    }}}
    free(i_list);
    return 0;
}
```

---

## 20. FFT Poisson solver (`src/psi_fft.c`)

All code in this section compiles only with `nvcc` (`#ifdef __NVCC__`).

### 20.1 Data structures (`src/psi_fft.h`)

```c
typedef enum {
  PSI_FFT_LAPLACIAN_ANALYTIC = 0,  /* continuous k² */
  PSI_FFT_LAPLACIAN_DISCRETE = 1   /* stencil eigenvalues */
} psi_fft_laplacian_t;

typedef enum {
  PSI_FFT_INFLUENCE_SIMPLE  = 0,  /* 1/(ε k²) */
  PSI_FFT_INFLUENCE_HOCKNEY = 1,  /* Hockney & Eastwood Eq. 8-22 */
  PSI_FFT_INFLUENCE_EWALD   = 2   /* exp(-k²/4α²)/(ε k² |P̂|²) */
} psi_fft_influence_t;

typedef enum {
  PSI_FFT_DECONV_NONE    = 0,
  PSI_FFT_DECONV_PESKIN  = 1,
  PSI_FFT_DECONV_BSPLINE4 = 2,
  PSI_FFT_DECONV_BSPLINE6 = 3,
  PSI_FFT_DECONV_KB4     = 4,
  PSI_FFT_DECONV_PESKIN6 = 5
} psi_fft_deconv_t;

typedef struct psi_solver_fft_s {
  psi_solver_t   super;           /* virtual table base */
  psi_t*         psi;
  int            is_initialised;
  int            nx, ny, nz;
  int            ntotal[3];
  double         lx, ly, lz;
  double         epsilon;
  double         beta;
  void*          plan_forward;    /* cufftHandle* */
  void*          plan_backward;   /* cufftHandle* */
  void*          rho_real_d;      /* cufftDoubleReal* on device */
  void*          rho_complex_d;   /* cufftDoubleComplex* on device */
  void*          psi_complex_d;   /* cufftDoubleComplex* on device */
  void*          psi_real_d;      /* cufftDoubleReal* on device */
  double*        laplacian_eigenval_d;  /* precomputed λ(k), DISCRETE only */
  double*        G_opt_d;         /* precomputed Ĝ(k), HOCKNEY or EWALD */
  psi_fft_laplacian_t laplacian_type;
  psi_fft_influence_t influence_type;
  double         shape_a;
  int            n_alias_max;
  double         ewald_alpha;
  double         ewald_rcut;
} psi_solver_fft_t;
```

### 20.2 CUDA kernel: `copy_rho_to_real_kernel`

Copies charge density from Ludwig's SOA field (with halos) to the contiguous FFT array.
`rho_elec = rho0 − rho1`, scaled by `beta`.

```c
__global__ void copy_rho_to_real_kernel(double* rho_real, const double* rho_field,
                                         double beta,
                                         int nx, int ny, int nz, int nhalo,
                                         int nk, int nsites) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= nx * ny * nz) return;

    int iz = idx % nz;
    int iy = (idx / nz) % ny;
    int ix = idx / (nz * ny);

    /* Ludwig SOA index: row-major with halos, Z fastest */
    int str_z = 1;
    int str_y = nz + 2 * nhalo;
    int str_x = str_y * (ny + 2 * nhalo);
    int ic = ix + 1, jc = iy + 1, kc = iz + 1;
    int ludwig_idx = str_x*(nhalo+ic-1) + str_y*(nhalo+jc-1) + str_z*(nhalo+kc-1);

    /* addr_rank1(nsites, nk, index, species) = nsites*species + index */
    double rho0 = rho_field[nsites * 0 + ludwig_idx];
    double rho1 = rho_field[nsites * 1 + ludwig_idx];
    rho_real[idx] = (rho0 - rho1) * beta;
}
```

### 20.3 CUDA kernel: `copy_psi_from_real_kernel`

Normalises by 1/N (cuFFT is unnormalized) and adds external field `psi_ext = −E0·r`.

```c
__global__ void copy_psi_from_real_kernel(double* psi_field, const double* psi_real,
                                           double e0x, double e0y, double e0z,
                                           int nx, int ny, int nz, int nhalo) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= nx * ny * nz) return;

    int iz = idx % nz, iy = (idx / nz) % ny, ix = idx / (nz * ny);

    int str_z = 1, str_y = nz + 2*nhalo, str_x = str_y*(ny + 2*nhalo);
    int ic = ix+1, jc = iy+1, kc = iz+1;
    int ludwig_idx = str_x*(nhalo+ic-1) + str_y*(nhalo+jc-1) + str_z*(nhalo+kc-1);

    double norm = 1.0 / (double)(nx * ny * nz);
    double psi_poisson = psi_real[idx] * norm;
    double psi_ext = -(e0x*(double)ic + e0y*(double)jc + e0z*(double)kc);
    psi_field[ludwig_idx] = psi_poisson + psi_ext;
}
```

### 20.4 CUDA kernel: `poisson_divide_kernel` (analytic k²)

```c
__global__ void poisson_divide_kernel(cufftDoubleComplex* psi_hat,
                                       const cufftDoubleComplex* rho_hat,
                                       double epsilon,
                                       double kx_factor, double ky_factor, double kz_factor,
                                       int nx, int ny, int nz_complex) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= nx * ny * nz_complex) return;

    int iz = idx % nz_complex;
    int iy = (idx / nz_complex) % ny;
    int ix = idx / (nz_complex * ny);

    double kx = (ix <= nx/2) ? (double)ix : (double)(ix - nx);
    double ky = (iy <= ny/2) ? (double)iy : (double)(iy - ny);
    double kz = (double)iz;

    kx *= kx_factor; ky *= ky_factor; kz *= kz_factor;
    double k_sq = kx*kx + ky*ky + kz*kz;

    if (k_sq < 1.0e-15) {
        psi_hat[idx].x = 0.0;
        psi_hat[idx].y = 0.0;
    } else {
        double factor = 1.0 / (epsilon * k_sq);
        psi_hat[idx].x = rho_hat[idx].x * factor;
        psi_hat[idx].y = rho_hat[idx].y * factor;
    }
}
```

### 20.5 CUDA kernel: `poisson_divide_discrete_kernel` (stencil eigenvalues)

```c
__global__ void poisson_divide_discrete_kernel(cufftDoubleComplex* psi_hat,
                                                const cufftDoubleComplex* rho_hat,
                                                const double* laplacian_eigenval,
                                                double epsilon,
                                                int nx, int ny, int nz_complex) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= nx * ny * nz_complex) return;

    double lambda = laplacian_eigenval[idx];
    if (lambda < 1.0e-15) {
        psi_hat[idx].x = 0.0;
        psi_hat[idx].y = 0.0;
    } else {
        double factor = 1.0 / (epsilon * lambda);
        psi_hat[idx].x = rho_hat[idx].x * factor;
        psi_hat[idx].y = rho_hat[idx].y * factor;
    }
}
```

### 20.6 CUDA kernel: `poisson_divide_hockney_kernel` (Ewald / Hockney influence)

```c
__global__ void poisson_divide_hockney_kernel(cufftDoubleComplex* psi_hat,
                                               const cufftDoubleComplex* rho_hat,
                                               const double* G_opt,
                                               int nx, int ny, int nz_complex) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= nx * ny * nz_complex) return;

    double G = G_opt[idx];
    psi_hat[idx].x = rho_hat[idx].x * G;
    psi_hat[idx].y = rho_hat[idx].y * G;
}
```

### 20.7 `psi_solver_fft_solve` — main solve routine

```c
int psi_solver_fft_solve(psi_solver_fft_t* solver, int ntimestep) {
    psi_t* psi = solver->psi;
    int nlocal[3], nhalo;
    cs_nlocal(psi->cs, nlocal);
    cs_nhalo(psi->cs, &nhalo);

    int nx = solver->nx, ny = solver->ny, nz = solver->nz;
    int n_real    = nx * ny * nz;
    int n_complex = nx * ny * (nz/2 + 1);

    /* 1. Ensure rho on device */
    field_memcpy(psi->rho, tdpMemcpyHostToDevice);

    /* Get device data pointers via offsetof trick */
    double* rho_data_d = NULL, *psi_data_d = NULL;
    size_t data_offset = offsetof(field_t, data);
    cudaMemcpy(&rho_data_d, (char*)(psi->rho->target) + data_offset,
               sizeof(double*), cudaMemcpyDeviceToHost);
    cudaMemcpy(&psi_data_d, (char*)(psi->psi->target) + data_offset,
               sizeof(double*), cudaMemcpyDeviceToHost);

    int nk, nsites = psi->nsites;
    psi_nk(psi, &nk);

    int threads = 256;
    int blocks = (n_real + threads - 1) / threads;

    /* 2. Copy rho -> device contiguous array */
    copy_rho_to_real_kernel<<<blocks, threads>>>(
        (double*)solver->rho_real_d, rho_data_d,
        solver->beta, nx, ny, nz, nhalo, nk, nsites);
    cudaDeviceSynchronize();

    /* 3. Forward FFT R2C */
    cufftExecD2Z(*((cufftHandle*)solver->plan_forward),
                 (cufftDoubleReal*)solver->rho_real_d,
                 (cufftDoubleComplex*)solver->rho_complex_d);

    /* 4. Divide by Laplacian (analytic, discrete, or influence function) */
    blocks = (n_complex + threads - 1) / threads;
    if ((solver->influence_type == PSI_FFT_INFLUENCE_HOCKNEY ||
         solver->influence_type == PSI_FFT_INFLUENCE_EWALD) && solver->G_opt_d) {
        poisson_divide_hockney_kernel<<<blocks, threads>>>(
            (cufftDoubleComplex*)solver->psi_complex_d,
            (cufftDoubleComplex*)solver->rho_complex_d,
            (double*)solver->G_opt_d, nx, ny, nz/2+1);
    } else if (solver->laplacian_type == PSI_FFT_LAPLACIAN_DISCRETE) {
        poisson_divide_discrete_kernel<<<blocks, threads>>>(
            (cufftDoubleComplex*)solver->psi_complex_d,
            (cufftDoubleComplex*)solver->rho_complex_d,
            (double*)solver->laplacian_eigenval_d,
            solver->epsilon, nx, ny, nz/2+1);
    } else {
        double kx_factor = 2.0*M_PI/solver->lx;
        double ky_factor = 2.0*M_PI/solver->ly;
        double kz_factor = 2.0*M_PI/solver->lz;
        poisson_divide_kernel<<<blocks, threads>>>(
            (cufftDoubleComplex*)solver->psi_complex_d,
            (cufftDoubleComplex*)solver->rho_complex_d,
            solver->epsilon, kx_factor, ky_factor, kz_factor,
            nx, ny, nz/2+1);
    }
    cudaDeviceSynchronize();

    /* 5. Inverse FFT C2R */
    cufftExecZ2D(*((cufftHandle*)solver->plan_backward),
                 (cufftDoubleComplex*)solver->psi_complex_d,
                 (cufftDoubleReal*)solver->psi_real_d);

    /* 6. Copy psi -> Ludwig field, normalise, add external field */
    blocks = (n_real + threads - 1) / threads;
    copy_psi_from_real_kernel<<<blocks, threads>>>(
        psi_data_d, (double*)solver->psi_real_d,
        psi->e0[X], psi->e0[Y], psi->e0[Z],
        nx, ny, nz, nhalo);
    cudaDeviceSynchronize();

    /* 7. Sync to host and update halos */
    field_memcpy(psi->psi, tdpMemcpyDeviceToHost);
    psi_halo_psi(psi);
    field_memcpy(psi->psi, tdpMemcpyHostToDevice);

    return 0;
}
```

### 20.8 `psi_solver_fft_set_influence_ewald`

Precomputes the Ewald influence function
`Ĝ(k) = exp(−k²/4α²) / (ε · k² · |P̂(k)|²)` for all k-vectors and uploads to GPU.

```c
int psi_solver_fft_set_influence_ewald(psi_solver_fft_t* solver,
                                        double alpha, double rcut,
                                        psi_fft_deconv_t deconv) {
    int nx = solver->nx, ny = solver->ny, nz = solver->nz;
    int nz_complex = nz/2 + 1;
    int n_complex  = nx * ny * nz_complex;
    double alpha_sq_4 = 4.0 * alpha * alpha;
    double kx_factor = 2.0*M_PI/solver->lx;
    double ky_factor = 2.0*M_PI/solver->ly;
    double kz_factor = 2.0*M_PI/solver->lz;

    double* G_ewald_h = (double*)malloc(n_complex * sizeof(double));

    for (int ix = 0; ix < nx; ix++) {
        int kx_idx = (ix <= nx/2) ? ix : ix - nx;
        double kx = kx_factor * kx_idx;
    for (int iy = 0; iy < ny; iy++) {
        int ky_idx = (iy <= ny/2) ? iy : iy - ny;
        double ky = ky_factor * ky_idx;
    for (int iz = 0; iz < nz_complex; iz++) {
        double kz = kz_factor * iz;
        int idx = ix*(ny*nz_complex) + iy*nz_complex + iz;

        if (ix == 0 && iy == 0 && iz == 0) { G_ewald_h[idx] = 0.0; continue; }

        double k_sq = kx*kx + ky*ky + kz*kz;

        /* Compute Fourier transform of the assignment kernel P̂(k) = Πᵢ Σ_m φ(m) cos(kᵢ·m) */
        double px = 0.0, py = 0.0, pz = 0.0;
        if (deconv == PSI_FFT_DECONV_BSPLINE4) {
            for (int m = -2; m <= 2; m++) {
                double w = d_bspline4((double)m);
                px += w * cos(kx * m);
                py += w * cos(ky * m);
                pz += w * cos(kz * m);
            }
        } else if (deconv == PSI_FFT_DECONV_BSPLINE6) {
            for (int m = -3; m <= 3; m++) {
                double w = d_bspline6((double)m);
                px += w * cos(kx * m); py += w * cos(ky * m); pz += w * cos(kz * m);
            }
        } else if (deconv == PSI_FFT_DECONV_PESKIN) {
            for (int m = -2; m <= 2; m++) {
                double w = d_peskin((double)m);
                px += w * cos(kx * m); py += w * cos(ky * m); pz += w * cos(kz * m);
            }
        } else if (deconv == PSI_FFT_DECONV_KB4) {
            for (int m = -2; m <= 2; m++) {
                double w = d_kb4((double)m);
                px += w * cos(kx * m); py += w * cos(ky * m); pz += w * cos(kz * m);
            }
        } else if (deconv == PSI_FFT_DECONV_PESKIN6) {
            for (int m = -3; m <= 3; m++) {
                double w = d_peskin6((double)m);
                px += w * cos(kx * m); py += w * cos(ky * m); pz += w * cos(kz * m);
            }
        } else {  /* PSI_FFT_DECONV_NONE */
            px = py = pz = 1.0;
        }

        double P_hat_sq = (px * py * pz) * (px * py * pz);
        if (P_hat_sq < 1.0e-10) { G_ewald_h[idx] = 0.0; continue; }

        G_ewald_h[idx] = exp(-k_sq / alpha_sq_4) / (solver->epsilon * k_sq * P_hat_sq);
    }}}

    if (solver->G_opt_d == NULL)
        cudaMalloc(&solver->G_opt_d, n_complex * sizeof(double));

    cudaMemcpy(solver->G_opt_d, G_ewald_h, n_complex * sizeof(double),
               cudaMemcpyHostToDevice);
    free(G_ewald_h);

    solver->influence_type = PSI_FFT_INFLUENCE_EWALD;
    solver->ewald_alpha = alpha;
    solver->ewald_rcut  = rcut;

    return 0;
}
```

### 20.9 `psi_solver_fft_create_internal` — allocation and cuFFT plans

```c
static int psi_solver_fft_create_internal(psi_t* psi, psi_fft_laplacian_t laplacian_type,
                                           psi_solver_fft_t** psolver) {
    /* Allocates psi_solver_fft_t, 4 CUDA arrays (rho_real, rho_complex,
       psi_complex, psi_real), creates cufftPlan3d forward (CUFFT_D2Z) and
       backward (CUFFT_Z2D) plans.
       If laplacian_type == PSI_FFT_LAPLACIAN_DISCRETE, precomputes
       lambda(k) = Σ_p wlaplacian[p] * cos(k · cv[p]) for all k and
       uploads to laplacian_eigenval_d. */

    /* n_real    = nx * ny * nz */
    /* n_complex = nx * ny * (nz/2 + 1) */
    /* cufftPlan3d(plan, nx, ny, nz, type) */
    /* ... full code in src/psi_fft.c lines 134-361 ... */
}
```

The full implementation is in [src/psi_fft.c:134](src/psi_fft.c#L134).

---

## 21. Summary of the adjoint constraint

For exact momentum conservation (Newton III), the gather operator must be the exact
adjoint of the scatter operator.  This means:

- Particle scatter: `ρ_node += q_p · φ(r_node − r_p)`
- Particle gather: `F_p    = Σ_node q_p · E(r_node) · φ(r_node − r_p)`
- Fluid scatter:   `ρ'_node = Σ_src ρ(r_src) · φ(r_node − r_src)` (source acts as "particle")
- Fluid gather:    `F(r_src) = ρ(r_src) · E_smooth(r_src)` where `E_smooth(r_src) = Σ_node E(r_node) · φ(r_node − r_src)`

The same kernel φ must be used for all four operations.  The `mesh_offset` shift in the
interlacing mode is applied consistently to both scatter and gather by using the exact
same `i_list` construction and `r_src = {si + mesh_offset, sj, sk}` formula in both
`subgrid_scatter_fluid_offset` and `psi_force_gradmu_e_ewald_offset`.
