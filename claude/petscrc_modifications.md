# Modificaciones para usar la configuración de `.petscrc` en PETSc/Ludwig

Este documento describe todas las modificaciones realizadas en `src/psi_petsc.c` para que
PETSc respete las opciones definidas en `microgel/local/.petscrc`, incluyendo vectores y
matrices en GPU (CUDA), el precondicionador BoomerAMG (hypre) y la tolerancia del solver.

Archivo de referencia: [microgel/local/.petscrc](../microgel/local/.petscrc)
Archivo modificado:    [src/psi_petsc.c](../src/psi_petsc.c)

---

## 1. Contenido de `.petscrc` y qué activa cada opción

```ini
# Vectores y matrices en GPU
-vec_type cuda
-mat_type aijcusparse

# Solver Krylov: Gradiente Conjugado estándar
-ksp_type cg

# Precondicionador: BoomerAMG de hypre (requiere petsc-cuda-hypre)
-pc_type hypre
-pc_hypre_type boomeramg

# Tolerancia relativa del solver
-ksp_rtol 5.0e-17
```

PETSc lee este archivo automáticamente si está en el **directorio de trabajo** al ejecutar
`Ludwig.exe` (donde se llama `PetscInitialize`), o en `$HOME` para efecto global.  Las
opciones del `.petscrc` tienen la misma prioridad que argumentos de línea de comandos.

---

## 2. Por qué el código original no respetaba `.petscrc`

En la versión upstream de Ludwig, `psi_solver_petsc_initialise` creaba el DMDA y el KSP
con tipos fijos en código (`VECSTANDARD`, `MATMPIAIJ`) y **no llamaba a
`KSPSetFromOptions`** después de crear el contexto.  Esto hacía que:

1. El tipo de vector y matriz quedara fijo en CPU aunque `.petscrc` pidiera `cuda`/
   `aijcusparse`.
2. El tipo de solver KSP y precondicionador del `.petscrc` fueran ignorados.
3. Al intentar usar hypre con matrices en GPU se producía el error:
   `HYPRE_MEMORY_DEVICE expects a device vector`.
4. Al copiar la solución de vuelta al host con `DMDAVecGetArray` sobre un `VECCUDA` se
   producía un segfault porque el vector no estaba vinculado a la CPU.

---

## 3. Modificaciones en `psi_solver_petsc_initialise` (líneas 257–298)

### 3.1 Detección dinámica de `-vec_type cuda` y selección de tipos GPU

**Código original (upstream):**
```c
PetscCall(DMSetVecType(solver->block->da, VECSTANDARD));
PetscCall(DMSetMatType(solver->block->da, MATMPIAIJ));
```

**Código modificado:**
```c
/*CHANGE INIT - 20260326 Use GPU vec/mat types when -vec_type cuda is set via .petscrc.
 * DMSetVecType/DMSetMatType must match the requested type BEFORE DMSetUp and
 * DMCreateGlobalVector, otherwise -vec_type cuda in .petscrc is ignored and
 * hypre PC fails with "HYPRE_MEMORY_DEVICE expects a device vector".
 * Original (CPU only):
 *   PetscCall(DMSetVecType(solver->block->da, VECSTANDARD));
 *   PetscCall(DMSetMatType(solver->block->da, MATMPIAIJ)); */
{
  PetscBool use_cuda = PETSC_FALSE;
  char vtype[64];
  PetscOptionsGetString(NULL, NULL, "-vec_type", vtype, sizeof(vtype), &use_cuda);
  if (use_cuda && !strcmp(vtype, "cuda")) {
    PetscCall(DMSetVecType(solver->block->da, VECCUDA));
    PetscCall(DMSetMatType(solver->block->da, MATAIJCUSPARSE));
  } else {
    PetscCall(DMSetVecType(solver->block->da, VECSTANDARD));
    PetscCall(DMSetMatType(solver->block->da, MATMPIAIJ));
  }
}
/*CHANGE END - 20260326 */
PetscCall(DMSetUp(solver->block->da));
```

**Por qué debe ir antes de `DMSetUp`:**
`DMSetUp` y `DMCreateGlobalVector` fijan el tipo de los vectores internos del DMDA.  Si se
llama a `DMSetUp` antes de `DMSetVecType`, los vectores ya quedan creados como
`VECSTANDARD` y el tipo no puede cambiarse después.  PETSc ignoraría silenciosamente la
opción `-vec_type cuda` del `.petscrc`.

**Por qué se lee la opción manualmente en lugar de dejarlo a PETSc:**
`DMSetFromOptions` aplica opciones generales de DM pero no propaga `-vec_type` al DMDA
cuando el `DM` ya tiene un tipo de vector por defecto.  La lectura explícita con
`PetscOptionsGetString` garantiza que el tipo coincida antes de que el DM sea configurado.

### 3.2 Activación de opciones KSP desde `.petscrc`

**Código original (upstream):**
```c
KSPCreate(PETSC_COMM_WORLD, &solver->block->ksp);
KSPSetOperators(solver->block->ksp, solver->block->a, solver->block->a);
KSPSetTolerances(solver->block->ksp, rtol, abstol, PETSC_DEFAULT, maxits);
/* KSPSetFromOptions ausente → .petscrc ignorado */
```

**Código modificado:**
```c
KSPCreate(PETSC_COMM_WORLD, &solver->block->ksp);
KSPSetOperators(solver->block->ksp, solver->block->a, solver->block->a);
KSPSetTolerances(solver->block->ksp, rtol, abstol, PETSC_DEFAULT, maxits);
/*CHANGE INIT - 20260326 Allow .petscrc / command-line options for KSP solver */
KSPSetFromOptions(solver->block->ksp);
/*CHANGE END - 20260326 */
```

`KSPSetFromOptions` aplica todas las opciones PETSc registradas (de `.petscrc`, variables
de entorno o línea de comandos) al contexto KSP **después** de las llamadas programáticas.
Las opciones del `.petscrc` **sobreescriben** los valores fijados por `KSPSetTolerances`,
por eso `-ksp_rtol 5.0e-17` en el `.petscrc` tiene efecto aunque el archivo de entrada de
Ludwig fije un valor distinto en `electrokinetics_solver_reltol`.

**Orden de precedencia resultante:**
```
KSPSetTolerances (valores del input de Ludwig)
    ↓  sobreescrito por
KSPSetFromOptions (valores de .petscrc / línea de comandos)
```

---

## 4. Modificación en `psi_solver_petsc_solve` — proyección del espacio nulo (líneas 405–416)

**Problema:** Con matrices `aijcusparse` en GPU, `MatSetNullSpace` no proyecta
automáticamente el modo constante del vector RHS y de la solución inicial.  El solver CG
recibe un sistema indefinido (contaminado por el modo nulo) y diverge con el error
`DIVERGED_INDEFINITE_MAT`.

**Código añadido:**
```c
/*CHANGE INIT - 20260326 Fix DIVERGED_INDEFINITE_MAT with aijcusparse:
 * MatSetNullSpace does not project the null space from the RHS/solution when
 * using GPU matrices. Explicit projection via MatNullSpaceRemove on b and x
 * ensures CG sees a consistent SPD system (no constant-mode contamination). */
{
  MatNullSpace nullsp;
  MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, NULL, &nullsp);
  MatNullSpaceRemove(nullsp, solver->block->b);
  MatNullSpaceRemove(nullsp, solver->block->x);
  MatNullSpaceDestroy(&nullsp);
}
/*CHANGE END - 20260326 */
```

Esta proyección explícita elimina la componente de valor medio del RHS y de la solución
inicial antes de cada `KSPSolve`, garantizando que el sistema sea simétrico definido
positivo (SPD) tal como requiere CG.  En CPU con `MATMPIAIJ` esta proyección la hacía
PETSc internamente al detectar el nullspace registrado; con `aijcusparse` esa detección
falla.

---

## 5. Modificación en `psi_solver_petsc_da_to_psi` — lectura en CPU del vector GPU (líneas 850–875)

**Problema:** `DMDAVecGetArray` es una operación sobre memoria del host (CPU).  Cuando el
vector `x` es de tipo `VECCUDA` su memoria vive en la GPU.  Llamar a `DMDAVecGetArray`
directamente produce un segfault o datos inválidos.

**Código original (upstream):**
```c
DMDAGetCorners(solver->block->da, &xs, &ys, &zs, &xw, &yw, &zw);
DMDAVecGetArray(solver->block->da, solver->block->x, &psi_3d);
/* ... copia de datos ... */
DMDAVecRestoreArray(solver->block->da, solver->block->x, &psi_3d);
```

**Código modificado:**
```c
DMDAGetCorners(solver->block->da, &xs, &ys, &zs, &xw, &yw, &zw);

/*CHANGE INIT - 20260326 If the solution vector is on GPU (VECCUDA), bind it to CPU
 * so that DMDAVecGetArray can access it from the host. */
VecBindToCPU(solver->block->x, PETSC_TRUE);
/*CHANGE END - 20260326 */

DMDAVecGetArray(solver->block->da, solver->block->x, &psi_3d);

/* ... copia de datos ... */

DMDAVecRestoreArray(solver->block->da, solver->block->x, &psi_3d);

/*CHANGE INIT - 20260326 Restore GPU binding after host read. */
VecBindToCPU(solver->block->x, PETSC_FALSE);
/*CHANGE END - 20260326 */
```

`VecBindToCPU(v, PETSC_TRUE)` fuerza a PETSc a mantener los datos en memoria del host y
sincroniza el contenido desde la GPU si es necesario.  `VecBindToCPU(v, PETSC_FALSE)`
restaura el comportamiento GPU para el siguiente `KSPSolve`.

---

## 6. Compatibilidad CUDA/C++ en asignación de structs (líneas 346–353, 962–969)

**Problema:** NVCC (compilador CUDA) en modo C++ no acepta *designated initializers* de C99
en estructuras:

```c
/* C99 válido, NVCC rechaza con -x cu */
MatStencil row = {.i = i, .j = j, .k = k};
```

**Código corregido:**
```c
/*CHANGE INIT - 20251119 CUDA C++ compatibility fix */
MatStencil row;
row.i = i;
row.j = j;
row.k = k;
/*CHANGE END*/
```

Esta corrección aplica en dos lugares:
- `psi_solver_petsc_matrix_set` (Laplaciano uniforme)
- `psi_solver_petsc_var_epsilon_matrix_set` (Laplaciano con ε variable)

---

## 7. Corrección de NULL pointer en `psi_solver_petsc_solve` (líneas 395–400)

**Problema:** La llamada original `assert(solver)` causa `SIGABRT` en tests que prueban el
manejo de puntero nulo.

**Código modificado:**
```c
/*CHANGE INIT - 20250112 Fix NULL pointer handling for test compatibility */
/* Original: assert(solver); */
if (solver == NULL) return -1;
/*CHANGE END - 20250112 */
```

---

## 8. Resumen de todas las modificaciones

| Función | Líneas aprox. | Cambio | Fecha |
|---------|---------------|--------|-------|
| `psi_solver_petsc_initialise` | 257–276 | Detección de `-vec_type cuda`, tipos GPU en DMDA antes de `DMSetUp` | 20260326 |
| `psi_solver_petsc_initialise` | 296–298 | Añadir `KSPSetFromOptions` para activar opciones de `.petscrc` | 20260326 |
| `psi_solver_petsc_solve` | 395–400 | NULL check en lugar de `assert` | 20250112 |
| `psi_solver_petsc_solve` | 405–416 | Proyección explícita del espacio nulo para GPU matrices | 20260326 |
| `psi_solver_petsc_matrix_set` | 346–353 | Compatibilidad CUDA/C++ en `MatStencil` | 20251119 |
| `psi_solver_petsc_da_to_psi` | 850–875 | `VecBindToCPU` antes/después de `DMDAVecGetArray` | 20260326 |
| `psi_solver_petsc_var_epsilon_matrix_set` | 962–969 | Compatibilidad CUDA/C++ en `MatStencil` | 20251119 |

---

## 9. Secuencia de inicialización completa con las modificaciones

```
PetscInitialize()          ← lee .petscrc del directorio de ejecución
    │
    ├─ DMDACreate3d(...)
    ├─ PetscOptionsGetString("-vec_type", ...)
    │   ├─ si "cuda":  DMSetVecType(VECCUDA)  + DMSetMatType(MATAIJCUSPARSE)
    │   └─ si no:      DMSetVecType(VECSTANDARD) + DMSetMatType(MATMPIAIJ)
    ├─ DMSetUp(da)            ← fija tipos de vectores
    ├─ DMCreateMatrix(da, &a) ← matriz en GPU si MATAIJCUSPARSE
    ├─ DMCreateGlobalVector(da, &x)  ← vector en GPU si VECCUDA
    ├─ VecDuplicate(x, &b)
    │
    ├─ KSPCreate(...)
    ├─ KSPSetOperators(ksp, a, a)
    ├─ KSPSetTolerances(ksp, rtol_input, ...)   ← valores del input de Ludwig
    ├─ KSPSetFromOptions(ksp)                   ← sobreescribe con .petscrc:
    │       -ksp_type cg
    │       -pc_type hypre / -pc_hypre_type boomeramg
    │       -ksp_rtol 5.0e-17
    └─ KSPSetUp(ksp)
```

```
Cada timestep (psi_solver_petsc_solve):
    │
    ├─ psi_solver_petsc_rhs_set(solver)     ← construye b en CPU
    ├─ psi_solver_petsc_psi_to_da(solver)   ← copia psi → x (guess inicial)
    ├─ MatNullSpaceRemove(nullsp, b)         ← proyecta modo nulo de b
    ├─ MatNullSpaceRemove(nullsp, x)         ← proyecta modo nulo de x
    ├─ KSPSetInitialGuessNonzero(PETSC_TRUE)
    ├─ KSPSolve(ksp, b, x)                  ← solve en GPU (CG + BoomerAMG)
    └─ psi_solver_petsc_da_to_psi(solver):
            VecBindToCPU(x, TRUE)           ← sincroniza GPU→CPU
            DMDAVecGetArray(da, x, &psi_3d) ← lectura en host
            ... copia x → psi->psi->data ...
            DMDAVecRestoreArray(...)
            VecBindToCPU(x, FALSE)          ← restaura modo GPU
```

---

## 10. Dependencias de compilación

Para que estas modificaciones compilen y funcionen se requiere:

| Requisito | Verificación |
|-----------|-------------|
| PETSc compilado con `--with-cuda` | `PETSC_INC = -I/usr/local/petsc-cuda-hypre/include` en `config.mk` |
| PETSc compilado con `--download-hypre` | `pc_type hypre` en `.petscrc` requiere hypre con soporte CUDA |
| CUDA disponible (nvcc) | `TARGET = nvcc` en `config.mk` |
| Flag `-DPETSC` activo | `HAVE_PETSC = true` en `config.mk` activa el bloque `#else` de `psi_petsc.c` |
| `.petscrc` en el directorio de ejecución | Copiar desde `microgel/local/.petscrc` al directorio donde se ejecuta `Ludwig.exe` |

La instalación activa en este sistema es `/usr/local/petsc-cuda-hypre` (PETSc con CUDA +
hypre con soporte GPU).  La instalación alternativa sin hypre es
`/usr/local/petsc-cuda` — con ella sólo funcionan los precondicionadores `jacobi` y
`none`.

---

## 11. Opciones comentadas en `.petscrc` y cuándo activarlas

| Opción | Cuándo activar |
|--------|----------------|
| `-ksp_type pipecg` | Alternativa a `cg` con mejor estabilidad numérica en GPU (menos pérdida de ortogonalidad); probar si CG diverge en sistemas grandes |
| `-pc_type jacobi` | Reemplaza hypre cuando PETSc se compila sin `--download-hypre` |
| `-ksp_monitor` | Diagnóstico: imprime residual en cada iteración; desactivar en producción |
| `-ksp_converged_reason` | Imprime la razón de convergencia al final de cada solve |
| `-ksp_view` | Imprime la configuración completa del solver al inicio |
| `-ksp_atol 1.0e-15` | Tolerancia absoluta; útil cuando el RHS es muy pequeño |
| `-ksp_max_it 10000` | Límite de iteraciones; el default de PETSc es 10000 |
| `-ksp_gmres_restart 100` | Sólo relevante si se cambia a `-ksp_type gmres` |

**Nota sobre `-ksp_rtol 5.0e-17`:** El input de Ludwig establece la tolerancia via
`electrokinetics_solver_reltol`.  Un valor de `1e-28` (típico en tests) es inalcanzable
en float64 y hace que CG pierda ortogonalidad.  El valor `5.0e-17` en `.petscrc`
sobreescribe ese valor a uno que CG+BoomerAMG puede alcanzar sin degradarse (~10 veces el
épsilon de máquina).
