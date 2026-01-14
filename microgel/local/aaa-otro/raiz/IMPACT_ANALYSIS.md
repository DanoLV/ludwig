# Análisis de Impacto del Cambio de Pesos de Stencil

## Resumen del Cambio Propuesto

Cambiar los factores de `wlaplacian` en:

| Archivo | Línea | Cambio | Factor |
|---------|-------|--------|--------|
| `src/stencil_d3q7.c` | 65 | **NO CAMBIAR** | -8.0 (correcto) |
| `src/stencil_d3q19.c` | 59 | `-36.0` → `-6.0` | **÷6** |
| `src/stencil_d3q27.c` | 63 | `-216.0` → `-6.0` | **÷36** |

**NO cambiar** `wgradients` (ya está correcto en +3.0 para D3Q19/D3Q27)

## Archivos que Usan Estos Pesos

Búsqueda en el código encontró **3 archivos** que usan `wlaplacian` y `wgradients`:

### 1. `src/psi_petsc.c` - Solver de Poisson

**Uso de `wlaplacian`**:

```c
// Línea 338 - Construcción de matriz del Laplaciano
v[p] = s->wlaplacian[p] * epsilon;

// Línea 1017 - Versión con epsilon variable
v[p] = s->wlaplacian[p] * epsilon0;
```

**Uso de `wgradients`**:

```c
// Líneas 1006-1008 - Cálculo de gradiente de epsilon
gradeps[X] += s->wgradients[p]*s->cv[p][X]*epsilon1;
gradeps[Y] += s->wgradients[p]*s->cv[p][Y]*epsilon1;
gradeps[Z] += s->wgradients[p]*s->cv[p][Z]*epsilon1;

// Líneas 1020-1022 - Términos adicionales en ecuación de Poisson generalizada
v[p] += s->wgradients[p]*s->cv[p][X]*gradeps[X];
v[p] += s->wgradients[p]*s->cv[p][Y]*gradeps[Y];
v[p] += s->wgradients[p]*s->cv[p][Z]*gradeps[Z];
```

**IMPACTO**:
- ✅ El cambio de `wlaplacian` **corrige** el operador Laplaciano
- ✅ El solver resolverá ∇²ψ = -ρ/ε correctamente
- ✅ Los términos con `wgradients` (grad epsilon) ya están correctos
- ⚠️ Soluciones ψ de D3Q19/D3Q27 serán **diferentes** (correctas ahora)

### 2. `src/psi_gradients.c` - Campo Eléctrico

**Uso de `wgradients`**:

```c
// Líneas 61-63 - Cálculo de E = -∇ψ
e[X] -= s->wgradients[p]*cx*psi0;
e[Y] -= s->wgradients[p]*cy*psi0;
e[Z] -= s->wgradients[p]*cz*psi0;
```

**IMPACTO**:
- ✅ **NO cambia** (wgradients no se modifica)
- ✅ Pero como ψ será correcto, E también será correcto
- 📊 Campos eléctricos de D3Q19/D3Q27 cambiarán de magnitud:
  - D3Q19: **×6** más grande
  - D3Q27: **×36** más grande

### 3. `src/psi_force.c` - Fuerza Electrostática

**Uso de `wgradients`**:

```c
// Línea 430 - Divergencia del tensor de stress
force[ia] -= s->wgradients[p] * pth[ia][ib] * s->cv[p][ib];
```

**IMPACTO**:
- ✅ **NO cambia** (wgradients no se modifica)
- ✅ Calcula ∇·P correctamente (P = tensor de stress)
- 📊 Fuerzas sobre el fluido cambiarán porque:
  - El tensor de stress P depende de E
  - E será diferente (correcto ahora)
  - Fuerzas escaladas por el mismo factor que E

## Operaciones Afectadas

### ✅ OPERACIONES QUE SE CORRIGEN

1. **Solver de Poisson** (`psi_petsc.c`)
   - Resuelve: ∇²ψ = -ρ/ε
   - **ANTES**: D3Q19 y D3Q27 sobre-escalaban ∇²ψ
   - **DESPUÉS**: Todos los stencils dan ∇²ψ consistente

2. **Campo Eléctrico** (`psi_gradients.c`)
   - Calcula: E = -∇ψ
   - **ANTES**: E sub-escalado porque ψ estaba mal
   - **DESPUÉS**: E correcto

3. **Fuerza Electrostática** (`psi_force.c`)
   - Calcula: F = ∇·P donde P depende de E
   - **ANTES**: F sub-escalada
   - **DESPUÉS**: F correcta

4. **Ecuación de Poisson Generalizada** (epsilon variable)
   - Resuelve: ∇·(ε∇ψ) = -ρ
   - **ANTES**: Inconsistente entre ∇² y ∇
   - **DESPUÉS**: Consistente

### ❌ OPERACIONES QUE NO SE AFECTAN

- **Lattice Boltzmann** (LB): No usa estos pesos
- **Advección-Difusión**: Independiente
- **Nernst-Planck**: Usa `psi_gradients.c` pero ya correcto
- **Boundary Conditions**: No afectadas

## Cambios en Resultados Numéricos

### D3Q7
- ✅ **SIN CAMBIOS** (ya estaba correcto)

### D3Q19
- ψ (potencial): **×6 más pequeño** (era 6× más grande antes)
- E (campo): **×6 más grande** (era 1/6 del correcto)
- F (fuerza): **×6 más grande** (era 1/6 del correcto)
- ρ_elec (densidad de carga): **NO cambia** (input)

### D3Q27
- ψ (potencial): **×36 más pequeño** (era 36× más grande antes)
- E (campo): **×36 más grande** (era 1/36 del correcto)
- F (fuerza): **×36 más grande** (era 1/36 del correcto)
- ρ_elec (densidad de carga): **NO cambia** (input)

## Verificación de Consistencia

Después del cambio, todos los stencils deberían dar resultados **consistentes**:

```
Para el mismo problema físico (misma ρ, mismo ε):

D3Q7:  ψ, E, F → valores correctos
D3Q19: ψ, E, F → valores consistentes con D3Q7
D3Q27: ψ, E, F → valores consistentes con D3Q7
```

La diferencia entre stencils será solo por **precisión numérica** (orden de truncamiento), NO por factores de escala.

## Impacto en Resultados Existentes

### ⚠️ ADVERTENCIA CRÍTICA

**Todos los resultados previos con D3Q19 o D3Q27 tienen errores sistemáticos**:

1. **Simulaciones con microgel**:
   - Campos eléctricos: 1/6 o 1/36 del correcto
   - Fuerzas electroforéticas: 1/6 o 1/36 del correcto
   - Movilidades: Incorrectas
   - Número de Dukhin: Incorrecto

2. **Publicaciones**:
   - Si reportaste valores absolutos de E o F: **INCORRECTOS**
   - Si reportaste relaciones/ratios: **Posiblemente correctos** (si todos los cálculos usaron el mismo stencil incorrecto)

3. **Comparaciones con teoría**:
   - Si comparaste con teoría analítica: **Descuadres explicados**

## Compatibilidad hacia Atrás

### Archivos de Datos
- **Formato**: Sin cambios
- **Valores**: ψ será diferente → **incompatibles**
- **Recomendación**: Re-ejecutar simulaciones

### Inputs
- **Archivos de configuración**: Sin cambios necesarios
- **Parámetros físicos**: Mantener los mismos

### Outputs
- **Formato**: Sin cambios
- **Valores**: Diferentes (correctos ahora)

## Pruebas Recomendadas

Después de aplicar el cambio:

### 1. Test de Consistencia entre Stencils
```bash
# Ejecutar el mismo caso con los 3 stencils
./ludwig input_d3q7
./ludwig input_d3q19
./ludwig input_d3q27

# Comparar campos eléctricos - deberían ser similares
python compare_stencils.py
```

### 2. Test de Ley de Gauss
```bash
# Verificar ∇·E = ρ/ε para cada stencil
python verify_gauss_law.py
```

### 3. Test de Casos Conocidos
```bash
# Comparar con soluciones analíticas
python verify_analytical.py
```

## Código de Prueba

He creado en `/home/bater/Sim/ludwig/`:

1. `verify_stencil_consistency.py` - Test matemático básico
2. `derive_consistent_weights.py` - Derivación de pesos correctos
3. `analyze_fd_theory.py` - Análisis de isotropía

Después del cambio, ejecutar:
```bash
python3 verify_stencil_consistency.py
```

Debería mostrar que todos los stencils dan ∇²ψ = -6.0

## Resumen Ejecutivo

### ✅ BENEFICIOS DEL CAMBIO

1. **Consistencia matemática**: ∇·∇ψ = ∇²ψ para todos los stencils
2. **Campos correctos**: E = -∇ψ da valores correctos
3. **Fuerzas correctas**: F = ∇·P con valores correctos
4. **Mezcla de stencils**: Puedes usar diferentes stencils en diferentes partes
5. **Validación**: Resultados comparables con teoría analítica

### ⚠️ COSTOS DEL CAMBIO

1. **Incompatibilidad**: Resultados previos tienen errores
2. **Re-ejecución**: Simulaciones pasadas necesitan repetirse
3. **Publicaciones**: Posibles correcciones necesarias

### 🎯 RECOMENDACIÓN

**APLICAR EL CAMBIO** porque:
- Los valores actuales son matemáticamente incorrectos
- El cambio es simple (2 líneas de código)
- Mejora fundamental la consistencia del código
- Permite usar D3Q19 y D3Q27 confiablemente

## Implementación

Ver archivos:
- `fix_stencil_weights.patch` - Patch para aplicar
- `SOLUTION_STENCIL_WEIGHTS.md` - Solución detallada
- `derive_consistent_weights.py` - Derivación matemática
