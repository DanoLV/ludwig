# Análisis de Consistencia de Stencils Electrocinéticos

## Problema Identificado

Los stencils D3Q7, D3Q19 y D3Q27 usan **factores de escala inconsistentes** entre el operador Laplaciano (`wlaplacian`) y el operador gradiente (`wgradients`), lo que resulta en campos eléctricos incorrectos cuando se mezclan diferentes stencils o tamaños de lattice LB.

## Verificación Matemática

### Test Ejecutado

Campo de prueba: `ψ(x,y,z) = x² + y² + z²`
- Laplaciano teórico: `∇²ψ = 6.0`
- Gradiente teórico en (1,1,1): `∇ψ = (2, 2, 2)`

### Resultados Obtenidos

| Stencil | ∇²ψ Calculado | Error Laplaciano | ∇ψ Calculado | Error Gradiente |
|---------|---------------|------------------|--------------|-----------------|
| D3Q7    | -6.0          | Signo incorrecto | (-2,-2,-2)   | ✓ Correcto      |
| D3Q19   | -36.0         | 6× sobre-escalado| (-2,-2,-2)   | ✓ Correcto      |
| D3Q27   | -216.0        | 36× sobre-escalado| (-2,-2,-2)  | ✓ Correcto      |

**Conclusión**: El gradiente es correcto en todos los casos, pero el Laplaciano tiene diferentes factores de escala.

## Análisis de Código Actual

### 1. Factores de Escala en Stencils

```c
// D3Q7 (stencil_d3q7.c:65-66)
s->wlaplacian[p] = -8.0*wv[p];   // wv[p] = {2/8, 1/8, ...}
s->wgradients[p] = +4.0*wv[p];   // Relación: 8/4 = 2:1

// D3Q19 (stencil_d3q19.c:59-60)
s->wlaplacian[p] = -36.0*wv[p];  // wv[p] = {12/36, 2/36, 1/36, ...}
s->wgradients[p] = +3.0*wv[p];   // Relación: 36/3 = 12:1

// D3Q27 (stencil_d3q27.c:63-64)
s->wlaplacian[p] = -216.0*wv[p]; // wv[p] = {64/216, 16/216, 4/216, 1/216, ...}
s->wgradients[p] = +3.0*wv[p];   // Relación: 216/3 = 72:1
```

### 2. Uso en el Solver de Poisson

```c
// psi_petsc.c:338 - Construcción de la matriz del Laplaciano
v[p] = s->wlaplacian[p] * epsilon;

// Esto resuelve: ∇²ψ = -ρ/ε
// PERO con el Laplaciano sobre-escalado, ψ queda sub-escalado
```

### 3. Uso en Cálculo de Campo Eléctrico

```c
// psi_gradients.c:61-63 - Cálculo del campo E = -∇ψ
e[X] -= s->wgradients[p]*cx*psi0;
e[Y] -= s->wgradients[p]*cy*psi0;
e[Z] -= s->wgradients[p]*cz*psi0;

// El gradiente está bien normalizado
// PERO actúa sobre un ψ que fue calculado con Laplaciano sobre-escalado
```

## Consecuencia del Problema

Para la ecuación de Poisson `∇²ψ = -ρ/ε`:

1. **D3Q7**: Laplaciano con factor -8 → ψ se resuelve con escala casi correcta
2. **D3Q19**: Laplaciano con factor -36 → ψ se resuelve 6× más pequeño
3. **D3Q27**: Laplaciano con factor -216 → ψ se resuelve 36× más pequeño

Cuando calculas `E = -∇ψ`:
- El operador gradiente está bien normalizado
- PERO actúa sobre un ψ que está sub-escalado
- Resultado: **E tiene la magnitud incorrecta**

Ejemplo numérico:
```
Para ρ/ε = 1.0:

D3Q7:  ∇²ψ ≈ -6ψ = -1.0  →  ψ ≈ 1/6    →  E = -∇ψ ≈ correcto
D3Q19: ∇²ψ ≈ -36ψ = -1.0 →  ψ ≈ 1/36   →  E = -∇ψ ≈ 1/6 del correcto
D3Q27: ∇²ψ ≈ -216ψ = -1.0 → ψ ≈ 1/216  →  E = -∇ψ ≈ 1/36 del correcto
```

## Intentos Previos de Corrección

### Evidencia en psi_petsc.c (líneas 467-486)

Ya identificaste este problema y comentaste:
```c
/*CHANGE INIT - 20251201 Correct scaling factor for external field RHS */
/*CHANGE UPDATED - 20251203 Removed RHS correction after Laplacian renormalization */
/* The external field contribution to RHS must be scaled for D3Q19 and D3Q27.
 * D3Q7 already works correctly in the original code (no correction needed).
 * For stencil D3Q19: correction = 36.0/6.0 = 6.0
 * For stencil D3Q27: correction = 216.0/6.0 = 36.0
 *
 * UPDATE 20251203: After renormalizing wlaplacian in stencil_d3q19.c and
 * stencil_d3q27.c from (-36, -216) to (-6, -6), the RHS no longer needs
 * this correction. All stencils now use the same normalization.
 */
```

**PERO**: Los archivos stencil_d3q19.c y stencil_d3q27.c **NO fueron modificados**
(los .old son idénticos a los actuales).

## Solución Propuesta

Hay **DOS enfoques posibles**:

### Opción A: Normalizar el Laplaciano (MÁS SIMPLE Y CORRECTA)

Cambiar los factores de `wlaplacian` para que todos los stencils den ∇²ψ = 6.0:

```c
// D3Q7: NO CAMBIAR (ya está casi correcto, solo signo)
s->wlaplacian[p] = -8.0*wv[p];   // Da -6.0 (correcto en magnitud)
s->wgradients[p] = +4.0*wv[p];   // Da 0.5 (correcto)

// D3Q19: CAMBIAR de -36.0 a -8.0
s->wlaplacian[p] = -8.0*wv[p];   // Nuevo: dará ≈-6.0
s->wgradients[p] = +3.0*wv[p];   // Mantener

// D3Q27: CAMBIAR de -216.0 a -8.0
s->wlaplacian[p] = -8.0*wv[p];   // Nuevo: dará ≈-6.0
s->wgradients[p] = +3.0*wv[p];   // Mantener
```

**PROBLEMA CON ESTA OPCIÓN**:
- Necesito verificar si -8.0 es el factor correcto para D3Q19 y D3Q27
- Los pesos LB son diferentes para cada stencil
- Podría no dar exactamente -6.0

### Opción B: Derivar Factores Correctos desde Teoría LB

Los factores actuales vienen de la teoría de Lattice Boltzmann. Necesito:

1. Entender **por qué** se usan estos factores específicos (-36, -216)
2. Verificar si hay una **normalización implícita** en la teoría LB
3. Determinar los factores correctos que garanticen consistencia

## Siguiente Paso Recomendado

Crear un **test numérico completo** que:

1. Construya la matriz del Laplaciano para cada stencil
2. Resuelva ∇²ψ = -ρ/ε con una ρ conocida
3. Calcule E = -∇ψ con el gradiente
4. Verifique que ∇·E = ρ/ε (consistencia de Poisson)

Este test nos dirá **exactamente** qué factores usar para cada stencil.

## Archivos a Modificar

1. `src/stencil_d3q19.c` - línea 59
2. `src/stencil_d3q27.c` - línea 63
3. Posiblemente `src/stencil_d3q7.c` - línea 65 (verificar signo)

## Referencias en el Código

- Solver Poisson: [psi_petsc.c:338](src/psi_petsc.c#L338)
- Campo eléctrico: [psi_gradients.c:61-63](src/psi_gradients.c#L61-L63)
- Comentario previo: [psi_petsc.c:467-486](src/psi_petsc.c#L467-L486)
- Pesos D3Q7: [stencil_d3q7.c:65-66](src/stencil_d3q7.c#L65-L66)
- Pesos D3Q19: [stencil_d3q19.c:59-60](src/stencil_d3q19.c#L59-L60)
- Pesos D3Q27: [stencil_d3q27.c:63-64](src/stencil_d3q27.c#L63-L64)
