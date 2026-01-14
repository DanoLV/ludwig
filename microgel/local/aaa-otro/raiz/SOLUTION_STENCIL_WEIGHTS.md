# Solución al Problema de Inconsistencia de Stencils

## Resumen Ejecutivo

**PROBLEMA**: Los stencils D3Q7, D3Q19 y D3Q27 producen campos eléctricos muy diferentes para la misma configuración física debido a normalizaciones inconsistentes entre el operador Laplaciano y el operador gradiente.

**CAUSA RAÍZ**: Los factores de escala actuales (-8, -36, -216) para `wlaplacian` están **sobre-escalando** el operador Laplaciano en D3Q19 y D3Q27, lo que resulta en potenciales ψ sub-escalados, y por ende campos eléctricos E incorrectos.

**SOLUCIÓN**: Normalizar todos los stencils para que usen el **mismo esquema de diferencias finitas**, ya sea:
- Opción A: Todos usan FD de 2º orden (solo vecinos cardinales)
- Opción B: Todos usan FD isotrópico de alto orden (incluyendo diagonales)

## Análisis Detallado

### 1. Verificación Matemática

Test con campo ψ(x,y,z) = x² + y² + z²  (∇²ψ debería dar -6.0):

| Stencil | ∇²ψ Actual | Factor Actual | ∇²ψ Esperado | Error |
|---------|------------|---------------|--------------|-------|
| D3Q7    | -6.0       | -8.0          | -6.0         | 0%    |
| D3Q19   | -36.0      | -36.0         | -6.0         | 500%  |
| D3Q27   | -216.0     | -216.0        | -6.0         | 3500% |

**El gradiente funciona correctamente en TODOS los casos**, pero actúa sobre un ψ incorrecto.

### 2. Análisis de Isotropía

Los pesos LB **SÍ son isotrópicos** (todos los puntos a la misma distancia tienen el mismo peso), PERO:

**D3Q7** (diferencias finitas de 2º orden estándar):
- Usa solo vecinos cardinales (d² = 1)
- wlaplacian = [-6.0, +1.0 × 6]
- wgradients = [0.0, +0.5 × 6]
- ✓ CORRECTO

**D3Q19** (diferencias finitas de alto orden):
- Usa cardinales (d² = 1) **Y** diagonales de cara (d² = 2)
- wlaplacian = [+24, -2.0 × 6, -1.0 × 12]
- wgradients = [0, 0.167 × 6, 0.083 × 12]
- ✗ SOBRE-ESCALADO

**D3Q27** (diferencias finitas de muy alto orden):
- Usa cardinales (d² = 1) + diagonales de cara (d² = 2) + diagonales de cubo (d² = 3)
- wlaplacian = [+152, -16.0 × 6, -4.0 × 12, -1.0 × 8]
- wgradients = [0, 0.222 × 6, 0.056 × 12, 0.014 × 8]
- ✗ SOBRE-ESCALADO

### 3. Origen del Problema

Los factores actuales:
```c
// D3Q7
s->wlaplacian[p] = -8.0  * wv[p];  // ✓ Correcto
s->wgradients[p] = +4.0  * wv[p];  // ✓ Correcto

// D3Q19
s->wlaplacian[p] = -36.0 * wv[p];  // ✗ Sobre-escalado 6×
s->wgradients[p] = +3.0  * wv[p];  // ✓ Correcto

// D3Q27
s->wlaplacian[p] = -216.0 * wv[p]; // ✗ Sobre-escalado 36×
s->wgradients[p] = +3.0  * wv[p];  // ✓ Correcto
```

La razón de estos factores viene de querer que `Σ w_lap * wv = -6` **usando todos los puntos**, pero esto crea una **inconsistencia dimensional**:

- El solver resuelve: `∇²ψ = -ρ/ε` usando `wlaplacian`
- Con D3Q27: `(-216 * wv) * ψ = -ρ/ε` → `ψ_resuelto = ψ_correcto / 36`
- El campo calcula: `E = -∇ψ` usando `wgradients` (correcto)
- Resultado: `E_calculado = E_correcto / 36`

## Soluciones Propuestas

### Opción A: Normalizar a FD de 2º Orden Estándar (RECOMENDADA)

**Ventaja**: Todos los stencils usan el mismo orden de aproximación → consistencia garantizada

**Cambios**:

```c
// stencil_d3q7.c - NO CAMBIAR (ya correcto)
s->wlaplacian[p] = -8.0*wv[p];
s->wgradients[p] = +4.0*wv[p];

// stencil_d3q19.c - CAMBIAR a usar solo vecinos cardinales
// Ignorar puntos diagonales (d² = 2)
for (int p = 0; p < s->npoints; p++) {
  int dist_sq = s->cv[p][X]*s->cv[p][X] +
                s->cv[p][Y]*s->cv[p][Y] +
                s->cv[p][Z]*s->cv[p][Z];

  if (dist_sq == 1) {  // Vecinos cardinales
    s->wlaplacian[p] = -8.0*wv[p];  // CAMBIO: -36.0 → -8.0
    s->wgradients[p] = +4.0*wv[p];  // CAMBIO: +3.0 → +4.0
  } else {  // Puntos diagonales
    s->wlaplacian[p] = 0.0;
    s->wgradients[p] = 0.0;
  }
}

// stencil_d3q27.c - CAMBIAR a usar solo vecinos cardinales
for (int p = 0; p < s->npoints; p++) {
  int dist_sq = s->cv[p][X]*s->cv[p][X] +
                s->cv[p][Y]*s->cv[p][Y] +
                s->cv[p][Z]*s->cv[p][Z];

  if (dist_sq == 1) {  // Vecinos cardinales
    s->wlaplacian[p] = -8.0*wv[p];  // CAMBIO: -216.0 → -8.0
    s->wgradients[p] = +4.0*wv[p];  // CAMBIO: +3.0 → +4.0
  } else {  // Puntos diagonales
    s->wlaplacian[p] = 0.0;
    s->wgradients[p] = 0.0;
  }
}
```

**Resultado**: Todos los stencils darán ∇²ψ = -6.0 y E = -∇ψ consistente.

### Opción B: Normalizar el Laplaciano (MÁS SIMPLE, POSIBLEMENTE INCORRECTA)

**ADVERTENCIA**: Esta opción cambia solo los factores pero mantiene el uso de puntos diagonales. Podría romper otras propiedades matemáticas.

```c
// stencil_d3q19.c
s->wlaplacian[p] = -8.0*wv[p];  // CAMBIO: -36.0 → -8.0
s->wgradients[p] = +3.0*wv[p];  // MANTENER

// stencil_d3q27.c
s->wlaplacian[p] = -8.0*wv[p];  // CAMBIO: -216.0 → -8.0
s->wgradients[p] = +3.0*wv[p];  // MANTENER
```

**PROBLEMA**: No garantiza que ∇²ψ = -6.0 porque los puntos diagonales también contribuyen con el factor -8.0, creando un sobre-conteo.

### Opción C: Derivar Pesos Isotrópicos de Alto Orden (CORRECTA PERO COMPLEJA)

Calcular pesos que usen **todos los puntos** del stencil pero de forma consistente:

**Requiere**:
1. Resolver sistema de ecuaciones desde expansión de Taylor hasta 4º-6º orden
2. Garantizar isotropía y consistencia
3. Validar numéricamente

**Ventaja**: Usa todo el stencil → potencialmente mayor precisión

**Desventaja**: Matemáticamente complejo, difícil de validar

## Recomendación Final

**IMPLEMENTAR OPCIÓN A**:

1. **Modificar** `stencil_d3q19.c` y `stencil_d3q27.c` para usar solo vecinos cardinales
2. **Eliminar** el código comentado de corrección RHS en `psi_petsc.c` (líneas 467-520)
3. **Probar** con los casos de test existentes para verificar consistencia
4. **Documentar** el cambio claramente

### Código de Implementación

Ver archivo: `fix_stencil_weights.patch`

## Validación

Después de aplicar los cambios:

1. Ejecutar `verify_stencil_consistency.py` para verificar matemáticamente
2. Comparar campos eléctricos entre D3Q7, D3Q19 y D3Q27 con misma ρ
3. Verificar que ∇·E = ρ/ε se cumpla para todos los stencils
4. Comparar energía electrostática total

## Referencias

- Análisis completo: `STENCIL_CONSISTENCY_ANALYSIS.md`
- Test de consistencia: `verify_stencil_consistency.py`
- Análisis de isotropía: `analyze_fd_theory.py`
- Cálculo de pesos óptimos: `calculate_correct_weights.py`

## Archivos a Modificar

1. **src/stencil_d3q19.c** - líneas 55-65
2. **src/stencil_d3q27.c** - líneas 55-70
3. **src/psi_petsc.c** - eliminar líneas 467-586 (código comentado)

## Impacto

**Cambios en resultados**:
- D3Q7: Sin cambios
- D3Q19: Campos eléctricos ~6× más grandes (correctos ahora)
- D3Q27: Campos eléctricos ~36× más grandes (correctos ahora)

**Compatibilidad**:
- Resultados previos con D3Q19/D3Q27 necesitarán re-ejecutarse
- Archivos de datos existentes quedan obsoletos
- Publicaciones basadas en D3Q19/D3Q27 podrían tener valores incorrectos
