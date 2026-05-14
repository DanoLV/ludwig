# Análisis de Corrección PN para FFT con Asignación Asimétrica de Carga

## 1. Situación del Sistema

En Ludwig, el sistema electrocinético tiene una **asimetría fundamental** en cómo se asignan las cargas al mesh:

### 1.1 Cargas Iónicas (en nodos)

Las densidades de carga iónica ρ₀, ρ₁ están definidas **directamente en los nodos** del mesh de Lattice Boltzmann. No hay función de asignación intermedia - las cargas están localizadas exactamente en los puntos de la red discreta.

```
q_ion(r_n) = (ρ₀(r_n) - ρ₁(r_n)) · ΔV
```

donde ΔV = 1 (volumen de celda unitario en unidades de red).

En términos de función de asignación, esto corresponde a una **delta de Kronecker** en el mesh discreto:
```
W_ions(k) = 1   (función de asignación = delta de Kronecker en el nodo)
```

**Nota**: Esto significa que el solver FFT está implícitamente optimizado para cargas "tipo ion" localizadas en nodos. Cualquier partícula subgrid con posición arbitraria está desajustada respecto a esta función de Green.

### 1.2 Partícula Subgrid

La carga de la partícula se distribuye a los nodos vecinos usando la **delta regularizada de Peskin**:

```
q_mesh(r_n) = q_p · d_peskin(r_p - r_n)
```

donde:
- `q_p = q0 - q1` es la carga neta de la partícula
- `r_p` es la posición de la partícula (generalmente NO coincide con un nodo)
- `d_peskin(r)` es la función delta de Peskin 3D

La función de Peskin tiene soporte compacto [-2, 2] en cada dimensión:
```
d_peskin(r) = d(rx) · d(ry) · d(rz)

donde d(x) = {
  (1/8)(3 - 2|x| + √(1 + 4|x| - 4x²))     si |x| ≤ 1
  (1/8)(5 - 2|x| - √(-7 + 12|x| - 4x²))   si 1 < |x| ≤ 2
  0                                         si |x| > 2
}
```

En el espacio de Fourier:
```
W_particle(k) = Ŵ_peskin(k) ≠ 1
```

---

## 2. El Solver FFT Actual

### 2.1 Ecuación que Resuelve

El solver FFT resuelve la ecuación de Poisson:
```
∇²ψ = -ρ_total / ε
```

donde `ρ_total` incluye tanto cargas iónicas como cargas de partículas distribuidas.

### 2.2 Función de Green del FFT

En el espacio de Fourier, el solver actual usa:
```
ψ̂(k) = ρ̂(k) / (ε · λ(k))
```

donde λ(k) es el eigenvalor del operador Laplaciano:
- **Analítico**: λ(k) = k² = kx² + ky² + kz²
- **Discreto**: λ(k) = Σ_p w_p · cos(k·c_p) (suma sobre stencil)

Llamemos a esta función de Green del solver:
```
G_solver(k) = 1 / (ε · λ(k))
```

La función de Green en espacio real:
```
G_solver(r) = (1/N³) · Σ_k [ cos(k·r) / (ε · λ(k)) ]
```

### 2.3 Función de Green Óptima según Hockney & Eastwood

Según Hockney & Eastwood (Ec. 8-26), la función de Green óptima para minimizar errores debería ser:
```
G'_opt(k) = Σ_m [ D*(k+k_m) · R̂(k+k_m) ] / |Σ_m W(k+k_m)|²
```

donde:
- **W(k)**: Función de asignación de carga (Peskin en nuestro caso para partículas)
- **D(k)**: Operador de diferenciación discreta para obtener el campo desde el potencial
- **R̂(k)**: "Shaping function" del potencial de referencia (para Coulomb: R̂(k) = 4π/k² en unidades Gaussianas)
- **k_m**: Vectores de aliasing del mesh (k_m = 2πm/L para enteros m)

**El solver actual usa G_solver(k) = 1/(ε·λ(k)), que:**
1. **NO considera W(k)** - la función de asignación de Peskin
2. **NO considera aliasing** - los términos k_m ≠ 0
3. **Asume implícitamente cargas perfectamente localizadas en nodos** (W = 1)

---

## 3. Análisis de Interacciones

### 3.1 Interacción Ion-Ion (Nodo-Nodo)

**Situación**: Dos cargas iónicas en nodos r_m y r_n.

**Lo que calcula el FFT**:
```
ψ_solver(r_n) = q_ion(r_m) · G_solver(r_n - r_m)
```

**Lo que debería ser (Coulomb)**:
```
ψ_Coulomb(r_n) = q_ion(r_m) · G_Coulomb(r_n - r_m)
     = q_ion(r_m) / (4πε|r_n - r_m|)
```

**Error**:
```
Δψ_ion-ion = q_ion · [G_Coulomb(r) - G_solver(r)]
           = q_ion · G_diff(r)
```

**Nota importante**: Este error **depende de la separación r** entre los nodos. Lo que es sistemático es que cualquier par de iones a la misma separación y en la misma dirección de la red experimenta el mismo error. Como todos los iones usan el mismo esquema discreto, este error afecta principalmente a propiedades colectivas (screening efectivo, coeficientes de difusión), pero **no introduce la asimetría partícula-ion** que es el foco de la corrección PN.

### 3.2 Interacción Partícula-Ion (Partícula-Nodo)

Esta es la interacción **más importante** para la fuerza sobre la partícula.

**Situación**: Partícula con carga q_p en posición r_p, ion con carga q_ion en nodo r_n.

#### 3.2.1 Potencial en el nodo debido a la partícula

**Lo que calcula el FFT**:

La carga de la partícula se distribuye primero con Peskin:
```
ψ_solver(r_n) = Σ_j [ q_p · w_j · G_solver(r_n - r_j) ]
```
donde:
- La suma es sobre nodos j en el soporte de Peskin (típicamente 4³ = 64 nodos)
- w_j = d_peskin(r_p - r_j) es el peso de Peskin

**Lo que debería ser (Coulomb de carga puntual)**:
```
ψ_Coulomb(r_n) = q_p · G_Coulomb(r_n - r_p)
               = q_p / (4πε|r_n - r_p|)
```

**Corrección del potencial** (para Nernst-Planck):
```
Δψ(r_n) = q_p · G_Coulomb(r_n - r_p) - Σ_j [ q_p · w_j · G_solver(r_n - r_j) ]

        = q_p · [ G_Coulomb(r_n - r_p) - Σ_j w_j · G_solver(r_n - r_j) ]
```

#### 3.2.2 Campo en la partícula debido a los iones

**Lo que calcula el FFT** (interpolado con Peskin):

El campo eléctrico en cada nodo es E_mesh(r_i) = -∇ψ(r_i).
La partícula siente el campo interpolado:
```
E_solver = Σ_i [ w_i · E_mesh(r_i) ]
         = Σ_i [ w_i · (-∇ψ_solver)(r_i) ]
```
donde w_i = d_peskin(r_p - r_i).

Expandiendo, el campo de un ion en r_n interpolado a la partícula:
```
E_solver_from_ion = Σ_i [ w_i · q_ion · (-∇G_solver)(r_i - r_n) ]
```

**Lo que debería ser (Coulomb directo)**:

El ion está en el nodo r_n (no distribuido), así que el campo real es:
```
E_Coulomb = q_ion · (-∇G_Coulomb)(r_p - r_n)
          = q_ion / (4πε) · (r_p - r_n) / |r_p - r_n|³
```

**Corrección del campo** (para fuerza en partícula):
```
E_corr_from_ion = E_Coulomb - E_solver_from_ion

                = q_ion · (-∇G_Coulomb)(r_p - r_n)
                  - Σ_i [ w_i · q_ion · (-∇G_solver)(r_i - r_n) ]

                = q_ion · [ (-∇G_Coulomb)(r_p - r_n)
                           - Σ_i w_i · (-∇G_solver)(r_i - r_n) ]
```

**Nota sobre G_solver y G_diff**: Definimos G_diff como la diferencia entre Coulomb y el solver:
```
G_diff(r) ≡ G_Coulomb(r) - G_solver(r)
```

Esto es una **definición operacional**, no una identidad exacta. Representa la regularización implícita que el kernel discreto hace de la singularidad de Coulomb en la escala de la celda.

Por lo tanto:
```
∇G_solver = ∇G_Coulomb - ∇G_diff
```

### 3.3 Self-Field de la Partícula

**Situación**: La partícula siente un campo espurio de su propia carga distribuida.

**Lo que calcula el FFT**:

1. La carga q_p se distribuye a nodos j con pesos w_j
2. El FFT calcula el potencial en todos los nodos
3. El campo se interpola de vuelta a la partícula con pesos w_i

```
E_self_solver = Σ_i [ w_i · Σ_j [ q_p · w_j · (-∇G_solver)(r_i - r_j) ] ]
              = q_p · Σ_i Σ_j [ w_i · w_j · (-∇G_solver)(r_i - r_j) ]
```

**Lo que debería ser**:

En teoría clásica, el self-field de una carga puntual es infinito. En simulaciones de partículas, se adopta la **convención física-numérica** de que la fuerza de auto-interacción es cero tras la regularización:
```
E_self_real = 0   (por convención/regularización)
```

**Corrección del self-field**:
```
E_self_corr = 0 - E_self_solver
            = -q_p · Σ_i Σ_j [ w_i · w_j · (-∇G_solver)(r_i - r_j) ]
```

**Nota**: Este término es **puramente numérico** - surge del procedimiento de asignación-solver-interpolación y no tiene contrapartida directa en la teoría continua. Es un artefacto del método que debe restarse.

---

## 4. Resumen de Correcciones Necesarias

### 4.1 Corrección del Potencial (para Nernst-Planck)

Para cada nodo r_n cercano a la partícula:
```
Δψ(r_n) = q_p · [ G_Coulomb(r_n - r_p) - Σ_j w_j · G_solver(r_n - r_j) ]
```

donde j recorre el soporte de Peskin alrededor de r_p.

### 4.2 Corrección del Campo (para fuerza en partícula)

**Parte 1: Corrección PN (partícula-nodo)**

Para cada nodo r_n con carga iónica q_ion:
```
E_PN_corr = Σ_nodos_cercanos q_ion · [ (-∇G_Coulomb)(r_p - r_n)
                                       - Σ_i w_i · (-∇G_solver)(r_i - r_n) ]
```

**Parte 2: Corrección del Self-Field**
```
E_self_corr = -q_p · Σ_i Σ_j [ w_i · w_j · (-∇G_solver)(r_i - r_j) ]
```

**Campo total corregido**:
```
E_total = E_mesh_interpolado + E_PN_corr + E_self_corr
```

### 4.3 Conservación de Momento (Newton III)

**Para la corrección PN:**

La fuerza sobre la partícula es:
```
F_particle = q_p · E_PN_corr
```

La reacción debe ir **directamente al nodo** (NO redistribuida con Peskin):
```
F_node = -q_p · E_PN_corr_from_this_node
```

**Justificación**: La interacción física es entre una partícula (tratada como puntual para la corrección) y un ion localizado exactamente en ese nodo. Si redistribuyéramos la reacción con Peskin, violaríamos el cierre local de momento aunque el momento total global se conserve.

**Para el self-field:**

No hay reacción externa porque es una auto-interacción. Solo se corrige el campo que ve la partícula sin aplicar fuerza al fluido.

---

## 5. Funciones de Green Necesarias

### 5.1 G_Coulomb (analítica)
```
G_Coulomb(r) = 1 / (4πε|r|)
∇G_Coulomb(r) = -r / (4πε|r|³)
```

### 5.2 G_solver (numérica, suma de Fourier)
```
G_solver(r) = (1/N³) · Σ_k [ cos(k·r) / (ε · λ(k)) ]
```

El gradiente:
```
∇G_solver(r) = -(1/N³) · Σ_k [ k · sin(k·r) / (ε · λ(k)) ]
```

### 5.3 G_diff = G_Coulomb - G_solver

Por definición:
```
G_diff(r) ≡ G_Coulomb(r) - G_solver(r)
```

Esta es la función que `psi_fft_pn_compute_G_diff_direct()` calcula actualmente.

Por lo tanto:
```
∇G_solver = ∇G_Coulomb - ∇G_diff
```

---

## 6. Diferencia con la Implementación Actual

### 6.1 Lo que hace `psi_fft_pn_correct_field()` actualmente

```c
E_corr = Σ_nodos [ q_node · grad(G_diff)(r_p - r_node) ]
```

Esto asume que el campo de corrección es simplemente `q_node · ∇G_diff` evaluado en la posición de la partícula.

### 6.2 Lo que debería hacer

```c
// Para cada nodo con carga iónica q_ion:
E_PN_corr_from_node = q_ion · [ (-∇G_Coulomb)(r_p - r_n)
                               - Σ_i w_i · (-∇G_solver)(r_i - r_n) ]
```

La diferencia clave es que el término del solver debe considerar la **interpolación de Peskin**: el campo que la partícula siente del solver es un promedio ponderado del campo en los nodos de interpolación, no el campo en la posición exacta de la partícula.

### 6.3 Lo que falta: Self-Field

La implementación actual NO resta el self-field:
```c
E_self = q_p · Σ_i Σ_j [ w_i · w_j · (-∇G_solver)(r_i - r_j) ]
```

Este término es crucial porque el `Esub` que calcula `subgrid_update_Esub()` incluye la contribución de la propia carga de la partícula distribuida en el mesh.

---

## 7. Pseudocódigo de la Corrección Completa

```c
int psi_fft_pn_correct_field_complete(psi_fft_pn_t * pn,
                                       colloids_info_t * cinfo,
                                       hydro_t * hydro) {

  for (cada partícula pc) {

    q_p = pc->s.q0 - pc->s.q1;
    r_p = pc->s.r;  // posición de partícula

    // Arrays para acumular correcciones
    E_PN_corr = {0, 0, 0};
    E_self = {0, 0, 0};

    // Calcular pesos de Peskin para esta partícula
    // Soporte: nodos con índices en [floor(r_p) - 1, floor(r_p) + 2] en cada dim
    for (i en soporte Peskin) {
      w[i] = d_peskin(r_p - r_i);
    }

    // ================================================
    // PARTE 1: SELF-FIELD (restar)
    // ================================================
    // E_self = q_p * Σ_i Σ_j [ w_i * w_j * (-∇G_solver)(r_i - r_j) ]
    //
    // Nota: ∇G_solver = ∇G_Coulomb - ∇G_diff
    //
    // Optimización: Como ∇G_solver es impar (∇G(-r) = -∇G(r)),
    // los términos (i,j) y (j,i) se cancelan parcialmente.
    // Se puede explotar esta simetría para reducir el costo.

    for (i en soporte Peskin) {      // nodos de interpolación
      for (j en soporte Peskin) {    // nodos donde está la carga

        r_ij = r_i - r_j;

        if (|r_ij| < tolerancia) {
          // Mismo nodo: ∇G_solver(0) = 0 por simetría
          continue;
        }

        // ∇G_solver = ∇G_Coulomb - ∇G_diff
        gradG_Coulomb = -r_ij / (4πε|r_ij|³);
        gradG_diff = compute_grad_G_diff_direct(r_ij);
        gradG_solver = gradG_Coulomb - gradG_diff;

        E_self += w[i] * w[j] * q_p * (-gradG_solver);
      }
    }

    // ================================================
    // PARTE 2: CORRECCIÓN PN (sumar)
    // ================================================
    // Para cada nodo con carga iónica

    for (n en nodos cercanos con |r_p - r_n| < cutoff) {

      q_ion = rho0[n] - rho1[n];
      if (|q_ion| < tolerance) continue;

      r_pn = r_p - r_n;  // vector nodo -> partícula

      // Campo Coulomb directo desde el ion a la partícula
      E_Coulomb_direct = q_ion * (-grad_G_Coulomb)(r_pn);
      //                = q_ion * r_pn / (4πε|r_pn|³)

      // Campo del solver interpolado desde el ion
      // (lo que la partícula "ve" a través del mesh)
      E_solver_interpolated = 0;
      for (i en soporte Peskin) {
        r_in = r_i - r_n;  // vector nodo_ion -> nodo_interpolación

        if (|r_in| < tolerancia) {
          // Ion en el mismo nodo que el de interpolación
          // ∇G_solver(0) = 0 por simetría
          continue;
        }

        gradG_Coulomb = -r_in / (4πε|r_in|³);
        gradG_diff = compute_grad_G_diff_direct(r_in);
        gradG_solver = gradG_Coulomb - gradG_diff;

        E_solver_interpolated += w[i] * q_ion * (-gradG_solver);
      }

      // Corrección para este ion
      E_corr_from_ion = E_Coulomb_direct - E_solver_interpolated;
      E_PN_corr += E_corr_from_ion;

      // Newton III: reacción DIRECTA al nodo (no redistribuir con Peskin)
      F_on_node = -q_p * E_corr_from_ion * (kT/e);
      hydro_f_local_add(hydro, index_n, F_on_node);
    }

    // ================================================
    // APLICAR CORRECCIONES
    // ================================================
    // Esub actualmente tiene: E_mesh interpolado (incluye self-field)
    // Corrección: E_correct = E_mesh - E_self + E_PN_corr

    pc->Esub[X] = pc->Esub[X] - E_self[X] + E_PN_corr[X];
    pc->Esub[Y] = pc->Esub[Y] - E_self[Y] + E_PN_corr[Y];
    pc->Esub[Z] = pc->Esub[Z] - E_self[Z] + E_PN_corr[Z];
  }

  return 0;
}
```

---

## 8. Notas de Implementación

### 8.1 Cálculo de ∇G_solver

Como G_solver = G_Coulomb - G_diff (por nuestra definición de G_diff), tenemos:
```
∇G_solver = ∇G_Coulomb - ∇G_diff
```

donde:
- `∇G_Coulomb(r) = -r / (4πε|r|³)` es analítico
- `∇G_diff` ya está implementado en `psi_fft_pn_compute_grad_G_diff_direct()`

### 8.2 Eficiencia y Costo Computacional

**Costo actual (método directo):**
- El self-field tiene O(64²) = O(4096) pares (i,j) por partícula
- La corrección PN tiene O(N_cutoff × 64) términos por partícula
- Cada evaluación de ∇G_diff es O(N³) con suma de Fourier directa

**Optimizaciones posibles:**

1. **Simetría del self-field**: Como ∇G_solver es impar, los términos (i,j) y (j,i) satisfacen:
   ```
   w_i·w_j·∇G_solver(r_i - r_j) + w_j·w_i·∇G_solver(r_j - r_i) = 0  si w_i = w_j
   ```
   Para Peskin, w_i ≠ w_j en general, pero se puede explotar la antisimetría de ∇G para reducir cálculos.

2. **Tabla precalculada**: Precalcular ∇G_solver (o ∇G_diff) en offsets enteros y para offsets fraccionarios típicos. Como los offsets (r_i - r_j) y (r_i - r_n) son combinaciones de enteros pequeños, una tabla de tamaño moderado cubriría todos los casos.

3. **Modelo analítico aproximado**: Usar una aproximación tipo erf(κr)/r² como en el documento original SHORT_RANGE_CORRECTIONS_ANALYSIS.md, calibrada para coincidir con G_diff en el rango de interés.

### 8.3 Verificación de la Implementación

Para validar:

1. **Partícula aislada sin iones**: E_total debería ser ≈ 0 (solo queda self-field, que se resta)

2. **Partícula en el centro de un nodo**: El self-field debería ser exactamente 0 por simetría

3. **Conservación de momento**: Verificar que Σ F_partículas + Σ F_fluido = 0

4. **Comparación con Coulomb directo**: Para una partícula y un ion, verificar que la fuerza total converge a Coulomb cuando el mesh se refina

---

## 9. Solución Implementada: Corrección del Potencial

### 9.1 Enfoque Elegido

En lugar de corregir directamente el campo E (fuerza), **corregimos el potencial ψ** cerca de las partículas. Este enfoque tiene ventajas fundamentales:

1. **Nernst-Planck usa ψ, no fuerzas**: El transporte iónico depende del gradiente del potencial electroquímico μ = μ_s + z·e·ψ. Aplicar fuerzas al fluido no afecta N-P.

2. **Self-field automáticamente eliminado**: El potencial de Coulomb de una carga puntual en r_p es:
   ```
   ψ_Coulomb(r_n) = q_p / (4πε|r_n - r_p|)
   ```
   Este potencial NO tiene componente de auto-interacción en r = r_p (o más bien, diverge pero físicamente no representa una fuerza sobre sí misma).

3. **Newton III implícito**: Si el ion en r_n ve el potencial correcto de la partícula, y la partícula ve el potencial correcto del ion, la reciprocidad está garantizada por el teorema de reciprocidad de Green.

4. **Esub calculado correctamente**: Después de corregir ψ, la interpolación estándar de Peskin para calcular Esub = Σ_i w_i E(r_i) = -Σ_i w_i ∇ψ(r_i) automáticamente da el campo correcto.

### 9.2 Fórmula Implementada

Para cada nodo r_n cercano a la partícula (dentro del cutoff):

```
Δψ(r_n) = q_p · [ G_Coulomb(r_n - r_p) - Σ_j w_j · G_solver(r_n - r_j) ]
```

donde:
- `G_Coulomb(r) = 1/(4πε|r|)` es la función de Green de Coulomb exacta
- `G_solver(r) = G_Coulomb(r) - G_diff(r)` es la función de Green del solver FFT
- `w_j = d_peskin(r_p - r_j)` son los pesos de Peskin para distribuir la carga
- La suma en j es sobre los 64 nodos del soporte de Peskin

**Interpretación física**: Restamos el potencial que el FFT calculó (de la carga distribuida) y sumamos el potencial de Coulomb correcto (de la carga puntual).

### 9.3 Caso Especial: G_solver(0)

Cuando el nodo fuente j coincide con el nodo objetivo n (r_n = r_j), necesitamos G_solver(0).

A diferencia de G_Coulomb(0) que diverge, G_solver(0) es **finito**:
```
G_solver(0) = (1/N³) · Σ_k [1 / (ε · λ(k))]    para k ≠ 0
```

Este es el "auto-potencial" o "self-energy" del mesh, y es la regularización implícita que el método FFT hace de la singularidad de Coulomb.

En la implementación, este valor se pre-calcula y cachea en `pn->G_solver_zero` durante la inicialización.

### 9.4 Orden de Operaciones en ludwig.c

El orden correcto en el timestep es:

1. `psi_solver_fft_solve()` - Resuelve Poisson con carga distribuida via Peskin
2. `subgrid_charge_from_particles_restore()` - Quita la carga de partículas de los nodos
3. **`psi_fft_pn_correction()`** - Corrige ψ cerca de partículas (NUEVO: antes de Esub)
4. `subgrid_update_Esub()` - Calcula Esub interpolando desde ψ corregido
5. `colloid_sums_halo()` - Sincroniza Esub entre procesos

### 9.5 Código Implementado

La función `psi_fft_pn_correct_potential()` en `psi_fft_pn.c` implementa esta corrección:

```c
/* Para cada partícula y cada nodo cercano: */
double G_Coulomb_pn = 1.0 / (4.0 * M_PI * epsilon * dist_pn);

/* Suma sobre nodos Peskin */
double G_solver_sum = 0.0;
for (j en soporte_Peskin) {
    w_j = d_peskin(r_p - r_j);
    if (j == n) {
        G_solver_jn = pn->G_solver_zero;  // valor cacheado
    } else {
        G_solver_jn = G_Coulomb(r_n - r_j) - G_diff(r_n - r_j);
    }
    G_solver_sum += w_j * G_solver_jn;
}

delta_psi = q_p * (G_Coulomb_pn - G_solver_sum);
psi->data[addr] += delta_psi;
```

### 9.6 Verificación

Para validar la corrección:

1. **Partícula centrada en nodo**: La corrección debería ser mínima (distribución Peskin más simétrica)

2. **Partícula entre nodos**: La corrección debería compensar el error de la distribución Peskin

3. **Esub de partícula aislada (sin iones)**: Debería ser ≈ 0 después de la corrección (sin self-field)

4. **Fuerza partícula-ion**: Debería converger a Coulomb directo q_p·q_ion/(4πε|r_p - r_ion|²)

5. **Conservación de momento**: Σ F_partículas = -Σ F_fluido (fuerzas via gradmu o divstress)

---

## 10. Referencias

1. Hockney, R. W., & Eastwood, J. W. (1988). Computer Simulation Using Particles. CRC Press. Capítulo 8.

2. Peskin, C. S. (2002). The immersed boundary method. Acta Numerica, 11, 479-517.

3. Documento interno: `/home/bater/Sim/ludwig/docs/SHORT_RANGE_CORRECTIONS_ANALYSIS.md`
