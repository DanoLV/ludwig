# Análisis de Correcciones de Corto Alcance para Partículas Subgrid en Ludwig

## Resumen Ejecutivo

Este documento analiza las alternativas para implementar correcciones de corto alcance
en las interacciones eléctricas de partículas subgrid en el código Ludwig. El objetivo
es mejorar la precisión del cálculo de fuerzas eléctricas sin perder la conservación
de masa y momento.

---

## 1. Contexto y Problema

### 1.1 Estado Actual del Código

El código Ludwig implementa electrocinética para partículas subgrid mediante:

1. **Distribución de carga (Peskin)**: La carga de cada partícula se distribuye a los 
   nodos vecinos usando la función delta regularizada de Peskin (`d_peskin()` en
   `src/subgrid.c`).

2. **Solver de Poisson** (PETSc): Resuelve ∇²ψ = -ρ_elec/ε para obtener el potencial
   eléctrico en todo el dominio.

3. **Nernst-Planck**: Transporta los iones siguiendo el esquema de Capuani et al. (2004)
   con flujos calculados en los enlaces entre nodos.

4. **Interpolación del campo**: El campo eléctrico sobre la partícula se calcula
   interpolando E_mesh desde los nodos vecinos con pesos Peskin.

### 1.2 El Problema Fundamental

Cuando se distribuye la carga de una partícula con la función de Peskin:

- El campo eléctrico calculado por el solver de Poisson es una **versión suavizada**
  del campo real.
- Los nodos cercanos a la partícula contienen el **self-field** (campo que la partícula
  genera sobre sí misma).
- La interacción con cargas iónicas cercanas está **suavizada** respecto a la
  interacción Coulombiana real.

Esto introduce errores sistemáticos en el cálculo de fuerzas, especialmente a
distancias cortas.

### 1.3 Sistema Fuera de Equilibrio

**IMPORTANTE**: El sistema electrocinético de Ludwig está **fuera de equilibrio**.
Hay flujo de fluido, campos eléctricos externos, gradientes de concentración iónica,
y las partículas se mueven. Esto tiene implicaciones fundamentales:

- **NO se puede usar el potencial de Yukawa/Debye-Hückel** para las correcciones.
  Yukawa asume que los iones están en equilibrio de Boltzmann alrededor de las cargas,
  lo cual no es cierto en nuestro sistema.

- La distribución de iones (ρ₀, ρ₁) es calculada dinámicamente por Nernst-Planck
  y **no sigue** la distribución de Boltzmann ρ ∝ exp(-zeψ/kT).

- Las correcciones de corto alcance deben usar **Coulomb puro** (1/r²), no Yukawa
  (exp(-κr)/r²), ya que el apantallamiento dinámico ya está incluido en la solución
  de Nernst-Planck + Poisson.

### 1.4 Base Teórica: Método P³M

Según Hockney & Eastwood ("Computer Simulation Using Particles", Cap. 8), el método
Particle-Particle-Particle-Mesh (P³M) divide la fuerza total en:

```
F_total = F_mesh (largo alcance) + F_PP (corto alcance)
```

Donde:
- **F_mesh**: Calculada eficientemente via FFT/solver de Poisson (ya implementado)
- **F_PP**: Corrección directa para interacciones de corto alcance

La corrección de corto alcance compensa el error introducido por:
1. La asignación de carga al mesh (charge assignment)
2. La interpolación del campo desde el mesh
3. La resolución finita del mesh

**Nota sobre el potencial**: En P³M clásico para sistemas en equilibrio se usa Ewald
con Yukawa. En nuestro caso fuera de equilibrio, la corrección PP debe ser Coulomb
puro, y el "apantallamiento" emerge naturalmente de la dinámica de Nernst-Planck.

### 1.5 Distinción Crítica: κ_Ewald vs κ_Debye-Hückel

Es fundamental distinguir entre dos usos diferentes de "κ":

1. **κ_E (Ewald splitting parameter)**: Parámetro **numérico** que controla cómo se
   divide la interacción entre corto y largo alcance. Se elige por eficiencia
   computacional, típicamente κ_E ≈ 2/Δx donde Δx es el tamaño de malla.

2. **κ_DH (Debye-Hückel)**: Parámetro **físico** de apantallamiento que asume
   equilibrio de Boltzmann. **NO debe usarse** en nuestro sistema fuera de equilibrio.

**El kernel de Ewald erfc(κ_E r)/r puede usarse** como función de splitting numérica
sin violar la física. El κ_E no representa apantallamiento físico, solo divide
matemáticamente la interacción para evitar doble conteo.

```
Corrección = F_Coulomb_real - F_reference_mesh
           = (1/r²) - erfc(κ_E r)/r²    para r < r_c
           = 0                           para r > r_c
```

Esto garantiza que:
- La corrección es exactamente cero más allá del radio de corte r_c
- No se altera la física de largo alcance que el mesh resuelve bien
- El apantallamiento físico sigue emergiendo de Nernst-Planck

---

## 2. Tipos de Interacciones que Requieren Corrección

### 2.1 Self-Field (Partícula consigo misma)

**Descripción**: La carga de la partícula, distribuida en los nodos vecinos con Peskin,
genera un campo eléctrico que actúa sobre la propia partícula.

**Magnitud**: Este es típicamente el efecto dominante a corta distancia.

**Estado en el código**: Parcialmente implementado en `subgrid_compute_self_field_single_particle()`.

```
E_self = Σ_nodos [ q_p × d_peskin(r_p - r_n) × E_from_distributed_charge(r_n → r_p) ]
```

### 2.2 Particle-Particle (PP)

**Descripción**: Interacción directa entre dos partículas subgrid cercanas. El campo
del mesh suaviza esta interacción.

**Corrección necesaria**:
```
F_PP_corr = F_Coulomb_real(r_12) - F_mesh_interpolated(r_12)
```

**Referencia en código**: El archivo `src/ewald.c` implementa sumas de Ewald para
dipolos magnéticos. La estructura de `ewald_real_space_sum()` puede adaptarse para
cargas.

### 2.3 Particle-Node (PN) - Interacción con Cargas Iónicas

**Descripción**: Interacción entre una partícula subgrid y la densidad de carga iónica
(ρ₀ - ρ₁) en los nodos cercanos de la red.

**El problema**:
- La carga iónica está localizada en nodos discretos
- El solver de Poisson calcula un campo suavizado
- La interacción real partícula-nodo es más "puntual" que lo que el mesh captura

**Corrección necesaria**:
```
F_PN_corr = Σ_nodos_cercanos [ F_Coulomb_directo(r_p, r_n, ρ_net)
                              - F_mesh_interpolated(r_p, r_n) ]
```

**IMPORTANTE - Evitar Doble Conteo**: La corrección PN debe usar el **kernel de Ewald**
(erfc), NO la fuerza de Coulomb total. El solver de Poisson ya calcula la parte de
largo alcance. Si usamos Coulomb total, contaríamos dos veces esa contribución:

```
F_PN_corr = F_shortrange_only = [1/r² - erfc(κ_E r)/r²] × q_p × Q_nodo
                              = erf(κ_E r)/r² × q_p × Q_nodo
```

Donde Q_nodo = ρ_net × Δx³ es la carga efectiva en el nodo.

---

## 3. Alternativas de Implementación

### 3.1 ALTERNATIVA A: Solo Corrección de Self-Field

**Descripción**: Implementar únicamente la corrección del self-field, asumiendo que
las interacciones PP y PN están razonablemente bien capturadas por el mesh.

**Implementación**:
```c
// En subgrid_update_Esub() o función nueva
for (cada partícula pc) {
    E_self = {0, 0, 0};

    for (cada nodo cercano n) {
        dr = d_peskin(r_p - r_n);
        q_at_node = q_p * dr;  // Carga asignada al nodo

        if (dist > epsilon) {
            // Campo Coulomb que esta carga genera en la posición de la partícula
            // E = (1/4πε) * q / r² * r_hat
            factor = prefactor / (dist * dist);
            E_self += q_at_node * factor * (r_n - r_p) / dist;
        }
    }

    // Restar el self-field del campo total
    pc->Esub -= E_self;
}
```

**Ventajas**:
- Implementación simple
- No afecta la conservación (solo modifica el campo sobre la partícula)
- Corrige el efecto dominante

**Desventajas**:
- No corrige interacciones PP ni PN
- Puede ser insuficiente si hay muchas partículas cercanas o alta concentración iónica

**Archivos a modificar**: `src/subgrid.c`

---

### 3.2 ALTERNATIVA B: Correcciones PP + Self-Field (Estilo Ewald)

**Descripción**: Implementar correcciones particle-particle siguiendo la estructura
del código de Ewald existente, más la corrección de self-field.

**Implementación**:
```c
// Nueva función: subgrid_shortrange_electric_sum()
int subgrid_shortrange_electric_sum(colloids_info_t* cinfo, psi_t* psi) {

    double rc = 2.0;  // Radio de corte para correcciones

    // 1. Corrección Self-Field para cada partícula
    for (cada partícula pc) {
        E_self = compute_self_field(pc, psi);
        pc->Esub -= E_self;
    }

    // 2. Corrección PP (estructura similar a ewald_real_space_sum)
    for (ic = 1; ic <= ncell[X]; ic++) {
        for (jc = 1; jc <= ncell[Y]; jc++) {
            for (kc = 1; kc <= ncell[Z]; kc++) {

                // Obtener partícula en esta celda
                colloids_info_cell_list_head(cinfo, ic, jc, kc, &p_c1);

                for (; p_c1; p_c1 = p_c1->next) {
                    if (p_c1->s.bc != COLLOID_BC_SUBGRID) continue;

                    // Loop sobre celdas vecinas
                    for (dx = -1; dx <= +1; dx++) {
                        for (dy = -1; dy <= +1; dy++) {
                            for (dz = -1; dz <= +1; dz++) {

                                colloids_info_cell_list_head(cinfo, ic+dx, jc+dy, kc+dz, &p_c2);

                                for (; p_c2; p_c2 = p_c2->next) {
                                    if (p_c2->s.bc != COLLOID_BC_SUBGRID) continue;
                                    if (p_c1->s.index >= p_c2->s.index) continue;

                                    // Calcular separación
                                    cs_minimum_distance(cs, p_c1->s.r, p_c2->s.r, r12);
                                    r = sqrt(r12[X]² + r12[Y]² + r12[Z]²);

                                    if (r < rc) {
                                        // Fuerza Coulomb directa (sin Yukawa - sistema fuera de equilibrio)
                                        // F = (1/4πε) * q1*q2 / r² * r_hat
                                        q1 = p_c1->s.q0 - p_c1->s.q1;
                                        q2 = p_c2->s.q0 - p_c2->s.q1;

                                        F_direct = q1 * q2 * prefactor / (r * r) * r12/r;

                                        // Restar lo que ya cuenta el mesh (aproximación)
                                        F_mesh = interpolate_mesh_force(p_c1, p_c2, psi);

                                        F_corr = F_direct - F_mesh;

                                        // Aplicar fuerzas (Newton III)
                                        p_c1->fex += F_corr;
                                        p_c2->fex -= F_corr;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return 0;
}
```

**Ventajas**:
- Corrige tanto self-field como interacciones PP
- Reutiliza la estructura de celdas existente para búsqueda eficiente de vecinos
- Patrón ya probado en el código de Ewald

**Desventajas**:
- No corrige interacciones PN con iones
- La estimación de F_mesh entre partículas es aproximada
- Requiere distribuir la fuerza de reacción al fluido para conservación

**Archivos a modificar**: `src/subgrid.c`, `src/subgrid.h`

---

### 3.3 ALTERNATIVA C: Corrección PN via Campo Eléctrico Efectivo

**Descripción**: Modificar el campo eléctrico que "ve" la partícula para incluir
correcciones de la interacción con cargas iónicas cercanas.

**Implementación**:
```c
int subgrid_update_Esub_with_PN_correction(colloids_info_t* cinfo, psi_t* psi) {

    for (cada partícula pc) {

        E_total = {0, 0, 0};
        q_p = pc->s.q0 - pc->s.q1;

        for (cada nodo cercano (i,j,k)) {

            index = cs_index(cs, i, j, k);
            r_pn = r_nodo - r_particula;
            dist = |r_pn|;
            dr = d_peskin(r_p[X]-i) * d_peskin(r_p[Y]-j) * d_peskin(r_p[Z]-k);

            // === Campo del mesh interpolado (método actual) ===
            psi_electric_field(psi, index, E_mesh);
            E_mesh_interp = E_mesh * dr;

            // === Corrección Self-Field (Coulomb puro) ===
            if (dist > 0.1) {
                q_at_node = q_p * dr;
                // E = (1/4πε) * q / r² * r_hat (sin screening Yukawa)
                factor = prefactor / (dist * dist);
                E_self = q_at_node * factor * (-r_pn / dist);  // Campo sobre partícula
            }

            // === Corrección PN (Coulomb puro) ===
            psi_rho(psi, index, 0, &rho0);
            psi_rho(psi, index, 1, &rho1);
            rho_net = rho0 - rho1;

            if (dist > 0.1 && |rho_net| > 1e-15) {
                // Campo Coulomb directo desde la carga iónica en el nodo
                // (sin screening - el apantallamiento emerge de Nernst-Planck)
                factor = prefactor / (dist * dist);
                E_PN_direct = rho_net * factor * (-r_pn / dist);

                // Lo que el mesh ya cuenta (aproximación)
                // El mesh "ve" la carga iónica distribuida, no puntual
                E_PN_mesh = E_mesh * dr;  // Ya incluido en E_mesh_interp

                // Corrección: diferencia entre directo y mesh
                // NOTA: Esta es una aproximación, ver discusión abajo
                E_PN_corr = E_PN_direct - E_PN_mesh;
            }

            // Acumular campo total corregido
            E_total += E_mesh_interp - E_self + E_PN_corr;
        }

        pc->Esub = E_total;
    }
}
```

**El problema sutil con E_PN_corr**:

El campo `E_mesh` en cada nodo es el resultado del solver de Poisson considerando
TODAS las cargas. No es trivial separar "qué parte del campo en el nodo n viene
de la carga en el nodo m".

Para una corrección PN exacta necesitaríamos la **función de Green del mesh discreto**,
que describe cómo una carga unitaria en un nodo contribuye al potencial/campo en
todos los demás nodos.

**Aproximaciones posibles**:

a) **Aproximación de campo local**: Asumir que E_mesh en un nodo representa principalmente
   el campo de cargas cercanas, y la corrección es pequeña.

b) **Función de Green analítica**: Usar la función de Green del Laplaciano continuo
   como aproximación:
   ```
   G(r) = 1/(4πε|r|)  (Coulomb puro - apropiado para sistema fuera de equilibrio)
   ```
   Nota: NO usar Debye-Hückel G(r) = exp(-κ|r|)/(4πε|r|) porque el sistema está
   fuera de equilibrio y el apantallamiento emerge de Nernst-Planck.

c) **Pre-cálculo de la función de Green del mesh**: Calcular numéricamente cómo
   responde el solver de Poisson a una carga delta en el origen. Esto se hace una
   vez y se tabula.

**Ventajas**:
- Aborda las tres fuentes de error (self, PP implícito, PN)
- Mantiene la estructura del cálculo de fuerzas (F = qE)
- No requiere modificar Nernst-Planck directamente

**Desventajas**:
- La corrección PN exacta es compleja
- Requiere aproximaciones cuya precisión es difícil de evaluar
- Puede introducir inconsistencias si las aproximaciones son pobres

**Archivos a modificar**: `src/subgrid.c`

---

### 3.4 ALTERNATIVA D: Corrección Completa P³M con Función de Referencia

**Descripción**: Implementación completa del esquema P³M de Hockney, definiendo
explícitamente una función de referencia R(r) que el mesh calcula bien.

**Teoría**:

En P³M, la fuerza total se escribe como:
```
F_total = F_mesh + F_PP
F_PP = F_real - F_reference
```

Donde:
- `F_real`: Fuerza Coulombiana pura (1/r²) - **NO Yukawa** ya que el sistema está fuera de equilibrio
- `F_reference`: Fuerza de una distribución de carga suavizada que el mesh resuelve exactamente

La función de referencia típica es la convolución de la carga puntual con la función
de asignación al mesh (Peskin en nuestro caso):
```
ρ_reference(r) = q × W(r)  donde W = d_peskin ⊗ d_peskin ⊗ d_peskin (3D)
```

**Implementación**:
```c
// Precálculo (una vez al inicio)
void compute_reference_force_table(double rc, double dr) {
    // Tabular F_reference(r) para r ∈ [0, rc]
    // F_reference = ∫ ρ_reference(r') × E_Coulomb(r-r') dr'
    // donde E_Coulomb = (1/4πε) / r² * r_hat (Coulomb puro)
    // Esta integral se puede hacer numéricamente o analíticamente para Peskin
}

// En tiempo de ejecución
int subgrid_pp_correction(colloids_info_t* cinfo, psi_t* psi) {

    for (cada par de partículas cercanas (p1, p2)) {
        r = distancia(p1, p2);

        if (r < rc) {
            q1 = p1->s.q0 - p1->s.q1;
            q2 = p2->s.q0 - p2->s.q1;

            // Fuerza Coulomb real (sin Yukawa - sistema fuera de equilibrio)
            F_real = q1 * q2 * coulomb_force(r);  // = prefactor / r²

            // Fuerza de referencia (tabulada)
            F_ref = q1 * q2 * interpolate_reference_force(r);

            // Corrección
            F_corr = F_real - F_ref;

            p1->fex += F_corr * r_hat;
            p2->fex -= F_corr * r_hat;
        }
    }

    // Similar para self-field y PN...
}
```

**Ventajas**:
- Método más riguroso y fundamentado teóricamente
- Errores bien caracterizados y controlables
- Escalable a sistemas grandes

**Desventajas**:
- Implementación más compleja
- Requiere precálculo de tablas
- La función de referencia para Peskin 3D requiere integración numérica

**Archivos a modificar**: `src/subgrid.c`, `src/subgrid.h`, nuevo archivo `src/p3m_reference.c`

---

### 3.5 ALTERNATIVA E: Método Híbrido para Conservación

**Descripción**: Combinar correcciones de fuerza sobre partículas con modificaciones
a los flujos de Nernst-Planck para mantener conservación estricta.

**El problema de conservación**:

Cuando agregamos una fuerza F_corr sobre una partícula, por Newton III debe haber
una fuerza de reacción -F_corr sobre "algo". Las opciones son:

1. **Sobre el fluido**: Distribuir -F_corr a los nodos con pesos Peskin y agregar
   a `hydro_f_local_add()`. Conserva momento total fluido+partículas.

2. **Sobre los iones**: Agregar -F_corr como fuerza sobre las especies iónicas.
   Pero Nernst-Planck usa flujos en enlaces, no fuerzas en nodos.

**Implementación para conservación de momento**:
```c
int subgrid_shortrange_corrections_conserving(colloids_info_t* cinfo,
                                               psi_t* psi,
                                               hydro_t* hydro) {

    for (cada partícula pc) {

        // Calcular todas las correcciones
        F_self_corr = compute_self_field_correction(pc, psi);
        F_PP_corr = compute_PP_corrections(pc, cinfo, psi);
        F_PN_corr = compute_PN_corrections(pc, psi);

        F_total_corr = F_self_corr + F_PP_corr + F_PN_corr;

        // Aplicar a la partícula
        pc->fex += F_total_corr;

        // Distribuir reacción al fluido (conserva momento)
        for (cada nodo cercano n) {
            dr = d_peskin(r_p - r_n);
            F_reaction = -F_total_corr * dr;
            hydro_f_local_add(hydro, index_n, F_reaction);
        }
    }
}
```

**Para conservación de masa iónica** (más complejo):

La fuerza sobre los iones debería traducirse en un flujo adicional en Nernst-Planck:
```c
// En nernst_planck_fluxes_d3qx(), agregar:
j_correction = (D / kT) * rho * F_PN_reaction
```

Pero esto requiere modificar el esquema de Capuani, que usa flujos simétricos en
enlaces para garantizar equilibrio Boltzmann.

**Ventajas**:
- Conservación estricta de momento
- Puede extenderse a conservación de masa

**Desventajas**:
- La parte de conservación de masa es compleja
- Puede romper las propiedades del esquema de Capuani
- Requiere modificar múltiples archivos

**Archivos a modificar**: `src/subgrid.c`, `src/nernst_planck.c` (opcional)

---

### 3.6 ALTERNATIVA F: Solo Fuerzas sobre Partículas (Enfoque Pragmático)

**Descripción**: Implementar correcciones solo para las fuerzas sobre las partículas,
sin modificar Nernst-Planck. La reacción va al fluido via Peskin.

**Justificación**:
- Las partículas subgrid son mucho más masivas que los iones
- El efecto sobre la dinámica iónica de la fuerza de reacción es pequeño
- Nernst-Planck ya incluye el campo eléctrico del mesh, que captura la física principal

**Implementación**:
```c
// Nueva función en subgrid.c
int subgrid_electric_shortrange_correction(colloids_info_t* cinfo,
                                            psi_t* psi,
                                            hydro_t* hydro) {

    double prefactor = compute_coulomb_prefactor(psi);  // = e²/(4πε)
    double kappa_E = 2.0;  // Parámetro de splitting ≈ 2/Δx
    double rc = 3.0 / kappa_E;  // Radio de corte

    for (cada partícula pc) {

        double F_corr[3] = {0, 0, 0};
        double q_p = pc->s.q0 - pc->s.q1;

        // 1. Self-field correction
        for (cada nodo cercano n) {
            dr = d_peskin(r_p - r_n);
            q_n = q_p * dr;  // Carga de la partícula asignada al nodo
            if (dist > 0.1) {
                // Usar kernel de corrección erf(κr)/r² (NO Coulomb total)
                shortrange_correction_force(dist, q_n, q_p, prefactor, kappa_E, F_self);
                F_corr -= F_self;  // Restar self-field
            }
        }

        // 2. PP correction con otras partículas
        for (cada partícula vecina pc2) {
            r12 = distancia(pc, pc2);
            if (r12 < rc && pc->index < pc2->index) {
                q2 = pc2->s.q0 - pc2->s.q1;
                // Kernel de corrección: erf(κr)/r² (solo corto alcance)
                shortrange_correction_force(r12, q_p, q2, prefactor, kappa_E, F_PP_corr);

                F_corr += F_PP_corr;
                pc2->fex -= F_PP_corr;  // Newton III
            }
        }

        // 3. PN correction con nodos cargados (NUEVO)
        for (cada nodo cercano k con carga iónica) {
            r_pk = distancia(pc, nodo_k);
            if (r_pk < rc && r_pk > 0.1) {
                Q_nodo = (rho0[k] - rho1[k]) * dx * dx * dx;
                // Kernel de corrección: erf(κr)/r² (evita doble conteo)
                shortrange_correction_force(r_pk, q_p, Q_nodo, prefactor, kappa_E, F_PN_corr);

                F_corr += F_PN_corr;
                // Newton III: reacción DIRECTA al nodo k (NO usar Peskin)
                hydro_f_local_add(hydro, index_k, -F_PN_corr);
            }
        }

        // Aplicar corrección total a la partícula
        pc->fex += F_corr;
    }

    return 0;
}
```

**Ventajas**:
- Corrige las tres fuentes de error: Self-field, PP, y PN
- Usa kernel de corrección erf(κr)/r² que evita doble conteo
- Conserva momento total (Newton III directo al nodo para PN)
- No modifica Nernst-Planck
- Radio de corte pequeño (rc ≈ 1.5Δx) hace el cálculo eficiente

**Desventajas**:
- Requiere iterar sobre nodos cercanos para corrección PN
- Pequeña inconsistencia en la dinámica iónica (reacción va al fluido, no a NP)

**Archivos a modificar**: `src/subgrid.c`, `src/subgrid.h`

---

## 4. Comparación de Alternativas

| Alternativa | Complejidad | Self-Field | PP | PN | Conservación | Precisión |
|-------------|-------------|------------|----|----|--------------|-----------|
| A: Solo Self | Baja | ✓ | ✗ | ✗ | ✓ | Media |
| B: PP + Self (Ewald) | Media | ✓ | ✓ | ✗ | ✓* | Media-Alta |
| C: Campo Efectivo | Media-Alta | ✓ | ~ | ✓ | ✓ | Media |
| D: P³M Completo | Alta | ✓ | ✓ | ✓ | ✓ | Alta |
| E: Híbrido Conserv. | Alta | ✓ | ✓ | ✓ | ✓✓ | Alta |
| F: Pragmático Híbrido | Media | ✓ | ✓ | ✓ | ✓* | Alta |

*Conserva momento (Newton III directo), no masa iónica estrictamente.

---

## 5. Funciones Auxiliares Necesarias

Independientemente de la alternativa elegida, se necesitan estas funciones:

### 5.1 Fuerza Coulomb Pura

**IMPORTANTE**: Usar Coulomb puro, NO Yukawa. El sistema está fuera de equilibrio
y el apantallamiento emerge dinámicamente de Nernst-Planck.

```c
void coulomb_force(double r, double q1, double q2,
                   double prefactor, double F[3], double r_hat[3]) {
    // F = prefactor * q1 * q2 / r² * r_hat
    // donde prefactor = e² / (4πε)
    // NO usar exp(-κr) - el sistema NO está en equilibrio de Debye-Hückel
    double magnitude = prefactor * q1 * q2 / (r * r);
    F[X] = magnitude * r_hat[X];
    F[Y] = magnitude * r_hat[Y];
    F[Z] = magnitude * r_hat[Z];
}
```

**¿Por qué no Yukawa?**

El potencial de Yukawa/Debye-Hückel:
```
φ(r) = (q/4πε) * exp(-κr) / r
```
asume que los iones están en equilibrio de Boltzmann: ρ± ∝ exp(∓zeψ/kT).

En nuestro sistema electrocinético:
- Hay campos eléctricos externos aplicados
- Hay flujo de fluido (advección de iones)
- Hay gradientes de concentración impuestos
- Las partículas se mueven

Por lo tanto, la distribución de iones NO sigue Boltzmann, y el "apantallamiento"
efectivo es calculado dinámicamente por el solver de Nernst-Planck + Poisson.
Usar Yukawa sería contar el apantallamiento dos veces.

### 5.2 Kernel de Ewald para Corrección (Alternativa más Rigurosa)

Para evitar doble conteo, la corrección de corto alcance debe usar el kernel de Ewald:

```c
// Fuerza de corrección usando splitting de Ewald
// F_corr = F_Coulomb - F_reference = (1/r² - erfc(κ_E*r)/r²) * prefactor * q1 * q2
void shortrange_correction_force(double r, double q1, double q2,
                                  double prefactor, double kappa_E,
                                  double F_corr[3], double r_hat[3]) {

    // erfc(x) = 1 - erf(x), función de error complementaria
    double x = kappa_E * r;
    double erfc_x = erfc(x);

    // Derivada de erfc(κr)/r da el factor de fuerza
    // F = -∇[erfc(κr)/r] = [erfc(κr)/r² + (2κ/√π)exp(-κ²r²)/r] * r_hat
    double exp_x2 = exp(-x * x);
    double two_kappa_sqrtpi = 2.0 * kappa_E / sqrt(PI);

    // Fuerza de corrección: Coulomb real menos lo que el mesh ya calcula
    // F_corr = (1 - erfc(κr))/r² + términos de derivada
    double factor = (1.0 - erfc_x) / (r * r) + two_kappa_sqrtpi * exp_x2 / r;

    double magnitude = prefactor * q1 * q2 * factor;
    F_corr[X] = magnitude * r_hat[X];
    F_corr[Y] = magnitude * r_hat[Y];
    F_corr[Z] = magnitude * r_hat[Z];
}
```

**Nota**: κ_E es el parámetro de splitting numérico (≈ 2/Δx), NO el κ de Debye-Hückel.
Este kernel garantiza que la corrección sea cero para r > r_c (donde erfc(κ_E*r_c) ≈ 1).

### 5.3 Prefactor de Coulomb
```c
double compute_coulomb_prefactor(psi_t* psi) {
    double eunit, epsilon;
    psi_unit_charge(psi, &eunit);
    psi_epsilon(psi, &epsilon);
    // prefactor = e² / (4πε₀ε_r) en unidades de energía/distancia
    return (eunit * eunit) / (4.0 * PI * epsilon);
}
```

### 5.3 Longitud de Bjerrum (para referencia)
```c
double compute_bjerrum_length(psi_t* psi) {
    // l_B = e² / (4πεkT) - distancia a la cual energía Coulomb = kT
    double eunit, epsilon, beta;
    psi_unit_charge(psi, &eunit);
    psi_epsilon(psi, &epsilon);
    psi_beta(psi, &beta);  // beta = 1/kT
    return beta * (eunit * eunit) / (4.0 * PI * epsilon);
}
```

**Nota sobre κ (kappa)**: Aunque se puede calcular un κ "instantáneo" de la densidad
iónica promedio, este NO debe usarse para las correcciones de corto alcance porque
el sistema no está en equilibrio. Solo sería útil como escala de longitud característica
para análisis.

---

## 6. Consideraciones de Implementación

### 6.1 Estructura de Celdas para Búsqueda de Vecinos

El código ya tiene una estructura de celdas en `colloids_info_t` usada por el código
de Ewald. Esta misma estructura sirve para buscar partículas vecinas eficientemente.

### 6.2 Radio de Corte

El radio de corte `rc` para las correcciones debe ser:
- Mayor que el rango de la función de Peskin (~2 unidades de red)
- Menor que el tamaño de celda para búsqueda eficiente
- Típicamente rc ≈ 2-4 unidades de red

### 6.3 Precisión Numérica

El código ya usa sumas de Klein (`klein_t`) para evitar errores de cancelación.
Las correcciones de corto alcance deben usar el mismo esquema.

### 6.4 Paralelización

Las correcciones PP requieren comunicación entre procesadores si las partículas
están en diferentes dominios. El código de Ewald ya maneja esto.

### 6.5 Tercera Ley de Newton y Conservación del Momento (CRÍTICO)

**Este es el punto donde el problema de velocidades espurias puede reaparecer.**

Cuando se aplica una fuerza de corrección F_corr a una partícula, por Newton III
debe haber una fuerza de reacción -F_corr sobre "algo". La pregunta es: ¿sobre qué?

#### Para correcciones Partícula-Partícula (PP):
```c
// Newton III se cumple automáticamente
p_c1->fex += F_PP_corr;
p_c2->fex -= F_PP_corr;  // Reacción exacta
```

#### Para correcciones Partícula-Nodo (PN):

**RECOMENDACIÓN CRÍTICA**: NO usar Peskin para distribuir la reacción de correcciones
de corto alcance. Aplicar la reacción **directamente al nodo específico**.

```c
// CORRECTO: Reacción directa al nodo (conserva momento exactamente)
hydro_f_local_add(hydro, index_nodo_k, -F_PN_corr);

// INCORRECTO: Distribuir con Peskin (introduce inconsistencia)
// for (cada nodo cercano n) {
//     dr = d_peskin(r_p - r_n);
//     hydro_f_local_add(hydro, index_n, -F_PN_corr * dr);  // NO HACER
// }
```

**Justificación**: La corrección de corto alcance es una interacción "punto a punto"
(partícula en r_p ↔ nodo en r_n). La reacción debe ir exactamente al nodo involucrado
para conservar momento. Usar Peskin para distribuir la reacción introduciría una
inconsistencia: la fuerza viene de un nodo específico pero la reacción se distribuye
a varios.

#### Resumen de conservación:

| Tipo de corrección | Fuerza sobre partícula | Reacción |
|-------------------|------------------------|----------|
| Self-field | -F_self sobre partícula | Ninguna (es auto-interacción) |
| PP (part-part) | +F_PP_corr sobre p1 | -F_PP_corr sobre p2 |
| PN (part-nodo) | +F_PN_corr sobre partícula | -F_PN_corr **directo** al nodo k |

### 6.6 Elección del Parámetro κ_E (Ewald Splitting)

El parámetro κ_E controla la "suavidad" de la transición entre corto y largo alcance:

- **κ_E pequeño**: Corrección de corto alcance tiene mayor rango, más cálculos PP
- **κ_E grande**: Corrección más localizada, pero más sensible a posición exacta

**Recomendación práctica**:
```c
double kappa_E = 2.0 / delta_x;  // donde delta_x = 1.0 (tamaño de malla)
double rc = 3.0 / kappa_E;       // radio de corte donde erfc(κ_E*rc) ≈ 0
```

Esto da rc ≈ 1.5 Δx, que es consistente con el rango de la función de Peskin.

---

## 7. Plan de Implementación Incremental: De Simple a Preciso

Este plan permite implementar las alternativas de forma progresiva, empezando por
la más simple (A) y avanzando hacia la más precisa (E). Cada nivel agrega
funcionalidad sobre el anterior.

```
Nivel 1 (Alt. A) → Nivel 2 (Alt. B) → Nivel 3 (Alt. F) → Nivel 4 (Alt. D) → Nivel 5 (Alt. E)
   Self-field        + PP corr.        + PN corr.         + Func. Ref.       + Conserv. Masa
   ⭐                 ⭐⭐               ⭐⭐⭐              ⭐⭐⭐⭐            ⭐⭐⭐⭐⭐
```

---

### Nivel 0: Infraestructura Común (Prerrequisito para todos)

**Objetivo**: Crear funciones auxiliares que serán usadas por todos los niveles.

**Archivos**: `src/subgrid.c`, `src/subgrid.h`

```c
/* ========== INFRAESTRUCTURA COMÚN ========== */

#include <math.h>

/* Parámetros configurables (pueden moverse a input) */
#define SHORTRANGE_KAPPA_E  2.0    /* Parámetro de splitting numérico */
#define SHORTRANGE_RC       1.5    /* Radio de corte */

/* Kernel de corrección: erf(κr)/r² */
static double shortrange_erf_kernel(double r, double kappa_E) {
    if (r < 1.0e-10) return 0.0;
    return erf(kappa_E * r) / (r * r);
}

/* Prefactor de Coulomb: e²/(4πε) */
static double shortrange_coulomb_prefactor(psi_t* psi) {
    double eunit, epsilon;
    psi_unit_charge(psi, &eunit);
    psi_epsilon(psi, &epsilon);
    return (eunit * eunit) / (4.0 * PI * epsilon);
}

/* Fuerza de corrección vectorial */
static void shortrange_force(double dist, double q1, double q2,
                              double prefactor, double kappa_E,
                              const double r_vec[3], double F[3]) {
    double kernel = shortrange_erf_kernel(dist, kappa_E);
    double mag = prefactor * q1 * q2 * kernel / dist;  /* F/r para obtener dirección */
    F[X] = mag * r_vec[X];
    F[Y] = mag * r_vec[Y];
    F[Z] = mag * r_vec[Z];
}
```

**Validación**:
- [ ] Compila sin errores
- [ ] `erf(2.0 * 1.5) ≈ 0.9999` (corrección casi completa en rc)

---

### Nivel 1: Solo Self-Field (Alternativa A)

**Corresponde a**: Alternativa A - Solo Corrección de Self-Field

**Qué corrige**: El campo que la partícula genera sobre sí misma via su carga distribuida.

**Dificultad**: ⭐ (Más simple)

**Precisión**: Media - Corrige el efecto dominante pero ignora PP y PN.

```c
/* ========== NIVEL 1: SELF-FIELD ========== */

int subgrid_shortrange_level1(colloids_info_t* cinfo, psi_t* psi) {
    /* Alternativa A: Solo corregir self-field */

    double prefactor = shortrange_coulomb_prefactor(psi);
    double kappa_E = SHORTRANGE_KAPPA_E;
    colloid_t* pc = NULL;

    colloids_info_all_head(cinfo, &pc);

    for (; pc; pc = pc->nextall) {
        if (pc->s.bc != COLLOID_BC_SUBGRID) continue;

        double E_self[3] = {0.0, 0.0, 0.0};
        double q_p = pc->s.q0 - pc->s.q1;

        /* Loop sobre nodos en soporte de Peskin */
        int i0 = (int) floor(pc->s.r[X]);
        int j0 = (int) floor(pc->s.r[Y]);
        int k0 = (int) floor(pc->s.r[Z]);

        for (int di = -1; di <= 2; di++) {
          for (int dj = -1; dj <= 2; dj++) {
            for (int dk = -1; dk <= 2; dk++) {

                double r_pn[3];
                r_pn[X] = (double)(i0 + di) - pc->s.r[X];
                r_pn[Y] = (double)(j0 + dj) - pc->s.r[Y];
                r_pn[Z] = (double)(k0 + dk) - pc->s.r[Z];

                double dist = sqrt(r_pn[X]*r_pn[X] + r_pn[Y]*r_pn[Y] + r_pn[Z]*r_pn[Z]);
                if (dist < 0.1) continue;

                /* Carga de partícula asignada a este nodo */
                double dr = d_peskin(r_pn[X]) * d_peskin(r_pn[Y]) * d_peskin(r_pn[Z]);
                double q_node = q_p * dr;

                /* Campo de corrección */
                double kernel = shortrange_erf_kernel(dist, kappa_E);
                double E_mag = prefactor * q_node * kernel;

                E_self[X] += E_mag * r_pn[X] / dist;
                E_self[Y] += E_mag * r_pn[Y] / dist;
                E_self[Z] += E_mag * r_pn[Z] / dist;
            }
          }
        }

        /* Restar self-field */
        pc->s.Esub[X] -= E_self[X];
        pc->s.Esub[Y] -= E_self[Y];
        pc->s.Esub[Z] -= E_self[Z];
    }

    return 0;
}
```

**Integración** (en ludwig.c después de `subgrid_update_Esub`):
```c
subgrid_shortrange_level1(ludwig->cinfo, ludwig->psi);
```

**Validación Nivel 1**:
- [ ] Partícula aislada sin campo externo: F → 0
- [ ] No introduce velocidades espurias nuevas
- [ ] Partícula en campo uniforme: F = qE exacto

**Cuándo pasar al Nivel 2**: Cuando Nivel 1 funcione pero necesites mejor precisión
en sistemas con múltiples partículas cercanas.

---

### Nivel 2: Self-Field + PP (Alternativa B)

**Corresponde a**: Alternativa B - PP + Self (estilo Ewald)

**Qué agrega**: Corrección partícula-partícula para pares cercanos.

**Dificultad**: ⭐⭐

**Precisión**: Media-Alta - Corrige interacciones entre partículas.

```c
/* ========== NIVEL 2: SELF-FIELD + PP ========== */

int subgrid_shortrange_level2(colloids_info_t* cinfo, psi_t* psi) {
    /* Alternativa B: Self-field + PP */

    /* Primero: corrección self-field (mismo que Nivel 1) */
    subgrid_shortrange_level1(cinfo, psi);

    /* Segundo: corrección PP */
    double prefactor = shortrange_coulomb_prefactor(psi);
    double kappa_E = SHORTRANGE_KAPPA_E;
    double rc = SHORTRANGE_RC;

    int ncell[3];
    colloids_info_ncell(cinfo, ncell);

    /* Loop sobre celdas (patrón de ewald_real_space_sum) */
    for (int ic = 1; ic <= ncell[X]; ic++) {
      for (int jc = 1; jc <= ncell[Y]; jc++) {
        for (int kc = 1; kc <= ncell[Z]; kc++) {

            colloid_t* pc1 = NULL;
            colloids_info_cell_list_head(cinfo, ic, jc, kc, &pc1);

            for (; pc1; pc1 = pc1->next) {
                if (pc1->s.bc != COLLOID_BC_SUBGRID) continue;
                double q1 = pc1->s.q0 - pc1->s.q1;

                /* Celdas vecinas */
                for (int dx = -1; dx <= +1; dx++) {
                  for (int dy = -1; dy <= +1; dy++) {
                    for (int dz = -1; dz <= +1; dz++) {

                        colloid_t* pc2 = NULL;
                        colloids_info_cell_list_head(cinfo, ic+dx, jc+dy, kc+dz, &pc2);

                        for (; pc2; pc2 = pc2->next) {
                            if (pc2->s.bc != COLLOID_BC_SUBGRID) continue;
                            if (pc1->s.index >= pc2->s.index) continue;

                            double q2 = pc2->s.q0 - pc2->s.q1;

                            double r12[3];
                            cs_minimum_distance(cinfo->cs, pc1->s.r, pc2->s.r, r12);
                            double dist = sqrt(r12[X]*r12[X] + r12[Y]*r12[Y] + r12[Z]*r12[Z]);

                            if (dist >= rc || dist < 0.1) continue;

                            /* Fuerza de corrección PP */
                            double F_corr[3];
                            shortrange_force(dist, q1, q2, prefactor, kappa_E, r12, F_corr);

                            /* Newton III entre partículas */
                            pc1->fex[X] += F_corr[X];  pc1->fex[Y] += F_corr[Y];  pc1->fex[Z] += F_corr[Z];
                            pc2->fex[X] -= F_corr[X];  pc2->fex[Y] -= F_corr[Y];  pc2->fex[Z] -= F_corr[Z];
                        }
                    }
                  }
                }
            }
        }
      }
    }

    return 0;
}
```

**Validación Nivel 2**:
- [ ] Todo de Nivel 1 sigue funcionando
- [ ] Dos partículas: F ∝ 1/r² a corta distancia
- [ ] Newton III: Σ F_partículas = 0

**Cuándo pasar al Nivel 3**: Cuando Nivel 2 funcione pero las velocidades espurias
sigan siendo problemáticas (indica que la interacción con iones es importante).

---

### Nivel 3: Self-Field + PP + PN (Alternativa F)

**Corresponde a**: Alternativa F - Enfoque Pragmático Híbrido

**Qué agrega**: Corrección partícula-nodo con iones. **Este es el nivel crítico
para eliminar velocidades espurias.**

**Dificultad**: ⭐⭐⭐

**Precisión**: Alta - Corrige las tres fuentes de error principales.

```c
/* ========== NIVEL 3: SELF-FIELD + PP + PN ========== */

int subgrid_shortrange_level3(colloids_info_t* cinfo, psi_t* psi, hydro_t* hydro) {
    /* Alternativa F: Self-field + PP + PN */

    /* Primero: correcciones de Nivel 2 */
    subgrid_shortrange_level2(cinfo, psi);

    /* Tercero: corrección PN (partícula-nodo) */
    double prefactor = shortrange_coulomb_prefactor(psi);
    double kappa_E = SHORTRANGE_KAPPA_E;
    double rc = SHORTRANGE_RC;
    int range = (int) ceil(rc) + 1;

    colloid_t* pc = NULL;
    colloids_info_all_head(cinfo, &pc);

    for (; pc; pc = pc->nextall) {
        if (pc->s.bc != COLLOID_BC_SUBGRID) continue;

        double q_p = pc->s.q0 - pc->s.q1;
        int i0 = (int) floor(pc->s.r[X]);
        int j0 = (int) floor(pc->s.r[Y]);
        int k0 = (int) floor(pc->s.r[Z]);

        for (int di = -range; di <= range; di++) {
          for (int dj = -range; dj <= range; dj++) {
            for (int dk = -range; dk <= range; dk++) {

                int i = i0 + di, j = j0 + dj, k = k0 + dk;

                double r_pn[3];
                r_pn[X] = (double)i - pc->s.r[X];
                r_pn[Y] = (double)j - pc->s.r[Y];
                r_pn[Z] = (double)k - pc->s.r[Z];

                double dist = sqrt(r_pn[X]*r_pn[X] + r_pn[Y]*r_pn[Y] + r_pn[Z]*r_pn[Z]);
                if (dist >= rc || dist < 0.1) continue;

                /* Carga iónica neta en el nodo */
                int index = cs_index(cinfo->cs, i, j, k);
                double rho0, rho1;
                psi_rho(psi, index, 0, &rho0);
                psi_rho(psi, index, 1, &rho1);
                double Q_node = rho0 - rho1;  /* × Δx³ = 1 */

                if (fabs(Q_node) < 1.0e-15) continue;

                /* Fuerza de corrección PN */
                double F_corr[3];
                shortrange_force(dist, q_p, Q_node, prefactor, kappa_E, r_pn, F_corr);

                /* Fuerza sobre partícula */
                pc->fex[X] += F_corr[X];
                pc->fex[Y] += F_corr[Y];
                pc->fex[Z] += F_corr[Z];

                /* CRÍTICO: Newton III - reacción DIRECTA al nodo k */
                double F_react[3] = {-F_corr[X], -F_corr[Y], -F_corr[Z]};
                hydro_f_local_add(hydro, index, F_react);
            }
          }
        }
    }

    return 0;
}
```

**Validación Nivel 3**:
- [ ] Todo de Niveles 1-2 sigue funcionando
- [ ] Conservación: Σ F_partículas + Σ F_fluido = 0
- [ ] **Velocidades espurias reducidas significativamente**

**Cuándo pasar al Nivel 4**: Cuando Nivel 3 funcione pero necesites mayor
precisión en la función de referencia del mesh.

---

### Nivel 4: Con Función de Referencia Tabulada (Alternativa D)

**Corresponde a**: Alternativa D - P³M Completo

**Qué agrega**: Función de referencia precalculada que representa exactamente
lo que el mesh calcula, en lugar de aproximación analítica.

**Dificultad**: ⭐⭐⭐⭐

**Precisión**: Muy Alta - Errores bien caracterizados.

```c
/* ========== NIVEL 4: P³M CON FUNCIÓN DE REFERENCIA ========== */

/* Tabla precalculada de la función de referencia */
typedef struct {
    int npoints;
    double dr;
    double* F_ref;  /* F_reference(r) tabulada */
} shortrange_table_t;

/* Crear tabla (llamar una vez al inicio) */
int shortrange_table_create(shortrange_table_t** table, double rc, int npoints) {

    *table = (shortrange_table_t*) malloc(sizeof(shortrange_table_t));
    (*table)->npoints = npoints;
    (*table)->dr = rc / npoints;
    (*table)->F_ref = (double*) malloc(npoints * sizeof(double));

    /* Calcular F_reference = convolución de Peskin × Peskin × Coulomb */
    /* Esto representa lo que el mesh "ve" de una carga puntual */
    for (int i = 0; i < npoints; i++) {
        double r = (i + 0.5) * (*table)->dr;
        /* Aproximación: usar la convolución numérica de Peskin */
        (*table)->F_ref[i] = compute_peskin_convolution_force(r);
    }

    return 0;
}

/* Interpolar de la tabla */
static double shortrange_table_lookup(shortrange_table_t* table, double r) {
    int i = (int)(r / table->dr);
    if (i >= table->npoints - 1) return 0.0;
    double t = (r - i * table->dr) / table->dr;
    return (1.0 - t) * table->F_ref[i] + t * table->F_ref[i+1];
}

int subgrid_shortrange_level4(colloids_info_t* cinfo, psi_t* psi,
                               hydro_t* hydro, shortrange_table_t* table) {
    /* Alternativa D: P³M con función de referencia */

    double prefactor = shortrange_coulomb_prefactor(psi);
    double rc = SHORTRANGE_RC;

    /* ... similar a Nivel 3, pero usando: */

    /* F_corr = F_Coulomb_real - F_reference_tabulada */
    double F_coulomb = prefactor * q1 * q2 / (dist * dist);
    double F_ref = shortrange_table_lookup(table, dist);
    double F_corr_mag = F_coulomb - F_ref;

    /* ... resto igual que Nivel 3 ... */

    return 0;
}
```

**Función auxiliar para calcular la convolución de Peskin**:
```c
/* Calcula numéricamente: F_ref(r) = ∫∫ d_peskin(r') d_peskin(r'') F_Coulomb(r-r'-r'') dr' dr'' */
static double compute_peskin_convolution_force(double r) {
    /* Integración numérica en 3D - costosa pero se hace una sola vez */
    double sum = 0.0;
    double h = 0.1;  /* Paso de integración */

    for (double x1 = -2; x1 <= 2; x1 += h) {
      for (double y1 = -2; y1 <= 2; y1 += h) {
        for (double z1 = -2; z1 <= 2; z1 += h) {
          double w1 = d_peskin(x1) * d_peskin(y1) * d_peskin(z1);
          if (w1 < 1e-10) continue;

          for (double x2 = -2; x2 <= 2; x2 += h) {
            for (double y2 = -2; y2 <= 2; y2 += h) {
              for (double z2 = -2; z2 <= 2; z2 += h) {
                double w2 = d_peskin(x2) * d_peskin(y2) * d_peskin(z2);
                if (w2 < 1e-10) continue;

                double dx = r - x1 - x2;  /* Asumiendo r en dirección x */
                double dy = -y1 - y2;
                double dz = -z1 - z2;
                double d = sqrt(dx*dx + dy*dy + dz*dz);

                if (d > 0.1) {
                    sum += w1 * w2 * (1.0 / (d * d)) * h * h * h * h * h * h;
                }
              }
            }
          }
        }
      }
    }

    return sum;
}
```

**Validación Nivel 4**:
- [ ] Todo de Niveles 1-3 sigue funcionando
- [ ] F_total converge a Coulomb analítico para r → 0
- [ ] Error sistemático reducido vs Nivel 3

**Cuándo pasar al Nivel 5**: Cuando Nivel 4 funcione pero necesites conservación
estricta de masa iónica además de momento.

---

### Nivel 5: Conservación Completa (Alternativa E)

**Corresponde a**: Alternativa E - Híbrido con Conservación

**Qué agrega**: Modificación de flujos en Nernst-Planck para conservar masa iónica.

**Dificultad**: ⭐⭐⭐⭐⭐ (Máxima)

**Precisión**: Máxima - Conservación estricta de momento Y masa.

```c
/* ========== NIVEL 5: CONSERVACIÓN COMPLETA ========== */

/* ADVERTENCIA: Este nivel requiere modificar nernst_planck.c */

/* La idea es agregar un término de flujo correctivo en NP:
 *
 * j_corr = (D/kT) × ρ_ion × F_PN_reaction
 *
 * Esto asegura que la fuerza de reacción sobre los iones
 * se traduzca en un flujo que conserva masa.
 */

/* Estructura para almacenar fuerzas de reacción por nodo */
typedef struct {
    int nsites;
    double* F_reaction;  /* [nsites][3] */
} shortrange_reaction_t;

int subgrid_shortrange_level5_forces(colloids_info_t* cinfo, psi_t* psi,
                                      hydro_t* hydro,
                                      shortrange_reaction_t* reaction) {
    /* Calcular fuerzas como en Nivel 4, pero guardar F_reaction por nodo */

    /* ... similar a Nivel 4 ... */

    /* En lugar de hydro_f_local_add, guardar para NP: */
    reaction->F_reaction[index*3 + X] += F_react[X];
    reaction->F_reaction[index*3 + Y] += F_react[Y];
    reaction->F_reaction[index*3 + Z] += F_react[Z];

    /* También agregar al fluido para conservar momento total */
    hydro_f_local_add(hydro, index, F_react);

    return 0;
}

/* Modificación en nernst_planck.c - agregar después de calcular flujos estándar */
int nernst_planck_shortrange_flux_correction(psi_t* psi,
                                              shortrange_reaction_t* reaction,
                                              double dt) {
    /* Agregar flujo correctivo: j_corr = (D/kT) × ρ × F_reaction */

    double diff[2], beta;
    psi_diffusivity(psi, diff);
    psi_beta(psi, &beta);

    for (int index = 0; index < reaction->nsites; index++) {
        double rho0, rho1;
        psi_rho(psi, index, 0, &rho0);
        psi_rho(psi, index, 1, &rho1);

        /* Flujo correctivo para cada especie */
        /* NOTA: Esto es aproximado - el esquema de Capuani usa flujos en enlaces */
        double mobility0 = diff[0] * beta;
        double mobility1 = diff[1] * beta;

        /* El flujo debería integrarse en los enlaces del esquema de Capuani */
        /* Esto requiere modificación más profunda de nernst_planck_fluxes_d3qx */
    }

    return 0;
}
```

**NOTA IMPORTANTE sobre Nivel 5**:

La implementación completa de conservación de masa iónica requiere:

1. Modificar `nernst_planck_fluxes_d3qx()` para incluir el flujo correctivo
2. Asegurar que el flujo correctivo sea simétrico en los enlaces (como Capuani)
3. Verificar que no rompa las propiedades de equilibrio del esquema

**Esto es complejo y puede afectar la estabilidad del esquema de Capuani.**

**Validación Nivel 5**:
- [ ] Todo de Niveles 1-4 sigue funcionando
- [ ] Conservación de momento: Σp = constante
- [ ] Conservación de masa iónica: Σρ_i = constante
- [ ] El esquema sigue siendo estable

---

### Resumen: Progresión de Niveles

| Nivel | Alternativa | Corrige | Conserva | Precisión | Complejidad |
|-------|-------------|---------|----------|-----------|-------------|
| 1 | A | Self-field | - | Media | ⭐ |
| 2 | B | + PP | Momento (PP) | Media-Alta | ⭐⭐ |
| 3 | F | + PN | + Momento (PN) | Alta | ⭐⭐⭐ |
| 4 | D | + Func. Ref. | + Momento | Muy Alta | ⭐⭐⭐⭐ |
| 5 | E | + Flujo NP | + Masa iónica | Máxima | ⭐⭐⭐⭐⭐ |

**Recomendación de uso**:

1. **Empezar siempre por Nivel 1** - Es el más simple y corrige el efecto dominante

2. **Nivel 3 es el objetivo práctico** - Balancea precisión y complejidad

3. **Nivel 4 solo si** la función de referencia analítica no es suficiente

4. **Nivel 5 solo si** la conservación de masa iónica es crítica para tu aplicación

---

### Selector de Nivel en Tiempo de Ejecución

Para facilitar pruebas y comparaciones:

```c
/* En subgrid.h */
typedef enum {
    SHORTRANGE_LEVEL_NONE = 0,
    SHORTRANGE_LEVEL_1_SELF = 1,
    SHORTRANGE_LEVEL_2_PP = 2,
    SHORTRANGE_LEVEL_3_PN = 3,
    SHORTRANGE_LEVEL_4_P3M = 4,
    SHORTRANGE_LEVEL_5_CONSERV = 5
} shortrange_level_t;

/* En subgrid.c */
int subgrid_shortrange_corrections(colloids_info_t* cinfo, psi_t* psi,
                                    hydro_t* hydro, shortrange_level_t level) {
    switch (level) {
        case SHORTRANGE_LEVEL_NONE:
            return 0;
        case SHORTRANGE_LEVEL_1_SELF:
            return subgrid_shortrange_level1(cinfo, psi);
        case SHORTRANGE_LEVEL_2_PP:
            return subgrid_shortrange_level2(cinfo, psi);
        case SHORTRANGE_LEVEL_3_PN:
            return subgrid_shortrange_level3(cinfo, psi, hydro);
        case SHORTRANGE_LEVEL_4_P3M:
            return subgrid_shortrange_level4(cinfo, psi, hydro, table);
        case SHORTRANGE_LEVEL_5_CONSERV:
            return subgrid_shortrange_level5(cinfo, psi, hydro, reaction);
        default:
            return -1;
    }
}

/* En input file:
 * shortrange_correction_level  3
 */
```

Esto permite cambiar entre niveles sin recompilar, facilitando comparaciones
y debugging.

---

## 8. Recomendación

Para una primera implementación, se recomienda la **Alternativa F (Enfoque Pragmático)**
por las siguientes razones:

1. **Balance complejidad/beneficio**: Corrige los efectos más importantes (self-field y PP)
   sin requerir modificaciones extensas.

2. **Reutilización de código**: Aprovecha la estructura de celdas y patrones del código
   de Ewald existente.

3. **Conservación de momento**: Garantiza conservación del momento total al distribuir
   la fuerza de reacción al fluido.

4. **Extensibilidad**: Puede mejorarse incrementalmente hacia las alternativas D o E
   si se requiere mayor precisión.

5. **Validación**: Es más fácil de validar comparando con casos analíticos conocidos.

### Pasos de implementación sugeridos:

1. **Funciones base**:
   - Implementar `coulomb_force()` (Coulomb puro)
   - Implementar `shortrange_correction_force()` con kernel erfc(κ_E r)
   - Elegir κ_E ≈ 2/Δx

2. **Corrección de self-field**:
   - Calcular E_self usando el kernel de corrección
   - Validar que E_self → 0 cuando partícula está en centro de nodo

3. **Correcciones PP (partícula-partícula)**:
   - Usar estructura de celdas existente (como ewald_real_space_sum)
   - Aplicar Newton III: F sobre p1, -F sobre p2

4. **Correcciones PN (partícula-nodo)**:
   - Tratar nodos como "partículas" de carga Q_nodo = ρ_net × Δx³
   - **CRÍTICO**: Reacción directa al nodo, NO distribuir con Peskin

5. **Validación de conservación**:
   - Verificar que Σ F = 0 (momento total)
   - Monitorear velocidades espurias

6. **Comparación con soluciones analíticas**:
   - Caso de una partícula aislada (self-field debe ser cero)
   - Dos partículas (recuperar ley de Coulomb)

### Puntos Críticos a Recordar:

| Aspecto | Correcto | Incorrecto |
|---------|----------|------------|
| Potencial físico | Coulomb 1/r | Yukawa exp(-κr)/r |
| Kernel de splitting | erfc(κ_E r) con κ_E numérico | exp(-κ_DH r) con κ físico |
| Reacción PN | Directa al nodo k | Distribuida con Peskin |
| κ_E | ≈ 2/Δx (numérico) | κ_DH de Debye (físico) |

---

## 8. Referencias

1. Hockney, R. W., & Eastwood, J. W. (1988). Computer Simulation Using Particles.
   CRC Press. Capítulo 8: The P³M Algorithm.

2. Nash, R. W., Adhikari, R., & Cates, M. E. (2008). Singular forces and pointlike
   colloids in lattice Boltzmann hydrodynamics. Physical Review E, 77(2), 026709.

3. Capuani, F., Pagonabarraga, I., & Frenkel, D. (2004). Discrete solution of the
   electrokinetic equations. The Journal of Chemical Physics, 121(2), 973-986.

4. Peskin, C. S. (2002). The immersed boundary method. Acta Numerica, 11, 479-517.

---

## Apéndice A: Código de Referencia - Ewald Real Space Sum

El siguiente código de `src/ewald.c` muestra el patrón para iterar sobre pares de
partículas cercanas:

```c
int ewald_real_space_sum(ewald_t * ewald) {

    colloids_info_ncell(ewald->cinfo, ncell);

    for (ic = 1; ic <= ncell[X]; ic++) {
        for (jc = 1; jc <= ncell[Y]; jc++) {
            for (kc = 1; kc <= ncell[Z]; kc++) {

                colloids_info_cell_list_head(ewald->cinfo, ic, jc, kc, &p_c1);

                for (; p_c1; p_c1 = p_c1->next) {

                    for (dx = -1; dx <= +1; dx++) {
                        for (dy = -1; dy <= +1; dy++) {
                            for (dz = -1; dz <= +1; dz++) {

                                colloids_info_cell_list_head(ewald->cinfo,
                                    ic+dx, jc+dy, kc+dz, &p_c2);

                                for (; p_c2; p_c2 = p_c2->next) {

                                    if (p_c1->s.index < p_c2->s.index) {
                                        // Calcular interacción...
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
```

---

## Apéndice B: Función Delta de Peskin

La función de Peskin implementada en `src/subgrid.c`:

```c
double d_peskin(double r) {
    double rmod = fabs(r);
    double delta = 0.0;

    if (rmod <= 1.0) {
        delta = 0.125 * (3.0 - 2.0*rmod + sqrt(1.0 + 4.0*rmod - 4.0*rmod*rmod));
    }
    else if (rmod <= 2.0) {
        delta = 0.125 * (5.0 - 2.0*rmod - sqrt(-7.0 + 12.0*rmod - 4.0*rmod*rmod));
    }

    return delta;
}
```

Propiedades:
- Soporte compacto: d(r) = 0 para |r| > 2
- Normalización: ∫ d(r) dr = 1
- Suavidad: C¹ continua
- Partición de unidad: Σ_n d(r-n) = 1 para todo r

