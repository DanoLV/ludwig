# Potencial de Exclusión de Volumen para Iones en Ludwig

## Problema Físico Actual

La redistribución de carga removida mediante el kernel de Peskin es **no-física**:

- Se remueve carga de nodos dentro del coloide
- Se redistribuye a través de un kernel Peskin 4³, esparciendo la carga por toda una región
- **Incorrecto:** Los iones solo deberían moverse **localmente** hacia nodos adyacentes
- **Dinámica violada:** No respeta la física de difusión/transporte electrocinético

## Solución: Potencial de Exclusión Suave

En lugar de remover y redistribuir discretamente, agregar un **potencial externo** $V_{ext}(\mathbf{r})$ que:

1. Actúa en la ecuación de Nernst-Planck como término en el potencial químico
2. Los iones se expulsan **naturalmente** a través de la dinámica de difusión
3. La redistribución ocurre **localmente** y es **físicamente consistente**

### Relación con la Accesibilidad S(r)

La función de accesibilidad cúbica $S(\mathbf{r})$ ya existe en el código:
- Define suavemente la transición sólido ↔ fluido en [a, a+ε]
- Es un **campo de fase geométrico**

La relación correcta es:
$$V_{ext}(\mathbf{r} - \mathbf{r}_p) = V_0(1 - S(\mathbf{r} - \mathbf{r}_p))$$

Esto convierte la geometría (S) en física (energía libre).

## Implementación en Nernst-Planck

### 1. Potencial Químico Modificado

En `nernst_planck.c`, el potencial químico de Capuani es:
$$\mu_k = kT \ln(\rho_k) + z_k e \psi$$

Debe extenderse a:
$$\mu_k = kT \ln(\rho_k) + z_k e \psi + V_0(1 - S(\mathbf{r} - \mathbf{r}_p))$$

**Clave:** El término $V_{ext}$ entra en el `exp[β μ^{ex}]` del flujo de Capuani.

### 2. Flujo de Iones

La ecuación de NP que ya usa Ludwig:
$$j_k = -D_k \rho_k \nabla(\beta \mu_k^{ex})$$

Con $V_{ext}$ incluido, automáticamente genera:
$$j_k = -D_k \rho_k \nabla(\beta(z_k e \psi + V_0(1-S)))$$

Los gradientes de $S$ actúan como **repulsión suave**.

## Energía Libre del Sistema

La energía libre total (consistente con Ludwig):

$$F = \int \sum_i \rho_i \left[ kT (\ln \rho_i - 1) + z_i e \psi + V_0(1-S(\mathbf{r} - \mathbf{r}_p)) \right] d\mathbf{r}$$

De esta energía derivan automáticamente:
- **Potencial químico:** $\mu_i = \frac{\delta F}{\delta \rho_i}$
- **Fuerza en el fluido:** $\mathbf{f}_{fluid} = -\sum_i \rho_i \nabla \mu_i$
- **Fuerza en la partícula:** $\mathbf{F}_p = -\frac{\partial F}{\partial \mathbf{r}_p}$

## Fuerza en la Partícula

### Forma General

$$\mathbf{F}_p = -\frac{\partial F}{\partial \mathbf{r}_p} = -V_0 \int \sum_i \rho_i(\mathbf{r}) \nabla_{\mathbf{r}} S(\mathbf{r} - \mathbf{r}_p) d\mathbf{r}$$

**Explícitamente en función de la partícula** (no del fluido):

Si $S$ es radial: $S = S(r)$ donde $r = |\mathbf{r} - \mathbf{r}_p|$

$$\mathbf{F}_p = -V_0 \int \sum_i \rho_i(\mathbf{r}) \frac{dS}{dr}\bigg|_{r=|\mathbf{r}-\mathbf{r}_p|} \frac{\mathbf{r} - \mathbf{r}_p}{|\mathbf{r} - \mathbf{r}_p|} d\mathbf{r}$$

### Interpretación Física

- La fuerza viene únicamente de la interfase (donde $\frac{dS}{dr} \neq 0$)
- Es equivalente a la **presión osmótica integrada sobre la capa de exclusión**
- Completamente física y necesaria para conservación de momento

### Implementación Discreta

```c
// Para cada partícula p
Fp = 0;
for (cada nodo i) {
    r_vec = node_pos[i] - particle_pos_p;
    r = |r_vec|;

    // Gradiente de S en dirección radial
    dS_dr = dS_dr_function(r, R, delta);

    // Densidad total de iones
    rho_tot = rho_plus[i] + rho_minus[i];

    // Contribución a la fuerza
    Fp += rho_tot * dS_dr * (r_vec / r) * dV;
}

Fp *= -V0;
```

**Clave:** La fuerza es estrictamente en función de $\mathbf{r}_p$ a través de $S(\mathbf{r} - \mathbf{r}_p)$.

## Caso de Múltiples Partículas

Para evitar acumulación infinita de potencial cuando partículas se acercan:

$$V_{ext}^{tot}(\mathbf{r}) = V_0 \left(1 - \prod_p S_p(\mathbf{r} - \mathbf{r}_p)\right)$$

**Interpretación:** Una única región sólida efectiva, aunque haya múltiples partículas.

### Fuerza con Múltiples Partículas

Para la partícula $p$:
$$\mathbf{F}_p = -V_0 \int \sum_i \rho_i(\mathbf{r}) \nabla_{\mathbf{r}_p} \prod_q S_q(\mathbf{r} - \mathbf{r}_q) d\mathbf{r}$$

Solo derivás respecto a ESA partícula.

## Parámetros

### Magnitud del Potencial: $V_0$

Controla cuán fuerte es la exclusión. En equilibrio:
$$\rho \propto e^{-\beta V_{ext}}$$

**Valores típicos:**
- $\beta V_0 \approx 5-8$: exclusión suave (~ρ reduce ×150)
- $\beta V_0 \approx 8-12$: exclusión fuerte (ρ ≈ 0)
- $\beta V_0 > 20$: inestable / gradientes explosivos

**Recomendación:** $\beta V_0 = 8-12$

### Validación

Criterio: $\rho_{inside} / \rho_{bulk} < 10^{-3}$

Si se cumple → $V_0$ es suficiente

### Ancho de Interfaz: $\delta$

Ya está definido en tu $S(r)$ actual.

**Rango estable:** $1 < \delta < 3$ Δx

**Óptimo:** $\delta \approx 1.5$ Δx

## Errores a Evitar

❌ **NO hacer:**

1. Multiplicar $D$ por $S(r)$ → rompe estructura variacional y equilibrio
2. Multiplicar $\rho$ directamente por $S(r)$ → viola conservación
3. Usar LJ directo → singular, incompatible con PNP continuo
4. Agregar la fuerza como $\mathbf{f} = -\rho \nabla V_{ext}$ sin derivarla de energía libre

✅ **SÍ hacer:**

1. Introducir $V_{ext}$ en $\mu_k$ (potencial químico)
2. Dejar que Capuani genere los flujos naturalmente
3. Derivar fuerzas desde energía libre (variacional)
4. Verificar conservación de momento: $\mathbf{F}_p = - \int \mathbf{f}_{fluid}$

## Extensiones: Modelos Más Realistas

Si querés ir más allá de $V_{ext}$ ad-hoc, podés usar:

### Modelo de Bikerman (Lattice Gas)

Exclusión de volumen automática con:
$$\mu_i^{ex} = kT \ln \left(\frac{1}{1 - v \sum_j \rho_j}\right)$$

donde $v$ es el volumen de un ion.

**Ventaja:** Basado en energía libre, sin parámetros arbitrarios.

### Modelo de Carnahan-Starling (Hard-Sphere)

Más preciso para altas concentraciones:
$$\mu_i^{ex} = kT\frac{\eta(8 - 9\eta + 3\eta^2)}{(1 - \eta)^3}$$

donde $\eta = v \sum_i \rho_i$.

## Resumen de la Implementación

| Paso | Qué | Dónde |
|------|-----|-------|
| 1 | Calcular $S(\mathbf{r} - \mathbf{r}_p)$ | `psi_exclusion.c` (ya existe) |
| 2 | Calcular $V_{ext} = V_0(1-S)$ | `nernst_planck.c` |
| 3 | Agregar a $\mu_k$ en cálculo de flux | `nernst_planck.c` (línea exp[β μ^{ex}]) |
| 4 | Calcular fuerza en partícula | Nueva función en `colloids.c` o loop de `ludwig.c` |
| 5 | Aplicar fuerza a partícula | Loop principal de dinámica |
| 6 | Verificar conservación | Tests de momento y carga |

## Referencias Conceptuales

- **Capuani et al. (2004):** Estructura variacional de NP en redes
- **Peskin (1977):** Interpolación suave (IB method)
- **Bikerman (1942):** Exclusión de volumen en electrolitos
- **Carnahan-Starling (1969):** EOS de esferas duras
- **Ludwig código:** Implementación de Capuani con Poisson solvido externamente

---

**Autor:** Discusión Claude-ChatGPT sobre correcciones a redistribución de carga en sistemas coloide-electrolito

**Fecha:** 2026-03-30

**Estado:** Conceptualmente completado, listo para implementación
