# Método Lattice Boltzmann con Multiple Relaxation Times (MRT)

## Autor: Documentación generada desde el código Ludwig
## Fecha: 2025-11-02

---

## Índice

1. [Introducción al método Lattice Boltzmann](#1-introducción-al-método-lattice-boltzmann)
2. [Fundamentos matemáticos](#2-fundamentos-matemáticos)
3. [El modelo D3Q19](#3-el-modelo-d3q19)
4. [Modos hidrodinámicos](#4-modos-hidrodinámicos)
5. [Paso de colisión MRT](#5-paso-de-colisión-mrt)
6. [Paso de propagación](#6-paso-de-propagación)
7. [Implementación en Ludwig](#7-implementación-en-ludwig)
8. [Referencias en el código](#8-referencias-en-el-código)

---

## 1. Introducción al método Lattice Boltzmann

El **método Lattice Boltzmann (LBM)** es un enfoque computacional para simular dinámica de fluidos basado en la teoría cinética de gases. En lugar de resolver directamente las ecuaciones de Navier-Stokes, el LBM modela el fluido como un conjunto de partículas ficticias que residen en los nodos de una red (lattice) y se mueven a lo largo de direcciones discretas.

### Ventajas del LBM

- Facilita el manejo de geometrías complejas
- Paralelización natural
- Inclusión sencilla de física multicomponente y multifase
- Condiciones de frontera relativamente simples

### El algoritmo básico

El LBM consiste en dos pasos fundamentales que se ejecutan alternadamente:

1. **Colisión**: Las distribuciones en cada sitio se relajan hacia el equilibrio
2. **Propagación** (streaming): Las distribuciones post-colisión se mueven a sitios vecinos

---

## 2. Fundamentos matemáticos

### 2.1 Ecuación de Boltzmann en red

La ecuación de evolución del Lattice Boltzmann es:

```
f_i(x + c_i Δt, t + Δt) = f_i(x, t) + Ω_i(x, t)
```

Donde:
- `f_i(x, t)` = función de distribución en la dirección `i`, posición `x` y tiempo `t`
- `c_i` = vector de velocidad discreta en la dirección `i`
- `Ω_i(x, t)` = operador de colisión
- `Δt` = paso temporal (usualmente Δt = 1 en unidades de red)

Esta ecuación se puede separar en:

**Paso de colisión:**
```
f̃_i(x, t) = f_i(x, t) + Ω_i(x, t)
```

**Paso de propagación:**
```
f_i(x + c_i, t + 1) = f̃_i(x, t)
```

### 2.2 Variables macroscópicas

Las cantidades hidrodinámicas se obtienen como momentos de las funciones de distribución:

**Densidad:**
```
ρ(x, t) = Σ_i f_i(x, t)
```

**Momentum:**
```
ρu_α(x, t) = Σ_i f_i(x, t) c_iα
```

**Tensor de esfuerzos:**
```
Π_αβ(x, t) = Σ_i f_i(x, t) c_iα c_iβ
```

### 2.3 Operador de colisión BGK

El operador BGK (Bhatnagar-Gross-Krook) más simple es:

```
Ω_i = -1/τ [f_i - f_i^eq]
```

Donde:
- `τ` = tiempo de relajación único
- `f_i^eq` = distribución de equilibrio

**Problema:** Un solo tiempo de relajación limita la estabilidad numérica.

### 2.4 Operador de colisión MRT

El **Multiple Relaxation Times (MRT)** usa diferentes tiempos de relajación para diferentes modos:

```
Ω = -M^(-1) · S · [m - m^eq]
```

Donde:
- `M` = matriz de transformación a modos
- `S` = matriz diagonal de tiempos de relajación
- `m` = vector de modos
- `m^eq` = vector de modos de equilibrio

**Ventaja:** Permite controlar independientemente la viscosidad de corte (shear) y volumétrica (bulk), mejorando estabilidad.

---

## 3. El modelo D3Q19

### 3.1 Descripción

El modelo **D3Q19** (3 Dimensiones, 19 Velocidades) es uno de los más utilizados para simulaciones 3D.

### 3.2 Vectores de velocidad

Los 19 vectores de velocidad discreta `c_i` son:

```c
// Velocidad 0: Partícula en reposo
c[0]  = { 0,  0,  0}

// Velocidades 1-6: Caras del cubo (peso 2/36)
c[3]  = { 1,  0,  0}
c[7]  = { 0,  1,  0}
c[9]  = { 0,  0,  1}
c[16] = {-1,  0,  0}
c[12] = { 0, -1,  0}
c[10] = { 0,  0, -1}

// Velocidades 7-18: Aristas del cubo (peso 1/36)
c[1]  = { 1,  1,  0}
c[2]  = { 1,  0,  1}
c[4]  = { 1,  0, -1}
c[5]  = { 1, -1,  0}
c[6]  = { 0,  1,  1}
c[8]  = { 0,  1, -1}
c[11] = { 0, -1,  1}
c[13] = { 0, -1, -1}
c[14] = {-1,  1,  0}
c[15] = {-1,  0,  1}
c[17] = {-1,  0, -1}
c[18] = {-1, -1,  0}
```

### 3.3 Pesos (weights)

Los pesos para la distribución de equilibrio son:

```
w[0]  = 12/36    (partícula en reposo)
w[i]  = 2/36     para velocidades en caras (i = 3, 7, 9, 10, 12, 16)
w[i]  = 1/36     para velocidades en aristas (resto)
```

### 3.4 Velocidad del sonido

La velocidad del sonido en la red es:

```
c_s² = 1/3  (en unidades de red)
```

### 3.5 Distribución de equilibrio

Para D3Q19, la distribución de equilibrio es:

```
f_i^eq = w_i ρ [1 + (c_i·u)/c_s² + (c_i·u)²/(2c_s⁴) - u²/(2c_s²)]
```

Donde `u` es la velocidad macroscópica del fluido.

---

## 4. Modos hidrodinámicos

### 4.1 Concepto de modos

En lugar de trabajar directamente con las 19 funciones de distribución `f_i`, el método MRT las transforma a un **espacio de modos** usando una matriz de proyección `M`.

### 4.2 Transformación a modos

```
m = M · f
```

Donde:
- `m` = vector de 19 modos
- `M` = matriz 19×19 de proyección
- `f` = vector de 19 distribuciones

### 4.3 Los 19 modos para D3Q19

Los modos se dividen en:

#### **Modos conservados (3 modos):**

```
m[0] = ρ                    (densidad)
m[1] = ρ u_x                (momentum en x)
m[2] = ρ u_y                (momentum en y)
m[3] = ρ u_z                (momentum en z)
```

Estos modos **NO se relajan** porque representan cantidades conservadas.

#### **Modos de esfuerzo (6 modos):**

```
m[4] = S_xx = ρ u_x² - ρ c_s²     (esfuerzo normal en x)
m[5] = S_xy = ρ u_x u_y           (esfuerzo de corte xy)
m[6] = S_xz = ρ u_x u_z           (esfuerzo de corte xz)
m[7] = S_yy = ρ u_y² - ρ c_s²     (esfuerzo normal en y)
m[8] = S_yz = ρ u_y u_z           (esfuerzo de corte yz)
m[9] = S_zz = ρ u_z² - ρ c_s²     (esfuerzo normal en z)
```

Estos modos determinan la **viscosidad** del fluido.

#### **Modos fantasma (ghost modes, 9 modos):**

```
m[10] = χ1       = (2c² - 3)(3c_z² - c²)
m[11] = j_χ1,x   = χ1 · c_x
m[12] = j_χ1,y   = χ1 · c_y
m[13] = j_χ1,z   = χ1 · c_z
m[14] = χ2       = (2c² - 3)(c_y² - c_x²)
m[15] = j_χ2,x   = χ2 · c_x
m[16] = j_χ2,y   = χ2 · c_y
m[17] = j_χ2,z   = χ2 · c_z
m[18] = χ3       = 3c⁴ - 6c² + 1
```

Estos modos **no tienen significado físico directo** pero afectan la estabilidad numérica.

### 4.4 Matriz de proyección M para D3Q19

La matriz `M` tiene 19 filas, una por cada modo. Cada fila `m` proyecta las 19 distribuciones según:

```c
// Código en lb_d3q19.c:107-150
for (int p = 0; p < 19; p++) {
    double rho  = 1.0;
    double cx   = model->cv[p][X];
    double cy   = model->cv[p][Y];
    double cz   = model->cv[p][Z];
    double sxx  = cx*cx - cs2;  // cs2 = 1/3
    double sxy  = cx*cy;
    double sxz  = cx*cz;
    double syy  = cy*cy - cs2;
    double syz  = cy*cz;
    double szz  = cz*cz - cs2;

    double c2   = cx*cx + cy*cy + cz*cz;
    double chi1 = (2.0*c2 - 3.0)*(3.0*cz*cz - c2);
    double chi2 = (2.0*c2 - 3.0)*(cy*cy - cx*cx);
    double chi3 = 3.0*c2*c2 - 6.0*c2 + 1;

    M[0][p] = rho;
    M[1][p] = cx;
    M[2][p] = cy;
    M[3][p] = cz;
    M[4][p] = sxx;
    M[5][p] = sxy;
    M[6][p] = sxz;
    M[7][p] = syy;
    M[8][p] = syz;
    M[9][p] = szz;
    M[10][p] = chi1;
    M[11][p] = chi1*cx;
    M[12][p] = chi1*cy;
    M[13][p] = chi1*cz;
    M[14][p] = chi2;
    M[15][p] = chi2*cx;
    M[16][p] = chi2*cy;
    M[17][p] = chi2*cz;
    M[18][p] = chi3;
}
```

### 4.5 Transformación inversa

Para volver de modos a distribuciones:

```
f = M^(-1) · m
```

La matriz inversa `M^(-1)` se denota como `mi` en el código.

---

## 5. Paso de colisión MRT

### 5.1 Algoritmo completo de `lb_collision_mrt1_site`

El paso de colisión MRT en Ludwig sigue estos pasos:

#### **PASO 1: Cargar distribuciones** (líneas 316-321)

```c
for (p = 0; p < NVEL; p++) {
    for_simd_v(iv, NSIMDVL) {
        fchunk[p * NSIMDVL + iv] = lb->f[...];
    }
}
```

Lee las funciones de distribución `f_i` del sitio actual.

#### **PASO 2: Transformar f → modos** (líneas 330-343)

```c
#ifdef _D3Q19_
    d3q19_f2mode_chunk(mode, fchunk);
#else
    for (m = 0; m < NVEL; m++) {
        mode[m] = 0.0;
        for (p = 0; p < NVEL; p++) {
            mode[m] += fchunk[p] * M[m][p];
        }
    }
#endif
```

**Ecuación:** `m = M · f`

Proyecta las distribuciones al espacio de modos.

#### **PASO 3: Extraer variables hidrodinámicas** (líneas 345-366)

```c
rho[iv] = mode[0];                          // Densidad
u[0][iv] = mode[1];                         // Momentum x
u[1][iv] = mode[2];                         // Momentum y
u[2][iv] = mode[3];                         // Momentum z
s[ia][ib][iv] = mode[1 + NDIM + m];        // Tensor esfuerzos
```

Extrae las cantidades físicas de los primeros modos.

#### **PASO 4: Calcular velocidad con fuerza** (líneas 368-376)

```c
rrho[iv] = 1.0 / rho[iv];
for (ia = 0; ia < 3; ia++) {
    u[ia][iv] = rrho[iv] * (u[ia][iv] + 0.5 * force[ia][iv]);
}
```

**Ecuación:** `u = (ρu + F/2) / ρ`

Corrige el momentum con la fuerza externa usando el **esquema de Guo**.

#### **PASO 5: Calcular tiempos de relajación** (líneas 378-398)

```c
if (have_visc_model == 0) {
    eta[iv] = eta_shear;           // Viscosidad fija
    eta_bulk[iv] = eta_bulk_const;
} else {
    eta[iv] = hydro->eta[index];   // Viscosidad variable (no-Newtoniano)
    eta_bulk[iv] = (eta_bulk_const / eta_shear) * eta[iv];
}

lb_relaxation_time_shear_v(lb, eta, rtau);
lb_relaxation_time_bulk_v(lb, eta, eta_bulk, rtau_bulk);
lb_relaxation_time_ghosts_v(lb, eta, rtau_ghost);
```

Los tiempos de relajación se relacionan con la viscosidad por:

```
η = ρ c_s² (τ - 0.5) Δt
```

Por tanto:

```
τ = η/(ρ c_s²) + 0.5
```

Y en el código se usa `rtau = 1/τ` (inverso del tiempo de relajación).

#### **PASO 6: Calcular tensor de esfuerzos de equilibrio** (líneas 407-441)

**Caso sin energía libre:**
```c
seq[ia][ib][iv] = rho[iv] * u[ia][iv] * u[ib][iv];
```

**Caso con energía libre (fe):**
```c
fe->func->str_symm_v(fe, index0, symm);  // Tensor simétrico de la energía libre
seq[ia][ib][iv] = rho[iv] * u[ia][iv] * u[ib][iv] + symm[ia][ib][iv];
```

**Ecuación general:**
```
S_αβ^eq = ρ u_α u_β + Π_αβ^fe
```

Donde `Π^fe` es la contribución del estrés químico/electrostático de la energía libre.

#### **PASO 7: Separar traza (presión) de parte traceless** (líneas 443-449)

```c
// Calcular trazas
tr_s[iv] = s[0][0][iv] + s[1][1][iv] + s[2][2][iv];
tr_seq[iv] = seq[0][0][iv] + seq[1][1][iv] + seq[2][2][iv];

// Formar partes traceless
for (ia = 0; ia < 3; ia++) {
    s[ia][ia][iv] -= (1.0/3.0) * tr_s[iv];
    seq[ia][ia][iv] -= (1.0/3.0) * tr_seq[iv];
}
```

**Por qué:** La viscosidad volumétrica actúa solo sobre la traza, mientras que la viscosidad de corte actúa sobre la parte sin traza.

Descomposición:
```
S_αβ = S̄_αβ + (1/3) δ_αβ Tr(S)
```

Donde `S̄_αβ` es la parte traceless (deviatoric).

#### **PASO 8: RELAJACIÓN - Núcleo del MRT** (líneas 451-468)

**Relajación de la traza (viscosidad bulk):**
```c
tr_s[iv] = tr_s[iv] - rtau_bulk[iv] * (tr_s[iv] - tr_seq[iv]);
```

**Ecuación:** `Tr(S)^new = Tr(S) - (1/τ_bulk) [Tr(S) - Tr(S)^eq]`

**Relajación de la parte traceless (viscosidad shear):**
```c
for (ia = 0; ia < 3; ia++) {
    for (ib = 0; ib < 3; ib++) {
        s[ia][ib][iv] -= rtau[iv] * (s[ia][ib][iv] - seq[ia][ib][iv]);
    }
}
```

**Ecuación:** `S̄_αβ^new = S̄_αβ - (1/τ_shear) [S̄_αβ - S̄_αβ^eq]`

**Reconstruir tensor completo:**
```c
s[ia][ib][iv] += delta[ia][ib] * (1.0/3.0) * tr_s[iv];
```

**Corrección por fuerza:**
```c
s[ia][ib][iv] += (2.0 - rtau[iv])
                  * (u[ia][iv]*force[ib][iv] + force[ia][iv]*u[ib][iv]);
```

**Ecuación de corrección:**
```
S_αβ += (2 - 1/τ) (u_α F_β + F_α u_β)
```

Este término asegura que la fuerza externa se incorpore correctamente en el estrés.

#### **PASO 9: Añadir fluctuaciones térmicas (opcional)** (líneas 470-511)

Si `noise` está activado:

```c
var = lb_fluctuations_var_eta(1.0/rtau[iv], kT);
var_bulk = lb_fluctuations_var_bulk(1.0/rtau_bulk[iv], kT);

lb_fluctuations_stress(noise, index, var, var_bulk, shat);
lb_fluctuations_ghosts(noise, index, var_ghost, ghat);
```

**Teorema de fluctuación-disipación:**

La varianza del ruido estocástico `ŝ` en el tensor de esfuerzos es:

```
⟨ŝ_αβ ŝ_γδ⟩ = 2 k_B T η (δ_αγ δ_βδ + δ_αδ δ_βγ) / (Δx³ Δt)
```

Esto simula las fluctuaciones hidrodinámicas presentes en fluidos reales a escala mesoscópica.

#### **PASO 10: Actualizar modos post-colisión** (líneas 513-538)

**Densidad (sin cambio):**
```c
// mode[0] no cambia
```

**Momentum (actualizar con fuerza):**
```c
for (ia = 0; ia < 3; ia++) {
    mode[1 + ia] += force[ia][iv];
}
```

**Tensor de esfuerzos (relajado + ruido):**
```c
m = 0;
for (ia = 0; ia < 3; ia++) {
    for (ib = ia; ib < 3; ib++) {
        mode[1 + 3 + m] = s[ia][ib][iv] + shat[ia][ib][iv];
        m++;
    }
}
```

**Modos fantasma (relajados a cero):**
```c
for (m = NHYDRO; m < NVEL; m++) {  // NHYDRO = 10
    mode[m] = mode[m] - rtau_ghost[m] * (mode[m] - 0.0) + ghat[m];
}
```

**Ecuación:** `m_ghost^new = m_ghost - (1/τ_ghost) m_ghost + ĝ`

Los modos fantasma se relajan hacia **equilibrio cero**.

#### **PASO 11: Transformar modos → distribuciones** (líneas 541-553)

```c
#ifdef _D3Q19_
    d3q19_mode2f_chunk(mode, fchunk);
#else
    for (p = 0; p < NVEL; p++) {
        ftmp[p] = 0.0;
        for (m = 0; m < NVEL; m++) {
            ftmp[p] += M_inv[p][m] * mode[m];
        }
        fchunk[p] = ftmp[p];
    }
#endif
```

**Ecuación:** `f^new = M^(-1) · m^new`

Proyecta los modos post-colisión de vuelta al espacio de distribuciones.

#### **PASO 12: Escribir resultados** (líneas 555-590)

```c
// Guardar distribuciones
lb->f[index, p] = fchunk[p];

// Guardar densidad
hydro->rho[index] = rho[iv];

// Guardar velocidad
hydro->u[index, ia] = u[ia][iv];
```

### 5.2 Resumen del operador de colisión MRT

En forma matricial compacta:

```
f^new = f + M^(-1) · S · (m^eq - m) + términos de fuerza + ruido
```

Donde `S` es la matriz diagonal:

```
S = diag(0, 0, 0, 0, -1/τ_shear, -1/τ_shear, ..., -1/τ_ghost, ...)
```

Los primeros 4 elementos son 0 porque densidad y momentum se conservan.

---

## 6. Paso de propagación

### 6.1 Concepto

Después de la colisión, las distribuciones se mueven a lo largo de sus velocidades discretas:

```
f_i(x + c_i, t+1) = f̃_i(x, t)
```

Este paso es puramente geométrico: cada distribución se copia al sitio vecino en la dirección de su velocidad.

### 6.2 Esquema Pull

Ludwig usa un esquema **pull** en lugar de push:

```c
// Código en propagation.c:124-132
for (int p = 0; p < NVEL; p++) {
    // Pull from neighbour
    icp = ic - cv[p][X];
    jcp = jc - cv[p][Y];
    kcp = kc - cv[p][Z];
    indexp = index(icp, jcp, kcp);

    fprime[index, p] = f[indexp, p];
}
```

**Interpretación:** En lugar de empujar (push) `f[index, p]` al vecino, jalamos (pull) desde el vecino con velocidad `-c_p`.

**Ventaja del pull:**
- Evita condiciones de carrera (race conditions) en paralelización
- Cada hilo/thread solo escribe en su propia celda

### 6.3 Swap de arrays

Después de la propagación, los arrays se intercambian:

```c
// propagation.c:92
lb_model_swapf(lb);
```

Esto intercambia los punteros `f` y `fprime` para evitar copias innecesarias:

```
temp = f;
f = fprime;
fprime = temp;
```

### 6.4 Propagación vectorizada

Para eficiencia, Ludwig usa una versión vectorizada que procesa múltiples sitios simultáneamente (SIMD):

```c
// propagation.c:150+
__global__ void lb_propagation_kernel(kernel_3d_v_t k3v, lb_t * lb)
```

---

## 7. Implementación en Ludwig

### 7.1 Estructura del ciclo temporal

El ciclo principal de simulación en Ludwig es:

```c
// ludwig.c (esquema simplificado)
for (timestep = 0; timestep < nsteps; timestep++) {

    // 1. Calcular fuerzas de energía libre
    fe_force_compute(fe, hydro);

    // 2. COLISIÓN
    lb_collision_mrt(lb, hydro, map, noise, fe);

    // 3. PROPAGACIÓN
    lb_propagation(lb);

    // 4. Aplicar condiciones de frontera
    boundary_conditions_apply(bc, lb);

    // 5. Actualizar otros campos (composición, etc.)
    advection_diffusion_step(phi);

    // 6. Salida de datos
    if (timestep % output_freq == 0) {
        write_output(lb, hydro);
    }
}
```

### 7.2 Jerarquía de llamadas para colisión

```
lb_collide()                              [collision.c:154]
  └─> lb_collision_mrt1()                 [collision.c:221] - Kernel CUDA
       └─> lb_collision_mrt1_site()       [collision.c:254] - Por sitio
            ├─> d3q19_f2mode_chunk()      Transformar f→modos
            ├─> lb_relaxation_time_*()     Calcular τ
            ├─> lb_fluctuations_*()        Añadir ruido
            └─> d3q19_mode2f_chunk()       Transformar modos→f
```

### 7.3 Jerarquía de llamadas para propagación

```
lb_propagation()                          [propagation.c:43]
  └─> lb_propagation_driver()             [propagation.c:58]
       ├─> lb_propagation_kernel()        [propagation.c:150+] - Kernel CUDA
       └─> lb_model_swapf()               Intercambiar f ↔ fprime
```

### 7.4 Optimizaciones clave

#### Vectorización SIMD

El código procesa `NSIMDVL` sitios simultáneamente (típicamente 1, 2, 4 u 8):

```c
#define NSIMDVL 1  // o 2, 4, 8 dependiendo del hardware

for_simd_v(iv, NSIMDVL) {
    // Operaciones vectorizadas
}
```

Esto permite usar instrucciones SIMD del procesador (SSE, AVX, etc.).

#### Kernels CUDA optimizados para D3Q19

Las funciones `d3q19_f2mode_chunk()` y `d3q19_mode2f_chunk()` están **altamente optimizadas** usando:
- Desenrollado de bucles (loop unrolling)
- Instrucciones FMA (Fused Multiply-Add)
- Implementación de Kahan para reducir errores de redondeo

#### Uso de memoria constante en GPU

```c
static __constant__ lb_collide_param_t lbp;
```

Los parámetros de colisión se almacenan en memoria constante de la GPU para acceso rápido.

---

## 8. Referencias en el código

### 8.1 Archivos principales

| Archivo | Descripción |
|---------|-------------|
| [src/collision.c](../src/collision.c) | Implementación de colisión MRT |
| [src/propagation.c](../src/propagation.c) | Implementación de propagación |
| [src/lb_d3q19.c](../src/lb_d3q19.c) | Definición modelo D3Q19 |
| [src/lb_d3q19.h](../src/lb_d3q19.h) | Vectores y pesos D3Q19 |
| [src/lb_model.c](../src/lb_model.c) | Funciones generales de modelos LB |

### 8.2 Funciones clave

#### Colisión
- `lb_collision_mrt1()` - [collision.c:221-232](../src/collision.c#L221-L232) - Kernel driver
- `lb_collision_mrt1_site()` - [collision.c:254-593](../src/collision.c#L254-L593) - Colisión por sitio

#### Propagación
- `lb_propagation()` - [propagation.c:43-50](../src/propagation.c#L43-L50) - Driver principal
- `lb_propagation_kernel()` - [propagation.c:150+](../src/propagation.c#L150) - Kernel CUDA

#### Modelo D3Q19
- `lb_d3q19_create()` - [lb_d3q19.c:30-89](../src/lb_d3q19.c#L30-L89) - Inicialización
- `lb_d3q19_matrix_ma()` - [lb_d3q19.c:107-153](../src/lb_d3q19.c#L107-L153) - Matriz de proyección

### 8.3 Variables y estructuras importantes

```c
typedef struct lb_t {
    double* f;           // Distribuciones actuales
    double* fprime;      // Distribuciones post-propagación
    lb_collide_param_t* param;  // Parámetros de colisión
    int nsite;           // Número de sitios
    int ndist;           // Número de distribuciones
} lb_t;

typedef struct hydro_t {
    double* rho;         // Densidad
    double* u;           // Velocidad [nsite][3]
    double* force;       // Fuerza [nsite][3]
    double* eta;         // Viscosidad dinámica (opcional)
} hydro_t;

typedef struct lb_collide_param_t {
    int8_t cv[NVEL][3];  // Vectores de velocidad
    double wv[NVEL];     // Pesos
    double** ma;         // Matriz M (proyección a modos)
    double** mi;         // Matriz M^(-1) (proyección inversa)
    double cs2;          // Velocidad del sonido al cuadrado
    double eta_shear;    // Viscosidad de corte
    double eta_bulk;     // Viscosidad volumétrica
    double kt;           // Temperatura (para fluctuaciones)
    int noise;           // Flag para fluctuaciones
} lb_collide_param_t;
```

---

## 9. Recuperación de ecuaciones de Navier-Stokes

### 9.1 Análisis de Chapman-Enskog

Mediante expansión de Chapman-Enskog multi-escala, se puede demostrar que el LBM recupera las ecuaciones de Navier-Stokes en el límite continuo.

**Ecuación de continuidad:**
```
∂ρ/∂t + ∇·(ρu) = 0
```

**Ecuación de momentum:**
```
∂(ρu)/∂t + ∇·(ρuu) = -∇p + ∇·(η[∇u + (∇u)^T]) + (η_bulk - 2η/3)∇(∇·u) + F
```

Donde:
- `p = ρ c_s²` (presión)
- `η = ρ c_s² (τ - 0.5) Δt` (viscosidad dinámica)
- `η_bulk` (viscosidad volumétrica)

### 9.2 Relación entre τ y viscosidad cinemática

La viscosidad cinemática `ν = η/ρ` se relaciona con el tiempo de relajación por:

```
ν = c_s² (τ - 0.5) Δt
```

Para D3Q19 con `c_s² = 1/3` y `Δt = 1`:

```
ν = (1/3)(τ - 0.5)
```

O equivalentemente:

```
τ = 3ν + 0.5
```

**Estabilidad:** El método es estable para `τ > 0.5`, lo que implica `ν > 0`.

### 9.3 Números adimensionales

**Número de Reynolds:**
```
Re = UL/ν = UL / [c_s²(τ - 0.5)]
```

Para simular flujos con alto Re, necesitamos `τ` cercano a 0.5, pero esto compromete la estabilidad. **El método MRT mejora esto** al usar diferentes `τ` para diferentes modos.

---

## 10. Ventajas del método MRT vs BGK

| Aspecto | BGK (Single Relaxation Time) | MRT (Multiple Relaxation Times) |
|---------|------------------------------|----------------------------------|
| **Estabilidad** | Limitada, especialmente en alto Re | Muy mejorada |
| **Control de viscosidades** | Solo viscosidad de corte | Shear y bulk independientes |
| **Flexibilidad** | Baja | Alta (ajuste fino de cada modo) |
| **Costo computacional** | Bajo | ~20% más alto (transformaciones M) |
| **Condiciones de frontera** | Más difíciles | Más robustas |
| **Aplicaciones multifase** | Problemático | Excelente |

---

## 11. Energía libre y acoplamiento multicomponente

### 11.1 Tensor de estrés químico

Cuando hay energía libre `fe` (por ejemplo, para fluidos multicomponente o con interfases), el tensor de estrés de equilibrio incluye una contribución extra:

```c
fe->func->str_symm_v(fe, index, symm);
seq[ia][ib] = rho * u[ia] * u[ib] + symm[ia][ib];
```

**Físicamente:**
```
Π_αβ^eq = ρ u_α u_β + Π_αβ^chemical
```

### 11.2 Ejemplos de energía libre

**Fluido binario (Cahn-Hilliard):**
```
Π_αβ^chemical = -κ (∂_α φ)(∂_β φ) + ...
```

**Electrostática:**
```
Π_αβ^chemical = ε₀ε_r [E_α E_β - (1/2)δ_αβ E²]
```

Esto permite acoplar:
- Hidrodinámica ↔ Composición
- Hidrodinámica ↔ Campo eléctrico
- Hidrodinámica ↔ Campo de fase

---

## 12. Fluctuaciones hidrodinámicas

### 12.1 Motivación física

A escala mesoscópica (micrones), las fluctuaciones térmicas son importantes. El teorema de fluctuación-disipación requiere ruido estocástico proporcional a `√(k_B T η)`.

### 12.2 Implementación

```c
if (lb->param->noise) {
    var = lb_fluctuations_var_eta(1.0/rtau, kT);
    lb_fluctuations_stress(noise, index, var, var_bulk, shat);

    // Añadir ruido al tensor de estrés
    mode[4..9] = s[ia][ib] + shat[ia][ib];

    // Añadir ruido a modos fantasma
    mode[10..18] += ghat[...];
}
```

**Varianza del ruido:**
```
var = 2 k_B T η / (Δx³ Δt)
```

Esto permite simular fenómenos como:
- Movimiento Browniano
- Difusión térmica
- Correlaciones hidrodinámicas de largo alcance

---

## 13. Conclusiones

El método Lattice Boltzmann con MRT implementado en Ludwig es un framework completo para:

1. **Hidrodinámica incompresible** a través de las ecuaciones de Navier-Stokes
2. **Flujos complejos** con geometrías irregulares
3. **Física multicomponente/multifase** mediante energías libres
4. **Acoplamiento electro-hidrodinámica** para sistemas coloidales cargados
5. **Fluctuaciones térmicas** para escala mesoscópica

### Flujo de datos resumido

```
┌─────────────────────────────────────────────────────┐
│  Inicialización: f_i = f_i^eq(ρ₀, u₀)             │
└────────────────┬────────────────────────────────────┘
                 │
        ┌────────▼────────┐
        │  LOOP TEMPORAL  │
        └────────┬────────┘
                 │
    ┌────────────▼─────────────┐
    │  1. COLISIÓN (MRT)       │
    │     f → m (proyección)   │
    │     Relax modos          │
    │     m → f (proyección⁻¹) │
    └────────────┬─────────────┘
                 │
    ┌────────────▼─────────────┐
    │  2. PROPAGACIÓN          │
    │     f(x) → f(x + c_i)    │
    └────────────┬─────────────┘
                 │
    ┌────────────▼─────────────┐
    │  3. CONDICIONES FRONTERA │
    └────────────┬─────────────┘
                 │
    ┌────────────▼─────────────┐
    │  4. OUTPUT (si necesario)│
    └────────────┬─────────────┘
                 │
                 └──────► Repetir
```

---

## Referencias bibliográficas

1. **d'Humières, D.** (1992). "Generalized Lattice-Boltzmann Equations". *Rarefied Gas Dynamics*.
2. **Lallemand, P. & Luo, L.S.** (2000). "Theory of the lattice Boltzmann method: Dispersion, dissipation, isotropy, Galilean invariance, and stability". *Physical Review E*, 61(6), 6546.
3. **Guo, Z., Zheng, C., & Shi, B.** (2002). "Discrete lattice effects on the forcing term in the lattice Boltzmann method". *Physical Review E*, 65(4), 046308.
4. **Adhikari, R., Stratford, K., Cates, M.E., & Wagner, A.J.** (2005). "Fluctuating lattice Boltzmann". *Europhysics Letters*, 71(3), 473.
5. **Stratford, K., et al.** Ludwig code documentation and papers.

---

**Documento generado a partir del código Ludwig**
**Repositorio:** `/home/bater/Sim/ludwig/`
**Autor de la documentación:** Asistente IA basado en análisis del código fuente
**Fecha:** 2025-11-02
