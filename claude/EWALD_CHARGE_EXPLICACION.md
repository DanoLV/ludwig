# Suma de Ewald para cargas puntuales: Teoría e Implementación GPU

**Archivo de referencia:** `src/ewald_charge.c`
**Función principal:** `ewald_charge_sum_full_gpu` (línea 3586)
**Referencia teórica:** Allen & Tildesley, *Computer Simulation of Liquids*, Cap. 5.5; Deserno & Holm, *J. Chem. Phys.* 109, 7694 (1998).

---

## Parte I: Fundamentos Teóricos

### 1.1 El Problema: Interacciones de Coulomb en Sistemas Periódicos

En un sistema con condiciones de contorno periódicas (PBC), el potencial electrostático en un punto **r** debido a todas las cargas del sistema y sus imágenes periódicas es:

```
φ(r) = (1/4πε) Σ_j Σ_n  q_j / |r - r_j - n·L|
```

donde:
- `q_j` es la carga de la partícula j
- `r_j` es su posición
- `n = (n_x, n_y, n_z)` indexa las imágenes periódicas
- `L` es el tamaño de la celda
- `ε` es la permitividad del medio

Esta suma **converge muy lentamente** porque el potencial de Coulomb `1/r` decae despacio. Incluso con truncamiento esférico, se necesitan muchas imágenes.

En esta implementación también hay densidad de carga continua en los nodos del lattice de Boltzmann:

```
φ(r) = (1/4πε) Σ_p q_p / |r - r_p|  +  (1/4πε) ∫ ρ_elec(r') / |r - r'| d³r'
```

donde el primer término es de coloides y el segundo de los nodos del lattice, tratados como cargas puntuales `q = ρ_elec · dV` con `dV = 1` en unidades de lattice.

---

### 1.2 El Truco de Ewald: Identidad de Splitting

La clave es la identidad matemática exacta:

```
1/r = erfc(α·r)/r  +  erf(α·r)/r
```

donde:
- `erfc(x) = 1 - erf(x) = (2/√π) ∫_x^∞ e^{-t²} dt` — función de error complementaria
- `α` es el **parámetro de splitting de Ewald** (unidades de 1/longitud)

**¿Por qué funciona este truco?**
- `erfc(αr)/r` → decae **exponencialmente** para r grande (parte de corto alcance → suma en espacio real)
- `erf(αr)/r` → es **suave** en el espacio real, lo que significa que su transformada de Fourier converge rápido (parte de largo alcance → suma en espacio Fourier)

El potencial total se separa entonces en tres contribuciones:

```
φ(r) = φ_real(r) + φ_fourier(r) - φ_self
```

---

### 1.3 Parte de Espacio Real

Para cada carga fuente j, la contribución de corto alcance al potencial en **r** es:

```
φ_real(r) = (1/4πε) Σ_j  q_j · erfc(α|r - r_j|) / |r - r_j|
```

Esta suma **converge rápidamente** para `r > rc`, donde el radio de corte típico es `rc ~ 5/α`. En unidades de lattice se trunca a `irc = ceil(rc)`.

La **fuerza** derivada de esta parte (gradiente negativo) es:

```
F_real(r) = q_i · (1/4πε) Σ_j  q_j · [erfc(αr_ij)/r_ij²  +  (2α/√π)·e^{-α²r_ij²}/r_ij] · r̂_ij
```

donde `r_ij = |r_i - r_j|` y `r̂_ij = (r_i - r_j)/r_ij`.

El primer término viene de `d/dr[erfc(αr)/r]` y el segundo del `d/dr[erfc(αr)]`:

```
d/dr [erfc(αr)/r] = -erfc(αr)/r²  -  (2α/√π)·e^{-α²r²}/r
```

En código ([ewald_charge.c:2672](../src/ewald_charge.c#L2672)):
```c
double ar = d_alpha * dist;
double r_inv = 1.0/dist;
double term1 = erfc(ar)/(dist*dist);
double term2 = (2.0*d_alpha*d_rpi)*exp(-ar*ar)*r_inv;   // d_rpi = 1/√π
double F_mag = q1*q_p*(term1+term2)/(4.0*M_PI*d_epsilon);
F_real[i] += F_mag * dr[i] * r_inv;
```

---

### 1.4 Parte de Espacio de Fourier

La contribución de largo alcance (`erf(αr)/r`) se evalúa en el espacio de Fourier. La transformada de Fourier del potencial de Coulomb atenuado es:

```
Ṽ(k) = (1/ε·V) · (4π/k²) · e^{-k²/(4α²)}
```

donde `V = L_x · L_y · L_z` es el volumen de la celda.

El potencial de la parte Fourier en el punto **r** es:

```
φ_fourier(r) = Σ_{k≠0}  G(k) · ρ̃(k) · e^{i·k·r}
```

donde el **propagador de Green en Fourier** es:

```
G(k) = (1/ε·V) · exp(-k²/(4α²)) / k²
     = b0 · exp(-k²/(4α²)) / k²
```

con `b0 = 1/(ε·V)`.

Y la **densidad de carga en Fourier** (factor de estructura) es:

```
ρ̃(k) = Σ_j  q_j · e^{i·k·r_j}
       = Σ_j  q_j · cos(k·r_j)  +  i · Σ_j  q_j · sin(k·r_j)
       = C(k) + i·S(k)
```

donde:
- `C(k) = Σ_j q_j · cos(k·r_j)` — parte real del factor de estructura
- `S(k) = Σ_j q_j · sin(k·r_j)` — parte imaginaria

Expandiendo `e^{ik·r} = cos(k·r) + i·sin(k·r)` y tomando solo la parte real del potencial:

```
φ_fourier(r) = Σ_{k≠0}  G(k) · [C(k)·cos(k·r) + S(k)·sin(k·r)]
```

En código ([ewald_charge.c:2453](../src/ewald_charge.c#L2453)):
```c
phi_fourier += factor * Gk[kn] * (Sk_cos[kn]*coskr + Sk_sin[kn]*sinkr);
```

El **campo eléctrico** es `E = -∇φ`, y usando `∇[cos(k·r)] = -k·sin(k·r)`, `∇[sin(k·r)] = k·cos(k·r)`:

```
E_fourier(r) = -∇φ_fourier = Σ_{k≠0}  G(k) · k · [S(k)·cos(k·r) - C(k)·sin(k·r)]
```

Y la **fuerza** sobre la carga `q` en **r** es `F = q · E`:

```
F_fourier(r) = q · Σ_{k≠0}  G(k) · k · [S(k)·cos(k·r) - C(k)·sin(k·r)]
```

En código ([ewald_charge.c:2594](../src/ewald_charge.c#L2594)):
```c
double im_part = Sk_sin[kn]*coskr - Sk_cos[kn]*sinkr;  // S·cos(kr) - C·sin(kr)
F_fourier[i] += factor * q1 * Gk[kn] * k[i] * im_part;
```

---

### 1.5 El Factor de Simetría (factor = 2)

Los vectores **k** se recorren solo con `kz ≥ 0` para evitar doble conteo (ya que `G(k) = G(-k)` y los términos trigonométricos son pares en cos y antisimétricos en sin). Para `kz > 0` se suma la contribución del vector `k` y su opuesto `-k` de una sola vez, dando el factor 2. Para `kz = 0` no hay duplicado:

```c
double factor = (kz_arr[kn] > 0) ? 2.0 : 1.0;
```

---

### 1.6 Los Vectores k Permitidos

Los vectores k del lattice recíproco de una caja ortorrómbica son:

```
k = (2π·n_x/L_x,  2π·n_y/L_y,  2π·n_z/L_z)
```

Solo se incluyen vectores con `0 < |k|² ≤ k_max`, donde `k_max` se determina por la convergencia deseada. En código:

```c
fkx = 2π/Lx;   fky = 2π/Ly;   fkz = 2π/Lz;
k[X] = fkx*kx;  k[Y] = fky*ky;  k[Z] = fkz*kz;
ksq = k[X]² + k[Y]² + k[Z]²;
if (ksq <= 0.0 || ksq > kmax_) continue;
Gk = b0 * exp(-r4alpha_sq * ksq) / ksq;   // r4alpha_sq = 1/(4α²)
```

---

### 1.7 Corrección de Energía Propia (Self-Energy)

La parte Fourier incluye la interacción de cada carga consigo misma (a través de las imágenes periódicas). Esto debe sustraerse:

```
φ_self = (α/√π) · q_i / (2πε)   (por partícula i)
```

En esta implementación la auto-energía no está explícitamente en `ewald_charge_sum_full_gpu` (se omite o se trata en otro lugar), pero es importante conceptualmente.

---

### 1.8 Corrección Dipolar (Deserno & Holm, 1998)

Al imponer condiciones de contorno periódicas con condiciones de contorno en vacío, hay una corrección adicional por el momento dipolar total del sistema. Según Deserno & Holm:

**Energía de corrección:**
```
E_dipole = 2π / ((1 + 2ε')·V·ε) · |M|²
```

**Campo eléctrico de corrección:**
```
E_dipole_field = -4π / ((1 + 2ε')·V·ε) · M
```

**Fuerza sobre carga q_i:**
```
F_dipole,i = q_i · E_dipole_field = -4π·q_i / ((1 + 2ε')·V·ε) · M
```

donde:
- `M = Σ_j q_j · r_j` — momento dipolar total del sistema
- `ε' = epsilon_prime` — constante dieléctrica del medio circundante
  - `ε' = 1` → condiciones de contorno en vacío
  - `ε' → ∞` → condiciones de contorno metálicas ("tinfoil") → corrección = 0
- `V` — volumen de la caja de simulación

En código ([ewald_charge.c:3853](../src/ewald_charge.c#L3853)):
```c
double dipole_prefactor = 4.0*π / ((1.0 + 2.0*epsilon_prime)*V);
E_dipole[i] = -dipole_prefactor * M[i];
```

---

## Parte II: Implementación GPU — Paso a Paso

### Flujo General

```
┌─────────────────────────────────────────────────────────────────┐
│              ewald_charge_sum_full_gpu                          │
│                                                                 │
│  [INIT]   Constantes → GPU (cudaMemcpyToSymbol)                │
│     ↓                                                           │
│  [1/6]   Pre-computar k-vectors y G(k) → GPU                  │
│     ↓                                                           │
│  [2/6]   KERNEL: S(k), C(k) en GPU (lattice + partículas)     │
│     ↓                                                           │
│  [3/6]   MPI reduce S(k), C(k)                                 │
│     ↓                                                           │
│  [DIPL]  Calcular M = Σ q·r → MPI reduce → E_dipole            │
│     ↓                                                           │
│  [4/6]   KERNEL: φ_fourier(r) en lattice (Fourier)             │
│          KERNEL: φ_real(r) partícula→nodo (real)              │
│          KERNEL: φ_real(r) nodo→nodo (real)                   │
│          [DIAG] φ_fourier + φ_real en arrays separados        │
│     ↓                                                           │
│  [5/6]   KERNEL: F_fourier + F_real en fluido                  │
│          KERNEL: E_fourier + E_real en fluido (diagnóstico)    │
│     ↓                                                           │
│  [6/6]   KERNEL: E_fourier + E_real en partículas → fex        │
│     ↓                                                           │
│  [CHECK] Verificación conservación de momento                  │
│     ↓                                                           │
│  [FREE]  Liberar memoria GPU y CPU                             │
└─────────────────────────────────────────────────────────────────┘
```

---

### Paso 0: Inicialización y Copia de Constantes a GPU

**Líneas:** [3586-3640](../src/ewald_charge.c#L3586)

Se calculan las constantes derivadas:

| Variable en código | Símbolo matemático | Significado |
|---|---|---|
| `fkx = 2π/Lx` | `Δk_x` | Paso en espacio k dirección x |
| `r4alpha_sq = 1/(4α²)` | `1/(4α²)` | Exponente del gaussiano en G(k) |
| `b0 = 1/(Lx·Ly·Lz·ε)` | `1/(V·ε)` | Prefactor del propagador Fourier |

Se copian a constantes de dispositivo con `cudaMemcpyToSymbol`:
- `d_alpha`, `d_epsilon`, `d_beta`, `d_eunit`, `d_rpi`
- `d_ewald_rc`, `d_fkx/y/z`, `d_r4alpha_sq`, `d_b0`, `d_kmax`
- `d_nk`, `d_nlocal`, `d_noffset`

---

### Paso 1: Pre-computar k-vectors y G(k)

**Líneas:** [3641-3687](../src/ewald_charge.c#L3641)

En CPU se itera sobre todos los índices de vectores **k** = `(kx, ky, kz)` con `kz ≥ 0`:

```c
for kz = 0..nk[Z]:
  for ky = -nk[Y]..nk[Y]:
    for kx = -nk[X]..nk[X]:
      ksq = (fkx*kx)² + (fky*ky)² + (fkz*kz)²
      if 0 < ksq ≤ kmax:
        kvec[kn] = (fkx*kx, fky*ky, fkz*kz)
        Gk[kn] = b0 * exp(-ksq/(4α²)) / ksq
        kz_arr[kn] = kz    // para saber si aplicar factor 2
        kn++
```

El propagador de Green es:

```
G(k) = exp(-|k|²/(4α²)) / (ε·V·|k|²)
```

Se copian `kvec_d`, `Gk_d`, `kz_arr_d` a GPU.

---

### Paso 2: Factor de Estructura S(k) y C(k) en GPU

**Líneas:** [3743-3780](../src/ewald_charge.c#L3743)

#### Kernel para el lattice: `ewald_structure_factor_lattice_kernel` ([2312](../src/ewald_charge.c#L2312))

Cada hilo GPU procesa **un nodo del lattice**. Para el nodo en posición global **r**:

```
q = ρ₀(r) - ρ₁(r)    // carga neta = densidad ión+ - densidad ión-
```

Para cada vector **k**:
```
S(k) += q · sin(k·r)   →   atomicAdd(&Sk_sin[kn], q*sin(k·r))
C(k) += q · cos(k·r)   →   atomicAdd(&Sk_cos[kn], q*cos(k·r))
```

Se usa `atomicAdd` porque múltiples hilos escriben en la misma dirección.

El índice de Ludwig (con halos) se calcula como:
```c
str_y = nlocal[Z] + 2*nhalo
str_x = str_y * (nlocal[Y] + 2*nhalo)
ludwig_idx = str_x*(nhalo + ic - 1) + str_y*(nhalo + jc - 1) + (nhalo + kc - 1)
```

La carga se lee en layout SOA (Structure of Arrays):
```c
rho0 = rho_data[nsites * 0 + ludwig_idx]   // ρ del ión positivo
rho1 = rho_data[nsites * 1 + ludwig_idx]   // ρ del ión negativo
q = rho0 - rho1
```

#### Kernel para partículas: `ewald_structure_factor_particle_kernel` ([2373](../src/ewald_charge.c#L2373))

Cada hilo procesa **una partícula coloidal**. La carga de la partícula p es `q = q0 - q1` (carga neta). Para cada **k**:

```
S(k) += q · sin(k·r_p)
C(k) += q · cos(k·r_p)
```

#### Reducción MPI

Como el dominio está distribuido entre procesos MPI, cada proceso calculó solo la contribución de sus nodos/partículas locales. Se necesita sumar global:

```c
MPI_Allreduce(MPI_IN_PLACE, Sk_sin_h, nk_actual, MPI_DOUBLE, MPI_SUM, comm);
MPI_Allreduce(MPI_IN_PLACE, Sk_cos_h, nk_actual, MPI_DOUBLE, MPI_SUM, comm);
```

---

### Paso 3: Corrección Dipolar

**Líneas:** [3782-3866](../src/ewald_charge.c#L3782)

Se calcula en CPU el momento dipolar total:

```
M_x = Σ_j  q_j · x_j
M_y = Σ_j  q_j · y_j
M_z = Σ_j  q_j · z_j
```

sumando contribuciones de coloides y del lattice. Se usa el acumulador de Kahan para alta precisión numérica.

Se hace MPI reduce del vector `(M_x, M_y, M_z, Q_total)`.

Luego se calcula el campo eléctrico de corrección dipolar:

```
prefactor = 4π / ((1 + 2ε') · V)
E_dip = -prefactor · M
```

Este campo se copia a la GPU como constante de dispositivo `d_E_dipole`.

---

### Paso 4: Potencial en los Nodos del Lattice

**Líneas:** [3990-4078](../src/ewald_charge.c#L3990)

El potencial total en cada nodo es la suma de tres contribuciones:

```
psi_data[r] = β·e · [φ_fourier(r) + φ_real,partículas(r) + φ_real,nodos(r)]
```

Los tres kernels se ejecutan secuencialmente; los kernels de espacio real **acumulan** sobre el resultado del kernel de Fourier con `atomicAdd` o `+=`.

#### Kernel Fourier: `ewald_potential_fourier_kernel` ([2418](../src/ewald_charge.c#L2418))

Cada hilo GPU procesa **un nodo del lattice**. Se calcula la parte Fourier del potencial en la posición global **r**:

```
φ_fourier(r) = Σ_{kn}  factor(kn) · G(k) · [C(k)·cos(k·r) + S(k)·sin(k·r)]
```

donde `factor = 2` si `kz > 0` (simetría), `factor = 1` si `kz = 0`.

Se agrega la corrección dipolar al potencial:

```
φ_dipole(r) = prefactor · (M · r)
```

(Nota: `E_dip = -∇φ_dip`, y `φ_dip = prefactor · (M · r)` reproduce `E_dip = -prefactor · M`.)

El potencial se almacena en formato reducido (escalado por β·e_unit):

```c
psi_data[ludwig_idx] = d_beta * d_eunit * (phi_fourier + phi_dipole);
```

#### Kernel Real (partícula→nodo): `ewald_potential_real_particle_kernel` ([2482](../src/ewald_charge.c#L2482))

Cada hilo procesa un nodo del lattice. Para cada partícula dentro del radio de corte `rc`:

```
φ_real(r) = (1/4πε) Σ_p  q_p · erfc(α·|r - r_p|) / |r - r_p|
```

```c
phi_real += q_p * erfc(d_alpha*dist) / (4.0*M_PI*d_epsilon*dist);
```

Se suma al potencial ya almacenado:
```c
psi_data[ludwig_idx] += d_beta * d_eunit * phi_real;
```

#### Kernel Real (nodo→nodo): `ewald_potential_real_lattice_kernel` ([2949](../src/ewald_charge.c#L2949))

Cada hilo procesa un nodo del lattice. Se agrega la contribución de espacio real de los **nodos vecinos del lattice** al nodo receptor, análogo a `ewald_efield_real_lattice_kernel`.

Para cada nodo vecino dentro del radio de corte `rc` (excluyendo auto-interacción):

```
φ_real(r) += (1/4πε) Σ_{r' ≠ r, |r-r'| < rc}  q_{r'} · erfc(α·|r - r'|) / |r - r'|
```

donde `q_{r'} = ρ₀(r') - ρ₁(r')` es la carga neta del nodo vecino.

```c
phi_real += q_node * erfc(d_alpha*dist) / (4.0*M_PI*d_epsilon*dist);
// ...
atomicAdd(&psi_data[ludwig_idx], d_beta * d_eunit * phi_real);
```

Se usa `atomicAdd` porque distintos hilos pueden escribir sobre el mismo `psi_data[ludwig_idx]` (en la práctica cada hilo escribe solo su propio nodo, pero el patrón `atomicAdd` es consistente con el kernel de campo eléctrico análogo).

**¿Por qué es necesario?** El factor de estructura `S(k)` incluye la contribución de todos los nodos del lattice, por lo que `φ_fourier` ya contiene la interacción nodo-nodo de largo alcance. Sin este kernel de espacio real, faltaría la parte de corto alcance de la interacción entre nodos del lattice, y el potencial total `φ = φ_fourier + φ_real` sería incorrecto.

#### Bloque diagnóstico: campos `psi_fourier` y `psi_real` ([4017-4078](../src/ewald_charge.c#L4017))

Si `psi->psi_fourier` o `psi->psi_real` están habilitados (mediante las opciones de campo en el input), se calculan las componentes por separado en arrays temporales de GPU y se almacenan como campos diagnóstico.

**Secuencia:**

1. Alocar dos arrays temporales en GPU: `psi_fourier_d` y `psi_real_d` (escalares, un valor por sitio).
2. **`psi_fourier_d`**: llenar con `ewald_potential_fourier_kernel` (incluye corrección dipolar).
3. **`psi_real_d`**:
   - Si hay partículas: `ewald_potential_real_particle_kernel` acumula la contribución partícula→nodo.
   - `ewald_potential_real_lattice_kernel` acumula la contribución nodo→nodo.
4. Copiar ambos arrays a CPU y usar `field_scalar_set` para escribir en los campos `psi->psi_fourier` y `psi->psi_real`.
5. Sincronizar de nuevo a GPU con `field_memcpy(..., tdpMemcpyHostToDevice)`.

**Invariante:**
```
psi_total[r] = psi_fourier[r] + psi_real[r]   (en unidades de β·e·φ)
```

Este bloque es **solo diagnóstico**: no afecta el potencial principal `psi_data_d` que se usa en las fuerzas.

---

### Paso 5: Fuerza en el Fluido

**Líneas:** [3890-4072](../src/ewald_charge.c#L3890)

#### Kernel Fourier de fuerza: `ewald_force_fourier_kernel` ([2544](../src/ewald_charge.c#L2544))

Cada hilo procesa un nodo con carga `q = ρ₀ - ρ₁`. La fuerza Fourier es:

```
F_fourier,i(r) = q · Σ_{kn}  factor · G(k) · k_i · [S(k)·cos(k·r) - C(k)·sin(k·r)]
```

La expresión `[S·cos(k·r) - C·sin(k·r)]` es la **parte imaginaria** de `ρ̃(k)·e^{ik·r}`. Esto surge porque:

```
∂/∂r_i [C(k)·cos(k·r) + S(k)·sin(k·r)] = k_i · [-C(k)·sin(k·r) + S(k)·cos(k·r)]
                                          = k_i · [S(k)·cos(k·r) - C(k)·sin(k·r)]
```

Se agrega la corrección dipolar a la fuerza:

```
F_dipole,i = q · E_dip,i
```

```c
F_fourier[i] += q1 * d_E_dipole[i];
```

Se acumula con `atomicAdd` en el array de fuerzas.

#### Kernel Real de fuerza: `ewald_force_real_particle_kernel` ([2622](../src/ewald_charge.c#L2622))

Para cada nodo con carga `q1` y cada partícula dentro del radio `rc`:

```
F_real,i(r) = (q1·q_p / 4πε) · [erfc(αr)/r²  +  (2α/√π)·e^{-α²r²}/r] · (r_i/r)
```

donde `r = |r_nodo - r_partícula|` y el vector apunta del nodo a la partícula.

```c
double ar = d_alpha * dist;
double r_inv = 1.0/dist;
double term1 = erfc(ar)/(dist*dist);
double term2 = (2.0*d_alpha*d_rpi)*exp(-ar*ar)*r_inv;
double F_mag = q1*q_p*(term1+term2)/(4.0*M_PI*d_epsilon);
F_real[i] += F_mag * dr[i] * r_inv;
```

#### Cálculo del Campo Eléctrico (Diagnóstico)

Si `ewald->psi->efield` existe, se calculan también los kernels de campo eléctrico:

- `ewald_efield_fourier_kernel` ([2702](../src/ewald_charge.c#L2702)) — igual que la fuerza Fourier pero **sin multiplicar por la carga del nodo receptor**:

  ```
  E_fourier,i(r) = -Σ_{kn}  factor · G(k) · k_i · [S(k)·cos(k·r) - C(k)·sin(k·r)]
  ```

  El signo negativo viene de `E = -∇φ`.

- `ewald_efield_real_particle_kernel` ([2776](../src/ewald_charge.c#L2776)) — campo de las partículas en los nodos del lattice
- `ewald_efield_real_lattice_kernel` ([2849](../src/ewald_charge.c#L2849)) — campo de nodos vecinos del lattice

---

### Paso 6: Fuerza en las Partículas

**Líneas:** [4085-4134](../src/ewald_charge.c#L4085)

#### Kernel Fourier: `ewald_particle_field_fourier_kernel` ([2938](../src/ewald_charge.c#L2938))

Cada hilo procesa **una partícula**. Se calcula el campo eléctrico Fourier en la posición de la partícula p:

```
E_fourier,i(r_p) = Σ_{kn}  factor · G(k) · k_i · [S(k)·cos(k·r_p) - C(k)·sin(k·r_p)]
```

Se almacena en `Esub` (campo eléctrico escalado: `β·e_unit·E`):

```c
E_fourier[i] += factor * d_beta * d_eunit * Gk[kn] * k[i] * im_part;
```

La fuerza sobre la partícula es `fex = q_p · E · (kT/e_unit)`:

```c
fex[i] = q_p * E_fourier[i] * kt / d_eunit;  // kt = 1/beta
```

#### Kernel Real: `ewald_particle_field_real_kernel` ([3006](../src/ewald_charge.c#L3006))

Para cada partícula p, se calcula la contribución de espacio real al campo que siente:

1. **De los nodos del lattice** — iterando sobre los nodos cercanos dentro del radio `rc`:

   ```
   E_real,i += q_nodo · β·e_unit · [erfc(αr)/r²  +  (2α/√π)·e^{-α²r²}/r] · (dr_i/r) / (4πε)
   ```

2. **De otras partículas** (p ≠ p2):

   ```
   E_real,i += q_p2 · β·e_unit · [erfc(αr)/r²  +  (2α/√π)·e^{-α²r²}/r] · (dr_i/r) / (4πε)
   ```

Se acumula sobre el campo Fourier y se actualiza la fuerza:

```c
Esub_data[3*p + i] += E_real[i];
fex_data[3*p + i] += q_p * E_real[i] * kt / d_eunit;
```

---

### Verificación de Conservación de Momento

**Líneas:** [4147-4165](../src/ewald_charge.c#L4147)

Por la tercera ley de Newton, la fuerza total sobre el sistema (fluido + partículas) debe ser cero:

```
F_total = F_fluido + F_partículas ≈ 0
```

Se calcula y reporta el residuo como diagnóstico:

```c
F_diff[i] = F_fluid[i] + F_particle[i]  // debe ser ~ 0
```

---

## Parte III: Campo Externo

**Función:** `ewald_charge_external_field` ([1152](../src/ewald_charge.c#L1152))

Se agrega un campo externo uniforme `E₀` a las partículas coloidales:

```
Esub[i] += E0[i]   // se suma directamente al campo eléctrico
```

Nota: la fuerza sobre la partícula se calcula en otro lugar como `fex = q · Esub · kT/e_unit`.

---

## Parte IV: Parámetros y Reglas de Elección

### El parámetro α

El parámetro `α` controla el balance entre las sumas en espacio real y Fourier:

- **α grande**: la parte real converge rápido (pocos términos), pero la parte Fourier necesita muchos vectores k.
- **α pequeño**: la parte Fourier converge rápido, pero la parte real necesita muchos términos.

El valor óptimo para eficiencia balanceada es:

```
α_opt ≈ 5 / (2 · rc)
```

donde `rc` es el radio de corte real. En el código:

```c
alpha_ = 5.0 / (2.0 * ewald_rc_);
```

### El radio de corte k_max

La condición de corte en el espacio de Fourier es:

```
k_max = (2·α·nk_max)²   (aproximado)
```

o más precisamente se escoge `nk_max` tal que los términos despreciados tengan error relativo < tolerancia.

### Escalado de Unidades

Los potenciales se almacenan en el campo `psi` en forma escalada:

```
psi_stored = β · e_unit · φ_físico
```

donde:
- `β = 1/(kT)` — inverso de la temperatura
- `e_unit` — unidad de carga en las unidades del código

Esto convierte el potencial electrostático en el potencial electroquímico adimensional usado en las ecuaciones de Nernst-Planck del código Ludwig.

---

## Parte V: Arquitectura GPU (CUDA)

### Organización de Hilos

Cada kernel lanza `N_nodos` hilos en bloques de 256:

```c
int threads = 256;
int blocks = (ntotal_nodes + threads - 1) / threads;
```

Cada hilo procesa **un nodo del lattice** (o una partícula en los kernels de partículas).

### Índices Ludwig con Halos

El lattice en Ludwig incluye celdas halo de ancho `nhalo` alrededor del dominio local. El índice de un nodo `(ic, jc, kc)` (1-based, sin halo) es:

```c
str_z = 1
str_y = nlocal[Z] + 2*nhalo
str_x = str_y * (nlocal[Y] + 2*nhalo)
ludwig_idx = str_x*(nhalo + ic - 1) + str_y*(nhalo + jc - 1) + str_z*(nhalo + kc - 1)
```

### Layout de Datos

Los campos vectoriales (fuerza, campo eléctrico) se almacenan en formato AOS (Array of Structures):

```
force_data[3*ludwig_idx + 0]  // Fx
force_data[3*ludwig_idx + 1]  // Fy
force_data[3*ludwig_idx + 2]  // Fz
```

Los campos escalares con múltiples especies usan SOA (Structure of Arrays):

```
rho_data[nsites * 0 + ludwig_idx]  // ρ de especie 0
rho_data[nsites * 1 + ludwig_idx]  // ρ de especie 1
```

---

## Parte VI: Resumen de Ecuaciones

### Potencial total en r

```
φ(r) = φ_real(r) + φ_fourier(r) + φ_dipole(r)

φ_real(r)    = (1/4πε) Σ_p  q_p · erfc(α|r - r_p|) / |r - r_p|          [partícula→nodo]
             + (1/4πε) Σ_{r'≠r, |r-r'|<rc}  q_{r'} · erfc(α|r-r'|) / |r-r'|  [nodo→nodo]

φ_fourier(r) = Σ_k  G(k) · [C(k)·cos(k·r) + S(k)·sin(k·r)]
               donde G(k) = exp(-k²/4α²) / (εVk²)
                     C(k) = Σ_j q_j·cos(k·r_j)
                     S(k) = Σ_j q_j·sin(k·r_j)

φ_dipole(r)  = [4π / (1+2ε')V] · (M·r)
               donde M = Σ_j q_j · r_j
```

### Fuerza total en carga q en r

```
F(r) = F_real(r) + F_fourier(r) + F_dipole(r)

F_real,i(r) = (q/4πε) Σ_j q_j · [erfc(αr_j)/r_j² + (2α/√π)e^{-α²r_j²}/r_j] · r̂_j,i

F_fourier,i(r) = q · Σ_k G(k) · k_i · [S(k)·cos(k·r) - C(k)·sin(k·r)]

F_dipole,i(r) = q · E_dip,i = -q · [4π/((1+2ε')V)] · M_i
```

### Campo eléctrico total en r (sin multiplicar por q)

```
E(r) = E_real(r) + E_fourier(r) + E_dipole

E_real,i(r) = (1/4πε) Σ_j q_j · [erfc(αr_j)/r_j² + (2α/√π)e^{-α²r_j²}/r_j] · r̂_j,i

E_fourier,i(r) = -Σ_k G(k) · k_i · [S(k)·cos(k·r) - C(k)·sin(k·r)]

E_dipole,i = -[4π/((1+2ε')V)] · M_i
```

---

*Documento actualizado el 2026-03-09. Referencia: src/ewald_charge.c, función ewald_charge_sum_full_gpu (línea 3677). Se agregó documentación de ewald_potential_real_lattice_kernel (línea 2949) y del bloque diagnóstico psi_fourier/psi_real (líneas 4017-4078).*
