# Análisis: separar sigma para fluido y partículas en `ewald_charge_sum_full_gaussian_gpu`

Fecha: 2026-05-14
Archivo principal: [src/ewald_charge.c](src/ewald_charge.c)

## 1. Objetivo

Actualmente `ewald_charge_sum_full_gaussian_gpu` usa un único parámetro `sigma` para **todas** las cargas (partículas coloidales + nodos del lattice). Queremos permitir dos anchuras gaussianas independientes:

- `sigma_p`: anchura de las cargas puntuales en las partículas (coloides).
- `sigma_f`: anchura de la "carga" continua del fluido en cada nodo del lattice.

Esto es necesario porque físicamente:
- la carga del coloide es discreta y vive sobre el centro del coloide → conviene un sigma propio (típicamente menor / ajustado al radio del coloide o al kernel subgrid).
- la "carga" del nodo es una representación regularizada de una densidad continua de iones → su anchura natural es el espaciado del lattice, distinto del coloide.

Con dos sigmas distintas la fórmula `erf(r/σ)/r` (interacción entre dos gaussianas idénticas) deja de ser válida en las interacciones cruzadas y hay que usar la convolución de dos gaussianas con anchuras distintas.

## 2. Marco matemático

Una gaussiana normalizada de anchura σ tiene la forma
$$g_\sigma(r) = \frac{1}{(\pi)^{3/2}\sigma^3} \exp(-r^2/\sigma^2).$$

El potencial de Coulomb generado por dicha gaussiana es
$$\phi_\sigma(r) = \frac{\mathrm{erf}(r/\sigma)}{4\pi\varepsilon\,r}.$$

La interacción **entre dos gaussianas** de anchuras σ_a y σ_b separadas distancia r se obtiene convolucionando ambas distribuciones y resulta en una gaussiana de anchura
$$\sigma_{ab} = \sqrt{\sigma_a^2 + \sigma_b^2}.$$
Por lo tanto el potencial cruzado es
$$\phi_{ab}(r) = \frac{\mathrm{erf}(r/\sigma_{ab})}{4\pi\varepsilon\,r}.$$

### 2.1 Caso de Ewald gaussiano con dos σ

El splitting de Ewald se construye introduciendo además la anchura algorítmica η (en el código `d_gauss_eta = alpha_`). El truco habitual es escribir

$$\frac{\mathrm{erf}(r/\sigma_{ab})}{r} = \underbrace{\frac{\mathrm{erf}(r/\sigma_{ab}) - \mathrm{erf}(\eta r)}{r}}_{\text{real-space (corto alcance)}}
+ \underbrace{\frac{\mathrm{erf}(\eta r)}{r}}_{\text{recíproco (largo alcance)}}.$$

- La parte real-space se calcula con un par `(σ_{ab}, η)` que depende de **qué dos especies interaccionan**.
- La parte recíproca (Fourier) usa solo η, y el factor de forma multiplicativo de cada especie es
  $$\tilde g_\sigma(k) = \exp(-k^2 \sigma^2 / 4).$$
  Es decir, al colocar la carga en la malla espectral se multiplica por `gaussian_fk = exp(-k²σ²/4)` con el **σ de esa especie** (partícula o nodo). El producto en Fourier de las dos gaussianas reproduce automáticamente la interacción cruzada con σ_{ab} = √(σ_p² + σ_f²), sin ningún esfuerzo adicional.

### 2.2 Self-energy (corrección de auto-interacción)

Cada especie tiene su propio término de self:
$$U_\text{self}(\sigma) = \frac{1}{4\pi\varepsilon}\cdot\frac{1}{\sigma\sqrt{2\pi}}\sum_i q_i^2.$$
Si se usan dos σ, hay que calcular **dos** self-energies, una con σ_p para las partículas y otra con σ_f para los nodos. Hoy no se usa explícitamente en el código que estamos modificando (no veo `self_table` aplicada aquí), pero es importante recordarlo si se conecta con el resto del solver.

### 2.3 Taylor en r→0

En el código aparecen dos límites Taylor para el caso `dist < 1e-5/1e-6`:
- Potencial: `phi_diff → (2/√π)·(1/σ - η)`
- Fuerza/derivada: `-dphi/dr → (4/(3√π))·(1/σ³ - η³)·r`

Cuando σ pasa a ser σ_{ab} = √(σ_a² + σ_b²), todos esos `σ` deben sustituirse por σ_{ab} **específico al par de especies** que interaccionan.

## 3. Mapeo de cambios por kernel

A continuación, para cada kernel involucrado, indico qué σ se necesita y qué hay que cambiar.

| Kernel | Interacción | σ que se necesita | Comentario |
|---|---|---|---|
| `ewald_structure_factor_lattice_kernel_gaussian` ([src/ewald_charge.c:8849](src/ewald_charge.c#L8849)) | nodo→Fourier | σ_f (vía `gaussian_fk`) | Hay que precomputar un `gaussian_fk_fluid[kn] = exp(-k²σ_f²/4)` y pasarlo aquí. |
| `ewald_structure_factor_particle_kernel_gaussian` ([src/ewald_charge.c:8905](src/ewald_charge.c#L8905)) | partícula→Fourier | σ_p (vía `gaussian_fk`) | Precomputar `gaussian_fk_particle[kn] = exp(-k²σ_p²/4)` y pasarlo. |
| `ewald_potential_real_particle_kernel_gaussian` ([src/ewald_charge.c:8954](src/ewald_charge.c#L8954)) | potencial real en **nodo** debido a una **partícula** | σ_pf = √(σ_p² + σ_f²) | Cambiar `d_gauss_sigma` por `d_gauss_sigma_pf`. |
| `ewald_potential_real_lattice_kernel_gaussian` ([src/ewald_charge.c:9015](src/ewald_charge.c#L9015)) | potencial real en **nodo** debido a otro **nodo** | σ_ff = σ_f·√2 | Usar `d_gauss_sigma_ff`. |
| `ewald_force_real_particle_kernel_gaussian` ([src/ewald_charge.c:9087](src/ewald_charge.c#L9087)) | fuerza real sobre **nodo** (q1=lattice) debida a **partícula** | σ_pf = √(σ_p² + σ_f²) | Usar `d_gauss_sigma_pf`. |
| `ewald_efield_real_particle_kernel_gaussian` ([src/ewald_charge.c:9171](src/ewald_charge.c#L9171)) | E real en **nodo** debido a **partícula** | σ_pf | Usar `d_gauss_sigma_pf`. |
| `ewald_efield_real_lattice_kernel_gaussian` ([src/ewald_charge.c:9252](src/ewald_charge.c#L9252)) | E real en **nodo** debido a otro **nodo** | σ_ff = σ_f·√2 | Usar `d_gauss_sigma_ff`. |
| `ewald_particle_field_real_kernel_gaussian` ([src/ewald_charge.c:9351](src/ewald_charge.c#L9351)) | E y fuerza real en **partícula** debido a (a) nodos y (b) otras **partículas** | σ_pf para (a), σ_pp = σ_p·√2 para (b) | El kernel tiene **dos loops** y cada uno necesita un σ distinto. |

### 3.1 Punto crítico: `ewald_particle_field_real_kernel_gaussian`

Este kernel mezcla dos tipos de fuentes (nodos y partículas) dentro del mismo kernel. Hay que pasarle **dos** sigmas (o equivalentemente σ_pf y σ_pp) y aplicar el que corresponde en cada bucle interno. Los límites Taylor también deben usar el σ correcto en cada caso. Es el único kernel donde el cambio no es una simple sustitución global de `d_gauss_sigma`.

## 4. Plan de implementación

### 4.1 Nuevas variables device-constant

En la sección de constantes (cerca de [src/ewald_charge.c:8951-8952](src/ewald_charge.c#L8951-L8952)):

```c
__device__ __constant__ double d_gauss_sigma_p;     /* sigma de partículas */
__device__ __constant__ double d_gauss_sigma_f;     /* sigma del fluido (nodos) */
__device__ __constant__ double d_gauss_sigma_pp;    /* σ_p · √2  para par partícula-partícula */
__device__ __constant__ double d_gauss_sigma_ff;    /* σ_f · √2  para par nodo-nodo */
__device__ __constant__ double d_gauss_sigma_pf;    /* √(σ_p² + σ_f²) para par partícula-nodo */
```

Mantener `d_gauss_eta` y `d_gauss_sigma` originales por compatibilidad si fuese necesario, o eliminar `d_gauss_sigma` del todo. Recomendación: **eliminarlo** para evitar bugs silenciosos.

### 4.2 Cambio de firma de la función pública

[src/ewald_charge.h:152](src/ewald_charge.h#L152):

```c
int ewald_charge_sum_full_gaussian_gpu(ewald_charge_t * ewald, FILE * fp,
                                       double sigma_particle,
                                       double sigma_fluid);
```

Y propagar el cambio al solver [src/psi_solver_ewald.c:141-147](src/psi_solver_ewald.c#L141-L147):
- añadir `solver->sigma_particle` y `solver->sigma_fluid` en `psi_solver_ewald_gaussian_t`
- ajustar `psi_solver_ewald_gaussian_create` para aceptar/almacenar ambos
- pasar ambos al kernel desde `psi_solver_ewald_gaussian_solve`

Quien construye el solver (mirar [src/ludwig.c:960](src/ludwig.c#L960) y entorno) debe leer dos valores del input file, p.ej. `ewald_gauss_sigma_particle` y `ewald_gauss_sigma_fluid`.

### 4.3 Cambios dentro de `ewald_charge_sum_full_gaussian_gpu`

En [src/ewald_charge.c:9515-9580](src/ewald_charge.c#L9515-L9580):

1. Calcular las tres sigmas efectivas:
   ```c
   double sigma_pp = sigma_p * sqrt(2.0);
   double sigma_ff = sigma_f * sqrt(2.0);
   double sigma_pf = sqrt(sigma_p*sigma_p + sigma_f*sigma_f);
   ```
2. Subir las cinco constantes a device con `cudaMemcpyToSymbol` (eliminar el `cudaMemcpyToSymbol(d_gauss_sigma, ...)` actual en línea 9578).

3. Precomputar **dos** arrays de form factor en Fourier:
   ```c
   double* gaussian_fk_p_h = malloc(nktot_ * sizeof(double));   // exp(-k²σ_p²/4)
   double* gaussian_fk_f_h = malloc(nktot_ * sizeof(double));   // exp(-k²σ_f²/4)
   ```
   y subir ambos a device (`gaussian_fk_p_d`, `gaussian_fk_f_d`). El array actual `gaussian_fk_d` queda obsoleto.

4. Pasar `gaussian_fk_p_d` al kernel de structure-factor de partículas y `gaussian_fk_f_d` al de lattice.

5. Liberar ambos arrays al final (sustituir el `free(gaussian_fk_h)` y `cudaFree(gaussian_fk_d)`).

### 4.4 Cambios en los kernels

Reemplazar dentro de cada kernel los usos de `d_gauss_sigma` según la tabla de §3:

- Kernels marcados con σ_pf → usar `d_gauss_sigma_pf` (inv_s3 e `u_sig` se reescriben con esa variable).
- Kernels marcados con σ_ff → usar `d_gauss_sigma_ff`.
- En `ewald_particle_field_real_kernel_gaussian`:
  - Loop "nodos" ([9382-9434](src/ewald_charge.c#L9382-L9434)) → `d_gauss_sigma_pf`.
  - Loop "otras partículas" ([9436-9476](src/ewald_charge.c#L9436-L9476)) → `d_gauss_sigma_pp`.

Conviene refactorizar la lógica `if (dist < taylor_thr) {...} else {...}` a una pequeña función `__device__` que tome `sigma_eff` y `eta` como argumentos y devuelva `dphi_diff` o `dphi_dr_diff`. Eso reduce la duplicación y hace los cambios menos propensos a errores.

Esbozo:
```c
__device__ inline double gauss_phi_diff(double r, double sig, double eta) {
  if (r < 1.0e-6) return (2.0/sqrt(M_PI))*(1.0/sig - eta);
  return (erf(r/sig) - erf(eta*r)) / r;
}
__device__ inline double gauss_minus_dphi_dr(double r, double sig, double eta) {
  if (r < 1.0e-5) {
    double inv_s3 = 1.0/(sig*sig*sig);
    double eta3 = eta*eta*eta;
    return (4.0/(3.0*sqrt(M_PI)))*(inv_s3 - eta3)*r;
  }
  double r2 = r*r, r_inv = 1.0/r;
  double u_s = r/sig, u_e = eta*r;
  double dg_s = (2.0/(sqrt(M_PI)*sig))*exp(-u_s*u_s)*r_inv - erf(u_s)/r2;
  double dg_e = (2.0*eta/sqrt(M_PI))*exp(-u_e*u_e)*r_inv - erf(u_e)/r2;
  return -(dg_s - dg_e);
}
```
Cada kernel pasaría a llamarlas con el `sig` adecuado.

## 5. Cuestiones abiertas

1. **Self-energy.** ¿Se está aplicando self correction en este código path? Si sí, hay que dividirla en aporte de partículas (σ_p) y aporte de nodos (σ_f). Revisar `ewald_self_table_t` y dónde se usa. En la función actual no aparece llamada al self, parece que se ignora — confirmar con el resto del solver.
2. **Cutoff real `ewald_rc_`.** El cutoff `ewald_rc_` actualmente es una constante global. Con dos sigmas, el "alcance efectivo" de la parte real-space depende de σ_{ab}: la parte real-space decae como `erfc(η r) − erfc(r/σ)`, así que para σ pequeñas la parte real-space domina hasta más lejos. Verificar que `ewald_rc_` es suficiente para σ_pp = max(σ_p√2, σ_f√2, σ_pf). Si no, dar opción a recalcular `irc` por canal.
3. **Compatibilidad con cuFFT path.** Si en algún momento se quiere usar `psi_fft_pn` con cargas gaussianas separadas, el deconvolution kernel también necesitará dos sigmas. No es necesario tocarlo ahora pero conviene dejarlo apuntado.
4. **Naming en input file.** Acordar las keys: propuesta `ewald_gauss_sigma_particle` / `ewald_gauss_sigma_fluid`, o como tag único `ewald_gauss_sigma = sigma_p sigma_f`. Mantener compatibilidad: si solo se da uno, asignar el mismo a los dos (equivalente al comportamiento actual).
5. **Tests de regresión.** Sanity check: si `sigma_p == sigma_f`, el resultado numérico debe coincidir con la versión single-sigma actual hasta precisión float-double. Añadir un test que compare ambos modos con la misma sigma.

## 6. Resumen de archivos a modificar

- [src/ewald_charge.c](src/ewald_charge.c) — kernels (§3) y función principal (§4.3).
- [src/ewald_charge.h](src/ewald_charge.h) — firma pública (§4.2).
- [src/psi_solver_ewald.c](src/psi_solver_ewald.c) — almacenar y pasar ambas sigmas; ajustar `psi_solver_ewald_gaussian_create` y `..._solve`.
- [src/psi_solver_ewald.h](src/psi_solver_ewald.h) — firma de create.
- [src/ludwig.c:960](src/ludwig.c#L960) y aledaños — lectura del input y construcción del solver.
- Input file de ejemplo en `microgel/local/input` — añadir las dos keys.

Todas las modificaciones deben seguir la convención del proyecto: envolverlas en bloques `/*CHANGE INIT - <fecha/desc> */ ... /*CHANGE END - <desc> */` y no borrar el código previo (dejarlo comentado dentro del bloque).
