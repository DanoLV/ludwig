# Cómo resuelve Ludwig el problema electrocinético con partículas subgrid (rama PETSc), paso a paso

Fecha: 2026-11-15. Pedido del usuario: el análisis paso a paso de cómo resuelve hoy Ludwig el
problema, con todas las ecuaciones, qué representan, qué suponen y qué simplifican. Se empieza por la
rama que resuelve Poisson con PETSc.

**Versión del código analizada:** rama `develop-LaValle`, commit `dc7f2059`, más el cambio sin
commitear en `src/psi_petsc.c`, que es solo un comentario. Los números de línea corresponden a ese
estado.

**Cómo se hizo:** leyendo el código. Ninguna afirmación de este documento se comprobó corriendo
simulaciones nuevas. Donde algo es una inferencia o no está verificado, se indica.

**Configuración cubierta:**

| Opción | Valor |
|---|---|
| Energía libre | `FE_ELECTRO` |
| Poisson | `electrokinetics_solver_type petsc`, con o sin paredes de potencial fijo |
| Interlacing (`interlacing_mode_`) | 0 |
| Ewald (`flag_ewald`) | 0 |
| Corrección `pm_sr` | apagada |
| Partículas | subgrid (`COLLOID_BC_SUBGRID`) cargadas |
| LB | D3Q19 (`config.mk`), relajación MRT `M10`, con o sin `lb_fluctuations` |
| Stencil de NP y de Poisson | D3Q27 (`electrokinetics_solver_stencil 27`) |
| Driver de NP | GPU (`nernst_planck_driver_d3qx_gpu`) |
| Advección de NP | `centred`, `flla` o `flla_mass` (§7.2) |

**Qué no se cubre:** FFT/Ewald, interlacing, `pm_sr`, partículas resueltas (bounce-back sobre
coloides), dos especies con sal (las ecuaciones son las mismas por especie) y el detalle de las
fuerzas de enlace y de pared. Nota: LONG-MASS-3E5 (la serie de Ingo) usó `solver_type fft`. Lo único
que cambia ahí es el paso 4 (Poisson); el resto de la cadena es la misma.

---

## 1. Unidades y símbolos

Todo está en unidades de red (LU): $\Delta x = \Delta t = 1$ y densidad de referencia del LB
$\rho_0 = 1$.

| Símbolo | Significado | En el código |
|---|---|---|
| $k_BT$ | temperatura | `physics_kt`, `temperature` |
| $\beta = 1/k_BT$ | | `psi->beta` |
| $e$ | carga unidad | `psi->e` (`electrokinetics_eunit`) |
| $\varepsilon$ | permitividad uniforme | `electrokinetics_epsilon` |
| $\ell_B = \beta e^2/(4\pi\varepsilon)$ | longitud de Bjerrum | (derivada) |
| $n_k(\mathbf x)$ | densidad numérica de la especie $k$ en el nodo $\mathbf x$ | `psi->rho` |
| $z_k$ | valencia | `psi->valency` |
| $\rho_{el}=e\sum_k z_k n_k$ | densidad de carga de los iones | `psi_rho_elec` |
| $\tilde\psi = \beta e\,\psi$ | potencial adimensional (lo que se guarda) | `psi->psi` |
| $\tilde{\mathbf E} = -\nabla\tilde\psi$ | campo adimensional | `psi_electric_field` |
| $\tilde{\mathbf E}_0$ | campo externo, en unidades de $\tilde\psi$ por LU | `electric_e0` (`psi->e0`) |
| $D_k^{in}$ | difusividad del input | `electrokinetics_d0/d1` |
| $A_0$ | factor del stencil de NP (§7.1) | |
| $D_k = A_0 D_k^{in}$ | difusividad física efectiva | |
| $f_i(\mathbf x,t)$ | distribuciones del LB, $i = 0..18$ | `lb->f` |
| $\mathbf c_i$, $w_i$ | velocidades y pesos del LB | |
| $c_s^2 = 1/3$ | velocidad del sonido al cuadrado | |
| $\eta$ | viscosidad de corte | `viscosity` |
| $\tau$ | tiempo de relajación | $\eta = \rho_0 c_s^2(\tau - 1/2)$ (`lb_relaxation_time_shear`) |
| $\mathbf c_p$, $w_p$ | vectores y pesos del stencil de NP/Poisson (D3Q27: 27 puntos) | `psi->stencil` |
| $\mathbf r_a$, $q_a$ | posición y carga (en unidades de $e$) de la partícula $a$ | `s.r`, `s.q0 - s.q1` |
| $a_h$, $a_L$ | radio hidrodinámico pedido y radio "de red" | `s.ah`, `s.al` |
| $\delta_h(\mathbf r)$ | kernel de reparto (producto de 1D) | `subgrid_kernel` |
| $\delta_{P4}$ | kernel de Peskin de 4 puntos | `d_peskin` |

---

## 2. El modelo continuo que se pretende resolver

Lo que sigue es una reconstrucción a partir del código: el código no escribe estas ecuaciones
continuas como tales.

**Iones (Nernst–Planck), sin ruido.** Para cada especie $k$:

$$
\partial_t n_k = \nabla\cdot\Big[ D_k\big(\nabla n_k + n_k\, z_k\nabla\tilde\psi\big) - \mathbf u\, n_k \Big].
$$

**Poisson (electrostática instantánea):**

$$
-\varepsilon\nabla^2\psi = \rho_{el} + e\sum_a q_a\,\delta_h(\mathbf x - \mathbf r_a)
\quad\Longleftrightarrow\quad
-\nabla^2\tilde\psi = 4\pi\ell_B\Big[\sum_k z_k n_k + \sum_a q_a\,\delta_h(\mathbf x-\mathbf r_a)\Big],
$$

con $\psi(\mathbf x + L\hat{\mathbf e}_\alpha) = \psi(\mathbf x) - E_{0,\alpha}L$ (campo externo
mediante un salto periódico). Con paredes se agrega la condición de Dirichlet
$\psi = \psi_{\rm wall}$ en los nodos de pared.

**Fluido (Navier–Stokes casi incompresible, resuelto con LB), con o sin tensor de tensiones
aleatorio:**

$$
\rho(\partial_t\mathbf u + \mathbf u\cdot\nabla\mathbf u) = -\nabla p + \eta\nabla^2\mathbf u
+ \mathbf f + \nabla\cdot\boldsymbol\sigma^{R},
$$

$$
\mathbf f = \rho_{el}\,\mathbf E + \sum_a \mathbf F_a\,\delta(\mathbf x - \mathbf r_a).
$$

**Partículas subgrid sobreamortiguadas (sin inercia):**

$$
\dot{\mathbf r}_a = \mathbf u(\mathbf r_a) + \mu_b\,\mathbf F_a + \boldsymbol\xi_a,
\qquad
\mu_b = \frac{1}{6\pi\eta}\Big(\frac{1}{a_h} - \frac{1}{a_L}\Big),
\qquad
\langle\xi_{a\alpha}\xi_{b\beta}\rangle = 2k_BT\mu_b\,\delta_{ab}\delta_{\alpha\beta}\delta(t-t'),
$$

$$
\mathbf F_a = q_a e\,\mathbf E(\mathbf r_a) + \mathbf F_a^{\rm enlaces/pares/pared}.
$$

Este es el modelo subgrid de Ludwig: Nash, Adhikari & Cates, *Phys. Rev. E* **77**, 026709 (2008).
No contrasté las ecuaciones contra ese artículo. La idea es que la movilidad total
$1/(6\pi\eta a_h)$ se reparte entre la respuesta del fluido a la fuerza repartida (que el modelo
supone igual a $1/(6\pi\eta a_L)$) y la parte "desnuda" $\mu_b$.

### 2.1 Supuestos físicos del modelo continuo

1. **Iones como campo continuo de campo medio.**
   - No hay correlaciones ion–ion ni tamaño iónico (no hay límite estérico: $\mu_{solv}=0$ en
     `fe_electro`).
   - Tampoco hay ruido propio de los iones, el "ruido desnudo" $\nabla\cdot(\sqrt{2D_kn_k}\,\mathbf W)$
     de la hidrodinámica fluctuante.
   - Con `lb_fluctuations on`, los iones solo reciben el ruido a través del término advectivo
     $\mathbf u n_k$.
2. **Electrostática instantánea (electrocuasiestática).** $\varepsilon$ es uniforme y no hay
   contraste dieléctrico. Las paredes, si existen, son conductores con potencial impuesto.
3. **Partícula puntual con la carga y la fuerza suavizadas por un kernel.** No hay volumen excluido
   entre la partícula y los iones: los iones pueden ocupar los nodos del kernel.
4. **Partícula sin inercia.** La fuerza que recibe se entrega íntegra al fluido.
5. **Fluido isotérmico, casi incompresible**, con error $O(\mathrm{Ma}^2)$ propio del LB.

---

## 3. Orden de las operaciones en un paso de tiempo

Bucle de `ludwig_run`, [src/ludwig.c:1517](../src/ludwig.c#L1517). Se denota con $n$ el paso que
empieza.

| # | Operación | Función | Línea en ludwig.c |
|---|---|---|---|
| 1 | Poner a cero la fuerza del fluido | `hydro_f_zero` | 1524 |
| 2 | Mover las partículas ($\mathbf r^{n} = \mathbf r^{n-1} + \Delta\mathbf r^{n-1}$) y calcular enlaces, pares y pared, repartidos al fluido | `ludwig_colloids_update` → `subgrid_force_from_particles` | 1528, 3968/4043 |
| 3 | Sumar la carga de las partículas a $n_k$ (temporalmente) | `subgrid_charge_from_particles` | 1661 |
| 4 | **Poisson con PETSc** | `psi_solver_petsc_solve` | 1666 |
| 5 | Restaurar $n_k$ (solo iones) | `memcpy(rho_saved)` | 1671 |
| 6 | Fuerza de los iones sobre el fluido, $\rho_{el}\mathbf E$ | `psi_force_gradmu_e` | 1855 |
| 7 | Campo en cada partícula, $\mathbf E_{sub}$ | `subgrid_update_Esub` | 1875 |
| 8 | Nernst–Planck (`multisteps` subpasos con $\tilde\psi$ fijo) | `nernst_planck_driver_d3qx_gpu` | 1836–2009 |
| 9 | Restar la media de $\tilde\psi$ | `psi_zero_mean` | 2027 |
| 10 | Fuerza eléctrica sobre cada partícula y reacción sobre el fluido | `subgrid_update_forces_electrokinetics` | 2076 |
| 11 | Colisión del LB (con fuerza y ruido) | `lb_collide` | 2216 |
| 12 | Velocidad de cada partícula, $\mathbf v = \mathbf u_I + \mu_b\mathbf F + \boldsymbol\xi$ | `subgrid_update` | 2252 |
| 13 | Bounce-back en paredes sólidas | `bounce_back_on_links`, `wall_bbl` | 2253–2255 |
| 14 | Guardar los flujos de masa del LB (para `flla` y `flla_mass`) | `nernst_planck_lattice_flux_compute` | 2272 |
| 15 | Propagación del LB | `lb_propagation` | 2278 |

**Consecuencias del orden (desfases de medio paso o un paso, propios del esquema):**
- NP (8) usa $\tilde\psi^n$ y $n_k^n$, pero la velocidad $\mathbf u$ y los flujos de masa del LB son
  los de la colisión del paso **anterior**.
- La velocidad de la partícula (12) usa la $\mathbf u$ de la colisión del **mismo** paso. Esa
  colisión ya incluye la fuerza de la propia partícula.
- La posición se actualiza al comienzo del paso **siguiente** (2).

---

## 4. Reparto de la carga de las partículas a la malla (paso 3)

[src/subgrid.c:447](../src/subgrid.c#L447):

$$
n_k(\mathbf x) \;\leftarrow\; n_k(\mathbf x) + q_{a,k}\,\delta_h(\mathbf x - \mathbf r_a),
\qquad
\delta_h(\mathbf r) = \phi(r_x)\,\phi(r_y)\,\phi(r_z),
$$

- $q_{a,0}$ = `s.q0` y $q_{a,1}$ = `s.q1`, de modo que la carga neta es $q_a = q_{a,0} - q_{a,1}$.
- $\phi$ es el kernel 1D elegido (`hann`, `peskin4`, `bspline4/6`, `kb4`, `gauss`, `peskin8`).
- La suma se acumula con suma compensada de Klein.

**Supuestos y simplificaciones:**
- La carga de la partícula se guarda en el mismo arreglo que los iones solo durante la resolución de
  Poisson. Después se restaura (paso 5), así que NP nunca transporta la carga de la partícula.
- El kernel es separable (producto de 1D). Con $\sum_{\mathbf x}\delta_h = 1$ conserva la carga
  exactamente. Sus momentos y su isotropía dependen del kernel.

---

## 5. Poisson con PETSc (paso 4)

Archivo [src/psi_petsc.c](../src/psi_petsc.c).

### 5.1 Operador discreto

Stencil D3Q27 ([src/stencil_d3q27.c:71](../src/stencil_d3q27.c#L71)), con pesos de red $w_p$
(8/27, 2/27, 1/54, 1/216):

$$
\nabla_h^2\tilde\psi(\mathbf x) = \sum_{p\neq 0} 6\,w_p\,\big[\tilde\psi(\mathbf x+\mathbf c_p) - \tilde\psi(\mathbf x)\big].
$$

- El factor 6 es la corrección de normalización del 2025-12-03 (antes −216); está publicada en el
  issue #389 y en el PR #411 de upstream.
- Reproduce $\nabla^2(x^2+y^2+z^2) = 6$ porque $\sum_p 6w_p\,c_{p\alpha}c_{p\beta} = 2\delta_{\alpha\beta}$.

**Sistema que se arma.** En `psi_solver_petsc_matrix_set` (l. 441) y `psi_solver_petsc_rhs_set`
(l. 648) se construye $A\mathbf x = \mathbf b$ con `wlaplacian[p]` $= -6w_p$. Eso es:

$$
\underbrace{-\varepsilon\,\nabla_h^2}_{A}\,\tilde\psi = \underbrace{\beta e\,\rho_{el}^{tot}}_{b}
\qquad\Longrightarrow\qquad
-\nabla_h^2\tilde\psi = \frac{\beta e^2}{\varepsilon}\sum_k z_k n_k^{tot} = 4\pi\ell_B\sum_k z_k n_k^{tot},
$$

donde $n_k^{tot}$ incluye la carga de las partículas del paso 3.

### 5.2 Caso periódico (sin paredes)

- **Espacio nulo.** La matriz es singular: la constante pertenece al núcleo. Se declara ese espacio
  nulo y además se proyecta explícitamente sobre $\mathbf b$ y $\mathbf x$ (l. 597–605; hacía falta
  con matrices en GPU). En consecuencia, la componente de carga neta de $\mathbf b$ se elimina en
  silencio. El sistema debe ser neutro: ver la electroneutralidad inicial, §10.
- **Campo externo.** El salto $\tilde\psi(\mathbf x+L) = \tilde\psi(\mathbf x) - \tilde E_0L$ entra
  en el lado derecho en las caras del dominio: $b \mathrel{+}= \varepsilon\tilde E_0 L$ en la cara
  $x=0$ y $-\varepsilon\tilde E_0L$ en la cara $x=L-1$ (l. 752–814).
  - Esa suma por cara solo es correcta porque los pesos con $c_\alpha=-1$ suman −1 (comentario
    CHANGE 20261028).

### 5.3 Con paredes de potencial fijo (Asta et al. 2019, ec. 15; CHANGE 20260926)

- **Nodo de pared (`MAP_BOUNDARY`):** se agrega la fila desacoplada $D\,x = D\,\tilde\psi_{\rm wall}$,
  con $D = \varepsilon\,w^{lap}_0$.
- **Nodo de fluido con un vecino de pared:** el enlace se **duplica**, porque la pared física está a
  medio paso de malla. La columna se elimina y su valor pasa al lado derecho, para que $A$ siga
  siendo simétrica (CG + BoomerAMG lo necesitan):

$$
b \mathrel{-}= 2\,w^{lap}_p\,\varepsilon\,\tilde\psi_{\rm wall}(\mathbf x+\mathbf c_p),
\qquad
A_{\mathbf x\mathbf x} = -\sum_p w^{lap}_p\,\varepsilon\,(1+\chi_p),
$$

  con $\chi_p = 1$ si el vecino es pared.
- **Espacio nulo:** con paredes no hay, y no se proyecta.
- **Valor de la pared** (`petsc_wall_value`, l. 258):

$$
\tilde\psi_{\rm wall}(\mathbf r) = \tilde\psi_{w,\rm lado} - \tilde{\mathbf E}_0\cdot(\mathbf r - \mathbf r_c).
$$

  La rampa existe solo si hay $\tilde{\mathbf E}_0$ tangencial. Se interpreta como un electrodo
  resistivo o un arreglo de electrodos (revisión del 2026-11-13): no es un conductor ideal único.
- **Carga inducida en cada pared:** se calcula a partir del salto de flujo en los enlaces pared–fluido
  (l. 1240) y se informa junto con la neutralidad total.

### 5.4 Resolución

- **Solver:** `KSPSolve`. El `.petscrc` actual pide CG con precondicionador BoomerAMG (hypre),
  vectores y matrices en GPU y `-ksp_rtol 5e-17`. Como `KSPSetFromOptions` se llama después de
  `KSPSetTolerances`, esos valores pisan las tolerancias del input.
- **Semilla:** el $\tilde\psi$ del paso anterior (`psi_to_da`, `KSPSetInitialGuessNonzero`).
- **Simplificación:** se resuelve **una vez por paso de LB**. El potencial no se actualiza en los
  subpasos de NP (§7.5).
- **Supuesto:** que la tolerancia alcanzada es suficiente. No la verifiqué aquí; el residuo se
  imprime cada `nfreq` pasos.

---

## 6. Campo eléctrico y fuerzas (pasos 6, 7 y 10)

### 6.1 Campo en los nodos

`psi_electric_field`, [src/psi_gradients.c](../src/psi_gradients.c). Gradiente centrado isotrópico,
con $\sum_p 3w_p\,c_{p\alpha}c_{p\beta}=\delta_{\alpha\beta}$:

$$
\tilde{\mathbf E}(\mathbf x) = -\sum_p 3w_p\,\mathbf c_p\,\tilde\psi(\mathbf x+\mathbf c_p).
$$

Los halos de $\tilde\psi$ incluyen el salto de $\tilde{\mathbf E}_0$ (`psi_halo_psijump`).

**Con paredes** se usa la forma en diferencias, duplicando los enlaces a la pared, y se lee la pared
en el punto medio del enlace (corrección $+\tilde{\mathbf E}_0\cdot\mathbf c/2$):

$$
\tilde{\mathbf E}(\mathbf x) = -\sum_p 3w_p\,\mathbf c_p\,\big[\tilde\psi(\mathbf x+\mathbf c_p)-\tilde\psi(\mathbf x)\big](1+\chi_p).
$$

En los nodos de pared, $\tilde{\mathbf E}=0$.

### 6.2 Fuerza de los iones sobre el fluido

`psi_force_gradmu_e`, [src/psi_force.c:115](../src/psi_force.c#L115):

$$
\mathbf f_{ion}(\mathbf x) = \rho_{el}(\mathbf x)\,\frac{k_BT}{e}\,\tilde{\mathbf E}(\mathbf x).
$$

- $\rho_{el}$ contiene **solo iones**: la carga de la partícula ya se restauró.
- $\tilde{\mathbf E}$ es el campo **total**: iones, partícula, campo externo y paredes.
- Si el nodo está ocupado por un coloide resuelto, la fuerza va al coloide. Con partículas subgrid,
  el nodo es fluido.
- **Simplificación:** la corrección de momento original (restar la fuerza media) está comentada.

### 6.3 Campo y fuerza sobre la partícula

`subgrid_update_Esub` ([src/subgrid.c:2241](../src/subgrid.c#L2241)) y
`subgrid_update_forces_electrokinetics` ([src/subgrid.c:1399](../src/subgrid.c#L1399)):

$$
\tilde{\mathbf E}_{sub,a} = \sum_{\mathbf x}\delta_h(\mathbf x-\mathbf r_a)\,\tilde{\mathbf E}(\mathbf x),
\qquad
\mathbf F_a^{el} = q_a\,\frac{k_BT}{e}\,\tilde{\mathbf E}_{sub,a}.
$$

La misma fuerza se reparte también sobre el fluido:

$$
\mathbf f(\mathbf x) \mathrel{+}= \mathbf F_a^{el}\,\delta_h(\mathbf x - \mathbf r_a),
$$

con el mismo kernel $\delta_h$ que la carga.

**Supuestos y simplificaciones:**
- **Entrega al fluido.** Como la partícula no tiene inercia, toda su fuerza se transmite al fluido.
  El signo es el mismo: el fluido es empujado como la partícula.
- **Conservación del momento.** Reparto y lectura usan el mismo kernel, y el gradiente discreto es
  antisimétrico. Por eso, en un sistema periódico neutro con $\mathbf E_0 = 0$, la fuerza eléctrica
  total sobre el fluido es

$$
\sum_{\mathbf x}\rho^{tot}(\mathbf x)\,\mathbf E(\mathbf x) = 0
$$

  a menos de redondeo. Es una inferencia algebraica que no comprobé numéricamente aquí. Con paredes o
  con $\mathbf E_0$, la fuerza neta no es cero: la absorben las paredes o el campo externo.
- **Autofuerza.** $\tilde{\mathbf E}_{sub}$ incluye el campo de la propia partícula. Con un kernel
  simétrico, ese campo no da autofuerza exacta cero fuera de los nodos: es la autofuerza de red
  ("pinning"). La corrección `pm_sr` que la ataca está apagada en esta configuración.
- **Halos.** Las copias de halo suman la fuerza solo en el dueño (correcciones de doble conteo del
  2026-09-26 y 2026-10-28).

---

## 7. Nernst–Planck (paso 8)

`nernst_planck_driver_d3qx_gpu`, [src/nernst_planck.c:1803](../src/nernst_planck.c#L1803). Es un
esquema de volúmenes finitos sobre los 26 enlaces del stencil D3Q27, con flujos por enlace
$J_{k,p}(\mathbf x)$ salientes de $\mathbf x$ hacia $\mathbf x+\mathbf c_p$.

### 7.1 Difusión y migración (esquema de Capuani)

Kernel `np_diffusive_fluxes_kernel` (l. 1308). Con $\mu_k(\mathbf x) = z_k\tilde\psi(\mathbf x)$,
nodo $i$ y vecino $j = \mathbf x+\mathbf c_p$:

$$
J^{dif}_{k,p}(i) = -\,D_k^{in}\;\frac{1+e^{\mu_i-\mu_j}}{2}\;\Big(n_{k,j}\,e^{\mu_j-\mu_i} - n_{k,i}\Big)\;\frac{1}{|\mathbf c_p|}.
$$

Referencia: Capuani, Pagonabarraga & Frenkel, *J. Chem. Phys.* **121**, 973 (2004). Ludwig lo
documenta como "Capuani scheme"; no contrasté la ecuación con el artículo.

**Qué representa:**
- Es la forma de Slotboom $-D\,e^{-\mu}\nabla(n\,e^{\mu})$ discretizada en el enlace.
- Se anula exactamente en el equilibrio de Boltzmann discreto, $n_j e^{\mu_j} = n_i e^{\mu_i}$. Por
  eso el estado estacionario sin flujo es exactamente PB discreta, como se verificó en
  NOFLUCT-D10/D30 (README de INGO-1PART-LBTERM-DIAG).

**Factor $A_0$:**
- La suma recorre todos los enlaces sin pesos de red. Desarrollando para campos suaves,
  $\sum_p (n_j-n_i)/|\mathbf c_p| \approx A_0\,\nabla^2 n$, con

$$
A_0 = \tfrac12\sum_p \frac{c_{p,x}^2}{|\mathbf c_p|} =
\begin{cases} 3.828 & \text{D3Q19} \\ 6.138 & \text{D3Q27}\end{cases}
$$

  (derivación propia; coincide con los valores de la memoria del proyecto).
- La difusividad física es entonces $D_k = A_0 D_k^{in}$. La movilidad de migración lleva el mismo
  factor, así que se cumple la relación de Einstein.

**Bordes:** los enlaces hacia nodos no fluidos (pared) se omiten, lo que da flujo nulo de difusión y
migración.

### 7.2 Advección

*(Sección reescrita el 2026-11-16 a pedido del usuario.)*

**Qué se quiere resolver.** El término advectivo del flujo de NP (§2), es decir, la ecuación de
Capuani et al. (2004, ec. 16):

$$
\partial_t n_k = -\nabla\cdot(n_k\,\mathbf u),
$$

donde $\mathbf u$ es la velocidad del fluido que calcula el LB. Capuani no la discretiza así: usa un
esquema de volumen de fluido (desplaza el contenido de cada celda según $\mathbf u$), que nunca se
implementó en Ludwig.

**Cómo se elige.** Con `electrokinetics_advection` en el input: `centred` (por defecto), `flla` o
`flla_mass`. Los tres escriben un flujo advectivo saliente $J^{adv}_{k,p}(i)$ por enlace. Ese flujo se
suma al difusivo (§7.1); luego se anulan los enlaces con un extremo no fluido, se antisimetriza y se
integra con Euler explícito (secuencia completa en §7.3). Solo actúa si hay hidrodinámica (`hydro` no
nulo).

**Datos de qué paso usan.** NP se ejecuta antes de la colisión (§3). Por eso los tres esquemas usan
$\mathbf u$ de la colisión del paso anterior. `flla` y `flla_mass` además usan el flujo de enlace que
se guardó al final del paso anterior (ver más abajo). En el primer paso de una corrida, o en el
primero después de un reinicio, ese flujo todavía no existe y no se aplica advección (driver,
l. 1896; variable `np_mflux_valid_`).

---

#### 7.2.1 `centred`: divergencia de red con pesos

Kernel `np_advective_fluxes_kernel` (l. 1248). Para el nodo $i$ y su vecino $j = i+\mathbf c_p$, en
cada uno de los 26 enlaces del stencil de NP:

$$
J^{adv}_{k,p}(i) = W_p\;\frac{(\mathbf u_i+\mathbf u_j)\cdot\mathbf c_p}{2}\;\frac{n_{k,i}+n_{k,j}}{2},
\qquad W_p = 2\,\texttt{wgradients}[p].
$$

Es decir: la velocidad normal al enlace por la concentración, ambas promediadas en el punto medio del
enlace, por un peso.

**De dónde salen los pesos.** `wgradients[p]` está definido en `stencil_d3q19.c:70` y
`stencil_d3q27.c:73` como $3w_p = w_p/c_s^2$. Aquí $w_p$ son los pesos de red de la familia DnQm del
mismo stencil y $c_s^2 = 1/3$. En D3Q7 vale $1/2$ (`stencil_d3q7.c:66`). Estos pesos cumplen

$$
\sum_p w_p\, c_{p\alpha} c_{p\beta} = c_s^2\,\delta_{\alpha\beta}
\quad\Longrightarrow\quad
\sum_p W_p\, c_{p\alpha} c_{p\beta} = 2\,\delta_{\alpha\beta}.
$$

En notación tensorial, $\sum_p W_p\,\mathbf c_p\mathbf c_p = 2\mathbf I$. Aquí $\mathbf c_p\mathbf c_p$
es el producto diádico (matriz $3\times3$ de componentes $c_{p\alpha}c_{p\beta}$) y $\mathbf I$ es la
**matriz identidad** $3\times 3$ ($I_{\alpha\beta}=\delta_{\alpha\beta}$). Que la suma sea
proporcional a $\mathbf I$ significa que la red no privilegia ninguna dirección a segundo orden
(isotropía).

**Por qué hacen falta.** Se aproxima el flujo del punto medio,
$J^{adv}_{k,p} \approx W_p\,(n_k\mathbf u)(\mathbf x+\mathbf c_p/2)\cdot\mathbf c_p$, y se desarrolla
por Taylor alrededor de $\mathbf x$:

$$
\sum_p J^{adv}_{k,p}
\approx \underbrace{\sum_p W_p\,c_{p\alpha}}_{=0\ (\text{enlaces opuestos})}(n_k u_\alpha)
+ \frac12\sum_p W_p\,c_{p\alpha}c_{p\beta}\;\partial_\beta(n_k u_\alpha) + O(\nabla^3)
= \nabla\cdot(n_k\mathbf u) + O(\nabla^3).
$$

El segundo término da exactamente la divergencia solo si $\tfrac12\sum_p W_p c_{p\alpha}c_{p\beta} =
\delta_{\alpha\beta}$, que es la condición de arriba. El factor 2 de $W_p$ compensa el $\tfrac12$ del
punto medio.

| Stencil NP | $W_p$ (enlaces $\lvert\mathbf c\rvert = 1,\sqrt2,\sqrt3$) | Sin pesos: $\sum_p J^{adv}_{k,p} \approx \big(\tfrac12\sum_p c_{px}^2\big)\,\nabla\cdot(n\mathbf u)$ |
|---|---|---|
| D3Q7 | 1, –, – | $1\times$ (correcto) |
| D3Q19 | $2/3,\ 1/6$, – | $5\times$ |
| D3Q27 | $4/9,\ 1/9,\ 1/36$ | $9\times$ |

**Historia.** Upstream (desde el commit `e4e9569a`, 2014, hasta master 0.23.1) suma los enlaces sin
pesos, lo que equivale a $W_p = 1$. Eso solo es correcto en D3Q7, el stencil por defecto, y con 19
o 27 puntos advecta los iones 5 o 9 veces más rápido que el fluido. Se corrigió en Ludwig local el
2026-11-12 (CHANGE 20261112, l. 1131 en CPU y l. 1843 en GPU). En D3Q7 el resultado no cambia bit a
bit.

**Referencia del operador.** Ramadugu, Thampi, Adhikari, Succi & Ansumali, *Lattice operators from
discrete hydrodynamics*, arXiv:1208.1009, ec. (10):
$\nabla\cdot\boldsymbol\phi = \frac1T\sum_i w_i\,\mathbf c_i\cdot\boldsymbol\phi(\mathbf r+\mathbf c_i)
+ O(\nabla^3)$, con $T = c_s^2 = 1/3$. Si se reemplaza el producto de promedios por el promedio del
producto (diferencia $O(c^2)$), la suma $\sum_p J^{adv}_{k,p}$ de arriba es esa ec. (10) aplicada a
$\boldsymbol\phi = n_k\mathbf u$. El paper trata solo el operador. Escribirlo como flujos de enlace
antisimétricos para la advección de NP es aplicación nuestra, no está en el paper.

**Propiedades:**
- Esquema centrado de segundo orden, sin limitador: no es monótono. En STEP-ADV (escalón advectado con
  $u = 0.01$ y $D\approx0$) dio sobre y subimpulsos de +27 %/−82 % y un retraso del frente de
  ~0.9 LU.
- Usa los 26 enlaces del stencil de NP, aunque el LB sea D3Q19.
- **Paredes y nodos sólidos.** El kernel calcula el flujo en todos los enlaces, también en los que
  van a un nodo sólido (pared o coloide BBL), donde $\mathbf u = 0$ y $n_k$ queda congelado en su valor
  inicial. Ese flujo, $W_p\,\tfrac12(\mathbf u_i\cdot\mathbf c_p)\,\tfrac12(n_{k,i}+n_{k,\text{sólido}})$,
  lo anula después `np_no_flux_boundary_kernel` (§7.3), de modo que la pared es impermeable a la
  advección, como a la difusión (§7.1).
  - Antes del 2026-11-16 ese paso no existía en el driver de GPU: el flujo salía del nodo fluido y no
    entraba en ningún nodo actualizado. Con LB térmico y paredes, el total de iones derivaba
    −4·10⁻⁵ (relativo) en 2000 pasos; con el cambio (CHANGE 20261116) se conserva a 10⁻¹⁵
    (TEST-NOFLUX-ADV).
- No conoce la densidad del LB: un campo de iones uniforme sigue uniforme aunque el fluido fluctúe en
  densidad. En NOPART con LB térmico, las fluctuaciones de iones quedaron prácticamente
  descorrelacionadas de las del fluido (corr. ≈ 0.03).

---

#### 7.2.2 `flla`: flujo de enlace del LB normalizado (DeVincentis & Thomas 2023)

Referencia: DeVincentis & Thomas, *Comput. Math. Appl.* **144**, 1–11 (2023). Los números de
ecuación de abajo son los que cita el comentario del código (l. 1505); no los volví a contrastar con el
artículo.

1. **Al final de cada paso de LB** (después de la colisión y del bounce-back, antes de la propagación;
   `ludwig.c:2268`), `np_lattice_mflux_kernel` guarda, para cada enlace $\alpha$ del LB (D3Q19):

   $$
   m_{ij} = \frac{f^*_{\bar\alpha}(\mathbf x_j)}{\rho_j} - \frac{f^*_\alpha(\mathbf x_i)}{\rho_i},
   \qquad j = i + \mathbf c_\alpha,
   $$

   donde $f^*$ son las distribuciones post-colisión y $\bar\alpha$ la dirección opuesta. La llamada
   solo lee las $f$: no modifica el LB.
2. **En el NP del paso siguiente** (`np_advective_fluxes_flla_kernel`, l. 1608), con $q = n_k$, se
   reconstruye el valor en la interfaz con un limitador (ec. 19):

   $$
   \bar q = q_{up} + \tfrac12\,\phi_{SB}(r)\,(q_{dn}-q_{up})\Big(1-\frac{|u_\alpha|}{|\mathbf c_\alpha|}\Big),
   \qquad r = \frac{q_{up}-q_{up-1}}{q_{dn}-q_{up}},
   $$

   - $u_\alpha = \tfrac12(\mathbf u_i+\mathbf u_j)\cdot\mathbf c_\alpha/|\mathbf c_\alpha|$ es la
     velocidad media del enlace.
   - El nodo *up* (aguas arriba) se elige por el **signo de $u_\alpha$**, no por el de $m_{ij}$.
     *dn* es el nodo aguas abajo y *up−1* el anterior a *up*.
   - $\phi_{SB}(r) = \max(0,\min(1,2r),\min(2,r))$ es el limitador superbee, elegido por el usuario.
   - El factor $(1-|u_\alpha|/|\mathbf c_\alpha|)$ es la corrección de Courant de tipo Lax–Wendroff.
3. El flujo saliente es $J^{adv}_{k,p}(i) = -\,m_{ij}\,\bar q$, así que
   $n_{k,i} \leftarrow n_{k,i} + \sum_\alpha m_{ij}\,\bar q$.

**Limitación medida.** Con $n_k$ uniforme, el cambio por paso es $n_k\sum_\alpha m_{ij}$. Eso no es
cero si la densidad del LB varía entre nodos: se cumple la descomposición exacta $S = D/\rho + R$
(TEST-HIDRO-dx2, `comprobar_flla_descomposicion.py`). Esas variaciones de densidad equilibran la fuerza
eléctrica y no dependen de $\eta$. Por eso, en el barrido de viscosidad, $\Delta\mu_-$ con `flla`
queda en un piso de ≈3.8·10⁻⁷ desde $\eta\approx0.5$ (12.6 veces el de `centred` en $\eta = 2$).

---

#### 7.2.3 `flla_mass`: los iones viajan con la masa del LB

Variante agregada el 2026-11-12 a pedido del usuario. Es el mismo kernel con `mode = 2`. La idea es
transportar los iones **con el mismo flujo de masa que mueve el fluido del LB**, de modo que la
concentración por unidad de masa $c_k = n_k/\rho$ se comporte como un escalar pasivo del fluido.

1. **Al final del paso de LB** (mismo punto que `flla`) se guardan, por nodo y enlace del LB,
   - el flujo de masa sin normalizar que el enlace lleva de $j$ a $i$ en la propagación:
     $J_{ij} = f^*_{\bar\alpha}(\mathbf x_j) - f^*_\alpha(\mathbf x_i)$;
   - la densidad del nodo antes de la propagación,
     $\rho^*_i = \sum_\alpha f^*_\alpha(\mathbf x_i)$. La colisión conserva la masa, así que es la
     densidad del paso que termina.

   Como la propagación hace $f_\alpha(\mathbf x,t+1) = f^*_\alpha(\mathbf x-\mathbf c_\alpha)$, se
   cumple **exactamente**
   $$\rho_i(t+1) = \rho^*_i + \sum_\alpha J_{ij}$$
   (Dubois & Lallemand 2008, Prop. 2, según el comentario del código; la identidad se verifica
   reindexando la suma).
2. **En el NP del paso siguiente:**
   - se calcula $c_k = n_k/\rho^*$ en los nodos *up*, *dn* y *up−1*;
   - se reconstruye $\bar c_k$ en la interfaz con el mismo superbee y la misma corrección de Courant de
     §7.2.2;
   - el flujo saliente es $J^{adv}_{k,p}(i) = -\,J_{ij}\,\bar c_k$, es decir
     $n_{k,i} \leftarrow n_{k,i} + \sum_\alpha J_{ij}\,\bar c_k$.
3. **Propiedad clave.** Si $c_k$ es uniforme, entonces $\bar c_k = c_k$ para cualquier limitador, y
   $n_{k,i}(t+1) = c_k(\rho^*_i + \sum_\alpha J_{ij}) = c_k\,\rho_i(t+1)$. Los iones siguen
   **exactamente** a la densidad del fluido, también con LB térmico. Comprobado en NOPART-MASS: la
   desviación estándar de los iones es 0.0173, igual a la del fluido, con correlación 1.000 y residuo
   0.0002. `flla` no tiene esta propiedad (§7.2.2).
4. **Diferencia con `flla`:** `flla` transporta $n_k$ con $m_{ij}$ (flujo dividido por la densidad de
   cada lado). `flla_mass` transporta $c_k = n_k/\rho^*$ con el flujo de masa sin normalizar. En el
   barrido de viscosidad, `flla_mass` baja con $\eta$ como `centred`, pero 1.39–1.44 veces por encima.
   Ese factor no está explicado.

**Restricciones de `flla` y `flla_mass`** (las verifica `nernst_planck_lattice_flux_compute`, l. 1712):
- Los enlaces del LB tienen que estar contenidos en el stencil de NP. Solo se usan esos enlaces: con
  stencil 27 y LB D3Q19, los 8 enlaces de esquina no llevan flujo advectivo; el arreglo de flujos se
  pone a cero en cada llamada.
- `nhalo` ≥ 2 (se lee el nodo *up−1* o *j+c*).
- Un solo proceso MPI en `flla_mass` (lee $\rho^*$ envolviendo periódicamente dentro del dominio
  local), igual que la antisimetrización (§7.3).
- Se saltan los enlaces hacia nodos no fluidos. Si *up−1* no es fluido, se usa $q_{up}$.
- Solo existen en el driver de GPU (l. 389).

---

#### 7.2.4 Resumen

| | `centred` | `flla` | `flla_mass` |
|---|---|---|---|
| Velocidad o flujo que usa | $\mathbf u$ en los nodos | $m_{ij}$ (LB, normalizado) | $J_{ij}$ (LB, masa) |
| Variable transportada | $n_k$ | $n_k$ | $c_k = n_k/\rho^*$ |
| Enlaces | 26 del stencil NP, con pesos $W_p$ | 18 del LB | 18 del LB |
| Limitador | no | superbee | superbee |
| Sigue la densidad del LB | no | parcialmente | exactamente |
| Nodos no fluidos | flujo anulado por `np_no_flux_boundary_kernel` (desde 2026-11-16) | los omite | los omite |
| Paso de los datos | $\mathbf u$ del paso anterior | $\mathbf u$ y flujo del paso anterior | ídem |

### 7.3 Secuencia de los flujos, antisimetrización y actualización

En cada subpaso, el driver de GPU arma el flujo total por enlace $J_{k,p}(i)$ (saliente de $i$ hacia
$i+\mathbf c_p$) en este orden:

1. **Advección** (§7.2): escribe $J^{adv}_{k,p}(i)$. El arreglo se pone a cero antes, en cada llamada.
2. **Difusión y migración** (§7.1, l. 1308): suma $J^{dif}_{k,p}(i)$. Omite los nodos y enlaces no
   fluidos.
3. **No-flujo** (`np_no_flux_boundary_kernel`, l. 1428; CHANGE 20261116): pone $J_{k,p}(i) = 0$ si
   $i$ o $i+\mathbf c_p$ no es fluido. Es la versión de GPU de `np_no_flux_boundary` del driver de
   CPU. Hace falta por la advección `centred`, que no mira el mapa (§7.2.1).
4. **Antisimetrización** (`np_antisym_kernel`, l. 1374): $J_p(i) \leftarrow \tfrac12[J_p(i) -
   J_{\bar p}(j)]$ y $J_{\bar p}(j) = -J_p(i)$. Lo que sale de $i$ por un enlace entra exactamente
   en $j$. La periodicidad se envuelve dentro del dominio local, así que es válida con un solo proceso
   MPI.
5. **Actualización** (`np_update_rho_kernel`, l. 1466), Euler explícito:

$$
n_{k,i}^{\,m+1} = n_{k,i}^{\,m} - \Delta t_{NP}\sum_p J_{k,p}(i),
\qquad
\Delta t_{NP} = 1/\texttt{multisteps}.
$$

  Solo se actualizan los nodos fluidos.

**Conservación.** La suma de $n_k$ sobre los nodos fluidos se conserva a redondeo si todo flujo que
sale de un nodo fluido entra en otro nodo fluido. La antisimetrización (4) lo asegura entre nodos
fluidos, y el paso de no-flujo (3) asegura que no haya flujo hacia los nodos sólidos, que no se
actualizan. Comprobado en TEST-NOFLUX-ADV (paredes + LB térmico, `centred`): $|\Delta N|/N \le
7\cdot10^{-15}$ en 2000 pasos.

### 7.4 Lo que NP no tiene

1. **Ruido propio de los iones.** No existe el término $\nabla\cdot(\sqrt{2D_kn_k}\,\mathbf W)$; el
   único ruido entra por $\mathbf u$.
2. **Potencial químico de solvatación o estérico** ($\mu_{solv}=0$ en el driver de GPU).
3. **Limitador o positividad en la difusión.** Solo `flla` y `flla_mass` tienen limitador, y solo en
   la advección.
4. **Control adaptativo de subpasos.** `maxacc` no se calcula en la GPU (se fija en 0). En
   consecuencia, `nernst_planck_adjust_multistep` nunca duplica los subpasos, y solo puede reducirlos
   si `diffacc > 0`. Esto es una inferencia de la lectura: no lo comprobé en una corrida.
*(Hasta el 2026-11-16 esta lista tenía un punto 5, "condición de no-flujo explícita en la GPU". Ya no
falta: ver §7.3, paso 3.)*

### 7.5 Subpasos

En el bucle de `multisteps` (ludwig.c:1836), **$\tilde\psi$ no se recalcula**: todos los subpasos
usan el potencial del paso 4.
- Con `multisteps 1` (el caso de LONG-MASS y H-microgel) no hay diferencia.
- Con más de un subpaso, la migración se integra con un campo congelado durante un paso de LB.

---

## 8. Media del potencial (paso 9)

`psi_zero_mean` ([src/psi.c:612](../src/psi.c#L612)): $\tilde\psi \leftarrow \tilde\psi - \langle\tilde\psi\rangle$
sobre todo el dominio, halos incluidos.

- Se llama siempre, también con paredes. Como todo lo dinámico usa diferencias de $\tilde\psi$ (NP,
  campo y carga inducida), no debería tener efecto dinámico.
- Pero el $\tilde\psi$ **escrito en los archivos** queda desplazado por una constante respecto del
  potencial de pared impuesto. Inferencia de la lectura, no comprobada.

---

## 9. Fluido: colisión MRT con fuerza y ruido (paso 11)

`lb_collision_mrt1_site`, [src/collision.c](../src/collision.c), D3Q19.

1. **Modos.** $m_a = \sum_i M_{ai}f_i$ (`d3q19_f2mode_chunk`), y vuelta con
   $f_i = \sum_a M^{-1}_{ia}m_a$ (`d3q19_mode2f_chunk_kahan`, con suma compensada; es el arreglo del
   issue #355).
2. **Velocidad con media fuerza:**

$$
\mathbf u = \frac{1}{\rho}\Big(\sum_i f_i\mathbf c_i + \tfrac12\mathbf F\Big),
\qquad
\mathbf F = \mathbf f_{ion} + \sum_a\big(\mathbf F_a^{el}+\mathbf F_a^{enl}\big)\delta + \mathbf F_{global}.
$$

3. **Relajación del tensor.**
   - La parte sin traza se relaja con $\omega = 1/\tau$ y la traza con $\omega_b$, hacia
     $s^{eq}_{\alpha\beta} = \rho u_\alpha u_\beta$.
   - El código agrega además una corrección de fuerza:

$$
s'_{\alpha\beta} = s_{\alpha\beta} - \omega\,(s_{\alpha\beta}-s^{eq}_{\alpha\beta}) + (2-\omega)\,(u_\alpha F_\beta + F_\alpha u_\beta) + \hat s_{\alpha\beta}.
$$

   - El comentario del código dice que esa corrección "assumes equal relaxation times". Es el
     forzado estilo Guo de Ludwig, pero **no verifiqué** la normalización de los modos que hace
     equivalente el factor $(2-\omega)$ a $(1-\omega/2)$ de Guo, Zheng & Shi, *PRE* **65**, 046308
     (2002).
4. **Momento:** $\mathbf j' = \mathbf j + \mathbf F$.
5. **Modos fantasma.** Con `M10`, $\omega_g = 1$ (se reinician a cero en cada paso), más ruido.
6. **Ruido** (`lb_fluctuations on`). Referencias: Adhikari, Stratford, Cates & Wagner, *Europhys.
   Lett.* **71**, 473 (2005); Dünweg, Schiller & Ladd, *PRE* **76**, 036704 (2007); no contrasté
   ecuaciones.
   - Desviación del tensor aleatorio sin traza, por componente fuera de la diagonal:

$$
\sigma_\eta = \sqrt{k_BT/c_s^2}\;\sqrt{1/9}\;\sqrt{\frac{2\tau-1}{\tau^2}}
$$

     ($\sqrt2$ veces eso en la diagonal).
   - Para la traza el factor es $\sqrt{2/9}$, y para los fantasmas, $\sqrt{k_BT\,r_{na}}$ por la
     misma raíz.
   - **Simplificación:** en estas expresiones no aparece la densidad local; equivale a usar
     $\rho=\rho_0=1$.
   - Los números aleatorios no son gaussianos: provienen de una tabla discreta de 8 valores con
     varianza unidad (Ladd, *Comput. Phys. Commun.* **180**, 2140 (2009)), con un generador por nodo
     (KISS99).

---

## 10. Velocidad y movimiento de la partícula (pasos 12 y 2)

`subgrid_update` ([src/subgrid.c:3556](../src/subgrid.c#L3556)):

$$
\mathbf v_a = \mathbf u_I(\mathbf r_a) + \mu_b\,\mathbf F_a^{ex} + \boldsymbol\xi_a,
\qquad
\mathbf u_I(\mathbf r_a) = \sum_{\mathbf x}\delta_{P4}(\mathbf x-\mathbf r_a)\,\mathbf u(\mathbf x),
$$

- $\mathbf F^{ex}$ = enlaces + pares + pared + fuerza eléctrica.
- $\xi_{a\alpha} = \sqrt{2k_BT\mu_b}\;g$, con $g$ gaussiana **truncada en $|g|<3$**, solo si
  `lb_fluctuations on`.

La posición se mueve con Euler explícito al comienzo del paso siguiente (`colloids.c:1078`):

$$
\mathbf r_a^{\,n+1} = \mathbf r_a^{\,n} + \mathbf v_a\,\Delta t.
$$

**Supuestos y simplificaciones:**
1. **Partícula sobreamortiguada**, sin masa ni memoria. Su movilidad es la del fluido local más
   $\mu_b$.
2. **Interpolación con Peskin 4 siempre.** La velocidad se interpola con $\delta_{P4}$
   (`drange_ = 1`), **cualquiera sea** `subgrid_kernel`. Las fuerzas de enlace también se reparten
   con $\delta_{P4}$. En cambio, la carga, el campo y la fuerza eléctrica usan el kernel elegido
   (Hann-8 en las corridas actuales).
   - La calibración de $a_L$ (la respuesta del fluido a la fuerza repartida, interpolada de vuelta)
     depende del par de kernels.
   - Con kernels distintos para repartir la fuerza eléctrica y para interpolar, el $a_L$ efectivo
     puede no ser el del input. **No lo verifiqué**: queda señalado como posible inconsistencia.
3. **Ruido truncado.** El truncamiento en 3σ reduce la varianza a 0.9733 de la nominal (varianza de
   la normal truncada, calculada con `scipy.stats.truncnorm`). Equivale a usar $D_b$ un 2.7 % menor.
4. **La $\mathbf u$ interpolada** es la de la colisión del mismo paso. Incluye la respuesta (con
   media fuerza) a la fuerza de la propia partícula.
5. **Control de salto.** La comprobación de `drmax` en `colloids.c:1074` es
   `s.dr[ia] > drmax`, **sin valor absoluto**: un salto grande negativo no detiene la corrida.
   Observación de la lectura; no lo corregí.

---

## 11. Condiciones iniciales relevantes

**Electroneutralidad** (`psi_electroneutral`, [src/psi.c:937](../src/psi.c#L937); se omite con
`electrokinetics_electroneutral no`, que solo está permitido con paredes):
- Se suma la carga total $Q = \sum z n + \sum_a q_a$.
- Se agrega la contracarga uniforme $|Q|/V_{f,ef}$ en todos los nodos fluidos de la especie de signo
  contrario. $V_{f,ef}$ excluye los nodos a menos de 0.1 LU de una partícula.

**Supuesto:** el sistema periódico debe ser exactamente neutro, porque el espacio nulo de Poisson
(§5.2) descarta la carga neta.

---

## 12. Resumen de supuestos y simplificaciones

| # | Supuesto o simplificación | Dónde | Efecto esperado |
|---|---|---|---|
| 1 | Iones de campo medio, sin ruido propio | NP, §7.4 | Con LB térmico, el estado estacionario no tiene por qué ser Boltzmann; la nube cerca de una partícula móvil sale por debajo de PB (INGO-1PART-LBTERM-DIAG) |
| 2 | Sin tamaño iónico ni límite estérico | $\mu_{solv}=0$ | Concentraciones sin cota a potenciales altos |
| 3 | $\varepsilon$ uniforme, electrostática instantánea | Poisson | Sin imágenes dieléctricas; las paredes solo pueden ser conductoras |
| 4 | Partícula puntual suavizada, sin volumen excluido con los iones | Kernel | Los iones pueden ocupar el centro de la partícula; el potencial cerca de $r\lesssim$ soporte depende del kernel |
| 5 | Partícula sin inercia; movilidad $= u_I + \mu_b F$ | `subgrid_update` | Válido para tiempos $\gg$ relajación de la velocidad |
| 6 | Kernels distintos: Peskin-4 para $\mathbf u$ y enlaces, kernel elegido para carga y fuerza eléctrica | §10 | Posible $a_L$ efectivo distinto del calibrado (no verificado) |
| 7 | Ruido de la partícula truncado en 3σ | §10 | $D_b$ efectivo × 0.973 |
| 8 | Ruido del LB con $\rho=1$ y deviates discretos | §9 | Correcto si $\rho\approx1$ |
| 9 | Un Poisson por paso de LB; $\tilde\psi$ fijo en los subpasos de NP | §5.4, §7.5 | Irrelevante con `multisteps 1` |
| 10 | NP explícito (Euler) | §7.3 | Límite de estabilidad en $D\Delta t$ (no analizado aquí) |
| 11 | NP usa la $\mathbf u$ y los flujos del paso anterior | §3 | Desfase de un paso entre fluido e iones |
| 12 | Advección centrada sin limitador (`centred`) | §7.2 | Oscilaciones a Péclet de celda alto |
| 13 | `flla_mass` y la antisimetrización: solo 1 proceso MPI | §7.2–7.3 | Restricción de paralelismo |
| 14 | Subpasos de NP no adaptativos en GPU | §7.4 | La precisión no se controla automáticamente |
| 15 | Autofuerza de red no corregida (`pm_sr` apagado) | §6.3 | Fuerza espuria dependiente de la posición dentro de la celda |
| 16 | Corrección de momento de la fuerza iónica comentada | §6.2 | Inofensivo si el sistema es neutro y periódico (inferencia) |
| 17 | `psi_zero_mean` también con paredes | §8 | Archivos de $\psi$ desplazados por una constante (inferencia) |
| 18 | `drmax` sin valor absoluto | §10 | Saltos negativos grandes no detectados |

## 13. Qué no está verificado en este documento

- La equivalencia de la corrección de fuerza $(2-\omega)$ con el forzado de Guo (§9).
- Las ecuaciones exactas de las referencias citadas (Nash 2008, Capuani 2004, Adhikari 2005,
  Dünweg 2007, Guo 2002): se citan como fuente del método según los comentarios del código y la
  documentación de Ludwig, no contrastadas.
- La inconsistencia de kernels en $a_L$ (punto 6), el desplazamiento de $\psi$ con paredes
  (punto 17) y que los subpasos no se adapten (punto 14): son inferencias de la lectura.
- Las fuerzas de enlace, pares y pared (Steele) entran en $\mathbf F^{ex}$, pero no se describen
  aquí.
