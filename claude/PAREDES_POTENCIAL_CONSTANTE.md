# Paredes a potencial constante en el solver PETSc (sistema no periódico)

**Estado:** implementado y validado el 2026-09-26 (camino de permitividad uniforme).
**Referencia:** A. J. Asta, I. Palaia, E. Trizac, M. Levesque, B. Rotenberg,
*Lattice Boltzmann Electrokinetics simulation of nanocapacitors*,
arXiv:1907.04732 (2019), sección II-C. PDF en la bibliografía:
`/mnt/datos/Github/Trabajo/Papers/Asta et al. - 2019 - Lattice Boltzmann Electrokinetics simulation of nanocapacitors.pdf`

**Archivos modificados:** `src/psi_petsc.c`, `src/psi_petsc.h`, `src/psi.h`,
`src/psi_gradients.c`, `src/ludwig.c`, `src/subgrid.c`. Todas las modificaciones están entre
marcas `/*CHANGE INIT - 20260926 ... */` y `/*CHANGE END - 20260926 */`.

**Pruebas:** `capacitor/local/WALL-DIRICHLET/` (entradas, salidas, scripts,
gráficos en `graficos/` e índice en `README.md`; `correr_tests.sh` las
reproduce todas).

---

## 0. Resumen

Hasta ahora el solver de Poisson de Ludwig con PETSc era estrictamente
periódico. Con estos cambios, los nodos de pared (`MAP_BOUNDARY`) pueden
funcionar como **electrodos metálicos a potencial fijo**, y el sistema deja de
ser periódico en la dirección normal a las paredes.

Qué quedó funcionando:

- Poisson con condición de Dirichlet en las paredes, con la pared física en el
  **plano medio** entre el último nodo fluido y el primer nodo sólido
  (segundo orden en el líquido).
- Cálculo de la **carga inducida** en cada electrodo en cada resolución.
- **Campo eléctrico** consistente con esa condición de borde en los nodos
  fluidos junto a la pared, que es el que usan las fuerzas sobre el fluido y
  sobre las partículas.
- Sin las claves nuevas en el input, el comportamiento es **bit a bit idéntico**
  al original.

Qué no está hecho o no está verificado: ver la sección 6.

---

## 1. El problema y la solución del paper

### 1.1 Dónde está la pared

En LB con *bounce-back*, la interfaz sólido–líquido hidrodinámica no está en
el último nodo sólido sino a mitad de camino entre el último nodo fluido y el
primer nodo sólido. Para que la electrostática sea consistente con la
hidrodinámica, el potencial del electrodo ψ_s tiene que imponerse en ese mismo
plano medio.

La forma ingenua (resolver Poisson en el líquido con el laplaciano de siempre
y fijar ψ_s en los nodos sólidos) pone la pared en el nodo sólido, a una
distancia h del último nodo fluido en lugar de h/2. Eso es de **primer orden**
en el espaciado de red (Fig. 1 del paper).

### 1.2 El laplaciano modificado

El laplaciano consistente con la red (ec. 12 del paper) es

```
∇²φ(r) = (2 / c_s² Δt²) Σ_i w_i [φ(r + c_i Δt) − φ(r)]
```

Con la función característica del sólido (ec. 14), χ_s(r) = 1 en un nodo
sólido y 0 en uno fluido, el paper lo reemplaza por (ec. 15)

```
∇²φ(r) = (2 / c_s² Δt²) Σ_i w_i [φ(r + c_i Δt) − φ(r)] · [1 + χ_s(r + c_i Δt) − χ_s(r)]
```

Para un nodo **fluido** el factor vale 1 en un enlace fluido–fluido y **2** en
un enlace hacia un nodo sólido. En palabras del paper: *"para cada enlace de
frontera se multiplica por 2 la diferencia"*. La razón es que el valor
conocido ψ_s está a h/2 y no a h, así que la pendiente sobre ese enlace es el
doble. Los nodos sólidos no se resuelven: se mantienen en ψ_s.

### 1.3 La carga inducida en el electrodo (ec. 16) — corrección al texto del paper

El paper calcula la carga del electrodo con Poisson al revés,

```
Q = Δx³ Σ_{r ∈ electrodo} ρ_el(r) = −(e Δx³ / 4π l_B) Σ_{r ∈ electrodo} ∇²φ(r)
```

y dice que ese laplaciano se calcula "con la ec. 15". Leída literalmente, con
r = nodo del electrodo (χ_s(r) = 1), la ec. 15 da **cero en todos lados**: en
un enlace sólido-sólido el factor vale 1 + 1 − 1 = 1 pero φ(r+cᵢ) − φ(r) =
ψ_s − ψ_s = 0; en el enlace hacia el fluido el factor vale 1 + 0 − 1 = **0**,
así que ese término también se anula aunque la diferencia de potencial no sea
cero. Esto contradice al propio texto del paper, que dice (pág. 5, debajo de
la ec. 16) que ese laplaciano *"vanishes everywhere inside the electrode
except at interfacial nodes"* — leída literalmente, la ec. 15 se anula
también en los nodos interfaciales.

**La resolución está en la propia ec. 15, no en una fórmula aparte.** El
factor de la ec. 15 para el enlace dirigido r → r+cᵢ es

```
F(r, r+cᵢ) = 1 + χ_s(r+cᵢ) − χ_s(r)
```

y no es simétrico en r ↔ r+cᵢ cuando un extremo es sólido y el otro fluido.
Si se invierte el signo del término χ_s —escribiendo F'(r, r+cᵢ) =
1 − χ_s(r+cᵢ) + χ_s(r) = 1 + χ_s(r) − χ_s(r+cᵢ)— se obtiene exactamente
F'(r, r+cᵢ) = F(r+cᵢ, r): **el mismo factor de la ec. 15, mirando el mismo
enlace desde el otro nodo.** No es una ecuación nueva, es la ec. 15 aplicada
con r = nodo sólido en vez de r = nodo fluido:

| r | vecino | ec. 15, F(r, r+cᵢ) = 1+χ_s(r+cᵢ)−χ_s(r) | ec. 15 "invertida", F'(r, r+cᵢ) = 1+χ_s(r)−χ_s(r+cᵢ) |
|---|---|---|---|
| fluido, vecino pared | χ_s(r)=0, χ_s(r+cᵢ)=1 | **2** (así se dobla el peso en la matriz, sección 1.2) | 0 |
| pared, vecino fluido | χ_s(r)=1, χ_s(r+cᵢ)=0 | 0 (por eso la ec. 15 literal da carga cero) | **2** |

Es decir: la ec. 15 tal cual está escrita en el paper es la que hace falta
para **armar la matriz de Poisson** (ecuaciones de los nodos fluidos, sección
1.2); evaluada con r = nodo del electrodo da cero y no sirve para la ec. 16.
La versión con el signo de χ_s invertido es la misma ecuación vista "desde el
sólido", y es la que hay que usar para calcular la carga de la ec. 16: en un
enlace sólido-sólido vale 1 (y el término es cero porque ambos nodos están al
mismo potencial ψ_s), y en el enlace hacia el fluido vale 2, reproduciendo
exactamente lo que dice el texto del paper (cero en el electrodo salvo en los
nodos interfaciales). Lo verificamos numéricamente
(`scripts/verificar_carga_pared.py`, capacitor de placas con stencil D3Q27):

| Laplaciano en los nodos del electrodo | carga / exacta (Gauss) | Q_inf + Q_sup + Q_iones |
|---|---|---|
| ec. 15 literal (factor sin invertir) | 0 | = Q_iones (no conserva) |
| ec. 12 sin modificar | 0.5 | Q_iones / 2 |
| **ec. 15 con el signo de χ_s invertido (implementada)** | **1.0000** | **5 × 10⁻¹⁵** |

En forma de flujo por el plano medio, que es como está implementado en
`psi_petsc.c`:

```
Q_electrodo = (ε / e β) Σ_{enlaces pared w → fluido f} 2 · wlaplacian[p] · (ψ_f − ψ_w)
```

(en unidades de `rho_elec` sumado sobre nodos, con los pesos `wlaplacian` de
Ludwig, que tienen el signo opuesto a los w_i del paper). El factor 2 de esta
fórmula es el mismo F'(r, r+cᵢ) = 2 de la tabla de arriba: la resolución de
Poisson no cambia (sigue usando la ec. 15 "desde el fluido" para armar la
matriz); solo el diagnóstico de carga de la ec. 16 necesita evaluarla "desde
el sólido".

### 1.4 El campo eléctrico junto a la pared

El mismo argumento vale para el gradiente. `psi_electric_field` usaba una
diferencia central que, en el primer nodo fluido, toma ψ_s como si estuviera a
distancia h. Para un perfil lineal eso da **exactamente 3/4 del campo real**,
y el error **no disminuye al refinar la red** (se queda en 25–28%): es un
error de formulación, no de resolución.

La corrección usa el mismo factor 2, pero **en forma de diferencias**:

```
E_a(r) = − Σ_p wgradients[p] · c_pa · [ψ(r + c_p) − ψ(r)] · (1 + χ(r + c_p))
```

La forma de diferencias es obligatoria. El código original calcula
`Σ_p wgradients[p] c_pa ψ(r + c_p)`, sin restar ψ(r); eso es equivalente solo
porque Σ_p wgradients[p] c_pa = 0, y esa identidad deja de valer en cuanto se
duplica un enlace (es la misma trampa que con la diagonal del laplaciano,
sección 3.2).

Esta fórmula es el promedio de los dos flujos de cara que usa el laplaciano de
la sección 1.2 (el del enlace hacia adentro y el de la pared, a h/2), así que es
el gradiente consistente con la ley de Gauss discreta que resuelve el solver.
Se probó también una interpolación cuadrática por ψ_s, ψ₀ y ψ₁, que en papel
es de segundo orden; aplicada a la solución discreta da 4% de error y
converge a primer orden, así que se descartó (ver sección 5.3).

### 1.5 Por qué "no periódico" sin cambiar el dominio de PETSc

El DMDA de PETSc sigue siendo periódico en las tres direcciones
(`DM_BOUNDARY_PERIODIC`), igual que el sistema de coordenadas de Ludwig. La no
periodicidad la aporta la pared: una fila de pared es una fila **desacoplada**
(una sola entrada, en la diagonal), y las filas de fluido no guardan columnas
de pared (sección 3.2). Así, ninguna ecuación conecta el fluido de un lado de
la pared con el del otro lado a través del borde periódico. Basta con que la
pared tenga al menos un nodo de espesor, que es el alcance del stencil de 27
puntos; las paredes de `porous_media_init wall_z` tienen exactamente uno
(z = 1 y z = N).

Esto simplificó el plan original, que preveía cambiar el tipo de borde del
DMDA y armar filas de Dirichlet en el borde del dominio.

### 1.6 Guía de lectura del código, paso a paso

Las secciones 1.1–1.5 son la teoría; esta es la misma idea contada como una
lectura guiada del código real, para poder explicarlo sin tener el archivo
abierto al lado.

**El punto de partida.** El laplaciano discreto que usa Ludwig (consistente
con el stencil de Boltzmann de red, D3Q27) es una suma sobre los 27 vecinos
(incluido `cv[0] = (0,0,0)`, el propio nodo):

```
∇²ψ(r) = Σ_p wlaplacian[p] · ψ(r + cv[p])
```

con una convención fija ([stencil_d3q27.c:71-77](../src/stencil_d3q27.c#L71-L77)):

```c
wlaplacian[p] = -6.0 * wv[p];      // p > 0, uno por cada vecino
wlaplacian[0] = -(suma de todos los wlaplacian[p], p > 0);   // el propio nodo
```

**La diagonal es menos la suma de los vecinos.** Eso no es arbitrario: es lo
que garantiza que si ψ es constante en todo el vecindario, ∇²ψ da cero. Es
también la primera cosa que se rompe en cuanto aparece una pared, y es la raíz
de casi todos los cambios de código.

**Paso 1 — la fila de un nodo fluido junto a la pared.** El bucle real en
`psi_solver_petsc_matrix_set` ([psi_petsc.c:457-477](../src/psi_petsc.c#L457-L477)):

```c
double diag = 0.0;
int np = 1;                    // la entrada 0 es la diagonal
col[0] = row;
for (int p = 1; p < s->npoints; p++) {
  int wall = petsc_wall_site(solver, ic + s->cv[p][X], jc + s->cv[p][Y], kc + s->cv[p][Z]);
  double coef = s->wlaplacian[p] * epsilon * (1.0 + wall);   // factor 1 o 2
  diag -= coef;
  if (!wall) {
    col[np].i = i + s->cv[p][X]; col[np].j = ...; col[np].k = ...;
    v[np] = coef;
    np += 1;
  }
}
v[0] = diag;
```

Para cada uno de los 26 vecinos: `petsc_wall_site` pregunta si es
`MAP_BOUNDARY` ([psi_petsc.c:230-236](../src/psi_petsc.c#L230-L236)) — es lo
único que pregunta, no le importa el eje ni la posición. El coeficiente de ese
vínculo es `wlaplacian[p]·ε·(1+wall)`: normal si es fluido, **doble** si es
pared. La diagonal se acumula como `−Σ coef` en vez de leerse de la tabla — ya
no puede ser la de siempre, porque una vez que algún vínculo vale el doble, la
cancelación automática de la sección 1.2 se rompe. Y si el vecino es pared,
**su columna directamente no se agrega** (`if (!wall)`): esa es la
"eliminación de columna" del paso 3.

**Paso 2 — la fila de un nodo de pared.** Es trivial
([psi_petsc.c:452-456](../src/psi_petsc.c#L452-L456)): una sola entrada, en la
diagonal, valor `wall_diag = ε·wlaplacian[0]` (del mismo orden que una
diagonal de fluido, para que el precondicionador no vea números dispares). Esa
fila dice "ψ = algo", y ese "algo" lo pone el lado derecho (paso 3).

**Paso 3 — por qué se elimina la columna, y adónde va.** Esto no sale de la
física sino de una restricción práctica: el `.petscrc` usa CG + BoomerAMG, y
CG solo funciona con matrices **simétricas**. Si la fila de un nodo fluido
guardara la columna hacia su vecino de pared, pero la fila de esa pared fuera
una identidad aislada sin la columna simétrica de vuelta, la matriz dejaría de
ser simétrica. La solución es no guardar esa columna — el término no
desaparece, se **traslada al lado derecho**, porque el valor de ψ en la pared
ya se conoce de antemano (es un dato, no una incógnita). Es el paso clásico de
separar A·x=b en incógnitas y datos conocidos. Lo hace
`psi_solver_petsc_rhs_set` ([psi_petsc.c:638-660](../src/psi_petsc.c#L638-L660)):

```c
if (petsc_wall_site(solver, ic, jc, kc)) {
  rho_3d[k][j][i] = solver->block->wall_diag * petsc_wall_value(solver, ic, jc, kc);
  continue;                                    // fila de pared: b = D·ψ_pared
}
rho_3d[k][j][i] = rho_elec * eunit * beta;      // fila de fluido, como siempre
for (cada vecino p)
  if (ese vecino es pared)
    rho_3d[k][j][i] -= 2.0 * wlaplacian[p] * epsilon * ψ_pared_del_vecino;
```

El resultado (verificado, `max|A−Aᵀ| = 0`) es matemáticamente equivalente a
haber dejado la columna — solo que ahora vive en `b` en lugar de en `A`.

**Paso 4 — por qué desaparece el espacio nulo.** Sin pared, el laplaciano
periódico tiene un modo nulo: sumarle una constante a ψ en todos lados no
cambia ∇²ψ, así que hay que sacárselo a la matriz (`MatSetNullSpace`) o CG no
converge (infinitas soluciones ψ, ψ+c, ψ+2c...). Con una pared a potencial
fijo esa ambigüedad se termina: ψ=0 en la pared ya elige una solución única.
Por eso ([psi_petsc.c:508-510](../src/psi_petsc.c#L508-L510)) con pared se
pasa `NULL` en vez del nullspace de modo constante.

**Paso 5 — cómo se prende todo (`psi_solver_petsc_wall_set`,
[psi_petsc.c:1095-1152](../src/psi_petsc.c#L1095-L1152)).** En orden: valida
(solver correcto, sin `e0`, con nodos de pared) → refresca el halo del mapa
(`map_halo`, porque el bucle de arriba pregunta por vecinos que pueden estar
en otra celda o al otro lado de la caja) → guarda eje, potenciales y
`wall_diag` → **copia el potencial de pared directamente en el campo ψ** en
esos nodos (así cualquier código que lea ψ antes de resolver ve algo
sensato) → `solver->psi->wall_map = map` (le avisa a `psi_t` para que
`psi_electric_field`, que vive en otro archivo, se entere) → llama a
`psi_solver_petsc_matrix_set` **una sola vez** (las paredes no se mueven, así
que la matriz no hace falta rehacerla en cada paso; solo el lado derecho
cambia, porque la carga sí cambia).

**Paso 6 — el mismo truco para el campo eléctrico.** E = −∇ψ tiene el mismo
problema que el laplaciano: tomar ψ de la pared como si estuviera a distancia
h en vez de h/2 da un gradiente mal calculado. La solución, en
`psi_electric_field` ([psi_gradients.c:69-92](../src/psi_gradients.c#L69-L92)),
es la misma idea con la misma trampa:

```c
double psi00 = psi->psi->data[...index...];
for (cada vecino p) {
  double fac = (vecino es pared) ? 2.0 : 1.0;
  e[X] -= wgradients[p] * cx * (psi_vecino - psi00) * fac;   // idem Y, Z
}
```

Tiene que ser en **forma de diferencia** (`psi_vecino − psi00`), no la suma
directa `Σ wgradients[p]·cx·ψ_vecino` que usaba el código original. En bulk
ambas dan lo mismo porque Σ wgradients[p]·cx = 0 por simetría del stencil —
pero esa cancelación se rompe en cuanto un vínculo se duplica, exactamente el
mismo problema que con la diagonal del laplaciano en el paso 1. Con la forma
de suma directa el campo sale al 75% junto a la pared y **no mejora al
refinar la malla** (sección 5.3): es un error de fórmula, no de resolución.
En un nodo de pared, la función devuelve directamente 0 — adentro del
conductor no hay campo.

**En una frase:** el laplaciano de Ludwig calcula ∇²ψ como una diferencia
pesada con cada vecino; la pared está a mitad de camino entre el último nodo
de fluido y el nodo sólido, así que cada vínculo hacia la pared vale el doble
de lo normal (Asta et al. 2019). Como eso duplica algunas columnas, hay que
(1) reconstruir la diagonal a mano en cada nodo, porque la cancelación
automática se rompe, (2) sacar esas columnas de la matriz para que siga siendo
simétrica —CG lo exige— y pasarlas al lado derecho, donde el valor de la pared
ya es un dato conocido, y (3) aplicar el mismo factor 2, en forma de
diferencia, al gradiente para que la fuerza eléctrica salga bien.

**Mapa de quién llama a quién:**

```
psi_solver_petsc_wall_set()                    ← una vez, al arrancar (ludwig_walls_attach, sección 3.5)
 ├─ map_halo(map)
 ├─ guarda wall_axis, wall_psi, wall_diag
 ├─ copia ψ_pared en los nodos de pared
 ├─ psi->wall_map = map                        ← para que psi_electric_field se entere
 └─ psi_solver_petsc_matrix_set()               ← arma la matriz UNA VEZ (factor 2, columnas eliminadas)

en cada paso de tiempo:
 psi_solver_petsc_solve()
  ├─ psi_solver_petsc_rhs_set()                 ← recalcula b (carga + columnas eliminadas)
  ├─ KSPSolve()                                 ← resuelve A x = b (A fija, b nuevo)
  ├─ psi_solver_petsc_wall_charge_compute()      ← carga inducida en cada pared (diagnóstico)
  └─ psi_electric_field()  (desde otro archivo)  ← usa psi->wall_map, mismo factor 2, en forma de diferencia
```

Las funciones auxiliares `petsc_wall_site`, `petsc_wall_side`,
`petsc_wall_value` ([psi_petsc.c:230-261](../src/psi_petsc.c#L230-L261)) son
las únicas que conocen la geometría (qué es pared, y qué potencial le toca);
todo lo demás —matriz, lado derecho, campo eléctrico— es genérico y solo le
pregunta a `petsc_wall_site` "¿sos pared?".

---

## 2. Cómo se usa

```
porous_media_init               wall_z        # paredes MAP_BOUNDARY en z=1 y z=N
electrokinetics_solver_type     petsc
electrokinetics_wall_axis       z             # x | y | z
electrokinetics_wall_potential  0.0_1.0e-6    # potencial pared inferior _ superior
```

- **Unidades del potencial:** las del campo ψ de Ludwig, que es el **potencial
  reducido** βeψ (en unidades de kT/e). Verificado en la prueba B: los iones
  en equilibrio cumplen ln(ρ₊/ρ₋) = −2.000000·ψ.
- **Qué pared tiene qué potencial:** un nodo de pared en la mitad inferior de
  la caja (a lo largo de `electrokinetics_wall_axis`) toma el primer valor, y
  uno en la mitad superior, el segundo.
- **Salida:** con la frecuencia `freq_psi_resid`, el log muestra

  ```
  Wall charge lower -1.985915490356800e-07 upper  1.985915488622381e-07
  Fluid charge       1.734723475976807e-16 total  3.051965588075651e-20
  ```

  La carga está en unidades de `rho_elec` sumado sobre nodos; `total` es
  Q_inf + Q_sup + Q_fluido, que tiene que dar cero (control de neutralidad).
- **Sin estas claves**, el solver es el periódico original.
- **Electroneutralidad (opcional con paredes):**
  `electrokinetics_electroneutral no` (por defecto `yes`) evita que
  `psi_electroneutral` agregue los contraiones que neutralizan las partículas.
  Con paredes a potencial fijo Poisson tiene solución única para cualquier
  carga neta y las paredes aportan la carga que compensa, así que una carga
  aislada entre placas es un problema bien planteado. Sin paredes se rechaza:
  los solvers periódicos restan la carga media, que es lo mismo que un fondo
  neutralizante.
- **Corrección de corto alcance:** `pm_sr_pair_ref measured` funciona con
  paredes. La tabla se calibra con el solver **periódico**, antes de enganchar
  las paredes (sección 3.5), así que es la tabla de bulk y el caché
  `pm_sr_meshref_*.dat` sigue siendo válido para corridas periódicas. No
  describe la interacción de corto alcance partícula–pared (imagen); eso es
  la fase 3.

Combinaciones rechazadas con un mensaje de error (no se ignoran en silencio):

| Combinación | Motivo |
|---|---|
| solver distinto de `petsc` | solo está implementado en PETSc |
| permitividad variable (`fe_electro_symmetric` con ε₁ ≠ ε₂) | no implementado en ese camino |
| `electric_e0` ≠ 0 | el término de campo externo del RHS supone periodicidad |
| mapa sin nodos `MAP_BOUNDARY` | no hay paredes |
| `electrokinetics_electroneutral no` sin paredes | la carga neta no tiene sentido en un sistema periódico |

---

## 3. Modificaciones por archivo

### 3.1 `src/psi_petsc.h` — interfaz pública

Dos funciones nuevas:

```c
int psi_solver_petsc_wall_set(psi_solver_petsc_t * solver, map_t * map,
                              int axis, const double psi_wall[2]);
int psi_solver_petsc_wall_charge(const psi_solver_petsc_t * solver,
                                 double q[2], double * q_fluid);
```

`wall_set` activa las paredes; `wall_charge` devuelve la carga de cada pared y
del fluido de la última resolución.

### 3.2 `src/psi_petsc.c` — el solver

**Estado interno** (`struct psi_solver_petsc_block_s`): puntero al mapa (NULL
= periódico), eje, potenciales de las dos paredes, diagonal de las filas de
pared, y las cargas de la última resolución. Como el bloque se aloca con
`calloc`, todo arranca en cero y el comportamiento por defecto es el original.

**Funciones auxiliares** (`petsc_wall_site`, `petsc_wall_side`,
`petsc_wall_value`): si un nodo (coordenadas locales, halo incluido) es pared,
a qué mitad de la caja pertenece y qué potencial le corresponde.

**Matriz** (`psi_solver_petsc_matrix_set`, camino de ε uniforme):

- *Fila de pared:* una sola entrada, `A[r,r] = D` con `D = ε · wlaplacian[0]`.
  Se escala como una diagonal de fluido para que el precondicionador vea filas
  comparables.
- *Fila de fluido:* para cada vecino p,
  `coef = wlaplacian[p] · ε · (1 + χ(vecino))`. Si el vecino es fluido, el
  coeficiente se guarda en su columna; si es pared, **la columna se elimina**
  (su término pasa al lado derecho). La diagonal se reconstruye nodo por nodo
  como `−Σ coef`: ya no se puede leer de `wlaplacian[0]`, porque esa tabla
  supone que todos los enlaces valen 1.
- *Por qué eliminar las columnas:* el `.petscrc` usa CG + BoomerAMG, que exigen
  una matriz **simétrica**. Con filas de pared identidad y filas de fluido que
  guardan la columna de pared, la matriz no es simétrica. Eliminando esas
  columnas, el bloque fluido–fluido queda simétrico y definido positivo.
- La matriz se arma **una sola vez**, en `wall_set` (las paredes no se mueven).
  Antes de rearmarla se llama a `MatZeroEntries`: la estructura de no-ceros
  está congelada (`MAT_NEW_NONZERO_LOCATIONS` falso) y una fila que ahora
  escribe menos entradas conservaría las viejas.
- *Nullspace:* con Dirichlet el operador ya no es singular, así que se quita
  (`MatSetNullSpace(a, NULL)`). Sin paredes se adjunta el de modo constante,
  como antes.

**Lado derecho** (`psi_solver_petsc_rhs_set`, se rearma en cada resolución):

- *Nodo de pared:* `b = D · ψ_pared` (la fila es `D x = b`).
- *Nodo fluido:* `b = rho_elec · e · β` como antes, y por cada vecino de pared
  se resta la columna eliminada: `b −= 2 · wlaplacian[p] · ε · ψ_pared`.
- Se acumula la carga total del fluido para el control de neutralidad.

**Resolución** (`psi_solver_petsc_solve`): la proyección explícita del
nullspace sobre `b` y `x` (`MatNullSpaceRemove`, agregada el 2026-03-26 para
CG en GPU) **no se hace con paredes**: no hay nullspace, y restar la media
desplazaría el potencial impuesto. Después de resolver se calcula la carga
inducida y se reporta.

**Carga inducida** (`psi_solver_petsc_wall_charge_compute`): suma de flujos por
el plano medio con la regla simétrica (sección 1.3), leyendo la solución desde
un vector local con fantasmas del DMDA y reduciendo con MPI. Se saltean los
enlaces pared–pared: a través del borde periódico unen los dos electrodos, y
eso no es un enlace físico.

**Activación** (`psi_solver_petsc_wall_set`): valida (ε uniforme, sin `e0`,
con nodos de pared), actualiza el halo del mapa, guarda los parámetros, copia
ψ_pared al campo ψ en los nodos de pared (valor inicial, y lo que ve cualquier
código que lea ψ antes de la primera resolución), le pasa el mapa a `psi` para
el campo eléctrico, rearma la matriz y vuelve a preparar el KSP.

Sin PETSc compilado, las dos funciones nuevas son stubs que devuelven −1.

### 3.3 `src/psi.h` — `psi_t`

Nuevo miembro al final del struct:

```c
map_t* wall_map;   /* MAP_BOUNDARY sites are walls; NULL = no walls */
```

Lo fija `psi_solver_petsc_wall_set`. `psi_t` se aloca con `calloc`, así que
sin paredes vale NULL.

### 3.4 `src/psi_gradients.c` — `psi_electric_field`

Si `psi->wall_map` es NULL, corre el código original sin cambios. Si no:

- en un nodo de pared devuelve **E = 0** (interior del conductor; antes
  devolvía un valor que mezclaba los dos electrodos a través del borde
  periódico);
- en un nodo fluido usa la fórmula de la sección 1.4.

Esta función es la única que alimenta el campo en todas las rutas que usa
PETSc: fuerza sobre el fluido (`psi_force.c`), fuerza sobre coloides
(`nernst_planck.c`), fuerza subgrid (`subgrid.c`) y la salida `efield`. El
solver FFT tiene su propio gradiente (`psi_fft_pn.c`), pero nunca se usa con
paredes.

### 3.5 `src/ludwig.c`

**Lectura y enganche de las paredes.** En `ludwig_rt`, después de
`psi_rt_init_rho`, solo se leen y validan `electrokinetics_wall_potential` y
`electrokinetics_wall_axis` (solver PETSc, eje válido). El enganche real
(`ludwig_walls_attach`, que llama a `psi_solver_petsc_wall_set`) ocurre en
`ludwig_run`, **justo antes del loop temporal y después de construir las
tablas de `pm_sr`**. Hay dos razones para ese orden:

- el solver de Poisson se crea en `free_energy_init_rt`, **antes** que el mapa
  (`map_init_rt`), y la inicialización de ψ todavía puede modificar el mapa;
- la referencia medida de `pm_sr` se calibra con el mismo objeto
  `ludwig->poisson`. Si las paredes ya estuvieran puestas, la tabla quedaría
  medida con paredes en lugar de en bulk periódico, y se guardaría en un caché
  (`pm_sr_meshref_k*_n*_L*_rc*_nr*.dat`) cuyo nombre no registra paredes, así
  que una corrida periódica posterior la reutilizaría contaminada.

Una versión anterior resolvía el segundo punto rechazando
`pm_sr_pair_ref measured` con paredes; con el enganche diferido ya no hace
falta.

**Electroneutralidad opcional.** En la inicialización (`ntstep == 0`),
`psi_electroneutral` se saltea con `electrokinetics_electroneutral no`, que
solo se acepta si hay paredes.

### 3.6 `src/subgrid.c`

**Corrección partícula–nodo junto a las paredes**
(`pm_sr_apply_force_correction`, nueva función `pm_sr_node_behind_wall`). Esa
corrección recorre los nodos hasta r_cut (~6) desde la partícula, más allá
del soporte del kernel. Con paredes saltea los nodos de pared y los nodos
fluidos alcanzados **a través del borde periódico** cuando ese borde es una
pared: su copia en el halo es el fluido del otro lado del electrodo. Con el
rango actual, que se recorta al dominio local `[1, nlocal]`, esos nodos del
otro lado ni siquiera se visitan, así que en la práctica el filtro actúa sobre
los nodos de pared; el chequeo del borde queda por si el rango se amplía.

**Arreglo de una fuerza duplicada (error previo, no relacionado con las
paredes)** en `subgrid_update_forces_electrokinetics`. El loop recorre también
las celdas halo, que contienen copias periódicas de las partículas de la celda
de coloides más externa. Para el reparto de la reacción sobre el fluido eso
está bien, porque cada copia reparte solo en sus nodos locales. Pero cada copia
sumaba la fuerza **completa** a su `fex`, y la suma de halos
(`colloid_sums.c`) se la agregaba a la partícula original: una partícula
cargada en una celda de cara recibía 2× la fuerza electrostática, en una
arista 4× y en una esquina 8×. Ahora solo la copia dueña (celdas 1…ncell)
acumula. Se detectó en la prueba C (fuerza·β/E = 2 exacto para z₀ ≤ 5, con
celdas de 4.86). Revisión de corridas anteriores: ninguna corrida de microgel
ni la batería de fuerzas entre 2 partículas tuvo partículas cargadas en celdas
de borde; sí las pruebas `CROSS-borde` (fuerzas absolutas duplicadas durante
toda la corrida, así que la continuidad al cruzar el borde, que era lo que se
medía, no se ve afectada), `old/` y las `PHASE0-*` (camino Ewald, que usa otra
función; sin verificar si le afecta).

---

## 4. Estructura del sistema lineal (resumen)

Para una fila de fluido r con vecinos p:

| | valor |
|---|---|
| A[r, vecino fluido] | ε · wlaplacian[p] |
| A[r, vecino pared] | 0 (columna eliminada) |
| A[r, r] | −Σ_p ε · wlaplacian[p] · (1 + χ_p) |
| b[r] | e β · rho_elec(r) − Σ_{p pared} 2 ε · wlaplacian[p] · ψ_p |

Para una fila de pared: A[r, r] = ε · wlaplacian[0], b[r] = ε · wlaplacian[0] · ψ_pared.

---

## 5. Validación

Todas las pruebas usan la configuración real de GPU (`.petscrc`: CG +
BoomerAMG, `aijcusparse`). Caja 4 × 4 × 34 con paredes en z = 1 y z = 34, así
que los planos medios quedan en z = 1.5 y z = 33.5 (L = 32).
Directorio: `capacitor/local/WALL-DIRICHLET/`.

### 5.1 Prueba A — capacitor sin iones (`A-sin-iones/`)

ψ_inf = 0, ψ_sup = 10⁻⁶, un paso. Solución exacta: perfil lineal entre planos
medios.

| Magnitud | Resultado |
|---|---|
| perfil ψ(z) | exacto, error relativo 2 × 10⁻¹⁴ |
| primer nodo fluido | 1.5625 × 10⁻⁸ (con la pared en el nodo sólido sería 3.03 × 10⁻⁸) |
| carga en cada pared | ± 5.000000000000 × 10⁻⁸ = ε·E·A/(eβ) exacto |
| convergencia de CG | 16 iteraciones, residuo 2.6 × 10⁻¹⁹ |
| campo E_z en el fluido | 1.000000000000 × exacto en **todos** los nodos (antes 0.75 junto a la pared) |
| campo tangencial | 10⁻²⁵ |
| campo en nodos de pared | 0 (antes 4.9 × 10⁻¹²) |

### 5.2 Prueba B — capacitor con electrolito (`B-debye-huckel/`)

λ_D = 4 (ρ = 3.125 × 10⁻³ por especie), D = 0.05, `hydrodynamics off`, 60000
pasos hasta el estacionario. Referencias: la solución de Debye-Hückel
continua (ec. 17 del paper) y la solución de Debye-Hückel **discreta** en la
misma red con la misma regla de pared. La segunda separa el error de
implementación del error de discretización. Script: `scripts/analizar_B.py`.

| Magnitud | Resultado |
|---|---|
| residuo de Poisson nodo a nodo (regla de pared) | 8.8 × 10⁻¹⁴ sobre escala 2.7 × 10⁻⁴ |
| equilibrio de Boltzmann | ln(ρ₊/ρ₋) = −2.000000·ψ, lineal a 10⁻¹⁴ |
| perfil vs DH discreta | 8 × 10⁻⁹ |
| perfil vs DH continua | 3.3 × 10⁻³ (error de red con κh = 0.25) |
| carga vs DH continua | 0.99229 |
| neutralidad Q_inf + Q_sup + Q_fluido | 3 × 10⁻²⁰ |
| campo junto a la pared vs DH continua | 1.00034 en ambas paredes (antes 0.719); 0 en los nodos de pared |

Con `hydrodynamics off` el campo no realimenta la dinámica, así que el
potencial y las densidades son idénticos a los de la corrida anterior a la
corrección del campo (diferencia 10⁻²¹); solo cambia el campo.

### 5.3 Orden de convergencia

Como Ludwig reproduce la DH discreta a 10⁻⁸, el orden del esquema se mide
refinando esa solución discreta con κL = 8 fijo
(`scripts/convergencia_campo.py` y el análisis del perfil).

| | error del perfil al refinar | carga / continua |
|---|---|---|
| pared en el plano medio (implementado) | segundo orden (razón 3.65 → 3.92) | 0.9923 → 0.9999 |
| pared en el nodo sólido (ingenuo) | primer orden (razón 1.82 → 1.95), 16× peor | 0.883 → 0.984 |

Campo en el nodo junto a la pared:

| Esquema | error relativo (h = 1 … 1/16) | razón al refinar |
|---|---|---|
| diferencia central (antes) | 28% → 25% | ≈ 1.0 (no converge) |
| **factor 2 en diferencias (implementado)** | 3.4 × 10⁻⁴ → 1.3 × 10⁻⁷ | ≈ 7 |
| interpolación cuadrática (descartada) | 4.2% → 0.26% | 2.0 |

### 5.4 Prueba C — partícula subgrid cargada entre paredes a tierra (`C-*`)

Una partícula subgrid con q = +1, sin sal y sin corrección `pm_sr` (PM puro),
en una caja 32 × 32 × 34 con las dos paredes a potencial 0, a distancias d de
la pared inferior entre 0.5 y 16 (un paso por posición). Teoría (detalle en el
`README.md` de la carpeta): carga inducida por Shockley–Ramo, y campo de
imagen como serie de Fourier en x,y con Dirichlet en z, incluyendo las
imágenes laterales de la caja periódica. El fondo de contraiones que agrega
`psi_electroneutral` cancela exactamente el término k = 0. Script:
`scripts/analizar_C.py`.

| | Hann-8 | Hann-6 |
|---|---|---|
| distancia mínima sin tocar la pared, d ≥ n/2 − ½ | 3.5 | 2.5 |
| carga inducida vs Shockley–Ramo | 1.2 × 10⁻¹⁴ | 1.5 × 10⁻¹⁴ |
| carga perdida en la pared | 3 × 10⁻¹⁴ | 3 × 10⁻¹⁴ |
| campo de imagen: error en la distancia mínima | 1.5% | 1.0% |
| campo de imagen: error mediano | 0.03% | 0.12% |
| partícula en el centro (d = 16), teoría 0 | 9 × 10⁻¹⁸ | 1 × 10⁻¹⁷ |

Cuando el soporte del kernel alcanza los nodos de pared, la carga que cae ahí
se pierde (la fila de pared impone el potencial y descarta el lado derecho):
hasta 16% de la carga con Hann-8 en d = 1.5, y el campo sale subestimado hasta
en un 70%. Esto confirma la distancia mínima d_min = n/2 − ½ prevista en el plan
(fase 2) y es la razón de la capa de sólido a distancia fija junto a los
electrodos en `util/porous_solid_cube`.

**Sin fondo de contraiones y con corrección** (Hann-8, figura
`graficos/C_comparacion.png`):

| batería | fondo | corrección | error del campo en d_min | mediana | carga vs teoría |
|---|---|---|---|---|---|
| `C-carga-imagen-hann08` | sí | no | 1.5% | 3.5 × 10⁻⁴ | 1 × 10⁻¹⁴ |
| `C-corr-hann08` | sí | sí | 1.5% | 3.5 × 10⁻⁴ | 1 × 10⁻¹⁴ |
| `C-sin-fondo-hann08` | no (carga neta +1) | no | 1.2% | 5 × 10⁻⁵ | 4 × 10⁻¹⁴ |
| `C-sin-fondo-corr-hann08` | no | sí | 1.2% | 5 × 10⁻⁵ | 4 × 10⁻¹⁴ |

Sin fondo la teoría conserva el término k = 0 y la carga inducida total es
−q; el error es menor porque el fondo discreto de contraiones agrega el suyo.
La corrección no cambia el resultado de una sola partícula: sin iones es
idéntico bit a bit y con el fondo difiere en 10⁻⁷, porque corrige pares
partícula–partícula y partícula–iones y la interacción con la imagen no está
en sus tablas.

El error residual de 1–1.5% en la distancia mínima es un efecto de corto
alcance entre la partícula y su imagen (a 2d ≈ 5–7, del orden del ancho del
kernel); a partir de d ≈ 5 cae por debajo de 0.3%. Es el tipo de error que
la corrección medida de corto alcance tendría que cubrir cerca de las
paredes (fase 3, pendiente).

### 5.5 Regresión periódica (`regresion-periodica/`)

Un paso de una corrida periódica anterior con PETSc (una partícula, Hann-8,
λ_D = 8, con corrección de corto alcance), con el binario nuevo: diferencia
**exactamente cero** en ψ, en las densidades iónicas y en el campo, antes y
después de modificar `psi_electric_field`, y otra vez después de diferir el
enganche de las paredes y de arreglar la fuerza duplicada (la partícula de esa
prueba está en el centro de la caja, fuera de las celdas de borde).

### 5.6 Calibración de `pm_sr` antes de las paredes (`C-corr-hann08/`)

El `run.log` de las corridas con corrección muestra el orden: primero
"measuring mesh reference (256 Poisson solves)" y "mesh reference written",
después "Constant-potential walls". El fondo de contraiones sobrevive a la
calibración (carga del fluido 3 × 10⁻¹⁴). La antigua prueba de rechazo quedó
obsoleta (`guarda-pmsr-obsoleta/NOTA.txt`).

---

## 6. Limitaciones y pendientes

**No verificado:**

- **MPI con más de un proceso.** El build de GPU exige una GPU por proceso y la
  máquina de desarrollo tiene una. El código está escrito para varios procesos
  (halo del mapa, fantasmas del DMDA, reducciones con MPI), pero hay que
  probarlo, por ejemplo en MareNostrum 5.
- **Hidrodinámica con paredes y electrocinética juntas.** Las pruebas usan
  `hydrodynamics off`. `porous_media_init wall_z` marca el mapa como medio
  poroso, y `wall_rt_init` arma el bounce-back para ese caso, pero la
  combinación no se probó.
- **Nernst-Planck cerca de la pared** (pendiente anotado). Los flujos iónicos
  saltean los vecinos que no son fluido (`nernst_planck.c`), así que NP nunca
  lee ψ de la pared, y el perfil de Boltzmann en equilibrio salió exacto. Lo
  que falta verificar es la **dinámica** de carga cerca de la pared.

**Limitación previa observada (sistemas periódicos, sin arreglar):** el rango
de nodos de la corrección partícula–nodo se recorta al dominio local y solo
la recorre la copia dueña de la partícula. Una partícula a menos de r_cut
(~6) del borde de la caja no recibe la corrección de los iones que están del
otro lado del borde periódico. No afecta al microgel (lejos de los bordes);
sí a sólidos que llenan la caja, como el del capacitor, cuando tengan carga.

**No implementado:**

- permitividad variable (`wall_set` devuelve −1);
- campo externo `electric_e0` con paredes (en Ludwig ese campo entra como un
  salto de potencial en el RHS, incompatible con electrodos metálicos; para
  electroósmosis tangencial habría que aplicarlo como fuerza de volumen);
- electrodos con potencial distinto a lo largo de la pared (hoy: mitad
  inferior / mitad superior de la caja);
- paredes que cambian durante la corrida (la matriz se arma una vez);
- corrección de corto alcance `pm_sr` cerca de paredes y partículas subgrid
  cuyo kernel alcanza la pared (fases 2–4 del plan).

**Validación pendiente sugerida:** la dinámica de carga del capacitor
(Fig. 4 del paper, carga en función del tiempo con tiempo característico
λ_D·L/2D) y el perfil electroosmótico estacionario (con el campo tangencial
como fuerza de volumen).

---

## 7. Cambios relacionados anteriores

No son parte de la condición de pared, pero son necesarios para usar PETSc en
estos sistemas:

- `src/psi_petsc.c`, 2026-09-21: ancho de stencil del DMDA = 1 en lugar de
  `nhalo`. Con `nhalo` (3–4 con kernels Hann/Peskin-6) se preasignaban 343–729
  entradas por fila para un operador de 27 puntos, y hypre en GPU abortaba por
  falta de memoria.
- `src/subgrid.c`, 2026-09-21: en `pm_sr_mesh_pair_esub_x` se eliminó una
  copia device→host posterior a la resolución que pisaba con ceros la solución
  de PETSc al medir la referencia de corto alcance.
- `util/porous_solid_cube.c`: opciones `--wall-layer`, `--wall-distance` y
  `--wall-layer-density` para generar el sólido poroso del capacitor con una
  capa de partículas a distancia mínima fija de cada pared, ancladas solo en
  la dirección normal (preparación de la fase 2: espacio libre de sólido junto
  al electrodo, sin renormalizar el kernel).
