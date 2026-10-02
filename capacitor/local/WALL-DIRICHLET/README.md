# Pruebas de las paredes a potencial constante (PETSc)

Todas las pruebas que verifican que las paredes a potencial fijo funcionan
correctamente, con sus entradas, salidas, scripts de análisis y gráficos.

La descripción completa de la implementación (física, ecuaciones, cambios en
el código y limitaciones) está en
[claude/PAREDES_POTENCIAL_CONSTANTE.md](../../../claude/PAREDES_POTENCIAL_CONSTANTE.md).

Para reproducir todo: `./correr_tests.sh` (o `./correr_tests.sh A B` para
algunas). Los gráficos se regeneran con `scripts/graficos_paredes.py`.

---

## Configuración común

Capacitor de placas: caja 4 × 4 × 34, paredes `MAP_BOUNDARY` en z = 1 y
z = 34 (`porous_media_init wall_z`), periódico en x e y. Los planos medios,
donde queda la pared física, están en z = 1.5 y z = 33.5, así que la distancia
entre electrodos es L = 32. Solver PETSc con la configuración real de GPU
(`.petscrc`: CG + BoomerAMG, `aijcusparse`). Potencial de las paredes 0 y
10⁻⁶ (unidades reducidas kT/e). Plantilla de entrada: `input.base`.

---

## Contenido

| Carpeta | Qué verifica | Resultado |
|---|---|---|
| `A-sin-iones/` | Poisson con paredes sin carga (1 paso): perfil lineal, carga inducida, campo | perfil exacto (2 × 10⁻¹⁴); carga = Gauss; campo exacto en todos los nodos |
| `B-debye-huckel/` | electrolito con λ_D = 4 hasta el estacionario (60000 pasos, sin hidrodinámica) | = Debye-Hückel discreta a 8 × 10⁻⁹; Boltzmann exacto; campo junto a la pared 1.00034 |
| `C-carga-imagen-hann08/`, `C-carga-imagen-hann06/` | una partícula subgrid q = +1 entre paredes a tierra (caja 32 × 32 × 34, sin sal), en función de la distancia d a la pared; 1 paso por posición | carga inducida = Shockley–Ramo a 10⁻¹⁴; campo de imagen vs serie de Fourier: ≤ 1.5% en el borde permitido, mediana 0.03% (Hann-8) y 0.1% (Hann-6); si el kernel toca la pared se pierde carga y el campo sale mal |
| `C-sin-fondo-hann08/` | lo mismo **sin** fondo de contraiones (`electrokinetics_electroneutral no`): carga neta +1, las paredes aportan −1 | carga = Shockley–Ramo sin fondo a 10⁻¹⁴ (total −1); campo: 1.2% en el borde permitido, mediana 5 × 10⁻⁵ |
| `C-corr-hann08/`, `C-sin-fondo-corr-hann08/` | las dos anteriores **con** corrección `pm_sr` (referencia medida, calibrada en periódico antes de enganchar las paredes) | sin iones, idéntico a sin corrección; con el fondo, difiere en 10⁻⁷. Con una sola partícula la corrección no tiene sobre qué actuar: corrige pares partícula–partícula y partícula–iones, no la imagen |
| `regresion-periodica/` | sin las claves nuevas el solver es el original (1 paso de una corrida vieja con PETSc) | diferencia exactamente 0 en ψ, iones y campo |
| `guarda-pmsr-obsoleta/` | verificaba que se rechazara `pm_sr_pair_ref measured` con paredes; ya no aplica (ver `NOTA.txt`) | — |
| `A-mpi2-no-ejecutada/` | la prueba A con 2 procesos MPI | **no se pudo correr**: el build de GPU exige una GPU por proceso (ver `NOTA.txt`) |
| `scripts/` | análisis y pruebas numéricas en Python | ver abajo |
| `graficos/` | figuras de todas las pruebas | ver abajo |

Las subcarpetas `previo/` dentro de A y B guardan el campo eléctrico
**antes** de corregir `psi_electric_field`, para los gráficos antes/después.

### Scripts

| Script | Qué hace |
|---|---|
| `verificar_carga_pared.py` | resuelve el capacitor en Python con el stencil D3Q27 y compara tres fórmulas para la carga inducida (ec. 16 del paper); también verifica que la matriz quede simétrica |
| `analizar_B.py` | sobre la salida de B: residuo de Poisson nodo a nodo, equilibrio de Boltzmann, perfil vs Debye-Hückel discreta y continua, campo junto a la pared |
| `convergencia_campo.py` | orden de convergencia del perfil y del campo junto a la pared, refinando la solución discreta de Debye-Hückel |
| `graficos_paredes.py` | genera las figuras de A, B, cargas y convergencia |
| `analizar_C.py` | prueba C: tabla contra la teoría y figura de una batería (`analizar_C.py <carpeta> <n>`); detecta solo si hubo fondo y corrección |
| `comparar_C.py` | error del campo de las cinco baterías C en una sola figura |

### Gráficos (`graficos/`)

| Figura | Contenido |
|---|---|
| `A_capacitor_sin_iones.png` | (a) ψ(z) contra el perfil exacto; (b) zoom junto a la pared, con la recta que daría una pared en el nodo sólido; (c) campo antes (0.75 junto a la pared) y después (1.000) |
| `B_debye_huckel.png` | (a) ψ(z) contra Debye-Hückel continua y discreta; (b) desviación de cada una (el error de implementación es 10⁻⁸; el de red, 10⁻³); (c) densidades iónicas contra Boltzmann; (d) campo antes (0.719) y después (1.00034) |
| `carga_inducida_reglas.png` | carga inducida y neutralidad con cada laplaciano posible en los nodos del electrodo: la ec. 15 literal da 0, la ec. 12 da 0.5, la regla simétrica da 1 y conserva la carga a 10⁻¹⁵ |
| `C_carga_imagen_hann08.png`, `C_carga_imagen_hann06.png`, `C_sin_fondo_hann08.png`, `C_corr_hann08.png`, `C_sin_fondo_corr_hann08.png` | (a) campo sobre la partícula contra la teoría de imágenes (y, como referencia, una pared sola sin imágenes laterales); (b) error relativo; (c) carga inducida en cada pared contra Shockley–Ramo; (d) error de carga y carga perdida en la pared |
| `C_comparacion.png` | error relativo del campo de imagen de las cinco baterías C juntas |
| `convergencia.png` | (a) error del perfil: segundo orden con la pared en el plano medio, primer orden con la pared en el nodo; (b) error del campo junto a la pared: la diferencia central no converge, la cuadrática es de primer orden, la implementada converge rápido |

---

## Prueba C: teoría

Caja periódica en x,y de lado ℓ = 32 (área A), paredes a tierra en los planos
medios separadas L = 32, carga q a distancia d de la pared inferior. En
unidades de Ludwig (potencial reducido, 4π l_B = β/ε = 10):

- **Carga inducida** (Shockley–Ramo): −q(L−d)/L en la inferior y −q·d/L en la
  superior, más la del fondo de contraiones que agrega `psi_electroneutral`
  (uniforme en el fluido salvo el nodo pegado a la partícula, que queda en
  cero; incluido exactamente). Es independiente del kernel mientras no toque
  la pared.
- **Campo de imagen** (Fourier en x,y, Dirichlet en z):
  E_z = −(4π l_B q/A)·[(L−2d)/(2L) + Σ_{k≠0} sinh(k(L−2d))/(2 sinh kL)],
  k = 2π|(m,n)|/ℓ. El fondo uniforme cancela exactamente el término k = 0.
  Para una carga de simetría esférica que no toca la pared la imagen actúa
  como sobre una carga puntual, así que la teoría puntual vale para el kernel.
- **Distancia mínima:** Hann-n tiene peso nulo para |x| ≥ n/2, así que no
  deposita carga en la pared si d ≥ n/2 − ½ (3.5 para Hann-8, 2.5 para
  Hann-6). Las posiciones más cercanas se corrieron a propósito para mostrar
  el efecto.

### Con y sin fondo de contraiones

Por defecto `psi_electroneutral` agrega contraiones uniformes que neutralizan
la partícula, así que el fluido es neutro y las cargas de las paredes suman
cero. Con paredes a potencial fijo eso no hace falta: Poisson tiene solución
única para cualquier carga neta, y las paredes aportan la carga que compensa.
`electrokinetics_electroneutral no` (solo permitido con paredes) deja la carga
neta. En ese caso la teoría conserva el término k = 0 y la carga inducida
total es −q. Es el problema más limpio (carga aislada entre placas a tierra),
y da el menor error: el fondo discreto de contraiones agrega su propio
pequeño error.

### Error encontrado de paso (no relacionado con las paredes)

En la primera versión de estas corridas la fuerza sobre la partícula daba
**exactamente el doble** de q·E/β para z₀ ≤ 5. Causa: en
`subgrid_update_forces_electrokinetics` (`src/subgrid.c`) cada copia halo
periódica de una partícula de la celda de coloides más externa sumaba la
fuerza completa a su `fex`, y la suma de halos se la agregaba a la original
(×2 en una cara, ×4 en una arista, ×8 en una esquina). Corregido: solo la
copia dueña acumula. Todas las baterías C se repitieron con el arreglo; como
los números de las primeras versiones venían del campo (`Esub`), que era
correcto, esos resultados no cambian. Ninguna corrida de microgel ni la
batería de fuerzas entre 2 partículas tuvo partículas cargadas en celdas de
borde.

## Pendiente

- Prueba con varios procesos MPI (necesita una GPU por proceso).
- Hidrodinámica activada junto con las paredes (todas estas pruebas usan
  `hydrodynamics off`).
- Dinámica de carga de Nernst-Planck cerca de la pared (el equilibrio sí está
  verificado).
- Prueba C con electrolito (imágenes apantalladas): necesita llegar al
  estacionario de los iones en cada posición y fijar en la teoría el modo
  k = 0 con la conservación de iones.
