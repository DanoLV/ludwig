# Diario de trabajo (Ludwig)

Formato de cada entrada: fecha y hora, objetivo, cambio, prueba, resultado (observación / interpretación),
pendiente. Las observaciones son lo medido; las interpretaciones van marcadas como tales.

---

## 2026-11-05 (mañana; registrada a posteriori a las 13:10) — Refinamiento con kernels C³ (Peskin-6 y Hann cos⁴)

- **Objetivo:** ver si el orden ~1 cerca del borde del soporte del kernel (Hann-8) se debe al salto de W''.
- **Cambio:** ninguno en el código en esta etapa (las corridas Δx = ½ ya estaban en cola). Análisis con
  `capacitor/local/WALL-DIRICHLET/scripts/grafico_C_refinamiento_peskin6.py` → `graficos/C_refinamiento_peskin6.png`.
- **Prueba:** series anidadas Δx = 2, 1, ½ para Hann-8 (C¹), Peskin-6 estirado (C³) y Hann cos⁴ (C³), con q = 1 y
  q = 0.25; deriva ≤ 3e-10 en todas.
- **Observación:** orden de L (parte lineal) entre Δx = 1 y ½ en las bandas justo afuera del soporte:
  Hann-8 1.39–1.46; cos⁴ 1.88–1.93; Peskin-6 1.80–1.86. El orden del residuo total d sigue en 1.1–1.6 en esas bandas
  con los tres kernels; ahí N (parte ∝ q²) es comparable o mayor que L a Δx = ½. cos⁴ da orden 1.45 en L adentro del
  soporte (z = −11, −5), Hann 1.90.
- **Interpretación:** la pérdida de orden de L en el borde es consistente con la falta de suavidad del Hann; la de d
  se debe a que N casi no cambia con Δx. Todo depende de la separación L/N, no validada.
- **Pendiente:** validar la separación con q = ½ (no autorizado); explicar el 1.45 de cos⁴ adentro del soporte.

## 2026-11-05 (registrada a posteriori a las 13:10) — Reglas de trabajo del asistente

- **Objetivo:** fijar reglas permanentes de cómo trabaja el asistente.
- **Cambio:** creado `~/.claude/CLAUDE.md` con las 9 reglas (alcance: todos los proyectos).
- **Prueba:** comprobado que el archivo existe; no probado en una sesión nueva.
- **Resultado (observación):** archivo creado.
- **Pendiente:** decidir si el alcance queda global o solo Ludwig.

## 2026-11-05 (registrada a posteriori a las 13:10) — ¿N es el piso del solver?

- **Objetivo:** evaluar si la parte N de la separación en q podría ser el piso de error del solver.
- **Cambio:** ninguno.
- **Prueba:** residuo de PETSc en los run.log (Hann-16, Δx = ½): 1.70e-16 (q = 1) y 4.25e-17 (q = 0.25), 21
  iteraciones; deriva temporal ≤ 3e-10.
- **Observación:** el residuo escala ∝ q (cociente 4.0); N (RMS por banda) va de 5e-7 a 3e-4.
- **Interpretación:** un error ∝ q va entero a L y se cancela en N; un error absoluto ε aparecería en N como −3.2 ε y
  haría falta ε ~ 1e-4. No hay indicio de eso.
- **Pendiente:** convertir el residuo de PETSc a unidades de ψ; validación con q = ½.

## 2026-11-05 13:10 — Carga desplazada en el centro del canal (en curso)

- **Objetivo:** potencial, campo y sus errores con la carga en el centro del capacitor y desplazada fracciones de
  celda hacia la placa superior (posiciones no de simetría).
- **Cambio:** nuevos `capacitor/local/WALL-DIRICHLET/crear_desplazada.py` (crea
  `C-electrolito-hann8-LD6-desplazada/z_{17.50,17.60,17.75,17.90,18.00}`: igual que `C-electrolito-hann8-LD6/z_17.50`
  salvo la z de la carga y la cabecera del input; Hann-8, Δx = 1) y `scripts/grafico_C_desplazada.py` (residuo de ψ y
  de E_z de Ludwig contra DH lineal en el eje, perfil y RMS por bandas de ancho 2). δ = 0 se repite con el ejecutable
  actual; la corrida vieja no se tocó.
- **Prueba:** diff de input y config contra la base (solo cambian cabecera y z); solución analítica con la carga
  desplazada contra cuadratura (acuerdo 1e-12) y simetría exacta para δ = 0.
- **Resultado:** corridas en marcha (3 de 5 terminadas OK a las 13:10).
- **Pendiente:** correr el análisis cuando terminen.

## 2026-11-05 13:10 — ¿La no linealidad debería ser despreciable en el borde del kernel?

- **Objetivo:** responder si N puede ser no linealidad en el borde del soporte, donde ψ es chico.
- **Cambio:** ninguno.
- **Prueba:** ψ analítico (DH, Hann-8, z0 = −8) por bandas, comparado con N (Hann-8, Δx = ½) de la separación en q.
- **Observación:** ψ(z0) = 0.317. Desde z = −15 hasta −3, N/ψ queda entre 4.4e-4 y 5.9e-4, mientras ψ²/6 (la
  no linealidad local relativa) va de 1.4 % en la carga a ~1e-3 en el borde. En el borde (z = −3) N = 4.1e-5 y ψ³/6 =
  1.4e-4.
- **Interpretación:** N sigue la forma de ψ, no la de ψ³ local: es consistente con una corrección no local (la
  respuesta apantallada a la fuente κ²ψ³/6 concentrada cerca de la carga, que se propaga como ψ). Es chica frente a ψ
  (~5e-4) pero comparable al error de discretización L en Δx = ½.
- **Pendiente:** contrastar N con una solución no lineal independiente; validación con q = ½.

## 2026-11-05 13:21 — Carga desplazada en el centro del canal: resultados

- **Objetivo:** potencial, campo y sus errores con la carga en z0 = δ (centro del canal + δ hacia la placa superior),
  δ = 0, 0.10, 0.25, 0.40, 0.50; Hann-8, Δx = 1.
- **Cambio:** ninguno nuevo (corridas y script de la entrada de las 13:10). Figura `graficos/C_desplazada.png`.
- **Prueba:** las 5 corridas terminaron OK. Con δ = 0 el RMS por bandas sale exactamente simétrico respecto del plano
  medio; con δ = 0.5 (carga sobre un nodo) d_E en ese nodo = 2.5e-7 (≈ 0 por simetría). NO medida: la deriva al
  estacionario (solo hay salidas en los pasos 1 y 5000, igual que la corrida base).
- **Observación:**
  - Dentro del soporte: max|d_ψ| = 1.15e-2 a 1.20e-2 (≈ 3.5 % del pico de ψ, 0.33); max|d_E| = 4.7e-3 a 5.1e-3
    (≈ 7 % del pico de E_z, 0.07). Varía menos de 10 % con δ; la forma del residuo se traslada con la carga.
  - Fuera del soporte: d_ψ ~ 1e-5–1e-4 y d_E ~ 1e-5. En las bandas pegadas al borde (z = −5, 5) el RMS de d_ψ cambia
    hasta un factor ~10 con δ (banda 5: 3.2e-5 con δ = 0.25, 3.4e-4 con δ = 0.5).
- **Interpretación:** el error grande está dentro del soporte y casi no depende de la posición fraccional. La
  sensibilidad a δ fuera del soporte aparece en el borde, donde cambia la posición del borde respecto de los nodos
  (con δ = 0.5 los bordes caen sobre nodos). Las bandas están fijas en el laboratorio, así que parte de esa variación
  es solo que la banda contiene otra distancia a la carga. El residuo incluye la parte no lineal (no separada).
- **Pendiente:** medir la deriva; separar L/N si hace falta (requiere corridas con q = 0.25); campo en la posición de
  la partícula (sugerido, no autorizado).

## 2026-11-05 13:54 — Verificación del estacionario de la carga desplazada (en curso)

- **Objetivo:** correr 5000 pasos más de C-electrolito-hann8-LD6-desplazada, verificar el estacionario y repetir el
  análisis si lo alcanzó; agregar el error RMS relativo por banda.
- **Cambio:**
  - Primero lancé por error las corridas desde cero hasta 10000 pasos (`crear_desplazada.py N=10000`,
    `C-electrolito-hann8-LD6-desplazada-N10000`), sin consultar; el usuario pidió continuar desde el paso 5000. Las
    detuve en el paso ~3000; el directorio queda sin usar.
  - Nuevo `crear_continuacion.py`: copia cada z_* (con las salidas del paso 5000) a
    `C-electrolito-hann8-LD6-desplazada-cont` y reanuda con N_start 5000, N_cycles 5000, freq_config 1000;
    `correr_cont.sh` no borra salidas previas. Las corridas originales no se tocaron.
  - `scripts/grafico_C_desplazada.py`: argumento `cont`; deriva D_psi y D_E en toda la red entre los dos últimos
    volcados; fila nueva e_b = RMS_b(d)/RMS_b(analítico).
- **Prueba:** reanudación de 1 paso en una copia en el scratchpad: Ludwig lee dist, rho/vel, psi/qsi y config.cds del
  paso 5000; max|psi(5001) − psi(5000)| = 2.1e-15 (max|psi| = 0.34); los archivos del paso 5000 quedan idénticos (cmp).
- **Criterio de estacionario (elegido por el usuario):** D < 1 % del RMS por banda más chico del residuo (para ψ y
  para E por separado). Heurístico.
- **Pendiente:** esperar las corridas y aplicar el criterio.

## 2026-11-05 14:05 — Fuerza espuria por el desplazamiento (provisorio, paso 5000)

- **Objetivo:** comparar la fuerza de Ludwig sobre la carga desplazada con la fuerza analítica distribuida por el
  kernel y con la fuerza analítica en la posición de la partícula.
- **Cambio:** nuevo `scripts/grafico_C_fuerza_desplazada.py` → `graficos/C_fuerza_desplazada.png`.
  F_L = Esub_z de `particle_Esub.csv` (Σ W(r_i − r0) E(r_i), subgrid_update_Esub); F_kernel = ∫ W(r − r0) E_DH,z d³r
  (DH lineal, x-y en Fourier, z por partes: ∫ W'(z − z0) Φ̄(z) dz); F_point = E_DH,z(r0).
- **Prueba:** F_L recalculado desde el volcado efield coincide con Esub a 7 cifras (unidades: Esub ya es eE/kT);
  δ = 0 da F = 0 (1e-18) en Ludwig y en ambas analíticas; F_kernel no cambia al duplicar los puntos de Gauss ni al
  bajar los modos laterales; forma integrada por partes vs derivada directa: 2e-14.
- **Observación (paso 5000, estacionario aún no verificado):**
  δ = 0.10 / 0.25 / 0.40 / 0.50: F_L = 1.273e-5 / 3.167e-5 / 5.017e-5 / 6.214e-5;
  F_L − F_kernel = 3.8e-7 / 7.8e-7 / 7.2e-7 / 2.9e-7 (3.1 / 2.5 / 1.4 / 0.5 % de F_kernel);
  F_L − F_point = 7.4e-7 / 1.7e-6 / 2.1e-6 / 2.1e-6.
- **Interpretación:** la referencia adecuada es F_kernel (misma definición que Ludwig); F_point difiere de F_kernel
  ~3 % por el ensanchamiento, no por error. F_L − F_kernel mezcla autofuerza de red, error de discretización de la
  fuerza de imagen y la parte no lineal; con estos datos no se separan. Con δ = 0.5 (nodo) sigue siendo distinta de 0.
- **Pendiente:** repetir con la continuación al paso 10000 si es estacionaria.

## 2026-11-05 14:20 — Carga desplazada: estacionario verificado; análisis repetido al paso 10000

- **Objetivo:** aplicar el criterio de estacionario a la continuación y repetir los análisis de campo y fuerza.
- **Cambio:** ninguno en el código. Figuras `graficos/C_desplazada_cont.png` y `graficos/C_fuerza_desplazada_cont.png`.
- **Prueba:** las 5 continuaciones terminaron OK (reanudadas desde el paso 5000).
- **Observación:**
  - Estacionario: D_psi = 2e-16 a 2e-15 y D_E = 7e-17 a 4e-16 (pasos 9000 → 10000), contra umbrales 1 % del RMS
    mínimo = 1.2e-7 a 1.9e-7 (ψ) y 2.4e-8 a 1.1e-7 (E): cumple en las 5 posiciones. Esub cambia ≤ 2.6e-14 entre los
    pasos 5000 y 10000.
  - Todos los números de campo y fuerza coinciden con los del paso 5000 en las cifras mostradas.
  - Error RMS relativo por banda: dentro del soporte ψ 3.0–3.5 %, E_z 4.8–7.9 %; fuera, ψ 0.02–0.4 % y E_z
    0.03–0.4 %; en las bandas junto a las placas ψ sube a 1.3 % (ψ_DH → 0 en las placas a tierra, el residuo no).
  - Fuerza: F_L − F_kernel = 3.8e-7 / 7.8e-7 / 7.2e-7 / 2.9e-7 para δ = 0.10 / 0.25 / 0.40 / 0.50 (3.1 / 2.5 / 1.4 /
    0.5 %), igual que al paso 5000.
- **Interpretación:** las corridas de 5000 pasos ya estaban en el estacionario; los resultados provisorios quedan
  confirmados.
- **Pendiente:** separar las contribuciones a F_L − F_kernel (autofuerza de red, discretización, no linealidad) si
  hace falta; borrar o no `C-electrolito-hann8-LD6-desplazada-N10000`.

## 2026-11-05 14:55 — Carga desplazada: celda completa (δ = 0.6–0.9) y errores nodo a nodo

- **Objetivo:** completar el desplazamiento de una celda cada 0.1 (δ = 0.60, 0.70, 0.80, 0.90) y ver cómo cambia la
  fuerza espuria; agregar al gráfico de campo el error absoluto y relativo en cada nodo.
- **Cambio:**
  - Borrado `C-electrolito-hann8-LD6-desplazada-N10000` (pedido por el usuario). La opción N= se mantiene.
  - Lancé por error δ = 0.60, 0.75, 0.90 (espaciado 0.15) con continuación automática a 10000; el usuario corrigió:
    paso 0.1 y no continuar si a 5000 ya está en equilibrio. Detenidas; 0.75 quedó apartada como
    `descartada_z_18.25` (sin borrar).
  - `crear_desplazada.py`: DELTAS con 0.60, 0.70, 0.80, 0.90; las de δ > 0.5 con freq_config 1000; no recrea
    corridas existentes. `crear_continuacion.py`: saltea continuaciones existentes.
  - `grafico_C_desplazada.py`: filas nuevas |d| y |d/analítico| en cada nodo; modo `celda` (δ ≤ 0.5 de -cont,
    δ > 0.5 de la base). `grafico_C_fuerza_desplazada.py`: mismo modo. Figuras `C_desplazada_celda.png`,
    `C_fuerza_desplazada_celda.png`.
- **Prueba:** 4 corridas OK. Estacionario (criterio del usuario) entre los pasos 4000 y 5000: D_psi = 1e-11, D_E =
  3e-12 a 5e-12, contra umbrales 1.3e-7 a 2.0e-7 y 3.7e-8 a 5.6e-8: cumple. Esub desde el volcado efield = Esub de
  Ludwig a 7 cifras.
- **Observación:** F_L − F_kernel (×1e-7) para δ = 0, 0.1, 0.25, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9:
  0, 3.8, 7.8, 7.2, 2.9, −1.4, −2.5, −1.0, +2.0. Relativo a F_kernel: ≤ 3.1 % (δ = 0.1), 0.18 % en δ = 0.9.
  F_L − F_point crece de 0 a 3.4e-6 en forma casi monótona (≈ 3 % de F). Dentro del soporte, max|d_ψ| 1.15–1.20e-2
  y max|d_E| 4.7–5.1e-3 para todos los δ.
- **Interpretación:** en una celda, F_L − F_kernel oscila con la posición respecto de los nodos (cambia de signo
  entre 0.5 y 0.9) en lugar de crecer al acercarse a la placa; no se anula en δ = 0.5, así que hay además una parte
  no periódica. Una celda (0.9 de 16 de distancia a la placa) no alcanza para saber cómo cambia cerca de la placa.
- **Pendiente:** posiciones más cerca de la placa (sugerido, no autorizado).

## 2026-11-05 15:57 — Carga desplazada: acercamiento a la placa (δ = 1 a 12)

- **Objetivo:** ver cómo cambia la fuerza espuria al acercar la carga a la placa superior.
- **Cambio:** `crear_desplazada.py`: 10 posiciones δ = 1.0, 2.2, 3.4, 4.7, 5.9, 7.1, 8.3, 9.6, 10.8, 12.0 (equiespaciadas,
  redondeadas a 0.1; δ = 12 deja el borde del soporte justo en la placa z = +16, a mitad del enlace). Modo `placa` en
  los dos scripts de análisis: `graficos/C_desplazada_placa.png`, `graficos/C_fuerza_desplazada_placa.png`.
- **Prueba:** 10 corridas OK; estacionario (criterio del usuario, pasos 4000 → 5000): D_psi = 2e-11 a 2e-10, D_E =
  5e-12 a 5e-11, umbrales ≥ 1.1e-7 (ψ) y ≥ 1.5e-8 (E): cumple. Esub desde efield = Esub a 7 cifras.
- **Observación:** distancia a la placa 15.0 → 4.0: F_kernel crece de 1.24e-4 a 9.19e-3. F_L − F_kernel = 6.0e-7,
  1.9e-6, 2.7e-6, 1.9e-6, 3.4e-6, 6.5e-6, 1.0e-5, 5.1e-6, 8.2e-6, 5.8e-5; relativo 0.48, 0.69, 0.58, 0.28, 0.35, 0.48,
  0.54, 0.17, 0.17, 0.63 %. F_point/F_kernel: 2.9 % por debajo lejos, 2.3 % en δ = 12. En las dos posiciones de punto
  medio (δ = 1 y 12) el relativo es 0.48 y 0.63 %. Campo dentro del soporte: max|d_ψ| 1.15–1.20e-2, max|d_E|
  4.7–5.1e-3, como en el centro.
- **Interpretación:** el error de la fuerza crece en valor absoluto con la fuerza, pero en relativo queda entre 0.2 y
  0.7 % sin tendencia clara al acercarse; la dispersión podría venir de la posición dentro de la celda (cambia de una
  posición a otra), que aquí no se controla por separado. La parte que oscilaba en la celda central (±5e-7) es
  chica frente a estas diferencias.
- **Pendiente:** si hace falta separar distancia y posición en la celda, barrer la celda a una distancia fija cerca de
  la placa (no autorizado).

## 2026-11-05 16:07 — Gráfico de fuerza con errores relativos

- **Objetivo:** agregar al gráfico de la fuerza el error relativo al valor analítico en la posición y al analítico
  distribuido por el kernel.
- **Cambio:** `scripts/grafico_C_fuerza_desplazada.py`: dos paneles nuevos, (F_L − F_point)/F_point y
  (F_L − F_kernel)/F_kernel en % (δ = 0 se omite: F = 0). Regeneradas `C_fuerza_desplazada_placa.png` y
  `C_fuerza_desplazada_celda.png` (mismos datos).
- **Prueba:** los valores de los paneles coinciden con la tabla calculada antes (p. ej. δ = 12: 2.95 % y 0.63 %).
- **Observación:** celda central: relativo a F_kernel baja de 3.1 % (δ = 0.1) a −0.3 % (δ = 0.7); relativo a F_point
  de 6.1 % a 2.6 %. Placa: relativo a F_kernel 0.17–0.69 %, relativo a F_point 2.9–3.7 %.
- **Interpretación:** en la celda central el relativo es grande porque la fuerza es chica (F → 0 en δ = 0) y la
  parte que depende de la posición en la celda (±5e-7) pesa más.
- **Pendiente:** ninguno nuevo.

## 2026-11-05 17:13 — Acercamiento a la placa: puntos medios (19 posiciones)

- **Objetivo:** agregar los puntos medios entre las posiciones del acercamiento a la placa.
- **Cambio:** `crear_desplazada.py`: δ = 1.6, 2.8, 4.0, 5.3, 6.5, 7.7, 9.0, 10.2, 11.4 (4.05 → 4.0 y 8.95 → 9.0,
  elección del usuario). Modo `placa` con las 19 posiciones; regenerados `C_desplazada_placa.png` y
  `C_fuerza_desplazada_placa.png`.
- **Prueba:** 9 corridas OK; estacionario 4000 → 5000 cumple en las 19 (ninguna "NOT steady"; D_psi ≤ 1.7e-10).
- **Observación:** (F_L − F_kernel)/F_kernel alterna entre ~0.2 y ~0.7 % de un punto al siguiente. Agrupado por
  posición en la celda (δ mod 1; 0 = punto medio, 0.5 = nodo): 0.0 → 0.48, 0.46, 0.38, 0.63 %; 0.2 → 0.69, 0.59;
  0.3 → 0.57, 0.54; 0.4 → 0.58, 0.64; 0.6 → 0.21, 0.17; 0.7 → 0.28, 0.21; 0.8 → 0.28, 0.17 % (δ creciente).
  (F_L − F_point)/F_point entre 2.9 y 3.7 %.
- **Interpretación:** la dispersión del error relativo de la fuerza viene sobre todo de la posición dentro de la celda
  (más alto entre 0.2 y 0.4, más bajo entre 0.6 y 0.8, el mismo patrón de signo que en la celda central); a posición
  en la celda fija cambia poco con la distancia a la placa (≤ 0.25 puntos porcentuales).
- **Pendiente:** ninguno nuevo.

## 2026-11-05 17:24 — Carga desplazada con Δx/2 (en curso)

- **Objetivo:** 4 posiciones con Δx = 1/2: extremos de la serie de acercamiento y dos intermedias lo más
  equiespaciadas posible entre las ya corridas: δ = 1.0, 4.7, 8.3, 12.0 (elegido por el usuario).
- **Cambio:** nuevo `crear_desplazada_ref.py` → `C-electrolito-hann16-LD12-ref2-desplazada/z_{35.50,42.90,50.10,57.50}`
  (Hann-16, 64×64×66, ε = 5000, ρ_el/8, N = 20000; misma receta que crear_refinado.py; psi cada 4000 pasos = mismo
  tiempo físico que 1000 con Δx = 1). Scripts de análisis: modo `ref2` (R = 2: lectura de la red fina, Esub y efield
  × R para pasar a campo físico).
- **Prueba:** input igual al de `C-electrolito-hann16-LD12-ref2/z_17.50` salvo comentarios, posición y freq_config;
  .petscrc idéntico. Recalculo de Σ W·E desde efield con Hann-16 en esa corrida existente: −1.74007191e-3 vs Esub·2 =
  −1.74007191e-3.
- **Criterio de estacionario elegido por mí (no consultado):** deriva entre los pasos 16000 y 20000 (equivalente en
  tiempo físico a 4000 → 5000 con Δx = 1), mismo umbral del 1 %.
- **Resultado:** corridas en marcha (2 a la vez, ~4–5 h).
- **Pendiente:** analizar cuando terminen.

## 2026-11-05 21:28 — Carga desplazada con Δx/2: δ = 1.0 y 4.7 (parcial)

- **Objetivo:** comparar la fuerza y el error de ψ/E con Δx/2 contra Δx = 1 en las mismas δ.
- **Cambio:** ninguno en código; corridos `grafico_C_desplazada.py ref2` y `grafico_C_fuerza_desplazada.py ref2`
  (las figuras `_ref2` tienen por ahora solo estas dos δ).
- **Prueba:** estacionario 16000 → 20000 cumple (D_psi ≤ 7e-11, D_E ≤ 2.4e-11, umbrales ~4e-8 y ~4e-9). F_L del csv
  coincide con Σ W·E recalculada desde efield. Duración real: ~4 h por corrida (no 4–5 h por el lote).
- **Observación:** (F_L − F_kernel)/F_kernel: δ = 1.0 0.48 % → 0.072 %; δ = 4.7 0.28 % → 0.073 % (cociente de
  diferencias 6.7 y 3.8). (F_L − F_point)/F_point: 3.43 → 3.01 %, 3.22 → 3.00 %. max|d_psi| 2.7e-3, max|d_E| 1.3e-3
  (Δx = 1: ~1.2e-2 y ~5e-3).
- **Interpretación:** F_L converge a F_kernel; el ~3 % respecto de F_point es la diferencia kernel vs puntual, no error
  de red. Los cocientes no dan un orden limpio: en δ = 4.7 la posición en la celda cambia (0.2 con Δx = 1, 0.9 con Δx/2).
- **Pendiente:** δ = 8.3 y 12.0 (~01:20); regenerar figuras _ref2 con las 4; el usuario no confirmó el criterio de
  estacionario.

## 2026-11-06 01:26 — Carga desplazada con Δx/2: las 4 posiciones

- **Objetivo:** completar la comparación Δx/2 vs Δx = 1 en δ = 1.0, 4.7, 8.3, 12.0.
- **Cambio:** ninguno en código; regeneradas `graficos/C_desplazada_ref2.png` y `graficos/C_fuerza_desplazada_ref2.png`.
- **Prueba:** las 4 corridas terminaron [OK] (01:23). Estacionario 16000 → 20000 cumple en las 4 (D_psi ≤ 1.0e-10,
  D_E ≤ 3.3e-11). F_L del csv = Σ W·E desde efield en las 4.
- **Observación:** (F_L − F_kernel)/F_kernel, Δx = 1 → Δx/2: δ = 1.0 0.48 → 0.072 %; 4.7 0.28 → 0.073 %;
  8.3 0.54 → 0.013 %; 12.0 0.63 → 0.043 %. (F_L − F_point)/F_point: 3.43 → 3.01; 3.22 → 3.00; 3.44 → 2.90;
  2.95 → 2.35 %. max|d_psi| 2.7–2.8e-3 y max|d_E| 1.3e-3 en las 4 (no cambia con la distancia a la placa).
- **Interpretación:** F_L converge a F_kernel en las 4; los cocientes (6.7, 3.8, 41, 14) no dan un orden único:
  mezclan refinamiento con posición en la celda (δ = 1 y 12 mantienen el punto medio y aun así dan 6.7 y 14).
- **Pendiente:** el usuario no confirmó el criterio de estacionario; figuras _ref2 no inspeccionadas visualmente.

## 2026-11-06 10:07 — Convergencia Δx = 1 → 1/2 de la carga desplazada

- **Objetivo:** gráfico de convergencia de ψ, E y F entre Hann-8 (Δx = 1) y Hann-16 (Δx/2) en δ = 1.0, 4.7, 8.3, 12.0,
  con F_kernel (límite continuo) como fuerza objetivo, y órdenes p.
- **Cambio:** nuevo `scripts/grafico_C_convergencia_desplazada.py` → `graficos/C_convergencia_desplazada_{campo,fuerza}.png`.
  Títulos de `grafico_C_desplazada.py` y `grafico_C_fuerza_desplazada.py` dependen del modo (KLAB); dpi 170 → 300 en
  ambos. ψ y E por bandas de ancho 2 (elección del usuario: los nodos z de las dos redes no coinciden).
- **Prueba:** F_L y F_kernel reproducen los valores ya impresos por grafico_C_fuerza_desplazada.py; figuras inspeccionadas.
- **Observación:** p de F_L − F_kernel = 2.73, 1.94, 5.37, 3.85. Bandas dentro del soporte: p ≈ 1.7–2.7 en ψ y E salvo
  puntos aislados; fuera, p ≈ 1.6–2.4 con valores sueltos de −2 a 4.7 en bandas donde el RMS pasa por un mínimo.
- **Interpretación:** consistente con orden ~2 para ψ y E; los p extremos son bandas donde el error cambia de signo
  (cociente sin significado de orden). El p de la fuerza no es un orden limpio (posición en la celda, una sola pareja).
- **Pendiente:** parte no lineal PB no separada; criterio de estacionario Δx/2 sin confirmar por el usuario.

## 2026-11-06 10:28 — F_point fuera del análisis de fuerza

- **Objetivo:** quitar F_point = E_DH,z(r0) de los análisis (decisión del usuario: no es la fuerza teórica sobre la
  carga distribuida).
- **Cambio:** `scripts/grafico_C_fuerza_desplazada.py`: sin F_point en la salida ni en la figura (3 paneles: F_L y
  F_kernel, F_L − F_kernel, relativo a F_kernel). Copia previa en el scratchpad de la sesión. grafico_C_convergencia_desplazada.py no lo usaba.
- **Prueba:** modo ref2: F_L, F_kernel y F_L − F_kernel idénticos a los anteriores; figura inspeccionada.
- **Pendiente:** regeneración de las figuras de fuerza (placa, celda, cont, base) y de campo ref2/placa a 300 dpi en curso.

## 2026-11-06 14:09 — Derivada analítica de φ_DH y Δx = 2 en la convergencia de la fuerza

- **Objetivo:** (1) calcular F_kernel y E_z,DH con la derivada analítica de φ_DH; (2) agregar Δx = 2 (Hann-4) al gráfico
  de convergencia de la fuerza.
- **Cambio:** `grafico_C_hoja_hann.py`: nuevas dpsit, dU, dphi (forma cerrada). `grafico_C_desplazada.py`: ez_an = −dphi
  (antes diferencia finita, queda como ez_fd). `grafico_C_fuerza_desplazada.py`: F_kernel = −∫W dΦ̄/dz directo (la forma
  por partes queda como F_kernel_partes). `crear_desplazada_ref.py` acepta r= (por defecto 2); `r=0.5` creó
  C-electrolito-hann4-LD3-ref0.5-desplazada (16×16×18, N = 2000, psi cada 250). `grafico_C_convergencia_desplazada.py`:
  Δx = 2 en la figura de fuerza (p por pareja de redes) y chequeo de estacionario de esas corridas.
- **Prueba:** dpsit vs cuadratura ≤ 4e-13 relativo, vs FD ≤ 9e-11; E_z vs FD ≤ 5.6e-11; F_kernel directo vs por partes
  ≤ 3.5e-17. Input Δx = 2 igual al de C-electrolito-hann4-LD3-ref0.5 salvo comentarios y freq_config; mismo ejecutable
  (md5). 4 corridas OK (3 min 46 s); estacionario 1750 → 2000 cumple en las 4.
- **Observación:** (F_L − F_kernel)/F_kernel con Δx = 2 / 1 / ½: δ = 1.0 1.93 / 0.48 / 0.072 %; 4.7 3.11 / 0.28 / 0.073 %;
  8.3 3.29 / 0.54 / 0.013 %; 12.0 5.81 / 0.63 / 0.043 %. p (2→1) = 2.01, 3.47, 2.61, 3.21; p (1→½) = 2.73, 1.94, 5.37, 3.85.
- **Interpretación:** convergencia de al menos orden ~2 en todas las δ; los p no son constantes (posición de la carga en la
  celda distinta en cada red: con Δx = 2, δ = 1 cae en un nodo y δ = 12 en un punto medio).
- **Pendiente:** modo cont de los scripts falla para δ > 0.5 (decidir); regeneración de figuras con dphi en curso.

## 2026-11-06 14:47 — Δx = 2 en todas las δ: por qué fluctúa el orden de la fuerza

- **Objetivo:** entender por qué p de F_L − F_kernel fluctúa entre posiciones.
- **Cambio:** `crear_desplazada_ref.py r=0.5` crea las 28 δ de Hann-8 (24 nuevas; comentarios del input con z a 3
  decimales y δ a 2). `grafico_C_convergencia_desplazada.py todos` → `graficos/C_convergencia_desplazada_fuerza_todos.png`
  (error relativo vs δ para Δx = 2, 1, ½ y p(2→1), p(1→½); nodos de cada red marcados). Δx = 1 con δ ≤ 0.5 usa -cont.
- **Prueba:** 24 corridas OK (bateria2.log); estacionario 1750 → 2000 cumple en las 28 Δx = 2. Figura inspeccionada.
- **Observación:** p(2→1) va de −0.3 a 5.9. Los p bajos (δ = 1.6: 0.9; 3.4: −0.3; 5.3: 0.4; 7.1: 1.2; 11.4: −0.1) están
  0.1–0.6 por encima (hacia la placa) de un nodo de Δx = 2; ahí el error de Δx = 2 es mínimo (en 1.6 cambia de signo:
  −0.39 %). Los máximos del error de Δx = 2 (2.2, 2.8, 4.7, 6.5, 8.3, 10.2, 10.8: 3–4 %) están por debajo de un nodo.
  Δx = 1 mostraba el mismo patrón respecto de sus nodos (alto en δ mod 1 = 0.2–0.4, bajo en 0.6–0.8). Cerca del centro
  (δ < 1) el error relativo crece (Δx = 2: 4–14 %, Δx = 1: 0.1–3 %) porque F_kernel → 0.
- **Interpretación (hipótesis):** el error tiene una parte que depende de la posición de la carga en la celda (período Δx,
  de signo opuesto a cada lado de un nodo) además de una parte suave en δ. A δ fijo cada red está en otra posición de su
  celda, así que p mezcla las dos; los p extremos aparecen donde el error de una red pasa cerca de cero.
- **Pendiente:** probar la hipótesis comparando a posición en la celda fija (propuesto, no hecho); modo cont.

## 2026-11-06 17:11 — Fuerza desplazada: δ = 1…4 cada 0.1 en Hann-4 y Hann-8
- **Objetivo:** entender por qué el orden observado p(Δx=2→1) de la fuerza fluctúa con δ.
- **Cambio:** 25 posiciones nuevas (δ = 1.1…3.9, paso 0.1) en C-electrolito-hann4-LD3-ref0.5-desplazada y C-electrolito-hann8-LD6-desplazada (crear_desplazada.py, crear_desplazada_ref.py). ALL en grafico_C_convergencia_desplazada.py pasa a 53 posiciones.
- **Prueba:** 50/50 corridas OK (bateria3.log, bateria6.log). Estacionariedad (variación entre los dos últimos volcados < 1 % del RMS mínimo por banda): todas OK; Δx=2 1750→2000, Δx=1 4000→5000 (script en scratchpad, steady_dx1.py).
- **Observación:** el error relativo oscila con período igual a la celda de cada red: 2 en Δx=2 (máx +3.99 % en δ=2.5, mín ~+0.5 % en δ=3.3–3.4, cambio de signo entre 1.6 y 1.7) y 1 en Δx=1 (máx ~0.7 % en δ≈1.3, 2.3, 3.3; mín ~0.25 % en δ≈1.7, 2.7, 3.7). Los p extremos coinciden con los mínimos de una red y los máximos de la otra (p=−2.3 en δ=1.1; p=3.9 en δ=2.7).
- **Interpretación:** con δ fijo, cada red está en otra fracción de celda, así que p(δ) mezcla el orden de convergencia con la dependencia del error respecto de la posición dentro de la celda. En fracción de celda, máximos y mínimos caen en la misma posición relativa en ambas redes.
- **Pendiente:** decidir cómo definir un orden que compare la misma fracción de celda o el promedio por celda; modo cont; etiqueta "NOT steady" del modo base.

## 2026-11-06 23:20 — Fuerza desplazada: grilla completa δ = 0…12 cada 0.1 en Hann-4 y Hann-8
- **Objetivo:** ver la dependencia del error de F_L con la posición en la celda en todo el canal.
- **Cambio:** 69 δ nuevas en cada red (crear_desplazada_ref.py r=0.5, crear_desplazada.py; δ = 0.2, 0.3 de Hann-8 con freq_config 1000). Figura "todos": 4 paneles (F_L sobre F_kernel; (F_L−F_kernel)/F_kernel con signo; |rel| log; p). ALL = 0…12 cada 0.1 más 0.25 (122 posiciones).
- **Prueba:** 69/69 OK en bateria4.log (Hann-4) y bateria7.log (Hann-8). Estacionariedad (< 1 % del RMS mínimo por banda): 69/69 en cada red (Δx=2 1750→2000, Δx=1 4000→5000).
- **Observación:** Δx=2 oscila con período 2: máximos +3.2…+4.8 % en δ = 2.5, 4.5, 6.5, 8.5, 10.5; mínimos +0.35…+0.7 % en δ ≈ 3.3–5.4–7.4–9.3 (y −1.2 % en 1.3). Δx=1 oscila con período 1: máximos +0.54…+0.88 % en δ = k.3; mínimos +0.10…+0.28 % en δ = k.7. p(2→1) oscila entre −0.9 y 5.4 fuera del centro, con período 2.
- **Interpretación:** en fracción de celda, máximos y mínimos caen en posiciones parecidas en ambas redes (mínimo ~0.2 celda después de un nodo; máximo ~0.2–0.25 celda antes del siguiente). p a δ fijo mezcla el orden con esa dependencia.
- **Pendiente:** definir un orden comparable (misma fracción de celda o promedio por celda) si el usuario lo decide; modo cont; etiqueta "NOT steady" del modo base.

## 2026-11-07 21:35 — Fuerza desplazada: análisis parcial Hann-16 (δ = 1.0…2.0 cada 0.1)
- **Objetivo:** ver la dependencia del error con la posición en la celda también en Δx = 1/2.
- **Cambio:** 20 δ nuevas (1.1…3.0) para r = 2 en crear_desplazada_ref.py, lanzadas 2 en paralelo (bateria2.log). grafico_C_convergencia_desplazada.py modo todos: saltea corridas Δx = 1/2 sin terminar (antes abortaba).
- **Prueba:** 10/10 OK hasta δ = 2.0. Estacionariedad Δx = 1/2 (16000→20000, < 1 % del RMS mínimo por banda): OK en las 14 terminadas (scratchpad steady_dxhalf.py).
- **Observación:** Δx = 1/2 positivo en 1.0–2.0, oscila con período ~0.5: máx 0.12 % (δ 1.1), 0.11 % (1.6); mín 0.030 % (1.4), 0.039 % (1.9). Nodos de Δx = 1/2 en δ = 1.25, 1.75. p(1→½) entre 0.53 (δ 1.7) y 4.8 (δ 1.3).
- **Interpretación:** mínimos ~0.1–0.3 celda después del nodo y máximos ~0.1–0.3 celda antes, como en Δx = 2 y 1 (resolución de 0.1 = 0.2 celda: localización gruesa). p a δ fijo vuelve a mezclar orden y posición en la celda.
- **Pendiente:** δ = 2.1–3.0 (~20 h), luego análisis completo de la celda 1–3.

## 2026-11-08 18:02 — Fuerza desplazada: Hann-16 completo en la celda δ = 1…3 cada 0.1
- **Objetivo:** cubrir con Δx = 1/2 la celda δ = 1–3 (decisión del usuario, frente a la grilla completa ~10 días).
- **Cambio:** ninguno nuevo; terminó el lote de 20 corridas (bateria2.log: 20 OK, 0 FALLO; fin 2026-11-08 17:58). Figura todos regenerada.
- **Prueba:** estacionariedad Δx = 1/2 (16000→20000, < 1 % del RMS mínimo por banda): 24/24 OK. Modo todos: ninguna "NOT steady", ninguna sin terminar.
- **Observación:** Δx = 1/2 positivo en toda la celda; período ~0.5; máximos 0.12, 0.11, 0.10, 0.09 % en δ = 1.1, 1.6, 2.1, 2.6; mínimos 0.030, 0.039, 0.044, 0.046 % en δ = 1.4, 1.9, 2.4, 2.9. p(1→½) entre 0.53 (δ 1.7) y 4.77 (δ 1.3); entre 2.5 y 3.1 en 11 de las 21 posiciones.
- **Interpretación:** los extremos de Δx = 1/2 caen a ~0.3 celda de sus nodos (δ = k.25, k.75), parecido a Δx = 2 y 1 (~0.2 celda); con muestreo de 0.2 celda la localización es gruesa. La fluctuación de p(1→½) tiene el mismo origen que la de p(2→1).
- **Pendiente:** si el usuario lo decide, un orden que compare la misma fracción de celda o el promedio por celda; modo cont; etiqueta "NOT steady" del modo base.

## 2026-11-08 19:36 — ψ y E_z puntuales en el centro de la carga desplazada, 3 redes
- **Objetivo:** figura como la de fuerzas para el potencial y el campo en r0 (decisión del usuario: valor puntual interpolado; paneles valor / error relativo con signo / p).
- **Cambio:** script nuevo scripts/grafico_C_carga_desplazada_puntual.py -> graficos/C_carga_desplazada_puntual.png. Ludwig: ψ y E_z del último paso en la línea x = y = carga, interpolación lineal en z en z0. Referencia: φ_DH(r0) y E_DH,z(r0) (DH lineal, fuente Hann continua, derivada cerrada).
- **Prueba:** corre sin errores; dos corridas con salida idéntica. No verifiqué la interpolación contra otra independiente.
- **Observación:** ψ: error positivo; Δx=2 de +6.5 % (δ par, mitad de celda) a +16.5 % (δ impar, nodo); Δx=1 de +1.65 % (δ entero) a +3.67 % (δ = k+0.5, nodo); Δx=1/2 +0.37…+0.70 %. p_ψ(2→1) 1.3–3.3, p_ψ(1→½) 1.55–3.3. E_z: en nodos de cada red el error es chico (Δx=2 en δ=1: +2.0 %; Δx=1 en δ=0.5, 1.5: +0.4 %); entre nodos llega a +2700 % / −950 % (Δx=2) y +1400 % (Δx=1) cerca del centro. p_E(2→1) con picos de 7–13 en δ = k+0.5.
- **Interpretación:** el error de E_z puntual entre nodos parece dominado por la interpolación lineal del campo propio de la carga (grande y de pendiente fuerte cerca de r0), no por el error de la solución en los nodos. Hipótesis, no verificada.
- **Pendiente:** decidir si interesa separar error de interpolación vs error en nodos.

## 2026-11-09 12:23 — Error de ψ y E en los nodos: capas de distancia y normas L2/L∞
- **Objetivo:** evaluar el error de ψ y E como campos (por separado de la fuerza), sin interpolar.
- **Decisiones del usuario:** capas de ancho 1 (parámetro w=), error absoluto dividido por el RMS de la referencia en la capa, δ = 1.0, 4.7, 8.3, 12; además normas L2 (√Σ Δx³ e²) e L∞.
- **Cambio:** script nuevo scripts/grafico_C_error_capas_desplazada.py -> graficos/C_error_capas_desplazada.png y _orden.png. φ_DH y E_DH en todos los nodos con FFT 2D de los modos plegados mod N (exacto en los nodos).
- **Prueba:** FFT vs HH.phi/HH.dphi en el eje ≤ 1e-15; fuera del eje φ ≤ 1e-15, E_x vs diferencia finita ≤ 5e-11; simetría x↔y 6e-17. No verifiqué otra cosa del cálculo de error.
- **Observación:** normas globales: ψ L2 rel 4.5–7.0 % / 1.0–1.6 % / 0.23–0.39 % (Δx = 2/1/½), E L2 rel ~15 % / ~4 % / ~1 %. p entre 1.84 y 2.23 para ψ y E, L2 y L∞, en los 4 δ y ambos pares de redes. Por capas: dentro del kernel (d < 4) p(1→½) ≈ 2 (ψ 2.1–2.3, E 1.5–2.0); lejos, ψ → 2 y E 1.7–2.0; caídas de p en capas donde el error relativo tiene un mínimo pronunciado (ψ en d ≈ 8–10).
- **Interpretación:** como campos, ψ y E convergen con orden ~2 de forma estable; la fluctuación fuerte del orden de la fuerza no aparece en estas normas. Los mínimos de error por capa parecen cruces de signo del error (no verificado).
- **Pendiente:** confirmar los cruces de signo; separar el error de la fuerza en solución vs cuadratura si el usuario lo pide.

## 2026-11-09 16:17 — Dos partículas: estacionario y error de ψ/E en nodos con Δx = 2 y 1; kernels Wendland; prueba Hann cos⁴ en el capacitor
- **Objetivo:** (1) métricas de error en nodos (capas, normas L2/L∞, orden p) para dos cargas en la caja periódica, Δx = 2 y 1; (2) ver si los kernels Wendland (Jiang et al. 2013, App. A) sirven como kernel subgrid; (3) probar el capacitor con Hann cos⁴ de σ ≈ Hann-8.
- **Cambio:** corridas SLIDES-LD8-hann08-corr_no-N10000 (Δx = 1, 10000 pasos) y SLIDES-LD8-hann04-corr_no-dx2 (Δx = 2, Hann-4, mismas posiciones físicas; crear_2part_dx2.py). microgel/local/scripts/error_capas_2part.py reescrito para las dos redes, κ fijo 0.124975 (el "Kappa calculado" del log da 0 en r = 16: psi_electroneutral anula las densidades del nodo central ocupado). capacitor/local/WALL-DIRICHLET/crear_desplazada_p4.py: copias de δ = 1 con cos⁴ (n = 10 en Δx = 1, n = 5 en Δx = 2).
- **Prueba:** estacionario 32/32 OK (Δx = 1: dψ(9000→10000) ≤ 1.1e-8 vs umbral ≥ 1.5e-8; Δx = 2: dψ(2750→3000) ≤ 1.4e-9 vs ≥ 1.1e-7). Truncación de la referencia 2N→3N: φ ≤ 2e-6, E ≤ 1.1e-5. Wendland y cos⁴: sumas S₀ = Σw, S₁ = Σ(j−c)w en la red.
- **Resultado (observación):** error relativo L2 de ψ 4.6–6.1 % (Δx = 2) y 1.1–1.4 % (Δx = 1); de E 12–15 % y 3.2–3.9 %. p(2→1) global: ψ L2 2.10–2.18, L∞ 1.89–2.10; E L2 1.91–1.98, L∞ 1.75–1.86, casi independiente de r. Por capas: p de ψ ≈ 2 hasta d ≈ 15, pico ~4 en d ≈ 20 y cae a 0.3–0.7 para d ≳ 25; p de E ≥ 1.5 en todas las capas. Wendland C2/C4/C6 cumplen las condiciones continuas del paper pero S₀ no es exacta (soporte 8: 2e-3/5e-4/2e-4); Hann cos⁴ n = 8: S₀ exacta, máx|S₁| = 3.4e-4. cos⁴ con n = 10: σ = 1.4145 vs 1.4461 del Hann-8 (−2.2 %; el σ igual exige n = 10.22 no entero).
- **Interpretación (no verificada):** la caída de p de ψ lejos de las cargas podría deberse a un error casi constante (p. ej. la constante de ψ fijada por la media nula) que domina donde ψ_DH es chico.
- **Pendiente:** resultado de las corridas cos⁴ en δ = 1 (Δx = 2 lista, Δx = 1 en curso); verificar la interpretación anterior si el usuario lo pide; items viejos (modo cont, etiqueta nan de grafico_C_desplazada.py, descomposición (b)).

## 2026-11-09 16:19 — Capacitor con Hann cos⁴ (n = 10 / 5), δ = 1: primeras corridas
- **Objetivo:** probar Hann cos⁴ con σ ≈ Hann-8 en δ = 1, Δx = 2 y 1, y compararlo con cos².
- **Cambio:** C-electrolito-hann10p4-LD6-desplazada/z_18.50 y C-electrolito-hann5p4-LD3-ref0.5-desplazada/z_10.00 (crear_desplazada_p4.py); mismo exe que cos².
- **Prueba:** las dos terminaron normalmente; F_L entre los últimos dos registros cambia < 2e-13 (relativo ~1e-9).
- **Resultado (observación):** F_L (físico): cos⁴ Δx = 2 1.267496e-4, Δx = 1 1.247857e-4; cos² Δx = 2 1.267656e-4, Δx = 1 1.249585e-4. Contra F_kernel de cos² (1.243621e-4): cos² +1.93 % / +0.48 % (p = 2.01), cos⁴ +1.92 % / +0.34 %. La diferencia Δx = 2 → 1 es 1.96e-6 (cos⁴) vs 1.81e-6 (cos²).
- **Interpretación:** la comparación de cos⁴ contra F_kernel de cos² no es su error: el blanco de cos⁴ es ∫W₄ E_DH[W₄], distinto. Hace falta la referencia DH con fuente cos⁴.
- **Pendiente:** decidir si se extiende la referencia analítica (grafico_C_hoja_hann.py) a cos⁴ para obtener el error y p de cos⁴.

## 2026-11-09 16:27 — Referencia DH para Hann cos⁴ y error de la fuerza en δ = 1
- **Objetivo:** evaluar el error de la fuerza con Hann cos⁴ (ancho 10) contra su propio blanco continuo, y compararlo con cos² (ancho 8).
- **Cambio:** scripts/grafico_C_hoja_hann_pot.py (W, Ŵ, ∫e^{bs}W y W′ para cos^p como suma de armónicos mω; reemplaza esas funciones en HH y FD en memoria, sin tocar los archivos originales); scripts/grafico_C_fuerza_p4.py (evaluación, solo texto).
- **Prueba:** W vs cos^p/(n C_p) ≤ 1e-16, ∫W = 1, Ŵ y J vs cuadratura ≤ 1e-15, psit y ⟨U⟩ vs cuadratura (12 cifras), F_kernel vs forma por partes ≤ 6e-14; con cos² el módulo generalizado reproduce F_kernel original (1.4e-14); F_L recalculado desde efield con pesos cos⁴ coincide con Esub (≤ 4e-15); estacionario OK en las 4 corridas.
- **Resultado (observación):** F_kernel cos⁴ = 1.242062e-4 (cos²: 1.243621e-4). Error (F_L − F_kernel)/F_kernel: cos² +1.933 % (Δx = 2), +0.480 % (Δx = 1), p = 2.01; cos⁴ +2.048 %, +0.467 %, p = 2.13.
- **Interpretación:** en δ = 1 los dos kernels dan errores casi iguales. δ = 1 es nodo de Δx = 2 y punto medio de Δx = 1: una sola posición no dice cómo se comporta cos⁴ dentro de la celda.
- **Pendiente:** decidir si se corren más δ con cos⁴ (p. ej. la celda δ = 1–3) y/o Δx = 1/2.

## 2026-11-09 17:13 — ψ y E en nodos con cos⁴ vs cos² (δ = 1), convergencia en el borde del kernel; bug de emparejamiento de capas
- **Objetivo:** ver el orden p(2→1) por capas cerca del borde del soporte, cos⁴ (ancho 10) vs cos² (ancho 8), δ = 1.
- **Cambio:** scripts/grafico_C_error_capas_p4.py (nuevo; graficos/C_error_capas_p4.png). Bug corregido: las capas de dos redes se emparejaban por índice; si una capa no tiene nodos en una red (δ = 1: Δx = 2 sin nodos en 1 ≤ d < 2; δ = 12: sin nodos en d < 1) los órdenes quedaban corridos. Nueva función shell_order (empareja por centro) en grafico_C_error_capas_desplazada.py, usada también en grafico_C_error_capas_p4.py y microgel/local/scripts/error_capas_2part.py. Figura vieja guardada en el scratchpad (C_error_capas_desplazada_orden.emparejado_por_indice.png); se regenera.
- **Prueba:** chequeo FFT vs fórmula cerrada en el eje ≤ 1.4e-15 (las dos referencias); lista de capas vacías por red: en dos partículas solo falta la capa más lejana de Δx = 2 (fuera del emparejamiento, sin efecto).
- **Resultado (observación):** global p(2→1): ψ L2 2.17 (cos²) / 2.14 (cos⁴), E L2 1.93 / 1.93. Por capas: ψ casi igual para los dos kernels desde d ≈ 5; p de ψ baja a 1.46 en d = 5.5 (ambos) y 1.6–1.8 en d = 7–9, vuelve a 2 para d ≳ 10. E: en d = 4.5 p = 0.17 (cos²) vs 1.03 (cos⁴); con Δx = 1 el error relativo de E en d = 6–8 es ~2× menor con cos⁴ (2.6e-4 vs 6.7e-4 en d = 6.5). Capa d = 0.5: no comparable (Δx = 2 tiene el nodo en d = 0, donde E = 0).
- **Interpretación (no verificada):** la caída de p de E justo fuera de la cara del soporte cos² podría deberse a su menor suavidad (C¹ vs C³); una sola δ y capas de ancho 1 con pocos nodos de Δx = 2 no alcanzan para afirmarlo.
- **Pendiente:** revisar los órdenes corregidos de δ = 1 y 12 del análisis anterior; posible Δx = 1/2 con cos⁴ para un segundo par de redes.

## 2026-11-09 18:03 — Dos partículas: fuerza con Δx = 2 y 1 (N10000), doble conteo en el halo en Δx = 2
- **Objetivo:** ver cómo quedó la fuerza entre dos partículas con las corridas nuevas.
- **Prueba:** cociente F/(kT Esub) por partícula en las 32 corridas; F física = kT Esub_x/dx de la partícula móvil.
- **Resultado (observación):** en Δx = 2 (exe viejo 1c30df4b) la partícula 2 (fija en x = 4.25, celda de coloides del borde) tiene F = 2 kT Esub exacto en las 16 corridas: es el bug del halo corregido el 2026-09-26; la partícula 1 está bien. Δx = 1: todo bien. Como las partículas están fijas, ψ, E y Esub no se ven afectados. Fuerza a 10000 vs 3000 pasos (Δx = 1): diferencia relativa 2e-5 a 1.2e-3. F(Δx = 2) vs F(Δx = 1): −14 % (r = 2), −12 % (r = 1), −5.5 % (r = 5); para r ≥ 6 los r pares dan +0.14 a +0.22 % y los impares −1.8 a −3.4 %. −F₂/F₁ = 1 exacto en Δx = 1; en Δx = 2 hasta 1 + 3e-5 (solo r pares).
- **Interpretación (no verificada):** la alternancia par/impar en Δx = 2 podría venir de la posición dentro de la celda (r par: las dos cargas con offset 0.25; r impar: 0.25 y 0.75), como la dependencia con la posición vista en el capacitor.
- **Pendiente:** elegir la referencia para el error de la fuerza (DH Gauss-Gauss del análisis viejo o DH periódica con Hann continuo como en ψ/E).

## 2026-11-09 18:10 — Dos partículas: error de la fuerza vs distancia, Δx = 2 y 1, referencia DH periódica con Hann continuo
- **Objetivo:** fuerza entre dos cargas en función de r contra la referencia elegida por el usuario (opción 2: DH periódica, cargas Hann continuas de ancho 8, fuerza = campo promediado con el kernel).
- **Cambio:** microgel/local/scripts/error_fuerza_2part.py (nuevo) -> FUERZA-2PART-MEAS/FUERZA_LD08.0_Hann-8_corr-no_dx2-1.png. F_ref/kT = (4π l_B/V) Σ_k kx Ŵ3²/(k²+κ²) sin(kx r); F_L = kT Esub_x/dx de la partícula 1 (Esub, por el bug del halo en la partícula 2 de Δx = 2).
- **Prueba:** truncación M = 200 → 400: 1.5e-11 relativo; F_ref vs −dU/dr por diferencias centrales: 5.8e-10.
- **Resultado (observación):** Δx = 1: e = −5.3 % (r = 1), −2.8 % (r = 4), −0.23 % (r = 8), −0.09 % (r = 16), monótono. Δx = 2: −16.7 % (r = 1), −18.0 % (r = 2), −7.0 % (r = 5); para r ≥ 8 alterna: r par −0.09 a +0.13 %, r impar −1.9 a −2.6 %. p(2→1) = 1.6–2.1 para r ≤ 6; para r ≥ 7 oscila entre −3.3 y 4.2 (no interpretable: el error de Δx = 2 con r par es del orden o menor que el de Δx = 1).
- **Interpretación (no verificada):** la alternancia de Δx = 2 sería la dependencia con la posición dentro de la celda (r par: ambas cargas con offset 0.25; r impar: 0.25 y 0.75).
- **Pendiente:** corrida cos⁴ Δx = 1/2 en δ = 1 (capacitor), a decidir.

## 2026-11-09 18:47 — ¿Cuándo se cumple el estacionario en las corridas del capacitor con Δx = 1/2?
- **Objetivo:** ver si las corridas Δx = 1/2 (Hann-16, 20000 pasos) cumplían el criterio antes, para estimar los pasos necesarios en dos partículas con Δx = 1/2.
- **Prueba:** scratchpad/steady_history.py: max|ψ(n) − ψ(n−1)| y max|E(n) − E(n−1)| entre volcados consecutivos (intervalo t = 1000 en pasos de Δx = 1) contra el umbral (1 % del menor RMS por banda de Ludwig − DH en el eje), δ = 1, 4.7, 8.3, 12, Δx = 2, 1, 1/2.
- **Resultado (observación):** primer volcado que cumple: Δx = 1/2 en paso 16000 (t = 4000) para δ = 1, 4.7, 8.3 y 20000 (t = 5000) para δ = 12; Δx = 1 en t = 3000 (δ = 1) o 4000; Δx = 2 en t = 3000–4000. La deriva cae ~80× cada t = 1000 en el capacitor; en las corridas de dos partículas Δx = 1 caía ~3.3× cada 1000 pasos y cumplían recién entre 8000 y 10000.
- **Interpretación:** en el capacitor Δx = 1/2 necesitó el mismo tiempo físico que Δx = 1 (o un volcado más); no sobraron pasos salvo uno o dos volcados. Para dos partículas con Δx = 1/2 eso indica ~36000–40000 pasos (no verificado: otro sistema, sin paredes, λ_D = 8).
- **Pendiente:** decidir si se corren las 8 corridas de dos partículas con Δx = 1/2 (~10 h c/u).

## 2026-11-09 19:04 — Dos partículas Δx = 1/2, solo r = 4 (corriendo)
- **Objetivo:** ver en cuántos pasos llega al estacionario una corrida de dos partículas con Δx = 1/2, antes de decidir el resto.
- **Cambio:** microgel/local/crear_2part_dxhalf.py -> SLIDES-LD8-hann16-corr_no-dxhalf/r_04: caja 96³, Hann-16, ε 5e3, ρ_el 9.6485e-5 (λ_D = 16 nodos), cargas en x = 15.5 y 23.5, y = z = 47.5 (misma posición física), N = 20000, ψ cada 2000; ejecutable ACTUAL (src/Ludwig.exe 4f0c20e7, pedido del usuario), no el de Δx = 1 y 2 (1c30df4b).
- **Prueba (en curso):** scratchpad/monitor_dxhalf.py evalúa cada volcado: max|ψ(n) − ψ(n−2000)| vs 1 % del menor RMS por capa del error de ψ (mismo criterio que Δx = 1 y 2) y la fuerza kT Esub_x/Δx de la partícula 1.
- **Resultado (observación):** arrancó bien (log: Hann orden 16, cos²); ~0.95 s/paso medido en los primeros 650 pasos -> ~5.3 h, fin ~00:10.
- **Pendiente:** resultado del monitor; decidir las demás separaciones.

## 2026-11-09 23:58 — Δx=½, dos partículas, r=4: corrida terminada
- Objetivo: fuerza a r=4 con Δx=½ (SLIDES-LD8-hann16-corr_no-dxhalf/r_04, exe actual 4f0c20e7), 20000 pasos, ψ cada 2000.
- Observación: terminó normalmente (23:57). F_L(20000) = 3.215082e-7, e = −0.902 % vs F_ref = 3.244338e-7. Cambio relativo de F 16000→18000 ≈ 2e-6 y 18000→20000 ≈ 2e-6. particle_force = kT·Esub exactamente, y las dos partículas son iguales y opuestas (sin doble conteo de halo).
- Observación: el criterio de ψ (max|Δψ| < 1 % del menor RMS de error por capa) NO se cumplió: en 20000 va 2.87e-6 contra un umbral de 2.21e-8 (130×). max|Δψ| baja unas 2.6× cada 2000 pasos; el umbral se estabilizó en ~2e-8 a partir de 16000.
- Errores a r=4: Δx=2 −10.66 %, Δx=1 −2.84 %, Δx=½ −0.90 %; p(2→1) = 1.91, p(1→½) = 1.66.
- Interpretación (no verificada): F ya convergió en la cifra que importa; el p<2 de 1→½ podría venir de una parte no lineal de PB o de un término de orden menor, no está separado.
- Pendiente: decidir si se corren otras r con Δx=½, y si el criterio de estacionario es demasiado exigente con Δx=½.

## 2026-11-10 09:04 — Δx=½, dos partículas, r=5: lanzada
- Objetivo: fuerza a r=5 con Δx=½, 12000 pasos (pedido del usuario), ψ cada 2000.
- Cambio: crear_2part_dxhalf.py acepta los argumentos opcionales "d N" (sin argumentos hace lo mismo que creó r_04; copia del original en el scratchpad). Creada SLIDES-LD8-hann16-corr_no-dxhalf/r_05 (cargas en x = 15.5 y 25.5, y = z = 47.5, exe 4f0c20e7). Log en bateria_r05.log (bateria.log de r_04 sin tocar).
- Prueba: diff del input con r_04 (solo cambian el encabezado, N_cycles y vel_io_freq) y config.cds (solo la línea 36). El run arrancó sin errores.
- Pendiente: con r=4, F a 12000 pasos ya estaba a −0.898 % (final −0.902 %); comparar r=5 cuando termine (~3 h).

## 2026-11-10 09:10 — Δx=½, dos partículas, r=6: lanzada en paralelo con r=5
- Objetivo: fuerza a r=6 con Δx=½, 12000 pasos (como r=5), en paralelo si había memoria.
- Observación: GPU de 8188 MiB; r_05 ocupaba 1810 MiB y con r_06 el total es ~3.6 GB. Creada r_06 con "crear_2part_dxhalf.py 6 12000" (cargas en x = 15.5 y 27.5; diff con r_05 solo en la línea 36 de config.cds). Lanzada directamente (no con correr.sh) para no tocar r_05.
- Pendiente: comparar F de r=5 y r=6 con F_ref y con Δx=1 y 2 cuando terminen.

## 2026-11-10 09:16 — Δx=½, dos partículas, r=7: tercera corrida en paralelo
- Objetivo: ver si correr 3 juntas es más rápido que 2 y después 1 (corte: más de ~32 pasos/min por corrida).
- Cambio: r_07 creada con "crear_2part_dxhalf.py 7 12000" (cargas en x = 15.5 y 29.5) y lanzada directamente.
- Observación: con 3 juntas, ~43 pasos/min cada una (en 150 s); antes, 2 juntas a 64–65 y 1 sola a ~63. GPU con 6068 MiB y 84 %. Throughput total ~130 pasos/min vs ~128 con dos: la GPU ya estaba casi saturada con dos.
- Decisión: se deja (43 > 32). Fin estimado de las tres ~13:45–13:55, contra ~15:20 con 2 + 1; r_05 y r_06 terminan más tarde que antes (~13:45 vs 12:15).

## 2026-11-10 12:38 — Energía libre y potencial químico, dos partículas r=4 (pasos 0–2)
- Objetivo: verificar que el sistema llega al mínimo de F y se queda ahí (antes de evaluar U(r)).
- Cambio: nuevo microgel/local/scripts/energia_2part.py. F/kT = Σ[Σ_a ρ_a(lnρ_a−1) + ½ρ_e ψ] (funcional de fe_electro.c; ψ ya en kT/e), ρ_fix con d_hann; F_exc = F − entropía de los mismos iones uniformes; μ_a/kT = lnρ_a + Z_aψ (como calc_mui.py), Δμ = máx−mín. Salida: FUERZA-2PART-MEAS/ENERGIA_r04.csv.
- Observación (Δx = 2, 1, ½):
  - N_± constante a 1e-16 relativo; carga total 0 a 1e-14 (N₋ = N₊ + 2: el init uniforme neutraliza las cargas fijas).
  - F nunca sube.
  - Último ΔF entre volcados: 3e-12, 7e-12 y 5e-9 kT (relativo a F_exc: 8e-12, 2e-11, 2e-8); en Δx=½ todavía baja ~3× por volcado.
  - Δμ se estanca en (2.6, 4.0)e-6 kT (Δx=2) y (1.9, 2.2)e-6 kT (Δx=1); en Δx=½ sigue bajando: (2.0, 3.2)e-5 a t=5000, con Δμ₊ no monótono (1.4e-5 a 14000, 3.2e-5 a 16000).
  - F_exc final: 0.34852, 0.33025, 0.32576 kT.
- Interpretación (no verificada): los tres llegan a un mínimo, Δx=½ todavía no del todo. El piso de Δμ ~2e-6 kT no está explicado: podrían ser la precisión del volcado, la tolerancia del solver o el esquema NP.
- Pendiente: el usuario decide la tolerancia de meseta; origen del piso de Δμ; paso 3 (U(r) por ΔF y por integración de la fuerza, errores relativo y absoluto).

## 2026-11-10 12:55 — ¿El piso de Δμ viene del desfase ψ(n) ← qsi(n−1)?
- Objetivo: ver si el retraso de un paso entre ψ y qsi explica que Δμ se estanque en ~2e-6 kT.
- Prueba (sin cambios de código; script en línea): cambio por paso de ψ y de ln ρ_± estimado como la diferencia entre volcados dividida por los pasos entre ellos (cota del cambio en un paso si cada nodo relaja monótonamente), contra Δμ.
- Observación: en el piso, el cambio de ψ por paso es 4.6e-12 (Δx=2, paso 3000) y 9.2e-12 (Δx=1, paso 10000), contra Δμ = (2.55, 3.99)e-6 y (1.94, 2.20)e-6. Es 5–6 órdenes menor. Además ψ sigue cambiando ~3× menos por volcado mientras Δμ queda fijo. Con Δx=½, ψ cambia 1.4e-9 por paso con Δμ = 2e-5, que todavía está bajando.
- Interpretación: el desfase de un paso no explica el piso (supone relajación monótona por nodo).
- Pendiente: origen del piso; localizar dónde está el máximo de |μ − ⟨μ⟩| (propuesto, no hecho).

## 2026-11-10 12:59 — Perfil de δμ alrededor de las cargas, r=4, tres resoluciones
- Objetivo: ver si el piso de Δμ es el mismo en todos los refinados y cómo cambia μ alrededor de la partícula.
- Cambio: nuevo microgel/local/scripts/mu_capas_2part.py. δμ = μ − ⟨μ⟩ (media de todos los nodos); capas de distancia física a la carga más cercana, ancho 1; media con signo, RMS y máx |δμ| por especie. Volcados a t=5000 en los tres Δx y el último de Δx=2 (t=12000) y Δx=1 (t=10000).
- Observación a t=5000: los tres perfiles son casi iguales. Hay un gradiente suave a escala de la caja, de −1.07e-5 (μ₊) cerca de las cargas a +1.04e-5 en el punto más lejano (μ₋: −1.65e-5 a +1.59e-5). Δμ total: (2.20, 3.41)e-5, (2.20, 3.35)e-5 y (2.01, 3.17)e-5 para Δx = 2, 1, ½.
- Observación en el piso (Δx = 2 y 1): la media por capa es ~0 (≲1e-7) pero RMS ≫ |media| (estructura angular). La RMS es máxima a d ≈ 2.5–4.5, dentro del soporte del kernel (radio 4): 7.1e-7 y 5.7e-7 (μ₊), 1.28e-6 y 6.8e-7 (μ₋). max|δμ| ≈ 1.5e-6 y 1.25e-6 (μ₊); decae a ~5e-8 lejos. En Δx=1 queda además un resto de modo de caja de ~7e-8 en el punto más lejano.
- Interpretación (no verificada): el Δμ de Δx=½ a t=5000 es el transitorio del modo lento de caja, igual que en los otros Δx a ese tiempo; no es su piso. El piso está localizado alrededor de las cargas y con estructura angular. Con dos refinados baja solo 1.2–1.9× (no 4×). Hipótesis no probada: flujo estacionario pequeño (fuerza ρ_eE discreta no exactamente un gradiente) que advecta iones.
- Pendiente: el piso de Δx=½ requiere t ≳ 8000–10000 (32000–40000 pasos, ~9–10 h); r_05–07 solo llegan a t=3000. Chequeo propuesto: max|u| en vel-*.

## 2026-11-10 13:05 — ¿Hay flujo en el estado final? (r=4)
- Objetivo: chequear la hipótesis de que el piso de Δμ viene de un flujo estacionario.
- Prueba (sin cambios; script en línea sobre vel-*): Δx=2 paso 3000: max|u| = 2.15e-7, RMS 1.34e-8; Δx=1: max|u| = 7.18e-8, RMS 3.87e-9, idéntico en los pasos 5000 y 10000.
- Observación: hay un flujo estacionario chico (no decae). Ludwig tiene la clave "hydrodynamics no" (hydro_rt.c:42).
- Interpretación (no verificada): el flujo podría mantener el piso de Δμ; habría que comprobarlo con una corrida sin hidrodinámica.
- Pendiente: autorización para una prueba Δx=2 con hydrodynamics no (y otra igual con hidro y el mismo exe).

## 2026-11-10 13:08 — Prueba "hydrodynamics no" con partículas subgrid: no corre
- Objetivo: ver si Ludwig corre sin hidrodinámica (para probar si el piso de Δμ viene del flujo).
- Cambio: microgel/local/TEST-HIDRO-dx2/{hidro_si,hidro_no}: copia de SLIDES-LD8-hann04-corr_no-dx2/r_04 con el exe actual 4f0c20e7; hidro_no agrega "hydrodynamics no". Solo se corrió hidro_no.
- Observación: run.log dice "Hydrodynamics: off" y aborta en el primer paso: subgrid.c:252 Assertion 'hydro' failed (subgrid_force_from_particles, llamado desde ludwig.c:3943/4018; también subgrid_update en ludwig.c:2241 recibe hydro). No hay volcados.
- Interpretación: el camino subgrid supone hidrodinámica; para apagarla habría que modificar código (alcance desconocido, puede haber más puntos).
- Pendiente: elegir entre modificar ludwig.c (proteger con if (ludwig->hydro)) o una prueba sin cambiar código (viscosidad mayor → flujo menor).

## 2026-11-10 13:14 — r_06 (Δx=½) murió por falta de memoria en la GPU (error mío)
- Observación: r_06 abortó en el paso 10775 de 12000 con "Line 1451 (nernst_planck.c): cudaErrorMemoryAllocation: out of memory", 2 s después de que lancé las dos pruebas de viscosidad (2 × 904 MiB) con r_05, r_06 y r_07 corriendo (3 × 1810 MiB). Quedan sus volcados hasta el paso 10000 (t=2500) y Esub hasta 10775. r_05 y r_07 siguen corriendo.
- Causa: lancé las pruebas sin revisar la memoria libre; NP pide memoria temporal durante el paso.
- Pendiente: decidir si se relanza r_06 (sin borrar la salida parcial: en un directorio nuevo) o si alcanza con el paso 10000.

## 2026-11-10 13:19 — Prueba de viscosidad (Δx=2, r=4): el piso de Δμ baja con el flujo
- Objetivo: ver si el piso de Δμ viene del flujo estacionario.
- Cambio: TEST-HIDRO-dx2/hidro_si (η = 0.05) y hidro_si_eta0.5 (η = η_bulk = 0.5, τ = 2.375), exe 4f0c20e7, 3000 pasos. Script TEST-HIDRO-dx2/scripts/analizar_visc.py.
- Observación: hidro_si reproduce la corrida original (exe viejo) en F y Δμ con todas las cifras impresas.
  - η 0.05 → 0.5: max|u| 2.15e-7 → 4.01e-8 (5.4×), RMS|u| 1.34e-8 → 1.77e-9 (7.6×).
  - Piso de Δμ (2.55, 3.99)e-6 → (3.97, 5.83)e-7 kT (6.4×, 6.9×); RMS de δμ por capa 6–9× menor, con la misma forma (máximo en d ≈ 2.5–3.5).
  - F y ΔF iguales en las dos corridas a ≥ 10 cifras.
  - En el paso 1, max|u| = 4.44e-8 en las dos (independiente de η).
- Interpretación (no verificada): el piso de Δμ lo mantiene el flujo estacionario (baja ~ como RMS|u|). Que u no baje 10× podría ser porque parte de u en el volcado es el término de media fuerza F/(2ρ), que no depende de η.
- Pendiente: la causa del flujo (fuerza discreta no gradiente); qué u usa NP para advectar.

## 2026-11-10 13:33 — Δx=½, r=5 terminada (y r=6 parcial)
- Observación r=5 (12000 pasos, t=3000): F_L = 2.50111e-7, e = −0.590 % vs F_ref = 2.51595e-7. Error en los pasos 6000/8000/10000/12000: −0.486/−0.558/−0.581/−0.590 %. Con Δx=1: −1.622 %; Δx=2: −7.029 %. p(2→1) = 2.12, p(1→½) = 1.46.
- Observación r=6 (cortada por falta de memoria en el paso 10775): e = −0.320 % (paso 10000: −0.313 %), todavía bajando; Δx=1 −0.632 %, Δx=2 −1.972 %. p(1→½) ≈ 0.98 con el valor no estacionario.
- Interpretación: r=5 con 12000 pasos todavía se mueve ~0.01 % por cada 2000 pasos; r=6 no está convergida, no usar su p.

## 2026-11-10 13:41 — Δx=½, r=7 terminada
- Observación: 12000 pasos, F_L = 1.28146e-7, e = −0.197 % vs F_ref = 1.28399e-7. Error en los pasos 6000/8000/10000/12000: +0.057/−0.119/−0.177/−0.197 %, todavía bajando (−0.020 en los últimos 2000 pasos). Δx=1: −0.237 %; Δx=2: −3.636 %. p(1→½) = 0.27.
- Interpretación: no está estacionaria; e(Δx=½) es casi igual a e(Δx=1). No se puede dar orden con este valor.

## 2026-11-10 13:59 — Continuación de Δx=½ r=5, 6, 7 hasta el paso 20000
- Objetivo: llevar las corridas existentes a 20000 pasos (pedido: continuar, no corridas nuevas).
- Cambio: microgel/local/crear_continuacion_dxhalf.py (receta de capacitor/.../crear_continuacion.py). SLIDES-LD8-hann16-corr_no-dxhalf-cont/r_05 (desde 12000), r_06 (desde 10000, su último volcado), r_07 (desde 12000): copia de input con N_start/N_cycles, exe 4f0c20e7, metadata y solo los archivos del paso inicial; CSV de proceced_data recortados al paso inicial (r_06 tenía filas hasta 10775). Originales sin tocar. correr_cont.sh → run_cont.log.
- Prueba: las tres dicen "Re-starting simulation at step N0". r_06: Esub_x del paso 10001 en la continuación es idéntico al de la corrida original en las 16 cifras.
- Pendiente: terminan ~17:00 (r_05, r_07) y ~17:45 (r_06); después, error de fuerza a 20000 y p(1→½).

## 2026-11-10 16:10 — Barrido de viscosidad Δx=2, r=4 (preparado, en espera)
- Objetivo: tendencia de u y del piso de Δμ con η.
- Cambio: TEST-HIDRO-dx2/hidro_si_eta{0.02,0.1,0.2,1.0,2.0} (copias de hidro_si, solo cambian viscosity y viscosity_bulk; exe 4f0c20e7). Valores elegidos por el usuario; con 0.05 y 0.5 ya corridas son 7 puntos.
- Ejecución: correr_visc.sh (log barrido_visc.log) espera a que terminen las continuaciones r_05/r_07 y corre de a una, solo con ≥2000 MiB libres en la GPU (decisión del usuario: no arriesgar las continuaciones).
- Pendiente: analizar con analizar_visc.py extendido a las 7 corridas cuando terminen.
- 2026-11-10 16:13 Cambio de ejecución (pedido del usuario): correr_visc.sh ahora lanza en paralelo las que entren, dejando ≥2100 MiB libres (reserva = lo libre con las 3 continuaciones). Lanza una si libre ≥ 2100 + 904 MiB y espera 60 s entre lanzamientos. Ahora no entra ninguna (libre ~2300).

## 2026-11-10 17:12 — Barrido de viscosidad Δx=2, r=4: resultado
- Corridas: 5 nuevas (η 0.02, 0.1, 0.2, 1, 2), todas "finished normally" (17:03–17:11, hasta 3 en paralelo, libre ≥ ~2600 MiB). Continuación r_06 sin interrumpir.
- Análisis: scripts/analizar_visc.py extendido a 7 corridas + tabla resumen (original en scratchpad).
- Observaciones (paso 3000): F idéntica (10 cifras) en las 7. Δμ+ 6.03e-6 → 1.89e-7, Δμ- 9.48e-6 → 2.56e-7, max|u| 4.91e-7 → 2.12e-8, RMS|u| 3.28e-8 → 8.1e-10 de η=0.02 a 2. Pendiente local d ln/d ln η: Δμ+ ≈ -0.94 (bajo η) → -0.48 (alto η); max|u| -0.90 → -0.41; RMS|u| -0.97 → -0.48. u igual entre 2750 y 3000 en todas; Δμ+ de η=1, 2 cambió 10–25% entre 2500 y 2750 y <1% después.
- Interpretación (no verificada): el piso escala ~1/η a baja η (lo sostiene el flujo viscoso) y se aplana a alta η, podría haber una parte de u independiente de η (p.ej. el término F/(2ρ)); ojo: max|u| a η=2 (2.1e-8) ya está por debajo del max|u| de paso 1 (4.4e-8), así que ese valor de paso 1 no es el piso.
- Pendiente: decidir si se investiga el origen de esa parte independiente de η.
- 2026-11-10 17:14 Gráfico Δμ± vs η (log-log, paso 3000): TEST-HIDRO-dx2/graficos/dmu_vs_viscosidad.png, script scripts/grafico_dmu_visc.py. Valores idénticos a la tabla de analizar_visc.py. Sin ajuste ni rectas de referencia.
- 2026-11-10 17:17 Gráfico de dos paneles (decisión del usuario: no doble eje), eje η compartido: Δμ± arriba, max|u| y RMS|u| abajo, paso 3000. TEST-HIDRO-dx2/graficos/dmu_u_vs_viscosidad.png, script scripts/grafico_dmu_u_visc.py. El gráfico anterior (dmu_vs_viscosidad.png) se conserva. Valores iguales a la tabla de analizar_visc.py. Observación: las 4 curvas se aplanan a η ≳ 0.5.

## 2026-11-10 17:35 — Continuaciones Δx=½ r=5,6,7 hasta 20000: resultado
- Las 3 "finished normally" con línea de reinicio (17:35). r_06 paso 10001 idéntico al original.
- Error de fuerza (partícula 1, F = 1e-5·Esub_x/0.5 vs F_ref) en 20000: r=5 −0.595 %, r=6 −0.334 %, r=7 −0.209 % (r=4: −0.902 %). Cambio 18000→20000 < 0.001 puntos en las tres (convergido en el tiempo a ese nivel).
- p(1→½) con Δx=1 (−1.622, −0.632, −0.237 %): r=5 1.45, r=6 0.92, r=7 0.18 (r=4 1.66).
- Interpretación (no verificada): el orden observado cae con r; en r=7 Δx=1 y Δx=½ casi coinciden → otro error domina (no de discretización espacial), p.ej. de la referencia o de resolución temporal/caja.

## 2026-11-10 19:25 — Δx=½, distancias 3 y 2 (pedido del usuario)
- Creados SLIDES-LD8-hann16-corr_no-dxhalf/r_03 (cargas x = 15.5 y 21.5) y r_02 (15.5 y 19.5) con crear_2part_dxhalf.py d 20000 (y = z = 47.5; exe 4f0c20e7; 20000 pasos como r_04). Input y config difieren de r_04 solo en cabecera y posición de una carga.
- r_03 lanzada con correr.sh 1 (bateria_r03.log; las otras r_* se saltean). r_02 lanzada directo (no con correr.sh, para no duplicar r_03 que aún no tiene psi). Las dos en paralelo, ~3.6 GB GPU.

## 2026-11-11 00:01 — Δx=½, distancias 1 y 8 (pedido del usuario)
- Creados r_01 (cargas x = 15.5 y 17.5) y r_08 (15.5 y 31.5) con crear_2part_dxhalf.py d 20000; exe 4f0c20e7; difieren de r_04 solo en cabecera y posición de una carga.
- Espera a que terminen r_02/r_03 (~01:45) y luego correr.sh 2 (bateria_r01_r08.log): corre r_01 y r_08 en paralelo; saltea las terminadas y no toca r_06 (salida parcial). ETA ~6.5 h más.

## 2026-11-11 00:13 — Δx=1, distancias intermedias 4.1…4.9 (pedido del usuario)
- Script crear_2part_dx1_intermedios.py → SLIDES-LD8-hann08-corr_no-N10000-d4a5/pos_12.10 … 12.90 (9 corridas). Copias de pos_12.00 (48³, λD=8, Hann-8, 10000 pasos) cambiando solo la x de la carga móvil (la otra en x=8; y=z=24).
- Ejecutable: el actual 4f0c20e7 (elección del usuario), NO el 1c30df4b de las corridas Δx=1 de d=4 y 5 (no se rehacen). Ojo al comparar con esas.
- Se lanzan solas al terminar r_01/r_08 (Δx=½), con correr.sh 3 (3 en paralelo; elección mía). Log bateria.log.

## 2026-11-11 01:47 — Δx=½ r=2 y 3: resultado
- Las dos "finished normally" (~01:46). r_01 y r_08 arrancaron después.
- F_ref (error_fuerza_2part.F_ref, M=400; M/2 cambia <1e-12): r=2 3.258525e-7, r=3 3.629434e-7.
- Error en 20000 (Δx=2 / 1 / ½): r=2 −18.021 / −4.785 / −1.376 %; r=3 −13.649 / −3.944 / −1.173 %. Cambio 18000→20000 < 0.001 puntos.
- p(2→1) / p(1→½): r=2 1.91 / 1.80; r=3 1.79 / 1.75 (r=4: 1.91 / 1.66).
- Observación: a r ≤ 4 el orden es ~1.7–1.9 en ambas parejas; cae a r ≥ 5. Δx=2 y 1 usan exe viejo, Δx=½ el actual.

## 2026-11-11 08:12 — Δx=½ r=1 y 8: resultado; arrancan las intermedias Δx=1
- r_01 y r_08 "finished normally" (08:10). F_ref: r=1 1.956777e-7, r=8 9.314642e-8.
- Error en 20000 (Δx=2 / 1 / ½): r=1 −16.662 / −5.339 / −1.506 %; r=8 −0.093 / −0.229 / −0.189 %. Cambio 18000→20000: r=1 < 0.001, r=8 −0.001 puntos.
- p(2→1) / p(1→½): r=1 1.64 / 1.83; r=8 −1.31 / 0.28 (en r=8 Δx=2 da MENOS error que Δx=1: no monotónico).
- Las 9 intermedias Δx=1 (d 4.1–4.9) arrancaron 08:11, de a 3 (3.1 GB GPU).
- 2026-11-11 08:51 Gráfico de fuerza con Δx=½ (pedido del usuario): error_fuerza_2part.py ahora además guarda FUERZA-2PART-MEAS/FUERZA_LD08.0_Hann-8_corr-no_dx2-1-half.png (puntos Δx=½ r=1–8, paso 20000, r=5–7 desde -cont) en paneles de fuerza y error; panel p sin cambios. El png anterior se regenera idéntico (mismo md5). Original del script en scratchpad.

## 2026-11-11 09:06 — Ludwig sin hidrodinámica con partículas subgrid (pedido del usuario)
- Cambios (CHANGE INIT/END 20261111):
  - subgrid.c subgrid_force_from_particles: si hydro == NULL, return después de lubricación + colloid_sums_halo (no se reparte la fuerza al fluido).
  - subgrid.c subgrid_update_forces_electrokinetics: hydro_memcpy y hydro_f_local_add protegidos con if (hydro); la fuerza sobre la partícula no cambia.
  - ludwig.c: stats_total_force_write (total_force.csv) solo si hay hydro (segundo assert que apareció).
- No protegidas (no se usan en estas corridas): variantes _ewald/_theory, pm_sr (subgrid.c ~5961), interlacing ≠ 0 en ludwig.c:1692.
- Exe: 36be1759 (solo subgrid.c; hidro_no abortó en stats_total_force_write) → 45dc94bc (definitivo). Backup de subgrid.c y del exe 4f0c20e7 en scratchpad. Las corridas d4a5 en curso usan su copia de 4f0c20e7.
- Prueba corta 300 pasos sin hidro (scratchpad): "finished normally".
- Pruebas completas: TEST-HIDRO-dx2/hidro_si_exe45dc94bc (debe reproducir hidro_si) e hidro_no_exe45dc94bc, correr_exe45dc.sh, en espera de memoria GPU.
- 2026-11-11 09:13 Resultado prueba sin hidro (exe 45dc94bc):
  - Obs: hidro_si con 36be1759 y con 45dc94bc idénticos byte a byte a hidro_si (psi, qsi, vel 3000 y particle_Esub.csv).
  - Obs: hidro_no_exe45dc94bc "finished normally"; N+ igual; F nunca sube, dF igual a las corridas con hidro (−3.2e-12 en 3000).
  - Obs: Δμ± sin hidro baja sin piso, ~×2.85 cada 250 pasos: 1.15e-6 (2000) → 1.76e-8 / 1.81e-8 (3000). Con η=2 se planta en 1.89e-7 / 2.56e-7 desde 2750. Hasta ~2000 las dos curvas coinciden.
  - Obs: Esub_x part. 1 a 3000: sin hidro 5.7970491e-2, η=2 5.7970494e-2, η=0.05 5.7970540e-2 (diferencia relativa ≤ 8.5e-7).
  - Interpretación: el piso de Δμ lo sostiene el flujo; sin flujo el sistema sigue relajando hacia μ uniforme. El efecto del flujo sobre la fuerza es < 1e-6 relativo.
  - Pendiente: con 3000 pasos no se ve aún el piso sin hidro (si lo hay); se podría continuar hidro_no.

## 2026-11-11 09:25 — Sin hidro: continuación 3000 → 6000 (pedido del usuario)
- TEST-HIDRO-dx2/hidro_no_exe45dc94bc-cont (copias del paso 3000, N_start 3000, N_cycles 3000; exe 45dc94bc). "Re-starting at step 3000" + "finished normally".
- Obs: Δμ± baja de forma exponencial sin piso, factor ~2.85 cada 250 pasos (≈ e-fold cada ~240 pasos) hasta 7.6e-14 / 7.8e-14 kT en 6000.
- Obs: ΔF entre volcados = 0 o ±2.3e-13 desde 3750 (= 1 ulp de F≈1049: redondeo). N± constantes.
- Obs: Esub_x part. 1 fijo en 5.7970491027e-2 desde 4250 (cambia <1e-12 relativo).
- Obs: residuo final de δμ suave a escala de caja (RMS 3.6e-14 cerca de las cargas, mínimo 7e-15 a d≈20–24, sube a 3.6e-14 en el extremo de la caja), sin estructura localizada en el kernel.
- Interpretación: sin flujo el sistema llega al mínimo de F (μ uniforme) hasta precisión de máquina en F; el piso con hidro se debía al flujo. Nota: la observación de "watcher" disparado por "Note assertions" (error mío, corregido).

## 2026-11-11 09:32 — Sin hidro: continuación 6000 → 9000 (pedido del usuario)
- hidro_no_exe45dc94bc-cont2 (copias del paso 6000; exe 45dc94bc). Reinicio en 6000 + "finished normally".
- Obs: Δμ± deja de bajar: 7.6e-14 (6000) → 4.44e-14 / 4.35e-14 desde 6250–7000 y constante hasta 9000.
- Obs: archivos psi y qsi idénticos byte a byte entre volcados desde 7500 (cambian hasta 7250). F idéntica en las 16 cifras escritas; Esub_x part. 1 = 5.7970491027126e-2 constante. N± constantes.
- Obs: residuo final de δμ ~1e-14–2e-14 RMS, suave a escala de caja.
- Interpretación: se alcanzó un punto fijo exacto del esquema en doble precisión; el piso ~4e-14 kT es redondeo (≈ 200 ε en μ; ε·|ln ρ| ≈ 1.2e-15 por nodo, más la acumulación en ψ). No hay piso físico/numérico sin flujo por encima de eso.

## 2026-11-11 09:37 — Sin hidro con Δx=1, Hann-8 (pedido del usuario)
- TEST-HIDRO-dx1/hidro_no_exe45dc94bc: copia de SLIDES-LD8-hann08-corr_no-N10000/pos_12.00 (r=4, 48³) con hydrodynamics no, exe 45dc94bc, N_cycles 30000 (elección mía: Δx=2 llegó al punto fijo en ~7500 pasos, ×4 por D igual por paso), freq_config 1000.
- Corre en paralelo con las intermedias d4a5 (3 × 1014 MiB); libre antes 4621 MiB.

## 2026-11-11 09:46 — Δx=1, d = 4.1…4.9: resultado
- Las 9 "finished normally" (08:11–09:46).
- Error de fuerza (partícula 1, paso 10000, vs F_ref): 4.0 −2.843 (exe viejo), 4.1 −2.738, 4.2 −2.631, 4.3 −2.522, 4.4 −2.413, 4.5 −2.302, 4.6 −2.187, 4.7 −2.067, 4.8 −1.936, 4.9 −1.790, 5.0 −1.622 % (exe viejo).
- Obs: monótono y suave entre 4 y 5, sin saltos en los extremos con exe distinto; el paso por 0.1 crece de −0.105 a −0.168 puntos. No se verificó la convergencia temporal de cada intermedia (solo el último paso).
- 2026-11-11 09:53 Gráfico de fuerza + intermedias Δx=1 (pedido del usuario): FUERZA-2PART-MEAS/FUERZA_LD08.0_Hann-8_corr-no_dx2-1-half-int.png (círculos vacíos, paneles fuerza y error), generado por error_fuerza_2part.py. Los dos png anteriores se regeneran idénticos (mismo md5).

## 2026-11-11 10:00 — Δx=2 y Δx=1, d = 8.2…9.8 cada 0.2 (pedido del usuario)
- Motivo (usuario): esperaba fluctuación de la fuerza con la posición dentro de la celda; Δx=2 muestra saltos par/impar en el error y Δx=1 no.
- crear_2part_d8a10.py → SLIDES-LD8-hann04-corr_no-dx2-d8a10/r_08.2…r_09.8 (desde r_08; carga móvil 4.25 + d/2, pasos de 0.1 en red) y SLIDES-LD8-hann08-corr_no-N10000-d8a10/pos_16.20…17.80 (desde pos_16.00; pasos de 0.2 en red). 18 corridas, exe 45dc94bc (= 4f0c20e7 con hidro). d=9.0 incluido (también rehecho con exe actual; sirve para comparar con el viejo).
- correr_d8a10.sh (log bateria_d8a10.log): paralelo con reserva ≥2100 MiB libres; Δx=2 primero.

## 2026-11-11 10:35 — Sin hidro Δx=1, Hann-8, r=4: resultado (30000 pasos)
- TEST-HIDRO-dx1/hidro_no_exe45dc94bc "finished normally".
- Obs: Δμ± baja exponencialmente, factor ~2.85 cada 1000 pasos (igual que Δx=2 cada 250 pasos: mismo tiempo físico), hasta ~1.8e-13 desde el paso ~23000 y constante hasta 30000. Con hidro el piso era (1.94, 2.20)e-6.
- Obs: ψ/qsi NO quedan idénticos byte a byte entre volcados (a diferencia de Δx=2), pero Δμ y Esub_x (3.1520927279710e-2, 14 cifras) constantes desde 22000–23000. ΔF ≤ 4.5e-13 (≤ 2 ulp de F≈1408). N± constantes.
- Obs: Esub_x part. 1 con hidro (exe viejo, paso 10000) 3.1520943714e-2 vs sin hidro 3.1520927473e-2 en 10000: diferencia relativa 5.2e-7.
- Obs: residuo final RMS δμ 2e-14–9e-14, suave a escala de caja.
- Interpretación: mismo resultado que Δx=2: sin flujo se llega al mínimo de F hasta redondeo; piso de redondeo ~4× mayor (8× más nodos). El flujo no es un artefacto de una sola resolución.
- 2026-11-11 10:52 Gráfico ψ y E de la corrida sin hidro Δx=1 (pedido del usuario; formato elegido: perfiles en el eje + mapas 2D, con referencia DH): TEST-HIDRO-dx1/graficos/psi_E_hidro_no_paso30000.png, script TEST-HIDRO-dx1/scripts/grafico_psi_E_hidro_no.py (usa dh_nodes de error_capas_2part.py). Obs: en el eje y=z=24, max|ψ_L − ψ_DH| = 1.16e-2 (máx ψ_DH 0.477, en el centro entre cargas, ψ_L por encima), max|E_x,L − E_x,DH| = 5.2e-3 (máx 0.0887, en los picos).
- 2026-11-11 11:04 Gráfico ψ/E sin hidro: referencia ahora evaluada en grilla fina (2401 puntos en x, suma de modos sobre la recta y=z=24, M=3N), línea fina (lw 1). Chequeo: en los nodos coincide con dh_nodes (constante uniforme a 4.8e-16, E_x a 7.5e-17). Mismo archivo png (sobrescrito, es un gráfico, no un resultado).
- 2026-11-11 11:26 Gráfico ψ/E sin hidro: agregada fila de error punto a punto en la recta y=z=24 (Ludwig − referencia en los nodos, con signo, absoluto). Obs: error de ψ ≈ +0.0115 casi constante entre x=8 y 12, cae a 0 hacia x≈5 y 15, pequeño negativo (~−6e-4) justo afuera; error de E_x antisimétrico: +0.0051 en x=6, −0.0051 en x=14 (2 nodos hacia afuera de cada carga), ≈0 entre las cargas y lejos (|x−X| > 4, fuera del soporte del kernel).
- 2026-11-11 11:32 Gráficos ψ/E con hidro (pedido del usuario): scripts/grafico_psi_E_2part.py {hann8_hidro, hann16_hidro} → FUERZA-2PART-MEAS/PSI_E_hann8_hidro.png (pos_12.00, exe viejo, paso 10000) y PSI_E_hann16_hidro.png (dxhalf r_04, paso 20000). Unidades físicas. Hann-16: la recta de nodos usada es y=z=23.75 (el eje y=z=24 no tiene nodos; elección del script).
  - Obs: Hann-8 con vs sin hidro: max|Δψ| 8.2e-8, max|ΔE| 2.6e-8 (eE/kT) → gráficos indistinguibles; max error en el eje igual (ψ 1.162e-2, E_x 5.156e-3).
  - Obs: Hann-16: max error ψ 2.25e-3 (Hann-8: 1.16e-2, ×5.2), E_x 1.35e-3 (×3.8). Misma forma: ψ meseta positiva entre cargas, E_x antisimétrico con picos ~2 unidades físicas afuera de cada carga. (Recta desplazada 0.25 del eje en Hann-16: comparación aproximada.)

## 2026-11-11 11:34 — d = 8…10 cada 0.2, Δx=2 y Δx=1: resultado
- 18 corridas OK (11:34). Error de fuerza (partícula 1, último paso) vs F_ref:
  d:    8.0    8.2    8.4    8.6    8.8    9.0    9.2    9.4    9.6    9.8   10.0
  Δx=2: −0.093 −0.141 −0.584 −1.273 −2.028 −2.632 −2.844 −2.397 −1.217 −0.382 −0.046 %
  Δx=1: −0.229 −0.002 −0.056 −0.237 −0.357 −0.204 −0.001 −0.053 −0.217 −0.324 −0.182 %
  (8.0, 10.0 exe viejo; 9.0 con exe nuevo da lo mismo que el viejo: −2.632 / −0.204.)
- Obs: oscilación con período = 1 celda de la red de cada resolución (2 físico en Δx=2, 1 en Δx=1). Amplitud pico a pico ~2.8 % (Δx=2) y ~0.36 % (Δx=1). Δx=2: mínimo |e| con la carga móvil a fracción ~0.25 de celda (x = 8.25, 9.25), máximo ~0.85. Δx=1: |e|≈0 con la carga a fracción 0.2 de nodo, máximo a 0.8.
- Interpretación: en Δx=1 todas las d enteras caen en la misma fase de celda (carga sobre nodo) → la oscilación no se veía; en Δx=2 las d enteras alternan fase 0.25/0.75 → saltos par/impar. El "orden observado" p de las d enteras mezcla fases distintas entre resoluciones y no es un orden de convergencia limpio.
- Pendiente: decidir si se grafica; si se estudia la dependencia con la posición de la carga fija (siempre en x=8 / 4.25).
- 2026-11-11 11:41 PSI_E_hann8/16_hidro: el "punto suelto" en x≈11 era el marcador de la leyenda sobre los datos (mi cambio a loc="upper right"). Leyenda movida afuera, encima de cada panel; regenerados y revisados los dos png.
- 2026-11-11 11:47 Etiqueta de la referencia unificada a "continuous Hann (width 8)" también en grafico_psi_E_hidro_no.py; la etiqueta más larga hizo que loc="best" pusiera la leyenda sobre la curva (confirma la causa en los otros dos) → leyenda afuera, encima del panel, como en los otros. Los tres gráficos ψ/E con el mismo formato de leyenda.
- 2026-11-11 11:49 Leyenda de Ludwig en los 3 gráficos ψ/E: ahora indica el kernel discreto, "Hann-N (N nodes = width 8, Δx)" (Hann-8 en Δx=1, Hann-16 en Δx=½), a pedido del usuario.
- 2026-11-11 11:53 Gráfico de fuerzas (pedido: el último, sin hacer uno nuevo): FUERZA_..._dx2-1-half-int.png ahora incluye d = 8.2…9.8 para Δx=1 (círculos vacíos, misma serie que 4.1–4.9) y Δx=2 (triángulos vacíos naranjas), en paneles de fuerza y error. Valores = tabla anterior. Los otros dos png del script se regeneran idénticos.
- 2026-11-11 11:55 Nuevo gráfico (pedido del usuario): FUERZA-2PART-MEAS/FUERZA_LD08.0_Hann-8_corr-no_abs-rel.png — fuerza / error absoluto F_L − F_ref (con signo) / error relativo, mismas series que -half-int, sin panel de orden. Generado al final de error_fuerza_2part.py; los otros 3 png se regeneran idénticos.
  - Obs: el error absoluto es máximo en r=2 para las tres resoluciones (Δx=2 −5.9e-8, Δx=1 −1.6e-8, Δx=½ −4.5e-9) y no en r=1 como el relativo.
- 2026-11-11 11:59 Nuevo gráfico zoom 8 ≤ r ≤ 10 (pedido del usuario): FUERZA-2PART-MEAS/FUERZA_LD08.0_Hann-8_corr-no_abs-rel_zoom8-10.png (fuerza / error absoluto / error relativo, solo puntos en la ventana, sin líneas que unan puntos, grilla vertical cada 0.2). Generado al final de error_fuerza_2part.py; los otros png sin cambios.
- 2026-11-11 12:13 Zoom 8–10 (pedido del usuario): una sola serie por red (enteros e intermedios juntos, mismo símbolo/color), sin Δx=½, sin aclaraciones de ejecutable. Mismo archivo png.
- 2026-11-11 12:14 Zoom 8–10: leyenda dentro del panel de fuerza (esquina superior derecha, libre de datos), a pedido del usuario.

## 2026-11-11 12:27 — Δx=2, r=4 fijo, corrimiento del par dentro de la celda (pedido del usuario)
- crear_2part_dx2_r4_shift.py → SLIDES-LD8-hann04-corr_no-dx2-r4-shift/d_0.0…d_2.0: copias de r_04 con las DOS cargas corridas en x δ físico = 0…2.0 cada 0.2 (δ/2 en red; elección del usuario: solo x, 11 puntos con control δ=0 ≡ δ=2). Exe 45dc94bc, 3000 pasos. correr.sh 4.
- Script de gráfico: scripts/grafico_fuerza_shift_r4_dx2.py (F en ambas cargas vs δ, con F_ref(4); error relativo).
- 2026-11-11 12:32 FUERZA_..._dx2-1-half.png sin el panel de orden (pedido del usuario): ahora figura propia de 2 paneles (fuerza, error relativo) en error_fuerza_2part.py; mismo archivo. Los demás png del script sin cambios (md5 iguales).
- 2026-11-11 12:37 Resultado corrimiento del par (Δx=2, r=4): 11 corridas OK. FUERZA-2PART-MEAS/FUERZA_shift_r4_dx2.png.
  - Obs: error relativo entre −11.60 % (δ=0.4–0.6, ambas cargas cerca de mitad de celda: x≈4.5/6.5 red) y −9.53 % (δ=1.4–1.6, ambas cerca de nodo: x≈5.0/7.0); δ=0: −10.66 % (igual a r_04 viejo). Pico a pico ≈ 2.1 puntos.
  - Obs: control δ=0 vs δ=2.0 idéntico (3e-16 relativo). F1 = −F2 a ≤ 5e-6 relativo. Simetría: F1(δ) = −F2(1.0−δ) (reflexión especular: valores de 0.0/1.0, 0.2/0.8, 0.4/0.6, 1.2/1.8, 1.4/1.6 iguales con cargas intercambiadas).
  - Interpretación: a r fijo, la fuerza depende de la posición del par en la celda con período de una celda; el valor de la tabla de r enteros (δ=0) es un punto de esa curva, no su promedio.
- 2026-11-11 12:41 FUERZA_shift_r4_dx2.png: 3er panel (F − ⟨F⟩)/⟨F⟩ en %, ⟨F⟩ = media simple sobre δ = 0…1.8 (un período; δ=2.0 excluido por ser la misma configuración que δ=0; elección mía). Obs: ⟨F⟩ = 2.9001e-7 (−10.61 % vs F_ref); variación entre −1.09 % (δ=0.4–0.6) y +1.21 % (δ=1.4–1.6); δ=0 y 1.0 ≈ −0.05 %.
- 2026-11-11 12:45 FUERZA_shift_r4_dx2.png: kT pasado al título ("k_BT = 1e-5 (lattice units)", título en 2 líneas); eje y "force (lattice units)". La fuerza NO está en unidades de kT: F = kT·Esub_x/Δx en unidades de red (F/kT = valor × 1e5). La etiqueta "(k_BT = 10^-5)" sigue en los gráficos de error_fuerza_2part.py (no cambiada).
- 2026-11-11 12:48 error_fuerza_2part.py: en los 5 gráficos de fuerza, eje y "F_x on charge 1 (lattice units)" y "k_BT = 1e-5 (lattice units)" en el título (pedido del usuario); título del zoom en 2 líneas.
- 2026-11-11 13:17 Error mío en grafico_psi_E_2part.py: los mapas usaban extent (0, L) en vez de (0.5, L+0.5) → mapa corrido medio nodo respecto de las cruces y flechas (por eso PSI_E_hann8_hidro se veía distinto de psi_E_hidro_no). Corregido; regenerados PSI_E_hann8_hidro.png y PSI_E_hann16_hidro.png. Los datos con/sin hidro difieren ≤ 8e-8; ahora los mapas coinciden a la vista.
- 2026-11-11 13:22 psi_E_hidro_no_paso30000.png: nueva fila (3ª) con la solución analítica (DH periódico, Hann continuo ancho 8) de cada carga por separado y la suma, ψ y E_x en la recta y=z=24 (pedido del usuario). Constante de ψ: mitad de la del total para cada carga (elección mía). Chequeo de superposición: |carga1 + carga2 − total| ≤ 2.4e-16. Obs: cada carga sola tiene ψ máx 0.323 y |E_x| máx 0.072, picos de E a ~2.3 de su carga.
- 2026-11-11 13:23 psi_E_hidro_no_paso30000.png: título general bajado (suptitle y=0.915) para quitar el espacio vacío sobre las leyendas.
- 2026-11-11 13:24 PSI_E_hann8_hidro.png y PSI_E_hann16_hidro.png: título general bajado (suptitle y=0.925) en grafico_psi_E_2part.py, como en el gráfico sin hidro.
- 2026-11-11 13:29 dmu_u_vs_viscosidad.png: título y ejes pasados a inglés (pedido del usuario; regla general: textos de gráficos en inglés, guardada en memoria). grafico_dmu_visc.py (dmu_vs_viscosidad.png) todavía tiene textos en español.
- 2026-11-11 13:31 dmu_vs_viscosidad.png: título y eje x pasados a inglés (grafico_dmu_visc.py).
- 2026-11-11 14:32 FUERZA_..._dx2-1-half.png (pedido del usuario): referencia más fina (lw 1); línea punteada en r = σ/2 = 0.723 (σ = 1.4460 = desvío estándar de la ventana Hann continua de ancho 8, σ² = n²(1/12 − 1/(2π²)); posición σ/2 elegida por el usuario); leyenda para las dos líneas verticales (r=8 y σ/2); el texto de la línea r=8 se sacó del eje x.
- 2026-11-11 14:36 Gráficos de fuerza vs r (error_fuerza_2part.py: dx2-1, -half, -half-int, abs-rel): referencia fina (lw 1) y líneas verticales con leyenda en r = σ/2 = 0.723 (punteada), r = σ = 1.446 (raya-punto; el usuario la define como separación mínima físicamente posible sin superposición real) y r = 8 (rayas). Función kernel_lines(). En el zoom 8–10 no se agregan (σ/2 y σ fuera de rango). Leyenda del panel de orden p arriba a la izquierda (tapaba puntos).
- 2026-11-11 14:39 Leyenda de la línea r = σ simplificada a "r = σ = 1.446: kernel std. dev." (pedido del usuario), en los 4 gráficos de fuerza vs r.
- 2026-11-11 14:44 Gráficos de fuerza vs r: líneas en múltiplos kσ dentro del ancho del kernel (kσ < 8 → k = 2…5: 2.892, 4.338, 5.784, 7.230), raya-punto gris claro, una sola entrada de leyenda (pedido del usuario). En kernel_lines(), aplica a los 4 gráficos.
- 2026-11-11 14:55 Zoom 8–10: puntos unidos con líneas (ordenados en r), a pedido del usuario.
- 2026-11-11 14:57 FUERZA_shift_r4_dx2.png: solo fuerza sobre la carga 1, puntos unidos con líneas; leyenda solo en el panel de fuerza (pedido del usuario). La carga 2 se sigue leyendo e imprimiendo en la tabla.
- 2026-11-11 15:13 Nuevo gráfico zoom 4 ≤ r ≤ 5, solo Hann-8 (Δx=1): FUERZA-2PART-MEAS/FUERZA_LD08.0_Hann-8_corr-no_abs-rel_zoom4-5.png (fuerza / error absoluto / relativo, puntos unidos). Obs: el error absoluto sube casi linealmente de −9.2e-9 (r=4) a −4.1e-9 (r=5), incrementos ~0.5e-9 por 0.1 sin oscilación visible; en 8–10 el mismo Δx=1 oscila con período 1 (amplitud ~0.26e-9 abs.). Interpretación pendiente.
- 2026-11-11 15:16 Zoom 4–5: serie Δx=1 en azul (#2a78d6) en vez de negro (pedido del usuario).

## 2026-11-11 16:59 — Teoría: iones explícitos vs implícitos en tamaño de polielectrolito lineal
- Objetivo: entender por qué las simulaciones con iones explícitos (colaborador, Alemania) dan Rg, R_ee y R_h menores que el modelo implícito, y cómo modelarlos vs N.
- Cambio: ninguno (solo explicación teórica). Fuente verificada por búsqueda: Stevens & Kremer, J. Phys. II France 6, 1607 (1996): en régimen diluido sin sal, las cadenas con DH están más estiradas que las de Coulomb con contraiones explícitos.
- Pendiente: definir qué es el "modelo implícito" (DH lineal par a par o campo medio PB/NP de Ludwig), parámetros (b, f, l_B, sal, concentración, rango de N) y si hay datos del colaborador para comparar.
- 2026-11-11 17:05 Corrección del usuario: el modelo implícito es Poisson–Boltzmann NO lineal (no DH). Consecuencia: la condensación de Manning ya está en campo medio (Fuoss–Katchalsky–Lifson, PNAS 37, 579 (1951)); la diferencia con iones explícitos debe venir de lo que PB omite. Fuente verificada (resumen): Deserno, Holm & May, Macromolecules 33, 199 (2000): MD da condensación consistentemente más fuerte que PB. Pendiente: parámetros (ξ, valencia, sal), cómo resuelve la red la zona cercana a la cadena.
- 2026-11-11 22:36 Lectura Netz & Orland, EPJE 1, 203 (2000) (12 pp, leído completo). Observación: corrección a un lazo (RPA) alrededor de PB; ion siente auto-energía (1/2)q² G_self(r) (ecs. 30, 43, 57-61), test particles con término Q²G_self/2q² atractivo (ec. 70); parámetro Σ=2πq³l_B²σ_s, inválido para Σ≳12; solo contraiones, sin sal (con sal requiere núcleo duro, p. 213); solo resuelto para pared plana. En Ludwig el lugar natural sería mu_solv en nernst_planck.c (hoy 0 para fe_electro; caminos GPU/PETSc asumen mu_solv=0, líneas ~1197-1307). Nada implementado.
- 2026-11-11 22:52 Usuario: iones monovalentes y con sal. Consecuencia: Netz-Orland (solo contraiones) no aplica directo. Estimación heurística (fórmula de pared plana aplicada a superficie de cilindro, NO validada): Σ ≈ l_B²/(a A); con l_B=0.7 nm, A=0.25 nm, a=0.5 nm da Σ≈4 → correlaciones no despreciables para cadenas muy cargadas aun monovalentes. Pendiente: ξ, A, a, κ, Δx, kernel del modelo implícito y del colaborador.
- 2026-11-11 22:57 Lectura Netz, EPJE 5, 557 (2001) (18 pp, completo). Observación: Ξ=2πq³l_B²σ_s (ec. 4) es el mismo Σ de Netz-Orland 2000; μ=1/(2πqσ_s l_B) (ec. 2). PB = orden cero de expansión en Ξ (válido Ξ<−1/(z̃ ln z̃) cerca, ec. 59; Ξ<z̃/ln z̃ lejos, ec. 60); SC = virial, ρ̃=e^{-z̃} (ec. 94), válido Ξ>z̃ cerca / Ξ>z̃² lejos (ecs. 97-98); hueco intermedio sin teoría (fig. 7). Atracción entre placas iguales en SC, MC umbral Ξ≈10 (p. 570). Solo contraiones, sin sal, geometría plana. Interpretación: no da corrección implementable en Ludwig; sirve como diagnóstico (calcular Ξ, μ, a⊥ vs Δx). Ejemplo ilustrativo (a=0.5, A=0.25 nm, NO datos del usuario): Ξ≈4, μ≈0.18 nm, PB válido a z ≳ 9μ ≈ 1.6 nm. Pendiente: parámetros reales.
- 2026-11-11 23:13 Lectura Kornyshev, J. Phys. Chem. B 111, 5545 (2007) (13 pp, completo). Observación: PB de red con volumen excluido (Bikerman/Eigen-Wicke/BAO): c± = c0 e^{∓u}/(1−γ+γ cosh u), γ=2c0/c_max (ecs. 5,7); Poisson-Fermi (ec. 9); C(u0) analítica (ec. 20) → GC si γ→0 (ec. 21), alas C≈C0/√(2γ|u0|) (ec. 22). Sigue siendo campo medio: sin correlaciones ni oscilaciones; cita que la corrección estérica MF no mejora frente a MC (Antypov-Barbosa-Holm PRE 71, 061106 (2005)). Interpretación: no explica cadenas más chicas (el efecto estérico limita la condensación → iría en sentido opuesto). Relevante para EDLC (secundamiento): da el mu_solv estérico μ=−kT ln(1−Σc/c_max) que falta en Ludwig y C(u0) ec. 20 como test analítico. Estimación heurística (a=0.7 nm, 0.1 M): γ≈0.04, saturación desde u0≈ln(2/γ)≈4 (~0.1 V). Nada implementado.
- 2026-11-11 23:20 Lectura Messina, "Electrostatics in soft matter", arXiv:0811.4578 (2008) (33 pp, completo; revisión). Observación: PB ≡ g_αβ=1 en jerarquía YBG (ecs. 12-13): ignora correlaciones laterales electrostáticas y estéricas. p. 8: para monovalentes acuosos PB e incluso DH concuerdan "astonishingly well" con simulaciones (refs 21-26, incl. Deserno-Holm-May). Sec. 2.3: discretitud de carga solo importa a T finita para σ baja y multivalentes. Sec. 3.1: con monovalentes GRANDES, volumen excluido + acople débil favorece correlaciones y sobrecarga (no es campo medio). Sec. 4: imágenes de corto alcance (~tamaño de ión). Sec. 5: adsorción de PE (no tamaño de cadena aislada). Interpretación: no da método implementable; sugiere que con monovalentes la diferencia del colaborador podría venir del tamaño de los iones explícitos o de diferencias de modelo/resolución, más que de correlaciones puras. Pendiente: parámetros (incl. diámetro de iones del colaborador).
- 2026-11-11 23:30 Búsqueda bibliográfica extensiones del acople fuerte (solo resúmenes/búsqueda, ningún texto completo leído). Encontrado: Naji-Kanduč-Forsman-Podgornik perspectiva JCP 2013 (arXiv 1307.5130); dressed counterions (arXiv 0912.5460, sal mono + multivalentes); Šamaj-Trizac Wigner SC (arXiv 1009.4640, 1112.3763), régimen fluido (arXiv 1804.08459), placas asimétricas PRE 110, 014609 (2024), paredes polarizables arXiv 2607.08851 (2026); cilindros Naji-Netz (cond-mat/0504447), aguja (arXiv 1212.0430), Sci Rep 2017. Afines: Agrawal-Wang teoría de fluctuaciones gaussiana renormalizada (arXiv 2212.00141; extendida a cepillos de polímero, resúmenes APS 2024/2026); Bazant-Storey-Kornyshev PRL 106, 046102 (2011), Poisson de 4º orden. Interpretación: SC estricto no aplica a monovalentes; lo más aplicable a cadenas con sal sería Agrawal-Wang; lo más implementable en Ludwig, BSK.
- 2026-11-11 23:38 Descarga a "/home/danolv/Desktop/papers a ordenar": 13 PDFs (arXiv/abiertos; JCP 2012 de FU Berlin devolvió HTML → reemplazado por arXiv 1208.3579; perspectiva 2013 y Messina ya estaban, no se tocaron). Leídos completos: Agrawal-Duan-Wang (arXiv 2212.00141), Bazant-Storey-Kornyshev PRL 2011, Naji-Netz cond-mat/0504447. Resto: resumen + conclusiones (perspectiva: también inicio sec. IV). Observaciones: Agrawal tabla 1: Ξ<1 PB sin sobrecarga, 1<Ξ<100 correlaciones, Ξ>100 + volumen excluido; fig. 3 monovalentes sin sobrecarga aun a σ alto; autoenergía local de corto alcance analítica (ec. 16) + parte larga numérica (ec. 13). BSK: Poisson 4º orden ε(l_c²∇²−1)∇²φ=ρ (ec. 2), BC extra n·∇∇²φ=0; μ con estérico (ec. 5) apto para NP. Naji-Netz: Ξ definido en superficie del cilindro; ξc=1 y exponentes de MF no cambian con correlaciones, pero densidad cerca del cilindro sube con Ξ. Kanduč 2012: simulación con iones implícitos (mono) reproduce explícitos, cuantitativo a sal alta. Dressed/Šamaj-Trizac/aguja/SciRep: multivalentes o solo contraiones, no aplicables a monovalentes con sal. Interpretación: SC no aplica al caso del colaborador; piezas útiles: calcular Ξ del cilindro real; ec. 16 de Agrawal como mu_solv local (parcial); BSK para EDLC (en FFT periódico sería cambiar el núcleo a 1/(k²(1+l_c²k²))). Nada implementado.
- 2026-11-11 23:46 Usuario: a futuro puede usar iones multivalentes. Reevaluación (sin implementar): Ξ∝q³ (×8 divalente, ×27 trivalente) → régimen intermedio/fuerte donde PB falla cualitativamente. Relevantes entonces: Agrawal-Wang (teoría unificada con mezclas mono+trivalente, validada vs experimento fig. 4b); dressed counterions/Kanduč 2012 (mono implícito + multi explícito ≈ todo explícito) → ruta híbrida posible en Ludwig: multivalentes como partículas subgrid cargadas + monovalentes como campos NP (sugerencia, requiere validar kernel vs μ y ruido térmico); Naji-Netz ecs. 2-3 (p_MF, p_SC alrededor de cilindro) como tests analíticos; Adžić-Podgornik si hay regulación de carga. Hsiao-Luijten PRL 97, 148301 (2006) (colapso/re-expansión de PE con sal, citado en Agrawal, no leído).
- 2026-11-11 23:51 Descargado y leído completo Hsiao & Luijten, PRL 97, 148301 (2006) (arXiv cond-mat/0610164) → "papers a ordenar". Observación: MC bead-spring (Stevens-Kremer), iones explícitos LJ, λ_B=3σ, N=16-96, Cm=0.008σ⁻³. Fig. 1: R_g² vs C_s: 1:1 baja lenta y monótona; Z:1 (Z≥2) colapso con mínimo en C_Z (neutralización, ∝ Cm) y re-expansión. Fig. 2 (4:1): ν(Cs) de >0.588 a 0.34±0.02 en C_Z y vuelve a ≈0.588. Fig. 3: sobrecarga Q(r) para 3:1 con iones de diámetro σ; con σ/2 casi no hay (pares de Bjerrum). Fig. 4: inversión de movilidad. Interpretación: tamaño iónico es parámetro crítico; ν(Cs) da directamente R_g∝N^ν para modelar vs N; en PB (sin correlaciones ni tamaño) no aparecería colapso/reexpansión con multivalentes (no verificado en Ludwig).

## 2026-11-12 11:55 — Test del colaborador (ESPResSo): polímero cargado con contraiones implícitos
- Objetivo: entender el test (charged_polymer.pdf) y ubicarlo en la teoría leída.
- Observaciones (del documento): monómeros MD q=+1, contraiones campo continuo deriva-difusión (EK), ℓB=2.84[x] (7.1 Å), grid 1[x], c=4.7e-5 (5 mM), z=−1, NO figura sal agregada; N=3..50; ν(Re): explícito 0.84 / continuo ≈1.05; ν(Rg): 0.81 / ≈0.98.
- Cálculos propios (no del doc, heurísticos): L=40..102 ✓ ec. 5; ξ=ℓB/A≈2.8–3.2 → ~65 % condensado (Manning); Ξ=2πℓB²σ_s≈8 (radio 1) a 17 (radio 0.5) → régimen de cruce 10–100 (Naji-Netz); μ≈0.17–0.34[x] < grid 1[x]; κ⁻¹≈24–42[x], comparable a L.
- Interpretación: aun monovalente, Ξ~10 → correlaciones plausibles; además capa condensada subresuelta por la grilla (mismo problema que en Ludwig). ν>1 es exponente efectivo en rango corto de N (b_eff crece con N), no asintótico.
- Pendiente: confirmar si hay sal; tamaño de iones en referencia explícita; cómo se acopla carga de monómeros a la grilla EK y si monómero–monómero va además por P3M (posible doble conteo).
- 2026-11-12 12:18 Esquema de acople del colaborador (ek_particle_coupling.pdf, anterior): partícula–partícula por P3M (Coulomb puntual); E_salt por FFT de ρ_salt, interpolado a partículas con Hann; sal siente E_tot de ρ_tot=ρ_salt+ρ_part(Hann). Observación: sin doble conteo; mismo kernel para repartir y recoger. Interpretación: asimetría — repulsión monómero–monómero sin suavizar, atracción monómero–contraión suavizada (Hann, σ=0.18·n [x]; n=3 en la figura → 0.54 [x], > μ≈0.2–0.3) → apantallamiento cercano subestimado → cadena más extendida, consistente con ν alto. Mismo esquema que mis pruebas Hann en Ludwig. Pendiente: n real, refinamiento.
- 2026-11-12 13:07 Búsqueda (3 consultas, solo resúmenes): ¿publicado refinar grilla+exclusión o híbrido contraiones explícitos+sal continua? No encontrado tal cual para cadenas. Cercanos: Kuron et al. JCP 2016 (partículas cargadas como fronteras móviles en EK de ESPResSo, iones excluidos de esferas resueltas); arXiv 2203.04603 (PB modificado para PE con iones de tamaño finito, monovalente); arXiv 1905.04960 (hinchamiento de geles: explícito → PB); Kanduč 2012 (multivalentes explícitos + monovalentes DH, MC sin grilla); Tischler et al. (EK termalizada en ESPResSo, no verificado). "No encontrado" ≠ inexistente.
- 2026-11-12 13:16 Observación: la EK termalizada es de Ingo Tischler (1er autor): Tischler, Weik, Kaufmann, Kuron, Weeber, Holm, J. Comput. Sci. (2022) 101770, doi 10.1016/j.jocs.2022.101770; ChemRxiv 10.26434/chemrxiv-2021-39nhv (v3) (Crossref). Descargados a "papers a ordenar": Kuron 2016 (arXiv 1607.04572), arXiv 2203.04603, arXiv 1905.04960. Fallaron (bloqueo Cloudflare/PMC, no se forzó): Tischler 2022 (ChemRxiv), Yin et al. JCP 2008 (doi 10.1063/1.2991179), Rathee et al. Polymers 2019 (doi 10.3390/polym11010183). No leídos aún.
- 2026-11-12 13:29 Lectura 6 papers (carpeta movida a ~/Descargas/papers a ordenar). Completo: Tischler et al. 2022 (ChemRxiv v3). Parcial (resumen, modelo, conclusiones): Vahid et al. arXiv 2203.04603 (SPB), Landsgesell et al. arXiv 1905.04960, Kuron et al. 2016, Yin et al. JCP 2008, Rathee et al. Polymers 2019. Observaciones: Tischler sec. IV: 3 acoples de partículas (puntual por interpolación: "accuracy and discretization artifacts strongly depend on the interpolation scheme and the grid spacing"; raspberry; frontera móvil); acople a EK termalizada = trabajo futuro; ruido recortado si n<1/a³ (ec. 44) → en su test (4.7e-5 iones/celda en bulk) el ruido se anularía. Kuron: Debye ≥ ~4 celdas, R=3–8 celdas; ×2 resolución = ×8 costo. SPB: PB + WCA blando rod–ión (1 parámetro) reproduce iones explícitos alrededor de varilla rígida λ=−3.6 e/nm (ξ≈2.6) hasta 1 M. Landsgesell: PB celda con varilla penetrable reproduce hinchamiento de gel MD explícito (ξ>1), límites multivalentes/alta carga. Rathee: DH implícito falla a alta sal mono y con divalentes. Yin: divalentes condensan preferentemente. Interpretación: opción 1 = refinar + exclusión tipo SPB; EK termalizada no resuelve el test.
- 2026-11-12 13:43 Creado /mnt/datos/Github/Trabajo/UB/Proyectos/Bibliografia_resumen.xlsx (hoja Papers: 29 artículos, columnas Proyecto/Importancia/Resumen/DOI/Año/Autores/Título; hoja Leyenda). DOIs, años y autores verificados con Crossref. Corrección: Agrawal-Duan-Wang publicado en J. Phys. Chem. B 2023 (10.1021/acs.jpcb.3c04739), no JCTC; Šamaj-dos Santos-Trizac ya publicado en JCP (10.1063/5.0343388). Importancia y etiquetas de proyecto = criterio propio, a revisar.

## 2026-11-12 14:54 — Ruido térmico en electrocinética (problema de Tischler) y qué haría falta en Ludwig
- Objetivo: entender por qué Ingo dice que la termalización falla a baja densidad y qué implicaría agregar ruido al NP de Ludwig.
- Fuente leída: Tischler et al. (chemrxiv v3), ecs. 40–45 (sección IV). Flujo fluctuante por enlace: j = sqrt(n̄_enlace)·sqrt(2D)·W, normalizado con 1/sqrt(a_i A0 Δt) (ec. 43). Si n<1/a³, la amplitud se multiplica por H(n)=n (ecs. 44–45) para evitar densidades negativas.
- Código revisado (sin cambios): nernst_planck.c. En GPU el flujo difusivo se calcula por enlace y después se antisimetriza (np_antisym_kernel, procesa cada par una vez desde el índice menor y con periodicidad local, o sea que solo vale para un rank). El NP es explícito con multisteps.
- Observación (cuenta): en el test de Ingo hay ~4.7e-5 contraiones por celda en el bulk. Estimación heurística: ~0.1 por celda en la capa condensada. Con el recorte de Tischler el ruido queda prácticamente apagado en todo el dominio.
- Interpretación: el problema no es numérico sino de descripción. Con N≪1 ion por celda, un campo continuo con ruido gaussiano no representa la estadística discreta. Para tener N≥1 hace falta a ≳ c^(-1/3) ≈ 28 [x], comparable con λD (24–42 [x]). Ruido gaussiano alrededor de PB daría a lo sumo la corrección a un lazo (Netz–Orland), insuficiente con Ξ~10.
- Corrección a una memoria: el NP de Ludwig es explícito, así que la ventaja frente a Péraud 2016 (límite Δt<λD²/D) que había anotado no está respaldada.
- Pendiente: nada implementado. Opciones discutidas: no agregar ruido, ruido con recorte, iones discretos o híbrido.

## 2026-11-12 15:17 — ¿La separación de Ingo (P3M partícula–partícula + EK para la sal) vale en PB no lineal?
- Separación: Poisson es lineal, así que φ_tot = φ_part + φ_sal es exacta. La no linealidad de PB está en la respuesta de los iones, y los iones sienten φ_tot. Con el mismo kernel para esparcir y recoger, el modelo es conservativo y no tiene autofuerza. Conclusión: sí resuelve un PB no lineal, pero de un modelo con la atracción monómero–ion suavizada (Coulomb convolucionado con Hann).
- Observación (Tischler ecs. 34 y 36): el flujo de deriva usa densidad promedio × ΔΦ. Con flujo cero, el cociente de densidades entre nodos es (1+x/2)/(1−x/2) en lugar de e^x (x = βzeΔΦ por enlace). Error: +1 % con x=0.5, +10 % con x=1, +56 % con x=1.5, y diverge con x=2.
- Ludwig: el flujo es D·½(1+b0)(ρ1·b1 − ρ0) con b = exp(Δμ), cuyo equilibrio es la distribución de Boltzmann exacta para cualquier ΔΦ.
- Estimación (cota superior, sin suavizado ni condensación): una línea de carga con ℓB=2.84 y b=1 da ΔΦ ≈ 3.9 kT entre r=1 y r=2.
- Interpretación: la discretización sesgaría la condensación hacia arriba y el suavizado hacia abajo. No sé cuál domina.
- No verificado: qué versión de ESPResSo y qué esquema EK usa Ingo (el EK nuevo, basado en waLBerla, puede diferir del paper).

## 2026-11-12 15:26 — ¿Es consistente un modelo con P3M suavizado con el mismo kernel y el mismo σ?
- Pregunta del usuario: si la interacción partícula–partícula usara el mismo kernel Hann y el mismo σ que el acople con EK, ¿el modelo sería consistente?
- Razonamiento: un modelo es consistente si cada tipo de par usa un potencial definido por una única imagen de la carga, U_ab = z_a z_b ℓB (K_a⊗K_b⊗1/r). "Mismo kernel en P3M" (esparcir y recoger) da partícula–partícula = K⊗K (ancho √2σ), mientras que partícula–ion = K una sola vez (el ion como valor de grilla). No es la misma forma.
- Opciones: (A) iones también como nubes K (carga = K*n, los iones sienten K*φ), así todos los pares son K⊗K; (B) definir explícitamente qué forma se quiere y usar la misma en la referencia explícita.
- Observación clave: con potenciales suavizados, la ν de Ingo (cargas puntuales) deja de ser la referencia. Habría que rehacer la simulación explícita con los mismos potenciales suaves.
- Sugerido, no hecho: comparar tres corridas (explícita puntual, explícita suave, continuo suave) para separar el efecto del suavizado del error de campo medio.
- Detalles: si K es radial y de soporte compacto, K⊗K⊗1/r = 1/r para r > 2·soporte (teorema de la cáscara), así que se puede implementar como P3M puntual más una corrección tabulada de corto alcance. Un kernel separable (producto 1D) no es radial: la forma es anisótropa y depende de la posición respecto de la grilla. Con grillas distintas, la igualdad vale solo hasta el error de discretización.

## 2026-11-12 15:48 — Respuesta de Ingo con resultados (con y sin advección, kernels CIC/Hann-3/Hann-5)
- Observaciones (leídas de sus figuras):
  - ν explícito en esta corrida: Re 0.86, Rg 0.79, Rh 0.39.
  - Con advección: CIC 1.01/0.91/0.44, Hann-3 1.00/0.88/0.39, Hann-5 1.01/0.90/0.44. Las curvas de desviación de los 3 kernels coinciden casi al píxel: en N−1≈32, Re +42 %, Rg +32 %, Rh +14 %.
  - Sin advección: CIC 0.93/0.84/0.41, Hann-3 0.92/0.82/0.36. En N−1≈32, Re +16–18 %, Rg +12–13 %, Rh +5 %.
  - Ingo confirma que con la densidad del test el ruido EK es prácticamente nulo. No probó refinar la grilla ni cambiar el P3M. No respondió por la sal ni por el tamaño de los iones explícitos.
- Interpretación (hipótesis):
  - (1) Apagar la advección reduce la discrepancia unas 2.5 veces. Si las dinámicas muestrearan la misma distribución de equilibrio, la estática no debería depender del transporte, así que algo en el acople advectivo rompe el equilibrio. Candidatos: difusión numérica del esquema upwind, o agitación por la velocidad térmica del LB sin ruido iónico que la compense. Ninguno verificado.
  - (2) Que el kernel no influya es consistente con que la condensación de Manning en PB dependa poco del radio. Eso debilita, sin refutarla, la hipótesis de subresolución del suavizado de la partícula.
  - (3) El residuo sin advección (0.93 contra 0.86) podría venir de correlaciones (Ξ≈8–17) o de la grilla de iones. Solo refinar la grilla lo distingue.
- Pendiente: decidir qué responderle a Ingo.

## 2026-11-12 16:15 — Fluctuaciones térmicas del LB frente al transporte EK, con los parámetros de Ingo
- Parámetros (PDF charged_polymer): kT=1, ρ=1, ν=1, a=1, δt=0.01 [t], D=6.515e-2 [x]²/[t], ℓB=2.84. El τ del LB no figura (supuse τ=δt para pasar a unidades de red).
- Aclaración del usuario: Ingo no probó cargas no puntuales en el P3M.
- Cuentas (equipartición por nodo; estimaciones de orden de magnitud):
  - u_rms = √(kT/ρa³) = 1 [x]/[t] por componente.
  - Péclet de celda u·a/D ≈ 15.
  - Velocidad de deriva EK v = D·|∇βeφ| = 0.02–0.2 → u_th/v ≈ 5–50.
  - Escala de difusión advectiva tipo Stokes–Einstein kT/(6πηa) ≈ 0.05, del orden de D.
  - Difusión numérica de un upwind ~⟨|u_x|⟩a/2 ≈ 0.4 ≈ 6 D (forma genérica de upwind de primer orden; no verifiqué el CTU de ESPResSo).
  - En unidades de red (τ=0.01): kT=1e-4, u=0.01, ν=0.01, D=6.5e-4.
- Ludwig (np_advective_fluxes_kernel): advección centrada, u_enlace·½(ρ0+ρ1). Sin difusión numérica, pero con Péclet de celda ≈ 15 > 2 un esquema centrado explícito puede oscilar o dar densidades negativas (criterio estándar, no probado aquí).
- Interpretación (hipótesis): las fluctuaciones térmicas del LB dominan el flujo iónico instantáneo. Si el esquema de advección agrega difusión sin la movilidad correspondiente, los iones quedan efectivamente "más calientes", hay menos condensación y la cadena se estira. Sería consistente con ν→1 con advección. No verificado.

## 2026-11-12 16:22 — ¿Las fluctuaciones del LB superan el límite drmax=0.8 LU de las partículas subgrid?
- Código: colloids.c:1074 corta la corrida si dr[ia] > drmax (0.8 por defecto). subgrid_update: dr = u_interpolada + drag·F + frand, con frand gaussiano truncado en 3σ.
- Cuentas con los parámetros de Ingo, suponiendo 1 paso LB = τ:
  - τ=0.01 [t]: u_th=0.01 LU/paso, patada aleatoria rms=0.036, fuerza crítica ≈1230 kT/[x] (WCA con r<0.70, U≈67 kT; FENE con r>1.47, U≈112 kT).
  - τ=0.1 [t]: u_th=0.1, patada 0.11, fuerza crítica ≈120 kT/[x] (WCA con r<0.83, U≈7 kT; FENE con r>1.25, U≈40 kT).
- Interpretación: con τ=0.01 el ruido térmico queda 20–80 veces por debajo del límite, y solo una configuración inicial mala lo dispara. Con τ=0.1 el WCA a ~7 kT sí es alcanzable por fluctuaciones. Es un riesgo real, aunque raro.
- Observación lateral: el chequeo es dr > drmax, sin valor absoluto; un desplazamiento negativo grande no lo detecta. No lo cambié.

## 2026-11-12 16:32 — Corrida de 1 partícula con los parámetros de Ingo (LB térmico, contraiones advectados)
- Objetivo: ver cuánto mueve el ruido térmico a una partícula subgrid por paso (frente a drmax=0.8) con los parámetros de Ingo.
- Setup: microgel/local/INGO-1PART-LBTERM.
  - Unidades: 1 LU = 0.01 [t]; kT=1e-4, η=0.01, ρ=1; D_ion=6.515e-4 (input D/A0, D3Q27); ℓB=2.84 (ε=280.2024).
  - Hann-8; ah=0.3079 y al=0.4951 (drag = 1/γ de Ingo).
  - Caja 27³; sin sal (contraiones por electroneutralidad, 5.08e-5 por nodo); posición inicial al azar en la celda central (14.085, 14.303, 13.563).
  - Decisiones del usuario: τ, kernel, movilidad, advección.
- Incidentes: el primer intento abortó en el paso 1 por falta de ./proceced_data/ (la creé). El segundo se cortó en el paso 100 por la aserción rho >= 0 (fe_electro.c:215).
- Observaciones:
  - En el paso 100 la densidad de contraiones tenía mínimo −1.50e-5 (−30 % de la media), máximo 1.52e-4 (3× la media); carga total conservada.
  - Desplazamiento de la partícula por paso: rms 0.036–0.042 LU por componente (la patada pura esperada es 0.036); máximo 0.100 LU.
  - Densidad del fluido: desviación estándar 0.0172 (la equipartición predice √(3kT) = 0.0173).
- Interpretación:
  - drmax no es un problema con este τ (8 veces de margen sobre el máximo observado en 100 pasos).
  - La densidad negativa es compatible con la advección centrada a Péclet de celda ≈15, pero no está aislada: falta una corrida sin hidrodinámica o sin fluctuaciones para separar la causa.
  - En una compilación con -DNDEBUG la corrida no se habría detenido y la densidad negativa habría pasado sin aviso.
- Pendiente: decidir si separar la causa de la densidad negativa.

## 2026-11-12 16:41 — Diagnóstico de la densidad iónica negativa (microgel/local/INGO-1PART-LBTERM-DIAG)
- Objetivo: encontrar por qué la densidad de contraiones se vuelve negativa en INGO-1PART-LBTERM.
- Corridas (un cambio cada una, estadísticas en cada paso):
  - BASE: densidad negativa en el paso 72.
  - NOHYDRO y NOFLUCT: 1000 pasos sin problema, resultados idénticos (n entre 5.07e-5 y 6.61e-5).
  - FIXED (partícula fija): idéntica a BASE.
- Observación: la causa es la advección de los iones por la velocidad térmica del LB. El movimiento de la partícula y la electromigración no intervienen.
- Código (nernst_planck.c:1213 y la versión CPU 1062): flujo por enlace ½(u0+u1)·c · ½(ρ0+ρ1), sumado sobre los 26 enlaces sin peso. Por Taylor da (M/2)∇·(ρu) con M = Σc_x² (D3Q27: factor 9, D3Q19: 5, D3Q7: 1). La difusión tiene el factor A0, que el input compensa; la advección no lo compensa.
- DUMP (paso 50): δn/n con desviación estándar 0.246 (mínimo −0.84), δρ/ρ del fluido 0.017, correlación 0.05. En el paso 60, 0.279 (crecimiento ~√t).
- DUMP7 (D3Q7, factor 1): desviación estándar 0.077 en el paso 50 y 0.085 en el 60, correlación 0.01–0.02.
- Interpretación:
  - En el continuo, sin difusión, n/ρ se conserva a lo largo de las trayectorias, así que δn/n debería seguir a δρ/ρ (~0.017, correlación ~1).
  - Las fluctuaciones iónicas no correlacionadas y crecientes indican que el operador advectivo de NP (velocidades nodales centradas) no es el transporte de masa del LB. Las velocidades térmicas de escala de grilla actúan como ruido multiplicativo sin restauración.
  - El factor 9 lo amplifica (unas 3 veces menos dispersión con D3Q7), pero no es la única causa.
  - Por linealidad, el efecto relativo no debería depender de la densidad media (no probado).
- Pendiente: decidir si corregir el factor de normalización de la advección y/o reformular el flujo advectivo para que use el flujo de masa del LB. No implementado.

## 2026-11-12 16:51 — Advección de Ludwig frente a Capuani 2004
- Fuente: Capuani, Pagonabarraga y Frenkel, JCP 121, 973 (2004), §III.A.2 (p. 976–977, ec. 16 y Fig. 1) y §VI.A (ecs. 30–32).
- Su esquema de advección no usa flujos por enlace. Desplaza la celda del nodo r en v(r)Δt y reparte la densidad entre las celdas que intersecta, en proporción al volumen (volume of fluid). Dicen explícitamente que lo hacen para evitar interpolar velocidades a los enlaces.
- Propiedades: conserva la masa, no da negativos si |v_α|Δt ≤ 1, y tiene difusión espuria D* = ½ Σ_α v_α(1−v_α) (ec. 31).
- Ludwig: el código upstream (K. Stratford; la versión previa al commit 20c21db33 lo llama "'Centred difference' advective fluxes ... TODO: belongs to Nernst Planck") usa flujos centrados por enlace. Antes de 2023 el stencil era de compilación, con D3Q7 por defecto (NP_D3Q18 y NP_D3Q26 como opción). Con D3Q7 el factor M/2 es 1, así que la hipótesis del usuario sobre la normalización es consistente con la historia. No puedo afirmar que nunca se probó con 19 o 27.
- Cuenta: con la velocidad térmica de Ingo (⟨|v_α|⟩ ≈ 0.008 LU/paso), Capuani da D* ≈ 0.004 por dirección ≈ 6 D_ion. En el equilibrio de deriva–difusión (D + D*), la temperatura efectiva de los iones sería ≈ 7T.
- Interpretación: ESPResSo usa este esquema (Tischler cita [1,38]). Sería un mecanismo concreto para el ν ≈ 1 de Ingo con advección. Copiar Capuani tal cual en Ludwig arreglaría los negativos pero introduciría el mismo efecto. No verificado.

## 2026-11-12 17:04 — Corrección de la normalización del flujo advectivo de NP (D3Q19/D3Q27)
- Cambio: src/nernst_planck.c (CPU np_advective_fluxes; GPU kernel np_advective_fluxes_kernel y driver). Cada enlace pesa W_p = 2·wgradients[p], con lo que ½ΣW c c = I, el mismo criterio que la corrección del laplaciano del 20251203. Valores: D3Q27 4/9, 1/9, 1/36; D3Q19 1/3, 1/6; D3Q7 1 (sin cambio). Marcado como CHANGE 20261112. El ejecutable anterior quedó copiado en el scratchpad (Ludwig.exe.pre_advnorm).
- Elección del peso: el usuario eligió 2·wgradients frente a 2/M uniforme (iguales a primer orden).
- Pruebas:
  - make build compiló sin errores.
  - Tests unitarios: abortan en test_lb_model.c:259 (pesos del modelo LB, archivo no tocado). No verifiqué si ya fallaban antes; los tests de NP no llegaron a correr.
  - DUMP7-W idéntico bit a bit a DUMP7 (control de regresión de D3Q7).
  - DUMP-W: desviación estándar de δn/n 0.024 (antes 0.246, unas 10 veces menos, compatible con el factor 9).
  - BASE-W: 1000 pasos sin negativos; mínimo de n/n̄ en 0.77 a los 1000 pasos.
- Observación: la desviación iónica sigue siendo ~1.4 veces la del fluido y no está correlacionada con él (0.06). El mínimo sigue bajando lentamente (0.88 en el paso 100, 0.77 en el 1000). No sé si satura.
- Interpretación: el factor 9 era la causa principal. Queda la inconsistencia con el transporte de masa del LB (opción c, pendiente).
- Impacto: toda corrida previa con hydrodynamics on y stencil 19/27 advectaba los iones 5 o 9 veces más rápido que el fluido. Solo afecta corridas con flujo; el equilibrio con el fluido en reposo no cambia.
- Pendiente: corrida larga para ver si el mínimo satura; opción a (Capuani, volume of fluid); opción c.

## 2026-11-12 17:16 — Prueba de 10⁴ pasos con la advección corregida (LONG-W)
- Setup: igual a BASE-W, cambiando solo la salida (estadísticas cada 10 pasos, posiciones cada 100, campos cada 2500). 6 minutos de GPU.
- Observaciones:
  - Sin densidades negativas en 10⁴ pasos.
  - El mínimo de n/n̄ satura en ~0.78 (promedios por bloque 0.787, 0.783, 0.779, 0.775; el más bajo 0.695).
  - Lejos de la partícula (r > 9), la desviación estándar iónica es 0.057 en el paso 2500 y 0.058 en el 10000 (estacionaria); la del fluido es 0.017 y la correlación 0.02–0.04.
  - Cerca de la partícula (r < 3), la nube de contraiones está +22 % (paso 2500) y +30 % (paso 10000) sobre la media, y sigue creciendo.
- Interpretación: el ruido espurio que genera el desacople entre la advección de NP y el transporte de masa del LB llega a un estado estacionario (~6 % rms) en el que la difusión equilibra la generación. Es del orden de una quinta parte de la señal de la nube a esta altura.
- Por qué no se vio antes: en las corridas previas lb_fluctuations estaba apagado (velocidades chicas). Pendiente comparar corridas viejas con flujo usando el ejecutable anterior y el nuevo.

## 2026-11-12 17:41 — Discusión: ¿deben los iones moverse igual que el fluido? Justificación de la opción c
- Objeción del usuario: los iones también se mueven por electromigración y difusión, así que no tienen por qué seguir al fluido. Además, la opción c no está ni en el paper ni en Ludwig.
- Respuesta:
  - Mi afirmación δn/n̄ = δρ/ρ̄ solo vale lejos de las cargas, donde n es uniforme y E≈0. Lo medí ahí (r > 9).
  - Control: NOFLUCT (electromigración y difusión, sin velocidades térmicas) mantiene el campo lejano uniforme a 0.13 % tras 1000 pasos (mínimo global 0.9987 n̄).
  - La restauración electrostática a escala de grilla es lenta: λD ≈ 24 LU con solo contraiones (heurístico) y λD²/D ≈ 8.6e5 pasos.
  - Conclusión: el ~6 % lejano viene de la advección térmica discreta, no de la electrodinámica.
- Opción c: no se implementa sin justificación. Faltan una búsqueda bibliográfica y tests que comparen a, b y c: campo uniforme bajo LB térmico, test de difusión espuria de Capuani (ecs. 30–31), traslación de un perfil, y perfil de equilibrio con y sin fluctuaciones.

## 2026-11-12 17:49 — Advección sola, sin partícula (NOPART-W)
- Setup: LONG-W sin partícula subgrid; sal uniforme n+ = n− = 5.0805e-5 (elección mía, informada al usuario); 2500 pasos; código con la advección corregida. El primer intento abortó porque un comentario superaba el largo de línea que acepta el parser.
- Observaciones:
  - La densidad de carga es exactamente 0 en todos los volcados (las dos especies evolucionan idénticas), así que nunca hay campo eléctrico.
  - Desviación estándar de δn/n̄: 0.054 (paso 500), 0.056 (1000), 0.057 (1500–2500); mínimo −0.24.
  - Fluido 0.017; correlación ion–fluido 0.015–0.030; número total de iones conservado.
- Interpretación: el ~6 % es idéntico al del campo lejano de LONG-W (0.057–0.058) y aparece sin carga ni campo, así que lo produce solo la advección discreta por las velocidades térmicas. La electrodinámica no interviene. En el continuo, un campo uniforme advectado no cambiaría (o seguiría a ρ con ~1.7 % y correlación 1).

## 2026-11-12 17:53 — Explicación: cómo se mide el nivel de ruido y de dónde sale el "ideal"
- Medida: δn_i = n_i/n̄ − 1 en cada nodo de un volcado; σ = desviación estándar sobre los 27³ nodos; r = coeficiente de Pearson entre los mapas δn y δρ. "Estacionario" significa que σ no cambia entre volcados.
- σ_ρ medido (0.0172–0.0174) coincide con la equipartición del LB, √(kT/(ρ c_s²)) = 0.0173 (no lo cito de una referencia verificada).
- "Ideal" (derivación propia, no es un estándar de la literatura): con c = n/ρ, sin carga y con D→0, ∂c/∂t + u·∇c = 0. Si c empieza uniforme, sigue uniforme y δn/n̄ = δρ/ρ̄ (σ = 0.017, r = 1). Con D > 0, la difusión de escala de grilla (~50–150 pasos) es mucho más lenta que el período del sonido (~2–4 pasos), así que n sigue a ρ casi adiabáticamente. Si el fluido fuera incompresible, el ideal sería 0.
- Hipótesis: los iones se advectan con la misma velocidad que la masa del LB; no hay campo eléctrico. La compresibilidad del LB es un artefacto, así que tanto 0 como 0.017 son referencias aceptables; 0.057 sin correlación no es ninguna de las dos.

## 2026-11-12 18:03 — Evaluación de la opción c (advección de iones con el flujo de masa del LB)
- Esquema: el flujo neto de masa por enlace es J_i(x) = f_i*(x) − f_ī*(x+c_i), con f* las poblaciones post-colisión. Flujo iónico = J_i · c_enlace, con c = n/ρ*. Si c es uniforme, n(t+1) = c·ρ(t+1) exactamente, para cualquier interpolación que devuelva c cuando c es uniforme.
- Con poblaciones brutas (Σ f_i* c en lugar de J neto) se agregaría una difusión espuria ~c_s²/2 por paso: descartado.
- Interpolación: centrada (sin difusión numérica, sin garantía de positividad) frente a donante (positiva si f* ≥ 0, pero con difusión numérica ~½|J|, el mismo problema que Capuani). Decisión pendiente del usuario.
- Ludwig: el orden del paso es NP (l.1959) → colisión (2205) → lb_halo (2219) → subgrid/BBL (2241) → propagación (2259). J se calcularía justo antes de la propagación, con halos válidos, y el NP del paso siguiente lo usaría junto con ρ* del mismo paso. LB D3Q19 (18 enlaces, subconjunto de los de NP D3Q27). Las paredes dan J neto = 0 por bounce-back.
- Bibliografía (solo títulos y resúmenes, no leídos): esquemas de escalar pasivo que usan las distribuciones del LB (Onishi et al.; Osmanlic & Körner; comparación de Küng et al. 2018; DeVincentis & Thomas 2023, conservativo con limitador). No encontré nada en electrocinética ni con LB térmico.

## 2026-11-12 18:24 — Lectura de Küng 2020, DeVincentis & Thomas 2023 y Dubois & Lallemand; nuevo análisis de la opción c
- Fuentes (PDF en ~/Descargas/papers a ordenar: main.pdf, main2.pdf, 2306.15291.pdf):
  - Küng, Osmanlic, Markl y Körner, CAMWA 79, 55 (2020), doi 10.1016/j.camwa.2018.01.017.
  - DeVincentis y Thomas, CAMWA 144, 1 (2023), doi 10.1016/j.camwa.2023.05.031.
  - Dubois y Lallemand, PCFD 8, 11 (2008), arXiv 2306.15291.
- Observaciones:
  - Dubois y Lallemand, Prop. 2 (ecs. 14, 17, 19): la actualización de densidad del LB es exactamente un volumen finito con flujo de masa ψ_j = f_j*(x) − f_σ(j)*(x_j), conservativo. Es el J de la opción c.
  - Küng, ecs. 14–15: el método de Osmanlic es la opción c con interpolación donante. Es de primer orden con difusión numérica D = ½u(1−u) (ec. 17, la misma forma que la ec. 31 de Capuani). El método de Onishi (ec. 13) usa poblaciones brutas con anti-difusión C: con C=0 es muy difusivo, con C=1 oscila.
  - DeVincentis y Thomas: FLLA = flujos de la red (ec. 21) + reconstrucción lineal + limitador TVD (ecs. 14–19), de segundo orden.
  - DeVincentis y Thomas §4.2 (flujo turbulento determinista, escalar uniforme): el esquema con velocidades nodales (ec. 20, análogo al de Ludwig) da errores >100 % porque ∇·u_LB ~1e-3 mientras ∂tρ_LB ~1e-6. FLLA no tiene ese error.
  - Ninguno trata difusión física con advección (DeVincentis lo deja como trabajo futuro), ni LB térmico, ni electrocinética.
- Interpretación: la opción c está publicada para escalares pasivos ("lattice advection": Osmanlic, FLLA), y DeVincentis documenta de forma independiente nuestro resultado NOPART. Lo nuevo sería combinarla con la difusión y electromigración de NP y con LB térmico.
- Requisitos de implementación: la reconstrucción con limitador necesita alcance 2 a lo largo del enlace; el halo actual es 4 (por el kernel subgrid).
- Decisión pendiente del usuario: el limitador (donor, LW, minmod, van Leer, superbee).
- Validación propuesta: reproducir D = ½u(1−u) del donor con la caída de amplitud de una gaussiana (Küng ec. 16), y después los tests NOPART, Capuani y nube de equilibrio.

## 2026-11-12 18:41 — Implementación de la advección FLLA (opción c, formulación de DeVincentis, superbee)
- Cambio (CHANGE 20261112):
  - nernst_planck.c: estado estático; kernel np_lattice_mflux_kernel (m = f_ᾱ(j)/ρ_j − f_α(i)/ρ_i, ec. 21, con f post-colisión); kernel np_advective_fluxes_flla_kernel (ecs. 13, 18, 19, superbee); funciones nernst_planck_advection_scheme_set y nernst_planck_lattice_flux_compute; fatal si se elige flla con el driver de CPU.
  - nernst_planck.h: declaraciones.
  - ludwig.c: clave electrokinetics_advection centred|flla (default centred); cálculo de m después de BBL y antes de la propagación (solo lectura del LB).
  - En el código queda anotada la alternativa (flujo de masa sin normalizar y c = n/ρ*) para probar a futuro.
  - Ejecutable anterior (con la corrección de pesos) copiado a scratchpad/Ludwig.exe.advnorm.
- Pruebas:
  - Compila.
  - Regresión: con centred, NOPART-REG es idéntico bit a bit a NOPART-W en el paso 500.
  - NOPART-FLLA (500 pasos, LB térmico, sal uniforme): desviación estándar iónica 0.057 (centred 0.054); correlación con ρ 0.26 (centred 0.03); pendiente 0.86; residuo 0.055; iones conservados; sin negativos (mínimo −0.22).
- Observación: FLLA aumenta la correlación con ρ pero no reduce el ruido.
- Interpretación (hipótesis, no verificada):
  - La normalización por ρ de cada nodo (ec. 21) hace que, aun con c uniforme, el flujo m·q̄ difiera del flujo de masa J·c en ~c·f·δρ (~1e-3 c por enlace y por paso con δρ ~ 1.7 %). Ese error es aleatorio y se acumula hasta que la difusión lo frena (~150 pasos), y da ~5 %.
  - En DeVincentis, ∂ρ/∂t ~ 1e-6 (flujo determinista), así que el efecto allí es despreciable. Con LB térmico no.
  - La correlación baja podría deberse además al desfase de un paso entre m (paso anterior) y el ρ volcado.
- Pendiente: probar la formulación alternativa (J sin normalizar, c = n/ρ*), que es exacta para c uniforme. Requiere autorización del usuario. Faltan también los tests de validación con flujo uniforme (Küng y DeVincentis §4.1).

## 2026-11-12 19:36 — Formulación alternativa flla_mass y test de desfase
- Cambio (nernst_planck.c, nernst_planck.h y comentario en ludwig.c; CHANGE 20261112):
  - Nueva opción electrokinetics_advection flla_mass: J = f_ᾱ(x_j) − f_α(x_i) sin normalizar (Dubois y Lallemand, ec. 17); c = n/ρ* reconstruida con superbee; ρ* guardado en el kernel de flujos y leído con envoltura periódica.
  - Fatal si hay más de un proceso MPI.
  - Los kernels nuevos usan bloques de 4×4×4: con 8×8×8 daban cudaErrorLaunchOutOfResources en la compilación de depuración. Agregué tdpPeekAtLastError después de los lanzamientos.
- Pruebas:
  - Compila.
  - Regresión: centred (NOPART-REG2) y flla (NOPART-FLLA2) idénticos bit a bit a las corridas previas.
  - NOPART-MASS (500 pasos): desviación estándar iónica 0.0173, igual a la del fluido; correlación 1.000; pendiente 1.000; residuo 2e-4; iones conservados; mínimo −0.069 (sin negativos).
  - Desfase (60 pasos, volcados en cada paso, s = 50): centred |r| < 0.07 en todo k ∈ [−2, 2]; flla, r máximo 0.31 (k = −1 y 0); flla_mass, r = 1.000 en k = 0.
- Observaciones:
  - Los volcados de n y ρ del mismo paso ya están alineados en el tiempo (flla_mass correlaciona en k = 0).
  - El desfase no explica la baja correlación de flla ni de centred.
- Interpretación: el ruido de flla se debía a la normalización por ρ de cada nodo (única diferencia con flla_mass), lo que confirma la hipótesis anterior. flla_mass hace que los iones sigan exactamente al fluido en un campo uniforme.
- No verificado: precisión con perfiles no uniformes (tests de flujo uniforme de Küng y DeVincentis), positividad con la nube de una partícula, comportamiento con MPI.

## 2026-11-12 19:45 — Partícula en la celda central con flla_mass (BASE-MASS)
- Setup: idéntico a BASE-W (1000 pasos, misma posición inicial y semilla), con electrokinetics_advection flla_mass y volcados en los pasos 500 y 1000. 1 minuto.
- Observaciones:
  - Sin negativos.
  - n/n̄: mínimo ~0.93–0.94 estable desde el paso 10 (BASE-W bajaba a 0.77); máximo 1.07 (paso 100), 1.14 (500), 1.20 (1000).
  - Campo lejano (r > 9): desviación estándar 0.0173, igual a la del fluido, correlación 1.000, en los pasos 500 y 1000.
  - Nube cerca de la partícula (r < 3): +6.1 % en el paso 500, +9.5 % en el 1000.
  - La trayectoria de la partícula es idéntica a la de BASE-W (rms por paso 0.0364, máximo 0.107, posición final igual).
- Interpretación: el mínimo de ~0.93 es el de las fluctuaciones de densidad del LB (unos 4σ con σ = 0.017). La nube crece de forma monótona. Que la trayectoria no cambie sugiere que la fuerza de la nube sobre la partícula es despreciable frente al ruido en estos 1000 pasos. Los máximos no se pueden comparar directamente con NOFLUCT (allí la partícula no se mueve) ni con BASE-W (incluye ruido espurio).
- MPI (sugerencia del usuario): agregar una capa al halo del LB cuando se usa flla_mass permitiría calcular ρ* en el halo. Alternativa: intercambiar solo el campo ρ*. Pendiente.
- Pendiente: un observable comparable de la nube (perfil radial en el marco de la partícula, promediado en el tiempo) frente a PB; tests de flujo uniforme.

## 2026-11-12 20:08 — Nube de contraiones con flla_mass, con y sin LB térmico (10⁴ pasos)
- Aclaración al usuario: la tabla anterior daba el mínimo y el máximo de n de los contraiones sobre todos los nodos (estadística de Ludwig), normalizados por n̄.
- Setup: LONG-MASS (LB térmico) y LONG-MASS-NOFLUCT (lb_fluctuations off; partícula quieta porque no hay patadas). Volcados cada 500 pasos. Análisis con cloud_profile.py (perfil radial en el marco de la partícula, imagen mínima), salida en cloud_profile_LONG-MASS.txt.
- Observaciones:
  - Exceso medio en r < 3: con fluctuaciones +6 % (paso 500) → meseta de ~+20–27 % desde el paso ~3000 (fluctúa). Sin fluctuaciones +7 % → +69 % a los 10⁴ pasos, todavía creciendo.
  - Perfil promedio entre los pasos 5000 y 10000: r < 1, 1.33 frente a 2.14; r 3–4, 1.13 frente a 1.22; r > 8, ≈ 1.00 en ambos.
  - Contraiones en exceso dentro de r < 6 a los 10⁴ pasos: 0.005 e (con fluctuaciones) y 0.009 e (sin), sobre 1 e en total.
  - La partícula se desplazó ~4.6 LU.
  - Campo lejano (r > 9, con fluctuaciones): correlación con el fluido 1.000 hasta el paso ~4000; baja a 0.91 a los 10⁴ (desviación iónica 0.019 frente a 0.017 del fluido).
  - Sin negativos; carga total conservada.
- Interpretación:
  - Ninguna de las dos corridas está en equilibrio (L²/D ≈ 1e6 pasos), así que no se pueden comparar con PB.
  - La meseta con fluctuaciones sería un balance dinámico entre la formación de la nube y el desplazamiento de la partícula. La partícula se mueve con la velocidad promediada por el kernel Hann-8 más su difusión; los iones, con el flujo local. No está demostrado.
  - La caída de la correlación en el campo lejano se atribuiría a restos de nube que la partícula deja atrás, no a ruido espurio: sin fluctuaciones, el campo lejano también desarrolla estructura (0.0009).
- Pendiente: llegar al equilibrio y compararlo con PB requiere ~1e6 pasos (~10 h de GPU estimadas); no autorizado.

## 2026-11-12 20:16 — Excel de bibliografía reestructurado y análisis de cuánto subir D
- Excel (/mnt/datos/Github/Trabajo/UB/Proyectos/Bibliografia_resumen.xlsx):
  - Respaldo previo: Bibliografia_resumen_respaldo_20261112.xlsx (incluye la edición del usuario de las 20:14).
  - Nuevo orden: Título | Importancia | Resumen | DOI | Año | Autores | 4 columnas de proyecto (TRUE/FALSE, con validación de lista, filtrables).
  - Agregados (verificados en Crossref): DeVincentis & Thomas 2023 (Alta), Küng et al. 2020 (Media; DOI de 2018), Dubois & Lallemand 2008 (Media, PCFD 8, 11).
  - Verificación: 32 filas, 10 columnas, 0 diferencias en los 29 previos frente al respaldo. No lo abrí en Excel ni en LibreOffice.
  - Las casillas nativas no se pueden crear con openpyxl; en Excel 365, Insertar > Casilla las muestra a partir de los valores TRUE/FALSE.
- D (estimaciones heurísticas):
  - Estabilidad difusiva explícita: D ≲ 0.32 LU (×490).
  - Relajación del modo más lento L²/(4π²D) ≈ 2.8e4/k pasos.
  - Péclet iónico térmico u_th/D = 15/k.
  - Deriva por paso D·|∇βeψ| = 0.002·k.
  - Si también se escala la movilidad de la partícula, la patada a 3σ supera drmax = 0.8 para k > 55.
- Interpretación: el equilibrio de la nube no depende de D, pero con LB térmico la forma estacionaria sí, porque cambia el Péclet iónico. Recomiendo ×10 para la forma en equilibrio sin fluctuaciones, sin cambiar la movilidad de la partícula.

## 2026-11-12 20:44 — Forma de la nube con D_ion ×10, flla_mass, con y sin LB térmico (1.5e4 pasos)
- Setup: LONG-MASS-D10 y LONG-MASS-NOFLUCT-D10; d = 1.0614504e-3 (A0 = 6.138, D_física = 6.515e-3 LU); movilidad de la partícula sin cambio. Análisis: cloud_profile.py con WIN0 = 10000, salida en cloud_profile_D10.txt.
- Observaciones:
  - Sin fluctuaciones: el exceso en r < 3 pasa de +1.244 (paso 10000) a +1.258 (15000); casi convergido.
  - Perfil promedio 10000–15000, sin / con fluctuaciones: r < 1, 3.17 / 2.71 (−14 %); 1–2, 2.65 / 2.36 (−11 %); 2–3, 2.06 / 1.93 (−6 %); 3–4, 1.63 / 1.57 (−4 %); 5–6, 1.23 / 1.21; r > 10 igual en ambos.
  - Contraiones en exceso dentro de r < 6: 0.022 e (sin) y ~0.019 e (con), sobre 1 e en total.
  - Con fluctuaciones, la partícula se desplazó ~8 LU y el exceso fluctúa entre +0.87 y +1.12.
  - Campo lejano (r > 9): la correlación ion–fluido baja a ~0.5, pero allí la desviación iónica incluye la cola de la nube (0.026 sin fluctuaciones), así que no es un indicador de ruido espurio en este caso.
- Interpretación:
  - Sin fluctuaciones, la nube llega a un equilibrio débil: no hay condensación para una carga puntual aislada con ℓB = 2.84 y una sola carga en la caja.
  - Con fluctuaciones, la nube cercana es 10–15 % menor. Es compatible con el desfase entre la nube y la partícula que se mueve (efecto de relajación): el tiempo de relajación de la nube a r ~ 3 es r²/D ≈ 1400 pasos, y el desplazamiento relativo en ese tiempo es ~2 LU.
  - No sé si la magnitud es física o numérica (por ejemplo, la difusión del limitador, o que la partícula use la velocidad promediada por el kernel y los iones la local).
  - D_ion/D_partícula = 10 acá; en el modelo de Ingo es 1, así que el efecto en su caso sería distinto.
- Pendiente: referencia PB de la nube de equilibrio (sin fluctuaciones) con la carga Hann-8 en la caja periódica; referencia para el caso con fluctuaciones (iones explícitos o teoría).

## 2026-11-12 21:08 — D_ion ×30 y referencia de equilibrio PB discreta
- Referencia (pb_reference.py, en la carpeta DIAG): PB discreto con la misma discretización que Ludwig (Hann-8 separable con p = 2, laplaciano discreto de 27 puntos con los pesos corregidos el 20251203, φ = βeψ, contraiones canónicos con N = 1, n ∝ e^φ). Iteración de Picard amortiguada: converge en 41 iteraciones a 1e-12.
- Validación: NOFLUCT-D10 en el paso 15000 coincide con PB nodo por nodo a 0.14 % de diferencia máxima; NOFLUCT-D30 en el paso 10000, a < 1e-4. Confirma también el orden de los volcados (x lento, z rápido) y que el estado sin fluctuaciones es el PB discreto exacto.
- Corridas D30: d = 3.1843512e-3 (A0 = 6.138, D_física = 1.95e-2 LU); 1e4 pasos; sin negativos. Factor elegido por mí dentro del rango seguro estimado (×30); informado al usuario.
- Observaciones (perfil promedio frente a PB): con fluctuaciones, r < 1: ×10 −14.7 %, ×30 −6.1 %; r 1–2: −11.3 % y −4.7 %; r 2–3: −6.7 % y −2.5 %; r > 6 dentro de ~0.5 %. El déficit del exceso baja ~2.5 veces al triplicar D.
- Interpretación: el déficit escala aproximadamente como 1/D_ion. Es compatible con dos causas que este test no separa:
  - (a) efecto de relajación por el movimiento relativo partícula–nube;
  - (b) difusión numérica del limitador sin movilidad asociada (temperatura efectiva de los iones mayor).
- Pendiente: LB térmico con la partícula fija para separar (a) de (b).

## 2026-11-12 21:39 — D×30 con LB térmico y la partícula fija (LONG-MASS-D30-FIXED)
- Setup: LONG-MASS-D30 con isfixedr = 1 (un solo cambio). 1e4 pasos, ~7 min.
- Observaciones (perfil promedio 5000–10000 frente a PB):
  - Partícula fija: r < 1, 3.099 (−2.5 %); 1–2, 2.629 (−1.3 %); 2–3, 2.067 (0.0 %); 3–4, 1.640 (+0.5 %).
  - Partícula libre (D30): −6.1 %, −4.7 %, −2.5 %, −1.3 %.
  - Promedio temporal nodo por nodo frente a PB: rms 0.6 % (compatible con el ruido de densidad del LB, 1.7 %/√11), máximo 3.2 %. Instantáneo: rms 1.8 %.
  - La desviación entre volcados en el bin r < 1 es 0.018, así que el −2.5 % parece significativo (correlación entre volcados no evaluada).
- Interpretación (hipótesis):
  - La mayor parte del déficit con la partícula libre se debe a su movimiento relativo a la nube (efecto de relajación): a r < 1, unos −3.6 de los −6.1 puntos.
  - Queda un residuo cerca de la partícula de ~−2.5 %. Candidatos: (i) "difusión de remolino" por advección con las velocidades térmicas, sin ruido iónico que la compense (escala kT/(6πηa) ≈ 5.3e-4 LU, que es el 2.7 % de D_ion×30, el mismo orden); (ii) difusión numérica del limitador. Que (i) sea un artefacto o física consistente depende de que el acople fuerza–fluido compense la movilidad (teoría de Donev); no verificado.
  - Predicción no probada: si es (i), el residuo con la partícula fija escalaría como 1/D (~−7 % con D×10).

## 2026-11-12 22:06 — Partícula fija con D×10, difusión numérica del limitador (escalón) y movimiento de la partícula libre
- LONG-MASS-D10-FIXED (1.5e4 pasos): frente a PB, r < 1 −8.4 %, 1–2 −5.9 %, 2–3 −2.2 %, > 3 dentro de 0.3 %. Con la partícula libre (D10): −14.7 / −11.3 / −6.7 %. El déficit del exceso a r < 1 es 0.123 (D10) frente a 0.037 (D30): cociente 3.4, frente al 3 que predice 1/D. La predicción previa era ~−7 % (observado −8.4 %).
- STEP-ADV (escalón con flujo uniforme u = 0.01, D ≈ 0, 2000 pasos; script analyse_step.py):
  - centred: sobrepaso +27 % y subpaso −82 % (casi negativo) en t = 2000; el frente se atrasa ~0.9 LU (dispersión).
  - flla y flla_mass idénticos (con ρ uniforme coinciden): sin sobrepasos (TVD), frente en la posición correcta (±0.08), σ de 0.75 (t = 250) a 1.01 (t = 2000), D_num ≈ 9e-5 LU (t ≥ 500) y decreciente. Es ~55 veces menor que la celda donante (Küng, ec. 17: 4.95e-3).
- Interpretación del residuo con la partícula fija:
  - La difusión numérica en flujo uniforme (D_num/D = 1.4 % con D10 y 0.5 % con D30) no alcanza a explicarlo.
  - La "difusión de remolino" kT/(6πηa) (De/D = 8.1 % y 2.7 %), traducida a una temperatura efectiva T(1 + δ), predice déficits de ~14 % y ~4.6 % frente a los observados 12.3 % y 3.7 %.
  - Favorece (i), pero el prefactor de De es heurístico (a = 1) y el comportamiento del limitador en flujo térmico (que invierte de signo) no está medido.
- Partícula libre con LB térmico (particle_motion.py):
  - Velocidad media compatible con 0 (|media| ≤ 2 errores estándar por bloques).
  - Autocorrelación de la velocidad 0.04 a lag 1.
  - D_partícula por MSD ≈ 7.0e-4 LU (D_bare 6.5e-4; el objetivo del modelo kT/(6πη a_h) sería 1.7e-3, así que la parte hidrodinámica casi se cancela con al = 0.495).
  - Desplazamiento neto 4.6 LU (1e4 pasos) y 7.7 (1.5e4) frente a √(6Dt) = 6.5 y 7.9.
  - La trayectoria es idéntica con D_ion ×1 y ×30: los iones no afectan a la partícula.
  - Respuesta: en equilibrio ⟨v⟩ = 0, pero la posición no está fija: es una caminata aleatoria con ⟨Δr⟩ = 0 y ⟨Δr²⟩ = 6Dt (solo el promedio de ensamble es nulo).

## 2026-11-12 23:21 — Corrida larga con parámetros de Ingo (LONG-MASS-3E5)
- **Objetivo:** estadística de la nube de contraiones con D_ion de Ingo (×1), LB térmico, partícula libre, advección flla_mass.
- **Cambio:** ninguno en el código. Nueva carpeta microgel/local/INGO-1PART-LBTERM-DIAG/LONG-MASS-3E5 = input de LONG-MASS con N_cycles 3e5 y salida cada 1000 pasos.
- **Prueba:** lanzada en segundo plano; arranca bien (2000 pasos en ~80 s, ~1 GB de GPU).
- **Resultado:** pendiente (~3.3 h).
- **Pendiente:** elegir la ventana de equilibración mirando la serie temporal del exceso (r<3) antes de promediar; comparar con PB y con D×10/D×30.

## 2026-11-13 02:22 — Resultado LONG-MASS-3E5 (D_ion ×1 de Ingo)
- **Prueba:** 3e5 pasos terminados normalmente (2.96 h). Perfil radial con cloud_profile.py → cloud_profile_3E5.txt.
- **Observación:** exceso r<3 por bloques de 2e4 pasos: 0.25 → 0.29 → entre 0.22 y 0.41 desde 4e4, sin tendencia. Perfil 1e5–3e5 (ventana PROVISIONAL): r<1 1.47 (PB 3.18), 1–2 1.40 (2.66), 2–3 1.30 (2.07), 3–4 1.20 (1.63). El exceso a r<1 es ~21% del de PB.
- **Interpretación (no verificada):** coherente con la tendencia D×30 → D×10 (déficit creciente al bajar D) y con la heurística de difusión de remolino (De/D ≈ 0.8 a D×1).
- **Pendiente:** el usuario elige la ventana de equilibración; las fluctuaciones lentas entre bloques indican tiempos de correlación largos (error estadístico aún no estimado).

## 2026-11-13 12:19 — Organización de git: local/ fuera, runs/ dentro
- **Objetivo:** que git contenga solo el código de Ludwig, los espejos <proyecto>/runs/ y las notas (claude/, diario/).
- **Cambio:** (1) runs/ actualizado con util/archive_runs.py (capacitor 11716 archivos, 88 MB; microgel 242793, 3.76 GB). (2) git rm -r --cached de capacitor/local (606) y microgel/local (206): fuera del índice, archivos intactos en disco. (3) .gitignore: las ~30 reglas microgel/local/* + microgel/runs/ + capacitor/ reemplazadas por /*/local/ con comentario explicativo; agregados __pycache__/ y *.pyc. Respaldo del .gitignore anterior en el scratchpad.
- **Prueba:** git check-ignore: local/ ignorado; runs/, claude/, diario/, src/ se suben. Sin commit.
- **Observación:** con las reglas actuales git vería 2.26 GB nuevos, 2.0 GB de ellos en microgel/runs/old; los config.cds.init.001-001 (2230) quedan ignorados por la regla global *001-001.
- **Pendiente:** decidir si microgel/runs/old va a git y si se exceptúan los config.cds.init; commit.

## 2026-11-13 13:58 — Registro de simulaciones en Excel por proyecto
- **Objetivo:** un Excel por proyecto con todas las simulaciones: descripción, parámetros, ID del archivo de coloides, fecha, ubicación (local/runs/borrada), una pestaña por serie científica (agrupación acordada con el usuario).
- **Cambio:** nuevo util/registro_simulaciones.py → capacitor/Registro_simulaciones.xlsx (661 sim., 7 pestañas, 413 archivos de coloides distintos) y microgel/Registro_simulaciones.xlsx (1202 sim., 10 pestañas, 262 coloides; 106 "solo runs", todas de old/).
- **Prueba:** fecha desde Start time del log en 660/661 y 1174/1202 (resto: mtime del input, indicado). Fila LONG-MASS-3E5 y coloide C261 (N=1, Σq0=1, ah=0.3079) coinciden con el input y config.cds.init. Edición manual de una descripción conservada al regenerar (y luego restaurada).
- **Observación:** 2 input de runs/old/ son binarios corruptos (marcados INPUT ILEGIBLE). Un input de 1 línea (NOPART-MASS_fallido) hace que la fila "Comunes" de la pestaña Ingo quede vacía (claves ausentes cuentan como distintas). Las cabeceras de input copiadas de otra corrida dan descripciones heredadas (p.ej. LONG-MASS-3E5 dice "LONG-MASS: ... 10^4 steps").
- **Pendiente:** revisión de las descripciones por el usuario; no abierto en Excel/LibreOffice.

## 2026-11-13 15:38 — Registro de simulaciones: descripciones con la idea de cada serie
- **Objetivo:** quitar las líneas "=====" de las descripciones y que cada descripción diga qué se intenta simular.
- **Cambio:** util/registro_propositos.py (nuevo, editable): IDEA de cada carpeta-serie con su fuente (README, cabecera input, sesiones previas, INFERIDO). util/registro_simulaciones.py: cabecera = solo bloque de comentarios inicial, sin líneas sin letras ni títulos de plantilla; Descripción = "IDEA: … / CASO: …"; CASO = línea del README con el nombre de la corrida > cabecera del input > nombre + archivo de coloides (N, q, posiciones, separaciones, fijas/libres); nueva columna "Fuente de la idea".
- **Prueba:** 0 descripciones con "===" en ambos Excel; muestras revisadas en todas las pestañas (p.ej. SLIDES pos_9 → 2 partículas fijas, r12 = 1.000; LONG-MASS-3E5 → línea del README).
- **Observación:** parte de las ideas son INFERIDAS del nombre/parámetros (marcadas en "Fuente de la idea"); 2 filas de microgel sin documentar (Varios).
- **Pendiente:** revisión del usuario de las ideas INFERIDAS.

## 2026-11-13 16:21 — Registro de simulaciones: columna Conclusión
- **Objetivo:** columna a la derecha de la descripción con qué se puede concluir de cada simulación; llenarla donde haya resultado registrado.
- **Cambio:** util/registro_propositos.py: lista CONCLUSIONES (por carpeta-serie, con fuente: README, diario, memoria, sesiones). util/registro_simulaciones.py: columna "Conclusión" (C), oculta "Conclusión (auto)"; ediciones manuales de Descripción y Conclusión se conservan al regenerar.
- **Prueba:** regenerados ambos Excel; conclusión presente en 588/661 (capacitor) y 473/1202 (microgel) filas; textos cotejados con las fuentes citadas (README WALL-DIRICHLET y F/G, diario 2026-11-05/09/11/12, memoria subgrid_halo y solido).
- **Observación:** sin conclusión registrada: la mayoría de SLIDES sin MEAS, FFT-Force, PETSC-LD8, PIN2GRID, SS-SCAN, RING, MGEL-*-neutro, C-dospotenciales, C-particula-libre, refinamientos -ref.
- **Pendiente:** el usuario completa las vacías; la de INGO-DIAG usa la ventana provisional de LONG-MASS-3E5.

## 2026-11-13 16:26 — Registro de cambios de código (code_changes/)
- **Objetivo:** un archivo que documente los cambios de código con descripción y archivos modificados.
- **Cambio:** carpeta nueva code_changes/: extraer_cambios.py (lee las marcas CHANGE INIT/END de src, target, tests, util, mpi_s, config) → CAMBIOS_CODIGO.md (índice por fecha + detalle con archivo:líneas); NOTAS_MANUALES.md (explicación del porqué de los cambios principales: 2026-11-12 advección NP, 2026-10-28 enlaces ×2, 2026-10-14 σ pm_sr, 2026-09-26 paredes Dirichlet + fuerza ×2, 2025-12-03 laplaciano), copiado al principio.
- **Prueba:** 527 bloques en 54 fechas; entrada 2026-11-12 cotejada con los cambios de esta sesión (ludwig.c, nernst_planck.c/.h); archivos de las notas manuales verificados con grep de las marcas.
- **Observación:** 8 bloques sin CHANGE END (marcados en el listado); cambios sin marca no aparecen.
- **Pendiente:** notas manuales para cambios anteriores a 2025-12 si el usuario las quiere.

## 2026-11-13 16:35 — code_changes: fechas de git y nota de suma compensada
- **Cambio:** extraer_cambios.py fecha las marcas CHANGE INIT sin fecha con el primer commit que introdujo esa línea (git log -S, señaladas "(fecha de git)"); NOTAS_MANUALES.md: entrada "Suma compensada (Kahan/Klein)" (f2mode/mode2f en collision.c: b5bf67b0 2025-10-30 D3Q19, 354f9e7b 2026-01-14 D3Q27 + lb_data/lb_d3q27; Klein en fuerzas subgrid, NP, stats_velocity, psi_fft_pn, ewald).
- **Prueba:** 186 marcas sin fecha fechadas; ninguna queda "sin fecha" ni "sin commitear"; bloques Kahan de collision.c ahora en 2026-01-14 y "Using Kahan version" en 2025-11-19.
- **Observación:** util_sum (kahan_t/klein_t) es del Ludwig original (Stratford 2021). El "por qué" del cambio Kahan no está registrado: queda "por completar".
- 2026-11-13 16:40 code_changes/NOTAS_MANUALES.md: agregado el porqué de la suma de Kahan en f2mode/mode2f (dato del usuario: con 1 partícula cargada + sal, sin campo externo, el redondeo acumulado aceleraba el fluido sin control en vez de conservar el momento). Pendiente: carpeta de la corrida y medida de la mejora.
- 2026-11-13 16:53 code_changes/NOTAS_MANUALES.md: evidencia de la suma de Kahan (microgel/runs/old/sim/particula cargada/*.png y Liquido solo/{Sin_correcion,Correcion_f2m_m2f}) e historia del issue ludwig-cf/ludwig#355. Obs: fluido solo sin corrección v_media ~1.7e-17 tras 3 pasos; con Kahan 0 exacto (inputs idénticos). Pendiente: medida cuantitativa con partícula cargada; la solución no está publicada en el issue.
- 2026-11-13 17:02 code_changes/NOTAS_MANUALES.md: laplaciano D3Q19/D3Q27 publicado en ludwig-cf/ludwig#389 (2026-01-12), incorporado al develop original por PR #411, cerrado 2026-05-07.
- 2026-11-13 17:15 util/registro_simulaciones.py: las corridas registradas cuya carpeta ya no existe en local ni en runs se conservan en el Excel con Ubicación 'BORRADA (no está en local ni en runs)' (valores tomados del registro anterior). Probado en copia simulando la desaparición de LONG-MASS-3E5. Hoy no hay ninguna (las 106 de old conservan su input).

## 2026-11-13 17:24 — git: ni local/ ni runs/; solo el registro y saved_data/
- **Cambio:** .gitignore: /*/local/ y /*/runs/ ignorados; Registro_simulaciones.xlsx y /*/saved_data/** se suben (saved_data exenta de *.csv/*.vtk/*001-001); .~lock.*# ignorado. git rm -r --cached capacitor/runs (438 archivos; quedan en disco). Nuevas capacitor/saved_data y microgel/saved_data con README.md. Cabecera de util/archive_runs.py actualizada.
- **Prueba:** git check-ignore sobre rutas de ejemplo (local, runs ignorados; xlsx, saved_data incl. .csv/.vtk/001-001, src, util, code_changes, diario se suben). Próximo commit: 0.9 MB nuevos + 1251 borrados del índice + 14 modificados.
- **Observación:** el Excel solo no alcanza para rehacer una corrida: no guarda el archivo de coloides (solo su ID/hash), ni el input completo (solo las claves que varían + comunes), ni el ejecutable.
- **Pendiente:** decidir qué pasa a saved_data (condiciones iniciales, README y scripts de WALL-DIRICHLET, evidencia del fix Kahan en old/sim).

## 2026-11-13 17:32 — saved_data: input y coloides base + deltas por corrida
- **Cambio:** util/guardar_bases.py (nuevo): por carpeta-serie guarda en <proyecto>/saved_data/<carpeta>/ input.base (input real de la corrida con menos diferencias), coloides_<k>.cds (contenido más frecuente; otro contenido de igual largo se guarda como cambios de línea si difiere en <= 40 líneas, elección práctica) y deltas.json. util/rehacer_corrida.py (nuevo): reconstruye input + config.cds.init.001-001. Registro: columna "Rehacer desde saved_data". README de saved_data ampliado.
- **Prueba:** reconstrucción de TODAS las corridas: coloides idéntico byte a byte, input con las mismas claves/valores: 0 fallas (661 capacitor, 1200 microgel; 2 inputs corruptos de old/ sin reconstruir). Ejemplo SLIDES-LD8-hann08-corr_yes/pos_12.00: delta = 1 línea (posición), cmp idéntico.
- **Observación:** saved_data 6.4 MB (capacitor; 5.1 MB es el sólido L48) + 4.2 MB (microgel). No incluye auxiliares del input ni ejecutable.

## 2026-11-13 17:40 — Evidencia Kahan en microgel/Graficos + documento de scripts
- **Cambio:** code_changes/NOTAS_MANUALES.md: evidencia de microgel/Graficos (en git, commit e93f2fb4 2025-11-10): Fluid/*-velocidades.png (sin corrección, v media ~±1e-17) vs *-momento_kahan*.png (0 exacto); Sin campo/: partícula cargada sin campo, v ~1e-16 y fluido ~1e-17 sin crecer (que sean con corrección se infiere de la fecha). util/documentar_scripts.py (nuevo) → SCRIPTS.md en la raíz (qué hace y cómo se usa cada script, desde su docstring/comentarios).
- **Prueba:** SCRIPTS.md generado: 10 herramientas generales (1 sin doc: util/extract.sh), 19 programas C de util, 64 scripts de capacitor (0 sin doc), microgel ~145 (sin doc: calculos*.py, plot*.py, copy_plots, extraer_posicion, lattice_error_*, 3 run.sh, freepos.py, aaa-otro/*). Duplicado de microgel/local/scripts corregido.
- **Observación:** los scripts de <proyecto>/local/ ya no van a git (local/ ignorado); SCRIPTS.md los nombra pero su código no queda versionado.
- **Pendiente:** decidir dónde versionar los scripts de proyecto; documentar los que no tienen docstring.

## 2026-11-13 17:45 — Scripts de proyecto movidos a <proyecto>/scripts/ (opción a)
- **Cambio:** 98 movimientos (218 archivos, 1.8 MB): carpetas local/scripts y local/<serie>/scripts → <proyecto>/scripts[/<serie>] con enlace simbólico en el lugar original; scripts sueltos de cada serie y de la raíz de local → <proyecto>/scripts/<serie>/ con enlace duro (22 lanzadores usan readlink -f "$0"). Lista en code_changes/scripts_movidos_20261113.json. README_ubicacion.md en cada scripts/. util/documentar_scripts.py lee <proyecto>/scripts/; SCRIPTS.md regenerado (243 scripts).
- **Prueba:** sin colisiones de nombres; readlink -f de lanzadores sigue dando su carpeta de serie; analizar_C.py corrido por la ruta vieja reproduce los números del README (campo 1.5e-2 máx, mediana 3.48e-4; carga 1.23e-14).
- **Error mío:** esa prueba regeneró graficos/C_carga_imagen_hann08.png (el script guarda la figura sin avisar). Resultado idéntico byte a byte a la copia de runs/ (mismo md5); fecha original restaurada con touch -r. No se perdió nada.
- **Pendiente:** commit; scripts sin docstring (30 en microgel, util/extract.sh).
- 2026-11-13 17:57 claude/petscrc_modifications.md: referencia del .petscrc cambiada a la raíz (./.petscrc, idéntico a microgel/local/.petscrc y a las 792 copias por corrida).

## 2026-11-13 18:59 — H-microgel-ruido: microgel de G con ruido térmico (prueba corta)
- **Objetivo:** microgel entre placas con ruido térmico; evaluar si el ruido es demasiado fuerte, con la escala de Ingo y el mismo problema físico que G.
- **Cambio:** capacitor/local/WALL-DIRICHLET/H-microgel-ruido/{neutro,superficie}-{steeleG,steele1e-5} (crear.py, correr.sh, analizar_ruido.py en capacitor/scripts/H-microgel-ruido). Escala C: kT 1e-4, η 0.01, ρ 1, a_h 0.308 (al 0.474 de G); ε 1e3, k 1e-3, ε_rep 1e-4 (lB/r0, √(kT/k)/r0, ε_rep/kT iguales a G). Sin sal. 5000 pasos.
- **Prueba:** las 4 terminaron normalmente (0.15–0.18 s/paso, 4 en paralelo).
- **Observación:** desplazamiento rms 0.0344/paso/componente (predicho √(2D_desnudo) = 0.0348), máx 0.151 = 19 % de drmax; ⟨U_enlaces⟩/kT = 17.3/17.0 ± 0.5 vs 17.5 (equipartición, aprox. armónica); fluido ⟨u_a²⟩/(kT/ρ) = 0.98–1.02; d mínimo a pared 8.8 → las dos paredes dan corridas idénticas byte a byte; Rg 4.64 → 5.07 (neutro) / 5.19 (superficie), aún relajando.
- **Interpretación:** con la escala C el ruido no es demasiado fuerte para los monómeros ni para el fluido. La estabilidad de la pared con ruido no quedó probada (ningún monómero llegó a ella).
- **Pendiente:** prueba de pared con un monómero libre cerca del plano; duración de la corrida larga.

## 2026-11-13 19:12 — H-microgel-ruido: prueba de pared con ruido (un monómero)
- **Cambio:** crear_pared.py / correr_pared.sh / analizar_pared.py (capacitor/scripts/H-microgel-ruido). Un monómero libre q = 0 y q = 1, paredes Steele ε_w = 6.2572e-19 (G ×10) y 1e-5 (σ_w = 1), 2000 pasos.
- **Prueba:** primer intento desde d = 4.0: Steele G aborta en el paso 21 (h llegó a 0.019, dr = 696/739); Steele 1e-5 aborta en el paso 1 por error mío (arranque en h = 0.5, dentro de su zona repulsiva) → renombradas _fallido / -d4. Repetidas desde d = 5.0 (h = 1.5).
- **Observación:** Steele G desde d = 5: q0 aborta en el paso 877 (dr = 1.3e4) tras pasar 772 pasos dentro del alcance con ⟨V⟩/kT = 41 (Boltzmann 0.005); q1 termina (h mín 0.040, máx dr 30 % drmax). Steele 1e-5: ambas terminan, h mín 0.768, V máx 2 kT, máx dr 17 % drmax, ⟨V⟩/kT 0.25 vs Boltzmann 0.30 (300 muestras). q1 ≈ q0 (imagen ~0.04 kT).
- **Interpretación:** confirma el criterio μ·V'' < 1: la Steele de G no sirve con ruido (aborta o muestrea mal); ε_w = 1e-5 sí.
- **Error mío:** posición inicial dentro de la zona repulsiva; sed -i rompió el enlace duro de correr_pared.sh (rehecho).
- **Pendiente:** corrida larga con ε_w = 1e-5.
- 2026-11-13 19:19 H-microgel-ruido: lanzada la corrida larga {neutro,superficie}-largo (crear_largo.py, correr_largo.sh): pared Steele ε_w = 1e-5, 50000 pasos, posiciones cada 100, volcados completos cada 10000 (reinicio). Misma configuración inicial y semilla que las pruebas cortas.
- 2026-11-13 19:21 H-microgel-ruido corrida larga: primer lanzamiento falló (neutro-largo abortó: comentario > 127 bytes partido por el lector de input, NKEY_LENGTH en src/runtime.c; error mío). superficie-largo detenida y sus 200 pasos parciales borrados; comentarios partidos en ambos inputs y en crear.py/crear_largo.py; relanzadas las dos, sin avisos de claves no usadas. Las pruebas cortas tenían la clave basura 'egaban' sin efecto.
- 2026-11-13 19:41 Memoria: idea de potencial de pared ψ(t) (idea_wall_potential_time_dependent.md). Nuevo PENDIENTES.xlsx en la raíz: 27 tareas ordenadas por prioridad y dificultad (valores propuestos por Claude, a revisar), con leyenda de criterios.
- 2026-11-13 20:06 Interpretación de la rampa de pared con e0 tangencial corregida (no es conductor ideal: electrodo resistivo / arreglo de electrodos; imagen de conductor; no modela paredes aislantes) en src/psi_petsc.c (solo comentario, no recompilado), claude/PAREDES_POTENCIAL_CONSTANTE.md §8 y README de F-campo-tangencial; nota en code_changes.
- 2026-11-13 21:15 H-microgel-ruido corrida larga: análisis con el script nuevo analizar_largo.py (capacitor/scripts/H-microgel-ruido/, hardlink en la carpeta de la serie; salida en logs_analisis_largo.txt; SCRIPTS.md regenerado).
  - **Observación:** las dos terminaron normalmente. Ningún monómero llegó a la pared (h ≥ 4.6). ⟨U_enlaces⟩/kT por bloque de 5000 pasos: 16.1–19.3 (ref. 17.5). Fluido: ⟨u²⟩/(kT/ρ) 0.99–1.02. Rg: 4.64 → ~5.1 en 10000 pasos, después oscila entre 4.2 y 5.0 sin tendencia. Neutro y superficie casi iguales bloque a bloque (misma semilla).
  - **Interpretación:** Rg aparentemente estacionario desde ~1e4 pasos; las fluctuaciones por bloque son mayores que la diferencia neutro–superficie, así que el efecto de la carga no es medible con esta estadística. Las corridas no son independientes.
  - **Pendiente:** ventana de equilibrio (decisión del usuario); ¿más pasos o semillas distintas? ¿cuánto conviene acercar el gel a la pared?
- 2026-11-13 21:40 Pregunta: ¿se puede incluir el ruido en PB para comparar con LONG-MASS-3E5? (respuesta conceptual, nada implementado)
  - **Observación (parámetros del input):** D_ion = D_desnudo = 6.515e-4, D_Stokes = 1.72e-3, con advección flla_mass. Medido: n/n̄(r<1) = 1.47; PB da 3.18.
  - **Interpretación / derivación propia, sin fuente verificada:**
    - Con NP determinista y la partícula browniana, en DH lineal la nube vista desde la partícula es DH con κ_ef² = κ²/(1 + D_p/D), es decir T_ef = T(1 + D_p/D).
    - En equilibrio consistente (iones con ruido y su correlación con la partícula) seguiría siendo PB: la diferencia viene de que NP es campo medio sin ruido.
    - Estimación grosera no lineal: 3.18^(T/T_ef) da 1.8 (D_p = D_desnudo) a 1.4 (D_p = D_Stokes), compatible en orden con 1.47. Es solo un indicador.
  - **Pendiente:** decidir si comparar con la fórmula lineal y qué D_p usar (difusión relativa partícula–iones, afectada por la advección correlacionada); buscar fuente (Démery & Dean, sin verificar).
- 2026-11-14 Pregunta conceptual (sin cambios): ¿la advección de NP da ruido a los iones? ¿PB depende de D?
  - **Respuesta:** la advección por el u fluctuante es ruido multiplicativo u·∇c. No genera fluctuaciones donde c es uniforme y no reemplaza al término de ruido propio ∇·(√(2Dc) W) de la hidrodinámica fluctuante (Donev et al. 2014; Péraud et al. 2016). La PB de equilibrio no depende de D, solo de kT, z y ε; D fija solo cuánto tarda la nube en seguir a la partícula.
  - **Propuesta (no lanzada):** variar D_ion con todo lo demás fijo. Predicción lineal: κ_ef² = κ²/(1 + D_p/D), así que con D mayor la nube debería acercarse a PB.
  - **Pendiente:** autorización, factores de D y duración.
- 2026-11-14 Pregunta: ¿qué es el ruido desnudo de NP y cómo se comporta el modelo frente a la teoría? (sin cambios; solo cuentas sobre perfiles ya existentes de INGO-1PART-LBTERM-DIAG)
  - **Observación (README de la serie):**
    - Sin fluctuaciones, el perfil coincide con PB discreta (D×10: 0.14 %; D×30: < 1e-4).
    - Con LB térmico y la partícula fija, n/n̄(r<1) da −8.4 % (D×10) y −2.5 % (D×30).
    - Con la partícula libre da −14.7 % (D×10), −6.1 % (D×30) y −54 % (D×1).
  - **Interpretación heurística:**
    - Ajusté una T_ef por corrida con n = n_PB^(T/T_ef), solo en el bin r < 1 (regla inventada, no derivada). Despejando D_extra = (T_ef/T − 1)·D:
      - partícula fija: ~5e-4 (lo atribuyo a la agitación advectiva, D_adv);
      - partícula libre: ~1.1e-3, que es D_adv más el D de la partícula (6.5–7e-4);
      - D×1: 1.3e-3.
    - Esto es consistente con T_ef = T(1 + (D_p + D_adv)/D). Pero la regla de exponente falla en r = 3–4 (el sistema es sin sal y la carga se conserva).
  - **Pendiente:** decidir si se ajusta T_ef con pb_reference.py usando lB_ef = lB·T/T_ef (solo Python, barato). Preguntarle a Ingo qué significa "ruido EK ≈ 0".
- 2026-11-16 Documento claude/MODELO_EK_PETSC_PASO_A_PASO.md: análisis paso a paso, leyendo el código, de cómo Ludwig resuelve partículas subgrid cargadas + NP + Poisson-PETSc + LB térmico. Incluye ecuaciones en LaTeX, el orden del paso de tiempo, y supuestos y simplificaciones (tabla de 18). No modifica código ni corre nada.
  - **Observaciones de la lectura (no verificadas en corridas):**
    - la velocidad de la partícula y los enlaces usan siempre Peskin-4, mientras que la carga y la fuerza eléctrica usan el kernel elegido;
    - el ruido de la partícula está truncado en 3σ (varianza × 0.973);
    - el ruido del LB no tiene ρ local;
    - drmax se compara sin valor absoluto;
    - en GPU, maxacc = 0, así que los subpasos de NP no se adaptan;
    - psi_zero_mean también se aplica con paredes;
    - la corrección de fuerza (2−ω) del LB no está contrastada con Guo.
  - **Pendiente:** que el usuario lo revise; analizar la rama FFT.
- 2026-11-16 Pregunta: ¿el arreglo de la advección (2026-11-12) cambia los gráficos de equilibrio químico (Δμ)? (respuesta sin correr nada)
  - **Observación:**
    - TEST-HIDRO-dx2 (hidro_si y barrido de η: graficos/dmu_vs_viscosidad.png, dmu_u_vs_viscosidad.png) y las Δμ de FUERZA-2PART-MEAS se hicieron el 10–11/11 con los exe 4f0c20e7/45dc94bc: advección centred sin peso, stencil 27, así que los iones se advectaban 9× más rápido.
    - Las corridas hidro_no no tienen advección: no cambian.
  - **Interpretación (heurística):**
    - El piso de Δμ con hidro debería bajar ~9× si es lineal en el flujo advectivo (el barrido en η lo sugiere).
    - F y la fuerza no deberían cambiar de forma visible (el efecto del flujo era < 1e-6 relativo).
    - La conclusión "sin flujo, μ uniforme a precisión de máquina" no cambia.
  - **Pendiente:** autorización para rehacer hidro_si (y el barrido) con el exe actual, en carpetas nuevas (~5 min cada una).
- 2026-11-16 00:45 TEST-HIDRO-dx2: barrido de viscosidad rehecho con el exe actual (221f3229, con el arreglo del peso de advección).
  - **Cambio:** carpetas advfix/hidro_si* creadas por scripts/crear_advfix.py (copias del input y del config.cds.init; solo cambia el exe). Corridas con correr_advfix.sh: 7/7 terminaron normalmente, ~6 min con 7 en paralelo. Gráfico graficos/dmu_u_vs_viscosidad_advfix.png (grafico_dmu_u_visc_advfix.py); el original no se tocó.
  - Error mío: el primer intento usó scripts/.. (symlink) y creó carpetas vacías en microgel/scripts/advfix; las borré con rmdir y ahora el script usa rutas absolutas.
  - **Observación (paso 3000):**
    - max|u| y RMS|u| idénticos antes y después (4 cifras).
    - Δμ₋ baja 10.1× (η = 0.02) a 8.4× (η = 2).
    - Δμ₊ baja 6.8× a 5.5×.
    - Después del arreglo Δμ₊ ≈ Δμ₋; antes, Δμ₋ ≈ 1.5 Δμ₊.
  - **Interpretación (no verificada):**
    - Se confirma que el piso viene de la advección por el flujo estacionario, y que el arreglo no cambia el flujo.
    - El factor no es 9 exacto para las dos especies; no lo explico.
    - En η ≥ 1 puede pesar el transitorio que todavía queda en el paso 3000 (la corrida sin hidro tenía 1.8e-8 a 3000).
  - **Pendiente:** si hace falta, extender a 6000 pasos para separar el piso del transitorio en η alto.
- 2026-11-16 01:05 TEST-HIDRO-dx2: el mismo barrido de viscosidad con advección flla (pedido del usuario).
  - **Cambio:** carpetas advfix_flla/hidro_si* creadas por crear_advfix_flla.py (exe 221f3229 + "electrokinetics_advection flla"). Corridas con correr_advfix_flla.sh: 7/7 terminaron normalmente. Gráfico de las tres variantes en graficos/dmu_u_vs_viscosidad_3adv.png (grafico_dmu_u_visc_3adv.py).
  - **Observación (paso 3000):**
    - |u| idéntico en las tres variantes.
    - Con flla, Δμ₋ = 1.5e-6 (η = 0.02) baja hasta ~3.8e-7 y se planta desde η ≈ 0.5.
    - Con flla, Δμ₊ tiene un mínimo de 1.1e-7 en η = 0.1 y sube a 2.9e-7 en η = 2.
    - flla queda por encima de centred en casi todo el rango (×1.6 a ×12), salvo Δμ₊ en η = 0.05 y 0.1.
  - **Interpretación (hipótesis, no verificada):**
    - El piso de flla no baja con η, así que no lo pone u.
    - Candidato: flla toma el flujo de masa post-colisión, Σ f* c = ρu + F/2, que en estado estacionario incluye la media fuerza F/(2ρ), independiente de η. centred usa u = (Σ f c + F/2)/ρ.
    - La corrida de η = 0.05 tenía max|u| = 4.4e-8 en el paso 1, independiente de η, que es del orden de ese término.
  - **Pendiente:** confirmar con el flujo m guardado frente a u; ver si flla_mass hace lo mismo.
- 2026-11-16 01:30 Comprobación de la hipótesis "flla advecta con u + F/(2ρ)" (pedido del usuario). Solo se leyeron archivos, sin corridas ni cambios de código.
  - **Método:** f* post-colisión del paso 3000, reconstruida exactamente desde dist-3000 (f*_a(x) = f_dump_a(x + c_a)). Scripts comprobar_flla_fuerza.py y comprobar_flla_descomposicion.py; salidas en advfix_flla/comprobar_flla_*_3000.txt.
  - **Observación:**
    - F/(2ρ) es independiente de η: RMS 8.5e-10, máximo 3.7e-8.
    - La fuente de flla para un campo uniforme, S = Σ_a m_a, tiene RMS de 5.2e-10 a 6.5e-10, casi independiente de η.
    - Descomposición exacta S = D/ρ + R: el flujo neto de masa D del LB es ~1e-16 (1e-13 en η = 0.02), y R, el término de diferencias de densidad del LB entre vecinos, da corr(S, R) = 1.0000.
    - Las variaciones de ρ_LB son independientes de η: RMS 1.6e-9, máximo 3.3e-8 a 5.2e-8.
    - S se correlaciona con −div(F/2): 0.90 en η = 0.02 y 0.84 en η = 2.
  - **Interpretación:**
    - La hipótesis F/(2ρ) queda REFUTADA como mecanismo, porque el flujo neto de masa es nulo.
    - La fuente viene de normalizar por la densidad local del LB en flla (modo 1): m = f̄/ρ_j − f/ρ_i. Es la limitación que ya documenta el código (flla conserva un campo uniforme solo hasta las diferencias de densidad del LB).
    - Esas diferencias de densidad equilibran la fuerza eléctrica (c_s²∇ρ ≈ F, presión), por eso no dependen de η, y R se correlaciona con div F. Esto último es inferencia: no comparé ∇ρ con F directamente.
    - flla_mass conserva c uniforme por construcción, así que no debería tener este piso. No verificado.
  - Corrijo la entrada de 01:05: el término F/(2ρ) no es la causa.
  - **Pendiente:** el barrido con flla_mass, si el usuario lo autoriza.
- 2026-11-16 01:50 TEST-HIDRO-dx2: el mismo barrido de viscosidad con advección flla_mass (pedido del usuario).
  - **Cambio:** carpetas advfix_flla_mass/hidro_si* creadas por crear_advfix_flla_mass.py y corridas con correr_advfix_flla_mass.sh: 7/7 terminaron normalmente. Gráfico de las cuatro variantes en graficos/dmu_u_vs_viscosidad_4adv.png (grafico_dmu_u_visc_4adv.py); tabla en advfix_flla_mass/tabla_4adv_3000.txt.
  - **Observación (paso 3000):**
    - flla_mass baja con η como centred, sin el piso de flla.
    - Queda ~1.4× por encima de centred (Δμ₊ ×1.38 a ×1.18; Δμ₋ ×1.41 a ×1.44).
    - |u| es idéntico en las 4 variantes.
  - **Interpretación:** es coherente con la explicación del piso de flla por la normalización con la densidad del LB, que flla_mass no tiene. No explico el factor ~1.4 respecto de centred; candidatos: difusión numérica del limitador, o el desfase de un paso en el flujo.
  - **Pendiente:** ninguna variante llega a μ uniforme con flujo; queda por evaluar si el paso 3000 ya es estacionario en η alto.
- 2026-11-16 02:10 TEST-HIDRO-dx2: continuación a 6000 pasos de las 28 corridas (4 variantes × 7 η) y curva de Δμ en estado estacionario (pedido del usuario).
  - **Cambio:**
    - Carpetas cont6000/<variante>/<corrida> creadas por scripts/crear_cont_visc.py (receta de crear_continuacion_dxhalf.py; mismo exe que la corrida de origen). Corridas con correr_cont_visc.sh 6000: 28/28 terminaron normalmente con "Re-starting simulation at step 3000".
    - En flla y flla_mass, el primer paso después del reinicio no tiene advección (no hay todavía flujo guardado).
  - **Criterio de estacionario** (heurístico, propuesto por mí): |Δμ(N) − Δμ(N−1000)|/Δμ(N) < 1 % para las dos especies.
  - **Observación:**
    - Todas cumplen con mucho margen: r(1000) va de 1e-8 a 1.7e-6 y r(1500) ≈ 8 × r(1000), así que sigue decayendo. Los archivos no son idénticos byte a byte.
    - Entre 3000 y 6000 casi nada cambió (< 0.5 %), salvo Δμ₊ de centred en η = 2: 3.43e-8 → 2.91e-8 (−15 %).
    - Salidas: cont6000/estacionario_6000.txt y estacionario_6000_detalle.txt; tabla en cont6000/tabla_4adv_6000.txt.
    - No hizo falta llegar a 9000.
  - **Gráfico:** graficos/dmu_u_vs_viscosidad_4adv_estacionario.png (grafico_dmu_u_visc_4adv_estacionario.py).
    - centred tiene el menor Δμ salvo Δμ₊ de flla en η = 0.05 y 0.1.
    - En estado estacionario, flla_mass queda ×1.39–1.44 sobre centred.
    - Al piso de flla le corresponde centred/flla = 0.08 en η = 2.
  - **Interpretación:** las conclusiones del paso 3000 se mantienen.

## 2026-11-16 — Excel de conversión de unidades Ingo → Ludwig
- Objetivo: planilla que haga la conversión de los parámetros de Ingo (Tabla 1, ec. 2 y 5 de su PDF) a LU, con explicación.
- Cambio: nuevo microgel/conversion_unidades_Ingo_Ludwig.xlsx (hojas Conversion y Explicacion; entradas en azul, fórmulas).
  Generado con un script en el scratchpad (no forma parte del proyecto).
- Prueba: LibreOffice no está instalado, así que no se recalculó con una planilla; las fórmulas se evaluaron con un
  intérprete propio en Python. Los 7 chequeos contra el input de INGO-1PART-LBTERM dan TRUE.
- Observación: el input usó μ_b = D_LU/kT_LU = 6.515; 1/γ_LU con γ = 15.35 da 6.51466 (diferencia relativa 5e-5,
  γ = 1/0.06515 = 15.349 redondeado) y daría a_h = 0.3078942 en vez de 0.3078881. La planilla usa D/kT y muestra 1/γ.
- Pendiente: abrir el archivo en Excel para confirmar que calcula igual que la evaluación en Python.

## 2026-11-16 12:59 — Lista de tareas en Excel (Proyectos/Tareas.xlsx)
- Objetivo: pasar el archivo UB/Tareas a un Excel con tareas hechas y pendientes, una pestaña por tipo/proyecto.
- Cambio: nuevo /mnt/datos/Github/Trabajo/UB/Proyectos/Tareas.xlsx (pestañas Administrativo, Común (Ludwig),
  Microgel (tesis), Stuttgart (Ingo), EDLC (Benjamin), Leyenda). Fuentes: archivo Tareas, PENDIENTES.xlsx y memorias
  de Claude. No se modificó ninguna fuente. Script en el scratchpad (no forma parte del proyecto).
- Prueba: se releyó el archivo con openpyxl (77 filas, 29 hechas); no se abrió en Excel/LibreOffice.
- Observación: el arreglo de la fuga de klein_t ya está en src/subgrid.c:1376 (la memoria decía "no aplicado").
- Interpretación: "Convertir unidades" de Ingo se marcó Hecho por inferencia (ya hay corridas con sus parámetros).
- Pendiente: revisar la asignación de tareas a pestañas y las prioridades vacías.

## 2026-11-16 14:41 — Tareas_v2.xlsx: subítems, color por tarea y formato condicional en Estado
- Objetivo: tareas con subítems y color distinto por tarea; color de la columna Estado según su valor.
- Cambio: nuevo Proyectos/Tareas_v2.xlsx (columnas Sub y Avance, colores fijos por tarea, sin formato condicional
  en las filas). Luego, formato condicional solo en la columna D (Pendiente naranja, Hecho verde, En curso amarillo,
  Descartado tachado gris), escrito en el XML en formato x14 como OnlyOffice. Respaldo previo: Tareas_v2_respaldo_1441.xlsx.
- Prueba: zip válido, XML bien formado, 4 reglas y 2 listas desplegables por pestaña, valores idénticos al respaldo.
- Observación: el formato condicional anterior (filas enteras) pisaba los colores manuales; se descartó.
- Pendiente: ver en OnlyOffice que se vean los colores y que la fórmula de Avance calcule (no evaluada).

## 2026-11-16 — Contraste de la interpretación de parámetros ESPResSo con Cimrák et al. 2014
- Objetivo: verificar con Cimrák et al. 2014 (Comput. Phys. Commun.; OIF en ESPResSo) la conversión Ingo -> Ludwig.
- Observación: el paper define el acoplamiento F = ξ(v − u) con u interpolada linealmente (trilineal) y la fuerza
  devuelta al fluido con los mismos pesos (ec. 3); ξ es el "friction" de lbfluid. Las partículas tienen masa y siguen
  m d²X/dt² = F (ec. 2, masa por defecto 1). Usa D3Q19. No define si "visc" es cinemática ni el tau del LB de Ingo.
- Interpretación: γ de Ingo = ξ desnudo -> μ_b = 1/γ en Ludwig es consistente. Falta en el mapeo: la inercia de la
  partícula de ESPResSo (m/γ = 0.065 [t] = 6.5 δt si m = 1), ausente en la partícula subgrid de Ludwig.
- Pendiente: confirmar con Ingo masa de los monómeros, tau y agrid del LB, versión de ESPResSo y si ν es cinemática.

## 2026-11-16 — Parámetros de Ingo contra la documentación de ESPResSo 5.0.1 (lb, ek, integration, API)
- Objetivo: revisar la conversión Ingo -> Ludwig suponiendo ESPResSo 5.0.1.
- Observación (doc): LBFluid(agrid, tau, density, kinematic_viscosity, kT) en unidades MD; viscosidad CINEMÁTICA.
  Termostato LB: m dv/dt = f − γ(v − u) + sqrt(2γkT) η, u interpolada lineal en 8 nodos (Ahlrichs-Dünweg).
  LB y MD no sincronizados (tau_LB puede ser varios δt); tau de EK = múltiplo entero de δt. Colisión TRT.
  EK: j = −D∇n − (D/kT) q n ∇Φ + n v_fl + sqrt(n) W; ΔΦ = −4π lB kT Σ q n; EKFFT(permittivity = ε0 εr);
  EKSpecies(diffusion, valency, advection, friction_coupling, kT, thermalized). Esquema de advección y
  acople carga-partícula -> grilla EK no documentados en esas páginas.
- Interpretación: k_BT, ρ, ν cinemático, D, ε = 1/(4π lB kT) (= 280.2 en LU) y μ_b = 1/γ son consistentes.
  Diferencias no incluidas en el mapeo: inercia (masa) de la partícula; colisión TRT vs MRT M10.
- Pendiente (preguntas a Ingo): tau del LB y de EK, agrid de EK vs LB, friction_coupling, thermalized,
  esquema de advección, cómo acopla la carga de los monómeros a EK (CIC/Hann no aparecen en la doc).

## 2026-11-16 — Justificación de los pesos advectivos de NP
- Pregunta: ¿el LB usa esos pesos?, ¿qué sustenta W_p = 2 wgradients[p] en la advección centrada de NP?
- Observación: wgradients[p] = 3 w_p = w_p/c_s² (stencil_d3q19.c:70, stencil_d3q27.c:73; D3Q7: 1/2).
  Ramadugu, Thampi, Adhikari, Succi, Ansumali, arXiv:1208.1009, ec. (10): ∇·φ = (1/T) Σ_i w_i c_i·φ(r+c_i) + O(∇³),
  T = 1/3 (= c_s²). Verificado en el texto del preprint; revista no verificada.
- Interpretación: el flujo de enlace de Ludwig con W_p, sumado sobre p, es esa ec. (10) aplicada a ρu (salvo el
  producto de promedios, O(c²)). El paper no trata advección de NP ni estabilidad; la aplicación es nuestra.
  El LB no tiene un término advectivo con pesos: los w_i entran por f_eq.

## 2026-11-16 — Reescritura de §7.2 (Advección) de claude/MODELO_EK_PETSC_PASO_A_PASO.md
- Pedido: mejorar la descripción; explicar flla_mass, por qué centred usa pesos y qué es I.
- Cambio: §7.2 reescrita en 7.2.1 centred (pesos W_p = 2 wgradients = 2 w_p/c_s², derivación por Taylor,
  tabla 1/5/9, I = identidad 3x3, ref. Ramadugu et al. ec. 10), 7.2.2 flla, 7.2.3 flla_mass
  (rho(t+1) = rho* + Σ J exacto, propiedad c uniforme), 7.2.4 tabla resumen. Basado en lectura de
  src/nernst_planck.c (l. 1246, 1461-1720, 1850-1875). Otras secciones sin tocar.
- Observaciones de la lectura: upwind en flla/flla_mass se elige por el signo de u, no de J; con stencil 27
  los 8 enlaces de esquina no llevan advección en flla*; el comentario CHANGE de l. 1123 llama M a Σc_x²
  con valores 1/5/9, que en realidad son ½Σc_x² (sin efecto en el código).
- Pendiente: §7.4 punto 3 dice "solo flla_mass tiene limitador" (flla también); no corregido.

## 2026-11-16 — Advección centred de NP con paredes (lectura de código + una corrida existente)
- Pregunta: ¿la advección centred omite los nodos de pared?
- Observación (código): en el driver de GPU, np_advective_fluxes_kernel no mira el mapa; la difusión sí omite
  enlaces a nodos no fluidos; np_antisym no mira el mapa; np_update_rho solo actualiza nodos fluidos;
  np_no_flux_boundary (pone a 0 los flujos fluido-sólido) existe solo en el driver de CPU.
- Observación (C-electrolito-hann4-LD3-ref0.5-desplazada/z_12.50, paso 2000, sin ruido): u = 0 en los planos de
  pared, n en la pared = valor inicial 1.111e-2 (no se actualiza); |u| en el fluido <= 2e-9; iones totales en el
  fluido conservados a ~1e-10 relativo.
- Interpretación: con centred, por un enlace fluido-pared pasa W_p ½(u_i·c_p) ½(n_i + n_pared) que sale del fluido
  y no entra a ningún lado (no conserva). Proporcional a la velocidad normal en la primera capa de fluido:
  despreciable en flujo paralelo sin ruido; relevante con LB térmico + paredes + iones. No medido en ese caso.
- Pendiente (no implementado): enmascarar los enlaces no fluidos en el kernel centred de GPU.

## 2026-11-16 — No-flujo en el driver de GPU de NP (CHANGE 20261116) + prueba TEST-NOFLUX-ADV
- Objetivo: que la advección centred no pierda iones por los enlaces fluido-sólido (pedido del usuario).
- Cambio: src/nernst_planck.c, nuevo np_no_flux_boundary_kernel (anula flujos de enlace con un extremo no
  fluido), llamado entre los kernels difusivo y de antisimetrización. Exe 221f3229 -> 269ec1b8.
  Registro: NOTAS_MANUALES.md + CAMBIOS_CODIGO.md. Documento del modelo §7.2.1, §7.2.4 y §7.4.5 actualizados.
- Prueba: capacitor/local/WALL-DIRICHLET/TEST-NOFLUX-ADV (16x16x18, paredes z = 1/18, electrolito 1.11e-2,
  LB térmico kT 1e-5, centred, 2000 pasos, ~40 s cada una) + control periódico sin paredes.
- Observación: total de iones en el fluido rel. al paso 1: antes entre +8.3e-6 y −4.1e-5 (−4.05e-5 en 2000);
  después |Δ| ≤ 7e-15. Sin paredes: idéntico byte a byte (21/21 qsi). RMS u_z primera capa ≈ 4e-3 en ambos.
  Primer intento sin paredes abortó (electroneutral no exige paredes): sinpared_fallido.
- Interpretación: la pérdida venía de los enlaces fluido-pared de la advección centred; corregida.
- Pendiente: no probado con coloides BBL ni MPI; no medido el efecto en perfiles cerca de la pared;
  identificar corridas anteriores afectadas (centred + paredes/BBL + iones + velocidad normal).

## 2026-11-16 — Registros de simulaciones actualizados + columna "Scripts"
- Pregunta previa: la versión CPU de la advección centred ya tiene los pesos (np_advective_fluxes, l. 1131,
  CHANGE 20261112: u *= 2 wgradients[p]). Observación de lectura de código; la CPU no se usa (driver GPU desde 2026-03-26).
- Cambio: util/registro_simulaciones.py: nueva columna "Scripts" (función run_scripts, regla heurística documentada en
  la hoja Leyenda); regenerados capacitor/Registro_simulaciones.xlsx (661 -> 680 sim.: + TEST-NOFLUX-ADV 5 y
  H-microgel-ruido 14, todas "solo local") y microgel/Registro_simulaciones.xlsx (1202 -> 1251: + 49 de
  TEST-HIDRO-dx2 advfix/advfix_flla/advfix_flla_mass/cont6000). Versión anterior de ambos en git (dc7f2059).
- Observación: Scripts específicos en 129 (microgel) y 121 (capacitor) corridas; "todos los de la serie" 149 / 557;
  "no identificado" 940 / 2 (series viejas sin carpeta scripts/<serie>/ o scripts que arman los nombres con formatos).
- Limitación: la regla busca nombres literales; no detecta scripts que construyen la ruta (p.ej. correr.sh con $1,
  crear.py con f-strings). No es un registro de lo efectivamente ejecutado.
- Pendiente: las corridas nuevas están "solo local" (falta util/archive_runs.py si se quiere espejar en runs/).

## 2026-11-16 — Documento del modelo: flujos de NP de centred corregidos
- Pedido: corregir en claude/MODELO_EK_PETSC_PASO_A_PASO.md la descripción de los flujos de NP de centred.
- Cambio (solo documento): §7.2 intro y §7.2.1 (nodos sólidos: el flujo se calcula y luego se anula; historia
  del arreglo), §7.3 reescrita como secuencia de 5 pasos (advección, difusión, no-flujo, antisimetría, Euler) +
  condición de conservación, §7.4 (punto 3: flla también tiene limitador; punto 5 retirado con nota),
  tabla de §1 (flla), paso 14 de §3, F_p -> J^adv, números de línea actualizados tras el CHANGE 20261116.
- No verificado: render del Markdown.

## 2026-11-16 — Corrida con parámetros de Ingo + sal (lambda_D = 2 l_B), advección centred, 3e4 pasos
- Objetivo: simulación equivalente a LONG-MASS-3E5 pero con sal tal que lambda_D sea comparable a l_B.
- Decisiones del usuario: lambda_D = 2 l_B = 5.68; base LONG-MASS-3E5 con centred; 3e4 pasos.
- Cambio: microgel/scripts/INGO-1PART-LBTERM-DIAG/crear_sal.py -> local/INGO-1PART-LBTERM-DIAG/LONG-SALT-LD2LB-CEN-3E4
  (rho_el 4.3426e-4 = 1/(8 pi l_B lambda_D^2), advection centred, N_cycles 3e4; resto idéntico). Exe 269ec1b8.
- Observación: terminó normalmente (17.7 min); log Debye 5.6799696, Bjerrum 2.84. N+, N- y Q = -1 conservados a
  8 decimales; sin negativos; n+ en [2.31e-4, 5.56e-4] y n- en [3.81e-4, 8.70e-4] en el paso 3e4 (medias 4.34e-4 /
  4.85e-4). Aviso "key not used (Line 72): #" = línea de comentario indentada heredada del input base (mismo aviso en
  LONG-MASS-3E5, línea 63); su texto ("total charge -1 ...") ahora queda debajo de la línea de la sal.
- Interpretación: ninguna todavía (sin análisis del perfil ni de fluctuaciones).
- Pendiente: decidir análisis (perfil radial vs PB con sal, fluctuaciones lejos de la partícula); registro microgel regenerado.

## 2026-11-16 — Explicación de cloud_profile.py y pb_reference.py (sin cambios)
- Observación (lectura): cloud_profile usa solo la especie 1 (qsi col. 1) y NBAR fijo = 5.0805e-5 (contraiones
  solos); pb_reference es PB discreta canónica de UNA especie (z = -1, N = 1) con la partícula fija en R0 inicial.
- Interpretación: para LONG-SALT hay que adaptar ambos (dos especies, n̄ por especie, N+ y N- fijos). No hecho.

## 2026-11-16 — Análisis de LONG-SALT-LD2LB-CEN-3E4 (perfil vs PB con sal) + CSV para ParaView
- Cambio: nuevos microgel/scripts/INGO-1PART-LBTERM-DIAG/pb_reference_sal.py y cloud_profile_sal.py (copias a
  2 especies de pb_reference.py / cloud_profile.py; originales sin tocar). CSV ParaView con util/extract_colloids
  en LONG-SALT-LD2LB-CEN-3E4/colloides_paraview (31 pasos).
- Elecciones mías (declaradas): PB resuelta en la posición de cada volcado (además de la inicial); ventana 15000-30000;
  tolerancia Picard 1e-12 (igual que el original).
- Observación: PB n-(r<1) 2.69 vs Ludwig 1.49; n+(r<1) 0.37 vs 0.65; carga iónica en r<6: PB -0.285, Ludwig -0.12..-0.16.
  ex-(r<3) Ludwig sube de 0.11 (1000) a ~0.3 hacia 1e4 y fluctúa 0.19-0.39 después. Ruido lejano 0.06 rel., corr. con
  fluido ~0.03 (como NOPART-W). PB a r0 vs PB por volcado: casi igual (dif. < 0.01 por bin).
- Interpretación (no verificada): la nube de Ludwig tiene ~1/3 del exceso de PB, como en LONG-MASS sin sal (r<1: 1.47 vs
  3.18). Candidatos: movimiento de la partícula (D_part = D_ion), ruido advectivo, o que no está estacionaria
  (tau = lambda_D^2/D ~ 4.7e4 pasos > 3e4). Además la PB con sal no se verificó contra una corrida sin ruido.
- Pendiente (propuesto, no corrido): gemela NOFLUCT con partícula fija (verifica PB con sal) y gemela térmica fija
  (separa movimiento de ruido); corrida más larga.

## 2026-11-16 — Serie SAL-SERIE armada (sin correr) + pendientes en PENDIENTES.xlsx
- Objetivo: ver cómo cambia la nube con D_ion, kT y dx/2, y si sigue a la partícula.
- Decisiones del usuario: kT menor con eps ∝ 1/kT (l_B fijo); D×10, D×30, kT/4, kT/16; dx/2 con escalado difusivo
  (dt/4) y 1.5e4 pasos físicos; gemelas fija-térmica en todos y fija-sin-ruido en la base.
- Cambio: crear_serie_sal.py + correr_serie_sal.sh (microgel/scripts/INGO-1PART-LBTERM-DIAG, enlaces en local);
  12 carpetas en SAL-SERIE/ con README. Verificado con diff que cada input/config solo cambia lo previsto.
  Primer intento abortó (línea de cabecera > 120 bytes en DXH-FIJA): carpetas recién creadas borradas y recreadas.
  PENDIENTES.xlsx: filas 28 (correr) y 29 (analizar); copia previa en el scratchpad.
- Derivación dx/2 (mía, aceptada por el usuario): kT_LU x2, eps/4, n_s/8, mu_b(LU)/2 con a_L fijo -> a_h 0.379664;
  el a_h físico no se puede conservar (mu_b < 0).
- Pendiente: correr; generalizar los scripts de análisis a L, l_B y kernel variables; definir la métrica de seguimiento.

## 2026-11-16 — Etiquetas de C_refinamiento_anidado_hann.png + ROOT de grafico_C_formulas.py
- Pregunta: qué magnitud muestra el gráfico. Respuesta (lectura del código): potencial electrostático psi = e phi/k_BT;
  arriba RMS por banda de psi_Ludwig - psi_analítico (carga Hann continua, capacitor con electrolito, eje rho = 0);
  abajo orden local log2(RMS_dx/RMS_dx/2).
- Cambio: título y ejes del gráfico (grafico_C_refinamiento_anidado.py); PNG regenerado (anterior en el scratchpad,
  mismos datos). grafico_C_formulas.py: ROOT absoluto (local/WALL-DIRICHLET); desde el traslado del 2026-11-13 apuntaba a
  capacitor/scripts y todos los scripts que lo importan fallaban al buscar las corridas (FileNotFoundError).
- Comprobación: el PNG nuevo tiene las mismas curvas que el viejo (comparación visual); los órdenes impresos coinciden con
  los del gráfico. No revisé los demás scripts que importan grafico_C_formulas (solo que ahora no deberían fallar por ROOT).

## 2026-11-16 — dmu_u_vs_viscosidad.png: agregados Δμ con pesos y Δμ sin hidrodinámica
- Pedido: agregar Δμ tras la corrección de pesos y Δμ sin hidro; ¿|u| es el máximo?
- Respuesta: el panel de abajo tiene max y RMS sobre los nodos de |u| (módulo de la velocidad por nodo), paso 3000.
- Cambio: grafico_dmu_u_visc.py (rutas absolutas; panel superior: viejo = guiones/vacío, con pesos (advfix) =
  sólido, sin hidro (hidro_no_exe45dc94bc, paso 3000) = línea horizontal punteada; texto con el valor al paso 9000).
  Figura anterior conservada como graficos/dmu_u_vs_viscosidad_v1_20261110.png.
- Observación (paso 3000): con pesos Δμ baja 6.8–10× (η 0.02: 6.0e-6/9.5e-6 -> 8.9e-7/9.4e-7; η 2: 1.9e-7/2.6e-7 ->
  3.4e-8/3.1e-8). Sin hidro: 1.76e-8/1.81e-8 en 3000 (todavía bajando), 4.4e-14 en 9000.
- Interpretación: con pesos, en η = 2 Δμ queda ~2× por encima del valor sin hidro en el mismo paso; el sin-hidro sigue
  bajando, así que la comparación al paso 3000 no es de estados estacionarios.

## 2026-11-16 — dmu_u_vs_viscosidad.png en estado estacionario (sin corridas nuevas)
- Pedido: correr más pasos para que los puntos sean estacionarios.
- Observación: no hizo falta correr: cont6000/old y cont6000/centred (corridas hoy) cumplen |Δμ(6000)-Δμ(5000)|/Δμ < 1e-4
  (estacionario_6000.txt); hidro_no_exe45dc94bc decae ×~2.85 cada 250 pasos hasta 6250 y queda constante en 4.44e-14
  hasta 9000.
- Cambio: grafico_dmu_u_visc.py usa paso 6000 (hidro) y 9000 (sin hidro); |u| de centred, old igual a 3.8e-6.
  Versión paso 3000 guardada como graficos/dmu_u_vs_viscosidad_v2_paso3000.png.
- Valores (6000): sin pesos Δμ± 6.03e-6/9.48e-6 (η 0.02) -> 1.90e-7/2.56e-7 (η 2); con pesos 8.87e-7/9.35e-7 -> 2.91e-8/3.05e-8.
- Interpretación: 4.4e-14 ≈ 200 ε_máquina, probablemente nivel de redondeo (no verificado). Con hidro, Δμ con pesos en η = 2
  sigue ~6 órdenes por encima de sin hidro: el piso con hidro no se explica por el redondeo.

## 2026-11-16 — dmu_u_vs_viscosidad.png: recuadro con los parámetros de la simulación
- Cambio: grafico_dmu_u_visc.py agrega a la derecha caja 24^3 (48 físico), 2 cargas fijas +1 Hann-4 en x = 4.25/6.25
  (separación 4 físico), sal 1:1 n = 6.175e-3/nodo, l_B 0.398 (0.796 físico), lambda_D 4 (8 físico), kappa l_B 0.10,
  kT 1e-5, eps 2e4, D input 0.01 (físico 0.0614, A0 6.138), rho 0.8, M10 sin ruido, FFT, stencil 27, centred.
  Fuente: hidro_si/input y run.log. Verificado con diff de inputs y cmp de config.cds.init que solo cambian viscosidad,
  ejecutable e hidrodinámica entre las corridas del gráfico. Título aclara r = separación.

## 2026-11-16 — C_convergencia_desplazada_fuerza_todos.png: curva y puntos más finos
- Cambio: grafico_C_convergencia_desplazada.py (figure_todos): F_kernel lw 2.5 -> 0.8 (gris oscuro, encima), marcadores
  ms 5-6 -> 2.5-3, líneas 1.0 -> 0.7. Mismos datos. PNG anterior en el scratchpad. Se corre con
  `python3 scripts/grafico_C_convergencia_desplazada.py todos` desde capacitor/local/WALL-DIRICHLET (~2-3 min).

## 2026-11-16 — C_convergencia_desplazada_fuerza_todos.png: recuadro con lambda_D y caja
- Observación (inputs + run.log de las 3 series *-desplazada): físico (dx = 1): caja 32x32 periódica x 32 entre placas
  (psi = 0), lambda_D = 6, l_B = 0.796, n = 1.389e-3, q = +1 fija, kernel ancho 8. En nodos: dx = 2/1/1/2 -> caja 16/32/64,
  lambda_D 3/6/12, l_B 0.398/0.796/1.59 (Ludwig imprime esos valores).
- Cambio: grafico_C_convergencia_desplazada.py (figure_todos) agrega el recuadro en el panel superior. Mismos datos.

## 2026-11-16 — particula_libre_pared_steele.png: sin panel (c), leyenda movida, texto al ancho de los paneles
- Cambio: grafico_pared_blanda.py: 2 paneles (se quita "wall charge"; el rango de carga se sigue imprimiendo), leyenda de (a)
  abajo a la izquierda, etiqueta de la línea 3.5 debajo de la línea, título y texto de fórmulas partidos en líneas para no
  exceder el ancho de los paneles, ROOT absoluto. Regenerado con STEELE=1 (PNG anterior en el scratchpad).
  particula_libre_pared_blanda.png (STEELE=0) no se regeneró: con el código nuevo también perdería el panel (c).
- Observación: la curva azul (esfera blanda) queda tapada por la violeta (Steele): trayectorias iguales dentro de 1.8e-3
  en 121 puntos (paradas en d = 3.6055 y 3.6066).

## 2026-11-16 — particula_libre_pared_steele.png: "Steele repulsion potential", sin esfera blanda, carga q = +1; estimación para q mayor
- Cambio: grafico_pared_blanda.py (STEELE=1): RUNS sin la corrida de esfera blanda; textos "Steele repulsion potential";
  título y texto con la carga q = +1 (q0 = 1 en config.cds, l_B 0.796, sin electrolito). PNG anterior en el scratchpad.
- Pregunta: ¿qué pasa con mayor carga? Estimación (NO simulada): sin electrolito Poisson es lineal -> atracción imagen ∝ q^2.
  Equilibrio F_Steele(h) = q^2 F_el: h = 0.100/0.088/0.082/0.075/0.069/0.066 para q = 1/2/3/5/8/10 (d > 3.5: sin pérdida de
  carga). Rigidez k = -dF/dh crece ~q^2.2. Movilidad efectiva medida mu ≈ 1.9e3 (v/F, pasos 31000-34000) -> mu k = 0.03 (q=1),
  1.1 (q=5), 3.1 (q=8): Euler explícito sobreamortiguado estable si mu k < 2 -> límite q ≈ 6-7 con eps_w fijo; con eps_w ∝ q^2
  (para parar en d = 3.6) mu k = 0.033 q^2 -> q ≈ 7-8. Hipótesis: atracción constante cerca de d = 3.6, mu sin cambio.

## 2026-11-16 — C_hoja_hann_closed_vs_ludwig.png y _2pot_: solo gráfico, fórmulas y parámetros
- Cambio: grafico_C_hoja_hann.py: se quitan el panel de residuos (Ludwig − fórmula; se siguen imprimiendo) y la nota de la
  derivación (Jackson 3.17, sumas |n|,|m| <= 96); quedan el panel de potencial (Ludwig, Hann-8, gaussiana), las fórmulas y la
  línea de parámetros. Regenerados los dos PNG (anteriores en el scratchpad). Las versiones *_ppt.png (LAYOUT=left) no se
  regeneraron.

## 2026-11-16 — C_hoja_hann_closed(_2pot)_vs_ludwig.png: parámetros arriba de las fórmulas, texto al ancho del gráfico
- Cambio: grafico_C_hoja_hann.py (layout normal): bloque de texto nuevo, una ecuación (o dos cortas) por línea, centrado en
  el eje; "Parameters" (3 líneas) arriba del título de las fórmulas; título del gráfico a 10 pt. Mismas fórmulas, sin
  cambios de contenido (se separaron las que iban juntas). Layout "left" (ppt) sin cambios. PNG anteriores en el scratchpad.

## 2026-11-16 — C_desplazada_ref2.png sin los paneles de RMS por banda
- Cambio: grafico_C_desplazada.py: en modo ref2 la figura tiene 4 filas (sin "band RMS" ni "RMS_b/RMS_b"); las bandas se
  siguen calculando e imprimiendo (las usa el chequeo de estado estacionario: ninguna corrida dio "NOT steady"). Los otros
  modos (cont, celda, placa, normal) no cambian. PNG anterior en el scratchpad.

## 2026-11-16 — C_convergencia_desplazada_fuerza_todos en "a -Mostrar" sin el panel de orden
- Cambio: grafico_C_convergencia_desplazada.py: argumento "sinorden" (con "todos") -> 3 paneles, sin "observed order p",
  en graficos/C_convergencia_desplazada_fuerza_todos_sin_orden.png; la figura de 4 paneles de graficos/ no cambia.
  Copiada sobre "a -Mostrar/C_convergencia_desplazada_fuerza_todos.png" (la anterior, idéntica a la de graficos/, en el
  scratchpad).

## 2026-11-17 — C_potencial_carga_vs_placas_infinitas.png (C-sin-fondo-hann08): referencia con kernel Hann-8
- Cambio: grafico_C_potencial_carga.py: nueva psi_hann (carga Hann-8 continua entre placas a tierra, caja periódica 32x32,
  sin fondo: Fourier en x,y con W^(kx)W^(ky) y convolución en z con W(s-d), G_k en forma de exponenciales estables, G_0
  lineal). Reemplaza en el panel (a) a "placas infinitas" y "periódica puntual". Textos en inglés. ROOT absoluto.
  PNG anterior en el scratchpad.
- Prueba: lejos del soporte psi_hann = psi_periodic (puntual) a ~1e-6; M = 60 vs 80: 1e-7.
- Observación: max|Ludwig - teoría Hann| fuera del soporte 0.9-1.0e-4 (d = 4, 8, 12), dentro 1.3e-2 (mismo orden que en
  C_hoja_hann con electrolito: 3e-4 / 1.2e-2).
- Interpretación: la diferencia dentro del soporte es probablemente el kernel discreto (Ludwig) vs continuo (teoría); no
  verificado aquí.

## 2026-11-18 — Pregunta conceptual: ¿la nube deformada por el ruido implica que faltan correcciones de muchos cuerpos / no linealidad?
- Respuesta (sin cálculos nuevos; datos de INGO-1PART-LBTERM-DIAG): el déficit de la nube depende de D_ion (fija D10 0.123 vs
  D30 0.037, ~1/D) y de si la partícula se mueve. En equilibrio la distribución estática (Boltzmann) no depende de D, de eta
  ni de la movilidad -> el estado estacionario de Ludwig con LB térmico NO es de equilibrio: el ruido entra solo por la
  advección y NP no tiene su propio ruido (no se cumple fluctuación-disipación para los iones).
- Interpretación: no es un efecto de muchos cuerpos (aparece con una sola partícula) ni una falla de las correcciones de corto
  alcance (que corrigen la discretización de una interacción estática); tampoco requiere no linealidad (en la imagen de
  T efectiva el efecto es lineal). Las correlaciones físicas más allá de PB necesitarían ruido iónico (Péraud 2016).
- Pruebas propuestas (no hechas): variar eta con partícula fija (en equilibrio no debería cambiar la nube); dos partículas
  fijas con LB térmico vs suma de una partícula (aditividad con ruido).

## 2026-11-18 — Aclaración: interacción entre VARIAS partículas, cada una con su nube
- Respuesta (conceptual, sin cálculos nuevos): en Ludwig la parte mediada por iones (muchos cuerpos, campo medio) sale de
  Poisson + NP en la malla; pm_sr solo corrige el error de malla de las interacciones directas partícula–partícula y
  partícula–nodo (estáticas, lineales, referencia de ancho del kernel). NBODY confirmó aditividad solo de la parte lineal
  sin sal. Con ruido, la nube deformada cambia la interacción mediada por iones de forma dependiente del estado (D, eta,
  movilidad, densidad): no se puede absorber en correcciones de par.
- Pendiente (propuesto, no hecho): NBODY con sal sin ruido (no aditividad de PB) y con ruido; F(r) de dos partículas fijas
  con/sin ruido vs PB.

## 2026-11-18 — Pregunta: MD con monómeros y contraiones explícitos + sal implícita
- Estimación (sin simulación; supuestos: partículas de 1 nm, fracción de volumen 0.1 hinchado, caja L = 2R + 10 lambda_D):
  R = 20 nm -> ~6e3 monómeros, 3e3-1.5e4 iones de sal (1-100 mM); R = 50 nm -> 1e5 monómeros, 1e4-1.6e5 iones;
  R = 100 nm -> 8e5 monómeros, 3e4-1e6 iones. Hasta ~10 mM los iones de sal explícitos son <= los monómeros: la sal
  implícita ahorra poco; pesa a 100 mM o en cajas grandes.
- Interpretación: el costo lo fijan los monómeros, no la sal. Sal implícita por DH: válida si rho_f/2c_s pequeño o c_s bajo;
  falla en el medio (partición de Donnan). Sal como campo PNP (Ludwig/Ingo invertido) capta Donnan y dinámica.

## 2026-11-18 — Pregunta: ¿capacitor con sólido poroso en Ludwig o MD explícita?
- Respuesta (conceptual, sin cálculos): depende del régimen. Nanoporos (< ~2-3 nm), 1 M / líquido iónico, voltajes de V:
  lambda_D < tamaño iónico, sin estérico ni correlaciones ni metalicidad -> MD con electrodos a potencial constante.
  Mesoporos/macroporos (poros >> lambda_D), sal diluida, psi ~ pocos kT/e, interés en carga/transporte (tiempos us-s) ->
  continuo (Ludwig) tiene sentido y MD no llega.
- Observación sobre el modelo actual: el sólido de partículas subgrid no excluye iones ni fluido y no es conductor; para un
  electrodo poroso haría falta sólido como nodos del mapa (exclusión + no-flujo) y Dirichlet en esos nodos (hoy solo planos).

## 2026-11-18 — Corrección: el sólido entre las placas NO es el electrodo
- Rehecha la respuesta: electrodos = paredes planas a potencial fijo (ya implementadas); el sólido es una matriz porosa
  entre ellas. Observación de código: Ludwig tiene porous_media_init (wall_*, circle_xy, square_xy, simple/fcc/bcc cubic)
  y medios porosos desde archivo, y psi_petsc.c tiene un solver de epsilon variable (var_epsilon), no probado con paredes.
- Interpretación: la red de partículas subgrid no excluye iones ni fluido (gel permeable, no sólido) y no tiene contraste
  dieléctrico. Ludwig sirve para poros >> lambda_D y sal diluida; MD para nanoporos/electrolito concentrado.

## 2026-11-18 — Resumen pedido: ¿el modelo no tiene sentido ni por costo ni por precisión?
- Respuesta: parcialmente. No tiene sentido usar partículas subgrid para representar objetos de escala molecular (monómeros,
  sólido rígido): sin volumen excluido, kernel más ancho que el objeto, costo dominado por las partículas (SOLIDO-L48: 2544
  partículas, ~3.4 s/paso). Sí tiene sentido el continuo PNP + LB a escala mesoscópica (poros/geles >> lambda_D, sal diluida,
  tiempos largos) con el sólido como mapa de nodos o el gel como campo continuo. Decisión pendiente del usuario: régimen físico.

## 2026-11-18 — Búsqueda bibliográfica: ¿el modelo continuo no sirve para el gel a resolución de monómero?
- Observación (abstracts verificados, textos completos no leídos): no encontré literatura que diga que no sirve en general.
  PB funciona para la hinchazón y los contraiones confinados de geles monovalentes (Claudio-Holm-Kremer 2009; Landsgesell
  et al. 2019). Falla en: condensación más fuerte que PB (Deserno-Holm-May 2000), correlaciones a alto acoplamiento
  (Netz-Orland 2000), apiñamiento a alto voltaje/concentración (Bazant et al. 2009).
- Interpretación: corrige mi afirmación previa de que el continuo es inválido a resolución de monómero; los argumentos
  sostenibles son costo/tamaño y los regímenes específicos. Memoria de referencia guardada.

## 2026-11-18 — Documento para el director: escala del modelo de micro/nanogel en confinamiento
- Pedido: fundamentos para discutir cambiar el modelo (MD) o la tesis.
- Cambio: página privada https://claude.ai/artifact/M4UBQVgAWMJRbu3zYzHHHv (fuente en el scratchpad). Contenido: tamaño físico
  del kernel vs Δx, costo Ludwig (medido) y MD (estimado), literatura (PB funciona para hinchazón; falla en condensación,
  correlaciones, apiñamiento), limitaciones medidas del modelo (ruido no de equilibrio, nube 1/3 con sal, pérdida de carga
  junto a pared), nanoconfinamiento, opciones A-D, lo ya hecho y preguntas para decidir. Números marcados medido/estimado.
- No verificado: textos completos de las referencias; la cifra de ESPResSo (ELC/ICC) no se citó con fuente.

## 2026-11-18 — Novedad para la tesis: búsqueda rápida de huecos (4 búsquedas, no sistemática)
- Observación: (1) ruido térmico consistente en LB-EK con iones continuos: Kuron et al. 2016 dicen que es no trivial;
  ESPResSo EK lo tiene opcional y amortiguado cerca de densidad 0; no encontré trabajo sobre el artefacto de no equilibrio
  que medimos. (2) Microgel iónico en rendija/confinamiento con electrolito: nada directo (cercanos: Alziyadi-Denton 2023
  celda, Del Monte 2019). (3) Translocación de micro/nanogel por nanoporo: simulación mecánica sin electrostática (Biswas-
  Chakrabarti 2026, capilar) y experimentos de pulso resistivo (2011); no encontré simulación con corriente iónica.
  (4) Gel entre electrodos: preprint 2024 de hidrogel polielectrolito entre electrodos (grano grueso).
- Interpretación: líneas candidatas con posible novedad; requieren búsqueda sistemática antes de comprometerse.

## 2026-11-18 — Pregunta conceptual: cómo se plantea el microgel en MD (solvente explícito o implícito)
- Objetivo: responder cómo se modela en MD un microgel con iones (qué es explícito).
- Cambio/prueba: ninguno (respuesta conceptual, sin código ni corridas).
- Resultado (interpretación, referencias de memoria sin verificar en esta sesión): el estándar es el "modelo primitivo":
  monómeros y contraiones explícitos (bead-spring tipo Kremer-Grest), solvente implícito (continuo dieléctrico, l_B) y
  termostato de Langevin; sal explícita (GCMC/reaction ensemble) o Debye-Hückel implícito. La hidrodinámica (LB acoplado,
  MPCD, DPD) solo se agrega si la pregunta es dinámica. Solvente atomístico: inviable a escala de microgel.
- Pendiente: decidir la pregunta física (equilibrio vs dinámica), porque de eso depende si hace falta hidrodinámica.
