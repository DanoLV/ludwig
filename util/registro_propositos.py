"""
Idea / propósito de cada carpeta-serie, para util/registro_simulaciones.py.
(2026-11-13). Editable a mano: es la fuente de la columna "Descripción".

Cada entrada: (proyecto, regex sobre el nombre del grupo, idea, fuente).
Se usa la PRIMERA que coincida (poner las más específicas antes).
Fuentes:
  README          = README de la carpeta (o README de WALL-DIRICHLET, tabla "Contenido")
  cabecera input  = comentarios al principio del input de las corridas
  sesiones        = conversaciones previas con Claude donde se armaron las corridas
  INFERIDO        = deducido del nombre de la carpeta y de los parámetros: VERIFICAR
"""

PROPOSITOS = [
    # ------------------------------------------------------------ capacitor
    ("capacitor", r"/A-sin-iones$",
     "Validar el solver de Poisson con paredes a potencial fijo sin iones (1 paso): perfil lineal, carga inducida y campo exactos.", "README"),
    ("capacitor", r"/A-mpi2",
     "Prueba A con 2 procesos MPI. No se pudo correr: el build de GPU exige una GPU por proceso (ver NOTA.txt).", "README"),
    ("capacitor", r"/B-debye-huckel",
     "Validar el electrolito entre paredes (λD = 4) hasta el estacionario, sin hidrodinámica, contra Debye-Hückel discreta y Boltzmann.", "README"),
    ("capacitor", r"/D-paper-asta",
     "Replicar el capacitor de la Fig. 3 de Asta et al. 2019 (L = 76, λD = 6, Δψ = 0.1): potencial y concentraciones en el estacionario.", "README"),
    ("capacitor", r"/E-paper-asta",
     "Dinámica de carga del capacitor de Asta (Fig. 4a): Q(t) con D dividido por A0, para comparar el tiempo de carga τ con la teoría.", "README"),
    ("capacitor", r"/regresion-periodica",
     "Regresión: sin las claves nuevas de pared el solver debe dar exactamente lo mismo que el original (1 paso de una corrida vieja con PETSc).", "README"),
    ("capacitor", r"/guarda-pmsr-obsoleta",
     "Verificaba que se rechazara pm_sr_pair_ref measured con paredes. Ya no aplica (ver NOTA.txt).", "README"),
    ("capacitor", r"/C-carga-imagen",
     "Una partícula subgrid q = +1 fija entre paredes a tierra, sin sal, en función de la distancia d a la pared (1 paso por posición): "
     "carga inducida vs fórmula de reciprocidad y campo de imagen vs serie de Fourier.", "README"),
    ("capacitor", r"/C-sin-fondo-corr",
     "Igual que C-sin-fondo (sin fondo de contraiones, carga neta +1) pero con la corrección pm_sr activa, para comprobar que no altera el caso de una partícula.", "README"),
    ("capacitor", r"/C-sin-fondo",
     "Igual que C-carga-imagen pero SIN fondo de contraiones (electroneutral no): carga neta +1 y las paredes aportan −1.", "README"),
    ("capacitor", r"/C-corr-hann",
     "Igual que C-carga-imagen pero con la corrección pm_sr activa (referencia medida), para comprobar que no altera el caso de una partícula.", "README"),
    ("capacitor", r"/C-dospotenciales",
     "Partícula q = +1 fija entre paredes a potenciales DISTINTOS (0.1 y 0.3), sin fondo de contraiones, barrido de d: campo aplicado + imagen.", "INFERIDO (parámetros vs C-carga-imagen)"),
    ("capacitor", r"/C-electrolito-dospotenciales",
     "Partícula q = +1 fija con electrolito 1:1 entre paredes a potenciales distintos (0.1 y 0.3), para distintos λD.", "INFERIDO (parámetros)"),
    ("capacitor", r"/C-electrolito-.*-desplazada-cont",
     "Continuación (más pasos) de corridas de C-electrolito-hann8-LD6-desplazada.", "INFERIDO (nombre)"),
    ("capacitor", r"/C-electrolito-.*desplazada",
     "Carga q = +1 fija con electrolito, desplazada una distancia δ desde el plano medio del canal hacia la placa superior (barrido de δ): "
     "fuerza neta cerca del centro, para una resolución de red / kernel dados.", "cabecera input"),
    ("capacitor", r"/C-electrolito-.*-ref[0-9.]+(-q0.25)?$",
     "Convergencia con la resolución de red: mismo sistema físico que C-electrolito-hann8-LD6/z_09.50 (L = 32, λD = 6, d = 8) con dx = 1/ref "
     "y un kernel del mismo ancho físico (Hann, Hann cos^4 'p4', Gauss o Peskin-6 estirado). Las '-q0.25' repiten con q = 0.25 (linealidad).", "cabecera input"),
    ("capacitor", r"/C-electrolito-hann[68]-LD[468]$",
     "Partícula q = +1 fija entre paredes a tierra con sal 1:1 (λD = 4, 6, 8), barrido de d: fuerza de imagen y cargas inducidas vs "
     "Debye-Hückel con número de iones conservado (5000 pasos, estacionario).", "README"),
    ("capacitor", r"/C-electrolito-hann8p4-LD6-ref1",
     "Como C-electrolito-hann8-LD6 con kernel Hann cos^4 (p4).", "cabecera input"),
    ("capacitor", r"/C-lineal",
     "Test de linealidad: igual que C-electrolito-hann8-LD4/LD8 con q = 0.25 y 2; F/q² y Q/q deben ser independientes de q si Debye-Hückel lineal es exacto.", "README"),
    ("capacitor", r"/C-particula-libre",
     "Partícula cargada LIBRE (isfixedr = 0) entre paredes a tierra, sin o con electrolito: ver su movimiento por la atracción de la carga imagen. "
     "Variantes: viscosidad, hasta la pared, distintos potenciales de pared (Steele, blanda, d = 3.5), equilibrado.", "cabecera input + INFERIDO (variantes por nombre)"),
    ("capacitor", r"/F-campo-tangencial",
     "Con placas a potencial constante y un campo externo e0 paralelo a las placas: ¿hay problemas cuando una partícula subgrid o un polímero cruza la frontera periódica en x o y?", "README"),
    ("capacitor", r"/G-microgel-relax",
     "Relajación de un microgel chico (20 monómeros, red al azar) entre placas a ψ = 0, neutro y cargado.", "README"),
    ("capacitor", r"^SOLIDO",
     "Sólido poroso neutro (L = 48): relajación a equilibrio mecánico de una red de puntos subgrid con enlaces armónicos (sin carga).", "cabecera input"),
    # ------------------------------------------------------------ microgel
    ("microgel", r"^INGO-1PART-LBTERM-DIAG",
     "Diagnosticar la densidad negativa de contraiones de INGO-1PART-LBTERM (advección por el LB térmico) y probar los esquemas de advección "
     "(peso corregido, flla, flla_mass), la nube de contraiones vs Poisson-Boltzmann y la estadística con los parámetros de Ingo.", "README"),
    ("microgel", r"^INGO-1PART-LBTERM",
     "Una carga subgrid q = +1 con los parámetros del test de Ingo Tischler (ESPResSo), LB térmico y contraiones continuos advectados: "
     "ver cuánto mueven las fluctuaciones térmicas a la partícula y a la nube.", "README"),
    ("microgel", r"^FFT-Force-.*-halfint",
     "Fuerza entre 2 partículas (FFT, Hann-4, sin corrección) en separaciones semienteras.", "INFERIDO (nombre)"),
    ("microgel", r"^FFT-Force-.*-rcutjump",
     "Fuerza entre 2 partículas (FFT, Hann-4) en separaciones alrededor del radio de corte: buscar un salto en r_cut.", "INFERIDO (nombre)"),
    ("microgel", r"^FFT-Force-.*corr_yes",
     "Fuerza entre 2 partículas q = +1 fijas vs separación (solver FFT, Hann-8) CON corrección pm_sr de σ = 1.446; la variante -rhoel usa otra concentración de sal.", "INFERIDO (nombre y parámetros)"),
    ("microgel", r"^FFT-Force-",
     "Fuerza entre 2 partículas q = +1 fijas vs separación con el solver FFT, sin corrección pm_sr, para un kernel Hann dado (primeras baterías de fuerza de pares).", "INFERIDO (nombre y parámetros)"),
    ("microgel", r"^PETSC-LD8",
     "Batería de fuerza de 2 partículas (fija en x = 8, móvil en 8 + r, r = 1…24) con el solver PETSc en GPU, λD = 8 exacto, Hann-8, con y sin corrección pm_sr.", "sesiones"),
    ("microgel", r"^SLIDES-DIR",
     "Fuerza de 2 partículas con la separación en otras direcciones de la red (no solo x), r = 1, 2, 3, 4 y 6 (r = 6 fuera del corte, como control).", "sesiones"),
    ("microgel", r"^SLIDES-SHIFT",
     "Batería de fuerza de 2 partículas con la partícula fija en x = 8.25 en vez de 8: cuánto cambia la fuerza según la posición respecto de la red.", "sesiones"),
    ("microgel", r"^SLIDES-NONINT",
     "Fuerza de 2 partículas en separaciones no enteras (la variante INTERP interpola la referencia de la corrección).", "INFERIDO (nombre)"),
    ("microgel", r"^SLIDES-NOSALT",
     "Batería de fuerza de 2 partículas sin sal añadida (λD efectivo ≈ 74), con y sin corrección pm_sr (MEAS = referencia medida, RCUT6 = r_cut 6).", "sesiones"),
    ("microgel", r"^SLIDES-.*(sigmafit|kerneltarget|-sigma_)",
     "Batería de fuerza de 2 partículas con la corrección pm_sr apuntando a un σ distinto del del kernel (σ ajustado o objetivo): "
     "ver si la corrección alcanza su propio objetivo.", "sesiones"),
    ("microgel", r"^SLIDES-.*MEAS",
     "Batería de fuerza de 2 partículas vs separación con la corrección pm_sr de referencia MEDIDA (MEAS; RCUT6 = r_cut 6; gauss08eq = kernel gaussiano "
     "de igual ancho que Hann-8): error de la fuerza respecto de la referencia.", "sesiones"),
    ("microgel", r"^SLIDES-.*-dx2-r4-shift",
     "Mismo sistema que SLIDES-LD8-hann04-corr_no-dx2/r_04 (dx = 2, separación 4) con las DOS cargas corridas en x: efecto de la posición respecto de la red.", "cabecera input"),
    ("microgel", r"^SLIDES-.*-(dx2|dxhalf)",
     "Fuerza de 2 partículas en el mismo sistema físico que SLIDES-LD8-hann08-corr_no-N10000 pero con otra resolución de red (dx = 2 ó 1/2, kernel de igual ancho físico): convergencia en dx.", "cabecera input"),
    ("microgel", r"^SLIDES-.*-N10000-d",
     "Separaciones intermedias (entre 4–5 u 8–10) de la batería N10000, moviendo solo la carga móvil.", "cabecera input"),
    ("microgel", r"^SLIDES-.*-N10000",
     "Repetición de SLIDES-LD8-hann08-corr_no con 10000 pasos (antes 3000) para llegar al estado estacionario.", "cabecera input"),
    ("microgel", r"^SLIDES-",
     "Fuerza entre 2 partículas q = +1 fijas en sal vs separación r (una fija en x = 8, la otra en 8 + r), para un λD y kernel dados, "
     "con o sin corrección pm_sr de corto alcance: comparar con la fuerza de referencia.", "sesiones + nombre"),
    ("microgel", r"^TEST-HIDRO",
     "¿Afecta la hidrodinámica a la fuerza de 2 partículas fijas? Misma corrida que SLIDES (r = 4) con hidrodinámica apagada/encendida, distintas viscosidades y ejecutables.", "cabecera input"),
    ("microgel", r"^SS-SCAN",
     "Barrido de la repulsión de núcleo blando entre 2 partículas liberadas (soft_sphere_epsilon, exponente ν = 6/12, distancia inicial): "
     "elegir parámetros que las mantengan separadas sin inestabilidad.", "INFERIDO (nombres y parámetros)"),
    ("microgel", r"^NBODY",
     "Aditividad de las fuerzas electrostáticas: combinaciones de 2 a 4 partículas (p12, p123, …) sin sal, fuerza en el paso 1, "
     "para comprobar que la fuerza total es la suma de las de pares (parte lineal).", "sesiones"),
    ("microgel", r"^PIN2GRID-LD8-corr_yes",
     "Como PIN2GRID pero con corrección pm_sr de σ = 0.40.", "INFERIDO (nombre)"),
    ("microgel", r"^PIN2GRID",
     "Una partícula cargada fija recorriendo una celda (x = 24.0 a 25.0, paso 0.1; caja 48, λD = 8): fuerza espuria sobre sí misma según su posición "
     "respecto de la red ('pin-to-grid'), para un kernel dado.", "sesiones"),
    ("microgel", r"^EFIELD-1PART-PETSC",
     "Campo eléctrico alrededor de una carga subgrid con y sin corrección pm_sr, con el solver PETSc.", "INFERIDO (nombre; gemela de EFIELD-1PART-LD8-hann08)"),
    ("microgel", r"^EFIELD-1PART",
     "Campo eléctrico alrededor de una sola carga subgrid (caja 48, λD = 8), con y sin corrección pm_sr; la gauss08eq es la gemela con kernel gaussiano.", "README"),
    ("microgel", r"^RING",
     "Anillo de 10 partículas cargadas: con y sin corrección pm_sr y con y sin campo eléctrico externo.", "INFERIDO (nombres de subcarpeta)"),
    ("microgel", r"^CROSS",
     "Anillo de 10 partículas cruzando la frontera periódica bajo campo externo: prueba del arreglo en el borde (con y sin arreglo; con E0 = 5 explotaron; control con campo chico).", "sesiones"),
    ("microgel", r"^MGEL-RUN",
     "Microgel de 127 monómeros (L = 48) con carga en todos o en la superficie, en sal (λD 4–8), Hann-6/8: estabilidad e hinchamiento (Rg).", "sesiones"),
    ("microgel", r"^MGEL-SURF-NOSALT|^MGEL-RED-SURF-NOSALT",
     "Microgel con los monómeros de superficie cargados y SIN sal: relajar hasta el equilibrio (Rg).", "sesiones"),
    ("microgel", r"^MGEL-SURF",
     "Microgel con los 65 monómeros de superficie a q = +1 desde la estructura equilibrada, para λD = 8, 6, 4 y 2.", "sesiones"),
    ("microgel", r"^MGEL-RED-neutro|^MGEL-REF-neutro",
     "Microgel (red) neutro: relajación de referencia sin carga.", "INFERIDO (nombre)"),
    ("microgel", r"^old/EWALD",
     "Pruebas antiguas del método de Ewald (kernel gaussiano dual; σ, α, r_c) para la fuerza sobre una carga.", "INFERIDO (nombre)"),
    ("microgel", r"^old/PM-FFT",
     "Pruebas antiguas del P3M (FFT o PETSc) con corrección de corto alcance para una partícula: tablas 3D, σ, halo, tamaño de caja.", "INFERIDO (nombre)"),
    ("microgel", r"^old/Pruebas",
     "Pruebas antiguas de resolución del PM y de FFT/Ewald con distintos kernels.", "INFERIDO (nombre)"),
    ("microgel", r"^old/dimero",
     "Pruebas antiguas de un dímero (2 partículas) con PM-FFT: interlacing y distancias.", "INFERIDO (nombre)"),
    ("microgel", r"^old/sim",
     "Simulaciones antiguas: líquido solo, pruebas de microgel y 'línea peskin' (una partícula con PETSc).", "INFERIDO (nombre)"),
    ("microgel", r"^aaa-otro",
     "Carpeta de scripts de visualización (carga, campo, velocidades) con un input suelto.", "README"),
    ("microgel", r"^\(raiz",
     "Input suelto en la raíz de microgel/local (plantilla).", "INFERIDO"),
]


# ---------------------------------------------------------------------------
# CONCLUSIONES (2026-11-13): solo donde hay un resultado registrado. Misma
# estructura (proyecto, regex del grupo, texto, fuente); se usa la primera que
# coincida. Grupo sin entrada -> columna vacía (completar a mano).
# Fuentes: README, diario (diario/diario.md, fecha), memoria (notas de Claude),
#          sesiones (conversaciones previas con Claude).
# ---------------------------------------------------------------------------
CONCLUSIONES = [
    # ------------------------------------------------------------ capacitor
    ("capacitor", r"/A-sin-iones$",
     "Perfil de ψ exacto (2e-14), carga inducida = ley de Gauss y campo exacto en todos los nodos: el solver con paredes Dirichlet es correcto sin iones.", "README"),
    ("capacitor", r"/A-mpi2",
     "No ejecutada: el build de GPU exige una GPU por proceso.", "README"),
    ("capacitor", r"/B-debye-huckel",
     "Coincide con Debye-Hückel discreta a 8e-9; Boltzmann exacto; campo junto a la pared 1.00034 del teórico.", "README"),
    ("capacitor", r"/D-paper-asta",
     "Potencial vs ec. (17) de Asta: diferencia 1.5e-3 de Δψ (error de red); Q∞/Q_DH = 0.9966. Reproduce la Fig. 3.", "README"),
    ("capacitor", r"/E-paper-asta",
     "τ = 4132.7 pasos vs 4145.2 teórico (0.3 %), con D dividido por A0: la dinámica de carga es correcta.", "README"),
    ("capacitor", r"/regresion-periodica",
     "Diferencia exactamente 0 en ψ, iones y campo: sin las claves de pared el solver es el original.", "README"),
    ("capacitor", r"/C-carga-imagen",
     "Carga inducida = fórmula de reciprocidad a 1e-14; campo de imagen vs serie de Fourier ≤ 1.5 % en el borde permitido (mediana 0.03 %). "
     "Esta prueba destapó el bug de fuerza ×2/×4/×8 en celdas de borde (corregido 2026-09-26).", "README + memoria"),
    ("capacitor", r"/C-sin-fondo-corr|/C-corr-hann",
     "Sin iones idéntico a sin corrección; con fondo difiere en 1e-7: con una sola partícula la corrección pm_sr no tiene sobre qué actuar.", "README"),
    ("capacitor", r"/C-sin-fondo",
     "Carga = reciprocidad sin fondo a 1e-14 (total −1); campo 1.2 % en el borde permitido, mediana 5e-5.", "README"),
    ("capacitor", r"/C-electrolito-hann5p4|/C-electrolito-hann10p4",
     "Hann cos⁴ en δ = 1: error de la fuerza vs su blanco continuo +2.05 % (Δx = 2) y +0.47 % (Δx = 1), orden p = 2.13; misma precisión que cos² (p = 2.01).", "diario 2026-11-09"),
    ("capacitor", r"/C-electrolito-hann(4-LD3-ref0.5|8-LD6|16-LD12-ref2)-desplazada$",
     "Dentro del soporte del kernel el error es ≈ 3.5 % del pico de ψ y ≈ 7 % del pico de E, casi independiente de δ; fuera del soporte 1e-5–1e-4. "
     "Fuerza en δ = 1 vs límite continuo: +1.93 % (Δx = 2), +0.48 % (Δx = 1), orden p = 2.01.", "diario 2026-11-05 y 2026-11-09"),
    ("capacitor", r"/C-electrolito-hann[68]-LD[468]$",
     "E_z teoría/Ludwig entre 0.983 y 0.998; Q_inf a < 2e-3. La desviación crece al bajar λD (1.05, 0.48, 0.34 % para λD = 4, 6, 8).", "README"),
    ("capacitor", r"/C-lineal",
     "F/q² y Q/q independientes de q a ~1e-3: la desviación de ~1 % con λD = 4 no es no linealidad de Debye-Hückel sino discretización de la red.", "README"),
    ("capacitor", r"/F-campo-tangencial",
     "Tras los arreglos 1–3, cruzar la frontera periódica no cambia nada: borde = ref a ≤ 5e-13 en trayectoria y ≤ 4e-14 en fuerza "
     "(antes del arreglo 3 el polímero difería 13 % en fuerza).", "README"),
    ("capacitor", r"/G-microgel-relax",
     "Neutro: equilibrio (Rg = 4.765, constante desde ~5000 pasos). Con carga en superficie: NO equilibrado a 20000 pasos (Rg crece sin frenar; "
     "la red tiene 19 modos blandos). El centro de masa no se mueve: la atracción por la imagen es despreciable a esa distancia.", "README"),
    ("capacitor", r"^SOLIDO",
     "El acople partícula-fluido explícito es inestable para k ≥ 4e-3 con Hann-4; la geometría de equilibrio depende solo de r0 = 3.5. ~3.4 s/paso con 2544 partículas.", "memoria"),
    # ------------------------------------------------------------ microgel
    ("microgel", r"^INGO-1PART-LBTERM-DIAG",
     "La densidad negativa venía de la advección por el LB térmico, amplificada ×9 (D3Q27) por un error de normalización del flujo advectivo (corregido 2026-11-12). "
     "flla_mass hace que los iones sigan exactamente la densidad del fluido. Sin fluctuaciones se llega al Poisson-Boltzmann discreto exacto. Con fluctuaciones la nube es "
     "más débil que PB y el déficit crece al bajar D_ion (r < 1: −2.5 % con D×30 y −8.4 % con D×10, partícula fija); con el D de Ingo (×1) el exceso a r < 1 es ~21 % del de PB "
     "(ventana provisional). La partícula libre difunde (⟨v⟩ ≈ 0, D ≈ 7e-4).", "diario 2026-11-12/13"),
    ("microgel", r"^INGO-1PART-LBTERM",
     "Abortó por densidad de contraiones negativa a los pocos pasos (diagnóstico en INGO-1PART-LBTERM-DIAG).", "diario 2026-11-12"),
    ("microgel", r"^SLIDES-.*(sigmafit|kerneltarget|-sigma_)",
     "La corrección pm_sr alcanza su propio objetivo σ con precisión (p. ej. σ = 0.723: error +1.22 % a r = 1); elegir mal σ apunta a la curva equivocada, no hace fallar la corrección.", "sesiones"),
    ("microgel", r"^SLIDES-.*MEAS",
     "Con corrección de referencia medida el error de la fuerza es 0.1–1.5 % a r = 1 (máximo 0.8–4 % en todo r); sin corrección −17 a −19 % a r = 1 con λD = 4.", "sesiones"),
    ("microgel", r"^SLIDES-.*-dx2-r4-shift",
     "A r fijo la fuerza depende de la posición del par en la celda con período de una celda: entre −11.6 % y −9.5 % respecto de la referencia (±1.1 % alrededor de la media −10.6 %).", "diario 2026-11-11"),
    ("microgel", r"^SLIDES-LD8-hann04-corr_no-dx2$",
     "Error de ψ en nodos (L2) 4.6–6.1 % y de E 12–15 % con Δx = 2 (vs 1.1–1.4 % y 3.2–3.9 % con Δx = 1): orden ≈ 2. "
     "OJO: corrida con el ejecutable viejo, la fuerza de la partícula 2 está duplicada por el bug del halo; usar kT·Esub/dx.", "diario 2026-11-09 + memoria"),
    ("microgel", r"^SLIDES-LD8-hann08-corr_no-N10000$",
     "Estacionario alcanzado. Error de ψ en nodos (L2) 1.1–1.4 % y de E 3.2–3.9 % con Δx = 1 (vs 4.6–6.1 % y 12–15 % con Δx = 2): orden ≈ 2.", "diario 2026-11-09"),
    ("microgel", r"^TEST-HIDRO",
     "La hidrodinámica cambia la fuerza en < 1e-6 relativo. Sin flujo el sistema relaja hasta el mínimo de F (redondeo); con flujo queda un piso de Δμ ~ 2e-6.", "diario 2026-11-11"),
    ("microgel", r"^NBODY",
     "Confirma la aditividad de las fuerzas en la parte lineal (sin sal, paso 1, σ = 0.4); no prueba la parte no lineal con iones.", "sesiones"),
    ("microgel", r"^EFIELD-1PART-LD8-(hann08|gauss08eq)",
     "Error del campo a r = 0–2: Hann-8 −23.7 % sin corrección y −1.1 % con; gaussiano de igual ancho −18.1 % y −0.8 %.", "sesiones / README"),
    ("microgel", r"^CROSS",
     "El cruce del borde funcionó con el arreglo (adv_arreglado); con E0 = 5 explotaron. Las fuerzas eléctricas están ×2 uniformes por el bug del halo "
     "(anillo siempre en celdas de borde), pero la conclusión de continuidad al cruzar se mantiene.", "memoria + sesiones"),
    ("microgel", r"^MGEL-RED-SURF-NOSALT",
     "Rg crece lentamente (τ ≈ 131000 pasos); proyección Rg∞ ≈ 7.27, equilibrio estimado cerca del paso 425000.", "sesiones"),
    ("microgel", r"^MGEL-RUN|^MGEL-SURF-hann06",
     "κR ≈ 1–4 (microgel real 10–160): régimen Donnan no lineal comparable a un microgel real con ~25 mM de carga fija en 1 mM de sal.", "memoria (análisis de escala)"),
    ("microgel", r"^old/",
     "Anteriores al arreglo del 2026-09-26: partículas cargadas en celdas de borde pueden tener la fuerza ×2/×4/×8; verificar particle_force.csv vs particle_Esub.csv antes de usarlas.", "memoria"),
]
