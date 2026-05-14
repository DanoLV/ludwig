# Scripts de análisis - Ludwig/Microgel

Descripción completa y opciones de todos los scripts en este directorio.

---

## Índice

- [Scripts de cálculo y estadísticas](#scripts-de-cálculo-y-estadísticas)
- [Scripts de análisis de campo eléctrico](#scripts-de-análisis-de-campo-eléctrico)
- [Scripts de visualización - velocidades](#scripts-de-visualización---velocidades)
- [Scripts de visualización - campo eléctrico y potencial](#scripts-de-visualización---campo-eléctrico-y-potencial)
- [Scripts de visualización - cargas y estructura](#scripts-de-visualización---cargas-y-estructura)
- [Scripts de organización y limpieza](#scripts-de-organización-y-limpieza)
- [Scripts de orquestación](#scripts-de-orquestación)
- [Referencia de archivos de entrada de Ludwig](#referencia-de-archivos-de-entrada-de-ludwig)

---

## Scripts de cálculo y estadísticas

### `calculos.py`

Calcula estadísticas del coloide desde archivos `colloids-NNNNNNNN.csv`: centro de masa, momentos de inercia, longitud de bond promedio, densidad y volumen (via convex hull).

```bash
./calculos.py -nciclos <num_ciclos> -npaso <paso> -o <archivo_salida>
```

| Flag | Requerido | Descripción |
|---|---|---|
| `-nciclos` | Sí | Número total de ciclos |
| `-npaso` | Sí | Frecuencia de salida (ciclos por paso) |
| `-o` | Sí | Nombre del archivo de salida CSV |

**Salida:** CSV con columnas `cycle, xcm, ycm, zcm, Ixcm, Iycm, Izcm, lbond, density, volume`

---

### `calculosvel.py`

Extrae datos de velocidad de la primera partícula desde archivos `colloids-NNNNNNNN.csv`.

```bash
./calculosvel.py -nciclos <num> -npaso <paso> -o <archivo_salida>
```

| Flag | Requerido | Descripción |
|---|---|---|
| `-nciclos` | Sí | Número total de ciclos |
| `-npaso` | Sí | Frecuencia de salida |
| `-o` | Sí | Nombre del archivo de salida CSV |

**Salida:** CSV con columnas `cycle, vx, vy, vz`

---

### `calculosvelfluid.py`

Combina velocidades de partícula (colloids CSV) y fluido (vel files), calculando también velocidades relativas.

```bash
./calculosvelfluid.py -nciclos <num> -ninicio <inicio> -npaso <paso> -o <salida> [--idir <dir>]
```

| Flag | Requerido | Default | Descripción |
|---|---|---|---|
| `-nciclos` | Sí | — | Número total de ciclos |
| `-ninicio` | Sí | — | Ciclo de inicio |
| `-npaso` | Sí | — | Frecuencia de salida |
| `-o` | Sí | — | Nombre del archivo de salida |
| `--idir` | No | `.` | Directorio con archivos colloids CSV |

**Salida:** CSV con columnas `cycle, x, y, z, vx, vy, vz, <vfx>, <vfy>, <vfz>, vxr, vyr, vzr`

---

### `calculosvelfluidonly.py`

Extrae velocidad del fluido únicamente (sin partículas) desde archivos vel. Para simulaciones sin coloide.

```bash
./calculosvelfluidonly.py -nciclos <num> -ninicio <inicio> -npaso <paso> -o <salida>
```

| Flag | Requerido | Descripción |
|---|---|---|
| `-nciclos` | Sí | Número total de ciclos |
| `-ninicio` | Sí | Ciclo de inicio |
| `-npaso` | Sí | Frecuencia de salida |
| `-o` | Sí | Nombre del archivo de salida |

**Salida:** CSV con columnas `cycle, vx, vy, vz, <vfx>, <vfy>, <vfz>, vxr, vyr, vzr`

---

### `extraer_posicion.py`

Extrae posición y velocidad de la primera partícula desde archivos `colloids-NNNNNNNN.csv`. Soporta procesamiento incremental (abre en modo append).

```bash
./extraer_posicion.py -nciclos <num> -ninicio <inicio> -npaso <paso> -o <archivo_salida>
```

| Flag | Requerido | Descripción |
|---|---|---|
| `-nciclos` | Sí | Número total de ciclos |
| `-ninicio` | Sí | Ciclo de inicio |
| `-npaso` | Sí | Frecuencia de salida |
| `-o` | Sí | Nombre del archivo de salida |

**Salida:** CSV con columnas `cycle, x, y, z, vx, vy, vz` (notación científica)

---

### `calc_mui.py`

Calcula el potencial químico μ_k = kT·ln(ρ_k) + z_k·e·ψ para cada especie iónica y su desviación del equilibrio (Δμ). Puede generar archivos theta-*.001-001 y plots de verificación de Boltzmann.

```bash
# Step único:
./calc_mui.py -n <nstep> -L <L> -kt <kT> [opciones]

# Serie temporal:
./calc_mui.py --n-start <inicio> --n-end <fin> --n-step <paso> -L <L> -kt <kT> [opciones]
```

| Flag | Requerido | Default | Descripción |
|---|---|---|---|
| `-n` / `--nstep` | No* | — | Step único a analizar |
| `--n-start` | No* | — | Step inicial de la serie |
| `--n-end` | No* | — | Step final de la serie |
| `--n-step` | No* | — | Incremento entre steps de la serie |
| `-L` / `--grid-size` | Sí | — | Tamaño de la grilla cúbica (L×L×L) |
| `-kt` / `--kt` | Sí | — | Energía térmica kT |
| `-d` / `--directory` | No | `.` | Directorio con archivos de Ludwig |
| `--z0` | No | `1.0` | Carga de la especie 0 (cationes, default: +1) |
| `--z1` | No | `-1.0` | Carga de la especie 1 (aniones, default: -1) |
| `--rho-min` | No | `1e-20` | Densidad mínima para evitar log(0) |
| `--theta` | No | False | Calcula θ_k = ρ_k·exp(z_k·ϕ) y escribe archivos theta-*.001-001 |
| `--theta-only` | No | False | Grafica θ_k por especie en ejes separados (implica `--theta`) |
| `--boltzmann-plot` | No | False | Scatter de exp(βz_k eψ) vs ρ_k por nodo y especie |
| `--ln-rho` | No | False | Scatter de ln(ρ_k) vs ϕ con línea de Boltzmann |
| `--csv` | No | `mui_series.csv` | Nombre del CSV de salida |
| `--plot` | No | `mui_delta_evolution.png` | Nombre del plot de salida |

*Se requiere `-n` solo, o los tres de `--n-start`, `--n-end`, `--n-step`.

---

### `qsi_subgrid_force.py`

Calcula las fuerzas sobre el coloide debidas a la distribución de carga iónica en la subgrilla. Soporta múltiples esquemas de distribución de carga en los subpuntos.

```bash
./qsi_subgrid_force.py -n <nstep> -L <L> --rcut <rcut> [opciones]
```

| Flag | Requerido | Default | Descripción |
|---|---|---|---|
| `-n` / `--nstep` | Sí | — | Número de paso temporal |
| `-L` / `--grid-size` | Sí | — | Tamaño cúbico N×N×N |
| `--rcut` | Sí | — | Radio de corte: solo voxeles dentro de esta distancia |
| `-q` / `--charge` | No | `1.0` | Carga del coloide |
| `-eps` / `--epsilon` | No | `1.0` | Permitividad dieléctrica |
| `--nsub` | No | `3` | Subdivisiones por arista de voxel (→ nsub³ subpartículas por voxel) |
| `--rmin-sub` | No | `0.0` | Radio mínimo: subpartículas más cercanas son excluidas y su carga redistribuida |
| `--rcut-sub` | No | False | Aplica rcut también a subpartículas individuales |
| `--gradient-dist` | No | False | Distribuye carga por gradiente lineal de ρ_net (diferencias centradas) |
| `--quadratic-dist` | No | False | Interpolación de Lagrange cuadrática, separable (3 nodos/eje) |
| `--quartic-dist` | No | False | Interpolación de Lagrange cuártica, separable (5 nodos/eje) |
| `--quadratic-full` | No | False | Producto tensorial cuadrático completo (3×3×3 = 27 nodos) |
| `--quartic-full` | No | False | Producto tensorial cuártico completo (5×5×5 = 125 nodos) |
| `--num-species` | No | `2` | Número de especies en archivo qsi |
| `--csv` | No | None | Nombre base para salida CSV |
| `--no-csv` | No | False | No exportar CSV aunque se especifique `--csv` |
| `--qsi-file` | No | None | Ruta explícita al archivo qsi (sobreescribe nombre estándar) |
| `--cds-file` | No | None | Ruta explícita al archivo config.cds (sobreescribe nombre estándar) |

---

## Scripts de análisis de campo eléctrico

### `analizar_campo_electrico.py`

Analiza la distribución del campo eléctrico subgrid desde archivos `particle_Esub.csv` en múltiples directorios. Extrae posiciones, calcula magnitudes medias y componentes, y genera plots comparativos.

```bash
./analizar_campo_electrico.py [directorio_base] [--no-show]
```

| Argumento | Requerido | Default | Descripción |
|---|---|---|---|
| `[directorio_base]` | No | `..` (padre del cwd) | Directorio base con subdirectorios de simulación |
| `--no-show` | No | False | Suprime la visualización interactiva de matplotlib |

**Entrada:** archivos `proceced_data/particle_Esub.csv` en subdirectorios con patrón `pos_X_Y_Z`
**Salida:** CSV con estadísticas y PNG comparativos

---

### `verificar_campo_electrico.py`

Verificación exhaustiva del campo eléctrico: compara todos los puntos simulados con el campo teórico Debye-Hückel, calculando error relativo punto a punto.

```bash
./verificar_campo_electrico.py -e <efield_file> -c <colloid.csv> -L <L> -q <carga> -eps <epsilon> [opciones]
```

| Flag | Requerido | Default | Descripción |
|---|---|---|---|
| `-e` / `--efield` | Sí | — | Archivo efield de Ludwig (ej: `efield-00000001.001-001`) |
| `-c` / `--colloid` | Sí | — | Archivo colloids (ej: `colloids-00000001.csv`) |
| `-L` / `--grid-size` | Sí | — | Tamaño de grilla (cubo L×L×L) |
| `-q` / `--charge` | Sí | — | Carga de la partícula en unidades Ludwig |
| `-eps` / `--epsilon` | Sí | — | Permitividad relativa ε_r |
| `-kt` / `--kt` | No | `1.0` | Energía térmica k_B·T (referencia) |
| `--kappa` | No | `0.0` | Parámetro de Debye κ = 1/λ_D (0 = Coulomb puro) |
| `-o` / `--output` | No | `efield_comparison.png` | Archivo de salida |
| `--rmin` | No | `0.5` | Radio mínimo de análisis |
| `--rmax` | No | L/2 | Radio máximo de análisis |
| `--max-points` | No | `10000` | Máximo de puntos a graficar |
| `--linear` | No | False | Usar escala lineal (default: logarítmica) |

**Salida:** PNG con 2 paneles: comparación de campo y métricas de error relativo

---

### `verificar_simetria_psi.py`

Verifica la simetría radial y planar del potencial eléctrico desde archivos psi, con métricas cuantitativas y comparación con teoría Debye-Hückel.

```bash
./verificar_simetria_psi.py -f <psi_file> -s <nx> <ny> <nz> (-c <cds_file> | -p <x> <y> <z>) [opciones]
```

| Flag | Requerido | Default | Descripción |
|---|---|---|---|
| `-f` / `--file` | Sí | — | Archivo psi a analizar |
| `-s` / `--size` | Sí | — | Tamaño de grilla (NX NY NZ) |
| `-c` / `--colloid-file` | No* | — | Archivo `config.cds*.001-001` o `colloids*.csv` para posición de partícula |
| `-p` / `--position` | No* | — | Posición de la partícula (X Y Z) |
| `--stencil` | No | `19` | Tipo de stencil para gradiente: `7` (D3Q7), `19` (D3Q19), `27` (D3Q27) |
| `--min-radius` | No | `0.01` | Radio mínimo para análisis radial |
| `--max-radius` | No | `15.0` | Radio máximo |
| `--num-radii` | No | `10` | Número de radios a analizar |
| `--shell-thickness` | No | `0.1` | Grosor de capa esférica (solo con `--use-shells`) |
| `--use-all-points` | No | True | Usa TODOS los puntos de grilla agrupados por bins de distancia (recomendado) |
| `--use-shells` | No | False | Método de capa esférica (alternativa a all-points) |
| `--distance-bin-size` | No | `0.05` | Tamaño de bin para agrupar distancias |
| `--planes` | No | `xy xz yz` | Planos a analizar |
| `--log-scale` | No | `none` | Escala logarítmica: `y`, `xy`, `symlog`, `symlog-xy`, `none` |
| `--symlog-linthresh` | No | `1e-6` | Umbral lineal para symlog (región [−linthresh, +linthresh]) |
| `--radial-only` | No | False | Solo análisis radial (omite planos, componentes y métricas) |
| `--efield-only` | No | False | En modo radial, solo grafica campo E y error (omite potencial); requiere `--radial-only` |
| `--save-csv` | No | False | Guarda perfiles radiales en archivos CSV |
| `--psi-offset-factor` | No | `1.1` | Factor para offset en escala log: offset = −ψ_min × factor |
| `--compare-theory` | No | False | Compara con modelo teórico Debye-Hückel |
| `--charge` | No | `1.0` | Carga de la partícula |
| `--epsilon` | No | `1.0` | Permitividad relativa |
| `--kt` | No | `1.0` | Energía térmica k_B·T |
| `--kappa` | No | `0.0` | Parámetro de Debye κ (0 = Coulomb) |
| `--no-fit-scale` | No | False | No ajusta el factor de escala teórico; usa q/(4πε·kt) directamente |
| `--theory-start-distance` | No | auto | Distancia mínima desde la que se grafica la curva teórica |
| `-o` / `--output-prefix` | No | None (pantalla) | Prefijo para archivos de salida |

*Se requiere exactamente uno de `-c` o `-p`.

---

### `compare_field_theory.py`

Compara el campo eléctrico simulado (desde archivo psi) con predicciones teóricas Debye-Hückel, en modo línea o plano.

```bash
./compare_field_theory.py -f <psi_file> -s <nx> <ny> <nz> -m <line|plane> --charge-pos <x> <y> <z> [opciones]
```

| Flag | Requerido | Default | Descripción |
|---|---|---|---|
| `-f` / `--file` | Sí | — | Archivo psi a leer |
| `-s` / `--size` | Sí | — | Tamaño de grilla (NX NY NZ) |
| `--charge-pos` | Sí | — | Posición de la carga puntual (X Y Z) |
| `-m` / `--mode` | Sí | — | Modo de comparación: `line` o `plane` |
| `--charge` | No | `1.0` | Magnitud de la carga |
| `--epsilon` | No | `1.0` | Permitividad del medio |
| `--kt` | No | `1.0` | Energía térmica k_B·T |
| `--kappa` | No | `0.0` | Parámetro κ (0 = Coulomb puro) |
| `--ionic-strength` | No | None | Fuerza iónica en mol/L (calcula κ automáticamente) |
| `--start` | No* | — | Punto de inicio de línea (X Y Z); requerido para mode=line |
| `--end` | No* | — | Punto final de línea (X Y Z); requerido para mode=line |
| `--num-points` | No | `100` | Número de puntos de interpolación a lo largo de la línea |
| `--show-nodes` | No | False | Muestra solo nodos de la red Ludwig sin interpolación |
| `-p` / `--plane` | No* | — | Plano: `xy`, `xz`, `yz`; requerido para mode=plane |
| `--position` | No | centro | Índice de posición del plano |
| `--external-field` | No | None | Componentes del campo externo a sustraer (EX EY EZ) |
| `-o` / `--output` | No | None (pantalla) | Ruta del archivo de salida |

---

## Scripts de visualización - velocidades

### `plot.py`

Grafica la evolución de la estructura del coloide: momentos de inercia normalizados, longitud de bond media y densidad vs tiempo. Sin argumentos; lee `datos.csv` directamente del directorio actual.

```bash
./plot.py
```

**Entrada:** `datos.csv` (columnas: `cycle, xcm, ycm, zcm, Ixcm, Iycm, Izcm, lbond, density, volume`)
**Salida:** Plot de 3 paneles

---

### `plotdatos.py`

Grafica posición o velocidad de partícula desde CSV con escalado automático de ejes.

```bash
./plotdatos.py -i <input.csv> [-t posicion|velocidad] [-o <salida.png>]
```

| Flag | Requerido | Default | Descripción |
|---|---|---|---|
| `-i` | Sí | — | Archivo CSV de entrada |
| `-t` / `--tipo` | No | `posicion` | Tipo de datos a graficar: `posicion` o `velocidad` |
| `-o` / `--output` | No | `{tipo}.png` | Nombre del archivo de salida |

**Entrada:** CSV con columnas `cycle, x, y, z, vx, vy, vz` (de `extraer_posicion.py` o `calculosvel.py`)

---

### `plotvel.py`

Grafica velocidades de partícula y/o fluido desde CSV con escalado automático en Y.

```bash
./plotvel.py -i <input.csv> [--out_dir <dir>] [--fluid-only]
```

| Flag | Requerido | Default | Descripción |
|---|---|---|---|
| `-i` | Sí | — | Archivo CSV de entrada |
| `--out_dir` | No | `.` | Directorio de salida |
| `--fluid-only` | No | False | Grafica solo velocidad media del fluido (1 panel en lugar de 2) |

**Salida:** `{out_dir}/velocidades.png`

---

### `plot_velocities_grouped.py`

Agrupa datos de velocidad (`datosfluid.csv`) por parámetro y genera plots combinados. Busca archivos recursivamente bajo el directorio especificado.

```bash
./plot_velocities_grouped.py -d <dir_padre> --group-by <param1> [param2 ...] [opciones]
```

| Flag | Requerido | Default | Descripción |
|---|---|---|---|
| `-d` / `--dir` | Sí | — | Directorio padre con datos de simulación |
| `--group-by` | Sí | — | Parámetro(s) por los que agrupar: `stencil`, `rhoel`, `pos`, `q`, `eps`, `kT`, `rho`, `eta`, `n`, `s`, `L`, `Lx`, `Ly`, `Lz` |
| `--plot-type` | No | `particle` | Tipo de datos: `both`, `particle`, `fluid` |
| `--out_dir` | No | `.` | Directorio de salida (si no se da `-o`) |
| `-o` / `--output` | No | auto | Ruta completa del archivo de salida |
| `--component` | No | `all` | Componente de velocidad: `all`, `x`, `y`, `z` |

---

### `plot_velocity_field.py`

Visualiza el campo de velocidad 3D desde archivos vel como cortes 2D con vectores quiver superpuestos.

```bash
# Desde archivo específico:
./plot_velocity_field.py -f <archivo_vel> [-s <nx> <ny> <nz>] [opciones]

# Desde directorio + timestep:
./plot_velocity_field.py -d <dir> -t <timestep> [-s <nx> <ny> <nz>] [opciones]
```

| Flag | Requerido | Default | Descripción |
|---|---|---|---|
| `-f` / `--file` | No* | — | Archivo de velocidades específico |
| `-d` / `--directory` | No* | — | Directorio con archivos de velocidad |
| `-t` / `--timestep` | No** | — | Número de paso temporal (requerido con `-d`) |
| `-s` / `--size` | No | auto (cubo) | Tamaño de grilla (NX NY NZ) |
| `--plane` | No | `all` | Plano a graficar: `xy`, `xz`, `yz`, `all` |
| `--index` | No | centro | Índice del plano de corte |
| `--skip` | No | `1` | Factor de submuestreo para flechas quiver |
| `--scale` | No | auto | Escala de las flechas |
| `--colormap` | No | `viridis` | Colormap |
| `--magnitude-only` | No | False | Solo magnitud en el plano, sin flechas |
| `--single-plot` | No | False | Un solo panel con magnitud y vectores (en lugar de 2 paneles) |
| `-o` / `--output` | No | `.` | Directorio de salida |
| `--prefix` | No | `vel` | Prefijo para nombres de archivos de salida |
| `--show` | No | False | Mostrar plots interactivamente |

*Se requiere exactamente uno de `-f` o `-d`.
**`-t` es requerido cuando se usa `-d`.

---

## Scripts de visualización - campo eléctrico y potencial

### `plot_efield.py`

Grafica el campo eléctrico desde archivos efield a lo largo de líneas o planos. Soporta comparación con teoría Debye-Hückel y visualización con vectores quiver.

```bash
./plot_efield.py -n <nstep> -L <L> --mode <line|plane> [opciones]
```

| Flag | Requerido | Default | Descripción |
|---|---|---|---|
| `-n` / `--nstep` | Sí | — | Número de paso temporal |
| `-L` / `--grid-size` | Sí | — | Tamaño de grilla (un valor L×L×L o tres NX NY NZ) |
| `--mode` | Sí | — | Modo: `line` (1D) o `plane` (2D) |
| `-d` / `--directory` | No | `.` | Directorio con archivos efield |
| `--field-type` | No | `total` | Tipo de archivo efield: `total`, `real`, `fourier` |
| `--colloid-file` | No | None | Archivo colloids-*.csv para leer posición de partícula |
| `--start` | No* | — | Punto de inicio de línea (X Y Z) |
| `--end` | No* | — | Punto final de línea (X Y Z) |
| `--num-points` | No | `200` | Puntos de interpolación a lo largo de la línea |
| `--show-nodes` | No | False | Muestra solo nodos de la red sin interpolación |
| `--plane` | No* | — | Plano de corte: `xy`, `xz`, `yz`; requerido para mode=plane |
| `--position` | No | centro | Coordenada física del plano de corte |
| `--quiver` | No | False | Superpone flechas del campo sobre el plano |
| `--quiver-step` | No | `2` | Submuestreo de flechas quiver |
| `--charge-pos` | No | None | Posición de la carga para comparación teórica (X Y Z) |
| `--charge` | No | `1.0` | Carga q |
| `--kappa` | No | `0.0` | Parámetro de Debye κ (0 = Coulomb) |
| `--epsilon` | No | `1.0` | Permitividad ε (informativo) |
| `--kt` | No | `1.0` | Energía térmica kT (informativo) |
| `-o` / `--output` | No | None (pantalla) | Archivo de salida |

---

### `plot_electric_field.py`

Grafica el campo eléctrico calculado como E = −∇ψ desde archivos psi, usando derivadas por stencil. Soporta planos 2D, líneas 1D, superficies 3D y sustracción de campo externo.

```bash
./plot_electric_field.py -f <psi_file> -s <nx> <ny> <nz> -m <mode> [opciones]
```

| Flag | Requerido | Default | Descripción |
|---|---|---|---|
| `-f` / `--file` | Sí | — | Archivo psi a leer |
| `-s` / `--size` | Sí | — | Tamaño de grilla (NX NY NZ) |
| `-m` / `--mode` | Sí | — | Modo: `plane` (2D), `line` (1D), `plane3d` (superficie 3D) |
| `-p` / `--plane` | No* | — | Plano: `xy`, `xz`, `yz`; requerido para plane/plane3d |
| `-pos` / `--position` | No | centro | Índice de posición del plano |
| `-v` / `--vectors` | No | False | Superpone vectores de campo (solo en mode=plane) |
| `--vector-stride` | No | `2` | Espaciado entre vectores |
| `--start` | No* | — | Inicio de línea (X Y Z); requerido para mode=line |
| `--end` | No* | — | Fin de línea (X Y Z); requerido para mode=line |
| `--num-points` | No | `100` | Puntos a lo largo de la línea |
| `-c` / `--component` | No | `magnitude` | Componente: `magnitude`, `x`, `y`, `z`, `psi`, `all` |
| `--colormap` | No | `viridis` | Colormap para plane3d |
| `--elevation` | No | `30` | Ángulo de elevación para vista 3D |
| `--azimuth` | No | `-60` | Ángulo azimutal para vista 3D |
| `--external-field` | No | None | Campo externo constante a sustraer (Ex Ey Ez) |
| `-o` / `--output` | No | None (pantalla) | Archivo de salida |

---

### `plot_componentes_ewald.py`

Grafica las componentes del campo eléctrico de la descomposición de Ewald (real, Fourier y total) con comparación con la teoría Debye-Hückel. Soporta filtrado por dirección y exportación CSV.

```bash
./plot_componentes_ewald.py -n <nstep> -L <L> -q <carga> -eps <epsilon> [opciones]
```

| Flag | Requerido | Default | Descripción |
|---|---|---|---|
| `-n` / `--nstep` | Sí | — | Número de paso temporal |
| `-L` / `--grid-size` | Sí | — | Tamaño de grilla (L×L×L) |
| `-q` / `--charge` | Sí | — | Carga de la partícula |
| `-eps` / `--epsilon` | Sí | — | Permitividad relativa ε_r |
| `-d` / `--directory` | No | `.` | Directorio con archivos de Ludwig |
| `-kt` / `--kt` | No | `1.0` | Energía térmica kT |
| `--kappa` | No | None | Parámetro de Debye κ (si no se da, se calcula desde `--rho-el`) |
| `--rho-el` | No | None | Densidad de electrolito para calcular κ |
| `--alpha` | No | `0.5` | Parámetro de separación de Ewald α |
| `--rc` | No | `5.0` | Radio de corte en espacio real rc |
| `--rmin` | No | `0.1` | Radio mínimo |
| `--rmax` | No | L/2 | Radio máximo |
| `--max-points` | No | `10000` | Máximo de puntos a graficar |
| `--linear` | No | False | Escala lineal (default: logarítmica) |
| `--log-log` | No | False | Escala log-log en ambos ejes |
| `--log-log-error` | No | False | Escala log en el eje Y del panel de error |
| `--xmin` | No | None | Valor mínimo del eje X |
| `--ymin` | No | None | Valor mínimo del eje Y |
| `--show-ewald-components` | No | auto | Muestra componentes real y Fourier de Ewald por separado |
| `--show-ewald-total` | No | auto | Muestra suma total numérica de Ewald (real+Fourier) |
| `--show-ewald-theory-total` | No | auto | Muestra suma total teórica de Ewald |
| `--show-dh` | No | auto | Muestra curva teórica Debye-Hückel |
| `--show-dh-periodic` | No | auto | Muestra DH con imágenes periódicas |
| `--dh-periodic-shells` | No | `1` | Número de capas de imágenes periódicas (1 = 26 vecinos a distancia ~L) |
| `--direction` | No | None | Filtra por dirección: `x`, `y`, `z`, `xy`, `xz`, `yz`, `xyz`, o `dx,dy,dz` |
| `--angle-tol` | No | `15.0` | Tolerancia angular en grados para filtrado por dirección |
| `--pointwise-error` | No | False | Error por nodo vs DH periódico en posición 3D exacta; requiere `--show-dh-periodic` y `--direction` |
| `--no-error-inf` | No | False | Oculta error vs DH infinito (serie morada) |
| `--no-error-dir` | No | False | Oculta error vs DH periódico direccional (serie naranja) |
| `-o` / `--output` | No | `ewald_components.png` | Archivo de salida |
| `--csv` | No | None | Exporta datos a CSV (nombre base; se agrega `-n_NNNNNN.csv`) |
| `--csv-points` | No | None | Exporta puntos graficados con errores a CSV |

---

### `plot_psi_ewald.py`

Grafica el potencial eléctrico ψ calculado via suma de Ewald comparado con teoría Debye-Hückel. Soporta filtrado por dirección, imágenes periódicas y exportación CSV.

```bash
./plot_psi_ewald.py -n <nstep> -L <L> [opciones]
```

| Flag | Requerido | Default | Descripción |
|---|---|---|---|
| `-n` / `--nstep` | Sí | — | Número de paso temporal |
| `-L` / `--grid-size` | Sí | — | Tamaño de grilla (L×L×L) |
| `-q` / `--charge` | No | `1.0` | Carga de la partícula |
| `-kt` / `--kt` | No | `1.0` | Energía térmica kT |
| `-eps` / `--epsilon` | No | `1.0` | Permitividad relativa |
| `--alpha` | No | `0.5` | Parámetro de separación de Ewald α |
| `--rc` | No | `5.0` | Radio de corte en espacio real |
| `--kappa` | No | None | Parámetro de Debye κ (se calcula desde `--rho-el` si no se da) |
| `--rho-el` | No | None | Densidad de electrolito para calcular κ |
| `-d` / `--directory` | No | `.` | Directorio con archivos de Ludwig |
| `--rmin` | No | `0.5` | Radio mínimo a graficar |
| `--rmax` | No | L/2 | Radio máximo |
| `--show-components` | No | False | Muestra psi_real y psi_fourier como series separadas |
| `--no-total` | No | False | Oculta el psi total de Ludwig |
| `--show-ewald-theory` | No | auto | Muestra curva total de teoría Ewald (erfc+erf)/r |
| `--show-dh` | No | auto | Muestra curva DH infinita |
| `--show-dh-periodic` | No | auto | Muestra curva DH periódica |
| `--dh-periodic-shells` | No | `1` | Número de capas de imágenes periódicas |
| `--pointwise-error` | No | False | Calcula error en posición exacta del nodo (serie verde) |
| `--no-error-inf` | No | False | Oculta error vs DH infinito (serie morada) |
| `--no-error-dir` | No | False | Oculta error vs DH periódico direccional (serie naranja) |
| `--linear` | No | False | Escala lineal Y (default: semilog) |
| `--log-log` | No | False | Escala log-log en ambos ejes |
| `--log-log-error` | No | False | Escala log Y en el panel de error |
| `--xmin` | No | None | Valor mínimo del eje X |
| `--ymin` | No | None | Valor mínimo del eje Y |
| `--subtract-background` | No | False | Sustrae potencial de fondo de contraiones (promedio en nodos lejanos) |
| `--bg-rmin` | No | `12.0` | Radio mínimo para calcular el fondo |
| `--rise-to-DHp` | No | False | Desplaza todos los valores de ψ para alinear el mínimo con la curva DH periódica |
| `--direction` | No | None | Filtra por dirección: `x`, `y`, `z`, `xy`, `xyz`, o `dx,dy,dz` |
| `--angle-tol` | No | `15.0` | Tolerancia angular en grados |
| `--output` | No | `psi_ewald` | Nombre base del archivo de salida |
| `--max-points` | No | `10000` | Máximo de puntos a graficar |
| `--csv` | No | None | Exporta datos a CSV (nombre base; se agrega `-n_NNNNNN.csv`) |

---

## Scripts de visualización - cargas y estructura

### `plot_charge_distribution.py`

Visualización completa de distribución de carga. Lee archivos qsi (densidades de carga positiva/negativa) y genera plots de carga neta, especies individuales y densidad total en distintos modos.

```bash
./plot_charge_distribution.py -f <qsi_file> -s <nx> <ny> <nz> -m <mode> [opciones]
```

| Flag | Requerido | Default | Descripción |
|---|---|---|---|
| `-f` / `--file` | Sí | — | Archivo qsi a leer |
| `-s` / `--size` | Sí | — | Tamaño de grilla (NX NY NZ) |
| `-m` / `--mode` | Sí | — | Modo: `plane` (2D), `line` (1D), `plane3d` (superficie 3D), `stats` (solo estadísticas) |
| `-p` / `--plane` | No* | — | Plano: `xy`, `xz`, `yz`; requerido para plane/plane3d |
| `-pos` / `--position` | No | centro | Índice de posición del plano |
| `--start` | No* | — | Punto de inicio de línea (X Y Z); requerido para mode=line |
| `--end` | No* | — | Punto final de línea (X Y Z); requerido para mode=line |
| `--num-points` | No | `100` | Número de puntos a lo largo de la línea |
| `-c` / `--component` | No | `net` | Componente: `net`, `species0`, `species1`, `total`, `all` |
| `--colormap` | No | `RdBu_r` | Colormap para plane3d |
| `--elevation` | No | `30` | Ángulo de elevación para vista 3D |
| `--azimuth` | No | `-60` | Ángulo azimutal para vista 3D |
| `--no-interp` | No | False | Usa pcolormesh (sin interpolación) en lugar de contourf |
| `--num-species` | No | `2` | Número de especies iónicas |
| `-o` / `--output` | No | None (pantalla) | Archivo de salida |

---

### `plot_qsi_ewald.py`

Grafica la distribución de carga iónica ρ_net calculada via Ewald, con comparación con teoría Debye-Hückel. Soporta promedios radiales, mapas angulares y mapas de plano.

```bash
./plot_qsi_ewald.py -n <nstep> -L <L> [opciones]
```

| Flag | Requerido | Default | Descripción |
|---|---|---|---|
| `-n` / `--nstep` | Sí | — | Número de paso temporal |
| `-L` / `--grid-size` | Sí | — | Tamaño de grilla (L×L×L) |
| `-q` / `--charge` | No | `1.0` | Carga de la partícula |
| `-kt` / `--kt` | No | `1.0` | Energía térmica kT |
| `-eps` / `--epsilon` | No | `1.0` | Permitividad relativa |
| `--kappa` | No | None | Parámetro de Debye κ (se calcula desde `--rho-el` si no se da) |
| `--rho-el` | No | None | Densidad de electrolito para calcular κ |
| `-d` / `--directory` | No | `.` | Directorio con archivos de Ludwig |
| `--rmin` | No | `0.5` | Radio mínimo a graficar |
| `--rmax` | No | L/2 | Radio máximo |
| `--particle-radius` | No | `0.0` | Radio de partícula (lu); si >0 agrega curva DH de radio finito |
| `--num-species` | No | `2` | Número de especies iónicas en archivo qsi |
| `--show-species0` | No | False | Muestra densidad de la especie 0 (cationes) |
| `--show-species1` | No | False | Muestra densidad de la especie 1 (aniones) |
| `--show-net` | No | auto | Muestra carga neta ρ_net = ρ₊ − ρ₋ |
| `--no-net` | No | False | Oculta la carga neta |
| `--show-dh` | No | auto | Muestra curva DH infinita para ρ_net |
| `--show-dh-periodic` | No | auto | Muestra curva DH periódica |
| `--show-dh-periodic-pts` | No | False | Muestra DH periódico evaluado en posición 3D exacta del nodo |
| `--dh-periodic-shells` | No | `1` | Número de capas de imágenes periódicas |
| `--pointwise-error` | No | False | Error por nodo en posición exacta |
| `--no-error-inf` | No | False | Oculta error vs DH infinito |
| `--no-error-dir` | No | False | Oculta error vs DH periódico direccional |
| `--direction` | No | None | Filtra por dirección: `x`, `y`, `z`, `xy`, `xyz`, o `dx,dy,dz` |
| `--angle-tol` | No | `15.0` | Tolerancia angular en grados |
| `--semilog` | No | False | Escala semilog Y |
| `--log-log` | No | False | Escala log-log en ambos ejes |
| `--log-log-error` | No | False | Escala log Y en panel de error |
| `--show-error` | No | False | Muestra panel de error relativo vs DH infinito |
| `--error-ymax` | No | None | Límite superior para el eje Y de error en % |
| `--xmin` | No | None | Valor mínimo del eje X |
| `--ymin` | No | None | Valor mínimo del eje Y |
| `--output` | No | `qsi_ewald` | Nombre base del archivo de salida |
| `--max-points` | No | `10000` | Máximo de puntos a graficar |
| `--csv` | No | None | Exporta datos a CSV |
| **Promedio radial** | | | |
| `--radial` | No | False | Grafica promedio radial (media ± std) en capas esféricas |
| `--radial-scatter` | No | False | Grafica puntos individuales coloreados por capa |
| `--shell-width` | No | `1.0` | Ancho de capa para promedio radial (lu) |
| `--csv-radial` | No | None | Exporta promedio radial a CSV |
| `--csv-scatter` | No | None | Exporta nodos de scatter a CSV |
| **Mapa angular** | | | |
| `--angular-map` | No | False | Genera mapa de calor (θ,φ) en una capa esférica |
| `--angular-r` | No | None | Radio central de la capa para el mapa angular (lu) |
| `--angular-dr` | No | `1.0` | Grosor de la capa angular (lu) |
| `--angular-n-theta` | No | `90` | Resolución en θ |
| `--angular-n-phi` | No | `180` | Resolución en φ |
| `--angular-abs` | No | False | Usa \|ρ_net\| en lugar de ρ_net en el mapa angular |
| **Mapa de plano** | | | |
| `--plane-map` | No | False | Genera mapa 2D de densidad de carga en un plano |
| `--plane` | No | `xy` | Normal del plano: `xy`, `xz`, `yz`, o `nx,ny,nz` |
| `--plane-dz` | No | `0.5` | Semiancho del corte en la dirección normal (lu) |
| `--plane-type` | No | `pcolor` | Tipo de plot: `pcolor` o `contour` |
| `--plane-contour-levels` | No | `20` | Número de niveles para contourf |
| `--plane-abs` | No | False | Usa \|ρ_net\| en el mapa de plano |

---

### `plot_qsi_subgrid.py`

Compara la distribución de carga subgrid con la teoría Debye-Hückel. Grafica nodos de malla y subpuntos coloreados por distancia al voxel padre.

```bash
./plot_qsi_subgrid.py -n <nstep> -L <L> --rcut <rcut> [opciones]
```

| Flag | Requerido | Default | Descripción |
|---|---|---|---|
| `-n` / `--nstep` | Sí | — | Número de paso temporal |
| `-L` / `--grid-size` | Sí | — | Tamaño de grilla (L×L×L) |
| `--rcut` | Sí | — | Radio de corte para voxeles |
| `-q` / `--charge` | No | `1.0` | Carga de la partícula |
| `-eps` / `--epsilon` | No | `1.0` | Permitividad |
| `--kappa` | No | None | Parámetro de Debye κ (se calcula desde `--rho-el` si no se da) |
| `--rho-el` | No | None | Densidad de electrolito para calcular κ |
| `-kt` / `--kt` | No | `1.0` | Energía térmica kT |
| `--num-species` | No | `2` | Número de especies iónicas |
| `-d` / `--directory` | No | `.` | Directorio con archivos de Ludwig |
| `--qsi-file` | No | None | Ruta explícita al archivo qsi |
| `--cds-file` | No | None | Ruta explícita al archivo config.cds |
| `--rmin` | No | `0.5` | Radio mínimo para nodos |
| `--rmax` | No | =rcut | Radio máximo para el plot |
| `--rmin-sub` | No | `0.0` | Radio de exclusión mínima para subpuntos |
| `--rcut-sub` | No | False | Aplica rcut también a subpuntos individuales |
| `--nsub` | No | `3` | Subdivisiones por arista (nsub³ subpuntos por voxel) |
| `--gradient-dist` | No | False | Distribuye carga por gradiente lineal de ρ_net |
| `--quadratic-dist` | No | False | Interpolación de Lagrange cuadrática separable (3 nodos/eje) |
| `--quartic-dist` | No | False | Interpolación de Lagrange cuártica separable (5 nodos/eje) |
| `--quadratic-full` | No | False | Producto tensorial cuadrático completo (3×3×3 = 27 nodos) |
| `--quartic-full` | No | False | Producto tensorial cuártico completo (5×5×5 = 125 nodos) |
| `--no-nodes` | No | False | No muestra nodos qsi |
| `--no-subs` | No | False | No muestra subpuntos |
| `--no-dh` | No | False | No muestra la curva teórica DH |
| `--log-y` | No | False | Escala logarítmica en Y (usa \|ρ\|) |
| `--output` | No | `qsi_subgrid` | Nombre base del archivo de salida |

---

### `plot_campo_electrico_rhoel.py`

Compara magnitudes y componentes medias del campo eléctrico entre simulaciones variando el parámetro `rhoel`. Lee archivos `campo_electrico_medio.csv` en subdirectorios bajo `../`.

```bash
./plot_campo_electrico_rhoel.py [linear|log] [symlog_threshold]
```

| Argumento posicional | Default | Descripción |
|---|---|---|
| `[scale_mode]` | `log` | Escala: `linear`/`lin`/`l` para lineal; cualquier otro valor para logarítmica |
| `[symlog_threshold]` | `1e-20` | Umbral symlog (float) para plots de componentes |

---

### `plot_mui_plane.py`

Grafica planos del potencial químico μ centrados en la posición de la partícula, con colormap divergente para valores positivos/negativos.

```bash
# Step único:
./plot_mui_plane.py -n <nstep> -L <L> [opciones]

# Serie temporal:
./plot_mui_plane.py --n-start <inicio> --n-end <fin> --n-step <paso> -L <L> [opciones]
```

| Flag | Requerido | Default | Descripción |
|---|---|---|---|
| `-n` / `--nstep` | No* | — | Step único a graficar |
| `--n-start` | No* | — | Step inicial de la serie |
| `--n-end` | No* | — | Step final de la serie |
| `--n-step` | No* | — | Incremento entre steps de la serie |
| `-L` / `--grid-size` | Sí | — | Tamaño de grilla cúbica L (L×L×L) |
| `-d` / `--directory` | No | `.` | Directorio con archivos de Ludwig |
| `--normal` | No | `z` | Vector normal del plano: `x`, `y`, `z`, o `nx,ny,nz`; puede repetirse |
| `--all-planes` | No | False | Genera los tres planos ortogonales (x, y, z) |
| `--output-prefix` | No | `mui_plane` | Prefijo para nombres de archivos de salida |

*Se requiere `-n` solo, o los tres de `--n-start`, `--n-end`, `--n-step`.

---

## Scripts de organización y limpieza

### `coloideacsv.sh`

Convierte en lote archivos de configuración `config.cds-NNNNNNNNN.001-001` al formato CSV. Llama a `./scripts/extract_colloids` por cada archivo en el rango.

```bash
./coloideacsv.sh -n <paso_final> -i <paso_inicial> -p <incremento> [-o <dir_salida>]
```

| Flag | Requerido | Default | Descripción |
|---|---|---|---|
| `-n` | Sí | — | Paso final (Ntotal) |
| `-i` | Sí | — | Paso inicial (Ninicio) |
| `-p` | Sí | — | Incremento de paso |
| `-o` | No | `colloid_data` | Directorio de salida |

---

### `delete_files.sh`

Elimina archivos intermedios de simulación conservando solo el último de cada tipo por directorio (para permitir reinicios). También elimina archivos `coll*.vtk` sin conservar ninguno.

```bash
./delete_files.sh
```

Sin argumentos. Ejecutar desde el directorio padre. Archivos que limpia (conserva el último): `dist-0*`, `rho-0*`, `psi-0*`, `qsi-0*`, `vel-0*`.

---

### `delete_files_all_leave_image.sh`

Similar a `delete_files.sh` pero también elimina todos los archivos `config.cds0*` salvo el último. Más agresivo que `delete_files.sh`.

```bash
./delete_files_all_leave_image.sh
```

Sin argumentos.

---

### `delete_files_dist_rho.sh`

Elimina archivos `dist-0*` y `rho-0*` conservando solo el último de cada tipo. Usa delimitadores nulos para manejo seguro de nombres de archivos.

```bash
./delete_files_dist_rho.sh
```

Sin argumentos.

---

### `organize_by_parameter.sh`

Organiza carpetas de simulación en subdirectorios agrupados por valor de parámetro extraído del nombre (ej: `rhoel`, `kappa`). **Mueve** las carpetas.

```bash
./organize_by_parameter.sh <directorio_base> <parametro> [prefijo]
```

| Posición | Requerido | Default | Descripción |
|---|---|---|---|
| `$1` (BASE_DIR) | Sí | — | Directorio base |
| `$2` (PARAM) | Sí | — | Nombre del parámetro a buscar en nombres de directorio (ej: `rhoel`) |
| `$3` (PREFIX) | No | vacío | Prefijo opcional para nombres de subdirectorios destino |

**Ejemplo:**
```bash
./organize_by_parameter.sh /path/to/sims rhoel PB_Sn
```

---

### `organize_by_parameter_dryrun.sh`

Versión de previsualización de `organize_by_parameter.sh`: muestra qué se movería sin ejecutar cambios.

```bash
./organize_by_parameter_dryrun.sh <directorio_base> <parametro> [prefijo]
```

Misma firma que `organize_by_parameter.sh`. Solo imprime las acciones, no las ejecuta.

---

### `copy_plots.py`

Copia recursivamente subdirectorios `plots/` desde origen a destino, excluyendo archivos `.py` y `.csv`.

```bash
./copy_plots.py -origen <dir_origen> -destino <dir_destino>
```

| Flag | Requerido | Descripción |
|---|---|---|
| `-origen` | Sí | Directorio fuente |
| `-destino` | Sí | Directorio destino |

---

### `graficos.sh`

Recolecta todos los archivos `velocidades.png` de subdirectorios y los copia al directorio `../graficos/` con el nombre del primer subdirectorio como prefijo.

```bash
./graficos.sh
```

Sin argumentos. Ejecutar desde el directorio padre de las simulaciones.

---

## Scripts de orquestación

### `make_plots.sh`

Orquesta la generación de múltiples plots comparando campo eléctrico teórico con simulaciones. Extrae parámetros de nombres de directorio y llama a `analizar_campo_electrico.py`, `plot_field_comparizon.sh`, `plot_verificar_simetria_radial.sh`, `plot_velocity_field.sh` y `plot_campo_electrico_rhoel.py`.

```bash
./make_plots.sh [-d <dim>] [-e <epsilon>] [-l <size>] [-o <outfile>] [-p <paso>] [-t <kT>] \
    [-x <Lx>] [-y <Ly>] [-z <Lz>]
```

| Flag | Requerido | Default | Descripción |
|---|---|---|---|
| `-d` | No | auto-detectado | Dimensión del sistema (2 o 3); auto-detectada desde la grilla |
| `-e` | No | desde nombre de dir | epsilon × kT electrostático |
| `-l` | No | max(Lx,Ly,Lz) | Tamaño del sistema L |
| `-o` | No | `single-peskin-Efield` | Nombre base del archivo de salida |
| `-p` | No | — | Incremento de pasos para salida (paso) |
| `-t` | No | desde nombre de dir | Energía térmica kT |
| `-x` | No | desde nombre de dir | Tamaño de grilla X (Lx) |
| `-y` | No | desde nombre de dir | Tamaño de grilla Y (Ly) |
| `-z` | No | desde nombre de dir | Tamaño de grilla Z (Lz) |

Los parámetros `epsilon`, `kappa`, `kT` y `rhoel` se extraen también de los nombres de subdirectorios via regex.

---

### `plot_field_comparizon.sh`

Compara el campo eléctrico teórico Debye-Hückel con la simulación a lo largo de líneas y diagonales. Llama a `compare_field_theory.py` por cada subdirectorio con patrón `pos_X_Y_Z`.

```bash
./plot_field_comparizon.sh -e <epsilon> -k <kappa> -p <paso> [opciones]
```

| Flag | Requerido | Default | Descripción |
|---|---|---|---|
| `-e` | Sí | — | epsilon × kT electrostático |
| `-k` | Sí | — | Longitud inversa de Debye (kappa) |
| `-p` | Sí | — | Incremento de pasos para salida |
| `-o` | No | — | Nombre base del archivo de salida |
| `-s` | No | — | Posición de inicio como `X_Y_Z` |
| `-f` | No | — | Posición final como `X_Y_Z` |
| `-t` | No | — | Energía térmica kT |
| `-x` | No | auto (32) | Tamaño de grilla X |
| `-y` | No | auto | Tamaño de grilla Y |
| `-z` | No | auto | Tamaño de grilla Z |

---

### `plot_velocity_field.sh`

Wrapper que aplica `plot_velocity_field.py` a todos los subdirectorios bajo `./*/`.

```bash
./plot_velocity_field.sh -p <paso> -x <Lx> -y <Ly> -z <Lz>
```

| Flag | Requerido | Descripción |
|---|---|---|
| `-p` | Sí | Paso temporal |
| `-x` | Sí | Tamaño de grilla X |
| `-y` | Sí | Tamaño de grilla Y |
| `-z` | Sí | Tamaño de grilla Z |

---

### `plot_verificar_simetria_radial.sh`

Ejecuta `verificar_simetria_psi.py` en múltiples subdirectorios (dos veces: escala lineal y log). Llama también a `plot_field_comparizon.sh`.

```bash
./plot_verificar_simetria_radial.sh -e <epsilon> -k <kappa> -p <paso> [opciones]
```

| Flag | Requerido | Default | Descripción |
|---|---|---|---|
| `-e` | Sí | — | Permitividad electrostática |
| `-k` | Sí | — | Longitud inversa de Debye (kappa) |
| `-p` | Sí | — | Incremento de pasos para salida |
| `-d` | No | — | Dimensión del sistema (2 o 3); afecta Lz (Lz=4 si dim=2) |
| `-l` | No | — | Tamaño del sistema L |
| `-o` | No | — | Nombre base del archivo de salida |
| `-q` | No | — | Flag para plot logarítmico |
| `-r` | No | — | Flag radial |
| `-s` | No | — | Posición de inicio |
| `-f` | No | — | Posición final |
| `-t` | No | — | Energía térmica kT |
| `-a` | No | — | Grosor de capa para análisis radial logarítmico |

---

## Referencia de archivos de entrada de Ludwig

| Archivo | Contenido |
|---|---|
| `colloids-NNNNNNNN.csv` | Posiciones, velocidades y propiedades de partículas por paso temporal |
| `vel-NNNNNNNNN.001-001` | Campo de velocidades del fluido en la grilla |
| `psi-NNNNNNNN.001-001` | Potencial eléctrico ψ en la grilla |
| `qsi-NNNNNNNN.001-001` | Densidades de carga iónica por especie en la grilla |
| `efield-NNNNNNNN.001-001` | Campo eléctrico en la grilla |
| `dist-NNNNNNNN.001-001` | Distribuciones de población LB |
| `rho-NNNNNNNN.001-001` | Densidad del fluido |
| `config.cds-NNNNNNNNN.001-001` | Configuración de coloides (formato binario) |
| `theta-NNNNNNNN.001-001` | θ_k = ρ_k·exp(z_k·ϕ) por especie (generado por `calc_mui.py`) |
