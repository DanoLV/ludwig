# Visualización de Distribución de Carga desde Archivos QSI de Ludwig

Este script permite analizar y visualizar la distribución de carga eléctrica a partir de los archivos de densidad de carga (qsi) generados por Ludwig.

## Descripción

Los archivos qsi contienen la densidad de carga de especies iónicas en cada punto de la malla. Este script proporciona:

- **Visualizaciones 2D y 3D** de la distribución de carga
- **Análisis de línea 1D** a lo largo de trayectorias específicas
- **Cálculo de carga total** del sistema
- **Estadísticas completas** por especie y carga neta

## Archivos

1. **`plot_charge_distribution.py`**: Script principal para graficar distribución de carga

## Requisitos

```bash
pip install numpy matplotlib scipy
```

## Conceptos

### Especies Iónicas

Los archivos qsi contienen 2 columnas con **densidades de carga** (siempre ≥ 0):
- **Columna 0**: Densidad de carga positiva ρ₊ (cationes, con carga +1)
- **Columna 1**: Densidad de carga negativa ρ₋ (aniones, con carga -1)

**IMPORTANTE**: Las densidades representan concentraciones y son **siempre no-negativas**.
Valores negativos indicarían un error en la simulación o en la lectura del archivo.

<!-- CHANGE INIT - Subgrid charge output file -->
### Carga de Coloides Subgrid

Los archivos **qsi_colloid** contienen 1 columna con **carga neta de coloides subgrid**:
- **Columna 0**: Carga neta interpolada de coloides ρ_colloid = q₀ - q₁
  - q₀: densidad positiva del coloide
  - q₁: densidad negativa del coloide
  - La carga se distribuye a la malla usando la función delta de Peskin

**NOTA**: Los valores pueden ser positivos, negativos o cero, dependiendo de la carga neta del coloide.
Este archivo se genera automáticamente junto con psi y qsi cuando hay coloides subgrid presentes.
<!-- CHANGE END - Subgrid charge output file -->

### Cantidades Derivadas

- **Carga neta en cada nodo**: ρ_net = ρ₊ - ρ₋
  - Positivo: exceso de cationes
  - Negativo: exceso de aniones
  - Cero: localmente neutro
- **Densidad iónica total**: ρ_total = ρ₊ + ρ₋ (concentración total de iones)
- **Carga total del sistema**:
  - Carga positiva total: ∑ρ₊ (sobre todo el dominio)
  - Carga negativa total: ∑ρ₋ (sobre todo el dominio)
  - Carga neta total: ∑(ρ₊ - ρ₋) ≈ 0 (electroneutralidad esperada)

## Uso Básico

### Sintaxis General

```bash
./plot_charge_distribution.py -f <archivo_qsi> -s <nx> <ny> <nz> -m <modo> [opciones]
```

**Argumentos obligatorios:**
- `-f, --file`: Archivo qsi a leer
- `-s, --size`: Tamaño de la malla (nx ny nz)
- `-m, --mode`: Modo de graficación (`plane`, `line`, `plane3d`, o `stats`)

### Modo Estadísticas (`stats`)

Calcula y muestra estadísticas de carga sin generar gráficos.

**Ejemplo:**

```bash
./plot_charge_distribution.py -f qsi-000010250.001-001 -s 32 32 32 -m stats
```

**Salida:**
```
======================================================================
ESTADÍSTICAS DE DISTRIBUCIÓN DE CARGA
======================================================================
Tamaño de malla: 32 × 32 × 32 = 32768 puntos

Densidad de Carga Positiva ρ₊ (Cationes, carga +1):
  Carga/Densidad total:    1.000000e+00    ← Carga positiva total
  Promedio:                3.051758e-05
  Desv. std:               2.087755e-03
  Mínimo:                  0.000000e+00    ← Debe ser ≥ 0
  Máximo:                  1.428571e-01

Densidad de Carga Negativa ρ₋ (Aniones, carga -1):
  Carga/Densidad total:    1.000000e+00    ← Carga negativa total
  Promedio:                3.051758e-05
  Desv. std:               4.460969e-07
  Mínimo:                  0.000000e+00    ← Debe ser ≥ 0
  Máximo:                  3.058189e-05

Carga Neta (ρ₊ - ρ₋):
  Carga/Densidad total:    1.665335e-16
  Promedio:                5.082198e-21
  Desv. std:               2.088202e-03
  Mínimo:                 -3.058189e-05
  Máximo:                  1.428571e-01

✓ Sistema electroneutro (carga neta total ≈ 0)    ← Verificación
======================================================================
```

**Interpretación:**
- Las dos especies tienen la **misma carga total** (1.0) → sistema balanceado
- Ambas densidades tienen **mínimo = 0** → correcta (no hay valores negativos)
- Carga neta total ≈ **10⁻¹⁶** → electroneutralidad verificada
- ✓ marca indica que el sistema está correctamente balanceado

### Modo Plano 2D (`plane`)

Visualiza la distribución de carga en un plano 2D cortando el dominio 3D.

**Opciones:**
- `-p, --plane`: Plano a graficar (`xy`, `xz`, o `yz`)
- `-pos, --position`: Posición del plano (índice). Por defecto: centro
- `-c, --component`: Componente a graficar:
  - `net`: Carga neta ρ₀ - ρ₁ (por defecto)
  - `species0`: Densidad de especie 0
  - `species1`: Densidad de especie 1
  - `total`: Densidad total |ρ₀| + |ρ₁|
- `-o, --output`: Archivo de salida PNG

**Ejemplos:**

```bash
# 1. Carga neta en plano XY central
./plot_charge_distribution.py -f qsi-000010250.001-001 -s 32 32 32 -m plane -p xy -c net

# 2. Especie 0 (cationes) en plano XZ en y=10
./plot_charge_distribution.py -f qsi-000010250.001-001 -s 32 32 32 -m plane -p xz -pos 10 -c species0

# 3. Especie 1 (aniones) en plano YZ
./plot_charge_distribution.py -f qsi-000010250.001-001 -s 32 32 32 -m plane -p yz -c species1

# 4. Densidad total en plano XY y guardar
./plot_charge_distribution.py -f qsi-000010250.001-001 -s 32 32 32 -m plane -p xy -c total -o carga_total.png
```

### Modo Superficie 3D (`plane3d`)

Visualiza la distribución de carga como una superficie tridimensional.

**Opciones:**
- `-p, --plane`: Plano a graficar (`xy`, `xz`, o `yz`)
- `-pos, --position`: Posición del plano (índice). Por defecto: centro
- `-c, --component`: Componente a graficar (`net`, `species0`, `species1`, `total`)
- `--colormap`: Mapa de colores (por defecto: `RdBu_r`)
  - `RdBu_r`: Rojo-azul (bueno para carga neta)
  - `viridis`, `plasma`: Buenos para densidades
  - `seismic`: Destacar regiones positivas/negativas
- `--elevation`: Ángulo de elevación en grados (por defecto: 30)
- `--azimuth`: Ángulo azimutal en grados (por defecto: -60)
- `-o, --output`: Archivo de salida PNG

**Ejemplos:**

```bash
# 1. Superficie 3D de carga neta en plano XY
./plot_charge_distribution.py -f qsi-000010250.001-001 -s 32 32 32 -m plane3d -p xy -c net

# 2. Superficie 3D de especie 0 con colormap viridis
./plot_charge_distribution.py -f qsi-000010250.001-001 -s 32 32 32 -m plane3d -p xz -c species0 --colormap viridis

# 3. Densidad total 3D con vista personalizada
./plot_charge_distribution.py -f qsi-000010250.001-001 -s 32 32 32 -m plane3d -p xy -c total --elevation 45 --azimuth -30

# 4. Vista superior de carga neta
./plot_charge_distribution.py -f qsi-000010250.001-001 -s 32 32 32 -m plane3d -p xy -c net --elevation 90 --azimuth 0
```

### Modo Línea 1D (`line`)

Visualiza la distribución de carga a lo largo de una línea en el espacio.

**Opciones:**
- `--start`: Punto inicial (x y z) en coordenadas de la malla
- `--end`: Punto final (x y z) en coordenadas de la malla
- `--num-points`: Número de puntos interpolados (por defecto: 100)
- `-c, --component`: Componente a graficar:
  - `all`: Todas las especies + carga neta en subgráficos separados (recomendado)
  - `net`: Solo carga neta
  - `species0`: Solo especie 0
  - `species1`: Solo especie 1
  - `total`: Solo densidad total
- `-o, --output`: Archivo de salida PNG

**Ejemplos:**

```bash
# 1. Todas las especies a lo largo del eje X
./plot_charge_distribution.py -f qsi-000010250.001-001 -s 32 32 32 -m line \
  --start 0 16 16 --end 31 16 16 -c all

# 2. Solo carga neta a lo largo de la diagonal
./plot_charge_distribution.py -f qsi-000010250.001-001 -s 32 32 32 -m line \
  --start 0 0 0 --end 31 31 31 -c net

# 3. Especie 0 a lo largo del eje Z
./plot_charge_distribution.py -f qsi-000010250.001-001 -s 32 32 32 -m line \
  --start 16 16 0 --end 16 16 31 -c species0

# 4. Densidad total a lo largo del eje Y con más puntos
./plot_charge_distribution.py -f qsi-000010250.001-001 -s 32 32 32 -m line \
  --start 16 0 16 --end 16 31 16 -c total --num-points 200
```

## Interpretación de los Gráficos

### Gráficos de Plano 2D

- **Mapa de colores**:
  - Para `net`: Rojo = carga positiva, Azul = carga negativa, Blanco = neutra
  - Para `species0/1`: Intensidad del color indica concentración de iones
  - Para `total`: Indica regiones de alta densidad iónica total

### Gráficos de Superficie 3D

- **Altura (eje Z)**: Representa el valor de la densidad de carga
- **Color**: Proporciona información adicional sobre la magnitud
- **Picos**: Regiones de alta concentración de carga
- **Valles**: Regiones de baja concentración

### Gráficos de Línea 1D

**Modo `all`**: Genera dos subgráficos
1. **Densidades por Especie**: Muestra cada especie con diferentes colores
2. **Carga Neta**: Muestra la diferencia entre especies

**Modos individuales**: Muestra solo la componente seleccionada

## Interpretación Física

### Densidades y Carga

**Densidades (ρ₊ y ρ₋):**
- Representan **concentraciones** de iones
- Son **siempre no-negativas** (≥ 0)
- Unidades típicas: carga por unidad de volumen
- Valores negativos indicarían error

**Carga neta en un nodo:**
```
ρ_net = ρ₊ - ρ₋

Si ρ_net > 0: exceso de cationes (región positiva)
Si ρ_net < 0: exceso de aniones (región negativa)
Si ρ_net = 0: localmente neutro
```

### Carga Total del Sistema

- **∑ρ₊**: Carga positiva total (suma de todos los cationes)
- **∑ρ₋**: Carga negativa total (suma de todos los aniones)
- **∑(ρ₊ - ρ₋)**: Carga neta total del sistema

**Electroneutralidad:**
- En un sistema cerrado: ∑ρ₊ ≈ ∑ρ₋
- Por tanto: ∑(ρ₊ - ρ₋) ≈ 0
- Desviaciones pequeñas (~10⁻¹⁶) son errores numéricos normales
- Desviaciones grandes pueden indicar:
  - Errores en la simulación
  - Condiciones de frontera no neutras
  - Partículas cargadas en el sistema

### Distribución Espacial

- **Regiones de alta densidad**: Pueden indicar:
  - Presencia de partículas cargadas (coloides)
  - Acumulación en capas dobles eléctricas
  - Efectos de confinamiento

- **Gradientes de carga**: Generan campos eléctricos (ver `plot_electric_field.py`)

## Casos de Uso Comunes

### 1. Verificar Electroneutralidad

```bash
./plot_charge_distribution.py -f qsi-000010250.001-001 -s 32 32 32 -m stats
```
Verificar que "Carga neta total" ≈ 0

### 2. Visualizar Capa Doble Eléctrica

```bash
# Vista 3D de carga neta alrededor de partícula
./plot_charge_distribution.py -f qsi-000010250.001-001 -s 32 32 32 -m plane3d -p xy -c net --colormap seismic

# Perfil radial de carga
./plot_charge_distribution.py -f qsi-000010250.001-001 -s 32 32 32 -m line \
  --start 16 16 16 --end 31 16 16 -c all
```

### 3. Comparar Distribución de Cationes y Aniones

```bash
# Cationes
./plot_charge_distribution.py -f qsi-000010250.001-001 -s 32 32 32 -m plane -p xy -c species0 -o cationes.png

# Aniones
./plot_charge_distribution.py -f qsi-000010250.001-001 -s 32 32 32 -m plane -p xy -c species1 -o aniones.png
```

### 4. Analizar Variación Temporal

Procesar múltiples timesteps y comparar la carga total:

```bash
for file in qsi-*.001-001; do
    echo "Archivo: $file"
    ./plot_charge_distribution.py -f $file -s 32 32 32 -m stats | grep "Carga total"
    echo ""
done
```

## Relación con Otros Scripts

- **`plot_electric_field.py`**: Los gradientes de carga generan campos eléctricos
  - Usar qsi para ver dónde está la carga
  - Usar psi/E para ver el campo resultante

- Secuencia de análisis recomendada:
  1. `plot_charge_distribution.py -m stats` → Verificar electroneutralidad
  2. `plot_charge_distribution.py -m plane3d` → Visualizar distribución
  3. `plot_electric_field.py -m plane3d` → Ver campo eléctrico resultante

## Consejos de Uso

1. **Identificar tamaño de malla**:
   ```bash
   grep "size" input
   ```

2. **Listar archivos qsi**:
   ```bash
   ls -lh qsi-*.001-001
   ```

3. **Verificar electroneutralidad**:
   La carga neta total debe ser muy pequeña (~10⁻¹⁶ o menor)

4. **Elegir colormap apropiado**:
   - `RdBu_r` o `seismic`: Para carga neta (distingue ±)
   - `viridis` o `plasma`: Para especies individuales
   - `coolwarm`: Alternativa para carga neta

5. **Resolución de gráficos**:
   Los archivos PNG se guardan a 300 DPI (alta calidad para publicaciones)

<!-- CHANGE INIT - Subgrid charge output file -->
## Visualización de Carga de Coloides Subgrid

Para visualizar la distribución de carga de coloides subgrid desde archivos qsi_colloid:

```bash
# Plano XY mostrando carga neta de coloides
python3 plot_charge_distribution.py -f qsi_colloid-000010000.001-001 -s 32 32 32 \
    -m plane -p xy -c species_0 -o colloid_charge_xy.png

# Línea a lo largo del eje Z
python3 plot_charge_distribution.py -f qsi_colloid-000010000.001-001 -s 32 32 32 \
    -m line --start 16 16 0 --end 16 16 32 -c species_0 -o colloid_charge_line.png

# Superficie 3D
python3 plot_charge_distribution.py -f qsi_colloid-000010000.001-001 -s 32 32 32 \
    -m plane3d -p xy --position 16 -c species_0 -o colloid_charge_3d.png
```

**Nota**: Para archivos qsi_colloid con una sola columna, usa `-c species_0` para visualizar la carga neta.
<!-- CHANGE END - Subgrid charge output file -->

## Solución de Problemas

**Error: "El archivo tiene X filas, pero se esperaban Y"**
- Verifica que el tamaño de malla especificado con `-s` sea correcto

**Carga total muy diferente de lo esperado**
- Revisa las unidades y la configuración de la simulación
- Verifica que el archivo qsi no esté corrupto

**Gráficos muestran valores muy pequeños o cero**
- Puede que no haya carga en esa región
- Prueba diferentes planos o posiciones
- Usa modo `stats` para ver el rango de valores

## Notas Técnicas

- La interpolación para líneas usa `scipy.interpolate.RegularGridInterpolator` con interpolación lineal
- Los gráficos 3D usan `plot_surface` de matplotlib con antialiasing
- Las estadísticas se calculan sobre todo el dominio 3D

## Referencias

- Documentación de Ludwig: http://ludwig.epcc.ed.ac.uk/
- Electrocinética en Ludwig: Ver documentación sobre `psi` y `electrokinetics`
- Teoría de doble capa eléctrica: DLVO theory

## Autor

Script generado para facilitar el análisis de simulaciones electrocinéticas con Ludwig.
