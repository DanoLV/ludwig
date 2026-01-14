# Graficación del Campo Eléctrico desde Archivos PSI de Ludwig

Este conjunto de scripts permite visualizar el campo eléctrico calculado a partir de los archivos de potencial eléctrico (psi) generados por Ludwig.

## Descripción

El campo eléctrico se calcula como **E = -∇ψ**, donde ψ es el potencial eléctrico. Los scripts leen los archivos psi, calculan el gradiente del potencial usando diferencias finitas, y generan visualizaciones en 2D (planos) o 1D (líneas).

## Archivos

1. **`plot_electric_field.py`**: Script principal para graficar un único archivo psi
2. **`batch_plot_electric_field.py`**: Script para procesar múltiples archivos psi de forma automática

## Requisitos

```bash
pip install numpy matplotlib scipy
```

## Uso Básico

### Script Individual: `plot_electric_field.py`

#### Sintaxis General

```bash
./plot_electric_field.py -f <archivo_psi> -s <nx> <ny> <nz> -m <modo> [opciones]
```

**Argumentos obligatorios:**
- `-f, --file`: Archivo psi a leer
- `-s, --size`: Tamaño de la malla (nx ny nz)
- `-m, --mode`: Modo de graficación (`plane`, `line`, o `plane3d`)

#### Modo Plano (`plane`)

Visualiza el campo eléctrico en un plano 2D cortando el dominio 3D.

**Opciones:**
- `-p, --plane`: Plano a graficar (`xy`, `xz`, o `yz`)
- `-pos, --position`: Posición del corte (índice). Por defecto: centro
- `-c, --component`: Componente a graficar:
  - `magnitude`: Magnitud |E| (por defecto)
  - `x`, `y`, `z`: Componentes individuales Ex, Ey, Ez
  - `psi`: Potencial eléctrico ψ
- `-v, --vectors`: Mostrar vectores del campo superpuestos
- `--vector-stride`: Espaciado entre vectores (por defecto: 2)
- `-o, --output`: Archivo de salida PNG

**Ejemplos:**

```bash
# 1. Graficar magnitud del campo en plano XY central
./plot_electric_field.py -f psi-000010050.001-001 -s 32 32 32 -m plane -p xy

# 2. Graficar componente Ez en plano XZ en z=10 con vectores
./plot_electric_field.py -f psi-000010050.001-001 -s 32 32 32 -m plane -p xz -pos 10 -c z -v

# 3. Graficar potencial eléctrico en plano YZ y guardar
./plot_electric_field.py -f psi-000010050.001-001 -s 32 32 32 -m plane -p yz -c psi -o potencial.png

# 4. Campo con vectores densos (stride=1)
./plot_electric_field.py -f psi-000010050.001-001 -s 32 32 32 -m plane -p xy -v --vector-stride 1
```

#### Modo Línea (`line`)

Visualiza el campo eléctrico y potencial a lo largo de una línea en el espacio 3D.

**Opciones:**
- `--start`: Punto inicial (x y z) en coordenadas de la malla
- `--end`: Punto final (x y z) en coordenadas de la malla
- `--num-points`: Número de puntos interpolados (por defecto: 100)
- `-c, --component`: Componente a graficar:
  - `all`: Todas las componentes (Ex, Ey, Ez, |E|) + potencial ψ en dos subgráficos (por defecto si no se especifica)
  - `magnitude`: Solo la magnitud |E|
  - `x`, `y`, `z`: Solo la componente Ex, Ey o Ez
  - `psi`: Solo el potencial eléctrico ψ
- `-o, --output`: Archivo de salida PNG

**Ejemplos:**

```bash
# 1. Todas las componentes a lo largo de la diagonal (gráfico completo con 2 subplots)
./plot_electric_field.py -f psi-000010050.001-001 -s 32 32 32 -m line --start 0 0 0 --end 31 31 31 -c all

# 2. Solo componente Ez a lo largo del eje X
./plot_electric_field.py -f psi-000010050.001-001 -s 32 32 32 -m line --start 0 16 16 --end 31 16 16 -c z

# 3. Solo magnitud del campo a lo largo del eje Z
./plot_electric_field.py -f psi-000010050.001-001 -s 32 32 32 -m line --start 16 16 0 --end 16 16 31 -c magnitude

# 4. Solo potencial eléctrico a lo largo del eje Y
./plot_electric_field.py -f psi-000010050.001-001 -s 32 32 32 -m line --start 16 0 16 --end 16 31 16 -c psi

# 5. Línea personalizada con más puntos
./plot_electric_field.py -f psi-000010050.001-001 -s 32 32 32 -m line \
  --start 5 10 8 --end 27 22 24 --num-points 200 -c all -o campo_linea.png
```

#### Modo Superficie 3D (`plane3d`)

Visualiza el potencial o campo eléctrico como una superficie tridimensional en un plano específico.

**Opciones:**
- `-p, --plane`: Plano a graficar (`xy`, `xz`, o `yz`)
- `-pos, --position`: Posición del plano (índice). Por defecto: centro
- `-c, --component`: Componente a graficar:
  - `magnitude`: Magnitud |E| (por defecto)
  - `x`, `y`, `z`: Componentes individuales Ex, Ey, Ez
  - `psi`: Potencial eléctrico ψ
- `--colormap`: Mapa de colores (por defecto: `viridis`)
  - Opciones: `viridis`, `plasma`, `inferno`, `magma`, `RdBu_r`, `coolwarm`, `seismic`, etc.
- `--elevation`: Ángulo de elevación de la vista en grados (por defecto: 30)
- `--azimuth`: Ángulo azimutal de la vista en grados (por defecto: -60)
- `-o, --output`: Archivo de salida PNG

**Ejemplos:**

```bash
# 1. Superficie 3D del potencial en plano XY central
./plot_electric_field.py -f psi-000010050.001-001 -s 32 32 32 -m plane3d -p xy -c psi

# 2. Magnitud del campo 3D con colormap plasma
./plot_electric_field.py -f psi-000010050.001-001 -s 32 32 32 -m plane3d -p xz -c magnitude --colormap plasma

# 3. Componente Ez en 3D con vista personalizada
./plot_electric_field.py -f psi-000010050.001-001 -s 32 32 32 -m plane3d -p xy -c z --elevation 45 --azimuth -30

# 4. Potencial en plano YZ con colormap RdBu_r (rojo-azul)
./plot_electric_field.py -f psi-000010050.001-001 -s 32 32 32 -m plane3d -p yz -c psi --colormap RdBu_r

# 5. Vista superior (elevation=90) del potencial
./plot_electric_field.py -f psi-000010050.001-001 -s 32 32 32 -m plane3d -p xy -c psi --elevation 90 --azimuth 0

# 6. Guardar superficie 3D en archivo
./plot_electric_field.py -f psi-000010050.001-001 -s 32 32 32 -m plane3d -p xy -c psi -o potencial_3d.png
```

### Script por Lotes: `batch_plot_electric_field.py`

Procesa múltiples archivos psi automáticamente.

#### Sintaxis General

```bash
./batch_plot_electric_field.py -d <directorio> -s <nx> <ny> <nz> -m <modo> [opciones]
```

**Argumentos obligatorios:**
- `-d, --directory`: Directorio con archivos psi (por defecto: `.`)
- `-s, --size`: Tamaño de la malla
- `-m, --mode`: Modo (`plane` o `line`)

**Opciones adicionales:**
- `--pattern`: Patrón de búsqueda (por defecto: `psi-*.001-001`)
- `--output-dir`: Directorio de salida (por defecto: `electric_field_plots`)
- `--max-files`: Procesar solo los primeros N archivos
- `--skip`: Procesar cada N archivos (útil para conjuntos grandes)

**Ejemplos:**

```bash
# 1. Procesar todos los archivos en plano XY con vectores
./batch_plot_electric_field.py -d . -s 32 32 32 -m plane -p xy -v

# 2. Procesar solo archivos con timestep 00001XXXX
./batch_plot_electric_field.py -d . -s 32 32 32 --pattern "psi-00001*.001-001" -m plane -p xz

# 3. Procesar cada 5 archivos para ahorrar tiempo
./batch_plot_electric_field.py -d . -s 32 32 32 -m plane -p xy --skip 5

# 4. Primeros 10 archivos, líneas a lo largo del eje X
./batch_plot_electric_field.py -d . -s 32 32 32 -m line \
  --start 0 16 16 --end 31 16 16 --max-files 10

# 5. Guardar en directorio personalizado
./batch_plot_electric_field.py -d . -s 32 32 32 -m plane -p yz --output-dir mis_graficos/

# 6. Generar superficies 3D del potencial para todos los archivos
./batch_plot_electric_field.py -d . -s 32 32 32 -m plane3d -p xy -c psi --output-dir graficos_3d/

# 7. Superficies 3D con colormap personalizado
./batch_plot_electric_field.py -d . -s 32 32 32 -m plane3d -p xz -c magnitude --colormap plasma --output-dir graficos_3d_plasma/
```

## Interpretación de los Gráficos

### Gráficos de Plano 2D

- **Mapa de colores**: Muestra la intensidad del componente seleccionado
  - Rojo: Valores positivos
  - Azul: Valores negativos
  - Blanco: Valores cercanos a cero

- **Vectores**: Si se activan con `-v`, muestran la dirección del campo eléctrico en el plano

- **Componentes**:
  - `magnitude`: Muestra la intensidad total del campo |E| = √(Ex² + Ey² + Ez²)
  - `x`, `y`, `z`: Muestran la componente direccional específica
  - `psi`: Muestra el potencial eléctrico (no el campo)

### Gráficos de Línea 1D

Dependiendo de la opción `-c, --component`:

**Modo `all` (por defecto)**: Genera dos subgráficos
1. **Campo Eléctrico**:
   - Líneas roja, verde, azul: Componentes Ex, Ey, Ez
   - Línea negra discontinua: Magnitud total |E|

2. **Potencial Eléctrico**:
   - Línea púrpura: Potencial ψ a lo largo de la línea

**Modos individuales** (`x`, `y`, `z`, `magnitude`, `psi`): Genera un único gráfico
- Muestra solo la componente seleccionada
- Útil para análisis enfocado de una componente específica
- Gráfico más simple y compacto

El eje X siempre muestra la distancia desde el punto inicial.

### Gráficos de Superficie 3D

Los gráficos 3D muestran el potencial o campo eléctrico como una superficie tridimensional:

- **Superficie**: La altura (eje Z) representa el valor del campo o potencial
- **Color**: El mapa de colores proporciona información adicional sobre la magnitud
- **Vista**: Los ángulos de elevación y azimut permiten ver la superficie desde diferentes perspectivas
  - `elevation`: 0° = vista lateral, 90° = vista superior
  - `azimuth`: controla la rotación alrededor del eje Z

**Mapas de colores recomendados:**
- `viridis`, `plasma`, `inferno`: Buenos para magnitudes (sin valores negativos)
- `RdBu_r`, `coolwarm`, `seismic`: Buenos para componentes con valores positivos y negativos
- `magma`: Bueno para destacar regiones de alto valor

**Ventajas del modo 3D:**
- Visualiza la "topografía" del potencial o campo
- Identifica fácilmente máximos, mínimos y gradientes
- Mejor percepción de variaciones espaciales complejas

## Estructura de Archivos PSI

Los archivos psi de Ludwig contienen:
- Una línea por cada punto de la malla
- Valores en notación científica (ej: `1.890742573739617e-03`)
- Orden de escritura: para cada (x, y), todos los valores de z

Para una malla de 32×32×32:
- Total de líneas: 32,768
- Tamaño típico: ~500-800 KB

## Consejos de Uso

1. **Identificar el tamaño de la malla**: Busca en tu archivo `input`:
   ```bash
   grep "size" input
   ```

2. **Listar archivos psi disponibles**:
   ```bash
   ls -lh psi-*.001-001
   ```

3. **Para animaciones**: Usa el script batch y luego combina las imágenes:
   ```bash
   # Procesar todos los archivos
   ./batch_plot_electric_field.py -d . -s 32 32 32 -m plane -p xy

   # Crear video con ffmpeg
   cd electric_field_plots
   ffmpeg -framerate 10 -pattern_type glob -i 'E_field_*.png' -c:v libx264 -pix_fmt yuv420p campo_electrico.mp4
   ```

4. **Comparar diferentes planos**: Ejecuta el script batch tres veces con diferentes planos:
   ```bash
   ./batch_plot_electric_field.py -d . -s 32 32 32 -m plane -p xy --output-dir graficos_xy
   ./batch_plot_electric_field.py -d . -s 32 32 32 -m plane -p xz --output-dir graficos_xz
   ./batch_plot_electric_field.py -d . -s 32 32 32 -m plane -p yz --output-dir graficos_yz
   ```

## Solución de Problemas

**Error: "El archivo tiene X valores, pero se esperaban Y"**
- Verifica que el tamaño de malla especificado con `-s` coincida con el usado en la simulación

**Error: "No se encuentra el archivo"**
- Usa rutas absolutas o asegúrate de estar en el directorio correcto
- Verifica que el archivo exista con `ls psi-*.001-001`

**Gráficos vacíos o con valores muy pequeños**
- Puede que no haya campo eléctrico significativo en esa región
- Prueba visualizar el potencial con `-c psi` para verificar

**Script batch no encuentra archivos**
- Verifica el patrón con `ls psi-*.001-001`
- Ajusta el `--pattern` si tus archivos tienen formato diferente

## Restar Campo Eléctrico Externo

### Motivación

En simulaciones con campo eléctrico aplicado externamente, el campo total es:
```
E_total = E_cargas + E_externo
```

Para visualizar **solo el campo generado por las cargas** (polarización, capa doble eléctrica, etc.), podemos restar el campo externo constante.

### Uso

Añadir el argumento `--external-field Ex Ey Ez`:

```bash
# Ejemplo: Campo externo E_ext = (0, 0, 1.0) en dirección Z
./plot_electric_field.py -f psi-000000500.001-001 -s 32 32 32 -m plane -p xy \
  --external-field 0 0 1.0
```

### Ejemplos

```bash
# 1. Restar campo externo en dirección Z
./plot_electric_field.py -f psi-000000500.001-001 -s 32 32 32 -m plane -p xy -c magnitude \
  --external-field 0 0 0.001

# 2. Restar campo externo en dirección X (electroforesis)
./plot_electric_field.py -f psi-000000500.001-001 -s 32 32 32 -m plane3d -p xz -c magnitude \
  --external-field 0.005 0 0

# 3. Comparar campo con y sin externo (línea)
# Con campo total:
./plot_electric_field.py -f psi-000000500.001-001 -s 32 32 32 -m line \
  --start 0 16 16 --end 31 16 16 -c all -o campo_total.png

# Solo campo de cargas:
./plot_electric_field.py -f psi-000000500.001-001 -s 32 32 32 -m line \
  --start 0 16 16 --end 31 16 16 -c all --external-field 0 0 1.0 -o campo_cargas.png

# 4. Visualizar polarización de partícula (restar campo uniforme)
./plot_electric_field.py -f psi-000000500.001-001 -s 32 32 32 -m plane3d -p xy \
  --external-field 0 0 1.0 -o polarizacion_3d.png
```

### Casos de Uso

**1. Análisis de Polarización:**
- Aplicar campo externo constante
- Restar ese campo para ver solo la respuesta de las cargas
- Visualizar dipolo inducido o redistribución de carga

**2. Capa Doble Eléctrica:**
- En electroforesis, restar el campo aplicado
- Ver solo el campo generado por la capa doble alrededor de la partícula

**3. Verificación Numérica:**
- Comparar campo total vs campo sin externo
- El campo restado debe mostrar claramente la contribución de las cargas

### Nota Importante

El campo externo se especifica en las **mismas unidades** que el campo calculado de ψ.
Consulta el archivo `input` de tu simulación para ver el valor del campo aplicado.

## Notas Técnicas

- El campo eléctrico se calcula usando `numpy.gradient()`, que implementa diferencias finitas centradas de segundo orden en el interior del dominio
- La resta del campo externo se hace componente por componente: `E_plot = E_total - E_ext`
- La interpolación para líneas usa `scipy.interpolate.RegularGridInterpolator` con interpolación lineal
- Los gráficos se guardan con resolución de 300 DPI por defecto

## Referencias

- Documentación de Ludwig: http://ludwig.epcc.ed.ac.uk/
- Electrocinética en Ludwig: Ver documentación sobre `psi` y `electrokinetics`

## Autor

Scripts generados para facilitar el análisis de simulaciones de Ludwig con efectos electrocinéticos.
