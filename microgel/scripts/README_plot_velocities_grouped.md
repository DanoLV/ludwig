# Script para Plotear Velocidades Agrupadas por Parámetro

## Descripción

`plot_velocities_grouped.py` es un script Python flexible que combina plots de velocidad de múltiples simulaciones, agrupándolas por el parámetro que elijas (stencil, rhoel, posición, etc.).

## Características

- **Agrupación flexible**: Elige por qué parámetro agrupar los datos
- **Búsqueda recursiva**: Encuentra automáticamente todos los archivos `datosfluid.csv` en la estructura de directorios
- **Extracción automática de parámetros**: Lee los parámetros de los nombres de los directorios
- **Opciones de visualización**: Elige entre velocidad de partícula, fluido, o ambas
- **Filtrado por componente**: Plotea todas las componentes (x,y,z) o solo una específica

## Requisitos

- Python 3 con las siguientes bibliotecas:
  - numpy
  - matplotlib

## Uso

### Sintaxis básica

```bash
./plot_velocities_grouped.py -d <directorio_padre> --group-by <parámetro> [opciones]
```

### Parámetros obligatorios

- `-d, --dir`: Directorio padre que contiene los datos (puede ser cualquier nivel de la jerarquía)
- `--group-by`: Parámetro(s) por el cual agrupar los datos. Puede especificar uno o varios parámetros separados por espacio
  - Un parámetro: `--group-by stencil`
  - Dos parámetros: `--group-by stencil rhoel`
  - Tres o más: `--group-by Lx stencil rhoel`

### Parámetros disponibles para agrupar

El script extrae automáticamente estos parámetros de los nombres de directorios:

- `stencil`: Tipo de stencil (7, 19, 27, etc.)
- `rhoel`: Valor de rhoel (0.000E+00, 1.000E-04, etc.)
- `n`: Número de pasos
- `s`: Parámetro s
- `L`: Tamaño del dominio
- `q`: Carga
- `pos`: Posición (x_y_z)
- `e`: Campo eléctrico (ex_ey_ez)
- `kT`: Temperatura
- `rho`: Densidad
- `eta`: Viscosidad
- `eps`: Epsilon (permitividad)
- `Lx`, `Ly`, `Lz`: Dimensiones del dominio

### Opciones adicionales

- `--plot-type`: Tipo de plot a generar
  - `particle`: Solo velocidad de partícula (por defecto)
  - `fluid`: Solo velocidad de fluido
  - `both`: Ambas velocidades en paneles separados

- `--component`: Componente de velocidad a plotear
  - `all`: Todas las componentes x, y, z (por defecto)
  - `x`: Solo componente x
  - `y`: Solo componente y
  - `z`: Solo componente z

- `--out_dir`: Directorio de salida (por defecto: directorio actual). Se usa cuando no se especifica `-o`

- `-o, --output`: Ruta completa del archivo de salida. Si se especifica, se ignora `--out_dir`
  - Ejemplo: `-o ./mis_plots/velocidades_custom.png`
  - El directorio se crea automáticamente si no existe

## Ejemplos

### 1. Agrupar por stencil y mostrar velocidad de partícula

```bash
./plot_velocities_grouped.py \
  -d "single-peskin-Lx_32_Ly_32_Lz_32-D3Q19-stencil_var/single-peskin-...-rhoel_0.000E+00" \
  --group-by stencil
```

Esto genera un plot comparando las velocidades de partícula para diferentes stencils (7, 19, 27).

### 2. Agrupar por rhoel y mostrar solo velocidad de fluido

```bash
./plot_velocities_grouped.py \
  -d "single-peskin-Lx_32_Ly_32_Lz_32-D3Q19-stencil_var" \
  --group-by rhoel \
  --plot-type fluid
```

Esto genera un plot comparando las velocidades de fluido para diferentes valores de rhoel.

### 3. Agrupar por stencil y mostrar ambas velocidades

```bash
./plot_velocities_grouped.py \
  -d "single-peskin-Lx_32_Ly_32_Lz_32-D3Q19-stencil_var" \
  --group-by stencil \
  --plot-type both
```

Genera un plot con dos paneles: uno para velocidad de partícula y otro para velocidad de fluido.

### 4. Agrupar por rhoel y mostrar solo componente x

```bash
./plot_velocities_grouped.py \
  -d "single-peskin-Lx_32_Ly_32_Lz_32-D3Q19-stencil_var" \
  --group-by rhoel \
  --component x
```

Muestra solo la componente x de la velocidad para diferentes valores de rhoel.

### 5. Especificar directorio de salida

```bash
./plot_velocities_grouped.py \
  -d "single-peskin-Lx_32_Ly_32_Lz_32-D3Q19-stencil_var" \
  --group-by stencil \
  --out_dir "./plots"
```

### 6. Especificar ruta completa del archivo de salida

```bash
./plot_velocities_grouped.py \
  -d "single-peskin-Lx_32_Ly_32_Lz_32-D3Q19-stencil_var" \
  --group-by rhoel \
  -o "./resultados/graficos/velocidades_rhoel_custom.png"
```

Este comando:
- Agrupa por rhoel
- Guarda el plot en `./resultados/graficos/velocidades_rhoel_custom.png`
- Crea el directorio `./resultados/graficos/` si no existe

### 7. Agrupar por múltiples parámetros (stencil + rhoel)

```bash
./plot_velocities_grouped.py \
  -d "single-peskin-Lx_32_Ly_32_Lz_32-D3Q19-stencil_var" \
  --group-by stencil rhoel \
  --component x
```

Este comando agrupa por la **combinación** de stencil y rhoel, creando una curva separada para cada combinación única (stencil=7 + rhoel=0.000E+00, stencil=7 + rhoel=1.000E-04, etc.).

### 8. Agrupar por tres parámetros

```bash
./plot_velocities_grouped.py \
  -d "single-peskin-Lx_32_Ly_32_Lz_32-D3Q19-stencil_var" \
  --group-by Lx stencil rhoel \
  --plot-type both \
  -o "./plots/multi_param.png"
```

Agrupa por Lx, stencil y rhoel simultáneamente. Útil cuando tienes múltiples dimensiones de sistema.

## Estructura de Directorios

El script funciona con cualquier estructura jerárquica. Ejemplos:

### Estructura típica de 2 niveles:

```
directorio_padre/
├── param1_value1/
│   ├── param2_value1/
│   │   └── proceced_data/
│   │       └── datosfluid.csv
│   └── param2_value2/
│       └── proceced_data/
│           └── datosfluid.csv
└── param1_value2/
    └── param2_value1/
        └── proceced_data/
            └── datosfluid.csv
```

### Ejemplo real:

```
single-peskin-Lx_32_Ly_32_Lz_32-D3Q19-stencil_var/
├── single-peskin-...-rhoel_0.000E+00/
│   ├── single-peskin-...-rhoel_0.000E+00-stencil_7/
│   │   └── proceced_data/datosfluid.csv
│   ├── single-peskin-...-rhoel_0.000E+00-stencil_19/
│   │   └── proceced_data/datosfluid.csv
│   └── single-peskin-...-rhoel_0.000E+00-stencil_27/
│       └── proceced_data/datosfluid.csv
├── single-peskin-...-rhoel_1.000E-04/
│   ├── single-peskin-...-rhoel_1.000E-04-stencil_7/
│   │   └── proceced_data/datosfluid.csv
│   └── ...
└── ...
```

## Archivo de Salida

El script genera un archivo PNG de alta resolución (300 dpi) con el nombre:
- Un parámetro: `velocidades_by_<parámetro>.png`
- Múltiples parámetros: `velocidades_by_<param1>_<param2>_<param3>.png`

Por ejemplo:
- `velocidades_by_stencil.png` (un parámetro)
- `velocidades_by_rhoel.png` (un parámetro)
- `velocidades_by_stencil_rhoel.png` (dos parámetros)
- `velocidades_by_Lx_stencil_rhoel.png` (tres parámetros)

## Formato de Datos

El script lee archivos CSV (`datosfluid.csv`) con el siguiente formato:

| Columna | Descripción |
|---------|-------------|
| 0 | cycle (tiempo/paso de simulación) |
| 4 | vx (velocidad x de partícula) |
| 5 | vy (velocidad y de partícula) |
| 6 | vz (velocidad z de partícula) |
| 7 | <vfx> (velocidad promedio x de fluido) |
| 8 | <vfy> (velocidad promedio y de fluido) |
| 9 | <vfz> (velocidad promedio z de fluido) |

## Cómo funciona

1. **Búsqueda recursiva**: El script busca todos los archivos `datosfluid.csv` en la estructura de directorios
2. **Extracción de parámetros**: Para cada archivo, extrae los parámetros de los nombres de los directorios padre usando expresiones regulares
3. **Agrupación**: Agrupa los archivos según el parámetro especificado en `--group-by`
4. **Visualización**: Genera un plot con diferentes colores para cada grupo

## Notas

- El script ordena automáticamente los grupos alfabéticamente
- Cada grupo se visualiza con un color diferente
- Las componentes x, y, z se distinguen por estilos de línea:
  - Sólida (`-`) para componente x
  - Discontinua (`--`) para componente y
  - Punto-raya (`-.`) para componente z
- Si hay múltiples archivos en un mismo grupo, todos se plotean con el mismo color
- Los límites de los ejes se ajustan automáticamente para el rango de datos visible

## Diferencias con `plot_all_velocities.py`

- `plot_all_velocities.py`: Script antiguo que solo agrupa por stencil
- `plot_velocities_grouped.py`: Script nuevo y flexible que permite agrupar por cualquier parámetro

Se recomienda usar `plot_velocities_grouped.py` para mayor flexibilidad.
