# Script para Combinar Plots de Velocidad

## Descripción

`plot_all_velocities.py` es un script Python que combina los plots de velocidad de múltiples subcarpetas con diferentes configuraciones de stencil en un solo gráfico comparativo.

## Características

- **Búsqueda automática**: Encuentra automáticamente todas las subcarpetas que contengan 'stencil_' en su nombre
- **Plots combinados**: Genera gráficos que comparan las velocidades de partícula y fluido de todos los stencils
- **Opciones flexibles**: Permite plotear solo la velocidad del fluido o ambas (partícula y fluido)
- **Visualización clara**: Usa diferentes colores y estilos de línea para distinguir entre diferentes stencils

## Requisitos

- Python 3 con las siguientes bibliotecas:
  - numpy
  - matplotlib

## Uso

### Sintaxis básica

```bash
./plot_all_velocities.py -d <directorio_padre> [opciones]
```

### Opciones

- `-d, --dir`: (Requerido) Directorio padre que contiene las subcarpetas stencil
- `--out_dir`: (Opcional) Directorio de salida para el plot. Por defecto: directorio actual
- `--fluid-only`: (Opcional) Plotear solo la velocidad promedio del fluido

### Ejemplos

#### 1. Generar plot completo (velocidad de partícula y fluido)

```bash
./plot_all_velocities.py -d "/home/bater/Sim/ludwig/microgel/local/single-peskin-Lx_32_Ly_32_Lz_32-D3Q19-stencil_var/single-peskin-Lx_32_Ly_32_Lz_32-D3Q19-stencil_var-n_2000-s_10-L_32-q_1.0-pos_16.20_16.50_16.50-e_0.0_0.0_0.0-kT_0.00001-rho_0.8-eta_0.05-eps_1.0e4-rhoel_0.000E+00"
```

Este comando generará un plot con dos paneles:
- Panel superior: Velocidad de la partícula (Vx, Vy, Vz) para todos los stencils
- Panel inferior: Velocidad del fluido (Vfx, Vfy, Vfz) para todos los stencils

#### 2. Generar plot solo de velocidad del fluido

```bash
./plot_all_velocities.py -d "/home/bater/Sim/ludwig/microgel/local/single-peskin-Lx_32_Ly_32_Lz_32-D3Q19-stencil_var/single-peskin-Lx_32_Ly_32_Lz_32-D3Q19-stencil_var-n_2000-s_10-L_32-q_1.0-pos_16.20_16.50_16.50-e_0.0_0.0_0.0-kT_0.00001-rho_0.8-eta_0.05-eps_1.0e4-rhoel_0.000E+00" --fluid-only
```

Este comando generará un plot de un solo panel con solo las velocidades del fluido.

#### 3. Especificar directorio de salida

```bash
./plot_all_velocities.py -d "/path/to/parent/folder" --out_dir "/path/to/output"
```

## Estructura de Directorios Esperada

El script espera la siguiente estructura de directorios:

```
directorio_padre/
├── nombre-largo-stencil_7/
│   └── proceced_data/
│       └── datosfluid.csv
├── nombre-largo-stencil_19/
│   └── proceced_data/
│       └── datosfluid.csv
└── nombre-largo-stencil_27/
    └── proceced_data/
        └── datosfluid.csv
```

El script:
1. Busca todas las subcarpetas que contengan 'stencil_' en su nombre
2. Extrae el número de stencil del nombre del directorio
3. Lee el archivo `proceced_data/datosfluid.csv` de cada subcarpeta
4. Genera un plot combinado comparando todos los stencils

## Archivo de Salida

El script genera un archivo PNG de alta resolución (300 dpi) con el nombre:
- `velocidades_combined.png`

## Formato de Datos

El script lee archivos CSV con el siguiente formato (columnas):
- 0: cycle (tiempo/paso de simulación)
- 4: vx (velocidad x de partícula)
- 5: vy (velocidad y de partícula)
- 6: vz (velocidad z de partícula)
- 7: <vfx> (velocidad promedio x de fluido)
- 8: <vfy> (velocidad promedio y de fluido)
- 9: <vfz> (velocidad promedio z de fluido)

## Notas

- El script ordena automáticamente los stencils por su número
- Cada stencil se visualiza con un color diferente
- Las componentes x, y, z se distinguen por estilos de línea (sólida, discontinua, punto-raya)
- Los límites de los ejes se ajustan automáticamente para el rango de datos visible
