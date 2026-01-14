# runmulti.sh - Script para ejecutar múltiples simulaciones

## Descripción

`runmulti.sh` es un script que permite ejecutar múltiples simulaciones de Ludwig modificando automáticamente los parámetros de la partícula coloidal en el archivo `config.cds.init.001-001`.

## Parámetros que se pueden modificar

El script modifica los siguientes parámetros:

### Parámetros del coloide (en config.cds.init.001-001):
1. **Carga (q0)**: Línea 51 del archivo config
2. **Posición X**: Primera coordenada en línea 40
3. **Posición Y**: Segunda coordenada en línea 40
4. **Posición Z**: Tercera coordenada en línea 40
5. **Parámetro AL**: Línea 58 del archivo config

### Parámetros del fluido/sistema (vía runbg.sh):
6. **Campo eléctrico**: Vector Ex_Ey_Ez
7. **Temperatura**: Valor de kT
8. **Densidad del fluido**: Valor de rho
9. **Viscosidad**: Valores de viscosidad de corte y volumétrica

## Uso básico

```bash
./runmulti.sh [OPCIONES_RUNBG] [OPCIONES_MULTI]
```

### Opciones específicas de runmulti.sh

**Parámetros del coloide:**
- `--charge VALUES`: Lista de valores de carga separados por comas
- `--position VALUES`: Lista de vectores de posición completos (formato x_y_z separados por comas)
- `--position-x VALUES`: Lista de valores de posición X separados por comas
- `--position-y VALUES`: Lista de valores de posición Y separados por comas
- `--position-z VALUES`: Lista de valores de posición Z separados por comas
- `--al VALUES`: Lista de valores del parámetro AL separados por comas

**Parámetros del fluido/sistema:**
- `--electric-field VALUES`: Lista de vectores de campo eléctrico (formato ex_ey_ez separados por comas)
- `--temperature VALUES`: Lista de valores de temperatura (kT) separados por comas
- `--rho VALUES`: Lista de valores de densidad del fluido separados por comas
- `--viscosity VALUES`: Lista de valores de viscosidad separados por comas

**Opciones de control:**
- `--base-dir DIR`: **OBLIGATORIO** - Directorio base para las simulaciones
- `--param-name NAME`: Nombre del parámetro para nombrar carpetas (opcional, se determina automáticamente)
- `--parallel`: Ejecutar simulaciones en paralelo
- `--max-parallel N`: Límite de simulaciones paralelas simultáneas

**Notas importantes**:
- No se puede usar `--position` junto con `--position-x`, `--position-y` o `--position-z`. Usar `--position` para posiciones completas o las opciones individuales para variar componentes específicas.
- Para variar parámetros del fluido, usar las opciones largas (ej: `--temperature` en lugar de `-k`)
- Todas las opciones de `runbg.sh` se pasan directamente al script. **Excepto** `-o/--output-dir` que se reemplaza por `--base-dir`, y `-e/-k/-r/-v` que se usan para variar parámetros.

### Opciones de runbg.sh disponibles

Todas estas opciones se pueden usar directamente en runmulti.sh:

- `-i/--initial-step`: Paso inicial (0 para nueva simulación)
- `-n/--nsteps`: Número de pasos de simulación
- `-s/--step-interval`: Intervalo de pasos para archivos de salida
- `-x/--size-x`: Tamaño del frame en dirección X
- `-y/--size-yz`: Tamaño del frame en direcciones Y y Z
- `-v/--viscosity`: Viscosidad del fluido
- `-r/--rho`: Densidad del fluido
- `-k/--temperature`: Temperatura (kT)
- `-t/--free-energy`: Modelo de energía libre (fe_electro, none)
- `-e/--electric-field`: Campo eléctrico externo (Ex_Ey_Ez)
- `-m/--mpi-procs`: Número de procesos MPI
- `-g/--grid-mpi`: Grid de descomposición MPI (NX_NY_NZ)
- `-c/--cores`: Threads OpenMP por proceso MPI
- `-w/--fluctuations`: Fluctuaciones LB (0=off, 1=on)
- `-q/--relaxation-scheme`: Esquema de relajación LB (M10, BGK)
- `-z/--solver`: Solver de electrocinética (petsc, sor)
- `-a/--angle-harmonic`: Parámetros angle harmonic
- `-b/--bond-harmonic`: Parámetros bond harmonic
- `-j/--gravity`: Vector de gravedad del coloide
- `-u/--single-monomer`: Modo monómero único (y/n)
- `-l/--fluid-only`: Graficar solo velocidad promedio del fluido (y/n)
- `-d/--delete-files`: Eliminar archivos previos antes de empezar (y/n)
- `-f/--freq-config`: Frecuencia de salida de configuración
- `-p/--plot-interval`: Intervalo de tiempo entre graficas periódicas

## Ejemplos

### Ejemplo 1: Variar solo la carga

```bash
./runmulti.sh \
  --charge 0.5,1.0,1.5,2.0 \
  --base-dir test-charge \
  -i 0 -n 10000 -s 500 \
  -x 32 -y 32 -v 0.5 \
  -t fe_electro -e 0.0_0.0_0.1
```

Esto creará 4 simulaciones con las siguientes carpetas:
- `test-charge_0.5`
- `test-charge_1.0`
- `test-charge_1.5`
- `test-charge_2.0`

### Ejemplo 2: Variar solo la posición en X

```bash
./runmulti.sh \
  --position-x 16.0,16.25,16.5,16.75,17.0 \
  --base-dir test-position \
  -i 0 -n 10000 -s 500 \
  -x 32 -y 32 -v 0.5
```

Esto creará 5 simulaciones con diferentes posiciones X:
- `test-position_16.0`
- `test-position_16.25`
- `test-position_16.5`
- `test-position_16.75`
- `test-position_17.0`

### Ejemplo 3: Variar múltiples parámetros (producto cartesiano)

```bash
./runmulti.sh \
  --charge 1.0,2.0 \
  --position-x 16.0,17.0 \
  --base-dir test-combined \
  -i 0 -n 10000 -s 500 \
  -x 32 -y 32 -v 0.5 \
  -t fe_electro
```

Esto creará 4 simulaciones (2 cargas × 2 posiciones):
- `test-combined_1.0_16.0`
- `test-combined_1.0_17.0`
- `test-combined_2.0_16.0`
- `test-combined_2.0_17.0`

### Ejemplo 4: Variar el parámetro AL

```bash
./runmulti.sh \
  --al 0.4,0.5,0.6,0.7 \
  --base-dir test-al \
  -i 0 -n 10000 -s 500 \
  -x 32 -y 32 -v 0.5 \
  -t fe_electro -e 0.0_0.0_0.0
```

### Ejemplo 5: Variar campo eléctrico

```bash
./runmulti.sh \
  --electric-field 0.001_0.0_0.0,0.005_0.0_0.0,0.01_0.0_0.0 \
  --base-dir test-efield \
  -i 0 -n 10000 -s 500 \
  -x 32 -y 32 -v 0.02
```

Esto creará 3 simulaciones con diferentes campos eléctricos en X:
- `test-efield_0.001_0.0_0.0`
- `test-efield_0.005_0.0_0.0`
- `test-efield_0.01_0.0_0.0`

### Ejemplo 6: Variar temperatura (kT)

```bash
./runmulti.sh \
  --temperature 0.0001,0.0005,0.001,0.005 \
  --base-dir test-temperature \
  -i 0 -n 10000 -s 500 \
  -x 32 -y 32 -v 0.02 -e 0.005_0.0_0.0
```

Esto creará 4 simulaciones con diferentes temperaturas:
- `test-temperature_0.0001`
- `test-temperature_0.0005`
- `test-temperature_0.001`
- `test-temperature_0.005`

### Ejemplo 7: Variar densidad y viscosidad (producto cartesiano)

```bash
./runmulti.sh \
  --rho 0.6,0.8,1.0 \
  --viscosity 0.01,0.02,0.05 \
  --base-dir test-fluid \
  -i 0 -n 10000 -s 500 \
  -x 32 -y 32
```

Esto creará 9 simulaciones (3 densidades × 3 viscosidades):
- `test-fluid_0.6_0.01`
- `test-fluid_0.6_0.02`
- `test-fluid_0.6_0.05`
- `test-fluid_0.8_0.01`
- ... (9 combinaciones totales)

### Ejemplo 8: Variar posición completa (X, Y, Z)

```bash
# Opción 1: Usando --position para vectores completos
# Útil cuando quieres posiciones específicas (ej: diagonal)
./runmulti.sh \
  --position 16.0_16.0_16.0,17.0_17.0_17.0,18.0_18.0_18.0 \
  --base-dir scan-diagonal \
  -i 0 -n 10000 -s 500 -x 32 -y 32 -v 0.5

# Opción 2: Variar solo X, Y y Z se mantienen iguales
./runmulti.sh \
  --position-x 16.0,16.5,17.0 \
  --base-dir scan-x \
  -i 0 -n 10000 -s 500 -x 32 -y 32 -v 0.5

# Opción 3: Variar X e Y independientemente (crea producto cartesiano)
./runmulti.sh \
  --position-x 16.0,17.0 \
  --position-y 16.0,17.0 \
  --base-dir scan-xy \
  -i 0 -n 10000 -s 500 -x 32 -y 32 -v 0.5
```

### Ejemplo 9: Variar múltiples parámetros físicos

```bash
./runmulti.sh \
  --electric-field 0.001_0.0_0.0,0.005_0.0_0.0 \
  --temperature 0.0001,0.0005 \
  --charge 0.5,1.0 \
  --base-dir scan-multiparametro \
  -i 0 -n 10000 -s 500 \
  -x 32 -y 32 -v 0.02
```

Esto creará 8 simulaciones (2 campos × 2 temperaturas × 2 cargas):
- `scan-multiparametro_0.5_0.001_0.0_0.0_0.0001`
- `scan-multiparametro_0.5_0.001_0.0_0.0_0.0005`
- `scan-multiparametro_0.5_0.005_0.0_0.0_0.0001`
- ... (8 combinaciones totales)

### Ejemplo 10: Ejecución en paralelo

```bash
# Ejecutar todas las simulaciones en paralelo (sin límite)
./runmulti.sh \
  --temperature 0.0001,0.0005,0.001 \
  --base-dir test-parallel \
  --parallel \
  -i 0 -n 10000 -s 500 -x 32 -y 32 -v 0.02

# Ejecutar con límite de 4 simulaciones simultáneas
./runmulti.sh \
  --electric-field 0.001_0.0_0.0,0.005_0.0_0.0,0.01_0.0_0.0 \
  --temperature 0.0001,0.0005 \
  --base-dir test-parallel-limited \
  --parallel --max-parallel 4 \
  -i 0 -n 10000 -s 500 -x 32 -y 32 -v 0.02
```

```bash
# Ejemplo de script wrapper para variar temperatura
for temp in 0.0005 0.001 0.002; do
  ./runmulti.sh \
    --charge 1.0 \
    --base-dir test-temp-${temp} \
    -i 0 -n 10000 -s 500 \
    -x 32 -y 32 -v 0.5 \
    -k $temp \
    -t fe_electro -e 0.0_0.0_0.1
done

# Ejemplo de script wrapper para variar densidad
for rho in 0.6 0.7 0.8 0.9; do
  ./runmulti.sh \
    --charge 1.0 \
    --base-dir test-rho-${rho} \
    -i 0 -n 10000 -s 500 \
    -x 32 -y 32 -v 0.5 \
    -r $rho \
    -t fe_electro
done

# Ejemplo de script wrapper para variar campo eléctrico
for efield in "0.0_0.0_0.0" "0.0_0.0_0.05" "0.0_0.0_0.1"; do
  # Reemplazar _ con - para nombre de carpeta
  dirname=$(echo $efield | tr '_' '-')
  ./runmulti.sh \
    --charge 1.0 \
    --base-dir test-efield-${dirname} \
    -i 0 -n 10000 -s 500 \
    -x 32 -y 32 -v 0.5 \
    -t fe_electro -e $efield
done
```

### Ejemplo 7: Usando con MPI

```bash
./runmulti.sh \
  --charge 1.0,1.5,2.0 \
  --base-dir test-mpi \
  -i 0 -n 100000 -s 1000 \
  -x 64 -y 64 -v 0.5 \
  -m 4 -g 2_2_1 -c 10 \
  -t fe_electro -e 0.0_0.0_0.1
```

### Ejemplo 8: Ejecución en paralelo

```bash
# Ejecutar TODAS las simulaciones en paralelo (simultáneamente)
./runmulti.sh \
  --charge 0.5,1.0,1.5,2.0 \
  --base-dir test-parallel \
  --parallel \
  -i 0 -n 10000 -s 500 \
  -x 32 -y 32 -v 0.5

# Ejecutar máximo 4 simulaciones en paralelo a la vez
./runmulti.sh \
  --charge 0.5,1.0,1.5,2.0,2.5,3.0 \
  --base-dir test-parallel-limited \
  --parallel --max-parallel 4 \
  -i 0 -n 10000 -s 500 \
  -x 32 -y 32 -v 0.5

# Ejemplo práctico: 9 simulaciones (3x3), máximo 3 en paralelo
./runmulti.sh \
  --charge 1.0,1.5,2.0 \
  --position-x 16.0,17.0,18.0 \
  --base-dir scan-parallel \
  --parallel --max-parallel 3 \
  -i 0 -n 50000 -s 1000 \
  -x 64 -y 64 -v 0.5 \
  -t fe_electro
```

## Comportamiento del script

1. **Backup automático**: El script crea un backup del archivo `config.cds.init.001-001` original
2. **Restauración**: Para cada simulación, restaura el archivo original antes de modificarlo
3. **Modo de ejecución**:
   - **Secuencial (default)**: Las simulaciones se ejecutan una tras otra
   - **Paralelo (--parallel)**: Todas las simulaciones se lanzan simultáneamente
   - **Paralelo limitado (--parallel --max-parallel N)**: Máximo N simulaciones simultáneas
4. **Limpieza**: Al finalizar todas las simulaciones, restaura el archivo original y elimina el backup

## Estructura de archivos generados

Cada simulación crea su propia carpeta con:
- Archivo `config.cds.init.001-001` modificado
- Archivo `input` con parámetros de la simulación
- Ejecutable `Ludwig.exe`
- Scripts de procesamiento y graficación
- Subcarpeta `graficos/` con scripts de visualización

## Notas importantes

1. El archivo `config.cds.init.001-001` **debe existir** en el directorio actual
2. La opción `--base-dir` es **obligatoria**
3. Al menos un parámetro variable debe especificarse
4. Si múltiples parámetros varían, se ejecutan **todas las combinaciones** (producto cartesiano)
5. El script preserva el archivo de configuración original

### Modo de ejecución

**Secuencial (default)**:
- Las simulaciones se ejecutan una después de otra
- Seguro para recursos limitados
- Fácil de seguir y debuggear

**Paralelo (--parallel)**:
- Todas las simulaciones se lanzan simultáneamente
- **ADVERTENCIA**: Requiere suficiente CPU, memoria y disco
- Cada simulación usa los recursos especificados (MPI, OpenMP, etc.)
- Ejemplo: 10 simulaciones con `-m 4` = 40 procesos MPI simultáneos

**Paralelo limitado (--parallel --max-parallel N)**:
- Máximo N simulaciones ejecutándose al mismo tiempo
- Cuando una termina, se lanza la siguiente
- Balance entre velocidad y recursos
- Recomendado para conjuntos grandes de simulaciones

### Consideraciones de recursos

Al usar `--parallel`, considera:
- **CPU**: Número de simulaciones × procesos MPI × threads OpenMP
- **Memoria**: Cada simulación carga el dominio completo en memoria
- **Disco**: Escritura simultánea de múltiples simulaciones
- **Ejemplo seguro**: `--max-parallel 2` o `--max-parallel 4` en máquina típica

## Debugging

Si algo sale mal, el script:
- Muestra mensajes de error descriptivos
- Restaura el archivo de configuración original
- Se detiene en la primera simulación que falle

Para ver qué comando exactamente se ejecuta, revisa la salida del script que muestra:
```
Running: ./runbg.sh [argumentos completos] -o [directorio]
```
