# Ludwig Microgel Simulation Runner

Este documento explica el uso del script `runbg_new.sh` para ejecutar simulaciones de microgeles con Ludwig usando soporte MPI.

## Descripción General

El script `runbg_new.sh` es una herramienta mejorada para ejecutar simulaciones de dinámica de fluidos con Ludwig. Incluye:

- Soporte completo para MPI (Message Passing Interface) con múltiples procesos
- Soporte para OpenMP con múltiples threads por proceso
- Generación automática y periódica de gráficos durante la simulación
- Opciones de línea de comandos largas y descriptivas
- Validación robusta de parámetros de entrada
- Gestión automática de archivos de configuración

## Requisitos

- Ludwig ejecutable compilado (`Ludwig.exe`)
- MPI instalado (`mpirun` o `mpiexec`)
- Python 3.x para scripts de post-procesamiento
- Archivo de configuración inicial (`config.cds.init.001-001`)
- Archivo de entrada de parámetros (`input`)
- Scripts auxiliares:
  - `runplot.sh` - Script de generación de gráficos
  - `del.sh` - Script de limpieza de archivos
  - Scripts Python de análisis (`calculosvel.py`, `plotvel.py`, etc.)

## Uso Básico

```bash
./runbg_new.sh [OPCIONES]
```

Para ver la ayuda completa:
```bash
./runbg_new.sh --help
```

## Parámetros Obligatorios

Estos parámetros son requeridos para todas las simulaciones:

| Opción Corta | Opción Larga | Descripción | Ejemplo |
|--------------|--------------|-------------|---------|
| `-i` | `--initial-step` | Paso inicial (0 para nueva simulación) | `-i 0` |
| `-n` | `--nsteps` | Número total de pasos a calcular | `-n 1000000` |
| `-o` | `--output-dir` | Directorio de salida para resultados | `-o resultados` |
| `-s` | `--step-interval` | Intervalo de pasos para archivos de salida | `-s 500` |

## Parámetros Opcionales

Los parámetros están organizados alfabéticamente para facilitar su búsqueda:

### Parámetros de Computación

| Opción Corta | Opción Larga | Descripción | Valor por Defecto | Ejemplo |
|--------------|--------------|-------------|-------------------|---------|
| `-c` | `--cores` | Threads OpenMP por proceso MPI | 1 | `-c 10` |
| `-m` | `--mpi-procs` | Número de procesos MPI | 1 | `-m 4` |
| `-g` | `--grid-mpi` | Descomposición de grilla MPI (NX_NY_NZ) | - | `-g 2_2_1` |

### Parámetros de Geometría

| Opción Corta | Opción Larga | Descripción | Valor por Defecto | Ejemplo |
|--------------|--------------|-------------|-------------------|---------|
| `-x` | `--size-x` | Tamaño del dominio en dirección X | - | `-x 26` |
| `-y` | `--size-yz` | Tamaño del dominio en direcciones Y y Z | - | `-y 26` |

### Parámetros Físicos del Fluido

| Opción Corta | Opción Larga | Descripción | Valor por Defecto | Ejemplo |
|--------------|--------------|-------------|-------------------|---------|
| `-v` | `--viscosity` | Viscosidad del fluido | - | `-v 0.5` |
| `-r` | `--rho` | Densidad del fluido (rho0) | - | `-r 1.0` |
| `-k` | `--temperature` | Temperatura (kT) para fluctuaciones térmicas | - | `-k 0.0005` |
| `-w` | `--fluctuations` | Fluctuaciones LB: 0=desactivadas, 1=activadas | 0 | `-w 1` |
| `-q` | `--relaxation-scheme` | Esquema de relajación LB: 'M10' o 'BGK' | M10 | `-q BGK` |
| `-t` | `--free-energy` | Modelo de energía libre | none | `-t fe_electro` |
| `-e` | `--electric-field` | Campo eléctrico externo (Ex_Ey_Ez) | - | `-e 0.0_0.0_0.1` |
| `-z` | `--solver` | Solver electrocinético | sor | `-z petsc` |

### Parámetros de Coloides e Interacciones

| Opción Corta | Opción Larga | Descripción | Valor por Defecto | Ejemplo |
|--------------|--------------|-------------|-------------------|---------|
| `-j` | `--gravity` | Vector de gravedad del coloide (Gx_Gy_Gz) | - | `-j 0.0_0.0_-0.000005` |
| `-b` | `--bond-harmonic` | Potencial armónico de enlace: on_k_r0 | - | `-b 1_5.0e-5_0.5` |
| `-a` | `--angle-harmonic` | Potencial armónico angular: on_k_theta0 | - | `-a 1_1.0e-04_2.0944` |

**Nota sobre potenciales:**
- **bond-harmonic**: Formato `on_k_r0` donde:
  - `on`: 1 para activar, 0 para desactivar
  - `k`: Constante del resorte (bond_harmonic_k)
  - `r0`: Longitud de equilibrio (bond_harmonic_r0)
- **angle-harmonic**: Formato `on_k_theta0` donde:
  - `on`: 1 para activar, 0 para desactivar
  - `k`: Constante del potencial (angle_harmonic_k)
  - `theta0`: Ángulo de equilibrio en radianes (angle_harmonic_theta0)

### Parámetros de Control

| Opción Corta | Opción Larga | Descripción | Valor por Defecto | Ejemplo |
|--------------|--------------|-------------|-------------------|---------|
| `-d` | `--delete-files` | Eliminar archivos previos (y/n) | n | `-d y` |
| `-f` | `--freq-config` | Frecuencia de salida de configuración | - | `-f 1000` |
| `-p` | `--plot-interval` | Intervalo de tiempo entre gráficos (segundos) | 300 | `-p 600` |
| `-u` | `--single-monomer` | Modo monómero único (y/n) | n | `-u y` |
| `-l` | `--fluid-only` | Graficar solo velocidad promedio del fluido (y/n) | n | `-l y` |

## Ejemplos de Uso

### Ejemplo 1: Nueva simulación básica

Simulación simple con un solo proceso MPI:

```bash
./runbg_new.sh \
  --initial-step 0 \
  --nsteps 1000000 \
  --step-interval 500 \
  --output-dir mi_simulacion \
  --size-x 26 \
  --size-yz 26 \
  --viscosity 0.5
```

### Ejemplo 2: Simulación con campo eléctrico

Simulación con energía libre electrocinética y campo eléctrico aplicado:

```bash
./runbg_new.sh \
  --initial-step 0 \
  --nsteps 2000000 \
  --step-interval 1000 \
  --output-dir sim_electro \
  --size-x 32 \
  --size-yz 32 \
  --viscosity 0.8 \
  --free-energy fe_electro \
  --electric-field 0.0_0.0_0.1 \
  --solver petsc \
  --rho 1.0
```

### Ejemplo 3: Simulación paralela con MPI

Simulación con 4 procesos MPI en una grilla 2×2×1, cada proceso con 8 threads OpenMP:

```bash
./runbg_new.sh \
  --initial-step 0 \
  --nsteps 5000000 \
  --step-interval 500 \
  --output-dir sim_paralela \
  --mpi-procs 4 \
  --grid-mpi 2_2_1 \
  --cores 8 \
  --size-x 64 \
  --size-yz 32 \
  --viscosity 0.5 \
  --free-energy fe_electro \
  --electric-field 0.0_0.0_0.05
```

### Ejemplo 4: Simulación con temperatura y fluctuaciones térmicas

Simulación con temperatura kT, fluctuaciones LB activadas y esquema de relajación BGK:

```bash
./runbg_new.sh \
  --initial-step 0 \
  --nsteps 1000000 \
  --step-interval 500 \
  --output-dir sim_thermal \
  --size-x 32 \
  --size-yz 32 \
  --viscosity 0.5 \
  --temperature 0.0005 \
  --fluctuations 1 \
  --relaxation-scheme BGK \
  --rho 0.8
```

**Nota:** Cuando se activan las fluctuaciones (`--fluctuations 1`), es recomendable especificar también una temperatura (`--temperature`) para controlar la magnitud de las fluctuaciones térmicas.

### Ejemplo 5: Simulación con potenciales moleculares (microgel)

Simulación con potenciales de enlace (bond) y ángulo (angle) armónicos, más gravedad:

```bash
./runbg_new.sh \
  --initial-step 0 \
  --nsteps 2000000 \
  --step-interval 1000 \
  --output-dir microgel_sim \
  --size-x 64 \
  --size-yz 64 \
  --viscosity 1e-04 \
  --temperature 0.0005 \
  --bond-harmonic 1_5.0e-5_0.5 \
  --angle-harmonic 1_1.0e-04_2.0944 \
  --gravity 0.0_0.0_-0.000005 \
  --free-energy fe_electro \
  --electric-field 0.0_0.0_0.0 \
  --rho 0.8
```

**Explicación de parámetros moleculares:**
- `--bond-harmonic 1_5.0e-5_0.5`: Activa (1) el potencial de enlace con k=5.0e-5 y r0=0.5
- `--angle-harmonic 1_1.0e-04_2.0944`: Activa (1) el potencial angular con k=1.0e-04 y θ0=2.0944 rad (120°)
- `--gravity 0.0_0.0_-0.000005`: Aplica gravedad en dirección -Z

### Ejemplo 6: Continuar simulación existente

Para continuar una simulación previamente iniciada:

```bash
./runbg_new.sh \
  --initial-step 1000000 \
  --nsteps 500000 \
  --step-interval 500 \
  --output-dir mi_simulacion
```

**Nota:** Al continuar una simulación, el directorio de salida debe existir y contener los archivos de configuración previos.

### Ejemplo 7: Usando opciones cortas

Ejemplo equivalente al Ejemplo 5 con opciones cortas:

```bash
./runbg_new.sh -i 0 -n 2000000 -s 1000 -o microgel_sim \
  -x 64 -y 64 -v 1e-04 -k 0.0005 \
  -b 1_5.0e-5_0.5 -a 1_1.0e-04_2.0944 \
  -j 0.0_0.0_-0.000005 -t fe_electro -e 0.0_0.0_0.0 -r 0.8
```

## Comportamiento del Script

### Nueva Simulación (initial-step = 0)

Cuando se inicia una nueva simulación:

1. **Validación:** Verifica que el directorio de salida NO exista
2. **Creación:** Crea el directorio de salida
3. **Preparación:** Copia archivos de configuración inicial
4. **MPI Setup:** Si se usan múltiples procesos MPI, crea configuraciones para cada proceso
5. **Archivos:** Copia todos los scripts necesarios (simulación, análisis, gráficos)
6. **Ejecución:** Lanza la simulación en segundo plano
7. **Monitoreo:** Genera gráficos periódicamente mientras la simulación corre

### Continuar Simulación (initial-step > 0)

Cuando se continúa una simulación existente:

1. **Validación:** Verifica que el directorio de salida EXISTA
2. **Verificación:** Comprueba que existan archivos de configuración previos
3. **Actualización:** Actualiza parámetros en el archivo `input`
4. **Ejecución:** Continúa la simulación desde el paso especificado

### Generación de Gráficos

Durante la ejecución:

- **Periódicos:** Se generan gráficos cada `--plot-interval` segundos (default: 300s)
- **Automático:** El script detecta nuevos archivos de salida y los procesa
- **Final:** Al terminar la simulación, se genera un gráfico final con todos los datos

## Archivos Generados

En el directorio de salida se generan:

- `output.txt` - Salida estándar de Ludwig
- `outputplot.txt` - Salida de los scripts de gráficos
- `config.cds*` - Archivos de configuración de partículas coloidales
- `vel-*` - Archivos de campo de velocidad
- `*.csv` - Datos procesados en formato CSV
- `*.png` / `*.pdf` - Gráficos generados

## Configuración de MPI

### Número de Procesos

El parámetro `--mpi-procs` especifica cuántos procesos MPI usar:

```bash
--mpi-procs 4  # Usa 4 procesos MPI
```

### Grilla de Descomposición

El parámetro `--grid-mpi` especifica cómo dividir el dominio espacial:

```bash
--grid-mpi 2_2_1  # Divide en 2×2×1 = 4 subdominios
--grid-mpi 2_1_1  # Divide en 2×1×1 = 2 subdominios
--grid-mpi 4_2_1  # Divide en 4×2×1 = 8 subdominios
```

**Importante:** El producto NX × NY × NZ debe ser igual a `--mpi-procs`

### Threads OpenMP

El parámetro `--cores` especifica threads por proceso MPI:

```bash
--cores 8  # Cada proceso MPI usa 8 threads OpenMP
```

**Total de threads = mpi-procs × cores**

Ejemplo: `--mpi-procs 4 --cores 8` usa 32 threads totales.

## Solución de Problemas

### Error: "Missing output directory"

**Causa:** No se especificó el parámetro obligatorio `-o` o `--output-dir`

**Solución:**
```bash
./runbg_new.sh -i 0 -n 1000 -s 100 -o mi_directorio
```

### Error: "Output directory already exists"

**Causa:** Intentando iniciar nueva simulación (initial-step=0) pero el directorio ya existe

**Soluciones:**
- Usar un nombre diferente para el directorio
- Eliminar el directorio existente
- Continuar la simulación con `--initial-step` > 0

### Error: "Output directory does not exist"

**Causa:** Intentando continuar simulación (initial-step>0) pero el directorio no existe

**Solución:** Verificar el nombre del directorio o iniciar nueva simulación con `--initial-step 0`

### Error: "No MPI command found"

**Causa:** No está instalado MPI (ni `mpirun` ni `mpiexec`)

**Solución:** Instalar MPI (OpenMPI o MPICH)
```bash
# Ubuntu/Debian
sudo apt-get install openmpi-bin libopenmpi-dev

# CentOS/RHEL
sudo yum install openmpi openmpi-devel
```

### La simulación corre pero no genera gráficos

**Posibles causas:**
- Scripts Python no tienen permisos de ejecución
- Faltan dependencias de Python (matplotlib, numpy, etc.)
- Error en `runplot.sh`

**Solución:** Revisar `outputplot.txt` para errores específicos

## Diferencias con la Versión Anterior (runbg.sh)

| Aspecto | Versión Antigua | Versión Nueva |
|---------|----------------|---------------|
| Opciones | Solo cortas (`-n`) | Cortas y largas (`-n`, `--nsteps`) |
| Orden | Aleatorio | Alfabético |
| Nombres | Ambiguos (`-a`, `-z`) | Descriptivos (`--freq-config`, `--solver`) |
| Ayuda | Solo comentarios | Función `--help` integrada |
| Validación | Básica | Mejorada con mensajes claros |
| Mensajes | Escasos | Informativos con resumen |
| Documentación | En comentarios | README completo |

## Mapeo de Opciones: Antigua → Nueva

| Antigua | Nueva Corta | Nueva Larga | Descripción |
|---------|-------------|-------------|-------------|
| `-a` | `-f` | `--freq-config` | Frecuencia de configuración |
| `-c` | `-c` | `--cores` | Número de cores |
| `-d` | `-d` | `--delete-files` | Eliminar archivos |
| `-e` | `-t` | `--free-energy` | Modelo de energía |
| `-f` | `-e` | `--electric-field` | Campo eléctrico |
| `-g` | `-z` | `--solver` | Tipo de solver |
| `-i` | `-i` | `--initial-step` | Paso inicial |
| `-l` | `-x` | `--size-x` | Tamaño en X |
| `-m` | `-m` | `--mpi-procs` | Procesos MPI |
| `-n` | `-n` | `--nsteps` | Número de pasos |
| `-o` | `-o` | `--output-dir` | Directorio de salida |
| `-p` | `-s` | `--step-interval` | Intervalo de pasos |
| `-q` | `-l` | `--fluid-only` | Solo fluido |
| `-r` | `-g` | `--grid-mpi` | Grilla MPI |
| `-s` | `-u` | `--single-monomer` | Monómero único |
| `-t` | `-p` | `--plot-interval` | Intervalo de gráficos |
| `-v` | `-v` | `--viscosity` | Viscosidad |
| `-y` | `-y` | `--size-yz` | Tamaño en Y-Z |
| `-z` | `-r` | `--rho` | Densidad |
| N/A | `-a` | `--angle-harmonic` | Potencial angular (NUEVO) |
| N/A | `-b` | `--bond-harmonic` | Potencial de enlace (NUEVO) |
| N/A | `-w` | `--fluctuations` | Fluctuaciones LB (NUEVO) |
| N/A | `-j` | `--gravity` | Gravedad del coloide (NUEVO) |
| N/A | `-k` | `--temperature` | Temperatura kT (NUEVO) |
| N/A | `-q` | `--relaxation-scheme` | Esquema de relajación LB (NUEVO) |

## Consejos de Rendimiento

1. **Optimización MPI:**
   - Para dominios grandes, usar más procesos MPI
   - La grilla debe ser balanceada (ej: 2_2_1 mejor que 4_1_1 para dominio cúbico)

2. **Optimización OpenMP:**
   - `cores` debe ser ≤ número de cores físicos por nodo
   - Para cluster: ajustar según cores disponibles por nodo

3. **Tamaño de Dominio:**
   - Cada subdominio MPI debe tener al menos ~10-20 nodos de grilla por dimensión
   - Ejemplo: dominio 64×64×32 con grilla 4_4_2 → subdominios de 16×16×16

4. **Intervalos de Salida:**
   - `step-interval` más grande → menos I/O, simulación más rápida
   - `plot-interval` ajustar según velocidad de generación de datos

## Referencias

- Documentación de Ludwig: https://ludwig.epcc.ed.ac.uk/
- Tutorial de MPI: https://mpitutorial.com/
- OpenMP: https://www.openmp.org/

## Autor y Licencia

Script desarrollado para simulaciones de microgeles con Ludwig.

Para reportar problemas o sugerencias, contactar al administrador del proyecto.
