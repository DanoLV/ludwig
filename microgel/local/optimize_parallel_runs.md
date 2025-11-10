# Optimización para correr 2 simulaciones en paralelo

## Problema identificado
- Ludwig usa GPU moderadamente (~30%)
- El cuello de botella NO es la GPU
- El problema es: I/O, memoria, y scheduler de WSL2

## Soluciones recomendadas

### Opción 1: Mejorar I/O (RECOMENDADO)
Reduce la frecuencia de escritura de archivos:

En tu comando runmulti.sh, aumenta los intervalos:
```bash
./runmulti.sh --parallel --max-parallel 2 \
    -s 1000 \          # Aumenta step-interval (era 500)
    -p 600 \           # Aumenta plot-interval (era 300)
    [otros parámetros]
```

### Opción 2: Usar nice para balancear prioridades
Modifica runbg.sh para que el segundo proceso tenga misma prioridad:
- Ambos procesos tendrán acceso equitativo a CPU/I/O

### Opción 3: Separar directorios de salida en diferentes discos
Si tienes múltiples discos, usa uno para cada simulación

### Opción 4: Modo SECUENCIAL (más confiable)
Si la diferencia de velocidad es inaceptable:
```bash
./runmulti.sh [parámetros sin --parallel]
```
Cada simulación corre óptimamente, pero una después de otra.

## Optimizaciones ya aplicadas en runbg.sh
✓ Process binding para MPI (--bind-to core)
✓ Variables de entorno GPU configuradas
✓ OMP_NUM_THREADS configurado

## Monitoreo
Para verificar durante la ejecución:
```bash
./check_gpu_usage.sh     # GPU
./monitor_performance.sh # CPU, memoria, I/O
```
