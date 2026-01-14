# Resumen: Inicialización y Actualización de la Velocidad del Fluido

Este documento describe dónde se inicializa y actualiza la velocidad del fluido en el código de Ludwig.

## INICIALIZACIÓN DE LA VELOCIDAD

### 1. Creación de la estructura hidrodinámica
- **Archivo:** `src/hydro.c` (líneas 49-114)
- **Función:** `hydro_create()`
- Crea el campo de velocidad `hydro->u` como parte de la estructura hidrodinámica

### 2. Lectura desde archivos de checkpoint
- **Archivo:** `src/ludwig.c` (líneas 345-349)
- **Función:** `hydro_io_read()`
- Lee la densidad y velocidad desde archivos guardados al reiniciar una simulación

### 3. Inicialización a cero antes de cada paso
- **Archivo:** `src/hydro.c` (líneas 217-235)
- **Función:** `hydro_u_zero()`
- **Llamada:** `src/ludwig.c` línea 813 en el bucle principal
- Pone la velocidad a cero antes de la colisión

## ACTUALIZACIÓN DE LA VELOCIDAD

La actualización ocurre principalmente durante la **etapa de colisión** del método de Lattice Boltzmann.

### Archivo principal: `src/collision.c`

**Función clave:** `lb_collision_mrt1_site()` (líneas 253-593)

### Pasos de actualización:

#### 1. Extrae velocidad de los momentos (líneas 349-352):
```c
for (ia = 0; ia < NDIM; ia++) {
  for_simd_v(iv, NSIMDVL) u[ia][iv] = mode[(1 + ia)*NSIMDVL+iv];
}
```

#### 2. Actualiza con la fuerza de cuerpo (líneas 370-376):
```c
u[ia][iv] = rrho[iv]*(u[ia][iv] + 0.5*force[ia][iv]);
```

**Esta es la ecuación clave:** `u_nueva = (u_modos + 0.5 × fuerza) / ρ`

#### 3. Guarda la velocidad actualizada (líneas 569-573):
```c
hydro->u->data[addr_rank1(hydro->nsite, 3, index0+iv, ia)] = u[ia][iv];
```

## SECUENCIA EN EL BUCLE PRINCIPAL

En `src/ludwig.c`, el bucle temporal ejecuta:

1. **Línea 813:** `hydro_u_zero()` - Pone velocidad a cero
2. **Líneas 824-825:** `lb_collide()` - **ACTUALIZA la velocidad** basándose en distribuciones y fuerzas
3. **Línea 847:** `lb_propagation()` - Propaga las distribuciones
4. **Línea 966:** `hydro_io_write()` - Escribe velocidad a archivo (opcional)

## RESUMEN

- La velocidad se **inicializa** en `src/hydro.c`
- La velocidad se **actualiza** cada paso temporal en `src/collision.c` (líneas 370-376) durante la colisión MRT
- La ecuación de actualización incorpora las fuerzas de cuerpo y los momentos de la función de distribución

## ARCHIVOS CLAVE

| Archivo | Propósito |
|---------|-----------|
| `src/hydro.h` | Definición de estructura hidrodinámica |
| `src/hydro.c` | Funciones core (creación, inicialización) |
| `src/hydro_rt.c` | Configuración runtime para hidrodinámica |
| `src/collision.c` | Colisión Lattice Boltzmann (actualización de velocidad) |
| `src/propagation.c` | Propagación de distribuciones después de colisión |
| `src/ludwig.c` | Bucle principal de simulación (orquesta todas las actualizaciones) |

---

**Fecha de generación:** 25 de octubre de 2025
**Código fuente:** Ludwig - Lattice Boltzmann Code
