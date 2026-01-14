# Cambios para conservación de momento y suma total de fuerzas

**Fecha**: 2025-11-12
**Autor**: Modificaciones para cálculo de momento y fuerzas en partículas subgrid

---

## Resumen

Se implementó un sistema completo para calcular y reportar:
1. **Conservación de momento del sistema** (incluyendo partículas subgrid con masa efectiva)
2. **Suma total de fuerzas del sistema** (fluido + coloides resueltos + coloides subgrid)

Ambos cálculos utilizan **Klein summation** para minimizar errores de redondeo numérico.

---

## 1. CAMBIOS PARA CONSERVACIÓN DE MOMENTO

### 1.1. Archivo: `src/stats_colloid.c`

#### Líneas 63-70: Masa efectiva para partículas subgrid

**Cambio principal**: Las partículas subgrid ahora tienen masa efectiva calculada como el volumen de fluido desplazado.

```c
/*CHANGE INIT - Subgrid momentum */
if (pc->s.bc == COLLOID_BC_SUBGRID) {
  /* Subgrid particles: use effective mass from displaced fluid volume */
  PI_DOUBLE(pi);
  double volume = (4.0/3.0) * pi * pow(pc->s.ah, 3.0);
  mass = rho0 * volume;
}
/*CHANGE END - Subgrid momentum */
```

**Antes**: `mass = 0.0` (sin inercia)
**Ahora**: `mass = ρ₀ × V_particula = ρ₀ × (4/3)π × ah³`

donde:
- `ρ₀` = densidad del fluido (colloids_info_rho0)
- `ah` = radio hidrodinámico (pc->s.ah)

#### Líneas 184-236: Nueva función `stats_colloid_momentum_subgrid()`

Función para calcular **solo** el momento de partículas subgrid (diagnóstico separado).

```c
int stats_colloid_momentum_subgrid(colloids_info_t * cinfo, double g[3])
```

**Retorna**: `g[3]` = momento total de todas las partículas subgrid

---

### 1.2. Archivo: `src/stats_colloid.h`

#### Líneas 19-21: Declaración nueva función

```c
/*CHANGE INIT - Subgrid momentum */
int stats_colloid_momentum_subgrid(colloids_info_t * cinfo, double g[3]);
/*CHANGE END - Subgrid momentum */
```

---

### 1.3. Archivo: `src/ludwig.c`

#### Líneas 1111-1113: Variables adicionales

```c
/*CHANGE INIT - Subgrid momentum */
int nsubgrid;
/*CHANGE END - Subgrid momentum */
```

#### Líneas 1118-1120: Variable momento subgrid

```c
/*CHANGE INIT - Subgrid momentum */
double gsubgrid[3];  /* Subgrid colloid momentum */
/*CHANGE END - Subgrid momentum */
```

#### Líneas 1144-1149: Cálculo momento subgrid

```c
/*CHANGE INIT - Subgrid momentum */
nsubgrid = ludwig->collinfo->nsubgrid;
if (nsubgrid > 0) {
  stats_colloid_momentum_subgrid(ludwig->collinfo, gsubgrid);
}
/*CHANGE END - Subgrid momentum */
```

#### Líneas 1174-1192: Reporte en pantalla (CORREGIDO)

**IMPORTANTE**: Se corrigió un bug donde `[colloids]` mostraba el momento de TODOS los coloides (resueltos + subgrid).

```c
/*CHANGE INIT - Subgrid momentum */
/* OLD CODE:
if (ncolloid > 0) {
  pe_info(pe, "[colloids] %14.7e %14.7e %14.7e\n", gc[X], gc[Y], gc[Z]);
}
*/
/* NEW CODE: Separate resolved and subgrid colloid momentum */
if (ncolloid > 0) {
  /* gc includes all colloids, subtract subgrid to show only resolved */
  double gc_resolved[3];
  gc_resolved[X] = gc[X] - gsubgrid[X];
  gc_resolved[Y] = gc[Y] - gsubgrid[Y];
  gc_resolved[Z] = gc[Z] - gsubgrid[Z];
  pe_info(pe, "[colloids] %14.7e %14.7e %14.7e\n", gc_resolved[X], gc_resolved[Y], gc_resolved[Z]);
}
if (nsubgrid > 0) {
  pe_info(pe, "[subgrid ] %14.7e %14.7e %14.7e\n", gsubgrid[X], gsubgrid[Y], gsubgrid[Z]);
}
/*CHANGE END - Subgrid momentum */
```

**Explicación del fix**:
- `stats_colloid_momentum()` calcula momento de **TODOS** los coloides → `gc`
- `stats_colloid_momentum_subgrid()` calcula momento de **SOLO subgrid** → `gsubgrid`
- Por lo tanto: `gc_resolved = gc - gsubgrid` (solo coloides resueltos)

**Output esperado**:
```
Momentum - x y z
[total   ]  1.234567e-10  2.345678e-10  3.456789e-10
[fluid   ]  1.000000e-10  2.000000e-10  3.000000e-10
[colloids]  0.000000e+00  0.000000e+00  0.000000e+00  <- SOLO resueltos
[subgrid ]  2.345670e-11  3.456780e-11  4.567890e-11  <- SOLO subgrid
[walls   ]  0.000000e+00  0.000000e+00  0.000000e+00
```

---

## 2. CAMBIOS PARA SUMA TOTAL DE FUERZAS

### 2.1. Archivo: `src/stats_velocity.c`

#### Líneas 25-28: Includes adicionales

```c
/*CHANGE INIT - Total force calculation */
#include "util_sum.h"
#include "colloids.h"
/*CHANGE END - Total force calculation */
```

#### Líneas 127-241: Función `stats_total_force()`

Calcula la suma total de fuerzas en el sistema usando **Klein summation** (compensación doble para minimizar errores de redondeo).

```c
int stats_total_force(hydro_t * hydro, map_t * map, colloids_info_t * cinfo,
                      double ffluid[3], double fcoll[3], double fsubgrid[3],
                      double ftotal[3])
```

**Componentes calculados**:
- **ffluid[3]**: Suma de fuerzas en todos los nodos de fluido (`hydro->force`)
- **fcoll[3]**: Suma de fuerzas en coloides resueltos (`force + fex`)
- **fsubgrid[3]**: Suma de fuerzas en coloides subgrid (`fex`)
- **ftotal[3]**: Suma total = ffluid + fcoll + fsubgrid

**Método numérico**: Klein summation (mejor que Kahan summation)
- Usa estructura `klein_t` con 3 componentes: sum, cs, ccs
- Reduce errores de redondeo en sumas largas
- MPI_Reduce con operación personalizada `klein_mpi_sum`

**Pseudocódigo**:
```
1. Inicializar Klein sums a cero para cada componente
2. Copiar fuerzas desde device a host
3. Loop sobre nodos fluidos:
   - Si status == MAP_FLUID:
     - Leer fuerza en el nodo
     - Klein_add_double a ffluid_local
4. Loop sobre coloides:
   - Si SUBGRID: Klein_add_double fex a fsubgrid_local
   - Si RESOLVED: Klein_add_double (force + fex) a fcoll_local
5. MPI_Reduce con klein_mpi_type y klein_mpi_sum
6. Extraer sumas finales con klein_sum()
```

#### Líneas 243-299: Función `stats_total_force_write()`

Escribe fuerzas a archivo de salida.

```c
int stats_total_force_write(hydro_t * hydro, map_t * map, colloids_info_t * cinfo,
                             int timestep, const char * filename)
```

**Características**:
- Solo el rank 0 escribe al archivo
- Modo append ("a")
- Escribe header si el archivo está vacío
- Formato: valores separados por espacios, alta precisión (%.15e)

**Formato del archivo `total_force.dat`**:
```
# Total force output
# timestep ftotal_x ftotal_y ftotal_z ffluid_x ffluid_y ffluid_z fcoll_x fcoll_y fcoll_z fsubgrid_x fsubgrid_y fsubgrid_z
1000 0.000000000000000e+00 0.000000000000000e+00 1.000000000000000e-04 0.000000000000000e+00 ...
2000 0.000000000000000e+00 0.000000000000000e+00 1.000000000000000e-04 0.000000000000000e+00 ...
```

---

### 2.2. Archivo: `src/stats_velocity.h`

#### Líneas 18-20: Include adicional

```c
/*CHANGE INIT - Total force calculation */
#include "colloids.h"
/*CHANGE END - Total force calculation */
```

#### Líneas 32-38: Declaraciones funciones

```c
/*CHANGE INIT - Total force calculation */
int stats_total_force(hydro_t * hydro, map_t * map, colloids_info_t * cinfo,
                      double ffluid[3], double fcoll[3], double fsubgrid[3],
                      double ftotal[3]);
int stats_total_force_write(hydro_t * hydro, map_t * map, colloids_info_t * cinfo,
                             int timestep, const char * filename);
/*CHANGE END - Total force calculation */
```

---

### 2.3. Archivo: `src/ludwig.c`

#### Líneas 1123-1125: Variables fuerzas

```c
/*CHANGE INIT - Total force calculation */
double ffluid[3], fcoll[3], fsubgrid[3], ftotal[3];
/*CHANGE END - Total force calculation */
```

#### Líneas 1181-1196: Reporte en pantalla

```c
/*CHANGE INIT - Total force calculation */
/* Calculate and report total forces */
stats_total_force(ludwig->hydro, ludwig->map, ludwig->collinfo,
                  ffluid, fcoll, fsubgrid, ftotal);

pe_info(pe, "\n");
pe_info(pe, "Total force - x y z\n");
pe_info(pe, "[total   ] %14.7e %14.7e %14.7e\n", ftotal[X], ftotal[Y], ftotal[Z]);
pe_info(pe, "[fluid   ] %14.7e %14.7e %14.7e\n", ffluid[X], ffluid[Y], ffluid[Z]);
if (ncolloid > 0) {
  pe_info(pe, "[colloids] %14.7e %14.7e %14.7e\n", fcoll[X], fcoll[Y], fcoll[Z]);
}
if (nsubgrid > 0) {
  pe_info(pe, "[subgrid ] %14.7e %14.7e %14.7e\n", fsubgrid[X], fsubgrid[Y], fsubgrid[Z]);
}
/*CHANGE END - Total force calculation */
```

**Output esperado**:
```
Total force - x y z
[total   ]  0.000000e+00  0.000000e+00  1.000000e-04
[fluid   ]  0.000000e+00  0.000000e+00  0.000000e+00
[colloids]  0.000000e+00  0.000000e+00  0.000000e+00
[subgrid ]  0.000000e+00  0.000000e+00  1.000000e-04
```

#### Líneas 916-920: Escritura automática a archivo

```c
/*CHANGE INIT - Total force calculation */
/* Write total force to file */
stats_total_force_write(ludwig->hydro, ludwig->map, ludwig->collinfo,
                        step, "total_force.dat");
/*CHANGE END - Total force calculation */
```

**Se ejecuta cuando**: Se guardan estados de coloides (`config.cds`)
- `is_config_step()` es verdadero, o
- `is_measurement_step()` es verdadero, o
- `is_colloid_io_step()` es verdadero

---

## 3. ARCHIVOS MODIFICADOS (RESUMEN)

| Archivo | Cambios |
|---------|---------|
| `src/stats_colloid.c` | Masa efectiva subgrid + función momentum_subgrid |
| `src/stats_colloid.h` | Declaración stats_colloid_momentum_subgrid() |
| `src/stats_velocity.c` | Funciones stats_total_force() y stats_total_force_write() |
| `src/stats_velocity.h` | Declaraciones funciones de fuerza + include colloids.h |
| `src/ludwig.c` | Integración en reportes + escritura automática |

---

## 4. INTERPRETACIÓN FÍSICA

### 4.1. Conservación de momento

**Momento total del sistema**:
```
p_total = p_fluid + p_colloids + p_subgrid + p_walls
```

donde:
- `p_fluid` = Σ(ρ × u) en nodos fluidos (de distribución LB)
- `p_colloids` = Σ(m × v) coloides resueltos (m calculado de geometría)
- `p_subgrid` = Σ(m_eff × v) coloides subgrid con `m_eff = ρ₀ × V_particula`
- `p_walls` = momento transferido a paredes (contabilidad)

**Verificación**: En ausencia de fuerzas externas, `p_total` debe ser constante.

### 4.2. Suma total de fuerzas

**Fuerza total del sistema**:
```
F_total = F_fluid + F_colloids + F_subgrid
```

donde:
- `F_fluid` = Σ f en nodos fluidos (body forces)
- `F_colloids` = Σ(force + fex) coloides resueltos
- `F_subgrid` = Σ fex coloides subgrid

**Verificación**:
- Sin fuerzas externas: `F_total ≈ 0`
- Con campo eléctrico: `F_total = E × Q_total` (carga total)
- 2ª ley de Newton: `F_total = dp_total/dt`

### 4.3. Relación entre momento y fuerza

```
dp/dt = F_total
```

Por lo tanto:
- Si `F_total ≠ 0` → el momento debe cambiar linealmente
- Si `F_total = 0` → el momento debe conservarse

---

## 5. EJEMPLO DE USO

### 5.1. Verificación en output estándar

Buscar en `stdout` o archivo de log:

```
Momentum - x y z
[total   ]  1.234567e-10  2.345678e-10  3.456789e-10
[fluid   ]  1.000000e-10  2.000000e-10  3.000000e-10
[colloids]  0.000000e+00  0.000000e+00  0.000000e+00
[subgrid ]  2.345670e-11  3.456780e-11  4.567890e-11
[walls   ]  0.000000e+00  0.000000e+00  0.000000e+00

Total force - x y z
[total   ]  0.000000e+00  0.000000e+00  1.000000e-04
[fluid   ]  0.000000e+00  0.000000e+00  0.000000e+00
[colloids]  0.000000e+00  0.000000e+00  0.000000e+00
[subgrid ]  0.000000e+00  0.000000e+00  1.000000e-04
```

### 5.2. Análisis del archivo `total_force.dat`

```python
import numpy as np
import matplotlib.pyplot as plt

# Leer datos
data = np.loadtxt('total_force.dat')
timestep = data[:, 0]
ftotal_z = data[:, 3]  # Componente z de fuerza total

# Graficar
plt.plot(timestep, ftotal_z)
plt.xlabel('Timestep')
plt.ylabel('Total Force Z [LB units]')
plt.title('Total Force Evolution')
plt.savefig('total_force_evolution.png')
```

### 5.3. Verificar conservación

```python
# Leer momento (del output estándar o archivo personalizado)
# Verificar que d(momento)/dt = fuerza_total

momentum = ...  # Leer de archivo
force = data[:, 1:4]  # ftotal_x, ftotal_y, ftotal_z

# Derivada numérica del momento
dt = timestep[1] - timestep[0]
dmom_dt = np.gradient(momentum, dt, axis=0)

# Comparar
residual = dmom_dt - force
print(f"Residual máximo: {np.max(np.abs(residual))}")
```

---

## 6. NOTAS IMPORTANTES

### 6.1. Precisión numérica

- Se usa **Klein summation** en lugar de sumas directas
- Esto reduce errores de redondeo de ~O(ε×N) a ~O(ε²×N)
- Especialmente importante para sumas con ~10⁶ elementos

### 6.2. Partículas subgrid vs resueltas

| Propiedad | Resueltas | Subgrid |
|-----------|-----------|---------|
| Masa | Calculada de geometría | Volumen fluido desplazado |
| Fuerzas usadas | `force + fex` | `fex` solamente |
| Método | BBL (Boundary Bounce-Back) | Nash (IBM-LB) |
| Inercia | Explícita | Implícita en fluido |

### 6.3. Interpretación de `fex` vs `force`

- **`fex`**: Fuerzas externas (ej: eléctricas, gravitatorias)
- **`force`**: Fuerza hidrodinámica del fluido sobre coloide resuelto
- **`fsub`**: En subgrid, es la velocidad interpolada (no fuerza)

Para subgrid:
```
v_new = fsub + drag × fex + thermal_noise
```

donde `drag = 1/(6πηa_h)` (Stokes)

---

## 7. DEBUGGING

Si los resultados son incorrectos:

### 7.1. Verificar momento

```bash
grep "Momentum" output.log | tail -20
```

Verificar:
- ¿El momento total es ~constante sin fuerzas externas?
- ¿La componente subgrid tiene magnitud razonable?

### 7.2. Verificar fuerzas

```bash
grep "Total force" output.log | tail -20
```

Verificar:
- ¿La fuerza total coincide con fuerzas externas aplicadas?
- ¿Las componentes individuales suman correctamente?

### 7.3. Verificar archivo

```bash
tail total_force.dat
wc -l total_force.dat
```

Verificar:
- ¿Se está escribiendo cada `colloid_io_step`?
- ¿Los valores tienen sentido físico?

---

## 8. REFERENCIAS

1. **Nash et al. (2007)**: "Point particle model for Lattice Boltzmann"
   Método IBM-LB para partículas subgrid

2. **Klein summation**: Doubly compensated summation algorithm
   Mejor que Kahan para sumas muy largas

3. **Código Ludwig original**: `src/subgrid.c`, `src/stats_colloid.c`

---

## 9. CHANGELOG

**2025-11-12**: Implementación inicial
- Masa efectiva para partículas subgrid
- Función stats_colloid_momentum_subgrid()
- Función stats_total_force() con Klein summation
- Escritura automática a total_force.dat
- Reporte en pantalla integrado

---

## 10. CONTACTO

Para preguntas o problemas con estos cambios, consultar el código fuente en:
- `src/stats_colloid.c` (líneas 63-70, 184-236)
- `src/stats_velocity.c` (líneas 127-299)
- `src/ludwig.c` (líneas 1106-1199, 916-920)

Todos los cambios están marcados con tags `/*CHANGE INIT - ...*/` y `/*CHANGE END - ...*/`
