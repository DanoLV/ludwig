# Changelog - runmulti.sh

## Versión 3.0 - Modo paralelo

### Nuevas características importantes

1. **Ejecución en paralelo**
   - Nueva opción `--parallel`: Ejecuta todas las simulaciones simultáneamente
   - Nueva opción `--max-parallel N`: Limita el número de simulaciones concurrentes
   - Cada simulación se ejecuta en su propio directorio independiente
   - Tracking de procesos en background con PIDs
   - Espera inteligente a que todas las simulaciones terminen
   - Reporte de simulaciones exitosas vs fallidas

2. **Gestión de recursos**
   - Control de concurrencia con `--max-parallel`
   - Monitoreo de procesos en background
   - Liberación automática de slots cuando terminan simulaciones

### Comportamiento

**Modo secuencial (default)**:
- Ejecuta una simulación tras otra
- Sale al primer error
- Comportamiento original preservado

**Modo paralelo (`--parallel`)**:
- Lanza todas las simulaciones en background
- Continúa aunque algunas fallen
- Espera a que todas terminen antes de finalizar
- Muestra reporte de éxitos/fallos

**Modo paralelo limitado (`--parallel --max-parallel N`)**:
- Mantiene máximo N simulaciones ejecutándose
- Cuando una termina, lanza la siguiente
- Balance entre velocidad y recursos

### Ejemplos

```bash
# Paralelo ilimitado (¡cuidado con recursos!)
./runmulti.sh --charge 1.0,2.0,3.0 --base-dir test --parallel \
  -i 0 -n 10000 -s 500 -x 32 -y 32 -v 0.5

# Paralelo limitado (recomendado)
./runmulti.sh --charge 1.0,1.5,2.0,2.5,3.0 --base-dir test \
  --parallel --max-parallel 3 \
  -i 0 -n 10000 -s 500 -x 32 -y 32 -v 0.5
```

## Versión 2.1 - Corrección de paso de parámetros

### Correcciones importantes

1. **Paso correcto de parámetros a runbg.sh**
   - Todos los parámetros de runbg.sh ahora se reconocen y pasan correctamente
   - Lista explícita de todas las opciones de runbg.sh en el parser
   - Validación de parámetros desconocidos

### Tests agregados

- `test_param_passing.sh`: Verifica que los parámetros se pasen correctamente

## Versión 2.0 - Actualización con opción --position

### Nuevas características

1. **Opción --position agregada**
   - Permite especificar vectores de posición completos en formato `x_y_z`
   - Ejemplo: `--position 16.0_16.0_16.0,17.0_17.0_17.0,18.0_18.0_18.0`
   - Útil para trayectorias específicas (diagonales, curvas, etc.)

2. **Validación mejorada**
   - El script ahora verifica que `--position` no se use junto con `--position-x/y/z`
   - Mensajes de error más claros

3. **Documentación actualizada**
   - README extendido con ejemplos de uso de `--position`
   - Ejemplos adicionales mostrando cómo variar parámetros del archivo `input` usando wrappers
   - Aclaración sobre cómo variar temperatura, densidad y campo eléctrico

### Parámetros disponibles

#### Modificación de config.cds.init.001-001:
- `--charge VALUES`: Valores de carga (q0)
- `--position VALUES`: Vectores de posición completos (x_y_z)
- `--position-x VALUES`: Solo componente X
- `--position-y VALUES`: Solo componente Y
- `--position-z VALUES`: Solo componente Z
- `--al VALUES`: Parámetro AL

#### Modificación del archivo input (vía runbg.sh):
Para variar estos parámetros, usar las opciones de `runbg.sh` directamente:
- `-k/--temperature`: Temperatura (kT)
- `-r/--rho`: Densidad del fluido
- `-e/--electric-field`: Campo eléctrico (Ex_Ey_Ez)
- `-v/--viscosity`: Viscosidad
- Y todas las demás opciones de runbg.sh

### Ejemplos de uso

```bash
# Variar posición completa (diagonal)
./runmulti.sh \
  --position 16.0_16.0_16.0,17.0_17.0_17.0,18.0_18.0_18.0 \
  --base-dir scan-diagonal \
  -i 0 -n 10000 -s 500 \
  -x 32 -y 32 -v 0.5

# Variar carga y posición X
./runmulti.sh \
  --charge 1.0,2.0 \
  --position-x 16.0,17.0 \
  --base-dir test-q-x \
  -i 0 -n 10000 -s 500 \
  -x 32 -y 32 -v 0.5

# Variar temperatura usando wrapper
for temp in 0.0005 0.001 0.002; do
  ./runmulti.sh \
    --charge 1.0 \
    --base-dir test-temp-${temp} \
    -i 0 -n 10000 -s 500 \
    -x 32 -y 32 -v 0.5 \
    -k $temp \
    -t fe_electro
done
```

### Archivos incluidos

1. `runmulti.sh` - Script principal
2. `README_runmulti.md` - Documentación completa
3. `ejemplo_runmulti.sh` - Ejemplos comentados
4. `test_config_modification.sh` - Test de modificación de config
5. `test_position_option.sh` - Test de opción --position
6. `CHANGELOG_runmulti.md` - Este archivo

### Notas técnicas

- El script mantiene un backup del archivo `config.cds.init.001-001` original
- Para cada simulación, restaura el config original antes de modificarlo
- Los loops están anidados para crear todas las combinaciones posibles (producto cartesiano)
- La opción `--position` tiene prioridad sobre `--position-x/y/z` individuales
- El script valida que no se mezclen ambos métodos de especificar posición

### Compatibilidad

- Totalmente compatible con todas las opciones de `runbg.sh`
- No requiere cambios en archivos existentes
- Preserva el formato del archivo `config.cds.init.001-001`
