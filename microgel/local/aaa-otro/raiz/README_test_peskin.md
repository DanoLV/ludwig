# Test Peskin Delta

Programa para calcular los deltas de Peskin para un punto dado y verificar que la suma de todos los deltas es 1.

## Descripción

Este programa calcula la función delta de Peskin de 4 puntos para una posición dada en el espacio 3D. La función delta de Peskin es una aproximación suave de la función delta de Dirac que se utiliza en métodos de frontera inmersa (Immersed Boundary Method).

### Función Delta de Peskin

La función delta de Peskin en 1D está definida como:

```
δ(r) = {
    1/8 * (3 - 2|r| + √(1 + 4|r| - 4r²))     si |r| ≤ 1
    1/8 * (5 - 2|r| - √(-7 + 12|r| - 4r²))   si 1 < |r| ≤ 2
    0                                         si |r| > 2
}
```

Para 3D, la función delta es el producto de las funciones 1D:
```
δ₃ᴰ(x,y,z) = δ(x) * δ(y) * δ(z)
```

### Propiedades

1. **Soporte compacto**: La función es cero fuera de un rango de 2 unidades de malla.
2. **Normalización**: La suma de todos los deltas en la malla es exactamente 1.
3. **Continuidad**: La función es continua en todo su dominio.

## Compilación

```bash
gcc -o test_peskin_delta test_peskin_delta.c -lm
```

O usando el Makefile:

```bash
make -f Makefile.test_peskin
```

## Uso

```bash
./test_peskin_delta <x> <y> <z>
```

Donde `<x>`, `<y>`, `<z>` son las coordenadas del punto en el espacio continuo.

### Ejemplos

```bash
# Punto en una posición arbitraria
./test_peskin_delta 5.3 5.7 5.2

# Punto en una posición simétrica
./test_peskin_delta 10.0 10.0 10.0

# Punto cerca de un nodo
./test_peskin_delta 8.1 8.1 8.1
```

## Salida

El programa muestra:

1. La posición del punto
2. El rango de nodos afectados (i_min, i_max, etc.)
3. Una tabla con los valores de delta para cada nodo:
   - `dx`, `dy`, `dz`: valores de delta 1D en cada dirección
   - `delta_3D`: producto de los tres deltas 1D
4. La suma total de todos los deltas
5. El error (debe ser ~0)
6. Verificación de que la suma es 1.0

### Ejemplo de salida

```
=======================================================
  Cálculo de deltas de Peskin para un punto
=======================================================

Posición del punto: (5.300000, 5.700000, 5.200000)

Rango de nodos afectados:
  x: [3, 7]
  y: [3, 7]
  z: [3, 7]

Deltas para cada nodo:
-------------------------------------------------------
  i   j   k  |   dx   |   dy   |   dz   | delta_3D
-------------------------------------------------------
  5   5   5 | 0.48508 | 0.48508 | 0.46956 | 0.11048742
  ...
-------------------------------------------------------

Suma total de deltas: 1.000000000000000
Error (debe ser ~0):  0.000000000000000e+00

✓ VERIFICACIÓN CORRECTA: La suma es 1.0
```

## Interpretación

- Cada nodo de la malla recibe una fracción del "peso" del punto.
- Los nodos más cercanos al punto tienen mayor peso.
- La suma de todos los pesos es exactamente 1, lo que garantiza conservación.
- El delta de Peskin de 4 puntos afecta a los nodos en un rango de 2 unidades en cada dirección.

## Referencias

- Peskin, C. S. (2002). "The immersed boundary method". Acta Numerica, 11, 479-517.
- Nash, R. W., Adhikari, R., Tailleur, J., & Cates, M. E. (2010). "Run-and-tumble particles with hydrodynamics: Sedimentation, trapping, and upstream swimming". Physical Review Letters, 104(25), 258101.

## Notas

- La coordenada del nodo se considera en el centro de la celda (i+0.5).
- Solo se muestran los deltas significativos (> 1e-15).
- El error numérico típico es del orden de 1e-16 (precisión de máquina).
