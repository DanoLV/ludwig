#!/usr/bin/env python3
"""
Script para determinar el offset correcto entre coordenadas del config.init
y las coordenadas del campo teórico.
"""

import numpy as np
import sys

# Simulación: sistema 64x64x64
# Posición en config.init: 32.20, 32.50, 32.50
# Posición en colloids CSV: 31.70, 32.00, 32.00

config_pos = np.array([32.20, 32.50, 32.50])
csv_pos = np.array([31.70, 32.00, 32.00])

print("="*60)
print("ANÁLISIS DE COORDENADAS EN LUDWIG")
print("="*60)
print(f"\nPosición en config.init: {config_pos}")
print(f"Posición en colloids CSV: {csv_pos}")
print(f"Diferencia: {config_pos - csv_pos}")

# En Ludwig, los índices de la red van de ic=1 a ic=N
# Veamos diferentes interpretaciones:

print("\n" + "="*60)
print("POSIBLES INTERPRETACIONES DEL SISTEMA DE COORDENADAS")
print("="*60)

print("\n1. Nodos en posiciones enteras (0, 1, 2, ..., N-1):")
print("   - config_pos 32.20 significa: nodo 32 está en x=32.0")
print("   - Para campo teórico en nodos (0.5, 1.5, ..., 63.5):")
print(f"     Posición de la carga: {config_pos} (sin cambio)")

print("\n2. Nodos en posiciones medio-enteras (0.5, 1.5, ..., N-0.5):")
print("   - config_pos 32.20 significa: nodo 32 está en x=32.5")
print("   - Para campo teórico en nodos (0.5, 1.5, ..., 63.5):")
print(f"     Posición de la carga: {config_pos} (sin cambio)")

print("\n3. Índices 1-based (ic=1 a ic=N) -> coordenadas (1.0 a N.0):")
print("   - config_pos 32.20 se refiere a la posición física 32.20")
print("   - Nodos físicos en (1.0, 2.0, ..., 64.0)")
print("   - Para campo teórico en nodos (0.5, 1.5, ..., 63.5):")
print(f"     Posición de la carga: {config_pos - 0.5}")

print("\n4. Sistema 0-based con offset implícito:")
print("   - config_pos 32.20 en sistema donde nodos están en (0, 1, ..., 63)")
print("   - Para campo teórico en (0.5, 1.5, ..., 63.5):")
print(f"     Posición de la carga: {config_pos + 0.5}")

print("\n" + "="*60)
print("VERIFICACIÓN CON CÓDIGO LUDWIG")
print("="*60)

print("""
En util.c (línea 492-532):
  - "Lattice sites are assumed to be at integer positions"
  - x0 = r0[X] - floor(r0[X])  // Extrae la parte fraccionaria
  - rsq = pow(1.0 * ic - x0, 2)  // ic son enteros

Esto sugiere:
  - Los NODOS están en posiciones ENTERAS: ic = 0, 1, 2, ..., N-1
  - Una partícula en r = 32.20 está a 0.20 del nodo ic=32
  - La parte fraccionaria x0 = 0.20

En psi_colloid.c (línea 64-67):
  - for (ic = 1; ic <= nlocal[X]; ic++)
  - Los índices van de 1 a N (1-based indexing)

CONCLUSIÓN:
  - Ludwig usa índices 1-based: ic ∈ [1, N]
  - El nodo ic tiene posición física: ic - 0.5 o ic ???

Necesito verificar empíricamente...
""")

print("\n" + "="*60)
print("RECOMENDACIÓN PARA VERIFICACIÓN EXPERIMENTAL")
print("="*60)
print("""
Para determinar el offset correcto:

1. Usar compare_field_theory.py con diferentes posiciones:

   a) Sin offset (asumiendo nodos en enteros):
      --charge-pos 32.20 32.50 32.50

   b) Con offset -0.5 (asumiendo 1-based a 0-based):
      --charge-pos 31.70 32.00 32.00

   c) Con offset +0.5:
      --charge-pos 32.70 33.00 33.00

2. Comparar el error relativo en cada caso.

3. El que tenga menor error es el correcto.
""")
