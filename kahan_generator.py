#!/usr/bin/env python3
"""
Script para generar el código de Kahan para los modos restantes de d3q19_f2mode_chunk
"""

# Modo 7: t3 positivos en 1,5,6,7,8,11,12,13,14,18 y -r3 negativos en 0,2,3,4,9,10,15,16,17
mode7_t3 = [1,5,6,7,8,11,12,13,14,18]
mode7_r3 = [0,2,3,4,9,10,15,16,17]

# Modo 8: c1 positivos en 6,13 y -c1 negativos en 8,11
mode8_c1_pos = [6,13]
mode8_c1_neg = [8,11]

# Modo 9: t3 positivos en 2,4,6,8,9,10,11,13,15,17 y -r3 negativos en 0,1,3,5,7,12,14,16,18
mode9_t3 = [2,4,6,8,9,10,11,13,15,17]
mode9_r3 = [0,1,3,5,7,12,14,16,18]

# Modo 10: c1 positivos en 2,3,4,6,7,8,11,12,13,15,16,17 y -c2 negativos en 1,5,9,10,14,18
mode10_c1 = [2,3,4,6,7,8,11,12,13,15,16,17]
mode10_c2 = [1,5,9,10,14,18]

# Modo 11: c1 positivos en 2,3,4, c2 positivo en 14,18 y -c2 negativos en 1,5, -c1 negativos en 15,16,17
mode11_c1_pos = [2,3,4]
mode11_c2_pos = [14,18]
mode11_c2_neg = [1,5]
mode11_c1_neg = [15,16,17]

# Modo 12: c1 positivos en 6,7,8, c2 positivos en 5,18, -c1 negativos en 11,12,13, -c2 negativos en 1,14
mode12_c1_pos = [6,7,8]
mode12_c2_pos = [5,18]
mode12_c1_neg = [11,12,13]
mode12_c2_neg = [1,14]

# Modo 13: c1 positivos en 2,6,11,15, c2 positivo en 10, -c1 negativos en 8,13,17, -c2 negativo en 9
mode13_c1_pos = [2,6,11,15]
mode13_c2_pos = [10]
mode13_c1_neg = [8,13,17]
mode13_c2_neg = [9]

# Modo 14: c1 positivos en 3,6,8,11,13, -c1 negativos en 2,4,7,12,15,17
mode14_c1_pos = [3,6,8,11,13]
mode14_c1_neg = [2,4,7,12,15,17]

# Modo 15: c1 positivos en 15,17, -c1 negativos en 2,4,16
mode15_c1_pos = [15,17]
mode15_c1_neg = [2,4,16]

# Modo 16: c1 positivos en 6,8,12, -c1 negativos en 7,11,13
mode16_c1_pos = [6,8,12]
mode16_c1_neg = [7,11,13]

# Modo 17: c1 positivos en 4,6,11,17, -c1 negativos en 2,8,13,15
mode17_c1_pos = [4,6,11,17]
mode17_c1_neg = [2,8,13,15]

# Modo 18: c1 positivos en 0,1,2,4,5,6,8,11,13,14,15,17,18, -c2 negativos en 3,7,9,10,12,16
mode18_c1 = [0,1,2,4,5,6,8,11,13,14,15,17,18]
mode18_c2 = [3,7,9,10,12,16]

print("Modos definidos correctamente")
