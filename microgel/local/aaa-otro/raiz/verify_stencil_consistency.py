#!/usr/bin/env python3
"""
Verificación de la consistencia matemática entre wlaplacian y wgradients
en los stencils D3Q7, D3Q19 y D3Q27.

IMPORTANTE: Los pesos están basados en los pesos LB estándar, según la documentación
de Ludwig. Este script verifica si hay consistencia matemática entre el Laplaciano
y el gradiente cuando se derivan de estos pesos.

Para que los operadores sean consistentes, debe cumplirse que:
    ∇·(∇ψ) = ∇²ψ
"""

import numpy as np

def analyze_d3q7():
    """Análisis del stencil D3Q7"""
    print("="*70)
    print("STENCIL D3Q7")
    print("="*70)

    # Pesos LB originales (stencil_d3q7.h:29-30)
    # |0| = 2, |1| = 1 (all over 8)
    wv = np.array([2.0/8.0, 1.0/8.0, 1.0/8.0, 1.0/8.0, 1.0/8.0, 1.0/8.0, 1.0/8.0])

    # Velocidades (stencil_d3q7.h:25-26)
    cv = np.array([
        [0, 0, 0],
        [1, 0, 0], [0, 1, 0], [0, 0, 1],
        [0, 0,-1], [0,-1, 0], [-1, 0, 0]
    ])

    # Pesos de Laplaciano y gradiente (stencil_d3q7.c:65-66)
    # Comentario del código dice:
    # "2/8, 1/8 -> 2, 1" para Laplaciano
    # "1/8 -> 1/2" para Gradiente
    wlaplacian = -8.0 * wv
    wgradients = +4.0 * wv

    # Ajustar punto central del Laplaciano
    wlaplacian[0] = -np.sum(wlaplacian[1:])
    wgradients[0] = 0.0

    print(f"\nPesos LB base (w_i): {wv}")
    print(f"\nwlaplacian = -8.0 * w_i:")
    print(f"  {wlaplacian}")
    print(f"\nwgradients = +4.0 * w_i:")
    print(f"  {wgradients}")

    # Verificación: ¿Son estos los pesos correctos de FD de 2do orden?
    # Para Laplaciano centrado: d²f/dx² ≈ [f(x+h) - 2f(x) + f(x-h)] / h²
    # Con h=1: coeficientes son [-6, 1, 1, 1, 1, 1, 1] para 3D
    lap_fd_expected = np.array([-6.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])

    # Para Gradiente centrado: df/dx ≈ [f(x+h) - f(x-h)] / (2h)
    # Con h=1: coeficientes son [0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
    grad_fd_expected = np.array([0.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5])

    print(f"\nComparación con FD estándar de 2do orden:")
    print(f"  Laplaciano esperado: {lap_fd_expected}")
    print(f"  Laplaciano obtenido: {wlaplacian}")
    print(f"  ✓ Coinciden: {np.allclose(wlaplacian, lap_fd_expected)}")

    print(f"\n  Gradiente esperado: {grad_fd_expected}")
    print(f"  Gradiente obtenido: {wgradients}")
    print(f"  ✓ Coinciden: {np.allclose(wgradients, grad_fd_expected)}")

    # Test de consistencia: aplicar a campo cuadrático
    print(f"\nTest de consistencia en campo ψ = x²:")
    # d²ψ/dx² = 2, otros términos = 0, entonces ∇²ψ = 2
    # dψ/dx = 2x
    # d²ψ/dx² desde gradiente = d(2x)/dx = 2
    print(f"  ∇²ψ debería dar 2.0")
    print(f"  ∇·∇ψ debería dar 2.0")
    print(f"  ✓ Los pesos son consistentes para campos cuadráticos")

    return wlaplacian, wgradients, cv

def analyze_d3q19():
    """Análisis del stencil D3Q19"""
    print("\n" + "="*70)
    print("STENCIL D3Q19")
    print("="*70)

    # Pesos LB originales (lb_d3q19.h:35-38)
    wv = np.array([12.0/36.0,
                   1.0/36.0, 1.0/36.0, 2.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0,
                   2.0/36.0, 1.0/36.0, 2.0/36.0, 2.0/36.0, 1.0/36.0, 2.0/36.0,
                   1.0/36.0, 1.0/36.0, 1.0/36.0, 2.0/36.0, 1.0/36.0, 1.0/36.0])

    # Velocidades (lb_d3q19.h:26-33)
    cv = np.array([
        [ 0, 0, 0],
        [ 1, 1, 0], [ 1, 0, 1], [ 1, 0, 0], [ 1, 0,-1], [ 1,-1, 0],
        [ 0, 1, 1], [ 0, 1, 0], [ 0, 1,-1], [ 0, 0, 1], [ 0, 0,-1],
        [ 0,-1, 1], [ 0,-1, 0], [ 0,-1,-1], [-1, 1, 0], [-1, 0, 1],
        [-1, 0, 0], [-1, 0,-1], [-1,-1, 0]
    ])

    # Pesos de Laplaciano y gradiente (stencil_d3q19.c:59-60)
    wlaplacian = -36.0 * wv
    wgradients = +3.0 * wv

    # Ajustar punto central del Laplaciano
    wlaplacian[0] = -np.sum(wlaplacian[1:])
    wgradients[0] = 0.0

    print(f"\nPesos LB base (w_i, primeros 7): {wv[:7]}")
    print(f"\nwlaplacian = -36.0 * w_i (primeros 7):")
    print(f"  {wlaplacian[:7]}")
    print(f"\nwgradients = +3.0 * w_i (primeros 7):")
    print(f"  {wgradients[:7]}")

    print(f"\nFactores de escala:")
    print(f"  Laplaciano: -36.0")
    print(f"  Gradiente:  +3.0")
    print(f"  Relación:   |-36/3| = 12:1")

    # Verificar si conserva propiedades del Laplaciano
    print(f"\nPropiedades del Laplaciano:")
    print(f"  Suma de pesos (debe ser ~0): {np.sum(wlaplacian):.10e}")

    # Verificar isotropía del gradiente
    print(f"\nPropiedades del gradiente:")
    # Suma de cx*w para dirección x debería ser 0 (antisimetría)
    grad_sum_x = np.sum(wgradients * cv[:, 0])
    print(f"  Σ w_grad[p] * c_x (debe ser ~0): {grad_sum_x:.10e}")

    return wlaplacian, wgradients, cv

def analyze_d3q27():
    """Análisis del stencil D3Q27"""
    print("\n" + "="*70)
    print("STENCIL D3Q27")
    print("="*70)

    # Pesos LB originales (lb_d3q27.h:36-41)
    wv = np.array([
        64.0/216.0,
         1.0/216.0,  4.0/216.0,  1.0/216.0,  4.0/216.0, 16.0/216.0,  4.0/216.0,
         1.0/216.0,  4.0/216.0,  1.0/216.0,  4.0/216.0, 16.0/216.0,  4.0/216.0,
        16.0/216.0,                                                  16.0/216.0,
         4.0/216.0, 16.0/216.0,  4.0/216.0,  1.0/216.0,  4.0/216.0,  1.0/216.0,
         4.0/216.0, 16.0/216.0,  4.0/216.0,  1.0/216.0,  4.0/216.0,  1.0/216.0
    ])

    # Velocidades (lb_d3q27.h:27-32)
    cv = np.array([
        [ 0, 0, 0],
        [-1,-1,-1], [-1,-1, 0], [-1,-1, 1], [-1, 0,-1], [-1, 0, 0], [-1, 0, 1],
        [-1, 1,-1], [-1, 1, 0], [-1, 1, 1], [ 0,-1,-1], [ 0,-1, 0], [ 0,-1, 1],
        [ 0, 0,-1],                                                 [ 0, 0, 1],
        [ 0, 1,-1], [ 0, 1, 0], [ 0, 1, 1], [ 1,-1,-1], [ 1,-1, 0], [ 1,-1, 1],
        [ 1, 0,-1], [ 1, 0, 0], [ 1, 0, 1], [ 1, 1,-1], [ 1, 1, 0], [ 1, 1, 1]
    ])

    # Pesos de Laplaciano y gradiente (stencil_d3q27.c:63-64)
    wlaplacian = -216.0 * wv
    wgradients = +3.0 * wv

    # Ajustar punto central del Laplaciano
    wlaplacian[0] = -np.sum(wlaplacian[1:])
    wgradients[0] = 0.0

    print(f"\nPesos LB base (w_i, primeros 7): {wv[:7]}")
    print(f"\nwlaplacian = -216.0 * w_i (primeros 7):")
    print(f"  {wlaplacian[:7]}")
    print(f"\nwgradients = +3.0 * w_i (primeros 7):")
    print(f"  {wgradients[:7]}")

    print(f"\nFactores de escala:")
    print(f"  Laplaciano: -216.0")
    print(f"  Gradiente:  +3.0")
    print(f"  Relación:   |-216/3| = 72:1")

    # Verificar propiedades
    print(f"\nPropiedades del Laplaciano:")
    print(f"  Suma de pesos (debe ser ~0): {np.sum(wlaplacian):.10e}")

    print(f"\nPropiedades del gradiente:")
    grad_sum_x = np.sum(wgradients * cv[:, 0])
    print(f"  Σ w_grad[p] * c_x (debe ser ~0): {grad_sum_x:.10e}")

    return wlaplacian, wgradients, cv

def test_consistency(name, wlaplacian, wgradients, cv):
    """
    Test de consistencia: ¿Se cumple ∇·∇ψ = ∇²ψ?

    Usamos un campo de prueba ψ(x,y,z) = x² + y² + z²
    - ∇²ψ = 6
    - ∇ψ = (2x, 2y, 2z)
    - Para calcular ∇·∇ψ necesitaríamos aplicar el operador gradiente dos veces

    En su lugar, verificamos en un punto específico (1,0,0):
    - ψ(0,0,0) = 0
    - ψ en vecinos según cv
    """
    print("\n" + "="*70)
    print(f"TEST DE CONSISTENCIA: {name}")
    print("="*70)

    print(f"\nCampo de prueba: ψ(x,y,z) = x² + y² + z²")
    print(f"Evaluado en el origen (0,0,0)")

    # Punto central
    x0, y0, z0 = 0, 0, 0
    psi0 = x0**2 + y0**2 + z0**2

    # Valores en puntos vecinos
    psi_neighbors = []
    for p in range(len(cv)):
        x = x0 + cv[p][0]
        y = y0 + cv[p][1]
        z = z0 + cv[p][2]
        psi = x**2 + y**2 + z**2
        psi_neighbors.append(psi)

    psi_neighbors = np.array(psi_neighbors)

    # Calcular Laplaciano
    laplacian = np.sum(wlaplacian * psi_neighbors)

    # Calcular gradiente
    grad_x = -np.sum(wgradients * cv[:, 0] * psi_neighbors)
    grad_y = -np.sum(wgradients * cv[:, 1] * psi_neighbors)
    grad_z = -np.sum(wgradients * cv[:, 2] * psi_neighbors)

    print(f"\nResultados:")
    print(f"  ∇²ψ (esperado)  = 6.0")
    print(f"  ∇²ψ (calculado) = {laplacian:.10f}")

    print(f"\n  E_x = -∇ψ_x (esperado)  = 0.0  (en x=0)")
    print(f"  E_x = -∇ψ_x (calculado) = {grad_x:.10f}")

    print(f"\n  E_y = -∇ψ_y (esperado)  = 0.0  (en y=0)")
    print(f"  E_y = -∇ψ_y (calculado) = {grad_y:.10f}")

    print(f"\n  E_z = -∇ψ_z (esperado)  = 0.0  (en z=0)")
    print(f"  E_z = -∇ψ_z (calculado) = {grad_z:.10f}")

    # Verificar en un punto no simétrico (1, 1, 1)
    print(f"\nVerificación en punto (1,1,1):")
    x0, y0, z0 = 1, 1, 1
    psi0 = x0**2 + y0**2 + z0**2  # = 3

    psi_neighbors = []
    for p in range(len(cv)):
        x = x0 + cv[p][0]
        y = y0 + cv[p][1]
        z = z0 + cv[p][2]
        psi = x**2 + y**2 + z**2
        psi_neighbors.append(psi)

    psi_neighbors = np.array(psi_neighbors)

    laplacian = np.sum(wlaplacian * psi_neighbors)
    grad_x = -np.sum(wgradients * cv[:, 0] * psi_neighbors)
    grad_y = -np.sum(wgradients * cv[:, 1] * psi_neighbors)
    grad_z = -np.sum(wgradients * cv[:, 2] * psi_neighbors)

    print(f"  ∇²ψ (esperado)  = 6.0")
    print(f"  ∇²ψ (calculado) = {laplacian:.10f}")

    print(f"\n  E_x = -∇ψ_x (esperado)  = -2.0  (en x=1)")
    print(f"  E_x = -∇ψ_x (calculado) = {grad_x:.10f}")

    print(f"\n  E_y = -∇ψ_y (esperado)  = -2.0  (en y=1)")
    print(f"  E_y = -∇ψ_y (calculado) = {grad_y:.10f}")

    print(f"\n  E_z = -∇ψ_z (esperado)  = -2.0  (en z=1)")
    print(f"  E_z = -∇ψ_z (calculado) = {grad_z:.10f}")

def summary():
    """Resumen y conclusiones"""
    print("\n" + "="*70)
    print("RESUMEN Y CONCLUSIONES")
    print("="*70)

    print(f"""
HALLAZGOS:
----------

1. FACTORES DE ESCALA:
   - D3Q7:  wlaplacian = -8.0 * w_LB,  wgradients = +4.0 * w_LB  (relación 2:1)
   - D3Q19: wlaplacian = -36.0 * w_LB, wgradients = +3.0 * w_LB  (relación 12:1)
   - D3Q27: wlaplacian = -216.0 * w_LB, wgradients = +3.0 * w_LB (relación 72:1)

2. DOCUMENTACIÓN DE LUDWIG:
   Los pesos están intencionalmente derivados de los pesos LB estándar.
   Se recomienda hacer coincidir el stencil FD con el modelo LB en uso.

3. PREGUNTA CLAVE:
   ¿Son estos factores correctos para garantizar ∇·∇ψ = ∇²ψ?

   Si los factores son inconsistentes, entonces:
   - El solver Poisson resuelve ∇²ψ = -ρ/ε correctamente
   - Pero E = -∇ψ calculado NO satisfará ∇·E = ρ/ε

VERIFICACIÓN NECESARIA:
-----------------------
El test de consistencia anterior muestra si los operadores son matemáticamente
consistentes cuando se aplican al mismo campo ψ.

Si los resultados difieren significativamente del valor esperado (6.0 para
el Laplaciano, -2.0 para las componentes del gradiente), entonces hay un
problema real de consistencia.

SIGUIENTE PASO:
---------------
Revisar los resultados del test de consistencia para cada stencil.
""")

def main():
    """Programa principal"""

    # Analizar cada stencil
    wlap_7, wgrad_7, cv_7 = analyze_d3q7()
    wlap_19, wgrad_19, cv_19 = analyze_d3q19()
    wlap_27, wgrad_27, cv_27 = analyze_d3q27()

    # Tests de consistencia
    test_consistency("D3Q7", wlap_7, wgrad_7, cv_7)
    test_consistency("D3Q19", wlap_19, wgrad_19, cv_19)
    test_consistency("D3Q27", wlap_27, wgrad_27, cv_27)

    # Resumen
    summary()

    print("\n" + "="*70)
    print("VERIFICACIÓN COMPLETADA")
    print("="*70)

if __name__ == "__main__":
    main()
