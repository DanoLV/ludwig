#!/usr/bin/env python3
"""
Cálculo de los pesos correctos para Laplaciano y Gradiente
en stencils D3Q7, D3Q19 y D3Q27.

Estrategia:
1. Para cada stencil, calcular los pesos que dan el Laplaciano exacto
2. Verificar que esos pesos también den el gradiente exacto
3. Determinar el factor de escala correcto desde los pesos LB base
"""

import numpy as np
from scipy.optimize import least_squares

def test_field_values(cv, point=(1, 1, 1)):
    """
    Valores del campo de prueba ψ(x,y,z) = x² + y² + z²
    en los puntos del stencil centrado en 'point'.
    """
    x0, y0, z0 = point
    psi_values = []

    for p in range(len(cv)):
        x = x0 + cv[p][0]
        y = y0 + cv[p][1]
        z = z0 + cv[p][2]
        psi = x**2 + y**2 + z**2
        psi_values.append(psi)

    return np.array(psi_values)

def compute_laplacian_weights(cv, target_laplacian=6.0):
    """
    Encuentra los pesos w tal que:
    Σ w[p] * ψ(x + c[p]) = target_laplacian

    para el campo de prueba ψ = x² + y² + z²
    """
    n = len(cv)

    # Punto central para la evaluación
    point = (0, 0, 0)  # Usar origen para simetría
    psi_vals = test_field_values(cv, point)

    # Restricciones:
    # 1. Suma de pesos = 0 (requerido para Laplaciano)
    # 2. w[0] = -Σ w[p] para p > 0

    # Construir matriz de restricciones
    # Para simplificar, usamos w[0] = -Σ w[1:]

    def residual(w_nonzero):
        """w_nonzero son los pesos w[1:], w[0] se calcula automáticamente"""
        w = np.zeros(n)
        w[1:] = w_nonzero
        w[0] = -np.sum(w_nonzero)

        # Calcular Laplaciano
        lap = np.dot(w, psi_vals)

        # Queremos que lap = -target_laplacian (signo negativo por convención)
        return lap + target_laplacian

    # Punto inicial: distribución uniforme
    w0 = np.ones(n-1) / (n-1)

    # Resolver
    result = least_squares(residual, w0)

    # Construir pesos completos
    w = np.zeros(n)
    w[1:] = result.x
    w[0] = -np.sum(result.x)

    return w

def compute_gradient_weights(cv):
    """
    Encuentra los pesos w_grad tal que:
    E_x = -Σ w_grad[p] * c_x[p] * ψ(x + c[p]) = -∂ψ/∂x

    Para ψ = x² + y² + z² en (1,1,1):
    ∂ψ/∂x = 2x = 2
    """
    n = len(cv)

    # Punto para evaluación (no usar origen porque gradiente es 0)
    point = (1, 1, 1)
    psi_vals = test_field_values(cv, point)

    # Gradiente esperado en (1,1,1): (2, 2, 2)
    target_grad = np.array([2.0, 2.0, 2.0])

    def residual(w_grad_nonzero):
        """w_grad_nonzero son los pesos w[1:], w[0] = 0 siempre"""
        w_grad = np.zeros(n)
        w_grad[1:] = w_grad_nonzero

        # Calcular gradiente
        grad = np.zeros(3)
        for p in range(n):
            grad[0] -= w_grad[p] * cv[p][0] * psi_vals[p]
            grad[1] -= w_grad[p] * cv[p][1] * psi_vals[p]
            grad[2] -= w_grad[p] * cv[p][2] * psi_vals[p]

        return grad - target_grad

    # Punto inicial
    w0 = np.ones(n-1) * 0.1

    # Resolver
    result = least_squares(residual, w0)

    # Construir pesos completos
    w_grad = np.zeros(n)
    w_grad[1:] = result.x

    return w_grad

def analyze_stencil(name, cv, wv):
    """
    Análisis completo de un stencil:
    - Calcular pesos óptimos de Laplaciano
    - Calcular pesos óptimos de Gradiente
    - Encontrar factores de escala desde pesos LB
    """
    print("\n" + "="*70)
    print(f"ANÁLISIS: {name}")
    print("="*70)

    n = len(cv)

    # Pesos LB base
    print(f"\nPesos LB base (primeros 7):")
    print(f"  {wv[:min(7, n)]}")

    # 1. Calcular pesos óptimos del Laplaciano
    print(f"\n1. LAPLACIANO ÓPTIMO")
    print(f"   Objetivo: Σ w_lap[p] * ψ(x+c[p]) = -6.0 para ψ = x²+y²+z²")

    w_lap_opt = compute_laplacian_weights(cv, target_laplacian=6.0)
    print(f"   Pesos óptimos (primeros 7):")
    print(f"   {w_lap_opt[:min(7, n)]}")

    # Verificar suma = 0
    print(f"   Suma de pesos (debe ser ~0): {np.sum(w_lap_opt):.10e}")

    # Test
    psi_test = test_field_values(cv, point=(0, 0, 0))
    lap_test = np.dot(w_lap_opt, psi_test)
    print(f"   Verificación: ∇²ψ = {lap_test:.10f} (esperado: -6.0)")

    # 2. Determinar factor de escala desde pesos LB
    print(f"\n2. FACTOR DE ESCALA PARA LAPLACIANO")

    # Para cada peso no central, calcular el ratio
    ratios_lap = []
    for p in range(1, min(7, n)):
        if abs(wv[p]) > 1e-10:
            ratio = w_lap_opt[p] / wv[p]
            ratios_lap.append(ratio)
            print(f"   w_lap[{p}] / wv[{p}] = {w_lap_opt[p]:.6f} / {wv[p]:.6f} = {ratio:.2f}")

    if ratios_lap:
        avg_ratio_lap = np.mean(ratios_lap)
        std_ratio_lap = np.std(ratios_lap)
        print(f"\n   Factor promedio: {avg_ratio_lap:.2f} ± {std_ratio_lap:.4f}")

        # Verificar si todos los ratios son consistentes
        if std_ratio_lap < 0.1:
            print(f"   ✓ Factor consistente: usar wlaplacian = {avg_ratio_lap:.1f} * wv")
        else:
            print(f"   ✗ Factores inconsistentes: usar pesos óptimos directamente")

    # 3. Calcular pesos óptimos del Gradiente
    print(f"\n3. GRADIENTE ÓPTIMO")
    print(f"   Objetivo: E_x = -Σ w_grad[p]*c_x[p]*ψ(x+c[p]) = -2.0 en (1,1,1)")

    w_grad_opt = compute_gradient_weights(cv)
    print(f"   Pesos óptimos (primeros 7):")
    print(f"   {w_grad_opt[:min(7, n)]}")

    # Test
    psi_test = test_field_values(cv, point=(1, 1, 1))
    grad_test = np.zeros(3)
    for p in range(n):
        grad_test[0] -= w_grad_opt[p] * cv[p][0] * psi_test[p]
        grad_test[1] -= w_grad_opt[p] * cv[p][1] * psi_test[p]
        grad_test[2] -= w_grad_opt[p] * cv[p][2] * psi_test[p]

    print(f"   Verificación: E = {grad_test} (esperado: [-2, -2, -2])")

    # 4. Determinar factor de escala desde pesos LB
    print(f"\n4. FACTOR DE ESCALA PARA GRADIENTE")

    ratios_grad = []
    for p in range(1, min(7, n)):
        if abs(wv[p]) > 1e-10:
            ratio = w_grad_opt[p] / wv[p]
            ratios_grad.append(ratio)
            print(f"   w_grad[{p}] / wv[{p}] = {w_grad_opt[p]:.6f} / {wv[p]:.6f} = {ratio:.2f}")

    if ratios_grad:
        avg_ratio_grad = np.mean(ratios_grad)
        std_ratio_grad = np.std(ratios_grad)
        print(f"\n   Factor promedio: {avg_ratio_grad:.2f} ± {std_ratio_grad:.4f}")

        if std_ratio_grad < 0.1:
            print(f"   ✓ Factor consistente: usar wgradients = {avg_ratio_grad:.1f} * wv")
        else:
            print(f"   ✗ Factores inconsistentes: usar pesos óptimos directamente")

    return {
        'w_lap_opt': w_lap_opt,
        'w_grad_opt': w_grad_opt,
        'lap_factor': avg_ratio_lap if ratios_lap else None,
        'grad_factor': avg_ratio_grad if ratios_grad else None,
    }

def main():
    """Análisis de los tres stencils"""

    # D3Q7
    cv_d3q7 = np.array([
        [0, 0, 0],
        [1, 0, 0], [0, 1, 0], [0, 0, 1],
        [0, 0,-1], [0,-1, 0], [-1, 0, 0]
    ])
    wv_d3q7 = np.array([2.0/8.0, 1.0/8.0, 1.0/8.0, 1.0/8.0, 1.0/8.0, 1.0/8.0, 1.0/8.0])

    results_d3q7 = analyze_stencil("D3Q7", cv_d3q7, wv_d3q7)

    # D3Q19
    cv_d3q19 = np.array([
        [ 0, 0, 0],
        [ 1, 1, 0], [ 1, 0, 1], [ 1, 0, 0], [ 1, 0,-1], [ 1,-1, 0],
        [ 0, 1, 1], [ 0, 1, 0], [ 0, 1,-1], [ 0, 0, 1], [ 0, 0,-1],
        [ 0,-1, 1], [ 0,-1, 0], [ 0,-1,-1], [-1, 1, 0], [-1, 0, 1],
        [-1, 0, 0], [-1, 0,-1], [-1,-1, 0]
    ])
    wv_d3q19 = np.array([12.0/36.0,
                         1.0/36.0, 1.0/36.0, 2.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0,
                         2.0/36.0, 1.0/36.0, 2.0/36.0, 2.0/36.0, 1.0/36.0, 2.0/36.0,
                         1.0/36.0, 1.0/36.0, 1.0/36.0, 2.0/36.0, 1.0/36.0, 1.0/36.0])

    results_d3q19 = analyze_stencil("D3Q19", cv_d3q19, wv_d3q19)

    # D3Q27
    cv_d3q27 = np.array([
        [ 0, 0, 0],
        [-1,-1,-1], [-1,-1, 0], [-1,-1, 1], [-1, 0,-1], [-1, 0, 0], [-1, 0, 1],
        [-1, 1,-1], [-1, 1, 0], [-1, 1, 1], [ 0,-1,-1], [ 0,-1, 0], [ 0,-1, 1],
        [ 0, 0,-1],                                                 [ 0, 0, 1],
        [ 0, 1,-1], [ 0, 1, 0], [ 0, 1, 1], [ 1,-1,-1], [ 1,-1, 0], [ 1,-1, 1],
        [ 1, 0,-1], [ 1, 0, 0], [ 1, 0, 1], [ 1, 1,-1], [ 1, 1, 0], [ 1, 1, 1]
    ])
    wv_d3q27 = np.array([
        64.0/216.0,
         1.0/216.0,  4.0/216.0,  1.0/216.0,  4.0/216.0, 16.0/216.0,  4.0/216.0,
         1.0/216.0,  4.0/216.0,  1.0/216.0,  4.0/216.0, 16.0/216.0,  4.0/216.0,
        16.0/216.0,                                                  16.0/216.0,
         4.0/216.0, 16.0/216.0,  4.0/216.0,  1.0/216.0,  4.0/216.0,  1.0/216.0,
         4.0/216.0, 16.0/216.0,  4.0/216.0,  1.0/216.0,  4.0/216.0,  1.0/216.0
    ])

    results_d3q27 = analyze_stencil("D3Q27", cv_d3q27, wv_d3q27)

    # Resumen
    print("\n" + "="*70)
    print("RESUMEN DE FACTORES ÓPTIMOS")
    print("="*70)

    print(f"\n{'Stencil':<10} {'Factor Laplaciano':<20} {'Factor Gradiente':<20}")
    print("-"*70)

    for name, results in [("D3Q7", results_d3q7), ("D3Q19", results_d3q19), ("D3Q27", results_d3q27)]:
        lap_f = results['lap_factor']
        grad_f = results['grad_factor']
        print(f"{name:<10} {lap_f:<20.2f} {grad_f:<20.2f}")

    print(f"\nRECOMENDACIÓN:")
    print(f"==============")
    print(f"\nUsar estos factores en los archivos stencil_*.c:")
    print(f"\n  stencil_d3q7.c:   wlaplacian = {results_d3q7['lap_factor']:.1f} * wv")
    print(f"                    wgradients = {results_d3q7['grad_factor']:.1f} * wv")
    print(f"\n  stencil_d3q19.c:  wlaplacian = {results_d3q19['lap_factor']:.1f} * wv")
    print(f"                    wgradients = {results_d3q19['grad_factor']:.1f} * wv")
    print(f"\n  stencil_d3q27.c:  wlaplacian = {results_d3q27['lap_factor']:.1f} * wv")
    print(f"                    wgradients = {results_d3q27['grad_factor']:.1f} * wv")

if __name__ == "__main__":
    main()
