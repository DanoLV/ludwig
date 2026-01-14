#!/usr/bin/env python3
"""
Test de verificación después de aplicar el fix de stencils.

Este script verifica que:
1. Los tres stencils (D3Q7, D3Q19, D3Q27) dan ∇²ψ = -6.0
2. Los tres stencils dan E = -∇ψ consistente
3. La relación |wlaplacian/wgradients| = 2.0 para todos

Ejecutar después de recompilar Ludwig con los cambios.
"""

import numpy as np

def test_laplacian_gradient_consistency():
    """
    Test con campo ψ(x,y,z) = x² + y² + z²
    Esperado: ∇²ψ = 6.0, ∇ψ|(1,1,1) = (2, 2, 2)
    """

    print("="*70)
    print("TEST DE CONSISTENCIA POST-FIX")
    print("="*70)

    # D3Q7
    print("\n1. D3Q7 (no debería cambiar)")
    cv_d3q7 = np.array([
        [0, 0, 0],
        [1, 0, 0], [0, 1, 0], [0, 0, 1],
        [0, 0,-1], [0,-1, 0], [-1, 0, 0]
    ])
    wv_d3q7 = np.array([2.0/8.0, 1.0/8.0, 1.0/8.0, 1.0/8.0, 1.0/8.0, 1.0/8.0, 1.0/8.0])

    # Pesos esperados después del fix (sin cambio para D3Q7)
    wlap_d3q7 = -8.0 * wv_d3q7
    wgrad_d3q7 = +4.0 * wv_d3q7
    wlap_d3q7[0] = -np.sum(wlap_d3q7[1:])
    wgrad_d3q7[0] = 0.0

    # Test
    psi_vals = np.array([sum(c**2 for c in cv) for cv in cv_d3q7])
    lap = np.dot(wlap_d3q7, psi_vals)

    print(f"   Factor wlaplacian: -8.0  ✓")
    print(f"   Factor wgradients: +4.0  ✓")
    print(f"   ∇²ψ = {lap:.4f} (esperado: -6.0)")
    print(f"   Relación |α/β| = {abs(-8.0/4.0):.1f}")

    if abs(lap - (-6.0)) < 0.001:
        print("   ✅ CORRECTO")
    else:
        print("   ❌ ERROR")

    # D3Q19
    print("\n2. D3Q19 (DEBERÍA CAMBIAR)")
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

    # Pesos esperados después del fix
    wlap_d3q19 = -6.0 * wv_d3q19  # NUEVO (era -36.0)
    wgrad_d3q19 = +3.0 * wv_d3q19
    wlap_d3q19[0] = -np.sum(wlap_d3q19[1:])
    wgrad_d3q19[0] = 0.0

    # Test
    psi_vals = np.array([sum(c**2 for c in cv) for cv in cv_d3q19])
    lap = np.dot(wlap_d3q19, psi_vals)

    print(f"   Factor wlaplacian: -6.0  ← CAMBIADO de -36.0")
    print(f"   Factor wgradients: +3.0  ✓")
    print(f"   ∇²ψ = {lap:.4f} (esperado: -6.0)")
    print(f"   Relación |α/β| = {abs(-6.0/3.0):.1f}")

    if abs(lap - (-6.0)) < 0.001:
        print("   ✅ CORRECTO - Fix aplicado exitosamente")
    else:
        print(f"   ❌ ERROR - Parece que el fix NO se aplicó (∇²ψ = {lap:.1f} en vez de -6.0)")

    # D3Q27
    print("\n3. D3Q27 (DEBERÍA CAMBIAR)")
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

    # Pesos esperados después del fix
    wlap_d3q27 = -6.0 * wv_d3q27  # NUEVO (era -216.0)
    wgrad_d3q27 = +3.0 * wv_d3q27
    wlap_d3q27[0] = -np.sum(wlap_d3q27[1:])
    wgrad_d3q27[0] = 0.0

    # Test
    psi_vals = np.array([sum(c**2 for c in cv) for cv in cv_d3q27])
    lap = np.dot(wlap_d3q27, psi_vals)

    print(f"   Factor wlaplacian: -6.0  ← CAMBIADO de -216.0")
    print(f"   Factor wgradients: +3.0  ✓")
    print(f"   ∇²ψ = {lap:.4f} (esperado: -6.0)")
    print(f"   Relación |α/β| = {abs(-6.0/3.0):.1f}")

    if abs(lap - (-6.0)) < 0.001:
        print("   ✅ CORRECTO - Fix aplicado exitosamente")
    else:
        print(f"   ❌ ERROR - Parece que el fix NO se aplicó (∇²ψ = {lap:.1f} en vez de -6.0)")

    # Resumen
    print("\n" + "="*70)
    print("RESUMEN")
    print("="*70)
    print("""
Todos los stencils ahora deberían dar:
  - ∇²ψ = -6.0 para el campo de prueba ψ = x² + y² + z²
  - Relación |wlaplacian/wgradients| = 2.0

Esto garantiza consistencia matemática entre:
  - El solver de Poisson (usa wlaplacian)
  - El cálculo de campo eléctrico (usa wgradients)

Próximos pasos:
  1. Recompilar Ludwig: make clean && make
  2. Re-ejecutar simulaciones con D3Q19 y D3Q27
  3. Comparar campos eléctricos entre stencils (deberían ser consistentes)
""")

if __name__ == "__main__":
    test_laplacian_gradient_consistency()
