#!/usr/bin/env python3
"""
Test D3Q27 matrix properties and symmetry
"""

import numpy as np
import sys

# Add path to import the generator
sys.path.insert(0, '/home/bater/Sim/ludwig')
from generate_d3q27_kahan import compute_transformation_matrix, compute_inverse_matrix, cv_d3q27, wv_d3q27

NVEL = 27

def test_symmetry():
    """Test if the velocity set has the required symmetry"""
    print("Testing D3Q27 velocity set symmetry...")

    # Sum of velocities weighted by lattice weights should be zero
    sum_x = np.sum(cv_d3q27[:, 0] * wv_d3q27)
    sum_y = np.sum(cv_d3q27[:, 1] * wv_d3q27)
    sum_z = np.sum(cv_d3q27[:, 2] * wv_d3q27)

    print(f"  Weighted velocity sum:")
    print(f"    X: {sum_x:.20e}")
    print(f"    Y: {sum_y:.20e}")
    print(f"    Z: {sum_z:.20e}")

    if abs(sum_x) < 1e-15 and abs(sum_y) < 1e-15 and abs(sum_z) < 1e-15:
        print("  ✓ Velocity symmetry OK")
    else:
        print("  ✗ Velocity symmetry FAILED")

def test_matrix_orthogonality():
    """Test orthogonality of the transformation matrices"""
    print("\nTesting matrix orthogonality...")

    M = compute_transformation_matrix()
    M_inv = compute_inverse_matrix(M)

    # Test: sum_p (M_inv[p,m] * M[m',p]) = delta[m, m']
    test_matrix = np.zeros((NVEL, NVEL))
    for m in range(NVEL):
        for mp in range(NVEL):
            for p in range(NVEL):
                test_matrix[m, mp] += M_inv[p, m] * M[mp, p]

    error = np.max(np.abs(test_matrix - np.eye(NVEL)))
    print(f"  Orthogonality error: {error:.2e}")

    if error < 1e-10:
        print("  ✓ Matrix orthogonality OK")
    else:
        print("  ✗ Matrix orthogonality FAILED")
        print(f"\n  Identity matrix diagonal: {np.diag(test_matrix)[:5]}")
        print(f"  Off-diagonal max: {np.max(np.abs(test_matrix - np.diag(np.diag(test_matrix)))):.2e}")

def test_momentum_conservation():
    """Test if momentum modes are correctly extracted"""
    print("\nTesting momentum extraction...")

    M = compute_transformation_matrix()
    M_inv = compute_inverse_matrix(M)

    # Create a distribution with zero velocity (equilibrium at rest)
    # f_i = w_i * rho
    rho = 1.0
    f = wv_d3q27 * rho

    # Calculate modes
    modes = np.zeros(NVEL)
    for m in range(NVEL):
        for p in range(NVEL):
            modes[m] += M[m, p] * f[p]

    print(f"  Mode 0 (density): {modes[0]:.15e} (should be {rho:.1f})")
    print(f"  Mode 1 (momentum x): {modes[1]:.15e} (should be 0)")
    print(f"  Mode 2 (momentum y): {modes[2]:.15e} (should be 0)")
    print(f"  Mode 3 (momentum z): {modes[3]:.15e} (should be 0)")

    # Reconstruct distribution
    f_reconstructed = np.zeros(NVEL)
    for p in range(NVEL):
        for m in range(NVEL):
            f_reconstructed[p] += M_inv[p, m] * modes[m]

    reconstruction_error = np.max(np.abs(f_reconstructed - f))
    print(f"  Reconstruction error: {reconstruction_error:.2e}")

    if abs(modes[1]) < 1e-14 and abs(modes[2]) < 1e-14 and abs(modes[3]) < 1e-14:
        print("  ✓ Momentum extraction OK")
    else:
        print("  ✗ Momentum extraction FAILED")
        print(f"    Expected: 0, Got: x={modes[1]:.3e}, y={modes[2]:.3e}, z={modes[3]:.3e}")

def test_mode_structure():
    """Verify the structure of specific modes"""
    print("\nTesting mode structure...")

    M = compute_transformation_matrix()

    # Mode 0 should be all ones
    mode0_correct = np.all(np.abs(M[0, :] - 1.0) < 1e-14)
    print(f"  Mode 0 (density): {'✓ OK' if mode0_correct else '✗ FAILED'}")

    # Modes 1-3 should match velocity components
    mode1_correct = np.allclose(M[1, :], cv_d3q27[:, 0])
    mode2_correct = np.allclose(M[2, :], cv_d3q27[:, 1])
    mode3_correct = np.allclose(M[3, :], cv_d3q27[:, 2])
    print(f"  Mode 1 (momentum x): {'✓ OK' if mode1_correct else '✗ FAILED'}")
    print(f"  Mode 2 (momentum y): {'✓ OK' if mode2_correct else '✗ FAILED'}")
    print(f"  Mode 3 (momentum z): {'✓ OK' if mode3_correct else '✗ FAILED'}")

if __name__ == "__main__":
    print("=" * 70)
    print("D3Q27 Matrix Properties Test")
    print("=" * 70)

    test_symmetry()
    test_matrix_orthogonality()
    test_momentum_conservation()
    test_mode_structure()

    print("\n" + "=" * 70)
    print("Tests completed")
    print("=" * 70)
