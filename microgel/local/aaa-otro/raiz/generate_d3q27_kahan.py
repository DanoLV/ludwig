#!/usr/bin/env python3
"""
Generator for D3Q27 Kahan/Klein summation functions

This script generates optimized d3q27_f2mode_chunk and d3q27_mode2f_chunk_kahan
functions with compensated summation to reduce roundoff errors.

Based on the D3Q19 implementation in collision.c

CHANGE INIT - Kahan summation for D3Q27
This file was auto-generated to create Kahan compensated summation versions
of the mode transformations for D3Q27, analogous to the D3Q19 implementation.
CHANGE END
"""

import numpy as np

# D3Q27 velocity vectors (from lb_d3q27.h)
cv_d3q27 = np.array([
    [ 0, 0, 0],  # 0
    [-1,-1,-1],  # 1
    [-1,-1, 0],  # 2
    [-1,-1, 1],  # 3
    [-1, 0,-1],  # 4
    [-1, 0, 0],  # 5
    [-1, 0, 1],  # 6
    [-1, 1,-1],  # 7
    [-1, 1, 0],  # 8
    [-1, 1, 1],  # 9
    [ 0,-1,-1],  # 10
    [ 0,-1, 0],  # 11
    [ 0,-1, 1],  # 12
    [ 0, 0,-1],  # 13
    [ 0, 0, 1],  # 14
    [ 0, 1,-1],  # 15
    [ 0, 1, 0],  # 16
    [ 0, 1, 1],  # 17
    [ 1,-1,-1],  # 18
    [ 1,-1, 0],  # 19
    [ 1,-1, 1],  # 20
    [ 1, 0,-1],  # 21
    [ 1, 0, 0],  # 22
    [ 1, 0, 1],  # 23
    [ 1, 1,-1],  # 24
    [ 1, 1, 0],  # 25
    [ 1, 1, 1],  # 26
])

# D3Q27 weights (from lb_d3q27.h)
wv_d3q27 = np.array([
    64.0/216.0,  # 0
     1.0/216.0, 4.0/216.0,  1.0/216.0,  # 1-3
     4.0/216.0, 16.0/216.0, 4.0/216.0,  # 4-6
     1.0/216.0, 4.0/216.0,  1.0/216.0,  # 7-9
     4.0/216.0, 16.0/216.0, 4.0/216.0,  # 10-12
    16.0/216.0, 16.0/216.0,              # 13-14
     4.0/216.0, 16.0/216.0, 4.0/216.0,  # 15-17
     1.0/216.0, 4.0/216.0,  1.0/216.0,  # 18-20
     4.0/216.0, 16.0/216.0, 4.0/216.0,  # 21-23
     1.0/216.0, 4.0/216.0,  1.0/216.0   # 24-26
])

NVEL = 27
cs2 = 1.0/3.0

def compute_transformation_matrix():
    """Compute the M matrix for D3Q27"""
    M = np.zeros((NVEL, NVEL))

    for p in range(NVEL):
        cx = cv_d3q27[p, 0]
        cy = cv_d3q27[p, 1]
        cz = cv_d3q27[p, 2]

        M[ 0, p] = 1.0
        M[ 1, p] = cx
        M[ 2, p] = cy
        M[ 3, p] = cz
        M[ 4, p] = cx*cx - cs2
        M[ 5, p] = cx*cy
        M[ 6, p] = cx*cz
        M[ 7, p] = cy*cy - cs2
        M[ 8, p] = cy*cz
        M[ 9, p] = cz*cz - cs2
        M[10, p] = 3.0*(cx*cx - cs2)*cy
        M[11, p] = 3.0*(cx*cx - cs2)*cz
        M[12, p] = 3.0*(cy*cy - cs2)*cz
        M[13, p] = 3.0*(cy*cy - cs2)*cx
        M[14, p] = 3.0*(cz*cz - cs2)*cx
        M[15, p] = 3.0*(cz*cz - cs2)*cy
        M[16, p] = cx*cy*cz
        M[17, p] = 9.0*(cx*cx - cs2)*(cy*cy - cs2)
        M[18, p] = 9.0*(cy*cy - cs2)*(cz*cz - cs2)
        M[19, p] = 9.0*(cz*cz - cs2)*(cx*cx - cs2)
        M[20, p] = 9.0*(cx*cx - cs2)*cy*cz
        M[21, p] = 9.0*(cy*cy - cs2)*cz*cx
        M[22, p] = 9.0*(cz*cz - cs2)*cx*cy
        M[23, p] = 9.0*(cx*cx - cs2)*(cy*cy - cs2)*cz
        M[24, p] = 9.0*(cy*cy - cs2)*(cz*cz - cs2)*cx
        M[25, p] = 9.0*(cz*cz - cs2)*(cx*cx - cs2)*cy
        M[26, p] = 27.0*(cx*cx - cs2)*(cy*cy - cs2)*(cz*cz - cs2)

    return M

def compute_inverse_matrix(M):
    """
    Compute M^{-1} using Ludwig's method:
    mi[p][m] = wv[p] * na[m] * ma[m][p]
    where na[m] = 1 / sum_p(wv[p] * ma[m][p]^2)

    This is the correct formula used in Ludwig for D3Q19.
    """
    # CHANGE INIT - Use Ludwig's formula correctly
    # Compute normalisers: na[m] = 1 / sum_p(wv[p] * ma[m][p]^2)
    na = np.zeros(NVEL)
    for m in range(NVEL):
        sum_val = 0.0
        for p in range(NVEL):
            sum_val += wv_d3q27[p] * M[m, p] * M[m, p]
        na[m] = 1.0 / sum_val

    # Compute inverse: mi[p][m] = wv[p] * na[m] * ma[m][p]
    M_inv = np.zeros((NVEL, NVEL))
    for p in range(NVEL):
        for m in range(NVEL):
            M_inv[p, m] = wv_d3q27[p] * na[m] * M[m, p]
    # CHANGE END

    return M_inv

def format_coefficient(val, tolerance=1e-14):
    """Format coefficient as C code constant"""
    if abs(val) < tolerance:
        return "c0"
    elif abs(val - 1.0) < tolerance:
        return "c1"
    elif abs(val + 1.0) < tolerance:
        return "-c1"
    else:
        # Return as explicit constant
        return f"{val:.18e}"

def generate_f2mode_chunk(M):
    """Generate d3q27_f2mode_chunk function with Kahan summation

    For momentum modes (m=1,2,3), sum positive contributions first, then negative,
    to match D3Q19's ordering which achieves exact zero.
    """

    code = []
    code.append("/* CHANGE INIT - D3Q27 Kahan summation for f2mode */")
    code.append("__device__ void d3q27_f2mode_chunk(double* mode, const double* __restrict__ fchunk)")
    code.append("{")
    code.append("  int m, iv;")
    code.append("")
    code.append("  /* Initialize all modes to zero */")
    code.append("  for (m = 0; m < NVEL; m++) {")
    code.append("    for_simd_v(iv, NSIMDVL) mode[m * NSIMDVL + iv] = 0.0;")
    code.append("  }")
    code.append("")

    # Mode 0 (density) - no Kahan needed (all coefficients are 1.0)
    code.append("  /* m=0 - density (no compensated summation needed) */")
    for p in range(NVEL):
        coeff = format_coefficient(M[0, p])
        code.append(f"  for_simd_v(iv, NSIMDVL) mode[0 * NSIMDVL + iv] += fchunk[{p} * NSIMDVL + iv] * {coeff};")
    code.append("")

    # Modes 1-26 with Kahan summation
    for m in range(1, NVEL):
        code.append(f"  /* m={m} */")
        code.append("  {")
        code.append(f"    volatile double mode{m}_sum[NSIMDVL];")
        code.append(f"    volatile double mode{m}_c[NSIMDVL];")
        code.append("    for_simd_v(iv, NSIMDVL) {")
        code.append(f"      mode{m}_sum[iv] = 0.0;")
        code.append(f"      mode{m}_c[iv] = 0.0;")
        code.append("    }")
        code.append("")

        # For momentum modes (1, 2, 3), group by sign to match D3Q19 ordering
        if m in [1, 2, 3]:
            # Momentum mode - sum positive first, then negative
            code.append(f"    /* Momentum mode {m} - sum positive then negative like D3Q19 */")

            # Collect positive and negative contributions
            positive_terms = []
            negative_terms = []
            for p in range(NVEL):
                coeff = M[m, p]
                if abs(coeff) > 1e-14:  # Skip zeros
                    if coeff > 0:
                        positive_terms.append((p, coeff))
                    else:
                        negative_terms.append((p, coeff))

            # Add positive contributions first
            for p, coeff in positive_terms:
                coeff_str = format_coefficient(coeff)
                code.append("    for_simd_v(iv, NSIMDVL) {")
                code.append(f"      volatile double val = fchunk[{p} * NSIMDVL + iv] * {coeff_str};")
                code.append(f"      volatile double y = val - mode{m}_c[iv];")
                code.append(f"      volatile double t = mode{m}_sum[iv] + y;")
                code.append(f"      mode{m}_c[iv] = (t - mode{m}_sum[iv]) - y;")
                code.append(f"      mode{m}_sum[iv] = t;")
                code.append("    }")

            # Then add negative contributions
            for p, coeff in negative_terms:
                coeff_str = format_coefficient(coeff)
                code.append("    for_simd_v(iv, NSIMDVL) {")
                code.append(f"      volatile double val = fchunk[{p} * NSIMDVL + iv] * {coeff_str};")
                code.append(f"      volatile double y = val - mode{m}_c[iv];")
                code.append(f"      volatile double t = mode{m}_sum[iv] + y;")
                code.append(f"      mode{m}_c[iv] = (t - mode{m}_sum[iv]) - y;")
                code.append(f"      mode{m}_sum[iv] = t;")
                code.append("    }")
        else:
            # Non-momentum modes - keep original order
            for p in range(NVEL):
                coeff = M[m, p]
                if abs(coeff) > 1e-14:  # Skip zeros
                    coeff_str = format_coefficient(coeff)
                    code.append("    for_simd_v(iv, NSIMDVL) {")
                    code.append(f"      volatile double val = fchunk[{p} * NSIMDVL + iv] * {coeff_str};")
                    code.append(f"      volatile double y = val - mode{m}_c[iv];")
                    code.append(f"      volatile double t = mode{m}_sum[iv] + y;")
                    code.append(f"      mode{m}_c[iv] = (t - mode{m}_sum[iv]) - y;")
                    code.append(f"      mode{m}_sum[iv] = t;")
                    code.append("    }")

        code.append("")
        code.append(f"    for_simd_v(iv, NSIMDVL) mode[{m} * NSIMDVL + iv] = mode{m}_sum[iv];")
        code.append("  }")
        code.append("")

    code.append("}")
    code.append("/* CHANGE END */")
    code.append("")

    return "\n".join(code)

def generate_mode2f_chunk_kahan(M_inv):
    """Generate d3q27_mode2f_chunk_kahan function"""

    code = []
    code.append("/* CHANGE INIT - D3Q27 Kahan summation for mode2f */")
    code.append("__device__ void d3q27_mode2f_chunk_kahan(double* mode, double* fchunk) {")
    code.append("")
    code.append("  double ftmp[NSIMDVL];")
    code.append("  double ftmp_c[NSIMDVL];  /* Kahan compensation */")
    code.append("  int iv;")
    code.append("")
    code.append("  /* Helper macro for Kahan summation */")
    code.append("#define KAHAN_ADD(coeff, mode_idx) \\")
    code.append("    for_simd_v(iv, NSIMDVL) { \\")
    code.append("      volatile double val = (coeff) * mode[(mode_idx) * NSIMDVL + iv]; \\")
    code.append("      volatile double y = val - ftmp_c[iv]; \\")
    code.append("      volatile double t = ftmp[iv] + y; \\")
    code.append("      ftmp_c[iv] = (t - ftmp[iv]) - y; \\")
    code.append("      ftmp[iv] = t; \\")
    code.append("    }")
    code.append("")

    # Generate code for each velocity
    for p in range(NVEL):
        code.append(f"  /* p={p} */")
        code.append("  for_simd_v(iv, NSIMDVL) { ftmp[iv] = 0.0; ftmp_c[iv] = 0.0; }")

        for m in range(NVEL):
            coeff = M_inv[p, m]
            if abs(coeff) > 1e-14:  # Skip zeros
                coeff_str = format_coefficient(coeff)
                code.append(f"  KAHAN_ADD({coeff_str}, {m});")

        code.append(f"  for_simd_v(iv, NSIMDVL) fchunk[{p} * NSIMDVL + iv] = ftmp[iv];")
        code.append("")

    code.append("#undef KAHAN_ADD")
    code.append("}")
    code.append("/* CHANGE END */")
    code.append("")

    return "\n".join(code)

def main():
    print("Generating D3Q27 transformation matrices...")
    M = compute_transformation_matrix()
    M_inv = compute_inverse_matrix(M)

    # Verify: The correct relation is M_inv @ M should give weighted identity
    # or equivalently, sum_p (mi[p][m] * ma[m'][p]) should equal delta[m, m']
    test_matrix = np.zeros((NVEL, NVEL))
    for m in range(NVEL):
        for mp in range(NVEL):
            for p in range(NVEL):
                test_matrix[m, mp] += M_inv[p, m] * M[mp, p]

    error = np.max(np.abs(test_matrix - np.eye(NVEL)))
    print(f"Matrix inversion error: {error:.2e}")

    if error > 1e-10:
        print(f"WARNING: Large inversion error! Check the formula.")
        # Print some diagnostic info
        print(f"M shape: {M.shape}")
        print(f"M_inv shape: {M_inv.shape}")
        print(f"M[0,:] = {M[0,:]}")
        print(f"M_inv[:,0] = {M_inv[:,0]}")

    print("\nGenerating d3q27_f2mode_chunk...")
    f2mode_code = generate_f2mode_chunk(M)

    print("Generating d3q27_mode2f_chunk_kahan...")
    mode2f_code = generate_mode2f_chunk_kahan(M_inv)

    # Write to file
    output_file = "/home/bater/Sim/ludwig/d3q27_kahan_functions.c"
    with open(output_file, 'w') as f:
        f.write("/* Auto-generated D3Q27 Kahan/Klein summation functions */\n")
        f.write("/* Generated by generate_d3q27_kahan.py */\n")
        f.write("/* CHANGE INIT - D3Q27 Kahan summation implementation */\n\n")
        f.write("#ifdef _D3Q27_\n")
        f.write("/* CHANGE INIT - D3Q27 constant definitions */\n")
        f.write("#define c0        0.0\n")
        f.write("#define c1        1.0\n")
        f.write("#define c2        2.0\n")
        f.write("/* CHANGE END */\n")
        f.write("#endif\n\n")
        f.write("/* Function declarations (add to collision.c header section): */\n")
        f.write("/* __device__ void d3q27_f2mode_chunk(double* mode, const double* __restrict__ fchunk); */\n")
        f.write("/* __device__ void d3q27_mode2f_chunk_kahan(double* mode, double* fchunk); */\n\n")
        f.write(f2mode_code)
        f.write("\n\n")
        f.write(mode2f_code)
        f.write("\n/* CHANGE END */\n")

    print(f"\nGenerated functions written to: {output_file}")
    print(f"Total lines: {f2mode_code.count(chr(10)) + mode2f_code.count(chr(10))}")

    # Also save matrices for reference
    np.savetxt("/home/bater/Sim/ludwig/d3q27_matrix_M.txt", M, fmt='%20.15f')
    np.savetxt("/home/bater/Sim/ludwig/d3q27_matrix_Minv.txt", M_inv, fmt='%20.15f')
    print("Matrices saved to d3q27_matrix_M.txt and d3q27_matrix_Minv.txt")

if __name__ == "__main__":
    main()
