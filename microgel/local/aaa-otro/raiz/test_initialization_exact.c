/*****************************************************************************
 *
 * test_initialization_exact.c
 *
 * Exactly replicate what happens in Ludwig initialization and momentum
 * calculation to find where the D3Q27 error comes from.
 *
 *****************************************************************************/

#include <stdio.h>
#include <math.h>

/* D3Q27 weights and velocities */
const double wv27[27] = {   64.0/216.0,
      1.0/216.0,  4.0/216.0, 1.0/216.0, 4.0/216.0, 16.0/216.0,  4.0/216.0,
      1.0/216.0,  4.0/216.0, 1.0/216.0, 4.0/216.0, 16.0/216.0,  4.0/216.0,
     16.0/216.0,                                               16.0/216.0,
      4.0/216.0, 16.0/216.0, 4.0/216.0, 1.0/216.0,  4.0/216.0,  1.0/216.0,
      4.0/216.0, 16.0/216.0, 4.0/216.0, 1.0/216.0,  4.0/216.0,  1.0/216.0};

const int cv27[27][3] = {
    { 0, 0, 0},
   {-1,-1,-1}, {-1,-1, 0}, {-1,-1, 1}, {-1, 0,-1}, {-1, 0, 0}, {-1, 0, 1},
   {-1, 1,-1}, {-1, 1, 0}, {-1, 1, 1}, { 0,-1,-1}, { 0,-1, 0}, { 0,-1, 1},
   { 0, 0,-1},                                                 { 0, 0, 1},
   { 0, 1,-1}, { 0, 1, 0}, { 0, 1, 1}, { 1,-1,-1}, { 1,-1, 0}, { 1,-1, 1},
   { 1, 0,-1}, { 1, 0, 0}, { 1, 0, 1}, { 1, 1,-1}, { 1, 1, 0}, { 1, 1, 1}
};

/* D3Q19 weights and velocities */
const double wv19[19] = {   12.0/36.0,
    1.0/36.0, 1.0/36.0, 2.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 2.0/36.0,
    1.0/36.0, 2.0/36.0, 2.0/36.0, 1.0/36.0, 2.0/36.0, 1.0/36.0, 1.0/36.0,
    1.0/36.0, 2.0/36.0, 1.0/36.0, 1.0/36.0};

const int cv19[19][3] = {
    { 0,  0,  0},
    { 1,  1,  0}, { 1,  0,  1}, { 1,  0,  0},
    { 1,  0, -1}, { 1, -1,  0}, { 0,  1,  1},
    { 0,  1,  0}, { 0,  1, -1}, { 0,  0,  1},
    { 0,  0, -1}, { 0, -1,  1}, { 0, -1,  0},
    { 0, -1, -1}, {-1,  1,  0}, {-1,  0,  1},
    {-1,  0,  0}, {-1,  0, -1}, {-1, -1,  0}
};

/* Kahan summation */
void kahan_add(volatile double *sum, volatile double *c, double val) {
    volatile double y = val - *c;
    volatile double t = *sum + y;
    *c = (t - *sum) - y;
    *sum = t;
}

/* Test exact replication of Ludwig's initialization and momentum calc */
void test_exact_replication(const char *name, int nvel,
                             const double *wv, const int cv[][3], double rho) {

    printf("\n=================================================================\n");
    printf("%s: Exact Replication Test (rho=%.1f)\n", name, rho);
    printf("=================================================================\n\n");

    /* Step 1: Initialize equilibrium distribution (lb_1st_moment_equilib_set) */
    double f[27];  /* Max size for D3Q27 */

    printf("Step 1: Initialize f_i = rho * w_i\n");
    printf("-----------------------------------------------------------------\n");
    for (int p = 0; p < nvel; p++) {
        f[p] = rho * wv[p];
    }
    printf("  Done (equilibrium at rest)\n\n");

    /* Step 2: Calculate momentum WITHOUT Kahan (naive summation) */
    printf("Step 2: Calculate momentum (naive summation)\n");
    printf("-----------------------------------------------------------------\n");
    double g_naive_x = 0.0, g_naive_y = 0.0, g_naive_z = 0.0;
    for (int p = 0; p < nvel; p++) {
        g_naive_x += f[p] * cv[p][0];
        g_naive_y += f[p] * cv[p][1];
        g_naive_z += f[p] * cv[p][2];
    }
    printf("  gx (naive): %.20e\n", g_naive_x);
    printf("  gy (naive): %.20e\n", g_naive_y);
    printf("  gz (naive): %.20e\n", g_naive_z);
    printf("\n");

    /* Step 3: Calculate momentum WITH Kahan (lb_1st_moment) */
    printf("Step 3: Calculate momentum (Kahan summation)\n");
    printf("-----------------------------------------------------------------\n");
    volatile double g_kahan_x = 0.0, g_kahan_y = 0.0, g_kahan_z = 0.0;
    volatile double c_x = 0.0, c_y = 0.0, c_z = 0.0;

    for (int p = 0; p < nvel; p++) {
        kahan_add(&g_kahan_x, &c_x, f[p] * cv[p][0]);
    }
    for (int p = 0; p < nvel; p++) {
        kahan_add(&g_kahan_y, &c_y, f[p] * cv[p][1]);
    }
    for (int p = 0; p < nvel; p++) {
        kahan_add(&g_kahan_z, &c_z, f[p] * cv[p][2]);
    }

    printf("  gx (Kahan): %.20e\n", g_kahan_x);
    printf("  gy (Kahan): %.20e\n", g_kahan_y);
    printf("  gz (Kahan): %.20e\n", g_kahan_z);
    printf("\n");

    /* Step 4: Difference */
    printf("Step 4: Difference (Kahan - naive)\n");
    printf("-----------------------------------------------------------------\n");
    printf("  Δgx: %.20e\n", g_kahan_x - g_naive_x);
    printf("  Δgy: %.20e\n", g_kahan_y - g_naive_y);
    printf("  Δgz: %.20e\n", g_kahan_z - g_naive_z);
    printf("\n");

    /* Step 5: Check if Kahan achieves exactly zero */
    printf("Step 5: Is momentum exactly zero with Kahan?\n");
    printf("-----------------------------------------------------------------\n");
    if (g_kahan_x == 0.0 && g_kahan_y == 0.0 && g_kahan_z == 0.0) {
        printf("  ✓ YES - Kahan summation achieves exactly zero!\n");
    } else {
        printf("  ✗ NO - Kahan summation does NOT achieve exactly zero\n");
        printf("  Maximum error: %.2e\n",
               fmax(fabs(g_kahan_x), fmax(fabs(g_kahan_y), fabs(g_kahan_z))));
    }
}

int main() {
    printf("=================================================================\n");
    printf("Exact Replication of Ludwig Initialization\n");
    printf("=================================================================\n");
    printf("\nThis test exactly replicates:\n");
    printf("1. lb_1st_moment_equilib_set: f[p] = rho * wv[p]\n");
    printf("2. lb_1st_moment: g[n] = sum_p (cv[p][n] * f[p]) with Kahan\n");
    printf("\n");

    /* Test with rho = 1.0 */
    test_exact_replication("D3Q19", 19, wv19, cv19, 1.0);
    test_exact_replication("D3Q27", 27, wv27, cv27, 1.0);

    /* Test with rho = 0.8 (from actual simulation) */
    test_exact_replication("D3Q19", 19, wv19, cv19, 0.8);
    test_exact_replication("D3Q27", 27, wv27, cv27, 0.8);

    printf("\n=================================================================\n");
    printf("Analysis Complete\n");
    printf("=================================================================\n");
    printf("\nKey Question: Why does D3Q19 achieve zero with Kahan but D3Q27 doesn't?\n");
    printf("\n");

    return 0;
}
