/*****************************************************************************
 *
 * Test weight roundoff errors in D3Q19 vs D3Q27
 *
 * This tests whether the weights themselves introduce systematic roundoff
 * that could explain the ~4e-18 error in D3Q27 but not D3Q19.
 *
 *****************************************************************************/

#include <stdio.h>
#include <math.h>

/* D3Q19 weights */
#define LB_WEIGHTS_D3Q19(wv) const double wv[19] = {   12.0/36.0, \
    1.0/36.0, 1.0/36.0, 2.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 2.0/36.0, \
    1.0/36.0, 2.0/36.0, 2.0/36.0, 1.0/36.0, 2.0/36.0, 1.0/36.0, 1.0/36.0, \
    1.0/36.0, 2.0/36.0, 1.0/36.0, 1.0/36.0}

/* D3Q27 weights */
#define LB_WEIGHTS_D3Q27(wv) const double wv[27] = {   64.0/216.0, \
      1.0/216.0,  4.0/216.0, 1.0/216.0, 4.0/216.0, 16.0/216.0,  4.0/216.0, \
      1.0/216.0,  4.0/216.0, 1.0/216.0, 4.0/216.0, 16.0/216.0,  4.0/216.0, \
     16.0/216.0,                                               16.0/216.0, \
      4.0/216.0, 16.0/216.0, 4.0/216.0, 1.0/216.0,  4.0/216.0,  1.0/216.0, \
      4.0/216.0, 16.0/216.0, 4.0/216.0, 1.0/216.0,  4.0/216.0,  1.0/216.0}

/* D3Q27 velocity vectors */
const int cv_d3q27[27][3] = {
    { 0, 0, 0},
   {-1,-1,-1}, {-1,-1, 0}, {-1,-1, 1}, {-1, 0,-1}, {-1, 0, 0}, {-1, 0, 1},
   {-1, 1,-1}, {-1, 1, 0}, {-1, 1, 1}, { 0,-1,-1}, { 0,-1, 0}, { 0,-1, 1},
   { 0, 0,-1},                                                 { 0, 0, 1},
   { 0, 1,-1}, { 0, 1, 0}, { 0, 1, 1}, { 1,-1,-1}, { 1,-1, 0}, { 1,-1, 1},
   { 1, 0,-1}, { 1, 0, 0}, { 1, 0, 1}, { 1, 1,-1}, { 1, 1, 0}, { 1, 1, 1}
};

/* D3Q19 velocity vectors */
const int cv_d3q19[19][3] = {
    { 0,  0,  0},
    { 1,  1,  0}, { 1,  0,  1}, { 1,  0,  0},
    { 1,  0, -1}, { 1, -1,  0}, { 0,  1,  1},
    { 0,  1,  0}, { 0,  1, -1}, { 0,  0,  1},
    { 0,  0, -1}, { 0, -1,  1}, { 0, -1,  0},
    { 0, -1, -1}, {-1,  1,  0}, {-1,  0,  1},
    {-1,  0,  0}, {-1,  0, -1}, {-1, -1,  0}
};

int main() {
    LB_WEIGHTS_D3Q19(wv19);
    LB_WEIGHTS_D3Q27(wv27);

    printf("=================================================================\n");
    printf("Testing Weight Roundoff Errors: D3Q19 vs D3Q27\n");
    printf("=================================================================\n\n");

    /* Test 1: Sum of weights should be exactly 1.0 */
    printf("Test 1: Sum of weights\n");
    printf("-----------------------------------------------------------------\n");

    double sum19 = 0.0;
    for (int p = 0; p < 19; p++) {
        sum19 += wv19[p];
    }

    double sum27 = 0.0;
    for (int p = 0; p < 27; p++) {
        sum27 += wv27[p];
    }

    printf("D3Q19 sum: %.20e (error: %.2e)\n", sum19, sum19 - 1.0);
    printf("D3Q27 sum: %.20e (error: %.2e)\n", sum27, sum27 - 1.0);
    printf("\n");

    /* Test 2: Weighted velocity sum for equilibrium at rest */
    printf("Test 2: Weighted velocity sum (rho=1.0, u=0)\n");
    printf("-----------------------------------------------------------------\n");
    printf("For equilibrium at rest, sum_i w_i * c_i should be exactly zero\n\n");

    /* D3Q19 */
    double sum19_x = 0.0, sum19_y = 0.0, sum19_z = 0.0;
    for (int p = 0; p < 19; p++) {
        sum19_x += wv19[p] * cv_d3q19[p][0];
        sum19_y += wv19[p] * cv_d3q19[p][1];
        sum19_z += wv19[p] * cv_d3q19[p][2];
    }

    printf("D3Q19 weighted velocity sum:\n");
    printf("  X: %.20e\n", sum19_x);
    printf("  Y: %.20e\n", sum19_y);
    printf("  Z: %.20e\n", sum19_z);
    printf("\n");

    /* D3Q27 */
    double sum27_x = 0.0, sum27_y = 0.0, sum27_z = 0.0;
    for (int p = 0; p < 27; p++) {
        sum27_x += wv27[p] * cv_d3q27[p][0];
        sum27_y += wv27[p] * cv_d3q27[p][1];
        sum27_z += wv27[p] * cv_d3q27[p][2];
    }

    printf("D3Q27 weighted velocity sum:\n");
    printf("  X: %.20e\n", sum27_x);
    printf("  Y: %.20e\n", sum27_y);
    printf("  Z: %.20e\n", sum27_z);
    printf("\n");

    /* Test 3: Check if 216.0 division introduces systematic error */
    printf("Test 3: Weight representation accuracy\n");
    printf("-----------------------------------------------------------------\n");
    printf("Checking if weights have systematic roundoff from division\n\n");

    /* Sample D3Q27 weights */
    double w1 = 1.0/216.0;
    double w4 = 4.0/216.0;
    double w16 = 16.0/216.0;
    double w64 = 64.0/216.0;

    /* Alternative: multiply then divide */
    double w1_alt = (1.0 * 36.0) / (216.0 * 36.0);
    double w4_alt = (4.0 * 36.0) / (216.0 * 36.0);

    printf("D3Q27 weight examples:\n");
    printf("  1/216 = %.20e\n", w1);
    printf("  4/216 = %.20e\n", w4);
    printf(" 16/216 = %.20e\n", w16);
    printf(" 64/216 = %.20e\n", w64);
    printf("\n");

    /* Test if 216 = 6^3 causes issues */
    double test_216 = 216.0;
    double test_36 = 36.0;
    printf("Division denominators:\n");
    printf("  216.0 = %.20e (binary representation may have roundoff)\n", test_216);
    printf("   36.0 = %.20e\n", test_36);
    printf("\n");

    /* Test 4: Equilibrium distribution for rho=0.8 */
    printf("Test 4: Equilibrium momentum for rho=0.8 (from simulation)\n");
    printf("-----------------------------------------------------------------\n");

    double rho = 0.8;

    /* D3Q19 equilibrium at rest */
    double g19_x = 0.0, g19_y = 0.0, g19_z = 0.0;
    for (int p = 0; p < 19; p++) {
        double f_eq = rho * wv19[p];
        g19_x += f_eq * cv_d3q19[p][0];
        g19_y += f_eq * cv_d3q19[p][1];
        g19_z += f_eq * cv_d3q19[p][2];
    }

    printf("D3Q19 momentum (rho=%.1f):\n", rho);
    printf("  gx: %.20e\n", g19_x);
    printf("  gy: %.20e\n", g19_y);
    printf("  gz: %.20e\n", g19_z);
    printf("\n");

    /* D3Q27 equilibrium at rest */
    double g27_x = 0.0, g27_y = 0.0, g27_z = 0.0;
    for (int p = 0; p < 27; p++) {
        double f_eq = rho * wv27[p];
        g27_x += f_eq * cv_d3q27[p][0];
        g27_y += f_eq * cv_d3q27[p][1];
        g27_z += f_eq * cv_d3q27[p][2];
    }

    printf("D3Q27 momentum (rho=%.1f):\n", rho);
    printf("  gx: %.20e\n", g27_x);
    printf("  gy: %.20e\n", g27_y);
    printf("  gz: %.20e\n", g27_z);
    printf("\n");

    printf("Difference (D3Q27 - D3Q19):\n");
    printf("  Δgx: %.20e\n", g27_x - g19_x);
    printf("  Δgy: %.20e\n", g27_y - g19_y);
    printf("  Δgz: %.20e\n", g27_z - g19_z);
    printf("\n");

    /* Compare with observed simulation error */
    printf("Observed simulation error: ~4.3e-18\n");
    printf("Predicted from equilibrium: %.2e\n", fabs(g27_x));

    if (fabs(g27_x) > 1e-20) {
        printf("\n*** FOUND IT: Equilibrium distribution has non-zero momentum! ***\n");
    }

    printf("\n=================================================================\n");
    printf("Analysis Complete\n");
    printf("=================================================================\n");

    return 0;
}
