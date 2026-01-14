#include <stdio.h>
#include <stdint.h>

/* D3Q27 model from lb_d3q27.h */
const int8_t cv[27][3] = {
    { 0, 0, 0},
    {-1,-1,-1}, {-1,-1, 0}, {-1,-1, 1}, {-1, 0,-1}, {-1, 0, 0}, {-1, 0, 1},
    {-1, 1,-1}, {-1, 1, 0}, {-1, 1, 1}, { 0,-1,-1}, { 0,-1, 0}, { 0,-1, 1},
    { 0, 0,-1},                                                 { 0, 0, 1},
    { 0, 1,-1}, { 0, 1, 0}, { 0, 1, 1}, { 1,-1,-1}, { 1,-1, 0}, { 1,-1, 1},
    { 1, 0,-1}, { 1, 0, 0}, { 1, 0, 1}, { 1, 1,-1}, { 1, 1, 0}, { 1, 1, 1}
};

const double wv[27] = {
    64.0/216.0,
     1.0/216.0,  4.0/216.0,  1.0/216.0,  4.0/216.0, 16.0/216.0,  4.0/216.0,
     1.0/216.0,  4.0/216.0,  1.0/216.0,  4.0/216.0, 16.0/216.0,  4.0/216.0,
    16.0/216.0,                                                  16.0/216.0,
     4.0/216.0, 16.0/216.0,  4.0/216.0,  1.0/216.0,  4.0/216.0,  1.0/216.0,
     4.0/216.0, 16.0/216.0,  4.0/216.0,  1.0/216.0,  4.0/216.0,  1.0/216.0
};

int main() {
    /* Test 1: Simple sum */
    double sum_x = 0.0, sum_y = 0.0, sum_z = 0.0;
    for (int p = 0; p < 27; p++) {
        sum_x += cv[p][0] * wv[p];
        sum_y += cv[p][1] * wv[p];
        sum_z += cv[p][2] * wv[p];
    }

    printf("Simple sum:\n");
    printf("  X: %.20e\n", sum_x);
    printf("  Y: %.20e\n", sum_y);
    printf("  Z: %.20e\n", sum_z);

    /* Test 2: Kahan sum */
    volatile double kahan_x = 0.0, kahan_y = 0.0, kahan_z = 0.0;
    volatile double c_x = 0.0, c_y = 0.0, c_z = 0.0;

    for (int p = 0; p < 27; p++) {
        volatile double val_x = cv[p][0] * wv[p];
        volatile double y_x = val_x - c_x;
        volatile double t_x = kahan_x + y_x;
        c_x = (t_x - kahan_x) - y_x;
        kahan_x = t_x;

        volatile double val_y = cv[p][1] * wv[p];
        volatile double y_y = val_y - c_y;
        volatile double t_y = kahan_y + y_y;
        c_y = (t_y - kahan_y) - y_y;
        kahan_y = t_y;

        volatile double val_z = cv[p][2] * wv[p];
        volatile double y_z = val_z - c_z;
        volatile double t_z = kahan_z + y_z;
        c_z = (t_z - kahan_z) - y_z;
        kahan_z = t_z;
    }

    printf("\nKahan sum:\n");
    printf("  X: %.20e\n", (double)kahan_x);
    printf("  Y: %.20e\n", (double)kahan_y);
    printf("  Z: %.20e\n", (double)kahan_z);

    return 0;
}
