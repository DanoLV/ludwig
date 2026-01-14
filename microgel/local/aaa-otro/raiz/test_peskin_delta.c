/*****************************************************************************
 *
 *  test_peskin_delta.c
 *
 *  Programa para calcular los deltas de Peskin para un punto dado
 *  y verificar que la suma de todos los deltas es 1.
 *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*****************************************************************************
 *
 *  d_peskin
 *
 *  Función delta de Peskin de 4 puntos.
 *  Basada en Peskin (2002).
 *
 *****************************************************************************/

double d_peskin(double r) {

    double rmod;
    double delta = 0.0;

    rmod = fabs(r);

    if (rmod <= 1.0) {
        delta = 0.125 * (3.0 - 2.0 * rmod + sqrt(1.0 + 4.0 * rmod - 4.0 * rmod * rmod));
    }
    else if (rmod <= 2.0) {
        delta = 0.125 * (5.0 - 2.0 * rmod - sqrt(-7.0 + 12.0 * rmod - 4.0 * rmod * rmod));
    }

    return delta;
}

/*****************************************************************************
 *
 *  main
 *
 *  Calcula los deltas de Peskin para un punto dado y verifica que suman 1.
 *
 *****************************************************************************/

int main(int argc, char** argv) {

    double r0[3];           /* Posición del punto (en coordenadas continuas) */
    int i_min, i_max;       /* Rango de nodos en x */
    int j_min, j_max;       /* Rango de nodos en y */
    int k_min, k_max;       /* Rango de nodos en z */
    int i, j, k;
    double r_vec[3];        /* Vector de distancia */
    double delta_x, delta_y, delta_z;
    double delta_3d;
    double sum_total = 0.0;
    double drange = 2.0;    /* Rango del delta de Peskin (4 puntos) */

    /* Leer la posición del punto desde argumentos de línea de comandos */
    if (argc != 4) {
        printf("Uso: %s <x> <y> <z>\n", argv[0]);
        printf("Ejemplo: %s 5.3 5.7 5.2\n", argv[0]);
        return 1;
    }

    r0[0] = atof(argv[1]);
    r0[1] = atof(argv[2]);
    r0[2] = atof(argv[3]);

    printf("\n");
    printf("=======================================================\n");
    printf("  Cálculo de deltas de Peskin para un punto\n");
    printf("=======================================================\n");
    printf("\n");
    printf("Posición del punto: (%.6f, %.6f, %.6f)\n", r0[0], r0[1], r0[2]);
    printf("\n");

    /* Determinar el rango de nodos afectados */
    i_min = (int)floor(r0[0] - drange);
    i_max = (int)floor(r0[0] + drange);
    j_min = (int)floor(r0[1] - drange);
    j_max = (int)floor(r0[1] + drange);
    k_min = (int)floor(r0[2] - drange);
    k_max = (int)floor(r0[2] + drange);

    printf("Rango de nodos afectados:\n");
    printf("  x: [%d, %d]\n", i_min, i_max);
    printf("  y: [%d, %d]\n", j_min, j_max);
    printf("  z: [%d, %d]\n", k_min, k_max);
    printf("\n");

    printf("Deltas para cada nodo:\n");
    printf("-------------------------------------------------------\n");
    printf("  i   j   k  |   dx   |   dy   |   dz   | delta_3D\n");
    printf("-------------------------------------------------------\n");

    /* Calcular deltas para cada nodo en el rango */
    for (i = i_min; i <= i_max; i++) {
        for (j = j_min; j <= j_max; j++) {
            for (k = k_min; k <= k_max; k++) {

                /* Vector de distancia desde el nodo al punto */
                r_vec[0] = r0[0] - (i + 0.5);  /* +0.5 porque el nodo está en el centro */
                r_vec[1] = r0[1] - (j + 0.5);
                r_vec[2] = r0[2] - (k + 0.5);

                /* Calcular deltas en cada dirección */
                delta_x = d_peskin(r_vec[0]);
                delta_y = d_peskin(r_vec[1]);
                delta_z = d_peskin(r_vec[2]);

                /* Delta 3D es el producto de los deltas 1D */
                delta_3d = delta_x * delta_y * delta_z;

                /* Acumular suma total */
                sum_total += delta_3d;

                /* Imprimir solo si el delta es significativo */
                if (delta_3d > 1.0e-15) {
                    printf("%3d %3d %3d | %7.5f | %7.5f | %7.5f | %10.8f\n",
                           i, j, k, delta_x, delta_y, delta_z, delta_3d);
                }
            }
        }
    }

    printf("-------------------------------------------------------\n");
    printf("\n");
    printf("Suma total de deltas: %.15f\n", sum_total);
    printf("Error (debe ser ~0):  %.15e\n", fabs(sum_total - 1.0));
    printf("\n");

    if (fabs(sum_total - 1.0) < 1.0e-10) {
        printf("✓ VERIFICACIÓN CORRECTA: La suma es 1.0\n");
    } else {
        printf("✗ ADVERTENCIA: La suma NO es 1.0 (diferencia = %.15e)\n",
               fabs(sum_total - 1.0));
    }
    printf("\n");

    return 0;
}
