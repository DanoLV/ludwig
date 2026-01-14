/*****************************************************************************
 *
 *  test_stencil_standalone.c
 *
 *  Standalone test for D3Q27 stencil - no dependencies on full Ludwig
 *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

/* Include only the stencil implementation */
#include "src/stencil_d3q27.c"
#include "src/stencils.c"
#include "src/lb_d3q27.c"

int main(int argc, char **argv) {

  stencil_t * s = NULL;
  int ifail = 0;

  printf("=============================================================\n");
  printf(" D3Q27 Stencil Analytical Tests (Standalone)\n");
  printf("=============================================================\n\n");

  /* Create stencil */
  ifail = stencil_d3q27_create(&s);
  if (ifail != 0) {
    printf("ERROR: Failed to create stencil\n");
    return -1;
  }

  printf("Stencil created successfully\n");
  printf("  ndim    = %d\n", s->ndim);
  printf("  npoints = %d\n", s->npoints);
  printf("  wlaplacian[0] = %.15f\n", s->wlaplacian[0]);
  printf("  wgradients[0] = %.15f\n\n", s->wgradients[0]);

  /* Test 1: Laplacian of f(x,y,z) = x² + y² + z² */
  printf("------------------------------------------------------\n");
  printf("Test 1: Laplacian of f(x,y,z) = x² + y² + z²\n");
  printf("        Expected: ∇²f = 6.0\n");
  printf("------------------------------------------------------\n");

  double x0 = 5.0, y0 = 3.0, z0 = -2.0;
  double f_center = x0*x0 + y0*y0 + z0*z0;
  double laplacian = s->wlaplacian[0] * f_center;

  for (int p = 1; p < s->npoints; p++) {
    double x = x0 + s->cv[p][0];
    double y = y0 + s->cv[p][1];
    double z = z0 + s->cv[p][2];
    double f = x*x + y*y + z*z;
    laplacian += s->wlaplacian[p] * f;
  }

  double expected_lap = 6.0;
  double error_lap = fabs(laplacian - expected_lap);

  printf("  Computed: %.15f\n", laplacian);
  printf("  Expected: %.15f\n", expected_lap);
  printf("  Error:    %.15e\n", error_lap);

  if (error_lap < 1.0e-10) {
    printf("  PASS ✓\n\n");
  } else {
    printf("  FAIL ✗\n\n");
    ifail = -1;
  }

  /* Test 2: Gradient of f(x,y,z) = 2x + 3y + 4z */
  printf("------------------------------------------------------\n");
  printf("Test 2: Gradient of f(x,y,z) = 2x + 3y + 4z\n");
  printf("        Expected: ∇f = (2, 3, 4)\n");
  printf("------------------------------------------------------\n");

  x0 = 1.0; y0 = -2.0; z0 = 3.0;
  double grad[3] = {0.0, 0.0, 0.0};

  for (int p = 0; p < s->npoints; p++) {
    double x = x0 + s->cv[p][0];
    double y = y0 + s->cv[p][1];
    double z = z0 + s->cv[p][2];
    double f = 2.0*x + 3.0*y + 4.0*z;

    for (int ia = 0; ia < 3; ia++) {
      grad[ia] += s->wgradients[p] * s->cv[p][ia] * f;
    }
  }

  double expected_grad[3] = {2.0, 3.0, 4.0};
  double error_grad[3];
  double max_error_grad = 0.0;

  printf("  Computed: (%.15f, %.15f, %.15f)\n", grad[0], grad[1], grad[2]);
  printf("  Expected: (%.15f, %.15f, %.15f)\n",
         expected_grad[0], expected_grad[1], expected_grad[2]);

  for (int ia = 0; ia < 3; ia++) {
    error_grad[ia] = fabs(grad[ia] - expected_grad[ia]);
    if (error_grad[ia] > max_error_grad) max_error_grad = error_grad[ia];
  }

  printf("  Error:    (%.15e, %.15e, %.15e)\n",
         error_grad[0], error_grad[1], error_grad[2]);

  if (max_error_grad < 1.0e-10) {
    printf("  PASS ✓\n\n");
  } else {
    printf("  FAIL ✗\n\n");
    ifail = -1;
  }

  /* Test 3: Consistency check - ratio of weights */
  printf("------------------------------------------------------\n");
  printf("Test 3: Weight ratio consistency\n");
  printf("        Expected: |wlaplacian/wgradients| = 2.0\n");
  printf("------------------------------------------------------\n");

  double ratio_ref = -1.0;
  int first_p = -1;

  for (int p = 1; p < s->npoints; p++) {
    if (s->wgradients[p] != 0.0) {
      double ratio = fabs(s->wlaplacian[p] / s->wgradients[p]);
      if (first_p < 0) {
        ratio_ref = ratio;
        first_p = p;
      } else {
        if (fabs(ratio - ratio_ref) > 1.0e-10) {
          printf("  ERROR: Inconsistent ratio at p=%d\n", p);
          printf("         ratio[%d] = %.15f\n", p, ratio);
          printf("         ratio_ref = %.15f\n", ratio_ref);
          ifail = -1;
        }
      }
    }
  }

  printf("  Reference ratio: %.15f\n", ratio_ref);
  printf("  Expected ratio:  %.15f\n", 2.0);
  printf("  Error:           %.15e\n", fabs(ratio_ref - 2.0));

  if (fabs(ratio_ref - 2.0) < 1.0e-10) {
    printf("  PASS ✓\n\n");
  } else {
    printf("  FAIL ✗\n\n");
    ifail = -1;
  }

  /* Clean up */
  stencil_free(&s);

  printf("=============================================================\n");
  if (ifail == 0) {
    printf(" ALL TESTS PASSED ✓\n");
  } else {
    printf(" SOME TESTS FAILED ✗\n");
  }
  printf("=============================================================\n");

  return ifail;
}
