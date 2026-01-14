/*****************************************************************************
 *
 *  test_stencil_d3q27_analytical.c
 *
 *  Comprehensive tests for D3Q27 stencil weights using analytical fields.
 *  Verifies that the stencil correctly computes Laplacian and gradients
 *  for known analytical functions.
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  (c) 2025 The University of Edinburgh
 *  Contributing authors:
 *  Daniel La Valle (University of Barcelona)
 *****************************************************************************/

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>

#include "pe.h"
#include "stencil_d3q27.h"

int test_stencil_d3q27_analytical_create(void);
int test_stencil_d3q27_analytical_laplacian_quadratic(void);
int test_stencil_d3q27_analytical_laplacian_polynomial(void);
int test_stencil_d3q27_analytical_gradient_linear(void);
int test_stencil_d3q27_analytical_gradient_quadratic(void);
int test_stencil_d3q27_analytical_consistency_check(void);

/*****************************************************************************
 *
 *  test_stencil_d3q27_analytical_suite
 *
 *****************************************************************************/

int test_stencil_d3q27_analytical_suite(void) {

  pe_t * pe = NULL;

  pe_create(MPI_COMM_WORLD, PE_QUIET, &pe);

  test_stencil_d3q27_analytical_create();
  test_stencil_d3q27_analytical_laplacian_quadratic();
  test_stencil_d3q27_analytical_laplacian_polynomial();
  test_stencil_d3q27_analytical_gradient_linear();
  test_stencil_d3q27_analytical_gradient_quadratic();
  test_stencil_d3q27_analytical_consistency_check();

  pe_info(pe, "%-9s %s\n", "PASS", __FILE__);
  pe_free(pe);

  return 0;
}

/*****************************************************************************
 *
 *  test_stencil_d3q27_analytical_create
 *
 *  Basic creation test with updated weight values
 *
 *****************************************************************************/

int test_stencil_d3q27_analytical_create(void) {

  int ifail = 0;
  stencil_t * s = NULL;

  ifail = stencil_d3q27_create(&s);
  assert(ifail == 0);
  assert(s);
  assert(s->ndim == 3);
  assert(s->npoints == 27);
  assert(s->cv);
  assert(s->wlaplacian);
  assert(s->wgradients);

  /* Updated expected values after -6.0 fix (changed from -216.0) */
  /* wlaplacian[0] should be approximately 152.0/36.0 = 4.222... */
  double expected_wlap0 = 152.0 / 36.0;
  double tolerance = 1.0e-10;

  if (fabs(s->wlaplacian[0] - expected_wlap0) > tolerance) {
    printf("ERROR: wlaplacian[0] = %.15f, expected = %.15f\n",
           s->wlaplacian[0], expected_wlap0);
    ifail = -1;
  }
  if (s->wgradients[0] != 0.0) ifail = -1;
  assert(ifail == 0);

  /* Print some diagnostic info */
  printf("Stencil D3Q27 weights:\n");
  printf("  wlaplacian[0] = %.15f (expected %.15f)\n",
         s->wlaplacian[0], expected_wlap0);
  printf("  wgradients[0] = %.15f\n", s->wgradients[0]);

  ifail = stencil_free(&s);
  assert(ifail == 0);
  assert(s == NULL);

  return ifail;
}

/*****************************************************************************
 *
 *  test_stencil_d3q27_analytical_laplacian_quadratic
 *
 *  Test Laplacian computation for f(x,y,z) = x² + y² + z²
 *  Analytical: ∇²f = 2 + 2 + 2 = 6
 *  Stencil convention: computes -∇²f = -6 (for Poisson eq: -∇²ψ = ρ)
 *
 *****************************************************************************/

int test_stencil_d3q27_analytical_laplacian_quadratic(void) {

  int ifail = 0;
  stencil_t * s = NULL;

  ifail = stencil_d3q27_create(&s);
  assert(ifail == 0);

  /* Test point at origin (but could be any point) */
  double x0 = 5.0, y0 = 3.0, z0 = -2.0;

  /* Field: f(x,y,z) = x² + y² + z² */
  /* Analytical Laplacian: ∇²f = 6 */
  double f_center = x0*x0 + y0*y0 + z0*z0;

  /* Compute stencil approximation to Laplacian */
  double laplacian = s->wlaplacian[0] * f_center;

  for (int p = 1; p < s->npoints; p++) {
    double x = x0 + s->cv[p][0];
    double y = y0 + s->cv[p][1];
    double z = z0 + s->cv[p][2];
    double f = x*x + y*y + z*z;
    laplacian += s->wlaplacian[p] * f;
  }

  /* Note: The stencil weights compute -∇²f by convention (for Poisson eq: -∇²ψ = ρ) */
  double expected = -6.0;
  double tolerance = 1.0e-10;

  printf("\nLaplacian test (quadratic field):\n");
  printf("  Field: f(x,y,z) = x² + y² + z²\n");
  printf("  Computed -∇²f = %.15f\n", laplacian);
  printf("  Expected -∇²f = %.15f (∇²f = 6.0)\n", expected);
  printf("  Error = %.15e\n", fabs(laplacian - expected));

  if (fabs(laplacian - expected) > tolerance) {
    printf("ERROR: Laplacian mismatch!\n");
    ifail = -1;
  }
  assert(ifail == 0);

  stencil_free(&s);
  return ifail;
}

/*****************************************************************************
 *
 *  test_stencil_d3q27_analytical_laplacian_polynomial
 *
 *  Test Laplacian for f(x,y,z) = 3x² + 2y² + 4z²
 *  Analytical: ∇²f = 6 + 4 + 8 = 18
 *  Stencil convention: computes -∇²f = -18
 *
 *****************************************************************************/

int test_stencil_d3q27_analytical_laplacian_polynomial(void) {

  int ifail = 0;
  stencil_t * s = NULL;

  ifail = stencil_d3q27_create(&s);
  assert(ifail == 0);

  double x0 = -1.0, y0 = 2.5, z0 = 0.5;

  /* Field: f(x,y,z) = 3x² + 2y² + 4z² */
  /* Analytical Laplacian: ∇²f = 6 + 4 + 8 = 18 */
  /* Stencil computes: -∇²f = -18 */
  double f_center = 3.0*x0*x0 + 2.0*y0*y0 + 4.0*z0*z0;

  double laplacian = s->wlaplacian[0] * f_center;

  for (int p = 1; p < s->npoints; p++) {
    double x = x0 + s->cv[p][0];
    double y = y0 + s->cv[p][1];
    double z = z0 + s->cv[p][2];
    double f = 3.0*x*x + 2.0*y*y + 4.0*z*z;
    laplacian += s->wlaplacian[p] * f;
  }

  double expected = -18.0;
  double tolerance = 1.0e-10;

  printf("\nLaplacian test (polynomial field):\n");
  printf("  Field: f(x,y,z) = 3x² + 2y² + 4z²\n");
  printf("  Computed -∇²f = %.15f\n", laplacian);
  printf("  Expected -∇²f = %.15f (∇²f = 18.0)\n", expected);
  printf("  Error = %.15e\n", fabs(laplacian - expected));

  if (fabs(laplacian - expected) > tolerance) {
    printf("ERROR: Laplacian mismatch!\n");
    ifail = -1;
  }
  assert(ifail == 0);

  stencil_free(&s);
  return ifail;
}

/*****************************************************************************
 *
 *  test_stencil_d3q27_analytical_gradient_linear
 *
 *  Test gradient computation for f(x,y,z) = 2x + 3y + 4z
 *  Analytical: ∇f = (2, 3, 4)
 *
 *****************************************************************************/

int test_stencil_d3q27_analytical_gradient_linear(void) {

  int ifail = 0;
  stencil_t * s = NULL;

  ifail = stencil_d3q27_create(&s);
  assert(ifail == 0);

  double x0 = 1.0, y0 = -2.0, z0 = 3.0;

  /* Field: f(x,y,z) = 2x + 3y + 4z */
  /* Analytical gradient: ∇f = (2, 3, 4) */

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

  double expected[3] = {2.0, 3.0, 4.0};
  double tolerance = 1.0e-10;

  printf("\nGradient test (linear field):\n");
  printf("  Field: f(x,y,z) = 2x + 3y + 4z\n");
  printf("  Computed gradient = (%.15f, %.15f, %.15f)\n", grad[0], grad[1], grad[2]);
  printf("  Expected gradient = (%.15f, %.15f, %.15f)\n", expected[0], expected[1], expected[2]);

  for (int ia = 0; ia < 3; ia++) {
    printf("  Error[%d] = %.15e\n", ia, fabs(grad[ia] - expected[ia]));
    if (fabs(grad[ia] - expected[ia]) > tolerance) {
      printf("ERROR: Gradient component %d mismatch!\n", ia);
      ifail = -1;
    }
  }
  assert(ifail == 0);

  stencil_free(&s);
  return ifail;
}

/*****************************************************************************
 *
 *  test_stencil_d3q27_analytical_gradient_quadratic
 *
 *  Test gradient for f(x,y,z) = x² + y² + z²
 *  Analytical: ∇f = (2x, 2y, 2z)
 *
 *****************************************************************************/

int test_stencil_d3q27_analytical_gradient_quadratic(void) {

  int ifail = 0;
  stencil_t * s = NULL;

  ifail = stencil_d3q27_create(&s);
  assert(ifail == 0);

  double x0 = 2.0, y0 = -1.5, z0 = 3.5;

  /* Field: f(x,y,z) = x² + y² + z² */
  /* Analytical gradient at (x0,y0,z0): ∇f = (2x0, 2y0, 2z0) */

  double grad[3] = {0.0, 0.0, 0.0};

  for (int p = 0; p < s->npoints; p++) {
    double x = x0 + s->cv[p][0];
    double y = y0 + s->cv[p][1];
    double z = z0 + s->cv[p][2];
    double f = x*x + y*y + z*z;

    for (int ia = 0; ia < 3; ia++) {
      grad[ia] += s->wgradients[p] * s->cv[p][ia] * f;
    }
  }

  double expected[3] = {2.0*x0, 2.0*y0, 2.0*z0};
  double tolerance = 1.0e-10;

  printf("\nGradient test (quadratic field):\n");
  printf("  Field: f(x,y,z) = x² + y² + z²\n");
  printf("  Test point: (%.15f, %.15f, %.15f)\n", x0, y0, z0);
  printf("  Computed gradient = (%.15f, %.15f, %.15f)\n", grad[0], grad[1], grad[2]);
  printf("  Expected gradient = (%.15f, %.15f, %.15f)\n", expected[0], expected[1], expected[2]);

  for (int ia = 0; ia < 3; ia++) {
    printf("  Error[%d] = %.15e\n", ia, fabs(grad[ia] - expected[ia]));
    if (fabs(grad[ia] - expected[ia]) > tolerance) {
      printf("ERROR: Gradient component %d mismatch!\n", ia);
      ifail = -1;
    }
  }
  assert(ifail == 0);

  stencil_free(&s);
  return ifail;
}

/*****************************************************************************
 *
 *  test_stencil_d3q27_analytical_consistency_check
 *
 *  Verify consistency relation: |wlaplacian[p]/wgradients[p]| should be
 *  constant for all p > 0 (similar to D3Q7 stencil where ratio = 2.0)
 *
 *****************************************************************************/

int test_stencil_d3q27_analytical_consistency_check(void) {

  int ifail = 0;
  stencil_t * s = NULL;

  ifail = stencil_d3q27_create(&s);
  assert(ifail == 0);

  printf("\nConsistency check (weight ratios):\n");
  printf("  Point | cv[x,y,z] | wlaplacian | wgradients | ratio\n");
  printf("  ------|-----------|------------|------------|-------\n");

  double ratio_ref = 0.0;
  int first_nonzero = -1;

  for (int p = 0; p < s->npoints; p++) {
    printf("  %4d  | [%2d,%2d,%2d] | %10.6f | %10.6f | ",
           p, s->cv[p][0], s->cv[p][1], s->cv[p][2],
           s->wlaplacian[p], s->wgradients[p]);

    if (p > 0 && s->wgradients[p] != 0.0) {
      double ratio = fabs(s->wlaplacian[p] / s->wgradients[p]);
      printf("%.6f", ratio);

      if (first_nonzero < 0) {
        ratio_ref = ratio;
        first_nonzero = p;
      } else {
        /* Check if ratio is consistent */
        double tolerance = 1.0e-10;
        if (fabs(ratio - ratio_ref) > tolerance) {
          printf(" <- INCONSISTENT!");
          ifail = -1;
        }
      }
    } else {
      printf("  N/A  ");
    }
    printf("\n");
  }

  printf("\nReference ratio (should be 2.0): %.15f\n", ratio_ref);

  double expected_ratio = 2.0;
  double tolerance = 1.0e-10;
  if (fabs(ratio_ref - expected_ratio) > tolerance) {
    printf("ERROR: Weight ratio = %.15f, expected = %.15f\n",
           ratio_ref, expected_ratio);
    ifail = -1;
  }

  assert(ifail == 0);

  stencil_free(&s);
  return ifail;
}
