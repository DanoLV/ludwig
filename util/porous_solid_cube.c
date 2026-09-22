/*****************************************************************************
 *
 *  porous_solid_cube.c
 *
 *  Produces a file containing a porous solid network suitable for input
 *  into the main code, using the same colloid configuration file format
 *  as microgel_poly_cross_density.c.
 *
 *  The cube side length, the mean bond distance, its tolerance ("delta")
 *  and the point density are specified on the command line. The cubic
 *  simulation volume is then filled with points at random subject to the
 *  constraint that no two points may lie closer than (lbond - delta).
 *  Bonds are then created between any pair of points whose separation
 *  falls within [lbond - delta, lbond + delta], up to a maximum number
 *  of bonds per point.
 *
 *  Periodicity: each of the x, y, z directions can independently be
 *  declared periodic (the default) or not. When periodic, both the
 *  minimum-separation rejection test and the bonding test use the
 *  minimum-image convention, so the network percolates seamlessly
 *  across that pair of faces -- matching how Ludwig itself resolves
 *  bonded colloids that straddle a periodic (or MPI rank) boundary,
 *  via cs_minimum_distance() and the colloid halo exchange.
 *
 *  Fixed faces: one or more cube faces (e.g. "z-", "x+") can be marked
 *  fixed via --fixed-face. Every point within --fixed-thickness of a
 *  marked face has isfixedr/isfixedv/isfixedw (and the per-axis
 *  variants) set, so it does not move during the simulation -- typically
 *  used to anchor a layer of points against a wall. A face that is
 *  marked fixed is, by default, also treated as non-periodic (this can
 *  be overridden explicitly with --periodic-x/-y/-z).
 *
 *  For compilation instructions see the Makefile.
 *
 *  $ make porous_solid_cube
 *
 *  $ ./porous_solid_cube
 *
 *  should produce a file config.cds.init.001-001 in the specified format.
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  Contributing authors
 *  Daniel La Valle (daniel.lavalle@ub.edu)
 *  (c) 2026- University of Barcelona
 *
 *****************************************************************************/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../src/colloid.h"
#include "../src/pe.h"
#include "../src/coords.h"
#include "../src/util.h"
#include "../src/util_fopen.h"

#include <time.h>
#include <float.h>
#include <string.h>
#include <getopt.h>
#include <errno.h>
#include <ctype.h>

#define M_PI 3.14159265358979323846
#define MAX_FIXED_FACES 6

/* A point in the network. Bonds are capped at NBOND_MAX (colloid.h). */

typedef struct {
  int id;                    /* 1-based id, matches colloid_state_t.index */
  double r[3];
  int neighbors[NBOND_MAX];
  int num_neighbors;
  int is_fixed;               /* 1 if anchored (no position/velocity update) */
} Point3D;

/* Uniform grid (cell list) over the cubic domain, used to avoid an
 * O(n^2) search for both the rejection-sampling placement step and
 * the bonding step. Cell size is (lbond + delta), so every pair of
 * points within bonding range are guaranteed to be found by scanning
 * only the 3x3x3 block of cells centred on a given cell (wrapped
 * around on any axis declared periodic). */

typedef struct {
  int* idx;
  int count;
  int capacity;
} Cell;

typedef struct {
  Cell* cells;
  int nx, ny, nz;
  double cellw[3];
  double lmin[3];
} Grid;

typedef struct {
  Point3D* points;
  int npoints;          /* number currently placed */
  int capacity;          /* allocated capacity */
  double side;            /* cube side length */
  double lmin[3];          /* domain origin (cs_lmin) */
  double lbond;            /* target mean bond length */
  double delta;            /* tolerance on bond length */
  double min_distance;      /* lbond - delta, clamped to >= 0 */
  int max_links;          /* maximum bonds per point */
  int periodic[3];          /* per-axis periodicity used for bonding */
} PorousSolid;

static void grid_create(Grid* grid, const double lmin[3], double side, double cellsize);
static void grid_free(Grid* grid);
static void grid_cell_of(const Grid* grid, const double r[3], int* ix, int* iy, int* iz);
static void grid_add(Grid* grid, int ix, int iy, int iz, int point_id);
static int neighborCellIndices(int base, int n, int periodic, int out[3]);

static double pairDistance(const PorousSolid* solid, const double ra[3], const double rb[3]);
static double randomDouble(double min, double max);
static int isValidPosition(const PorousSolid* solid, const Grid* grid, const double r[3]);
static void generateRandomPoints(PorousSolid* solid, Grid* grid, int target_npoints);
static void createBonds(PorousSolid* solid, Grid* grid);
static int alreadyBonded(const Point3D* p, int other_id);
static void markFixedFaces(PorousSolid* solid, const int fixed_axis[], const int fixed_side[],
                            int n_fixed_faces, double thickness);
static void printNetworkStats(const PorousSolid* solid);

enum format { ASCII, BINARY };

static void colloid_init_write_file(const int nc, const colloid_state_t* pc, const char* filename, const int form);
static void colloid_init_state(colloid_state_t* state, const PorousSolid* solid,
                                double a0, double ah, double q0, double epsilon,
                                double sa, double saf, int bc, int shape,
                                double drmax, double al);

static int parseFixedFaceArg(const char* arg, int axis_out[], int side_out[], int max_faces);

static int clopt(int argc, char** argv,
                  int* side, double* lbond, double* delta, double* density,
                  int* max_links, unsigned int* seed,
                  double* irad, double* hrad, double* charge, double* permittivity,
                  double* sa, double* saf, double* drmax, double* offset,
                  int* periodic_x, int* periodic_y, int* periodic_z,
                  const char** fixed_face_arg, double* fixed_thickness);

/*****************************************************************************
 *
 *  main
 *
 *****************************************************************************/

int main(int argc, char** argv) {

  int side = 32;              /* Cube side length (also sets ntotal) */

  double lbond = 1.0;         /* Mean bond distance */
  double delta = 0.2;         /* Bond length tolerance */
  double density = 0.3;       /* Points per unit volume */
  int max_links = NBOND_MAX;  /* Max bonds per point */
  unsigned int seed = (unsigned int)time(NULL);

  double a0 = 0.1;
  double ah = 0.1;
  double al = 1.5;
  double drmax = 0.8;
  double q0 = 0.0;
  double epsilon = 0.0;
  double sa = 0.0;
  double saf = 0.0;

  /* -1 means "not set on the command line" and is resolved below */
  int periodic_x = -1, periodic_y = -1, periodic_z = -1;
  const char* fixed_face_arg = NULL;
  double fixed_thickness = -1.0; /* resolved to lbond if unset */

  int bc = COLLOID_BC_SUBGRID;
  int shape = COLLOID_SHAPE_SPHERE;
  int file_format = ASCII;

  colloid_state_t* state;
  pe_t* pe;
  cs_t* cs;
  int ntotal[3];

  MPI_Init(&argc, &argv);
  pe_create(MPI_COMM_WORLD, PE_QUIET, &pe);
  assert(pe_mpi_size(pe) == 1);

  clopt(argc, argv, &side, &lbond, &delta, &density, &max_links, &seed,
        &a0, &ah, &q0, &epsilon, &sa, &saf, &drmax, &al,
        &periodic_x, &periodic_y, &periodic_z, &fixed_face_arg, &fixed_thickness);

  if (max_links > NBOND_MAX) max_links = NBOND_MAX;
  if (max_links < 0) max_links = 0;
  if (delta < 0.0) delta = 0.0;
  if (delta > lbond) delta = lbond; /* keep min_distance >= 0 */
  if (fixed_thickness < 0.0) fixed_thickness = lbond;

  int fixed_axis[MAX_FIXED_FACES];
  int fixed_side[MAX_FIXED_FACES];
  int n_fixed_faces = 0;
  if (fixed_face_arg != NULL) {
    n_fixed_faces = parseFixedFaceArg(fixed_face_arg, fixed_axis, fixed_side, MAX_FIXED_FACES);
  }

  /* A face marked fixed is, unless the user overrode it explicitly,
   * treated as non-periodic: a wall face is generally not periodic. */
  for (int k = 0; k < n_fixed_faces; k++) {
    int axis = fixed_axis[k];
    if (axis == X && periodic_x == -1) periodic_x = 0;
    if (axis == Y && periodic_y == -1) periodic_y = 0;
    if (axis == Z && periodic_z == -1) periodic_z = 0;
  }
  if (periodic_x == -1) periodic_x = 1;
  if (periodic_y == -1) periodic_y = 1;
  if (periodic_z == -1) periodic_z = 1;

  ntotal[X] = side;
  ntotal[Y] = side;
  ntotal[Z] = side;

  int periodic[3];
  periodic[X] = periodic_x;
  periodic[Y] = periodic_y;
  periodic[Z] = periodic_z;

  cs_create(pe, &cs);
  cs_ntotal_set(cs, ntotal);
  cs_periodicity_set(cs, periodic);
  cs_init(cs);

  srand(seed);

  PorousSolid solid;
  memset(&solid, 0, sizeof(PorousSolid));
  cs_lmin(cs, solid.lmin);
  solid.side = (double)side;
  solid.lbond = lbond;
  solid.delta = delta;
  solid.min_distance = lbond - delta;
  solid.max_links = max_links;
  solid.periodic[X] = periodic_x;
  solid.periodic[Y] = periodic_y;
  solid.periodic[Z] = periodic_z;

  int target_npoints = (int)llround(density * pow(solid.side, 3.0));
  if (target_npoints < 0) target_npoints = 0;

  printf("Generating porous solid with parameters:\n");
  printf("  Random seed: %u\n", seed);
  printf("  Cube side: %d\n", side);
  printf("  Periodicity (x,y,z): %d,%d,%d\n", periodic_x, periodic_y, periodic_z);
  printf("  Mean bond distance: %.4f\n", lbond);
  printf("  Bond delta: %.4f (accepted range [%.4f, %.4f])\n",
         delta, lbond - delta, lbond + delta);
  printf("  Point density: %.4f\n", density);
  printf("  Target number of points: %d\n", target_npoints);
  printf("  Max bonds per point: %d\n", max_links);
  if (n_fixed_faces > 0) {
    printf("  Fixed faces: %s (thickness %.4f)\n", fixed_face_arg, fixed_thickness);
  }

  solid.capacity = target_npoints;
  solid.points = (Point3D*)calloc(solid.capacity > 0 ? solid.capacity : 1, sizeof(Point3D));
  assert(solid.points != NULL);

  Grid grid;
  grid_create(&grid, solid.lmin, solid.side, lbond + delta);

  generateRandomPoints(&solid, &grid, target_npoints);
  markFixedFaces(&solid, fixed_axis, fixed_side, n_fixed_faces, fixed_thickness);
  createBonds(&solid, &grid);
  printNetworkStats(&solid);

  grid_free(&grid);

  state = (colloid_state_t*)calloc(solid.npoints, sizeof(colloid_state_t));
  assert(state != NULL);

  colloid_init_state(state, &solid, a0, ah, q0, epsilon, sa, saf, bc, shape, drmax, al);

  colloid_init_write_file(solid.npoints, state, "config.cds.init.001-001", file_format);

  free(state);
  free(solid.points);

  cs_free(cs);
  pe_free(pe);
  MPI_Finalize();

  return 0;
}

/*****************************************************************************
 *
 *  grid_create / grid_free / grid_cell_of / grid_add / neighborCellIndices
 *
 *****************************************************************************/

static void grid_create(Grid* grid, const double lmin[3], double side, double cellsize) {

  if (cellsize < DBL_EPSILON) cellsize = DBL_EPSILON;

  grid->nx = (int)(side / cellsize);
  grid->ny = grid->nx;
  grid->nz = grid->nx;
  if (grid->nx < 1) grid->nx = 1;
  if (grid->ny < 1) grid->ny = 1;
  if (grid->nz < 1) grid->nz = 1;

  grid->cellw[X] = side / grid->nx;
  grid->cellw[Y] = side / grid->ny;
  grid->cellw[Z] = side / grid->nz;

  grid->lmin[X] = lmin[X];
  grid->lmin[Y] = lmin[Y];
  grid->lmin[Z] = lmin[Z];

  int ncells = grid->nx * grid->ny * grid->nz;
  grid->cells = (Cell*)calloc(ncells, sizeof(Cell));
  assert(grid->cells != NULL);
}

static void grid_free(Grid* grid) {
  int ncells = grid->nx * grid->ny * grid->nz;
  for (int i = 0; i < ncells; i++) {
    free(grid->cells[i].idx);
  }
  free(grid->cells);
  grid->cells = NULL;
}

static void grid_cell_of(const Grid* grid, const double r[3], int* ix, int* iy, int* iz) {

  *ix = (int)((r[X] - grid->lmin[X]) / grid->cellw[X]);
  *iy = (int)((r[Y] - grid->lmin[Y]) / grid->cellw[Y]);
  *iz = (int)((r[Z] - grid->lmin[Z]) / grid->cellw[Z]);

  if (*ix < 0) *ix = 0;
  if (*iy < 0) *iy = 0;
  if (*iz < 0) *iz = 0;
  if (*ix >= grid->nx) *ix = grid->nx - 1;
  if (*iy >= grid->ny) *iy = grid->ny - 1;
  if (*iz >= grid->nz) *iz = grid->nz - 1;
}

static void grid_add(Grid* grid, int ix, int iy, int iz, int point_id) {

  Cell* cell = &grid->cells[(ix * grid->ny + iy) * grid->nz + iz];

  if (cell->count >= cell->capacity) {
    int new_capacity = cell->capacity == 0 ? 4 : cell->capacity * 2;
    int* new_idx = (int*)realloc(cell->idx, new_capacity * sizeof(int));
    assert(new_idx != NULL);
    cell->idx = new_idx;
    cell->capacity = new_capacity;
  }

  cell->idx[cell->count++] = point_id;
}

/* Candidate cell indices (-1, 0, +1 offset from base) along one axis.
 * When periodic, offsets wrap around and duplicates (which occur when
 * the grid has fewer than 3 cells along this axis) are removed. When
 * not periodic, out-of-range offsets are simply omitted. Returns the
 * number of indices written to out[]. */

static int neighborCellIndices(int base, int n, int periodic, int out[3]) {

  int count = 0;

  for (int d = -1; d <= 1; d++) {
    int c = base + d;
    if (periodic) {
      c = ((c % n) + n) % n;
    }
    else {
      if (c < 0 || c >= n) continue;
    }

    int dup = 0;
    for (int k = 0; k < count; k++) {
      if (out[k] == c) { dup = 1; break; }
    }
    if (!dup) out[count++] = c;
  }

  return count;
}

/*****************************************************************************
 *
 *  pairDistance
 *
 *  Separation between two positions, using the minimum-image convention
 *  on any axis declared periodic -- matching cs_minimum_distance(), so
 *  a bond spanning a periodic face is measured the same way Ludwig
 *  itself will measure it when computing the bond force at run time.
 *
 *****************************************************************************/

static double pairDistance(const PorousSolid* solid, const double ra[3], const double rb[3]) {

  double d[3];

  for (int ia = 0; ia < 3; ia++) {
    d[ia] = ra[ia] - rb[ia];
    if (solid->periodic[ia]) {
      if (d[ia] >  0.5 * solid->side) d[ia] -= solid->side;
      if (d[ia] < -0.5 * solid->side) d[ia] += solid->side;
    }
  }

  return sqrt(d[X] * d[X] + d[Y] * d[Y] + d[Z] * d[Z]);
}

static double randomDouble(double min, double max) {
  return min + ((double)rand() / RAND_MAX) * (max - min);
}

/*****************************************************************************
 *
 *  isValidPosition
 *
 *  A trial position is valid if it is not closer than min_distance to
 *  any point already placed (minimum-image on periodic axes). Only the
 *  neighbouring grid cells need to be examined.
 *
 *****************************************************************************/

static int isValidPosition(const PorousSolid* solid, const Grid* grid, const double r[3]) {

  int ix, iy, iz;
  grid_cell_of(grid, r, &ix, &iy, &iz);

  int cxs[3], cys[3], czs[3];
  int ncx = neighborCellIndices(ix, grid->nx, solid->periodic[X], cxs);
  int ncy = neighborCellIndices(iy, grid->ny, solid->periodic[Y], cys);
  int ncz = neighborCellIndices(iz, grid->nz, solid->periodic[Z], czs);

  for (int a = 0; a < ncx; a++) {
    for (int b = 0; b < ncy; b++) {
      for (int c = 0; c < ncz; c++) {
        const Cell* cell = &grid->cells[(cxs[a] * grid->ny + cys[b]) * grid->nz + czs[c]];
        for (int k = 0; k < cell->count; k++) {
          double dist = pairDistance(solid, r, solid->points[cell->idx[k]].r);
          if (dist < solid->min_distance) return 0;
        }
      }
    }
  }

  return 1;
}

/*****************************************************************************
 *
 *  generateRandomPoints
 *
 *  Rejection sampling: fill the cube with points at random subject to
 *  the minimum separation constraint (lbond - delta).
 *
 *****************************************************************************/

static void generateRandomPoints(PorousSolid* solid, Grid* grid, int target_npoints) {

  int max_attempts_per_point = 2000;
  long max_total_attempts = (long)max_attempts_per_point * (target_npoints > 0 ? target_npoints : 1);
  long attempts = 0;

  printf("\nFilling cube with random points...\n");

  while (solid->npoints < target_npoints && attempts < max_total_attempts) {

    double r[3];
    r[X] = solid->lmin[X] + randomDouble(0.0, solid->side);
    r[Y] = solid->lmin[Y] + randomDouble(0.0, solid->side);
    r[Z] = solid->lmin[Z] + randomDouble(0.0, solid->side);

    if (isValidPosition(solid, grid, r)) {
      int n = solid->npoints;
      solid->points[n].id = n + 1;
      solid->points[n].r[X] = r[X];
      solid->points[n].r[Y] = r[Y];
      solid->points[n].r[Z] = r[Z];
      solid->points[n].num_neighbors = 0;
      solid->points[n].is_fixed = 0;

      int ix, iy, iz;
      grid_cell_of(grid, r, &ix, &iy, &iz);
      grid_add(grid, ix, iy, iz, n);

      solid->npoints++;
    }

    attempts++;
  }

  printf("Placed %d of %d target points after %ld attempts\n",
         solid->npoints, target_npoints, attempts);

  if (solid->npoints < target_npoints) {
    printf("Warning: could not reach the target point count -- the minimum "
           "separation (%.4f) may be too large for the requested density.\n",
           solid->min_distance);
  }
}

/*****************************************************************************
 *
 *  markFixedFaces
 *
 *  Mark every point within `thickness` of one of the given faces as
 *  fixed (isfixedr/isfixedv/isfixedw -- see colloid_init_state), i.e.
 *  anchored in place for the whole simulation.
 *
 *****************************************************************************/

static void markFixedFaces(PorousSolid* solid, const int fixed_axis[], const int fixed_side[],
                            int n_fixed_faces, double thickness) {

  if (n_fixed_faces == 0) return;

  for (int k = 0; k < n_fixed_faces; k++) {
    int axis = fixed_axis[k];
    int side = fixed_side[k];
    double face_coord = (side < 0) ? solid->lmin[axis] : solid->lmin[axis] + solid->side;
    int nmarked = 0;

    for (int i = 0; i < solid->npoints; i++) {
      double dist_to_face = fabs(solid->points[i].r[axis] - face_coord);
      if (dist_to_face <= thickness) {
        solid->points[i].is_fixed = 1;
        nmarked++;
      }
    }

    printf("Marked %d points as fixed near the %c%s face\n",
           nmarked, "xyz"[axis], side < 0 ? "-" : "+");
  }
}

/*****************************************************************************
 *
 *  alreadyBonded
 *
 *****************************************************************************/

static int alreadyBonded(const Point3D* p, int other_id) {
  for (int k = 0; k < p->num_neighbors; k++) {
    if (p->neighbors[k] == other_id) return 1;
  }
  return 0;
}

/*****************************************************************************
 *
 *  createBonds
 *
 *  Bond every pair of points whose separation lies within
 *  [lbond - delta, lbond + delta] (minimum-image on periodic axes), up
 *  to max_links bonds per point.
 *
 *****************************************************************************/

static void createBonds(PorousSolid* solid, Grid* grid) {

  int bonds_added = 0;
  double dmax = solid->lbond + solid->delta;

  printf("\nCreating bonds (target length %.4f +/- %.4f)...\n", solid->lbond, solid->delta);

  for (int i = 0; i < solid->npoints; i++) {

    Point3D* pi = &solid->points[i];
    if (pi->num_neighbors >= solid->max_links) continue;

    int ix, iy, iz;
    grid_cell_of(grid, pi->r, &ix, &iy, &iz);

    int cxs[3], cys[3], czs[3];
    int ncx = neighborCellIndices(ix, grid->nx, solid->periodic[X], cxs);
    int ncy = neighborCellIndices(iy, grid->ny, solid->periodic[Y], cys);
    int ncz = neighborCellIndices(iz, grid->nz, solid->periodic[Z], czs);

    for (int a = 0; a < ncx && pi->num_neighbors < solid->max_links; a++) {
      for (int b = 0; b < ncy && pi->num_neighbors < solid->max_links; b++) {
        for (int c = 0; c < ncz && pi->num_neighbors < solid->max_links; c++) {

          const Cell* cell = &grid->cells[(cxs[a] * grid->ny + cys[b]) * grid->nz + czs[c]];
          for (int k = 0; k < cell->count && pi->num_neighbors < solid->max_links; k++) {
            int j = cell->idx[k];
            if (j <= i) continue; /* each pair considered once, from lower index */

            Point3D* pj = &solid->points[j];
            if (pj->num_neighbors >= solid->max_links) continue;
            if (alreadyBonded(pi, pj->id)) continue;

            double dist = pairDistance(solid, pi->r, pj->r);
            if (dist >= solid->min_distance && dist <= dmax) {
              pi->neighbors[pi->num_neighbors++] = pj->id;
              pj->neighbors[pj->num_neighbors++] = pi->id;
              bonds_added++;
            }
          }
        }
      }
    }
  }

  printf("Added %d bonds\n", bonds_added);
}

/*****************************************************************************
 *
 *  printNetworkStats
 *
 *****************************************************************************/

static void printNetworkStats(const PorousSolid* solid) {

  int min_bonds = NBOND_MAX + 1;
  int max_bonds = 0;
  double avg_bonds = 0.0;
  int no_bonds = 0;
  int n_fixed = 0;

  double min_dist = DBL_MAX;
  double max_dist = 0.0;
  double avg_dist = 0.0;
  int nbond_pairs = 0;

  for (int i = 0; i < solid->npoints; i++) {
    int num = solid->points[i].num_neighbors;
    if (num < min_bonds) min_bonds = num;
    if (num > max_bonds) max_bonds = num;
    avg_bonds += num;
    if (num == 0) no_bonds++;
    if (solid->points[i].is_fixed) n_fixed++;

    for (int k = 0; k < num; k++) {
      int other_id = solid->points[i].neighbors[k];
      if (other_id > solid->points[i].id) {
        const Point3D* other = &solid->points[other_id - 1];
        double dist = pairDistance(solid, solid->points[i].r, other->r);
        if (dist < min_dist) min_dist = dist;
        if (dist > max_dist) max_dist = dist;
        avg_dist += dist;
        nbond_pairs++;
      }
    }
  }

  if (solid->npoints > 0) avg_bonds /= solid->npoints;
  if (nbond_pairs > 0) avg_dist /= nbond_pairs;

  printf("\nNetwork statistics:\n");
  printf("  Total points: %d\n", solid->npoints);
  printf("  Fixed points: %d\n", n_fixed);
  printf("  Cube side: %.4f\n", solid->side);
  printf("  Bond length target: %.4f (delta: %.4f)\n", solid->lbond, solid->delta);
  printf("  Min bonds per point: %d\n", solid->npoints > 0 ? min_bonds : 0);
  printf("  Max bonds per point: %d\n", max_bonds);
  printf("  Avg bonds per point: %.4f\n", avg_bonds);
  printf("  Points with no bonds: %d\n", no_bonds);
  printf("  Total bonds: %d\n", nbond_pairs);
  if (nbond_pairs > 0) {
    printf("  Min bonded distance: %.6f\n", min_dist);
    printf("  Max bonded distance: %.6f\n", max_dist);
    printf("  Avg bonded distance: %.6f\n", avg_dist);
  }
}

/*****************************************************************************
 *
 *  colloid_init_state
 *
 *****************************************************************************/

static void colloid_init_state(colloid_state_t* state, const PorousSolid* solid,
                                double a0, double ah, double q0, double epsilon,
                                double sa, double saf, int bc, int shape,
                                double drmax, double al) {

  for (int j = 0; j < solid->npoints; j++) {
    const Point3D* p = &solid->points[j];

    state[j].index = p->id;
    state[j].rebuild = 1;
    state[j].a0 = a0;
    state[j].ah = ah;
    state[j].q0 = q0;
    state[j].q1 = 0.0;
    state[j].epsilon = epsilon;
    state[j].sa = sa;
    state[j].saf = saf;
    state[j].b1 = 0.0;
    state[j].b2 = 0.0;
    state[j].m[X] = 0.0;
    state[j].m[Y] = 0.0;
    state[j].m[Z] = 0.0;
    state[j].bc = bc;
    state[j].shape = shape;
    state[j].al = al;
    state[j].rng = p->id;
    state[j].r[X] = p->r[X];
    state[j].r[Y] = p->r[Y];
    state[j].r[Z] = p->r[Z];
    state[j].nbonds = p->num_neighbors;

    for (int i = 0; i < state[j].nbonds; i++) {
      state[j].bond[i] = p->neighbors[i];
    }

    if (state[j].nbonds > 0) {
      state[j].nangles = (int)(0.5 * (double)(state[j].nbonds * (state[j].nbonds - 1)));
    }
    else {
      state[j].nangles = 0;
    }

    if (p->is_fixed) {
      state[j].isfixedr = 1;
      state[j].isfixedv = 1;
      state[j].isfixedw = 1;
      state[j].isfixedrxyz[X] = 1;
      state[j].isfixedrxyz[Y] = 1;
      state[j].isfixedrxyz[Z] = 1;
      state[j].isfixedvxyz[X] = 1;
      state[j].isfixedvxyz[Y] = 1;
      state[j].isfixedvxyz[Z] = 1;
    }
  }
}

/*****************************************************************************
 *
 *  colloid_init_write_file
 *
 *****************************************************************************/

static void colloid_init_write_file(const int nc, const colloid_state_t* pc,
                                     const char* filename, const int form) {
  int n;
  FILE* fp;

  fp = util_fopen(filename, "w");
  if (fp == NULL) {
    printf("Could not open %s\n", filename);
    exit(0);
  }

  if (form == BINARY) {
    fwrite(&nc, sizeof(int), 1, fp);
  }
  else {
    fprintf(fp, "%22d\n", nc);
  }

  for (n = 0; n < nc; n++) {
    if (form == BINARY) {
      colloid_state_write_binary(pc + n, fp);
    }
    else {
      colloid_state_write_ascii(pc + n, fp);
    }
  }

  if (ferror(fp)) {
    perror("perror: ");
    printf("Error reported on write to %s\n", filename);
  }

  fclose(fp);
}

/*****************************************************************************
 *
 *  parseFixedFaceArg
 *
 *  Parses a comma-separated list of faces such as "z-,z+" or "x-" into
 *  axis (X/Y/Z) and side (-1/+1) pairs. Returns the number parsed.
 *
 *****************************************************************************/

static int parseFixedFaceArg(const char* arg, int axis_out[], int side_out[], int max_faces) {

  if (arg == NULL || strcmp(arg, "none") == 0 || arg[0] == '\0') return 0;

  char* buf = strdup(arg);
  assert(buf != NULL);

  int n = 0;
  char* token = strtok(buf, ",");
  while (token != NULL && n < max_faces) {
    /* Trim leading/trailing whitespace */
    while (*token == ' ') token++;
    size_t len = strlen(token);
    while (len > 0 && token[len - 1] == ' ') token[--len] = '\0';

    int axis = -1;
    int side = 0;
    if (len == 2) {
      char c0 = (char)tolower((unsigned char)token[0]);
      if (c0 == 'x') axis = X;
      else if (c0 == 'y') axis = Y;
      else if (c0 == 'z') axis = Z;

      if (token[1] == '-') side = -1;
      else if (token[1] == '+') side = +1;
    }

    if (axis == -1 || side == 0) {
      fprintf(stderr, "Error: invalid --fixed-face token '%s' (expected e.g. x-, y+, z-)\n", token);
      free(buf);
      exit(1);
    }

    axis_out[n] = axis;
    side_out[n] = side;
    n++;

    token = strtok(NULL, ",");
  }

  free(buf);
  return n;
}

/*****************************************************************************
 *
 *  clopt
 *
 *****************************************************************************/

static int clopt(int argc, char** argv,
                  int* side, double* lbond, double* delta, double* density,
                  int* max_links, unsigned int* seed,
                  double* irad, double* hrad, double* charge, double* permittivity,
                  double* sa, double* saf, double* drmax, double* offset,
                  int* periodic_x, int* periodic_y, int* periodic_z,
                  const char** fixed_face_arg, double* fixed_thickness) {
  int c;
  char* eptr;

  static struct option long_options[] = {
    {"side",             required_argument, 0, 'L'}, /* Cube side length (also grid size) */
    {"lbond",            required_argument, 0, 'l'}, /* Mean bond length */
    {"delta",            required_argument, 0, 'e'}, /* Bond length tolerance */
    {"density",          required_argument, 0, 'd'}, /* Point density */
    {"max_links",        required_argument, 0, 'm'}, /* Max bonds per point */
    {"seed",             required_argument, 0, 'S'}, /* Random seed */
    {"irad",             required_argument, 0, 'i'}, /* Input radius */
    {"hrad",             required_argument, 0, 'h'}, /* Hydrodynamic radius */
    {"charge",           required_argument, 0, 'c'}, /* Charge q0 */
    {"permittivity",     required_argument, 0, 'p'}, /* Dielectric permittivity */
    {"sa",                required_argument, 0, 's'}, /* Surface area */
    {"saf",               required_argument, 0, 'f'}, /* Surface area to fluid */
    {"drmax",             required_argument, 0, 'r'}, /* Max displacement */
    {"offset",            required_argument, 0, 'o'}, /* Subgrid offset parameter */
    {"periodic-x",        required_argument, 0, 'x'}, /* 0 or 1: periodic bonding in x */
    {"periodic-y",        required_argument, 0, 'y'}, /* 0 or 1: periodic bonding in y */
    {"periodic-z",        required_argument, 0, 'z'}, /* 0 or 1: periodic bonding in z */
    {"fixed-face",        required_argument, 0, 'F'}, /* e.g. "z-" or "x-,x+" or "none" */
    {"fixed-thickness",   required_argument, 0, 'T'}, /* layer thickness for fixed faces */
    {0, 0, 0, 0}
  };

  while (1) {
    int option_index = 0;
    c = getopt_long(argc, argv, "L:l:e:d:m:S:i:h:c:p:s:f:r:o:x:y:z:F:T:", long_options, &option_index);
    if (c == -1) break;

    switch (c) {
    case 'L': *side = atoi(optarg); break;
    case 'l': *lbond = strtod(optarg, &eptr); break;
    case 'e': *delta = strtod(optarg, &eptr); break;
    case 'd': *density = strtod(optarg, &eptr); break;
    case 'm': *max_links = atoi(optarg); break;
    case 'S': *seed = (unsigned int)strtoul(optarg, &eptr, 10); break;
    case 'i': *irad = strtod(optarg, &eptr); break;
    case 'h': *hrad = strtod(optarg, &eptr); break;
    case 'c': *charge = strtod(optarg, &eptr); break;
    case 'p': *permittivity = strtod(optarg, &eptr); break;
    case 's': *sa = strtod(optarg, &eptr); break;
    case 'f': *saf = strtod(optarg, &eptr); break;
    case 'r': *drmax = strtod(optarg, &eptr); break;
    case 'o': *offset = strtod(optarg, &eptr); break;
    case 'x': *periodic_x = atoi(optarg); break;
    case 'y': *periodic_y = atoi(optarg); break;
    case 'z': *periodic_z = atoi(optarg); break;
    case 'F': *fixed_face_arg = optarg; break;
    case 'T': *fixed_thickness = strtod(optarg, &eptr); break;
    default: abort();
    }
  }

  if (optind < argc) {
    printf("non-option ARGV-elements: ");
    while (optind < argc) printf("%s ", argv[optind++]);
    putchar('\n');
  }

  return 0;
}
