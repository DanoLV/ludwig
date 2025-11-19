/*****************************************************************************
 *
 *  microgel_poly_cross_density.c
 *
 *  Produces a file containing microgel polymer information suitable for input 
 *  into the main code. The microgel is a network of monomers distributed within 
 *  a spherical volume. Monomer positions are initialized using a Fibonacci 
 *  approach to distribute them evenly over the sphere. Then, linear chains are 
 *  grown inward from the surface monomers, with an average chain length adjusted 
 *  to match the target density. Afterward, crosslinks are created between chains 
 *  according to the specified crosslinking density.
 *
 *  For compilation instructions see the Makefile.
 *
 *  $ make microgel_poly_cross_density
 *
 *  $ ./microgel_poly_cross_density
 *
 *  should produce a file config.cds.init.001-001 in the specified format.
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  Contributing authors
 *  Daniel La Valle (daniel.lavalle@ub.edu)
 *  (c) 2025- University of Barcelona
 *
 *  From original program
 *  Kai Qi (kai.qi@epfl.ch)
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *  (c) 2012-2023 The University of Edinburgh
 *  (c) 2020- Swiss Federal Institute of Technology Lausanne
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

#define M_PI 3.14159265358979323846
#define MAX_BONDS 6
#define MIN_BONDS 1

/* Structure to represent a point in 3D space */
typedef struct
{
  int id;                 /* Unique identifier for each point */
  double x, y, z;         /* 3D coordinates */
  int* neighbors;         /* Array of neighbor IDs */
  int num_neighbors;      /* Number of neighbors */
  int neighbors_capacity; /* Capacity of neighbors array */
  int is_surface;         /* Flag to indicate if point is on surface */
} Point3D;

/* Structure to represent a pair of distance and index */
typedef struct
{
  double distance;
  int index;
} DistanceIndexPair;

/* Global variables for the microgel polymer */
typedef struct
{
  Point3D* monomers;         /* Array of monomers */
  int nmon;                  /* Number of monomers on surface */
  int interior_nmon;         /* Number of monomers in interior */
  int total_nmon;            /* Total number of monomers */
  double radius;             /* Sphere radius */
  double target_distance;    /* Target distance between nearest monomers */
  double distance_threshold; /* Threshold for determining neighbors */
  double min_distance;       /* Minimum allowed distance between monomers */
  double density;            /* Density of interior monomers */
  double avg_chain_length;    // Longitud media de las cadenas lineales
  double crosslink_density;   // Densidad de entrecruzamiento (probabilidad por monomero)
} MicrogelPolymer;

typedef struct {
  int monomer_index;  // ID o índice del monómero actual
  int* neighbors;     // Lista de índices de monómeros cercanos
  int num_neighbors;  // Cantidad de vecinos cercanos
} NearbyMonomerList;

/* Point Management Functions */
void initPoint3D(Point3D* point, double x, double y, double z, int id);
void freePoint3D(Point3D* point);
void addNeighbor(Point3D* point, int neighbor_id);
void normalizePoint(Point3D* point, double r);

/* Distance Calculation Functions */
double euclideanDistance(const Point3D* a, const Point3D* b);
double greatCircleDistance(const Point3D* a, const Point3D* b, double r);
int compareDistancePairs(const void* a, const void* b);

/* Polymer Management Functions */
void initMicrogelPolymer(MicrogelPolymer* polymer, int nmon, double target_distance, double distance_threshold);
void freeMicrogelPolymer(MicrogelPolymer* polymer);
void generateInitialMonomers(MicrogelPolymer* polymer);
void generatePolymerChains(MicrogelPolymer* polymer);
int areBonded(MicrogelPolymer* polymer, int idx1, int idx2);
void createCrosslinks(MicrogelPolymer* polymer);

/* Optimization Functions */
double calculateEnergyFunction(MicrogelPolymer* polymer);
double randomDouble(double min, double max);
void optimizeDistribution(MicrogelPolymer* polymer, int iterations);

/* Neighbor Finding and Analysis Functions */
void findNearestNeighbors(MicrogelPolymer* polymer);
// NearbyMonomerList* findNearbyMonomers(const MicrogelPolymer* polymer, double factor)
void verifyDistribution(const MicrogelPolymer* polymer);
void printConnectivity(const MicrogelPolymer* polymer);
void exportConnectivityData(const MicrogelPolymer* polymer, const char* filename);
int clopt(int argc, char** argv,
          int* gsize,
          int* nmon,
          double* density,
          double* lbond,
          double* distt,
          double* charge,
          double* permittivity,
          double* sa,
          double* saf,
          double* irad,
          double* hrad,
          double* drmax,
          double* crosslink_density,
          int* max_links,
          double* offset);

/* New function prototypes */
void initExtendedMicrogelPolymer(MicrogelPolymer* polymer,
                                 int surface_nmon,
                                 double target_distance,
                                 double distance_threshold,
                                 double density,
                                 double avg_chain_len,  // Nuevo parametro
                                 double crosslink_den); // Nuevo parametro
void generateInteriorMonomers(MicrogelPolymer* polymer);
int isValidPosition(const MicrogelPolymer* polymer, double x, double y, double z, int current_count);
void createFullPolymerNetwork(MicrogelPolymer* polymer);
void exportFullNetworkData(const MicrogelPolymer* polymer, const char* filename);
void printNetworkStats(const MicrogelPolymer* polymer);
void createBidirectionalBond(MicrogelPolymer* polymer, int monomer_idx_1, int monomer_idx_2);
void fillRemainingBonds(MicrogelPolymer* polymer);
/* mio 1 fin*/

enum format
{
  ASCII,
  BINARY
};

void colloid_init_write_file(const int nc, const colloid_state_t* pc,
                             const int form);

void colloid_init_state(double offset, colloid_state_t* state, MicrogelPolymer* polymer,
                        double a0,
                        double ah,
                        double q0,
                        double epsilon,
                        double sa,
                        double saf,
                        double q1,
                        double b1,
                        double b2,
                        int bc,
                        int shape,
                        double drmax,
                        double al);

/*****************************************************************************
 *
 *  main
 *
 *  You need to set the system parameters found directly below.
 *
 *****************************************************************************/

int main(int argc, char** argv)
{
  int gsize = 32; // Default grid size
  int ntotal[3] = { gsize, gsize, gsize }; /* Total system size (cf. input) */
  int periodic[3] = { 1, 1, 1 };  /* 0 = wall, 1 = periodic */
  int file_format = ASCII;

  int bc = COLLOID_BC_SUBGRID;
  int shape = COLLOID_SHAPE_SPHERE;

  int surface_monomers = 100; // Default surface monomers
  double lbond = 1.0;   // Bond length. Default value 
  double distt = 0.2;   // Distance threshold percentage
  double density = 0.5; // Density of interior monomers
  int    max_links = 6;

  double a0 = 0.1;      // Input radius. Default value  
  double ah = 0.1;      // Hydrodynamic radius. Default value  
  double al = 1.5;      // Offset parameter for subgrid particle. Default value  

  double drmax = 0.8;     // Maximum movement per time step 
  double q0 = 0.0;      // positive charge 
  double q1 = 0.0;      // negative charge 
  double epsilon = 0.0; // Dielectric permittivity
  double sa = 0.0;      // surface area
  double saf = 0.0;     // surface area to fluid
  double b1 = 0.00;
  double b2 = 0.00;

  // --- NUEVOS PARAMETROS POR DEFECTO ---
  double avg_chain_length = 5.0;  // Longitud media de las cadenas
  double crosslink_density = 0.1; // Probabilidad de entrecruzamiento (0.0 a 1.0)
  // ------------------------------------

  colloid_state_t* state;
  pe_t* pe;
  cs_t* cs;

  MPI_Init(&argc, &argv);

  pe_create(MPI_COMM_WORLD, PE_QUIET, &pe);

  /* This program is intended to be serial */
  assert(pe_mpi_size(pe) == 1);

  cs_create(pe, &cs);
  cs_ntotal_set(cs, ntotal);
  cs_periodicity_set(cs, periodic);
  cs_init(cs);

  // Initialize values from command line

  // clopt(argc, argv, &gsize, &nmon, &lbond, &distt, &a0, &ah, &al);
  clopt(argc, argv, &gsize, &surface_monomers, &density, &lbond, &distt, &q0, &epsilon, &sa, &saf, &a0, &ah, &drmax, &crosslink_density,
          &max_links, &al);

  /* Initialize random seed */
  srand((unsigned int)time(NULL));

  /* Initialize the microgel polymer */
  MicrogelPolymer polymer;
  // initMicrogelPolymer(&polymer, nmon, lbond, distt);
  initExtendedMicrogelPolymer(&polymer, surface_monomers, lbond,
                              distt, density, avg_chain_length, crosslink_density);

  printf("Generating polymer network with parameters:\n");
  printf("  Surface monomers: %d\n", polymer.nmon);
  printf("  Bond length: %.2f\n", polymer.target_distance);
  printf("  Distance threshold: %.2f\n", polymer.distance_threshold);
  printf("  Interior density: %.2f monomers/unit³\n", polymer.density);
  printf("  Sphere radius: %.4f\n", polymer.radius);

  /* Generate surface monomers using the Fibonacci sphere method */
  printf("\nDistributing surface monomers...\n");
  generateInitialMonomers(&polymer);

  /* Generate interior monomers (now via chains) and create bonds */
  printf("Creando red polimerica completa (cadenas y entrecruzamientos)...\n");
  createFullPolymerNetwork(&polymer); // Esta funcion orquesta la generacion de cadenas y entrecruzamientos

  /* Output the full network */
  printNetworkStats(&polymer);

  // Set values on colloid state array
  state = (colloid_state_t*)calloc(polymer.total_nmon, sizeof(colloid_state_t));
  assert(state != NULL);

  colloid_init_state((double)gsize / 2, state, &polymer, a0,
                      ah,
                      q0,
                      q1,
                      epsilon,
                      sa,
                      saf,
                      b1,
                      b2,
                      bc,
                      shape,
                      drmax,
                      al);

  /* Free memory */
  freeMicrogelPolymer(&polymer);

  /* mio 1 fin*/

  /* Write out */
  colloid_init_write_file(polymer.total_nmon, state, file_format);

  free(state);

  cs_free(cs);
  pe_free(pe);
  MPI_Finalize();

  return 0;
}

/****************************************************************************
 *
 *  colloid_init_write_file
 *
 ****************************************************************************/

void colloid_init_write_file(const int nc, const colloid_state_t* pc,
                             const int form)
{
  int n;
  const char* filename = "config.cds.init.001-001";
  FILE* fp;

  fp = util_fopen(filename, "w");
  if (fp == NULL)
  {
    printf("Could not open %s\n", filename);
    exit(0);
  }

  if (form == BINARY)
  {
    fwrite(&nc, sizeof(int), 1, fp);
  }
  else
  {
    fprintf(fp, "%22d\n", nc);
  }

  for (n = 0; n < nc; n++)
  {
    if (form == BINARY)
    {
      colloid_state_write_binary(pc + n, fp);
    }
    else
    {
      colloid_state_write_ascii(pc + n, fp);
    }
  }

  if (ferror(fp))
  {
    perror("perror: ");
    printf("Error reported on write to %s\n", filename);
  }

  fclose(fp);

  return;
}

/********************************************************************************/
/* Added functions for microgel                                                 */
/********************************************************************************/

/* Function to initialize a Point3D structure */
void initPoint3D(Point3D* point, double x, double y, double z, int id)
{
  point->x = x;
  point->y = y;
  point->z = z;
  point->id = id;
  point->neighbors = NULL;
  point->num_neighbors = 0;
  point->neighbors_capacity = 0;
}

/* Function to free memory allocated for a Point3D structure */
void freePoint3D(Point3D* point)
{
  if (point->neighbors != NULL)
  {
    free(point->neighbors);
    point->neighbors = NULL;
  }
  point->num_neighbors = 0;
  point->neighbors_capacity = 0;
}

/* Function to add a neighbor to a Point3D */
void addNeighbor(Point3D* point, int neighbor_id)
{
  /* Reallocate memory if needed */
  if (point->num_neighbors >= point->neighbors_capacity)
  {
    int new_capacity = point->neighbors_capacity == 0 ? 4 : point->neighbors_capacity * 2;
    int* new_neighbors = (int*)malloc(new_capacity * sizeof(int));
    if (new_neighbors == NULL)
    {
      printf("Error: Memory allocation failed.\n");
      exit(1);
    }

    /* Copy existing neighbors */
    for (int i = 0; i < point->num_neighbors; i++)
    {
      new_neighbors[i] = point->neighbors[i];
    }

    /* Free old array and update pointers */
    if (point->neighbors != NULL)
    {
      free(point->neighbors);
    }
    point->neighbors = new_neighbors;
    point->neighbors_capacity = new_capacity;
  }

  /* Add the new neighbor */
  point->neighbors[point->num_neighbors] = neighbor_id;
  point->num_neighbors++;
}

/* Function to normalize a point to lie on a sphere of radius r */
void normalizePoint(Point3D* point, double r)
{
  double len = sqrt(point->x * point->x + point->y * point->y + point->z * point->z);
  if (len > 0)
  {
    point->x = r * point->x / len;
    point->y = r * point->y / len;
    point->z = r * point->z / len;
  }
}

/* Function to calculate Euclidean distance between two points */
double euclideanDistance(const Point3D* a, const Point3D* b)
{
  return sqrt(pow(a->x - b->x, 2) + pow(a->y - b->y, 2) + pow(a->z - b->z, 2));
}

/* Function to calculate great circle distance between two points on a sphere */
double greatCircleDistance(const Point3D* a, const Point3D* b, double r)
{
  double dot = a->x * b->x + a->y * b->y + a->z * b->z;
  /* Ensure dot product is within valid range for acos */
  if (dot > 1.0)
    dot = 1.0;
  if (dot < -1.0)
    dot = -1.0;
  return r * acos(dot);
}

/* Comparison function for qsort */
int compareDistancePairs(const void* a, const void* b)
{
  const DistanceIndexPair* pair_a = (const DistanceIndexPair*)a;
  const DistanceIndexPair* pair_b = (const DistanceIndexPair*)b;

  if (pair_a->distance < pair_b->distance)
    return -1;
  if (pair_a->distance > pair_b->distance)
    return 1;
  return 0;
}

/* Function to initialize the MicrogelPolymer structure */
void initMicrogelPolymer(MicrogelPolymer* polymer, int nmon, double target_distance, double distance_threshold)
{
  polymer->nmon = nmon;
  polymer->target_distance = target_distance;
  polymer->distance_threshold = distance_threshold;

  /* Calculate appropriate sphere radius based on the number of points and target distance */
  polymer->radius = target_distance * sqrt(nmon / (4.0 * M_PI));

  /* Allocate memory for monomers */
  polymer->monomers = (Point3D*)malloc(nmon * sizeof(Point3D));
  if (polymer->monomers == NULL)
  {
    printf("Error: Memory allocation failed.\n");
    exit(1);
  }

  /* Initialize all monomers */
  for (int i = 0; i < nmon; i++)
  {
    initPoint3D(&polymer->monomers[i], 0, 0, 0, i);
  }
}

/* Function to free memory allocated for the MicrogelPolymer structure */
void freeMicrogelPolymer(MicrogelPolymer* polymer)
{
  if (polymer->monomers != NULL)
  {
    for (int i = 0; i < polymer->nmon; i++)
    {
      freePoint3D(&polymer->monomers[i]);
    }
    free(polymer->monomers);
    polymer->monomers = NULL;
  }
}

/* Function to generate initial monomers using the Fibonacci sphere method */
void generateInitialMonomers(MicrogelPolymer* polymer)
{
  const double phi = M_PI * (3.0 - sqrt(5.0)); /* Golden angle in radians */

  for (int i = 0; i < polymer->nmon; i++)
  {
    double y = 1.0 - (i / (double)(polymer->nmon - 1)) * 2.0; /* y goes from 1 to -1 */
    double radius_at_y = sqrt(1.0 - y * y);                           /* radius at y */

    double theta = phi * i; /* Golden angle increment */

    double x = cos(theta) * radius_at_y;
    double z = sin(theta) * radius_at_y;

    initPoint3D(&polymer->monomers[i], x, y, z, i + 1);
    normalizePoint(&polymer->monomers[i], polymer->radius);
  }
}

/* Helper function for optimization */
double calculateEnergyFunction(MicrogelPolymer* polymer)
{
  double energy = 0.0;
  for (int i = 0; i < polymer->nmon; i++)
  {
    for (int j = i + 1; j < polymer->nmon; j++)
    {
      /* We still use great circle distance for optimization */
      double dist = greatCircleDistance(&polymer->monomers[i], &polymer->monomers[j], polymer->radius);
      /* Penalize distances that differ from target */
      energy += pow(dist - polymer->target_distance, 2);
    }
  }
  return energy;
}

/* Optimize point distribution using a simple relaxation algorithm */
void optimizeDistribution(MicrogelPolymer* polymer, int iterations)
{
  /* Seed the random number generator */
  srand((unsigned int)time(NULL));

  /* Create a copy of the current points for working */
  Point3D* best_points = (Point3D*)malloc(polymer->nmon * sizeof(Point3D));
  Point3D* new_points = (Point3D*)malloc(polymer->nmon * sizeof(Point3D));
  if (best_points == NULL || new_points == NULL)
  {
    printf("Error: Memory allocation failed.\n");
    exit(1);
  }

  for (int i = 0; i < polymer->nmon; i++)
  {
    initPoint3D(&best_points[i], polymer->monomers[i].x, polymer->monomers[i].y, polymer->monomers[i].z, polymer->monomers[i].id);
  }

  double best_energy = calculateEnergyFunction(polymer);
  double temperature = 1.0;
  double cooling_rate = 0.995;

  for (int iter = 0; iter < iterations; iter++)
  {
    /* Create a copy of the current points */
    for (int i = 0; i < polymer->nmon; i++)
    {
      initPoint3D(&new_points[i], polymer->monomers[i].x, polymer->monomers[i].y, polymer->monomers[i].z, polymer->monomers[i].id);
    }

    /* Perturb points slightly */
    for (int i = 0; i < polymer->nmon; i++)
    {
      new_points[i].x += randomDouble(-0.1, 0.1) * temperature;
      new_points[i].y += randomDouble(-0.1, 0.1) * temperature;
      new_points[i].z += randomDouble(-0.1, 0.1) * temperature;
      normalizePoint(&new_points[i], polymer->radius);
    }

    /* Swap points to calculate energy */
    Point3D* temp = polymer->monomers;
    polymer->monomers = new_points;
    new_points = temp;

    /* Calculate energy of the new configuration */
    double new_energy = calculateEnergyFunction(polymer);

    /* Accept or reject the new configuration */
    if (new_energy < best_energy)
    {
      /* Update best configuration */
      best_energy = new_energy;
      for (int i = 0; i < polymer->nmon; i++)
      {
        best_points[i].x = polymer->monomers[i].x;
        best_points[i].y = polymer->monomers[i].y;
        best_points[i].z = polymer->monomers[i].z;
      }
    }
    else
    {
      /* Probability to accept worse solution decreases with temperature */
      double accept_probability = exp((best_energy - new_energy) / temperature);
      if (randomDouble(0.0, 1.0) < accept_probability)
      {
        /* Accept worse solution */
        best_energy = new_energy;
        for (int i = 0; i < polymer->nmon; i++)
        {
          best_points[i].x = polymer->monomers[i].x;
          best_points[i].y = polymer->monomers[i].y;
          best_points[i].z = polymer->monomers[i].z;
        }
      }
      else
      {
        /* Revert to previous state */
        temp = polymer->monomers;
        polymer->monomers = new_points;
        new_points = temp;
      }
    }

    /* Cool down */
    temperature *= cooling_rate;
  }

  /* Copy best configuration back to monomers */
  for (int i = 0; i < polymer->nmon; i++)
  {
    polymer->monomers[i].x = best_points[i].x;
    polymer->monomers[i].y = best_points[i].y;
    polymer->monomers[i].z = best_points[i].z;
  }

  /* Free temporary arrays */
  for (int i = 0; i < polymer->nmon; i++)
  {
    freePoint3D(&best_points[i]);
    freePoint3D(&new_points[i]);
  }
  free(best_points);
  free(new_points);
}

/* Find nearest neighbors for each monomer using Euclidean distance */
void findNearestNeighbors(MicrogelPolymer* polymer)
{
  /* Clear any existing neighbor information */
  for (int i = 0; i < polymer->nmon; i++)
  {
    if (polymer->monomers[i].neighbors != NULL)
    {
      free(polymer->monomers[i].neighbors);
      polymer->monomers[i].neighbors = NULL;
    }
    polymer->monomers[i].num_neighbors = 0;
    polymer->monomers[i].neighbors_capacity = 0;
  }

  /* Allocate memory for distance-index pairs */
  DistanceIndexPair* distances = (DistanceIndexPair*)malloc((polymer->nmon - 1) * sizeof(DistanceIndexPair));
  if (distances == NULL)
  {
    printf("Error: Memory allocation failed.\n");
    exit(1);
  }

  /* Calculate straight-line distances and find neighbors */
  for (int i = 0; i < polymer->nmon; i++)
  {
    /* Store distances to all other monomers */
    int pair_idx = 0;
    for (int j = 0; j < polymer->nmon; j++)
    {
      if (i != j)
      {
        /* Use Euclidean distance */
        distances[pair_idx].distance = euclideanDistance(&polymer->monomers[i], &polymer->monomers[j]);
        distances[pair_idx].index = j;
        pair_idx++;
      }
    }

    /* Sort by distance */
    qsort(distances, polymer->nmon - 1, sizeof(DistanceIndexPair), compareDistancePairs);

    /* Add neighbors that are within the threshold of the target distance */
    for (int k = 0; k < polymer->nmon - 1; k++)
    {
      if (distances[k].distance <= polymer->target_distance * (1.0 + polymer->distance_threshold))
      {
        addNeighbor(&polymer->monomers[i], polymer->monomers[distances[k].index].id);
      }
    }
  }

  /* Free memory */
  free(distances);
}

/* Verify the distribution's quality by checking distances */
void verifyDistribution(const MicrogelPolymer* polymer)
{
  double min_dist = DBL_MAX;
  double max_dist = 0.0;
  double avg_dist = 0.0;
  int count = 0;

  /* Calculate statistics for neighbor distances */
  double min_neighbor_dist = DBL_MAX;
  double max_neighbor_dist = 0.0;
  double avg_neighbor_dist = 0.0;
  int neighbor_count = 0;

  for (int i = 0; i < polymer->total_nmon; i++)
  {
    for (int j = i + 1; j < polymer->total_nmon; j++)
    {
      /* Using Euclidean distance */
      double dist = euclideanDistance(&polymer->monomers[i], &polymer->monomers[j]);
      if (dist < min_dist)
        min_dist = dist;
      if (dist > max_dist)
        max_dist = dist;
      avg_dist += dist;
      count++;
    }

    /* Check actual neighbors */
    for (int k = 0; k < polymer->monomers[i].num_neighbors; k++)
    {
      int neighbor_id = polymer->monomers[i].neighbors[k];

      /* Find the index corresponding to the neighbor ID */
      int j = -1;
      for (int m = 0; m < polymer->total_nmon; m++)
      {
        if (polymer->monomers[m].id == neighbor_id)
        {
          j = m;
          break;
        }
      }

      if (j != -1)
      {
        double dist = euclideanDistance(&polymer->monomers[i], &polymer->monomers[j]);
        if (dist < min_neighbor_dist)
          min_neighbor_dist = dist;
        if (dist > max_neighbor_dist)
          max_neighbor_dist = dist;
        avg_neighbor_dist += dist;
        neighbor_count++;
      }
    }
  }

  avg_dist /= count;
  avg_neighbor_dist = neighbor_count > 0 ? avg_neighbor_dist / neighbor_count : 0;

  /* Count average connections per monomer */
  double avg_connections = 0.0;
  for (int i = 0; i < polymer->total_nmon; i++)
  {
    avg_connections += polymer->monomers[i].num_neighbors;
  }
  avg_connections /= polymer->total_nmon;

  printf("Distribution quality metrics:\n");
  printf("  Target distance (ldist): %.4f\n", polymer->target_distance);
  printf("  Minimum Euclidean distance between any monomers: %.4f\n", min_dist);
  printf("  Maximum Euclidean distance between any monomers: %.4f\n", max_dist);
  printf("  Average Euclidean distance between all monomers: %.4f\n", avg_dist);
  printf("  Sphere radius: %.4f\n", polymer->radius);
  printf("\nConnectivity metrics:\n");
  printf("  Minimum Euclidean distance between connected monomers: %.8f\n", min_neighbor_dist);
  printf("  Maximum Euclidean distance between connected monomers: %.8f\n", max_neighbor_dist);
  printf("  Average Euclidean distance between connected monomers: %.8f\n", avg_neighbor_dist);
  printf("  Average connections per monomer: %.8f\n", avg_connections);
}

/* Print connectivity for each monomer */
void printConnectivity(const MicrogelPolymer* polymer)
{
  printf("\nMonomer connectivity:\n");
  for (int i = 0; i < polymer->total_nmon; i++)
  {
    printf("Monomer %d connects to: ", polymer->monomers[i].id);
    if (polymer->monomers[i].num_neighbors == 0)
    {
      printf("none");
    }
    else
    {
      for (int j = 0; j < polymer->monomers[i].num_neighbors; j++)
      {
        printf("%d", polymer->monomers[i].neighbors[j]);
        if (j < polymer->monomers[i].num_neighbors - 1)
        {
          printf(", ");
        }
      }
    }
    printf("\n");
  }
}

/* Export connectivity data for each monomer */
void exportConnectivityData(const MicrogelPolymer* polymer, const char* filename)
{
  FILE* outfile = fopen(filename, "w");
  if (outfile == NULL)
  {
    printf("Failed to open file: %s\n", filename);
    return;
  }

  fprintf(outfile, "# Microgel Polymer Connectivity Data\n");
  fprintf(outfile, "# Number of monomers: %d\n", polymer->nmon);
  fprintf(outfile, "# Target distance: %.4f\n", polymer->target_distance);
  fprintf(outfile, "# Sphere radius: %.4f\n\n", polymer->radius);

  fprintf(outfile, "# Format: Monomer_ID x y z [Neighbor_IDs...]\n");

  for (int i = 0; i < polymer->total_nmon; i++)
  {
    fprintf(outfile, "%d %.6f %.6f %.6f",
            polymer->monomers[i].id,
            polymer->monomers[i].x,
            polymer->monomers[i].y,
            polymer->monomers[i].z);

    for (int j = 0; j < polymer->monomers[i].num_neighbors; j++)
    {
      fprintf(outfile, " %d", polymer->monomers[i].neighbors[j]);
    }
    fprintf(outfile, "\n");
  }

  fclose(outfile);
  printf("Connectivity data exported to %s\n", filename);
}

int clopt(int argc, char** argv,
          int* gsize,
          int* nmon,
          double* density,
          double* lbond,
          double* distt,
          double* charge,
          double* epsilon,
          double* sa,
          double* saf,
          double* irad,
          double* hrad,
          double* drmax,
          double* crosslink_density,
          int* max_links,
          double* offset)
{
  int c;
  char* eptr;

  /* Structure to represent command line option */
  static struct option long_options[] =
  {
    {"nmon",              required_argument,       0, 'n'}, // Number of monomers
    {"lbond",             required_argument,       0, 'l'}, // Bong length
    {"distt",             required_argument,       0, 't'}, // Distance threshold
    {"gsize",             required_argument,       0, 'g'}, // Grid size (x=y=z)
    {"irad",              required_argument,       0, 'i'}, // Input radius 
    {"hrad",              required_argument,       0, 'h'}, // Hydrodynamic radius 
    {"offset",            required_argument,       0, 'o'}, // Offset parameter for subgrid particle 
    {"idensity",          required_argument,       0, 'd'}, // Interior density 
    {"charge",            required_argument,       0, 'c'}, // Interior density 
    {"permittivity",      required_argument,       0, 'p'}, // Dielectric permittivity
    {"sa",                required_argument,       0, 's'}, // Surface area
    {"saf",               required_argument,       0, 'f'}, // Surface area to fluid
    {"drmax",             required_argument,       0, 'r'}, // Max displacement for colloids
    {"crosslink_density", required_argument,  0, 'x'}, // crosslink density 
    {"max_links",         required_argument,          0, 'm'}, // Max links 
    {0, 0, 0, 0}
  };

  while (1)
  {
    /* getopt_long stores the option index here. */
    int option_index = 0;

    c = getopt_long(argc, argv, "abc:d:f:",
                     long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1)
      break;

    switch (c)
    {
    case 'n':
      *nmon = atoi(optarg);
      break;

    case 'l':
      *lbond = strtod(optarg, &eptr);
      break;

    case 't':
      *distt = strtod(optarg, &eptr);
      break;

    case 'g':
      *gsize = strtod(optarg, &eptr);
      break;

    case 'i':
      *irad = strtod(optarg, &eptr);
      break;

    case 'h':
      *hrad = strtod(optarg, &eptr);
      break;

    case 'o':
      *offset = strtod(optarg, &eptr);
      break;

    case 'd':
      *density = strtod(optarg, &eptr);
      break;

    case 'c':
      *charge = strtod(optarg, &eptr);
      break;

    case 'p':
      *epsilon = strtod(optarg, &eptr);
      break;

    case 's':
      *sa = strtod(optarg, &eptr);
      break;

    case 'f':
      *saf = strtod(optarg, &eptr);
      break;

    case 'r':
      *drmax = strtod(optarg, &eptr);
      break;

    case 'x':
      *crosslink_density = strtod(optarg, &eptr);
      break;

    case 'm':
      *max_links = atoi(optarg);
      break;

    default:
      abort();
    }
  }

  /* Print any remaining command line arguments (not options). */
  if (optind < argc)
  {
    printf("non-option ARGV-elements: ");
    while (optind < argc)
      printf("%s ", argv[optind++]);
    putchar('\n');
  }

  return 0;
}

void colloid_init_state(double offset, colloid_state_t* state, MicrogelPolymer* polymer,
                        double a0,
                        double ah,
                        double q0,
                        double q1,
                        double epsilon,
                        double sa,
                        double saf,
                        double b1,
                        double b2,
                        int bc,
                        int shape,
                        double drmax,
                        double al)
{
  // int a;
  for (int j = 0; j < polymer->total_nmon; j++)
  {
    state[j].index = polymer->monomers[j].id;
    state[j].rebuild = 1;
    state[j].a0 = a0;
    state[j].ah = ah;
    state[j].q0 = q0;
    state[j].q1 = q1;
    state[j].epsilon = epsilon;
    state[j].sa = sa;
    state[j].saf = saf;
    state[j].b1 = b1;
    state[j].b2 = b2;
    state[j].m[X] = 0.0;
    state[j].m[Y] = 0.0;
    state[j].m[Z] = 0.0;
    state[j].bc = bc;
    state[j].shape = shape;
    state[j].al = al;
    state[j].rng = polymer->monomers[j].id;
    state[j].r[X] = polymer->monomers[j].x + offset;
    state[j].r[Y] = polymer->monomers[j].y + offset;
    state[j].r[Z] = polymer->monomers[j].z + offset;
    state[j].nbonds = polymer->monomers[j].num_neighbors;
    for (int i = 0; i < state[j].nbonds; i++)
    {
      state[j].bond[i] = polymer->monomers[j].neighbors[i];
    }
    if (state[j].nbonds > 0) {
      // a = (int)(0.5 * (double)(state[j].nbonds * (state[j].nbonds - 1)));
      // printf("angulos: %d, enlaces: %d", a, state[j].nbonds);
      state[j].nangles = (int)(0.5 * (double)(state[j].nbonds * (state[j].nbonds - 1)));
    }
    else
      state[j].nangles = 0;
  }
}

/* Initialize the microgel polymer with extended functionality */
void initExtendedMicrogelPolymer(MicrogelPolymer* polymer,
                                 int surface_monomers,
                                 double bond_length,
                                 double dist_threshold,
                                 double target_density,
                                 double avg_chain_len, // Nuevo parametro
                                 double crosslink_den) // Nuevo parametro
{
  assert(polymer);

  polymer->nmon = surface_monomers;
  polymer->interior_nmon = 0; // Se llenara mas tarde por las cadenas
  polymer->total_nmon = 0;    // Se llenara mas tarde

  polymer->target_distance = bond_length;
  polymer->distance_threshold = dist_threshold;
  polymer->min_distance = polymer->target_distance * (1.0 - polymer->distance_threshold);
  polymer->density = target_density;

  // Calculo del radio: Ajustado para la nueva logica de densidad
  // El radio se calculara en base a la densidad total deseada y el numero de monomeros
  // Inicialmente, se puede estimar basado en los monomeros superficiales
  // o dejarlo para que se ajuste dinamicamente.
  // Por ahora, lo dejaremos como una estimacion o lo pasamos como parametro si se prefiere.
  // Para este ejemplo, calculamos un radio inicial basado en una densidad 'ideal'
  // que se ajustará al total de monomeros generados.
  // polymer->radius = pow(((double)polymer->nmon / polymer->density) * (3.0 / (4.0 * M_PI)), 1.0 / 3.0);

  /* Calculate appropriate sphere radius based on the number of surface points and target distance */
  polymer->radius = polymer->target_distance * sqrt(polymer->nmon / (4.0 * M_PI)) * 1.1;
  printf("Initial radius estimate: %.4f\n", polymer->radius);


  // Asignacion de memoria para monomeros (ajusta el tamaño si es necesario)
  // Es una estimacion, ya que el numero total exacto de monomeros interiores
  // dependera de la generacion de cadenas. Podemos sobredimensionar.
  int estimated_total_monomers = (int)(polymer->density * (4.0 / 3.0) * M_PI * pow(polymer->radius, 3) * 2.0); // Doble de la estimacion inicial
  if (estimated_total_monomers < polymer->nmon + 100) estimated_total_monomers = polymer->nmon + 100; // Asegurar espacio minimo

  polymer->monomers = (Point3D*)malloc(estimated_total_monomers * sizeof(Point3D));
  if (polymer->monomers == NULL) {
    fprintf(stderr, "Error: No se pudo asignar memoria para los monomeros.\n");
    exit(1);
  }

  // Asignacion de los nuevos parametros
  polymer->avg_chain_length = (estimated_total_monomers - polymer->nmon) / polymer->nmon; //LONGITUD MEDIA ASUMIENDO QUE SALE UNA CADENA DE CAdA MONOMERO DE SUPERFICIE
  polymer->crosslink_density = crosslink_den;

  // Inicializar todos los monomeros a un estado por defecto (opcional)
  for (int i = 0; i < estimated_total_monomers; i++) {
    initPoint3D(&polymer->monomers[i], 0.0, 0.0, 0.0, -1); // ID -1 para no inicializado
  }

  return;
}

/* Generate random monomers inside the sphere based on density */
void generateInteriorMonomers(MicrogelPolymer* polymer)
{
  int count = polymer->nmon;  /* Start from after surface monomers */
  int max_attempts = 1000;

  /* Calculate expected number of interior monomers */
  double sphere_volume = (4.0 / 3.0) * M_PI * pow(polymer->radius, 3);
  int target_interior = (int)(polymer->density * sphere_volume) - polymer->nmon;

  printf("Target interior monomers: %d\n", target_interior);

  /* Generate interior monomers using rejection sampling */
  int interior_added = 0;
  int total_attempts = 0;

  while (interior_added < target_interior && total_attempts < max_attempts * target_interior) {
    /* Generate random point in cube then reject if outside sphere */
    double x = (2.0 * rand() / RAND_MAX - 1.0) * polymer->radius;
    double y = (2.0 * rand() / RAND_MAX - 1.0) * polymer->radius;
    double z = (2.0 * rand() / RAND_MAX - 1.0) * polymer->radius;

    double dist_from_center = sqrt(x * x + y * y + z * z);

    /* Check if point is inside sphere but not too close to surface */
    if (dist_from_center < polymer->radius * 0.95) {
      /* Check if point is valid (not too close to other monomers) */
      if (isValidPosition(polymer, x, y, z, count)) {
        polymer->monomers[count].x = x;
        polymer->monomers[count].y = y;
        polymer->monomers[count].z = z;
        polymer->monomers[count].id = count + 1;
        polymer->monomers[count].is_surface = 0;
        count++;
        interior_added++;
      }
    }

    total_attempts++;
  }

  polymer->interior_nmon = interior_added;
  polymer->total_nmon = polymer->nmon + interior_added;

  printf("Added %d interior monomers after %d attempts\n",
  interior_added, total_attempts);

  /* Mark surface monomers */
  for (int i = 0; i < polymer->nmon; i++) {
    polymer->monomers[i].is_surface = 1;
  }
}

/* Check if a position is valid (not too close to existing monomers) */
int isValidPosition(const MicrogelPolymer* polymer, double x, double y, double z, int current_count)
{
  Point3D temp_point;
  initPoint3D(&temp_point, x, y, z, -1);

  for (int i = 0; i < current_count; i++) {
    double dist = euclideanDistance(&temp_point, &polymer->monomers[i]);

    if (dist < polymer->min_distance) {
      freePoint3D(&temp_point);
      return 0;  /* Too close to an existing monomer */
    }
  }

  freePoint3D(&temp_point);
  return 1;  /* Position is valid */
}

/* Function to create a bidirectional bond between two monomers */
void createBidirectionalBond(MicrogelPolymer* polymer, int monomer_idx_1, int monomer_idx_2)
{
  /* Check if distance is within threshold before creating bond */
  double distance = euclideanDistance(&polymer->monomers[monomer_idx_1], &polymer->monomers[monomer_idx_2]);

  /* Only create bond if within distance threshold */
  if (distance > polymer->target_distance * (1.0 + polymer->distance_threshold)) {
    return; /* Distance too large, don't create bond */
  }

  /* Only create bond if it doesn't already exist and both monomers have room for more bonds */
  if (polymer->monomers[monomer_idx_1].num_neighbors < MAX_BONDS &&
      polymer->monomers[monomer_idx_2].num_neighbors < MAX_BONDS) {

    /* Check if bond already exists in first monomer */
    int bond_exists = 0;
    for (int b = 0; b < polymer->monomers[monomer_idx_1].num_neighbors; b++) {
      if (polymer->monomers[monomer_idx_1].neighbors[b] == polymer->monomers[monomer_idx_2].id) {
        bond_exists = 1;
        break;
      }
    }

    if (!bond_exists) {
      /* Add bidirectional bonds */
      addNeighbor(&polymer->monomers[monomer_idx_1], polymer->monomers[monomer_idx_2].id);
      addNeighbor(&polymer->monomers[monomer_idx_2], polymer->monomers[monomer_idx_1].id);
    }
  }
}

/* Function to fill remaining bond capacity up to MAX_BONDS */
void fillRemainingBonds(MicrogelPolymer* polymer)
{
  printf("Filling remaining bond capacity...\n");
  int bonds_added = 0;

  /* Process monomers in order of increasing current bonds */
  for (int num_bonds = MIN_BONDS; num_bonds < MAX_BONDS; num_bonds++) {
    for (int i = 0; i < polymer->total_nmon; i++) {
      /* Only process monomers with exactly num_bonds bonds */
      if (polymer->monomers[i].num_neighbors == num_bonds) {
        /* Find closest monomers that have room for more bonds */
        DistanceIndexPair* distances = (DistanceIndexPair*)
          malloc((polymer->total_nmon - 1) * sizeof(DistanceIndexPair));

        if (distances == NULL) {
          printf("Error: Memory allocation failed for distances.\n");
          exit(1);
        }

        int dist_idx = 0;

        /* Calculate distances to all other monomers */
        for (int j = 0; j < polymer->total_nmon; j++) {
          if (i != j) {
            double dist = euclideanDistance(&polymer->monomers[i], &polymer->monomers[j]);
            /* Only include monomers within distance threshold */
            if (dist <= polymer->target_distance * (1.0 + polymer->distance_threshold)) {
              distances[dist_idx].distance = dist;
              distances[dist_idx].index = j;
              dist_idx++;
            }
          }
        }

        /* Sort by distance, but only consider the valid ones we've added */
        if (dist_idx > 0) {
          qsort(distances, dist_idx, sizeof(DistanceIndexPair), compareDistancePairs);

          /* Add bonds until maximum is reached */
          for (int k = 0; k < dist_idx && polymer->monomers[i].num_neighbors < MAX_BONDS; k++) {
            int j = distances[k].index;

            /* Only add bond if the other monomer has room for more bonds */
            if (polymer->monomers[j].num_neighbors < MAX_BONDS) {
              int before = polymer->monomers[i].num_neighbors + polymer->monomers[j].num_neighbors;
              createBidirectionalBond(polymer, i, j);
              int after = polymer->monomers[i].num_neighbors + polymer->monomers[j].num_neighbors;

              if (after > before) {
                bonds_added++;
              }
            }
          }
        }

        free(distances);
      }
    }
  }

  printf("Added %d additional bonds to fill capacity\n", bonds_added);
}

/* Modified function to ensure minimum number of bonds */
void ensureMinimumBonds(MicrogelPolymer* polymer)
{
  printf("Ensuring minimum number of bonds for all monomers...\n");
  int bonds_added = 0;

  for (int i = 0; i < polymer->total_nmon; i++) {
    if (polymer->monomers[i].num_neighbors < MIN_BONDS) {
      /* Find closest monomers that have room for more bonds */
      DistanceIndexPair* distances = (DistanceIndexPair*)
        malloc((polymer->total_nmon - 1) * sizeof(DistanceIndexPair));

      if (distances == NULL) {
        printf("Error: Memory allocation failed for distances.\n");
        exit(1);
      }

      int dist_idx = 0;

      /* Calculate distances to all other monomers */
      for (int j = 0; j < polymer->total_nmon; j++) {
        if (i != j) {
          double dist = euclideanDistance(&polymer->monomers[i], &polymer->monomers[j]);
          /* Filter by distance threshold */
          if (dist <= polymer->target_distance * (1.0 + polymer->distance_threshold)) {
            distances[dist_idx].distance = dist;
            distances[dist_idx].index = j;
            dist_idx++;
          }
        }
      }

      /* Sort by distance */
      if (dist_idx > 0) {
        qsort(distances, dist_idx, sizeof(DistanceIndexPair), compareDistancePairs);

        /* Add bonds until minimum is reached */
        for (int k = 0; k < dist_idx && polymer->monomers[i].num_neighbors < MIN_BONDS; k++) {
          int j = distances[k].index;

          /* Only add bond if the other monomer has room for more bonds */
          if (polymer->monomers[j].num_neighbors < MAX_BONDS) {
            int before = polymer->monomers[i].num_neighbors;
            createBidirectionalBond(polymer, i, j);
            int after = polymer->monomers[i].num_neighbors;

            if (after > before) {
              bonds_added++;
            }
          }
        }
      }

      free(distances);

      /* If still not enough bonds, warn about it */
      if (polymer->monomers[i].num_neighbors < MIN_BONDS) {
        printf("Warning: Monomer %d has only %d bonds, less than minimum %d\n",
               polymer->monomers[i].id, polymer->monomers[i].num_neighbors, MIN_BONDS);
      }
    }
  }

  printf("Added %d bonds to ensure minimum connectivity\n", bonds_added);
}

/* createFullPolymerNetwork function */
void createFullPolymerNetwork(MicrogelPolymer* polymer)
{
  assert(polymer);

  // 1. Los monomeros superficiales ya han sido generados por generateInitialMonomers
  // (Esta funcion se llamara antes de createFullPolymerNetwork en main)

  // 2. Generar cadenas lineales desde la superficie hacia el interior
  generatePolymerChains(polymer);

  // 3. Asegurar que ningun monomero quede suelto sin pertenecer a una cadena
  //    La logica de generatePolymerChains y createCrosslinks ya deberia asegurar
  //    que cada monomero tenga al menos un enlace de cadena o entrecruzamiento.
  //    Sin embargo, `ensureMinimumBonds` se puede usar para una verificacion final
  //    y añadir enlaces si algun monomero todavia tiene 0.
  //    Asegurate de que ensureMinimumBonds busque enlaces con CUALQUIER monomero cercano.
  printf("Verificando enlaces minimos y llenando capacidad restante...\n");
  ensureMinimumBonds(polymer); // Ajusta esta funcion para que trabaje con todos los monomeros
  // fillRemainingBonds(polymer); // Esta funcion llena la capacidad restante de enlaces

  // 4. Crear enlaces cruzados
  createCrosslinks(polymer);

  /* Print network statistics */
  printNetworkStats(polymer);

  return;
}

/* Export full network data to files */
void exportFullNetworkData(const MicrogelPolymer* polymer, const char* filename)
{
  /* Export the network data using the existing function */
  exportConnectivityData(polymer, filename);

  /* Also create CSV file for easier visualization */
  char csv_filename[256];
  sprintf(csv_filename, "%s.csv", filename);
  FILE* file = fopen(csv_filename, "w");

  if (file != NULL) {
    fprintf(file, "ID,X,Y,Z,Surface,NumBonds\n");
    for (int i = 0; i < polymer->total_nmon; i++) {
      fprintf(file, "%d,%.6f,%.6f,%.6f,%d,%d\n",
      polymer->monomers[i].id,
      polymer->monomers[i].x,
      polymer->monomers[i].y,
      polymer->monomers[i].z,
      polymer->monomers[i].is_surface,
      polymer->monomers[i].num_neighbors);
    }
    fclose(file);
    printf("CSV data exported to %s\n", csv_filename);
  }

  /* Create bond pairs file for visualization */
  char bonds_filename[256];
  sprintf(bonds_filename, "%s.bonds", filename);
  file = fopen(bonds_filename, "w");

  if (file != NULL) {
    fprintf(file, "# Bond pairs (ID1 ID2)\n");
    for (int i = 0; i < polymer->total_nmon; i++) {
      for (int j = 0; j < polymer->monomers[i].num_neighbors; j++) {
        int bond_id = polymer->monomers[i].neighbors[j];
        /* Only output each bond once (when ID1 < ID2) */
        if (polymer->monomers[i].id < bond_id) {
          fprintf(file, "%d %d\n", polymer->monomers[i].id, bond_id);
        }
      }
    }
    fclose(file);
    printf("Bond pairs exported to %s\n", bonds_filename);
  }
}

/* Print statistics about the network */
void printNetworkStats(const MicrogelPolymer* polymer)
{
  int min_bonds = MAX_BONDS + 1;
  int max_bonds = 0;
  double avg_bonds = 0.0;
  int no_bonds = 0;

  for (int i = 0; i < polymer->total_nmon; i++) {
    int num = polymer->monomers[i].num_neighbors;

    if (num < min_bonds) min_bonds = num;
    if (num > max_bonds) max_bonds = num;
    avg_bonds += num;

    if (num == 0) no_bonds++;
  }

  avg_bonds /= polymer->total_nmon;

  printf("\nNetwork statistics:\n");
  printf("  Total monomers: %d\n", polymer->total_nmon);
  printf("  Surface monomers: %d\n", polymer->nmon);
  printf("  Interior monomers: %d\n", polymer->interior_nmon);
  printf("  Sphere radius: %.4f\n", polymer->radius);
  printf("  Bond length: %.4f (tolerance: %.4f)\n",
  polymer->target_distance, polymer->distance_threshold);
  printf("  Min bonds per monomer: %d\n", min_bonds);
  printf("  Max bonds per monomer: %d\n", max_bonds);
  printf("  Avg bonds per monomer: %.2f\n", avg_bonds);
  printf("  Monomers with no bonds: %d\n", no_bonds);

  /* Verify the distribution quality */
  verifyDistribution(polymer);

  // /* Print connectivity information */
  // printf("\nSample of monomer connectivity:\n");
  // int sample_size = polymer->total_nmon < 10 ? polymer->total_nmon : 10;
  // for (int i = 0; i < sample_size; i++) {
  //   printf("Monomer %d connects to: ", polymer->monomers[i].id);
  //   if (polymer->monomers[i].num_neighbors == 0) {
  //     printf("none");
  //   }
  //   else {
  //     for (int j = 0; j < polymer->monomers[i].num_neighbors; j++) {
  //       printf("%d", polymer->monomers[i].neighbors[j]);
  //       if (j < polymer->monomers[i].num_neighbors - 1) {
  //         printf(", ");
  //       }
  //     }
  //   }
  //   printf("\n");
  // }
}

// Genera un numero aleatorio double entre min y max
double randomDouble(double min, double max) {
  return min + ((double)rand() / RAND_MAX) * (max - min);
}

// generatePolymerChains 
void generatePolymerChains(MicrogelPolymer* polymer) {
  printf("Generando cadenas polimericas lineales desde la superficie hacia el interior...\n");

  // Calculamos la cantidad total de monomeros interiores necesaria para cumplir la densidad
  double sphere_volume = (4.0 / 3.0) * M_PI * pow(polymer->radius, 3);
  int target_total_monomers = (int)(polymer->density * sphere_volume);
  printf("target_total_monomers: %i\n", target_total_monomers);
  // Los monomeros superficiales ya estan contados
  int interior_monomers_to_add = target_total_monomers - polymer->nmon;

  if (interior_monomers_to_add < 0) {
    printf("Advertencia: La densidad objetivo implica menos monomeros totales que los superficiales. No se añadiran cadenas interiores.\n");
    interior_monomers_to_add = 0;
  }

  int current_monomer_id = polymer->nmon; // IDs de los monomeros interiores comienzan despues de los superficiales
  int monomers_added_in_chains = 0;
  int max_placement_attempts_per_monomer = 30000; // Intentos para colocar un monomero
  int max_chain_extension_attempts = 5000; // Intentos para extender una cadena si se atasca

  // Iterar a traves de los monomeros superficiales para iniciar cadenas
  for (int i = 0; i < polymer->nmon && monomers_added_in_chains < interior_monomers_to_add; i++) {
    Point3D* current_chain_end = &polymer->monomers[i]; // El monomero superficial es el inicio de la cadena

    // Calcular una longitud de cadena aleatoria alrededor de la media
    // Usamos una distribucion uniforme simple para variabilidad
    int chain_length = (int)round(randomDouble(polymer->avg_chain_length * 0.5, polymer->avg_chain_length * 1.5));
    if (chain_length < 1) chain_length = 1; // Longitud minima de cadena

    int chain_extension_attempts_count = 0;

    for (int k = 0; k < chain_length && monomers_added_in_chains < interior_monomers_to_add; k++) {
      double new_x, new_y, new_z;
      int placement_attempts = 0;
      int placed = 0;

      while (placement_attempts < max_placement_attempts_per_monomer && !placed) {
        // Generar una direccion aleatoria con tendencia hacia el centro
        // Esta es una heuristica; podria necesitar mas refinamiento para formas complejas o densidades muy altas.
        // double dir_x = -current_chain_end->x + randomDouble(-1.0, 1.0);
        // double dir_y = -current_chain_end->y + randomDouble(-1.0, 1.0);
        // double dir_z = -current_chain_end->z + randomDouble(-1.0, 1.0);
        double dir_x = randomDouble(-1.0, 1.0);
        double dir_y = randomDouble(-1.0, 1.0);
        double dir_z = randomDouble(-1.0, 1.0);
        double len_dir = sqrt(dir_x * dir_x + dir_y * dir_y + dir_z * dir_z);
        if (len_dir == 0) { // Evitar division por cero, generar nueva direccion si es nula
          dir_x = randomDouble(-1.0, 1.0);
          dir_y = randomDouble(-1.0, 1.0);
          dir_z = randomDouble(-1.0, 1.0);
          len_dir = sqrt(dir_x * dir_x + dir_y * dir_y + dir_z * dir_z);
        }

        // Calcular la nueva posicion a 'target_distance' del monomero anterior
        // new_x = current_chain_end->x + (dir_x / len_dir) * polymer->target_distance;
        // new_y = current_chain_end->y + (dir_y / len_dir) * polymer->target_distance;
        // new_z = current_chain_end->z + (dir_z / len_dir) * polymer->target_distance;
        double target_dist = randomDouble(polymer->min_distance, polymer->target_distance *(1+polymer->distance_threshold));
        new_x = current_chain_end->x + (dir_x / len_dir) * target_dist;
        new_y = current_chain_end->y + (dir_y / len_dir) * target_dist;
        new_z = current_chain_end->z + (dir_z / len_dir) * target_dist;

        // Verificar que la nueva posicion este dentro de la esfera y no se superponga
        double dist_from_center = sqrt(new_x * new_x + new_y * new_y + new_z * new_z);

        if (dist_from_center < polymer->radius &&
            isValidPosition(polymer, new_x, new_y, new_z, current_monomer_id)) { // current_monomer_id es el ID del nuevo monomero
          initPoint3D(&polymer->monomers[current_monomer_id], new_x, new_y, new_z, current_monomer_id + 1);
          polymer->monomers[current_monomer_id].is_surface = 0; // Es un monomero interior

          // Crear el enlace entre el monomero actual y el anterior en la cadena
          createBidirectionalBond(polymer, current_chain_end->id - 1, current_monomer_id);

          current_chain_end = &polymer->monomers[current_monomer_id];
          current_monomer_id++;
          monomers_added_in_chains++;
          placed = 1;
          chain_extension_attempts_count = 0; // Resetear intentos al extender con exito
        }
        placement_attempts++;
      } // Fin while placement_attempts

      if (!placed) {
        // Si no se pudo colocar el monomero, intentar con una nueva direccion para la cadena
        chain_extension_attempts_count++;
        if (chain_extension_attempts_count >= max_chain_extension_attempts) {
          // printf("Info: No se pudo extender la cadena desde el monomero %d despues de %d intentos. Terminando cadena.\n", current_chain_end->id, max_chain_extension_attempts);
          break; // Terminar la cadena si se atasca demasiado
        }
        // Si no se pudo colocar, el bucle k sigue, y current_chain_end no cambia,
        // por lo que intentara extender desde el mismo punto con otra direccion aleatoria.
      }
    } // Fin for k (longitud de cadena)
  } // Fin for i (monomeros superficiales)

  polymer->interior_nmon = monomers_added_in_chains;
  polymer->total_nmon = polymer->nmon + polymer->interior_nmon;
  printf("Generados %d monomeros interiores en cadenas. Total de monomeros: %d.\n", polymer->interior_nmon, polymer->total_nmon);

  // Si aún quedan monómeros por añadir para la densidad objetivo, se podrían
  // implementar estrategias adicionales aquí (ej. iniciar nuevas cadenas
  // desde el interior si hay espacio, o ajustar el radio).
  // Por simplicidad, este ejemplo se basa en las cadenas que se forman desde la superficie.
  if (polymer->total_nmon < target_total_monomers) {
    printf("Advertencia: No se alcanzo la densidad objetivo. Faltan %d monomeros para un total de %d.\n",
           target_total_monomers - polymer->total_nmon, target_total_monomers);
  }
  else if (polymer->total_nmon > target_total_monomers) {
    printf("Advertencia: Se supero la densidad objetivo. Se generaron %d monomeros, se esperaban %d.\n",
           polymer->total_nmon, target_total_monomers);
  }
}

// createCrosslinks 
void createCrosslinks(MicrogelPolymer* polymer) {
  printf("Creando enlaces cruzados (cross-links)...\n");
  int crosslinks_added = 0;
  int max_crosslink_attempts_per_monomer = 10; // Intentos para encontrar un compañero de entrecruzamiento

  // Iterar a traves de todos los monomeros (excluyendo los superficiales si se desea,
  // pero para entrecruzamientos, cualquier monomero interior puede ser un objetivo).
  // Comenzamos desde los monomeros superficiales (indice 0) para considerar todos los monomeros.
  for (int i = 0; i < polymer->total_nmon; i++) {
    // Probabilidad de entrecruzamiento para este monomero
    if (((double)rand() / RAND_MAX) < polymer->crosslink_density) {
      int attempts = 0;
      int crosslink_found = 0;
      while (attempts < max_crosslink_attempts_per_monomer && !crosslink_found) {
        // Seleccionar un monomero aleatorio como posible compañero de entrecruzamiento
        int target_j = rand() % polymer->total_nmon;

        // Evitar auto-enlaces y enlaces con monomeros ya unidos, y asegurar capacidad de enlaces
        if (i != target_j && !areBonded(polymer, i, target_j) &&
            polymer->monomers[i].num_neighbors < MAX_BONDS &&
            polymer->monomers[target_j].num_neighbors < MAX_BONDS) {

          double dist = euclideanDistance(&polymer->monomers[i], &polymer->monomers[target_j]);
          // Verificar si la distancia es adecuada para un enlace
          if (dist <= polymer->target_distance * (1.0 + polymer->distance_threshold) &&
              dist >= polymer->target_distance * (1.0 - polymer->distance_threshold)) {
            createBidirectionalBond(polymer, i, target_j);
            crosslinks_added++;
            crosslink_found = 1;
          }
        }
        attempts++;
      }
    }
  }
  printf("Anadidos %d enlaces cruzados.\n", crosslinks_added);
}

// --- Funcion auxiliar para verificar si dos monomeros ya estan unidos ---
int areBonded(MicrogelPolymer* polymer, int idx1, int idx2) {
  for (int i = 0; i < polymer->monomers[idx1].num_neighbors; i++) {
    if (polymer->monomers[idx1].neighbors[i] == polymer->monomers[idx2].id) {
      return 1;
    }
  }
  return 0;
}

NearbyMonomerList* findNearbyMonomers(const MicrogelPolymer* polymer, double factor) {
  int N = polymer->total_nmon;
  double max_dist = polymer->target_distance * factor;

  /*CHANGE INIT - 20251119 CUDA C++ compatibility fix */
  /* Original: NearbyMonomerList* result = malloc(N * sizeof(NearbyMonomerList)); */
  /* C++ requires explicit cast from void* */
  NearbyMonomerList* result = (NearbyMonomerList*)malloc(N * sizeof(NearbyMonomerList));
  /*CHANGE END*/
  if (!result) {
    fprintf(stderr, "Error: No se pudo asignar memoria para NearbyMonomerList\n");
    exit(1);
  }

  for (int i = 0; i < N; i++) {
    result[i].monomer_index = i;
    result[i].num_neighbors = 0;

    // Primero contar cuántos vecinos hay dentro del rango
    for (int j = 0; j < N; j++) {
      if (i == j) continue;
      double dist = euclideanDistance(&polymer->monomers[i], &polymer->monomers[j]);
      if (dist <= max_dist) {
        result[i].num_neighbors++;
      }
    }

    // Asignar memoria para los vecinos
    /*CHANGE INIT - 20251119 CUDA C++ compatibility fix */
    /* Original: result[i].neighbors = malloc(result[i].num_neighbors * sizeof(int)); */
    /* C++ requires explicit cast from void* */
    result[i].neighbors = (int*)malloc(result[i].num_neighbors * sizeof(int));
    /*CHANGE END*/
    if (!result[i].neighbors) {
      fprintf(stderr, "Error: No se pudo asignar memoria para vecinos de monomero %d\n", i);
      exit(1);
    }

    // Volver a recorrer y guardar los índices
    int idx = 0;
    for (int j = 0; j < N; j++) {
      if (i == j) continue;
      double dist = euclideanDistance(&polymer->monomers[i], &polymer->monomers[j]);
      if (dist <= max_dist) {
        result[i].neighbors[idx++] = j;
      }
    }
  }

  return result;
}