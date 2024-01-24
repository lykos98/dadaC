#pragma once
#include "heap.h"
#include "kdtreeV2.h"
#include "kdtree.h"
#include "vptree.h"
#include "vptreeV2.h"
#include <math.h>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <string.h>

#define DTHR 23.92812698
#define PI_F 3.1415926f
#define ARRAY_INCREMENT 500
#define DA_DTYPE idx_t
#define NOBORDER MY_SIZE_MAX

#define VERBOSE_TRUE 1
#define VERBOSE_FALSE 0

/**********************************
 * DATA STRUCTURES FOR CLUSTERING *
 **********************************/

struct lu_dynamic_array {
  idx_t *data;
  idx_t size;
  idx_t count;
};

struct datapoint_info {
  float_t g;
  heap ngbh;
  idx_t array_idx;
  float_t log_rho;
  float_t log_rho_c;
  float_t log_rho_err;
  idx_t kstar;
  int is_center;
  int cluster_idx;
};

struct border_t {
  float_t density;
  float_t error;
  idx_t idx;
};

struct sparse_border_t {
  idx_t i;
  idx_t j;
  idx_t idx;
  float_t density;
  float_t error;
};

struct adj_list_t {
  idx_t count;
  idx_t size;
  struct sparse_border_t* data;
};


struct clusters {
  int use_sparse_borders;
  struct adj_list_t *sparse_borders;
  struct lu_dynamic_array centers;
  struct border_t **borders;
  struct border_t *__borders_data;
  idx_t n;
};

struct merge_t {
  idx_t source;
  idx_t target;
  float_t density;
};

typedef struct datapoint_info datapoint_info;
typedef struct lu_dynamic_array lu_dynamic_array;
typedef struct clusters clusters;
typedef struct merge_t merge_t;
typedef struct border_t border_t;
typedef struct sparse_border_t sparse_border_t;
typedef struct adj_list_t adj_list_t; 

void lu_dynamic_array_allocate(lu_dynamic_array *a);
void lu_dynamic_array_pushBack(lu_dynamic_array *a, idx_t p);
//void Clusters_allocate(clusters *c);
void clusters_allocate(clusters *c, int s);
void clusters_free(clusters *c);

int cmp(const void *a, const void *b);
float_t avg(const float_t *x, const idx_t n);
float_t mEst2(float_t *x, float_t *y, idx_t n);
float_t mEst(float_t *x, float_t *y, idx_t n);
float_t id_estimate(datapoint_info *particles, idx_t n, float_t fraction,int verbose);

datapoint_info* ngbh_search_vptree(void* data, size_t n, size_t byteSize, size_t dims, size_t k, float_t (*metric)(void *, void *), int verbose);
datapoint_info* ngbh_search_vptree_v2(void* data, size_t n, size_t byteSize, size_t dims, size_t k, float_t (*metric)(void *, void *), int verbose);
datapoint_info* ngbh_search_kdtree(float_t* data, size_t n, size_t ndims, size_t k, int verbose);
datapoint_info* ngbh_search_kdtree_v2(float_t* data, size_t n, size_t ndims, size_t k, int verbose);
datapoint_info* ngbh_search_bruteforce(void* data, size_t n, size_t byteSize, size_t dims, size_t k, float_t (*metric)(void *, void *), int verbose);
void compute_density_kstarnn(datapoint_info *particles, const float_t d,
                const idx_t points, int verbose);

void PAk(datapoint_info* dpInfo, const float_t d, const idx_t points, int verbose);
int cmpPP(const void *p1, const void *p2);
void compute_correction(datapoint_info *particles, idx_t n, float_t Z);

clusters Heuristic1(datapoint_info *dp, idx_t n, int verbose);
void Heuristic2(clusters *cluster, datapoint_info *particles, int verbose);
void Heuristic3(clusters* cluster, datapoint_info* dp_info, float_t Z, int halo, int verbose);
void free_datapoint_array(datapoint_info* d, size_t n);

float_t* eucMetricPx3(void*, void*);
float_t eud(void*, void*);
float_t eud_sq(void*, void*);
float_t eudOpt(void*, void*);
