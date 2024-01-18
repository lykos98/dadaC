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

/**********************************
 * DATA STRUCTURES FOR CLUSTERING *
 **********************************/
struct Node {
  idx_t data;
  struct Node *next;
};

struct LinkedList {
  idx_t count;
  struct Node *head;
};

struct lu_dynamicArray {
  idx_t *data;
  idx_t size;
  idx_t count;
};

struct Datapoint_info {
  float_t g;
  Heap ngbh;
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

struct SparseBorder_t {
  idx_t i;
  idx_t j;
  idx_t idx;
  float_t density;
  float_t error;
};

struct AdjList_t {
  idx_t count;
  idx_t size;
  struct SparseBorder_t* data;
};

struct aClusters {
  struct lu_dynamicArray centers;
  float_t **border_density;
  float_t **border_err;
  idx_t **border_idx;
  float_t *__border_density_data;
  float_t *__border_err_data;
  idx_t *__border_idx_data;
  idx_t n;
};

struct Clusters {
  int UseSparseBorders;
  struct AdjList_t *SparseBorders;
  struct lu_dynamicArray centers;
  struct border_t **borders;
  struct border_t *__borders_data;
  idx_t n;
};

struct merge_t {
  idx_t source;
  idx_t target;
  float_t density;
};

typedef struct Datapoint_info Datapoint_info;
typedef struct lu_dynamicArray lu_dynamicArray;
typedef struct Clusters Clusters;
typedef struct Node Node;
typedef struct LinkedList LinkedList;
typedef struct merge_t merge_t;
typedef struct border_t border_t;
typedef struct SparseBorder_t SparseBorder_t;
typedef struct AdjList_t AdjList_t; 

void LinkedList_Insert(LinkedList *L, Node *n);
void DynamicArray_allocate(lu_dynamicArray *a);
void DynamicArray_pushBack(lu_dynamicArray *a, idx_t p);
//void Clusters_allocate(Clusters *c);
void Clusters_allocate(Clusters *c, int s);
void Clusters_free(Clusters *c);

int cmp(const void *a, const void *b);
float_t avg(const float_t *x, const idx_t n);
float_t mEst2(float_t *x, float_t *y, idx_t n);
float_t mEst(float_t *x, float_t *y, idx_t n);
float_t idEstimate(Datapoint_info *particles, idx_t n, float_t fraction);

Datapoint_info* NgbhSearch_vptree(void* data, size_t n, size_t byteSize, size_t dims, size_t k, float_t (*metric)(void *, void *));
Datapoint_info* NgbhSearch_vptree_V2(void* data, size_t n, size_t byteSize, size_t dims, size_t k, float_t (*metric)(void *, void *));
Datapoint_info* NgbhSearch_kdtree(float_t* data, size_t n, size_t ndims, size_t k);
Datapoint_info* NgbhSearch_kdtree_V2(float_t* data, size_t n, size_t ndims, size_t k);
Datapoint_info* NgbhSearch_bruteforce(void* data, size_t n, size_t byteSize, size_t dims, size_t k, float_t (*metric)(void *, void *));
void computeRho(Datapoint_info *particles, const float_t d,
                const idx_t points);

void PAk(Datapoint_info* dpInfo, const float_t d, const idx_t points);
int cmpPP(const void *p1, const void *p2);
void computeCorrection(Datapoint_info *particles, idx_t n, float_t Z);
void KNN_search_kdtree(Datapoint_info *particles, float_t *data, kd_node *root,
                idx_t n, idx_t k);

Clusters Heuristic1(Datapoint_info *dp, idx_t n);
void Heuristic2(Clusters *cluster, Datapoint_info *particles);
void Heuristic3(Clusters *cluster, Datapoint_info *particles, float_t Z,int halo);
void freeDatapointArray(Datapoint_info* d, size_t n);

float_t* eucMetricPx3(void*, void*);
float_t eud(void*, void*);
float_t eud_sq(void*, void*);
float_t eudOpt(void*, void*);
