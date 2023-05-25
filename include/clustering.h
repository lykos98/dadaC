#include "kdtree.h"
#include<stdlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <stdint.h>

#define DTHR 23.92812698
#define PI_F 3.1415926f
#define ARRAY_INCREMENT 500
#define DA_DTYPE size_t
#define NOBORDER SIZE_MAX 

/**********************************
 * DATA STRUCTURES FOR CLUSTERING *
 **********************************/
struct Node {
    size_t data;
    struct Node* next;
    
};

struct LinkedList {
    size_t count;
    struct Node* head; 
};


struct lu_dynamicArray {
    size_t * data;
    size_t size;
    size_t count;
};

struct Datapoint_info {
    FLOAT_TYPE g;
    Heap ngbh;
    size_t array_idx;
    FLOAT_TYPE log_rho;
    FLOAT_TYPE log_rho_c;
    FLOAT_TYPE log_rho_err;
    size_t kstar;
    int is_center;
    int cluster_idx;
};

struct border_t
{
    FLOAT_TYPE density;
    FLOAT_TYPE error;
    size_t idx;        
};

struct aClusters {
    struct lu_dynamicArray centers;
    FLOAT_TYPE** border_density;
    FLOAT_TYPE** border_err;
    size_t** border_idx;
    FLOAT_TYPE* __border_density_data;
    FLOAT_TYPE* __border_err_data;
    size_t* __border_idx_data;
    size_t n;
    
};

struct Clusters {
    struct lu_dynamicArray centers;
    struct border_t ** borders;
    struct border_t * __borders_data; 
    size_t n;
};


struct merge_t {
  size_t source;
  size_t target;
  FLOAT_TYPE density;
};


typedef struct Datapoint_info Datapoint_info;
typedef struct lu_dynamicArray lu_dynamicArray;
typedef struct Clusters Clusters;
typedef struct Node Node;
typedef struct LinkedList LinkedList;
typedef struct merge_t merge_t;
typedef struct border_t border_t;

void LinkedList_Insert(LinkedList* L, Node* n);
void DynamicArray_allocate(lu_dynamicArray * a);
void DynamicArray_pushBack(lu_dynamicArray * a, size_t p);
void Clusters_allocate(Clusters * c);
void Clusters_free(Clusters * c);

int cmp(const void * a, const void * b);
FLOAT_TYPE avg(const FLOAT_TYPE * x, const size_t n);
FLOAT_TYPE mEst2(FLOAT_TYPE * x, FLOAT_TYPE *y, size_t n);
FLOAT_TYPE mEst(FLOAT_TYPE * x, FLOAT_TYPE *y, size_t n);
FLOAT_TYPE idEstimate(Datapoint_info* particles, size_t n);
void computeRho(Datapoint_info* particles, const FLOAT_TYPE d, const size_t points);
void computeRhoOpt(Datapoint_info* particles, kd_node* root, FLOAT_TYPE* data, const FLOAT_TYPE d, const size_t points);
int cmpPP(const void* p1, const void *p2);
void computeCorrection(Datapoint_info* particles, size_t n, FLOAT_TYPE Z);
void KNN_search(Datapoint_info * particles, FLOAT_TYPE * data, kd_node* root, size_t n, size_t k);


Clusters Heuristic1(Datapoint_info* particles, FLOAT_TYPE* data, size_t n);
void Heuristic2(Clusters* cluster, Datapoint_info* particles);
void Heuristic3(Clusters* cluster, Datapoint_info* particles, FLOAT_TYPE Z, int halo);
