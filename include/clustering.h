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


struct Clusters {
    struct lu_dynamicArray centers;
    struct LinkedList* clusters;   
    struct Node* _LLnodes;
    FLOAT_TYPE** border_density;
    FLOAT_TYPE** border_err;
    size_t** border_idx;
    FLOAT_TYPE* __border_density_data;
    FLOAT_TYPE* __border_err_data;
    size_t* __border_idx_data;
    
};

typedef struct Datapoint_info Datapoint_info;
typedef struct lu_dynamicArray lu_dynamicArray;
typedef struct Clusters Clusters;
typedef struct Node Node;
typedef struct LinkedList LinkedList;

void LinkedList_Insert(LinkedList* L, Node* n);
void DynamicArray_allocate(lu_dynamicArray * a);
void DynamicArray_pushBack(lu_dynamicArray * a, size_t p);
void Clusters_allocate(Clusters * c);
void Clusters_free(Clusters * c);

void mergeClusters(Clusters * cc, size_t i, size_t j);
int cmp(const void * a, const void * b);
FLOAT_TYPE avg(const FLOAT_TYPE * x, const size_t n);
FLOAT_TYPE mEst2(FLOAT_TYPE * x, FLOAT_TYPE *y, size_t n);
FLOAT_TYPE mEst(FLOAT_TYPE * x, FLOAT_TYPE *y, size_t n);
FLOAT_TYPE idEstimate(Datapoint_info* particles, size_t n);
void computeRho(Datapoint_info* particles, const FLOAT_TYPE d, const size_t points);
int cmpPP(const void* p1, const void *p2);
void calculateCorrection(Datapoint_info* particles, size_t n, FLOAT_TYPE Z);
void KNN_search(Datapoint_info * particles, FLOAT_TYPE * data, kd_node* root, size_t n, size_t k);


Clusters Heuristic1(Datapoint_info* particles, FLOAT_TYPE* data, size_t n);
void Heuristic2(Clusters* cluster, Datapoint_info* particles);
void Heuristic3(Clusters* cluster, Datapoint_info* particles, FLOAT_TYPE Z, int halo);
