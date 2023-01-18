#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#define T float
#define DATA_DIMS 0 
#define FLOAT_TYPE double 

struct heap_node
{
   FLOAT_TYPE value;
   size_t array_idx;
} ;

struct Heap
{
   size_t N; 
   size_t count;
   struct heap_node* data;
   
} ;

struct SimpleHeap
{
   size_t N; 
   T * data;
} ;

struct kd_node
{   
   int level;
   int split_var;
   FLOAT_TYPE * data;
   size_t array_idx;
   struct kd_node* parent;
   struct kd_node* lch;
   struct kd_node* rch;
};

typedef struct kd_node kd_node;
typedef struct heap_node heap_node;
typedef struct Heap Heap;
typedef struct SimpleHeap SimpleHeap;

void swap(T* a, T* b);

FLOAT_TYPE euclidean_distance(FLOAT_TYPE* p1, FLOAT_TYPE* p2);

void swapHeapNode(heap_node* a, heap_node* b);

void swap_kd_node_ptrs(kd_node **x, kd_node **y);

void allocateSimpleHeap(SimpleHeap* H, size_t n);

void allocateHeap(Heap* H, size_t n);

void initSimpleHeap(SimpleHeap* H);

void initHeap(Heap* H);

void freeSimpleHeap(SimpleHeap * H);

void freeHeap(Heap * H);

void heapifyMaxSimpleHeap(SimpleHeap* H, size_t node);

void heapifyMaxHeap(Heap* H, size_t node);


void setRootMaxSimpleHeap(SimpleHeap * H, T val);

void setRootMaxHeap(Heap * H, FLOAT_TYPE val, size_t array_idx);

void insertMaxHeap(Heap * H, FLOAT_TYPE val, size_t array_idx);

/**
 * 
 * KDtree implementation 
 * 
 * 
*/

void initializeKDnodes(kd_node * node_array, FLOAT_TYPE* d, size_t n );

void initializePTRS(kd_node** node_ptr_array, kd_node* node_array, size_t n );

int cmpKDnodes(kd_node* a, kd_node* b, int var);

void printKDnode(kd_node* node);

// Standard Lomuto partition function

int partition(kd_node** arr, int low, int high, int split_var);

// Implementation of QuickSelect
int medianOfNodes(kd_node** a, int left, int right, int split_var);

kd_node* make_tree(kd_node** t, int start, int end, kd_node* parent, int level);

inline FLOAT_TYPE hyper_plane_dist(FLOAT_TYPE* p1, FLOAT_TYPE* p2, int var);

inline int hyper_plane_side(FLOAT_TYPE* p1, FLOAT_TYPE* p2, int var);

void KNN_sub_tree_search(FLOAT_TYPE* point, kd_node* kdtree_root, Heap * H);

void HeapSort(Heap* H);

Heap KNN(FLOAT_TYPE* point, kd_node* kdtree_root, int maxk);