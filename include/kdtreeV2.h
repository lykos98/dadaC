#pragma once
#include "heap.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <stdint.h>
#define T double
#define DATA_DIMS 0 

#ifdef USE_FLOAT32
	#define FLOAT_TYPE float
#else
	#define FLOAT_TYPE double 
#endif

#ifdef USE_INT32
	#define MY_SIZE_MAX UINT32_MAX
	#define idx_t uint32_t
#else
	#define MY_SIZE_MAX UINT64_MAX
	#define idx_t uint64_t
#endif

struct kdTreeNodeV2List
{
	size_t count;
	struct kdNodeV2** data;
};


struct kdNodeV2
{   
   FLOAT_TYPE * data;
   int level;
   int split_var;
   int isLeaf;
   idx_t array_idx;
   struct kdNodeV2* parent;
   struct kdNodeV2* lch;
   struct kdNodeV2* rch;
   struct kdTreeNodeV2List nodeList;
};

typedef struct kdNodeV2 kdNodeV2;
typedef struct heap_node heap_node;
typedef struct Heap Heap;
typedef struct SimpleHeap SimpleHeap;

void initializeKDnodesV2(kdNodeV2* node_array, FLOAT_TYPE* d, idx_t n );

// Standard Lomuto partition function

Heap KNN_kdTreeV2(FLOAT_TYPE* point, kdNodeV2* kdtree_root, int maxk);

kdNodeV2 * build_tree_kdTreeV2(kdNodeV2* kd_ptrs, size_t n, size_t dimensions);
