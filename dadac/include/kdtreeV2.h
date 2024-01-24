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

struct kdnode_v2_list
{
	size_t count;
	struct kdnode_v2** data;
};


struct kdnode_v2
{   
   FLOAT_TYPE * data;
   int level;
   int split_var;
   int is_leaf;
   idx_t array_idx;
   struct kdnode_v2* parent;
   struct kdnode_v2* lch;
   struct kdnode_v2* rch;
   struct kdnode_v2_list node_list;
};

typedef struct kdnode_v2 kdnode_v2;

void initialize_kdnodes_v2(kdnode_v2* node_array, FLOAT_TYPE* d, idx_t n );

// Standard Lomuto partition function

heap knn_kdtree_v2(FLOAT_TYPE* point, kdnode_v2* kdtree_root, int maxk);

kdnode_v2 * build_tree_kdtree_v2(kdnode_v2* kd_ptrs, size_t n, size_t dimensions);
