#pragma once
#include "heap.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <stdint.h>

#ifdef USE_FLOAT32
	#define float_t float
#else
	#define float_t double 
#endif

#ifdef USE_INT32
	#define MY_SIZE_MAX UINT32_MAX
	#define idx_t uint32_t
#else
	#define MY_SIZE_MAX UINT64_MAX
	#define idx_t uint64_t
#endif

struct vpnode_v2_list
{
	size_t count;
	#ifdef VOPT
		idx_t* indexes;
		void* start_ptr;
		void* end_ptr;
	#else
		struct vpnode_v2** data;
	#endif
};


struct vpnode_v2
{
   void * data;
   struct vpnode_v2_list node_list;
   idx_t array_idx;
   struct vpnode_v2* outside;
   struct vpnode_v2* inside;
   struct vpnode_v2* parent;
   float_t mu;
   float_t __dist;
   int is_leaf;
   idx_t __bytesize;
};


typedef struct vpnode_v2 vpnode_v2;
typedef struct vpnode_v2_list vpnode_v2_list;


void initialize_vpnode_v2_ptrs(vpnode_v2** pointersArray, vpnode_v2* nodeArray, idx_t n);
void initialize_vpnode_v2_array(vpnode_v2* noV2deArray, void* data, idx_t n, idx_t bytesPerElement);
vpnode_v2* build_vptree_v2(vpnode_v2* t, int start, int end, vpnode_v2* parent, float_t (*metric)(void*, void*));
void knn_sub_vptree_v2_search(void* point, vpnode_v2* root, heap * H, float_t (*metric)(void*,void*));
heap knn_vptree_v2(void* point, vpnode_v2* root, int maxk, float_t (*metric)(void*, void*));

