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

struct vpTreeNodeV2List
{
	size_t count;
	#ifdef VOPT
		idx_t* indexes;
		void* start_ptr;
		void* end_ptr;
	#else
		struct vpTreeNodeV2** data;
	#endif
};


struct vpTreeNodeV2
{
   void * data;
   struct vpTreeNodeV2List nodeList;
   idx_t array_idx;
   struct vpTreeNodeV2* outside;
   struct vpTreeNodeV2* inside;
   struct vpTreeNodeV2* parent;
   float_t mu;
   float_t __dist;
   int isLeaf;
   idx_t __bytesize;
};


typedef struct vpTreeNodeV2 vpTreeNodeV2;
typedef struct vpTreeNodeV2List vpTreeNodeV2List;


void initialize_vpTreeNodes_pointers_V2(vpTreeNodeV2** pointersArray, vpTreeNodeV2* nodeArray, idx_t n);
void initialize_vpTreeNode_array_V2(vpTreeNodeV2* noV2deArray, void* data, idx_t n, idx_t bytesPerElement);
vpTreeNodeV2* build_vpTree_V2(vpTreeNodeV2* t, int start, int end, vpTreeNodeV2* parent, float_t (*metric)(void*, void*));
void KNN_sub_vpTree_search_V2(void* point, vpTreeNodeV2* root, Heap * H, float_t (*metric)(void*,void*));
Heap KNN_vpTree_V2(void* point, vpTreeNodeV2* root, int maxk, float_t (*metric)(void*, void*));

