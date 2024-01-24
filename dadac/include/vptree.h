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


struct vpnode
{
   void * data;
   idx_t array_idx;
   float_t mu;
   float_t __dist;
   struct vpnode* parent;
   struct vpnode* outside;
   struct vpnode* inside;

};

#ifdef ITERATIVE_VPTREE

struct stackNode
{
	int side;
	struct vpnode* node;
	float_t current_distance;
	float_t mu;
};

struct stack_vpnode
{
	struct stackNode* data;
	size_t count;
	size_t size;
};

typedef struct stackNode stackNode;
typedef struct stack_vpnode stack_vpnode;

#define DEFAULT_STACK_SIZE 200

#endif



typedef struct vpnode vpnode;


void initialize_vpnode_ptrs(vpnode** pointersArray, vpnode* nodeArray, idx_t n);
void initialize_vpnode_array(vpnode* nodeArray, void* data, idx_t n, idx_t bytesPerElement);
vpnode* build_vptree(vpnode** t, int start, int end, vpnode* parent, float_t (*metric)(void*, void*));
void knn_sub_vptree_search(void* point, vpnode* root, heap * H, float_t (*metric)(void*,void*));

#ifdef ITERATIVE_VPTREE
	heap knn_vptree(void* point, vpnode* root, int maxk, stack_vpnode* s, float_t (*metric)(void*, void*));
	void stackInit(stack_vpnode* s);
#else
	heap knn_vptree(void* point, vpnode* root, int maxk, float_t (*metric)(void*, void*));
#endif 

