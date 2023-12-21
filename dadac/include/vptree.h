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


struct vpTreeNode
{
   void * data;
   idx_t array_idx;
   float_t mu;
   float_t __dist;
   struct vpTreeNode* parent;
   struct vpTreeNode* outside;
   struct vpTreeNode* inside;

};

#ifdef ITERATIVE_VPTREE

struct stackNode
{
	int side;
	struct vpTreeNode* node;
	float_t current_distance;
	float_t mu;
};

struct stack_vpTreeNodes
{
	struct stackNode* data;
	size_t count;
	size_t size;
};

typedef struct stackNode stackNode;
typedef struct stack_vpTreeNodes stack_vpTreeNodes;

#define DEFAULT_STACK_SIZE 200

#endif



typedef struct vpTreeNode vpTreeNode;


void initialize_vpTreeNodes_pointers(vpTreeNode** pointersArray, vpTreeNode* nodeArray, idx_t n);
void initialize_vpTreeNode_array(vpTreeNode* nodeArray, void* data, idx_t n, idx_t bytesPerElement);
vpTreeNode* build_vpTree(vpTreeNode** t, int start, int end, vpTreeNode* parent, float_t (*metric)(void*, void*));
void KNN_sub_vpTree_search(void* point, vpTreeNode* root, Heap * H, float_t (*metric)(void*,void*));

#ifdef ITERATIVE_VPTREE
	Heap KNN_vpTree(void* point, vpTreeNode* root, int maxk, stack_vpTreeNodes* s, float_t (*metric)(void*, void*));
	void stackInit(stack_vpTreeNodes* s);
#else
	Heap KNN_vpTree(void* point, vpTreeNode* root, int maxk, float_t (*metric)(void*, void*));
#endif 

