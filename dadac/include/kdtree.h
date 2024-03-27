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



struct kdnode
{   
   int level;
   int split_var;
   FLOAT_TYPE * data;
   idx_t array_idx;
   struct kdnode* parent;
   struct kdnode* lch;
   struct kdnode* rch;
};

typedef struct kdnode kdnode;
typedef struct heap_node heap_node;

void swap(T* a, T* b);

FLOAT_TYPE euclidean_distance(FLOAT_TYPE* p1, FLOAT_TYPE* p2);

void swap_heap_node(heap_node* a, heap_node* b);

void swap_kdnode_ptrs(kdnode **x, kdnode **y);


/**
 * 
 * KDtree implementation 
 * 
 * 
*/

void initialize_kdnodes(kdnode * node_array, FLOAT_TYPE* d, idx_t n );

void initialize_kdnode_ptrs(kdnode** node_ptr_array, kdnode* node_array, idx_t n );

int cmp_kdnodes(kdnode* a, kdnode* b, int var);

void print_kdnode(kdnode* node);

// Standard Lomuto partition function

int partition(kdnode** arr, int low, int high, int split_var);

// Implementation of QuickSelect
int median_of_nodes(kdnode** a, int left, int right, int split_var);

kdnode* make_tree(kdnode** t, int start, int end, kdnode* parent, int level);

//inline FLOAT_TYPE hyper_plane_dist(FLOAT_TYPE* p1, FLOAT_TYPE* p2, int var);

//inline int hyper_plane_side(FLOAT_TYPE* p1, FLOAT_TYPE* p2, int var);

void knn_sub_tree_search(FLOAT_TYPE* point, kdnode* kdtree_root, heap * H);


heap knn(FLOAT_TYPE* point, kdnode* kdtree_root, int maxk);

kdnode * build_tree(kdnode** kd_ptrs, size_t n, size_t dimensions);
