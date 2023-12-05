#pragma once
#include "kdtree.h"
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

#define HEAP_LCH(x) (2*x + 1)
#define HEAP_RCH(x) (2*x + 2)
#define HEAP_PARENT(x) (x-1)/2  

#define HP_LEFT_SIDE 0
#define HP_RIGHT_SIDE 1


struct heap_node
{
   FLOAT_TYPE value;
   idx_t array_idx;
} ;

struct Heap
{
   idx_t N; 
   idx_t count;
   struct heap_node* data;
   
} ;

struct SimpleHeap
{
   idx_t N; 
   T * data;
} ;

typedef struct SimpleHeap SimpleHeap;
typedef struct Heap Heap;
typedef struct heap_node heap_node;

void allocateSimpleHeap(SimpleHeap* H, idx_t n);

void allocateHeap(Heap* H, idx_t n);

void initSimpleHeap(SimpleHeap* H);

void initHeap(Heap* H);

void freeSimpleHeap(SimpleHeap * H);

void freeHeap(Heap * H);

void heapifyMaxSimpleHeap(SimpleHeap* H, idx_t node);

void heapifyMaxHeap(Heap* H, idx_t node);

void setRootMaxSimpleHeap(SimpleHeap * H, T val);

void setRootMaxHeap(Heap * H, FLOAT_TYPE val, idx_t array_idx);

void insertMaxHeap(Heap * H, FLOAT_TYPE val, idx_t array_idx);

void HeapSort(Heap* H);
void insertMaxHeap_InsertionSort(Heap * H,const FLOAT_TYPE val,const idx_t array_idx);


int cmpHeapNodes(const void* a, const void* b);
