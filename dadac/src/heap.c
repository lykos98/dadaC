#include "../include/heap.h"
#include <stdint.h>
#include <time.h>

void swapHeapNode(heap_node* a, heap_node* b){
    heap_node tmp;
    memcpy(&tmp,a,sizeof(heap_node));
    memcpy(a,b,sizeof(heap_node));
    memcpy(b,&tmp,sizeof(heap_node));
    return;
}

void allocateSimpleHeap(SimpleHeap* H, idx_t n){
    H -> data = (T*)malloc(n*sizeof(T));
    H -> N = n;
    return;
}


void allocateHeap(Heap* H, idx_t n){
	H -> data = (heap_node*)malloc(n*sizeof(heap_node));
    H -> N = n;
    H -> count = 0;
    return;
}

void initSimpleHeap(SimpleHeap* H){
    for(idx_t i = 0; i < H -> N; ++i)
    {
        H -> data[i] = 1e12;
    }
    return;
}

void initHeap(Heap* H){
    for(idx_t i = 0; i < H -> N; ++i)
    {
        H -> data[i].value = 0.;
        H -> data[i].array_idx = MY_SIZE_MAX; 
    }
    return;
}

void freeSimpleHeap(SimpleHeap * H){ free(H -> data);}
void freeHeap(Heap * H){ free(H -> data);}

void heapifyMaxSimpleHeap(SimpleHeap* H, idx_t node){
    idx_t largest = node; 
    /*
    Found gratest between children of node and boundcheck if the node is a leaf 
    */
    if(HEAP_LCH(node) < H -> N){
        if(H -> data[HEAP_LCH(node)] > H -> data[largest] ) largest = HEAP_LCH(node);
    }
    if(HEAP_RCH(node) < H -> N){
        if(H -> data[HEAP_RCH(node)] > H -> data[largest] ) largest = HEAP_RCH(node);
    }
    if(largest == node){
        return;
    }
    else{
        swap(H -> data + node, H -> data + largest);
        heapifyMaxSimpleHeap(H, largest);
    }
}

void heapifyMaxHeap(Heap* H, idx_t node){
	idx_t nn = node;
    idx_t largest = nn; 
    /*
    Found gratest between children of node and boundcheck if the node is a leaf 
    */
	while(1)
	{
		largest = 	(HEAP_LCH(nn) < H -> N) && 
					(H -> data[HEAP_LCH(nn)].value > H -> data[largest].value ) ? HEAP_LCH(nn) : largest;

		largest = 	(HEAP_RCH(nn) < H -> N) && 
					(H -> data[HEAP_RCH(nn)].value > H -> data[largest].value ) ? HEAP_RCH(nn) : largest;
		if(largest != nn) 
		{
			swapHeapNode(H -> data + nn, H -> data + largest);
			nn = largest;
		}
		else
		{
			break;
		}
	}

    //if(HEAP_LCH(node) < H -> N){
    //    //if(H -> data[HEAP_LCH(node)].value > H -> data[largest].value ) largest = HEAP_LCH(node);
	//	largest = (H -> data[HEAP_LCH(nn)].value > H -> data[largest].value ) ? HEAP_LCH(nn) : largest;
    //}
    //if(HEAP_RCH(node) < H -> N){
    //    //if(H -> data[HEAP_RCH(node)].value > H -> data[largest].value ) largest = HEAP_RCH(node);
	//	largest = (H -> data[HEAP_RCH(nn)].value > H -> data[largest].value ) ? HEAP_RCH(nn) : largest;
    //}
    //if(largest == node){
    //    return;
    //}
    //else{
    //    swapHeapNode(H -> data + node, H -> data + largest);
    //    heapifyMaxHeap(H, largest);
    //}
}

void setRootMaxSimpleHeap(SimpleHeap * H, T val){
    H -> data[0] = val;
    heapifyMaxSimpleHeap(H,0);
    return;
}

void setRootMaxHeap(Heap * H, FLOAT_TYPE val, idx_t array_idx){
    H -> data[0].value = val;
    H -> data[0].array_idx = array_idx;
    heapifyMaxHeap(H,0);
    return;
}

void insertMaxHeap_InsertionSort(Heap * H,const FLOAT_TYPE val,const idx_t array_idx){
	heap_node tmpNode = {.value = val, .array_idx = array_idx};
	if(H -> count < H -> N)
	{
		idx_t idx = H -> count;
		H -> data[idx] = tmpNode;
		++(H -> count);
		while(idx >= 1)
		{
			if(H -> data[idx].value < H -> data[idx - 1].value)
			{
				swapHeapNode((H -> data) + idx, (H -> data) + idx - 1);
				idx--;
			}
			else
			{
				break;
			}
		}

	}
	else
	{
		if(H -> data[H -> count - 1].value > val)
		{
			idx_t idx = H -> count - 1;
			H -> data[idx] = tmpNode;
			while(idx >= 1)
			{
				if(H -> data[idx].value < H -> data[idx - 1].value)
				{
					swapHeapNode(H -> data + idx, H -> data + idx - 1);
					idx--;
				}
				else
				{
					break;
				}
			}
		}
	}
	return;
}

void insertMaxHeap(Heap * H,const FLOAT_TYPE val,const idx_t array_idx){
	int c1 = H -> count < H -> N;
	int c2 = (val < H -> data[0].value) && (!c1);
	int ctot = c1 + 2*c2;
	switch (ctot) {
		case 1:
			{
				idx_t node = H->count;
				++(H -> count);
				H -> data[node].value = val;
				H -> data[node].array_idx = array_idx;
				/*
				* Push up the node through the heap 
				*/
				while(node && H -> data[node].value > H -> data[HEAP_PARENT(node)].value)
				{
					swapHeapNode(H -> data + node, H -> data + HEAP_PARENT(node));
					node = HEAP_PARENT(node);
					//if(node == 0) break;
			}
			}
		break;
		case 2: 
			{
				setRootMaxHeap(H,val,array_idx);
			}
			break;
		default:
			break;
	}
    //if(H -> count < H -> N){
    //    idx_t node = H->count;
    //    ++(H -> count);
    //    H -> data[node].value = val;
    //    H -> data[node].array_idx = array_idx;
    //    /*
    //    * Push up the node through the heap 
    //    */
    //    while(node && H -> data[node].value > H -> data[HEAP_PARENT(node)].value)
    //    {
    //        swapHeapNode(H -> data + node, H -> data + HEAP_PARENT(node));
    //        node = HEAP_PARENT(node);
    //        //if(node == 0) break;
    //    }
    //    return;
    //}
    //else if (val < H -> data[0].value)
    //{
    //    setRootMaxHeap(H,val,array_idx);
    //    return;
    //}
    //else
    //{
    //    return;
    //}
	
    
}

#ifdef USE_FLOAT32
	#define EPS 5.96e-08
#else
	#define EPS 2.11e-16
#endif

int cmpHeapNodes(const void* a, const void* b)
{
	const heap_node* aa = (const heap_node*)a;
	const heap_node* bb = (const heap_node*)b;
	int val =  (aa -> value > bb -> value) - (aa -> value < bb -> value); 
	//return vv; 
	return val;


}


void HeapSort(Heap* H){
    idx_t n = H -> N;
	qsort(H -> data, n, sizeof(heap_node),cmpHeapNodes);
	//for(idx_t i= (H -> N) - 1; i > 0; --i)
    //{
    //    swapHeapNode(H -> data, H -> data + i);
    //    H -> N = i;
    //    heapifyMaxHeap(H,0);
    //}
    //H -> N = n;
}
