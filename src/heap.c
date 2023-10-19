#include "../include/heap.h"
#include <stdint.h>
#include <time.h>

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
    idx_t largest = node; 
    /*
    Found gratest between children of node and boundcheck if the node is a leaf 
    */
    if(HEAP_LCH(node) < H -> N){
        if(H -> data[HEAP_LCH(node)].value > H -> data[largest].value ) largest = HEAP_LCH(node);
    }
    if(HEAP_RCH(node) < H -> N){
        if(H -> data[HEAP_RCH(node)].value > H -> data[largest].value ) largest = HEAP_RCH(node);
    }
    if(largest == node){
        return;
    }
    else{
        swapHeapNode(H -> data + node, H -> data + largest);
        heapifyMaxHeap(H, largest);
    }
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

void insertMaxHeap(Heap * H, FLOAT_TYPE val, idx_t array_idx){
    if (H -> count == 0)
    {
        ++(H -> count);
        H -> data[0].value = val;
        H -> data[0].array_idx = array_idx;
    }
    else if(H -> count < H -> N){
        idx_t node = H->count;
        ++(H -> count);
        H -> data[node].value = val;
        H -> data[node].array_idx = array_idx;
        /*
        * Push up the node through the heap 
        */
        while(H -> data[node].value > H -> data[HEAP_PARENT(node)].value)
        {
            swapHeapNode(H -> data + node, H -> data + HEAP_PARENT(node));
            node = HEAP_PARENT(node);
            if(node == 0) break;
        }
        return;
    }
    else if (val < H -> data[0].value)
    {
        setRootMaxHeap(H,val,array_idx);
        return;
    }
    else
    {
        return;
    }
    
}

void HeapSort(Heap* H){
    idx_t n = H -> N;
    for(idx_t i= (H -> N) - 1; i > 0; --i)
    {
        swapHeapNode(H -> data, H -> data + i);
        H -> N = i;
        heapifyMaxHeap(H,0);
    }
    H -> N = n;
}
