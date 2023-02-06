#include "../include/kdtree.h"
#include <time.h>

#define HEAP_LCH(x) (2*x + 1)
#define HEAP_RCH(x) (2*x + 2)
#define HEAP_PARENT(x) (x-1)/2

#define HP_LEFT_SIDE 0
#define HP_RIGHT_SIDE 1

extern unsigned int data_dims;

int c = 0;
void swap(T* a, T* b){
    T tmp;
    memcpy(&tmp,a,sizeof(T));
    memcpy(a,b,sizeof(T));
    memcpy(b,&tmp,sizeof(T));
    return;
}

FLOAT_TYPE euclidean_distance(FLOAT_TYPE* p1, FLOAT_TYPE* p2){
    FLOAT_TYPE d = 0;
    for(int i = 0; i<data_dims; ++i){
        FLOAT_TYPE dd = (p1[i] - p2[i]);
        d += dd*dd;
    }
    return d;
}

void swapHeapNode(heap_node* a, heap_node* b){
    heap_node tmp;
    memcpy(&tmp,a,sizeof(heap_node));
    memcpy(a,b,sizeof(heap_node));
    memcpy(b,&tmp,sizeof(heap_node));
    return;
}

void swap_kd_node_ptrs(kd_node **x, kd_node **y) {
    kd_node* tmp;
    tmp = *x;
    *x = *y;
    *y = tmp;
    //memcpy(&tmp,  x, sizeof(tmp));
    //memcpy(x, y, sizeof(tmp));
    //memcpy(y, &tmp,  sizeof(tmp));
}

void allocateSimpleHeap(SimpleHeap* H, size_t n){
    H -> data = (T*)malloc(n*sizeof(T));
    H -> N = n;
    return;
}

void allocateHeap(Heap* H, size_t n){
    H -> data = (heap_node*)malloc(n*sizeof(heap_node));
    H -> N = n;
    H -> count = 0;
    return;
}

void initSimpleHeap(SimpleHeap* H){
    for(size_t i = 0; i < H -> N; ++i)
    {
        H -> data[i] = 1e12;
    }
    return;
}

void initHeap(Heap* H){
    for(size_t i = 0; i < H -> N; ++i)
    {
        H -> data[i].value = 0.;
        H -> data[i].array_idx = ULONG_MAX; 
    }
    return;
}

void freeSimpleHeap(SimpleHeap * H){ free(H -> data);}
void freeHeap(Heap * H){ free(H -> data);}

void heapifyMaxSimpleHeap(SimpleHeap* H, size_t node){
    size_t largest = node; 
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

void heapifyMaxHeap(Heap* H, size_t node){
    size_t largest = node; 
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

void setRootMaxHeap(Heap * H, FLOAT_TYPE val, size_t array_idx){
    H -> data[0].value = val;
    H -> data[0].array_idx = array_idx;
    heapifyMaxHeap(H,0);
    return;
}

void insertMaxHeap(Heap * H, FLOAT_TYPE val, size_t array_idx){
    if (H -> count == 0)
    {
        ++(H -> count);
        H -> data[0].value = val;
        H -> data[0].array_idx = array_idx;
    }
    else if(H -> count < H -> N){
        size_t node = H->count;
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

/**
 * 
 * KDtree implementation 
 * 
 * 
*/

void initializeKDnodes(kd_node * node_array, FLOAT_TYPE* d, size_t n )
{
    for(size_t i = 0; i < n; ++i)
    {
        node_array[i].data = d + (i*data_dims);
        node_array[i].array_idx = i;
        node_array[i].lch = NULL;
        node_array[i].rch = NULL;
        node_array[i].parent = NULL;
        node_array[i].level = -1;
        node_array[i].split_var = -1;
    }
}

void initializePTRS(kd_node** node_ptr_array, kd_node* node_array, size_t n )
{
    for(size_t i = 0; i < n; ++i)
    {
        node_ptr_array[i] = node_array + i;
    }
}

int cmpKDnodes(kd_node* a, kd_node* b, int var){
    
    FLOAT_TYPE res = a->data[var] - b->data[var];
    return (res > 0);
}

void printKDnode(kd_node* node)
{
    printf("Node %p:\n",node);
    printf("\t array_idx: %ld\n", node -> array_idx);
    printf("\t data: ");
    for(int i=0; i<data_dims; ++i) printf(" %f ",node->data[i]);
    printf("\n");
    printf("\t parent: %p\n",node -> parent);
    printf("\t lch: %p\n",node -> lch);
    printf("\t rch: %p\n",node -> rch);
    printf("\t level: %d\n", node -> level);
    printf("\t split_var: %d\n", node -> split_var);
    printf("\n");
}

// Standard Lomuto partition function

int partition(kd_node** arr, int low, int high, int split_var)
{
    kd_node* pivot = arr[high];
    
    int i = (low - 1);
    for (int j = low; j <= high - 1; j++) {
        if (!cmpKDnodes(arr[j],pivot,split_var)) {
            i++;
            swap_kd_node_ptrs(arr + i, arr + j);
        }
    }
    swap_kd_node_ptrs(arr + i + 1, arr + high);
    return (i + 1);
}

// Implementation of QuickSelect
int medianOfNodes(kd_node** a, int left, int right, int split_var)
{
    //printf("----------\n");
    int k = left + ((right - left + 1)/2); 
    int c = right - left + 1;
    //if(c < 20){
    //    v = split_var;
    //    qsort(a + left, c, sizeof(kd_node*),cmpKDN);
    //    return k;

    //}

    if(left == right) return left;
    if(left == (right - 1)){
        if(cmpKDnodes(a[left],a[right],split_var)) {swap_kd_node_ptrs(a + left, a + right);}
        return right;
    }
    //if(c == 3){
    //    printKDnode(a[left]);
    //    printKDnode(a[left + 1]);
    //    printKDnode(a[left + 2]);
    //}
    while (left <= right) {
 
        // Partition a[left..right] around a pivot
        // and find the position of the pivot
        int pivotIndex = partition(a, left, right,split_var);
        //printf("%d %d %d %d\n",left, right, k, pivotIndex);

        // If pivot itself is the k-th smallest element
        if (pivotIndex == k)
            return pivotIndex;
 
        // If there are more than k-1 elements on
        // left of pivot, then k-th smallest must be
        // on left side.
        else if (pivotIndex > k)
            right = pivotIndex - 1;
 
        // Else k-th smallest is on right side.
        else
            left = pivotIndex + 1;
    }
    //printf("Nope\n");
    return -1;
}

kd_node* make_tree(kd_node** t, int start, int end, kd_node* parent, int level)
{
    kd_node *n = NULL;
    int split_var = level % data_dims; 
    //printf("%d \n",level);

    int median_idx = -1;
    if ((end - start) < 0) return 0;
    if (end  == start) {
        n = t[start];
        n -> split_var = split_var;
        n->parent = parent;
        n->level = level;
        n -> lch = NULL;
        n -> rch = NULL;
        return n;
    }
    median_idx = medianOfNodes(t, start, end, split_var);
    //printf("%d median idx\n", median_idx);
    if(median_idx > -1){
        n = t[median_idx];
        n->lch  = make_tree(t, start, median_idx - 1, n, level + 1);
        n->rch = make_tree(t, median_idx + 1, end, n, level + 1);
        n -> split_var = split_var;
        n->parent = parent;
        n->level = level;
    }
    return n;
}

inline FLOAT_TYPE hyper_plane_dist(FLOAT_TYPE* p1, FLOAT_TYPE* p2, int var)
{
    return p1[var] - p2[var];
}

inline int hyper_plane_side(FLOAT_TYPE* p1, FLOAT_TYPE* p2, int var)
{
    return p1[var] > p2[var];
}

void KNN_sub_tree_search(FLOAT_TYPE* point, kd_node* kdtree_root, Heap * H)
{
    int split_var = kdtree_root -> split_var;
    FLOAT_TYPE current_distance = euclidean_distance(point, kdtree_root -> data);
    FLOAT_TYPE hp_distance = hyper_plane_dist(point, kdtree_root -> data, kdtree_root -> split_var);
    insertMaxHeap(H, current_distance, kdtree_root -> array_idx);

    int side = hp_distance > 0.f;

    //if(kdtree_root -> lch) KNN_sub_tree_search(point, kdtree_root -> lch, H);
    //if(kdtree_root -> rch) KNN_sub_tree_search(point, kdtree_root -> rch, H);
    switch (side)
    {
        case HP_LEFT_SIDE:
            if(kdtree_root -> lch) KNN_sub_tree_search(point, kdtree_root -> lch, H);
            break;
        
        case HP_RIGHT_SIDE:
            if(kdtree_root -> rch) KNN_sub_tree_search(point, kdtree_root -> rch, H);
            break;

        default:
            break;
    }
    FLOAT_TYPE max_d = H -> data[0].value;
    int c   = max_d > (hp_distance * hp_distance);
    //if(!c) printf("%f %f\n",max_d, hp_distance*hp_distance);
    if(c || (H -> count) < (H -> N))
    {

        switch (side)
        {
            case HP_LEFT_SIDE:
                if(kdtree_root -> rch) KNN_sub_tree_search(point, kdtree_root -> rch, H);
                break;
            
            case HP_RIGHT_SIDE:
                if(kdtree_root -> lch) KNN_sub_tree_search(point, kdtree_root -> lch, H);
                break;

            default:
                break;
        }
    }
    return;
}


void HeapSort(Heap* H){
    size_t n = H -> N;
    for(size_t i= (H -> N) - 1; i > 0; --i)
    {
        swapHeapNode(H -> data, H -> data + i);
        H -> N = i;
        heapifyMaxHeap(H,0);
    }
    H -> N = n;
}


Heap KNN(FLOAT_TYPE* point, kd_node* kdtree_root, int maxk)
{
    Heap H;
    allocateHeap(&H,maxk);
    initHeap(&H);
    KNN_sub_tree_search(point, kdtree_root,&H);
    HeapSort(&H);
    return H;
}

kd_node * build_tree(kd_node** kd_ptrs, size_t n)
{
    

    /*************************************************
     * Wrapper for make_tree function.               *
     * Simplifies interfaces and takes time measures *
     *************************************************/
    
    
    struct timespec start_tot, finish_tot;
    double elapsed_tot;
    printf("Building the KDtree:\n");
    clock_gettime(CLOCK_MONOTONIC, &start_tot);

    kd_node* root = make_tree(kd_ptrs, 0, n-1, NULL ,0);

    clock_gettime(CLOCK_MONOTONIC, &finish_tot);
    elapsed_tot = (finish_tot.tv_sec - start_tot.tv_sec);
    elapsed_tot += (finish_tot.tv_nsec - start_tot.tv_nsec) / 1000000000.0;
    printf("\tTotal time: %.3lfs\n\n", elapsed_tot);
    return root;

}
