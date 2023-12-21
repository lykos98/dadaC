#include "../include/kdtree.h"
#include "../include/heap.h"
#include <math.h>
#include <stdint.h>
#include <time.h>


extern unsigned int data_dims;

int c = 0;
void swap(T* a, T* b){
    T tmp;
    memcpy(&tmp,a,sizeof(T));
    memcpy(a,b,sizeof(T));
    memcpy(b,&tmp,sizeof(T));
    return;
}

FLOAT_TYPE euclidean_distance(FLOAT_TYPE* restrict p1, FLOAT_TYPE* restrict p2){
    register FLOAT_TYPE d = 0;
    for(unsigned int i = 0; i<data_dims; ++i){
        register FLOAT_TYPE dd = (p1[i] - p2[i]);
        d += dd*dd;
    }
	return d;
	//return sqrt(d);
}

#ifdef USE_FLOAT32 
	typedef float v4f __attribute__ ((vector_size (16)));
#else
	typedef double v4f __attribute__ ((vector_size (32)));
#endif

FLOAT_TYPE euclidean_distance2(FLOAT_TYPE* restrict u, FLOAT_TYPE* restrict v)
{
    register float_t s;
    uint32_t i = 0;
    // manually unrolled loop, might be vectorized
    register v4f acc = {0., 0., 0., 0.};
    for (; i + 4 <= data_dims; i += 4) {
        register v4f _u = {u[i], u[i + 1], u[i + 2], u[i + 3]};
        register v4f _v = {v[i], v[i + 1], v[i + 2], v[i + 3]};
        register v4f _diff = _u - _v;
		acc = _diff*_diff;
    }

    s = acc[0] + acc[1] + acc[2] + acc[3];
    if (i < data_dims) {
        for(; i<data_dims; ++i) {
            float_t d = u[i] - v[i];
            s += d * d;
        }
    }
    return s;
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


/**
 * 
 * KDtree implementation 
 * 
 * 
*/

void initializeKDnodes(kd_node * node_array, FLOAT_TYPE* d, idx_t n )
{
    for(idx_t i = 0; i < n; ++i)
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

void initializePTRS(kd_node** node_ptr_array, kd_node* node_array, idx_t n )
{
    for(idx_t i = 0; i < n; ++i)
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
    printf("\t array_idx: %lu\n", (uint64_t)(node -> array_idx));
    printf("\t data: ");
    for(unsigned int i=0; i<data_dims; ++i) printf(" %f ",node->data[i]);
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
    //:w
    //int c = right - left + 1;
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
		#pragma omp parallel
		{
			#pragma omp single
			{
				#pragma omp task shared(n)
				n->lch  = make_tree(t, start, median_idx - 1, n, level + 1);
				#pragma omp task shared(n)
				n->rch = make_tree(t, median_idx + 1, end, n, level + 1);
			}
		}
        n -> split_var = split_var;
        n->parent = parent;
        n->level = level;
    }
	/*
    if(median_idx > -1){
        n = t[median_idx];
        n->lch  = make_tree(t, start, median_idx - 1, n, level + 1);
        n->rch = make_tree(t, median_idx + 1, end, n, level + 1);
        n -> split_var = split_var;
        n->parent = parent;
        n->level = level;
    }
	*/
    return n;
}

static inline FLOAT_TYPE hyper_plane_dist(FLOAT_TYPE* p1, FLOAT_TYPE* p2, int var)
{
    return p1[var] - p2[var];
}

inline int hyper_plane_side(FLOAT_TYPE* p1, FLOAT_TYPE* p2, int var)
{
    return p1[var] > p2[var];
}

void KNN_sub_tree_search(FLOAT_TYPE* point, kd_node* kdtree_root, Heap * H)
{
    //int split_var = kdtree_root -> split_var;
    FLOAT_TYPE current_distance = euclidean_distance(point, kdtree_root -> data);
    FLOAT_TYPE hp_distance = hyper_plane_dist(point, kdtree_root -> data, kdtree_root -> split_var);
    insertMaxHeap(H, current_distance, kdtree_root -> array_idx);
	__builtin_prefetch(kdtree_root -> lch, 0, 3);
	__builtin_prefetch(kdtree_root -> rch, 0, 3);

    int side = hp_distance > 0.f;

    //if(kdtree_root -> lch) KNN_sub_tree_search(point, kdtree_root -> lch, H);
    //if(kdtree_root -> rch) KNN_sub_tree_search(point, kdtree_root -> rch, H);
    switch (side)
    {
        case HP_LEFT_SIDE:
            if(kdtree_root -> lch)
			{
				__builtin_prefetch(kdtree_root -> lch -> data);
				KNN_sub_tree_search(point, kdtree_root -> lch, H);
			}
            break;
        
        case HP_RIGHT_SIDE:
			if(kdtree_root -> rch)
			{
				__builtin_prefetch(kdtree_root -> rch -> data);
				KNN_sub_tree_search(point, kdtree_root -> rch, H);
			}
            break;

        default:
            break;
    }
    FLOAT_TYPE max_d = H -> data[0].value;
    int c   = max_d > (hp_distance * hp_distance);
    //int c   = max_d > fabs(hp_distance);
    //if(!c) printf("%f %f\n",max_d, hp_distance*hp_distance);
    //if(c || (H -> count) < (H -> N))
    if(c)
    {

        switch (side)
        {
            case HP_LEFT_SIDE:
                if(kdtree_root -> rch) 
				{
					__builtin_prefetch(kdtree_root -> rch -> data);
					KNN_sub_tree_search(point, kdtree_root -> rch, H);
				}
                break;
            
            case HP_RIGHT_SIDE:
                if(kdtree_root -> lch) 
				{
					__builtin_prefetch(kdtree_root -> lch -> data);
					KNN_sub_tree_search(point, kdtree_root -> lch, H);
				}
                break;

            default:
                break;
        }
    }
    return;
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

kd_node * build_tree(kd_node** kd_ptrs, size_t n, size_t dimensions )
{
	
    

    /*************************************************
     * Wrapper for make_tree function.               *
     * Simplifies interfaces and takes time measures *
     *************************************************/
    
   	data_dims = dimensions; 

    kd_node* root = make_tree(kd_ptrs, 0, n-1, NULL ,0);

    return root;

}
