#include "../include/kdtreeV2.h"
#include "../include/heap.h"
#include <math.h>
#include <stdint.h>
#include <time.h>

#define DEFAULT_LEAF_SIZE 30 


extern unsigned int data_dims;

FLOAT_TYPE eud_kdTreeV2(FLOAT_TYPE* restrict p1, FLOAT_TYPE* restrict p2){
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

FLOAT_TYPE eud_kdTreeV2_2(FLOAT_TYPE* restrict u, FLOAT_TYPE* restrict v)
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

FLOAT_TYPE* swapMem_kdv2;

void swap_kdNodeV2(kdNodeV2 *x, kdNodeV2 *y) {
    kdNodeV2 tmp;
    tmp = *x;
    *x = *y;
    *y = tmp;

	#ifdef SWMEM
		memcpy(swapMem_kdv2, x -> data, sizeof(FLOAT_TYPE)*data_dims);
		memcpy(x -> data, y -> data, sizeof(FLOAT_TYPE)*data_dims);
		memcpy(y -> data, swapMem_kdv2, sizeof(FLOAT_TYPE)*data_dims);
		FLOAT_TYPE* tmpPtr = x -> data;
		x -> data = y -> data;
		y -> data = tmpPtr;
	#endif
}


/**
 * 
 * KDtree implementation 
 * 
 * 
*/

void initializeKDnodesV2(kdNodeV2* node_array, FLOAT_TYPE* d, idx_t n )
{
    for(idx_t i = 0; i < n; ++i)
    {
        node_array[i].data = d + (i*data_dims);
        node_array[i].array_idx = i;
        node_array[i].lch = NULL;
        node_array[i].rch = NULL;
        node_array[i].parent = NULL;
        node_array[i].level = -1;
        node_array[i].isLeaf = 0;
        node_array[i].split_var = -1;
		node_array[i].nodeList.data = NULL;
		node_array[i].nodeList.count = 0;
    }
} 
/*
void initializePTRS(kdNodeV2** node_ptr_array, kdNodeV2* node_array, idx_t n )
{
    for(idx_t i = 0; i < n; ++i)
    {
        node_ptr_array[i] = node_array + i;
    }
}
*/

int cmpKDnodesV2(kdNodeV2* a, kdNodeV2* b, int var){
    
    FLOAT_TYPE res = a->data[var] - b->data[var];
    return (res > 0);
}

void printKDnodeV2(kdNodeV2* node)
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

int partition_kdNodeV2(kdNodeV2* arr, int low, int high, int split_var)
{
    kdNodeV2 pivot = arr[high];
    
    int i = (low - 1);
    for (int j = low; j <= high - 1; j++) {
        if (!cmpKDnodesV2(arr + j,&pivot,split_var)) {
            i++;
            swap_kdNodeV2(arr + i, arr + j);
        }
    }
    swap_kdNodeV2(arr + i + 1, arr + high);
    return (i + 1);
}

// Implementation of QuickSelect
int medianOfNodes_kdNodeV2(kdNodeV2* a, int left, int right, int split_var)
{
    //printf("----------\n");
    int k = left + ((right - left + 1)/2); 

    if(left == right) return left;
    if(left == (right - 1)){
        if(cmpKDnodesV2(a + left,a + right,split_var)) {swap_kdNodeV2(a + left, a + right);}
        return right;
    }
    while (left <= right) {
 
        // Partition a[left..right] around a pivot
        // and find the position of the pivot
        int pivotIndex = partition_kdNodeV2(a, left, right,split_var);
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
    return -1;
}

kdNodeV2* make_tree_kdNodeV2(kdNodeV2* t, int start, int end, kdNodeV2* parent, int level)
{
    kdNodeV2 *n = NULL;
    int split_var = level % data_dims; 
	FLOAT_TYPE max_diff = -999999.;
	for(unsigned int v = 0; v < data_dims; ++v)
	{
		FLOAT_TYPE max_v = -9999999.;
		FLOAT_TYPE min_v = 9999999.;
		for(int i = start; i <= end; ++i)
		{
			max_v = t[i].data[v] > max_v ? t[i].data[v] : max_v;
			min_v = t[i].data[v] < min_v ? t[i].data[v] : min_v;
		}
		if((max_v - min_v) > max_diff)
		{
			max_diff = max_v - min_v;
			split_var = v;
		}
	}
	
	#ifdef SWMEM	
		if(parent == NULL)
		{
			swapMem_kdv2 = (FLOAT_TYPE*)malloc(sizeof(FLOAT_TYPE)*data_dims);
		}
	#endif
	
	
	
	if(end - start < DEFAULT_LEAF_SIZE)
	{
		n =  t + start;
		n -> isLeaf = 1;
        n -> parent = parent;
        n -> lch = NULL;
        n -> rch = NULL;
		size_t j = 0;
		n -> nodeList.count = (size_t)(end - start + 1);
		n -> nodeList.data = (kdNodeV2**)malloc(n -> nodeList.count * sizeof(kdNodeV2*));
		for(int i = start; i <= end; ++i){
			t[i].parent = n;
			t[i].isLeaf = 1;
			t[i].lch = NULL;
			t[i].rch = NULL;
			n -> nodeList.data[j] = t + i;
			++j;
		}
		return n;
		
	}
	
	


    int median_idx = -1;
	
    //if ((end - start) < 0) return 0;
    //if (end  == start) {
    //    n = t + start;
    //    n -> split_var = split_var;
    //    n->parent = parent;
    //    n->level = level;
    //    n -> lch = NULL;
    //    n -> rch = NULL;
    //    return n;
    //}

    median_idx = medianOfNodes_kdNodeV2(t, start, end, split_var);
    //printf("%d median idx\n", median_idx);
    if(median_idx > -1){
		swap_kdNodeV2(t + start, t + median_idx);
        //n = t + median_idx;
        n = t + start;
				//n->lch  = make_tree_kdNodeV2(t, start, median_idx - 1, n, level + 1);
				n->lch  = make_tree_kdNodeV2(t, start + 1, median_idx, n, level + 1);
				//n->rch = make_tree_kdNodeV2(t, median_idx + 1, end, n, level + 1);
				n->rch = make_tree_kdNodeV2(t, median_idx + 1, end, n, level + 1);
        n -> split_var = split_var;
        n->parent = parent;
        n->level = level;
    }
	
	#ifdef SWMEM
		if(parent == NULL)
		{
			swapMem_kdv2 = malloc(sizeof(FLOAT_TYPE)*data_dims);
		}
	#endif

    return n;
}

static inline FLOAT_TYPE hyper_plane_dist(FLOAT_TYPE* p1, FLOAT_TYPE* p2, int var)
{
    return p1[var] - p2[var];
}

static inline int hyper_plane_side(FLOAT_TYPE* p1, FLOAT_TYPE* p2, int var)
{
    return p1[var] > p2[var];
}

void KNN_sub_tree_search_kdTreeV2(FLOAT_TYPE* point, kdNodeV2* root, Heap * H)
{
	if(root -> isLeaf)
	{
		for(size_t i = 0; i < root -> nodeList.count; ++i)
		{
			kdNodeV2* n = root -> nodeList.data[i];
			__builtin_prefetch(root -> nodeList.data + i + 1, 0, 3);
			FLOAT_TYPE distance = eud_kdTreeV2(point, n -> data);
			insertMaxHeap(H, distance,n -> array_idx);
		}
		return;
	}

    FLOAT_TYPE current_distance = eud_kdTreeV2(point, root -> data);
    FLOAT_TYPE hp_distance = hyper_plane_dist(point, root -> data, root -> split_var);
    insertMaxHeap(H, current_distance, root -> array_idx);
	__builtin_prefetch(root -> lch, 0, 3);
	__builtin_prefetch(root -> rch, 0, 3);

    int side = hp_distance > 0.f;

    switch (side)
    {
        case HP_LEFT_SIDE:
            if(root -> lch)
			{
				KNN_sub_tree_search_kdTreeV2(point, root -> lch, H);
			}
            break;
        
        case HP_RIGHT_SIDE:
			if(root -> rch)
			{
				KNN_sub_tree_search_kdTreeV2(point, root -> rch, H);
			}
            break;

        default:
            break;
    }
    FLOAT_TYPE max_d = H -> data[0].value;
    int c   = max_d > (hp_distance * hp_distance);

    //if(c || (H -> count) < (H -> N))
    if(c)
    {

        switch (side)
        {
            case HP_LEFT_SIDE:
                if(root -> rch) 
				{
					KNN_sub_tree_search_kdTreeV2(point, root -> rch, H);
				}
                break;
            
            case HP_RIGHT_SIDE:
                if(root -> lch) 
				{
					KNN_sub_tree_search_kdTreeV2(point, root -> lch, H);
				}
                break;

            default:
                break;
        }
    }
    return;
}




Heap KNN_kdTreeV2(FLOAT_TYPE* point, kdNodeV2* kdtree_root, int maxk)
{
    Heap H;
    allocateHeap(&H,maxk);
    initHeap(&H);
    KNN_sub_tree_search_kdTreeV2(point, kdtree_root,&H);
    HeapSort(&H);
    return H;
}

kdNodeV2 * build_tree_kdTreeV2(kdNodeV2* kd_ptrs, size_t n, size_t dimensions )
{
    /*************************************************
     * Wrapper for make_tree function.               *
     * Simplifies interfaces and takes time measures *
     *************************************************/
   	data_dims = dimensions; 
    kdNodeV2* root = make_tree_kdNodeV2(kd_ptrs, 0, n-1, NULL ,0);
    return root;

}
