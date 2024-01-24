#include "../include/kdtreeV2.h"
#include "../include/heap.h"
#include <math.h>
#include <stdint.h>
#include <time.h>

#define DEFAULT_LEAF_SIZE 30 


extern unsigned int data_dims;

FLOAT_TYPE eud_kdtree_v2(FLOAT_TYPE* restrict p1, FLOAT_TYPE* restrict p2){
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

FLOAT_TYPE eud_kdtree_v2_2(FLOAT_TYPE* restrict u, FLOAT_TYPE* restrict v)
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

void swap_kdnode_v2(kdnode_v2 *x, kdnode_v2 *y) {
    kdnode_v2 tmp;
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

void initialize_kdnodes_v2(kdnode_v2* node_array, FLOAT_TYPE* d, idx_t n )
{
    for(idx_t i = 0; i < n; ++i)
    {
        node_array[i].data = d + (i*data_dims);
        node_array[i].array_idx = i;
        node_array[i].lch = NULL;
        node_array[i].rch = NULL;
        node_array[i].parent = NULL;
        node_array[i].level = -1;
        node_array[i].is_leaf = 0;
        node_array[i].split_var = -1;
		node_array[i].node_list.data = NULL;
		node_array[i].node_list.count = 0;
    }
} 
/*
void initializePTRS(kdnode_v2** node_ptr_array, kdnode_v2* node_array, idx_t n )
{
    for(idx_t i = 0; i < n; ++i)
    {
        node_ptr_array[i] = node_array + i;
    }
}
*/

int cmpKDnodesV2(kdnode_v2* a, kdnode_v2* b, int var){
    
    FLOAT_TYPE res = a->data[var] - b->data[var];
    return (res > 0);
}

void printKDnodeV2(kdnode_v2* node)
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

int partition_kdnode_v2(kdnode_v2* arr, int low, int high, int split_var)
{
    kdnode_v2 pivot = arr[high];
    
    int i = (low - 1);
    for (int j = low; j <= high - 1; j++) {
        if (!cmpKDnodesV2(arr + j,&pivot,split_var)) {
            i++;
            swap_kdnode_v2(arr + i, arr + j);
        }
    }
    swap_kdnode_v2(arr + i + 1, arr + high);
    return (i + 1);
}

// Implementation of QuickSelect
int median_of_nodes_kdnode_v2(kdnode_v2* a, int left, int right, int split_var)
{
    //printf("----------\n");
    int k = left + ((right - left + 1)/2); 

    if(left == right) return left;
    if(left == (right - 1)){
        if(cmpKDnodesV2(a + left,a + right,split_var)) {swap_kdnode_v2(a + left, a + right);}
        return right;
    }
    while (left <= right) {
 
        // Partition a[left..right] around a pivot
        // and find the position of the pivot
        int pivotIndex = partition_kdnode_v2(a, left, right,split_var);
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

kdnode_v2* make_tree_kdnode_v2(kdnode_v2* t, int start, int end, kdnode_v2* parent, int level)
{
    kdnode_v2 *n = NULL;
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
		n -> is_leaf = 1;
        n -> parent = parent;
        n -> lch = NULL;
        n -> rch = NULL;
		size_t j = 0;
		n -> node_list.count = (size_t)(end - start + 1);
		n -> node_list.data = (kdnode_v2**)malloc(n -> node_list.count * sizeof(kdnode_v2*));
		for(int i = start; i <= end; ++i){
			t[i].parent = n;
			t[i].is_leaf = 1;
			t[i].lch = NULL;
			t[i].rch = NULL;
			n -> node_list.data[j] = t + i;
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

    median_idx = median_of_nodes_kdnode_v2(t, start, end, split_var);
    //printf("%d median idx\n", median_idx);
    if(median_idx > -1){
		swap_kdnode_v2(t + start, t + median_idx);
        //n = t + median_idx;
        n = t + start;
				//n->lch  = make_tree_kdnode_v2(t, start, median_idx - 1, n, level + 1);
				n->lch  = make_tree_kdnode_v2(t, start + 1, median_idx, n, level + 1);
				//n->rch = make_tree_kdnode_v2(t, median_idx + 1, end, n, level + 1);
				n->rch = make_tree_kdnode_v2(t, median_idx + 1, end, n, level + 1);
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

void knn_sub_tree_search_kdtree_v2(FLOAT_TYPE* point, kdnode_v2* root, heap * H)
{
	if(root -> is_leaf)
	{
		for(size_t i = 0; i < root -> node_list.count; ++i)
		{
			kdnode_v2* n = root -> node_list.data[i];
			__builtin_prefetch(root -> node_list.data + i + 1, 0, 3);
			FLOAT_TYPE distance = eud_kdtree_v2(point, n -> data);
			insert_max_heap(H, distance,n -> array_idx);
		}
		return;
	}

    FLOAT_TYPE current_distance = eud_kdtree_v2(point, root -> data);
    FLOAT_TYPE hp_distance = hyper_plane_dist(point, root -> data, root -> split_var);
    insert_max_heap(H, current_distance, root -> array_idx);
	__builtin_prefetch(root -> lch, 0, 3);
	__builtin_prefetch(root -> rch, 0, 3);

    int side = hp_distance > 0.f;

    switch (side)
    {
        case HP_LEFT_SIDE:
            if(root -> lch)
			{
				knn_sub_tree_search_kdtree_v2(point, root -> lch, H);
			}
            break;
        
        case HP_RIGHT_SIDE:
			if(root -> rch)
			{
				knn_sub_tree_search_kdtree_v2(point, root -> rch, H);
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
					knn_sub_tree_search_kdtree_v2(point, root -> rch, H);
				}
                break;
            
            case HP_RIGHT_SIDE:
                if(root -> lch) 
				{
					knn_sub_tree_search_kdtree_v2(point, root -> lch, H);
				}
                break;

            default:
                break;
        }
    }
    return;
}




heap knn_kdtree_v2(FLOAT_TYPE* point, kdnode_v2* kdtree_root, int maxk)
{
    heap H;
    allocate_heap(&H,maxk);
    init_heap(&H);
    knn_sub_tree_search_kdtree_v2(point, kdtree_root,&H);
    heap_sort(&H);
    return H;
}

kdnode_v2 * build_tree_kdtree_v2(kdnode_v2* kd_ptrs, size_t n, size_t dimensions )
{
    /*************************************************
     * Wrapper for make_tree function.               *
     * Simplifies interfaces and takes time measures *
     *************************************************/
   	data_dims = dimensions; 
    kdnode_v2* root = make_tree_kdnode_v2(kd_ptrs, 0, n-1, NULL ,0);
    return root;

}
