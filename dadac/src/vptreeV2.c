#include <stdint.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <stdint.h>
#include  "../include/vptreeV2.h"

#define DEFAULT_LEAF_SIZE 30 

#ifdef VOPT
	#define SWMEM
#endif

/*
 * 
 * KDtree implementation 
 * 
 * 
 */


void* swapMem;

void swap_vpnode_v2_ptrs(vpnode_v2 *x, vpnode_v2 *y) {
    vpnode_v2 tmp;
    tmp = *x;
    *x = *y;
    *y = tmp;
	#ifdef SWMEM
		memcpy(swapMem, x -> data, x -> __bytesize);
		memcpy(x -> data, y -> data, x -> __bytesize);
		memcpy(y -> data, swapMem, x -> __bytesize);
		
		void* tmpPtr = x -> data;
		x -> data = y -> data;
		y -> data = tmpPtr;
	#endif
	//vpnode_v2 tmpNode = *(*x);
	//*(*x) = *(*y);
	//*(*y)  = tmpNode;
	 
}

/*
 * 
 * KDtree implementation 
 * 
 * 
 */

void initialize_vpnode_v2_array(vpnode_v2* node_array, void* data, idx_t n, idx_t bytes_per_elements)
{
    for(idx_t i = 0; i < n; ++i)
    {
        node_array[i].data = data + (i*bytes_per_elements);
        node_array[i].array_idx = i;
        node_array[i].inside = NULL;
        node_array[i].outside = NULL;
        node_array[i].parent = NULL;
        node_array[i].mu = 0;
        node_array[i].__dist = 0;
        node_array[i].is_leaf = 0;
        node_array[i].node_list.data = NULL;
        node_array[i].node_list.count = 0;
		node_array[i].__bytesize = bytes_per_elements;
		
    }

}

void initialize_vpnode_v2_ptrs(vpnode_v2** pointers_array, vpnode_v2* node_array, idx_t n)
{
    for(idx_t i = 0; i < n; ++i) pointers_array[i] = node_array + i;

}

int cmp_vpnode_v2(vpnode_v2* point, vpnode_v2* pivot)
{
    float_t res = point->__dist - pivot->__dist;
    return (res > 0);
}

// Standard Lomuto partition function

int partition_vpnode_v2(vpnode_v2* arr, int low, int high)
{
    vpnode_v2 pivot = arr[high];
    
    int i = (low - 1);
    for (int j = low; j <= high - 1; j++) {
        if (!cmp_vpnode_v2(arr + j,&pivot)) {
            i++;
            swap_vpnode_v2_ptrs(arr + i, arr + j);
        }
    }
    swap_vpnode_v2_ptrs(arr + i + 1, arr + high);
    return (i + 1);

}

int median_of_vpnode_v2(vpnode_v2* a, int left, int right)
{
    //printf("----------\n");
    int k = left + ((right - left + 1)/2); 

    if(left == right) return left;
    if(left == (right - 1)){
        if(cmp_vpnode_v2(a + left,a + right)) {swap_vpnode_v2_ptrs(a + left, a + right);}
        return right;
    }
    while (left <= right) {
 
        // Partition a[left..right] around a pivot
        // and find the position of the pivot
        int pivotIndex = partition_vpnode_v2(a, left, right);
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


vpnode_v2* build_vptree_v2(vpnode_v2* t, int start, int end, vpnode_v2* parent, float_t (*metric)(void*, void*))
{
	#ifdef SWMEM
		if(parent == NULL)
		{
			swapMem = malloc((t + start) ->__bytesize);
		}
	#endif

    vpnode_v2 *n = NULL;
    //printf("%d \n",level);

    int median_idx = -1;
    if ((end - start) < 0) return NULL;

	if(end - start < DEFAULT_LEAF_SIZE)
	//if(end - start < -1)
	{
		n =  t + start;
		n -> is_leaf = 1;
        n -> parent = parent;
        n -> inside = NULL;
        n -> outside = NULL;
		n -> mu = 0;
		size_t j = 0;
		n -> node_list.count = (size_t)(end - start + 1);
		#ifdef VOPT
			n -> node_list.indexes = (idx_t*)malloc(n -> node_list.count * sizeof(idx_t));
			n -> node_list.start_ptr = t[start].data;
			n -> node_list.end_ptr 	= t[end].data;
		#else
			n -> node_list.data = (vpnode_v2**)malloc(n -> node_list.count * sizeof(vpnode_v2*));
		#endif
		for(int i = start; i <= end; ++i){
			t[i].parent = n;
			t[i].is_leaf = 1;
			t[i].inside = NULL;
			t[i].outside = NULL;
			#ifdef VOPT
				n -> node_list.indexes[j] = t[i].array_idx;
			#else
				n -> node_list.data[j] = t + i;
			#endif

			++j;
		}
		return n;
		
	}

	//int vp_idx = start + (end - start + 1)/2; 
	int vp_idx = start; 

	//compute distances
	for(int i = start; i <= end; ++i) t[i].__dist = metric(t[vp_idx].data,t[i].data);
	
	//now swap the vp_idx in the first place to retrieve it later
	//swap_vptreeNode_ptrs(t + vp_idx, t + start);
	//compute the median, as byproduct obtain the array partitioned on inside and outside BUT, with the median in the median place
    median_idx = median_of_vpnode_v2(t, start, end);
	//printf("Median: %d %.4lf\n\n", median_idx - start, t[median_idx]. __dist );
	//for(int i = start; i <= end; ++i) 
	//
	float_t mu = t[median_idx].__dist;
	//now swap start with the median
	
    if(median_idx > -1){
        n = t + vp_idx;
		n->mu = mu;
		n->inside  = build_vptree_v2(t, start + 1, median_idx, n, metric);
		n->outside = build_vptree_v2(t, median_idx + 1, end, n, metric);
				//#pragma omp taskwait

        n->parent = parent;
    }

	#ifdef SWMEM
	if(parent == NULL)
	{
		free(swapMem);
	}
	#endif


    return n;

}


void knn_sub_vptree_search_v2(void* point, vpnode_v2* root, heap * H, float_t (*metric)(void*,void*))
{
	
	#ifdef VOPT	
	if(root -> is_leaf)
	{
		
		for(idx_t i = 0; i < root -> node_list.count; ++i)
		{

			idx_t j = root -> node_list.indexes[i];
			float_t distance = metric(point, root -> node_list.start_ptr + (i*(root -> __bytesize)));
			insert_max_heap(H, distance,j);
		}
		return;
	}
	#else
		if(root -> is_leaf)
		{
			for(size_t i = 0; i < root -> node_list.count; ++i)
			{
				__builtin_prefetch(root -> node_list.data + i + 1, 0, 3);
				vpnode_v2* n = root -> node_list.data[i];
				float_t distance = metric(point, n -> data);
				insert_max_heap(H, distance,n -> array_idx);
			}
			return;
		}
	#endif


    float_t current_distance = metric(point, root -> data);
    insert_max_heap(H, current_distance, root -> array_idx);
	#define INSIDE 0
	#define OUTSIDE 1
	int side = (current_distance > root -> mu);

    switch (side)
    {
        case INSIDE:
            knn_sub_vptree_search_v2(point, root -> inside, H, metric);
            //if(root -> inside) knn_sub_vptree_search_v2(point, root -> inside, H, metric);
            break;
        case OUTSIDE:
            knn_sub_vptree_search_v2(point, root -> outside, H, metric);
            //if(root -> outside) knn_sub_vptree_search_v2(point, root -> outside, H, metric);
            break;

        default:
            break;
    }
    float_t tau = H -> data[0].value;
	//retrieve the maximum distance found so far
	//int heapNotFull = (H -> count) < (H -> N);
	switch (side)
	{
		case INSIDE:
			//if 	( ((current_distance + tau) > root->mu  || heapNotFull)) knn_sub_vptree_search_v2(point, root -> outside, H, metric);
			if 	( (current_distance + tau) > root->mu  ) knn_sub_vptree_search_v2(point, root -> outside, H, metric);
			//if 	( root -> outside && ((current_distance + tau) > root->mu  || heapNotFull)) knn_sub_vptree_search_v2(point, root -> outside, H, metric);
			break;

		case OUTSIDE:
			//if 	( root -> inside && ((current_distance - tau) < root->mu  || heapNotFull)) knn_sub_vptree_search(point, root -> inside, H, metric);
			//if 	( (current_distance  < (root->mu + tau)  || heapNotFull)) knn_sub_vptree_search_v2(point, root -> inside, H, metric);
			if 	( (current_distance  < (root->mu + tau) )) knn_sub_vptree_search_v2(point, root -> inside, H, metric);
			//if 	( root -> inside && (current_distance  < (root->mu + tau)  || heapNotFull)) knn_sub_vptree_search_v2(point, root -> inside, H, metric);
			break;
		
		default:
			break;
	}
    return;
	

}

heap knn_vptree_v2(void* point, vpnode_v2* root, int maxk, float_t (*metric)(void*, void*))
{
    heap H;
    allocate_heap(&H,maxk);
    init_heap(&H);
    knn_sub_vptree_search_v2(point, root,&H,metric);
    heap_sort(&H);
	for(size_t i = 0; i < H.count; ++i) H.data[i].value = H.data[i].value*H.data[i].value;
    return H;
}


