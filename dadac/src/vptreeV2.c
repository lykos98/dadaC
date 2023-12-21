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

/**
 * 
 * KDtree implementation 
 * 
 * 
*/


void* swapMem;

void swap_vpTreeNode_ptrs_V2(vpTreeNodeV2 *x, vpTreeNodeV2 *y) {
    vpTreeNodeV2 tmp;
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
	//vpTreeNodeV2 tmpNode = *(*x);
	//*(*x) = *(*y);
	//*(*y)  = tmpNode;
	 
}

/**
 * 
 * KDtree implementation 
 * 
 * 
*/

void initialize_vpTreeNode_array_V2(vpTreeNodeV2* nodeArray, void* data, idx_t n, idx_t bytesPerElement)
{
    for(idx_t i = 0; i < n; ++i)
    {
        nodeArray[i].data = data + (i*bytesPerElement);
        nodeArray[i].array_idx = i;
        nodeArray[i].inside = NULL;
        nodeArray[i].outside = NULL;
        nodeArray[i].parent = NULL;
        nodeArray[i].mu = 0;
        nodeArray[i].__dist = 0;
        nodeArray[i].isLeaf = 0;
        nodeArray[i].nodeList.data = NULL;
        nodeArray[i].nodeList.count = 0;
		nodeArray[i].__bytesize = bytesPerElement;
		
    }

}

void initialize_vpTreeNodes_pointers_V2(vpTreeNodeV2** pointersArray, vpTreeNodeV2* nodeArray, idx_t n)
{
    for(idx_t i = 0; i < n; ++i) pointersArray[i] = nodeArray + i;

}

int cmp_vpTreeNodes_V2(vpTreeNodeV2* point, vpTreeNodeV2* pivot)
{
    float_t res = point->__dist - pivot->__dist;
    return (res > 0);
}

// Standard Lomuto partition function

int partition_vpTreeNodes_V2(vpTreeNodeV2* arr, int low, int high)
{
    vpTreeNodeV2 pivot = arr[high];
    
    int i = (low - 1);
    for (int j = low; j <= high - 1; j++) {
        if (!cmp_vpTreeNodes_V2(arr + j,&pivot)) {
            i++;
            swap_vpTreeNode_ptrs_V2(arr + i, arr + j);
        }
    }
    swap_vpTreeNode_ptrs_V2(arr + i + 1, arr + high);
    return (i + 1);

}

int median_of_vpTreeNodes_V2(vpTreeNodeV2* a, int left, int right)
{
    //printf("----------\n");
    int k = left + ((right - left + 1)/2); 

    if(left == right) return left;
    if(left == (right - 1)){
        if(cmp_vpTreeNodes_V2(a + left,a + right)) {swap_vpTreeNode_ptrs_V2(a + left, a + right);}
        return right;
    }
    while (left <= right) {
 
        // Partition a[left..right] around a pivot
        // and find the position of the pivot
        int pivotIndex = partition_vpTreeNodes_V2(a, left, right);
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


vpTreeNodeV2* build_vpTree_V2(vpTreeNodeV2* t, int start, int end, vpTreeNodeV2* parent, float_t (*metric)(void*, void*))
{
	#ifdef SWMEM
		if(parent == NULL)
		{
			swapMem = malloc((t + start) ->__bytesize);
		}
	#endif

    vpTreeNodeV2 *n = NULL;
    //printf("%d \n",level);

    int median_idx = -1;
    if ((end - start) < 0) return NULL;

	if(end - start < DEFAULT_LEAF_SIZE)
	//if(end - start < -1)
	{
		n =  t + start;
		n -> isLeaf = 1;
        n -> parent = parent;
        n -> inside = NULL;
        n -> outside = NULL;
		n -> mu = 0;
		size_t j = 0;
		n -> nodeList.count = (size_t)(end - start + 1);
		#ifdef VOPT
			n -> nodeList.indexes = (idx_t*)malloc(n -> nodeList.count * sizeof(idx_t));
			n -> nodeList.start_ptr = t[start].data;
			n -> nodeList.end_ptr 	= t[end].data;
		#else
			n -> nodeList.data = (vpTreeNodeV2**)malloc(n -> nodeList.count * sizeof(vpTreeNodeV2*));
		#endif
		for(int i = start; i <= end; ++i){
			t[i].parent = n;
			t[i].isLeaf = 1;
			t[i].inside = NULL;
			t[i].outside = NULL;
			#ifdef VOPT
				n -> nodeList.indexes[j] = t[i].array_idx;
			#else
				n -> nodeList.data[j] = t + i;
			#endif

			++j;
		}
		return n;
		
	}

	//int vpIdx = start + (end - start + 1)/2; 
	int vpIdx = start; 

	//compute distances
	for(int i = start; i <= end; ++i) t[i].__dist = metric(t[vpIdx].data,t[i].data);

	//for(int i = start; i <= end; ++i) 
	//{
	//	printf("%.5lf ", t[i] -> __dist );
	//}
	//printf("----\n\n");
	
	//now swap the vpIdx in the first place to retrieve it later
	//swap_vpTreeNode_ptrs(t + vpIdx, t + start);
	//compute the median, as byproduct obtain the array partitioned on inside and outside BUT, with the median in the median place
    median_idx = median_of_vpTreeNodes_V2(t, start, end);
	//printf("Median: %d %.4lf\n\n", median_idx - start, t[median_idx]. __dist );
	//for(int i = start; i <= end; ++i) 
	//{
	//	printf("%.5lf ", t[i].__dist );
	//}
	//printf("*****\n\n");
	float_t mu = t[median_idx].__dist;
	//now swap start with the median
	
    if(median_idx > -1){
        n = t + vpIdx;
		n->mu = mu;
		n->inside  = build_vpTree_V2(t, start + 1, median_idx, n, metric);
		n->outside = build_vpTree_V2(t, median_idx + 1, end, n, metric);
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


void KNN_sub_vpTree_search_V2(void* point, vpTreeNodeV2* root, Heap * H, float_t (*metric)(void*,void*))
{
	
	#ifdef VOPT	
	if(root -> isLeaf)
	{
		
		for(idx_t i = 0; i < root -> nodeList.count; ++i)
		{

			idx_t j = root -> nodeList.indexes[i];
			float_t distance = metric(point, root -> nodeList.start_ptr + (i*(root -> __bytesize)));
			insertMaxHeap(H, distance,j);
		}
		return;
	}
	#else
		if(root -> isLeaf)
		{
			for(size_t i = 0; i < root -> nodeList.count; ++i)
			{
				__builtin_prefetch(root -> nodeList.data + i + 1, 0, 3);
				vpTreeNodeV2* n = root -> nodeList.data[i];
				float_t distance = metric(point, n -> data);
				insertMaxHeap(H, distance,n -> array_idx);
			}
			return;
		}
	#endif


    float_t current_distance = metric(point, root -> data);
    insertMaxHeap(H, current_distance, root -> array_idx);
	#define INSIDE 0
	#define OUTSIDE 1
	int side = (current_distance > root -> mu);

    switch (side)
    {
        case INSIDE:
            KNN_sub_vpTree_search_V2(point, root -> inside, H, metric);
            //if(root -> inside) KNN_sub_vpTree_search_V2(point, root -> inside, H, metric);
            break;
        case OUTSIDE:
            KNN_sub_vpTree_search_V2(point, root -> outside, H, metric);
            //if(root -> outside) KNN_sub_vpTree_search_V2(point, root -> outside, H, metric);
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
			//if 	( ((current_distance + tau) > root->mu  || heapNotFull)) KNN_sub_vpTree_search_V2(point, root -> outside, H, metric);
			if 	( (current_distance + tau) > root->mu  ) KNN_sub_vpTree_search_V2(point, root -> outside, H, metric);
			//if 	( root -> outside && ((current_distance + tau) > root->mu  || heapNotFull)) KNN_sub_vpTree_search_V2(point, root -> outside, H, metric);
			break;

		case OUTSIDE:
			//if 	( root -> inside && ((current_distance - tau) < root->mu  || heapNotFull)) KNN_sub_vpTree_search(point, root -> inside, H, metric);
			//if 	( (current_distance  < (root->mu + tau)  || heapNotFull)) KNN_sub_vpTree_search_V2(point, root -> inside, H, metric);
			if 	( (current_distance  < (root->mu + tau) )) KNN_sub_vpTree_search_V2(point, root -> inside, H, metric);
			//if 	( root -> inside && (current_distance  < (root->mu + tau)  || heapNotFull)) KNN_sub_vpTree_search_V2(point, root -> inside, H, metric);
			break;
		
		default:
			break;
	}
    return;
	

}

Heap KNN_vpTree_V2(void* point, vpTreeNodeV2* root, int maxk, float_t (*metric)(void*, void*))
{
    Heap H;
    allocateHeap(&H,maxk);
    initHeap(&H);
    KNN_sub_vpTree_search_V2(point, root,&H,metric);
    HeapSort(&H);
	for(size_t i = 0; i < H.count; ++i) H.data[i].value = H.data[i].value*H.data[i].value;
    return H;
}


