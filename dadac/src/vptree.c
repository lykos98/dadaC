#include <stdint.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <stdint.h>
#include  "../include/dadac.h"



/**
 * 
 * KDtree implementation 
 * 
 * 
*/



// Standard Lomuto partition function


// Implementation of QuickSelect


//inline float_t hyper_plane_dist(float_t* p1, float_t* p2, int var);

//inline int hyper_plane_side(float_t* p1, float_t* p2, int var);


void swap_vpnode_ptrs(vpnode **x, vpnode **y) {
    vpnode* tmp;
    tmp = *x;
    *x = *y;
    *y = tmp;
}

/**
 * 
 * KDtree implementation 
 * 
 * 
*/

void initialize_vpnode_array(vpnode* nodeArray, void* data, idx_t n, idx_t bytesPerElement)
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
		
    }

}

void initialize_vpnode_ptrs(vpnode** pointersArray, vpnode* nodeArray, idx_t n)
{
    for(idx_t i = 0; i < n; ++i) pointersArray[i] = nodeArray + i;

}

int cmp_vpnodes(vpnode* point, vpnode* pivot)
{
    float_t res = point->__dist - pivot->__dist;
    return (res > 0);
}

// Standard Lomuto partition function

int partition_vpnodes(vpnode** arr, int low, int high)
{
    vpnode* pivot = arr[high];
    
    int i = (low - 1);
    for (int j = low; j <= high - 1; j++) {
        if (!cmp_vpnodes(arr[j],pivot)) {
            i++;
            swap_vpnode_ptrs(arr + i, arr + j);
        }
    }
    swap_vpnode_ptrs(arr + i + 1, arr + high);
    return (i + 1);

}

int median_of_vpnodes(vpnode** a, int left, int right)
{
    //printf("----------\n");
    int k = left + ((right - left + 1)/2); 
    if(left == right) return left;
    if(left == (right - 1)){
        if(cmp_vpnodes(a[left],a[right])) {swap_vpnode_ptrs(a + left, a + right);}
        return right;
    }
    while (left <= right) {
 
        // Partition a[left..right] around a pivot
        // and find the position of the pivot
        int pivotIndex = partition_vpnodes(a, left, right);
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



vpnode* build_vptree(vpnode** t, int start, int end, vpnode* parent, float_t (*metric)(void*, void*))
{
    vpnode *n = NULL;
    //printf("%d \n",level);

    int median_idx = -1;
    if ((end - start) < 0) return 0;
    if (end  == start) {
        n = t[start];
        n -> parent = parent;
        n -> inside = NULL;
        n -> outside = NULL;
		n -> mu = 0;
        return n;
    }

	int vpIdx = start; 

	//compute distances
	for(int i = start; i <= end; ++i) t[i] -> __dist = metric(t[vpIdx] -> data,t[i] -> data);

	//now swap the vpIdx in the first place to retrieve it later
	//swap_vpnode_ptrs(t + vpIdx, t + start);
	//compute the median, as byproduct obtain the array partitioned on inside and outside BUT, with the median in the median place
    median_idx = median_of_vpnodes(t, start, end);
	float_t mu = t[median_idx]->__dist;
	//now swap start with the median
	
    if(median_idx > -1){
        n = t[vpIdx];
		n->mu = mu;
		#pragma omp parallel
		{
			#pragma omp single
			{
				#pragma omp task shared(n)
				n->inside  = build_vptree(t, start + 1, median_idx, n, metric);
				#pragma omp task shared(n)
				n->outside = build_vptree(t, median_idx + 1, end, n, metric);
				//#pragma omp taskwait
			}

		}
        n->parent = parent;
    }

    return n;

}



#define INSIDE 0
#define OUTSIDE 1

#ifdef ITERATIVE_VPTREE
/* 
 * Experimental feauture, not working so well
 * Iterative implementation of the vp tree not being not very fast
 * or sensibly (?) faster than recursive implementation
 */

const stackNode nodeNULL = { .node = NULL, .side = -1, .mu = 0, .current_distance = 0};

void stackInit(stack_vpnodes* s)
{
	s -> count = 0;
	s -> size  = DEFAULT_STACK_SIZE;
	s -> data  = (stackNode*)malloc(DEFAULT_STACK_SIZE*sizeof(stackNode));
}

void stackPush(stack_vpnodes* s, stackNode n)
{
	if(s -> count < s -> size)	
	{
		size_t idx = s -> count;
		s -> data[idx] = n;
		s -> count ++;
	}
	else
	{
		size_t new_size = s -> size + DEFAULT_STACK_SIZE;
		s -> data = realloc(s -> data, new_size * sizeof(stackNode));
		s -> size = new_size;
		size_t idx = s -> count;
		s -> data[idx] = n;
		s -> count ++;
	}
}

stackNode stackPop(stack_vpnodes* s)
{
	if(s -> count == 0)
	{
		return nodeNULL;
	}
	else
	{
		s -> count--;
		return s -> data[s -> count];
	}
}

void stackReset(stack_vpnodes* s)
{
	s -> count = 0;
	return;
}

void knn_sub_vptree_search_iterative(void* point, vpnode* root, heap * H, stack_vpnodes* s, float_t (*metric)(void*,void*))
{
    //int split_var = kdtree_root -> split_var;
	vpnode* n = root;
    //float_t current_distance = metric(point, root -> data);
    //insert_max_heap(H, current_distance, root -> array_idx);
	stackReset(s);

	while(s -> count > 0 || n != NULL)
	{
		switch(n != NULL)
		{
			case 1:
				{
					float_t current_distance = metric(point, n -> data);
					insert_max_heap(H, current_distance, n -> array_idx);
					int side = (current_distance > n -> mu);
					switch (side)
					{
						case INSIDE:
						{
							if(n -> outside)
							{
								stackNode sn = {.node = n -> outside, .side = side, .mu = n -> mu, .current_distance = current_distance};	
								stackPush(s, sn);
							}
							n = n -> inside;
							break;
						}
						case OUTSIDE:
						{
							if(n -> inside)
							{
								stackNode sn = {.node = n -> inside, .side = side, .mu = n -> mu, .current_distance = current_distance};	
								stackPush(s, sn);
							}
							n  = n -> outside;
							break;
						}
						default:
							break;
					}
					
				}
				break;
			
			case 0:
				{
					stackNode sn = stackPop(s);
					float_t tau = H -> data[0].value;
					//retrieve the maximum distance found so far
					int heapNotFull = (H -> count) < (H -> N);
					switch (sn.side)
					{
						case INSIDE:
							// the node MUST have the CHILD, then or the condition holds or the heap is not full yet 
							//if ((sn.current_distance + tau) > sn.mu  || heapNotFull) 
							//{
							//	n = sn.node;
							//}
							n = ((sn.current_distance + tau) > sn.mu  || heapNotFull) ? sn.node : NULL;
							break;

						case OUTSIDE:
							//if 	( root -> inside && ((current_distance - tau) < root->mu  || heapNotFull)) knn_sub_vptree_search(point, root -> inside, H, metric);
							//if 	(sn.current_distance  < (sn.mu + tau)  || heapNotFull) 
							//{
							//	n = sn.node;
							//}
							n = (sn.current_distance  < (sn.mu + tau)  || heapNotFull) ? sn.node : NULL; 
							break;
						
						default:
							break;
					}
				}
				break;

			default:
				break;
		}
	}
    return;
	

}

#else

void knn_sub_vptree_search(void* point, vpnode* root, heap * H, float_t (*metric)(void*,void*))
{
    float_t current_distance = metric(point, root -> data);
    insert_max_heap(H, current_distance, root -> array_idx);
	#define INSIDE 0
	#define OUTSIDE 1
	int side = (current_distance > root -> mu);

    switch (side)
    {
        case INSIDE:
            if(root -> inside) knn_sub_vptree_search(point, root -> inside, H, metric);
            break;
        case OUTSIDE:
            if(root -> outside) knn_sub_vptree_search(point, root -> outside, H, metric);
            break;

        default:
            break;
    }
    float_t tau = H -> data[0].value;
	//retrieve the maximum distance found so far
	int heapNotFull = (H -> count) < (H -> N);
	switch (side)
	{
		case INSIDE:
			// the node MUST have the CHILD, then or the condition holds or the heap is not full yet 
			if 	( root -> outside && ((current_distance + tau) > root->mu  || heapNotFull)) knn_sub_vptree_search(point, root -> outside, H, metric);
			break;

		case OUTSIDE:
			//if 	( root -> inside && ((current_distance - tau) < root->mu  || heapNotFull)) knn_sub_vptree_search(point, root -> inside, H, metric);
			if 	( root -> inside && (current_distance  < (root->mu + tau)  || heapNotFull)) knn_sub_vptree_search(point, root -> inside, H, metric);
			break;
		
		default:
			break;
	}
    return;
	

}

#endif



#ifdef ITERATIVE_VPTREE 
	heap knn_vptree(void* point, vpnode* root, int maxk, stack_vpnodes* s, float_t (*metric)(void*, void*))
#else
	heap knn_vptree(void* point, vpnode* root, int maxk, float_t (*metric)(void*, void*))
#endif
{
    heap H;
	allocate_heap(&H, maxk);
    init_heap(&H);
	#ifdef ITERATIVE_VPTREE
		knn_sub_vptree_search_iterative(point, root,&H,s,metric);
	#else
		knn_sub_vptree_search(point, root,&H,metric);
	#endif
    heap_sort(&H);
	for(size_t i = 0; i < H.count; ++i) H.data[i].value = H.data[i].value*H.data[i].value;
    return H;
}


