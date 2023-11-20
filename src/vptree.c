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


void swap_vpTreeNode_ptrs(vpTreeNode **x, vpTreeNode **y) {
    vpTreeNode* tmp;
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

void initialize_vpTreeNode_array(vpTreeNode* nodeArray, void* data, idx_t n, idx_t bytesPerElement)
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

void initialize_vpTreeNodes_pointers(vpTreeNode** pointersArray, vpTreeNode* nodeArray, idx_t n)
{
    for(idx_t i = 0; i < n; ++i) pointersArray[i] = nodeArray + i;

}

int cmp_vpTreeNodes(vpTreeNode* point, vpTreeNode* pivot)
{
    float_t res = point->__dist - pivot->__dist;
    return (res > 0);
}

// Standard Lomuto partition function

int partition_vpTreeNodes(vpTreeNode** arr, int low, int high)
{
    vpTreeNode* pivot = arr[high];
    
    int i = (low - 1);
    for (int j = low; j <= high - 1; j++) {
        if (!cmp_vpTreeNodes(arr[j],pivot)) {
            i++;
            swap_vpTreeNode_ptrs(arr + i, arr + j);
        }
    }
    swap_vpTreeNode_ptrs(arr + i + 1, arr + high);
    return (i + 1);

}

int median_of_vpTreeNodes(vpTreeNode** a, int left, int right)
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
        if(cmp_vpTreeNodes(a[left],a[right])) {swap_vpTreeNode_ptrs(a + left, a + right);}
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
        int pivotIndex = partition_vpTreeNodes(a, left, right);
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



vpTreeNode* build_vpTree(vpTreeNode** t, int start, int end, vpTreeNode* parent, float_t (*metric)(void*, void*))
{
    vpTreeNode *n = NULL;
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

	//int vpIdx = start + (end - start + 1)/2; 
	int vpIdx = start; 

	//compute distances
	for(int i = start; i <= end; ++i) t[i] -> __dist = metric(t[vpIdx] -> data,t[i] -> data);
	//for(int i = start; i <= end; ++i) 
	//{
	//	printf("%.5lf ", t[i] -> __dist );
	//}
	//printf("----\n\n");
	//now swap the vpIdx in the first place to retrieve it later
	//swap_vpTreeNode_ptrs(t + vpIdx, t + start);
	//compute the median, as byproduct obtain the array partitioned on inside and outside BUT, with the median in the median place
    median_idx = median_of_vpTreeNodes(t, start, end);
	//printf("Median: %d %.4lf\n\n", median_idx - start, t[median_idx] -> __dist );
	//for(int i = start; i <= end; ++i) 
	//{
	//	printf("%.5lf ", t[i] -> __dist );
	//}
	//printf("*****\n\n");
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
				n->inside  = build_vpTree(t, start + 1, median_idx, n, metric);
				#pragma omp task shared(n)
				n->outside = build_vpTree(t, median_idx + 1, end, n, metric);
				//#pragma omp taskwait
			}

		}
        n->parent = parent;
    }

	/*
    if(median_idx > -1){
        n = t[vpIdx];
		n->mu = mu;
		n->inside  = build_vpTree(t, start + 1, median_idx, n, metric);
		n->outside = build_vpTree(t, median_idx + 1, end, n, metric);
        n->parent = parent;
    }
	*/
    return n;

}



#define INSIDE 0
#define OUTSIDE 1

#ifdef ITERATIVE_VPTREE

const stackNode nodeNULL = { .node = NULL, .side = -1, .mu = 0, .current_distance = 0};

void stackInit(stack_vpTreeNodes* s)
{
	s -> count = 0;
	s -> size  = DEFAULT_STACK_SIZE;
	s -> data  = (stackNode*)malloc(DEFAULT_STACK_SIZE*sizeof(stackNode));
}

void stackPush(stack_vpTreeNodes* s, stackNode n)
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

stackNode stackPop(stack_vpTreeNodes* s)
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

void stackReset(stack_vpTreeNodes* s)
{
	s -> count = 0;
	return;
}

void KNN_sub_vpTree_search_iterative(void* point, vpTreeNode* root, Heap * H, stack_vpTreeNodes* s, float_t (*metric)(void*,void*))
{
    //int split_var = kdtree_root -> split_var;
	vpTreeNode* n = root;
    //float_t current_distance = metric(point, root -> data);
    //insertMaxHeap(H, current_distance, root -> array_idx);
	stackReset(s);

	while(s -> count > 0 || n != NULL)
	{
		switch(n != NULL)
		{
			case 1:
				{
					float_t current_distance = metric(point, n -> data);
					insertMaxHeap(H, current_distance, n -> array_idx);
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
							// the node MUST have the CHILD, then or the condition holds or the Heap is not full yet 
							//if ((sn.current_distance + tau) > sn.mu  || heapNotFull) 
							//{
							//	n = sn.node;
							//}
							n = ((sn.current_distance + tau) > sn.mu  || heapNotFull) ? sn.node : NULL;
							break;

						case OUTSIDE:
							//if 	( root -> inside && ((current_distance - tau) < root->mu  || heapNotFull)) KNN_sub_vpTree_search(point, root -> inside, H, metric);
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

void KNN_sub_vpTree_search(void* point, vpTreeNode* root, Heap * H, float_t (*metric)(void*,void*))
{
    //int split_var = kdtree_root -> split_var;
    float_t current_distance = metric(point, root -> data);
    insertMaxHeap(H, current_distance, root -> array_idx);
	#define INSIDE 0
	#define OUTSIDE 1
	int side = (current_distance > root -> mu);

    switch (side)
    {
        case INSIDE:
            if(root -> inside) KNN_sub_vpTree_search(point, root -> inside, H, metric);
            break;
        case OUTSIDE:
            if(root -> outside) KNN_sub_vpTree_search(point, root -> outside, H, metric);
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
			// the node MUST have the CHILD, then or the condition holds or the Heap is not full yet 
			if 	( root -> outside && ((current_distance + tau) > root->mu  || heapNotFull)) KNN_sub_vpTree_search(point, root -> outside, H, metric);
			break;

		case OUTSIDE:
			//if 	( root -> inside && ((current_distance - tau) < root->mu  || heapNotFull)) KNN_sub_vpTree_search(point, root -> inside, H, metric);
			if 	( root -> inside && (current_distance  < (root->mu + tau)  || heapNotFull)) KNN_sub_vpTree_search(point, root -> inside, H, metric);
			break;
		
		default:
			break;
	}
    return;
	

}

#endif



#ifdef ITERATIVE_VPTREE 
Heap KNN_vpTree(void* point, vpTreeNode* root, int maxk, stack_vpTreeNodes* s, float_t (*metric)(void*, void*))
#else
Heap KNN_vpTree(void* point, vpTreeNode* root, int maxk, float_t (*metric)(void*, void*))
#endif
{
    Heap H;
	allocateHeap(&H, maxk);
    initHeap(&H);
	#ifdef ITERATIVE_VPTREE
		KNN_sub_vpTree_search_iterative(point, root,&H,s,metric);
	#else
		KNN_sub_vpTree_search(point, root,&H,metric);
	#endif
    HeapSort(&H);
	for(size_t i = 0; i < H.count; ++i) H.data[i].value = H.data[i].value*H.data[i].value;
    return H;
}


