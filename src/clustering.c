//#include "../include/read_fof_snapshot.h"
#include "../include/clustering.h"
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#define MAX_SERIAL_MERGING 40000
#define MAX_N_NGBH 1000
#define PREALLOC_BORDERS 10

extern unsigned int data_dims;
idx_t Npart;
const border_t border_null = {.density = -1.0, .error = 0, .idx = NOBORDER};
const SparseBorder_t SparseBorder_null = {.density = -1.0, .error = 0, .idx = NOBORDER, .i = NOBORDER, .j = NOBORDER};

void LinkedList_Insert(LinkedList* L, Node* n)
{
    ++(L -> count);
    n -> next = L -> head;
    L -> head = n;
}

/*****************************
 * Clusters object functions *
 *****************************/

void dummy_Clusters_allocate(Clusters * c, int s)
{
    /*************************************
     * allocate additional resources and *
     * pointers for Clusters object      *
     *************************************/
    if(c -> centers.count == 0)
    {
        printf("Provide a valid cluster centers list\n");
        return;
    }

    idx_t nclus = c -> centers.count;
    
    if(s)
    {
	    //printf("Using sparse implementation\n");
	    c -> UseSparseBorders = 1;
	    c -> SparseBorders = (AdjList_t*)malloc(nclus*sizeof(AdjList_t));
	    for(idx_t i = 0; i < nclus; ++i)
	    {
		    c -> SparseBorders[i].count = 0;
		    c -> SparseBorders[i].size  = PREALLOC_BORDERS;
		    c -> SparseBorders[i].data  = (SparseBorder_t*)malloc(PREALLOC_BORDERS*sizeof(SparseBorder_t));
	    }

    }
    else
    {
	    //printf("Using dense implementation\n");
	    c -> UseSparseBorders = 0;
	    c -> __borders_data         = (border_t*)malloc(nclus*nclus*sizeof(border_t)); 
	    c -> borders                = (border_t**)malloc(nclus*sizeof(border_t*));

	    #pragma omp parallel for
	    for(idx_t i = 0; i < nclus; ++i)
	    {
		c -> borders[i]         = c -> __borders_data + i*nclus;
		for(idx_t j = 0; j < nclus; ++j)
		{
		    c -> borders[i][j] = border_null;
		}
	    }
    }
}

void Clusters_allocate(Clusters * c)
{
    /*************************************
     * allocate additional resources and *
     * pointers for Clusters object      *
     *************************************/
    if(c -> centers.count == 0)
    {
        printf("Provide a valid cluster centers list\n");
        return;
    }

    idx_t nclus = c -> centers.count;
    
    if(nclus > 10)
    {
	    c -> UseSparseBorders = 1;
	    c -> SparseBorders = (AdjList_t*)malloc(nclus*sizeof(AdjList_t));
	    for(idx_t i = 0; i < nclus; ++i)
	    {
		    c -> SparseBorders[i].count = 0;
		    c -> SparseBorders[i].size  = PREALLOC_BORDERS;
		    c -> SparseBorders[i].data  = (SparseBorder_t*)malloc(PREALLOC_BORDERS*sizeof(SparseBorder_t));
	    }

    }
    else
    {
	    c -> UseSparseBorders = 0;
	    c -> __borders_data         = (border_t*)malloc(nclus*nclus*sizeof(border_t)); 
	    c -> borders                = (border_t**)malloc(nclus*sizeof(border_t*));

	    #pragma omp parallel for
	    for(idx_t i = 0; i < nclus; ++i)
	    {
		c -> borders[i]         = c -> __borders_data + i*nclus;
		for(idx_t j = 0; j < nclus; ++j)
		{
		    c -> borders[i][j] = border_null;
		}
	    }
    }
}

void AdjList_Insert(AdjList_t* l, SparseBorder_t b)
{
	if(l -> count < l -> size)
	{
		l -> data[l -> count] = b;
		l -> count++;
	}
	else
	{
		l -> size += PREALLOC_BORDERS; 
		l -> data = realloc( l -> data, sizeof(SparseBorder_t) * ( l -> size));
		l -> data[l -> count] = b;
		l -> count++;
	}
}

void AdjList_reset(AdjList_t* l)
{
	free(l -> data);
	l -> count = 0;
	l -> size  = 0;
	l -> data  = NULL;
}

void Clusters_Reset(Clusters * c)
{
	if(c -> UseSparseBorders)
	{
		for(idx_t i = 0; i < c -> centers.count; ++i)
		{
			AdjList_reset((c -> SparseBorders) + i);
		
		}
		free(c -> SparseBorders);
		c -> SparseBorders = NULL;
	}
	else
	{
		free(c -> __borders_data);
		free(c -> borders);
	}
    free(c -> centers.data);
}

void Clusters_free(Clusters * c)
{

    Clusters_Reset(c);
}


void SparseBorder_Insert(Clusters *c, SparseBorder_t b)
{
	idx_t i = b.i;
	AdjList_t l = c -> SparseBorders[i];
	int check = 1;
	for(idx_t k = 0; k < l.count; ++k)
	{
		SparseBorder_t p = l.data[k];
		if(p.i == b.i && p.j == b.j)
		{
			if( b.density > p.density)
			{
				l.data[k] = b;
			}
			check = 0;
		}
	}
	if(check) AdjList_Insert(c -> SparseBorders + i, b);
	return;
}

SparseBorder_t SparseBorder_get(Clusters* c, idx_t i, idx_t j)
{
	SparseBorder_t b = SparseBorder_null;
	AdjList_t l = c -> SparseBorders[i];
	for(idx_t el = 0; el < l.count; ++el)
	{
		SparseBorder_t candidate = l.data[el];
		if(candidate.i == i && candidate.j == j)
		{
			b = candidate;
		}
	}

	return b;
}

/*****************
 * Dyanmic Array *
 *****************/

void DynamicArray_allocate(lu_dynamicArray * a)
{
    a -> data = (idx_t*)malloc(ARRAY_INCREMENT*sizeof(idx_t));
    a -> count = 0;
    a -> size = ARRAY_INCREMENT;
}

void DynamicArray_pushBack(lu_dynamicArray * a, idx_t p)
{
    if(a -> count < a -> size)
    {
        a -> data[a -> count] =  p;
        a -> count += 1;
    }
    else{
        a -> size += ARRAY_INCREMENT;
        a -> data = realloc(a -> data, a -> size * sizeof(idx_t));
        a -> data[a -> count] =  p;
        a -> count += 1;
    }
}

void DynamicArray_Reset(lu_dynamicArray * a){
    a -> count = 0;
}

void DynamicArray_Reserve(lu_dynamicArray * a, idx_t n)
{
    a -> data = realloc(a -> data, n*sizeof(idx_t));
    a -> size = n;
}

void DynamicArray_Init(lu_dynamicArray * a)
{
    a -> data = NULL;
    a -> count = 0;
    a -> size = 0;
}


/*******************
 * Clustering part *
 *******************/

void KNN_search(Datapoint_info * particles, FLOAT_TYPE * data, kd_node* root, idx_t n, idx_t k)
{
    struct timespec start_tot, finish_tot;
    double elapsed_tot;
    printf("KNN search:\n");
    clock_gettime(CLOCK_MONOTONIC, &start_tot);


    #pragma omp parallel for
    for(int p = 0; p < n; ++p)
    {
        particles[p].ngbh = KNN(data + data_dims*p, root, k);
        particles[p].array_idx = p;
    }

    clock_gettime(CLOCK_MONOTONIC, &finish_tot);
    elapsed_tot = (finish_tot.tv_sec - start_tot.tv_sec);
    elapsed_tot += (finish_tot.tv_nsec - start_tot.tv_nsec) / 1000000000.0;
    printf("\tTotal time: %.3lfs\n\n", elapsed_tot);
    return;

}

int cmp(const void * a, const void * b){
    FLOAT_TYPE aa = *((FLOAT_TYPE*)a);
    FLOAT_TYPE bb = *((FLOAT_TYPE*)b);
    return 2*(aa > bb ) - 1; 
}



FLOAT_TYPE avg(const FLOAT_TYPE * x, const idx_t n)
{
    FLOAT_TYPE f = 0;
    for(idx_t i = 0; i < n; ++i)
    {
        f += x[i];
    }
    return f/(FLOAT_TYPE)n;
}



FLOAT_TYPE mEst2(FLOAT_TYPE * x, FLOAT_TYPE *y, idx_t n)
{

    /********************************************
     * Estimate the m coefficient of a straight *
     * line passing through the origin          *
     * params:                                  *
     * - x: x values of the points              *
     * - y: y values of the points              *
     * - n: size of the arrays                  *
     ********************************************/
     

    FLOAT_TYPE x_avg, y_avg;
    FLOAT_TYPE num = 0;
    FLOAT_TYPE den = 0;
    FLOAT_TYPE dd;
    for(idx_t i = 0; i < n; ++i)
    {
        FLOAT_TYPE xx = x[i];
        FLOAT_TYPE yy = y[i];

        dd = xx;
        num += dd*yy;
        den += dd*dd;

    }
  
    return num/den;
}
FLOAT_TYPE mEst(FLOAT_TYPE * x, FLOAT_TYPE *y, idx_t n)
{
    FLOAT_TYPE x_avg, y_avg;
    x_avg = avg(x,n);
    y_avg = avg(y,n);
    FLOAT_TYPE num = 0;
    FLOAT_TYPE den = 0;
    FLOAT_TYPE dd;
    for(idx_t i = 0; i < n - 1; ++i)
    {
        FLOAT_TYPE xx = x[i];
        FLOAT_TYPE yy = y[i];

        dd = (xx - x_avg);
        num += dd*(yy - y_avg);
        den += dd*dd;

    }
  
    return num/den;
}

FLOAT_TYPE idEstimate(Datapoint_info* particles, idx_t n)
{

    /*********************************************************************************************
     * Estimation of the intrinsic dimension of a dataset                                        *
     * args:                                                                                     *
     * - particles: array of structs                                                             *
     * - n: number of particles                                                                  *
     * Estimates the id via 2NN method. Computation of the log ratio of the                      *
     * distances of the first 2 neighbors of each point. Then compute the empirical distribution *
     * of these log ratii                                                                        *
     * The intrinsic dimension is found by fitting with a straight line passing through the      *
     * origin                                                                                    *
     *********************************************************************************************/

    struct timespec start_tot, finish_tot;
    double elapsed_tot;

    printf("ID estimation:\n");
    clock_gettime(CLOCK_MONOTONIC, &start_tot);

    FLOAT_TYPE fraction = 0.9;
    FLOAT_TYPE* r = (FLOAT_TYPE*)malloc(n*sizeof(FLOAT_TYPE));
    FLOAT_TYPE* Pemp = (FLOAT_TYPE*)malloc(n*sizeof(FLOAT_TYPE));

    for(int i = 0; i < n; ++i)
    {
        r[i] = 0.5 * log(particles[i].ngbh.data[2].value/particles[i].ngbh.data[1].value);
        Pemp[i] = -log(1 - (FLOAT_TYPE)(i + 1)/(FLOAT_TYPE)n);
    }
    qsort(r,n,sizeof(FLOAT_TYPE),cmp);

    idx_t Neff = (idx_t)(n*fraction);

    FLOAT_TYPE d = mEst2(r,Pemp,Neff); 
    free(r);
    free(Pemp);

    clock_gettime(CLOCK_MONOTONIC, &finish_tot);
    elapsed_tot = (finish_tot.tv_sec - start_tot.tv_sec);
    elapsed_tot += (finish_tot.tv_nsec - start_tot.tv_nsec) / 1000000000.0;
    printf("\tID value: %.6lf\n", d);
    printf("\tTotal time: %.3lfs\n\n", elapsed_tot);

    return d;

}

//void computeRhoOpt(Datapoint_info* particles, kd_node* root, FLOAT_TYPE* data,  const FLOAT_TYPE d, const idx_t points){
//
//    /****************************************************
//     * Point density computation:                       *
//     * args:                                            *
//     * -   paricles: array of structs                   *
//     * -   d       : intrinsic dimension of the dataset *
//     * -   points  : number of points in the dataset    *
//     ****************************************************/
//
//    struct timespec start_tot, finish_tot;
//    double elapsed_tot;
//
//    printf("Density and k* estimation:\n");
//    clock_gettime(CLOCK_MONOTONIC, &start_tot);
//
//    idx_t kMAX = particles[0].ngbh.N;   
//
//    FLOAT_TYPE omega = 0.;  
//    if(sizeof(FLOAT_TYPE) == sizeof(float)){ omega = powf(PI_F,d/2)/tgammaf(d/2.0f + 1.0f);}  
//    else{omega = pow(M_PI,d/2.)/tgamma(d/2.0 + 1.0);}
//
//    //printf("Omega d %f\n", omega);
//
//    #pragma omp parallel for
//    for(idx_t i = 0; i < points; ++i)
//    {
//
//        idx_t j = 4;
//        idx_t k;
//        FLOAT_TYPE dL = 0.0;
//        FLOAT_TYPE vvi, vvj, vp;
//        int cont = 0;
//        while(1)
//        {
//            idx_t kMAX = particles[i].ngbh.N;   
//            while(j < kMAX && dL < DTHR)
//            {
//                idx_t ksel = j - 1;
//                vvi = omega * pow(particles[i].ngbh.data[ksel].value,d/2.);
//                idx_t jj = particles[i].ngbh.data[j].array_idx;
//                vvj = omega * pow(particles[jj].ngbh.data[ksel].value,d/2.);
//                vp = (vvi + vvj)*(vvi + vvj);
//                dL = -2.0 * ksel * log(4.*vvi*vvj/vp);
//                j = j + 1;
//            }
//            if(j == MAX_N_NGBH)
//            {
//                k = j - 1;
//                break;
//
//            }
//            else if(j == kMAX)
//            {
//                //reset ngbh and recalculate using a new size
//                j = j - 1;
//                freeHeap(&particles[i].ngbh);
//                kMAX = 1.2*particles[i].ngbh.N;
//                kMAX = kMAX > MAX_N_NGBH ? MAX_N_NGBH : kMAX;
//                particles[i].ngbh = KNN(data + i*data_dims, root, kMAX );
//            }
//            else
//            {
//                k = j - 2;
//                break;
//            }
//        }
//
//        particles[i].kstar = k;
//        particles[i].log_rho = log((FLOAT_TYPE)(k)/vvi/((FLOAT_TYPE)(points)));
//        //particles[i].log_rho = log((FLOAT_TYPE)(k)) - log(vvi) -log((FLOAT_TYPE)(points));
//        particles[i].log_rho_err =   1.0/sqrt((FLOAT_TYPE)k); //(FLOAT_TYPE)(-Q_rsqrt((float)k));
//        particles[i].g = particles[i].log_rho - particles[i].log_rho_err;
//    }
//
//    clock_gettime(CLOCK_MONOTONIC, &finish_tot);
//    elapsed_tot = (finish_tot.tv_sec - start_tot.tv_sec);
//    elapsed_tot += (finish_tot.tv_nsec - start_tot.tv_nsec) / 1000000000.0;
//    printf("\tTotal time: %.3lfs\n\n", elapsed_tot);
//
//    return;
//
//
//}

void computeRho(Datapoint_info* particles, const FLOAT_TYPE d, const idx_t points){

    /****************************************************
     * Point density computation:                       *
     * args:                                            *
     * -   paricles: array of structs                   *
     * -   d       : intrinsic dimension of the dataset *
     * -   points  : number of points in the dataset    *
     ****************************************************/

    struct timespec start_tot, finish_tot;
    double elapsed_tot;

    printf("Density and k* estimation:\n");
    clock_gettime(CLOCK_MONOTONIC, &start_tot);

    idx_t kMAX = particles[0].ngbh.N;   

    FLOAT_TYPE omega = 0.;  
    if(sizeof(FLOAT_TYPE) == sizeof(float)){ omega = powf(PI_F,d/2)/tgammaf(d/2.0f + 1.0f);}  
    else{omega = pow(M_PI,d/2.)/tgamma(d/2.0 + 1.0);}

    //printf("Omega d %f\n", omega);

    #pragma omp parallel for
    for(idx_t i = 0; i < points; ++i)
    {

        idx_t j = 4;
        idx_t k;
        FLOAT_TYPE dL = 0.0;
        FLOAT_TYPE vvi, vvj, vp;
        while(j < kMAX && dL < DTHR)
        {
            idx_t ksel = j - 1;
            vvi = omega * pow(particles[i].ngbh.data[ksel].value,d/2.);
            idx_t jj = particles[i].ngbh.data[j].array_idx;
            vvj = omega * pow(particles[jj].ngbh.data[ksel].value,d/2.);
            vp = (vvi + vvj)*(vvi + vvj);
            dL = -2.0 * ksel * log(4.*vvi*vvj/vp);
            j = j + 1;
        }
        if(j == kMAX)
        {
            k = j - 1;
        }
        else
        {
            k = j - 2;
        }
        particles[i].kstar = k;
        particles[i].log_rho = log((FLOAT_TYPE)(k)/vvi/((FLOAT_TYPE)(points)));
        //particles[i].log_rho = log((FLOAT_TYPE)(k)) - log(vvi) -log((FLOAT_TYPE)(points));
        particles[i].log_rho_err =   1.0/sqrt((FLOAT_TYPE)k); //(FLOAT_TYPE)(-Q_rsqrt((float)k));
        particles[i].g = particles[i].log_rho - particles[i].log_rho_err;
    }

    clock_gettime(CLOCK_MONOTONIC, &finish_tot);
    elapsed_tot = (finish_tot.tv_sec - start_tot.tv_sec);
    elapsed_tot += (finish_tot.tv_nsec - start_tot.tv_nsec) / 1000000000.0;
    printf("\tTotal time: %.3lfs\n\n", elapsed_tot);

    return;


}

int cmpPP(const void* p1, const void *p2)
{
    /***********************************************
     * Utility function to perform quicksort then  *
     * when clustering assignment is performed     *
     ***********************************************/
    Datapoint_info* pp1 = *(Datapoint_info**)p1;
    Datapoint_info* pp2 = *(Datapoint_info**)p2;
    return 2*(pp1 -> g < pp2 -> g) - 1;
}
void computeCorrection(Datapoint_info* particles, idx_t n, FLOAT_TYPE Z)
{
    /*****************************************************************************
     * Utility function, find the minimum value of the density of the datapoints *
     * and shift them up in order to further work with values greater than 0     *
     *****************************************************************************/
    FLOAT_TYPE min_log_rho = 999999.9;
    

    #pragma omp parallel
    {
        FLOAT_TYPE thread_min_log_rho = 9999999.;
        #pragma omp for
        for(idx_t i = 0; i < n; ++i)
        {
            FLOAT_TYPE tmp = particles[i].log_rho - Z*particles[i].log_rho_err;
            if(tmp < thread_min_log_rho){
                thread_min_log_rho = tmp;
            }
        }
        #pragma omp critical
        if(thread_min_log_rho < min_log_rho) min_log_rho = thread_min_log_rho;

        #pragma omp barrier 
        #pragma omp for
        for(idx_t i = 0; i < n; ++i)
        {
            particles[i].log_rho_c = particles[i].log_rho - min_log_rho + 1;
            particles[i].g = particles[i].log_rho_c - particles[i].log_rho_err;
        }
    }
    //printf("%lf\n",min_log_rho);
}

Clusters Heuristic1(Datapoint_info* particles, FLOAT_TYPE* data, idx_t n)
{
    /**************************************************************
     * Heurisitc 1, from paper of Errico, Facco, Laio & Rodriguez *
     * ( https://doi.org/10.1016/j.ins.2021.01.010 )              *
     *                                                            *
     * args:                                                      *
     * - particles: array of Datapoint structures                 *
     * - data: pointer to the dataset                             *
     * - n: number of Datapoints                                  *
     **************************************************************/

    struct timespec start_tot, finish_tot;
    double elapsed_tot;

    printf("H1: Preliminary cluster assignment\n");
    clock_gettime(CLOCK_MONOTONIC, &start_tot);

    idx_t ncenters = 0;
    idx_t putativeCenters = n;
    lu_dynamicArray allCenters, removedCenters, actualCenters, max_rho;
    DynamicArray_allocate(&allCenters);
    DynamicArray_allocate(&removedCenters);
    DynamicArray_allocate(&actualCenters);
    DynamicArray_allocate(&max_rho);

    Datapoint_info** particles_ptrs = (Datapoint_info**)malloc(n*sizeof(Datapoint_info*));

    struct timespec start, finish;
    double elapsed;


    #ifdef VERBOSE
        clock_gettime(CLOCK_MONOTONIC, &start);
    #endif

    for(idx_t i = 0; i < n; ++i)
    {   
        /*

        Find the centers of the clusters as the points of higher density in their neighborhoods
        A point is tagged as a putative center if it is the point of higer density of its neighborhood 
        
        */

        particles_ptrs[i] = particles + i;
        idx_t maxk = particles[i].kstar + 1;
        FLOAT_TYPE gi = particles[i].g;
        particles[i].is_center = 1;
        particles[i].cluster_idx = -1;
        //printf("%lf\n",p -> g);
        Heap i_ngbh = particles[i].ngbh;
        for(idx_t k = 1; k < maxk; ++k)
        {
            idx_t ngbh_index = i_ngbh.data[k].array_idx;
            FLOAT_TYPE gj = particles[ngbh_index].g;
            if(gj > gi){
                particles[i].is_center = 0;
                break;
            }
        }
        if(particles[i].is_center){
                DynamicArray_pushBack(&allCenters, i);
        }


    }

    #ifdef VERBOSE
        clock_gettime(CLOCK_MONOTONIC, &finish);
        elapsed = (finish.tv_sec - start.tv_sec);
        elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
        printf("\tFinding putative centers: %.3lfs\n",elapsed);
        clock_gettime(CLOCK_MONOTONIC, &start);
    #endif

    idx_t * to_remove = (idx_t*)malloc(allCenters.count*sizeof(idx_t));
    for(idx_t c = 0; c < allCenters.count; ++c) {to_remove[c] = MY_SIZE_MAX;}

    qsort(particles_ptrs, n, sizeof(Datapoint_info*), cmpPP);
    /*****************************************************************************************************
     * /!\ This part is VERY time consuming, complexity depends on the number of center previously found *
     * /!\ It is actually faster when using only a few threads                                           *
     *****************************************************************************************************/
    //#pragma omp parallel for num_threads(4)
    
    //#pragma omp parallel for 
    //for(idx_t p = 0; p < allCenters.count; ++p)
    //{   
    //    /*

    //    Check if the center spotted in the previous part belongs to the neighborhood
    //    of a point of higher density. If this is true remove it from the actual centers

    //    */
    //    idx_t i = allCenters.data[p];
    //    int e = 0;
    //    FLOAT_TYPE gi = particles[i].g;
    //    idx_t i_arrIdx = particles[i].array_idx;
    //    idx_t mr = MY_SIZE_MAX;
    //    FLOAT_TYPE max_g = -99999.0;
    //    for(idx_t j = 0; j < n; ++j)
    //    {
    //        //retrive the particle pointed by the pointer
    //        Datapoint_info pp = *(particles_ptrs[j]);
    //        idx_t kMAXj = pp.kstar;
    //        Heap j_ngbh = pp.ngbh;
    //        FLOAT_TYPE gj = pp.g;
    //        //e = 0;
    //        //check if there is point in which the point i is a neighbor with grater g
    //            //if gj > gi check the neighborhood
    //        //preliminarity check, if i is more distant than k* neighbor break
    //        FLOAT_TYPE dk = j_ngbh.data[kMAXj + 1].value;
    //        FLOAT_TYPE di = euclidean_distance(data + (i*data_dims), data + (pp.array_idx*data_dims));
    //        int stop = 0;
    //        if(dk > di)
    //        {
    //            mr = pp.array_idx;
    //            //found a neighborhood with higher g, break
    //            stop = 1;
    //            //max_g = gj;
    //        } 
    //        if(stop == 1|| i == pp.array_idx)
    //        {
    //            break;
    //        }
    //        
    //    }
    //    to_remove[p] = mr;
    //}
    //
    

    #pragma omp parallel
    {
            
        idx_t * to_remove_private = (idx_t*)malloc(allCenters.count*sizeof(idx_t));
    	for(idx_t c = 0; c < allCenters.count; ++c) {to_remove_private[c] = MY_SIZE_MAX;}
        #pragma omp for
        for(idx_t p = 0; p < n; ++p)
        {
        	Datapoint_info pp = *(particles_ptrs[p]);
        	for(idx_t j = 1; j < pp.kstar + 1; ++j)
        	{
        		idx_t jidx = pp.ngbh.data[j].array_idx;
        		if(particles[jidx].is_center && pp.g > particles[jidx].g)
        		{
        			//particles[jidx].is_center = 0;
        			for(idx_t c = 0; c < allCenters.count; ++c){
        				if(allCenters.data[c] == jidx)
					{

						if(to_remove_private[c] != MY_SIZE_MAX)
						{
							to_remove_private[c] = pp.g > 	particles[to_remove_private[c]].g  ? 
											pp.array_idx : to_remove_private[c];
						}
						else
						{
							to_remove_private[c] = pp.array_idx;
						}
					}
        			}
        		}
        	}
        }
        #pragma omp critical
        {
        	for(idx_t c = 0; c < allCenters.count; ++c)
        	{
        		if(to_remove_private[c] != MY_SIZE_MAX)
			{
				if(to_remove[c] != MY_SIZE_MAX)
				{
					to_remove[c] = particles[to_remove_private[c]].g > particles[to_remove[c]].g ?
						       to_remove_private[c] : to_remove[c];
				}
				else
				{
					to_remove[c] = to_remove_private[c];
				}
			}
        	}
        }

            free(to_remove_private);
    }

	

    for(idx_t p = 0; p < allCenters.count; ++p)
    {
        idx_t i = allCenters.data[p];
        int e = 0;
        FLOAT_TYPE gi = particles[i].g;
        idx_t mr = to_remove[p];
        if(mr != MY_SIZE_MAX)
        {
            //if(particles[mr].g > gi) e = 1;
	    e = 1;
        }
        switch (e)
        {
            case 1:
                {
                    DynamicArray_pushBack(&removedCenters,i);
                    particles[i].is_center = 0;
                    //for(idx_t c = 0; c < removedCenters.count - 1; ++c)
                    //{
                    //    if(mr == removedCenters.data[c])
                    //    {
                    //        mr = max_rho.data[c];
                    //    }
                    //}
                    DynamicArray_pushBack(&max_rho,mr);
                    
                }
                break;
            case 0:
                {
                    DynamicArray_pushBack(&actualCenters,i);
                    particles[i].cluster_idx = actualCenters.count - 1;
                }
                break;
            default:
                break;
        }
    }

    #ifdef VERBOSE
        clock_gettime(CLOCK_MONOTONIC, &finish);
        elapsed = (finish.tv_sec - start.tv_sec);
        elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
        printf("\tFinding actual centers:   %.3lfs\n",elapsed);

        clock_gettime(CLOCK_MONOTONIC, &start);
    #endif

    free(to_remove);

    idx_t nclusters = 0;


    /*****************************************************************************
     * Sort all the particles based on g and then perform the cluster assignment *
     * in asceding order                                                         *
     * UPDATE: particles already sorted                                          *
     *****************************************************************************/
                                                                                

    //qsort(particles_ptrs, n, sizeof(Datapoint_info*), cmpPP);

    for(idx_t i = 0; i < n; ++i)
    {   
        Datapoint_info* p = particles_ptrs[i];
        idx_t ele = p -> array_idx;
        //fprintf(f,"%lu\n",ele);
        if(!(p -> is_center))
        {
            int cluster = -1;
            idx_t k = 0;
            idx_t p_idx;
            idx_t max_k = p -> ngbh.N;
            //assign each particle at the same cluster as the nearest particle of higher density
            while( k < max_k - 1 && cluster == -1)
            {
                ++k;
                p_idx = p -> ngbh.data[k].array_idx;
                cluster = particles[p_idx].cluster_idx; 
            }

            //
            if(cluster == -1)
            {
                FLOAT_TYPE gmax = -99999.;               
                idx_t gm_index;
                for(idx_t k = 0; k < max_k; ++k)
                {
                    idx_t ngbh_index = p -> ngbh.data[k].array_idx;
                    for(idx_t m = 0; m < removedCenters.count; ++m)
                    {
                        FLOAT_TYPE gcand = particles[max_rho.data[m]].g;
                        if(ngbh_index == removedCenters.data[m] && gcand > gmax)
                        {   
                            //printf("%lu -- %lu\n", ele, m);
                            gmax = gcand;
                            gm_index = max_rho.data[m];
                        }
                    }
                }

                cluster = particles[gm_index].cluster_idx;

            }
            p -> cluster_idx = cluster;

        }
    }

    #ifdef VERBOSE
        clock_gettime(CLOCK_MONOTONIC, &finish);
        elapsed = (finish.tv_sec - start.tv_sec);
        elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
        printf("\tTentative clustering:     %.3lfs\n",elapsed);

        clock_gettime(CLOCK_MONOTONIC, &start);
    #endif

    free(particles_ptrs);
    free(max_rho.data);
    free(removedCenters.data);
    free(allCenters.data);


    Clusters c_all;
    c_all.centers = actualCenters;


    #ifdef VERBOSE
        clock_gettime(CLOCK_MONOTONIC, &finish);
        elapsed = (finish.tv_sec - start.tv_sec);
        elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
        printf("\tFinalizing clustering:    %.3lfs\n",elapsed);
        printf("\n");
    #endif

    clock_gettime(CLOCK_MONOTONIC, &finish_tot);
    elapsed_tot = (finish_tot.tv_sec - start_tot.tv_sec);
    elapsed_tot += (finish_tot.tv_nsec - start_tot.tv_nsec) / 1000000000.0;


    printf("\tFound %lu clusters\n",(uint64_t)actualCenters.count);
    printf("\tTotal time: %.3lfs\n\n", elapsed_tot);

    c_all.n = n;
    return c_all;
}

void Heuristic2(Clusters* cluster, Datapoint_info* particles)
{

    #define borders cluster->borders

    struct timespec start_tot, finish_tot;
    double elapsed_tot;
    idx_t n = cluster -> n;

    printf("H2: Finding border points\n");
    clock_gettime(CLOCK_MONOTONIC, &start_tot);


    idx_t nclus = cluster->centers.count; 
    idx_t max_k = particles[0].ngbh.N;

    for(idx_t i = 0; i < n; ++i)
    {
            idx_t pp = NOBORDER;
            /*loop over n neighbors*/
            int c = particles[i].cluster_idx;
            if(!particles[i].is_center)
            {
                for(idx_t k = 1; k < particles[i].kstar + 1; ++k)
                {
                    /*index of the kth ngbh of n*/
                    idx_t j = particles[i].ngbh.data[k].array_idx;
                    pp = NOBORDER;
                    /*Loop over kn neigbhours to find if n is the nearest*/
                    /*if cluster of the particle in nbhg is c then check is neighborhood*/                                                
                    if(particles[j].cluster_idx != c)
                    {
                        pp = j;
                        break;
                    }

                }
            }

            if(pp != NOBORDER)
            {
                for(idx_t k = 1; k < max_k; ++k)
                {
                    idx_t pp_ngbh_idx = particles[pp].ngbh.data[k].array_idx;
                    if(pp_ngbh_idx == i)
                    {
                        break;
                    }
                    if(particles[pp_ngbh_idx].cluster_idx == c)
                    {
                        pp = NOBORDER;
                        break;
                    }
                }
            }
                            /*if it is the maximum one add it to the cluster*/
            if(pp != NOBORDER)
            {
		int ppc = particles[pp].cluster_idx;
		if(cluster -> UseSparseBorders)
		{
			//insert one and symmetric one
			SparseBorder_t b = {.i = c, .j = ppc, .idx = i, .density = particles[i].g, .error = particles[i].log_rho_err}; 
			SparseBorder_Insert(cluster, b);
			//get symmetric border
			SparseBorder_t bsym = {.i = ppc, .j = c, .idx = i, .density = particles[i].g, .error = particles[i].log_rho_err}; 
			SparseBorder_Insert(cluster, bsym);

		}
		else
		{
			if(particles[i].g > borders[c][ppc].density)
			{
			    borders[c][ppc].density = particles[i].g;
			    borders[ppc][c].density = particles[i].g;
			    borders[c][ppc].idx = i;
			    borders[ppc][c].idx = i;
			}
		}
            }

    }


    if(cluster -> UseSparseBorders)
    {
	    for(idx_t c = 0; c < nclus; ++c)
	    {
		    for(idx_t el = 0; el < cluster -> SparseBorders[c].count; ++el)
		    {
			    //fix border density, write log rho c
			    idx_t idx = cluster -> SparseBorders[c].data[el].idx; 
			    cluster -> SparseBorders[c].data[el].density = particles[idx].log_rho_c;
		    }
	    }

    }
    else
    {
	    for(idx_t i = 0; i < nclus - 1; ++i)
	    {
		for(idx_t j = i + 1; j < nclus; ++j)
		{
		    idx_t p = borders[i][j].idx;
		    if(p != NOBORDER)
		    {   

			borders[i][j].density = particles[p].log_rho_c;
			borders[j][i].density = particles[p].log_rho_c;

			borders[i][j].error = particles[p].log_rho_err;
			borders[j][i].error = particles[p].log_rho_err;
		    }
		}
	    }

	    for(idx_t i = 0; i < nclus; ++i)
	    {
		borders[i][i].density = -1.0;
		borders[i][i].error = 0.0;
	    }
    }

    clock_gettime(CLOCK_MONOTONIC, &finish_tot);
    elapsed_tot = (finish_tot.tv_sec - start_tot.tv_sec);
    elapsed_tot += (finish_tot.tv_nsec - start_tot.tv_nsec) / 1000000000.0;
    printf("\tTotal time: %.3lfs\n\n", elapsed_tot);

    return;
    #undef borders
   }



void Merge_A_into_B(idx_t* who_amI, idx_t cluster_A, idx_t cluster_B, idx_t n)
{
    #pragma omp parallel if(n > MAX_SERIAL_MERGING)
    {
	    idx_t tmp;
	    #pragma omp for
	    for(idx_t i = 0; i < n; ++i)
	    {   
		//substitute occurencies of b with a 
		tmp = who_amI[i] == cluster_A ? cluster_B : who_amI[i];
		who_amI[i] = tmp;
	    }
    }
    return;
}


int compare_merging_density( const void *A, const void *B)
{
  FLOAT_TYPE DensA = ((merge_t*)A)->density;
  FLOAT_TYPE DensB = ((merge_t*)B)->density;

  return - ( DensA > DensB) + (DensA < DensB);
}


inline int is_a_merging( FLOAT_TYPE dens1, FLOAT_TYPE dens1_err,
			 FLOAT_TYPE dens2, FLOAT_TYPE dens2_err,
			 FLOAT_TYPE dens_border, FLOAT_TYPE dens_border_err,
			 FLOAT_TYPE Z)
/*
 * dens1 : the density of the particle that is the center of the first cluster
 * dens2 : the density of the particle that is the center of the second cluster
 * dens_border : the density of the border btw the cluster 1 and the cluster 2
 * *_err : the errors on the densities
 * Z     : the desired accuracy
 */
{
  /* in the original code it was:
   *
  FLOAT_TYPE a1 = particles[cluster->centers.data[i]].log_rho_c - border_density[i][j];
  FLOAT_TYPE a2 = particles[cluster->centers.data[j]].log_rho_c - border_density[i][j];
  
  FLOAT_TYPE e1 = Z*(particles[cluster->centers.data[i]].log_rho_err + border_err[i][j]);
  FLOAT_TYPE e2 = Z*(particles[cluster->centers.data[j]].log_rho_err + border_err[i][j]);
  */

  FLOAT_TYPE a1 = dens1 - dens_border;
  FLOAT_TYPE a2 = dens2 - dens_border;

  FLOAT_TYPE e1 = Z*(dens1_err + dens_border_err);
  FLOAT_TYPE e2 = Z*(dens2_err + dens_border_err);

  return (a1 < e1 || a2 < e2);
}


inline int merging_roles( FLOAT_TYPE dens1, FLOAT_TYPE dens1_err,
			  FLOAT_TYPE dens2, FLOAT_TYPE dens2_err,
			  FLOAT_TYPE dens_border, FLOAT_TYPE dens_border_err )
{
      
  FLOAT_TYPE c1 = (dens1 - dens_border) / (dens1_err + dens_border_err); 
  FLOAT_TYPE c2 = (dens2 - dens_border) / (dens2_err + dens_border_err);
  //printf("%.10lf %.10lf %d\n",c1,c2, c1 > c2);
  
  return ( c1 < c2 );     // if 1, this signal to swap 1 and 2
}

void fix_borders_A_into_B(idx_t A, idx_t B, border_t** borders, idx_t n)
{
   #pragma omp parallel for if(n > MAX_SERIAL_MERGING)
   for(idx_t i = 0; i < n; ++i) 
   {
        if(borders[A][i].idx != NOBORDER )
        {
            if(borders[B][i].idx != NOBORDER)
            {
                int mb = (borders[A][i].density > borders[B][i].density); 

                borders[B][i] = mb ? borders[A][i] : borders[B][i];
                borders[i][B] = borders[B][i];
            }
            else
            {
                borders[B][i] = borders[A][i];
                borders[i][B] = borders[B][i];
            }
        } 
        borders[A][i] = border_null;
        borders[i][A] = border_null;
   }
}

void inline Delete_adjlist_element(Clusters * c, const idx_t list_idx, const idx_t el)
{
	//swap last element with 
	idx_t count = c -> SparseBorders[list_idx].count;
	c -> SparseBorders[list_idx].data[el] = c -> SparseBorders[list_idx].data[count-1];
	c -> SparseBorders[list_idx].data[count-1] = SparseBorder_null;
	c -> SparseBorders[list_idx].count -= 1;
}

void fix_SparseBorders_A_into_B(idx_t s,idx_t t,Clusters* c)
{
	idx_t nclus = c -> centers.count;	
	//delete border trg -> src
	for(idx_t el = 0; el < c -> SparseBorders[t].count; ++el)
	{
		SparseBorder_t b = c -> SparseBorders[t].data[el];
		if(b.i == t && b.j == s)
		{
			//delete the border src trg
			Delete_adjlist_element(c, t, el);
		}
	}
	//find the border and delete it, other insert them in correct place
		
	for(idx_t el = 0; el < c -> SparseBorders[s].count; ++el)
	{
		SparseBorder_t b = c -> SparseBorders[s].data[el];
		if(b.j != t)
		{
			//insert these borders as trg -> j and j -> trg
			b.i = t;
			SparseBorder_Insert(c, b);
			SparseBorder_t bsym = b;
			bsym.i = b.j;
			bsym.j = b.i;
			SparseBorder_Insert(c, bsym);
		}
	}

	//clean up all borders
	//delete the src list
	AdjList_reset((c->SparseBorders) + s);
	//delete all borders containing src
	for(idx_t i = 0; i < nclus; ++i)
	{
		for(idx_t el = 0; el < c -> SparseBorders[i].count; ++el)
		{
			SparseBorder_t b = c -> SparseBorders[i].data[el];
			if(b.j == s)
			{
				//delete the border src trg
				Delete_adjlist_element(c, i, el);
			}
		}
			
	}


}

void Heuristic3_sparse(Clusters* cluster, Datapoint_info* particles, FLOAT_TYPE Z, int halo)
{
  printf("Using sparse implementation\n");
  #define borders cluster->borders

  struct timespec start_tot, finish_tot;
  double elapsed_tot;

  struct timespec start, finish;
  double elapsed;

  printf("H3: Merging clusters\n");
  clock_gettime(CLOCK_MONOTONIC, &start_tot);
  #ifdef VERBOSE
 	 clock_gettime(CLOCK_MONOTONIC, &start); 
  #endif

  idx_t nclus                 = cluster -> centers.count;  
  idx_t *  surviving_clusters = (idx_t*)malloc(nclus*sizeof(idx_t));
  for(idx_t i = 0; i < nclus; ++i)
    { 
        surviving_clusters[i] = i; 
    }

  idx_t   merge_count        = 0;
  idx_t   merging_table_size = 1000;
  merge_t *merging_table      = (merge_t*)malloc(sizeof(merge_t)*merging_table_size);
  
  /*Find clusters to be merged*/
  for(idx_t i = 0; i < nclus - 1; ++i)   
  {
    idx_t count = cluster -> SparseBorders[i].count;
    for(idx_t el = 0; el < count; ++el)   
    {
	      SparseBorder_t b = cluster -> SparseBorders[i].data[el];
	      if( b.j > b.i)
	      {
		      FLOAT_TYPE dens1           = particles[cluster->centers.data[b.i]].log_rho_c;
		      FLOAT_TYPE dens1_err       = particles[cluster->centers.data[b.i]].log_rho_err;
		      FLOAT_TYPE dens2           = particles[cluster->centers.data[b.j]].log_rho_c;
		      FLOAT_TYPE dens2_err       = particles[cluster->centers.data[b.j]].log_rho_err;
		      FLOAT_TYPE dens_border     = b.density;
		      FLOAT_TYPE dens_border_err = b.error;
	      
		      if ( is_a_merging( dens1, dens1_err, dens2, dens2_err, dens_border, dens_border_err, Z ) )
			{
			  
			  if ( merge_count == merging_table_size ) {
			    merging_table_size *= 1.1;
			    merging_table = (merge_t*)realloc( merging_table, sizeof(merge_t) * merging_table_size ); }

			  idx_t src = b.j;
			  idx_t trg = b.i;

			  merging_table[merge_count].source = src;
			  merging_table[merge_count].target = trg;
			  merging_table[merge_count].density = b.density;
			  ++merge_count;
			}
	      }

	}
            
  }

  qsort( (void*)merging_table, merge_count, sizeof(merge_t), compare_merging_density);

  #ifdef VERBOSE
	clock_gettime(CLOCK_MONOTONIC, &finish); 
	elapsed = (finish.tv_sec - start.tv_sec);
	elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
	printf("\tFinding merges:   %.3lfs\n", elapsed);
	clock_gettime(CLOCK_MONOTONIC, &start); 
  #endif
  
  
    for( idx_t m = 0; m < merge_count; m++ )
    {
      
        #define src surviving_clusters[merging_table[m].source]
        #define trg surviving_clusters[merging_table[m].target]
        //printf("Found: %lu, %lu which now is %lu, %lu\n",merging_table[m].source, merging_table[m].target, src,trg);

        int re_check = ( (src != merging_table[m].source) || (trg != merging_table[m].target) );
		idx_t new_src = (src < trg) ? src : trg;
		idx_t new_trg = (src < trg) ? trg : src;

                //pick who am I

                FLOAT_TYPE dens1           = particles[cluster->centers.data[new_src]].log_rho_c;
                FLOAT_TYPE dens1_err       = particles[cluster->centers.data[new_src]].log_rho_err;
                FLOAT_TYPE dens2           = particles[cluster->centers.data[new_trg]].log_rho_c;
                FLOAT_TYPE dens2_err       = particles[cluster->centers.data[new_trg]].log_rho_err;

		//borders get
		SparseBorder_t b 	   = SparseBorder_get(cluster, new_src, new_trg);
                FLOAT_TYPE dens_border     = b.density;
                FLOAT_TYPE dens_border_err = b.error;

                int i_have_to_merge = is_a_merging(dens1,dens1_err,dens2,dens2_err,dens_border,dens_border_err,Z);            
                switch (i_have_to_merge && src != trg)
                {
                case 1:
                    {
                        int side = merging_roles(dens1,dens1_err,dens2,dens2_err,dens_border,dens_border_err);
                        if(!side)
                        {
                            idx_t tmp;
                            tmp = new_src;
                            new_src = new_trg;
                            new_trg = tmp;
                        }

                        //borders[new_src][new_trg] = border_null;
                        //borders[new_trg][new_src] = border_null;
                        //printf("Merging %lu into %lu\n",new_src,new_trg);
                        fix_SparseBorders_A_into_B(new_src,new_trg,cluster);
                        Merge_A_into_B ( surviving_clusters, new_src, new_trg, nclus );	  
                    }
                    break;
                
                default:
                    break;
                }
        
        #undef src
        #undef trg
    }

  #ifdef VERBOSE
	clock_gettime(CLOCK_MONOTONIC, &finish); 
	elapsed = (finish.tv_sec - start.tv_sec);
	elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
	printf("\tCluster merging:  %.3lfs\n", elapsed);
	clock_gettime(CLOCK_MONOTONIC, &start); 
  #endif
  
    /*Finalize clustering*/
    /*Acutally copying */
    lu_dynamicArray tmp_centers;
    lu_dynamicArray tmp_cluster_idx;


    DynamicArray_Init(&tmp_centers);
    DynamicArray_Init(&tmp_cluster_idx);

    DynamicArray_Reserve(&tmp_centers, nclus);
    DynamicArray_Reserve(&tmp_cluster_idx, nclus);

    idx_t final_cluster_count = 0;

    idx_t* old_to_new = (idx_t*)malloc(nclus*sizeof(idx_t));
    idx_t incremental_k = 0;
    for(idx_t i = 0; i < nclus; ++i)
    {
        
        if(surviving_clusters[i] == i){
            DynamicArray_pushBack(&tmp_centers, cluster->centers.data[i]);
            DynamicArray_pushBack(&tmp_cluster_idx, i);
            old_to_new[i] = incremental_k;
            ++incremental_k;
            ++final_cluster_count;
        }
    }

    //fill the rest of old_to_new
    for(idx_t i = 0; i < nclus; ++i)
    {
        if(surviving_clusters[i] != i){
            idx_t cidx_to_copy_from = surviving_clusters[i];
            old_to_new[i] = old_to_new[cidx_to_copy_from];
        }
    }

    /*allocate auxiliary pointers to store results of the finalization of the procedure*/

    AdjList_t* tmp_borders      = (AdjList_t*)malloc(final_cluster_count*sizeof(AdjList_t));

    //initialize temporary borders
    for(idx_t i = 0; i < final_cluster_count; ++i)
    {
	    tmp_borders[i].count = 0;
	    tmp_borders[i].size  = PREALLOC_BORDERS;
	    tmp_borders[i].data  = (SparseBorder_t*)malloc(PREALLOC_BORDERS*sizeof(SparseBorder_t));
    }

    /*initialize all pointers*/

    /*Fix cluster assignment*/
    #pragma omp parallel for
    for(idx_t i = 0; i < cluster -> n; ++i)
    {
        particles[i].is_center = 0;
        int old_cidx = particles[i].cluster_idx;
        particles[i].cluster_idx = old_to_new[old_cidx];
    }

    
    #pragma omp parallel for
    for(idx_t c = 0; c < final_cluster_count; ++c)
    {
        idx_t c_idx = tmp_cluster_idx.data[c];
	for(idx_t el = 0; el < cluster -> SparseBorders[c_idx].count; ++el)
	{
		//retrieve border
		SparseBorder_t b = cluster -> SparseBorders[c_idx].data[el];
		//change idexes of clusters
		b.i = old_to_new[b.i];
		b.j = old_to_new[b.j];

		AdjList_Insert(tmp_borders + c, b);
	}
    }

    Clusters_Reset(cluster);
    /*pay attention to the defined borders*/
    /*copy into members*/
    cluster -> SparseBorders = tmp_borders;


    cluster -> centers = tmp_centers;
    /**
     * Fix center assignment
    */
    for(int i = 0; i < cluster -> centers.count; ++i)
    {
        int idx = cluster -> centers.data[i];
        particles[idx].is_center = 1;
    }
    /*Halo*/
    switch (halo)
    {
    case 1:
	{
		FLOAT_TYPE* max_border_den_array = (FLOAT_TYPE*)malloc(final_cluster_count*sizeof(FLOAT_TYPE));
		#pragma omp parallel
		{
		    #pragma omp for
		    for(idx_t c = 0; c < final_cluster_count; ++c)
		    {
			FLOAT_TYPE max_border_den = -2.;
			for(idx_t el = 0; el < cluster -> SparseBorders[c].count; ++el)
			{
			    SparseBorder_t b = cluster -> SparseBorders[c].data[el];
			    if(b.density > max_border_den)
			    {
				max_border_den = b.density;
			    }
			}
			max_border_den_array[c] = max_border_den;
		    }

		    #pragma omp barrier

		    #pragma omp for
		    for(idx_t i = 0; i < cluster -> n; ++i)
		    {
			int cidx = particles[i].cluster_idx;
			int halo_flag = particles[i].log_rho_c < max_border_den_array[cidx]; 
			particles[i].cluster_idx = halo_flag ? -1 : cidx;
		    }
		}
		free(max_border_den_array);
	}
        break;
    
    default:
        break;
    }    

    /*free memory and put the correct arrays into place*/
    free(tmp_cluster_idx.data);
    free(merging_table);
    //free(ipos.data);
    //free(jpos.data);
    free(surviving_clusters);
    free(old_to_new);

  #ifdef VERBOSE
	clock_gettime(CLOCK_MONOTONIC, &finish); 
	elapsed = (finish.tv_sec - start.tv_sec);
	elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
	printf("\tFinal operations: %.3lfs\n\n", elapsed);
  #endif

    clock_gettime(CLOCK_MONOTONIC, &finish_tot);
    elapsed_tot = (finish_tot.tv_sec - start_tot.tv_sec);
    elapsed_tot += (finish_tot.tv_nsec - start_tot.tv_nsec) / 1000000000.0;
    printf("\tFound %lu possible merges\n",(uint64_t)merge_count);
    printf("\tSurviving clusters %lu\n",(uint64_t)final_cluster_count);
    printf("\tTotal time: %.3lfs\n\n", elapsed_tot);

  #undef  borders  
}


void Heuristic3_dense(Clusters* cluster, Datapoint_info* particles, FLOAT_TYPE Z, int halo)
{
  printf("Using dense implementation\n");
  #define borders cluster->borders

  struct timespec start_tot, finish_tot;
  double elapsed_tot;

  struct timespec start, finish;
  double elapsed;

  printf("H3: Merging clusters\n");
  clock_gettime(CLOCK_MONOTONIC, &start_tot);
  #ifdef VERBOSE
 	 clock_gettime(CLOCK_MONOTONIC, &start); 
  #endif

  idx_t nclus              = cluster -> centers.count;  
  idx_t *  surviving_clusters = (idx_t*)malloc(nclus*sizeof(idx_t));
  for(idx_t i = 0; i < nclus; ++i)
    { 
        surviving_clusters[i] = i; 
    }

  idx_t   merge_count        = 0;
  idx_t   merging_table_size = 1000;
  merge_t *merging_table      = (merge_t*)malloc(sizeof(merge_t)*merging_table_size);
  
  /*Find clusters to be merged*/
  for(idx_t i = 0; i < nclus - 1; ++i)   
    for(idx_t j = i + 1; j < nclus; ++j)   
    {
	switch(borders[i][j].idx != NOBORDER)
	{
                    
	  case 1:		
	    {
	      FLOAT_TYPE dens1           = particles[cluster->centers.data[i]].log_rho_c;
	      FLOAT_TYPE dens1_err       = particles[cluster->centers.data[i]].log_rho_err;
	      FLOAT_TYPE dens2           = particles[cluster->centers.data[j]].log_rho_c;
	      FLOAT_TYPE dens2_err       = particles[cluster->centers.data[j]].log_rho_err;
	      FLOAT_TYPE dens_border     = borders[i][j].density;
	      FLOAT_TYPE dens_border_err = borders[i][j].error;
	      
	    if ( is_a_merging( dens1, dens1_err, dens2, dens2_err, dens_border, dens_border_err, Z ) )
		{
		  
		  if ( merge_count == merging_table_size ) {
		    merging_table_size *= 1.1;
		    merging_table = (merge_t*)realloc( merging_table, sizeof(merge_t) * merging_table_size ); }

		  //int swap = merging_roles( dens1, dens1_err, dens2, dens2_err, dens_border, dens_border_err);
		  idx_t src = j;
		  idx_t trg = i;
		  //switch ( swap )
		  //  {
		  //  case 0: { src = j; trg = i;} break;
		  //  case 1: { src = i; trg = j;} break;
		  //  }

		  merging_table[merge_count].source = src;
		  merging_table[merge_count].target = trg;
		  merging_table[merge_count].density = borders[src][trg].density;
          ++merge_count;
		}
	      break;
	    }
        default:
	    {
	      break;
	    }
            
	  }
      }

  qsort( (void*)merging_table, merge_count, sizeof(merge_t), compare_merging_density);
  #ifdef VERBOSE
	clock_gettime(CLOCK_MONOTONIC, &finish); 
	elapsed = (finish.tv_sec - start.tv_sec);
	elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
	printf("\tFinding merges:   %.3lfs\n", elapsed);
	clock_gettime(CLOCK_MONOTONIC, &start); 
  #endif
  
  
    for( idx_t m = 0; m < merge_count; m++ )
    {
      
        #define src surviving_clusters[merging_table[m].source]
        #define trg surviving_clusters[merging_table[m].target]
        //printf("Found: %lu, %lu which now is %lu, %lu\n",merging_table[m].source, merging_table[m].target, src,trg);

        int re_check = ( (src != merging_table[m].source) || (trg != merging_table[m].target) );
		idx_t new_src = (src < trg) ? src : trg;
		idx_t new_trg = (src < trg) ? trg : src;

                //pick who am I

                FLOAT_TYPE dens1           = particles[cluster->centers.data[new_src]].log_rho_c;
                FLOAT_TYPE dens1_err       = particles[cluster->centers.data[new_src]].log_rho_err;
                FLOAT_TYPE dens2           = particles[cluster->centers.data[new_trg]].log_rho_c;
                FLOAT_TYPE dens2_err       = particles[cluster->centers.data[new_trg]].log_rho_err;

                FLOAT_TYPE dens_border     = borders[new_src][new_trg].density;
                FLOAT_TYPE dens_border_err = borders[new_src][new_trg].error;

                int i_have_to_merge = is_a_merging(dens1,dens1_err,dens2,dens2_err,dens_border,dens_border_err,Z);            
                switch (i_have_to_merge && src != trg)
                {
                case 1:
                    {
                        int side = merging_roles(dens1,dens1_err,dens2,dens2_err,dens_border,dens_border_err);
                        if(!side)
                        {
                            idx_t tmp;
                            tmp = new_src;
                            new_src = new_trg;
                            new_trg = tmp;
                        }

                        borders[new_src][new_trg] = border_null;
                        borders[new_trg][new_src] = border_null;
                        //printf("Merging %lu into %lu\n",new_src,new_trg);
                        fix_borders_A_into_B(new_src,new_trg,borders,nclus);
                        Merge_A_into_B ( surviving_clusters, new_src, new_trg, nclus );	  
                    }
                    break;
                
                default:
                    break;
                }
        
        #undef src
        #undef trg
    }

  #ifdef VERBOSE
	clock_gettime(CLOCK_MONOTONIC, &finish); 
	elapsed = (finish.tv_sec - start.tv_sec);
	elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
	printf("\tCluster merging:  %.3lfs\n", elapsed);
	clock_gettime(CLOCK_MONOTONIC, &start); 
  #endif
  
    /*Finalize clustering*/
    /*Acutally copying */
    lu_dynamicArray tmp_centers;
    lu_dynamicArray tmp_cluster_idx;


    DynamicArray_Init(&tmp_centers);
    DynamicArray_Init(&tmp_cluster_idx);

    DynamicArray_Reserve(&tmp_centers, nclus);
    DynamicArray_Reserve(&tmp_cluster_idx, nclus);

    idx_t final_cluster_count = 0;

    idx_t* old_to_new = (idx_t*)malloc(nclus*sizeof(idx_t));
    idx_t incremental_k = 0;
    for(idx_t i = 0; i < nclus; ++i)
    {
        
        if(surviving_clusters[i] == i){
            DynamicArray_pushBack(&tmp_centers, cluster->centers.data[i]);
            DynamicArray_pushBack(&tmp_cluster_idx, i);
            old_to_new[i] = incremental_k;
            ++incremental_k;
            ++final_cluster_count;
        }
    }

    //fill the rest of old_to_new
    for(idx_t i = 0; i < nclus; ++i)
    {
        if(surviving_clusters[i] != i){
            idx_t cidx_to_copy_from = surviving_clusters[i];
            old_to_new[i] = old_to_new[cidx_to_copy_from];
        }
    }

    /*allocate auxiliary pointers to store results of the finalization of the procedure*/

    border_t** tmp_borders      = (border_t**)malloc(final_cluster_count*sizeof(border_t*));
    border_t*  tmp_borders_data = (border_t*)malloc(final_cluster_count*final_cluster_count*sizeof(border_t));

    /*initialize all pointers*/
    for(idx_t i = 0; i < final_cluster_count; ++i)
    {
        tmp_borders[i] = tmp_borders_data + i*final_cluster_count;
    }

    /*Fix cluster assignment*/
    #pragma omp parallel for
    for(idx_t i = 0; i < cluster -> n; ++i)
    {
        particles[i].is_center = 0;
        int old_cidx = particles[i].cluster_idx;
        particles[i].cluster_idx = old_to_new[old_cidx];
    }

    
    #pragma omp parallel for
    for(idx_t c = 0; c < final_cluster_count; ++c)
    {
        idx_t c_idx = tmp_cluster_idx.data[c];
        for(idx_t d = c; d < final_cluster_count; ++d)
        {
            idx_t c_jdx = tmp_cluster_idx.data[d];
            tmp_borders[c][d].density = borders[c_idx][c_jdx].density;
            tmp_borders[d][c].density = borders[c_idx][c_jdx].density;

            tmp_borders[c][d].idx = borders[c_idx][c_jdx].idx;
            tmp_borders[d][c].idx = borders[c_idx][c_jdx].idx;


            tmp_borders[c][d].error = borders[c_idx][c_jdx].error;
            tmp_borders[d][c].error = borders[c_idx][c_jdx].error;
        } 
    }

    Clusters_Reset(cluster);
    /*pay attention to the defined borders*/
    /*copy into members*/
    borders = tmp_borders;

    cluster -> __borders_data = tmp_borders_data;

    cluster -> centers = tmp_centers;
    /**
     * Fix center assignment
    */
    for(int i = 0; i < cluster -> centers.count; ++i)
    {
        int idx = cluster -> centers.data[i];
        particles[idx].is_center = 1;
    }
    /*Halo*/
    switch (halo)
    {
    case 1:
	{
		FLOAT_TYPE* max_border_den_array = (FLOAT_TYPE*)malloc(final_cluster_count*sizeof(FLOAT_TYPE));
		#pragma omp parallel
		{
		    #pragma omp for
		    for(idx_t c = 0; c < final_cluster_count; ++c)
		    {
			FLOAT_TYPE max_border_den = -2.;
			for(idx_t d = 0; d < final_cluster_count; ++d)
			{
			    if(tmp_borders[c][d].density > max_border_den)
			    {
				max_border_den = tmp_borders[c][d].density;
			    }
			}
			max_border_den_array[c] = max_border_den;
		    }

		    #pragma omp barrier

		    #pragma omp for
		    for(idx_t i = 0; i < cluster -> n; ++i)
		    {
			int cidx = particles[i].cluster_idx;
			int halo_flag = particles[i].log_rho_c < max_border_den_array[cidx]; 
			particles[i].cluster_idx = halo_flag ? -1 : cidx;
		    }
		}
		free(max_border_den_array);
	}
        break;
    
    default:
        break;
    }    

    /*free memory and put the correct arrays into place*/
    free(tmp_cluster_idx.data);
    free(merging_table);
    //free(ipos.data);
    //free(jpos.data);
    free(surviving_clusters);
    free(old_to_new);

  #ifdef VERBOSE
	clock_gettime(CLOCK_MONOTONIC, &finish); 
	elapsed = (finish.tv_sec - start.tv_sec);
	elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
	printf("\tFinal operations: %.3lfs\n\n", elapsed);
  #endif

    clock_gettime(CLOCK_MONOTONIC, &finish_tot);
    elapsed_tot = (finish_tot.tv_sec - start_tot.tv_sec);
    elapsed_tot += (finish_tot.tv_nsec - start_tot.tv_nsec) / 1000000000.0;
    printf("\tFound %lu possible merges\n", (uint64_t)merge_count);
    printf("\tSurviving clusters %lu\n", (uint64_t)final_cluster_count);
    printf("\tTotal time: %.3lfs\n\n", elapsed_tot);

  #undef  borders  
}


void Heuristic3(Clusters* cluster, Datapoint_info* particles, FLOAT_TYPE Z, int halo)
{
	if(cluster -> UseSparseBorders)
	{
		Heuristic3_sparse(cluster, particles,  Z,  halo);
	}
	else
	{
		Heuristic3_dense(cluster, particles,  Z,  halo);
	}
}
