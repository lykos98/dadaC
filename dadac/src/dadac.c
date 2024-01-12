//#include "../include/read_fof_snapshot.h"
#include "../include/dadac.h"
#include <math.h>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "./metrics.c"

#define MIN(x,y) x < y ? x : y
#define MAX(x,y) x > y ? x : y
#define DEPS 2.220446049250313e-16


#define DEFAULT_SLICE 10000 

#ifdef USE_BLAS
	#include <cblas.h>
	#include <sys/sysinfo.h>
#endif

#define MAX_SERIAL_MERGING 40000
#define MAX_N_NGBH 1000
#define PREALLOC_BORDERS 10

#ifdef VOPT
	#define SWMEM
#endif


unsigned int data_dims;
idx_t Npart;
const border_t border_null = {.density = -1.0, .error = 0, .idx = NOBORDER};
const SparseBorder_t SparseBorder_null = {.density = -1.0, .error = 0, .idx = NOBORDER, .i = NOBORDER, .j = NOBORDER};

void LinkedList_Insert(LinkedList* L, Node* n)
{
    ++(L -> count);
    n -> next = L -> head;
    L -> head = n;
}

int blas_are_in_use()
{
	#ifdef USE_BLAS
		return 1;
	#else
		return 0;
	#endif
}

/*****************************
 * Clusters object functions *
 *****************************/

void Clusters_allocate(Clusters * c, int s)
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

void old_Clusters_allocate(Clusters * c)
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
void KNN_search_kdTreeV2(Datapoint_info * dpInfo,kdNodeV2* kdNodeArray, kdNodeV2* root, idx_t n, idx_t k)
{
    struct timespec start_tot, finish_tot;
    double elapsed_tot;
    printf("KNN search:\n");
    clock_gettime(CLOCK_MONOTONIC, &start_tot);
	
	#ifdef PROGRESS_BAR
		idx_t progress_count = 0;
		idx_t step = n/100;
		printf("Progress 0/%lu -> 0%%\r",(uint64_t)n);
		fflush(stdout);
	#endif
    
    #pragma omp parallel
    {

	    #pragma omp for schedule(dynamic)
	    for(int p = 0; p < n; ++p)
	    {
			idx_t idx = kdNodeArray[p].array_idx;
			dpInfo[idx].ngbh = KNN_kdTreeV2(kdNodeArray[p].data, root, k);
			dpInfo[idx].cluster_idx = -1;
			dpInfo[idx].is_center = 0;
			dpInfo[idx].array_idx = idx;

			#ifdef PROGRESS_BAR
				idx_t aa;

				#pragma omp atomic capture
				aa = ++progress_count;

				if(aa % step == 0 )
				{
					printf("Progress %lu/%lu -> %u%%\r",(uint64_t)aa, (uint64_t)n, (uint32_t)((100*aa)/n) );
					fflush(stdout);
				}
			#endif
		}
    }

	

	#ifdef PROGRESS_BAR
		printf("Progress %lu/%lu -> 100%%\n",(uint64_t)n, (uint64_t)n);
	#endif
    //printf("Progress %lu/%lu\n",(uint64_t)progress_count, (uint64_t)n);

    clock_gettime(CLOCK_MONOTONIC, &finish_tot);
    elapsed_tot = (finish_tot.tv_sec - start_tot.tv_sec);
    elapsed_tot += (finish_tot.tv_nsec - start_tot.tv_nsec) / 1000000000.0;
    printf("\tTotal time: %.3lfs\n\n", elapsed_tot);
    return;

}


void KNN_search_kdtree(Datapoint_info * dpInfo, FLOAT_TYPE * data, kd_node* root, idx_t n, idx_t k)
{
    struct timespec start_tot, finish_tot;
    double elapsed_tot;
    printf("KNN search:\n");
    clock_gettime(CLOCK_MONOTONIC, &start_tot);
	
	#ifdef PROGRESS_BAR
		idx_t progress_count = 0;
		idx_t step = n/100;
		printf("Progress 0/%lu -> 0%%\r",(uint64_t)n);
		fflush(stdout);
	#endif

	/*		
	heap_node* preallocatedHeaps = (heap_node*)malloc(k*n*sizeof(heap_node));
    #pragma omp parallel
    {

	    #pragma omp for schedule(dynamic)
	    for(int p = 0; p < n; ++p)
	    {
		Heap H;
		H.count = 0;
		H.N = k;
		H.data = preallocatedHeaps + (k*p);

		KNN_sub_tree_search(data + data_dims*p, root, &H);
		HeapSort(&H);
		dpInfo[p].ngbh = H;
		dpInfo[p].array_idx = p;

		#ifdef PROGRESS_BAR
			idx_t aa;

			#pragma omp atomic capture
			aa = ++progress_count;

			if(aa % step == 0 )
			{
				printf("Progress %lu/%lu -> %u%%\r",(uint64_t)aa, (uint64_t)n, (uint32_t)((100*aa)/n) );
				fflush(stdout);
			}
		#endif
		}
    }
	*/
    
    #pragma omp parallel
    {

	    #pragma omp for schedule(dynamic)
	    for(int p = 0; p < n; ++p)
	    {
		dpInfo[p].ngbh = KNN(data + data_dims*p, root, k);
		dpInfo[p].array_idx = p;

		#ifdef PROGRESS_BAR
			idx_t aa;

			#pragma omp atomic capture
			aa = ++progress_count;

			if(aa % step == 0 )
			{
				printf("Progress %lu/%lu -> %u%%\r",(uint64_t)aa, (uint64_t)n, (uint32_t)((100*aa)/n) );
				fflush(stdout);
			}
		#endif
		}
    }

	

	#ifdef PROGRESS_BAR
		printf("Progress %lu/%lu -> 100%%\n",(uint64_t)n, (uint64_t)n);
	#endif
    //printf("Progress %lu/%lu\n",(uint64_t)progress_count, (uint64_t)n);

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
     

    //FLOAT_TYPE x_avg, y_avg;
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

FLOAT_TYPE idEstimate(Datapoint_info* dpInfo, idx_t n,FLOAT_TYPE fraction)
{

    /*********************************************************************************************
     * Estimation of the intrinsic dimension of a dataset                                        *
     * args:                                                                                     *
     * - dpInfo: array of structs                                                             *
     * - n: number of dpInfo                                                                  *
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

    //FLOAT_TYPE fraction = 0.7;
    FLOAT_TYPE* r = (FLOAT_TYPE*)malloc(n*sizeof(FLOAT_TYPE));
    FLOAT_TYPE* Pemp = (FLOAT_TYPE*)malloc(n*sizeof(FLOAT_TYPE));

    for(idx_t i = 0; i < n; ++i)
    {
        r[i] = 0.5 * log(dpInfo[i].ngbh.data[2].value/dpInfo[i].ngbh.data[1].value);
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

void computeRho(Datapoint_info* dpInfo, const FLOAT_TYPE d, const idx_t points){

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

    idx_t kMAX = dpInfo[0].ngbh.N - 1;   

    FLOAT_TYPE omega = 0.;  
    if(sizeof(FLOAT_TYPE) == sizeof(float)){ omega = powf(PI_F,d/2)/tgammaf(d/2.0f + 1.0f);}  
    else{omega = pow(M_PI,d/2.)/tgamma(d/2.0 + 1.0);}

    //printf("Omega d %f\n", omega);

    #pragma omp parallel for
    for(idx_t i = 0; i < points; ++i)
    {

        idx_t j = 4;
        idx_t k;
        FLOAT_TYPE dL  = 0.;
        FLOAT_TYPE vvi = 0.;
		FLOAT_TYPE vvj = 0.;
		FLOAT_TYPE vp  = 0.;
        while(j < kMAX  && dL < DTHR)
        {
            idx_t ksel = j - 1;
            vvi = omega * pow(dpInfo[i].ngbh.data[ksel].value,d/2.);

            idx_t jj = dpInfo[i].ngbh.data[j].array_idx;

            vvj = omega * pow(dpInfo[jj].ngbh.data[ksel].value,d/2.);

            vp = (vvi + vvj)*(vvi + vvj);
            dL = -2.0 * ksel * log(4.*vvi*vvj/vp);
            j = j + 1;
        }
        if(j == kMAX)
        {
            k = j - 1;
            vvi = omega * pow(dpInfo[i].ngbh.data[k].value,d/2.);
        }
        else
        {
            k = j - 2;
        }
        dpInfo[i].kstar = k;
        dpInfo[i].log_rho = log((FLOAT_TYPE)(k)/vvi/((FLOAT_TYPE)(points)));
        //dpInfo[i].log_rho = log((FLOAT_TYPE)(k)) - log(vvi) -log((FLOAT_TYPE)(points));
        dpInfo[i].log_rho_err =   1.0/sqrt((FLOAT_TYPE)k); //(FLOAT_TYPE)(-Q_rsqrt((float)k));
        dpInfo[i].g = dpInfo[i].log_rho - dpInfo[i].log_rho_err;
    }

    clock_gettime(CLOCK_MONOTONIC, &finish_tot);
    elapsed_tot = (finish_tot.tv_sec - start_tot.tv_sec);
    elapsed_tot += (finish_tot.tv_nsec - start_tot.tv_nsec) / 1000000000.0;
    printf("\tTotal time: %.3lfs\n\n", elapsed_tot);

    return;


}

void PAk(Datapoint_info* dpInfo, const FLOAT_TYPE d, const idx_t points)
{
    float_t omega = pow(M_PI,d/2.)/tgamma(d/2.0 + 1.0);
	computeRho(dpInfo,d,points);

	#pragma omp parallel
	{
		float_t* vi = (float_t*)malloc(dpInfo[0].ngbh.count*sizeof(float_t));
		#pragma omp for
		for(idx_t i = 0; i < points; ++i)
		{
			vi[0] = 0.;
			idx_t kstar = dpInfo[i].kstar;
			for(idx_t k = 0; k < kstar; ++k)
			{
				//vi[k] = omega * pow(dpInfo[i].ngbh.data[k].value/dpInfo[i].ngbh.data[k+1].value,d/2.);
				vi[k] = omega * (pow(dpInfo[i].ngbh.data[k+1].value,d/2.) - pow(dpInfo[i].ngbh.data[k].value,d/2.));
			}

			float_t f = dpInfo[i].log_rho + log((double)points);
			//if(i < 3)
			//{
			//	printf("ff %lf \n",f);
			//	printf("vi %lf %lf %lf\n",vi[0],vi[1],vi[2]);
			//}
			float_t a = 0.;

			float_t H00 = 0.;
			float_t H01 = 0.;
			float_t H10 = 0.;
			float_t H11 = 0.;

			float_t gf = 0.;
			float_t ga = 0.;

			float_t alpha = 0.;

			//1st iteration
			//
			gf = kstar;
			ga = (kstar + 1)*kstar/2.;
			for(idx_t k = 0; k < kstar; ++k)
			{
				idx_t l = k + 1;
				float_t dg = vi[k]*exp(f+a*l); 
				gf -= dg;
				ga -= l*dg;

				H00 -= dg;
				H10 -= dg*l;
				H11 -= l*l*dg;
			}

			H01 = H10;

			float_t dethess = (H00*H11 - H10*H01);
			float_t H00_inv = 0.;
			float_t H01_inv = 0.;
			float_t H10_inv = 0.;
			float_t H11_inv = 0.;
			float_t detinv 	= 0.;
			if(dethess > DEPS)
			{
				 detinv = 1/dethess;
				 H00_inv = +detinv*H11;
				 H11_inv = +detinv*H00;
				 H01_inv = -detinv*H01;
				 H10_inv = -detinv*H10;
			}

			float_t conv_condition = MAX(fabs(gf),fabs(ga)) > 1e-3;
			idx_t niter = 0;
			alpha = 0.1;
			while(conv_condition && niter < 10000)
			{
				if(dethess > DEPS)
				{
					f -= alpha*(gf*H00_inv + ga*H01_inv);	
					a -= alpha*(gf*H10_inv + ga*H11_inv);	
				}
				else
				{
					f = ( -a*kstar*(kstar+1)/2  - (gf-kstar) ) / kstar;
                  	a = ( -f*kstar - (gf-kstar)) / (kstar*(kstar+1)/2.);

				}
				H00 = 0.;
				H01 = 0.;
				H10 = 0.;
				H11 = 0.;

				gf = 0.;
				ga = 0.;

				gf = kstar;
				ga = (kstar + 1)*kstar/2.;
				for(idx_t k = 0; k < kstar; ++k)
				{
					idx_t l = k + 1;
					float_t dg = vi[k]*exp(f + a*l); 
					gf -= dg;
					ga -= l*dg;

					H00 -= dg;
					H10 -= l*dg;
					H11 -= l*l*dg;
				}

				H01 = H10;

				dethess = (H00*H11 - H10*H01);
				if(dethess > DEPS)
				{
					 detinv = 1/dethess;
					 H00_inv = +detinv*H11;
					 H11_inv = +detinv*H00;
					 H01_inv = -detinv*H01;
					 H10_inv = -detinv*H10;
				}

				if(fabs(gf) < DEPS || fabs(ga) < DEPS)
				{
					conv_condition = MAX(fabs(gf),fabs(ga)) > 1e-3;
				}
				else
				{
					conv_condition = MAX(fabs(gf/f),fabs(ga/a)) > 1e-3;
				}

				niter++;
			}
			dpInfo[i].log_rho =  f - log((double)points);
			dpInfo[i].log_rho_err = sqrt((float_t)(4*kstar + 2)/(float_t)((kstar-1)*kstar)); 
		}
	free(vi);

	}
		
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

void computeCorrection(Datapoint_info* dpInfo, idx_t n, FLOAT_TYPE Z)
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
            FLOAT_TYPE tmp = dpInfo[i].log_rho - Z*dpInfo[i].log_rho_err;
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
            dpInfo[i].log_rho_c = dpInfo[i].log_rho - min_log_rho + 1;
            dpInfo[i].g = dpInfo[i].log_rho_c - dpInfo[i].log_rho_err;
        }
    }
    //printf("%lf\n",min_log_rho);
}

Clusters Heuristic1(Datapoint_info* dpInfo, idx_t n)
{
    /**************************************************************
     * Heurisitc 1, from paper of Errico, Facco, Laio & Rodriguez *
     * ( https://doi.org/10.1016/j.ins.2021.01.010 )              *
     *                                                            *
     * args:                                                      *
     * - dpInfo: array of Datapoint structures                 *
     * - data: pointer to the dataset                             *
     * - n: number of Datapoints                                  *
     **************************************************************/

    struct timespec start_tot, finish_tot;
    double elapsed_tot;

    printf("H1: Preliminary cluster assignment\n");
    clock_gettime(CLOCK_MONOTONIC, &start_tot);

    //idx_t ncenters = 0;
    //idx_t putativeCenters = n;
    lu_dynamicArray allCenters, removedCenters, actualCenters, max_rho;
    DynamicArray_allocate(&allCenters);
    DynamicArray_allocate(&removedCenters);
    DynamicArray_allocate(&actualCenters);
    DynamicArray_allocate(&max_rho);

    Datapoint_info** dpInfo_ptrs = (Datapoint_info**)malloc(n*sizeof(Datapoint_info*));

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

        dpInfo_ptrs[i] = dpInfo + i;
        idx_t maxk = dpInfo[i].kstar + 1;
        FLOAT_TYPE gi = dpInfo[i].g;
        dpInfo[i].is_center = 1;
        dpInfo[i].cluster_idx = -1;
        //printf("%lf\n",p -> g);
        Heap i_ngbh = dpInfo[i].ngbh;
        for(idx_t k = 1; k < maxk; ++k)
        {
            idx_t ngbh_index = i_ngbh.data[k].array_idx;
            FLOAT_TYPE gj = dpInfo[ngbh_index].g;
            if(gj > gi){
                dpInfo[i].is_center = 0;
                break;
            }
        }
        if(dpInfo[i].is_center){
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
	qsort(dpInfo_ptrs, n, sizeof(Datapoint_info*), cmpPP);


	#ifdef NO_OPT_H1
	
    idx_t * to_remove = (idx_t*)malloc(allCenters.count*sizeof(idx_t));
    for(idx_t c = 0; c < allCenters.count; ++c) {to_remove[c] = MY_SIZE_MAX;}



    #pragma omp parallel
    {
            
        idx_t * to_remove_private = (idx_t*)malloc(allCenters.count*sizeof(idx_t));
    	for(idx_t c = 0; c < allCenters.count; ++c) {to_remove_private[c] = MY_SIZE_MAX;}

        #pragma omp for
        for(idx_t p = 0; p < n; ++p)
        {
        	Datapoint_info pp = *(dpInfo_ptrs[p]);
        	for(idx_t j = 1; j < pp.kstar + 1; ++j)
        	{
        		idx_t jidx = pp.ngbh.data[j].array_idx;
        		if(dpInfo[jidx].is_center && pp.g > dpInfo[jidx].g)
        		{
        			//dpInfo[jidx].is_center = 0;
        			for(idx_t c = 0; c < allCenters.count; ++c){
        				if(allCenters.data[c] == jidx)
					{

						if(to_remove_private[c] != MY_SIZE_MAX)
						{
							to_remove_private[c] = pp.g > 	dpInfo[to_remove_private[c]].g  ? pp.array_idx : to_remove_private[c];
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
					to_remove[c] = dpInfo[to_remove_private[c]].g > dpInfo[to_remove[c]].g ?
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
        //FLOAT_TYPE gi = dpInfo[i].g;
        idx_t mr = to_remove[p];
        if(mr != MY_SIZE_MAX)
        {
            //if(dpInfo[mr].g > gi) e = 1;
	    e = 1;
        }
        switch (e)
        {
            case 1:
                {
                    DynamicArray_pushBack(&removedCenters,i);
                    dpInfo[i].is_center = 0;
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
                    dpInfo[i].cluster_idx = actualCenters.count - 1;
                }
                break;
            default:
                break;
        }
    }
	free(to_remove);

	#else	
		
	idx_t* to_remove_mask = (idx_t*)malloc(n*sizeof(idx_t));
    for(idx_t p = 0; p < n; ++p) {to_remove_mask[p] = MY_SIZE_MAX;}

	
    #pragma omp parallel shared(to_remove_mask)
    {
        #pragma omp for
        for(idx_t p = 0; p < n; ++p)
        {
        	Datapoint_info pp = *(dpInfo_ptrs[p]);
			int flag = 0;
			idx_t ppp = 0;
			
        	for(idx_t j = 1; j < pp.kstar + 1; ++j)
			{
				idx_t jidx = pp.ngbh.data[j].array_idx; 
				if(dpInfo[jidx].is_center && pp.g > dpInfo[jidx].g)
				{
					
					#pragma omp critical 
					{
						ppp = to_remove_mask[jidx];
						flag = ppp != MY_SIZE_MAX;							
						to_remove_mask[jidx] = flag ? (pp.g > dpInfo[ppp].g ? pp.array_idx : ppp) : pp.array_idx; 
					}
					
					//#pragma omp atomic read 
					//ppp = to_remove_mask[jidx];

					//flag = ppp != MY_SIZE_MAX;							
					//
					//#pragma omp atomic write
					//to_remove_mask[jidx] = flag ? (pp.g > dpInfo[ppp].g ? pp.array_idx : ppp) : pp.array_idx; 
				}
			}
		}
	}
    
    

    for(idx_t p = 0; p < allCenters.count; ++p)
    {
        idx_t i = allCenters.data[p];
        int e = 0;
        //FLOAT_TYPE gi = dpInfo[i].g;
        idx_t mr = to_remove_mask[i];
        if(mr != MY_SIZE_MAX)
        {
            //if(dpInfo[mr].g > gi) e = 1;
			e = 1;
        }
        switch (e)
        {
            case 1:
                {
                    DynamicArray_pushBack(&removedCenters,i);
                    dpInfo[i].is_center = 0;
                    for(idx_t c = 0; c < removedCenters.count - 1; ++c)
                    {
                        if(mr == removedCenters.data[c])
                        {
                            mr = max_rho.data[c];
                        }
                    }
                    DynamicArray_pushBack(&max_rho,mr);
                    
                }
                break;
            case 0:
                {
                    DynamicArray_pushBack(&actualCenters,i);
                    dpInfo[i].cluster_idx = actualCenters.count - 1;
                }
                break;
            default:
                break;
        }
    }

	free(to_remove_mask);

	#endif

    #ifdef VERBOSE
        clock_gettime(CLOCK_MONOTONIC, &finish);
        elapsed = (finish.tv_sec - start.tv_sec);
        elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
        printf("\tFinding actual centers:   %.3lfs\n",elapsed);

        clock_gettime(CLOCK_MONOTONIC, &start);
    #endif


    //idx_t nclusters = 0;


    /*****************************************************************************
     * Sort all the dpInfo based on g and then perform the cluster assignment *
     * in asceding order                                                         *
     * UPDATE: dpInfo already sorted                                          *
     *****************************************************************************/
                                                                                

    //qsort(dpInfo_ptrs, n, sizeof(Datapoint_info*), cmpPP);

    for(idx_t i = 0; i < n; ++i)
    {   
        Datapoint_info* p = dpInfo_ptrs[i];
        //idx_t ele = p -> array_idx;
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
                cluster = dpInfo[p_idx].cluster_idx; 
            }

            //
            if(cluster == -1)
            {
                FLOAT_TYPE gmax = -99999.;               
                idx_t gm_index = 0;
                for(idx_t k = 0; k < max_k; ++k)
                {
                    idx_t ngbh_index = p -> ngbh.data[k].array_idx;
                    for(idx_t m = 0; m < removedCenters.count; ++m)
                    {
                        FLOAT_TYPE gcand = dpInfo[max_rho.data[m]].g;
                        if(ngbh_index == removedCenters.data[m] && gcand > gmax)
                        {   
                            //printf("%lu -- %lu\n", ele, m);
                            gmax = gcand;
                            gm_index = max_rho.data[m];
                        }
                    }
                }

                cluster = dpInfo[gm_index].cluster_idx;

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

    free(dpInfo_ptrs);
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

void Heuristic2(Clusters* cluster, Datapoint_info* dpInfo)
{

    #define borders cluster->borders

    struct timespec start_tot, finish_tot;
    double elapsed_tot;
    idx_t n = cluster -> n;

    printf("H2: Finding border points\n");
    clock_gettime(CLOCK_MONOTONIC, &start_tot);


    idx_t nclus = cluster->centers.count; 
    idx_t max_k = dpInfo[0].ngbh.N;

    for(idx_t i = 0; i < n; ++i)
    {
            idx_t pp = NOBORDER;
            /*loop over n neighbors*/
            int c = dpInfo[i].cluster_idx;
            if(!dpInfo[i].is_center)
            {
                for(idx_t k = 1; k < dpInfo[i].kstar + 1; ++k)
                {
                    /*index of the kth ngbh of n*/
                    idx_t j = dpInfo[i].ngbh.data[k].array_idx;
                    pp = NOBORDER;
                    /*Loop over kn neigbhours to find if n is the nearest*/
                    /*if cluster of the particle in nbhg is c then check is neighborhood*/                                                
                    if(dpInfo[j].cluster_idx != c)
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
                    idx_t pp_ngbh_idx = dpInfo[pp].ngbh.data[k].array_idx;
                    if(pp_ngbh_idx == i)
                    {
                        break;
                    }
                    if(dpInfo[pp_ngbh_idx].cluster_idx == c)
                    {
                        pp = NOBORDER;
                        break;
                    }
                }
            }
                            /*if it is the maximum one add it to the cluster*/
            if(pp != NOBORDER)
            {
		int ppc = dpInfo[pp].cluster_idx;
		if(cluster -> UseSparseBorders)
		{
			//insert one and symmetric one
			SparseBorder_t b = {.i = c, .j = ppc, .idx = i, .density = dpInfo[i].g, .error = dpInfo[i].log_rho_err}; 
			SparseBorder_Insert(cluster, b);
			//get symmetric border
			SparseBorder_t bsym = {.i = ppc, .j = c, .idx = i, .density = dpInfo[i].g, .error = dpInfo[i].log_rho_err}; 
			SparseBorder_Insert(cluster, bsym);

		}
		else
		{
			if(dpInfo[i].g > borders[c][ppc].density)
			{
			    borders[c][ppc].density = dpInfo[i].g;
			    borders[ppc][c].density = dpInfo[i].g;
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
			    cluster -> SparseBorders[c].data[el].density = dpInfo[idx].log_rho_c;
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

			borders[i][j].density = dpInfo[p].log_rho_c;
			borders[j][i].density = dpInfo[p].log_rho_c;

			borders[i][j].error = dpInfo[p].log_rho_err;
			borders[j][i].error = dpInfo[p].log_rho_err;
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


static inline int is_a_merging( FLOAT_TYPE dens1, FLOAT_TYPE dens1_err,
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
  FLOAT_TYPE a1 = dpInfo[cluster->centers.data[i]].log_rho_c - border_density[i][j];
  FLOAT_TYPE a2 = dpInfo[cluster->centers.data[j]].log_rho_c - border_density[i][j];
  
  FLOAT_TYPE e1 = Z*(dpInfo[cluster->centers.data[i]].log_rho_err + border_err[i][j]);
  FLOAT_TYPE e2 = Z*(dpInfo[cluster->centers.data[j]].log_rho_err + border_err[i][j]);
  */

  FLOAT_TYPE a1 = dens1 - dens_border;
  FLOAT_TYPE a2 = dens2 - dens_border;

  FLOAT_TYPE e1 = Z*(dens1_err + dens_border_err);
  FLOAT_TYPE e2 = Z*(dens2_err + dens_border_err);

  return (a1 < e1 || a2 < e2);
}


static inline int merging_roles( FLOAT_TYPE dens1, FLOAT_TYPE dens1_err,
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

static inline void Delete_adjlist_element(Clusters * c, const idx_t list_idx, const idx_t el)
{
	//swap last element with 
	idx_t count = c -> SparseBorders[list_idx].count;
	c -> SparseBorders[list_idx].data[el] = c -> SparseBorders[list_idx].data[count-1];
	c -> SparseBorders[list_idx].data[count-1] = SparseBorder_null;
	c -> SparseBorders[list_idx].count -= 1;
}

void fix_SparseBorders_A_into_B(idx_t s,idx_t t,Clusters* c)
{
	//delete border trg -> src
	
	//idx_t nclus = c -> centers.count;
	
	{
		{
			for(idx_t el = 0; el < c -> SparseBorders[t].count; ++el)
			{
				SparseBorder_t b = c -> SparseBorders[t].data[el];
				if(b.i == t && b.j == s)
				{
					//delete the border src trg
					Delete_adjlist_element(c, t, el);
				}
			}
		}
		//find the border and delete it, other insert them in correct place
		for(idx_t el = 0; el < c -> SparseBorders[s].count; ++el)
		{
			SparseBorder_t b = c -> SparseBorders[s].data[el];
		//	idx_t ii = b.i;
			if(b.j != t)
			{
				//insert these borders as trg -> j and j -> trg
				b.i = t;
				SparseBorder_Insert(c, b);
				SparseBorder_t bsym = b;
				bsym.i = b.j;
				bsym.j = b.i;
				SparseBorder_Insert(c, bsym);
				for(idx_t dl = 0; dl < c -> SparseBorders[b.j].count; ++dl)
				{
					SparseBorder_t b_del = c -> SparseBorders[b.j].data[dl];
					if(b_del.j == s)
					{
						//delete the border src trg
						Delete_adjlist_element(c, b.j, dl);
					}
				}
						
			}
		}
		//clean up all borders
		//delete the src list
		{
			AdjList_reset((c->SparseBorders) + s);
		}
		//delete all borders containing src
	//	for(idx_t i = 0; i < nclus; ++i)
	//	{
	//		for(idx_t el = 0; el < c -> SparseBorders[i].count; ++el)
	//		{
	//			SparseBorder_t b = c -> SparseBorders[i].data[el];
	//			if(b.j == s)
	//			{
	//				//delete the border src trg
	//				Delete_adjlist_element(c, i, el);
	//			}
	//		}
	//			
	//	}
	}


}

void Heuristic3_sparse(Clusters* cluster, Datapoint_info* dpInfo, FLOAT_TYPE Z, int halo)
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
		      FLOAT_TYPE dens1           = dpInfo[cluster->centers.data[b.i]].log_rho_c;
		      FLOAT_TYPE dens1_err       = dpInfo[cluster->centers.data[b.i]].log_rho_err;
		      FLOAT_TYPE dens2           = dpInfo[cluster->centers.data[b.j]].log_rho_c;
		      FLOAT_TYPE dens2_err       = dpInfo[cluster->centers.data[b.j]].log_rho_err;
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

        //int re_check = ( (src != merging_table[m].source) || (trg != merging_table[m].target) );
	//if(re_check)
	{
		idx_t new_src = (src < trg) ? src : trg;
		idx_t new_trg = (src < trg) ? trg : src;

                //pick who am I

                FLOAT_TYPE dens1           = dpInfo[cluster->centers.data[new_src]].log_rho_c;
                FLOAT_TYPE dens1_err       = dpInfo[cluster->centers.data[new_src]].log_rho_err;
                FLOAT_TYPE dens2           = dpInfo[cluster->centers.data[new_trg]].log_rho_c;
                FLOAT_TYPE dens2_err       = dpInfo[cluster->centers.data[new_trg]].log_rho_err;

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
        dpInfo[i].is_center = 0;
        int old_cidx = dpInfo[i].cluster_idx;
        dpInfo[i].cluster_idx = old_to_new[old_cidx];
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
    for(idx_t i = 0; i < cluster -> centers.count; ++i)
    {
        int idx = cluster -> centers.data[i];
        dpInfo[idx].is_center = 1;
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
				int cidx = dpInfo[i].cluster_idx;
				int halo_flag = dpInfo[i].log_rho_c < max_border_den_array[cidx]; 
				//int halo_flag = max_border_den_array[cidx] > dpInfo[i].log_rho_c  ; 
				dpInfo[i].cluster_idx = halo_flag ? -1 : cidx;
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


void Heuristic3_dense(Clusters* cluster, Datapoint_info* dpInfo, FLOAT_TYPE Z, int halo)
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
	      FLOAT_TYPE dens1           = dpInfo[cluster->centers.data[i]].log_rho_c;
	      FLOAT_TYPE dens1_err       = dpInfo[cluster->centers.data[i]].log_rho_err;
	      FLOAT_TYPE dens2           = dpInfo[cluster->centers.data[j]].log_rho_c;
	      FLOAT_TYPE dens2_err       = dpInfo[cluster->centers.data[j]].log_rho_err;
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

        //int re_check = ( (src != merging_table[m].source) || (trg != merging_table[m].target) );
	//if(re_check)
	{
		idx_t new_src = (src < trg) ? src : trg;
		idx_t new_trg = (src < trg) ? trg : src;

                //pick who am I

                FLOAT_TYPE dens1           = dpInfo[cluster->centers.data[new_src]].log_rho_c;
                FLOAT_TYPE dens1_err       = dpInfo[cluster->centers.data[new_src]].log_rho_err;
                FLOAT_TYPE dens2           = dpInfo[cluster->centers.data[new_trg]].log_rho_c;
                FLOAT_TYPE dens2_err       = dpInfo[cluster->centers.data[new_trg]].log_rho_err;

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
        dpInfo[i].is_center = 0;
        int old_cidx = dpInfo[i].cluster_idx;
        dpInfo[i].cluster_idx = old_to_new[old_cidx];
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
    for(idx_t i = 0; i < cluster -> centers.count; ++i)
    {
        int idx = cluster -> centers.data[i];
        dpInfo[idx].is_center = 1;
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
				int cidx = dpInfo[i].cluster_idx;
				int halo_flag = dpInfo[i].log_rho_c < max_border_den_array[cidx];  
				//int halo_flag =  (max_border_den_array[cidx] > dpInfo[i].log_rho_c) - (max_border_den_array[cidx] < dpInfo[i].log_rho_c); 
				//halo_flag = halo_flag > 0;
				dpInfo[i].cluster_idx = halo_flag ? -1 : cidx;
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


void Heuristic3(Clusters* cluster, Datapoint_info* dpInfo, FLOAT_TYPE Z, int halo)
{
	if(cluster -> UseSparseBorders)
	{
		Heuristic3_sparse(cluster, dpInfo,  Z,  halo);
	}
	else
	{
		Heuristic3_dense(cluster, dpInfo,  Z,  halo);
	}
}

void computeLevel(kd_node* root, idx_t prev_lvl)
{
	idx_t curr_lvl = prev_lvl + 1;
	root -> level = curr_lvl;
	if(root -> lch) computeLevel(root -> lch, curr_lvl);
	if(root -> rch) computeLevel(root -> rch, curr_lvl);
	return;
}

Heap KNN_bruteforce(void* point, void* data, size_t n, idx_t elementSize, idx_t k, float_t (*metric)(void*, void*))
{
    Heap H;
    allocateHeap(&H,k);
    initHeap(&H);
	for(idx_t j = 0; j < n; ++j)
	{
		float_t distance = metric(point, data + elementSize*j);
		insertMaxHeap_InsertionSort(&H, distance, j);
		//insertMaxHeap(&H, distance, j);
	}
	//HeapSort(&H);
	//for(size_t i = 0; i < H.count; ++i) H.data[i].value = H.data[i].value*H.data[i].value;
    return H;
}

int partition_heapNode(heap_node *array, int left, int right, int pivotIndex) {
    float_t pivotValue = array[pivotIndex].value;
    int storeIndex = left;
    int i;
    /* Move pivot to end */
    swapHeapNode(array + pivotIndex, array + right);
    for(i=left; i < right; i = i + 1 ){
        if(array[i].value < pivotValue){
    		swapHeapNode(array + storeIndex, array + i);
            storeIndex += 1;
        }
    }
    /* Move pivot to its final place */
    swapHeapNode(array + storeIndex , array + right);

    return storeIndex;
}

int qselect_heapNode(heap_node *array, int left, int right, int n) {
    int pivotIndex;
    if(left == right){
        return left;
    }
    pivotIndex = left + (rand() % (right-left + 1)); /* random int left <= x <= right */
    pivotIndex = partition_heapNode(array, left, right, pivotIndex);
    /* The pivot is in its final sorted position */
    if(n == pivotIndex){
        return pivotIndex;
    }else if(n < pivotIndex){
        return qselect_heapNode(array, left, pivotIndex-1, n);
    }else{
        return qselect_heapNode(array, pivotIndex+1, right, n);
    }
}

int quickselect_heapNode(heap_node *array, int array_size, int k){
    return qselect_heapNode(array, 0, array_size-1, k-1);
}

void KNN_BruteForce(Datapoint_info* points, void* data,size_t n, size_t byteSize, size_t dims, size_t k, float_t (*metric)(void*, void*))
{
	/*
	 * TODO: fix when threads are more than points
	 * */
	#pragma omp parallel
	{

		idx_t slice_len = n / omp_get_num_threads();
		heap_node* pvt_working_mem = (heap_node*)malloc(sizeof(heap_node)*n);

		#pragma omp for 
		for(idx_t slice = 0; slice < n; slice += slice_len)
		{
			for(idx_t p = slice; p < slice + slice_len; ++p)
			{
				for(idx_t j = 0; j < n; ++j)
				{
					pvt_working_mem[j].value = metric(data + p*dims*byteSize, data + j*dims*byteSize);
					pvt_working_mem[j].array_idx = j; 
				}
				//printf("Thread %d got slice %lu to %lu\n ", omp_get_thread_num(), slice, slice + slice_len);

				//points[p].ngbh = KNN_bruteforce(data + p*dims*byteSize, data,  n, byteSize*dims, k, metric);
				//
				Heap H = points[p].ngbh;

				quickselect_heapNode(pvt_working_mem, n,  k + 1);
				qsort(pvt_working_mem, k + 1, sizeof(heap_node), cmpHeapNodes);
				memcpy(H.data, pvt_working_mem, k*sizeof(heap_node));
			

				#ifdef PROGRESS_BAR
					idx_t aa;

					#pragma omp atomic capture
					aa = ++progress_count;

					if(aa % step == 0 )
					{
						printf("Progress %lu/%lu -> %u%%\r",(uint64_t)aa, (uint64_t)n, (uint32_t)((100*aa)/n) );
						fflush(stdout);
					}
				#endif
			}
		}

		idx_t remainder = omp_get_num_threads()*slice_len;
		#pragma omp for 
		for(idx_t p = remainder; p < n; ++p)
		{

			//points[p].ngbh = KNN_bruteforce(data + p*dims*byteSize, data,  n, byteSize*dims, k, metric);
			for(idx_t j = 0; j < n; ++j)
			{
				pvt_working_mem[j].value = metric(data + p*dims*byteSize, data + j*dims*byteSize);
				pvt_working_mem[j].array_idx = j; 
			}

			Heap H = points[p].ngbh;
			quickselect_heapNode(pvt_working_mem, n,  k + 1);
			qsort(pvt_working_mem, k, sizeof(heap_node), cmpHeapNodes);
			memcpy(H.data, pvt_working_mem, k*sizeof(heap_node));
		
			#ifdef PROGRESS_BAR
				idx_t aa;

				#pragma omp atomic capture
				aa = ++progress_count;

				if(aa % step == 0 )
				{
					printf("Progress %lu/%lu -> %u%%\r",(uint64_t)aa, (uint64_t)n, (uint32_t)((100*aa)/n) );
					fflush(stdout);
				}
			#endif
		}

		free(pvt_working_mem);

	}

}

#ifdef USE_BLAS
	void __handle_w_blas_d(Datapoint_info* points, double* data,size_t n, size_t dims, size_t k)
	{
		struct sysinfo info;
		sysinfo(&info);

		size_t slice_size = (size_t)((float_t)(info.freeram)*0.9/(float_t)(8*n));
		slice_size = MIN(n,slice_size);
		printf("%lu slice size \n",slice_size);
		float_t* norms_sq = (float_t*)malloc(n*sizeof(double));
		#pragma omp parallel for
		for(size_t i = 0; i < n; ++i)
		{
			norms_sq[i] = 0.;
			for(size_t j = 0; j < dims; ++j) norms_sq[i] += data[i*dims+j]*data[i*dims+j];
		}

		heap_node* working_mem 	= NULL;

		#pragma omp parallel
		{
			#pragma omp master
			{
				working_mem = (heap_node*)malloc(omp_get_num_threads()*sizeof(heap_node)*n); 
			}
		}

		float_t* middle_terms 	= (float_t*)malloc(slice_size*sizeof(float_t)*n); 

		size_t s = 0;
		for(s = 0; s < n; s+=slice_size)
		{
			size_t lower_idx = s; 
			size_t upper_idx = MIN(n,s + slice_size);
				#pragma omp parallel
				{
					#pragma omp for
					for(size_t i = lower_idx; i < upper_idx; ++i)
						for(size_t j = 0; j < n; ++j)
						{
							middle_terms[(i%slice_size)*n + j] = norms_sq[i] + norms_sq[j];
							//middle_terms[(i%DEFAULT_SLICE)*n + j] = 0;
						}
				}

				float_t* A = data + dims*lower_idx;
				cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
									upper_idx - lower_idx,n,dims, -2.0,
									A, dims, data, dims, 1.0, middle_terms,n);

				//reduce the arrays
				#pragma omp parallel 
				{
					heap_node* my_working_mem = working_mem + omp_get_thread_num()*n;
					#pragma omp for
					for(size_t i = lower_idx; i < upper_idx; ++i)
					{
						//for(size_t j = 0; j < n; ++j)
						//{
						////	my_working_mem[j].value = middle_terms[(i%DEFAULT_SLICE)*n + j];
						////	my_working_mem[j].array_idx = j;
						//	insertMaxHeap(&(points[i].ngbh), middle_terms[(i%DEFAULT_SLICE)*n + j],j);
						//}
						//HeapSort(&(points[i].ngbh));
						//copy into the working mem
						for(size_t j = 0; j < n; ++j)
						{
							my_working_mem[j].value = middle_terms[(i%slice_size)*n + j];
							my_working_mem[j].array_idx = j;
						}

						Heap H = points[i].ngbh;
						quickselect_heapNode(my_working_mem, n,  k + 1);
						//printf("%lu %lf \t ",my_working_mem[k].array_idx, my_working_mem[k].value);
						qsort(my_working_mem, k, sizeof(heap_node), cmpHeapNodes);
						//printf("%lu %lf\n",my_working_mem[k].array_idx, my_working_mem[k].value);
						memcpy(H.data, my_working_mem, k*sizeof(heap_node));
					}
				}
		}

		free(working_mem);
		free(middle_terms);

		free(norms_sq);	

	}
#endif

Datapoint_info* NgbhSearch_bruteforce(void* data, size_t n, size_t byteSize, size_t dims, size_t k, float_t (*metric)(void *, void *))
{
	/*
	 *
	 * /!\ Extremely inefficient implementation for euclidean metric
	 *     currently working on optimized one with the use of blas
	 *
	 * */

	METRICS_DATADIMS = dims;

	#ifdef VERBOSE
		printf("Brute-forcing kNN computation:\n");
		if(!blas_are_in_use())
		{
			printf("/!\\ Naive implementation, terribly unoptimized and slow, do not use unless forced to \n");
			printf("/!\\ If you want to use euclidean metric follow the instruction to compile and link against\n");
			printf("/!\\ a CBLAS implementation to drammatically speed up computation\n");
			if(!metric)
			{
				printf("\n");
				printf("/!\\ Pass a valid function as a metric! Exiting\n");
				exit(1);
			}

		}
		else
		{
			if(metric)	
			{
				printf("/!\\ CBLAS accelerated implementation, if you want to use euclidean metric \n");
				printf("/!\\ pass NULL to the function pointer of the metric \n");
			}
		}
	#endif
    Datapoint_info* points = (Datapoint_info*)malloc(n*sizeof(Datapoint_info));

    struct timespec start_tot, finish_tot;
    double elapsed_tot;
    printf("KNN search:\n");
    clock_gettime(CLOCK_MONOTONIC, &start_tot);
	
	#ifdef PROGRESS_BAR
		idx_t progress_count = 0;
		idx_t step = n/100;
		printf("Progress 0/%lu -> 0%%\r",(uint64_t)n);
		fflush(stdout);
	#endif
	#pragma omp parallel for
	for(idx_t j = 0; j < n; ++j)
	{
		Heap H;
		allocateHeap(&H, k);
		initHeap(&H);
		points[j].ngbh = H;
		points[j].array_idx = j;

	}

	if(metric)
	{
		KNN_BruteForce(points, data, n, byteSize, dims, k, metric);
	}
	else
	{
		#ifdef USE_BLAS
			__handle_w_blas_d(points, data, n, dims, k);
		#else
			KNN_BruteForce(points, data, n, byteSize, dims, k, eud);
		#endif

	}






	#ifdef PROGRESS_BAR
		printf("Progress %lu/%lu -> 100%%\n",(uint64_t)n, (uint64_t)n);
	#endif
    //printf("Progress %lu/%lu\n",(uint64_t)progress_count, (uint64_t)n);

    clock_gettime(CLOCK_MONOTONIC, &finish_tot);
    elapsed_tot = (finish_tot.tv_sec - start_tot.tv_sec);
    elapsed_tot += (finish_tot.tv_nsec - start_tot.tv_nsec) / 1000000000.0;
    printf("\tTotal time: %.3lfs\n\n", elapsed_tot);
	return points;

}

Datapoint_info* NgbhSearch_kdtree(FLOAT_TYPE* data, size_t n, size_t ndims, size_t k)
{
    struct timespec start, finish;
    double elapsed;


	data_dims = (unsigned int)ndims;

	#ifdef VERBOSE
		printf("Building the KDtree v1:\n");
		clock_gettime(CLOCK_MONOTONIC, &start);
	#endif

    kd_node* kd_node_array = (kd_node*)malloc(n*sizeof(kd_node));
    kd_node** kd_ptrs = (kd_node**)malloc(n*sizeof(kd_node*));

    initializeKDnodes(kd_node_array,data,n);
    initializePTRS(kd_ptrs, kd_node_array,n);

    kd_node* root = build_tree(kd_ptrs, n, ndims);

    //printf("The root of the tree is\n");
    //printKDnode(root);
	#ifdef VERBOSE
		clock_gettime(CLOCK_MONOTONIC, &finish);
		elapsed = (finish.tv_sec - start.tv_sec);
		elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
		printf("\tTotal time: %.3lfs\n\n", elapsed);
	#endif


    Datapoint_info* points = (Datapoint_info*)malloc(n*sizeof(Datapoint_info));

    /**************
     * KNN search *
     **************/
	//int lvls[25] = {0};
	//computeLevel(root, 0);
	//for(int i = 0; i < n; ++i) lvls[kd_ptrs[i]->level]++;
	//for(int i = 0; i < 25; ++i) printf("lvl %d counts %d\n", i, lvls[i]);
    KNN_search_kdtree(points,data, root, n, k);

    free(kd_ptrs);
    free(kd_node_array);

	return points;
}

Datapoint_info* NgbhSearch_kdtree_V2(FLOAT_TYPE* data, size_t n, size_t ndims, size_t k)
{
    struct timespec start, finish;
    double elapsed;


	data_dims = (unsigned int)ndims;

	#ifdef VERBOSE
		printf("Building the KDtree v2:\n");
		clock_gettime(CLOCK_MONOTONIC, &start);
	#endif

	#ifdef SWMEM
		FLOAT_TYPE* dummy_data = malloc(sizeof(FLOAT_TYPE)*ndims*n);
		memcpy(dummy_data,data,sizeof(FLOAT_TYPE)*ndims*n);
		data = dummy_data;
	#endif

    kdNodeV2* kdNode_array = (kdNodeV2*)malloc(n*sizeof(kdNodeV2));

    initializeKDnodesV2(kdNode_array,data,n);

    kdNodeV2* root = build_tree_kdTreeV2(kdNode_array, n, ndims);

    //printf("The root of the tree is\n");
    //printKDnode(root);
	#ifdef VERBOSE
		clock_gettime(CLOCK_MONOTONIC, &finish);
		elapsed = (finish.tv_sec - start.tv_sec);
		elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
		printf("\tTotal time: %.3lfs\n\n", elapsed);
	#endif


    Datapoint_info* points = (Datapoint_info*)malloc(n*sizeof(Datapoint_info));

    /**************
     * KNN search *
     **************/
	//int lvls[25] = {0};
	//computeLevel(root, 0);
	//for(int i = 0; i < n; ++i) lvls[kd_ptrs[i]->level]++;
	//for(int i = 0; i < 25; ++i) printf("lvl %d counts %d\n", i, lvls[i]);
	KNN_search_kdTreeV2(points,kdNode_array,root,n,k);

	#ifdef SWMEM
		free(data);
	#endif
	
	for(idx_t i = 0; i < n; ++i) if(kdNode_array[i].nodeList.data) free(kdNode_array[i].nodeList.data);
	free(kdNode_array);

	return points;
}

void freeDatapointArray(Datapoint_info* d, size_t n)
{
    for (idx_t i = 0; i < n; ++i)
    {        
        freeHeap(&d[i].ngbh);
    }
    free(d);
}

int FloatAndUintSize()
{
	int v = 0;
	int vf = sizeof(FLOAT_TYPE) == 8 ? 1 : 0; 
	int vi = sizeof(idx_t) == 8 ? 1 : 0; 
	v = vf + vi*2;
	return v;
}


void KNN_search_vpTree(Datapoint_info* dpInfo, vpTreeNode* vpNodeArray,vpTreeNode* root,idx_t k,size_t n, float_t (*metric)(void*, void*))
{	
    struct timespec start_tot, finish_tot;
    double elapsed_tot;
    printf("KNN search:\n");
    clock_gettime(CLOCK_MONOTONIC, &start_tot);
	

	#ifdef PROGRESS_BAR
		idx_t progress_count = 0;
		idx_t step = n/100;
		printf("Progress 0/%lu -> 0%%\r",(uint64_t)n);
	#endif
    fflush(stdout);
    
    #pragma omp parallel
    {

		#ifdef ITERATIVE_VPTREE
			stack_vpTreeNodes stack;
			stackInit(&stack);
		#endif

	    #pragma omp for schedule(dynamic)
		for(size_t i = 0; i < n; ++i) 
		{
			#ifdef ITERATIVE_VPTREE
				dpInfo[i].ngbh = KNN_vpTree(vpNodeArray[i].data, root,  k, &stack, metric);
			#else
				dpInfo[i].ngbh = KNN_vpTree(vpNodeArray[i].data, root,  k, metric);
			#endif 

			dpInfo[i].cluster_idx = -1;
			dpInfo[i].is_center = 0;
			dpInfo[i].array_idx = i;
			
			#ifdef PROGRESS_BAR
				idx_t aa;
				#pragma omp atomic capture
				aa = ++progress_count;

				if(aa % step == 0 )
				{
					printf("Progress %lu/%lu -> %u%%\r",(uint64_t)aa, (uint64_t)n, (uint32_t)((100*aa)/n) );
					fflush(stdout);
				}
			#endif
		}

		#ifdef ITERATIVE_VPTREE
			free(stack.data);
		#endif
	}
	

    
	
	#ifdef PROGRESS_BAR
		printf("Progress %lu/%lu -> 100%%\n",(uint64_t)n, (uint64_t)n);
	#endif
    //printf("Progress %lu/%lu\n",(uint64_t)progress_count, (uint64_t)n);

    clock_gettime(CLOCK_MONOTONIC, &finish_tot);
    elapsed_tot = (finish_tot.tv_sec - start_tot.tv_sec);
    elapsed_tot += (finish_tot.tv_nsec - start_tot.tv_nsec) / 1000000000.0;
    printf("\tTotal time: %.3lfs\n\n", elapsed_tot);
    return;
}


Datapoint_info* NgbhSearch_vptree(void* data, size_t n, size_t byteSize, size_t dims, size_t k, float_t (*metric)(void *, void *))
{

    struct timespec start, finish;
    double elapsed;
	METRICS_DATADIMS = (uint32_t)dims;
	
	#ifdef VERBOSE
		clock_gettime(CLOCK_MONOTONIC, &start);
		printf("Building the vp tree v1\n");
	#endif

	vpTreeNode* vpNodeArray = (vpTreeNode*)malloc(n*sizeof(vpTreeNode));
	vpTreeNode** vpPtrArray = (vpTreeNode**)malloc(n*sizeof(vpTreeNode*));
	initialize_vpTreeNode_array(vpNodeArray, data, n, byteSize*dims);
	initialize_vpTreeNodes_pointers(vpPtrArray, vpNodeArray, n);

	vpTreeNode* root = build_vpTree(vpPtrArray, 0, n-1, NULL, metric);
	#ifdef VERBOSE
		clock_gettime(CLOCK_MONOTONIC, &finish);
		elapsed = (finish.tv_sec - start.tv_sec);
		elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    	printf("\tTotal time: %.3lfs\n\n", elapsed);
	#endif

    Datapoint_info* points = (Datapoint_info*)malloc(n*sizeof(Datapoint_info));

    /**************
     * KNN search *
     **************/

	
	KNN_search_vpTree(points, vpNodeArray, root, k, n, metric);
	

	//printf("NODE STRAMBO: %lf %p %lu\n", vpNodeArray[1516].mu, vpNodeArray[1516].inside, vpNodeArray[1516].parent -> array_idx);
    free(vpPtrArray);
    free(vpNodeArray);
	return points;
}

void KNN_search_vpTree_V2(Datapoint_info* dpInfo, vpTreeNodeV2* vpNodeArray,vpTreeNodeV2* root,idx_t k,size_t n, float_t (*metric)(void*, void*))
{	
    struct timespec start_tot, finish_tot;
    double elapsed_tot;
    printf("KNN search:\n");
    clock_gettime(CLOCK_MONOTONIC, &start_tot);

	#ifdef PROGRESS_BAR
		idx_t progress_count = 0;
		idx_t step = n/100;
		printf("Progress 0/%lu -> 0%%\r",(uint64_t)n);
		fflush(stdout);
	#endif
    
    #pragma omp parallel
    {

	    #pragma omp for schedule(dynamic)
		for(size_t i = 0; i < n; ++i) 
		{
			idx_t idx = vpNodeArray[i].array_idx;
			dpInfo[idx].ngbh = KNN_vpTree_V2(vpNodeArray[i].data, root, k, metric);
			dpInfo[idx].cluster_idx = -1;
			dpInfo[idx].is_center = 0;
			dpInfo[idx].array_idx = idx;
			
			#ifdef PROGRESS_BAR
				idx_t aa;
				#pragma omp atomic capture
				aa = ++progress_count;

				if(aa % step == 0 )
				{
					printf("Progress %lu/%lu -> %u%%\r",(uint64_t)aa, (uint64_t)n, (uint32_t)((100*aa)/n) );
					fflush(stdout);
				}
			#endif
		}
	}
	

    
	
	#ifdef PROGRESS_BAR
		printf("Progress %lu/%lu -> 100%%\n",(uint64_t)n, (uint64_t)n);
	#endif
    //printf("Progress %lu/%lu\n",(uint64_t)progress_count, (uint64_t)n);

    clock_gettime(CLOCK_MONOTONIC, &finish_tot);
    elapsed_tot = (finish_tot.tv_sec - start_tot.tv_sec);
    elapsed_tot += (finish_tot.tv_nsec - start_tot.tv_nsec) / 1000000000.0;
    printf("\tTotal time: %.3lfs\n\n", elapsed_tot);
    return;
}

Datapoint_info* NgbhSearch_vptree_V2(void* data, size_t n, size_t byteSize, size_t dims, size_t k, float_t (*metric)(void *, void *))
{

    struct timespec start, finish;
    double elapsed;
	METRICS_DATADIMS = (uint32_t)dims;
	
	#ifdef VERBOSE
		clock_gettime(CLOCK_MONOTONIC, &start);
		printf("Building the vp tree v2\n");
	#endif

	#ifdef SWMEM
		void* dummy_data = malloc(byteSize*dims*n);
		memcpy(dummy_data,data,byteSize*dims*n);
		data = dummy_data;
	#endif

	vpTreeNodeV2* vpNodeArray = (vpTreeNodeV2*)malloc(n*sizeof(vpTreeNodeV2));
	initialize_vpTreeNode_array_V2(vpNodeArray, data, n, byteSize*dims);


	//vpTreeNodeV2* root = build_vpTree_V2(vpPtrArray, 0, n-1, NULL, metric);
	vpTreeNodeV2* root = build_vpTree_V2(vpNodeArray, 0, n-1, NULL, metric);
	#ifdef VERBOSE
		clock_gettime(CLOCK_MONOTONIC, &finish);
		elapsed = (finish.tv_sec - start.tv_sec);
		elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    	printf("\tTotal time: %.3lfs\n\n", elapsed);
	#endif

    Datapoint_info* points = (Datapoint_info*)malloc(n*sizeof(Datapoint_info));

    /**************
     * KNN search *
     **************/



	
	KNN_search_vpTree_V2(points, vpNodeArray, root, k, n, metric);
	
	#ifdef SWMEM
		free(data);
	#endif

	//printf("NODE STRAMBO: %lf %p %lu\n", vpNodeArray[1516].mu, vpNodeArray[1516].inside, vpNodeArray[1516].parent -> array_idx);
	#ifdef VOPT
		for(idx_t i = 0; i < n; ++i) if(vpNodeArray[i].nodeList.indexes) free(vpNodeArray[i].nodeList.indexes);
	#else
		for(idx_t i = 0; i < n; ++i) if(vpNodeArray[i].nodeList.data) free(vpNodeArray[i].nodeList.data);
	#endif

    free(vpNodeArray);
	return points;
}

void setRhoErrK(Datapoint_info* points, FLOAT_TYPE* rho, FLOAT_TYPE* rhoErr, idx_t* k, size_t n)
{
	for(size_t i = 0; i < n; ++i)
	{
		points[i].log_rho = rho[i];
		points[i].log_rho_err = rhoErr[i];
		points[i].g = points[i].log_rho - points[i].log_rho_err;
		points[i].kstar = k[i];
	}
	return;
}

void computeAvg(Datapoint_info* p, FLOAT_TYPE *va, FLOAT_TYPE* ve, FLOAT_TYPE* vals, FLOAT_TYPE* verr, size_t k, size_t n)
{
	#pragma omp parallel for
	for(size_t i = 0; i < n; ++i)
	{
		Heap H = p[i].ngbh;
		float_t v_acc = 0.;
		float_t err_acc = 0.;
		for(size_t j = 0; j < k; ++j)
		{
			idx_t ngbhIdx =  H.data[j].array_idx; 
			v_acc += vals[ngbhIdx]; 
			err_acc += verr[ngbhIdx]; 
		}
		va[i] = v_acc/(float_t)k;
		ve[i] = err_acc/(float_t)k;
	}
		
}


Datapoint_info* allocDatapoints(idx_t n)
{
	Datapoint_info* points = (Datapoint_info*)malloc(sizeof(Datapoint_info)*n);
	for(idx_t p = 0; p < n; ++p)
	{
		points[p].array_idx = p;
		points[p].log_rho = 0.;
		points[p].log_rho_c = 0.;
		points[p].log_rho_err = 0.;
		points[p].g = 0.;
		points[p].kstar = 0;
	}
	return points;
}

void importNeighborsAndDistances(Datapoint_info* points, idx_t* indeces, float_t* distances, idx_t n, idx_t k)
{
	for(idx_t p = 0; p < n; ++p)
	{
		Heap H;
		allocateHeap(&H, k);
		initHeap(&H);
		for(idx_t j = 0; j < k; ++j) 
		{
	//		if(j < 10 && p < 10) printf("%lu %lu -> %lu;  \t",p,k, indeces[p*k + j]);
			H.data[j].value = distances[p*k + j]*distances[p*k + j];
			H.data[j].array_idx = indeces[p*k + j];
		}
	//	if(p < 10) printf("\n");
		H.count = k;
		points[p].ngbh = H;
	}
}

void importDensity(Datapoint_info* points, idx_t* kstar, float_t* density, float_t* density_err, idx_t n)
{
	for(idx_t p = 0; p < n; ++p)
	{
		points[p].log_rho 	  = density[p];
		points[p].log_rho_err = density_err[p];
		points[p].g 		  = density[p] - density_err[p]; 
		points[p].kstar		  = kstar[p]; 
	}
}

void exportDensity(Datapoint_info* points, idx_t* kstar, float_t* density, float_t* density_err, idx_t n)
{
	for(idx_t p = 0; p < n; ++p)
	{
		density[p] 		= 	points[p].log_rho;
		density_err[p]  =  	points[p].log_rho_err;
		kstar[p]        =  	points[p].kstar;		  
	}
}

void exportNeighborsAndDistances(Datapoint_info* points, idx_t* dist_indices, float_t* dists, idx_t n, idx_t k)
{
	for(idx_t i = 0; i < n; ++i)
	{
		for(idx_t j = 0; j < k; ++j)
		{
			dists[i*k + j] 		 	= points[i].ngbh.data[j].value;
			dist_indices[i*k + j] 	= points[i].ngbh.data[j].array_idx;
		}
	}
}

void exportClusterAssignment(Datapoint_info* points, int* labels, idx_t n)
{
	for(idx_t i = 0; i < n; ++i) labels[i] = points[i].cluster_idx;
}

void exportBorders(Clusters* clusters, int* border_idx, float_t* border_den, float_t* border_err)
{
	idx_t nclus = clusters -> centers.count; 
	if(clusters->UseSparseBorders)		
	{
		for(idx_t i = 0; i < nclus; ++i)		
			for(idx_t el = 0; el < clusters -> SparseBorders[i].count; ++el)
			{
				idx_t j = clusters -> SparseBorders[i].data[el].i;
				idx_t p = i*nclus + j;
				border_idx[p] = (int)j;  
				border_den[p] = clusters -> SparseBorders[i].data[el].density; 
				border_err[p] = clusters -> SparseBorders[i].data[el].error;
			}
	}
	else
	{
		for(idx_t i = 0; i < nclus; ++i)
			for(idx_t j = 0; j < nclus; ++j)
			{
				idx_t p = i*nclus + j;
				border_idx[p] = (int)clusters -> __borders_data[p].idx;
				border_den[p] = clusters -> __borders_data[p].density;
				border_err[p] = clusters -> __borders_data[p].error;
			}
	}
}
