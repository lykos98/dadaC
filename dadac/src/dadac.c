//#include "../include/read_fof_snapshot.h"
#include "../include/dadac.h"
#include <math.h>
#include <omp.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "./metrics.c"

#define MIN(x,y) x < y ? x : y
#define MAX(x,y) x > y ? x : y
#define DEPS 2.220446049250313e-16

int verbose = VERBOSE_TRUE;


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
const sparse_border_t sparse_border_null = {.density = -1.0, .error = 0, .idx = NOBORDER, .i = NOBORDER, .j = NOBORDER};


int blas_are_in_use()
{
	/*
	 * Helper function for retrieving if the code is compiled with
	 * a BLAS implementation.
	 * Usefull only for integration in python
	 */
	#ifdef USE_BLAS
		return 1;
	#else
		return 0;
	#endif
}

/* Clustering part */

void clusters_allocate(clusters * c, int s)
{
	/*
	 * Helper function for handling allocation of resources 
	 */ 
    if(c -> centers.count == 0)
    {
        printf("Provide a valid cluster centers list\n");
        return;
    }

    idx_t nclus = c -> centers.count;
    
    if(s)
    {
	    //printf("Using sparse implementation\n");
	    c -> use_sparse_borders = 1;
	    c -> sparse_borders = (adj_list_t*)malloc(nclus*sizeof(adj_list_t));
	    for(idx_t i = 0; i < nclus; ++i)
	    {
		    c -> sparse_borders[i].count = 0;
		    c -> sparse_borders[i].size  = PREALLOC_BORDERS;
		    c -> sparse_borders[i].data  = (sparse_border_t*)malloc(PREALLOC_BORDERS*sizeof(sparse_border_t));
	    }

    }
    else
    {
	    //printf("Using dense implementation\n");
	    c -> use_sparse_borders = 0;
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


void adj_list_insert(adj_list_t* l, sparse_border_t b)
{
	/*
	 * Handling of sparse border implementation as an adjecency list
	 */
	if(l -> count < l -> size)
	{
		l -> data[l -> count] = b;
		l -> count++;
	}
	else
	{
		l -> size += PREALLOC_BORDERS; 
		l -> data = realloc( l -> data, sizeof(sparse_border_t) * ( l -> size));
		l -> data[l -> count] = b;
		l -> count++;
	}
}

void adj_list_reset(adj_list_t* l)
{
	/*
	 * Handling of sparse border implementation as an adjecency list
	 */
	if(l -> data) free(l -> data);
	l -> count = 0;
	l -> size  = 0;
	l -> data  = NULL;
}

void clusters_reset(clusters * c)
{
	/* 
	 * Handling reset of clusters object 
	 */
	if(c -> use_sparse_borders)
	{
		for(idx_t i = 0; i < c -> centers.count; ++i)
		{
			adj_list_reset((c -> sparse_borders) + i);
		
		}
		free(c -> sparse_borders);
		c -> sparse_borders = NULL;
	}
	else
	{
		if(c -> __borders_data)  free(c -> __borders_data);
		if(c -> borders) free(c -> borders);
	}
	if(c -> centers.data) free(c -> centers.data);
}

void clusters_free(clusters * c)
{
	/*
	 * Free cluster object
	 */
    clusters_reset(c);
}


void sparse_border_insert(clusters *c, sparse_border_t b)
{
	/*
	 * Insert a border element in the sparse implementation
	 */

	idx_t i = b.i;
	adj_list_t l = c -> sparse_borders[i];
	int check = 1;
	for(idx_t k = 0; k < l.count; ++k)
	{
		sparse_border_t p = l.data[k];
		if(p.i == b.i && p.j == b.j)
		{
			if( b.density > p.density)
			{
				l.data[k] = b;
			}
			check = 0;
		}
	}
	if(check) adj_list_insert(c -> sparse_borders + i, b);
	return;
}

sparse_border_t sparse_border_get(clusters* c, idx_t i, idx_t j)
{
	/*
	 * Get a border element in the sparse implementation
	 * - i,j: cluster to search for borders
	 * return border_null if not found
	 */

	sparse_border_t b = sparse_border_null;
	adj_list_t l = c -> sparse_borders[i];
	for(idx_t el = 0; el < l.count; ++el)
	{
		sparse_border_t candidate = l.data[el];
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

void lu_dynamic_array_allocate(lu_dynamic_array * a)
{
    a -> data = (idx_t*)malloc(ARRAY_INCREMENT*sizeof(idx_t));
    a -> count = 0;
    a -> size = ARRAY_INCREMENT;
}

void lu_dynamic_array_pushBack(lu_dynamic_array * a, idx_t p)
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

void lu_dynamic_array_Reset(lu_dynamic_array * a)
{
    a -> count = 0;
}

void lu_dynamic_array_reserve(lu_dynamic_array * a, idx_t n)
{
    a -> data = realloc(a -> data, n*sizeof(idx_t));
    a -> size = n;
}

void lu_dynamic_array_init(lu_dynamic_array * a)
{
    a -> data = NULL;
    a -> count = 0;
    a -> size = 0;
}

int cmp(const void * a, const void * b)
{
    float_t aa = *((float_t*)a);
    float_t bb = *((float_t*)b);
    return (aa > bb ) - (aa < bb); 
}



float_t avg(const float_t * x, const idx_t n)
{
    float_t f = 0;
    for(idx_t i = 0; i < n; ++i)
    {
        f += x[i];
    }
    return f/(float_t)n;
}

/*
 * Density estimation
 */

float_t mEst2(float_t * x, float_t *y, idx_t n)
{

    /*
     * Estimate the m coefficient of a straight 
     * line passing through the origin          
     * params:                                  
     * - x: x values of the points              
     * - y: y values of the points              
     * - n: size of the arrays                  
     */
     
    float_t num = 0;
    float_t den = 0;
    float_t dd;
    for(idx_t i = 0; i < n; ++i)
    {
        float_t xx = x[i];
        float_t yy = y[i];

        dd = xx;
        num += dd*yy;
        den += dd*dd;

    }
  
    return num/den;
}
float_t mEst(float_t * x, float_t *y, idx_t n)
{
    float_t x_avg, y_avg;
    x_avg = avg(x,n);
    y_avg = avg(y,n);
    float_t num = 0;
    float_t den = 0;
    float_t dd;
    for(idx_t i = 0; i < n - 1; ++i)
    {
        float_t xx = x[i];
        float_t yy = y[i];

        dd = (xx - x_avg);
        num += dd*(yy - y_avg);
        den += dd*dd;

    }
  
    return num/den;
}

float_t id_estimate(datapoint_info* dp_info, idx_t n,float_t fraction, int verbose)
{

    /*
     * Estimation of the intrinsic dimension of a dataset                                       
     * args:                                                                                    
     * - dp_info: array of structs                                                             
     * - n: number of dp_info                                                                  
     * Estimates the id via 2NN method. Computation of the log ratio of the                      
     * distances of the first 2 neighbors of each point. Then compute the empirical distribution 
     * of these log ratii                                                                        
     * The intrinsic dimension is found by fitting with a straight line passing through the      
     * origin                                                                                    
     */

    struct timespec start_tot, finish_tot;
    double elapsed_tot;

	if(verbose)
	{
		printf("ID estimation:\n");
		clock_gettime(CLOCK_MONOTONIC, &start_tot);
	}

    //float_t fraction = 0.7;
    float_t* r = (float_t*)malloc(n*sizeof(float_t));
    float_t* Pemp = (float_t*)malloc(n*sizeof(float_t));

    for(idx_t i = 0; i < n; ++i)
    {
        r[i] = 0.5 * log(dp_info[i].ngbh.data[2].value/dp_info[i].ngbh.data[1].value);
        Pemp[i] = -log(1 - (float_t)(i + 1)/(float_t)n);
    }
    qsort(r,n,sizeof(float_t),cmp);

    idx_t Neff = (idx_t)(n*fraction);

    float_t d = mEst2(r,Pemp,Neff); 
    free(r);
    free(Pemp);

	if(verbose)
	{
		clock_gettime(CLOCK_MONOTONIC, &finish_tot);
		elapsed_tot = (finish_tot.tv_sec - start_tot.tv_sec);
		elapsed_tot += (finish_tot.tv_nsec - start_tot.tv_nsec) / 1000000000.0;
		printf("\tID value: %.6lf\n", d);
		printf("\tTotal time: %.3lfs\n\n", elapsed_tot);
	}

    return d;

}

void compute_density_kstarnn(datapoint_info* dp_info, const float_t d, const idx_t points, int verbose){

    /*
     * Point density computation:                       
     * args:                                            
     * - paricles: array of structs                   
     * - d       : intrinsic dimension of the dataset 
     * - points  : number of points in the dataset    
     */

    struct timespec start_tot, finish_tot;
    double elapsed_tot;

	if(verbose)
	{
		printf("Density and k* estimation:\n");
		clock_gettime(CLOCK_MONOTONIC, &start_tot);
	}

    idx_t kMAX = dp_info[0].ngbh.N - 1;   

    float_t omega = 0.;  
    if(sizeof(float_t) == sizeof(float)){ omega = powf(PI_F,d/2)/tgammaf(d/2.0f + 1.0f);}  
    else{omega = pow(M_PI,d/2.)/tgamma(d/2.0 + 1.0);}

    #pragma omp parallel for
    for(idx_t i = 0; i < points; ++i)
    {

        idx_t j = 4;
        idx_t k;
        float_t dL  = 0.;
        float_t vvi = 0.;
		float_t vvj = 0.;
		float_t vp  = 0.;
        while(j < kMAX  && dL < DTHR)
        {
            idx_t ksel = j - 1;
            vvi = omega * pow(dp_info[i].ngbh.data[ksel].value,d/2.);

            idx_t jj = dp_info[i].ngbh.data[j].array_idx;

            vvj = omega * pow(dp_info[jj].ngbh.data[ksel].value,d/2.);

            vp = (vvi + vvj)*(vvi + vvj);
            dL = -2.0 * ksel * log(4.*vvi*vvj/vp);
            j = j + 1;
        }
        if(j == kMAX)
        {
            k = j - 1;
            vvi = omega * pow(dp_info[i].ngbh.data[k].value,d/2.);
        }
        else
        {
            k = j - 2;
        }
        dp_info[i].kstar = k;
        dp_info[i].log_rho = log((float_t)(k)/vvi/((float_t)(points)));
        //dp_info[i].log_rho = log((float_t)(k)) - log(vvi) -log((float_t)(points));
        dp_info[i].log_rho_err =   1.0/sqrt((float_t)k); //(float_t)(-Q_rsqrt((float)k));
        dp_info[i].g = dp_info[i].log_rho - dp_info[i].log_rho_err;
    }

	if(verbose)
	{
		clock_gettime(CLOCK_MONOTONIC, &finish_tot);
		elapsed_tot = (finish_tot.tv_sec - start_tot.tv_sec);
		elapsed_tot += (finish_tot.tv_nsec - start_tot.tv_nsec) / 1000000000.0;
		printf("\tTotal time: %.3lfs\n\n", elapsed_tot);
	}

    return;


}

void PAk(datapoint_info* dp_info, const float_t d, const idx_t n, int verbose)
{
	/*
	 * /!\ Work in progress
	 * Implementation of the PAk estimator (kstarnn + linear correction)
	 * from the paper https://pubs.acs.org/doi/10.1021/acs.jctc.7b00916 
	 *
	 * args:
	 * - datapoint_info* dp_info : array of objects to store density to 
	 * - const float_t d 		: intrinsic dimension estimator  
	 * - const idx_t n 	  	  	: number of points in the dataset
	 */

    float_t omega = pow(M_PI,d/2.)/tgamma(d/2.0 + 1.0);
	compute_density_kstarnn(dp_info,d,n,verbose);

	#pragma omp parallel
	{
		float_t* vi = (float_t*)malloc(dp_info[0].ngbh.count*sizeof(float_t));
		#pragma omp for
		for(idx_t i = 0; i < n; ++i)
		{
			vi[0] = 0.;
			idx_t kstar = dp_info[i].kstar;
			for(idx_t k = 0; k < kstar; ++k)
			{
				//vi[k] = omega * pow(dp_info[i].ngbh.data[k].value/dp_info[i].ngbh.data[k+1].value,d/2.);
				vi[k] = omega * (pow(dp_info[i].ngbh.data[k+1].value,d/2.) - pow(dp_info[i].ngbh.data[k].value,d/2.));
			}

			float_t f = dp_info[i].log_rho + log((double)n);
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
			dp_info[i].log_rho =  f - log((double)n);
			dp_info[i].log_rho_err = sqrt((float_t)(4*kstar + 2)/(float_t)((kstar-1)*kstar)); 
		}
	free(vi);

	}
		
}


int cmpPP(const void* p1, const void *p2)
{
    /*
     * Utility function to perform quicksort   
     * when clustering assignment is performed    
     */
    datapoint_info* pp1 = *(datapoint_info**)p1;
    datapoint_info* pp2 = *(datapoint_info**)p2;
    return 2*(pp1 -> g < pp2 -> g) - 1;
}

void compute_correction(datapoint_info* dp_info, idx_t n, float_t Z)
{
    /*
     * Utility function, find the minimum value of the density of the datapoints
     * and shift them up in order to further work with values greater than 0     
     */
    float_t min_log_rho = 999999.9;
    

    #pragma omp parallel
    {
        float_t thread_min_log_rho = 9999999.;
        #pragma omp for
        for(idx_t i = 0; i < n; ++i)
        {
            float_t tmp = dp_info[i].log_rho - Z*dp_info[i].log_rho_err;
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
            dp_info[i].log_rho_c = dp_info[i].log_rho - min_log_rho + 1;
            dp_info[i].g = dp_info[i].log_rho_c - dp_info[i].log_rho_err;
        }
    }
    //printf("%lf\n",min_log_rho);
}

clusters Heuristic1(datapoint_info* dp_info, idx_t n, int verbose)
{
    /*
     * Heurisitc 1, from paper of Errico, Facco, Laio & Rodriguez 
     * ( https://doi.org/10.1016/j.ins.2021.01.010 )              
     *                                                            
     * args:                                                      
     * - dp_info: array of Datapoint structures                 
     * - data: pointer to the dataset                             
     * - n: number of Datapoints                                  
     */

    struct timespec start_tot, finish_tot;
    double elapsed_tot;

	if(verbose)
	{
		printf("H1: Preliminary cluster assignment\n");
		clock_gettime(CLOCK_MONOTONIC, &start_tot);
	}

    //idx_t ncenters = 0;
    //idx_t putativeCenters = n;
    lu_dynamic_array all_centers, removed_centers, actual_centers, max_rho;
    lu_dynamic_array_allocate(&all_centers);
    lu_dynamic_array_allocate(&removed_centers);
    lu_dynamic_array_allocate(&actual_centers);
    lu_dynamic_array_allocate(&max_rho);

    datapoint_info** dp_info_ptrs = (datapoint_info**)malloc(n*sizeof(datapoint_info*));

    struct timespec start, finish;
    double elapsed;


    if(verbose)
	{
		clock_gettime(CLOCK_MONOTONIC, &start);
	}


    for(idx_t i = 0; i < n; ++i)
    {   
        /*
         * Find the centers of the clusters as the points of higher density in their neighborhoods
         * A point is tagged as a putative center if it is the point of higer density of its neighborhood 
         */

        dp_info_ptrs[i] = dp_info + i;
        idx_t maxk = dp_info[i].kstar + 1;
        float_t gi = dp_info[i].g;
        dp_info[i].is_center = 1;
        dp_info[i].cluster_idx = -1;
        //printf("%lf\n",p -> g);
        heap i_ngbh = dp_info[i].ngbh;
        for(idx_t k = 1; k < maxk; ++k)
        {
            idx_t ngbh_index = i_ngbh.data[k].array_idx;
            float_t gj = dp_info[ngbh_index].g;
            if(gj > gi){
                dp_info[i].is_center = 0;
                break;
            }
        }
        if(dp_info[i].is_center){
                lu_dynamic_array_pushBack(&all_centers, i);
        }


    }

    if(verbose)
	{
		clock_gettime(CLOCK_MONOTONIC, &finish);
		elapsed = (finish.tv_sec - start.tv_sec);
		elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
		printf("\tFinding putative centers: %.3lfs\n",elapsed);
		clock_gettime(CLOCK_MONOTONIC, &start);
	}

	qsort(dp_info_ptrs, n, sizeof(datapoint_info*), cmpPP);


	#ifdef NO_OPT_H1
	
	/* 
	 * debug thing left here, implemented an aggressive optimization
	 * in future to be removed
	 */
	
    idx_t * to_remove = (idx_t*)malloc(all_centers.count*sizeof(idx_t));
    for(idx_t c = 0; c < all_centers.count; ++c) {to_remove[c] = MY_SIZE_MAX;}



    #pragma omp parallel
    {
            
        idx_t * to_remove_private = (idx_t*)malloc(all_centers.count*sizeof(idx_t));
    	for(idx_t c = 0; c < all_centers.count; ++c) {to_remove_private[c] = MY_SIZE_MAX;}

        #pragma omp for
        for(idx_t p = 0; p < n; ++p)
        {
        	datapoint_info pp = *(dp_info_ptrs[p]);
        	for(idx_t j = 1; j < pp.kstar + 1; ++j)
        	{
        		idx_t jidx = pp.ngbh.data[j].array_idx;
        		if(dp_info[jidx].is_center && pp.g > dp_info[jidx].g)
        		{
        			//dp_info[jidx].is_center = 0;
        			for(idx_t c = 0; c < all_centers.count; ++c){
        				if(all_centers.data[c] == jidx)
					{

						if(to_remove_private[c] != MY_SIZE_MAX)
						{
							to_remove_private[c] = pp.g > 	dp_info[to_remove_private[c]].g  ? pp.array_idx : to_remove_private[c];
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
        	for(idx_t c = 0; c < all_centers.count; ++c)
        	{
        		if(to_remove_private[c] != MY_SIZE_MAX)
			{
				if(to_remove[c] != MY_SIZE_MAX)
				{
					to_remove[c] = dp_info[to_remove_private[c]].g > dp_info[to_remove[c]].g ?
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

	

    for(idx_t p = 0; p < all_centers.count; ++p)
    {
        idx_t i = all_centers.data[p];
        int e = 0;
        //float_t gi = dp_info[i].g;
        idx_t mr = to_remove[p];
        if(mr != MY_SIZE_MAX)
        {
            //if(dp_info[mr].g > gi) e = 1;
	    e = 1;
        }
        switch (e)
        {
            case 1:
                {
                    lu_dynamic_array_pushBack(&removed_centers,i);
                    dp_info[i].is_center = 0;
                    //for(idx_t c = 0; c < removed_centers.count - 1; ++c)
                    //{
                    //    if(mr == removed_centers.data[c])
                    //    {
                    //        mr = max_rho.data[c];
                    //    }
                    //}
                    lu_dynamic_array_pushBack(&max_rho,mr);
                    
                }
                break;
            case 0:
                {
                    lu_dynamic_array_pushBack(&actual_centers,i);
                    dp_info[i].cluster_idx = actual_centers.count - 1;
                }
                break;
            default:
                break;
        }
    }
	free(to_remove);

	#else	

	/* 
	 * optimized version
	 *
	 * Generate a mask that keeps track of the point has been eliminating the 
	 * point considered. Each thread updates this mask, then after the procedure
	 * ends, center, removed centers, and max_rho arrays are populated
	 */
		
	idx_t* to_remove_mask = (idx_t*)malloc(n*sizeof(idx_t));
    for(idx_t p = 0; p < n; ++p) {to_remove_mask[p] = MY_SIZE_MAX;}

	
    #pragma omp parallel shared(to_remove_mask)
    {
		#pragma omp for
		for(idx_t p = 0; p < n; ++p)
		{
			datapoint_info pp = *(dp_info_ptrs[p]);
			int flag = 0;
			idx_t ppp = 0;

			for(idx_t j = 1; j < pp.kstar + 1; ++j)
			{
				idx_t jidx = pp.ngbh.data[j].array_idx; 
				if(dp_info[jidx].is_center && pp.g > dp_info[jidx].g)
				{

					#pragma omp critical 
					{
						ppp = to_remove_mask[jidx];
						flag = ppp != MY_SIZE_MAX;							
						to_remove_mask[jidx] = flag ? (pp.g > dp_info[ppp].g ? pp.array_idx : ppp) : pp.array_idx; 
					}


					/* Alternative
					 * not sure if it works

					#pragma omp atomic read 
					ppp = to_remove_mask[jidx];

					flag = ppp != MY_SIZE_MAX;							

					#pragma omp atomic write
					to_remove_mask[jidx] = flag ? (pp.g > dp_info[ppp].g ? pp.array_idx : ppp) : pp.array_idx; 

					*/
				}
			}
		}
	}

	/* populate the usual arrays */
    for(idx_t p = 0; p < all_centers.count; ++p)
    {
        idx_t i = all_centers.data[p];
        int e = 0;
        //float_t gi = dp_info[i].g;
        idx_t mr = to_remove_mask[i];
        if(mr != MY_SIZE_MAX)
        {
            //if(dp_info[mr].g > gi) e = 1;
			e = 1;
        }
        switch (e)
        {
            case 1:
                {
                    lu_dynamic_array_pushBack(&removed_centers,i);
                    dp_info[i].is_center = 0;
                    for(idx_t c = 0; c < removed_centers.count - 1; ++c)
                    {
                        if(mr == removed_centers.data[c])
                        {
                            mr = max_rho.data[c];
                        }
                    }
                    lu_dynamic_array_pushBack(&max_rho,mr);
                    
                }
                break;

            case 0:
                {
                    lu_dynamic_array_pushBack(&actual_centers,i);
                    dp_info[i].cluster_idx = actual_centers.count - 1;
                }
                break;

            default:
                break;
        }
    }

	free(to_remove_mask);

	#endif

    if(verbose)
	{
		clock_gettime(CLOCK_MONOTONIC, &finish);
		elapsed = (finish.tv_sec - start.tv_sec);
		elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
		printf("\tFinding actual centers:   %.3lfs\n",elapsed);
		clock_gettime(CLOCK_MONOTONIC, &start);
	}


    /*
     * Sort all the dp_info based on g and then perform the cluster assignment
     * in asceding order                                                     
     * UPDATE: dp_info already sorted                                          
     */
                                                                                
    for(idx_t i = 0; i < n; ++i)
    {   
        datapoint_info* p = dp_info_ptrs[i];
        if(!(p -> is_center))
        {
            int cluster = -1;
            idx_t k = 0;
            idx_t p_idx;
            idx_t max_k = p -> ngbh.N;
            /*assign each particle at the same cluster as the nearest particle of higher density*/
            while( k < max_k - 1 && cluster == -1)
            {
                ++k;
                p_idx = p -> ngbh.data[k].array_idx;
                cluster = dp_info[p_idx].cluster_idx; 
            }

            //
            if(cluster == -1)
            {
                float_t gmax = -99999.;               
                idx_t gm_index = 0;
                for(idx_t k = 0; k < max_k; ++k)
                {
                    idx_t ngbh_index = p -> ngbh.data[k].array_idx;
                    for(idx_t m = 0; m < removed_centers.count; ++m)
                    {
                        float_t gcand = dp_info[max_rho.data[m]].g;
                        if(ngbh_index == removed_centers.data[m] && gcand > gmax)
                        {   
                            //printf("%lu -- %lu\n", ele, m);
                            gmax = gcand;
                            gm_index = max_rho.data[m];
                        }
                    }
                }

                cluster = dp_info[gm_index].cluster_idx;

            }
            p -> cluster_idx = cluster;

        }
    }

    if(verbose)
	{
		clock_gettime(CLOCK_MONOTONIC, &finish);
		elapsed = (finish.tv_sec - start.tv_sec);
		elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
		printf("\tTentative clustering:     %.3lfs\n",elapsed);
		clock_gettime(CLOCK_MONOTONIC, &start);
	}

    free(dp_info_ptrs);
    free(max_rho.data);
    free(removed_centers.data);
    free(all_centers.data);


    clusters c_all;
    c_all.centers = actual_centers;


    if(verbose)
	{
		clock_gettime(CLOCK_MONOTONIC, &finish);
		elapsed = (finish.tv_sec - start.tv_sec);
		elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
		printf("\tFinalizing clustering:    %.3lfs\n",elapsed);
		printf("\n");
	}

    clock_gettime(CLOCK_MONOTONIC, &finish_tot);
    elapsed_tot = (finish_tot.tv_sec - start_tot.tv_sec);
    elapsed_tot += (finish_tot.tv_nsec - start_tot.tv_nsec) / 1000000000.0;


	if(verbose)
	{
		printf("\tFound %lu clusters\n",(uint64_t)actual_centers.count);
		printf("\tTotal time: %.3lfs\n\n", elapsed_tot);
	}

    c_all.n = n;
    return c_all;
}

void Heuristic2(clusters* cluster, datapoint_info* dp_info, int verbose)
{
    /*
     * Heurisitc 2, from paper of Errico, Facco, Laio & Rodriguez 
     * ( https://doi.org/10.1016/j.ins.2021.01.010 )              
     *                                                            
     * args:                                                      
     * - clusters* cluster 		: object containing information about the clusters (borders, centers)
	 * - datapoint_info* dp_info : array of Datapoint structures                 
	 *
     */

    #define borders cluster->borders

    struct timespec start_tot, finish_tot;
    double elapsed_tot;
    idx_t n = cluster -> n;

	if(verbose)
	{
		printf("H2: Finding border points\n");
		clock_gettime(CLOCK_MONOTONIC, &start_tot);
	}


    idx_t nclus = cluster->centers.count; 
    idx_t max_k = dp_info[0].ngbh.N;

    for(idx_t i = 0; i < n; ++i)
    {
		idx_t pp = NOBORDER;
		/*loop over n neighbors*/
		int c = dp_info[i].cluster_idx;
		if(!dp_info[i].is_center)
		{
			for(idx_t k = 1; k < dp_info[i].kstar + 1; ++k)
			{
				/*index of the kth ngbh of n*/
				idx_t j = dp_info[i].ngbh.data[k].array_idx;
				pp = NOBORDER;
				/*Loop over kn neigbhours to find if n is the nearest*/
				/*if cluster of the particle in nbhg is c then check is neighborhood*/                                                
				if(dp_info[j].cluster_idx != c)
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
				idx_t pp_ngbh_idx = dp_info[pp].ngbh.data[k].array_idx;
				if(pp_ngbh_idx == i)
				{
					break;
				}
				if(dp_info[pp_ngbh_idx].cluster_idx == c)
				{
					pp = NOBORDER;
					break;
				}
			}
		}

		/*if it is the maximum one add it to the cluster*/
		if(pp != NOBORDER)
		{
			int ppc = dp_info[pp].cluster_idx;
			if(cluster -> use_sparse_borders)
			{
				//insert one and symmetric one
				sparse_border_t b = {.i = c, .j = ppc, .idx = i, .density = dp_info[i].g, .error = dp_info[i].log_rho_err}; 
				sparse_border_insert(cluster, b);
				//get symmetric border
				sparse_border_t bsym = {.i = ppc, .j = c, .idx = i, .density = dp_info[i].g, .error = dp_info[i].log_rho_err}; 
				sparse_border_insert(cluster, bsym);

			}
			else
			{
				if(dp_info[i].g > borders[c][ppc].density)
				{
					borders[c][ppc].density = dp_info[i].g;
					borders[ppc][c].density = dp_info[i].g;
					borders[c][ppc].idx = i;
					borders[ppc][c].idx = i;
				}
			}
		}

    }


    if(cluster -> use_sparse_borders)
    {
	    for(idx_t c = 0; c < nclus; ++c)
	    {
		    for(idx_t el = 0; el < cluster -> sparse_borders[c].count; ++el)
		    {
			    //fix border density, write log rho c
			    idx_t idx = cluster -> sparse_borders[c].data[el].idx; 
			    cluster -> sparse_borders[c].data[el].density = dp_info[idx].log_rho_c;
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

			borders[i][j].density = dp_info[p].log_rho_c;
			borders[j][i].density = dp_info[p].log_rho_c;

			borders[i][j].error = dp_info[p].log_rho_err;
			borders[j][i].error = dp_info[p].log_rho_err;
		    }
		}
	    }

	    for(idx_t i = 0; i < nclus; ++i)
	    {
		borders[i][i].density = -1.0;
		borders[i][i].error = 0.0;
	    }
    }

	if(verbose)
	{
		clock_gettime(CLOCK_MONOTONIC, &finish_tot);
		elapsed_tot = (finish_tot.tv_sec - start_tot.tv_sec);
		elapsed_tot += (finish_tot.tv_nsec - start_tot.tv_nsec) / 1000000000.0;
		printf("\tTotal time: %.3lfs\n\n", elapsed_tot);
	}

    return;
    #undef borders
   }



void merge_A_into_B(idx_t* who_amI, idx_t cluster_A, idx_t cluster_B, idx_t n)
{
	/*
	 * Utility function
	 * Performs correctino of the labels in the array that
	 * keep tracks of what cluster is after a merging
	 */

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
	/*
	 * Utility function 
	 * comparision between two merging
	 * elements
	 */
	
	float_t DensA = ((merge_t*)A)->density;
	float_t DensB = ((merge_t*)B)->density;

	return - ( DensA > DensB) + (DensA < DensB);
}


static inline int is_a_merging(  float_t dens1, float_t dens1_err,
								 float_t dens2, float_t dens2_err,
								 float_t dens_border, float_t dens_border_err,
								 float_t Z)
{
	/*
	 * dens1	   : the density of the particle that is the center of the first cluster
	 * dens2	   : the density of the particle that is the center of the second cluster
	 * dens_border : the density of the border btw the cluster 1 and the cluster 2
	 * border_err  : the errors on the densities
	 * Z     	   : the desired accuracy
	 */

	/* in the original code it was:
	 *
	 float_t a1 = dp_info[cluster->centers.data[i]].log_rho_c - border_density[i][j];
	 float_t a2 = dp_info[cluster->centers.data[j]].log_rho_c - border_density[i][j];

	 float_t e1 = Z*(dp_info[cluster->centers.data[i]].log_rho_err + border_err[i][j]);
	 float_t e2 = Z*(dp_info[cluster->centers.data[j]].log_rho_err + border_err[i][j]);
	 */

	float_t a1 = dens1 - dens_border;
	float_t a2 = dens2 - dens_border;

	float_t e1 = Z*(dens1_err + dens_border_err);
	float_t e2 = Z*(dens2_err + dens_border_err);

	return (a1 < e1 || a2 < e2);
}


static inline int merging_roles(  float_t dens1, float_t dens1_err,
								  float_t dens2, float_t dens2_err,
								  float_t dens_border, float_t dens_border_err )
{
	/*
	 * Utility function 
	 * Retrieve if cluster 1 is merged into 2 or
	 * vice versa
	 */
      
	float_t c1 = (dens1 - dens_border) / (dens1_err + dens_border_err); 
	float_t c2 = (dens2 - dens_border) / (dens2_err + dens_border_err);
	//printf("%.10lf %.10lf %d\n",c1,c2, c1 > c2);

	return ( c1 < c2 );     // if 1, this signal to swap 1 and 2
}

void fix_borders_A_into_B(idx_t A, idx_t B, border_t** borders, idx_t n)
{
	/*
	 * Dense border implementation
	 * - idx_t A 			: cluster A the one which goes into the other 
	 * - idx_t B 			: cluster B the one that recieves the merging
	 * - border_t** borders : whole border matrix
	 * - idx_t n 			: number of clusters
	 *
	 */

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

static inline void delete_adj_list_element(clusters * c, const idx_t list_idx, const idx_t el)
{
	/*
	 * Utility function
	 * Deletes an element into an adjecency list,
	 * representing the borders in the cluster topology
	 */

	idx_t count = c -> sparse_borders[list_idx].count;
	c -> sparse_borders[list_idx].data[el] = c -> sparse_borders[list_idx].data[count-1];
	c -> sparse_borders[list_idx].data[count-1] = sparse_border_null;
	c -> sparse_borders[list_idx].count -= 1;
}

void fix_sparse_borders_A_into_B(idx_t s,idx_t t,clusters* c)
{
	/*
	 * Handle borders after two clusters are merged
	 * - idx_t s 	 : source cluster, the one has to be merged
	 * - idx_t t 	 : target cluster, the one recieves the merge
	 * - clusters* c : object containing all the data 
	 *
	 * When s goes into t all the clusters which had a border with s now they must have
	 * a border with t. If t already has a border like that, the border with higher 
	 * density is kept
	 */
	
	{
		{
			for(idx_t el = 0; el < c -> sparse_borders[t].count; ++el)
			{
				sparse_border_t b = c -> sparse_borders[t].data[el];
				if(b.i == t && b.j == s)
				{
					//delete the border src trg
					delete_adj_list_element(c, t, el);
				}
			}
		}
		//find the border and delete it, other insert them in correct place
		for(idx_t el = 0; el < c -> sparse_borders[s].count; ++el)
		{
			sparse_border_t b = c -> sparse_borders[s].data[el];
		//	idx_t ii = b.i;
			if(b.j != t)
			{
				//insert these borders as trg -> j and j -> trg
				b.i = t;
				sparse_border_insert(c, b);
				sparse_border_t bsym = b;
				bsym.i = b.j;
				bsym.j = b.i;
				sparse_border_insert(c, bsym);
				for(idx_t dl = 0; dl < c -> sparse_borders[b.j].count; ++dl)
				{
					sparse_border_t b_del = c -> sparse_borders[b.j].data[dl];
					if(b_del.j == s)
					{
						//delete the border src trg
						delete_adj_list_element(c, b.j, dl);
					}
				}
						
			}
		}
		//clean up all borders
		//delete the src list
		{
			adj_list_reset((c->sparse_borders) + s);
		}
		//delete all borders containing src
	//	for(idx_t i = 0; i < nclus; ++i)
	//	{
	//		for(idx_t el = 0; el < c -> sparse_borders[i].count; ++el)
	//		{
	//			sparse_border_t b = c -> sparse_borders[i].data[el];
	//			if(b.j == s)
	//			{
	//				//delete the border src trg
	//				delete_adj_list_element(c, i, el);
	//			}
	//		}
	//			
	//	}
	}


}

void Heuristic3_sparse(clusters* cluster, datapoint_info* dp_info, float_t Z, int halo, int verbose)
{

    /*
     * Heurisitc 3, from paper of Errico, Facco, Laio & Rodriguez 
     * ( https://doi.org/10.1016/j.ins.2021.01.010 )              
	 * Dense implementation, makes use of a dense matrix to store
	 * borders between clusters, so it is more performant when the number of clusters is low
     *                                                            
     * args:                                                      
	 * - clusters* cluster 			: cluster object storing border info and cluster centers                 
     * - datapoint_info* dp_info 	: array of Datapoint structures                             
	 * - float_t Z 					: parameter to assess density peak significance
     * - halo 					    : flag to set if you want to compute also the halo points                               
     */

	#define borders cluster->borders

	struct timespec start_tot, finish_tot;
	double elapsed_tot;

	struct timespec start, finish;
	double elapsed;

	clock_gettime(CLOCK_MONOTONIC, &start_tot);

	if(verbose)
	{
		printf("H3: Merging clusters\n");
		printf("Using sparse implementation\n");
		clock_gettime(CLOCK_MONOTONIC, &start); 
	}

	idx_t nclus                 = cluster -> centers.count;  
	idx_t *  surviving_clusters = (idx_t*)malloc(nclus*sizeof(idx_t));
	for(idx_t i = 0; i < nclus; ++i)
	{ 
		surviving_clusters[i] = i; 
	}

	idx_t   merge_count        = 0;
	idx_t   merging_table_size = 1000;
	merge_t *merging_table      = (merge_t*)malloc(sizeof(merge_t)*merging_table_size);
  
	/*
	 * Find clusters to be merged
	 * Loop over borders and find which ones will generate a merge,
	 * store them later in the merging table
	 **/
	for(idx_t i = 0; i < nclus - 1; ++i)   
	{
		idx_t count = cluster -> sparse_borders[i].count;
		for(idx_t el = 0; el < count; ++el)   
		{
			sparse_border_t b = cluster -> sparse_borders[i].data[el];
			if( b.j > b.i)
			{
				float_t dens1           = dp_info[cluster->centers.data[b.i]].log_rho_c;
				float_t dens1_err       = dp_info[cluster->centers.data[b.i]].log_rho_err;
				float_t dens2           = dp_info[cluster->centers.data[b.j]].log_rho_c;
				float_t dens2_err       = dp_info[cluster->centers.data[b.j]].log_rho_err;
				float_t dens_border     = b.density;
				float_t dens_border_err = b.error;

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

	if(verbose)
	{
		clock_gettime(CLOCK_MONOTONIC, &finish); 
		elapsed = (finish.tv_sec - start.tv_sec);
		elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
		printf("\tFinding merges:   %.3lfs\n", elapsed);
		clock_gettime(CLOCK_MONOTONIC, &start); 
	}
  
  
    for( idx_t m = 0; m < merge_count; m++ )
    {
      
        #define src surviving_clusters[merging_table[m].source]
        #define trg surviving_clusters[merging_table[m].target]

        //int re_check = ( (src != merging_table[m].source) || (trg != merging_table[m].target) );
		//if(re_check)
		{
			/* 
			 * Enforce a that in case of symmetric merging condition the lowest idx cluster 
			 * is merged into the higher idx cluster, only to preserve compatibility with 
			 * original ADP implementation
			 *
			 * Process each element in the merging table
			 */
			idx_t new_src = (src < trg) ? src : trg;
			idx_t new_trg = (src < trg) ? trg : src;

			/*
			 * pick who am I and retrieve all needed data from the 
			 * border matrices
			 */

			float_t dens1           = dp_info[cluster->centers.data[new_src]].log_rho_c;
			float_t dens1_err       = dp_info[cluster->centers.data[new_src]].log_rho_err;
			float_t dens2           = dp_info[cluster->centers.data[new_trg]].log_rho_c;
			float_t dens2_err       = dp_info[cluster->centers.data[new_trg]].log_rho_err;

			//borders get
			sparse_border_t b 	   	= sparse_border_get(cluster, new_src, new_trg);
			float_t dens_border     = b.density;
			float_t dens_border_err = b.error;

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

						
					/* 
					 * Perform the actual meriging,
					 * first  -> fix the borders, delete old ones and spawn new one in the correct position
					 * second -> update the surviving_clusters buffer
					 */
					fix_sparse_borders_A_into_B(new_src, new_trg, cluster);
					merge_A_into_B(surviving_clusters, new_src, new_trg, nclus );	  
				}
				break;
			
			default:
				break;
			}
		}
        
        #undef src
        #undef trg
    }

	if(verbose)
	{
		clock_gettime(CLOCK_MONOTONIC, &finish); 
		elapsed = (finish.tv_sec - start.tv_sec);
		elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
		printf("\tCluster merging:  %.3lfs\n", elapsed);
		clock_gettime(CLOCK_MONOTONIC, &start); 
	}
  
    /*Finalize clustering*/
    /*Acutally copying */
    lu_dynamic_array tmp_centers;
    lu_dynamic_array tmp_cluster_idx;


    lu_dynamic_array_init(&tmp_centers);
    lu_dynamic_array_init(&tmp_cluster_idx);

    lu_dynamic_array_reserve(&tmp_centers, nclus);
    lu_dynamic_array_reserve(&tmp_cluster_idx, nclus);

    idx_t final_cluster_count = 0;

    idx_t* old_to_new = (idx_t*)malloc(nclus*sizeof(idx_t));
    idx_t incremental_k = 0;
    for(idx_t i = 0; i < nclus; ++i)
    {
        
        if(surviving_clusters[i] == i){
            lu_dynamic_array_pushBack(&tmp_centers, cluster->centers.data[i]);
            lu_dynamic_array_pushBack(&tmp_cluster_idx, i);
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

    adj_list_t* tmp_borders = (adj_list_t*)malloc(final_cluster_count*sizeof(adj_list_t));

    //initialize temporary borders
    for(idx_t i = 0; i < final_cluster_count; ++i)
    {
	    tmp_borders[i].count = 0;
	    tmp_borders[i].size  = PREALLOC_BORDERS;
	    tmp_borders[i].data  = (sparse_border_t*)malloc(PREALLOC_BORDERS*sizeof(sparse_border_t));
    }

    /*initialize all pointers*/

    /*Fix cluster assignment*/
    #pragma omp parallel for
    for(idx_t i = 0; i < cluster -> n; ++i)
    {
        dp_info[i].is_center = 0;
        int old_cidx = dp_info[i].cluster_idx;
        dp_info[i].cluster_idx = old_to_new[old_cidx];
    }

    
    #pragma omp parallel for
    for(idx_t c = 0; c < final_cluster_count; ++c)
    {
        idx_t c_idx = tmp_cluster_idx.data[c];
		for(idx_t el = 0; el < cluster -> sparse_borders[c_idx].count; ++el)
		{
			//retrieve border
			sparse_border_t b = cluster -> sparse_borders[c_idx].data[el];
			//change idexes of clusters
			b.i = old_to_new[b.i];
			b.j = old_to_new[b.j];

			adj_list_insert(tmp_borders + c, b);
		}
    }

    clusters_reset(cluster);
    /*pay attention to the defined borders*/
    /*copy into members*/
    cluster -> sparse_borders = tmp_borders;


    cluster -> centers = tmp_centers;
    /*
     * Fix center assignment
     */
    for(idx_t i = 0; i < cluster -> centers.count; ++i)
    {
        int idx = cluster -> centers.data[i];
        dp_info[idx].is_center = 1;
    }

    /*Halo*/
    switch (halo)
    {
		case 1:
		{
			float_t* max_border_den_array = (float_t*)malloc(final_cluster_count*sizeof(float_t));
			#pragma omp parallel
			{
				#pragma omp for
				for(idx_t c = 0; c < final_cluster_count; ++c)
				{
					float_t max_border_den = -2.;
					for(idx_t el = 0; el < cluster -> sparse_borders[c].count; ++el)
					{
						sparse_border_t b = cluster -> sparse_borders[c].data[el];
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
					int cidx = dp_info[i].cluster_idx;
					int halo_flag = dp_info[i].log_rho_c < max_border_den_array[cidx]; 
					//int halo_flag = max_border_den_array[cidx] > dp_info[i].log_rho_c  ; 
					dp_info[i].cluster_idx = halo_flag ? -1 : cidx;
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

	if(verbose)
	{
		clock_gettime(CLOCK_MONOTONIC, &finish); 
		elapsed = (finish.tv_sec - start.tv_sec);
		elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
		printf("\tFinal operations: %.3lfs\n\n", elapsed);

		clock_gettime(CLOCK_MONOTONIC, &finish_tot);
		elapsed_tot = (finish_tot.tv_sec - start_tot.tv_sec);
		elapsed_tot += (finish_tot.tv_nsec - start_tot.tv_nsec) / 1000000000.0;
		printf("\tFound %lu possible merges\n",(uint64_t)merge_count);
		printf("\tSurviving clusters %lu\n",(uint64_t)final_cluster_count);
		printf("\tTotal time: %.3lfs\n\n", elapsed_tot);
	}


  #undef  borders  
}


void Heuristic3_dense(clusters* cluster, datapoint_info* dp_info, float_t Z, int halo, int verbose)
{
    /*
     * Heurisitc 3, from paper of Errico, Facco, Laio & Rodriguez 
     * ( https://doi.org/10.1016/j.ins.2021.01.010 )              
	 * Dense implementation, makes use of a dense matrix to store
	 * borders between clusters, so it is more performant when the number of clusters is low
     *                                                            
     * args:                                                      
	 * - clusters* cluster 			: cluster object storing border info and cluster centers                 
     * - datapoint_info* dp_info 	: array of Datapoint structures                             
	 * - float_t Z 					: parameter to assess density peak significance
     * - halo 					    : flag to set if you want to compute also the halo points                               
     */

	#define borders cluster->borders

	struct timespec start_tot, finish_tot;
	double elapsed_tot;

	struct timespec start, finish;
	double elapsed;

	clock_gettime(CLOCK_MONOTONIC, &start_tot);

	if(verbose)
	{
		printf("H3: Merging clusters\n");
		printf("Using dense implementation\n");
		clock_gettime(CLOCK_MONOTONIC, &start); 
	}

	idx_t nclus              	= cluster -> centers.count;  
	idx_t *  surviving_clusters = (idx_t*)malloc(nclus*sizeof(idx_t));

	for(idx_t i = 0; i < nclus; ++i)
	{ 
		surviving_clusters[i] = i; 
	}

	idx_t   merge_count        = 0;
	idx_t   merging_table_size = 1000;
	merge_t *merging_table     = (merge_t*)malloc(sizeof(merge_t)*merging_table_size);
  
	/*
	 * Find clusters to be merged
	 * Loop over borders and find which ones will generate a merge,
	 * store them later in the merging table
	 */
	for(idx_t i = 0; i < nclus - 1; ++i)   
		for(idx_t j = i + 1; j < nclus; ++j)   
		{
			switch(borders[i][j].idx != NOBORDER)
			{
						
				case 1:		
				{
					float_t dens1           = dp_info[cluster->centers.data[i]].log_rho_c;
					float_t dens1_err       = dp_info[cluster->centers.data[i]].log_rho_err;
					float_t dens2           = dp_info[cluster->centers.data[j]].log_rho_c;
					float_t dens2_err       = dp_info[cluster->centers.data[j]].log_rho_err;
					float_t dens_border     = borders[i][j].density;
					float_t dens_border_err = borders[i][j].error;

					if ( is_a_merging( dens1, dens1_err, dens2, dens2_err, dens_border, dens_border_err, Z ) )
					{

						if ( merge_count == merging_table_size ) {
						merging_table_size *= 1.1;
						merging_table = (merge_t*)realloc( merging_table, sizeof(merge_t) * merging_table_size ); }

						idx_t src = j;
						idx_t trg = i;

						//int swap = merging_roles( dens1, dens1_err, dens2, dens2_err, dens_border, dens_border_err);
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

	/*
	 * sort the merging table
	 */

	qsort( (void*)merging_table, merge_count, sizeof(merge_t), compare_merging_density);

	if(verbose)
	{
		clock_gettime(CLOCK_MONOTONIC, &finish); 
		elapsed = (finish.tv_sec - start.tv_sec);
		elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
		printf("\tFinding merges:   %.3lfs\n", elapsed);
		clock_gettime(CLOCK_MONOTONIC, &start); 
	}
  
  
    for( idx_t m = 0; m < merge_count; m++ )
    {
        #define src surviving_clusters[merging_table[m].source]
        #define trg surviving_clusters[merging_table[m].target]

		/* Not needed anymore */

        //int re_check = ( (src != merging_table[m].source) || (trg != merging_table[m].target) );
		//if(re_check)
		{
			/* 
			 * Enforce a that in case of symmetric merging condition the lowest idx cluster 
			 * is merged into the higher idx cluster, only to preserve compatibility with 
			 * original ADP implementation
			 *
			 * Process each element in the merging table
			 */

			idx_t new_src = (src < trg) ? src : trg;
			idx_t new_trg = (src < trg) ? trg : src;

			/*
			 * pick who am I and retrieve all needed data from the 
			 * border matrices
			 */

			float_t dens1           = dp_info[cluster->centers.data[new_src]].log_rho_c;
			float_t dens1_err       = dp_info[cluster->centers.data[new_src]].log_rho_err;
			float_t dens2           = dp_info[cluster->centers.data[new_trg]].log_rho_c;
			float_t dens2_err       = dp_info[cluster->centers.data[new_trg]].log_rho_err;

			float_t dens_border     = borders[new_src][new_trg].density;
			float_t dens_border_err = borders[new_src][new_trg].error;

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
						
						/* 
						 * Perform the actual meriging,
						 * first  -> fix the borders, delete old ones and spawn new one in the correct position
						 * second -> update the surviving_clusters buffer
						 */
						fix_borders_A_into_B(new_src,new_trg,borders,nclus);
						merge_A_into_B ( surviving_clusters, new_src, new_trg, nclus );	  
					}
					break;
			
			default:
				break;
			}
		}
        
        #undef src
        #undef trg
    }

	if(verbose)
	{
		clock_gettime(CLOCK_MONOTONIC, &finish); 
		elapsed = (finish.tv_sec - start.tv_sec);
		elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
		printf("\tCluster merging:  %.3lfs\n", elapsed);
		clock_gettime(CLOCK_MONOTONIC, &start); 
	}
  
    /*Finalize clustering*/
    /*Acutally copying old matrices into new */
    lu_dynamic_array tmp_centers;
    lu_dynamic_array tmp_cluster_idx;


    lu_dynamic_array_init(&tmp_centers);
    lu_dynamic_array_init(&tmp_cluster_idx);

    lu_dynamic_array_reserve(&tmp_centers, nclus);
    lu_dynamic_array_reserve(&tmp_cluster_idx, nclus);

    idx_t final_cluster_count = 0;

    idx_t* old_to_new = (idx_t*)malloc(nclus*sizeof(idx_t));
    idx_t incremental_k = 0;

	/*
	 * Map old cluster labels to new
	 */

    for(idx_t i = 0; i < nclus; ++i)
    {
        
        if(surviving_clusters[i] == i){
            lu_dynamic_array_pushBack(&tmp_centers, cluster->centers.data[i]);
            lu_dynamic_array_pushBack(&tmp_cluster_idx, i);
            old_to_new[i] = incremental_k;
            ++incremental_k;
            ++final_cluster_count;
        }
    }

    /*fill the rest of old_to_new*/
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
        dp_info[i].is_center = 0;
        int old_cidx = dp_info[i].cluster_idx;
        dp_info[i].cluster_idx = old_to_new[old_cidx];
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

    clusters_reset(cluster);
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
        dp_info[idx].is_center = 1;
    }
    /*Halo*/
    switch (halo)
    {
		case 1:
		{
			float_t* max_border_den_array = (float_t*)malloc(final_cluster_count*sizeof(float_t));
			#pragma omp parallel
			{
				#pragma omp for
				for(idx_t c = 0; c < final_cluster_count; ++c)
				{
				float_t max_border_den = -2.;
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
					int cidx = dp_info[i].cluster_idx;
					int halo_flag = dp_info[i].log_rho_c < max_border_den_array[cidx];  
					//int halo_flag =  (max_border_den_array[cidx] > dp_info[i].log_rho_c) - (max_border_den_array[cidx] < dp_info[i].log_rho_c); 
					//halo_flag = halo_flag > 0;
					dp_info[i].cluster_idx = halo_flag ? -1 : cidx;
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

  if(verbose)
  {
	  clock_gettime(CLOCK_MONOTONIC, &finish); 
	  elapsed = (finish.tv_sec - start.tv_sec);
	  elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
	  printf("\tFinal operations: %.3lfs\n\n", elapsed);

	  clock_gettime(CLOCK_MONOTONIC, &finish_tot);
	  elapsed_tot = (finish_tot.tv_sec - start_tot.tv_sec);
	  elapsed_tot += (finish_tot.tv_nsec - start_tot.tv_nsec) / 1000000000.0;
	  printf("\tFound %lu possible merges\n", (uint64_t)merge_count);
	  printf("\tSurviving clusters %lu\n", (uint64_t)final_cluster_count);
	  printf("\tTotal time: %.3lfs\n\n", elapsed_tot);
  }


  #undef  borders  
}


void Heuristic3(clusters* cluster, datapoint_info* dp_info, float_t Z, int halo, int verbose)
{
	/*
	 * Wapper for dense and sparse implementation
	 */
	if(cluster -> use_sparse_borders)
	{
		Heuristic3_sparse(cluster, dp_info,  Z,  halo, verbose);
	}
	else
	{
		Heuristic3_dense(cluster, dp_info,  Z,  halo, verbose);
	}
}

void compute_level(kdnode* root, idx_t prev_lvl)
{
	idx_t curr_lvl = prev_lvl + 1;
	root -> level = curr_lvl;
	if(root -> lch) compute_level(root -> lch, curr_lvl);
	if(root -> rch) compute_level(root -> rch, curr_lvl);
	return;
}

/*
 * Other Helper functions
 */

void knn_search_kdtree_v2(datapoint_info * dp_info,kdnode_v2* kdNodeArray, kdnode_v2* root, idx_t n, idx_t k, int verbose)
{
	/*
	 * Helper function for knn serch using KDtree
	 *
	 * - datapoint_info * dp_info 	: array of datapoints information structs to store the neighborhood on 
	 * - float_t * data 			: array of data 
	 * - kdnode_v2* kdNodeArray  	: array of the nodes of the KDtree
	 * - kdnode_v2* root 			: root of the KDtree 
	 * - idx_t n 					: number of elements in the dataset 
	 * - idx_t k 					: number of neighbors to retrieve
	 */
    struct timespec start_tot, finish_tot;
    double elapsed_tot;

	if(verbose)
	{
		printf("knn search:\n");
		clock_gettime(CLOCK_MONOTONIC, &start_tot);
	}
	
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
			dp_info[idx].ngbh = knn_kdtree_v2(kdNodeArray[p].data, root, k);
			dp_info[idx].cluster_idx = -1;
			dp_info[idx].is_center = 0;
			dp_info[idx].array_idx = idx;

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

	if(verbose)
	{
		clock_gettime(CLOCK_MONOTONIC, &finish_tot);
		elapsed_tot = (finish_tot.tv_sec - start_tot.tv_sec);
		elapsed_tot += (finish_tot.tv_nsec - start_tot.tv_nsec) / 1000000000.0;
		printf("\tTotal time: %.3lfs\n\n", elapsed_tot);
	}
    return;

}


void knn_search_kdtree(datapoint_info * dp_info, float_t * data, kdnode* root, idx_t n, idx_t k, int verbose)
{
	/*
	 * Helper function for knn serch using KDtree
	 *
	 * - datapoint_info* dp_info : array of datapoints information structs to store the neighborhood on 
	 * - float_t * data 		: array of data 
	 * - kdnode* root 			: root of the KDtree 
	 * - idx_t n 				: number of elements in the dataset 
	 * - idx_t k 				: number of neighbors to retrieve
	 */
    struct timespec start_tot, finish_tot;
    double elapsed_tot;

	if(verbose)
	{
		printf("knn search:\n");
		clock_gettime(CLOCK_MONOTONIC, &start_tot);
	}
	
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
		dp_info[p].ngbh = knn(data + data_dims*p, root, k);
		dp_info[p].array_idx = p;

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

/*
 * Standard partition function 
 * and quick select method used for knn_bruteforce
 */

int partition_heap_node(heap_node *array, int left, int right, int pivot_index) {
    float_t pivot_value = array[pivot_index].value;
    int storeIndex = left;
    int i;
    /* Move pivot to end */
    swap_heap_node(array + pivot_index, array + right);
    for(i=left; i < right; i = i + 1 ){
        if(array[i].value < pivot_value){
    		swap_heap_node(array + storeIndex, array + i);
            storeIndex += 1;
        }
    }
    /* Move pivot to its final place */
    swap_heap_node(array + storeIndex , array + right);

    return storeIndex;
}

int qselect_heap_node(heap_node *array, int left, int right, int n) {
    int pivot_index;
    if(left == right){
        return left;
    }
    pivot_index = left + (rand() % (right-left + 1)); /* random int left <= x <= right */
    pivot_index = partition_heap_node(array, left, right, pivot_index);
    /* The pivot is in its final sorted position */
    if(n == pivot_index){
        return pivot_index;
    }else if(n < pivot_index){
        return qselect_heap_node(array, left, pivot_index-1, n);
    }else{
        return qselect_heap_node(array, pivot_index+1, right, n);
    }
}

int quickselect_heap_node(heap_node *array, int array_size, int k){
    return qselect_heap_node(array, 0, array_size-1, k-1);
}


void knn_BruteForce(datapoint_info* points, void* data,size_t n, size_t byteSize, size_t dims, size_t k, float_t (*metric)(void*, void*))
{
	/*
	 * /!\ work in progress
	 * knn function using brute force method, possibly using BLAS for euclidean
	 * metrics
	 *
	 * - datapoint_info* points : datapoints structs to store ngbh to 
	 * - void* data 			: dataset pointer
	 * - size_t byteSize 		: bytesize of each element in the vector of a row of the dataset
	 * - size_t dims 			: number of elements per row
	 * - size_t k 				: number of neighbors to retrieve
	 * - float_t (*metric)(void*, void*)) : callable metric function
	 *
	 */
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
				heap H = points[p].ngbh;

				quickselect_heap_node(pvt_working_mem, n,  k + 1);
				qsort(pvt_working_mem, k + 1, sizeof(heap_node), cmp_heap_nodes);
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

			for(idx_t j = 0; j < n; ++j)
			{
				pvt_working_mem[j].value = metric(data + p*dims*byteSize, data + j*dims*byteSize);
				pvt_working_mem[j].array_idx = j; 
			}

			heap H = points[p].ngbh;
			quickselect_heap_node(pvt_working_mem, n,  k + 1);
			qsort(pvt_working_mem, k, sizeof(heap_node), cmp_heap_nodes);
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
	void __handle_w_blas_d(datapoint_info* points, double* data,size_t n, size_t dims, size_t k)
	{
		/*
		 * A way to speed up (drammatically) euclidean metric computation, when dealing with high dimensional data
		 * The trick is to compute ||xi - xj||^2 as ||xi||^2 - 2*<xi|xj> + ||xj||^2 
		 * ||xi|| and ||xj|| can be computed upfront, then <xi|xj> can be computed as the product 
		 * of the matrix  X^T @ X, these terms are called _MIDDLE TERMS_
		 *
		 * Actually a slice of the dataset, Y, is used against the whole dataset and then the "all-to-all" distances are 
		 * reduced to the first k ones 
		 */
		struct sysinfo info;
		sysinfo(&info);

		/* 
		 * allocate a piece of memory to perform middle term computation
		 * estimate how much memory is free and take up to 90% for that
		 */

		size_t slice_size = (size_t)((float_t)(info.freeram)*0.9/(float_t)(8*n));
		slice_size = MIN(n,slice_size);

		printf("%lu slice size \n",slice_size);
		float_t* norms_sq = (float_t*)malloc(n*sizeof(double));

		float_t* middle_terms 	= (float_t*)malloc(slice_size*sizeof(float_t)*n); 

		/* compute the norms of the vectors */

		#pragma omp parallel for
		for(size_t i = 0; i < n; ++i)
		{
			norms_sq[i] = 0.;
			for(size_t j = 0; j < dims; ++j) norms_sq[i] += data[i*dims+j]*data[i*dims+j];
		}

		heap_node* working_mem 	= NULL;

		/* 
		 * allocate for each thread a piece of working memory to perform the reduction
		 * this memory should contain distance and index of the kth neighbor
		 * Recycle the heap for that purpose
		 */

		#pragma omp parallel
		{
			#pragma omp master
			{
				working_mem = (heap_node*)malloc(omp_get_num_threads()*sizeof(heap_node)*n); 
			}
		}


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

			#pragma omp parallel 
			{
				heap_node* my_working_mem = working_mem + omp_get_thread_num()*n;
				#pragma omp for
				for(size_t i = lower_idx; i < upper_idx; ++i)
				{
					/*alternative version, trying to understand what is faster*/

					//for(size_t j = 0; j < n; ++j)
					//{
					// 		my_working_mem[j].value = middle_terms[(i%DEFAULT_SLICE)*n + j];
					//		my_working_mem[j].array_idx = j;
					//		insertMaxheap(&(points[i].ngbh), middle_terms[(i%DEFAULT_SLICE)*n + j],j);
					//}
					//heapSort(&(points[i].ngbh));

					/* copy into the working mem */
					for(size_t j = 0; j < n; ++j)
					{
						my_working_mem[j].value = middle_terms[(i%slice_size)*n + j];
						my_working_mem[j].array_idx = j;
					}

					heap H = points[i].ngbh;
					quickselect_heap_node(my_working_mem, n,  k + 1);
					qsort(my_working_mem, k, sizeof(heap_node), cmp_heap_nodes);
					memcpy(H.data, my_working_mem, k*sizeof(heap_node));
				}
			}
		}

		free(working_mem);
		free(middle_terms);

		free(norms_sq);	

	}
#endif

datapoint_info* ngbh_search_bruteforce(void* data, size_t n, size_t byteSize, size_t dims, size_t k, float_t (*metric)(void *, void *), int verbose)
{
	/*
	 * Neighborhood search using bruteforce method
	 *
	 * - void* data 	 : the dataset, in general a collection of elements					 
	 * - size_t n 		 : number of points in the dataset
	 * - size_t byteSize : number of bytes for each element in the vector representing a datapooint 
	 * - size_t dims 	 : lenght of each datapoint in terms of elements 
	 * - size_t k 		 : number of neighbors to search for 
	 * - float_t (*metric)(void *, void *)) : distance function between two elements in the dataset MUST satisy triangle inequality 
	 *
	 */
	METRICS_DATADIMS = dims;

	if(verbose)
	{
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
	}
    datapoint_info* points = (datapoint_info*)malloc(n*sizeof(datapoint_info));

    struct timespec start_tot, finish_tot;
    double elapsed_tot;
    printf("knn search:\n");
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
		heap H;
		allocate_heap(&H, k);
		init_heap(&H);
		points[j].ngbh = H;
		points[j].array_idx = j;

	}

	if(metric)
	{
		knn_BruteForce(points, data, n, byteSize, dims, k, metric);
	}
	else
	{
		#ifdef USE_BLAS
			__handle_w_blas_d(points, data, n, dims, k);
		#else
			knn_BruteForce(points, data, n, byteSize, dims, k, eud);
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

datapoint_info* ngbh_search_kdtree(float_t* data, size_t n, size_t ndims, size_t k, int verbose)
{
	/*
	 * Neighborhood search using a KDtree
	 *
	 * - float_t* data 	: dataset as a matrix of float or double
	 * - size_t n 		: number of elements in the dataset 
	 * - size_t ndims 	: lenght of each vector in the dataset
	 * - size_t k 		: number of neighbors to retrieve
	 *
	 */
    struct timespec start, finish;
    double elapsed;


	data_dims = (unsigned int)ndims;

	if(verbose)
	{
		printf("Building the KDtree v1:\n");
		clock_gettime(CLOCK_MONOTONIC, &start);
	}

    kdnode* kdnode_array = (kdnode*)malloc(n*sizeof(kdnode));
    kdnode** kd_ptrs = (kdnode**)malloc(n*sizeof(kdnode*));

    initialize_kdnodes(kdnode_array,data,n);
    initialize_kdnode_ptrs(kd_ptrs, kdnode_array,n);

    kdnode* root = build_tree(kd_ptrs, n, ndims);

	if(verbose)
	{
		clock_gettime(CLOCK_MONOTONIC, &finish);
		elapsed = (finish.tv_sec - start.tv_sec);
		elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
		printf("\tTotal time: %.3lfs\n\n", elapsed);
	}


    datapoint_info* points = (datapoint_info*)malloc(n*sizeof(datapoint_info));

    /**************
     * knn search *
     **************/
    knn_search_kdtree(points,data, root, n, k, verbose);

    free(kd_ptrs);
    free(kdnode_array);

	return points;
}

datapoint_info* ngbh_search_kdtree_v2(float_t* data, size_t n, size_t ndims, size_t k, int verbose)
{
	/*
	 * Neighborhood search using a KDtree
	 *
	 * - float_t* data 	: dataset as a matrix of float or double
	 * - size_t n 		: number of elements in the dataset 
	 * - size_t ndims 	: lenght of each vector in the dataset
	 * - size_t k 		: number of neighbors to retrieve
	 *
	 */
    struct timespec start, finish;
    double elapsed;


	data_dims = (unsigned int)ndims;

	if(verbose)
	{
		printf("Building the KDtree v2:\n");
		clock_gettime(CLOCK_MONOTONIC, &start);
	}

	#ifdef SWMEM
		float_t* dummy_data = malloc(sizeof(float_t)*ndims*n);
		memcpy(dummy_data,data,sizeof(float_t)*ndims*n);
		data = dummy_data;
	#endif

    kdnode_v2* kdnode_array = (kdnode_v2*)malloc(n*sizeof(kdnode_v2));

    initialize_kdnodes_v2(kdnode_array,data,n);

    kdnode_v2* root = build_tree_kdtree_v2(kdnode_array, n, ndims);

	if(verbose)
	{
		clock_gettime(CLOCK_MONOTONIC, &finish);
		elapsed = (finish.tv_sec - start.tv_sec);
		elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
		printf("\tTotal time: %.3lfs\n\n", elapsed);
	}


    datapoint_info* points = (datapoint_info*)malloc(n*sizeof(datapoint_info));

    /**************
     * knn search *
     **************/
	knn_search_kdtree_v2(points,kdnode_array,root,n,k,verbose);

	#ifdef SWMEM
		free(data);
	#endif
	
	for(idx_t i = 0; i < n; ++i) if(kdnode_array[i].node_list.data) free(kdnode_array[i].node_list.data);
	free(kdnode_array);

	return points;
}

void free_datapoint_array(datapoint_info* d, size_t n)
{
    for (idx_t i = 0; i < n; ++i)
    {        
        free_heap(&d[i].ngbh);
    }
    free(d);
}

int float_and_uint_size()
{
	int v = 0;
	int vf = sizeof(float_t) == 8 ? 1 : 0; 
	int vi = sizeof(idx_t) == 8 ? 1 : 0; 
	v = vf + vi*2;
	return v;
}


void knn_search_vptree(datapoint_info* dp_info, vpnode* vpnode_array,vpnode* root,idx_t k,size_t n, float_t (*metric)(void*, void*), int verbose)
{	
	/*
	 * Helper function for performing knn search
	 *
	 * - datapoint_info* dp_info 	: array of datapoints to store the neighborhood on 
	 * - vpnode_v2* vpnode_array 	: array of nodes in the tree
	 * - vpnode_v2* root 		: root of the tree
	 * - idx_t k 					: number of neighbors to retrieve
	 * - size_t n 					: number of elements in the dataset 
	 * - float_t (*metric)(void*, void*)) : distance function
	 */
    struct timespec start_tot, finish_tot;
    double elapsed_tot;

	if(verbose)
	{
		printf("knn search:\n");
		clock_gettime(CLOCK_MONOTONIC, &start_tot);
	}
	

	#ifdef PROGRESS_BAR
		idx_t progress_count = 0;
		idx_t step = n/100;
		printf("Progress 0/%lu -> 0%%\r",(uint64_t)n);
	#endif
    fflush(stdout);
    
    #pragma omp parallel
    {

		#ifdef ITERATIVE_VPTREE
			stack_vpnodes stack;
			stackInit(&stack);
		#endif

	    #pragma omp for schedule(dynamic)
		for(size_t i = 0; i < n; ++i) 
		{
			#ifdef ITERATIVE_VPTREE
				dp_info[i].ngbh = knn_vptree(vpnode_array[i].data, root,  k, &stack, metric);
			#else
				dp_info[i].ngbh = knn_vptree(vpnode_array[i].data, root,  k, metric);
			#endif 

			dp_info[i].cluster_idx = -1;
			dp_info[i].is_center = 0;
			dp_info[i].array_idx = i;
			
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

	if(verbose)
	{
		clock_gettime(CLOCK_MONOTONIC, &finish_tot);
		elapsed_tot = (finish_tot.tv_sec - start_tot.tv_sec);
		elapsed_tot += (finish_tot.tv_nsec - start_tot.tv_nsec) / 1000000000.0;
		printf("\tTotal time: %.3lfs\n\n", elapsed_tot);
	}
    return;
}


datapoint_info* ngbh_search_vptree(void* data, size_t n, size_t byteSize, size_t dims, size_t k, float_t (*metric)(void *, void *), int verbose)
{
	/*
	 * Neighborhood search using vantage-point tree
	 *
	 * - void* data 	 : the dataset, in general a collection of elements					 
	 * - size_t n 		 : number of points in the dataset
	 * - size_t byteSize : number of bytes for each element in the vector representing a datapooint 
	 * - size_t dims 	 : lenght of each datapoint in terms of elements 
	 * - size_t k 		 : number of neighbors to search for 
	 * - float_t (*metric)(void *, void *)) : distance function between two elements in the dataset MUST satisy triangle inequality 
	 */
    struct timespec start, finish;
    double elapsed;
	METRICS_DATADIMS = (uint32_t)dims;
	
	if(verbose)
	{
		clock_gettime(CLOCK_MONOTONIC, &start);
		printf("Building the vp tree v1\n");
	}

	vpnode* vpnode_array = (vpnode*)malloc(n*sizeof(vpnode));
	vpnode** vpPtrArray = (vpnode**)malloc(n*sizeof(vpnode*));
	initialize_vpnode_array(vpnode_array, data, n, byteSize*dims);
	initialize_vpnode_ptrs(vpPtrArray, vpnode_array, n);

	vpnode* root = build_vptree(vpPtrArray, 0, n-1, NULL, metric);
	if(verbose)
	{
		clock_gettime(CLOCK_MONOTONIC, &finish);
		elapsed = (finish.tv_sec - start.tv_sec);
		elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
		printf("\tTotal time: %.3lfs\n\n", elapsed);
	}

    datapoint_info* points = (datapoint_info*)malloc(n*sizeof(datapoint_info));

    /**************
     * knn search *
     **************/

	
	knn_search_vptree(points, vpnode_array, root, k, n, metric, verbose);
	

    free(vpPtrArray);
    free(vpnode_array);
	return points;
}

void knn_search_vptree_v2(datapoint_info* dp_info, vpnode_v2* vpnode_array,vpnode_v2* root,idx_t k,size_t n, float_t (*metric)(void*, void*), int verbose)
{	
	/*
	 * Helper function for performing knn search
	 *
	 * - datapoint_info* dp_info 	: array of datapoints to store the neighborhood on 
	 * - vpnode_v2* vpnode_array 	: array of nodes in the tree
	 * - vpnode_v2* root 		: root of the tree
	 * - idx_t k 					: number of neighbors to retrieve
	 * - size_t n 					: number of elements in the dataset 
	 * - float_t (*metric)(void*, void*)) : distance function
	 */
    struct timespec start_tot, finish_tot;
    double elapsed_tot;

	if(verbose)
	{
		printf("knn search:\n");
		clock_gettime(CLOCK_MONOTONIC, &start_tot);
	}

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
			idx_t idx = vpnode_array[i].array_idx;
			dp_info[idx].ngbh = knn_vptree_v2(vpnode_array[i].data, root, k, metric);
			dp_info[idx].cluster_idx = -1;
			dp_info[idx].is_center = 0;
			dp_info[idx].array_idx = idx;
			
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

	if(verbose)
	{
		clock_gettime(CLOCK_MONOTONIC, &finish_tot);
		elapsed_tot = (finish_tot.tv_sec - start_tot.tv_sec);
		elapsed_tot += (finish_tot.tv_nsec - start_tot.tv_nsec) / 1000000000.0;
		printf("\tTotal time: %.3lfs\n\n", elapsed_tot);
	}
    return;
}

datapoint_info* ngbh_search_vptree_v2(void* data, size_t n, size_t byteSize, size_t dims, size_t k, float_t (*metric)(void *, void *), int verbose)
{
	/*
	 * Neighborhood search using vantage-point tree
	 *
	 * - void* data 	 : the dataset, in general a collection of elements					 
	 * - size_t n 		 : number of points in the dataset
	 * - size_t byteSize : number of bytes for each element in the vector representing a datapooint 
	 * - size_t dims 	 : lenght of each datapoint in terms of elements 
	 * - size_t k 		 : number of neighbors to search for 
	 * - float_t (*metric)(void *, void *)) : distance function between two elements in the dataset MUST satisy triangle inequality 
	 */

    struct timespec start, finish;
    double elapsed;
	METRICS_DATADIMS = (uint32_t)dims;
	
	if(verbose)
	{
		clock_gettime(CLOCK_MONOTONIC, &start);
		printf("Building the vp tree v2\n");
	}

	#ifdef SWMEM
		void* dummy_data = malloc(byteSize*dims*n);
		memcpy(dummy_data,data,byteSize*dims*n);
		data = dummy_data;
	#endif

	vpnode_v2* vpnode_array = (vpnode_v2*)malloc(n*sizeof(vpnode_v2));
	initialize_vpnode_v2_array(vpnode_array, data, n, byteSize*dims);


	//vpnode_v2* root = build_vptree_v2(vpPtrArray, 0, n-1, NULL, metric);
	vpnode_v2* root = build_vptree_v2(vpnode_array, 0, n-1, NULL, metric);
	if(verbose)
	{
		clock_gettime(CLOCK_MONOTONIC, &finish);
		elapsed = (finish.tv_sec - start.tv_sec);
		elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
		printf("\tTotal time: %.3lfs\n\n", elapsed);
	}

    datapoint_info* points = (datapoint_info*)malloc(n*sizeof(datapoint_info));

	knn_search_vptree_v2(points, vpnode_array, root, k, n, metric, verbose);
	
	#ifdef SWMEM
		free(data);
	#endif

	#ifdef VOPT
		for(idx_t i = 0; i < n; ++i) if(vpnode_array[i].node_list.indexes) free(vpnode_array[i].node_list.indexes);
	#else
		for(idx_t i = 0; i < n; ++i) if(vpnode_array[i].node_list.data) free(vpnode_array[i].node_list.data);
	#endif

    free(vpnode_array);
	return points;
}

void set_rho_err_k(datapoint_info* points, float_t* rho, float_t* rho_err, idx_t* k, size_t n)
{
	for(size_t i = 0; i < n; ++i)
	{
		points[i].log_rho = rho[i];
		points[i].log_rho_err = rho_err[i];
		points[i].g = points[i].log_rho - points[i].log_rho_err;
		points[i].kstar = k[i];
	}
	return;
}

void compute_avg(datapoint_info* p, float_t *va, float_t* ve, float_t* vals, float_t* verr, size_t k, size_t n)
{
	#pragma omp parallel for
	for(size_t i = 0; i < n; ++i)
	{
		heap H = p[i].ngbh;
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

/*
 * Interface to python
 */

datapoint_info* alloc_datapoints(idx_t n)
{
	datapoint_info* points = (datapoint_info*)malloc(sizeof(datapoint_info)*n);
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

void import_neighbors_and_distances(datapoint_info* points, idx_t* indeces, float_t* distances, idx_t n, idx_t k)
{
	for(idx_t p = 0; p < n; ++p)
	{
		heap H;
		allocate_heap(&H, k);
		init_heap(&H);
		for(idx_t j = 0; j < k; ++j) 
		{
			H.data[j].value = distances[p*k + j]*distances[p*k + j];
			H.data[j].array_idx = indeces[p*k + j];
		}
		H.count = k;
		points[p].ngbh = H;
	}
}

void import_density(datapoint_info* points, idx_t* kstar, float_t* density, float_t* density_err, idx_t n)
{
	for(idx_t p = 0; p < n; ++p)
	{
		points[p].log_rho 	  = density[p];
		points[p].log_rho_err = density_err[p];
		points[p].g 		  = density[p] - density_err[p]; 
		points[p].kstar		  = kstar[p]; 
	}
}

void export_density(datapoint_info* points, idx_t* kstar, float_t* density, float_t* density_err, idx_t n)
{
	for(idx_t p = 0; p < n; ++p)
	{
		density[p] 		= 	points[p].log_rho;
		density_err[p]  =  	points[p].log_rho_err;
		kstar[p]        =  	points[p].kstar;		  
	}
}

void export_neighbors_and_distances(datapoint_info* points, idx_t* dist_indices, float_t* dists, idx_t n, idx_t k)
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

void export_cluster_assignment(datapoint_info* points, int* labels, idx_t n)
{
	for(idx_t i = 0; i < n; ++i) labels[i] = points[i].cluster_idx;
}

void export_borders(clusters* clusters, int* border_idx, float_t* border_den, float_t* border_err)
{
	idx_t nclus = clusters -> centers.count; 
	if(clusters->use_sparse_borders)		
	{
		for(idx_t i = 0; i < nclus; ++i)		
			for(idx_t el = 0; el < clusters -> sparse_borders[i].count; ++el)
			{
				idx_t j = clusters -> sparse_borders[i].data[el].i;
				idx_t p = i*nclus + j;
				border_idx[p] = (int)j;  
				border_den[p] = clusters -> sparse_borders[i].data[el].density; 
				border_err[p] = clusters -> sparse_borders[i].data[el].error;
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

void reset_datapoints(datapoint_info* dp, size_t n)
{
	for(size_t i = 0; i < n; ++i)	
	{
		dp[i].cluster_idx = -1;
	}
}

void set_verbose_output(int s)
{
	verbose = s;
}












