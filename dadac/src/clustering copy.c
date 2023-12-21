//#include "../include/read_fof_snapshot.h"
#include "../include/clustering.h"
#include <time.h>
extern unsigned int data_dims;
size_t Npart;
const border_t border_null = {.density = -1.0, .error = 0, .idx = NOBORDER};

void LinkedList_Insert(LinkedList* L, Node* n)
{
    ++(L -> count);
    n -> next = L -> head;
    L -> head = n;
}

/*****************************
 * Clusters object functions *
 *****************************/


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

    size_t nclus = c -> centers.count;
    
    c -> __borders_data         = (border_t*)malloc(nclus*nclus*sizeof(border_t)); 
    c -> borders                = (border_t**)malloc(nclus*sizeof(border_t*));
    for(size_t i = 0; i < nclus; ++i)
    {
        c -> borders[i]         = c -> __borders_data + i*nclus;
        for(size_t j = 0; j < nclus; ++j)
        {
            c -> borders[i][j] = border_null;
        }
    }
}

void Clusters_Reset(Clusters * c)
{
    free(c -> centers.data);
    free(c -> __borders_data);
    free(c -> borders);
}

void Clusters_free(Clusters * c)
{

    Clusters_Reset(c);
}

/*****************
 * Dyanmic Array *
 *****************/

void DynamicArray_allocate(lu_dynamicArray * a)
{
    a -> data = (size_t*)malloc(ARRAY_INCREMENT*sizeof(size_t));
    a -> count = 0;
    a -> size = ARRAY_INCREMENT;
}

void DynamicArray_pushBack(lu_dynamicArray * a, size_t p)
{
    if(a -> count < a -> size)
    {
        a -> data[a -> count] =  p;
        a -> count += 1;
    }
    else{
        a -> size += ARRAY_INCREMENT;
        a -> data = realloc(a -> data, a -> size * sizeof(size_t));
        a -> data[a -> count] =  p;
        a -> count += 1;
    }
}

void DynamicArray_Reset(lu_dynamicArray * a){
    a -> count = 0;
}

void DynamicArray_Reserve(lu_dynamicArray * a, size_t n)
{
    a -> data = realloc(a -> data, n*sizeof(size_t));
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

void KNN_search(Datapoint_info * particles, FLOAT_TYPE * data, kd_node* root, size_t n, size_t k)
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



FLOAT_TYPE avg(const FLOAT_TYPE * x, const size_t n)
{
    FLOAT_TYPE f = 0;
    for(size_t i = 0; i < n; ++i)
    {
        f += x[i];
    }
    return f/(FLOAT_TYPE)n;
}



FLOAT_TYPE mEst2(FLOAT_TYPE * x, FLOAT_TYPE *y, size_t n)
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
    for(size_t i = 0; i < n; ++i)
    {
        FLOAT_TYPE xx = x[i];
        FLOAT_TYPE yy = y[i];

        dd = xx;
        num += dd*yy;
        den += dd*dd;

    }
  
    return num/den;
}
FLOAT_TYPE mEst(FLOAT_TYPE * x, FLOAT_TYPE *y, size_t n)
{
    FLOAT_TYPE x_avg, y_avg;
    x_avg = avg(x,n);
    y_avg = avg(y,n);
    FLOAT_TYPE num = 0;
    FLOAT_TYPE den = 0;
    FLOAT_TYPE dd;
    for(size_t i = 0; i < n - 1; ++i)
    {
        FLOAT_TYPE xx = x[i];
        FLOAT_TYPE yy = y[i];

        dd = (xx - x_avg);
        num += dd*(yy - y_avg);
        den += dd*dd;

    }
  
    return num/den;
}

FLOAT_TYPE idEstimate(Datapoint_info* particles, size_t n)
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

    size_t Neff = (size_t)(n*fraction);

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

void computeRho(Datapoint_info* particles, const FLOAT_TYPE d, const size_t points){

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

    size_t kMAX = particles[0].ngbh.N;   

    FLOAT_TYPE omega = 0.;  
    if(sizeof(FLOAT_TYPE) == sizeof(float)){ omega = powf(PI_F,d/2)/tgammaf(d/2.0f + 1.0f);}  
    else{omega = pow(M_PI,d/2.)/tgamma(d/2.0 + 1.0);}

    //printf("Omega d %f\n", omega);

    #pragma omp parallel for
    for(size_t i = 0; i < points; ++i)
    {

        size_t j = 4;
        size_t k;
        FLOAT_TYPE dL = 0.0;
        FLOAT_TYPE vvi, vvj, vp;
        while(j < kMAX && dL < DTHR)
        {
            size_t ksel = j - 1;
            vvi = omega * pow(particles[i].ngbh.data[ksel].value,d/2.);
            size_t jj = particles[i].ngbh.data[j].array_idx;
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
void computeCorrection(Datapoint_info* particles, size_t n, FLOAT_TYPE Z)
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
        for(size_t i = 0; i < n; ++i)
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
        for(size_t i = 0; i < n; ++i)
        {
            particles[i].log_rho_c = particles[i].log_rho - min_log_rho + 1;
            particles[i].g = particles[i].log_rho_c - particles[i].log_rho_err;
        }
    }
    //printf("%lf\n",min_log_rho);
}

Clusters Heuristic1(Datapoint_info* particles, FLOAT_TYPE* data, size_t n)
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

    size_t ncenters = 0;
    size_t putativeCenters = n;
    size_t max_k = particles[0].ngbh.N;
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

    for(size_t i = 0; i < n; ++i)
    {   
        /*

        Find the centers of the clusters as the points of higher density in their neighborhoods
        A point is tagged as a putative center if it is the point of higer density of its neighborhood 
        
        */

        particles_ptrs[i] = particles + i;
        size_t maxk = particles[i].kstar + 1;
        FLOAT_TYPE gi = particles[i].g;
        particles[i].is_center = 1;
        particles[i].cluster_idx = -1;
        //printf("%lf\n",p -> g);
        Heap i_ngbh = particles[i].ngbh;
        for(size_t k = 1; k < maxk; ++k)
        {
            size_t ngbh_index = i_ngbh.data[k].array_idx;
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

    size_t * to_remove = (size_t*)malloc(allCenters.count*sizeof(size_t));

    qsort(particles_ptrs, n, sizeof(Datapoint_info*), cmpPP);
    /*****************************************************************************************************
     * /!\ This part is VERY time consuming, complexity depends on the number of center previously found *
     * /!\ It is actually faster when using only a few threads                                           *
     *****************************************************************************************************/
    //#pragma omp parallel for num_threads(4)
    #pragma omp parallel for 
    for(size_t p = 0; p < allCenters.count; ++p)
    {   
        /*

        Check if the center spotted in the previous part belongs to the neighborhood
        of a point of higher density. If this is true remove it from the actual centers

        */
        size_t i = allCenters.data[p];
        int e = 0;
        FLOAT_TYPE gi = particles[i].g;
        size_t i_arrIdx = particles[i].array_idx;
        size_t mr = SIZE_MAX;
        FLOAT_TYPE max_g = -99999.0;
        for(size_t j = 0; j < n; ++j)
        {
            //retrive the particle pointed by the pointer
            Datapoint_info pp = *(particles_ptrs[j]);
            size_t kMAXj = pp.kstar;
            Heap j_ngbh = pp.ngbh;
            FLOAT_TYPE gj = pp.g;
            //e = 0;
            //check if there is point in which the point i is a neighbor with grater g
                //if gj > gi check the neighborhood
            //preliminarity check, if i is more distant than k* neighbor break
            FLOAT_TYPE dk = j_ngbh.data[kMAXj + 1].value;
            FLOAT_TYPE di = euclidean_distance(data + (i*data_dims), data + (pp.array_idx*data_dims));
            int stop = 0;
            if(dk > di)
            {
                mr = pp.array_idx;
                //found a neighborhood with higher g, break
                stop = 1;
                //max_g = gj;
            } 
            if(stop == 1|| i == pp.array_idx)
            {
                break;
            }
            
        }
        to_remove[p] = mr;
    }
    
    //#pragma omp parallel for 
    //for(size_t p = 0; p < allCenters.count; ++p)
    //{   
    //    /*

    //    Check if the center spotted in the previous part belongs to the neighborhood
    //    of a point of higher density. If this is true remove it from the actual centers

    //    */
    //    size_t i = allCenters.data[p];
    //    int e = 0;
    //    FLOAT_TYPE gi = particles[i].g;
    //    size_t i_arrIdx = particles[i].array_idx;
    //    size_t mr = SIZE_MAX;
    //    FLOAT_TYPE max_g = -99999.0;
    //    for(size_t j = 0; j < n; ++j)
    //    {
    //        size_t kMAXj = particles[j].kstar;
    //        Heap j_ngbh = particles[j].ngbh;
    //        FLOAT_TYPE gj = particles[j].g;
    //        //e = 0;
    //        //check if there is point in which the point i is a neighbor with grater g
    //            //if gj > gi check the neighborhood
    //        //preliminarity check, if i is more distant than k* neighbor break
    //        FLOAT_TYPE dk = j_ngbh.data[kMAXj + 1].value;
    //        FLOAT_TYPE di = euclidean_distance(data + (i*data_dims), data + (j*data_dims));
    //        if(dk > di && gj > gi)
    //        {
    //            for(size_t k = 1; k < kMAXj + 1; ++k )
    //            {
    //                if(j_ngbh.data[k].array_idx == i_arrIdx )
    //                {
    //                    if(gj > max_g)
    //                    {
    //                            mr = j;
    //                            max_g = gj;
    //                    }
    //                    break;
    //                }
    //            }
    //            
    //        } 
    //        
    //    }
    //    to_remove[p] = mr;
    //}
    

    for(size_t p = 0; p < allCenters.count; ++p)
    {
        size_t i = allCenters.data[p];
        int e = 0;
        FLOAT_TYPE gi = particles[i].g;
        size_t mr = to_remove[p];
        if(mr != SIZE_MAX)
        {
            if(particles[mr].g > gi) e = 1;
        }
        switch (e)
        {
            case 1:
                {
                    DynamicArray_pushBack(&removedCenters,i);
                    particles[i].is_center = 0;
                    //for(size_t c = 0; c < removedCenters.count - 1; ++c)
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

    size_t nclusters = 0;


    /*****************************************************************************
     * Sort all the particles based on g and then perform the cluster assignment *
     * in asceding order                                                         *
     * UPDATE: particles already sorted                                          *
     *****************************************************************************/
                                                                                

    //qsort(particles_ptrs, n, sizeof(Datapoint_info*), cmpPP);

    for(size_t i = 0; i < n; ++i)
    {   
        Datapoint_info* p = particles_ptrs[i];
        size_t ele = p -> array_idx;
        //fprintf(f,"%lu\n",ele);

        if(!(p -> is_center))
        {
            int cluster = -1;
            size_t k = 0;
            size_t p_idx;
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
                size_t gm_index;
                for(size_t k = 0; k < max_k; ++k)
                {
                    size_t ngbh_index = p -> ngbh.data[k].array_idx;
                    for(size_t m = 0; m < removedCenters.count; ++m)
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
    //FILE* f = fopen("nope7.dat","w");
    //for(int i = 0; i < allCenters.count; ++i)
    //{
    //    fprintf(f,"%lu\n",allCenters.data[i]);
    //}
    //fclose(f);

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


    printf("\tFound %ld clusters\n",actualCenters.count);
    printf("\tTotal time: %.3lfs\n\n", elapsed_tot);

    c_all.n = n;
    return c_all;
}

void Heuristic2(Clusters* cluster, Datapoint_info* particles)
{

    #define borders cluster->borders

    struct timespec start_tot, finish_tot;
    double elapsed_tot;
    size_t n = cluster -> n;

    printf("H2: Finding border points\n");
    clock_gettime(CLOCK_MONOTONIC, &start_tot);


    size_t nclus = cluster->centers.count; 
    size_t max_k = particles[0].ngbh.N;

    for(size_t i = 0; i < n; ++i)
    {
            /*retrive the index of the particle in the cluster c*/
            size_t pp = NOBORDER;
            /*loop over n neighbors*/
            int c = particles[i].cluster_idx;
            if(!particles[i].is_center)
            {
                for(size_t k = 1; k < particles[i].kstar + 1; ++k)
                {
                    /*index of the kth ngbh of n*/
                    size_t j = particles[i].ngbh.data[k].array_idx;
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
                for(size_t k = 1; k < max_k; ++k)
                {
                    size_t pp_ngbh_idx = particles[pp].ngbh.data[k].array_idx;
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
                if(particles[i].g > borders[c][ppc].density)
                {
                    borders[c][ppc].density = particles[i].g;
                    borders[ppc][c].density = particles[i].g;
                    borders[c][ppc].idx = i;
                    borders[ppc][c].idx = i;
                }
            }

            /*step into the linked list*/
            //printf("%p %d %d\n", n -> next, particles[n -> data].cluster_idx, c);
    }


    //for(int c = 0; c < nclus; ++c)
    //{
    //    Node* n = cluster -> clusters[c].head;
    //    while(n)
    //    {
    //        /*retrive the index of the particle in the cluster c*/
    //        size_t i = n -> data;
    //        size_t pp = NOBORDER;
    //        /*loop over n neighbors*/
    //        if(!particles[i].is_center)
    //        {
    //            for(size_t k = 1; k < particles[i].kstar + 1; ++k)
    //            {
    //                /*index of the kth ngbh of n*/
    //                size_t j = particles[i].ngbh.data[k].array_idx;
    //                pp = NOBORDER;
    //                /*Loop over kn neigbhours to find if n is the nearest*/
    //                /*if cluster of the particle in nbhg is c then check is neighborhood*/                                                
    //                if(particles[j].cluster_idx != c)
    //                {
    //                    pp = j;
    //                    break;
    //                }

    //            }
    //        }

    //        if(pp != NOBORDER)
    //        {
    //            for(size_t k = 1; k < max_k; ++k)
    //            {
    //                size_t pp_ngbh_idx = particles[pp].ngbh.data[k].array_idx;
    //                if(pp_ngbh_idx == i)
    //                {
    //                    break;
    //                }
    //                if(particles[pp_ngbh_idx].cluster_idx == c)
    //                {
    //                    pp = NOBORDER;
    //                    break;
    //                }
    //            }
    //        }
    //                        /*if it is the maximum one add it to the cluster*/
    //        if(pp != NOBORDER)
    //        {
    //            int ppc = particles[pp].cluster_idx;
    //            if(particles[i].g > border_density[c][ppc])
    //            {
    //                border_density[c][ppc] = particles[i].g;
    //                border_density[ppc][c] = particles[i].g;
    //                border_idx[c][ppc] = i;
    //                border_idx[ppc][c] = i;
    //            }
    //        }

    //        /*step into the linked list*/
    //        //printf("%p %d %d\n", n -> next, particles[n -> data].cluster_idx, c);
    //        n = n -> next;
    //    }


    //}
    for(size_t i = 0; i < nclus - 1; ++i)
    {
        for(size_t j = i + 1; j < nclus; ++j)
        {
            size_t p = borders[i][j].idx;
            if(p != NOBORDER)
            {   

                borders[i][j].density = particles[p].log_rho_c;
                borders[j][i].density = particles[p].log_rho_c;

                borders[i][j].error = particles[p].log_rho_err;
                borders[j][i].error = particles[p].log_rho_err;
            }
        }
    }

    for(size_t i = 0; i < nclus; ++i)
    {
        borders[i][i].density = -1.0;
        borders[i][i].error = 0.0;
    }

    clock_gettime(CLOCK_MONOTONIC, &finish_tot);
    elapsed_tot = (finish_tot.tv_sec - start_tot.tv_sec);
    elapsed_tot += (finish_tot.tv_nsec - start_tot.tv_nsec) / 1000000000.0;
    printf("\tTotal time: %.3lfs\n\n", elapsed_tot);

    return;
    #undef borders
   }


//void mergeClusters(Clusters * cc, size_t i, size_t j)
//{
//    LinkedList * li = (cc->clusters) + i;
//    LinkedList * lj = (cc->clusters) + j;
//    Node * n = li -> head;
//    /*find the tail node of li*/
//    while(n -> next)
//    {
//        n = n -> next;
//    }
//    
//    n -> next = lj -> head;
//    /*Reset list lj*/
//    lj -> head = NULL;
//    lj -> count = 0;
//}

void Merge_A_into_B(size_t* who_amI, size_t cluster_A, size_t cluster_B, size_t n)
{
    size_t tmp;
    for(size_t i = 0; i < n; ++i)
    {   
        //substitute occurencies of b with a 
        tmp = who_amI[i] == cluster_A ? cluster_B : who_amI[i];
        who_amI[i] = tmp;
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
      
  FLOAT_TYPE c1 = (dens1 - dens_border) /
	(dens1_err + dens_border_err); 
  FLOAT_TYPE c2 = (dens2 - dens_border) /
    (dens2_err + dens_border_err);
  
  return ( c1 < c2 );     // if 1, this signal to swap 1 and 2
}

void fix_borders_A_into_B(size_t A, size_t B, border_t** borders, size_t n)
{
   for(size_t i = 0; i < n; ++i) 
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

void Heuristic3(Clusters* cluster, Datapoint_info* particles, FLOAT_TYPE Z, int halo)
{
  #define borders cluster->borders

  struct timespec start_tot, finish_tot;
  double elapsed_tot;

  printf("H3: Merging clusters\n");
  clock_gettime(CLOCK_MONOTONIC, &start_tot);
  

  size_t nclus              = cluster -> centers.count;  
  size_t *  surviving_clusters = (size_t*)malloc(nclus*sizeof(size_t));
  for(size_t i = 0; i < nclus; ++i)
    { 
        surviving_clusters[i] = i; 
    }

  size_t   merge_count        = 0;
  size_t   merging_table_size = 1000;
  merge_t *merging_table      = (merge_t*)malloc(sizeof(merge_t)*merging_table_size);
  
  /*Find clusters to be merged*/
  for(size_t i = 0; i < nclus - 1; ++i)   
    for(size_t j = i + 1; j < nclus; ++j)   
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

		  int swap = merging_roles( dens1, dens1_err, dens2, dens2_err, dens_border, dens_border_err);
		  size_t src;
		  size_t trg;
		  switch ( swap )
		    {
		    case 0: { src = j; trg = i;} break;
		    case 1: { src = i; trg = j;} break;
		    }

		  merging_table[merge_count].source = src;
		  merging_table[merge_count].target = trg;
		  merging_table[merge_count].density = borders[trg][src].density;
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
  
    for( size_t m = 0; m < merge_count; m++ )
    {
 //    #define BORDER_NULL( S, T ) { border_density[(S)][(T)] = -1.0; \
 //     border_density[(T)][(S)] = -1.0;				\
 //     border_err[(S)][(T)]     = 0;				\
 //     border_err[(T)][(S)]     = 0;				\
 //     border_idx[(S)][(T)]     = NOBORDER;			\
 //     border_idx[(T)][(S)]     = NOBORDER; }

      
        #define src surviving_clusters[merging_table[m].source]
        #define trg surviving_clusters[merging_table[m].target]
        //printf("Found: %lu, %lu which now is %lu, %lu\n",merging_table[m].source, merging_table[m].target, src,trg);

        int re_check = ( (src != merging_table[m].source) || (trg != merging_table[m].target) );

        switch ( re_check )
        {
            case 0: 
            {
                printf("Merging %lu into %lu\n",src,trg);
                borders[src][trg] = border_null;
                borders[trg][src] = border_null;
                fix_borders_A_into_B(src,trg,borders,nclus);
                Merge_A_into_B (surviving_clusters, src, trg, nclus );	  
            } break;
        
            case 1: 
            {
                //pick who am I
                size_t new_src             = surviving_clusters[src];
                size_t new_trg             = surviving_clusters[trg];

                FLOAT_TYPE dens1           = particles[cluster->centers.data[src]].log_rho_c;
                FLOAT_TYPE dens1_err       = particles[cluster->centers.data[src]].log_rho_err;
                
                FLOAT_TYPE dens2           = particles[cluster->centers.data[trg]].log_rho_c;
                FLOAT_TYPE dens2_err       = particles[cluster->centers.data[trg]].log_rho_err;

                FLOAT_TYPE dens_border     = borders[src][trg].density;
                FLOAT_TYPE dens_border_err = borders[src][trg].error;

                int i_have_to_merge = is_a_merging(dens1,dens1_err,dens2,dens2_err,dens_border,dens_border_err,Z);            
                switch (i_have_to_merge && src != trg)
                {
                case 1:
                    {
                        size_t new_src = src;
                        size_t new_trg = trg;
                        int side = merging_roles(dens1,dens1_err,dens2,dens2_err,dens_border,dens_border_err);
                        if(!side)
                        {
                            size_t tmp;
                            tmp = new_src;
                            new_src = new_trg;
                            new_trg = tmp;
                        }

                        borders[new_src][new_trg] = border_null;
                        borders[new_trg][new_src] = border_null;
                        fix_borders_A_into_B(new_src,new_trg,borders,nclus);
                        Merge_A_into_B ( surviving_clusters, new_src, new_trg, nclus );	  
                        printf("Merging %lu into %lu\n",new_src,new_trg);
                    }
                    break;
                
                default:
                    break;
                }
            }
        }
        
        #undef src
        #undef trg
    }


  
        /*Merge clusters*/
    //    if(check)
    //    {
    //        //printf("Hey!\n");
    //        size_t barrier_index = 0;
    //        FLOAT_TYPE max_border = 0;
    //        for(size_t k = 0; k < ipos.count; ++k)
    //        {
    //            if(border_density[ipos.data[k]][jpos.data[k]] > max_border)
    //            {
    //                max_border = border_density[ipos.data[k]][jpos.data[k]];
    //                barrier_index = k;
    //            }
    //        }

    //        size_t imod = ipos.data[barrier_index];
    //        size_t jmod = jpos.data[barrier_index];
    //        //printf("%lu %lu %lu\n",imod,jmod,barrier_index);
    //        size_t ci = cluster->centers.data[imod];
    //        size_t cj = cluster->centers.data[jmod];
    //        

    //        FLOAT_TYPE c1 = (particles[ci].log_rho_c - border_density[imod][jmod]) / (particles[ci].log_rho_err + border_err[imod][jmod]); 
    //        FLOAT_TYPE c2 = (particles[cj].log_rho_c - border_density[imod][jmod]) / (particles[cj].log_rho_err + border_err[imod][jmod]); 

    //        if(c1 < c2) 
    //        {
    //            size_t tmp = jmod;
    //            jmod = imod;
    //            imod = tmp;
    //        }
    //        
    //        surviving_clusters[jmod] = 0;

    //        border_density[imod][jmod] = -1.0;
    //        border_density[jmod][imod] = -1.0;

    //        border_err[imod][jmod] = 0;
    //        border_err[jmod][imod] = 0;

    //        border_idx[imod][jmod] = NOBORDER;
    //        border_idx[jmod][imod] = NOBORDER;
    //        
    //        Merge_A_into_B(who_amI, jmod, imod, nclus);
    //        //mergeClusters(cluster, imod, jmod);

    //        /*Recompute borders*/
    //        for(size_t i = 0; i < nclus; ++i)
    //        {
    //            if(i != imod && i != jmod)
    //            {
    //                if(border_density[imod][i] < border_density[jmod][i])
    //                {
    //                    border_idx[imod][i] = border_idx[jmod][i];
    //                    border_idx[i][imod] = border_idx[i][jmod];

    //                    border_density[imod][i] = border_density[jmod][i];
    //                    border_density[i][imod] = border_density[i][jmod];

    //                    border_err[imod][i] = border_err[jmod][i];
    //                    border_err[i][imod] = border_err[i][jmod];
    //                }

    //                border_density[jmod][i] = -1.0;
    //                border_density[i][jmod] = -1.0;
    //                border_err[jmod][i] = 0;
    //                border_err[i][jmod] = 0;
    //                border_idx[jmod][i] = NOBORDER;
    //                border_idx[i][jmod] = NOBORDER;

    //            }
    //        }
    //    }
    //}

    /*Finalize clustering*/
    /*Acutally copying */
    lu_dynamicArray tmp_centers;
    lu_dynamicArray tmp_cluster_idx;


    DynamicArray_Init(&tmp_centers);
    DynamicArray_Init(&tmp_cluster_idx);

    DynamicArray_Reserve(&tmp_centers, nclus);
    DynamicArray_Reserve(&tmp_cluster_idx, nclus);

    size_t final_cluster_count = 0;

    size_t* old_to_new = (size_t*)malloc(nclus*sizeof(size_t));
    size_t incremental_k = 0;
    for(size_t i = 0; i < nclus; ++i)
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
    for(size_t i = 0; i < nclus; ++i)
    {
        if(surviving_clusters[i] != i){
            size_t cidx_to_copy_from = surviving_clusters[i];
            old_to_new[i] = old_to_new[cidx_to_copy_from];
        }
    }

    /*allocate auxiliary pointers to store results of the finalization of the procedure*/

    border_t** tmp_borders      = (border_t**)malloc(final_cluster_count*sizeof(border_t*));
    border_t*  tmp_borders_data = (border_t*)malloc(final_cluster_count*final_cluster_count*sizeof(border_t));

    /*initialize all pointers*/
    for(size_t i = 0; i < final_cluster_count; ++i)
    {
        tmp_borders[i] = tmp_borders_data + i*final_cluster_count;
    }

    /*Fix cluster assignment*/
    #pragma omp parallel for
    for(size_t i = 0; i < cluster -> n; ++i)
    {
        particles[i].is_center = 0;
        int old_cidx = particles[i].cluster_idx;
        particles[i].cluster_idx = old_to_new[old_cidx];
    }

    
    #pragma omp parallel for
    for(size_t c = 0; c < final_cluster_count; ++c)
    {
        size_t c_idx = tmp_cluster_idx.data[c];
        for(size_t d = c; d < final_cluster_count; ++d)
        {
            size_t c_jdx = tmp_cluster_idx.data[d];
            tmp_borders[c][d].density = borders[c_idx][c_jdx].density;
            tmp_borders[d][c].density = borders[c_idx][c_jdx].density;

            tmp_borders[c][d].idx = borders[c_idx][c_jdx].idx;
            tmp_borders[d][c].idx = borders[c_idx][c_jdx].idx;


            tmp_borders[c][d].error = borders[c_idx][c_jdx].error;
            tmp_borders[d][c].error = borders[c_idx][c_jdx].error;
        } 
    }

    Clusters_Reset(cluster);
    /*pay attention to the define borders*/
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
        FLOAT_TYPE* max_border_den_array = (FLOAT_TYPE*)malloc(final_cluster_count*sizeof(FLOAT_TYPE));
        #pragma omp parallel
        {
            #pragma omp for
            for(size_t c = 0; c < final_cluster_count; ++c)
            {
                FLOAT_TYPE max_border_den = -2.;
                for(size_t d = 0; d < final_cluster_count; ++d)
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
            for(size_t i = 0; i < cluster -> n; ++i)
            {
                int cidx = particles[i].cluster_idx;
                int halo_flag = particles[i].log_rho_c < max_border_den_array[cidx]; 
                particles[i].cluster_idx = halo_flag ? -1 : cidx;
            }
        }
        break;
    
    default:
        break;
    }    

    /*free memory and put the correct arrays into place*/
    free(tmp_cluster_idx.data);
    //free(ipos.data);
    //free(jpos.data);
    free(surviving_clusters);
    free(old_to_new);

    clock_gettime(CLOCK_MONOTONIC, &finish_tot);
    elapsed_tot = (finish_tot.tv_sec - start_tot.tv_sec);
    elapsed_tot += (finish_tot.tv_nsec - start_tot.tv_nsec) / 1000000000.0;
    printf("\tSurviving clusters %lu\n", final_cluster_count);
    printf("\tTotal time: %.3lfs\n\n", elapsed_tot);

  #undef  borders  
}


