#include "../include/read_fof_snapshot.h"
#include "../include/clustering.h"
#include <time.h>
extern unsigned int data_dims;

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

    c -> __border_idx_data      = (size_t*)malloc(nclus*nclus*sizeof(size_t));
    c -> border_idx             = (size_t**)malloc(nclus*sizeof(size_t*));
    c -> __border_density_data  = (FLOAT_TYPE*)malloc(nclus*nclus*sizeof(FLOAT_TYPE));
    c -> border_density         = (FLOAT_TYPE**)malloc(nclus*sizeof(FLOAT_TYPE*));
    c -> __border_err_data      = (FLOAT_TYPE*)malloc(nclus*nclus*sizeof(FLOAT_TYPE));
    c -> border_err             = (FLOAT_TYPE**)malloc(nclus*sizeof(FLOAT_TYPE*));
    for(size_t i = 0; i < nclus; ++i)
    {
        c -> border_idx[i]      = c -> __border_idx_data + i*nclus;
        c -> border_density[i]  = c -> __border_density_data + i*nclus;
        c -> border_err[i]      = c -> __border_err_data + i*nclus;
        for(size_t j = 0; j < nclus; ++j)
        {
            c -> border_err[i][j]       = 0.0;
            c -> border_idx[i][j]       = NOBORDER;
            c -> border_density[i][j]   = 0.0;
        }
    }
}

void Clusters_Reset(Clusters * c)
{
    free(c -> centers.data);
    free(c -> border_err);
    free(c -> border_density);
    free(c -> border_idx);
    free(c -> __border_err_data);
    free(c -> __border_density_data);
    free(c -> __border_idx_data);
    free(c -> clusters);
}

void Clusters_free(Clusters * c)
{

    Clusters_Reset(c);
    free(c -> _LLnodes);
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

/*Quake inverse square root, just an experiment
**nothing more, VERY fast, but not precise and works
 *only with float*/
//float Q_rsqrt( float number )
//{
//        long i;
//        float x2, y;
//        const float threehalfs = 1.5F;
// 
//        x2 = number * 0.5F;
//        y  = number;
//        i  = * ( long * ) &y; 
//        i  = 0x5f3759df - ( i >> 1 );
//        y  = * ( float * ) &i;
//        y  = y * ( threehalfs - ( x2 * y * y ) );   
// 
//        return y;
//}

/*******************
 * Clustering part *
 *******************/

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

    return;


}

int cmpPP(const void* p1, const void *p2)
{
/***********************************************
 * Utility function to perform quicksort then, *
 * when clustering assignment is performed     *
 ***********************************************/
    Datapoint_info* pp1 = *(Datapoint_info**)p1;
    Datapoint_info* pp2 = *(Datapoint_info**)p2;
    return 2*(pp1 -> g < pp2 -> g) - 1;
}
void calculateCorrection(Datapoint_info* particles, size_t n, FLOAT_TYPE Z)
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


    clock_gettime(CLOCK_MONOTONIC, &start);
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

    clock_gettime(CLOCK_MONOTONIC, &finish);
    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;

    printf("%.3lfs p1\n",elapsed);
    
    clock_gettime(CLOCK_MONOTONIC, &start);

    size_t * to_remove = (size_t*)malloc(allCenters.count*sizeof(size_t));

    #pragma omp parallel for
    for(size_t p = 0; p < allCenters.count; ++p)
    {   
        /*

        Check if the center spotted in the previous part belongs to the neighborhood
        of a point of higher density. If this is true remove it from the actual centers

        */
        size_t i = allCenters.data[p];
        int e = 0;
        //FLOAT_TYPE gi = particles[i].g;
        size_t i_arrIdx = particles[i].array_idx;
        size_t mr = SIZE_MAX;
        FLOAT_TYPE max_g = -99999.0;
        for(size_t j = 0; j < n; ++j)
        {
            //e = 0;
            Heap j_ngbh = particles[j].ngbh;
            size_t kMAXj = particles[j].kstar;
            FLOAT_TYPE gj = particles[j].g;
            //check if there is point in which the point i is a neighbor with grater g
                //if gj > gi check the neighborhood
            //preliminarity check, if i is more distant than k* neighbor break
            FLOAT_TYPE dk = j_ngbh.data[kMAXj + 1].value;
            FLOAT_TYPE di = euclidean_distance(data + (i*data_dims), data + (j*data_dims));
            if(dk < di){
                continue;
            }
            else
            {
                for(size_t k = 1; k < kMAXj + 1; ++k )
                {
                    if(j_ngbh.data[k].array_idx == i_arrIdx )
                    {
                        if(gj > max_g)
                        {
                                mr = j;
                                max_g = gj;
                        }
                        break;
                    }
                }
            } 
            
        }
        to_remove[p] = mr;
        //if(mr != SIZE_MAX)
        //{
        //    if(particles[mr].g > gi) e = 1;
        //}
        //if(e)
        //{
        //            DynamicArray_pushBack(&removedCenters,i);
        //            particles[i].is_center = 0;
        //            for(size_t c = 0; c < removedCenters.count - 1; ++c)
        //            {
        //                if(mr == removedCenters.data[c])
        //                {
        //                    mr = max_rho.data[c];
        //                }
        //            }
        //            DynamicArray_pushBack(&max_rho,mr);
        //}
        //else
        //{
        //            DynamicArray_pushBack(&actualCenters,i);
        //            particles[i].cluster_idx = actualCenters.count - 1;
        //}
        //switch (e)
        //{
        //    case 1:
        //        #pragma omp critical
        //        {
        //            DynamicArray_pushBack(&removedCenters,i);
        //            particles[i].is_center = 0;
        //            //for(size_t c = 0; c < removedCenters.count - 1; ++c)
        //            //{
        //            //    if(mr == removedCenters.data[c])
        //            //    {
        //            //        mr = max_rho.data[c];
        //            //    }
        //            //}
        //            DynamicArray_pushBack(&max_rho,mr);
        //            
        //        }
        //        break;
        //    case 0:
        //        #pragma omp critical
        //        {
        //            DynamicArray_pushBack(&actualCenters,i);
        //            particles[i].cluster_idx = actualCenters.count - 1;
        //        }
        //        break;
        //    default:
        //        break;
        //}
    }
    
    clock_gettime(CLOCK_MONOTONIC, &finish);
    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;

    printf("%.3lfs p2\n",elapsed);

    clock_gettime(CLOCK_MONOTONIC, &start);
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
    clock_gettime(CLOCK_MONOTONIC, &finish);
    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;

    printf("%.3lfs p2\n",elapsed);

    printf("Found %lu centers\n", actualCenters.count);
    size_t nclusters = 0;
    qsort(particles_ptrs, n, sizeof(Datapoint_info*), cmpPP);

    clock_gettime(CLOCK_MONOTONIC, &start);
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

    clock_gettime(CLOCK_MONOTONIC, &finish);
    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;

    printf("%.3lfs p3\n",elapsed);
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


    /**
     * Create the clusters object in order to have a more usefull division of the particles
    */
    Clusters c_all;
    c_all.centers = actualCenters;
    c_all._LLnodes = (Node*)malloc(n*sizeof(Node));
    c_all.clusters = (LinkedList*)malloc(actualCenters.count*sizeof(LinkedList));

    for(size_t i = 0; i < actualCenters.count; ++i)
    {
        c_all.clusters[i].count = 0;
        c_all.clusters[i].head = NULL;
    }

    /**
     * Initialize nodes and put them in the correct cluster
    */
    for(size_t i = 0; i < n; ++i)
    {
        c_all._LLnodes[i].data = i;
        c_all._LLnodes[i].next = NULL;
        int cluster_idx = particles[i].cluster_idx;
        /* Insert particle i in his cluster*/
        LinkedList_Insert(c_all.clusters + cluster_idx, c_all._LLnodes + i);

    }
    //printf("created lisist\n");
    return c_all;
}

void Heuristic2(Clusters* cluster, Datapoint_info* particles)
{
    size_t**     border_idx     = cluster->border_idx;
    FLOAT_TYPE** border_density = cluster->border_density;
    FLOAT_TYPE** border_err     = cluster->border_err;

    size_t nclus = cluster->centers.count; 
    size_t max_k = particles[0].ngbh.N;
    for(int c = 0; c < nclus; ++c)
    {
        Node* n = cluster -> clusters[c].head;
        while(n)
        {
            /*retrive the index of the particle in the cluster c*/
            size_t i = n -> data;
            size_t pp = NOBORDER;
            /*loop over n neighbours*/
            if(!particles[i].is_center)
            {
                for(size_t k = 1; k < particles[i].kstar + 1; ++k)
                {
                    /*index of the kth ngbh of n*/
                    size_t j = particles[i].ngbh.data[k].array_idx;
                    pp = NOBORDER;
                    /*Loop over kn neigbhours to find if n is the nearest*/
                    /*if cluster of the particle in nbhg is c then check is neighbourhood*/                                                
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
                if(particles[i].g > border_density[c][ppc])
                {
                    border_density[c][ppc] = particles[i].g;
                    border_density[ppc][c] = particles[i].g;
                    border_idx[c][ppc] = i;
                    border_idx[ppc][c] = i;
                }
            }

            /*step into the linked list*/
            //printf("%p %d %d\n", n -> next, particles[n -> data].cluster_idx, c);
            n = n -> next;
        }


    }
    for(size_t i = 0; i < nclus - 1; ++i)
    {
        for(size_t j = i + 1; j < nclus; ++j)
        {
            size_t p = border_idx[i][j];
            if(p != NOBORDER)
            {   

                border_density[i][j] = particles[p].log_rho_c;
                border_density[j][i] = particles[p].log_rho_c;

                border_err[i][j] = particles[p].log_rho_err;
                border_err[j][i] = particles[p].log_rho_err;
            }
        }
    }

    for(size_t i = 0; i < nclus; ++i)
    {
        border_density[i][i] = -1.0;
        border_err[i][i] = 0.0;
    }
    return;
   }


void mergeClusters(Clusters * cc, size_t i, size_t j)
{
    LinkedList * li = (cc->clusters) + i;
    LinkedList * lj = (cc->clusters) + j;
    Node * n = li -> head;
    /*find the tail node of li*/
    while(n -> next)
    {
        n = n -> next;
    }
    
    n -> next = lj -> head;
    /*Reset list lj*/
    lj -> head = NULL;
    lj -> count = 0;
}

void Heuristic3(Clusters* cluster, Datapoint_info* particles, FLOAT_TYPE Z, int halo)
{
    size_t**     border_idx     = cluster->border_idx;
    FLOAT_TYPE** border_density = cluster->border_density;
    FLOAT_TYPE** border_err     = cluster->border_err;

    int check = 1;
    size_t nclus = cluster -> centers.count;  
    int * surviving_clusters = (int*)malloc(nclus*sizeof(int));
    for(size_t i = 0; i < nclus; ++i)
    {
        surviving_clusters[i] = 1;
    }

    lu_dynamicArray jpos, ipos;
    DynamicArray_Init(&ipos);
    DynamicArray_Init(&jpos);
    DynamicArray_Reserve(&ipos,nclus);
    DynamicArray_Reserve(&jpos,nclus);

    printf("%ld\n",nclus);

    while(check)
    {
        check = 0;
        DynamicArray_Reset(&ipos);
        DynamicArray_Reset(&jpos);
        /*Find clusters to be merged*/
        for(size_t i = 0; i < nclus - 1; ++i)   
        {
            for(size_t j = i + 1; j < nclus; ++j)   
            {
                FLOAT_TYPE a1 = particles[cluster->centers.data[i]].log_rho_c - border_density[i][j];
                FLOAT_TYPE a2 = particles[cluster->centers.data[j]].log_rho_c - border_density[i][j];

                FLOAT_TYPE e1 = Z*(particles[cluster->centers.data[i]].log_rho_err + border_err[i][j]);
                FLOAT_TYPE e2 = Z*(particles[cluster->centers.data[j]].log_rho_err + border_err[i][j]);

                //printf("%lf %lf %lf %lf\n", a1,a2,e1,e2);
                //printf("%lf %lf\n",particles[cluster->centers.data[i]].log_rho_c, border_density[i][j]);

                
                if( a1 < e1 || a2 < e2)
                {
                    //printf("%lu %lu\n", cluster -> centers.data[i], cluster -> centers.data[j]);
                    DynamicArray_pushBack(&ipos, i);
                    DynamicArray_pushBack(&jpos, j);
                    check = 1;
                }
            }
        }

        /*Merge clusters*/
        if(check)
        {
            //printf("Hey!\n");
            size_t barrier_index = 0;
            FLOAT_TYPE max_border = 0;
            for(size_t k = 0; k < ipos.count; ++k)
            {
                if(border_density[ipos.data[k]][jpos.data[k]] > max_border)
                {
                    max_border = border_density[ipos.data[k]][jpos.data[k]];
                    barrier_index = k;
                }
            }

            size_t imod = ipos.data[barrier_index];
            size_t jmod = jpos.data[barrier_index];
            //printf("%lu %lu %lu\n",imod,jmod,barrier_index);
            size_t ci = cluster->centers.data[imod];
            size_t cj = cluster->centers.data[jmod];
            

            FLOAT_TYPE c1 = (particles[ci].log_rho_c - border_density[imod][jmod]) / (particles[ci].log_rho_err + border_err[imod][jmod]); 
            FLOAT_TYPE c2 = (particles[cj].log_rho_c - border_density[imod][jmod]) / (particles[cj].log_rho_err + border_err[imod][jmod]); 

            if(c1 < c2) 
            {
                size_t tmp = jmod;
                jmod = imod;
                imod = tmp;
            }
            
            surviving_clusters[jmod] = 0;

            border_density[imod][jmod] = -1.0;
            border_density[jmod][imod] = -1.0;

            border_err[imod][jmod] = 0;
            border_err[jmod][imod] = 0;

            mergeClusters(cluster, imod, jmod);

            /*Recompute borders*/
            for(size_t i = 0; i < nclus; ++i)
            {
                if(i != imod && i != jmod)
                {
                    if(border_density[imod][i] < border_density[jmod][i])
                    {
                        border_idx[imod][i] = border_idx[jmod][i];
                        border_idx[i][imod] = border_idx[i][jmod];

                        border_density[imod][i] = border_density[jmod][i];
                        border_density[i][imod] = border_density[i][jmod];

                        border_err[imod][i] = border_err[jmod][i];
                        border_err[i][imod] = border_err[i][jmod];
                    }

                    border_density[jmod][i] = -1.0;
                    border_density[i][jmod] = -1.0;
                    border_err[jmod][i] = 0;
                    border_err[i][jmod] = 0;

                }
            }
        }
    }

    /*Finalize clustering*/
    /*Acutally copying */
    lu_dynamicArray tmp_centers;
    lu_dynamicArray tmp_cluster_idx;


    DynamicArray_Init(&tmp_centers);
    DynamicArray_Init(&tmp_cluster_idx);

    DynamicArray_Reserve(&tmp_centers, nclus);
    DynamicArray_Reserve(&tmp_cluster_idx, nclus);

    size_t final_cluster_count = 0;
    for(size_t i = 0; i < nclus; ++i)
    {
        final_cluster_count += surviving_clusters[i];
        if(surviving_clusters[i]){
            DynamicArray_pushBack(&tmp_centers, cluster->centers.data[i]);
            DynamicArray_pushBack(&tmp_cluster_idx, i);
        }
    }

    printf("Surviving clusters %lu\n", final_cluster_count);
    /*allocate auxiliary pointers to store results of the finalization of the procedure*/

    FLOAT_TYPE** tmp_border_density = (FLOAT_TYPE**)malloc(final_cluster_count*sizeof(FLOAT_TYPE*));
    FLOAT_TYPE** tmp_border_err     = (FLOAT_TYPE**)malloc(final_cluster_count*sizeof(FLOAT_TYPE*));
    size_t**     tmp_border_idx     = (size_t**)malloc(final_cluster_count*sizeof(size_t*));

    FLOAT_TYPE* tmp_border_density_data = (FLOAT_TYPE*)malloc(final_cluster_count*final_cluster_count*sizeof(FLOAT_TYPE));
    FLOAT_TYPE* tmp_border_err_data     = (FLOAT_TYPE*)malloc(final_cluster_count*final_cluster_count*sizeof(FLOAT_TYPE)); 
    size_t* tmp_border_idx_data         = (size_t*)malloc(final_cluster_count*final_cluster_count*sizeof(size_t*));

    /*initialize all pointers*/
    for(size_t i = 0; i < final_cluster_count; ++i)
    {
        tmp_border_density[i]   = tmp_border_density_data + i*final_cluster_count;
        tmp_border_err[i]       = tmp_border_err_data + i*final_cluster_count;
        tmp_border_idx[i]       = tmp_border_idx_data + i*final_cluster_count;
    }
    LinkedList * tmp_cluster_list = (LinkedList*)malloc((final_cluster_count+1)*sizeof(LinkedList));

    /*Fix cluster assignment*/
    #pragma omp parallel for
    for(size_t c = 0; c < final_cluster_count; ++c)
    {
        size_t c_idx = tmp_cluster_idx.data[c];
        /*Put the correct ll head in tmp_cluster_list*/
        tmp_cluster_list[c] = (cluster -> clusters[c_idx]);
        /*walk the linked list and fix cluster assignment*/
        Node * n = tmp_cluster_list[c].head;
        while(n){
            particles[n -> data].cluster_idx = c;
            n = (n -> next);
        }
        /*Fix border matrices*/
        for(size_t d = c; d < final_cluster_count; ++d)
        {
            size_t c_jdx = tmp_cluster_idx.data[d];
            tmp_border_density[c][d] = border_density[c_idx][c_jdx];
            tmp_border_density[d][c] = border_density[c_idx][c_jdx];

            tmp_border_idx[c][d] = border_idx[c_idx][c_jdx];
            tmp_border_idx[d][c] = border_idx[c_idx][c_jdx];


            tmp_border_err[c][d] = border_err[c_idx][c_jdx];
            tmp_border_err[d][c] = border_err[c_idx][c_jdx];
        } 


    }

    Clusters_Reset(cluster);
    cluster -> border_density = tmp_border_density;
    cluster -> border_idx = tmp_border_idx;
    cluster -> border_err = tmp_border_err;

    cluster -> __border_density_data = tmp_border_density_data;
    cluster -> __border_idx_data = tmp_border_idx_data;
    cluster -> __border_err_data = tmp_border_err_data;

    cluster -> clusters = tmp_cluster_list;
    cluster -> centers = tmp_centers;
    /*Halo*/
    switch (halo)
    {
    case 1:
        for(size_t c = 0; c < final_cluster_count; ++c)
        {
            FLOAT_TYPE max_border_den = -2.;
            for(size_t d = 0; d < final_cluster_count; ++d)
            {
                if(tmp_border_density[c][d] > max_border_den)
                {
                    max_border_den = tmp_border_density[c][d];
                }
            }
            Node * n = tmp_cluster_list[c].head;
            /*if particle n has density lower than the one on the border put it into the halo*/
            while(n)
            {
                if(particles[n -> data].log_rho_c < max_border_den) particles[n -> data].cluster_idx = -1;
                n = n -> next;
            }
        }
        break;
    
    default:
        break;
    }    

    /*free memory and put the correct arrays into place*/
    free(tmp_cluster_idx.data);
    free(ipos.data);
    free(jpos.data);
    free(surviving_clusters);
}


