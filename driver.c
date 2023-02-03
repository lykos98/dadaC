#include "./include/clustering.h"
#include <stdio.h>
#include <time.h>

unsigned int data_dims = 5;

int main(int argc, char** argv){

    int print_results = 0;

    /*************************************************************
     * TODO: Boilerplate code for reading the datafile from argv *
     *************************************************************/
    if(argc < 3 )
    {
        printf("USAGE: ./driver [INPUT_FILE] [OUTPUT_FILE]");
        printf("\nThe program gives as output the cluster assignment of each datapoint\n");
        return;
    }
    
    FILE* f = fopen(argv[1],"r");
    if(!f)
    {
        printf("Nope\n");
        exit(1);
    }
    fseek(f,0,SEEK_END);
    size_t n = ftell(f);
    rewind(f);

    n = n/(sizeof(float)*data_dims);
    printf("Reading %lu particles\n",n);


    FLOAT_TYPE* data = (FLOAT_TYPE*)malloc(data_dims*n*sizeof(FLOAT_TYPE));
    float* df = (float*)malloc(data_dims*n*sizeof(float));
    fread(df,sizeof(float),data_dims*n,f);
    fclose(f);

    for(size_t i = 0; i < n*data_dims; ++i) data[i] = (FLOAT_TYPE)(df[i]);

    free(df);

    kd_node* kd_node_array = (kd_node*)malloc(n*sizeof(kd_node));
    kd_node** kd_ptrs = (kd_node**)malloc(n*sizeof(kd_node*));

    initializeKDnodes(kd_node_array,data,n);
    initializePTRS(kd_ptrs, kd_node_array,n);

    data_dims = 5;

    struct timespec start, finish;
    double elapsed;
    /*Making the tree*/

    clock_gettime(CLOCK_MONOTONIC, &start);

    kd_node* root = make_tree(kd_ptrs, 0, n-1, NULL ,0);

    clock_gettime(CLOCK_MONOTONIC, &finish);

    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;

    printf("%.3lfs for building the KDtree\n",elapsed);

    printf("The root of the tree is\n");
    printKDnode(root);

    int k = 1001;

    Datapoint_info* particles = (Datapoint_info*)malloc(n*sizeof(Datapoint_info));

    /**************
     * KNN search *
     **************/

    clock_gettime(CLOCK_MONOTONIC, &start);
    #pragma omp parallel for
    for(int p = 0; p < n; ++p)
    {
        particles[p].ngbh = KNN(data + data_dims*p, root, k);
        particles[p].array_idx = p;
    }

    clock_gettime(CLOCK_MONOTONIC, &finish);

    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;

    printf("%.3lfs KNN search\n",elapsed);

    /********************************
     * Intrinsic Dimension estimate *
     ********************************/

    clock_gettime(CLOCK_MONOTONIC, &start);

    double id = idEstimate(particles, n);

    printf("Instrinsic dimension %lf\n",id);

    clock_gettime(CLOCK_MONOTONIC, &finish);

    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;

    printf("%.3lfs id estimation\n",elapsed);
    /***********************
     * Density computation *
     ***********************/
    clock_gettime(CLOCK_MONOTONIC, &start);
    computeRho(particles,id,n);
    calculateCorrection(particles,n,1.96);

    clock_gettime(CLOCK_MONOTONIC, &finish);

    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;

    printf("%.3lfs Density and kstar computation\n",elapsed);
    /********************
     * First clustering *
     ********************/

    clock_gettime(CLOCK_MONOTONIC, &start);
    Clusters c = Heuristic1(particles, data, n);

    clock_gettime(CLOCK_MONOTONIC, &finish);

    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;

    printf("%.3lfs H1\n",elapsed);
    /***************************************************************************************
     * Allocate borders and other things to store clustering info                          *
     * Then Find borders between clusters and then merge clusters using peaks significance *
     ***************************************************************************************/

    clock_gettime(CLOCK_MONOTONIC, &start);
    Clusters_allocate(&c);  

    

    Heuristic2(&c, particles);

    clock_gettime(CLOCK_MONOTONIC, &finish);

    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;

    printf("%.3lfs H2\n",elapsed);
    
    clock_gettime(CLOCK_MONOTONIC, &start);
    Heuristic3(&c, particles, 1.96, 0);
    
    clock_gettime(CLOCK_MONOTONIC, &finish);

    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;

    printf("%.3lfs H3\n",elapsed);

    clock_gettime(CLOCK_MONOTONIC, &start);
    f = fopen(argv[2],"w");
    for(size_t i = 0; i < n; ++i)
    {
        fprintf(f,"%lu\t",particles[i].kstar);
        fprintf(f,"%d\t", particles[i].cluster_idx);
        fprintf(f,"%.12lf\t",particles[i].log_rho);
        //fprintf(f,"%.12lf\t",particles[i].log_rho_c);
        //fprintf(f,"%.12lf\t",particles[i].log_rho_err);
        //fprintf(f,"%.12lf\t",particles[i].g);
        fprintf(f,"%d\t",particles[i].is_center);
        fprintf(f,"\n");
    }
    fclose(f);

    clock_gettime(CLOCK_MONOTONIC, &finish);

    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;

    printf("%.3lfs writing results \n",elapsed);

    /*******************
     * Free all memory *
     *******************/

    for (size_t i = 0; i < n; ++i)
    {
        freeHeap(&particles[i].ngbh);
    }




    free(particles);
    free(kd_ptrs);
    free(kd_node_array);
    free(data);
    Clusters_free(&c);

}