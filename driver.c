#include "./include/clustering.h"
#include <stdio.h>
#include <time.h>

unsigned int data_dims;

void write_border_idx(const char * fname, Clusters * c)
{
    FILE * f = fopen(fname, "w");
    for(int i = 0; i < c -> centers.count; ++i)
    {
        for(int j = 0; j < c -> centers.count; ++j)
        {
            int a = c -> border_idx[i][j] == NOBORDER ? -1 : c -> border_idx[i][j];
            fprintf(f,"%d ",a);
        }
        fprintf(f,"\n");
    }
    fclose(f);
}

void write_point_info(const char * fname, Datapoint_info * particles, size_t n)
{
    FILE * f = fopen(fname,"w");
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
}

int main(int argc, char** argv){

    int print_results = 0;
    double Z;
    int halo;
    char aux_fname[80];

    data_dims = 5;

    /***********************************************************************
     * TODO: Make a function to perform KNN search, fix verbose and timing *
     ***********************************************************************/
    if(argc < 5 )
    {
        printf("USAGE: ./driver [INPUT_FILE] [OUTPUT_FILE] [Z] [HALO]");
        printf("\nThe program gives as output the cluster assignment of each datapoint\n");
        return;
    }
    else
    {
        Z = atof(argv[3]);
        halo = atoi(argv[4]);
        if(halo != 0 && halo != 1){
            printf("Insert valid halo identifier: 0 do not assign halo, 1 assign particles to the halo\n");
            return;
        }
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

    struct timespec start, finish;
    double elapsed;


    kd_node* root = build_tree(kd_ptrs, n);

    clock_gettime(CLOCK_MONOTONIC, &finish);

    printf("The root of the tree is\n");
    printKDnode(root);

    int k = 1001;

    Datapoint_info* particles = (Datapoint_info*)malloc(n*sizeof(Datapoint_info));

    /**************
     * KNN search *
     **************/

    KNN_search(particles,data, root, n, k);

    /********************************
     * Intrinsic Dimension estimate *
     ********************************/

    double id = idEstimate(particles, n);

    /***********************
     * Density computation *
     ***********************/
    computeRho(particles,id,n);
    calculateCorrection(particles,n,Z);

    /********************
     * First clustering *
     ********************/

    Clusters c = Heuristic1(particles, data, n);

    /***************************************************************************************
     * Allocate borders and other things to store clustering info                          *
     * Then Find borders between clusters and then merge clusters using peaks significance *
     ***************************************************************************************/
    Clusters_allocate(&c);  

    // sprintf(aux_fname, "%s_int", argv[2]);
    // write_point_info(aux_fname,particles,n);

    Heuristic2(&c, particles);

    //sprintf(aux_fname, "%s_bord_int", argv[2]);
    //write_border_idx(aux_fname,&c);
    
    Heuristic3(&c, particles, Z, halo);

    sprintf(aux_fname, "%s_bord", argv[2]);
    write_border_idx(aux_fname,&c);

    clock_gettime(CLOCK_MONOTONIC, &start);
    write_point_info(argv[2],particles,n);
    clock_gettime(CLOCK_MONOTONIC, &finish);

    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;

    printf("Writing results: %.3lf\n",elapsed);
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