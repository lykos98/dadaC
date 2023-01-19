#include "./include/clustering.h"
#include <stdio.h>

unsigned int data_dims = 5;

int main(int argc, char** argv){
    
    FILE* f = fopen("/home/francesco/Desktop/dssc/tirocinio/datafiles/0001_std","r");
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

    kd_node* root = make_tree(kd_ptrs, 0, n-1, NULL ,0);

    printf("The root of the tree is\n");
    printKDnode(root);

    int k = 1001;


    Datapoint_info* particles = (Datapoint_info*)malloc(n*sizeof(Datapoint_info));

    /**************
     * KNN search *
     **************/

    #pragma omp parallel for
    for(int p = 0; p < n; ++p)
    {
        particles[p].ngbh = KNN(data + data_dims*p, root, k);
        particles[p].array_idx = p;
    }

//    f = fopen("/home/francesco/Desktop/dssc/tirocinio/datafiles/ngbh","w");
//    for(size_t i = 0; i < n; ++i)
//    {
//        for(size_t j = 0; j < k; ++j) fprintf(f,"%lu\t",particles[i].ngbh.data[j].array_idx);
//
//        fprintf(f,"\n");
//    }
//    fclose(f);
//
    double id = idEstimate(particles, n);
    printf("Instrinsic dimension %lf\n",id);

    computeRho(particles,id,n);
    calculateCorrection(particles,n,1.96);
    Clusters c = Heuristic1(particles, data, n);

    Clusters_allocate(&c);  

    Heuristic2(&c, particles);

    //size_t nclus = c.centers.count; 
    //f = fopen("/home/francesco/Desktop/dssc/tirocinio/datafiles/brdrs","w");
    //for(size_t i = 0; i < nclus; ++i)
    //{
    //    for(size_t j = 0; j < nclus; ++j) c.border_idx[i][j] == NOBORDER ? fprintf(f,"%d\t",-1) : fprintf(f,"%d\t",(int)c.border_idx[i][j]);
    //    fprintf(f,"\n");
    //}
    //fclose(f);
    
    //f = fopen("/home/francesco/Desktop/dssc/tirocinio/datafiles/bden","w");
    //for(size_t i = 0; i < nclus; ++i)
    //{
    //    for(size_t j = 0; j < nclus; ++j) fprintf(f,"%lf\t",c.border_density[i][j]);
    //    fprintf(f,"\n");
    //}
    //fclose(f);

    Heuristic3(&c, particles, 1.96, 0);


    //f = fopen("/home/francesco/Desktop/dssc/tirocinio/datafiles/nope","w");
    //for(size_t i = 0; i < n; ++i)
    //{
    //    fprintf(f,"%lu\t",particles[i].kstar);
    //    fprintf(f,"%d\t", particles[i].cluster_idx);
    //    fprintf(f,"%.12lf\t",particles[i].log_rho);
    //    fprintf(f,"%.12lf\t",particles[i].log_rho_c);
    //    fprintf(f,"%.12lf\t",particles[i].log_rho_err);
    //    fprintf(f,"%.12lf\t",particles[i].g);
    //    fprintf(f,"\n");
    //}
    //fclose(f);




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