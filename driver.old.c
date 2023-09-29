#include "./include/clustering.h"
#include <stdio.h>
#include <time.h>

unsigned int data_dims;

void write_border_idx(const char * fname, Clusters * c)
{
    FILE * f = fopen(fname, "w");

    if(c -> UseSparseBorders)
    {
	    for(idx_t i = 0; i < c -> centers.count; ++i)
	    {
		for(idx_t el = 0; el < c ->SparseBorders[i].count; ++el)
		{
		    int a = c -> SparseBorders[i].data[el].idx; 
		    fprintf(f, "%d ",a);
		}
		fprintf(f,"\n");
	    }
    }
    else
    {
	    for(idx_t i = 0; i < c -> centers.count; ++i)
	    {
		for(idx_t j = 0; j < c -> centers.count; ++j)
		{
		    //int a = c -> borders[i][j].idx == NOBORDER ? -1 : c -> borders[i][j].idx;
		    //fprintf(f,"%d ",a);
		    int a = c -> borders[i][j].idx; 
		    if(c -> borders[i][j].idx != NOBORDER) fprintf(f, "%d ",a);
		}
		fprintf(f,"\n");
	    }
    }
    fclose(f);
}

void write_point_info(const char * fname, Datapoint_info * particles, idx_t n)
{
    FILE * f = fopen(fname,"w");
    for(idx_t i = 0; i < n; ++i)
    {
        fprintf(f,"%lu\t",(uint64_t)particles[i].kstar);
        fprintf(f,"%d\t", particles[i].cluster_idx);
	#ifdef USE_FLOAT32
		fprintf(f,"%.6f\t",particles[i].log_rho);
	#else
		fprintf(f,"%.11lf\t",particles[i].log_rho);
	#endif
        //fprintf(f,"%.12lf\t",particles[i].log_rho_c);
        //fprintf(f,"%.12lf\t",particles[i].log_rho_err);
        //fprintf(f,"%.12lf\t",particles[i].g);
        fprintf(f,"%d\t",particles[i].is_center);
        fprintf(f,"\n");
    }
    fclose(f);
}

int main(int argc, char** argv){

    double Z;
    int halo;
    char aux_fname[80];
    int k;
    int UseSparseBorders;

    struct timespec start, finish;
    double elapsed;
    struct timespec start_tot, finish_tot;
    double elapsed_tot;
    //Start timer
    clock_gettime(CLOCK_MONOTONIC, &start_tot);

    data_dims = 2;

    /***********************************************************************
     * TODO: Make a function to perform KNN search, fix verbose and timing *
     ***********************************************************************/

    /* default args */
    Z		     = 2;
    halo	     = 1;
    k		     = 1001;
    UseSparseBorders = 1;


    /* parsing */


    if(argc < 6 )
    {
		printf("USAGE: ./driver [INPUT_FILE] [OUTPUT_FILE] [Z] [HALO] [k] [s]\n");
		printf("\tINPUT_FILE : input file, file path\n");
		printf("\tOUTPUT_FILE: output file, file path\n");
		printf("\tZ	     : Z value, float\n");
		printf("\tHALO	     : Assign halo, bool 0/1 \n");
		printf("\tk	     : Number of neighbors to use, int (>0) \n");
		printf("\ts	     : Use sparse borders implementation, s/d [sparse/dense]\n");
		printf("\nThe program gives as output the cluster assignment of each datapoint\n");
		return 0;
    }
    else
    {
        Z = atof(argv[3]);
        halo = atoi(argv[4]);
	if(argc ==  6)
	{
		k = atoi(argv[5]);
	}
	else
	{
		k = 1001;
	}

        if(halo != 0 && halo != 1){
            printf("Insert valid halo identifier: 0 do not assign halo, 1 assign particles to the halo\n");
            return 0;
        }
    }

    printf("Using: \n\t Z    	   : %.2lf \n\t Halo      : %s \n\t Neighbors : %d \n\n", Z, halo ? "yes" : "no", k);
    
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
    printf("Reading %lu particles\n",(uint64_t)n);


    FLOAT_TYPE* data = (FLOAT_TYPE*)malloc(data_dims*n*sizeof(FLOAT_TYPE));
    float* df = (float*)malloc(data_dims*n*sizeof(float));
    size_t fff = fread(df,sizeof(float),data_dims*n,f);
    printf("Read %luB\n",fff);
    fclose(f);

    for(idx_t i = 0; i < n*data_dims; ++i) data[i] = (FLOAT_TYPE)(df[i]);

    free(df);



	Datapoint_info* particles = NgbhSearch(data, n, data_dims, k); 
    /********************************
     * Intrinsic Dimension estimate *
     ********************************/

    double id = idEstimate(particles, n);

    /***********************
     * Density computation *
     ***********************/
    computeRho(particles,id,n);
    computeCorrection(particles,n,Z);

    /********************
     * First clustering *
     ********************/

    Clusters c = Heuristic1(particles, data, n);

    /***************************************************************************************
     * Allocate borders and other things to store clustering info                          *
     * Then Find borders between clusters and then merge clusters using peaks significance *
     ***************************************************************************************/
   // Clusters_allocate(&c);  
    Clusters_allocate(&c, argv[6][0] == 's');  

    // sprintf(aux_fname, "%s_int", argv[2]);
    // write_point_info(aux_fname,particles,n);

    Heuristic2(&c, particles);

    //sprintf(aux_fname, "%s_bord_int", argv[2]);
    //write_border_idx(aux_fname,&c);

    c.n = n;
    
    Heuristic3(&c, particles, Z, halo);

    //sprintf(aux_fname, "%s_bord", argv[2]);
    //write_border_idx(aux_fname,&c);

    clock_gettime(CLOCK_MONOTONIC, &start);
    write_point_info(argv[2],particles,n);
    clock_gettime(CLOCK_MONOTONIC, &finish);

    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;

    printf("Writing results: %.3lf\n",elapsed);
    /*******************
     * Free all memory *
     *******************/

	freeDatapointArray(particles,n);
    free(data);
    Clusters_free(&c);

    clock_gettime(CLOCK_MONOTONIC, &finish_tot);
    elapsed_tot = (finish_tot.tv_sec - start_tot.tv_sec);
    elapsed_tot += (finish_tot.tv_nsec - start_tot.tv_nsec) / 1000000000.0;
    printf("ELAPSED time (measured by driver): %.3lfs\n\n", elapsed_tot);

    return 0;
}
