#include "./include/dadac.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <time.h>

struct Options 
{
	char* inputFile;
	char* outputFile;
    double Z;
    int halo;
    int k;
    int UseSparseBorders;
    int FileInFloat32;
    unsigned int data_dims;

};

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

void printHelp()
{
	printf("USAGE: ./driver i=[INPUT_FILE] o=[OUTPUT_FILE] d=[d] t=[t] z=[Z] h=[HALO] k=[k] s=[s] t=[t]\n");
	printf("\tINPUT_FILE : input file, file path\n");
	printf("\tOUTPUT_FILE: output file, file path\n");
	printf("\td	     : Lenght of the data vectors (int) \n");
	printf("\tZ	     : Z value, float\n");
	printf("\tHALO	     : Assign halo, y/n [yes/no]\n");
	printf("\tk	     : Number of neighbors to use, int (>0) \n");
	printf("\ts	     : Use sparse borders implementation, y/n [sparse/dense]\n");
	printf("\tt	     : Input binary is in Float32, y/n [float/double]\n");
	printf("\nThe program gives as output the cluster assignment of each datapoint\n");
	return;
}

void checkEq(char c)
{
	if(c != '='){
		printf("Wrongly formatted input\n");
		printHelp();
		exit(1);
	};
	return;
}

struct Options Parser(int argc, char** argv)
{
	struct Options opt = {.inputFile=NULL, .outputFile=NULL,.FileInFloat32 = 1,.k=1001, .halo = 1, .Z = 2, .data_dims = 0, .UseSparseBorders = 0 };
	if(argc < 2)
	{
		printHelp();
		exit(1);
	}
	for(int i = 1; i < argc; ++i)	
	{
		checkEq(argv[i][1]);
		switch(argv[i][0])
		{
			case 'i':
				opt.inputFile = argv[i] + 2;
				break;
			case 'o':
				opt.outputFile = argv[i] + 2;
				break;
			case 't':
				opt.FileInFloat32 = argv[i][2] == 'y' ? 1 : 0;
				break;
			case 'd':
				opt.data_dims = atoi(argv[i] + 2);
				break;
			case 'k':
				opt.k = atoi(argv[i] + 2);
				break;
			case 'h':
				opt.halo = argv[i][2] == 'y' ? 1 : 0;
				break;
			case 's':
				opt.UseSparseBorders = argv[i][2] == 'y' ? 1 : 0;
				break;
			case 'z':
				opt.Z = atof(argv[i] + 2);
				break;
			default:
				printHelp();
				exit(1);
				break;
		}
	}
	if(!(opt.inputFile) || !(opt.outputFile)){
		printf("Please provide input and output file paths\n");
		printHelp();
		exit(1);
	}
	if(opt.data_dims == 0)
	{
		printf("Please specify the lenght of each vector\n");
		printHelp();
		exit(1);
	}
	return opt;
}

float_t* generateRandomMatrix(idx_t ncols,idx_t nrows)
{
	float_t* mat = (float_t*)malloc(sizeof(float_t)*ncols*nrows);
	for(idx_t i = 0; i < nrows*ncols; ++i) mat[i] = ((float_t)rand())/((float_t)RAND_MAX);
	return mat;
}


int main(){

    //char aux_fname[80];

    //struct timespec start, finish;
    //double elapsed;
    //struct timespec start_tot, finish_tot;
    //double elapsed_tot;
	
	

	size_t n;
	float_t* data;
	uint32_t data_dims, k;

	//n = 35947;
	n = 20000;
	data_dims = 100;
	k = 600;
	//FILE* f = fopen("Bunny.txt","r");
	//data = (float_t*)malloc(n*data_dims*sizeof(float_t));
	//for(int i = 0; i < n; ++i)
	//{
	//	float x, y, z;
	//	fscanf(f,"%f %f %f\n",&x, &y, &z);
	//	data[i*data_dims] = (float_t)x;
	//	data[i*data_dims + 1] = (float_t)y;
	//	data[i*data_dims + 2] = (float_t)z;
	//}
	data = generateRandomMatrix(n, data_dims);

//	n = 2000;
//	data_dims = 2;
//	METRICS_DATADIMS = data_dims;
//	FILE* f = fopen("dummy.dat","r");
//	float* dd = (float*)malloc(4*n*data_dims);
//	data = (float_t*)malloc(sizeof(float_t)*n*data_dims);
//	fread(dd, 4, n*data_dims, f);
//	fclose(f);
//	for(int i=0; i < n*data_dims; ++i) data[i] = (float_t)dd[i];



	
    //Start timer

	Datapoint_info* particles = NgbhSearch_kdtree_V2(data, n, data_dims, k); 
	Datapoint_info* pp = NgbhSearch_bruteforce(data,n ,sizeof(float_t),data_dims, k, NULL); 

	
	for(idx_t idx=0; idx < n; ++idx)
	for(size_t i = 0; i < k; ++i)
	{
		if(particles[idx].ngbh.data[i].array_idx != pp[idx].ngbh.data[i].array_idx ) {
			printf("particle %lu --> ngbh %lu --> got m1 %lu %lf -- m2 %lu %lf\n", idx, i, 
													particles[idx].ngbh.data[i].array_idx,
													particles[idx].ngbh.data[i].value,
													pp[idx].ngbh.data[i].array_idx,
													pp[idx].ngbh.data[i].value);
		}
	}


	free(particles);
	free(pp);
	free(data);
    return 0;
}
