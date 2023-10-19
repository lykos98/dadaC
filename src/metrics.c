#include "../include/dadac.h"
#include <stdint.h>

uint32_t METRICS_DATADIMS;


float_t eud(void* a, void* b)
{
	float_t* aa = (float_t*)a;
	float_t* bb = (float_t*)b;
	float_t acc = 0;
	for(uint32_t i=0; i < METRICS_DATADIMS; ++i) 
	{
		float_t dd = (aa[i] - bb[i]);
		acc += dd*dd;
	}
	return sqrt(acc); 
	//return acc; 
}
