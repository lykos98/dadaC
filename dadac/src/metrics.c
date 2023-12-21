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

float_t eud_sq2(void* a, void* b)
{
	float_t* aa = (float_t*)a;
	float_t* bb = (float_t*)b;
	float_t acc = 0;
	for(uint32_t i=0; i < METRICS_DATADIMS; ++i) 
	{
		float_t dd = (aa[i] - bb[i]);
		acc += dd*dd;
	}
	return acc; 
	//return acc; 
}

float_t eudOpt(void* uvoid, void* vvoid)
{
    float_t s;
    uint32_t i = 0;
	float_t* u = (float_t*)uvoid;
	float_t* v = (float_t*)vvoid;
    // manually unrolled loop, might be vectorized
    float_t acc[4] = {0., 0., 0., 0.};
    for (; i + 4 <= METRICS_DATADIMS; i += 4) {
        float_t _u[4] = {u[i], u[i + 1], u[i + 2], u[i + 3]};
        float_t _v[4] = {v[i], v[i + 1], v[i + 2], v[i + 3]};
        float_t diff[4] = {_u[0] - _v[0],
                               _u[1] - _v[1],
                               _u[2] - _v[2],
                               _u[3] - _v[3]};
        acc[0] += diff[0] * diff[0];
        acc[1] += diff[1] * diff[1];
        acc[2] += diff[2] * diff[2];
        acc[3] += diff[3] * diff[3];
    }
    s = acc[0] + acc[1] + acc[2] + acc[3];
    if (i < METRICS_DATADIMS) {
        for(; i<METRICS_DATADIMS; ++i) {
            float_t d = u[i] - v[i];
            s += d * d;
        }
    }
    return sqrt(s);
}

float_t eud_sq(void* uvoid, void* vvoid)
{
    float_t s;
    uint32_t i = 0;
	float_t* u = (float_t*)uvoid;
	float_t* v = (float_t*)vvoid;
    // manually unrolled loop, might be vectorized
    float_t acc[4] = {0., 0., 0., 0.};
    for (; i + 4 <= METRICS_DATADIMS; i += 4) {
        float_t _u[4] = {u[i], u[i + 1], u[i + 2], u[i + 3]};
        float_t _v[4] = {v[i], v[i + 1], v[i + 2], v[i + 3]};
        float_t diff[4] = {_u[0] - _v[0],
                               _u[1] - _v[1],
                               _u[2] - _v[2],
                               _u[3] - _v[3]};
        acc[0] += diff[0] * diff[0];
        acc[1] += diff[1] * diff[1];
        acc[2] += diff[2] * diff[2];
        acc[3] += diff[3] * diff[3];
    }
    s = acc[0] + acc[1] + acc[2] + acc[3];
    if (i < METRICS_DATADIMS) {
        for(; i<METRICS_DATADIMS; ++i) {
            float_t d = u[i] - v[i];
            s += d * d;
        }
    }
    return s;
}
