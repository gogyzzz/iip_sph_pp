#ifndef MOTHER_H
#define MOTHER_H

#include "iip_blas_lv1.h"
#include "iip_matrix.h"
#include "iip_wav.h"

void init()
{

#if USE_CUDA
	cudaDeviceProp prop;
	cublasCreate(&handle);
	cudaGetDeviceProperties(&prop,0);
	max_thread = prop.maxThreadsDim[0];
	max_block = prop.maxGridSize[1];
#endif
}

void finit()
{
#if USE_CUDA
	cublasDestroy(handle);
#endif
}

#endif
