#ifndef MOTHER_H
#define MOTHER_H

#include "iip_blas_lv1.h"
#include "iip_matrix.h"

#if USE_CUDA
//cublasHandle_t handle;
//UINT max_thread;
//UINT max_block;

class cublas_setting
{
	cudaDeviceProp prop;
	
	public:
	cublas_setting()
	{
	cublasCreate(&handle);
/*
Assume that we always use first GPU device
What if, internal gpu and external gpu, 2 gpus exist?
*/
	cudaGetDeviceProperties(&prop,0);

	max_thread = prop.maxThreadsDim[0];
	max_block = prop.maxGridSize[1];

//printf("%d %d\n",max_thread,max_block);
	}
	~cublas_setting()
	{
	cublasDestroy(handle);
	}
};

static cublas_setting aa;

#endif



#endif
