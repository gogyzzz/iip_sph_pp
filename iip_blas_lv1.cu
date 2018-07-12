#include "iip_blas_lv1.h"

#if DEBUG
#define CUDA_CALL(x) \
{ \
	const cudaError_t a = (x); \
	if(a != cudaSuccess) { \
		printf("\nCuda Error: %s (err_num=%d) at line:%d\n", cudaGetErrorString(a), a, __LINE__); \
		cudaDeviceReset(); assert(0); \
	} \
}
#else
#define CUDA_CALL(x) {(x);}
#endif


/*
void copy(MAT* src, SINT src_increment, MAT* des, SINT des_increment) // for real mat
{
#if DEBUG
  printf("%s\n",__func__);
#endif
	UINT mat_size = src->d0 * src->d1 * src->d2;
	printf("max_thread : %d\n",max_thread);
	UINT num_block = (UINT)(len/(UINT)max_thread)+1;
	
	cu_copy<<<num_block,max_thread>>>(mat->data,len-1,val,max_thread);
	CUDA_CALL(cudaThreadSynchronize())
}

__global__ void cu_copy(DTYPE* data, UINT len,DTYPE val,UINT size_block)
{
	ITER idx = threadIdx.x + blockIdx.x * size_block;
	if(idx > len)
		return;
  data[idx]= val;
}
*/
