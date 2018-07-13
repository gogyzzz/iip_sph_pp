#include "iip_blas_lv1.h"

#if DEBUG
#define CUDA_CALL(x) \
{ \
	const cudaError_t a = (x); \
	if(a != cudaSuccess) { \
		printf("\nCuda Error: %s (err_num=%d) at line:%d\n", cudaGetErrorString(a), a, __LINE__); \
		cudaDeviceReset(); exit(0); \
	} \
}
#else
#define CUDA_CALL(x) {(x);}
#endif


/*** AXPY ***/
void axpy(DTYPE alpha, MAT *x, MAT *y) // for real mat
{
#if DEBUG
	printf("%s\n",__func__);
#endif
	UINT mat_size = x->d0 * x->d1 * x->d2;
	printf("max_thread : %d\n",max_thread);
	UINT num_block = (UINT)(mat_size/(UINT)max_thread)+1;
	
	cu_axpy<<<num_block,max_thread>>>(alpha, x->data, 1, y->data, 1, mat_size-1, max_thread);
	CUDA_CALL(cudaThreadSynchronize())
}

__global__ void cu_axpy(DTYPE alpha, DTYPE *X, UINT INCX, DTYPE *Y, UINT INCY, UINT len, UINT block_size)
{
	ITER idx = threadIdx.x + blockIdx.x * block_size;
	if(idx > len)
		return;
	Y[idx * INCY] += X[idx * INCX] * alpha;
}


void caxpy(CTYPE alpha, CMAT *x, CMAT *y) // for real mat
{
#if DEBUG
	printf("%s\n",__func__);
#endif
	UINT mat_size = x->d0 * x->d1 * x->d2;
	printf("max_thread : %d\n",max_thread);
	UINT num_block = (UINT)(mat_size/(UINT)max_thread)+1;
	
	cu_caxpy<<<num_block,max_thread>>>(alpha, x->data, 1, y->data, 1, mat_size-1, max_thread);
	CUDA_CALL(cudaThreadSynchronize())
}

__global__ void cu_caxpy(CTYPE alpha, CTYPE *X, UINT INCX, CTYPE *Y, UINT INCY, UINT len, UINT block_size)
{
	ITER idx = threadIdx.x + blockIdx.x * block_size;
	if(idx > len)
		return;
	Y[idx * INCY].re += X[idx * INCX].re * alpha.re;
	Y[idx * INCY].im += X[idx * INCX].im * alpha.im;
}


/*** COPY ***/
void copy(MAT *src, MAT *des)
{
#if DEBUG
	printf("%s\n",__func__);
#endif
	UINT mat_size = src->d0 * src->d1 * src->d2;
	printf("max_thread : %d\n",max_thread);
	UINT num_block = (UINT)(mat_size/(UINT)max_thread)+1;
	
	cu_copy<<<num_block,max_thread>>>(src->data, 1, des->data, 1, mat_size-1, max_thread);
	CUDA_CALL(cudaThreadSynchronize())
}

__global__ void cu_copy(DTYPE *SRC, UINT INC_SRC, DTYPE *DES, UINT INC_DES, UINT len, UINT block_size)
{
	ITER idx = threadIdx.x + blockIdx.x * block_size;
	if(idx > len)
		return;
	DES[idx * INC_DES] = SRC[idx * INC_SRC];
}

void ccopy(CMAT *src, CMAT *des)
{
#if DEBUG
	printf("%s\n",__func__);
#endif
	UINT mat_size = src->d0 * src->d1 * src->d2;
	printf("max_thread : %d\n",max_thread);
	UINT num_block = (UINT)(mat_size/(UINT)max_thread)+1;
	
	cu_ccopy<<<num_block,max_thread>>>(src->data, 1, des->data, 1, mat_size-1, max_thread);
	CUDA_CALL(cudaThreadSynchronize())
}

__global__ void cu_ccopy(CTYPE *SRC, UINT INC_SRC, CTYPE *DES, UINT INC_DES, UINT len, UINT block_size)
{
	ITER idx = threadIdx.x + blockIdx.x * block_size;
	if(idx > len)
		return;
	DES[idx * INC_DES].re = SRC[idx * INC_SRC].re;
	DES[idx * INC_DES].im = SRC[idx * INC_SRC].im;
}



/*** SUM ***/
DTYPE asum(MAT *mat, UINT inc)
{
#if DEBUG
	printf("%s\n",__func__);
#endif
	UINT mat_size = mat->d0 * mat->d1 * mat->d2;
	printf("max_thread : %d\n",max_thread);
	UINT num_block = (UINT)(mat_size/(UINT)max_thread)+1;

	DTYPE sum=0;
	
	cu_asum<<<num_block,max_thread>>>(mat->data, 1, mat_size-1, max_thread, &sum);
	CUDA_CALL(cudaThreadSynchronize())

	return sum;
}

__global__ void cu_asum(DTYPE *data, UINT inc, UINT len, UINT block_size, DTYPE *sum)
{
	ITER idx = threadIdx.x + blockIdx.x * block_size;
	__shared__ DTYPE s_cnt;

	if(idx > len)
		return;

	if(threadIdx.x == 0){
		s_cnt = 0;
	}

	__syncthreads();
	atomicAdd(&s_cnt, 1);
	__syncthreads();

	if(threadIdx.x == 0){
		//atomicAdd(sum, s_cnt);
	}
}

DTYPE casum(CMAT *mat, UINT inc)
{
#if DEBUG
	printf("%s\n",__func__);
#endif
	UINT mat_size = mat->d0 * mat->d1 * mat->d2;
	printf("max_thread : %d\n",max_thread);
	UINT num_block = (UINT)(mat_size/(UINT)max_thread)+1;

	CTYPE sum;
	sum.re = 0;
	sum.im = 0;

	cu_casum<<<num_block,max_thread>>>(mat->data, 1, mat_size-1, max_thread, &sum);
	CUDA_CALL(cudaThreadSynchronize())

	return sum.re + sum.im;
}

__global__ void cu_casum(CTYPE *data, UINT inc, UINT len, UINT block_size, CTYPE *sum)
{
	ITER idx = threadIdx.x + blockIdx.x * block_size;
	__shared__ int s_cnt_re;
	__shared__ int s_cnt_im;

	if(idx > len)
		return;

	if(threadIdx.x == 0){
		s_cnt_re = 0;
		s_cnt_im = 0;
	}
	
	__syncthreads();
	//atomicAdd(&s_cnt_re, data[idx * inc].re);
	//atomicAdd(&s_cnt_im, data[idx * inc].im);
	__syncthreads();

	if(threadIdx.x == 0){
		/*atomicAdd(&(*sum).re, s_cnt_re);
		atomicAdd(&(*sum).im, s_cnt_im);*/
	}
}