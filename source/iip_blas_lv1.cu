/*
 * ===========================================================
 *           Copyright (c) 2018, __IIPLAB__
 *                All rights reserved.
 * 
 * This Source Code Form is subject to the terms of
 * the Mozilla Public License, v. 2.0. 
 * If a copy of the MPL was not distributed with this file,
 *  You can obtain one at http://mozilla.org/MPL/2.0/.
 * ===========================================================
 */
#include "iip_blas_lv1.h"
#include "cublas_v2.h"

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

/**** Atomic ADD for double(Not Work) ****/
/*
__device__ double atomicAdd_(double* address, double val) {
	unsigned long long int* address_as_ull = (unsigned long long int*)address;
	unsigned long long int old = *address_as_ull, assumed;
	do {
		assumed = old;
		old = atomicCAS(address_as_ull, assumed, __double_as_longlong(val + __longlong_as_double(assumed)));
		// Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN) 
	} 
	while (assumed != old);
	return __longlong_as_double(old);
}*/



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

	cublasHandle_t handle;

	if (cublasCreate(&handle) != CUBLAS_STATUS_SUCCESS){
		printf("Error to initialize cublas!\n");
		return 0;
	}

	DTYPE sum = 0;

	//cu_asum<<<num_block,max_thread>>>(mat->data, 1, mat_size-1, max_thread, &sum);
	//CUDA_CALL(cudaThreadSynchronize())

#if NTYPE == 1
	cublasDasum(handle, mat_size, mat->data, inc, &sum);
#elif NTYPE == 0
	cublasSasum(handle, mat_size, mat->data, inc, &sum);
#endif
	cublasDestroy(handle);
	return sum;

	//return val_for_sum;
}
/*
__global__ void cu_asum(DTYPE *data, UINT inc, UINT len, UINT block_size, DTYPE *sum)
{
	ITER idx = threadIdx.x + blockIdx.x * block_size;
	ITER idx_thr = threadIdx.x;
	__shared__ DTYPE s_cnt;

	if(idx > len)
		return;

	if(threadIdx.x == 0){
		s_cnt = 0;
		cnt = 0;
	}

	__syncthreads();
	printf("thr%d is waiting for cnt == %d.\n", idx, idx_thr);
	//wait for other threads
	while (cnt != idx_thr){
		printf("thr%d is waiting for cnt == %d.\n", idx, idx_thr);
	}

	printf("thr%d starting!\n", idx);
	s_cnt += data[idx*inc];
	if (blockIdx.x == 0){
		cnt++;
		printf("cnt increased : %d\n", cnt);
	}
	__syncthreads();

	printf("data[%d] = %.2f\n", idx*inc, data[idx*inc]);
	printf("s_cnt(%d)_before = %.2f\n", idx*inc, s_cnt);
	s_cnt = atomicAdd_(&s_cnt, 10.0);
	printf("s_cnt(%d)_after = %.2f\n", idx*inc, s_cnt);
	__syncthreads();

	if(threadIdx.x == 0){
		printf("running %d\n", threadIdx.x);
		atomicAdd_(&val_for_sum, s_cnt);
		printf("end %d\n", threadIdx.x);
	}


}
*/

DTYPE casum(CMAT *mat, UINT inc)
{
#if DEBUG
	printf("%s\n",__func__);
#endif
	UINT mat_size = mat->d0 * mat->d1 * mat->d2;
	printf("max_thread : %d\n",max_thread);
	UINT num_block = (UINT)(mat_size / (UINT)max_thread) + 1;

	cublasHandle_t handle;

	if (cublasCreate(&handle) != CUBLAS_STATUS_SUCCESS){
		printf("Error to initialize cublas!\n");
		return 0;
	}

	DTYPE sum = 0;
	/*sum.re = 0;
	sum.im = 0;

	cu_casum<<<num_block,max_thread>>>(mat->data, 1, mat_size-1, max_thread, &sum);
	CUDA_CALL(cudaThreadSynchronize())

	return sum.re + sum.im;*/


#if NTYPE == 1
	cublasDzasum(handle, mat_size, (cuDoubleComplex*)(mat->data), inc, &sum);
#elif NTYPE == 0
	cublasScasum(handle, mat_size, (cuComplex*)(mat->data), inc, &sum);
#endif
	cublasDestroy(handle);
	return sum;
}
/*
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
		//atomicAdd(&(*sum).re, s_cnt_re);
		//atomicAdd(&(*sum).im, s_cnt_im);
	}
}*/


/*** DOT ***/
DTYPE dot(MAT *src_x, UINT x_increment, MAT *src_y, UINT y_increment)
{
#if DEBUG
	printf("%s\n",__func__);
#endif
	UINT mat_size = src_x->d0 * src_x->d1 * src_x->d2;
	printf("max_thread : %d\n",max_thread);
	UINT num_block = (UINT)(mat_size/(UINT)max_thread)+1;

	DTYPE result = 0;

	MAT* temp = (MAT*)malloc(sizeof(MAT));
	temp->ndim = 2;
	temp->d0 = src_x->d0;
	temp->d1 = src_x->d1;
	temp->d2 = src_x->d2;

	cudaMalloc((void**)&(temp->data), sizeof(DTYPE) * src_x->d0 * src_x->d1 * src_x->d2);
	cudaMemset(temp->data, 0, sizeof(DTYPE)*src_x->d0*src_x->d1*src_x->d2);
	
	cu_dot<<<num_block,max_thread>>>(temp->data, src_x->data, 1, src_y->data, 1, mat_size-1, max_thread);
	CUDA_CALL(cudaThreadSynchronize())

	result = asum(temp, 1);

	cudaFree(temp->data);
	free(temp);

	return result;
}
__global__ void cu_dot(DTYPE *result, DTYPE *src_x, UINT x_inc, DTYPE *src_y, UINT y_inc, UINT len, UINT block_size)
{
	ITER idx = threadIdx.x + blockIdx.x * block_size;
	if(idx > len)
		return;
		
	result[idx] = src_x[idx * x_inc] * src_y[idx * y_inc];
}