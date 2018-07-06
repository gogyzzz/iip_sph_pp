//#include "mother.h"
#include "iip_matrix.h"

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
cublasHandle_t handle;
UINT max_thread;
UINT max_block;


/**** allocMAT  ****/
MAT* alloc_MAT_1d(UINT d0)
{
#if DEBUG
printf("%s\n",__func__);
#endif
	MAT* mat = (MAT*)malloc(sizeof(MAT));
	mat->ndim = 0;
	mat->d0 = d0;
	mat->d1 = 1;
	mat->d2 = 1;
	cudaMalloc((void**)&(mat->data),sizeof(DTYPE) * d0);

	return mat;
}
MAT* alloc_MAT_2d(UINT d0,UINT d1)
{
#if DEBUG
printf("%s\n",__func__);
#endif
	MAT* mat = (MAT*)malloc(sizeof(MAT));
	mat->ndim = 1;
	mat->d0 = d0;
	mat->d1 = d1;
	mat->d2 = 1;
	cudaMalloc((void**)&(mat->data),sizeof(DTYPE) * d0* d1);

	return mat;

}
MAT* alloc_MAT_3d(UINT d0,UINT d1,UINT d2)
{
#if DEBUG
printf("%s\n",__func__);
#endif
	MAT* mat = (MAT*)malloc(sizeof(MAT));
	mat->ndim = 2;
	mat->d0 = d0;
	mat->d1 = d1;
	mat->d2 = d2;
	cudaMalloc((void**)&(mat->data),sizeof(DTYPE) * d0* d1*d2);
	
	return mat;
}

CMAT* calloc_MAT_1d(UINT d0)
{
#if DEBUG
printf("%s\n",__func__);
#endif
	CMAT* mat = (CMAT*)malloc(sizeof(CMAT));
	mat->ndim = 0;
	mat->d0 = d0;
	mat->d1 = 1;
	mat->d2 = 1;
	cudaMalloc((void**)&(mat->data),sizeof(CTYPE) * d0);

	return mat;
}
CMAT* calloc_MAT_2d(UINT d0,UINT d1)
{
#if DEBUG
printf("%s\n",__func__);
#endif
	CMAT* mat = (CMAT*)malloc(sizeof(CMAT));
	mat->ndim = 1;
	mat->d0 = d0;
	mat->d1 = d1;
	mat->d2 = 1;
	cudaMalloc((void**)&(mat->data),sizeof(CTYPE) * d0* d1);

	return mat;
}
CMAT* calloc_MAT_3d(UINT d0,UINT d1,UINT d2)
{
#if DEBUG
printf("%s\n",__func__);
#endif
	CMAT* mat = (CMAT*)malloc(sizeof(CMAT));
	mat->ndim = 2;
	mat->d0 = d0;
	mat->d1 = d1;
	mat->d2 = d2;
	cudaMalloc((void**)&(mat->data),sizeof(CTYPE) * d0* d1* d2);

	return mat;
}
/**** zeros  ****/

MAT* zeros_1d(UINT d0)
{
#if DEBUG
printf("%s\n",__func__);
#endif
	MAT* mat = (MAT*)malloc(sizeof(MAT));
	mat->ndim = 0;
	mat->d0 = d0;
	mat->d1 = 1;
	mat->d2 = 1;
	cudaMalloc((void**)&(mat->data),sizeof(DTYPE) * d0);
    cudaMemset(mat->data,0,sizeof(DTYPE)*d0  );
	return mat;
}

MAT* zeros_2d(UINT d0,UINT d1)
{
#if DEBUG
printf("%s\n",__func__);
#endif
	MAT* mat = (MAT*)malloc(sizeof(MAT));
	mat->ndim = 1;
	mat->d0 = d0;
	mat->d1 = d1;
	mat->d2 = 1;
	
	cudaMalloc((void**)&(mat->data),sizeof(DTYPE) * d0 *d1);
    cudaMemset(mat->data,0,sizeof(DTYPE)*d0*d1);

	return mat;
}
MAT* zeros_3d(UINT d0,UINT d1,UINT d2)
{
#if DEBUG
printf("%s\n",__func__);
#endif
	MAT* mat = (MAT*)malloc(sizeof(MAT));
	mat->ndim = 2;
	mat->d0 = d0;
	mat->d1 = d1;
	mat->d2 = d2;

	cudaMalloc((void**)&(mat->data),sizeof(DTYPE) * d0 *d1 *d2);
    cudaMemset(mat->data,0,sizeof(DTYPE)*d0*d1*d2);
	return mat;
}


CMAT* czeros_1d(UINT d0)
{
#if DEBUG
printf("%s\n",__func__);
#endif
	CMAT* mat = (CMAT*)malloc(sizeof(CMAT));
	mat->ndim = 0;
	mat->d0 = d0;
	mat->d1 = 1;
	mat->d2 = 1;
	
	cudaMalloc((void**)&(mat->data),sizeof(CTYPE) * d0);
    cudaMemset(mat->data,0,sizeof(CTYPE)*d0);
	return mat;
}

CMAT* czeros_2d(UINT d0,UINT d1)
{
#if DEBUG
printf("%s\n",__func__);
#endif
	CMAT* mat = (CMAT*)malloc(sizeof(CMAT));
	mat->ndim = 1;
	mat->d0 = d0;
	mat->d1 = d1;
	mat->d2 = 1;
	
	cudaMalloc((void**)&(mat->data),sizeof(CTYPE) * d0 *d1);
    cudaMemset(mat->data,0,sizeof(CTYPE)*d0*d1);
	return mat;
}
CMAT* czeros_3d(UINT d0,UINT d1,UINT d2)
{
#if DEBUG
printf("%s\n",__func__);
#endif
	CMAT* mat = (CMAT*)malloc(sizeof(CMAT));
	mat->ndim = 2;
	mat->d0 = d0;
	mat->d1 = d1;
	mat->d2 = d2;
	
	cudaMalloc((void**)&(mat->data),sizeof(CTYPE) * d0 *d1*d2);
    cudaMemset(mat->data,0,sizeof(CTYPE)*d0*d1*d2);
	return mat;
}

/**** set ****/

void set_1d(MAT*mat, UINT idx0, DTYPE val)
{
#if DEBUG
printf("%s\n",__func__);
#endif
	cu_set<<<1,1>>>(mat->data,idx0,val);
	cudaThreadSynchronize();
}
void set_2d(MAT*mat, UINT idx0, UINT idx1, DTYPE val)
{
#if DEBUG
printf("%s\n",__func__);
#endif
	cu_set<<<1,1>>>(mat->data,idx0 + (mat->d0)*idx1,val);
	cudaThreadSynchronize();
}
void set_3d(MAT*mat, UINT idx0, UINT idx1, UINT idx2, DTYPE val)
{
#if DEBUG
printf("%s\n",__func__);
#endif
	cu_set<<<1,1>>>(mat->data,idx0 + (mat->d0)*idx1 + (mat->d0)*(mat->d1)*idx2,val);
	cudaThreadSynchronize();
}

__global__ void cu_set(DTYPE*data,UINT idx,DTYPE val)
{
	data[idx]=val;
}

void cset_1d(CMAT*mat, UINT idx0, DTYPE re,DTYPE im)
{
#if DEBUG
printf("%s\n",__func__);
#endif
	cu_cset<<<1,1>>>(mat->data,idx0,re,im);
	cudaThreadSynchronize();
}
void cset_2d(CMAT*mat, UINT idx0, UINT idx1, DTYPE re,DTYPE im)
{
#if DEBUG
printf("%s\n",__func__);
#endif
	cu_cset<<<1,1>>>(mat->data,idx0 + (mat->d0)*idx1,re,im);
	cudaThreadSynchronize();
}
void cset_3d(CMAT*mat, UINT idx0, UINT idx1, UINT idx2, DTYPE re,DTYPE im)
{
#if DEBUG
printf("%s\n",__func__);
#endif
	cu_cset<<<1,1>>>(mat->data,idx0 + (mat->d0)*idx1 + (mat->d0)*(mat->d1)*idx2,re,im);
	cudaThreadSynchronize();
}

__global__ void cu_cset(CTYPE*data,UINT idx,DTYPE re,DTYPE im)
{
	data[idx].re=re;
	data[idx].im=im;
}


/**** fill ****/

void fill(MAT*mat, DTYPE val) // for real mat
{
#if DEBUG
  printf("%s\n",__func__);
	printf("max thread : %d\n",max_thread);
#endif
  UINT len = (mat->d0) * (mat->d1) * (mat->d2);
	UINT num_block = (UINT)(len/(UINT)max_thread)+1;
	cu_fill<<<num_block,max_thread>>>(mat->data,len-1,val,max_thread);
	cudaThreadSynchronize();
}

__global__ void cu_fill(DTYPE* data, UINT len,DTYPE val,UINT size_block)
{
	ITER idx = threadIdx.x + blockIdx.x * size_block;
	if(idx > len)
		return;
  data[idx]= val;
}

void cfill(CMAT*mat, DTYPE re, DTYPE im) // for complex mat
{
#if DEBUG
  printf("%s\n",__func__);
#endif
  UINT len = (mat->d0) * (mat->d1) * (mat->d2);

	UINT num_block = (len/max_thread)+1;
	cu_cfill<<<num_block,max_thread>>>(mat->data,len-1,re,im,max_thread);
	cudaThreadSynchronize();
}
__global__ void cu_cfill(CTYPE* data, UINT len,DTYPE re,DTYPE im,UINT size_block)
{
	ITER idx = threadIdx.x + blockIdx.x * size_block;
	if(idx > len)
		return;
  data[idx].re = re;
	data[idx].im = im; 
}

/**** get ****/

DTYPE get_1d(MAT*mat, UINT idx0)
{
#if DEBUG
printf("%s\n",__func__);
#endif
  DTYPE ret;
	cudaMemcpy(&ret,&(mat->data[idx0]),sizeof(DTYPE),cudaMemcpyDeviceToHost);
	return ret;
}

DTYPE get_2d(MAT*mat, UINT idx0, UINT idx1)
{
#if DEBUG
printf("%s\n",__func__);
#endif
  DTYPE ret;
	cudaMemcpy(&ret,&(mat->data[idx0 + (mat->d0)*idx1]),sizeof(DTYPE),cudaMemcpyDeviceToHost);
	return ret;
}
DTYPE get_3d(MAT*mat, UINT idx0, UINT idx1, UINT idx2)
{
#if DEBUG
printf("%s\n",__func__);
#endif
  DTYPE ret;
	cudaMemcpy(&ret,&(mat->data[idx0 + (mat->d0)*idx1 + (mat->d0)*(mat->d1)*idx2]),sizeof(DTYPE),cudaMemcpyDeviceToHost);
	return ret;
}

CTYPE cget_1d(CMAT*mat, UINT idx0)
{
#if DEBUG
printf("%s\n",__func__);
#endif
  CTYPE ret;
	cudaMemcpy(&ret,&(mat->data[idx0]),sizeof(CTYPE),cudaMemcpyDeviceToHost);
	return ret;
}

CTYPE cget_2d(MAT*mat, UINT idx0, UINT idx1)
{
#if DEBUG
printf("%s\n",__func__);
#endif
  CTYPE ret;
	cudaMemcpy(&ret,&(mat->data[idx0 + (mat->d0)*idx1]),sizeof(CTYPE),cudaMemcpyDeviceToHost);
	return ret;
}
CTYPE cget_3d(MAT*mat, UINT idx0, UINT idx1, UINT idx2)
{
#if DEBUG
printf("%s\n",__func__);
#endif
  CTYPE ret;
	cudaMemcpy(&ret,&(mat->data[idx0 + (mat->d0)*idx1 + (mat->d0)*(mat->d1)*idx2]),sizeof(CTYPE),cudaMemcpyDeviceToHost);
	return ret;
}
/**** submat ****/

void submat_1d(MAT* mat, MAT* submat, ITER d0_st, ITER d0_ed)
{
	UINT num_thread;
#if DEBUG
  printf("%s\n",__func__);
#endif
   if (mat->d0 < submat->d0) {printf("error in \nmat->d0 < submat->d0\n");return;}
  if ( d0_st == -1 ) d0_st = 0;
  if ( d0_ed == -1 ) d0_ed = mat->d0;
   
  num_thread = (d0_ed - d0_st +1) > max_thread ? max_thread : (d0_ed - d0_st +1);
 cu_submat<<<1,num_thread>>>(mat->data,submat->data,d0_st,d0_ed,0,0,mat->d0,1,submat->d0,1,num_thread); 
	cudaThreadSynchronize();
}


void submat_2d(MAT* mat, MAT* submat,
    ITER d0_st, ITER d0_ed,
    ITER d1_st, ITER d1_ed)
{
	UINT num_thread;
#if DEBUG
  printf("%s\n",__func__);
#endif
    if (mat->d0 < submat->d0 || mat->d1 < submat->d1) {printf("error in \nmat->d0 < submat->d0 || mat->d1 < submat->d1 \n");return;}
  if ( d0_st == -1 ) d0_st = 0;
  if ( d0_ed == -1 ) d0_ed = mat->d0;
  if ( d1_st == -1 ) d1_st = 0;
  if ( d1_ed == -1 ) d1_ed = mat->d1;
   
  num_thread = (d0_ed - d0_st +1) > max_thread ? max_thread : (d0_ed - d0_st +1);

 cu_submat<<<d1_ed - d1_st +1,num_thread>>>(mat->data,submat->data,d0_st,d0_ed,d1_st,0,mat->d0,mat->d1,submat->d0,submat->d1,num_thread); 
	cudaThreadSynchronize();
}


void submat_3d(MAT* mat, MAT* submat, 
    ITER d0_st, ITER d0_ed,
    ITER d1_st, ITER d1_ed,
    ITER d2_st, ITER d2_ed)
{
	UINT num_thread;
	dim3 Dg;
#if DEBUG
  printf("%s\n",__func__);
#endif
  if (mat->d0 < submat->d0 || mat->d1 < submat->d1 || mat->d2 < submat->d2) {printf("error in \nmat->d0 < submat->d0 || mat->d1 < submat->d1 || mat->d2 < submat->d2\n");return;}
  if ( d0_st == -1 ) d0_st = 0;
  if ( d0_ed == -1 ) d0_ed = mat->d0;
  if ( d1_st == -1 ) d1_st = 0;
  if ( d1_ed == -1 ) d1_ed = mat->d1;
  if ( d2_st == -1 ) d2_st = 0;
  if ( d2_ed == -1 ) d2_ed = mat->d2;
   
	Dg = dim3(d1_ed - d1_st +1,d2_ed - d2_st +1,1);
  num_thread = (d0_ed - d0_st +1) > max_thread ? max_thread : (d0_ed - d0_st +1);
 cu_submat<<<Dg,num_thread>>>(mat->data,submat->data,d0_st,d0_ed,d1_st,d2_st,mat->d0,mat->d1,submat->d0,submat->d1,num_thread); 
	cudaThreadSynchronize();
}

__global__ void cu_submat(DTYPE*data,DTYPE*sub_data,
													ITER d0_st, ITER d0_ed,
													ITER d1_st,	ITER d2_st,
													UINT mat_d0, UINT mat_d1,
													UINT sub_d0, UINT sub_d1,
													UINT num_thread
)
{
	ITER i,j,k;
  
	k = blockIdx.y;
	j = blockIdx.x;
	for (i = threadIdx.x; d0_st+i < d0_ed; i+=num_thread)
			{ 				
        sub_data[i + j*(sub_d0) + k*(sub_d0*sub_d1)]
         = data[(d0_st+i) + (d1_st+j)*(mat_d0) + (d2_st+k)*(mat_d0*mat_d1)];
			}
}

void csubmat_1d(CMAT* mat, CMAT* csubmat, ITER d0_st, ITER d0_ed)
{
	UINT num_thread;
#if DEBUG
  printf("%s\n",__func__);
#endif
   if (mat->d0 < csubmat->d0) {printf("error in \nmat->d0 < csubmat->d0\n");return;}
  if ( d0_st == -1 ) d0_st = 0;
  if ( d0_ed == -1 ) d0_ed = mat->d0;
   
  num_thread = (d0_ed - d0_st +1) > max_thread ? max_thread : (d0_ed - d0_st +1);
 cu_csubmat<<<1,num_thread>>>(mat->data,csubmat->data,d0_st,d0_ed,0,0,mat->d0,1,csubmat->d0,1,num_thread); 
	cudaThreadSynchronize();
}


void csubmat_2d(CMAT* mat, CMAT* csubmat,
    ITER d0_st, ITER d0_ed,
    ITER d1_st, ITER d1_ed)
{
	UINT num_thread;
#if DEBUG
  printf("%s\n",__func__);
#endif
    if (mat->d0 < csubmat->d0 || mat->d1 < csubmat->d1) {printf("error in \nmat->d0 < csubmat->d0 || mat->d1 < csubmat->d1 \n");return;}
  if ( d0_st == -1 ) d0_st = 0;
  if ( d0_ed == -1 ) d0_ed = mat->d0;
  if ( d1_st == -1 ) d1_st = 0;
  if ( d1_ed == -1 ) d1_ed = mat->d1;
   
  num_thread = (d0_ed - d0_st +1) > max_thread ? max_thread : (d0_ed - d0_st +1);

 cu_csubmat<<<d1_ed - d1_st +1,num_thread>>>(mat->data,csubmat->data,d0_st,d0_ed,d1_st,0,mat->d0,mat->d1,csubmat->d0,csubmat->d1,num_thread); 
	cudaThreadSynchronize();
}


void csubmat_3d(CMAT* mat, CMAT* csubmat, 
    ITER d0_st, ITER d0_ed,
    ITER d1_st, ITER d1_ed,
    ITER d2_st, ITER d2_ed)
{
	UINT num_thread;
	dim3 Dg;
#if DEBUG
  printf("%s\n",__func__);
#endif
  if (mat->d0 < csubmat->d0 || mat->d1 < csubmat->d1 || mat->d2 < csubmat->d2) {printf("error in \nmat->d0 < csubmat->d0 || mat->d1 < csubmat->d1 || mat->d2 < csubmat->d2\n");return;}
  if ( d0_st == -1 ) d0_st = 0;
  if ( d0_ed == -1 ) d0_ed = mat->d0;
  if ( d1_st == -1 ) d1_st = 0;
  if ( d1_ed == -1 ) d1_ed = mat->d1;
  if ( d2_st == -1 ) d2_st = 0;
  if ( d2_ed == -1 ) d2_ed = mat->d2;
   
	Dg = dim3(d1_ed - d1_st +1,d2_ed - d2_st +1,1);
  num_thread = (d0_ed - d0_st +1) > max_thread ? max_thread : (d0_ed - d0_st +1);
 cu_csubmat<<<Dg,num_thread>>>(mat->data,csubmat->data,d0_st,d0_ed,d1_st,d2_st,mat->d0,mat->d1,csubmat->d0,csubmat->d1,num_thread); 
	cudaThreadSynchronize();
}

__global__ void cu_csubmat(CTYPE*data,CTYPE*sub_data,
													ITER d0_st, ITER d0_ed,
													ITER d1_st,	ITER d2_st,
													UINT mat_d0, UINT mat_d1,
													UINT sub_d0, UINT sub_d1,
													UINT num_thread
)
{
	ITER i,j,k;
  
	k = blockIdx.y;
	j = blockIdx.x;
	for (i = threadIdx.x; d0_st+i < d0_ed; i+=num_thread)
			{ 				
        sub_data[i + j*(sub_d0) + k*(sub_d0*sub_d1)].re
         = data[(d0_st+i) + (d1_st+j)*(mat_d0) + (d2_st+k)*(mat_d0*mat_d1)].re;
        sub_data[i + j*(sub_d0) + k*(sub_d0*sub_d1)].im
         = data[(d0_st+i) + (d1_st+j)*(mat_d0) + (d2_st+k)*(mat_d0*mat_d1)].im;
			}
}


/**** miscellaneous  ****/
void free_MAT(MAT *mat)
{
#if DEBUG
printf("%s\n",__func__);
#endif
	cudaFree(mat->data);
	free(mat);
}

void free_CMAT(CMAT *mat)
{
#if DEBUG
printf("%s\n",__func__);
#endif
	cudaFree(mat->data);
	free(mat);
}


void print_MAT(MAT* mat)
{
#if DEBUG
printf("%s\n",__func__);
#endif

  cu_print_MAT<<<1,1>>>(mat->data,mat->d0,mat->d1,mat->d2);
	cudaDeviceSynchronize();
}

__global__ void cu_print_MAT(DTYPE* data,UINT d0,UINT d1, UINT d2)
{
	ITER i,j,k;
	
	for(k = 0; k < d2; k++)
	{
		for(i=0; i < d0 ; i++)	
		{
			for(j=0; j < d1;j++)
				printf("%.3lf ",data[k*(d1)*(d0) + j*(d0) + i]);
			printf("\n");
		}
	printf("\n");		
	}
}

void print_CMAT(CMAT* mat)
{
#if DEBUG
printf("%s\n",__func__);
#endif

	cu_print_CMAT<<<1,1>>>(mat->data,mat->d0,mat->d1,mat->d2);
  cudaDeviceSynchronize();
}

__global__ void cu_print_CMAT(CTYPE*data,UINT d0,UINT d1,UINT d2)
{
ITER i,j,k;
	
	for(k = 0; k < d2; k++)
	{
		for(i=0; i < d0 ; i++)	
		{
			for(j=0; j < d1;j++)
			{printf("%.3lf ",data[k*(d1)*(d0) + j*(d0) + i].re);
				printf("%.3lf|",data[k*(d1)*(d0) + j*(d0) + i].im);
			}printf("\n");
		}
	printf("\n");		
	}
}


