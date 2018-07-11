#include "iip_blas_lv2.h"

/****  dbmv ****/

void gemv(int32_t transA, DTYPE alpha, MAT* A, MAT*X, DTYPE beta, MAT*Y )
{
UINT m,n,lda;
#if DEBUG
printf("%s\n",__func__);
#endif

if(X->ndim != 0 || Y->ndim != 0 )
	{
		printf("A : %u  X : %u Y : %u\n",A->ndim,X->ndim,Y->ndim);
		printf("Use Vector for Vector operation\n");
		return;
	}
if(A->ndim != 1)
	{
		printf("Use 2D-Matrix for BLAS operation\n");
		return;
	}

if(transA == NoTran)
{
	m = A->d1;
	n = A->d0;
}
else
{
	m = A->d0;
	n = A->d1;
}
lda = A->d1;

#if DEBUG
printf("m : %u n: %u lda : %u\nalpha : %lf beta : %lf\n",m,n,lda,alpha,beta);
#endif

#if USE_CBLAS
	#if NTYPE==0
		cblas_sgemv(CblasColMajor,transA, m , n, alpha, A->data, lda, X->data,1,beta,Y->data,1);
	#else
		cblas_dgemv(CblasColMajor,transA, m , n,  alpha, A->data, lda, X->data,1,beta,Y->data,1);
		
	#endif
#else
		mp_gemv(transA, m , n, alpha, A->data, lda, X->data,1,beta,Y->data,1);
#endif
}

void mp_gemv(int32_t tranA, UINT m, UINT n, DTYPE alpha, DTYPE*A, UINT lda, DTYPE* X, SINT incx, DTYPE beta,DTYPE*Y, SINT incy)
{
	ITER i,j;
	DTYPE temp;

#if DEBUG
printf("%s\n",__func__);
#endif

	if(tranA==Tran)
	{
#pragma omp parallel for shared(A,X,Y) private(temp,j,i)					
	    for(j=0;j<lda;j++)
			{
				temp=0;
				for(i=0;i<m;i++)
				{
						temp += A[j + i*n]*X[i*incx];		
				}
				temp *=alpha;
				temp +=beta*Y[j*incy];
				Y[j*incy] = temp;
			}

	}
	else if(tranA==NoTran)
	{
#pragma omp parallel for shared(A,X,Y) private(temp,j,i)					
		for(j=0;j<lda;j++)
		{
			temp = 0;
			for(i=0;i<n;i++)
			{
				temp += A[i + n*j]*X[i*incx];
			}
				temp*=alpha;
				temp+= beta*Y[j*incy];
			  Y[j*incy] = temp;
		}	
	}
	else
	{
		printf("ERROR : Transpose argument is invalid\n");
		return;
	}


}
