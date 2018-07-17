#include "iip_blas_lv2.h"

/****  gemv ****/

void gemv(cublasOperation_t transA, DTYPE alpha, MAT* A, MAT*X, DTYPE beta, MAT*Y )
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


m = A->d1;
n = A->d0;
lda = m;

#if DEBUG
printf("trans : %d m : %u n: %u lda : %u\nalpha : %lf beta : %lf\n",transA,m,n,lda,alpha,beta);
#endif

	#if NTYPE==0
		cublasSgemv(handle,transA, m , n, &alpha,(A->data), lda, X->data,1,&beta,Y->data,1);
	#else
		cublasDgemv(handle,transA, m , n,  &alpha, (A->data), lda,(X->data),1,&beta,Y->data,1);
	#endif
}

void cgemv(cublasOperation_t transA, CTYPE alpha, CMAT* A, CMAT*X, CTYPE beta, CMAT*Y )
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


m = A->d1;
n = A->d0;
lda = m;

#if DEBUG
printf("trans : %d m : %u n: %u lda : %u\nalpha : %lf|%lf beta : %lf|%lf\n",transA,m,n,lda,alpha.re,alpha.im,beta.re,beta.im);
#endif

#if NTYPE==0
	cublasCgemv(handle,transA, m , n, CU_CX(&alpha), CU_CX(A->data), lda, CU_CX(X->data),1,CU_CX(&beta),CU_CX(Y->data),1);
#else
	cublasZgemv(handle,transA, m , n,  CU_CX(&alpha), CU_CX(A->data), lda,CU_CX(X->data),1,CU_CX(&beta),CU_CX(Y->data),1);
		
#endif
//		mp_cgemv(transA, m , n, alpha, A->data, lda, X->data,1,beta,Y->data,1);
}

