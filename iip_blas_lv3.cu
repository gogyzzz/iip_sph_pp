#include "iip_blas_lv3.h"

/*
 **  cblas_?gemm(layout,transA,transB,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc)
 **
 **    layout :   i) --->CblasRowMajor
 **    			   [0][1]  =  {0,1,2,3}
 **                	   [2][3]
 **
 **             ii)  |  [0][2] = {0,1,2,3}
 **                  |  [1][3]
 **                 \_/ CblasColMajor
 **
 **   
 **   C := alpha * op(A)*op(B) + beta*C
 **
 **     op(X) =  	  i) X      when transX = CblasNoTrans
 **
 **     		 	 ii) X**T     ''        = CblasTrans
 **
 **     			iii) X**H     ''        = CblasConjTrans
 **
 **      m  = the number of rows of op(A)
 **      n  = the number of columns of op(B) and C 
 **      k  = the number of columns of op(A) and rows of op(B)
 **
 **
 **      lda : the first dimension of A
 **      ldb : 		''	     B
 **      ldc :		''	     C
 **
 ** 		-the first dimension : the number of columns when CblasRowMajor
 **		                 	   	''    rows   when CblasColMajor
*/

void gemm(cublasOperation_t transA,cublasOperation_t transB, DTYPE alpha, MAT*A,MAT*B, DTYPE beta, MAT*C)
{

UINT m,n,k;
UINT lda,ldb,ldc;

#if DEBUG
printf("%s\n",__func__);
#endif

if(transA == NoTran)
{
m = A->d0;
k = A->d1;

lda = m;
ldc = m;
}
else
{
m = A->d1;
k = A->d0;

lda = k;
ldc = m;
}

if(transB == NoTran)
{
n = B->d1;

ldb = B->d0;
}
else
{
n = B->d0;

ldb = B ->d0;
}



	#if NTYPE == 0
	cublasSgemm(handle,transA,transB,m,n,k,&alpha,A->data,lda,B->data,ldb,&beta,C->data,ldc);
	#else
	cublasDgemm(handle,transA,transB,m,n,k,&alpha,A->data,lda,B->data,ldb,&beta,C->data,ldc);
	#endif
}


void cgemm(cublasOperation_t transA,cublasOperation_t transB, CTYPE alpha, CMAT*A,CMAT*B, CTYPE beta, CMAT*C)
{

UINT m,n,k;
UINT lda,ldb,ldc;

#if DEBUG
printf("%s\n",__func__);
#endif

if(transA == NoTran)
{
m = A->d0;
k = A->d1;

lda = m;
ldc = m;
}
else
{
m = A->d1;
k = A->d0;

lda = k;
ldc = m;
}

if(transB == NoTran)
{
n = B->d1;

ldb = B->d0;
}
else
{
n = B->d0;

ldb = B ->d0;
}


	#if NTYPE == 0
	cublasCgemm(handle,transA,transB,m,n,k,CU_CX(&alpha),  CU_CX(A->data),lda,CU_CX(B->data),ldb,CU_CX(&beta),CU_CX(C->data),ldc);
	#else
	cublasZgemm(handle,transA,transB,m,n,k,CU_CX(&alpha),CU_CX(A->data),lda,CU_CX(B->data),ldb,CU_CX(&beta),CU_CX(C->data),ldc);
	#endif
}
