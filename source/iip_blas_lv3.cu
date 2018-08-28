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

void gemm_mat(cublasOperation_t transA,cublasOperation_t transB, DTYPE alpha, MAT*A,MAT*B, DTYPE beta, MAT*C)
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

if((transA == CTran)||(transB == CTran ) )
{
	printf("ERROR : can't conjugate transpose real number matirx\n");
	return;
}


	#if NTYPE == 0
	cublasSgemm(handle,transA,transB,m,n,k,&alpha,A->data,lda,B->data,ldb,&beta,C->data,ldc);
	#else
	cublasDgemm(handle,transA,transB,m,n,k,&alpha,A->data,lda,B->data,ldb,&beta,C->data,ldc);
	#endif
}


void gemm_cmat(cublasOperation_t transA,cublasOperation_t transB, CTYPE alpha, CMAT*A,CMAT*B, CTYPE beta, CMAT*C)
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


/**** REAL ****/
void aABpbC(DTYPE alpha,MAT*A,MAT*B,DTYPE beta,MAT*C)
{
	gemm(NoTran,NoTran,alpha,A,B,beta,C);
}
void aABtpbC(DTYPE alpha,MAT*A,MAT*B,DTYPE beta,MAT*C)
{
	gemm(NoTran,Tran,alpha,A,B,beta,C);
}
void aAtBpbC(DTYPE alpha,MAT*A,MAT*B,DTYPE beta,MAT*C)
{
	gemm(Tran,NoTran,alpha,A,B,beta,C);
}
void aAtBtpbC(DTYPE alpha,MAT*A,MAT*B,DTYPE beta,MAT*C)
{
	gemm(Tran,Tran,alpha,A,B,beta,C);
}

void matmul(MAT*A,MAT*B,MAT*C)
{
	gemm(NoTran,NoTran,1,A,B,0,C);
}

/**** COMPLEX ****/
void caABpbC(CTYPE alpha,CMAT*A,CMAT*B,CTYPE beta,CMAT*C)
{
	cgemm(NoTran,NoTran,alpha,A,B,beta,C);
}
void caABtpbC(CTYPE alpha,CMAT*A,CMAT*B,CTYPE beta,CMAT*C)
{
	cgemm(NoTran,Tran,alpha,A,B,beta,C);
}
void caABjpbC(CTYPE alpha,CMAT*A,CMAT*B,CTYPE beta,CMAT*C)
{
	cgemm(NoTran,CTran,alpha,A,B,beta,C);
}


void caAtBpbC(CTYPE alpha,CMAT*A,CMAT*B,CTYPE beta,CMAT*C)
{
	cgemm(Tran,NoTran,alpha,A,B,beta,C);
}
void caAtBtpbC(CTYPE alpha,CMAT*A,CMAT*B,CTYPE beta,CMAT*C)
{
	cgemm(Tran,Tran,alpha,A,B,beta,C);
}
void caAtBhpbC(CTYPE alpha,CMAT*A,CMAT*B,CTYPE beta,CMAT*C)
{
	cgemm(Tran,CTran,alpha,A,B,beta,C);
}

void caAhBpbC(CTYPE alpha,CMAT*A,CMAT*B,CTYPE beta,CMAT*C)
{
	cgemm(CTran,NoTran,alpha,A,B,beta,C);
}
void caAhBtpbC(CTYPE alpha,CMAT*A,CMAT*B,CTYPE beta,CMAT*C)
{
	cgemm(CTran,Tran,alpha,A,B,beta,C);
}
void caAhBhpbC(CTYPE alpha,CMAT*A,CMAT*B,CTYPE beta,CMAT*C)
{
	cgemm(CTran,CTran,alpha,A,B,beta,C);
}

void cmatmul(CMAT*A,CMAT*B,CMAT*C)
{
CTYPE one_zero;
CTYPE zero_zero;
one_zero.re=1;
one_zero.im=0;
zero_zero.re=0;
zero_zero.im=0;
	cgemm(NoTran,NoTran,alpha,A,B,beta,C);
}
