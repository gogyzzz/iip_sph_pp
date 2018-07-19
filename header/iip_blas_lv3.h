#ifndef BLAS_LV3_H
#define BLAS_LV3_H
#include "iip_type.h"

/**** REAL ****/

void aABpbC(DTYPE alpha, MAT*A,MAT*B,DTYPE beta, MAT*C);
void aABtpbC(DTYPE alpha, MAT*A,MAT*B,DTYPE beta, MAT*C);

void aAtBpbC(DTYPE alpha, MAT*A,MAT*B,DTYPE beta, MAT*C);
void aAtBtpbC(DTYPE alpha, MAT*A,MAT*B,DTYPE beta, MAT*C);

void matmul(MAT* A,MAT* B,MAT* C);

/**** COMPLEX ****/

void caABpbC(CTYPE alpha, CMAT*A,CMAT*B,CTYPE beta,CMAT*C);
void caABtpbC(CTYPE alpha, CMAT*A,CMAT*B,CTYPE beta,CMAT*C);
void caABhpbC(CTYPE alpha, CMAT*A,CMAT*B,CTYPE beta,CMAT*C);

void caAtBpbC(CTYPE alpha, CMAT*A,CMAT*B,CTYPE beta,CMAT*C);
void caAtBtpbC(CTYPE alpha, CMAT*A,CMAT*B,CTYPE beta,CMAT*C);
void caAtBhpbC(CTYPE alpha, CMAT*A,CMAT*B,CTYPE beta,CMAT*C);

void caAhBpbC(CTYPE alpha, CMAT*A,CMAT*B,CTYPE beta,CMAT*C);
void caAhBtpbC(CTYPE alpha, CMAT*A,CMAT*B,CTYPE beta,CMAT*C);
void caAhBhpbC(CTYPE alpha, CMAT*A,CMAT*B,CTYPE beta,CMAT*C);

void cmatmul(CMAT* A,CMAT* B,CMAT* C);

	#ifndef USE_CUDA
void gemm(char transA, char transB, DTYPE alpha, MAT* A,MAT* B,DTYPE beta,MAT* C);
void mp_gemm(char transA, char transB, UINT m,UINT n,UINT k,DTYPE alpha, DTYPE* A,UINT lda,DTYPE* B,UINT ldb,DTYPE beta,DTYPE* C,UINT ldc);

void cgemm(char transA, char transB, CTYPE alpha, CMAT* A,CMAT* B,CTYPE beta,CMAT* C);
void mp_cgemm(char transA, char transB, UINT m,UINT n,UINT k,CTYPE alpha, CTYPE* A,UINT lda,CTYPE* B,UINT ldb,CTYPE beta,CTYPE* C,UINT ldc);
#else

void gemm(cublasOperation_t transA, cublasOperation_t transB, DTYPE alpha, MAT* A,MAT* B,DTYPE beta,MAT* C);

void cgemm(cublasOperation_t transA,cublasOperation_t transB,  CTYPE alpha, CMAT* A,CMAT* B,CTYPE beta,CMAT* C);
	#endif




#endif
