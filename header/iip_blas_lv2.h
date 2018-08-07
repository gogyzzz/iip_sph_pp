#ifndef BLAS_LV2_H
#define BLAS_LV2_H
#include "iip_type.h"

/*
 * TODO
 *
 * gemv
 *
 * */

#ifndef USE_CUDA
void gemv(char transA, DTYPE alpha, MAT* A, MAT* X, DTYPE beta, MAT* Y);
void omp_gemv(char transA, UINT m, UINT n, DTYPE alpha, DTYPE* A, UINT lda,
             DTYPE* X, SINT incx, DTYPE beta, DTYPE* Y, SINT incy);

void cgemv(char transA, CTYPE alpha, CMAT* A, CMAT* X, CTYPE beta, CMAT* Y);
void omp_cgemv(char transA, UINT m, UINT n, CTYPE alpha, CTYPE* A, UINT lda,
              CTYPE* X, SINT incx, CTYPE beta, CTYPE* Y, SINT incy);

#else
void gemv(cublasOperation_t transA, DTYPE alpha, MAT* A, MAT* X, DTYPE beta,
          MAT* Y);
void cgemv(cublasOperation_t transA, CTYPE alpha, CMAT* A, CMAT* X, CTYPE beta,
           CMAT* Y);

#endif

#endif
