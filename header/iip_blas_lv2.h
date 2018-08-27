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
#ifndef BLAS_LV2_H
#define BLAS_LV2_H
#include "iip_type.h"

/*  matrix * vector operation
 *  Since we don't have vector type
 *  we recommed to use gemm of blas_lv3.
 *
 *
 *  */

#ifndef USE_CUDA
void gemv_mat(char transA, DTYPE alpha, MAT* A, MAT* X, DTYPE beta, MAT* Y);
void omp_gemv(char transA, UINT m, UINT n, DTYPE alpha, DTYPE* A, UINT lda,
              DTYPE* X, SINT incx, DTYPE beta, DTYPE* Y, SINT incy);

void gemv_cmat(char transA, CTYPE alpha, CMAT* A, CMAT* X, CTYPE beta, CMAT* Y);
void omp_cgemv(char transA, UINT m, UINT n, CTYPE alpha, CTYPE* A, UINT lda,
               CTYPE* X, SINT incx, CTYPE beta, CTYPE* Y, SINT incy);

#else
void gemv_mat(cublasOperation_t transA, DTYPE alpha, MAT* A, MAT* X, DTYPE beta,
          MAT* Y);
void gemv_cmat(cublasOperation_t transA, CTYPE alpha, CMAT* A, CMAT* X, CTYPE beta,
           CMAT* Y);

#endif

#endif
