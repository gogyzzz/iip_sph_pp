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
#ifndef BLAS_LV3_H
#define BLAS_LV3_H
#include "iip_type.h"

/**** REAL ****/
/* C = alpha * A * B + beta * C   and its derivations
 * t is transposed
 * h is hermitian transposed
 * */
void aABpbC(DTYPE alpha, MAT* A, MAT* B, DTYPE beta, MAT* C);
void aABtpbC(DTYPE alpha, MAT* A, MAT* B, DTYPE beta, MAT* C);

void aAtBpbC(DTYPE alpha, MAT* A, MAT* B, DTYPE beta, MAT* C);
void aAtBtpbC(DTYPE alpha, MAT* A, MAT* B, DTYPE beta, MAT* C);

/* C = A*B   */
void matmul(MAT* A, MAT* B, MAT* C);

/**** COMPLEX ****/

void caABpbC(CTYPE alpha, CMAT* A, CMAT* B, CTYPE beta, CMAT* C);
void caABtpbC(CTYPE alpha, CMAT* A, CMAT* B, CTYPE beta, CMAT* C);
void caABhpbC(CTYPE alpha, CMAT* A, CMAT* B, CTYPE beta, CMAT* C);

void caAtBpbC(CTYPE alpha, CMAT* A, CMAT* B, CTYPE beta, CMAT* C);
void caAtBtpbC(CTYPE alpha, CMAT* A, CMAT* B, CTYPE beta, CMAT* C);
void caAtBhpbC(CTYPE alpha, CMAT* A, CMAT* B, CTYPE beta, CMAT* C);

void caAhBpbC(CTYPE alpha, CMAT* A, CMAT* B, CTYPE beta, CMAT* C);
void caAhBtpbC(CTYPE alpha, CMAT* A, CMAT* B, CTYPE beta, CMAT* C);
void caAhBhpbC(CTYPE alpha, CMAT* A, CMAT* B, CTYPE beta, CMAT* C);

void cmatmul(CMAT* A, CMAT* B, CMAT* C);

#ifndef USE_CUDA
/* C = alhpa * op(A) * op(B) + beta * C
 * trans options :
 * NoTran : No transposed
 * Tran   : Trasposed
 * CTran  : Hermitian transposed
 *
 * ex) gemm_mat(NoTran,NoTran,1,A,B,0,C)
 *     is equal to matmul(A,B,C)
 *
 *     gemm_cmat(CTran,Tran,alpha,A,B,beta,C)
 *     is equal to caAhBtpbC(alpha,A,B,beta,C)
 *
 *
 * */
void gemm_mat(char transA, char transB, DTYPE alpha, MAT* A, MAT* B, DTYPE beta,
          MAT* C);
void omp_gemm(char transA, char transB, UINT m, UINT n, UINT k, DTYPE alpha,
              DTYPE* A, UINT lda, DTYPE* B, UINT ldb, DTYPE beta, DTYPE* C,
              UINT ldc);

void gemm_cmat(char transA, char transB, CTYPE alpha, CMAT* A, CMAT* B, CTYPE beta,
           CMAT* C);
void omp_cgemm(char transA, char transB, UINT m, UINT n, UINT k, CTYPE alpha,
               CTYPE* A, UINT lda, CTYPE* B, UINT ldb, CTYPE beta, CTYPE* C,
               UINT ldc);
#else

void gemm_mat(cublasOperation_t transA, cublasOperation_t transB, DTYPE alpha,
          MAT* A, MAT* B, DTYPE beta, MAT* C);

void gemm_cmat(cublasOperation_t transA, cublasOperation_t transB, CTYPE alpha,
           CMAT* A, CMAT* B, CTYPE beta, CMAT* C);
#endif

#endif
