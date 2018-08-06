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

void gemm(char transA, char transB, DTYPE alpha, MAT* A, MAT* B, DTYPE beta,
          MAT* C) {
  UINT m, n, k;
  UINT lda, ldb, ldc;
  ITER i;
#if DEBUG
  printf("%s\n", __func__);
#endif

  if (transA == NoTran) {
    m = A->d0;
    k = A->d1;

    lda = m;
    ldc = m;
  } else {
    m = A->d1;
    k = A->d0;

    lda = k;
    ldc = m;
  }

  if (transB == NoTran) {
    n = B->d1;

    ldb = B->d0;
  } else {
    n = B->d0;

    ldb = B->d0;
  }

  if ((transA == CTran) || (transB == CTran)) {
    printf("ERROR : can't conjugate transpose real number matrix\n");
    return;
  }

  /** BATCH OPERATION **/
  if (A->d2 == B->d2) {
    for (i = 0; i < A->d2; i++) {
#if USE_CBLAS
#if NTYPE == 0
      cblas_sgemm(CblasColMajor, transA, transB, m, n, k, alpha,
                  &(A->data[i * (A->d0 * A->d1)]), lda,
                  &(B->data[i * B->d0 * B->d1]), ldb, beta,
                  &(C->data[i * C->d0 * C->d1]), ldc);
#else
      cblas_dgemm(CblasColMajor, transA, transB, m, n, k, alpha,
                  &(A->data[i * A->d0 * A->d1]), lda,
                  &(B->data[i * B->d0 * B->d1]), ldb, beta,
                  &(C->data[i * C->d0 * C->d1]), ldc);
#endif
#else
      mp_gemm(transA, transB, m, n, k, alpha, &(A->data[i * A->d0 * A->d1]),
              lda, &(B->data[i * B->d0 * B->d1]), ldb, beta,
              &(C->data[i * C->d0 * C->d1]), ldc);
#endif
    }
  } else if (A->d2 == 1 && B->d2 != 1) {
    if (C->d2 != B->d2) ASSERT(DIM_INVAL)
    for (i = 0; i < B->d2; i++) {
#if USE_CBLAS
#if NTYPE == 0
      cblas_sgemm(CblasColMajor, transA, transB, m, n, k, alpha, A->data, lda,
                  &(B->data[i * B->d0 * B->d1]), ldb, beta,
                  &(C->data[i * C->d0 * C->d1]), ldc);
#else
      cblas_dgemm(CblasColMajor, transA, transB, m, n, k, alpha, A->data, lda,
                  &(B->data[i * B->d0 * B->d1]), ldb, beta,
                  &(C->data[i * C->d0 * C->d1]), ldc);
#endif
#else
      mp_gemm(transA, transB, m, n, k, alpha, A->data, lda,
              &(B->data[i * B->d0 * B->d1]), ldb, beta,
              &(C->data[i * C->d0 * C->d1]), ldc);
#endif
    }
  } else if (A->d2 != 1 && B->d2 == 1) {
    if (C->d2 != A->d2) ASSERT(DIM_INVAL)
    for (i = 0; i < A->d2; i++) {
#if USE_CBLAS
#if NTYPE == 0
      cblas_sgemm(CblasColMajor, transA, transB, m, n, k, alpha,
                  &(A->data[i * A->d0 * A->d1]), lda, B->data, ldb, beta,
                  &(C->data[i * C->d0 * C->d1]), ldc);
#else
      cblas_dgemm(CblasColMajor, transA, transB, m, n, k, alpha,
                  &(A->data[i * A->d0 * A->d1]), lda, B->data, ldb, beta,
                  &(C->data[i * C->d0 * C->d1]), ldc);
#endif
#else
      mp_gemm(transA, transB, m, n, k, alpha, &(A->data[i * A->d0 * A->d1]),
              lda, B->data, ldb, beta, &(C->data[i * C->d0 * C->d1]), ldc);
#endif
    }
  } else
    ASSERT(DIM_INVAL)
}

void mp_gemm(char transA, char transB, UINT m, UINT n, UINT k, DTYPE alpha,
             DTYPE* A, UINT lda, DTYPE* B, UINT ldb, DTYPE beta, DTYPE* C,
             UINT ldc) {
  ITER i, j, l;
  DTYPE temp;
#if DEBUG
  printf("%s\n", __func__);
#endif

  if ((transA == NoTran) && (transB == NoTran)) {
#pragma omp parallel for shared(A, B, C) private(temp, i, j, l)
    for (l = 0; l < m; l++) {
      for (j = 0; j < n; j++) {
        temp = 0;
        for (i = 0; i < k; i++) {
          temp += A[i * m + l] * B[i + j * k];
        }
        C[l + m * j] *= beta;

        temp *= alpha;
        C[l + m * j] += temp;
      }
    }
  }

  if ((transA != NoTran) && (transB == NoTran)) {
#pragma omp parallel for shared(A, B, C) private(temp, i, j, l)
    for (l = 0; l < m; l++) {
      for (j = 0; j < n; j++) {
        temp = 0;
        for (i = 0; i < k; i++) {
          temp += A[l * k + i] * B[i + j * k];
        }
        C[l + m * j] *= beta;

        temp *= alpha;
        C[l + m * j] += temp;
      }
    }
  }
  if ((transA != NoTran) && (transB != NoTran)) {
#pragma omp parallel for shared(A, B, C) private(temp, i, j, l)
    for (l = 0; l < m; l++) {
      for (j = 0; j < n; j++) {
        temp = 0;
        for (i = 0; i < k; i++) {
          temp += A[l * k + i] * B[i * n + j];
        }
        C[l + m * j] *= beta;

        temp *= alpha;
        C[l + m * j] += temp;
      }
    }
  }
  if ((transA == NoTran) && (transB != NoTran)) {
#pragma omp parallel for shared(A, B, C) private(temp, i, j, l)
    for (l = 0; l < m; l++) {
      for (j = 0; j < n; j++) {
        temp = 0;
        for (i = 0; i < k; i++) {
          temp += A[i * m + l] * B[i * n + j];
        }
        C[l + m * j] *= beta;

        temp *= alpha;
        C[l + m * j] += temp;
      }
    }
  }
}

void cgemm(char transA, char transB, CTYPE alpha, CMAT* A, CMAT* B, CTYPE beta,
           CMAT* C) {
  UINT m, n, k;
  UINT lda, ldb, ldc;
  ITER i;
#if DEBUG
  printf("%s\n", __func__);
#endif

  if (transA == NoTran) {
    m = A->d0;
    k = A->d1;

    lda = m;
    ldc = m;
  } else {
    m = A->d1;
    k = A->d0;

    lda = k;
    ldc = m;
  }

  if (transB == NoTran) {
    n = B->d1;

    ldb = B->d0;
  } else {
    n = B->d0;

    ldb = B->d0;
  }

  if (A->d2 == B->d2) {
    for (i = 0; i < A->d2; i++) {
#if USE_CBLAS
#if NTYPE == 0
      cblas_cgemm(CblasColMajor, transA, transB, m, n, k, &alpha,
                  &(A->data[i * (A->d0 * A->d1)]), lda,
                  &(B->data[i * B->d0 * B->d1]), ldb, &beta,
                  &(C->data[i * C->d0 * C->d1]), ldc);
#else
      cblas_zgemm(CblasColMajor, transA, transB, m, n, k, &alpha,
                  &(A->data[i * A->d0 * A->d1]), lda,
                  &(B->data[i * B->d0 * B->d1]), ldb, &beta,
                  &(C->data[i * C->d0 * C->d1]), ldc);
#endif
#else
      mp_cgemm(transA, transB, m, n, k, alpha, &(A->data[i * A->d0 * A->d1]),
               lda, &(B->data[i * B->d0 * B->d1]), ldb, beta,
               &(C->data[i * C->d0 * C->d1]), ldc);
#endif
    }
  } else if (A->d2 == 1 && B->d2 != 1) {
    for (i = 0; i < B->d2; i++) {
      if (C->d2 != B->d2) ASSERT(DIM_INVAL)
#if USE_CBLAS
#if NTYPE == 0
      cblas_cgemm(CblasColMajor, transA, transB, m, n, k, &alpha, A->data, lda,
                  &(B->data[i * B->d0 * B->d1]), ldb, &beta,
                  &(C->data[i * C->d0 * C->d1]), ldc);
#else
      cblas_zgemm(CblasColMajor, transA, transB, m, n, k, &alpha, A->data, lda,
                  &(B->data[i * B->d0 * B->d1]), ldb, &beta,
                  &(C->data[i * C->d0 * C->d1]), ldc);
#endif
#else
      mp_cgemm(transA, transB, m, n, k, alpha, A->data, lda,
               &(B->data[i * B->d0 * B->d1]), ldb, beta,
               &(C->data[i * C->d0 * C->d1]), ldc);
#endif
    }
  } else if (A->d2 != 1 && B->d2 == 1) {
    for (i = 0; i < A->d2; i++) {
      if (C->d2 != A->d2) ASSERT(DIM_INVAL)
#if USE_CBLAS
#if NTYPE == 0
      cblas_cgemm(CblasColMajor, transA, transB, m, n, k, &alpha,
                  &(A->data[i * A->d0 * A->d1]), lda, B->data, ldb, &beta,
                  &(C->data[i * C->d0 * C->d1]), ldc);
#else
      cblas_zgemm(CblasColMajor, transA, transB, m, n, k, &alpha,
                  &(A->data[i * A->d0 * A->d1]), lda, B->data, ldb, &beta,
                  &(C->data[i * C->d0 * C->d1]), ldc);
#endif
#else
      mp_cgemm(transA, transB, m, n, k, alpha, &(A->data[i * A->d0 * A->d1]),
               lda, B->data, ldb, beta, &(C->data[i * C->d0 * C->d1]), ldc);
#endif
    }
  } else
    ASSERT(DIM_INVAL)
}

void mp_cgemm(char transA, char transB, UINT m, UINT n, UINT k, CTYPE alpha,
              CTYPE* A, UINT lda, CTYPE* B, UINT ldb, CTYPE beta, CTYPE* C,
              UINT ldc) {
  ITER i, j, l;
  ITER adx, bdx;
  CTYPE temp;
  DTYPE temp2;
  DTYPE o1, o2;
  DTYPE t_r, t_i;
#if DEBUG
  printf("%s\n", __func__);
#endif

  if ((transA == NoTran)) {
    if (transB == NoTran) {
#pragma omp parallel for shared(A, B, C) private(temp, i, j, l, temp2)
      for (l = 0; l < m; l++) {
        for (j = 0; j < n; j++) {
          temp.im = 0;
          temp.re = 0;
          for (i = 0; i < k; i++) {
            CXADD_mul(temp, A[i * m + l], B[i + j * k]);
          }
          CXMUL(C[l + m * j], beta, temp2);
          CXMUL(temp, alpha, temp2);
          CXADD(C[l + m * j], temp);
        }
      }

    } else if (transB == Tran) {
#pragma omp parallel for shared(A, B, C) private(temp, i, j, l, temp2)
      for (l = 0; l < m; l++) {
        for (j = 0; j < n; j++) {
          temp.re = 0;
          temp.im = 0;
          for (i = 0; i < k; i++) {
            CXADD_mul(temp, A[i * m + l], B[i * n + j])
          }
          CXMUL(C[l + m * j], beta, temp2) CXMUL(temp, alpha, temp2)
              CXADD(C[l + m * j], temp)
        }
      }
    } else  // tranB == CTran
    {
#pragma omp parallel for shared(A, B, C) private(temp, i, j, l, temp2)
      for (l = 0; l < m; l++) {
        for (j = 0; j < n; j++) {
          temp.re = 0;
          temp.im = 0;
          for (i = 0; i < k; i++) {
            temp.re += A[i * m + l].re * B[i * n + j].re;
            temp.re += A[i * m + l].im * B[i * n + j].im;
            temp.im += A[i * m + l].im * (B[i * n + j].re);
            temp.im -= A[i * m + l].re * (B[i * n + j].im);
          }
          CXMUL(C[l + m * j], beta, temp2) CXMUL(temp, alpha, temp2)
              CXADD(C[l + m * j], temp)
        }
      }
    }
  } else if (transA == Tran) {
    if (transB == NoTran) {
#pragma omp parallel for shared(A, B, C) private(temp, i, j, l, temp2)
      for (l = 0; l < m; l++) {
        for (j = 0; j < n; j++) {
          temp.re = 0;
          temp.im = 0;
          for (i = 0; i < k; i++) {
            CXADD_mul(temp, A[i + l * k], B[i + j * k])
          }
          CXMUL(C[l + m * j], beta, temp2) CXMUL(temp, alpha, temp2)
              CXADD(C[l + m * j], temp)
        }
      }

    } else if (transB == Tran) {
#pragma omp parallel for shared(A, B, C) private(temp, i, j, l, temp2)
      for (l = 0; l < m; l++) {
        for (j = 0; j < n; j++) {
          temp.re = 0;
          temp.im = 0;
          for (i = 0; i < k; i++) {
            CXADD_mul(temp, A[i + l * k], B[i * n + j])
          }
          CXMUL(C[l + m * j], beta, temp2) CXMUL(temp, alpha, temp2)
              CXADD(C[l + m * j], temp)
        }
      }
    } else  // transB == CTran
    {
#pragma omp parallel for shared(A, B, C) private(temp, i, j, l, temp2, adx, bdx)
      for (l = 0; l < m; l++) {
        for (j = 0; j < n; j++) {
          temp.re = 0;
          temp.im = 0;
          for (i = 0; i < k; i++) {
            adx = i + l * k;
            bdx = i * n + j;
            temp.re += A[adx].re * B[bdx].re - A[adx].im * (-B[bdx].im);
            temp.im += A[adx].re * (-B[bdx].im) + A[adx].im * (B[bdx].re);
          }
          CXMUL(C[l + m * j], beta, temp2) CXMUL(temp, alpha, temp2)
              CXADD(C[l + m * j], temp)
        }
      }
    }
  } else  // transA == CTran
  {
    if (transB == NoTran) {
#pragma omp parallel for shared(A, B, C) private(temp, i, j, l, temp2, adx, bdx)
      for (l = 0; l < m; l++) {
        for (j = 0; j < n; j++) {
          temp.re = 0;
          temp.im = 0;
          for (i = 0; i < k; i++) {
            adx = i + l * k;
            bdx = i + j * k;
            temp.re += A[adx].re * B[bdx].re + A[adx].im * (B[bdx].im);
            temp.im += A[adx].re * B[bdx].im - A[adx].im * (B[bdx].re);
          }
          CXMUL(C[l + m * j], beta, temp2) CXMUL(temp, alpha, temp2)
              CXADD(C[l + m * j], temp)
        }
      }

    } else if (transB == Tran) {
#pragma omp parallel for shared(A, B, C) private(i, j, l, adx, bdx, temp, temp2)
      for (l = 0; l < m; l++) {
        for (j = 0; j < n; j++) {
          temp.re = 0;
          temp.im = 0;
          for (i = 0; i < k; i++) {
            adx = i + l * k;
            bdx = i * n + j;
            temp.re += A[adx].re * B[bdx].re + (A[adx].im) * B[bdx].im;
            temp.im += A[adx].re * B[bdx].im - A[adx].im * (B[bdx].re);
          }
          CXMUL(C[l + m * j], beta, temp2) CXMUL(temp, alpha, temp2)
              CXADD(C[l + m * j], temp)
        }
      }
    } else  // transB == CTran
    {
#pragma omp parallel for shared(A, B, C) private(i, j, l, adx, bdx, temp, temp2)
      for (l = 0; l < m; l++) {
        for (j = 0; j < n; j++) {
          temp.re = 0;
          temp.im = 0;
          for (i = 0; i < k; i++) {
            adx = i + l * k;
            bdx = i * n + j;
            temp.re += A[adx].re * B[bdx].re - A[adx].im * B[bdx].im;
            temp.im -= A[adx].re * B[bdx].im + A[adx].im * B[bdx].re;
          }
          CXMUL(C[l + m * j], beta, temp2) CXMUL(temp, alpha, temp2)
              CXADD(C[l + m * j], temp)
        }
      }
    }
  }
}

/**** REAL ****/

void aABpbC(DTYPE alpha, MAT* A, MAT* B, DTYPE beta, MAT* C) {
#if DEBUG
  printf("%s\n", __func__);
#endif
  gemm(NoTran, NoTran, alpha, A, B, beta, C);
}
void aABtpbC(DTYPE alpha, MAT* A, MAT* B, DTYPE beta, MAT* C) {
#if DEBUG
  printf("%s\n", __func__);
#endif
  gemm(NoTran, Tran, alpha, A, B, beta, C);
}

void aAtBpbC(DTYPE alpha, MAT* A, MAT* B, DTYPE beta, MAT* C) {
#if DEBUG
  printf("%s\n", __func__);
#endif
  gemm(Tran, NoTran, alpha, A, B, beta, C);
}
void aAtBtpbC(DTYPE alpha, MAT* A, MAT* B, DTYPE beta, MAT* C) {
#if DEBUG
  printf("%s\n", __func__);
#endif
  gemm(Tran, Tran, alpha, A, B, beta, C);
}

void matmul(MAT* A, MAT* B, MAT* C) {
#if DEBUG
  printf("%s\n", __func__);
#endif
  gemm(NoTran, NoTran, 1, A, B, 0, C);
}

/**** COMPLEX ****/

void caABpbC(CTYPE alpha, CMAT* A, CMAT* B, CTYPE beta, CMAT* C) {
#if DEBUG
  printf("%s\n", __func__);
#endif
  cgemm(NoTran, NoTran, alpha, A, B, beta, C);
}
void caABtpbC(CTYPE alpha, CMAT* A, CMAT* B, CTYPE beta, CMAT* C) {
#if DEBUG
  printf("%s\n", __func__);
#endif
  cgemm(NoTran, Tran, alpha, A, B, beta, C);
}
void caABhpbC(CTYPE alpha, CMAT* A, CMAT* B, CTYPE beta, CMAT* C) {
#if DEBUG
  printf("%s\n", __func__);
#endif
  cgemm(NoTran, CTran, alpha, A, B, beta, C);
}

void caAtBpbC(CTYPE alpha, CMAT* A, CMAT* B, CTYPE beta, CMAT* C) {
#if DEBUG
  printf("%s\n", __func__);
#endif
  cgemm(Tran, NoTran, alpha, A, B, beta, C);
}
void caAtBtpbC(CTYPE alpha, CMAT* A, CMAT* B, CTYPE beta, CMAT* C) {
#if DEBUG
  printf("%s\n", __func__);
#endif
  cgemm(Tran, Tran, alpha, A, B, beta, C);
}
void caAtBhpbC(CTYPE alpha, CMAT* A, CMAT* B, CTYPE beta, CMAT* C) {
#if DEBUG
  printf("%s\n", __func__);
#endif
  cgemm(Tran, CTran, alpha, A, B, beta, C);
}

void caAhBpbC(CTYPE alpha, CMAT* A, CMAT* B, CTYPE beta, CMAT* C) {
#if DEBUG
  printf("%s\n", __func__);
#endif
  cgemm(CTran, NoTran, alpha, A, B, beta, C);
}
void caAhBtpbC(CTYPE alpha, CMAT* A, CMAT* B, CTYPE beta, CMAT* C) {
#if DEBUG
  printf("%s\n", __func__);
#endif
  cgemm(CTran, Tran, alpha, A, B, beta, C);
}
void caAhBhpbC(CTYPE alpha, CMAT* A, CMAT* B, CTYPE beta, CMAT* C) {
#if DEBUG
  printf("%s\n", __func__);
#endif
  cgemm(CTran, CTran, alpha, A, B, beta, C);
}

void cmatmul(CMAT* A, CMAT* B, CMAT* C) {
#if DEBUG
  printf("%s\n", __func__);
#endif
  CTYPE one_zero;
  CTYPE zero_zero;
  one_zero.re = 1;
  one_zero.im = 0;
  zero_zero.re = 0;
  zero_zero.im = 0;
  cgemm(NoTran, NoTran, one_zero, A, B, zero_zero, C);
}
