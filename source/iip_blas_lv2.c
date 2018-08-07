#include "iip_blas_lv2.h"

/****  gemv ****/

void gemv(char transA, DTYPE alpha, MAT *A, MAT *X, DTYPE beta, MAT *Y) {
  UINT m, n, lda;
#if DEBUG
  printf("%s\n", __func__);
#endif

  if (X->ndim != 0 || Y->ndim != 0) {
    printf("A : %u  X : %u Y : %u\n", A->ndim, X->ndim, Y->ndim);
    printf("Use Vector for Vector operation\n");
    return;
  }
  if (A->ndim != 1) {
    printf("Use 2D-Matrix for BLAS operation\n");
    return;
  }

  m = A->d1;
  n = A->d0;
  lda = m;

#if DEBUG
  printf("trans : %d m : %u n: %u lda : %u\nalpha : %lf beta : %lf\n", transA,
         m, n, lda, alpha, beta);
#endif

#if USE_CBLAS
#if NTYPE == 0
  cblas_sgemv(CblasColMajor, transA, m, n, alpha, A->data, lda, X->data, 1,
              beta, Y->data, 1);
#else
  cblas_dgemv(CblasColMajor, transA, m, n, alpha, A->data, lda, X->data, 1,
              beta, Y->data, 1);

#endif
#else
  omp_gemv(transA, m, n, alpha, A->data, lda, X->data, 1, beta, Y->data, 1);
#endif
}

void omp_gemv(char tranA, UINT m, UINT n, DTYPE alpha, DTYPE *A, UINT lda,
             DTYPE *X, SINT incx, DTYPE beta, DTYPE *Y, SINT incy) {
  ITER i, j;
  DTYPE temp;

#if DEBUG
  printf("%s\n", __func__);
#endif

  if (tranA == Tran) {
#pragma omp parallel for shared(A, X, Y) private(temp, j, i)
    for (j = 0; j < n; j++) {
      temp = 0;
      for (i = 0; i < m; i++) {
        temp += A[j + i * n] * X[i * incx];
      }
      temp *= alpha;
      temp += beta * Y[j * incy];
      Y[j * incy] = temp;
    }
  } else if (tranA == NoTran) {
#pragma omp parallel for shared(A, X, Y) private(temp, j, i)
    for (j = 0; j < lda; j++) {
      temp = 0;
      for (i = 0; i < n; i++) {
        temp += A[i + n * j] * X[i * incx];
      }
      temp *= alpha;
      temp += beta * Y[j * incy];
      Y[j * incy] = temp;
    }
  } else {
    printf("ERROR : Transpose argument is invalid\n");
    return;
  }
}

void cgemv(char transA, CTYPE alpha, CMAT *A, CMAT *X, CTYPE beta, CMAT *Y) {
  UINT m, n, lda;
#if DEBUG
  printf("%s\n", __func__);
#endif

  if (X->ndim != 0 || Y->ndim != 0) {
    printf("A : %u  X : %u Y : %u\n", A->ndim, X->ndim, Y->ndim);
    printf("Use Vector for Vector operation\n");
    return;
  }
  if (A->ndim != 1) {
    printf("Use 2D-Matrix for BLAS operation\n");
    return;
  }

  m = A->d1;
  n = A->d0;
  lda = m;

#if DEBUG
  printf("trans : %d m : %u n: %u lda : %u\nalpha : %lf|%lf beta : %lf|%lf\n",
         transA, m, n, lda, alpha.re, alpha.im, beta.re, beta.im);
#endif

#if USE_CBLAS
#if NTYPE == 0
  cblas_cgemv(CblasColMajor, transA, m, n, alpha, A->data, lda, X->data, 1,
              beta, Y->data, 1);
#else
  cblas_zgemv(CblasColMajor, transA, m, n, &alpha, A->data, lda, X->data, 1,
              &beta, Y->data, 1);

#endif
#else
  omp_cgemv(transA, m, n, alpha, A->data, lda, X->data, 1, beta, Y->data, 1);
#endif
}

void omp_cgemv(char tranA, UINT m, UINT n, CTYPE alpha, CTYPE *A, UINT lda,
              CTYPE *X, SINT incx, CTYPE beta, CTYPE *Y, SINT incy) {
  ITER i, j;
  CTYPE temp;
  DTYPE temp2;
#if DEBUG
  printf("%s\n", __func__);
#endif

  if (tranA == Tran) {
#pragma omp parallel for shared(A, X, Y) private(temp, j, i)
    for (j = 0; j < n; j++) {
      temp.re = 0;
      temp.im = 0;
      for (i = 0; i < m; i++) {
        CXADD_mul(temp, A[j + i * n], X[i * incx]);
      }
      CXMUL(temp, alpha, temp2);
      CXADD_mul(temp, beta, Y[j * incy]);
      Y[j * incy] = temp;
    }
  } else if (tranA == NoTran) {
#pragma omp parallel for shared(A, X, Y) private(temp, j, i)
    for (j = 0; j < lda; j++) {
      temp.re = 0;
      temp.im = 0;
      for (i = 0; i < n; i++) {
        CXADD_mul(temp, A[i + n * j], X[i * incx]);
      }
      CXMUL(temp, alpha, temp2);
      CXADD_mul(temp, beta, Y[j * incy]);
      Y[j * incy] = temp;
    }
  } else {
    printf("ERROR : Transpose argument is invalid\n");
    return;
  }
}
