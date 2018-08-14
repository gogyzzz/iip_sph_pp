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
#include "iip_math.h"
/**** SQUARE ROOT ****/
void sqrt_mat(MAT* mat) {
  sqrt_mat_inc(mat->d0 * mat->d1 * mat->d2, mat->data, 1);
}

void sqrt_mat_inc(UINT size, DTYPE* X, ITER inc) {
  ITER i;
#if DEBUG
  printf("%s\n", __func__);
#endif

//#pragma omp parallel for shared(X) private(i)
  for (i = 0; i < size; i += inc) {
#if NTYPE == 0
    X[i] = sqrtf(fabs(X[i]));
#elif NTYPE == 1
    X[i] = X[i] < 0 ? sqrt(-X[i]) : sqrt(X[i]);
#endif
  }
}

void sqrt_cmat(CMAT* mat) {
  sqrt_cmat_inc(mat->d0 * mat->d1 * mat->d2, mat->data, 1);
}

void sqrt_cmat_inc(UINT size, CTYPE* X, ITER incx) {
  ITER i;
#if DEBUG
  printf("%s\n", __func__);
#endif

#pragma omp parallel for shared(X) private(i)
  for (i = 0; i < size; i += incx) {
#if OS_WIN
#if NTYPE == 0
    *(_Fcomplex*)(&X[i]) = csqrtf(*(_Fcomplex*)(&X[i]));
#elif NTYPE == 1
    *(_Dcomplex*)(&X[i]) = csqrt(*(_Dcomplex*)(&X[i]));
#endif
#elif OS_UNIX
#if NTYPE == 0
    *(complex float*)(&X[i]) = csqrtf(*(complex float*)(&X[i]));
#elif NTYPE == 1
    *(complex double*)(&X[i]) = csqrt(*(complex double*)(&X[i]));
#endif
#endif
  }
}

/**** POWER ****/

void pow_mat(MAT* mat, DTYPE n) {
  pow_mat_inc(mat->d0 * mat->d1 * mat->d2, mat->data, n, 1);
}

void pow_mat_inc(UINT size, DTYPE* X, DTYPE n, ITER incx) {
  ITER i;
#if DEBUG
  printf("%s\n", __func__);
#endif

#pragma omp parallel for shared(X) private(i)
  for (i = 0; i < size; i += incx) {
#if NTYPE == 0
    X[i] = powf(X[i], n);
#elif NTYPE == 1
    X[i] = pow(X[i], n);
#endif
  }
}
void pow_cmat(CMAT* mat, DTYPE n) {
  pow_cmat_inc(mat->d0 * mat->d1 * mat->d2, mat->data, n, 1);
}

void pow_cmat_inc(UINT size, CTYPE* X, DTYPE n, ITER incx) {
  ITER i;
  CTYPE temp;
  temp.re = n;
  temp.im = 0;
#if DEBUG
  printf("%s\n", __func__);
#endif

#pragma omp parallel for shared(X) private(i)
  for (i = 0; i < size; i += incx) {
#if OS_WIN
#if NTYPE == 0
    *(_Fcomplex*)(&X[i]) = cpowf(*(_Fcomplex*)(&X[i]), *(_Fcomplex*)(&temp));
#elif NTYPE == 1
    *(_Dcomplex*)(&X[i]) = cpow(*(_Dcomplex*)(&X[i]), *(_Dcomplex*)(&temp));
#endif
#elif OS_UNIX
#if NTYPE == 0
    *(complex float*)(&X[i]) =
        cpowf(*(complex float*)(&X[i]), *(complex float*)(&temp));
#elif NTYPE == 1
    *(complex double*)(&X[i]) =
        cpow(*(complex double*)(&X[i]), *(complex double*)(&temp));
#endif
#endif
  }
}

void cpow_cmat(CMAT* mat, CTYPE n) {
  cpow_cmat_inc(mat->d0 * mat->d1 * mat->d2, mat->data, n, 1);
}

void cpow_cmat_inc(UINT size, CTYPE* X, CTYPE n, ITER incx) {
  ITER i;
#if DEBUG
  printf("%s\n", __func__);
#endif

#pragma omp parallel for shared(X) private(i)
  for (i = 0; i < size; i += incx) {
#if OS_WIN
#if NTYPE == 0
    *(_Fcomplex*)(&X[i]) = cpowf(*(_Fcomplex*)(&X[i]), *(_Fcomplex*)(&n));
#elif NTYPE == 1
    *(_Dcomplex*)(&X[i]) = cpow(*(_Dcomplex*)(&X[i]), *(_Dcomplex*)(&n));
#endif
#elif OS_UNIX
#if NTYPE == 0
    *(complex float*)(&X[i]) =
        cpowf(*(complex float*)(&X[i]), *(complex float*)(&n));
#elif NTYPE == 1
    *(complex double*)(&X[i]) =
        cpow(*(complex double*)(&X[i]), *(complex double*)(&n));

#endif
#endif
  }
}

/**** Uniform distribution ****/

void randu(MAT* mat, DTYPE a, DTYPE b) {
  randu_inc(mat->d0 * mat->d1 * mat->d2, mat->data, a, b, 1);
}

void randu_inc(UINT size, DTYPE* X, DTYPE a, DTYPE b, ITER incx) {
  ITER i;
#if DEBUG
  printf("%s\n", __func__);
#endif
  srand(get_micro_sec());

#pragma omp parallel for shared(X) private(i)
  for (i = 0; i < size; i += incx) {
#if NTYPE == 0
    X[i] = ((float)rand() / RAND_MAX) * (b - a) + a;
#elif NTYPE == 1
    X[i] = ((double)rand() / RAND_MAX) * (b - a) + a;
#endif
  }
}

void crandu(CMAT* mat, DTYPE ra, DTYPE rb, DTYPE ia, DTYPE ib) {
  crandu_inc(mat->d0 * mat->d1 * mat->d2, mat->data, ra, rb, ia, ib, 1);
}

void crandu_inc(UINT size, CTYPE* X, DTYPE ra, DTYPE rb, DTYPE ia, DTYPE ib,
                ITER incx) {
  ITER i;
#if DEBUG
  printf("%s\n", __func__);
#endif
  srand(get_micro_sec());
#pragma omp parallel for shared(X) private(i)
  for (i = 0; i < size; i += incx) {
#if NTYPE == 0
    X[i].re = ((float)rand() / RAND_MAX) * (rb - ra) + ra;
    X[i].im = ((float)rand() / RAND_MAX) * (ib - ia) + ia;
#elif NTYPE == 1
    X[i].re = ((double)rand() / RAND_MAX) * (rb - ra) + ra;
    X[i].im = ((double)rand() / RAND_MAX) * (ib - ia) + ia;
#endif
  }
}

void randn(MAT* mat, DTYPE mean, DTYPE std) {
  randn_inc(mat->d0 * mat->d1 * mat->d2, mat->data, mean, std, 1);
}
/*
 * Using Box-Muller Transform, But I heard ziggurat is fastest
 * */
void randn_inc(UINT size, DTYPE* X, DTYPE mean, DTYPE std, ITER incx) {
  ITER i;
  DTYPE u, v;
  DTYPE s;
#if DEBUG
  printf("%s\n", __func__);
#endif

  srand(get_micro_sec());
#pragma omp parallel for shared(X) private(i, u, v, s)
  for (i = 0; i < size; i += incx) {
    do {
      u = (2.0 * rand() / RAND_MAX) - 1;
      v = (2.0 * rand() / RAND_MAX) - 1;
      s = u * u + v * v;
    }
    // check if s is in unit circle or is zero
    while (s >= 1.0 || s == 0);

// You can use u,v both
#if NTYPE == 0
    X[i] = std * (u * sqrtf(-2.0 * logf(s) / s)) + mean;
#elif NTYPE == 1
    X[i] = std * (u * sqrt(-2.0 * log(s) / s)) + mean;
#endif
  }
}

void crandn(CMAT* mat, CTYPE mean, CTYPE std) {
  crandn_inc(mat->d0 * mat->d1 * mat->d2, mat->data, mean, std, 1);
}
void crandn_inc(UINT size, CTYPE* X, CTYPE mean, CTYPE std, ITER incx) {
  ITER i;
  DTYPE u, v;
  DTYPE s;
#if DEBUG
  printf("%s\n", __func__);
#endif

  srand(get_micro_sec());
#pragma omp parallel for shared(X) private(i, u, v, s)
  for (i = 0; i < size; i += incx) {
    do {
      u = (2.0 * rand() / RAND_MAX) - 1;
      v = (2.0 * rand() / RAND_MAX) - 1;
      s = u * u + v * v;
    }
    // check if s is in unit circle or is zero
    while (s >= 1.0 || s == 0);

// You can use u,v both
#if NTYPE == 0
    X[i].re = std.re * (u * sqrtf(-2.0 * logf(s) / s)) + mean.re;
    X[i].im = std.im * (v * sqrtf(-2.0 * logf(s) / s)) + mean.im;
#elif NTYPE == 1
    X[i].re = std.re * (u * sqrt(-2.0 * log(s) / s)) + mean.re;
    X[i].im = std.im * (u * sqrt(-2.0 * log(s) / s)) + mean.im;
#endif
  }
}

/**** round ****/
void round_mat(MAT* mat) {
  round_mat_inc(mat->d0 * mat->d1 * mat->d2, mat->data, 1);
}
void round_mat_inc(UINT size, DTYPE* X, ITER incx) {
  ITER i;
#if DEBUG
  printf("%s\n", __func__);
#endif
#pragma omp parallel for shared(X) private(i)
  for (i = 0; i < size; i += incx) {
#if NTYPE == 0
    X[i] = roundf(X[i]);
#elif NTYPE == 1
    X[i] = round(X[i]);
#endif
  }
}

void round_cmat(CMAT* mat) {
  round_cmat_inc(mat->d0 * mat->d1 * mat->d2, mat->data, 1);
}
void round_cmat_inc(UINT size, CTYPE* X, ITER incx) {
  ITER i;
#if DEBUG
  printf("%s\n", __func__);
#endif
#pragma omp parallel for shared(X) private(i)
  for (i = 0; i < size; i += incx) {
#if NTYPE == 0
    X[i].re = roundf(X[i].re);
    X[i].im = roundf(X[i].im);
#elif NTYPE == 1
    X[i].re = round(X[i].re);
    X[i].im = round(X[i].im);
#endif
  }
}
/**** floor ****/
void floor_mat(MAT* mat) {
  floor_mat_inc(mat->d0 * mat->d1 * mat->d2, mat->data, 1);
}
void floor_mat_inc(UINT size, DTYPE* X, ITER incx) {
  ITER i;
#if DEBUG
  printf("%s\n", __func__);
#endif
#pragma omp parallel for shared(X) private(i)
  for (i = 0; i < size; i += incx) {
#if NTYPE == 0
    X[i] = floorf(X[i]);
#elif NTYPE == 1
    X[i] = floor(X[i]);
#endif
  }
}

void floor_cmat(CMAT* mat) {
  floor_cmat_inc(mat->d0 * mat->d1 * mat->d2, mat->data, 1);
}
void floor_cmat_inc(UINT size, CTYPE* X, ITER incx) {
  ITER i;
#if DEBUG
  printf("%s\n", __func__);
#endif
#pragma omp parallel for shared(X) private(i)
  for (i = 0; i < size; i += incx) {
#if NTYPE == 0
    X[i].re = floorf(X[i].re);
    X[i].im = floorf(X[i].im);
#elif NTYPE == 1
    X[i].re = floor(X[i].re);
    X[i].im = floor(X[i].im);
#endif
  }
}
/**** ceil ****/
void ceil_mat(MAT* mat) {
  ceil_mat_inc(mat->d0 * mat->d1 * mat->d2, mat->data, 1);
}
void ceil_mat_inc(UINT size, DTYPE* X, ITER incx) {
  ITER i;
#if DEBUG
  printf("%s\n", __func__);
#endif
#pragma omp parallel for shared(X) private(i)
  for (i = 0; i < size; i += incx) {
#if NTYPE == 0
    X[i] = ceilf(X[i]);
#elif NTYPE == 1
    X[i] = ceil(X[i]);
#endif
  }
}

void ceil_cmat(CMAT* mat) {
  ceil_cmat_inc(mat->d0 * mat->d1 * mat->d2, mat->data, 1);
}
void ceil_cmat_inc(UINT size, CTYPE* X, ITER incx) {
  ITER i;
#if DEBUG
  printf("%s\n", __func__);
#endif
#pragma omp parallel for shared(X) private(i)
  for (i = 0; i < size; i += incx) {
#if NTYPE == 0
    X[i].re = ceilf(X[i].re);
    X[i].im = ceilf(X[i].im);
#elif NTYPE == 1
    X[i].re = ceil(X[i].re);
    X[i].im = ceil(X[i].im);
#endif
  }
}
/**** log ****/
void log_mat(MAT* mat) {
  log_mat_inc(mat->d0 * mat->d1 * mat->d2, mat->data, 1);
}
void log_mat_inc(UINT size, DTYPE* X, ITER incx) {
  ITER i;
#if DEBUG
  printf("%s\n", __func__);
#endif
#pragma omp parallel for shared(X) private(i)
  for (i = 0; i < size; i += incx) {
#if NTYPE == 0
    X[i] = logf(X[i] < 0 ? -X[i] : X[i]);
#elif NTYPE == 1
    X[i] = log(X[i] < 0 ? -X[i] : X[i]);
#endif
  }
}

void log_cmat(CMAT* mat) {
  log_cmat_inc(mat->d0 * mat->d1 * mat->d2, mat->data, 1);
}
void log_cmat_inc(UINT size, CTYPE* X, ITER incx) {
  ITER i;
#if DEBUG
  printf("%s\n", __func__);
#endif
#pragma omp parallel for shared(X) private(i)
  for (i = 0; i < size; i += incx) {
#if NTYPE == 0
    CXF(X[i]) = clogf(CXF(X[i]));
#elif NTYPE == 1
    CXD(X[i]) = clog(CXD(X[i]));
#endif
  }
}

void log2_mat(MAT* mat) {
  log2_mat_inc(mat->d0 * mat->d1 * mat->d2, mat->data, 1);
}
void log2_mat_inc(UINT size, DTYPE* X, ITER incx) {
  ITER i;
#if DEBUG
  printf("%s\n", __func__);
#endif
#pragma omp parallel for shared(X) private(i)
  for (i = 0; i < size; i += incx) {
#if NTYPE == 0
    X[i] = log2f(X[i] < 0 ? -X[i] : X[i]);
#elif NTYPE == 1
    X[i] = log2(X[i] < 0 ? -X[i] : X[i]);
#endif
  }
}

void log2_cmat(CMAT* mat) {
  log2_cmat_inc(mat->d0 * mat->d1 * mat->d2, mat->data, 1);
}
void log2_cmat_inc(UINT size, CTYPE* X, ITER incx) {
  ITER i;
#if DEBUG
  printf("%s\n", __func__);
#endif
#pragma omp parallel for shared(X) private(i)
  for (i = 0; i < size; i += incx) {
/*	There is no clog2
 *	But clog2(cx) = clog(cx)/log(2)
 *
 * */
#if NTYPE == 0
    CXF(X[i]) = clogf(CXF(X[i]));
    X[i].re = X[i].re / logf(2.);
    X[i].im = X[i].im / logf(2.);
#elif NTYPE == 1
    CXD(X[i]) = clog(CXD(X[i]));
    X[i].re = X[i].re / log(2.);
    X[i].im = X[i].im / log(2.);
#endif
  }
}

void log10_mat(MAT* mat) {
  log10_mat_inc(mat->d0 * mat->d1 * mat->d2, mat->data, 1);
}
void log10_mat_inc(UINT size, DTYPE* X, ITER incx) {
  ITER i;
#if DEBUG
  printf("%s\n", __func__);
#endif
#pragma omp parallel for shared(X) private(i)
  for (i = 0; i < size; i += incx) {
#if NTYPE == 0
    X[i] = log10f(X[i] < 0? -X[i] : X[i]);
#elif NTYPE == 1
    X[i] = log10(X[i] < 0 ? -X[i] : X[i]);
#endif
  }
}

void log10_cmat(CMAT* mat) {
  log10_cmat_inc(mat->d0 * mat->d1 * mat->d2, mat->data, 1);
}
void log10_cmat_inc(UINT size, CTYPE* X, ITER incx) {
  ITER i;
#if DEBUG
  printf("%s\n", __func__);
#endif
#pragma omp parallel for shared(X) private(i)
  for (i = 0; i < size; i += incx) {
#if NTYPE == 0
    CXF(X[i]) = clog10f(CXF(X[i]));
#elif NTYPE == 1
    CXD(X[i]) = clog10(CXD(X[i]));
#endif
  }
}

/**** log with base ****/
void logb_mat(MAT*mat, UINT base){
  logb_mat_inc(mat->d0*mat->d1*mat->d2,mat->data,base,1);
}

void logb_mat_inc(UINT size,DTYPE*X , UINT base, ITER incx){
  ITER i;
#if DEBUG
  printf("%s\n", __func__);
#endif
#pragma omp parallel for shared(X) private(i)
  for (i = 0; i < size; i += incx) {
#if NTYPE == 0
    X[i] = logf(X[i] < 0 ? -X[i] : X[i]) / logf(base);
#elif NTYPE == 1
    X[i] = log(X[i] < 0 ? -X[i] : X[i]) / log(base);
#endif
  }
 
}

void clogb_mat(CMAT*mat, UINT base){
  clogb_mat_inc(mat->d0*mat->d1*mat->d2,mat->data,base,1);
}

void clogb_mat_inc(UINT size,CTYPE*X , UINT base, ITER incx){
  ITER i;
  DTYPE temp;
#if DEBUG
  printf("%s\n", __func__);
#endif
#pragma omp parallel for shared(X) private(i)
  for (i = 0; i < size; i += incx) {
#if NTYPE == 0
    CXF(X[i]) = clogf(CXF(X[i]));
    temp = logf(base);
    X[i].re /= temp;
    X[i].im /= temp;
#elif NTYPE == 1
    CXD(X[i]) = clog(CXD(X[i]));
    temp = log(base);
    X[i].re /= temp;
    X[i].im /= temp;
#endif
  }
}

/**** exp ****/
void exp_mat(MAT* mat) {
  exp_mat_inc(mat->d0 * mat->d1 * mat->d2, mat->data, 1);
}
void exp_mat_inc(UINT size, DTYPE* X, ITER incx) {
  ITER i;
#if DEBUG
  printf("%s\n", __func__);
#endif
#pragma omp parallel for shared(X) private(i)
  for (i = 0; i < size; i += incx) {
#if NTYPE == 0
    X[i] = expf(X[i]);
#elif NTYPE == 1
    X[i] = exp(X[i]);
#endif
  }
}

void exp_cmat(CMAT* mat) {
  exp_cmat_inc(mat->d0 * mat->d1 * mat->d2, mat->data, 1);
}
void exp_cmat_inc(UINT size, CTYPE* X, ITER incx) {
  ITER i;
#if DEBUG
  printf("%s\n", __func__);
#endif
#pragma omp parallel for shared(X) private(i)
  for (i = 0; i < size; i += incx) {
#if NTYPE == 0
    CXF(X[i]) = cexpf(CXF(X[i]));
#elif NTYPE == 1
    CXD(X[i]) = cexp(CXD(X[i]));
#endif
  }
}
/**** abs ****/
void abs_mat(MAT* mat) {
  abs_mat_inc(mat->d0 * mat->d1 * mat->d2, mat->data, 1);
}
void abs_mat_inc(UINT size, DTYPE* X, ITER incx) {
  ITER i;
#if DEBUG
  printf("%s\n", __func__);
#endif
#pragma omp parallel for shared(X) private(i)
  for (i = 0; i < size; i += incx) {
#if NTYPE == 0
    X[i] = fabsf(X[i]);
#elif NTYPE == 1
    X[i] = fabs(X[i]);
#endif
  }
}

void abs_cmat(CMAT* mat) {
  abs_cmat_inc(mat->d0 * mat->d1 * mat->d2, mat->data, 1);
}
void abs_cmat_inc(UINT size, CTYPE* X, ITER incx) {
  ITER i;
#if DEBUG
  printf("%s\n", __func__);
#endif
#pragma omp parallel for shared(X) private(i)
  for (i = 0; i < size; i += incx) {
#if NTYPE == 0
    X[i].re = cabsf(CXF(X[i]));
    X[i].im = 0;
#elif NTYPE == 1
    X[i].re = cabs(CXD(X[i]));
    X[i].im = 0;
#endif
  }
}

/**** max ****/
DTYPE max_mat(MAT* mat, DIM* dim) {
  ITER i;
  DTYPE max = mat->data[0];
  UINT max_idx = 0;
#if DEBUG
  printf("%s\n", __func__);
#endif
  for (i = 0; i < mat->d0 * mat->d1 * mat->d2; i++) {
    if (mat->data[i] > max) {
      max = mat->data[i];
      max_idx = i;
    }
  }
  dim->d2 = max_idx / (mat->d0 * mat->d1);
  dim->d1 = (max_idx % (mat->d0 * mat->d1)) / mat->d0;
  dim->d0 = max_idx % mat->d0;
  return max;
}

CTYPE max_cmat(CMAT* mat, DIM* dim) {
  ITER i;
  CTYPE max = mat->data[0];
  UINT max_idx = 0;
#if DEBUG
  printf("%s\n", __func__);
#endif
  for (i = 0; i < mat->d0 * mat->d1 * mat->d2; i++) {
    if (mat->data[i].re >= max.re) {
      if (mat->data[i].re == max.re) {
        if (mat->data[i].im > max.im) {
          max = mat->data[i];
          max_idx = i;
        }
      } else {
        max = mat->data[i];
        max_idx = i;
      }
    }
  }
  dim->d2 = max_idx / (mat->d0 * mat->d1);
  dim->d1 = (max_idx % (mat->d0 * mat->d1)) / mat->d0;
  dim->d0 = max_idx % mat->d0;
  return max;
}

/**** absolute max ****/
DTYPE amax_mat(MAT* mat, DIM* dim) {
  ITER i;
  UINT max_idx = 0;
#if NTYPE == 0
  DTYPE max = fabsf(mat->data[0]);
#if DEBUG
  printf("%s\n", __func__);
#endif
  for (i = 0; i < mat->d0 * mat->d1 * mat->d2; i++) {
    if (fabsf(mat->data[i]) > max) {
      max = fabsf(mat->data[i]);
      max_idx = i;
    }
  }
#elif NTYPE == 1
  DTYPE max = fabs(mat->data[0]);
#if DEBUG
  printf("%s\n", __func__);
#endif
  for (i = 0; i < mat->d0 * mat->d1 * mat->d2; i++) {
    if (fabs(mat->data[i]) > max) {
      max = fabs(mat->data[i]);
      max_idx = i;
    }
  }
#endif
  dim->d2 = max_idx / (mat->d0 * mat->d1);
  dim->d1 = (max_idx % (mat->d0 * mat->d1)) / mat->d0;
  dim->d0 = max_idx % mat->d0;
  return max;
}

DTYPE amax_cmat(CMAT* mat, DIM* dim) {
  ITER i;
  UINT max_idx = 0;
  DTYPE max;
#if DEBUG
  printf("%s\n", __func__);
#endif
#if NTYPE == 0
  max = cabsf(CXF(mat->data[0]));
  for (i = 0; i < mat->d0 * mat->d1 * mat->d2; i++) {
    if (cabsf(CXF(mat->data[i])) > max) {
      max = cabsf(CXF(mat->data[i]));
      max_idx = i;
    }
  }
#elif NTYPE == 1
  max = cabs(CXD(mat->data[0]));
  for (i = 0; i < mat->d0 * mat->d1 * mat->d2; i++) {
    if (cabs(CXD(mat->data[i])) > max) {
      max = cabs(CXD(mat->data[i]));
      max_idx = i;
    }
  }
#endif
  dim->d2 = max_idx / (mat->d0 * mat->d1);
  dim->d1 = (max_idx % (mat->d0 * mat->d1)) / mat->d0;
  dim->d0 = max_idx % mat->d0;
  return max;
}

/**** min ****/
DTYPE min_mat(MAT* mat, DIM* dim) {
  ITER i;
  DTYPE min = mat->data[0];
  UINT min_idx = 0;
#if DEBUG
  printf("%s\n", __func__);
#endif
  for (i = 0; i < mat->d0 * mat->d1 * mat->d2; i++) {
    if (mat->data[i] < min) {
      min = mat->data[i];
      min_idx = i;
    }
  }
  dim->d2 = min_idx / (mat->d0 * mat->d1);
  dim->d1 = (min_idx % (mat->d0 * mat->d1)) / mat->d0;
  dim->d0 = min_idx % mat->d0;
  return min;
}

CTYPE min_cmat(CMAT* mat, DIM* dim) {
  ITER i;
  CTYPE min = mat->data[0];
  UINT min_idx = 0;
#if DEBUG
  printf("%s\n", __func__);
#endif
  for (i = 0; i < mat->d0 * mat->d1 * mat->d2; i++) {
    if (mat->data[i].re <= min.re) {
      if (mat->data[i].re == min.re) {
        if (mat->data[i].im < min.im) {
          min = mat->data[i];
          min_idx = i;
        }
      } else {
        min = mat->data[i];
        min_idx = i;
      }
    }
  }
  dim->d2 = min_idx / (mat->d0 * mat->d1);
  dim->d1 = (min_idx % (mat->d0 * mat->d1)) / mat->d0;
  dim->d0 = min_idx % mat->d0;
  return min;
}

/**** absolute min ****/
DTYPE amin_mat(MAT* mat, DIM* dim) {
  ITER i;
  UINT min_idx = 0;
#if NTYPE == 0
  DTYPE min = absf(mat->data[0]);
#if DEBUG
  printf("%s\n", __func__);
#endif
  for (i = 0; i < mat->d0 * mat->d1 * mat->d2; i++) {
    if (absf(mat->data[i]) < min) {
      min = absf(mat->data[i]);
      min_idx = i;
    }
  }
#elif NTYPE == 1
  DTYPE min = abs(mat->data[0]);
#if DEBUG
  printf("%s\n", __func__);
#endif
  for (i = 0; i < mat->d0 * mat->d1 * mat->d2; i++) {
    if (abs(mat->data[i]) < min) {
      min = abs(mat->data[i]);
      min_idx = i;
    }
  }
#endif
  dim->d2 = min_idx / (mat->d0 * mat->d1);
  dim->d1 = (min_idx % (mat->d0 * mat->d1)) / mat->d0;
  dim->d0 = min_idx % mat->d0;
  return min;
}

DTYPE amin_cmat(CMAT* mat, DIM* dim) {
  ITER i;
  UINT min_idx = 0;
  DTYPE min;
#if DEBUG
  printf("%s\n", __func__);
#endif
#if NTYPE == 0
  min = cabsf(CXF(mat->data[0]));
  for (i = 0; i < mat->d0 * mat->d1 * mat->d2; i++) {
    if (cabsf(CXF(mat->data[i])) < min) {
      min = cabsf(CXF(mat->data[i]));
      min_idx = i;
    }
  }
#elif NTYPE == 1
  min = cabs(CXD(mat->data[0]));
  for (i = 0; i < mat->d0 * mat->d1 * mat->d2; i++) {
    if (cabs(CXD(mat->data[i])) < min) {
      min = cabs(CXD(mat->data[i]));
      min_idx = i;
    }
  }
#endif
  dim->d2 = min_idx / (mat->d0 * mat->d1);
  dim->d1 = (min_idx % (mat->d0 * mat->d1)) / mat->d0;
  dim->d0 = min_idx % mat->d0;
  return min;
}

/**** accumulated sum ****/
void sum_mat(MAT*src, MAT*des, UINT axis){
ITER i,j,k;
DTYPE temp;
UINT d0,d1,d2;
#if DEBUG
  printf("%s\n", __func__);
#endif
if(src->d2 != des->d2) ASSERT_DIM_INVALID()
d0 = src->d0;
d1 = src->d1;
d2 = src->d2;
  if(axis == 0){ // d0 d1 d2 -> d0 1 d2
    if(d0 != des->d0) ASSERT_DIM_INVALID()
    if(des->d1 != 1) ASSERT_DIM_INVALID()
  
    for(k=0;k<d2;k++){
#pragma omp parallel for shared(des) private(i,j,temp)
      for(j=0;j<d0;j++){
        temp = 0;
        for(i=0;i<d1;i++){
          temp += src->data[k*d0*d1 + i*d0 + j];
        }  
        des->data[k*d0*d1 + j] = temp;
      }
    }

  }else if(axis == 1){ // d0 d1 d2 -> 1 d1 d2
    if(des->d0 != 1) ASSERT_DIM_INVALID()
    if(d1 != des->d1) ASSERT_DIM_INVALID()
    
    for(k=0;k<d2;k++){
#pragma omp parallel for shared(des) private(i,j,temp)
      for(j=0;j<d1;j++){
        temp = 0;
        for(i=0;i<d0;i++){
          temp += src->data[k*d0*d1 + j*d0 + i];
        }  
        des->data[k*d0*d1 + j] = temp;
      }
    } 
  }
  else ASSERT_ARG_INVALID()
}

void sum_cmat(CMAT*src, CMAT*des, UINT axis){
ITER i,j,k;
CTYPE temp;
UINT d0,d1,d2;
#if DEBUG
  printf("%s\n", __func__);
#endif
if(src->d2 != des->d2) ASSERT_DIM_INVALID()
d0 = src->d0;
d1 = src->d1;
d2 = src->d2;
  if(axis == 0){ // d0 d1 d2 -> d0 1 d2
    if(d0 != des->d0) ASSERT_DIM_INVALID()
    if(des->d1 != 1) ASSERT_DIM_INVALID()
  
    for(k=0;k<d2;k++){
#pragma omp parallel for shared(des) private(i,j,temp)
      for(j=0;j<d0;j++){
        temp.re = 0;
        temp.im = 0;
        for(i=0;i<d1;i++){
          temp.re += src->data[k*d0*d1 + i*d0 + j].re;
          temp.im += src->data[k*d0*d1 + i*d0 + j].im;
        }  
        des->data[k*d0*d1 + j] = temp;
      }
    }

  }else if(axis == 1){ // d0 d1 d2 -> 1 d1 d2
    if(des->d0 != 1) ASSERT_DIM_INVALID()
    if(d1 != des->d1) ASSERT_DIM_INVALID()
    
    for(k=0;k<d2;k++){
#pragma omp parallel for shared(des) private(i,j,temp)
      for(j=0;j<d1;j++){
        temp.re = 0;
        temp.im = 0;
        for(i=0;i<d0;i++){
          temp.re += src->data[k*d0*d1 + j*d0 + i].re;
          temp.im += src->data[k*d0*d1 + j*d0 + i].im;
        }  
        des->data[k*d0*d1 + j] = temp;
      }
    } 
  }
  else ASSERT_ARG_INVALID()
}

/**** misc****/
CTYPE CX(DTYPE r, DTYPE i) {
  CTYPE t;
  t.re = r;
  t.im = i;
  return t;
}
