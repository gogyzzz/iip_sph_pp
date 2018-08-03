#ifndef BLAS_LV1_H
#define BLAS_LV1_H

//#include "mother.h"
#include "iip_type.h"

/*********************
**** iip_blas_lv1 ****
**********************/

#define ABS_CTYPE(x) ((DTYPE)fabs((double)x.re) + (DTYPE)fabs((double)x.im))

void mp_axpy(UINT, DTYPE, DTYPE *, UINT, DTYPE *, UINT);
void axpy_inc(UINT size, DTYPE alpha, DTYPE *X, ITER incx, DTYPE *Y, ITER incy);
void axpy(DTYPE alpha, MAT *x, MAT *y);

void mp_caxpy(UINT, CTYPE, CTYPE *, UINT, CTYPE *, UINT);
void caxpy(CTYPE alpha, CMAT *x, CMAT *y);
void caxpy_inc(UINT size,CTYPE alpha,CTYPE *X, ITER incx,CTYPE *Y, ITER incy);

#if USE_CUDA
__global__ void cu_axpy(DTYPE alpha, DTYPE *X, UINT INCX, DTYPE *Y, UINT INCY,
                        UINT len, UINT block_size);
__global__ void cu_caxpy(CTYPE alpha, CTYPE *X, UINT INCX, CTYPE *Y, UINT INCY,
                         UINT len, UINT block_size);
#endif

void mp_copy(UINT N, DTYPE *src, SINT src_inc, DTYPE *des, SINT des_inc);
void copy(MAT *src, MAT *des);

void mp_ccopy(UINT N, CTYPE *src, SINT src_inc, CTYPE *des, SINT des_inc);
void ccopy(CMAT *src, CMAT *des);

#if USE_CUDA
__global__ void cu_copy(DTYPE *SRC, UINT INC_SRC, DTYPE *DES, UINT INC_DES,
                        UINT len, UINT block_size);
__global__ void cu_ccopy(CTYPE *SRC, UINT INC_SRC, CTYPE *DES, UINT INC_DES,
                         UINT len, UINT block_size);
#endif

DTYPE mp_asum(UINT N, DTYPE *data, UINT inc);
DTYPE asum(MAT *mat, UINT inc);

DTYPE mp_casum(UINT N, CTYPE *data, UINT inc);
DTYPE casum(CMAT *mat, UINT inc);

#if USE_CUDA
__global__ void cu_asum(DTYPE *data, UINT inc, UINT len, UINT block_size,
                        DTYPE *sum);
__global__ void cu_casum(CTYPE *data, UINT inc, UINT len, UINT block_size,
                         CTYPE *sum);
#endif

DTYPE mp_dot(UINT N, DTYPE *src_x, UINT x_inc, DTYPE *src_y, UINT y_inc);
DTYPE dot(MAT *src_x, UINT x_increment, MAT *src_y, UINT y_increment);

CTYPE mp_cdot(UINT N, CTYPE *src_x, UINT x_inc, DTYPE *src_y, UINT y_inc);
CTYPE cdot(CMAT *src_x, UINT x_increment, MAT *src_y, UINT y_increment);

CTYPE mp_udot(UINT N, CTYPE *src_x, UINT x_inc, CTYPE *src_y, UINT y_inc);
CTYPE udot(CMAT *src_x, UINT x_increment, CMAT *src_y, UINT y_increment);

#if USE_CUDA
__global__ void cu_dot(DTYPE *result, DTYPE *src_x, UINT x_inc, DTYPE *src_y,
                       UINT y_inc, UINT len, UINT block_size);
#endif

void swap(MAT *src_x, MAT *src_y);
void swap_inc(MAT *src_x, UINT x_inc, MAT *src_y, UINT y_inc);
void mp_swap(UINT N, DTYPE *src_x, UINT x_inc, DTYPE *src_y, UINT y_inc);

void cswap(CMAT *src_x, CMAT *src_y);
void cswap_inc(CMAT *src_x, UINT x_inc, CMAT *src_y, UINT y_inc);
void mp_cswap(UINT N, CTYPE *src_x, UINT x_inc, CTYPE *src_y, UINT y_inc);

UINT amax(MAT *src);
UINT amax_inc(MAT *src, UINT inc);
UINT mp_amax(UINT N, DTYPE *src, UINT inc);

UINT camax(CMAT *src);
UINT camax_inc(CMAT *src, UINT inc);
UINT mp_camax(UINT N, CTYPE *src, UINT inc);

UINT amin(MAT *src);
UINT amin_inc(MAT *src, UINT inc);
UINT mp_amin(UINT N, DTYPE *src, UINT inc);

UINT camin(CMAT *src);
UINT camin_inc(CMAT *src, UINT inc);
UINT mp_camin(UINT N, CTYPE *src, UINT inc);

DTYPE cabs1(CTYPE val);



void rotg(DTYPE *a, DTYPE *b, DTYPE *c, DTYPE *s);
void mp_rotg(DTYPE *a, DTYPE *b, DTYPE *c, DTYPE *s);
void crotg(CTYPE *a, CTYPE *b, DTYPE *c, CTYPE *s);
void mp_crotg(CTYPE *a, CTYPE *b, DTYPE *c, CTYPE *s);

DTYPE nrm2(MAT *src);
DTYPE nrm2_inc(MAT *src, UINT inc);
DTYPE mp_nrm2(UINT N, DTYPE *data, UINT inc);

DTYPE cnrm2(CMAT *src);
DTYPE cnrm2_inc(CMAT *src, UINT inc);
DTYPE mp_cnrm2(UINT N, CTYPE *data, UINT inc);

void rot(MAT *src_x, MAT *src_y, DTYPE c, DTYPE s);
void rot_inc(MAT *src_x, UINT x_inc, MAT *src_y, UINT y_inc, DTYPE c, DTYPE s);
void mp_rot(UINT N, DTYPE *src_x, UINT x_inc, DTYPE *src_y, UINT y_inc, DTYPE c,
            DTYPE s);

void crot(CMAT *src_x, CMAT *src_y, DTYPE c, DTYPE s);
void crot_inc(CMAT *src_x, UINT x_inc, CMAT *src_y, UINT y_inc, DTYPE c,
              DTYPE s);
void mp_crot(UINT N, CTYPE *src_x, UINT x_inc, CTYPE *src_y, UINT y_inc,
             DTYPE c, DTYPE s);

/** 실수행렬 * 실수**/
void scal(DTYPE alpha, MAT *mat);
void scal_inc(UINT size, DTYPE alpha, DTYPE *X, UINT incx);
void mp_scal(UINT size, DTYPE alpha, DTYPE *X, UINT incx);
/** 복소수행렬 * 실수 **/
void cscal(DTYPE alpha, CMAT *mat);
void cscal_inc(UINT size, DTYPE alpha, CTYPE *X, UINT incx);
void mp_cscal(UINT size, DTYPE alpha, CTYPE *X, UINT incx);
/** 복소수행렬 * 복소수 **/
void uscal(CTYPE alpha, CMAT *mat);
void uscal_inc(UINT size, CTYPE alpha, CTYPE *X, UINT incx);
void mp_uscal(UINT size, CTYPE alpha, CTYPE *X, UINT incx);

/** 열 스케일링 **/
void col_scal(DTYPE alpha, MAT *X, UINT idx);
void col_cscal(DTYPE alpha, CMAT *X, UINT idx);
void col_uscal(CTYPE alpha, CMAT *X, UINT idx);
/** 행 스케일링 **/
void row_scal(DTYPE alpha, MAT *X, UINT idx);
void row_cscal(DTYPE alpha, CMAT *X, UINT idx);
void row_uscal(CTYPE alpha, CMAT *X, UINT idx);

/** 실수행렬 + 실수 **/
void add(DTYPE alpha, MAT *mat);
void add_inc(UINT size, DTYPE alpha, DTYPE *X, UINT incx);
/** 복소수행렬 + 실수 **/
void cadd(DTYPE alpha, CMAT *mat);
void cadd_inc(UINT size, DTYPE alpha, CTYPE *X, UINT incx);

/** 복소수행렬 + 복소수 **/
void uadd(CTYPE alpha, CMAT *mat);
void uadd_inc(UINT size, CTYPE alpha, CTYPE *X, UINT incx);

/** 열 더하기 **/
void col_add(DTYPE alpha, MAT *X, UINT idx);
void col_cadd(DTYPE alpha, CMAT *X, UINT idx);
void col_uadd(CTYPE alpha, CMAT *X, UINT idx);
/** 행 더하기 **/
void row_add(DTYPE alpha, MAT *X, UINT idx);
void row_cadd(DTYPE alpha, CMAT *X, UINT idx);
void row_uadd(CTYPE alpha, CMAT *X, UINT idx);

#endif
