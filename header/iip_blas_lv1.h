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
#ifndef BLAS_LV1_H
#define BLAS_LV1_H

//#include "mother.h"
#include "iip_type.h"

/* GENERAL RULE
 * {func} :  DTYPE matrix MAT operation
 * c{func} : {func} with CTYPE matrix CMAT
 * {func}_inc : {func} with DTYPE array with increment and size
 * c{func}_inc : c{func} with CTYPE array
 *
 * + If there is {func} about two mat
 * c{func} : {func} with CMAT and MAT
 * u{func} : {func} with CMAT and CMAT
 *
 * + in {func}_inc
 * size : how many elements?
 * increment(inc) : interval of array
 * ex)
 * MAT A(4,3)
 * {func}_inc(UINT size,DTYPE*X,ITER incx)
 *
 * {func}_inc(12,A->data,1) : for all elements, equal to {func}(A)
 * {func}_inc(4,A->data,1) : for first column
 * {func}_inc(4,&(A->data[4]),1) : for second column
 * {func}_inc(3,A,4) : for first row
 * {func}_inc(3,&(A->data[1]),4) : for second row
 *
 * By adjusting argument, row or column operation is available.
 * */

/*********************
**** iip_blas_lv1 ****
**********************/

#define ABS_CTYPE(x) ((DTYPE)fabs((double)x.re) + (DTYPE)fabs((double)x.im))

/* Y = AX + Y */
void axpy(DTYPE alpha, MAT *x, MAT *y);
void axpy_inc(UINT size, DTYPE alpha, DTYPE *X, ITER incx, DTYPE *Y, ITER incy);
void omp_axpy(UINT N, DTYPE alpha, DTYPE *X, UINT INCX, DTYPE *Y, UINT INCY);

/* axpy about complex number. */
void caxpy(CTYPE alpha, CMAT *x, CMAT *y);
void caxpy_inc(UINT size, CTYPE alpha, CTYPE *X, ITER incx, CTYPE *Y,
               ITER incy);
void omp_caxpy(UINT, CTYPE, CTYPE *, UINT, CTYPE *, UINT);

#if USE_CUDA
__global__ void cu_axpy(DTYPE alpha, DTYPE *X, UINT INCX, DTYPE *Y, UINT INCY,
                        UINT len, UINT block_size);
__global__ void cu_caxpy(CTYPE alpha, CTYPE *X, UINT INCX, CTYPE *Y, UINT INCY,
                         UINT len, UINT block_size);
#endif

/* Copy MAT from src to des */

void copy_mat(MAT *src, MAT *des);
void copy_mat_inc(UINT size,DTYPE*X,ITER incx, DTYPE *Y, ITER incy);
void omp_copy_mat(UINT N, DTYPE *src, SINT src_inc, DTYPE *des, SINT des_inc);

/* Copy CMAT from src to des */
void ccopy_mat(CMAT *src, CMAT *des);
void ccopy_mat_inc(UINT size, CTYPE  *X, ITER incx, CTYPE *Y,ITER incxy);
void omp_ccopy_mat(UINT N, CTYPE *src, SINT src_inc, CTYPE *des, SINT des_inc);

#if USE_CUDA
__global__ void cu_copy(DTYPE *SRC, UINT INC_SRC, DTYPE *DES, UINT INC_DES,
                        UINT len, UINT block_size);
__global__ void cu_ccopy(CTYPE *SRC, UINT INC_SRC, CTYPE *DES, UINT INC_DES,
                         UINT len, UINT block_size);
#endif

/* Sum magnitude of the MAT elements */
DTYPE asum(MAT *mat, UINT inc);
DTYPE omp_asum(UINT N, DTYPE *data, UINT inc);

/* Sum magnitude of the MAT elements *
 * (a + bi)'s magnitude = |a| + |b|  */
DTYPE casum(CMAT *mat, UINT inc);
DTYPE omp_casum(UINT N, CTYPE *data, UINT inc);

#if USE_CUDA
__global__ void cu_asum(DTYPE *data, UINT inc, UINT len, UINT block_size,
                        DTYPE *sum);
__global__ void cu_casum(CTYPE *data, UINT inc, UINT len, UINT block_size,
                         CTYPE *sum);
#endif

/* (real vector)-(real vector) dot product. */
DTYPE dot(MAT *src_x, UINT x_increment, MAT *src_y, UINT y_increment);
DTYPE omp_dot(UINT N, DTYPE *src_x, UINT x_inc, DTYPE *src_y, UINT y_inc);

/* (complex vector)-(real vector) dot product. */
CTYPE cdot(CMAT *src_x, UINT x_increment, MAT *src_y, UINT y_increment);
CTYPE omp_cdot(UINT N, CTYPE *src_x, UINT x_inc, DTYPE *src_y, UINT y_inc);

/* (complex vector)-(complex vector) dot product. */
CTYPE udot(CMAT *src_x, UINT x_increment, CMAT *src_y, UINT y_increment);
CTYPE omp_udot(UINT N, CTYPE *src_x, UINT x_inc, CTYPE *src_y, UINT y_inc);

#if USE_CUDA
__global__ void cu_dot(DTYPE *result, DTYPE *src_x, UINT x_inc, DTYPE *src_y,
                       UINT y_inc, UINT len, UINT block_size);
#endif

/* Swaps a real vector with another real vector. */
void swap(MAT *src_x, MAT *src_y);
void swap_inc(UINT N, DTYPE *src_x, UINT x_inc, DTYPE *src_y, UINT y_inc);
void omp_swap(UINT N, DTYPE *src_x, UINT x_inc, DTYPE *src_y, UINT y_inc);

/* Swaps a complex vector with another complex vector. */
void cswap(CMAT *src_x, CMAT *src_y);
void cswap_inc(UINT N, CTYPE *src_x, UINT x_inc, CTYPE *src_y, UINT y_inc);
void omp_cswap(UINT N, CTYPE *src_x, UINT x_inc, CTYPE *src_y, UINT y_inc);

/* column swap*/
void col_swap(MAT *mat, UINT a, UINT b);
void col_cswap(CMAT *mat, UINT a, UINT b);

/* row swap*/
void row_swap(MAT *mat, UINT a, UINT b);
void row_cswap(CMAT *mat, UINT a, UINT b);
/* Index of the maximum absolute value element of a real vector */
UINT amax(MAT *src);
UINT amax_inc(MAT *src, UINT inc);
UINT omp_amax(UINT N, DTYPE *src, UINT inc);

/* Index of the maximum absolute value element of a complex vector */
UINT camax(CMAT *src);
UINT camax_inc(CMAT *src, UINT inc);
UINT omp_camax(UINT N, CTYPE *src, UINT inc);

/* Index of the minimum absolute value element of a real vector */
UINT amin(MAT *src);
UINT amin_inc(MAT *src, UINT inc);
UINT omp_amin(UINT N, DTYPE *src, UINT inc);

/* Index of the minimum absolute value element of a complex vector */
UINT camin(CMAT *src);
UINT camin_inc(CMAT *src, UINT inc);
UINT omp_camin(UINT N, CTYPE *src, UINT inc);

/* Compute the absolute value of a complex number */
DTYPE cabs1(CTYPE val);

/* Generate Givens rotation of real points */
void rotg(DTYPE *a, DTYPE *b, DTYPE *c, DTYPE *s);
void omp_rotg(DTYPE *a, DTYPE *b, DTYPE *c, DTYPE *s);

/* Generate Givens rotation of complex points */
void crotg(CTYPE *a, CTYPE *b, DTYPE *c, CTYPE *s);
void omp_crotg(CTYPE *a, CTYPE *b, DTYPE *c, CTYPE *s);

/* Vector 2-norm (Euclidean norm)(Real Vector) */
DTYPE nrm2(MAT *src);
DTYPE nrm2_inc(MAT *src, UINT inc);
DTYPE omp_nrm2(UINT N, DTYPE *data, UINT inc);

/* Vector 2-norm (Euclidean norm)(Complex Vector) */
DTYPE cnrm2(CMAT *src);
DTYPE cnrm2_inc(CMAT *src, UINT inc);
DTYPE omp_cnrm2(UINT N, CTYPE *data, UINT inc);

/* Plane rotation of real points */
void rot(MAT *src_x, MAT *src_y, DTYPE c, DTYPE s);
void rot_inc(MAT *src_x, UINT x_inc, MAT *src_y, UINT y_inc, DTYPE c, DTYPE s);
void omp_rot(UINT N, DTYPE *src_x, UINT x_inc, DTYPE *src_y, UINT y_inc,
             DTYPE c, DTYPE s);

/* Plane rotation of complex points */
void crot(CMAT *src_x, CMAT *src_y, DTYPE c, DTYPE s);
void crot_inc(CMAT *src_x, UINT x_inc, CMAT *src_y, UINT y_inc, DTYPE c,
              DTYPE s);
void omp_crot(UINT N, CTYPE *src_x, UINT x_inc, CTYPE *src_y, UINT y_inc,
              DTYPE c, DTYPE s);

/* (real MAT) * (real number) */
void scal(DTYPE alpha, MAT *mat);
void scal_inc(UINT size, DTYPE alpha, DTYPE *X, UINT incx);
void omp_scal(UINT size, DTYPE alpha, DTYPE *X, UINT incx);

/* (complex MAT) * (real number) */
void cscal(DTYPE alpha, CMAT *mat);
void cscal_inc(UINT size, DTYPE alpha, CTYPE *X, UINT incx);
void omp_cscal(UINT size, DTYPE alpha, CTYPE *X, UINT incx);

/* (complex MAT) * (complex number) */
void uscal(CTYPE alpha, CMAT *mat);
void uscal_inc(UINT size, CTYPE alpha, CTYPE *X, UINT incx);
void omp_uscal(UINT size, CTYPE alpha, CTYPE *X, UINT incx);

/* column scaling */
void col_scal(DTYPE alpha, MAT *X, UINT idx);
void col_cscal(DTYPE alpha, CMAT *X, UINT idx);
void col_uscal(CTYPE alpha, CMAT *X, UINT idx);

/* row scaling */
void row_scal(DTYPE alpha, MAT *X, UINT idx);
void row_cscal(DTYPE alpha, CMAT *X, UINT idx);
void row_uscal(CTYPE alpha, CMAT *X, UINT idx);

/* (real MAT) + (real number) */
void add(DTYPE alpha, MAT *mat);
void add_inc(UINT size, DTYPE alpha, DTYPE *X, UINT incx);

/* (complex MAT) + (real number) */
void cadd(DTYPE alpha, CMAT *mat);
void cadd_inc(UINT size, DTYPE alpha, CTYPE *X, UINT incx);

/* (complex MAT) + (complex number) */
void uadd(CTYPE alpha, CMAT *mat);
void uadd_inc(UINT size, CTYPE alpha, CTYPE *X, UINT incx);

/* column add */
void col_add(DTYPE alpha, MAT *X, UINT idx);
void col_cadd(DTYPE alpha, CMAT *X, UINT idx);
void col_uadd(CTYPE alpha, CMAT *X, UINT idx);

/* row add */
void row_add(DTYPE alpha, MAT *X, UINT idx);
void row_cadd(DTYPE alpha, CMAT *X, UINT idx);
void row_uadd(CTYPE alpha, CMAT *X, UINT idx);

#endif
