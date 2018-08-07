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
#ifndef IIP_MATH_H
#define IIP_MATH_H

#include "iip_time.h"
#include "iip_type.h"

/* GENERAL RULE
 * {func} :  DTYPE matrix MAT operation
 * c{func} : {func} with CTYPE matrix CMAT 
 * {func}_inc : {func} with DTYPE array with increment and size
 * c{func}_inc : c{func} with CTYPE array 
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

/**** SQUARE ROOT ****/
void sqrt_mat(MAT* mat);
void sqrt_mat_inc(UINT size, DTYPE* X, ITER incx);

void sqrt_cmat(CMAT* mat);
void sqrt_cmat_inc(UINT size, CTYPE* X, ITER incx);

/**** POWER ****/
void pow_mat(MAT* mat, DTYPE n);
void pow_mat_inc(UINT size, DTYPE* X, DTYPE n, ITER incx);

void pow_cmat(CMAT* mat, DTYPE n);
void pow_cmat_inc(UINT size, CTYPE* X, DTYPE n, ITER incx);

void cpow_cmat(CMAT* mat, CTYPE n);
void cpow_cmat_inc(UINT size, CTYPE* X, CTYPE n, ITER incx);

/* Uniform distribution of random number */
void randu(MAT* mat, DTYPE a, DTYPE b);
void randu_inc(UINT size, DTYPE* X, DTYPE a, DTYPE b, ITER incx);

void crandu(CMAT* mat, DTYPE ra, DTYPE rb, DTYPE ia, DTYPE ib);
void crandu_inc(UINT size, CTYPE* X, DTYPE ra, DTYPE rb, DTYPE ia, DTYPE ib,
                ITER incx);

/* Generate Normal Distribution Matrix Using Box-Muller Transform */
void randn(MAT* mat, DTYPE mean, DTYPE std);
void randn_inc(UINT size, DTYPE* X, DTYPE mean, DTYPE std, ITER incx);

void crandn(CMAT* mat, CTYPE mean, CTYPE var);
void crandn_inc(UINT size, CTYPE* X, CTYPE mean, CTYPE std, ITER incx);

/**** round ****/
void round_mat(MAT* mat);
void round_mat_inc(UINT size, DTYPE* X, ITER incx);

void round_cmat(CMAT* mat);
void round_cmat_inc(UINT size, CTYPE* X, ITER incx);

/**** floor ****/
void floor_mat(MAT* mat);
void floor_mat_inc(UINT size, DTYPE* X, ITER incx);

void floor_cmat(CMAT* mat);
void floor_cmat_inc(UINT size, CTYPE* X, ITER incx);

/**** ceil ****/
void ceil_mat(MAT* mat);
void ceil_mat_inc(UINT size, DTYPE* X, ITER incx);

void ceil_cmat(CMAT* mat);
void ceil_cmat_inc(UINT size, CTYPE* X, ITER incx);

/**** log ****/
void log_mat(MAT* mat);
void log_mat_inc(UINT size, DTYPE* X, ITER incx);

void log_cmat(CMAT* mat);
void log_cmat_inc(UINT size, CTYPE* X, ITER incx);

void log10_mat(MAT* mat);
void log10_mat_inc(UINT size, DTYPE* X, ITER incx);

void log10_cmat(CMAT* mat);
void log10_cmat_inc(UINT size, CTYPE* X, ITER incx);

void log2_mat(MAT* mat);
void log2_mat_inc(UINT size, DTYPE* X, ITER incx);

void log2_cmat(CMAT* mat);
void log2_cmat_inc(UINT size, CTYPE* X, ITER incx);

/**** exp ****/
void exp_mat(MAT* mat);
void exp_mat_inc(UINT size, DTYPE* X, ITER incx);

void exp_cmat(CMAT* mat);
void exp_cmat_inc(UINT size, CTYPE* X, ITER incx);

/**** abs ****/
void abs_mat(MAT* mat);
void abs_mat_inc(UINT size, DTYPE* X, ITER incx);

void abs_cmat(CMAT* mat);
void abs_cmat_inc(UINT size, CTYPE* X, ITER incx);

/**** max ****/
DTYPE max_mat(MAT* mat, DIM* dim);
CTYPE max_cmat(CMAT* mat, DIM* dim);

DTYPE amax_mat(MAT* mat, DIM* dim);
DTYPE amax_cmat(CMAT* mat, DIM* dim);
/**** min ****/
DTYPE min_mat(MAT* mat, DIM* dim);
CTYPE min_cmat(CMAT* mat, DIM* dim);

DTYPE amin_mat(MAT* mat, DIM* dim);
DTYPE amin_cmat(CMAT* mat, DIM* dim);

/**** return complex type : r + ii ****/
CTYPE CX(DTYPE r, DTYPE i);
#endif
