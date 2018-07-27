#ifndef IIP_MATH_H
#define IIP_MATH_H

#include "iip_time.h"
#include "iip_type.h"

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

// Uniform distribution of random number
void randu(MAT* mat, DTYPE a, DTYPE b);
void randu_inc(UINT size, DTYPE* X, DTYPE a, DTYPE b, ITER incx);

void crandu(CMAT* mat, DTYPE ra, DTYPE rb, DTYPE ia, DTYPE ib);
void crandu_inc(UINT size, CTYPE* X, DTYPE ra, DTYPE rb, DTYPE ia, DTYPE ib,
                ITER incx);

// Generate Normal Distribution Matrix Using Box-Muller Transform
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
CTYPE amax_cmat(CMAT* mat, DIM* dim);
/**** min ****/
DTYPE min_mat(MAT* mat, DIM* dim);
CTYPE min_cmat(CMAT* mat, DIM* dim);

DTYPE amin_mat(MAT* mat, DIM* dim);
CTYPE amin_cmat(CMAT* mat, DIM* dim);

#endif
