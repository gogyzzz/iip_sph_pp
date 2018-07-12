#ifndef BLAS_LV1_H
#define BLAS_LV1_H

//#include "mother.h"
#include "iip_type.h"

/*********************
**** iip_blas_lv1 ****
**********************/

/* TODO
 * + common
 * asum
 * nrm2
 * rot
 * scal
 * swap
 * i?amax
 * i?amin
 * -------copy
 * -------axpy
 * + s,d
 * dot
 * rotm
 * rotmg
 * rotg
 * + c,z
 * dotc
 * dotu
 *
 * */
void mp_axpy(UINT,DTYPE,DTYPE*,UINT,DTYPE*,UINT);
void axpy(DTYPE,MAT*,MAT*);

void mp_caxpy(UINT,CTYPE,CTYPE*,UINT,CTYPE*,UINT);
void caxpy(CTYPE,CMAT*,CMAT*);

void mp_copy(UINT N, DTYPE* src, SINT src_inc, DTYPE* des, SINT des_inc);
void copy(MAT* src, MAT* des);

void mp_ccopy(UINT N, CTYPE* src, SINT src_inc, CTYPE* des, SINT des_inc);
void ccopy(CMAT* src, CMAT* des);

DTYPE mp_asum(UINT N, DTYPE* data, UINT inc);
DTYPE asum(UINT N, MAT *mat, UINT inc);

DTYPE mp_casum(UINT N, CTYPE* data, UINT inc);
DTYPE casum(UINT N, CMAT *mat, UINT inc);

DTYPE mp_dot(UINT N, DTYPE* src_x, UINT x_inc, DTYPE* src_y, UINT y_inc);
DTYPE dot(UINT N, MAT* src_x, UINT x_increment, MAT *src_y, UINT y_increment);

DTYPE mp_sdot(UINT N, DTYPE* src_x, UINT x_inc, DTYPE* src_y, UINT y_inc);
DTYPE sdot(UINT N, MAT* src_x, UINT x_increment, MAT *src_y, UINT y_increment);

DTYPE mp_cdot(UINT N, CTYPE* src_x, UINT x_inc, CTYPE* src_y, UINT y_inc);
DTYPE cdot(UINT N, CMAT* src_x, UINT x_increment, CMAT *src_y, UINT y_increment);

#endif
