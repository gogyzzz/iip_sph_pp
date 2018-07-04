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
 * copy
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

void mp_caxpy(UINT,DTYPE,CTYPE*,UINT,CTYPE*,UINT);
void caxpy(DTYPE,CMAT*,CMAT*);

void mp_copy(UINT N, DTYPE* src, SINT src_inc, DTYPE* des, SINT des_inc);
void copy(MAT* src, SINT src_increment, MAT* des, SINT des_increment);

void mp_ccopy(UINT N, DTYPE* src, SINT src_inc, DTYPE* des, SINT des_inc);
void ccopy(CMAT* src, SINT src_increment, CMAT* des, SINT des_increment);


#endif