#ifndef BLAS_LV2_H
#define BLAS_LV2_H
#include "iip_type.h"

/*
 * TODO
 *
 * gemv
 *
 * */

void gemv(int32_t,DTYPE,MAT*,MAT*,DTYPE,MAT*);
void mp_gemv(int32_t,UINT,UINT,DTYPE,DTYPE*,UINT,DTYPE*,SINT,DTYPE,DTYPE*,SINT );


void cgemv(int32_t, CTYPE,CMAT*,CMAT*,CTYPE,CMAT*);
void mp_cgemv(int32_t, UINT,UINT,CTYPE,CTYPE*,UINT,CTYPE*,SINT,CTYPE,CTYPE*,SINT);


#endif
