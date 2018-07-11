#ifndef BLAS_LV2_H
#define BLAS_LV2_H
#include "iip_type.h"

/*
 * TODO
 *
 * gbmv
 *
 * */

void gemv(int,DTYPE,MAT*,MAT*,DTYPE,MAT*);
void mp_gemv(int32_t,UINT,UINT,DTYPE,DTYPE*,UINT,DTYPE*,SINT,DTYPE,DTYPE*,SINT );


void cgemv();
void mp_cgemv();


#endif
