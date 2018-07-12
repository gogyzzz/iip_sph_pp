#ifndef BLAS_LV2_H
#define BLAS_LV2_H
#include "iip_type.h"

/*
 * TODO
 *
 * gemv
 *
 * */

#ifndef USE_CUDA
void gemv(char,DTYPE,MAT*,MAT*,DTYPE,MAT*);
void mp_gemv(char,UINT,UINT,DTYPE,DTYPE*,UINT,DTYPE*,SINT,DTYPE,DTYPE*,SINT );

void cgemv(char, CTYPE,CMAT*,CMAT*,CTYPE,CMAT*);
void mp_cgemv(char, UINT,UINT,CTYPE,CTYPE*,UINT,CTYPE*,SINT,CTYPE,CTYPE*,SINT);

#else
void gemv(cublasOperation_t,DTYPE,MAT*,MAT*,DTYPE,MAT*);
void mp_gemv(cublasOperation_t,UINT,UINT,DTYPE,DTYPE*,UINT,DTYPE*,SINT,DTYPE,DTYPE*,SINT );

void cgemv(cublasOperation_t, CTYPE,CMAT*,CMAT*,CTYPE,CMAT*);
void mp_cgemv(cublasOperation_t, UINT,UINT,CTYPE,CTYPE*,UINT,CTYPE*,SINT,CTYPE,CTYPE*,SINT);
#endif

#endif
