#ifndef BLAS_LV3_H
#define BLAS_LV3_H
#include "iip_type.h"

	#ifndef USE_CUDA
void gemm(char, char, DTYPE, MAT*,MAT*,DTYPE,MAT*);
void mp_gemm(char, char, UINT,UINT,UINT,DTYPE, DTYPE*,UINT,DTYPE*,UINT,DTYPE,DTYPE*,UINT);

void cgemm(char, char, CTYPE, CMAT*,CMAT*,CTYPE,CMAT*);
void mp_cgemm(char, char, UINT,UINT,UINT,CTYPE, CTYPE*,UINT,CTYPE*,UINT,CTYPE,CTYPE*,UINT);
	#else





	#endif




#endif
