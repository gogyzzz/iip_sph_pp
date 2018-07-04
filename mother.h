#ifndef MOTHER_H
#define MOTHER_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>

#include "iip_blas_lv1.h"
#include "iip_matrix.h"

#if USE_CUDA
extern cublasHandle_t handle;

class cublas_setting
{
	public:
	cublas_setting(){cublasCreate(&handle);}
	~cublas_setting(){cublasDestroy(&handle);}
}

static cublas_setting;
#endif



#endif



#endif
