#ifndef IIP_TYPE_H
#define IIP_TYPE_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>

#define DEBUG 1

/***********************************
* ÀÌ ºÎºÐÀº Á÷Á¢ ÇØÁÖ¼¼¿ä
************************************/

#define DTYPE double
/*
* If DTYPE = float  ->  set NTYPE = 0
* If DTYPE = double ->  set NTYPE = 1
*
* */
#define NTYPE 1
/************************************
*¿©±â±îÁö ¼³Á¤ÇØÁÖ¼¼¿ä
*********************************** */

/*
INDEX
iip_matrix
iip_blas_lv1

*/

#define UINT uint32_t
#define SINT int32_t
#define ITER long

/*
* #define AA BB
* str(AA) -> 'AA'
* xstr(AA) -> 'BB'
* */
#define str(x) #x
#define xstr(x) str(x)

typedef struct MAT
{
	DTYPE* data;
	UINT ndim;
	UINT d0;
	UINT d1;
	UINT d2;
}MAT;

typedef struct CTYPE
{
	DTYPE re;
	DTYPE im;
}CTYPE;

typedef struct CMAT
{
	CTYPE* data;
	UINT ndim;
	UINT d0;
	UINT d1;
	UINT d2;
}CMAT;


/*
¿À¹ö·Îµù
#define ¿À¹ö·Îµù¸ÅÅ©·Î(_x(ÇÔ¼öÁß °¡Àå ÀÎÀÚ°¡ ÀûÀº¼öÀÇ ÀÎÀÚ¼ö¸¸Å­), ÇÔ¼ö ¼ö¸¸Å­, ...)¸Ç ¿À¸¥ÂÊ °Í
#define ¿À¹ö·ÎµùÇÒÇÔ¼ö(...) ¿À¹ö·Îµù¸ÅÅ©·Î(__VA_ARGS__, ÇÔ¼öµé)(__VA_ARGS__)

__VA_ARGS__´Â ÇÔ¼öÀÇ ÀÎÀÚÀÇ ¸ÅÅ©·Î

EX)

void f1(arg1,arg2)
void f2(arg1,arg2,arg3)

ÀÎÀÚ°¡ °¡ÀåÀûÀº ÇÔ¼öÀÇ ÀÎÀÚ°¡ 2°³ÀÌ±â ¶§¹®¿¡
_x,_xx
ÇÔ¼ö°¡ 2°³´Ï±î
_2,_1

#define o_f(_x,_xx, _2, _1) _1
#define f(...) o_f(__VA_ARGS__, f2, f1)(__VA_ARGS__)

*/

/***************************
**** pre main for CUDA ****
***************************/

#if USE_CBLAS 
#include "cblas.h"
#endif
#if USE_CUDA
#include "cublas_v2.h"
#endif

static cublasHandle_t handle;
static UINT max_thread;
static UINT max_block;


#endif
