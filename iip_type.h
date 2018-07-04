#ifndef IIP_TYPE_H
#define IIP_TYPE_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>

#define DEBUG 0

/***********************************
* 이 부분은 직접 해주세요
************************************/

#define DTYPE double
/*
* If DTYPE = float  ->  set NTYPE = 0
* If DTYPE = double ->  set NTYPE = 1
*
* */
#define NTYPE 1
/************************************
*여기까지 설정해주세요
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
오버로딩
#define 오버로딩매크로(_x(함수중 가장 인자가 적은수의 인자수만큼), 함수 수만큼, ...)맨 오른쪽 것
#define 오버로딩할함수(...) 오버로딩매크로(__VA_ARGS__, 함수들)(__VA_ARGS__)

__VA_ARGS__는 함수의 인자의 매크로

EX)

void f1(arg1,arg2)
void f2(arg1,arg2,arg3)

인자가 가장적은 함수의 인자가 2개이기 때문에
_x,_xx
함수가 2개니까
_2,_1

#define o_f(_x,_xx, _2, _1) _1
#define f(...) o_f(__VA_ARGS__, f2, f1)(__VA_ARGS__)

*/

/***************************
**** pre main for CUDA ****
***************************/

/*
//This only works in GNU
void __attribute__ ((constructor)) premain()
{


}

*/

#if USE_CBLAS 
#include "cblas.h"
#endif
#if BLAS_CU
#include "cublas_v2.h"
#endif



#endif