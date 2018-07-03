#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>

#define DEBUG 0

#if USE_CBLAS 
#include "cblas.h"
#endif
#if BLAS_CU
#include "cublas_v2.h"
#endif


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


/********************
 **** iip_matrix ****
 ********************/

void fill(MAT*, DTYPE);
void cfill(CMAT*, DTYPE,DTYPE);

/*** allocMAT ***/
#define alloc_MAT_load(_x,_3,_2,_1,...) _1
#define alloc_MAT(...) allocMAT_load(__VA_ARGS__, alloc_MAT_3d,alloc_MAT_2d,alloc_MAT_1d)(__VA_ARGS__) 

MAT* alloc_MAT_1d(UINT);
MAT* alloc_MAT_2d(UINT,UINT);
MAT* alloc_MAT_3d(UINT,UINT,UINT);

#define calloc_MAT_load(_x,_3,_2,_1,...) _1
#define calloc_MAT(...) callocMAT_load(__VA_ARGS__, calloc_MAT_3d,calloc_MAT_2d,calloc_MAT_1d)(__VA_ARGS__) 
CMAT* calloc_MAT_1d(UINT);
CMAT* calloc_MAT_2d(UINT,UINT);
CMAT* calloc_MAT_3d(UINT,UINT,UINT);


/**** zeros  ****/
#define zeros_load(_x,_3,_2,_1,...) _1
#define zeros(...) zeros_load(__VA_ARGS__, zeros_3d,zeros_2d,zeros_1d)(__VA_ARGS__) 

MAT* zeros_1d(UINT);
MAT* zeros_2d(UINT,UINT);
MAT* zeros_3d(UINT,UINT,UINT);

#define czeros_load(_x,_3,_2,_1,...) _1
#define czeros(...) czeros_load(__VA_ARGS__, czeros_3d,czeros_2d,czeros_1d)(__VA_ARGS__) 
CMAT* czeros_1d(UINT);
CMAT* czeros_2d(UINT,UINT);
CMAT* czeros_3d(UINT,UINT,UINT);

/**** set overloading  ****/
#define set_load(_x,_xx,_xxx,_3,_2,_1,...) _1
#define set(...) set_load(__VA_ARGS__, set_3d,set_2d,set_1d)(__VA_ARGS__)
void set_1d(MAT*,UINT,DTYPE);
void set_2d(MAT*,UINT,UINT,DTYPE);
void set_3d(MAT*,UINT,UINT,UINT,DTYPE);

#define cset_load(_x,_xx,_xxx,_xxxx,_3,_2,_1,...) _1
#define cset(...) cset_load(__VA_ARGS__, cset_3d,cset_2d,cset_1d)(__VA_ARGS__)
void cset_1d(CMAT*,UINT,DTYPE,DTYPE);
void cset_2d(CMAT*,UINT,UINT,DTYPE,DTYPE);
void cset_3d(CMAT*,UINT,UINT,UINT,DTYPE,DTYPE);

#if USE_CUDA
__global__ void cu_set(DTYPE*,UINT,UINT,UINT,DTYPE);
__global__ void cu_set(CTYPE*,UINT,UINT,UINT,DTYPE,DTYPE);
#endif


/****  get overloadnig ****/
#define get_load(_x,_xx,_3,_2,_1,...) _1
#define get(...) get_load(__VA_ARGS__, get_3d,get_2d,get_1d)(__VA_ARGS__)
DTYPE get_1d(MAT*,UINT);
DTYPE get_2d(MAT*,UINT,UINT);
DTYPE get_3d(MAT*,UINT,UINT,UINT);

#define cget_load(_x,_xx,_3,_2,_1,...) _1
#define cget(...) cget_load(__VA_ARGS__, cget_3d,cget_2d,cget_1d)(__VA_ARGS__)
CTYPE cget_1d(CMAT*,UINT);
CTYPE cget_2d(CMAT*,UINT,UINT);
CTYPE cget_3d(CMAT*,UINT,UINT,UINT);

/**** submat overloading ****/

//since arg +=2, pend _x for each function
#define submat_load(_x1,_x2,_x3,_x4,_3,_x5,_2,_x6,_1,...) _1
#define submat(...) submat_load(__VA_ARGS__, submat_3d,_,submat_2d,_,submat_1d)(__VA_ARGS__)
void submat_1d(MAT*, MAT*, ITER,ITER);
void submat_2d(MAT*, MAT*, ITER,ITER, ITER,ITER);
void submat_3d(MAT*, MAT*, ITER,ITER, ITER,ITER, ITER,ITER);

#define csubmat_load(_x1,_x2,_x3,_x4,_3,_x5,_2,_x6,_1,...) _1
#define csubmat(...) csubmat_load(__VA_ARGS__, csubmat_3d,_,csubmat_2d,_,csubmat_1d)(__VA_ARGS__)
void csubmat_1d(CMAT*, CMAT*, ITER,ITER);
void csubmat_2d(CMAT*, CMAT*, ITER,ITER, ITER,ITER);
void csubmat_3d(CMAT*, CMAT*, ITER,ITER, ITER,ITER, ITER,ITER);

/****  miscellaneous ****/ 
void free_MAT(MAT*);
void free_CMAT(CMAT*);
void print_MAT(MAT*);
void print_CMAT(CMAT*);

/********************
 **** iip_blas_lv1 ****
 ********************/

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



/***************************
 **** pre main for CUDA ****
 ***************************/

//This only works in GNU
void __attribute__ ((constructor)) premain()
{


}

