#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define DEBUG 1

#if USE_CBLAS 

#include "cblas.h"
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


#define CTYPE double
/*
 * If CTYPE = float  ->  set NTCYPE = 0
 * If CTYPE = double ->  set NCTYPE = 1
 *
 * */
#define NCTYPE 1

/************************************
 *여기까지 설정해주세요
 *********************************** */
#define UINT int
#define ITER long

#define str(x) #x
#define xstr(x) str(x)

typedef struct MAT
{
	DTYPE* data;
	/*
	 *	if ndim = 0  --> 1D
	 *			= 1  --> 2D
	 *			= 2  --> 3D
	 * */
	UINT ndim;
	UINT d0;
	UINT d1;
	UINT d2;
}MAT;

typedef struct COMPLEX
{
	CTYPE re;
	CTYPE im;
}COMPLEX;

typedef struct CMAT
{
	COMPLEX* data;
	/*
	 *	if ndim = 0  --> 1D
	 *			= 1  --> 2D
	 *			= 2  --> 3D
	 * */
	UINT ndim;
	UINT d0;
	UINT d1;
	UINT d2;
}CMAT;

/********************
 **** iip_math.c ****
 ********************/

/**** zeros overloading  ****/
#define zeros_load(_x,_3,_2,_1,...) _1
#define zeros(...) zeros_load(__VA_ARGS__, zeros_3d,zeros_2d,zeros_1d)(__VA_ARGS__) 

MAT* zeros_1d(UINT);
MAT* zeros_2d(UINT,UINT);
MAT* zeros_3d(UINT,UINT,UINT);

/**** set overloading  ****/
#define set_load(_x,_xx,_xxx,_3,_2,_1,...) _1
#define set(...) set_load(__VA_ARGS__, set_3d,set_2d,set_1d)(__VA_ARGS__)
void set_1d(MAT*,UINT,UINT);
void set_2d(MAT*,UINT,UINT,UINT);
void set_3d(MAT*,UINT,UINT,UINT,UINT);

/****  get overloadnig ****/
#define get_load(_x,_xx,_3,_2,_1,...) _1
#define get(...) get_load(__VA_ARGS__, get_3d,get_2d,get_1d)(__VA_ARGS__)
DTYPE get_1d(MAT*,UINT);
DTYPE get_2d(MAT*,UINT,UINT);
DTYPE get_3d(MAT*,UINT,UINT,UINT);

/****  miscellaneous ****/ 
void free_MAT(MAT*);
void print_MAT(MAT*);

/********************
 **** iip_blas_lv1.c ****
 ********************/
void mp_axpy(UINT,DTYPE,DTYPE*,UINT,DTYPE*,UINT);
void axpy(DTYPE,MAT*,MAT*);

