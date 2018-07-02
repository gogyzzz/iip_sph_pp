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

typedef struct CTYPE
{
	DTYPE re;
	DTYPE im;
}CTYPE;

typedef struct CMAT
{
	CTYPE* data;
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
void set_1d(MAT*,UINT,DTYPE);
void set_2d(MAT*,UINT,UINT,DTYPE);
void set_3d(MAT*,UINT,UINT,UINT,DTYPE);

void fill(MAT*, DTYPE);

/****  get overloadnig ****/
#define get_load(_x,_xx,_3,_2,_1,...) _1
#define get(...) get_load(__VA_ARGS__, get_3d,get_2d,get_1d)(__VA_ARGS__)
DTYPE get_1d(MAT*,UINT);
DTYPE get_2d(MAT*,UINT,UINT);
DTYPE get_3d(MAT*,UINT,UINT,UINT);

#define submat_load(_x,_xx,_7,_6,_5,_4,_3,_2,_1,...) _1
#define submat(...) submat_load(__VA_ARGS__, submat_3d,_,submat_2d,_,submat_1d)(__VA_ARGS__)
void submat_1d(MAT*, MAT*, int,int);
void submat_2d(MAT*, MAT*, int,int, int,int);
void submat_3d(MAT*, MAT*, int,int, int,int, int,int);

/****  miscellaneous ****/ 
void free_MAT(MAT*);
void print_MAT(MAT*);

/********************
 **** iip_blas_lv1.c ****
 ********************/
void mp_axpy(UINT,DTYPE,DTYPE*,UINT,DTYPE*,UINT);
void axpy(DTYPE,MAT*,MAT*);

