#ifndef IIP_TYPE_H
#define IIP_TYPE_H

#include <complex.h>
#include <math.h>
#include <memory.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
// random
#include <time.h>
// DBL_MAX, FLT_MAX
#include <float.h>

#if USE_CUDA
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#endif

#if OS_WIN
#include <windows.h>
#define __func__ __FUNCTION__
#endif
#define DEBUG 0

/***********************************
* 이 부분은 직접 해주세요
*********************************** */
#define DTYPE double
/*
* If DTYPE = float  ->  set NTYPE = 0
* If DTYPE = double ->  set NTYPE = 1
* */
#define NTYPE 1
/************************************
*여기까지 설정해주세요
*********************************** */

#define MAX_CHAR 256
#define FZERO 1e-10

#if NTYPE == 0
#define INF FLT_MAX
#elif NTYPE == 1
#define INF DBL_MAX
#endif

#define UINT uint32_t
#define SINT int32_t
#define ITER long

#if USE_CUDA
#define NoTran CUBLAS_OP_N
#define Tran CUBLAS_OP_T
#define CTran CUBLAS_OP_C
#else
#define NoTran 111
#define Tran 112
#define CTran 113
#endif

/*
* #define AA BB
* str(AA) -> 'AA'
* xstr(AA) -> 'BB'
* */
#define str(x) #x
#define xstr(x) str(x)

/**** STRUCT ****/
typedef struct MAT {
  DTYPE* data;
  UINT ndim;
  UINT d0;
  UINT d1;
  UINT d2;
} MAT;

typedef struct CTYPE {
  DTYPE re;
  DTYPE im;
} CTYPE;

typedef struct CMAT {
  CTYPE* data;
  UINT ndim;
  UINT d0;
  UINT d1;
  UINT d2;
} CMAT;

typedef struct RANGE {
  UINT s0, e0;  // d0 range
  UINT s1, e1;  // d1 range
  UINT s2, e2;  // d2 range
} RANGE;

typedef struct DIM {
  UINT d0;
  UINT d1;
  UINT d2;
} DIM;

/*
오버로딩
#define 오버로딩매크로(_x(함수중 가장 인자가 적은수의 인자수만큼), 함수 수만큼,
...)맨 오른쪽 것
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

/**** LIBRARY SETTING ****/

#if USE_OPEN
#include "cblas.h"
#endif

#if USE_MKL
#include "mkl.h"
#endif

#if USE_CUDA
#include "cublas_v2.h"
extern cublasHandle_t handle;
extern UINT max_thread;
extern UINT max_block;

// CAST CTYPE POINTER TO CUDA_COMPLEX TYPE POINTER
#if NTYPE == 0
#define CU_CX(x) (cuComplex*)(void*)(x)

#elif NTYPE == 1
#define CU_CX(x) (cuDoubleComplex*)(void*)(x)

#endif
#endif

/*******************************
**** MACRO FOR TPYE CASITNG*****
********************************/
#if OS_UNIX
#define CXF(X) (*(complex float*)(&X))
#define CXD(X) (*(complex double*)(&X))
#elif OS_WIN
#define CXF(X) (*(_Fcomplex*)(&X))
#define CXD(X) (*(_Dcomplex*)(&X))
#endif

/*************************************
 **** MACRO FUNCTIONS ****
 *************************************/
// Y*=X
#define cxmul(Y, X, T)                \
  {                                   \
    T = Y.re;                         \
    Y.re = Y.re * X.re - Y.im * X.im; \
    Y.im = T * X.im + Y.im * X.re;    \
  }
// Y+=X
#define cxadd(Y, X)     \
  {                     \
    Y.re = Y.re + X.re; \
    Y.im = Y.im + X.im; \
  }
// Y+=A*B
#define cxadd_mul(Y, A, B)                               \
  {                                                      \
    (Y.re) = (Y.re) + (A.re) * (B.re) - (A.im) * (B.im); \
    (Y.im) = (Y.im) + (A.re) * (B.im) + (A.im) * (B.re); \
  }
// Y = A*B
#define cxemul(Y, A, B)                       \
  {                                           \
    Y.re = (A.re) * (B.re) - (A.im) * (B.im); \
    Y.im = (A.re) * (B.im) + (A.im) + (B.re); \
  }

// Y = A/B
#define cxediv(Y, A, B)                             \
  {                                                 \
    Y.re = ((A.re) * (B.re) + (A.im) * (B.im)) /    \
           (((B.re) * (B.re)) + ((B.im) * (B.im))); \
    Y.im = ((A.im) * (B.re) - (A.re) * (B.im)) /    \
           (((B.re) * (B.re)) + ((B.im) * (B.im))); \
  }

#define cxeadd(Y, A, B) \
  {                     \
    Y.re = A.re + B.re; \
    Y.im = A.im + B.im; \
  }

#define SWAP(x, y, t) \
  {                   \
    t = y;            \
    y = x;            \
    x = t;            \
  }

/**** ASSERT MACRO ****/
#define DIM_INVAL 1
#define ARG_INVAL 2
#define NO_FILE 3
#define NOT_SQUARE 4
#define DET_FAIL 5

/**** ASSERT FUNCTION ****/
#define ASSERT(x)                                                             \
  {                                                                           \
    if (x == 1) {                                                             \
      printf("WARNING(%s) : invalid dimension\n", __func__);                  \
      return;                                                                 \
    } else if (x == 2) {                                                      \
      printf("WARNING(%s) : invalid argument\n", __func__);                   \
      return;                                                                 \
    } else if (x == 3) {                                                      \
      printf("ERROR(%s)   : no such file exits \n", __func__);                \
      exit(-1);                                                               \
    } else if (x == 4) {                                                      \
      printf("WARNING(%s) : This function requires SQAURE MATRIX \n",         \
             __func__);                                                       \
      return;                                                                 \
    } else if (x == 5) {                                                      \
      printf("WARNING(%s) : Matrix is singular, |det| < FEZRO \n", __func__); \
      return;                                                                 \
    }                                                                         \
  }

/*****************************
 **** MEMORY MANAGER *********
 *****************************/
#define MAX_MEM_PAGE 20
#define MEM_PAGE_BASIC_SIZE 256
#define LOG_ALLOC_UNIT 8  // 4 for float, 8 for double

typedef struct USED_MEM {
  unsigned long long int frag_idx;
  unsigned long long int size;
} USED_MEM;

//메모리풀
static void* memory_pool[MAX_MEM_PAGE];
//
static USED_MEM* memory_log[MAX_MEM_PAGE];

//해당 페이지에 몇 덩어리가 할당 되어있나
static unsigned long long int log_cnt[MAX_MEM_PAGE];
static unsigned int pool_cnt;
signed long long int page_alloc_isable(int page_idx, unsigned long long int require_size);
void* iip_malloc(unsigned long long int size);
void iip_free(void* ptr);
void init(UINT mem_pool_size);
void finit();
#endif
