/*
 * ===========================================================
 *           Copyright (c) 2018, __IIPLAB__
 *                All rights reserved.
 *
 * This Source Code Form is subject to the terms of
 * the Mozilla Public License, v. 2.0.
 * If a copy of the MPL was not distributed with this file,
 *  You can obtain one at http://mozilla.org/MPL/2.0/.
 * ===========================================================
 */
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

/*If set by 1, every function call will print its name*/
#define DEBUG 1

/***********************************
* Set This Section Manually
*********************************** */
#define DTYPE double
/*
* If DTYPE = float  ->  set NTYPE = 0
* If DTYPE = double ->  set NTYPE = 1
* */
#define NTYPE 1
/************************************
*********************************** */

#define MAX_CHAR 256
#define FZERO 1e-0

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

/**** LIBRARY SETTING ****/

#if USE_CUDA
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#endif

#if OS_WIN
#include <windows.h>
#endif

#if USE_OPEN
#include "cblas.h"
#endif

#if USE_MKL
#include "mkl.h"
#endif

#if USE_CBLAS
#include "lapacke.h"
#endif

#if USE_CUDA
#include "cublas_v2.h"
extern cublasHandle_t handle;
extern UINT max_thread;
extern UINT max_block;

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

/**** MACRO FOR MACRO ****/
#if OS_WIN
#define __func__ __FUNCTION__
#endif
/*
   overloading

#define overloading_macro(_x(the number of argument of function with least
arguments), many as the number of functions,
...) rightest one
#define function_to_overload(...) overloading_macro(__VA_ARGS__,
functions)(__VA_ARGS__)

__VA_ARGS__is MACRO for arguments of function

EX)

void f1(arg1,arg2)
void f2(arg1,arg2,arg3)

function with leasts has 2 arguments
_x,_xx
2 function to overload
_2,_1

#define o_f(_x,_xx, _2, _1) _1
#define f(...) o_f(__VA_ARGS__, f2, f1)(__VA_ARGS__)

*/

// CAST CTYPE POINTER TO CUDA_COMPLEX TYPE POINTER
#if NTYPE == 0
#define CU_CX(x) (cuComplex*)(void*)(x)

#elif NTYPE == 1
#define CU_CX(x) (cuDoubleComplex*)(void*)(x)

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
#define CXMUL(Y, X, T)                \
  {                                   \
    T = Y.re;                         \
    Y.re = Y.re * X.re - Y.im * X.im; \
    Y.im = T * X.im + Y.im * X.re;    \
  }
// Y+=X
#define CXADD(Y, X)     \
  {                     \
    Y.re = Y.re + X.re; \
    Y.im = Y.im + X.im; \
  }
// Y+=A*B
#define CXADD_mul(Y, A, B)                               \
  {                                                      \
    (Y.re) = (Y.re) + (A.re) * (B.re) - (A.im) * (B.im); \
    (Y.im) = (Y.im) + (A.re) * (B.im) + (A.im) * (B.re); \
  }
// Y = A*B
#define CXEMUL(Y, A, B)                       \
  {                                           \
    Y.re = (A.re) * (B.re) - (A.im) * (B.im); \
    Y.im = (A.re) * (B.im) + (A.im) + (B.re); \
  }

// Y = A/B
#define CXEDIV(Y, A, B)                                 \
  {                                                     \
    ASSERT(B.re != 0. || B.im != 0., "Divide by zero.\n") \
    Y.re = ((A.re) * (B.re) + (A.im) * (B.im)) /        \
           (((B.re) * (B.re)) + ((B.im) * (B.im)));     \
    Y.im = ((A.im) * (B.re) - (A.re) * (B.im)) /        \
           (((B.re) * (B.re)) + ((B.im) * (B.im)));     \
  }

#define CXEADD(Y, A, B) \
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

/****************
**** ASSERT *****
*****************/

/**** ASSERT FUNCTION ****/
extern char str_assert[MAX_CHAR];
/*
void func (MAT*A, MAT*B){
        //check dimension of A and B
        //sprintf(str, " *** [iip_sh_pp] [%s] [A's dim : %d, B's dim : %d]\n",
__func__, A->d, B->d);
        ASSERT_DIM(A, B)
        ASSERT_NULL(mat)
        ASSERT_EQUAL(A,B)
}*/

#define ASSERT(x, msg)                                    \
  {                                                       \
    if ((x) == 0) {                                       \
      printf(" *** [iip_sph_pp] [%s] %s", __func__, msg); \
      getchar();                                          \
      exit(-1);                                           \
    }                                                     \
  }

#define ASSERT_DIM_INVALID()                     \
  {                                              \
    sprintf(str_assert, "Invalid dimension.\n"); \
    ASSERT(0, str_assert);                       \
  }

#define ASSERT_DIM_EQUAL(A, B)                                                 \
  {                                                                            \
    if (A->d0 != B->d0 || A->d1 != B->d1) {                                    \
      sprintf(str_assert, " %d x %d | %d x %d\n", A->d0, A->d1, B->d0, B->d1); \
      ASSERT(0, str_assert)                                                    \
    }                                                                          \
  }

#define ASSERT_MUL(A, B, C)                                                    \
  {                                                                            \
    if (A->d1 != B->d0 || A->d0 != C->d0 || B->d1 != C->d1) {                  \
      sprintf(str_assert, "(%d * %d) X (%d * %d) = (%d * %d)\n", A->d0, A->d1, \
              B->d0, B->d1, C->d0, C->d1);                                     \
      ASSERT(0, str_assert)                                                    \
    }                                                                          \
  }

#define ASSERT_ARG_INVALID()                     \
  {                                              \
    sprintf(str_assert, "Invalid arguments.\n"); \
    ASSERT(0, str_assert);                       \
  }

#define ASSERT_NULL(f)                        \
  {                                           \
    if (f == 0) {                             \
      sprintf(str_assert, "NULL pointer.\n"); \
      ASSERT(0, str_assert)                   \
    }                                         \
  }

#define ASSERT_FILE(f, file_name)                             \
  {                                                           \
    if (f == 0) {                                             \
      sprintf(str_assert, "Invalid file (%s).\n", file_name); \
      ASSERT(0, str_assert)                                   \
    }                                                         \
  }

#define ASSERT_EQUAL(A, B)                                   \
  {                                                          \
    if (A == B) {                                            \
      sprintf(str_assert, "Parameters must be different\n"); \
      ASSERT(0, str_assert)                                  \
    }                                                        \
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

static void* memory_pool[MAX_MEM_PAGE];
static USED_MEM* memory_log[MAX_MEM_PAGE];

static unsigned long long int log_cnt[MAX_MEM_PAGE];
static unsigned int pool_cnt;
signed long long int page_alloc_isable(int page_idx,
                                       unsigned long long int require_size);
void* mpalloc(unsigned long long int size);
void mpfree(void* ptr);
void init(UINT mem_pool_size);
void finit();
#endif
