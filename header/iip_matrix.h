#ifndef IIP_MATRIX_H
#define IIP_MATRIX_H

//#include "mother.h"
#include "iip_type.h"

/********************
 **** iip_matrix ****
 ********************/

void fill(MAT*, DTYPE);
void cfill(CMAT*, DTYPE, DTYPE);

#if USE_CUDA
__global__ void cu_fill(DTYPE*, UINT, DTYPE, UINT);
__global__ void cu_cfill(CTYPE*, UINT, DTYPE, DTYPE, UINT);
#endif

/*** allocMAT ***/
#define alloc_MAT_load(_x, _3, _2, _1, ...) _1
#define alloc_MAT_load_(args_list) alloc_MAT_load args_list
#define alloc_MAT(...) \
  alloc_MAT_load_(     \
      (__VA_ARGS__, alloc_MAT_3d, alloc_MAT_2d, alloc_MAT_1d)(__VA_ARGS__))

MAT* alloc_MAT_1d(UINT);
MAT* alloc_MAT_2d(UINT, UINT);
MAT* alloc_MAT_3d(UINT, UINT, UINT);

#define alloc_CMAT_load(_x, _3, _2, _1, ...) _1
#define alloc_CMAT_load_(args_list) alloc_CMAT_load args_list
#define alloc_CMAT(...) \
  alloc_CMAT_load_(     \
      (__VA_ARGS__, alloc_CMAT_3d, alloc_CMAT_2d, alloc_CMAT_1d)(__VA_ARGS__))
CMAT* alloc_CMAT_1d(UINT);
CMAT* alloc_CMAT_2d(UINT, UINT);
CMAT* alloc_CMAT_3d(UINT, UINT, UINT);

/**** allocate MAT in memory pool : mem_MAT ***/
#define mem_MAT_load(_x, _3, _2, _1, ...) _1
#define mem_MAT_load_(args_list) mem_MAT_load args_list
#define mem_MAT(...) \
  mem_MAT_load_((__VA_ARGS__, mem_MAT_3d, mem_MAT_2d, mem_MAT_1d)(__VA_ARGS__))
MAT* mem_MAT_1d(UINT d0);
MAT* mem_MAT_2d(UINT d0, UINT d1);
MAT* mem_MAT_3d(UINT d0, UINT d1, UINT d2);

#define mem_CMAT_load(_x, _3, _2, _1, ...) _1
#define mem_CMAT_load_(args_list) mem_CMAT_load args_list
#define mem_CMAT(...) \
  mem_CMAT_load_(     \
      (__VA_ARGS__, mem_CMAT_3d, mem_CMAT_2d, mem_CMAT_1d)(__VA_ARGS__))
CMAT* mem_CMAT_1d(UINT d0);
CMAT* mem_CMAT_2d(UINT d0, UINT d1);
CMAT* mem_CMAT_3d(UINT d0, UINT d1, UINT d2);

/**** zeros  ****/
#define zeros_load(_x, _3, _2, _1, ...) _1
#define zeros_load_(args_list) zeros_load args_list
#define zeros(...) \
  zeros_load_((__VA_ARGS__, zeros_3d, zeros_2d, zeros_1d)(__VA_ARGS__))

MAT* zeros_1d(UINT);
MAT* zeros_2d(UINT, UINT);
MAT* zeros_3d(UINT, UINT, UINT);

#define czeros_load(_x, _3, _2, _1, ...) _1
#define czeros_load_(args_list) czeros_load args_list
#define czeros(...) \
  czeros_load_((__VA_ARGS__, czeros_3d, czeros_2d, czeros_1d)(__VA_ARGS__))
CMAT* czeros_1d(UINT);
CMAT* czeros_2d(UINT, UINT);
CMAT* czeros_3d(UINT, UINT, UINT);

/**** set overloading  ****/
#define set_load(_x, _xx, _xxx, _3, _2, _1, ...) _1
#define set_load_(args_list) set_load args_list
#define set(...) set_load_((__VA_ARGS__, set_3d, set_2d, set_1d)(__VA_ARGS__))
void set_1d(MAT*, UINT, DTYPE);
void set_2d(MAT*, UINT, UINT, DTYPE);
void set_3d(MAT*, UINT, UINT, UINT, DTYPE);

#define cset_load(_x, _xx, _xxx, _xxxx, _3, _2, _1, ...) _1
#define cset_load_(args_list) cset_load args_list
#define cset(...) \
  cset_load_((__VA_ARGS__, cset_3d, cset_2d, cset_1d)(__VA_ARGS__))
void cset_1d(CMAT*, UINT, DTYPE, DTYPE);
void cset_2d(CMAT*, UINT, UINT, DTYPE, DTYPE);
void cset_3d(CMAT*, UINT, UINT, UINT, DTYPE, DTYPE);

#if USE_CUDA
__global__ void cu_set(DTYPE*, UINT, DTYPE);
__global__ void cu_cset(CTYPE*, UINT, DTYPE, DTYPE);
#endif

/****  get overloadnig ****/
#define get_load(_x, _xx, _3, _2, _1, ...) _1
#define get_load_(args_list) get_load args_list
#define get(...) get_load_((__VA_ARGS__, get_3d, get_2d, get_1d)(__VA_ARGS__))
DTYPE get_1d(MAT*, UINT);
DTYPE get_2d(MAT*, UINT, UINT);
DTYPE get_3d(MAT*, UINT, UINT, UINT);

#define cget_load(_x, _xx, _3, _2, _1, ...) _1
#define cget_load_(args_list) cget_load args_list
#define cget(...) \
  cget_load_((__VA_ARGS__, cget_3d, cget_2d, cget_1d)(__VA_ARGS__))
CTYPE cget_1d(CMAT*, UINT);
CTYPE cget_2d(CMAT*, UINT, UINT);
CTYPE cget_3d(CMAT*, UINT, UINT, UINT);

/**** submat overloading ****/

#if USE_CUDA
__global__ void cu_submat(DTYPE*, DTYPE*, ITER, ITER, ITER, ITER, UINT, UINT,
                          UINT, UINT, UINT);
__global__ void cu_csubmat(CTYPE*, CTYPE*, ITER, ITER, ITER, ITER, UINT, UINT,
                           UINT, UINT, UINT);
#endif
// since arg +=2, pend _x for each function
#define submat_load(_x1, _x2, _x3, _x4, _3, _x5, _2, _x6, _1, ...) _1
#define submat_load_(args_list) submat_load args_list
#define submat(...) \
  submat_load_(     \
      (__VA_ARGS__, submat_3d, _, submat_2d, _, submat_1d)(__VA_ARGS__))
void submat_1d(MAT*, MAT*, ITER, ITER);
void submat_2d(MAT*, MAT*, ITER, ITER, ITER, ITER);
void submat_3d(MAT*, MAT*, ITER, ITER, ITER, ITER, ITER, ITER);

#define csubmat_load(_x1, _x2, _x3, _x4, _3, _x5, _2, _x6, _1, ...) _1
#define csubmat_load_(args_list) csubmat_load args_list
#define csubmat(...) \
  csubmat_load_(     \
      (__VA_ARGS__, csubmat_3d, _, csubmat_2d, _, csubmat_1d)(__VA_ARGS__))
void csubmat_1d(CMAT*, CMAT*, ITER, ITER);
void csubmat_2d(CMAT*, CMAT*, ITER, ITER, ITER, ITER);
void csubmat_3d(CMAT*, CMAT*, ITER, ITER, ITER, ITER, ITER, ITER);

/**** mem_submat overloading ****/

#define mem_submat_load(_x2, _x3, _x4, _3, _x5, _2, _x6, _1, ...) _1
#define mem_submat_load_(args_list) mem_submat_load args_list
#define mem_submat(...)                                              \
  mem_submat_load_((__VA_ARGS__, mem_submat_3d, _, mem_submat_2d, _, \
                    mem_submat_1d)(__VA_ARGS__))
MAT* mem_submat_1d(MAT* src, ITER s0, ITER e0);
MAT* mem_submat_2d(MAT* src, ITER s0, ITER e0, ITER s1, ITER e1);
MAT* mem_submat_3d(MAT* src, ITER s0, ITER e0, ITER s1, ITER e1, ITER s2,
                   ITER e2);

#define mem_csubmat_load(_x2, _x3, _x4, _3, _x5, _2, _x6, _1, ...) _1
#define mem_csubmat_load_(args_list) mem_csubmat_load args_list
#define mem_csubmat(...)                                                \
  mem_csubmat_load_((__VA_ARGS__, mem_csubmat_3d, _, mem_csubmat_2d, _, \
                     mem_csubmat_1d)(__VA_ARGS__))
CMAT* mem_csubmat_1d(CMAT* src, ITER s0, ITER e0);
CMAT* mem_csubmat_2d(CMAT* src, ITER s0, ITER e0, ITER s1, ITER e1);
CMAT* mem_csubmat_3d(CMAT* src, ITER s0, ITER e0, ITER s1, ITER e1, ITER s2,
                     ITER e2);

/**** DIM allocator ***/
DIM* new_dim();

/**** element operation by DIM ***/
DTYPE getbydim(MAT* mat, DIM* dim);
CTYPE cgetbydim(CMAT* mat, DIM* dim);

void setbydim(MAT* mat, DIM* dim, DTYPE val);
void csetbydim(CMAT* mat, DIM* dim, CTYPE val);

/**** add elements - broadcasting operation ****/
void add_elements(MAT* A, MAT* B, MAT* C);
void cadd_elements(CMAT* A, CMAT* B, CMAT* C);

/**** multiply elements - broadcasting operation ****/
void mul_elements(MAT* A, MAT* B, MAT* C);
void cmul_elements(CMAT* A, CMAT* B, CMAT* C);

/**** divide elements - broadcasting operation ****/
void div_elements(MAT* A, MAT* B, MAT* C);
void cdiv_elements(CMAT* A, CMAT* B, CMAT* C);

/**** inverse elements ***/
void inv_elements(MAT* mat);
void inv_elements_inc(UINT size, DTYPE* X, ITER incx);

void cinv_elelments(CMAT* mat);
void cinv_elelments_inc(UINT size, CTYPE* X, ITER incx);

/**** repeat matrix****/
void repmat(MAT* mat, DIM* dim);
void crepmat(CMAT* mat, DIM* dim);

/**** reshape matrix ****/
void reshape(MAT* mat, DIM* dim);
void creshape(CMAT* mat, DIM* dim);

/**** shift dimension ****/
void shiftdim(MAT* mat, SINT n);
void cshiftdim(CMAT* mat, SINT n);

/**** permute matrix ****/
void permute(MAT* mat, UINT seq);
void cpermute(CMAT* mat, UINT seq);

/**** transpose  ****/
MAT* create_trans(MAT* mat);
CMAT* create_ctrans(CMAT* mat);

void trans(MAT* mat);  // with memory pool
void ctrans(CMAT* mat);
/**** hermitian ****/
CMAT* create_hermit(CMAT* mat);

void hermit(CMAT* mat);
/**** Identity Matrix ****/
void id_MAT(MAT* mat);

/**** Inverse Matrix WIP ****/
void invers(MAT* mat);

/****  free memory of MAT ****/
void free_MAT(MAT*);
void free_CMAT(CMAT*);

/****  print MAT ****/
void print_MAT(MAT*);
void print_CMAT(CMAT*);

/**** free_MAT in memory pool****/
void free_mem_MAT(MAT* mat_in_memory_pool);
void free_mem_CMAT(CMAT* mat_in_memory_pool);

#if USE_CUDA
__global__ void cu_print_MAT(DTYPE*, UINT, UINT, UINT);
__global__ void cu_print_CMAT(CTYPE*, UINT, UINT, UINT);
#endif

#endif
