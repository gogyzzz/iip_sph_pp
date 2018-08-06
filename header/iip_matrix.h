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
#define alloc_mat_load(_x, _3, _2, _1, ...) _1
#define alloc_mat_load_(args_list) alloc_mat_load args_list
#define alloc_mat(...) \
  alloc_mat_load_(     \
      (__VA_ARGS__, alloc_mat_3d, alloc_mat_2d, alloc_mat_1d)(__VA_ARGS__))

MAT* alloc_mat_1d(UINT);
MAT* alloc_mat_2d(UINT, UINT);
MAT* alloc_mat_3d(UINT, UINT, UINT);

#define alloc_cmat_load(_x, _3, _2, _1, ...) _1
#define alloc_cmat_load_(args_list) alloc_cmat_load args_list
#define alloc_cmat(...) \
  alloc_cmat_load_(     \
      (__VA_ARGS__, alloc_cmat_3d, alloc_cmat_2d, alloc_cmat_1d)(__VA_ARGS__))
CMAT* alloc_cmat_1d(UINT);
CMAT* alloc_cmat_2d(UINT, UINT);
CMAT* alloc_cmat_3d(UINT, UINT, UINT);

/**** allocate MAT in memory pool : mpalloc_mat ***/
#define mpalloc_mat_load(_x, _3, _2, _1, ...) _1
#define mpalloc_mat_load_(args_list) mpalloc_mat_load args_list
#define mpalloc_mat(...) \
  mpalloc_mat_load_((__VA_ARGS__, mpalloc_mat_3d, mpalloc_mat_2d, mpalloc_mat_1d)(__VA_ARGS__))
MAT* mpalloc_mat_1d(UINT d0);
MAT* mpalloc_mat_2d(UINT d0, UINT d1);
MAT* mpalloc_mat_3d(UINT d0, UINT d1, UINT d2);

#define mpalloc_cmat_load(_x, _3, _2, _1, ...) _1
#define mpalloc_cmat_load_(args_list) mpalloc_cmat_load args_list
#define mpalloc_cmat(...) \
  mpalloc_cmat_load_(     \
      (__VA_ARGS__, mpalloc_cmat_3d, mpalloc_cmat_2d, mpalloc_cmat_1d)(__VA_ARGS__))
CMAT* mpalloc_cmat_1d(UINT d0);
CMAT* mpalloc_cmat_2d(UINT d0, UINT d1);
CMAT* mpalloc_cmat_3d(UINT d0, UINT d1, UINT d2);

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

/**** mpsubmat overloading ****/

#define mpsubmat_load(_x2, _x3, _x4, _3, _x5, _2, _x6, _1, ...) _1
#define mpsubmat_load_(args_list) mpsubmat_load args_list
#define mpsubmat(...)                                              \
  mpsubmat_load_((__VA_ARGS__, mpsubmat_3d, _, mpsubmat_2d, _, \
                    mpsubmat_1d)(__VA_ARGS__))
MAT* mpsubmat_1d(MAT* src, ITER s0, ITER e0);
MAT* mpsubmat_2d(MAT* src, ITER s0, ITER e0, ITER s1, ITER e1);
MAT* mpsubmat_3d(MAT* src, ITER s0, ITER e0, ITER s1, ITER e1, ITER s2,
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
DTYPE get_by_dim(MAT* mat, DIM* dim);
CTYPE cget_by_dim(CMAT* mat, DIM* dim);

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

void cinv_elements(CMAT* mat);
void cinv_elements_inc(UINT size, CTYPE* X, ITER incx);

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

/**** get Diagonal of Matrix ****/
/* ex
 * mat = | 1 3 |
 *       | 2 4 |
 *
 *       | 5 7 |
 *       | 6 8 |
 * =================
 * dia = | 1 |
 *       | 4 |
 *
 *       | 5 |
 *       | 8 |
 */
void diagonal(MAT* mat, MAT* dia);
void cdiagonal(CMAT* mat, CMAT* dia);

/**** get trace of Matrix ****/
/* ex
 * mat = | 1 3 |
 *       | 2 4 |
 *
 *       | 5 7 |
 *       | 6 8 |
 * =================
 * tr =  | 5  |(= 1 + 4 )
 *
 *       | 13 |(= 5 + 8 )
 */
void trace(MAT* mat, MAT* tr);
void ctrace(CMAT* mat, CMAT* tr);

/**** hermitian ****/
CMAT* create_hermit(CMAT* mat);

void hermit(CMAT* mat);

/**** set Identity Matrix ****/
void ident_mat(MAT* mat);

/**** Inverse Matrix WIP ****/
void invers(MAT* mat);

/****  free memory of MAT ****/
void free_mat(MAT*);
void free_cmat(CMAT*);

/****  print MAT ****/
void print_mat(MAT*);
void print_cmat(CMAT*);

/**** free_mat in memory pool****/
void free_mpalloc_mat(MAT* mat_in_memory_pool);
void free_mpalloc_cmat(CMAT* mat_in_memory_pool);

#if USE_CUDA
__global__ void cu_print_mat(DTYPE*, UINT, UINT, UINT);
__global__ void cu_print_cmat(CTYPE*, UINT, UINT, UINT);
#endif

#endif
