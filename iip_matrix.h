#ifndef IIP_MATRIX_H
#define IIP_MATRIX_H

//#include "mother.h"
#include "iip_type.h"

/********************
 **** iip_matrix ****
 ********************/

void fill(MAT*, DTYPE);
void cfill(CMAT*, DTYPE,DTYPE);

/*** allocMAT ***/
#define alloc_MAT_load(_x,_3,_2,_1,...) _1
#define alloc_MAT_load_(args_list) alloc_MAT_load args_list
#define alloc_MAT(...) alloc_MAT_load_((__VA_ARGS__, alloc_MAT_3d,alloc_MAT_2d,alloc_MAT_1d)(__VA_ARGS__))

MAT* alloc_MAT_1d(UINT);
MAT* alloc_MAT_2d(UINT,UINT);
MAT* alloc_MAT_3d(UINT,UINT,UINT);

#define calloc_MAT_load(_x,_3,_2,_1,...) _1
#define calloc_MAT_load_(args_list) calloc_MAT_load args_list
#define calloc_MAT(...) calloc_MAT_load_((__VA_ARGS__, calloc_MAT_3d,calloc_MAT_2d,calloc_MAT_1d)(__VA_ARGS__))
CMAT* calloc_MAT_1d(UINT);
CMAT* calloc_MAT_2d(UINT,UINT);
CMAT* calloc_MAT_3d(UINT,UINT,UINT);


/**** zeros  ****/
#define zeros_load(_x,_3,_2,_1,...) _1
#define zeros_load_(args_list) zeros_load args_list
#define zeros(...) zeros_load_((__VA_ARGS__, zeros_3d,zeros_2d,zeros_1d)(__VA_ARGS__))

MAT* zeros_1d(UINT);
MAT* zeros_2d(UINT,UINT);
MAT* zeros_3d(UINT,UINT,UINT);

#define czeros_load(_x,_3,_2,_1,...) _1
#define czeros_load_(args_list) czeros_load args_list
#define czeros(...) czeros_load_((__VA_ARGS__, czeros_3d,czeros_2d,czeros_1d)(__VA_ARGS__))
CMAT* czeros_1d(UINT);
CMAT* czeros_2d(UINT,UINT);
CMAT* czeros_3d(UINT,UINT,UINT);

/**** set overloading  ****/
#define set_load(_x,_xx,_xxx,_3,_2,_1,...) _1
#define set_load_(args_list) set_load args_list
#define set(...) set_load_((__VA_ARGS__, set_3d,set_2d,set_1d)(__VA_ARGS__))
void set_1d(MAT*,UINT,DTYPE);
void set_2d(MAT*,UINT,UINT,DTYPE);
void set_3d(MAT*,UINT,UINT,UINT,DTYPE);

#define cset_load(_x,_xx,_xxx,_xxxx,_3,_2,_1,...) _1
#define cset_load_(args_list) cset_load args_list
#define cset(...) cset_load_((__VA_ARGS__, cset_3d,cset_2d,cset_1d)(__VA_ARGS__))
void cset_1d(CMAT*,UINT,DTYPE,DTYPE);
void cset_2d(CMAT*,UINT,UINT,DTYPE,DTYPE);
void cset_3d(CMAT*,UINT,UINT,UINT,DTYPE,DTYPE);

#if USE_CUDA
__global__ void cu_set(DTYPE*,UINT,UINT,UINT,DTYPE);
__global__ void cu_set(CTYPE*,UINT,UINT,UINT,DTYPE,DTYPE);
#endif


/****  get overloadnig ****/
#define get_load(_x,_xx,_3,_2,_1,...) _1
#define get_load_(args_list) get_load args_list
#define get(...) get_load_((__VA_ARGS__, get_3d,get_2d,get_1d)(__VA_ARGS__))
DTYPE get_1d(MAT*,UINT);
DTYPE get_2d(MAT*,UINT,UINT);
DTYPE get_3d(MAT*,UINT,UINT,UINT);

#define cget_load(_x,_xx,_3,_2,_1,...) _1
#define cget_load_(args_list) cget_load args_list
#define cget(...) cget_load_((__VA_ARGS__, cget_3d,cget_2d,cget_1d)(__VA_ARGS__))
CTYPE cget_1d(CMAT*,UINT);
CTYPE cget_2d(CMAT*,UINT,UINT);
CTYPE cget_3d(CMAT*,UINT,UINT,UINT);

/**** submat overloading ****/

//since arg +=2, pend _x for each function
#define submat_load(_x1,_x2,_x3,_x4,_3,_x5,_2,_x6,_1,...) _1
#define submat_load_(args_list) submat_load args_list
#define submat(...) submat_load_((__VA_ARGS__, submat_3d,_,submat_2d,_,submat_1d)(__VA_ARGS__))
void submat_1d(MAT*, MAT*, ITER,ITER);
void submat_2d(MAT*, MAT*, ITER,ITER, ITER,ITER);
void submat_3d(MAT*, MAT*, ITER,ITER, ITER,ITER, ITER,ITER);

#define csubmat_load(_x1,_x2,_x3,_x4,_3,_x5,_2,_x6,_1,...) _1
#define csubmat_load_(args_list) csubmat_load args_list
#define csubmat(...) csubmat_load_((__VA_ARGS__, csubmat_3d,_,csubmat_2d,_,csubmat_1d)(__VA_ARGS__))
void csubmat_1d(CMAT*, CMAT*, ITER,ITER);
void csubmat_2d(CMAT*, CMAT*, ITER,ITER, ITER,ITER);
void csubmat_3d(CMAT*, CMAT*, ITER,ITER, ITER,ITER, ITER,ITER);

/****  miscellaneous ****/ 
void free_MAT(MAT*);
void free_CMAT(CMAT*);
void print_MAT(MAT*);
void print_CMAT(CMAT*);


#endif