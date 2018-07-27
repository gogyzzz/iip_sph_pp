//#include "mother.h"
#include "iip_matrix.h"

/**** alloc_MAT ****/

MAT *alloc_MAT_1d(UINT d0) {
  MAT *mat;
#if DEBUG
  printf("%s\n", __func__);
#endif

  mat = (MAT *)malloc(sizeof(MAT));
  mat->ndim = 0;
  mat->d0 = d0;
  mat->d1 = 1;
  mat->d2 = 1;

  mat->data = (DTYPE *)malloc(sizeof(DTYPE) * d0);

  return mat;
}
MAT *alloc_MAT_2d(UINT d0, UINT d1) {
  MAT *mat;
#if DEBUG
  printf("%s\n", __func__);
#endif

  mat = (MAT *)malloc(sizeof(MAT));
  mat->ndim = 1;
  mat->d0 = d0;
  mat->d1 = d1;
  mat->d2 = 1;

  mat->data = (DTYPE *)malloc(sizeof(DTYPE) * d0 * d1);

  return mat;
}
MAT *alloc_MAT_3d(UINT d0, UINT d1, UINT d2) {
  MAT *mat;
#if DEBUG
  printf("%s\n", __func__);
#endif

  mat = (MAT *)malloc(sizeof(MAT));
  mat->ndim = 2;
  mat->d0 = d0;
  mat->d1 = d1;
  mat->d2 = d2;

  mat->data = (DTYPE *)malloc(sizeof(DTYPE) * d0 * d1 * d2);

  return mat;
}

CMAT *alloc_CMAT_1d(UINT d0) {
  CMAT *mat;
#if DEBUG
  printf("%s\n", __func__);
#endif

  mat = (CMAT *)malloc(sizeof(CMAT));
  mat->ndim = 0;
  mat->d0 = d0;
  mat->d1 = 1;
  mat->d2 = 1;

  mat->data = (CTYPE *)malloc(sizeof(CTYPE) * d0);

  return mat;
}
CMAT *alloc_CMAT_2d(UINT d0, UINT d1) {
  CMAT *mat;
#if DEBUG
  printf("%s\n", __func__);
#endif

  mat = (CMAT *)malloc(sizeof(CMAT));
  mat->ndim = 1;
  mat->d0 = d0;
  mat->d1 = d1;
  mat->d2 = 1;

  mat->data = (CTYPE *)malloc(sizeof(CTYPE) * d0 * d1);

  return mat;
}
CMAT *alloc_CMAT_3d(UINT d0, UINT d1, UINT d2) {
  CMAT *mat;
#if DEBUG
  printf("%s\n", __func__);
#endif

  mat = (CMAT *)malloc(sizeof(CMAT));
  mat->ndim = 2;
  mat->d0 = d0;
  mat->d1 = d1;
  mat->d2 = d2;

  mat->data = (CTYPE *)malloc(sizeof(CTYPE) * d0 * d1 * d2);

  return mat;
}

/**** allocate MAT in memory pool : mem_MAT ***/
MAT *mem_MAT_1d(UINT d0) {
  MAT *mat;
#if DEBUG
  printf("%s\n", __func__);
#endif
  mat = iip_malloc(sizeof(MAT));
  mat->ndim = 0;
  mat->d0 = d0;
  mat->d1 = 1;
  mat->d2 = 1;
  mat->data = iip_malloc(sizeof(DTYPE) * d0);
  return mat;
}
MAT *mem_MAT_2d(UINT d0, UINT d1) {
  MAT *mat;
#if DEBUG
  printf("%s\n", __func__);
#endif
  mat = iip_malloc(sizeof(MAT));
  mat->ndim = 1;
  mat->d0 = d0;
  mat->d1 = d1;
  mat->d2 = 1;
  mat->data = iip_malloc(sizeof(DTYPE) * d0 * d1);
  return mat;
}
MAT *mem_MAT_3d(UINT d0, UINT d1, UINT d2) {
  MAT *mat;
#if DEBUG
  printf("%s\n", __func__);
#endif
  mat = iip_malloc(sizeof(MAT));
  mat->ndim = 2;
  mat->d0 = d0;
  mat->d1 = d1;
  mat->d2 = d2;
  mat->data = iip_malloc(sizeof(DTYPE) * d0 * d1 * d2);
  return mat;
}

CMAT *mem_CMAT_1d(UINT d0) {
  CMAT *mat;
#if DEBUG
  printf("%s\n", __func__);
#endif
  mat = iip_malloc(sizeof(CMAT));
  mat->ndim = 0;
  mat->d0 = d0;
  mat->d1 = 1;
  mat->d2 = 1;
  mat->data = iip_malloc(sizeof(CTYPE) * d0);
  return mat;
}
CMAT *mem_CMAT_2d(UINT d0, UINT d1) {
  CMAT *mat;
#if DEBUG
  printf("%s\n", __func__);
#endif
  mat = iip_malloc(sizeof(CMAT));
  mat->ndim = 1;
  mat->d0 = d0;
  mat->d1 = d1;
  mat->d2 = 1;
  mat->data = iip_malloc(sizeof(CTYPE) * d0 * d1);
  return mat;
}
CMAT *mem_CMAT_3d(UINT d0, UINT d1, UINT d2) {
  CMAT *mat;
#if DEBUG
  printf("%s\n", __func__);
#endif
  mat = iip_malloc(sizeof(CMAT));
  mat->ndim = 2;
  mat->d0 = d0;
  mat->d1 = d1;
  mat->d2 = d2;
  mat->data = iip_malloc(sizeof(CTYPE) * d0 * d1 * d2);
  return mat;
}
/**** zeros  ****/

MAT *zeros_1d(UINT d0) {
  MAT *mat;
#if DEBUG
  printf("%s\n", __func__);
#endif
  mat = (MAT *)malloc(sizeof(MAT));
  mat->ndim = 0;
  mat->d0 = d0;
  mat->d1 = 1;
  mat->d2 = 1;

  mat->data = (DTYPE *)malloc(sizeof(DTYPE) * d0);
  memset(mat->data, 0, sizeof(DTYPE) * d0);

  return mat;
}

MAT *zeros_2d(UINT d0, UINT d1) {
#if DEBUG
  printf("%s\n", __func__);
#endif
  MAT *mat = (MAT *)malloc(sizeof(MAT));
  mat->ndim = 1;
  mat->d0 = d0;
  mat->d1 = d1;
  mat->d2 = 1;

  mat->data = (DTYPE *)malloc(sizeof(DTYPE) * d0 * d1);
  memset(mat->data, 0, sizeof(DTYPE) * d0 * d1);

  return mat;
}
MAT *zeros_3d(UINT d0, UINT d1, UINT d2) {
#if DEBUG
  printf("%s\n", __func__);
#endif
  MAT *mat = (MAT *)malloc(sizeof(MAT));
  mat->ndim = 2;
  mat->d0 = d0;
  mat->d1 = d1;
  mat->d2 = d2;

  mat->data = (DTYPE *)malloc(sizeof(DTYPE) * d0 * d1 * d2);
  memset(mat->data, 0, sizeof(DTYPE) * d0 * d1 * d2);

  return mat;
}

CMAT *czeros_1d(UINT d0) {
#if DEBUG
  printf("%s\n", __func__);
#endif
  CMAT *mat = (CMAT *)malloc(sizeof(CMAT));
  mat->ndim = 0;
  mat->d0 = d0;
  mat->d1 = 1;
  mat->d2 = 1;

  mat->data = (CTYPE *)malloc(sizeof(CTYPE) * d0);
  memset(mat->data, 0, sizeof(CTYPE) * d0);

  return mat;
}

CMAT *czeros_2d(UINT d0, UINT d1) {
#if DEBUG
  printf("%s\n", __func__);
#endif
  CMAT *mat = (CMAT *)malloc(sizeof(CMAT));
  mat->ndim = 1;
  mat->d0 = d0;
  mat->d1 = d1;
  mat->d2 = 1;

  mat->data = (CTYPE *)malloc(sizeof(CTYPE) * d0 * d1);
  memset(mat->data, 0, sizeof(CTYPE) * d0 * d1);

  return mat;
}
CMAT *czeros_3d(UINT d0, UINT d1, UINT d2) {
#if DEBUG
  printf("%s\n", __func__);
#endif
  CMAT *mat = (CMAT *)malloc(sizeof(CMAT));
  mat->ndim = 2;
  mat->d0 = d0;
  mat->d1 = d1;
  mat->d2 = d2;

  mat->data = (CTYPE *)malloc(sizeof(CTYPE) * d0 * d1 * d2);
  memset(mat->data, 0, sizeof(CTYPE) * d0 * d1 * d2);

  return mat;
}

/**** set ****/

void set_1d(MAT *mat, UINT idx0, DTYPE val) {
#if DEBUG
  printf("%s\n", __func__);
#endif
  mat->data[idx0] = val;
}
void set_2d(MAT *mat, UINT idx0, UINT idx1, DTYPE val) {
#if DEBUG
  printf("%s\n", __func__);
#endif
  mat->data[idx0 + (mat->d0) * idx1] = val;
  printf("i0: %d, i1: %d, val: %lf\n", idx0, idx1, val);
}
void set_3d(MAT *mat, UINT idx0, UINT idx1, UINT idx2, DTYPE val) {
#if DEBUG
  printf("%s\n", __func__);
#endif
  mat->data[idx0 + (mat->d0) * idx1 + (mat->d0) * (mat->d1) * idx2] = val;
}

void cset_1d(CMAT *mat, UINT idx0, DTYPE re, DTYPE im) {
#if DEBUG
  printf("%s\n", __func__);
#endif
  mat->data[idx0].re = re;
  mat->data[idx0].im = im;
}
void cset_2d(CMAT *mat, UINT idx0, UINT idx1, DTYPE re, DTYPE im) {
#if DEBUG
  printf("%s\n", __func__);
#endif
  mat->data[idx0 + (mat->d0) * idx1].re = re;
  mat->data[idx0 + (mat->d0) * idx1].im = im;
}
void cset_3d(CMAT *mat, UINT idx0, UINT idx1, UINT idx2, DTYPE re, DTYPE im) {
#if DEBUG
  printf("%s\n", __func__);
#endif
  mat->data[idx0 + (mat->d0) * idx1 + (mat->d0) * (mat->d1) * idx2].re = re;
  mat->data[idx0 + (mat->d0) * idx1 + (mat->d0) * (mat->d1) * idx2].im = im;
}

/**** fill ****/

void fill(MAT *mat, DTYPE val)  // for real mat
{
  ITER i = 0;
  UINT len = mat->d0 * mat->d1 * mat->d2;
#if DEBUG
  printf("%s\n", __func__);
#endif
  for (i = 0; i < len; i++) mat->data[i] = val;
}

void cfill(CMAT *cmat, DTYPE re, DTYPE im)  // for complex mat
{
  ITER i = 0;
  UINT len = cmat->d0 * cmat->d1 * cmat->d2;
#if DEBUG
  printf("%s\n", __func__);
#endif
  for (i = 0; i < len; i++) {
    cmat->data[i].re = re;
    cmat->data[i].im = im;
  }
}

/**** get ****/

DTYPE get_1d(MAT *mat, UINT idx0) {
#if DEBUG
  printf("%s\n", __func__);
#endif
  return mat->data[idx0];
}

DTYPE get_2d(MAT *mat, UINT idx0, UINT idx1) {
#if DEBUG
  printf("%s\n", __func__);
#endif
  return mat->data[idx0 + (mat->d0) * idx1];
}
DTYPE get_3d(MAT *mat, UINT idx0, UINT idx1, UINT idx2) {
#if DEBUG
  printf("%s\n", __func__);
#endif
  return mat->data[idx0 + (mat->d0) * idx1 + (mat->d0) * (mat->d1) * idx0];
}

CTYPE cget_1d(CMAT *mat, UINT idx0) {
#if DEBUG
  printf("%s\n", __func__);
#endif
  return mat->data[idx0];
}

CTYPE cget_2d(CMAT *mat, UINT idx0, UINT idx1) {
#if DEBUG
  printf("%s\n", __func__);
#endif
  return mat->data[idx0 + (mat->d0) * idx1];
}
CTYPE cget_3d(CMAT *mat, UINT idx0, UINT idx1, UINT idx2) {
#if DEBUG
  printf("%s\n", __func__);
#endif
  return mat->data[idx0 + (mat->d0) * idx1 + (mat->d0) * (mat->d1) * idx0];
}

/**** submat ****/

void submat_1d(MAT *mat, MAT *submat, ITER d0_st, ITER d0_ed) {
#if DEBUG
  printf("%s\n", __func__);
#endif
  submat_3d(mat, submat, d0_st, d0_ed, -1, -1, -1, -1);
}

void submat_2d(MAT *mat, MAT *submat, ITER d0_st, ITER d0_ed, ITER d1_st,
               ITER d1_ed) {
#if DEBUG
  printf("%s\n", __func__);
#endif
  submat_3d(mat, submat, d0_st, d0_ed, d1_st, d1_ed, -1, -1);
}

void submat_3d(MAT *mat, MAT *submat, ITER d0_st, ITER d0_ed, ITER d1_st,
               ITER d1_ed, ITER d2_st, ITER d2_ed) {
  ITER i, j, k;

#if DEBUG
  printf("%s\n", __func__);
#endif
  if (mat->d0 < submat->d0 || mat->d1 < submat->d1 || mat->d2 < submat->d2) {
    printf(
        "error in \nmat->d0 < submat->d0 || mat->d1 < submat->d1 || mat->d2 < "
        "submat->d2\n");
    return;
  }

  if (d0_st == -1) d0_st = 0;
  if (d0_ed == -1) d0_ed = mat->d0;
  if (d1_st == -1) d1_st = 0;
  if (d1_ed == -1) d1_ed = mat->d1;
  if (d2_st == -1) d2_st = 0;
  if (d2_ed == -1) d2_ed = mat->d2;
#pragma omp parallel for shared(submat, mat) private(i, j, k)
  for (i = 0; i < d0_ed - d0_st; i++)
    for (j = 0; d1_st + j < d1_ed; j++)
      for (k = 0; d2_st + k < d2_ed; k++) {
        submat->data[i + j * (submat->d0) + k * (submat->d0 * submat->d1)] =
            mat->data[(d0_st + i) + (d1_st + j) * (mat->d0) +
                      (d2_st + k) * (mat->d0 * mat->d1)];
      }
}

void csubmat_1d(CMAT *mat, CMAT *submat, ITER d0_st, ITER d0_ed) {
#if DEBUG
  printf("%s\n", __func__);
#endif
  csubmat_3d(mat, submat, d0_st, d0_ed, -1, -1, -1, -1);
}

void csubmat_2d(CMAT *mat, CMAT *submat, ITER d0_st, ITER d0_ed, ITER d1_st,
                ITER d1_ed) {
#if DEBUG
  printf("%s\n", __func__);
#endif
  csubmat_3d(mat, submat, d0_st, d0_ed, d1_st, d1_ed, -1, -1);
}

void csubmat_3d(CMAT *mat, CMAT *submat, ITER d0_st, ITER d0_ed, ITER d1_st,
                ITER d1_ed, ITER d2_st, ITER d2_ed) {
  ITER i, j, k;

#if DEBUG
  printf("%s\n", __func__);
#endif
  if (mat->d0 < submat->d0 || mat->d1 < submat->d1 || mat->d2 < submat->d2) {
    printf(
        "error in \nmat->d0 < submat->d0 || mat->d1 < submat->d1 || mat->d2 < "
        "submat->d2\n");
    return;
  }
  if (d0_st == -1) d0_st = 0;
  if (d0_ed == -1) d0_ed = mat->d0;
  if (d1_st == -1) d1_st = 0;
  if (d1_ed == -1) d1_ed = mat->d1;
  if (d2_st == -1) d2_st = 0;
  if (d2_ed == -1) d2_ed = mat->d2;

#pragma omp parallel for shared(submat, mat) private(i, j, k)
  for (i = 0; i < d0_ed - d0_st; i++)
    for (j = 0; d1_st + j < d1_ed; j++)
      for (k = 0; d2_st + k < d2_ed; k++) {
        submat->data[i + j * (submat->d0) + k * (submat->d0 * submat->d1)].re =
            mat->data[(d0_st + i) + (d1_st + j) * (mat->d0) +
                      (d2_st + k) * (mat->d0 * mat->d1)]
                .re;
        submat->data[i + j * (submat->d0) + k * (submat->d0 * submat->d1)].im =
            mat->data[(d0_st + i) + (d1_st + j) * (mat->d0) +
                      (d2_st + k) * (mat->d0 * mat->d1)]
                .im;
      }
}

/**** mem_submat ****/

MAT *mem_submat_1d(MAT *mat, ITER d0_st, ITER d0_ed) {
  return mem_submat_3d(mat, d0_st, d0_ed, -1, -1, -1, -1);
}
MAT *mem_submat_2d(MAT *mat, ITER d0_st, ITER d0_ed, ITER d1_st, ITER d1_ed) {
  return mem_submat_3d(mat, d0_st, d0_ed, d1_st, d1_ed, -1, -1);
}
MAT *mem_submat_3d(MAT *mat, ITER d0_st, ITER d0_ed, ITER d1_st, ITER d1_ed,
                   ITER d2_st, ITER d2_ed) {
  ITER i, j, k;
  MAT *submat;

#if DEBUG
  printf("%s\n", __func__);
#endif

  if (d0_st == -1) d0_st = 0;
  if (d0_ed == -1) d0_ed = mat->d0;
  if (d1_st == -1) d1_st = 0;
  if (d1_ed == -1) d1_ed = mat->d1;
  if (d2_st == -1) d2_st = 0;
  if (d2_ed == -1) d2_ed = mat->d2;

  submat = mem_MAT((d0_ed - d0_st), (d1_ed - d1_st), (d2_ed - d2_st));
  if (mat->d0 < submat->d0 || mat->d1 < submat->d1 || mat->d2 < submat->d2) {
    printf(
        "ERROR-FALSE : mat->d0 < submat->d0 || mat->d1 < submat->d1 || mat->d2 "
        "< submat->d2\n");
    printf("%u %u %u %u %u %u\n", mat->d0, submat->d0, mat->d1, submat->d1,
           mat->d2, submat->d2);
    iip_free(submat);
    return;
  }

#pragma omp parallel for shared(submat, mat) private(i, j, k)
  for (i = 0; i < d0_ed - d0_st; i++)
    for (j = 0; d1_st + j < d1_ed; j++)
      for (k = 0; d2_st + k < d2_ed; k++) {
        submat->data[i + j * (submat->d0) + k * (submat->d0 * submat->d1)] =
            mat->data[(d0_st + i) + (d1_st + j) * (mat->d0) +
                      (d2_st + k) * (mat->d0 * mat->d1)];
      }
  return submat;
}

CMAT *mem_csubmat_1d(CMAT *mat, ITER d0_st, ITER d0_ed) {
  return mem_csubmat_3d(mat, d0_st, d0_ed, -1, -1, -1, -1);
}
CMAT *mem_csubmat_2d(CMAT *Cmat, ITER d0_st, ITER d0_ed, ITER d1_st,
                     ITER d1_ed) {
  return mem_csubmat_3d(Cmat, d0_st, d0_ed, d1_st, d1_ed, -1, -1);
}
CMAT *mem_csubmat_3d(CMAT *mat, ITER d0_st, ITER d0_ed, ITER d1_st, ITER d1_ed,
                     ITER d2_st, ITER d2_ed) {
  ITER i, j, k;
  CMAT *submat;

#if DEBUG
  printf("%s\n", __func__);
#endif

  if (d0_st == -1) d0_st = 0;
  if (d0_ed == -1) d0_ed = mat->d0;
  if (d1_st == -1) d1_st = 0;
  if (d1_ed == -1) d1_ed = mat->d1;
  if (d2_st == -1) d2_st = 0;
  if (d2_ed == -1) d2_ed = mat->d2;

  submat = mem_CMAT((d0_ed - d0_st), (d1_ed - d1_st), (d2_ed - d2_st));
  if (mat->d0 < submat->d0 || mat->d1 < submat->d1 || mat->d2 < submat->d2) {
    printf(
        "ERROR-FALSE : mat->d0 < submat->d0 || mat->d1 < submat->d1 || mat->d2 "
        "< submat->d2\n");
    printf("%u %u %u %u %u %u\n", mat->d0, submat->d0, mat->d1, submat->d1,
           mat->d2, submat->d2);
    iip_free(submat);
    return;
  }

#pragma omp parallel for shared(submat, mat) private(i, j, k)
  for (i = 0; i < d0_ed - d0_st; i++)
    for (j = 0; d1_st + j < d1_ed; j++)
      for (k = 0; d2_st + k < d2_ed; k++) {
        submat->data[i + j * (submat->d0) + k * (submat->d0 * submat->d1)].re =
            mat->data[(d0_st + i) + (d1_st + j) * (mat->d0) +
                      (d2_st + k) * (mat->d0 * mat->d1)]
                .re;
        submat->data[i + j * (submat->d0) + k * (submat->d0 * submat->d1)].im =
            mat->data[(d0_st + i) + (d1_st + j) * (mat->d0) +
                      (d2_st + k) * (mat->d0 * mat->d1)]
                .im;
      }
  return submat;
}

/**** DIM allocator ****/
DIM *new_dim() {
  DIM *dim = (DIM *)malloc(sizeof(DIM));
#if DEBUG
  printf("%s\n", __func__);
#endif
  return dim;
}

/**** element operation by DIM ****/
DTYPE getbydim(MAT *mat, DIM *dim) {
  return mat->data[dim->d2 * mat->d0 * mat->d1 + dim->d1 * mat->d1 + dim->d0];
}
CTYPE cgetbydim(CMAT *mat, DIM *dim) {
  return mat->data[dim->d2 * mat->d0 * mat->d1 + dim->d1 * mat->d1 + dim->d0];
}

void setbydim(MAT *mat, DIM *dim, DTYPE val) {
  mat->data[dim->d2 * mat->d0 * mat->d1 + dim->d1 * mat->d1 + dim->d0] = val;
}
void csetbydim(CMAT *mat, DIM *dim, CTYPE val) {
  mat->data[dim->d2 * mat->d0 * mat->d1 + dim->d1 * mat->d1 + dim->d0] = val;
}

/**** add elements - broadcasting ****/
void add_elements(MAT *A, MAT *B, MAT *C) {
  ITER i, j;
  UINT a0, a1, b0, b1;
  UINT a2, b2, c2;
  UINT ia, ib, ic;
#if DEBUG
  printf("%s\n", __func__);
#endif

  a2 = A->d2;
  b2 = B->d2;
  c2 = C->d2;
  a0 = A->d0;
  a1 = A->d1;
  b0 = B->d0;
  b1 = B->d1;

  ic = C->d0 * C->d1;
  /*1 to 1*/
  if (a2 == b2) {
    if (c2 != a2) ASSERT(DIM_INVAL)
    ia = a0 * a1;
    ib = b0 * b1;
  }
  /*1 to b2*/
  else if (a2 == 1 && b2 != 1) {
    if (b2 != c2) ASSERT(DIM_INVAL)
    ia = 0;
    ib = b0 * b1;
  }
  /* 1 to a2 */
  else if (a2 != 1 && b2 == 1) {
    if (a2 != c2) ASSERT(DIM_INVAL)
    ib = 0;
    ia = a0 * a1;
  }
  /* unmatching*/
  else
    ASSERT(DIM_INVAL)

  /* 12 broadcasting cases..
   * a0 a1 b0 b1 | c0 c1
   * 1  1  1  b1 | 1  b1
   * 1  1  b0 1  | b0 1
   * 1  a1 1  1  | 1  a1
   * a0 1  1  1  | a0 1
   *
   * 1  1  b0 b1 | b0 b1
   * a0 a1 1  1  | a0 a1
   * 1  a1 b0 1  | b0 a1
   * a0 1  1  b1 | a0 b1
   *
   * 1  a1 b0 b1 | b0 a1 (a1==b1)
   * a0 1  b0 b1 | b0 b1 (a0==b0)
   * a0 a1 b0  1 | b0 a1 (a0==b0)
   * a0 a1 1  b1 | a0 b1 (a1==b1)
   */
  // 1 1 1 b1
  if (((a0 == 1) && (a1 == 1)) && ((b0 == 1) && (b1 != 1))) {
    if ((C->d0 != 1) || (C->d1 != b1)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < b1; i++) {
        C->data[ic * j + i] = A->data[ia * j + 0] + B->data[ib * j + i];
      }
  }
  // 1 1 b0 1
  else if (((a0 == 1) && (a1 == 1)) && ((b0 != 1) && (b1 == 1))) {
    if ((C->d0 != b0) || (C->d1 != 1)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < b1; i++) {
        C->data[ic * j + i] = A->data[ia * j + 0] + B->data[ib * j + i];
      }
  }
  // 1 a1 1 1
  else if (((a0 == 1) && (a1 != 1)) && ((b0 == 1) && (b1 == 1))) {
    if ((C->d0 != 1) || (C->d1 != a1)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < b1; i++) {
        C->data[ic * j + i] = A->data[ia * j + i] + B->data[ib * j + 0];
      }
  }
  // a0 1 1 1
  else if (((a0 != 1) && (a1 == 1)) && ((b0 == 1) && (b1 == 1))) {
    if ((C->d0 != a0) || (C->d1 != 1)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < b1; i++) {
        C->data[ic * j + i] = A->data[ia * j + i] + B->data[ib * j + 0];
      }
  }
  // 1 1 b0 b1
  else if (((a0 == 1) && (a1 == 1)) && ((b0 != 1) && (b1 != 1))) {
    if ((C->d0 != b0) || (C->d1 != b1)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < b0 * b1; i++) {
        C->data[ic * j + i] = A->data[ia * j + 0] + B->data[ib * j + i];
      }
  }
  // a0 a1 1 1
  else if (((a0 != 1) && (a1 != 1)) && ((b0 == 1) && (b1 == 1))) {
    if ((C->d0 != a0) || (C->d1 != a1)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < a0 * a1; i++) {
        C->data[ic * j + i] = A->data[ia * j + i] + B->data[ib * j + 0];
      }
  }
  // 1 a1 b0 1
  else if (((a0 == 1) && (a1 != 1)) && ((b0 != 1) && (b1 == 1))) {
    if ((C->d0 != b0) || (C->d1 != a1)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < b0 * a1; i++) {
        C->data[ic * j + i] =
            A->data[ia * j + i / a1] + B->data[ib * j + i % b0];
      }
  }
  // a0 1 1 b1
  else if (((a0 != 1) && (a1 == 1)) && ((b0 == 1) && (b1 != 1))) {
    if ((C->d0 != a0) || (C->d1 != b1)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < a0 * b1; i++) {
        C->data[ic * j + i] =
            A->data[ia * j + i % a0] + B->data[ib * j + i / b1];
      }
  }
  // 1  a1 b0 b1 | b0 a1 (a1==b1)
  else if (((a0 == 1) && (a1 != 1)) && ((b0 != 1) && (b1 != 1))) {
    if (((C->d0 != b0) || (C->d1 != a1)) || (a1 != b1)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < b0 * a1; i++) {
        C->data[ic * j + i] = A->data[ia * j + i / b1] + B->data[ib * j + i];
      }
  }
  // a0 1  b0 b1 | b0 b1 (a0==b0)
  else if (((a0 != 1) && (a1 == 1)) && ((b0 != 1) && (b1 != 1))) {
    if (((C->d0 != b0) || (C->d1 != b1)) || (a0 != b0)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < b0 * b1; i++) {
        C->data[ic * j + i] = A->data[ia * j + i % b0] + B->data[ib * j + i];
      }
  }
  // a0 a1 b0  1 | b0 a1 (a0==b0)
  else if (((a0 != 1) && (a1 != 1)) && ((b0 != 1) && (b1 == 1))) {
    if (((C->d0 != b0) || (C->d1 != a1)) || (a0 != b0)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < b0 * a1; i++) {
        C->data[ic * j + i] = A->data[ia * j + i] + B->data[ib * j + i % a1];
      }
  }
  // a0 a1 1  b1 | a0 b1 (a1==b1)
  else if (((a0 != 1) && (a1 != 1)) && ((b0 == 1) && (b1 != 1))) {
    if (((C->d0 != a0) || (C->d1 != b1)) || (a1 != b1)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < a0 * b1; i++) {
        C->data[ic * j + i] = A->data[ia * j + i] + B->data[ib * j + i / a0];
      }
  } else
    ASSERT(DIM_INVAL)
}

void cadd_elements(CMAT *A, CMAT *B, CMAT *C) {
  ITER i, j;
  UINT a0, a1, b0, b1;
  UINT a2, b2, c2;
  UINT ia, ib, ic;
#if DEBUG
  printf("%s\n", __func__);
#endif

  a2 = A->d2;
  b2 = B->d2;
  c2 = C->d2;
  a0 = A->d0;
  a1 = A->d1;
  b0 = B->d0;
  b1 = B->d1;

  ic = C->d0 * C->d1;
  /*1 to 1*/
  if (a2 == b2) {
    if (c2 != a2) ASSERT(DIM_INVAL)
    ia = a0 * a1;
    ib = b0 * b1;
  }
  /*1 to b2*/
  else if (a2 == 1 && b2 != 1) {
    if (b2 != c2) ASSERT(DIM_INVAL)
    ia = 0;
    ib = b0 * b1;
  }
  /* 1 to a2 */
  else if (a2 != 1 && b2 == 1) {
    if (a2 != c2) ASSERT(DIM_INVAL)
    ib = 0;
    ia = a0 * a1;
  }
  /* unmatching*/
  else
    ASSERT(DIM_INVAL)

  /* 12 broadcasting cases..
   * a0 a1 b0 b1 | c0 c1
   * 1  1  1  b1 | 1  b1
   * 1  1  b0 1  | b0 1
   * 1  a1 1  1  | 1  a1
   * a0 1  1  1  | a0 1
   *
   * 1  1  b0 b1 | b0 b1
   * a0 a1 1  1  | a0 a1
   * 1  a1 b0 1  | b0 a1
   * a0 1  1  b1 | a0 b1
   *
   * 1  a1 b0 b1 | b0 a1 (a1==b1)
   * a0 1  b0 b1 | b0 b1 (a0==b0)
   * a0 a1 b0  1 | b0 a1 (a0==b0)
   * a0 a1 1  b1 | a0 b1 (a1==b1)
   */
  // 1 1 1 b1
  if (((a0 == 1) && (a1 == 1)) && ((b0 == 1) && (b1 != 1))) {
    if ((C->d0 != 1) || (C->d1 != b1)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < b1; i++) {
        cxeadd(C->data[ic * j + i], A->data[ia * j + 0], B->data[ib * j + i])
      }
  }
  // 1 1 b0 1
  else if (((a0 == 1) && (a1 == 1)) && ((b0 != 1) && (b1 == 1))) {
    if ((C->d0 != b0) || (C->d1 != 1)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < b1; i++) {
        cxeadd(C->data[ic * j + i], A->data[ia * j + 0], B->data[ib * j + i])
      }
  }
  // 1 a1 1 1
  else if (((a0 == 1) && (a1 != 1)) && ((b0 == 1) && (b1 == 1))) {
    if ((C->d0 != 1) || (C->d1 != a1)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < b1; i++) {
        cxeadd(C->data[ic * j + i], A->data[ia * j + i], B->data[ib * j + 0])
      }
  }
  // a0 1 1 1
  else if (((a0 != 1) && (a1 == 1)) && ((b0 == 1) && (b1 == 1))) {
    if ((C->d0 != a0) || (C->d1 != 1)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < b1; i++) {
        cxeadd(C->data[ic * j + i], A->data[ia * j + i], B->data[ib * j + 0])
      }
  }
  // 1 1 b0 b1
  else if (((a0 == 1) && (a1 == 1)) && ((b0 != 1) && (b1 != 1))) {
    if ((C->d0 != b0) || (C->d1 != b1)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < b0 * b1; i++) {
        cxeadd(C->data[ic * j + i], A->data[ia * j + 0], B->data[ib * j + i])
      }
  }
  // a0 a1 1 1
  else if (((a0 != 1) && (a1 != 1)) && ((b0 == 1) && (b1 == 1))) {
    if ((C->d0 != a0) || (C->d1 != a1)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < a0 * a1; i++) {
        cxeadd(C->data[ic * j + i], A->data[ia * j + i], B->data[ib * j + 0])
      }
  }
  // 1 a1 b0 1
  else if (((a0 == 1) && (a1 != 1)) && ((b0 != 1) && (b1 == 1))) {
    if ((C->d0 != b0) || (C->d1 != a1)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < b0 * a1; i++) {
        cxeadd(C->data[ic * j + i], A->data[ia * j + i / a1],
               B->data[ib * j + i % b0])
      }
  }
  // a0 1 1 b1
  else if (((a0 != 1) && (a1 == 1)) && ((b0 == 1) && (b1 != 1))) {
    if ((C->d0 != a0) || (C->d1 != b1)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < a0 * b1; i++) {
        cxeadd(C->data[ic * j + i], A->data[ia * j + i % a0],
               B->data[ib * j + i / b1])
      }
  }
  // 1  a1 b0 b1 | b0 a1 (a1==b1)
  else if (((a0 == 1) && (a1 != 1)) && ((b0 != 1) && (b1 != 1))) {
    if (((C->d0 != b0) || (C->d1 != a1)) || (a1 != b1)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < b0 * a1; i++) {
        cxeadd(C->data[ic * j + i], A->data[ia * j + i / b1],
               B->data[ib * j + i])
      }
  }
  // a0 1  b0 b1 | b0 b1 (a0==b0)
  else if (((a0 != 1) && (a1 == 1)) && ((b0 != 1) && (b1 != 1))) {
    if (((C->d0 != b0) || (C->d1 != b1)) || (a0 != b0)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < b0 * b1; i++) {
        cxeadd(C->data[ic * j + i], A->data[ia * j + i % b0],
               B->data[ib * j + i])
      }
  }
  // a0 a1 b0  1 | b0 a1 (a0==b0)
  else if (((a0 != 1) && (a1 != 1)) && ((b0 != 1) && (b1 == 1))) {
    if (((C->d0 != b0) || (C->d1 != a1)) || (a0 != b0)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < b0 * a1; i++) {
        cxeadd(C->data[ic * j + i], A->data[ia * j + i],
               B->data[ib * j + i % a1])
      }
  }
  // a0 a1 1  b1 | a0 b1 (a1==b1)
  else if (((a0 != 1) && (a1 != 1)) && ((b0 == 1) && (b1 != 1))) {
    if (((C->d0 != a0) || (C->d1 != b1)) || (a1 != b1)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < a0 * b1; i++) {
        cxeadd(C->data[ic * j + i], A->data[ia * j + i],
               B->data[ib * j + i / a0])
      }
  } else
    ASSERT(DIM_INVAL)
}

/**** mul_elements ****/
void mul_elements(MAT *A, MAT *B, MAT *C) {
  ITER i, j;
  UINT a0, a1, b0, b1;
  UINT a2, b2, c2;
  UINT ia, ib, ic;
#if DEBUG
  printf("%s\n", __func__);
#endif

  a2 = A->d2;
  b2 = B->d2;
  c2 = C->d2;
  a0 = A->d0;
  a1 = A->d1;
  b0 = B->d0;
  b1 = B->d1;

  ic = C->d0 * C->d1;
  /*1 to 1*/
  if (a2 == b2) {
    if (c2 != a2) ASSERT(DIM_INVAL)
    ia = a0 * a1;
    ib = b0 * b1;
  }
  /*1 to b2*/
  else if (a2 == 1 && b2 != 1) {
    if (b2 != c2) ASSERT(DIM_INVAL)
    ia = 0;
    ib = b0 * b1;
  }
  /* 1 to a2 */
  else if (a2 != 1 && b2 == 1) {
    if (a2 != c2) ASSERT(DIM_INVAL)
    ib = 0;
    ia = a0 * a1;
  }
  /* unmatching*/
  else
    ASSERT(DIM_INVAL)

  /* 12 broadcasting cases..
   * a0 a1 b0 b1 | c0 c1
   * 1  1  1  b1 | 1  b1
   * 1  1  b0 1  | b0 1
   * 1  a1 1  1  | 1  a1
   * a0 1  1  1  | a0 1
   *
   * 1  1  b0 b1 | b0 b1
   * a0 a1 1  1  | a0 a1
   * 1  a1 b0 1  | b0 a1
   * a0 1  1  b1 | a0 b1
   *
   * 1  a1 b0 b1 | b0 a1 (a1==b1)
   * a0 1  b0 b1 | b0 b1 (a0==b0)
   * a0 a1 b0  1 | b0 a1 (a0==b0)
   * a0 a1 1  b1 | a0 b1 (a1==b1)
   */
  // 1 1 1 b1

  if (((a0 == 1) && (a1 == 1)) && ((b0 == 1) && (b1 != 1))) {
    if ((C->d0 != 1) || (C->d1 != b1)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < b1; i++)
        C->data[ic * j + i] =
            B->data[ib * j + ib * j + i] * A->data[ia * j + ia * j + 0];
  }
  // 1 1 b0 1
  else if (((a0 == 1) && (a1 == 1)) && ((b0 != 1) && (b1 == 1))) {
    if ((C->d0 != b0) || (C->d1 != 1)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < b1; i++)
        C->data[ic * j + i] =
            B->data[ib * j + ib * j + i] * A->data[ia * j + ia * j + 0];
  }
  // 1 a1 1 1
  else if (((a0 == 1) && (a1 != 1)) && ((b0 == 1) && (b1 == 1))) {
    if ((C->d0 != 1) || (C->d1 != a1)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < b1; i++)
        C->data[ic * j + i] =
            A->data[ia * j + ia * j + i] * B->data[ib * j + ib * j + 0];
  }
  // a0 1 1 1
  else if (((a0 != 1) && (a1 == 1)) && ((b0 == 1) && (b1 == 1))) {
    if ((C->d0 != a0) || (C->d1 != 1)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < b1; i++)
        C->data[ic * j + i] = A->data[ia * j + i] * B->data[ib * j + 0];
  }
  // 1 1 b0 b1
  else if (((a0 == 1) && (a1 == 1)) && ((b0 != 1) && (b1 != 1))) {
    if ((C->d0 != b0) || (C->d1 != b1)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < b0 * b1; i++)
        C->data[ic * j + i] = B->data[ib * j + i] * A->data[ia * j + 0];
  }
  // a0 a1 1 1
  else if (((a0 != 1) && (a1 != 1)) && ((b0 == 1) && (b1 == 1))) {
    if ((C->d0 != a0) || (C->d1 != a1)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < a0 * a1; i++)
        C->data[ic * j + i] = A->data[ia * j + i] * B->data[ib * j + 0];
  }
  // 1 a1 b0 1
  else if (((a0 == 1) && (a1 != 1)) && ((b0 != 1) && (b1 == 1))) {
    if ((C->d0 != b0) || (C->d1 != a1)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < b0 * a1; i++)
        C->data[ic * j + i] =
            B->data[ib * j + i % b0] * A->data[ia * j + i / a1];
  }
  // a0 1 1 b1
  else if (((a0 != 1) && (a1 == 1)) && ((b0 == 1) && (b1 != 1))) {
    if ((C->d0 != a0) || (C->d1 != b1)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < a0 * b1; i++)
        C->data[ic * j + i] =
            B->data[ib * j + i / b1] * A->data[ia * j + i % a0];
  }
  // 1  a1 b0 b1 | b0 a1 (a1==b1)
  else if (((a0 == 1) && (a1 != 1)) && ((b0 != 1) && (b1 != 1))) {
    if (((C->d0 != b0) || (C->d1 != a1)) || (a1 != b1)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < b0 * a1; i++)
        C->data[ic * j + i] = A->data[ia * j + i / b1] * B->data[ib * j + i];
  }
  // a0 1  b0 b1 | b0 b1 (a0==b0)
  else if (((a0 != 1) && (a1 == 1)) && ((b0 != 1) && (b1 != 1))) {
    if (((C->d0 != b0) || (C->d1 != b1)) || (a0 != b0)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < b0 * b1; i++)
        C->data[ic * j + i] = A->data[ia * j + i % b0] * B->data[ib * j + i];
  }
  // a0 a1 b0  1 | b0 a1 (a0==b0)
  else if (((a0 != 1) && (a1 != 1)) && ((b0 != 1) && (b1 == 1))) {
    if (((C->d0 != b0) || (C->d1 != a1)) || (a0 != b0)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < b0 * a1; i++)
        C->data[ic * j + i] = B->data[ib * j + i % a1] * A->data[ia * j + i];
  }
  // a0 a1 1  b1 | a0 b1 (a1==b1)
  else if (((a0 != 1) && (a1 != 1)) && ((b0 == 1) && (b1 != 1))) {
    if (((C->d0 != a0) || (C->d1 != b1)) || (a1 != b1)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < a0 * b1; i++)
        C->data[ic * j + i] = B->data[ib * j + i / a0] * A->data[ia * j + i];
  } else
    ASSERT(DIM_INVAL)
}

void cmul_elements(CMAT *A, CMAT *B, CMAT *C) {
  ITER i, j;
  UINT a0, a1, b0, b1;
  UINT a2, b2, c2;
  UINT ia, ib, ic;
#if DEBUG
  printf("%s\n", __func__);
#endif

  a2 = A->d2;
  b2 = B->d2;
  c2 = C->d2;
  a0 = A->d0;
  a1 = A->d1;
  b0 = B->d0;
  b1 = B->d1;

  ic = C->d0 * C->d1;
  /*1 to 1*/
  if (a2 == b2) {
    if (c2 != a2) ASSERT(DIM_INVAL)
    ia = a0 * a1;
    ib = b0 * b1;
  }
  /*1 to b2*/
  else if (a2 == 1 && b2 != 1) {
    if (b2 != c2) ASSERT(DIM_INVAL)
    ia = 0;
    ib = b0 * b1;
  }
  /* 1 to a2 */
  else if (a2 != 1 && b2 == 1) {
    if (a2 != c2) ASSERT(DIM_INVAL)
    ib = 0;
    ia = a0 * a1;
  }
  /* unmatching*/
  else
    ASSERT(DIM_INVAL)

  // 1 1 1 b1
  if (((a0 == 1) && (a1 == 1)) && ((b0 == 1) && (b1 != 1))) {
    if ((C->d0 != 1) || (C->d1 != b1)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < b1; i++)
        cxemul(C->data[i], A->data[ia * j + 0], B->data[ib * j + i])
  }
  // 1 1 b0 1
  else if (((a0 == 1) && (a1 == 1)) && ((b0 != 1) && (b1 == 1))) {
    if ((C->d0 != b0) || (C->d1 != 1)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < b1; i++)
        cxemul(C->data[ic * j + i], A->data[ia * j + 0], B->data[ib * j + i])
  }
  // 1 a1 1 1
  else if (((a0 == 1) && (a1 != 1)) && ((b0 == 1) && (b1 == 1))) {
    if ((C->d0 != 1) || (C->d1 != a1)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < b1; i++)
        cxemul(C->data[ic * j + i], A->data[ia * j + i], B->data[ib * j + 0])
  }
  // a0 1 1 1
  else if (((a0 != 1) && (a1 == 1)) && ((b0 == 1) && (b1 == 1))) {
    if ((C->d0 != a0) || (C->d1 != 1)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < b1; i++)
        cxemul(C->data[ic * j + i], A->data[ia * j + i], B->data[ib * j + 0])
  }
  // 1 1 b0 b1
  else if (((a0 == 1) && (a1 == 1)) && ((b0 != 1) && (b1 != 1))) {
    if ((C->d0 != b0) || (C->d1 != b1)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < b0 * b1; i++)
        cxemul(C->data[ic * j + i], A->data[ia * j + 0], B->data[ib * j + i])
  }
  // a0 a1 1 1
  else if (((a0 != 1) && (a1 != 1)) && ((b0 == 1) && (b1 == 1))) {
    if ((C->d0 != a0) || (C->d1 != a1)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < a0 * a1; i++)
        cxemul(C->data[ic * j + i], A->data[ia * j + i], B->data[ib * j + 0])
  }
  // 1 a1 b0 1
  else if (((a0 == 1) && (a1 != 1)) && ((b0 != 1) && (b1 == 1))) {
    if ((C->d0 != b0) || (C->d1 != a1)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < b0 * a1; i++)
        cxemul(C->data[ic * j + i], A->data[ia * j + i / a1],
               B->data[ib * j + i % b0])
  }
  // a0 1 1 b1
  else if (((a0 != 1) && (a1 == 1)) && ((b0 == 1) && (b1 != 1))) {
    if ((C->d0 != a0) || (C->d1 != b1)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < a0 * b1; i++)
        cxemul(C->data[ic * j + i], A->data[ia * j + i % a0],
               B->data[ib * j + i / b1])
  }
  // 1  a1 b0 b1 | b0 a1 (a1==b1)
  else if (((a0 == 1) && (a1 != 1)) && ((b0 != 1) && (b1 != 1))) {
    if (((C->d0 != b0) || (C->d1 != a1)) || (a1 != b1)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < b0 * a1; i++)
        cxemul(C->data[ic * j + i], A->data[ia * j + i / b1],
               B->data[ib * j + i])
  }
  // a0 1  b0 b1 | b0 b1 (a0==b0)
  else if (((a0 != 1) && (a1 == 1)) && ((b0 != 1) && (b1 != 1))) {
    if (((C->d0 != b0) || (C->d1 != b1)) || (a0 != b0)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < b0 * b1; i++)
        cxemul(C->data[ic * j + i], A->data[ia * j + i % b0],
               B->data[ib * j + i])
  }
  // a0 a1 b0  1 | b0 a1 (a0==b0)
  else if (((a0 != 1) && (a1 != 1)) && ((b0 != 1) && (b1 == 1))) {
    if (((C->d0 != b0) || (C->d1 != a1)) || (a0 != b0)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < b0 * a1; i++)
        cxemul(C->data[ic * j + i], A->data[ia * j + i],
               B->data[ib * j + i % a1])
  }
  // a0 a1 1  b1 | a0 b1 (a1==b1)
  else if (((a0 != 1) && (a1 != 1)) && ((b0 == 1) && (b1 != 1))) {
    if (((C->d0 != a0) || (C->d1 != b1)) || (a1 != b1)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < a0 * b1; i++)
        cxemul(C->data[ic * j + i], A->data[ia * j + i],
               B->data[ib * j + i / a0])
  } else
    ASSERT(DIM_INVAL)
}

/**** div_elements - broadcasting operation ****/
void div_elements(MAT *A, MAT *B, MAT *C) {
  ITER i, j;
  UINT a0, a1, b0, b1;
  UINT a2, b2, c2;
  UINT ia, ib, ic;
#if DEBUG
  printf("%s\n", __func__);
#endif

  a2 = A->d2;
  b2 = B->d2;
  c2 = C->d2;
  a0 = A->d0;
  a1 = A->d1;
  b0 = B->d0;
  b1 = B->d1;

  ic = C->d0 * C->d1;
  /*1 to 1*/
  if (a2 == b2) {
    if (c2 != a2) ASSERT(DIM_INVAL)
    ia = a0 * a1;
    ib = b0 * b1;
  }
  /*1 to b2*/
  else if (a2 == 1 && b2 != 1) {
    if (b2 != c2) ASSERT(DIM_INVAL)
    ia = 0;
    ib = b0 * b1;
  }
  /* 1 to a2 */
  else if (a2 != 1 && b2 == 1) {
    if (a2 != c2) ASSERT(DIM_INVAL)
    ib = 0;
    ia = a0 * a1;
  }
  /* unmatching*/
  else
    ASSERT(DIM_INVAL)

  /* 12 broadcasting cases..
   * a0 a1 b0 b1 | c0 c1
   * 1  1  1  b1 | 1  b1
   * 1  1  b0 1  | b0 1
   * 1  a1 1  1  | 1  a1
   * a0 1  1  1  | a0 1
   *
   * 1  1  b0 b1 | b0 b1
   * a0 a1 1  1  | a0 a1
   * 1  a1 b0 1  | b0 a1
   * a0 1  1  b1 | a0 b1
   *
   * 1  a1 b0 b1 | b0 a1 (a1==b1)
   * a0 1  b0 b1 | b0 b1 (a0==b0)
   * a0 a1 b0  1 | b0 a1 (a0==b0)
   * a0 a1 1  b1 | a0 b1 (a1==b1)
   */
  // 1 1 1 b1
  if (((a0 == 1) && (a1 == 1)) && ((b0 == 1) && (b1 != 1))) {
    if ((C->d0 != 1) || (C->d1 != b1)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < b1; i++) {
        C->data[ic * j + i] = A->data[ia * j + 0] / B->data[ib * j + i];
      }
  }
  // 1 1 b0 1
  else if (((a0 == 1) && (a1 == 1)) && ((b0 != 1) && (b1 == 1))) {
    if ((C->d0 != b0) || (C->d1 != 1)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < b1; i++) {
        C->data[ic * j + i] = A->data[ia * j + 0] / B->data[ib * j + i];
      }
  }
  // 1 a1 1 1
  else if (((a0 == 1) && (a1 != 1)) && ((b0 == 1) && (b1 == 1))) {
    if ((C->d0 != 1) || (C->d1 != a1)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < b1; i++) {
        C->data[ic * j + i] = A->data[ia * j + i] / B->data[ib * j + 0];
      }
  }
  // a0 1 1 1
  else if (((a0 != 1) && (a1 == 1)) && ((b0 == 1) && (b1 == 1))) {
    if ((C->d0 != a0) || (C->d1 != 1)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < b1; i++) {
        C->data[ic * j + i] = A->data[ia * j + i] / B->data[ib * j + 0];
      }
  }
  // 1 1 b0 b1
  else if (((a0 == 1) && (a1 == 1)) && ((b0 != 1) && (b1 != 1))) {
    if ((C->d0 != b0) || (C->d1 != b1)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < b0 * b1; i++) {
        C->data[ic * j + i] = A->data[ia * j + 0] / B->data[ib * j + i];
      }
  }
  // a0 a1 1 1
  else if (((a0 != 1) && (a1 != 1)) && ((b0 == 1) && (b1 == 1))) {
    if ((C->d0 != a0) || (C->d1 != a1)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < a0 * a1; i++) {
        C->data[ic * j + i] = A->data[ia * j + i] / B->data[ib * j + 0];
      }
  }
  // 1 a1 b0 1
  else if (((a0 == 1) && (a1 != 1)) && ((b0 != 1) && (b1 == 1))) {
    if ((C->d0 != b0) || (C->d1 != a1)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < b0 * a1; i++) {
        C->data[ic * j + i] =
            A->data[ia * j + i / a1] / B->data[ib * j + i % b0];
      }
  }
  // a0 1 1 b1
  else if (((a0 != 1) && (a1 == 1)) && ((b0 == 1) && (b1 != 1))) {
    if ((C->d0 != a0) || (C->d1 != b1)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < a0 * b1; i++) {
        C->data[ic * j + i] =
            A->data[ia * j + i % a0] / B->data[ib * j + i / b1];
      }
  }
  // 1  a1 b0 b1 | b0 a1 (a1==b1)
  else if (((a0 == 1) && (a1 != 1)) && ((b0 != 1) && (b1 != 1))) {
    if (((C->d0 != b0) || (C->d1 != a1)) || (a1 != b1)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < b0 * a1; i++) {
        C->data[ic * j + i] = A->data[ia * j + i / b1] / B->data[ib * j + i];
      }
  }
  // a0 1  b0 b1 | b0 b1 (a0==b0)
  else if (((a0 != 1) && (a1 == 1)) && ((b0 != 1) && (b1 != 1))) {
    if (((C->d0 != b0) || (C->d1 != b1)) || (a0 != b0)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < b0 * b1; i++) {
        C->data[ic * j + i] = A->data[ia * j + i % b0] / B->data[ib * j + i];
      }
  }
  // a0 a1 b0  1 | b0 a1 (a0==b0)
  else if (((a0 != 1) && (a1 != 1)) && ((b0 != 1) && (b1 == 1))) {
    if (((C->d0 != b0) || (C->d1 != a1)) || (a0 != b0)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < b0 * a1; i++) {
        C->data[ic * j + i] = A->data[ia * j + i] / B->data[ib * j + i % a1];
      }
  }
  // a0 a1 1  b1 | a0 b1 (a1==b1)
  else if (((a0 != 1) && (a1 != 1)) && ((b0 == 1) && (b1 != 1))) {
    if (((C->d0 != a0) || (C->d1 != b1)) || (a1 != b1)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < a0 * b1; i++) {
        C->data[ic * j + i] = A->data[ia * j + i] / B->data[ib * j + i / a0];
      }
  } else
    ASSERT(DIM_INVAL)
}

void cdiv_elements(CMAT *A, CMAT *B, CMAT *C) {
  ITER i, j;
  UINT a0, a1, b0, b1;
  UINT a2, b2, c2;
  UINT ia, ib, ic;
#if DEBUG
  printf("%s\n", __func__);
#endif

  a2 = A->d2;
  b2 = B->d2;
  c2 = C->d2;
  a0 = A->d0;
  a1 = A->d1;
  b0 = B->d0;
  b1 = B->d1;

  ic = C->d0 * C->d1;
  /*1 to 1*/
  if (a2 == b2) {
    if (c2 != a2) ASSERT(DIM_INVAL)
    ia = a0 * a1;
    ib = b0 * b1;
  }
  /*1 to b2*/
  else if (a2 == 1 && b2 != 1) {
    if (b2 != c2) ASSERT(DIM_INVAL)
    ia = 0;
    ib = b0 * b1;
  }
  /* 1 to a2 */
  else if (a2 != 1 && b2 == 1) {
    if (a2 != c2) ASSERT(DIM_INVAL)
    ib = 0;
    ia = a0 * a1;
  }
  /* unmatching*/
  else
    ASSERT(DIM_INVAL)

  /* 12 broadcasting cases..
   * a0 a1 b0 b1 | c0 c1
   * 1  1  1  b1 | 1  b1
   * 1  1  b0 1  | b0 1
   * 1  a1 1  1  | 1  a1
   * a0 1  1  1  | a0 1
   *
   * 1  1  b0 b1 | b0 b1
   * a0 a1 1  1  | a0 a1
   * 1  a1 b0 1  | b0 a1
   * a0 1  1  b1 | a0 b1
   *
   * 1  a1 b0 b1 | b0 a1 (a1==b1)
   * a0 1  b0 b1 | b0 b1 (a0==b0)
   * a0 a1 b0  1 | b0 a1 (a0==b0)
   * a0 a1 1  b1 | a0 b1 (a1==b1)
   */
  // 1 1 1 b1
  if (((a0 == 1) && (a1 == 1)) && ((b0 == 1) && (b1 != 1))) {
    if ((C->d0 != 1) || (C->d1 != b1)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < b1; i++) {
        cxediv(C->data[ic * j + i], A->data[ia * j + 0], B->data[ib * j + i])
      }
  }
  // 1 1 b0 1
  else if (((a0 == 1) && (a1 == 1)) && ((b0 != 1) && (b1 == 1))) {
    if ((C->d0 != b0) || (C->d1 != 1)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < b1; i++) {
        cxediv(C->data[ic * j + i], A->data[ia * j + 0], B->data[ib * j + i])
      }
  }
  // 1 a1 1 1
  else if (((a0 == 1) && (a1 != 1)) && ((b0 == 1) && (b1 == 1))) {
    if ((C->d0 != 1) || (C->d1 != a1)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < b1; i++) {
        cxediv(C->data[ic * j + i], A->data[ia * j + i], B->data[ib * j + 0])
      }
  }
  // a0 1 1 1
  else if (((a0 != 1) && (a1 == 1)) && ((b0 == 1) && (b1 == 1))) {
    if ((C->d0 != a0) || (C->d1 != 1)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < b1; i++) {
        cxediv(C->data[ic * j + i], A->data[ia * j + i], B->data[ib * j + 0])
      }
  }
  // 1 1 b0 b1
  else if (((a0 == 1) && (a1 == 1)) && ((b0 != 1) && (b1 != 1))) {
    if ((C->d0 != b0) || (C->d1 != b1)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < b0 * b1; i++) {
        cxediv(C->data[ic * j + i], A->data[ia * j + 0], B->data[ib * j + i])
      }
  }
  // a0 a1 1 1
  else if (((a0 != 1) && (a1 != 1)) && ((b0 == 1) && (b1 == 1))) {
    if ((C->d0 != a0) || (C->d1 != a1)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < a0 * a1; i++) {
        cxediv(C->data[ic * j + i], A->data[ia * j + i], B->data[ib * j + 0])
      }
  }
  // 1 a1 b0 1
  else if (((a0 == 1) && (a1 != 1)) && ((b0 != 1) && (b1 == 1))) {
    if ((C->d0 != b0) || (C->d1 != a1)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < b0 * a1; i++) {
        cxediv(C->data[ic * j + i], A->data[ia * j + i / a1],
               B->data[ib * j + i % b0])
      }
  }
  // a0 1 1 b1
  else if (((a0 != 1) && (a1 == 1)) && ((b0 == 1) && (b1 != 1))) {
    if ((C->d0 != a0) || (C->d1 != b1)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < a0 * b1; i++) {
        cxediv(C->data[ic * j + i], A->data[ia * j + i % a0],
               B->data[ib * j + i / b1])
      }
  }
  // 1  a1 b0 b1 | b0 a1 (a1==b1)
  else if (((a0 == 1) && (a1 != 1)) && ((b0 != 1) && (b1 != 1))) {
    if (((C->d0 != b0) || (C->d1 != a1)) || (a1 != b1)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < b0 * a1; i++) {
        cxediv(C->data[ic * j + i], A->data[ia * j + i / b1],
               B->data[ib * j + i])
      }
  }
  // a0 1  b0 b1 | b0 b1 (a0==b0)
  else if (((a0 != 1) && (a1 == 1)) && ((b0 != 1) && (b1 != 1))) {
    if (((C->d0 != b0) || (C->d1 != b1)) || (a0 != b0)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < b0 * b1; i++) {
        cxediv(C->data[ic * j + i], A->data[ia * j + i % b0],
               B->data[ib * j + i])
      }
  }
  // a0 a1 b0  1 | b0 a1 (a0==b0)
  else if (((a0 != 1) && (a1 != 1)) && ((b0 != 1) && (b1 == 1))) {
    if (((C->d0 != b0) || (C->d1 != a1)) || (a0 != b0)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < b0 * a1; i++) {
        cxediv(C->data[ic * j + i], A->data[ia * j + i],
               B->data[ib * j + i % a1])
      }
  }
  // a0 a1 1  b1 | a0 b1 (a1==b1)
  else if (((a0 != 1) && (a1 != 1)) && ((b0 == 1) && (b1 != 1))) {
    if (((C->d0 != a0) || (C->d1 != b1)) || (a1 != b1)) ASSERT(DIM_INVAL)
    for (j = 0; j < c2; j++)
#pragma omp parallel for shared(C, B, A) private(i)
      for (i = 0; i < a0 * b1; i++) {
        cxediv(C->data[ic * j + i], A->data[ia * j + i],
               B->data[ib * j + i / a0])
      }
  } else
    ASSERT(DIM_INVAL)
}

/**** inverse elements ****/
void inv_elements(MAT *mat) {
  inv_elements_inc(mat->d0 * mat->d1 * mat->d2, mat->data, 1);
}

void inv_elements_inc(UINT size, DTYPE *X, ITER incx) {
  ITER i;
#if DEBUG
  printf("%s\n", __func__);
#endif
#pragma omp parallel for shared(X) private(i)
  for (i = 0; i < size; i += incx) X[i] = 1 / X[i];
}

void cinv_elements(CMAT *mat) {
  cinv_elements_inc(mat->d0 * mat->d1 * mat->d2, mat->data, 1);
}

void cinv_elements_inc(UINT size, CTYPE *X, ITER incx) {
  ITER i;
#if DEBUG
  printf("%s\n", __func__);
#endif
#pragma omp parallel for shared(X) private(i)
  for (i = 0; i < size; i += incx) {
    X[i].re = X[i].re / ((X[i].re * X[i].re) + (X[i].im * X[i].im));
    X[i].im = -X[i].im / ((X[i].re * X[i].re) + (X[i].im * X[i].im));
  }
}

/**** repeat matrix ****/
void repmat(MAT *mat, DIM *dim) {
  MAT *t_mat;
  ITER i, j, k, l;
  UINT size, idx;
  UINT v, h, b;  // vertical, horizontal, blockal
#if DEBUG
  printf("%s\n", __func__);
#endif

  t_mat = mem_MAT(mat->d0, mat->d1, mat->d2);
  copy(mat, t_mat);
  free_MAT(mat);
  mat =
      alloc_MAT(t_mat->d0 * dim->d0, t_mat->d1 * dim->d1, t_mat->d2 * dim->d2);

  size = t_mat->d0 * t_mat->d1;

#pragma omp parallel for shared(mat, t_mat) private(i, l, k, j)
  for (i = 0; i < size; i++) {
    for (l = 0; l < dim->d2; l++)
      for (k = 0; k < dim->d1; k++)
        for (j = 0; j < dim->d0; j++) {
          idx = (l * mat->d0 * mat->d1) + (k * t_mat->d1 * mat->d0) +
                (j * t_mat->d0) + (i / t_mat->d0) * mat->d0 + i % t_mat->d0;
          mat->data[idx] = t_mat->data[i];
        }
  }

  free_mem_MAT(t_mat);
}

void crepmat(CMAT *mat, DIM *dim) {
  CMAT *t_mat;
  ITER i, j, k, l;
  UINT size, idx;
  UINT v, h, b;  // vertical, horizontal, blockal
#if DEBUG
  printf("%s\n", __func__);
#endif

  t_mat = mem_CMAT(mat->d0, mat->d1, mat->d2);
  ccopy(mat, t_mat);
  free_CMAT(mat);
  mat =
      alloc_CMAT(t_mat->d0 * dim->d0, t_mat->d1 * dim->d1, t_mat->d2 * dim->d2);

  size = t_mat->d0 * t_mat->d1;

#pragma omp parallel for shared(mat, t_mat) private(i, l, k, j)
  for (i = 0; i < size; i++) {
    for (l = 0; l < dim->d2; l++)
      for (k = 0; k < dim->d1; k++)
        for (j = 0; j < dim->d0; j++) {
          idx = (l * mat->d0 * mat->d1) + (k * t_mat->d1 * mat->d0) +
                (j * t_mat->d0) + (i / t_mat->d0) * mat->d0 + i % t_mat->d0;
          mat->data[idx].re = t_mat->data[i].re;
          mat->data[idx].im = t_mat->data[i].im;
        }
  }

  free_mem_CMAT(t_mat);
}
/**** reshape matrix ****/
void reshape(MAT *mat, DIM *dim) {
  UINT ori, neo;
#if DEBUG
  printf("%s\n", __func__);
#endif
  ori = mat->d0 * mat->d1 * mat->d2;
  neo = dim->d0 * dim->d1 * dim->d2;

  if (ori != neo) ASSERT(DIM_INVAL)

  mat->d0 = dim->d0;
  mat->d1 = dim->d1;
  mat->d2 = dim->d2;
}
void creshape(CMAT *mat, DIM *dim) {
  UINT ori, neo;
#if DEBUG
  printf("%s\n", __func__);
#endif
  ori = mat->d0 * mat->d1 * mat->d2;
  neo = dim->d0 * dim->d1 * dim->d2;

  if (ori != neo) ASSERT(DIM_INVAL)

  mat->d0 = dim->d0;
  mat->d1 = dim->d1;
  mat->d2 = dim->d2;
}

/**** shift dimension ****/

void shiftdim(MAT *mat, SINT n) {
  /*어차피 3차원만있어서 빙글빙글 **/
  if (n < 0) {
    n = -n % 2 + 1;
  }
  n = n % 2;
  if (n == 1)
    permute(mat, 231);
  else
    permute(mat, 312);
}

void cshiftdim(CMAT *mat, SINT n) {
  /*어차피 3차원만있어서 빙글빙글 **/
  if (n < 0) {
    n = -n % 2 + 1;
  }
  n = n % 2;
  if (n == 1)
    cpermute(mat, 231);
  else
    cpermute(mat, 312);
}
/**** permutate ****/

/*
 * d2 = i/d0*d1
 * d1 = (i%(d0*d1))/d0
 * d0 = i % d0
 *
 * */
void permute(MAT *mat, UINT seq) {
  ITER i;
  MAT *t;
  UINT d0d1;
#if DEBUG
  printf("%s\n", __func__);
#endif

  /** Nothing to do**/
  if (seq == 123) return;

  if (seq == 132) {
    t = mem_MAT(mat->d0, mat->d1, mat->d2);
    copy(mat, t);
    d0d1 = t->d0 * t->d1;
    mat->d1 = t->d2;
    mat->d2 = t->d1;
#pragma omp parallel for shared(mat, t) private(i)
    for (i = 0; i < mat->d0 * mat->d1 * mat->d2; i++) {
      mat->data[((i % d0d1) / t->d0) * mat->d0 * mat->d1 +
                (i / d0d1) * mat->d0 + i % t->d0] = t->data[i];
    }
    free_mem_MAT(t);
  } else if (seq == 213) {
    t = mem_MAT(mat->d0, mat->d1, mat->d2);
    copy(mat, t);
    d0d1 = t->d0 * t->d1;
    mat->d0 = t->d1;
    mat->d1 = t->d0;
#pragma omp parallel for shared(mat, t) private(i)
    for (i = 0; i < mat->d0 * mat->d1 * mat->d2; i++) {
      mat->data[((i / d0d1)) * mat->d0 * mat->d1 + (i % t->d0) * mat->d0 +
                (i % d0d1) / t->d0] = t->data[i];
    }
    free_mem_MAT(t);
  } else if (seq == 231) {
    t = mem_MAT(mat->d0, mat->d1, mat->d2);
    copy(mat, t);
    d0d1 = t->d0 * t->d1;
    mat->d0 = t->d1;
    mat->d1 = t->d2;
    mat->d2 = t->d0;
#pragma omp parallel for shared(mat, t) private(i)
    for (i = 0; i < mat->d0 * mat->d1 * mat->d2; i++) {
      mat->data[(i % t->d0) * mat->d0 * mat->d1 + (i / d0d1) * mat->d0 +
                (i % d0d1) / t->d0] = t->data[i];
    }
    free_mem_MAT(t);
  } else if (seq == 312) {
    t = mem_MAT(mat->d0, mat->d1, mat->d2);
    copy(mat, t);
    d0d1 = t->d0 * t->d1;
    mat->d0 = t->d2;
    mat->d2 = t->d0;
#pragma omp parallel for shared(mat, t) private(i)
    for (i = 0; i < mat->d0 * mat->d1 * mat->d2; i++) {
      mat->data[(i % d0d1) / t->d0 * mat->d0 * mat->d1 + (i % t->d0) * mat->d0 +
                i / d0d1] = t->data[i];
    }
    free_mem_MAT(t);
  } else if (seq == 321) {
    t = mem_MAT(mat->d0, mat->d1, mat->d2);
    copy(mat, t);
    d0d1 = t->d0 * t->d1;
    mat->d0 = t->d2;
    mat->d2 = t->d0;
#pragma omp parallel for shared(mat, t) private(i)
    for (i = 0; i < mat->d0 * mat->d1 * mat->d2; i++) {
      mat->data[(i % t->d0) * mat->d0 * mat->d1 + (i % d0d1) / t->d0 * mat->d0 +
                i / d0d1] = t->data[i];
    }
    free_mem_MAT(t);
  } else
    ASSERT(ARG_INVAL)
}

void cpermute(CMAT *mat, UINT seq) {
  ITER i;
  CMAT *t;
  UINT d0d1;
#if DEBUG
  printf("%s\n", __func__);
#endif

  /** Nothing to do**/
  if (seq == 123) return;

  if (seq == 132) {
    t = mem_CMAT(mat->d0, mat->d1, mat->d2);
    ccopy(mat, t);
    d0d1 = t->d0 * t->d1;
    mat->d1 = t->d2;
    mat->d2 = t->d1;
#pragma omp parallel for shared(mat, t) private(i)
    for (i = 0; i < mat->d0 * mat->d1 * mat->d2; i++) {
      mat->data[((i % d0d1) / t->d0) * mat->d0 * mat->d1 +
                (i / d0d1) * mat->d0 + i % t->d0]
          .re = t->data[i].re;
      mat->data[((i % d0d1) / t->d0) * mat->d0 * mat->d1 +
                (i / d0d1) * mat->d0 + i % t->d0]
          .im = t->data[i].im;
    }
    free_mem_CMAT(t);
  } else if (seq == 213) {
    t = mem_CMAT(mat->d0, mat->d1, mat->d2);
    ccopy(mat, t);
    d0d1 = t->d0 * t->d1;
    mat->d0 = t->d1;
    mat->d1 = t->d0;
#pragma omp parallel for shared(mat, t) private(i)
    for (i = 0; i < mat->d0 * mat->d1 * mat->d2; i++) {
      mat->data[((i / d0d1)) * mat->d0 * mat->d1 + (i % t->d0) * mat->d0 +
                (i % d0d1) / t->d0]
          .re = t->data[i].re;
      mat->data[((i / d0d1)) * mat->d0 * mat->d1 + (i % t->d0) * mat->d0 +
                (i % d0d1) / t->d0]
          .im = t->data[i].im;
    }
    free_mem_CMAT(t);
  } else if (seq == 231) {
    t = mem_CMAT(mat->d0, mat->d1, mat->d2);
    ccopy(mat, t);
    d0d1 = t->d0 * t->d1;
    mat->d0 = t->d1;
    mat->d1 = t->d2;
    mat->d2 = t->d0;
#pragma omp parallel for shared(mat, t) private(i)
    for (i = 0; i < mat->d0 * mat->d1 * mat->d2; i++) {
      mat->data[(i % t->d0) * mat->d0 * mat->d1 + (i / d0d1) * mat->d0 +
                (i % d0d1) / t->d0]
          .re = t->data[i].re;
      mat->data[(i % t->d0) * mat->d0 * mat->d1 + (i / d0d1) * mat->d0 +
                (i % d0d1) / t->d0]
          .im = t->data[i].im;
    }
    free_mem_CMAT(t);
  } else if (seq == 312) {
    t = mem_CMAT(mat->d0, mat->d1, mat->d2);
    ccopy(mat, t);
    d0d1 = t->d0 * t->d1;
    mat->d0 = t->d2;
    mat->d2 = t->d0;
#pragma omp parallel for shared(mat, t) private(i)
    for (i = 0; i < mat->d0 * mat->d1 * mat->d2; i++) {
      mat->data[(i % d0d1) / t->d0 * mat->d0 * mat->d1 + (i % t->d0) * mat->d0 +
                i / d0d1]
          .re = t->data[i].re;
      mat->data[(i % d0d1) / t->d0 * mat->d0 * mat->d1 + (i % t->d0) * mat->d0 +
                i / d0d1]
          .im = t->data[i].im;
    }
    free_mem_CMAT(t);
  } else if (seq == 321) {
    t = mem_CMAT(mat->d0, mat->d1, mat->d2);
    ccopy(mat, t);
    d0d1 = t->d0 * t->d1;
    mat->d0 = t->d2;
    mat->d2 = t->d0;
#pragma omp parallel for shared(mat, t) private(i)
    for (i = 0; i < mat->d0 * mat->d1 * mat->d2; i++) {
      mat->data[(i % t->d0) * mat->d0 * mat->d1 + (i % d0d1) / t->d0 * mat->d0 +
                i / d0d1]
          .re = t->data[i].re;
      mat->data[(i % t->d0) * mat->d0 * mat->d1 + (i % d0d1) / t->d0 * mat->d0 +
                i / d0d1]
          .im = t->data[i].im;
    }
    free_mem_CMAT(t);
  } else
    ASSERT(ARG_INVAL)
}
/**** transpose ****/
MAT *create_trans(MAT *mat) {
  ITER i, j;
  UINT d0 = mat->d0;
  UINT d1 = mat->d1;
  UINT d2 = mat->d2;

  MAT *t_mat;
#if DEBUG
  printf("%s\n", __func__);
#endif
  if (mat->ndim == 0) {
    t_mat = alloc_MAT(d1, d0);
#pragma omp parallel for shared(t_mat, mat) private(i)
    for (i = 0; i < d0; i++) {
      t_mat->data[i] = mat->data[i];
    }
  } else if (mat->ndim == 1) {
    t_mat = alloc_MAT(d1, d0);
#pragma omp parallel for shared(t_mat, mat) private(i)
    for (i = 0; i < d0 * d1; i++) {
      t_mat->data[i / d0 + i % d0 * d1] = mat->data[i];
    }
  } else {
    t_mat = alloc_MAT(d1, d0, d2);
//일단은 3차원 배열을 2차원 배열의 batch로 여기고 구현함
#pragma omp parallel for shared(t_mat, mat) private(i, j)
    for (j = 0; j < d2; j++) {
      for (i = 0; i < d0 * d1; i++) {
        t_mat->data[j * d0 * d1 + i / d0 + i % d0 * d1] =
            mat->data[j * d0 * d1 + i];
      }
    }
  }
  return t_mat;
}

CMAT *create_ctrans(CMAT *mat) {
  ITER i, j;
  UINT d0 = mat->d0;
  UINT d1 = mat->d1;
  UINT d2 = mat->d2;

  CMAT *t_mat;
#if DEBUG
  printf("%s\n", __func__);
#endif
  if (mat->ndim == 0) {
    t_mat = alloc_CMAT(d1, d0);
#pragma omp parallel for shared(t_mat, mat) private(i)
    for (i = 0; i < d0; i++) {
      t_mat->data[i].re = mat->data[i].re;
      t_mat->data[i].im = mat->data[i].im;
    }
  } else if (mat->ndim == 1) {
    t_mat = alloc_CMAT(d1, d0);
#pragma omp parallel for shared(t_mat, mat) private(i)
    for (i = 0; i < d0 * d1; i++) {
      t_mat->data[i / d0 + i % d0 * d1].re = mat->data[i].re;
      t_mat->data[i / d0 + i % d0 * d1].im = mat->data[i].im;
    }
  } else {
    t_mat = alloc_CMAT(d1, d0, d2);
//일단은 3차원 배열을 2차원 배열의 batch로 여기고 구현함
#pragma omp parallel for shared(t_mat, mat) private(i, j)
    for (j = 0; j < d2; j++) {
      for (i = 0; i < d0 * d1; i++) {
        t_mat->data[j * d0 * d1 + i / d0 + i % d0 * d1].re =
            mat->data[j * d0 * d1 + i].re;
        t_mat->data[j * d0 * d1 + i / d0 + i % d0 * d1].im =
            mat->data[j * d0 * d1 + i].im;
      }
    }
  }
  return t_mat;
}

void trans(MAT *mat) {
  MAT *temp;
  ITER i, j, k;
  UINT d0, d1, d2;

#if DEBUG
  printf("%s\n", __func__);
#endif
  d0 = mat->d0;
  d1 = mat->d1;
  d2 = mat->d2;
  temp = mem_MAT(mat->d0, mat->d1, mat->d2);
  if (mat->ndim == 0) {
#pragma omp parallel for shared(temp, mat) private(i)
    for (i = 0; i < d0; i++) {
      temp->data[i] = mat->data[i];
    }
    mat->ndim = 1;
    mat->d0 = d1;
    mat->d1 = d0;
  } else if (mat->ndim == 1) {
#pragma omp parallel for shared(temp, mat) private(i)
    for (i = 0; i < d0 * d1; i++) {
      temp->data[i / d0 + i % d0 * d1] = mat->data[i];
    }
    mat->d0 = d1;
    mat->d1 = d0;
  } else {
// TODO 구현 덜됨
//일단은 3차원 배열을 2차원 배열의 batch로 여기고 구현함
#pragma omp parallel for shared(temp, mat) private(i, j)
    for (j = 0; j < d2; j++) {
      for (i = 0; i < d0 * d1; i++) {
        temp->data[j * d0 * d1 + i / d0 + i % d0 * d1] =
            mat->data[j * d0 * d1 + i];
      }
    }
    mat->d0 = d1;
    mat->d1 = d0;
  }
  copy(temp, mat);
  free_mem_MAT(temp);
}

void ctrans(CMAT *mat) {
  CMAT *temp;
  ITER i, j, k;
  UINT d0, d1, d2;

#if DEBUG
  printf("%s\n", __func__);
#endif
  d0 = mat->d0;
  d1 = mat->d1;
  d2 = mat->d2;
  temp = mem_CMAT(mat->d0, mat->d1, mat->d2);
  if (mat->ndim == 0) {
#pragma omp parallel for shared(temp, mat) private(i)
    for (i = 0; i < d0; i++) {
      temp->data[i].re = mat->data[i].re;
      temp->data[i].im = mat->data[i].im;
    }
    mat->ndim = 1;
    mat->d0 = d1;
    mat->d1 = d0;
  } else if (mat->ndim == 1) {
#pragma omp parallel for shared(temp, mat) private(i)
    for (i = 0; i < d0 * d1; i++) {
      temp->data[i / d0 + i % d0 * d1].re = mat->data[i].re;
      temp->data[i / d0 + i % d0 * d1].im = mat->data[i].im;
    }
    mat->d0 = d1;
    mat->d1 = d0;
  } else {
//일단은 3차원 배열을 2차원 배열의 batch로 여기고 구현함
#pragma omp parallel for shared(temp, mat) private(i, j)
    for (j = 0; j < d2; j++) {
      for (i = 0; i < d0 * d1; i++) {
        temp->data[j * d0 * d1 + i / d0 + i % d0 * d1].re =
            mat->data[j * d0 * d1 + i].re;
        temp->data[j * d0 * d1 + i / d0 + i % d0 * d1].im =
            mat->data[j * d0 * d1 + i].im;
      }
    }
    mat->d0 = d1;
    mat->d1 = d0;
  }
  ccopy(temp, mat);
  free_mem_CMAT(temp);
}

/**** hermitial ****/

CMAT *create_hermit(CMAT *mat) {
  ITER i, j;
  UINT d0 = mat->d0;
  UINT d1 = mat->d1;
  UINT d2 = mat->d2;

  CMAT *t_mat;
#if DEBUG
  printf("%s\n", __func__);
#endif
  if (mat->ndim == 0) {
    t_mat = alloc_CMAT(d1, d0);
#pragma omp parallel for shared(t_mat, mat) private(i)
    for (i = 0; i < d0; i++) {
      t_mat->data[i].re = mat->data[i].re;
      t_mat->data[i].im = -(mat->data[i].im);
    }
  } else if (mat->ndim == 1) {
    t_mat = alloc_CMAT(d1, d0);
#pragma omp parallel for shared(t_mat, mat) private(i)
    for (i = 0; i < d0 * d1; i++) {
      t_mat->data[i / d0 + i % d0 * d1].re = mat->data[i].re;
      t_mat->data[i / d0 + i % d0 * d1].im = -(mat->data[i].im);
    }
  } else {
    t_mat = alloc_CMAT(d1, d0, d2);
//일단은 3차원 배열을 2차원 배열의 batch로 여기고 구현함
#pragma omp parallel for shared(t_mat, mat) private(i, j)
    for (j = 0; j < d2; j++) {
      for (i = 0; i < d0 * d1; i++) {
        t_mat->data[j * d0 * d1 + i / d0 + i % d0 * d1].re =
            mat->data[j * d0 * d1 + i].re;
        t_mat->data[j * d0 * d1 + i / d0 + i % d0 * d1].im =
            -(mat->data[j * d0 * d1 + i].im);
      }
    }
  }
  return t_mat;
}

void hermit(CMAT *mat) {
  CMAT *temp;
  ITER i, j, k;
  UINT d0, d1, d2;

#if DEBUG
  printf("%s\n", __func__);
#endif
  d0 = mat->d0;
  d1 = mat->d1;
  d2 = mat->d2;
  temp = mem_CMAT(mat->d0, mat->d1, mat->d2);
  if (mat->ndim == 0) {
#pragma omp parallel for shared(temp, mat) private(i)
    for (i = 0; i < d0; i++) {
      temp->data[i].re = mat->data[i].re;
      temp->data[i].im = -mat->data[i].im;
    }
    mat->ndim = 1;
    mat->d0 = d1;
    mat->d1 = d0;
  } else if (mat->ndim == 1) {
#pragma omp parallel for shared(temp, mat) private(i)
    for (i = 0; i < d0 * d1; i++) {
      temp->data[i / d0 + i % d0 * d1].re = mat->data[i].re;
      temp->data[i / d0 + i % d0 * d1].im = -mat->data[i].im;
    }
    mat->d0 = d1;
    mat->d1 = d0;
  } else {
//일단은 3차원 배열을 2차원 배열의 batch로 여기고 구현함
#pragma omp parallel for shared(temp, mat) private(i, j)
    for (j = 0; j < d2; j++) {
      for (i = 0; i < d0 * d1; i++) {
        temp->data[j * d0 * d1 + i / d0 + i % d0 * d1].re =
            mat->data[j * d0 * d1 + i].re;
        temp->data[j * d0 * d1 + i / d0 + i % d0 * d1].im =
            -mat->data[j * d0 * d1 + i].im;
      }
    }
    mat->d0 = d1;
    mat->d1 = d0;
  }
  ccopy(temp, mat);
  free_mem_CMAT(temp);
}

/**** Identity Matrix****/
void id_MAT(MAT *mat) {
  ITER i, j, k;
  UINT d0, d2;

  if (mat->d0 != mat->d1) {
    printf("ERROR : This matrix is not square matrix.\n");
    return;
  }

  if (mat->ndim == 1) {
    d0 = mat->d0;

#pragma omp parallel for shared(mat) private(i, j)
    for (i = 0; i < d0; i++) {
      for (j = 0; j < d0; j++) {
        if (i == j) {
          mat->data[i * d0 + j] = 1;
        } else
          mat->data[i * d0 + j] = 0;
      }
    }
  } else if (mat->ndim == 2) {
    d0 = mat->d0;
    d2 = mat->d2;
#pragma omp parallel for shared(mat) private(i, j, k)
    for (i = 0; i < d2; i++) {
      for (j = 0; j < d0; j++) {
        for (k = 0; k < d0; k++) {
          if (i == j) {
            mat->data[i * d0 * d0 + j * d0 + k] = 1;
          } else
            mat->data[i * d0 * d0 + j * d0 + k] = 0;
        }
      }
    }
  }
}

/**** Inverser Matrix ****/

/**** miscellaneous  ****/
void free_MAT(MAT *mat) {
#if DEBUG
  printf("%s\n", __func__);
#endif
  free(mat->data);
  free(mat);
}

void free_CMAT(CMAT *mat) {
#if DEBUG
  printf("%s\n", __func__);
#endif
  free(mat->data);
  free(mat);
}

void print_MAT(MAT *mat) {
  ITER k, j, i;
#if DEBUG
  printf("%s\n", __func__);
#endif
  for (k = 0; k < mat->d2; k++) {
    for (i = 0; i < mat->d0; i++) {
      for (j = 0; j < mat->d1; j++)
        printf("%+3.3lf ",
               mat->data[k * (mat->d1) * (mat->d0) + j * (mat->d0) + i]);
      printf("\n");
    }
    printf("\n");
  }
}

void print_CMAT(CMAT *mat) {
  ITER k, j, i;
#if DEBUG
  printf("%s\n", __func__);
#endif
  for (k = 0; k < mat->d2; k++) {
    for (i = 0; i < mat->d0; i++) {
      for (j = 0; j < mat->d1; j++) {
        printf("%+3.3lf ",
               mat->data[k * (mat->d1) * (mat->d0) + j * (mat->d0) + i].re);
        printf("%+3.3lf| ",
               mat->data[k * (mat->d1) * (mat->d0) + j * (mat->d0) + i].im);
      }
      printf("\n");
    }
    printf("\n");
  }
}

/**** free_mem_MAT ***/
void free_mem_MAT(MAT *mat) {
  iip_free(mat->data);
  iip_free(mat);
}
void free_mem_CMAT(CMAT *mat) {
  iip_free(mat->data);
  iip_free(mat);
}
