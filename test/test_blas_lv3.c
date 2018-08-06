#include "header_for_test.h"
#include "mother.h"

#define _print 1

#define _A 0  // N N
#define _B 0  // T N
#define _C 0  // T T
#define _D 0  // N T

#define _E 0  // N N
#define _F 0  // N T
#define _G 0  // N C
#define _H 0  // T N
#define _I 0  // T T
#define _J 0  // T C
#define _K 0  // C M
#define _L 0  // C T
#define _M 1  // C C

#define colA 10
#define rowA 10
#define rowB 10
int main() {
  MAT *A, *B, *C;
  CMAT *CA, *CB, *CC;
  CTYPE alpha, beta;
  ITER i;

  stopwatch(0);

#if USE_CUDA
  init();
#endif

  alpha.re = 1;
  alpha.im = 1;
  beta.re = 0;
  beta.im = 0;
/**** REAL MATRIX ****/

#if _A
  A = alloc_mat(colA, rowA);
  B = alloc_mat(rowA, rowB);
  C = alloc_mat(colA, rowB);

  fill(A, 1.0);
  fill(B, 3.0);
  fill(C, -0.1);
  for (i = 0; i < A->d0 * A->d1; i++) {
    A->data[i] = i;
  }
  for (i = 0; i < B->d0 * B->d1; i++) {
    B->data[i] = i;
  }

#if _print
  print_mat(A);
  print_mat(B);
  print_mat(C);
#endif

  gemm(NoTran, NoTran, 1, A, B, 0, C);

#if _print
  print_mat(C);
#endif

  free_mat(A);
  free_mat(B);
  free_mat(C);
#endif

#if _B
  A = alloc_mat(rowA, colA);
  B = alloc_mat(rowA, rowB);
  C = alloc_mat(colA, rowB);

  fill(A, 1.0);
  fill(B, 3.0);
  fill(C, -0.1);
  for (i = 0; i < A->d0 * A->d1; i++) {
    A->data[i] = i;
  }
  for (i = 0; i < B->d0 * B->d1; i++) {
    B->data[i] = i;
  }

#if _print
  print_mat(A);
  print_mat(B);
  print_mat(C);
#endif

  gemm(Tran, NoTran, 1, A, B, 0, C);

#if _print
  print_mat(C);
#endif
  free_mat(A);
  free_mat(B);
  free_mat(C);
#endif

#if _C
  A = alloc_mat(rowA, colA);
  B = alloc_mat(rowB, rowA);
  C = alloc_mat(colA, rowB);

  fill(A, 1.0);
  fill(B, 3.0);
  fill(C, -0.1);
  for (i = 0; i < A->d0 * A->d1; i++) {
    A->data[i] = i;
  }
  for (i = 0; i < B->d0 * B->d1; i++) {
    B->data[i] = i;
  }

#if _print
  print_mat(A);
  print_mat(B);
  print_mat(C);
#endif

  gemm(Tran, Tran, 1, A, B, 0, C);

#if _print
  print_mat(C);
#endif
  free_mat(A);
  free_mat(B);
  free_mat(C);
#endif

#if _D
  A = alloc_mat(colA, rowA);
  B = alloc_mat(rowB, rowA);
  C = alloc_mat(colA, rowB);

  fill(A, 1.0);
  fill(B, 3.0);
  fill(C, -0.1);
  for (i = 0; i < A->d0 * A->d1; i++) {
    A->data[i] = i;
  }
  for (i = 0; i < B->d0 * B->d1; i++) {
    B->data[i] = i;
  }

#if _print
  print_mat(A);
  print_mat(B);
  print_mat(C);
#endif

  gemm(NoTran, Tran, 1, A, B, 0, C);

#if _print
  print_mat(C);
#endif
  free_mat(A);
  free_mat(B);
  free_mat(C);
#endif

/**** COMPLEX MATRIX ****/

#if _E
  CA = alloc_cmat(colA, rowA);
  CB = alloc_cmat(rowA, rowB);
  CC = alloc_cmat(colA, rowB);

  cfill(CA, 1.0, 1);
  cfill(CB, 3.0, 3);
  cfill(CC, -0.1, -0.1);

  for (i = 0; i < CA->d0 * CA->d1; i++) {
    CA->data[i].re = i;
    CA->data[i].im = i;
  }
  for (i = 0; i < CB->d0 * CB->d1; i++) {
    CB->data[i].re = i;
    CB->data[i].im = i;
  }
#if _print
  print_cmat(CA);
  print_cmat(CB);
  print_cmat(CC);
#endif

  cgemm(NoTran, NoTran, alpha, CA, CB, beta, CC);

#if _print
  print_cmat(CC);
#endif
  free_cmat(CA);
  free_cmat(CB);
  free_cmat(CC);
#endif
#if _F
  CA = alloc_cmat(colA, rowA);
  CB = alloc_cmat(rowB, rowA);
  CC = alloc_cmat(colA, rowB);

  cfill(CA, 1.0, 1);
  cfill(CB, 3.0, 3);
  cfill(CC, -0.1, -0.1);
  for (i = 0; i < CA->d0 * CA->d1; i++) {
    CA->data[i].re = i;
    CA->data[i].im = i;
  }
  for (i = 0; i < CB->d0 * CB->d1; i++) {
    CB->data[i].re = i;
    CB->data[i].im = i;
  }

#if _print
  print_cmat(CA);
  print_cmat(CB);
  print_cmat(CC);
#endif

  cgemm(NoTran, Tran, alpha, CA, CB, beta, CC);

#if _print
  print_cmat(CC);
#endif
  free_cmat(CA);
  free_cmat(CB);
  free_cmat(CC);
#endif

#if _G
  CA = alloc_cmat(colA, rowA);
  CB = alloc_cmat(rowB, rowA);
  CC = alloc_cmat(colA, rowB);

  cfill(CA, 1.0, 1);
  cfill(CB, 3.0, -3);
  cfill(CC, -0.1, -0.1);
  for (i = 0; i < CA->d0 * CA->d1; i++) {
    CA->data[i].re = i;
    CA->data[i].im = i;
  }
  for (i = 0; i < CB->d0 * CB->d1; i++) {
    CB->data[i].re = i;
    CB->data[i].im = i;
  }

#if _print
  print_cmat(CA);
  print_cmat(CB);
  print_cmat(CC);
#endif

  cgemm(NoTran, CTran, alpha, CA, CB, beta, CC);

#if _print
  print_cmat(CC);
#endif
  free_cmat(CA);
  free_cmat(CB);
  free_cmat(CC);
#endif
#if _H
  CA = alloc_cmat(rowA, colA);
  CB = alloc_cmat(rowA, rowB);
  CC = alloc_cmat(colA, rowB);

  cfill(CA, 1.0, 1);
  cfill(CB, 3.0, 3);
  cfill(CC, -0.1, -0.1);
  for (i = 0; i < CA->d0 * CA->d1; i++) {
    CA->data[i].re = i;
    CA->data[i].im = i;
  }
  for (i = 0; i < CB->d0 * CB->d1; i++) {
    CB->data[i].re = i;
    CB->data[i].im = i;
  }
#if _print
  print_cmat(CA);
  print_cmat(CB);
  print_cmat(CC);
#endif

  cgemm(Tran, NoTran, alpha, CA, CB, beta, CC);

#if _print
  print_cmat(CC);
#endif
  free_cmat(CA);
  free_cmat(CB);
  free_cmat(CC);
#endif
#if _I
  CA = alloc_cmat(rowA, colA);
  CB = alloc_cmat(rowB, rowA);
  CC = alloc_cmat(colA, rowB);

  cfill(CA, 1.0, 1);
  cfill(CB, 3.0, 3);
  cfill(CC, -0.1, -0.1);
  for (i = 0; i < CA->d0 * CA->d1; i++) {
    CA->data[i].re = i;
    CA->data[i].im = i;
  }
  for (i = 0; i < CB->d0 * CB->d1; i++) {
    CB->data[i].re = i;
    CB->data[i].im = i;
  }
#if _print
  print_cmat(CA);
  print_cmat(CB);
  print_cmat(CC);
#endif

  cgemm(Tran, Tran, alpha, CA, CB, beta, CC);

#if _print
  print_cmat(CC);
#endif
  free_cmat(CA);
  free_cmat(CB);
  free_cmat(CC);
#endif
#if _J
  CA = alloc_cmat(rowA, colA);
  CB = alloc_cmat(rowB, rowA);
  CC = alloc_cmat(colA, rowB);

  cfill(CA, 1.0, 1);
  cfill(CB, 3.0, -3);
  cfill(CC, -0.1, -0.1);
  for (i = 0; i < CA->d0 * CA->d1; i++) {
    CA->data[i].re = i;
    CA->data[i].im = i;
  }
  for (i = 0; i < CB->d0 * CB->d1; i++) {
    CB->data[i].re = i;
    CB->data[i].im = i;
  }

#if _print
  print_cmat(CA);
  print_cmat(CB);
  print_cmat(CC);
#endif

  cgemm(Tran, CTran, alpha, CA, CB, beta, CC);

#if _print
  print_cmat(CC);
#endif
  free_cmat(CA);
  free_cmat(CB);
  free_cmat(CC);
#endif
#if _K
  CA = alloc_cmat(rowA, colA);
  CB = alloc_cmat(rowA, rowB);
  CC = alloc_cmat(colA, rowB);

  cfill(CA, 1.0, -1);
  cfill(CB, 1.0, 1);
  cfill(CC, -0.1, -0.1);
  for (i = 0; i < CA->d0 * CA->d1; i++) {
    CA->data[i].re = i;
    CA->data[i].im = i;
  }
  for (i = 0; i < CB->d0 * CB->d1; i++) {
    CB->data[i].re = i;
    CB->data[i].im = i;
  }
#if _print
  print_cmat(CA);
  print_cmat(CB);
  print_cmat(CC);
#endif

  cgemm(CTran, NoTran, alpha, CA, CB, beta, CC);

#if _print
  print_cmat(CC);
#endif
  free_cmat(CA);
  free_cmat(CB);
  free_cmat(CC);
#endif
#if _L
  CA = alloc_cmat(rowA, colA);
  CB = alloc_cmat(rowB, rowA);
  CC = alloc_cmat(colA, rowB);

  cfill(CA, 1.0, -1);
  cfill(CB, 1.0, 1);
  cfill(CC, -0.1, -0.1);
  for (i = 0; i < CA->d0 * CA->d1; i++) {
    CA->data[i].re = i;
    CA->data[i].im = i;
  }
  for (i = 0; i < CB->d0 * CB->d1; i++) {
    CB->data[i].re = i;
    CB->data[i].im = i;
  }

#if _print
  print_cmat(CA);
  print_cmat(CB);
  print_cmat(CC);
#endif

  cgemm(CTran, Tran, alpha, CA, CB, beta, CC);

#if _print
  print_cmat(CC);
#endif
  free_cmat(CA);
  free_cmat(CB);
  free_cmat(CC);
#endif
#if _M
  CA = alloc_cmat(rowA, colA);
  CB = alloc_cmat(rowB, rowA);
  CC = alloc_cmat(colA, rowB);

  cfill(CA, 1.0, -1);
  cfill(CB, 1.0, -1);
  cfill(CC, -0.1, -0.1);

  for (i = 0; i < CA->d0 * CA->d1; i++) {
    CA->data[i].re = i;
    CA->data[i].im = i;
  }
  for (i = 0; i < CB->d0 * CB->d1; i++) {
    CB->data[i].re = i;
    CB->data[i].im = i;
  }

#if _print
  print_cmat(CA);
  print_cmat(CB);
  print_cmat(CC);
#endif

  cgemm(CTran, CTran, alpha, CA, CB, beta, CC);

#if _print
  print_cmat(CC);
#endif
  free_cmat(CA);
  free_cmat(CB);
  free_cmat(CC);
#endif
#if USE_CUDA
  finit();
#endif
  stopwatch(1);

  return 0;
}
