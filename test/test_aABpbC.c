#include "mother.h"

#define colA 2
#define rowA 6
#define rowB 10

int main() {
  CTYPE alpha, beta;
  MAT *A, *TA, *B, *TB, *C;
  CMAT *CA, *TCA, *CB, *TCB, *CC;
  int i;

  A = alloc_mat(colA, rowA);
  TA = alloc_mat(rowA, colA);
  B = alloc_mat(rowA, rowB);
  TB = alloc_mat(rowB, rowA);
  C = alloc_mat(colA, rowB);

  CA = alloc_cmat(colA, rowA);
  TCA = alloc_cmat(rowA, colA);
  CB = alloc_cmat(rowA, rowB);
  TCB = alloc_cmat(rowB, rowA);
  CC = alloc_cmat(colA, rowB);

  for (i = 0; i < colA * rowA; i++) A->data[i] = 1;
  for (i = 0; i < colA * rowA; i++) TA->data[i] = 1;
  for (i = 0; i < rowB * rowA; i++) B->data[i] = 1;
  for (i = 0; i < rowB * rowA; i++) TB->data[i] = 1;
  for (i = 0; i < colA * rowB; i++) C->data[i] = 1;

  for (i = 0; i < colA * rowA; i++) {
    CA->data[i].re = 1;
    CA->data[i].im = 2;
  }
  for (i = 0; i < colA * rowA; i++) {
    TCA->data[i].re = 1;
    TCA->data[i].im = 2;
  }
  for (i = 0; i < rowB * rowA; i++) {
    CB->data[i].re = 1;
    CB->data[i].im = 2;
  }
  for (i = 0; i < rowB * rowA; i++) {
    TCB->data[i].re = 1;
    TCB->data[i].im = 2;
  }
  for (i = 0; i < colA * rowB; i++) {
    CC->data[i].re = 1;
    CC->data[i].im = 2;
  }

  /**** REAL ****/

  aABpbC(1, A, B, 0, C);
  print_mat(C);

  aABtpbC(1, A, TB, 0, C);
  print_mat(C);

  aAtBpbC(1, TA, B, 0, C);
  print_mat(C);

  aAtBpbC(1, TA, B, 0, C);
  print_mat(C);

  matmul(A, B, C);
  print_mat(C);

  /**** COMPLEX ****/

  alpha.re = 1;
  beta.re = 0;
  alpha.im = 0;
  beta.im = 0;

  caABpbC(alpha, CA, CB, beta, CC);
  print_cmat(CC);
  caABtpbC(alpha, CA, TCB, beta, CC);
  print_cmat(CC);
  caABhpbC(alpha, CA, TCB, beta, CC);
  print_cmat(CC);

  caAtBpbC(alpha, TCA, CB, beta, CC);
  print_cmat(CC);
  caAtBtpbC(alpha, TCA, TCB, beta, CC);
  print_cmat(CC);
  caAtBhpbC(alpha, TCA, TCB, beta, CC);
  print_cmat(CC);

  caAhBpbC(alpha, TCA, CB, beta, CC);
  print_cmat(CC);
  caAhBtpbC(alpha, TCA, TCB, beta, CC);
  print_cmat(CC);
  caAhBhpbC(alpha, TCA, TCB, beta, CC);
  print_cmat(CC);

  cmatmul(CA, CB, CC);
  print_cmat(CC);

  free_mat(A);
  free_mat(B);
  free_mat(TA);
  free_mat(TB);
  free_mat(C);

  free_cmat(CA);
  free_cmat(CB);
  free_cmat(TCA);
  free_cmat(TCB);
  free_cmat(CC);

  return 0;
}
