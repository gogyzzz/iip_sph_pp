#include "header_for_test.h"
#include "mother.h"

int main() {
  MAT *A, *X, *Y;
  MAT *TA, *TX, *TY;

  stopwatch(0);

  A = zeros(5000, 4000);
  X = alloc_mat(5000);
  Y = alloc_mat(4000);

  fill(A, 2);
  fill(X, 1);
  fill(Y, -1);

  gemv(NoTran, 10, A, X, 3, Y);

  free_mat(A);
  free_mat(X);
  free_mat(Y);

  /**** TRANSPOSE  ****/

  TA = alloc_mat(4000, 5000);
  TX = alloc_mat(5000);
  TY = alloc_mat(4000);

  fill(TA, 2);
  fill(TX, 1);
  fill(TY, -1);

  gemv(Tran, 10, TA, TX, 3, TY);

  free_mat(TA);
  free_mat(TX);
  free_mat(TY);

  stopwatch(1);
  return 0;
}
