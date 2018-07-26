#include "header_for_test.h"
#include "mother.h"

int main() {
  MAT *A, *X, *Y;
  MAT *TA, *TX, *TY;

  stopwatch(0);

  A = zeros(5000, 4000);
  X = alloc_MAT(5000);
  Y = alloc_MAT(4000);

  fill(A, 2);
  fill(X, 1);
  fill(Y, -1);

  gemv(NoTran, 10, A, X, 3, Y);

  free_MAT(A);
  free_MAT(X);
  free_MAT(Y);

  /**** TRANSPOSE  ****/

  TA = alloc_MAT(4000, 5000);
  TX = alloc_MAT(5000);
  TY = alloc_MAT(4000);

  fill(TA, 2);
  fill(TX, 1);
  fill(TY, -1);

  gemv(Tran, 10, TA, TX, 3, TY);

  free_MAT(TA);
  free_MAT(TX);
  free_MAT(TY);

  stopwatch(1);
  return 0;
}
