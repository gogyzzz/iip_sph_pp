#include "mother.h"

#define d0 10
#define d1 10
#define d2 10

int main() {
  MAT *A, *B, *C;

  A = alloc_mat(d0, d1, d2);
  B = alloc_mat(d0, d1, d2);
  C = alloc_mat(d0, d1, d2);

  read_mat("../test/d_10_10_10.bin", A);
  stopwatch(0);
  randn(B, 0, 100);
  stopwatch(1);

  stopwatch(0);
  matmul(A, B, C);
  stopwatch(1);

  stopwatch(0);
  invert(C, A);
  stopwatch(1);
  stopwatch(0);
  matmul(A, C, B);
  stopwatch(1);
  print_mat(B);
  free_mat(A);
  free_mat(B);
  free_mat(C);
  return 0;
}
