#include "mother.h"

int main() {
  MAT *A;
  CMAT *CA;

  A = alloc_MAT(20, 5);
  randu(A, -10., 10.);
  print_MAT(A);

  randn(A, 0., 1.);
  print_MAT(A);

  free_MAT(A);

  CA = alloc_CMAT(20, 5);
  crandn(CA, 5., -5., 10., 100.);

  print_CMAT(CA);
  free_CMAT(CA);

  return 0;
}
