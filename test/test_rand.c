#include "mother.h"

int main() {
  MAT *A;
  CMAT *CA;

  A = alloc_mat(20, 5);
  randu(A, -10., 10.);
  print_mat(A);

  randn(A, 0., 1.);
  print_mat(A);

  free_mat(A);

  CA = alloc_cmat(20, 5);
  crandn(CA, 5., -5., 10., 100.);

  print_cmat(CA);
  free_cmat(CA);

  return 0;
}
