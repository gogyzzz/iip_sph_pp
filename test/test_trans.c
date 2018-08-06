#include "mother.h"

int main() {
  MAT *A;
  CMAT *CA;
  ITER i;

  A = alloc_mat(4);
  for (i = 0; i < A->d0 * A->d1 * A->d2; i++) A->data[i] = i;
  print_mat(A);
  trans(A);

  print_mat(A);
  free_mat(A);

  A = alloc_mat(4, 3);
  for (i = 0; i < A->d0 * A->d1 * A->d2; i++) A->data[i] = i;
  print_mat(A);
  trans(A);
  print_mat(A);

  free_mat(A);

  CA = alloc_cmat(5, 3);
  for (i = 0; i < 15; i++) {
    CA->data[i].re = 1;
    CA->data[i].im = 5;
  }
  print_cmat(CA);

  ctrans(CA);
  print_cmat(CA);

  hermit(CA);
  print_cmat(CA);

  free_cmat(CA);

  return 0;
}
