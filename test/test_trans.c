#include "mother.h"

int main() {
  MAT *A;
  CMAT *CA;
  ITER i;

  A = alloc_MAT(4);
  for (i = 0; i < A->d0 * A->d1 * A->d2; i++) A->data[i] = i;
  print_MAT(A);
  trans(A);

  print_MAT(A);
  free_MAT(A);

  A = alloc_MAT(4, 3);
  for (i = 0; i < A->d0 * A->d1 * A->d2; i++) A->data[i] = i;
  print_MAT(A);
  trans(A);
  print_MAT(A);

  free_MAT(A);

  CA = alloc_CMAT(5, 3);
  for (i = 0; i < 15; i++) {
    CA->data[i].re = 1;
    CA->data[i].im = 5;
  }
  print_CMAT(CA);

  ctrans(CA);
  print_CMAT(CA);

  hermit(CA);
  print_CMAT(CA);

  free_CMAT(CA);

  return 0;
}
