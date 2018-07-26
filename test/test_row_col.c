#include "mother.h"

int main() {
  MAT* A;
  ITER i;
  A = alloc_MAT(3, 5);

  for (i = 0; i < A->d0 * A->d1; i++) A->data[i] = 1;
  print_MAT(A);

  col_scal(-40, A, 1);
  row_scal(10, A, 2);

  print_MAT(A);

  add(10, A);

  print_MAT(A);

  row_add(-20, A, 1);
  print_MAT(A);

  free_MAT(A);

  return 0;
}
