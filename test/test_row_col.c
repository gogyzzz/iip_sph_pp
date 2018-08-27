#include "mother.h"

int main() {
  MAT* A;
  ITER i;
  A = alloc_mat(3, 5);

  for (i = 0; i < A->d0 * A->d1; i++) A->data[i] = 1;
  print_mat(A);

  col_scal(-40, A, 1);
  row_scal(10, A, 2);

  print_mat(A);

  add_mat(10, A);

  print_mat(A);

  row_add(-20, A, 1);
  print_mat(A);

  free_mat(A);

  return 0;
}
