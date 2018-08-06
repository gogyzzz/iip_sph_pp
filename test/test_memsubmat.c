#include "mother.h"

#define m 3
#define n 4
#define k 5

int main() {
  MAT* A;
  MAT* subA;

  A = mpalloc_mat(m, n);

  fill(A, 1);

  row_add(-10, A, 1);
  col_scal(3, A, 1);

  print_mat(A);

  subA = mpsubmat(A, -1, -1, 0, 3);

  print_mat(subA);

  free_mpalloc_mat(subA);
  free_mpalloc_mat(A);
  return 0;
}
