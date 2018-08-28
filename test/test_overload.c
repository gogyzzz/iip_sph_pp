#include "mother.h"

int main() {
  MAT *A, *B, *C;

  A = zeros(4);
  B = zeros(4, 3);
  C = alloc_mat(4, 3, 2);
  fill(C, 5);

  free_mat(A);
  free_mat(B);
  free_mat(C);

  return;
}
