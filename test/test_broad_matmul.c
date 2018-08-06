#include "mother.h"

#define colA 4
#define rowA 3
#define rowB 3

int main() {
  MAT* A;
  MAT* B;
  MAT* C;
  ITER i;
  A = zeros(colA, rowA, 4);
  B = zeros(rowA, rowB, 4);
  C = zeros(colA, rowB, 4);

  for (i = 0; i < colA * rowA * 4; i++) A->data[i] = i;

  for (i = 0; i < rowB * rowA * 4; i++) B->data[i] = i;

  print_mat(A);
  print_mat(B);
  matmul(A, B, C);
  print_mat(C);

  free_mat(A);
  free_mat(B);
  free_mat(C);
  return 0;
}
