#include "mother.h"
#define seq 321

int main() {
  MAT* A;
  ITER i;

  printf("permute seq : %d\n", seq);

  init(0);
  A = zeros(4, 3, 2);
  for (i = 0; i < 24; i++) A->data[i] = i;
  print_mat(A);
  permute(A, seq);
  print_mat(A);
  free_mat(A);

  A = zeros(1, 5, 2);
  for (i = 0; i < 10; i++) A->data[i] = i;
  print_mat(A);
  permute(A, seq);
  print_mat(A);
  free_mat(A);

  A = zeros(5, 1, 2);
  for (i = 0; i < 10; i++) A->data[i] = i;
  print_mat(A);
  permute(A, seq);
  print_mat(A);
  free_mat(A);

  A = zeros(50, 50, 50);
  free_mat(A);

  finit();
  return 0;
}
