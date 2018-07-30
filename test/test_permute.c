#include "mother.h"
#define seq 321

int main() {
  MAT* A;
  ITER i;

  printf("permute seq : %d\n", seq);

  init();
  A = zeros(4, 3, 2);
  for (i = 0; i < 24; i++) A->data[i] = i;
  print_MAT(A);
  permute(A, seq);
  print_MAT(A);
  free_MAT(A);

  A = zeros(1, 5, 2);
  for (i = 0; i < 10; i++) A->data[i] = i;
  print_MAT(A);
  permute(A, seq);
  print_MAT(A);
  free_MAT(A);

  A = zeros(5, 1, 2);
  for (i = 0; i < 10; i++) A->data[i] = i;
  print_MAT(A);
  permute(A, seq);
  print_MAT(A);
  free_MAT(A);

  finit();
  return 0;
}
