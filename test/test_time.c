#include "mother.h"

#define d0 7
#define d2 4

int main() {
  MAT *A, *B, *C;
  ITER i;
  long long total = 0;
  char file[MAX_CHAR];
  A = zeros(d0, d0, d2);
  B = zeros(d0, d0, d2);

  sprintf(file, "../test_data/d_%d_%d_4.bin", d0, d0);
  printf("reading : %s\n", file);
  read_mat(file, A);
  for (i = 0; i < 5; i++) {
    // add(1,A);
    invert(A, B);
  }
  for (i = 0; i < 50; i++) {
    stopwatch(0);
    // add(1,A);
    invert_nbyn(A->data, B->data, d0);
    total += stopwatch(1);
  }
  print_mat(B);
  printf("%lf\n", (double)total / (50 * 1000));
  free_mat(A);
  free_mat(B);

  return 0;
}
