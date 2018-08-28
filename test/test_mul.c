#include "mother.h"

int main() {
  MAT *A, *B, *C;
  ITER i;
  long long total = 0;

  A = alloc_mat(1024, 1024, 4);
  B = alloc_mat(1024, 1024, 4);
  C = alloc_mat(1024, 1024, 4);

  read_mat("../test_data/d_1024_1024_4.bin", A);
  read_mat("../test_data/d_1024_1024_4.bin", B);

  for (i = 0; i < 10; i++) matmul(A, B, C);

  for (i = 0; i < 50; i++) {
    stopwatch(0);
    matmul(A, B, C);
    total += stopwatch(1);
  }

  printf("%lf\n", (double)total / 50000);
}
