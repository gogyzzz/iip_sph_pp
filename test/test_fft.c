#include <fftw3.h>
#include "mother.h"

#define a0 0  // FFT
#define a1 0  // Complex
#define a2 0  // Half
#define a3 1  // fftw3
#define a4 0  // ooura vs fftw3

#define N 4
int main() {
  MAT *A;
  CMAT *B;
  CTYPE t;
  char file[MAX_CHAR];
  long long t1 = 0, t2 = 0;
  ITER i;

#if a3
  fftw_plan p;
#endif

#if a4
  fftw_plan p;
#endif
  // Init
  t.re = 0;
  t.im = 0;
  A = zeros(N);
  sprintf(file, "../test_data/d_%d_1_1.bin", N);
  read_mat(file, A);
#if a0

  print_mat(A);
  B = czeros(N);
  fft(A, B);

  print_cmat(B);
  fill(A, 0);

  ifft(B, A);

  print_mat(A);

#endif

#if a1
  B = czeros(N);
  fft(A, B);
  print_cmat(B);
  cfft(B, A);
  print_mat(A);
  cfill(B, t);
  cifft(A, B);
  print_cmat(B);
#endif

#if a2
  print_mat(A);
  B = czeros(N / 2 + 1, N);
  hfft(A, B);
  print_cmat(B);
  fill(A, 0);
  hifft(B, A);
  print_mat(A);
#endif

#if a3
  fftw_complex *in, *out;
  out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);
  print_mat(A);
  B = czeros(N);
  p = fftw_plan_dft_r2c_1d(N, A->data, B->data, 0);

  printf("%d\n", sizeof(fftw_complex));

  fftw_execute(p); /* repeat as needed */
  print_cmat(B);
  fftw_execute(p); /* repeat as needed */
  print_cmat(B);
  fftw_execute(p); /* repeat as needed */

  print_cmat(B);

  fftw_destroy_plan(p);
#endif

#if a4
  printf("==== %d ====\n", N);
  B = czeros(N / 2 + 1);

  for (i = 0; i < 50; i++) hfft(A, B);

  t1 = 0;
  for (i = 0; i < 500; i++) {
    stopwatch(0);
    hfft(A, B);
    t1 += stopwatch(1);
  }

  printf("OOURA : %lf\n", ((double)t1) / (500));
  // print_cmat(B);

  t2 = 0;
  stopwatch(0);
  p = fftw_plan_dft_r2c_1d(N, A->data, B->data, 0);
  t2 = stopwatch(1);

  printf("PLAN  : %lf\n", (double)t2);

  t1 = 0;
  for (i = 0; i < 500; i++) {
    stopwatch(0);
    fftw_execute(p);
    t1 += stopwatch(1);
  }
  printf("EXEC  : %lf\n", ((double)t1) / (500));
  printf("FFTW  : %lf\n", (double)t2 + (double)t1 / (500));
  // print_cmat(B);

  fftw_destroy_plan(p);

#endif

  free_mat(A);
  free_cmat(B);
  return 0;
}
