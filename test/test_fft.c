#define test0 0  // FFT
#define test1 0  // Complex
#define test2 0  // Half
#define test3 0  // fftw3
#define test4 0  // ooura vs fftw3
#define test5 1  // ooura vs fftw

#if !test5
#include <fftw3.h>
#endif

#if test5
#include <fftw.h>
#endif

#include "mother.h"

#define N 4
int main() {
  MAT *A;
  CMAT *A_, *B, *B_;
  CTYPE t;
  char file[MAX_CHAR];
  long long t1 = 0, t2 = 0;
  ITER i;

#if test3
  fftw_plan p;
#endif

#if test4
  fftw_plan p;
#endif

#if test5
  fftw_plan p;
#endif

  // Init
  t.re = 0;
  t.im = 0;
  A = zeros(N);
  //sprintf(file, "../test_data/d_%d_1_1.bin", N);
  //read_mat(file, A);
  randn(A,0,100);
  copy_mat_inc(N, A->data, 1, A_->data, 2);
#if test0

  print_mat(A);
  B = czeros(N);
  fft(A, B);

  print_cmat(B);
  fill(A, 0);

  ifft(B, A);

  print_mat(A);

#endif

#if test1
  B = czeros(N);
  fft(A, B);
  print_cmat(B);
  cfft(B, A);
  print_mat(A);
  cfill(B, t);
  cifft(A, B);
  print_cmat(B);
#endif

#if test2
  print_mat(A);
  B = czeros(N / 2 + 1, N);
  hfft(A, B);
  print_cmat(B);
  fill(A, 0);
  hifft(B, A);
  print_mat(A);
#endif

#if test3
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

#if test4
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

#if test5
  printf("==== %d ====\n", N);
  B = czeros(N / 2 + 1);
  B_ = czeros(N / 2 + 1);

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
  p = fftw_create_plan(N, FFTW_FORWARD, FFTW_MEASURE | FFTW_OUT_OF_PLACE);
  t2 = stopwatch(1);

  printf("PLAN  : %lf\n", (double)t2);

  t1 = 0;
  for (i = 0; i < 500; i++) {
	  stopwatch(0);
	  fftw_one(p, A_, B_);
	  t1 += stopwatch(1);
  }
  printf("EXEC  : %lf\n", ((double)t1) / (500));
  printf("FFTW  : %lf\n", (double)t2 + (double)t1 / (500));
  // print_cmat(B);

  fftw_destroy_plan(p);
#endif
  free_mat(A);
  free_cmat(A_);
  free_cmat(B);
  free_cmat(B_);
  return 0;
}
