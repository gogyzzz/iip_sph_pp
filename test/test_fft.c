#include "mother.h"

#define a0 0 //FFT
#define a1 0 //Complex
#define a2 1 //Half

#define N 8
int main()
{
 MAT*A;
 CMAT*B;
 CTYPE t;
 char file[MAX_CHAR];
 long long t1=0,t2=0;
 ITER i;
//Init
 t.re = 0;
 t.im = 0;
 A = zeros(N);
 sprintf(file,"../test_data/d_%d_1_1.bin",N);
 read_mat(file,A);
#if a0

// print_mat(A);
 B = czeros(N);
 for(i=0;i<50;i++)fft(A,B);
 for(i=0;i<500;i++){
   stopwatch(0);
   fft(A,B);
   t1 +=stopwatch(1);
 }

// print_cmat(B);
 fill(A,0);

 for(i=0;i<50;i++){
   stopwatch(0);
   ifft(B,A);
   t2 +=stopwatch(1);
 }

// print_mat(A);

 printf("fft  : %lf\n",((double)t1)/500000);
 printf("ifft : %lf\n",((double)t2)/50000);
#endif

#if a1
 B = czeros(N);
fft(A,B);
print_cmat(B);
cfft(B,A);
print_mat(A);
cfill(B, t);
cifft(A,B);
print_cmat(B);
#endif

#if a2
print_mat(A);
 B = czeros(N/2 + 1);
hfft(A,B);
print_cmat(B);
fill(A,0);
hifft(B,A);
print_mat(A);
#endif

free_mat(A);
 free_cmat(B);
 return 0;
}
