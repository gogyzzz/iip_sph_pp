#include "mother.h"

#define N 8

int main()
{
 mkl_handle* hand,ihand;
 MAT*A,*B;
 CMAT*CA,*CB;
 ITER i;

 init(1024);
 hand = fft_handle(N);

 ihand = ifft_handle(N);
 A = zeros(N);
 B = zeros(N);
 CA = czeros(N);
 CB = czeros(N);

 randn(A,0.,100.);
 print_mat(A);
 
 printf("==== MKL ====\n");
 
 mkl_fft(hand,A,CA);
 print_cmat(CA);
 fill(A,0);
 
 mkl_ifft(ihand,CA,A);
 print_mat(A);

 printf("==== Ooura ====\n");

 fft(A,CA);
 print_cmat(CA);

 ifft(CA,A);
 print_mat(A);

 free_handle(hand);
 free_handle(ihand);
 free_mat(A);
 free_cmat(CA);

 finit();
  return 0;
}
