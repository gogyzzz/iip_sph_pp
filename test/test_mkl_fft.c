#include "mother.h"

#define N 65536

int main(){
 mkl_handle* handle;
 MAT*in;
 CMAT*out;
 char file[MAX_CHAR];
 double*a; 
 int* ip;
 double* w;
 ITER i,j;
 long long t=0;

  init(1024);

 in = zeros(N);
 out = czeros((N/2)+1);

 sprintf(file,"../test_data/d_%d_1_1.bin",N); 
 read_mat(file,in);
 stopwatch(0);
 handle = create_handle(N);
 t = stopwatch(1);
 printf("%lld\n",t);
 t=0;


 for(j=0;j<50;j++){
   mkl_fft(handle,in->data,out->data);
 }

 for(j=0;j<500;j++){
   stopwatch(0);
   mkl_fft(handle,in->data,out->data);
   t+=stopwatch(1);
 }
 printf("%lld\n",t/500);

  a = mpalloc(sizeof(double)*N);
  ip = mpalloc(sizeof(int)*((int)(sqrt(N/2))+1));
  w = mpalloc(sizeof(double)*(N/2)); 
  ip[0]=0;
  
  t=0;
  for(j=0;j<500;j++){
   for(i=0;i<N;i++){
     a[i] = in->data[i];
    }
   stopwatch(0);
   rdft(N,1,a,ip,w);
   for(i=0;i<(N/2);i++){
     out->data[i].re = a[2*i];
     out->data[i].im =- a[2*i+1];
   }
   t +=stopwatch(1);
  }
 printf("%lld\n",t/500);

  free_handle(handle);
  finit();
  free_mat(in);
  free_cmat(out);

  return 0;
}
