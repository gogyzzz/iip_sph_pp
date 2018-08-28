#include "mother.h"
#include "fftw3.h"

#define loop 100

int main(){
 mkl_handle* handle;
 MAT*in;
 CMAT*out;
 char file[MAX_CHAR];
 double*a; 
 int* ip;
 double* w;
 ITER i,j,k;
 fftw_plan p;
 long long t=0;
 UINT N;
 CTYPE cz={0.,0.};

 init(1024);

 printf("N | MKL_HANDLE | MKL_FFT | FFTW3 PLAN | FFTW3_FFT | OOURA\n");
 printf("--- | --- | ---|--- | --- | ---\n");

 for(N = 65536;N>=32;N=N/2){

 in = zeros(N);
 out = czeros((N/2)+1);

 sprintf(file,"../test_data/d_%d_1_1.bin",N); 
 read_mat(file,in);
 printf("%6d | ",N);
 stopwatch(0);
 handle = fft_handle(in->d0);
 t = stopwatch(1);
 printf("%6lld |",t);
 t=0;

 for(j=0;j<loop;j++){
   mkl_hfft(handle,in,out);
 }
 cfill(out,cz);
 for(j=0;j<loop*5;j++){
   stopwatch(0);
   mkl_fft(handle,in,out);
   t+=stopwatch(1);
 }
 printf("%6lf |",(double)t/(loop*5));
//  if(N==32){printf("== mkl ==\n");print_cmat(out);}

  stopwatch(0);
	p = fftw_plan_dft_r2c_1d(N, in->data, out->data,0);
	t = stopwatch(1);
	printf("%6lld |", t);
	t = 0;
	
 cfill(out,cz);
	for (i = 0; i < loop*5; i++) {
		stopwatch(0);
		fftw_execute(p);
		t += stopwatch(1);
	}
	printf("%6lf |",(double)t/(loop*5));
  //if(N==32){printf("== fftw ==\n");print_cmat(out);}


  a = mpalloc(sizeof(double)*N);
  ip = mpalloc(sizeof(int)*((int)(sqrt(N/2))+1));
  w = mpalloc(sizeof(double)*(N/2)); 
  ip[0]=0;
  
  t=0;
 
 cfill(out,cz);
  for(j=0;j<loop*5;j++){
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
  printf("%6lf \n",(double)t/(loop*5));
 // if(N==32){printf("== ooura ==\n");print_cmat(out);}
 
  mpfree(w);
  mpfree(ip);
  mpfree(a);
  free_mat(in);
  free_cmat(out);
  fftw_destroy_plan(p);
  free_handle(handle);
 }
  finit();

  return 0;
}
