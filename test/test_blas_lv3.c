#include "mother.h"

#define _print 1
	
#define _A 0  //N N
#define _B 0  //T N
#define _C 0  //T T
#define _D 0  //N T

#define _E 1  //N N
#define _F 1  //T N
#define _G 1  //T T
#define _H 1  //N T


int main()
{
	MAT*A,*B,*C;
	CMAT* CA, *CB, *CC;
  CTYPE alpha,beta;

#if USE_CUDA
init();
#endif

	alpha.re=1;
	alpha.im=1;
	beta.re=0;
	beta.im=0;
/**** REAL MATRIX ****/

#if _A
	A = alloc_MAT(3,4);
	B = alloc_MAT(4,5);
	C = alloc_MAT(3,5);

  fill(A,1.0); 
  fill(B,3.0); 
  fill(C,-0.1); 

	#if _print
	print_MAT(A);	
	print_MAT(B);	
	print_MAT(C);	
	#endif

	gemm(NoTran,NoTran,1,A,B,0,C);

	#if _print
	print_MAT(C);	
	#endif
  free_MAT(A);
  free_MAT(B);
  free_MAT(C);
#endif

#if _B
	A = alloc_MAT(4,3);
	B = alloc_MAT(4,5);
	C = alloc_MAT(3,5);

  fill(A,1.0); 
  fill(B,3.0); 
  fill(C,-0.1); 

	#if _print
	print_MAT(A);	
	print_MAT(B);	
	print_MAT(C);	
	#endif

	gemm(Tran,NoTran,1,A,B,0,C);

	#if _print
	print_MAT(C);	
	#endif
  free_MAT(A);
  free_MAT(B);
  free_MAT(C);
#endif

#if _C
	A = alloc_MAT(4,3);
	B = alloc_MAT(5,4);
	C = alloc_MAT(3,5);

  fill(A,1.0); 
  fill(B,3.0); 
  fill(C,-0.1); 

	#if _print
	print_MAT(A);	
	print_MAT(B);	
	print_MAT(C);	
	#endif

	gemm(Tran,Tran,1,A,B,0,C);

	#if _print
	print_MAT(C);	
	#endif
  free_MAT(A);
  free_MAT(B);
  free_MAT(C);
#endif

#if _D
	A = alloc_MAT(3,4);
	B = alloc_MAT(5,4);
	C = alloc_MAT(3,5);

  fill(A,1.0); 
  fill(B,3.0); 
  fill(C,-0.1); 

	#if _print
	print_MAT(A);	
	print_MAT(B);	
	print_MAT(C);	
	#endif

	gemm(NoTran,Tran,1,A,B,0,C);

	#if _print
	print_MAT(C);	
	#endif
  free_MAT(A);
  free_MAT(B);
  free_MAT(C);
#endif

/**** COMPLEX MATRIX ****/

#if _E
	CA = alloc_CMAT(3,4);
	CB = alloc_CMAT(4,5);
	CC = alloc_CMAT(3,5);

  cfill(CA,1.0,0); 
  cfill(CB,3.0,0); 
  cfill(CC,-0.1,-0.1); 

	#if _print
	print_CMAT(CA);	
	print_CMAT(CB);	
	print_CMAT(CC);	
	#endif

	cgemm(NoTran,NoTran,alpha,CA,CB,beta,CC);

	#if _print
	print_CMAT(CC);	
	#endif
  free_CMAT(CA);
  free_CMAT(CB);
  free_CMAT(CC);
#endif

#if _F
	CA = alloc_CMAT(4,3);
	CB = alloc_CMAT(4,5);
	CC = alloc_CMAT(3,5);

  cfill(CA,1.0,0); 
  cfill(CB,3.0,0); 
  cfill(CC,-0.1,-0.1); 

	#if _print
	print_CMAT(CA);	
	print_CMAT(CB);	
	print_CMAT(CC);	
	#endif

	cgemm(Tran,NoTran,alpha,CA,CB,beta,CC);

	#if _print
	print_CMAT(CC);	
	#endif
  free_CMAT(CA);
  free_CMAT(CB);
  free_CMAT(CC);
#endif

#if _G
	CA = alloc_CMAT(4,3);
	CB = alloc_CMAT(5,4);
	CC = alloc_CMAT(3,5);

  cfill(CA,1.0,0); 
  cfill(CB,3.0,0); 
  cfill(CC,-0.1,-0.1); 

	#if _print
	print_CMAT(CA);	
	print_CMAT(CB);	
	print_CMAT(CC);	
	#endif

	cgemm(Tran,Tran,alpha,CA,CB,beta,CC);

	#if _print
	print_CMAT(CC);	
	#endif
  free_CMAT(CA);
  free_CMAT(CB);
  free_CMAT(CC);
#endif

#if _H
	CA = alloc_CMAT(3,4);
	CB = alloc_CMAT(5,4);
	CC = alloc_CMAT(3,5);

  cfill(CA,1.0,1.11); 
  cfill(CB,3.0,-2.33); 
  cfill(CC,-0.1,-0.1); 

	#if _print
	print_CMAT(CA);	
	print_CMAT(CB);	
	print_CMAT(CC);	
	#endif

	cgemm(NoTran,Tran,alpha,CA,CB,beta,CC);

	#if _print
	print_CMAT(CC);	
	#endif
  free_CMAT(CA);
  free_CMAT(CB);
  free_CMAT(CC);
#endif

#if USE_CUDA
finit();
#endif

	return 0;
}
