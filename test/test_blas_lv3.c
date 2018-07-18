#include "mother.h"

#define _print 1 
	
#define _A 0  //N N
#define _B 0  //T N
#define _C 0  //T T
#define _D 0  //N T

#define _E 0  //N N
#define _F 0  //N T
#define _G 1  //N C
#define _H 0  //T N
#define _I 0  //T T
#define _J 0 //T C
#define _K 0  //C M 
#define _L 0  //C T
#define _M 0  //C C

#define colA 3
#define rowA 4 
#define rowB 5
int main()
{
	MAT*A,*B,*C;
	CMAT* CA, *CB, *CC;
  CTYPE alpha,beta;
  ITER i;
#if USE_CUDA
init();
#endif

	alpha.re=1;
	alpha.im=1;
	beta.re=0;
	beta.im=0;
/**** REAL MATRIX ****/

#if _A
	A = alloc_MAT(colA,rowA);
	B = alloc_MAT(rowA,rowB);
	C = alloc_MAT(colA,rowB);

  fill(A,1.0); 
  fill(B,3.0); 
  fill(C,-0.1); 
for(i=0;i<A->d0*A->d1;i++){A->data[i] = i;}
for(i=0;i<B->d0*B->d1;i++){B->data[i] = i;}

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
	A = alloc_MAT(rowA,colA);
	B = alloc_MAT(rowA,rowB);
	C = alloc_MAT(colA,rowB);

  fill(A,1.0); 
  fill(B,3.0); 
  fill(C,-0.1); 
for(i=0;i<A->d0*A->d1;i++){A->data[i] = i;}
for(i=0;i<B->d0*B->d1;i++){B->data[i] = i;}

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
	A = alloc_MAT(rowA,colA);
	B = alloc_MAT(rowB,rowA);
	C = alloc_MAT(colA,rowB);

  fill(A,1.0); 
  fill(B,3.0); 
  fill(C,-0.1); 
for(i=0;i<A->d0*A->d1;i++){A->data[i] = i;}
for(i=0;i<B->d0*B->d1;i++){B->data[i] = i;}

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
	A = alloc_MAT(colA,rowA);
	B = alloc_MAT(rowB,rowA);
	C = alloc_MAT(colA,rowB);

  fill(A,1.0); 
  fill(B,3.0); 
  fill(C,-0.1); 
for(i=0;i<A->d0*A->d1;i++){A->data[i] = i;}
for(i=0;i<B->d0*B->d1;i++){B->data[i] = i;}

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
	CA = alloc_CMAT(colA,rowA);
	CB = alloc_CMAT(rowA,rowB);
	CC = alloc_CMAT(colA,rowB);

  cfill(CA,1.0,1); 
  cfill(CB,3.0,3); 
  cfill(CC,-0.1,-0.1); 

for(i=0;i<CA->d0*CA->d1;i++){CA->data[i].re = i;CA->data[i].im=i;}
for(i=0;i<CB->d0*CB->d1;i++){CB->data[i].re = i;CB->data[i].im=i;}
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
	CA = alloc_CMAT(colA,rowA);
	CB = alloc_CMAT(rowB,rowA);
	CC = alloc_CMAT(colA,rowB);

  cfill(CA,1.0,1); 
  cfill(CB,3.0,3); 
  cfill(CC,-0.1,-0.1); 
for(i=0;i<CA->d0*CA->d1;i++){CA->data[i].re = i;CA->data[i].im=i;}
for(i=0;i<CB->d0*CB->d1;i++){CB->data[i].re = i;CB->data[i].im=i;}

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

#if _G
	CA = alloc_CMAT(colA,rowA);
	CB = alloc_CMAT(rowB,rowA);
	CC = alloc_CMAT(colA,rowB);

  cfill(CA,1.0,1); 
  cfill(CB,3.0,-3); 
  cfill(CC,-0.1,-0.1); 
for(i=0;i<CA->d0*CA->d1;i++){CA->data[i].re = i;CA->data[i].im=i;}
for(i=0;i<CB->d0*CB->d1;i++){CB->data[i].re = i;CB->data[i].im=i;}

	#if _print
	print_CMAT(CA);	
	print_CMAT(CB);	
	print_CMAT(CC);	
	#endif

	cgemm(NoTran,CTran,alpha,CA,CB,beta,CC);

	#if _print
	print_CMAT(CC);	
	#endif
  free_CMAT(CA);
  free_CMAT(CB);
  free_CMAT(CC);
#endif
#if _H
	CA = alloc_CMAT(rowA,colA);
	CB = alloc_CMAT(rowA,rowB);
	CC = alloc_CMAT(colA,rowB);

  cfill(CA,1.0,1); 
  cfill(CB,3.0,3); 
  cfill(CC,-0.1,-0.1); 
for(i=0;i<CA->d0*CA->d1;i++){CA->data[i].re = i;CA->data[i].im=i;}
for(i=0;i<CB->d0*CB->d1;i++){CB->data[i].re = i;CB->data[i].im=i;}
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
#if _I
	CA = alloc_CMAT(rowA,colA);
	CB = alloc_CMAT(rowB,rowA);
	CC = alloc_CMAT(colA,rowB);

  cfill(CA,1.0,1); 
  cfill(CB,3.0,3); 
  cfill(CC,-0.1,-0.1); 
for(i=0;i<CA->d0*CA->d1;i++){CA->data[i].re = i;CA->data[i].im=i;}
for(i=0;i<CB->d0*CB->d1;i++){CB->data[i].re = i;CB->data[i].im=i;}
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
#if _J
	CA = alloc_CMAT(rowA,colA);
	CB = alloc_CMAT(rowB,rowA);
	CC = alloc_CMAT(colA,rowB);

  cfill(CA,1.0,1); 
  cfill(CB,3.0,-3); 
  cfill(CC,-0.1,-0.1); 
for(i=0;i<CA->d0*CA->d1;i++){CA->data[i].re = i;CA->data[i].im=i;}
for(i=0;i<CB->d0*CB->d1;i++){CB->data[i].re = i;CB->data[i].im=i;}

	#if _print
	print_CMAT(CA);	
	print_CMAT(CB);	
	print_CMAT(CC);	
	#endif

	cgemm(Tran,CTran,alpha,CA,CB,beta,CC);

	#if _print
	print_CMAT(CC);	
	#endif
  free_CMAT(CA);
  free_CMAT(CB);
  free_CMAT(CC);
#endif
#if _K
	CA = alloc_CMAT(rowA,colA);
	CB = alloc_CMAT(rowA,rowB);
	CC = alloc_CMAT(colA,rowB);

  cfill(CA,1.0,-1); 
  cfill(CB,1.0,1); 
  cfill(CC,-0.1,-0.1); 
for(i=0;i<CA->d0*CA->d1;i++){CA->data[i].re = i;CA->data[i].im=i;}
for(i=0;i<CB->d0*CB->d1;i++){CB->data[i].re = i;CB->data[i].im=i;}
	#if _print
	print_CMAT(CA);	
	print_CMAT(CB);	
	print_CMAT(CC);	
	#endif

	cgemm(CTran,NoTran,alpha,CA,CB,beta,CC);

	#if _print
	print_CMAT(CC);	
	#endif
  free_CMAT(CA);
  free_CMAT(CB);
  free_CMAT(CC);
#endif
#if _L
	CA = alloc_CMAT(rowA,colA);
	CB = alloc_CMAT(rowB,rowA);
	CC = alloc_CMAT(colA,rowB);
  
  cfill(CA,1.0,-1); 
  cfill(CB,1.0,1); 
  cfill(CC,-0.1,-0.1); 
for(i=0;i<CA->d0*CA->d1;i++){CA->data[i].re = i;CA->data[i].im=i;}
for(i=0;i<CB->d0*CB->d1;i++){CB->data[i].re = i;CB->data[i].im=i;}

	#if _print
	print_CMAT(CA);	
	print_CMAT(CB);	
	print_CMAT(CC);	
	#endif

	cgemm(CTran,Tran,alpha,CA,CB,beta,CC);

	#if _print
	print_CMAT(CC);	
	#endif
  free_CMAT(CA);
  free_CMAT(CB);
  free_CMAT(CC);
#endif
#if _M
	CA = alloc_CMAT(rowA,colA);
	CB = alloc_CMAT(rowB,rowA);
	CC = alloc_CMAT(colA,rowB);
  
  cfill(CA,1.0,-1); 
  cfill(CB,1.0,-1); 
  cfill(CC,-0.1,-0.1); 

for(i=0;i<CA->d0*CA->d1;i++){CA->data[i].re = i;CA->data[i].im=i;}
for(i=0;i<CB->d0*CB->d1;i++){CB->data[i].re = i;CB->data[i].im=i;}

	#if _print
	print_CMAT(CA);	
	print_CMAT(CB);	
	print_CMAT(CC);	
	#endif

	cgemm(CTran,CTran,alpha,CA,CB,beta,CC);

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
