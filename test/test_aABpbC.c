#include "mother.h"

#define colA 2
#define rowA 6
#define rowB 10

int main()
{
	CTYPE alpha,beta;
	MAT *A,*TA,*B,*TB,*C;
	CMAT *CA,*TCA,*CB,*TCB,*CC;
  int i;
	

A =	alloc_MAT(colA,rowA);
TA=	alloc_MAT(rowA,colA);
B =	alloc_MAT(rowA,rowB);
TB=	alloc_MAT(rowB,rowA);
C =	alloc_MAT(colA,rowB);

CA =	alloc_CMAT(colA,rowA);
TCA=	alloc_CMAT(rowA,colA);
CB =	alloc_CMAT(rowA,rowB);
TCB=	alloc_CMAT(rowB,rowA);
CC =	alloc_CMAT(colA,rowB);

	for(i=0;i<colA*rowA;i++)A->data[i]=1;
	for(i=0;i<colA*rowA;i++)TA->data[i]=1;
	for(i=0;i<rowB*rowA;i++)B->data[i]=1;
	for(i=0;i<rowB*rowA;i++)TB->data[i]=1;
	for(i=0;i<colA*rowB;i++)C->data[i]=1;

	for(i=0;i<colA*rowA;i++){CA->data[i].re=1;CA->data[i].im = 2;}
	for(i=0;i<colA*rowA;i++){TCA->data[i].re=1;TCA->data[i].im = 2;}
	for(i=0;i<rowB*rowA;i++){CB->data[i].re=1;CB->data[i].im = 2;}
	for(i=0;i<rowB*rowA;i++){TCB->data[i].re=1;TCB->data[i].im = 2;}
	for(i=0;i<colA*rowB;i++){CC->data[i].re=1;CC->data[i].im = 2;}

/**** REAL ****/

	aABpbC(1,A,B,0,C);
  print_MAT(C);

	aABtpbC(1,A,TB,0,C);
  print_MAT(C);

	aAtBpbC(1,TA,B,0,C);
  print_MAT(C);

	aAtBpbC(1,TA,B,0,C);
  print_MAT(C);

  matmul(A,B,C);
	print_MAT(C);

/**** COMPLEX ****/

	alpha.re=1;
	beta.re=0;
	alpha.im=0;
	beta.im=0;

	caABpbC(alpha,CA,CB,beta,CC);
  print_CMAT(CC);
	caABtpbC(alpha,CA,TCB,beta,CC);
  print_CMAT(CC);
	caABhpbC(alpha,CA,TCB,beta,CC);
  print_CMAT(CC);

	caAtBpbC(alpha,TCA,CB,beta,CC);
  print_CMAT(CC);
	caAtBtpbC(alpha,TCA,TCB,beta,CC);
  print_CMAT(CC);
	caAtBhpbC(alpha,TCA,TCB,beta,CC);
  print_CMAT(CC);

	caAhBpbC(alpha,TCA,CB,beta,CC);
  print_CMAT(CC);
	caAhBtpbC(alpha,TCA,TCB,beta,CC);
  print_CMAT(CC);
	caAhBhpbC(alpha,TCA,TCB,beta,CC);
  print_CMAT(CC);

	cmatmul(CA,CB,CC);
	print_CMAT(CC);

	free_MAT(A);
	free_MAT(B);
	free_MAT(TA);
	free_MAT(TB);
	free_MAT(C);

	free_CMAT(CA);
	free_CMAT(CB);
	free_CMAT(TCA);
	free_CMAT(TCB);
	free_CMAT(CC);


	return 0;
}
