#include "mother.h"

int main()
{
	MAT *A,*B;
	CMAT *CA,*CB;
	CTYPE ctemp={5,-7};  

	init();

	printf("===== init =====\nA=0,B=10\n");

	A  = zeros(2,3);
	CA = czeros(2,3);
	
	B  = alloc_MAT(2,3);
	CB = alloc_CMAT(2,3);
	
	printf("===== A,B ,CA,CB  =====\n");

	fill(B,10);
	cfill(CB,10,10);

	print_MAT(A);
	print_MAT(B);
	print_CMAT(CA);
	print_CMAT(CB);
	printf("===== axpy(5,B,A) =====\n");

	axpy(5,B,A);
	caxpy(ctemp,CB,CA);

	printf("===== print A,CA =====\n");
	print_MAT(A);
	print_CMAT(CA);
	
	printf("===== copy(A,B) =====\n");

	copy(A,B);
	ccopy(CA,CB);

	printf("===== print B,CB =====\n");
	print_MAT(B);
	print_CMAT(CB);
	

	printf("===== asum(B,1) =====\n");
	printf("asum(B,1) = %.2f\n", asum(B, 1));
	printf("casum(CB,1) = %.2f\n", casum(CB, 1));
	
	
	free_MAT(A);
	free_MAT(B);

	free_CMAT(CA);
	free_CMAT(CB); 

	finit();

	return 0;
}
