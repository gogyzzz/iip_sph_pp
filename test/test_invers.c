#include "mother.h"

int main()
{
	MAT *A,*B,*C;
	ITER i,j;
	A = alloc_MAT(4,4);
	B = alloc_MAT(4,4);
	C = alloc_MAT(4,4);

	for(i=0;i<16;i++)
		A->data[i]=i;
		
	print_MAT(A);

	copy(A,B);
	invers(B);
	print_MAT(B);

	matmul(A,B,C);
	print_MAT(C);


	return 0;
}
