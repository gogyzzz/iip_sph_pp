#include "mother.h"

int main()
{
	MAT * A;
	ITER i;
	A = alloc_MAT(3,5);

	for(i=0;i<A->d0 * A-> d1; i++)A->data[i]=1;
	print_MAT(A);

	col_scal(A,1,-40);
	row_scal(A,2,-10);

	print_MAT(A);

	free_MAT(A);


	return 0;
}
