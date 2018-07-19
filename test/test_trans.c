#include "mother.h"

int main()
{
	MAT*A;
	ITER i;

	A=alloc_MAT(4);
	for(i=0;i<A->d0 * A->d1 * A->d2; i++)
		A->data[i] = i;
	print_MAT(A);
	trans(A);
	print_MAT(A);
	free_MAT(A);

	A=alloc_MAT(4,3);
	for(i=0;i<A->d0 * A->d1 * A->d2; i++)
		A->data[i] = i;
	print_MAT(A);
	trans(A);
	print_MAT(A);
	free_MAT(A);
/*
	A=alloc_MAT(4,3,2);
	for(i=0;i<A->d0 * A->d1 * A->d2; i++)
		A->data[i] = i;
	print_MAT(A);
	trans(A);
	print_MAT(A);
	free_MAT(A);
	*/
	return 0;
}
