#include "mother.h"

#define m 3
#define n 4
#define k 5

int main()
{
	MAT*A;
	MAT*subA;

	A = mem_MAT(m,n);

	fill(A,1);

	row_add(-10,A,1);
	col_scal(3,A,1);

	print_MAT(A);

	subA = mem_submat(A,-1,-1,0,3);

	print_MAT(subA);

	free_mem_MAT(subA);
	free_mem_MAT(A);
	return 0;
}
