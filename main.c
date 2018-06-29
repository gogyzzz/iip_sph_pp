#include "mother.h"

int main()
{
	MAT* a1;
	MAT* a2;

	a1 = zeros(4);
	a2 = zeros(4,4);

	print_MAT(a1);
	print_MAT(a2);

	free_MAT(a1);
	free_MAT(a2);

	return 0;
}
