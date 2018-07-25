#include "mother.h"
#include <stdio.h>

int main()
{
	MAT*A;
	CMAT*CA;
	ITER i;

	A = zeros(5,5);
	CA = czeros(5,5);

	for(i=0;i<25;i++)
		A->data[i]=i+1;
	for(i=0;i<25;i++)
		{CA->data[i].re=i+1;CA->data[i].im = i+1;}

	inv_elements(A);
	cinv_elements(CA);

	print_MAT(A);
	free_MAT(A);

	print_CMAT(CA);
	free_CMAT(CA);
	return 0;
}

