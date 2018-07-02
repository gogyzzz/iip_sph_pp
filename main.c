#include "mother.h"

int main()
{
	//printf("%s | %s\n",str(DTYPE),xstr(DTYPE));	

	
	MAT* a1;
    MAT* a2;
	MAT* a3;

    a1 = zeros(4);
    a2 = zeros(3,4);
	a3 = zeros(2,3,4);

	printf("a1\n");
    print_MAT(a1);

	printf("a2\n");
    print_MAT(a2);

	printf("a3\n");
	print_MAT(a3);

	set(a1,2,1);
	set(a2,1,2,1);
	set(a3,0,1,2,1);

	printf("a1\n");
    print_MAT(a1);

	printf("a2\n");
    print_MAT(a2);

	printf("a3\n");
	print_MAT(a3);

	MAT * b1;

	b1 = zeros(4);
	for(int i=0;i<b1->d0;i++)
		set(b1,i,-123);
	printf("b1\n");
	print_MAT(b1);
	printf("axpy(200,a1,b1)\n");
	axpy(999,a1,b1);
	print_MAT(b1);


	free_MAT(a1);
    free_MAT(a2);
	free_MAT(a3);

	free_MAT(b1);
	return 0;
}




