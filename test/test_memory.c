#include "mother.h"

#define _print 0

#define m 400
#define n 500
#define k 600

int main()
{
	MAT*A,*B,*C;
	ITER i,j;
	init();

	for(j=0;j<200;j++)
	{
	A = iip_malloc(sizeof(MAT));
	B = iip_malloc(sizeof(MAT));
	C= iip_malloc(sizeof(MAT));
	
	A->data = iip_malloc(sizeof(DTYPE)*m*n);
	B->data = iip_malloc(sizeof(DTYPE)*n*k);
	C->data = iip_malloc(sizeof(DTYPE)*m*k);

	A->ndim = 1;
	B->ndim = 1;
	C->ndim = 1;

	A->d2 = 1;
	B->d2 = 1;
	C->d2 = 1;

	A->d1 = n;
	A->d0 = m;
	B->d1 = k;
	B->d0 = n;
	C->d1 = k;
	C->d0 = m;

	fill(A,1);
	fill(B,2);
	fill(C,0);
#if _print
	print_MAT(A);
	print_MAT(B);
	print_MAT(C);
#endif
	printf("[%d X %d] = [%d X %d] * [%d X %d]\n",m,k,m,n,n,k);
	matmul(A,B,C);
#if _print
	print_MAT(A);
	print_MAT(B);
	print_MAT(C);
#endif
	iip_free(A->data);
	iip_free(A);
	iip_free(C->data);
	iip_free(B->data);
	iip_free(C);
	iip_free(B);

	A = iip_malloc(sizeof(MAT));
	A->data = iip_malloc(sizeof(DTYPE)*20*100);
	B = iip_malloc(sizeof(MAT));
	C = iip_malloc(sizeof(MAT));
	C->data = iip_malloc(sizeof(DTYPE)*24*1);
	B->data = iip_malloc(sizeof(DTYPE)*10*133);
	
	iip_free(B->data);
	iip_free(A->data);
	iip_free(A);
	iip_free(B);
	iip_free(C->data);
	iip_free(C);
//	A = mem_MAT(10);

//	free_mem_MAT(10);
	}
	finit();
	return 0 ;
}
