#include "mother.h"
#include "header_for_test.h"

#define MAT_SIZE 500
#define BATCH_SIZE 50
int main()
{

	MAT **A,**B,**C;
  ITER i;

A = (MAT**)malloc(sizeof(MAT*)*BATCH_SIZE);
B = (MAT**)malloc(sizeof(MAT*)*BATCH_SIZE);
C = (MAT**)malloc(sizeof(MAT*)*BATCH_SIZE);

for(i=0;i<BATCH_SIZE;i++)
{
	A[i] = zeros(MAT_SIZE,MAT_SIZE);
	B[i] = zeros(MAT_SIZE,MAT_SIZE);
	C[i] = zeros(MAT_SIZE,MAT_SIZE);
}	
for(i=0;i<BATCH_SIZE;i++)
{	
	fill(A[i],10.10);
	fill(B[i],-10.10);
	fill(C[i],0);
}
	printf("JUST %d batch gemm (%d,%d) : ",BATCH_SIZE,MAT_SIZE,MAT_SIZE);
	stopwatch(0);
	for(i=0;i<BATCH_SIZE;i++)
		gemm(NoTran,NoTran,1,A[i],B[i],0,C[i]);
	stopwatch(1);

for(i=0;i<BATCH_SIZE;i++)
{	
	fill(A[i],10.10);
	fill(B[i],-10.10);
	fill(C[i],0);
}

	printf("OMP  %d batch gemm (%d,%d) : ",BATCH_SIZE,MAT_SIZE,MAT_SIZE);
	stopwatch(0);
#pragma omp parallel for shared(A,B,C) private(i)
	for(i=0;i<BATCH_SIZE;i++)
		gemm(NoTran,NoTran,1,A[i],B[i],0,C[i]);
	stopwatch(1);

	for(i=0;i<BATCH_SIZE;i++)
{
	free_MAT(A[i]);
	free_MAT(B[i]);
	free_MAT(C[i]);
}
	free(A);
	free(B);
	free(C);
	return 0;
}
