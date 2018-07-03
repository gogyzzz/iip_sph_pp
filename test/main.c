#include "mother.h"

#define PART 2
/*
 * 1: Dmatrix
 * 2: Cmatrix
 * 3: 
 *
 *
 * */

int main()
{
	//printf("%s | %s\n",str(DTYPE),xstr(DTYPE));	

  ITER i,j,k;
	
#if PART == 1	
	MAT* a1;
  MAT* a2;
	MAT* a3;
  MAT* a2sub;

  a1 = zeros(4);
  a2 = zeros(3,4);
	a3 = zeros(2,3,4);
  a2sub = zeros(1, 2);

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
	for(i=0;i<b1->d0;i++)
		set(b1,i,-123);
	printf("b1\n");
	print_MAT(b1);
	printf("axpy(200,a1,b1)\n");
	axpy(999,a1,b1);
	print_MAT(b1);

  printf("  fill b1 0.5 and print\n");
  fill(b1, 0.5);
  print_MAT(b1); // [0.5, 0.5, 0.5, 0.5]' expected

  printf("  get submat\n");
  for(i=0; i<3; i++)
    for(j=0; j<4; j++)
      set(a2, i, j, (DTYPE)((i+1) + a2->d0*j));
  print_MAT(a2);

  // expected below
  // 1.000000 4.000000 7.000000 10.000000
  // 2.000000 5.000000 8.000000 11.000000
  // 3.000000 6.000000 9.000000 12.00000

  submat(a2, a2sub,  1,2,  1,3,  -1,-1);
  print_MAT(a2sub); // [5.0, 8.0] expected


	free_MAT(a1);
  free_MAT(a2);
	free_MAT(a3);
  free_MAT(a2sub);

	free_MAT(b1);

#elif PART==2 
	
	CMAT* a1;
  CMAT* a2;
	CMAT* a3;
  CMAT* a2sub;

  a1 = czeros(4);
  a2 = czeros(3,4);
	a3 = czeros(2,3,4);
  a2sub = czeros(1, 2);

	printf("a1\n");
  print_CMAT(a1);

	printf("a2\n");
  print_CMAT(a2);

	printf("a3\n");
	print_CMAT(a3);

	cset(a1,2,1,1);
	cset(a2,1,2,1,1);
	cset(a3,0,1,2,1,1);

	printf("a1\n");
  print_CMAT(a1);

	printf("a2\n");
  print_CMAT(a2);

	printf("a3\n");
	print_CMAT(a3);

	CMAT * b1;

	b1 = czeros(4);
	for(i=0;i<b1->d0;i++)
		cset(b1,i,-123,-123);
	printf("b1\n");
	print_CMAT(b1);
	printf("caxpy(200,a1,b1)\n");
	caxpy(999,a1,b1);
	print_CMAT(b1);

  printf("  fill b1 0.5 and print\n");
  cfill(b1, 0.5,0.5);
  print_CMAT(b1); // [0.5, 0.5, 0.5, 0.5]' expected

  printf("  get csubmat\n");
  for(i=0; i<3; i++)
    for(j=0; j<4; j++)
      cset(a2, i, j, (DTYPE)((i+1) + a2->d0*j),(DTYPE)((i+1) + a2->d0*j) );
  print_CMAT(a2);

  csubmat(a2, a2sub,  1,2,  1,3,  -1,-1);
  print_CMAT(a2sub); // [5.0, 8.0] expected

	free_CMAT(a1);
  free_CMAT(a2);
	free_CMAT(a3);
  free_CMAT(a2sub);

	free_CMAT(b1);
#endif	

	return 0;
}




