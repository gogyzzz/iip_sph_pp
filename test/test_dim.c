#include "mother.h"

int main()
{
  MAT*A;
  ITER i;
  DIM* d;
  
  d=(DIM*)malloc(sizeof(DIM));
  printf("==== repmat ====\n");
  d->d0 = 3;
  d->d1= 2;
  d->d2 = 2;
  A= zeros(1,5);
  for(i=0;i<5;i++)
   A->data[i] = i+1;
  print_MAT(A);
  repmat(A,d);
  print_MAT(A);
  free_MAT(A);

  printf("==== shiftdim(+1) ====\n");
  A= zeros(2,1,5);
  for(i=0;i<10;i++)
   A->data[i] = i*10;
  print_MAT(A);
  shiftdim(A,1);
  print_MAT(A);
  free_MAT(A);

  printf("==== shiftdim(+1) ====\n");
  A= zeros(3,3,2);
  for(i=0;i<18;i++)
   A->data[i] = i;
  print_MAT(A);
  shiftdim(A,1);
  print_MAT(A);
  free_MAT(A);

}
