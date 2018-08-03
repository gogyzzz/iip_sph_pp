#include "mother.h"

int main()
{
  CMAT *CA;
  CMAT *CB;
  CMAT *CC;

  CA = czeros(3,3);
  CB = czeros(3,3);
  CC = czeros(3,3);

  CA->data[0] = CX(4.,2.);
  CA->data[1] = CX(0.,2.);
  CA->data[2] = CX(2.,2.);

  CA->data[3] = CX(3,-4.);
  CA->data[4] = CX(4.,0.);
  CA->data[5] = CX(5.,-1.);

  CA->data[6] = CX(0,1);
  CA->data[7] = CX(3,-1);
  CA->data[8] = CX(0,0);
  print_CMAT(CA);
  cinvert(CA,CB);
  cmatmul(CA,CB,CC);
  print_CMAT(CB);
  print_CMAT(CC);

  free_CMAT(CA);
  free_CMAT(CB);
  free_CMAT(CC);

  CA = czeros(4,4);
  CB = czeros(4,4);
  CC = czeros(4,4);

  CA->data[0] = CX(4.,2.);
  CA->data[1] = CX(0.,2.);
  CA->data[2] = CX(2.,2.);
  CA->data[3] = CX(1.,0.);

  CA->data[4] = CX(3,-4.);
  CA->data[5] = CX(4.,0);
  CA->data[6] = CX(5.,-1.);
  CA->data[7] = CX(0.,2.);

  CA->data[8] = CX(0,1);
  CA->data[9] = CX(3,-1);
  CA->data[10] = CX(0,0);
  CA->data[11] = CX(3,0);

  CA->data[12] = CX(2,2);
  CA->data[13] = CX(0,-5);
  CA->data[14] = CX(3,1);
  CA->data[15] = CX(0,4);
  print_CMAT(CA);
  cinvert(CA,CB);
  //cinvert_nbyn(CA,CB);
  cmatmul(CA,CB,CC);
  print_CMAT(CB);
  print_CMAT(CC);

  free_CMAT(CA);
  free_CMAT(CB);
  free_CMAT(CC);




  return 0;
}


