#include "mother.h"

int main() {
  MAT* A;
  CMAT* CA;
  MAT *dia, *tr;
  CMAT *cdia, *ctr;
  MAT* NV;
  ITER i;

  A = alloc_MAT(2, 2, 2);
  CA = czeros(2, 2, 2);

  dia = alloc_MAT(2, 1, 2);
  cdia = alloc_CMAT(2, 1, 2);

  tr = alloc_MAT(1, 1, 2);
  ctr = alloc_CMAT(1, 1, 2);

  NV = alloc_MAT(3, 2, 2);

  for (i = 0; i < 8; i++) A->data[i] = i + 1;
  for (i = 0; i < 8; i++) CA->data[i].re = i + 1;

  print_MAT(A);
  diagonal(A, dia);
  diagonal(A, NV);
  trace(A, NV);
  trace(A, tr);

  print_MAT(dia);
  print_MAT(tr);

  crandu(CA, 0., 10., 0., 10.);

  print_CMAT(CA);
  cdiagonal(CA, cdia);
  ctrace(CA, ctr);
  print_CMAT(cdia);
  print_CMAT(ctr);
  return 0;
}
