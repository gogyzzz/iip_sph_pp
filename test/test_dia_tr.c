#include "mother.h"

int main() {
  MAT* A;
  CMAT* CA;
  MAT *dia, *tr;
  CMAT *cdia, *ctr;
  MAT* NV;
  ITER i;

  A = alloc_mat(2, 2, 2);
  CA = czeros(2, 2, 2);

  dia = alloc_mat(2, 1, 2);
  cdia = alloc_cmat(2, 1, 2);

  tr = alloc_mat(1, 1, 2);
  ctr = alloc_cmat(1, 1, 2);

  NV = alloc_mat(3, 2, 2);

  for (i = 0; i < 8; i++) A->data[i] = i + 1;
  for (i = 0; i < 8; i++) CA->data[i].re = i + 1;

  print_mat(A);
  diagonal(A, dia);
  diagonal(A, NV);
  trace(A, NV);
  trace(A, tr);

  print_mat(dia);
  print_mat(tr);

  crandu(CA, 0., 10., 0., 10.);

  print_cmat(CA);
  cdiagonal(CA, cdia);
  ctrace(CA, ctr);
  print_cmat(cdia);
  print_cmat(ctr);
  return 0;
}
