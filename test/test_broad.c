#include "mother.h"

int main() {
  MAT *A, *B, *C;
  CMAT *CA, *CB, *CC;
  ITER i;
  int cnt = 0;
  printf("==== %d =====\n", cnt++);
  A = zeros(1, 1);
  B = zeros(1, 3);
  C = zeros(1, 3);
  fill(A, 2);
  for (i = 0; i < 3; i++) B->data[i] = i + 1;
  mul_elements(A, B, C);
  print_MAT(A);
  print_MAT(B);
  print_MAT(C);
  free_MAT(A);
  free_MAT(B), free_MAT(C);

  printf("==== %d =====\n", cnt++);
  A = zeros(1, 1);
  B = zeros(3, 3);
  C = zeros(3, 3);
  fill(A, 2);
  for (i = 0; i < 9; i++) B->data[i] = i + 1;
  mul_elements(A, B, C);
  print_MAT(A);
  print_MAT(B);
  print_MAT(C);
  free_MAT(A);
  free_MAT(B), free_MAT(C);

  printf("==== %d =====\n", cnt++);
  A = zeros(1, 3);
  B = zeros(3, 1);
  C = zeros(3, 3);
  for (i = 0; i < 3; i++) B->data[i] = i + 1;
  for (i = 0; i < 3; i++) A->data[i] = i + 1;
  mul_elements(A, B, C);
  print_MAT(A);
  print_MAT(B);
  print_MAT(C);
  free_MAT(A);
  free_MAT(B), free_MAT(C);

  printf("==== %d =====\n", cnt++);
  A = zeros(3, 1);
  B = zeros(1, 3);
  C = zeros(3, 3);
  for (i = 0; i < 3; i++) B->data[i] = i + 1;
  for (i = 0; i < 3; i++) A->data[i] = i + 1;
  mul_elements(A, B, C);
  print_MAT(A);
  print_MAT(B);
  print_MAT(C);
  free_MAT(A);
  free_MAT(B), free_MAT(C);

  printf("==== %d =====\n", cnt++);
  A = zeros(3, 3);
  B = zeros(1, 3);
  C = zeros(3, 3);
  for (i = 0; i < 3; i++) B->data[i] = i + 1;
  for (i = 0; i < 9; i++) A->data[i] = i + 1;
  mul_elements(A, B, C);
  print_MAT(A);
  print_MAT(B);
  print_MAT(C);
  free_MAT(A);
  free_MAT(B), free_MAT(C);

  printf("==== %d =====\n", cnt++);
  A = zeros(3, 1);
  B = zeros(3, 3);
  C = zeros(3, 3);
  for (i = 0; i < 9; i++) B->data[i] = i + 1;
  for (i = 0; i < 3; i++) A->data[i] = i + 1;
  mul_elements(A, B, C);
  print_MAT(A);
  print_MAT(B);
  print_MAT(C);
  free_MAT(A);
  free_MAT(B), free_MAT(C);

  printf("==== %d =====\n", cnt++);
  A = zeros(3, 1,3);
  B = zeros(3, 3);
  C = zeros(3, 3,3);
  for (i = 0; i < 9; i++) B->data[i] = i + 1;
  for (i = 0; i < 9; i++) A->data[i] = i + 1;
  mul_elements(A, B, C);
  print_MAT(A);
  print_MAT(B);
  print_MAT(C);
  free_MAT(A);
  free_MAT(B), free_MAT(C);


  /**** COMPLEX ****/

  printf("==== %d =====\n", cnt++);
  CA = czeros(1, 3);
  CB = czeros(3, 1);
  CC = czeros(3, 3);
  for (i = 0; i < 3; i++) {
    CB->data[i].re = i + 1;
    CB->data[i].im = i + 3;
  };
  for (i = 0; i < 3; i++) {
    CA->data[i].re = i + 1;
    CA->data[i].im = i + 3;
  };
  cmul_elements(CA, CB, CC);
  print_CMAT(CA);
  print_CMAT(CB);
  print_CMAT(CC);
  free_CMAT(CA);
  free_CMAT(CB), free_CMAT(CC);

  printf("==== %d =====\n", cnt++);
  CA = czeros(3, 3);
  CB = czeros(1, 3);
  CC = czeros(3, 3);
  for (i = 0; i < 3; i++) {
    CB->data[i].re = i + 1;
    CB->data[i].im = i + 3;
  };
  for (i = 0; i < 9; i++) {
    CA->data[i].re = i + 1;
    CA->data[i].im = i + 3;
  };
  cmul_elements(CA, CB, CC);
  print_CMAT(CA);
  print_CMAT(CB);
  print_CMAT(CC);
  free_CMAT(CA);
  free_CMAT(CB), free_CMAT(CC);

   printf("==== %d =====\n", cnt++);
  CA = czeros(4, 1,2);
  CB = czeros(4, 4);
  CC = czeros(4, 4,2);
  for (i = 0; i < 16; i++) {
    CB->data[i].re = i + 1;
    CB->data[i].im = i + 3;
  };
  for (i = 0; i < 8; i++) {
    CA->data[i].re = i + 1;
    CA->data[i].im = i + 3;
  };
  cadd_elements(CA, CB, CC);
  print_CMAT(CA);
  print_CMAT(CB);
  print_CMAT(CC);
  free_CMAT(CA);
  free_CMAT(CB), free_CMAT(CC);


  return 0;
}
