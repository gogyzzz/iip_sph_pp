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
  print_mat(A);
  print_mat(B);
  print_mat(C);
  free_mat(A);
  free_mat(B), free_mat(C);

  printf("==== %d =====\n", cnt++);
  A = zeros(1, 1);
  B = zeros(3, 3);
  C = zeros(3, 3);
  fill(A, 2);
  for (i = 0; i < 9; i++) B->data[i] = i + 1;
  mul_elements(A, B, C);
  print_mat(A);
  print_mat(B);
  print_mat(C);
  free_mat(A);
  free_mat(B), free_mat(C);

  printf("==== %d =====\n", cnt++);
  A = zeros(1, 3);
  B = zeros(3, 1);
  C = zeros(3, 3);
  for (i = 0; i < 3; i++) B->data[i] = i + 1;
  for (i = 0; i < 3; i++) A->data[i] = i + 1;
  mul_elements(A, B, C);
  print_mat(A);
  print_mat(B);
  print_mat(C);
  free_mat(A);
  free_mat(B), free_mat(C);

  printf("==== %d =====\n", cnt++);
  A = zeros(3, 1);
  B = zeros(1, 3);
  C = zeros(3, 3);
  for (i = 0; i < 3; i++) B->data[i] = i + 1;
  for (i = 0; i < 3; i++) A->data[i] = i + 1;
  mul_elements(A, B, C);
  print_mat(A);
  print_mat(B);
  print_mat(C);
  free_mat(A);
  free_mat(B), free_mat(C);

  printf("==== %d =====\n", cnt++);
  A = zeros(3, 3);
  B = zeros(1, 3);
  C = zeros(3, 3);
  for (i = 0; i < 3; i++) B->data[i] = i + 1;
  for (i = 0; i < 9; i++) A->data[i] = i + 1;
  mul_elements(A, B, C);
  print_mat(A);
  print_mat(B);
  print_mat(C);
  free_mat(A);
  free_mat(B), free_mat(C);

  printf("==== %d =====\n", cnt++);
  A = zeros(3, 1);
  B = zeros(3, 3);
  C = zeros(3, 3);
  for (i = 0; i < 9; i++) B->data[i] = i + 1;
  for (i = 0; i < 3; i++) A->data[i] = i + 1;
  mul_elements(A, B, C);
  print_mat(A);
  print_mat(B);
  print_mat(C);
  free_mat(A);
  free_mat(B), free_mat(C);

  printf("==== %d =====\n", cnt++);
  A = zeros(3, 1, 3);
  B = zeros(3, 3);
  C = zeros(3, 3, 3);
  for (i = 0; i < 9; i++) B->data[i] = i + 1;
  for (i = 0; i < 9; i++) A->data[i] = i + 1;
  mul_elements(A, B, C);
  print_mat(A);
  print_mat(B);
  print_mat(C);
  free_mat(A);
  free_mat(B), free_mat(C);

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
  print_cmat(CA);
  print_cmat(CB);
  print_cmat(CC);
  free_cmat(CA);
  free_cmat(CB), free_cmat(CC);

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
  print_cmat(CA);
  print_cmat(CB);
  print_cmat(CC);
  free_cmat(CA);
  free_cmat(CB), free_cmat(CC);

  printf("==== %d =====\n", cnt++);
  CA = czeros(4, 1, 2);
  CB = czeros(4, 4);
  CC = czeros(4, 4, 2);
  for (i = 0; i < 16; i++) {
    CB->data[i].re = i + 1;
    CB->data[i].im = i + 3;
  };
  for (i = 0; i < 8; i++) {
    CA->data[i].re = i + 1;
    CA->data[i].im = i + 3;
  };
  cadd_elements(CA, CB, CC);
  print_cmat(CA);
  print_cmat(CB);
  print_cmat(CC);
  free_cmat(CA);
  free_cmat(CB), free_cmat(CC);

  return 0;
}
