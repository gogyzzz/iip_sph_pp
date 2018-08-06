#include "mother.h"

int main() {
  MAT* A;
  CMAT* CA;

  A = zeros(5, 4);
  CA = czeros(5, 4);

  printf("==== round ====\n");
  randu(A, -10., 10.);
  crandu(CA, -10., 10., -10., 10.);
  print_mat(A);
  print_cmat(CA);
  round_mat(A);
  round_cmat(CA);
  print_mat(A);
  print_cmat(CA);

  printf("==== floor ====\n");
  randu(A, -10., 10.);
  crandu(CA, -10., 10., -10., 10.);
  print_mat(A);
  print_cmat(CA);
  floor_mat(A);
  floor_cmat(CA);
  print_mat(A);
  print_cmat(CA);

  printf("==== ceil ====\n");
  randu(A, -10., 10.);
  crandu(CA, -10., 10., -10., 10.);
  print_mat(A);
  print_cmat(CA);
  ceil_mat(A);
  ceil_cmat(CA);
  print_mat(A);
  print_cmat(CA);

  printf("==== log ====\n");
  randu(A, -10., 10.);
  crandu(CA, -10., 10., -10., 10.);
  print_mat(A);
  print_cmat(CA);
  log_mat(A);
  log_cmat(CA);
  print_mat(A);
  print_cmat(CA);

  printf("==== log2 ====\n");
  randu(A, -10., 10.);
  crandu(CA, -10., 10., -10., 10.);
  print_mat(A);
  print_cmat(CA);
  log2_mat(A);
  log2_cmat(CA);
  print_mat(A);
  print_cmat(CA);

  printf("==== log10 ====\n");
  randu(A, -10., 10.);
  crandu(CA, -10., 10., -10., 10.);
  print_mat(A);
  print_cmat(CA);
  log10_mat(A);
  log10_cmat(CA);
  print_mat(A);
  print_cmat(CA);

  printf("==== exp ====\n");
  randu(A, -10., 10.);
  crandu(CA, -10., 10., -10., 10.);
  print_mat(A);
  print_cmat(CA);
  exp_mat(A);
  exp_cmat(CA);
  print_mat(A);
  print_cmat(CA);

  return 0;
}
