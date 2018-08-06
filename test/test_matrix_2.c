#include "mother.h"

#define _A 1  // real number
#define _B 1  // complex number

#define _print 1
int main() {
  MAT *A, *TA;
  CMAT *CA, *TCA, *HCA;

/**** REAL ****/
#if _A
  // 1D
  A = alloc_mat(8);
  fill(A, 8);
  TA = trans(A);
#if _print
  print_mat(A);
  print_mat(TA);
#endif
  free_mat(A);
  free_mat(TA);
  // 2D
  A = alloc_mat(5, 2);
  fill(A, 25);
  TA = trans(A);
#if _print
  print_mat(A);
  print_mat(TA);
#endif
  free_mat(A);
  free_mat(TA);
  // 3D as batched 2D
  A = alloc_mat(2, 3, 4);
  fill(A, 1);
  TA = trans(A);
#if _print
  print_mat(A);
  print_mat(TA);
#endif
  free_mat(A);
  free_mat(TA);
#endif

/**** COMPLEX ****/
#if _B
  // 1D
  CA = alloc_cmat(4);
  cfill(CA, 4, -4);
  TCA = ctrans(CA);
  HCA = hermit(CA);
#if _print
  print_cmat(CA);
  print_cmat(TCA);
  print_cmat(HCA);
#endif
  free_cmat(CA);
  free_cmat(TCA);
  free_cmat(HCA);
  // 2D
  CA = alloc_cmat(2, 4);
  cfill(CA, 2, -4);
  TCA = ctrans(CA);
  HCA = hermit(CA);
#if _print
  print_cmat(CA);
  print_cmat(TCA);
  print_cmat(HCA);
#endif
  free_cmat(CA);
  free_cmat(TCA);
  free_cmat(HCA);
  // 3D as batched 2D
  CA = alloc_cmat(2, 3, 4);
  cfill(CA, 2, -3);
  TCA = ctrans(CA);
  HCA = hermit(CA);
#if _print
  print_cmat(CA);
  print_cmat(TCA);
  print_cmat(HCA);
#endif
  free_cmat(CA);
  free_cmat(TCA);
  free_cmat(HCA);

#endif

  return 0;
}
