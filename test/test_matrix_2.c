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
  A = alloc_MAT(8);
  fill(A, 8);
  TA = trans(A);
#if _print
  print_MAT(A);
  print_MAT(TA);
#endif
  free_MAT(A);
  free_MAT(TA);
  // 2D
  A = alloc_MAT(5, 2);
  fill(A, 25);
  TA = trans(A);
#if _print
  print_MAT(A);
  print_MAT(TA);
#endif
  free_MAT(A);
  free_MAT(TA);
  // 3D as batched 2D
  A = alloc_MAT(2, 3, 4);
  fill(A, 1);
  TA = trans(A);
#if _print
  print_MAT(A);
  print_MAT(TA);
#endif
  free_MAT(A);
  free_MAT(TA);
#endif

/**** COMPLEX ****/
#if _B
  // 1D
  CA = alloc_CMAT(4);
  cfill(CA, 4, -4);
  TCA = ctrans(CA);
  HCA = hermit(CA);
#if _print
  print_CMAT(CA);
  print_CMAT(TCA);
  print_CMAT(HCA);
#endif
  free_CMAT(CA);
  free_CMAT(TCA);
  free_CMAT(HCA);
  // 2D
  CA = alloc_CMAT(2, 4);
  cfill(CA, 2, -4);
  TCA = ctrans(CA);
  HCA = hermit(CA);
#if _print
  print_CMAT(CA);
  print_CMAT(TCA);
  print_CMAT(HCA);
#endif
  free_CMAT(CA);
  free_CMAT(TCA);
  free_CMAT(HCA);
  // 3D as batched 2D
  CA = alloc_CMAT(2, 3, 4);
  cfill(CA, 2, -3);
  TCA = ctrans(CA);
  HCA = hermit(CA);
#if _print
  print_CMAT(CA);
  print_CMAT(TCA);
  print_CMAT(HCA);
#endif
  free_CMAT(CA);
  free_CMAT(TCA);
  free_CMAT(HCA);

#endif

  return 0;
}
