#include "mother.h"
//#include "header_for_test.h"

// real
#define _A 0  // N
#define _B 0  // T
// complex
#define _C 0  // N
#define _D 1  // T
#define _E 1  // H

int main() {
  MAT *A, *X, *Y;
  MAT *TA, *TX, *TY;

  CMAT *CA, *CX, *CY;
  CMAT *TCA, *TCX, *TCY;
  CTYPE calpha, cbeta;
  calpha.re = 5;
  calpha.im = 0;
  cbeta.re = 0;
  cbeta.im = 0;
#if USE_CUDA
  init();
#endif

//	stopwatch(0);

#if _A
  A = zeros(5, 4);
  X = alloc_mat(5);
  Y = alloc_mat(4);

  fill(A, 2);
  fill(X, 1);
  fill(Y, -1);

  print_mat(A);
  print_mat(X);
  print_mat(Y);

  gemv_mat(NoTran, 10, A, X, 3, Y);

  print_mat(A);
  print_mat(X);
  print_mat(Y);

  free_mat(A);
  free_mat(X);
  free_mat(Y);

#endif
/**** TRANSPOSE  ****/
#if _B
  TA = zeros(4, 5);
  TX = alloc_mat(5);
  TY = alloc_mat(4);

  fill(TA, 2);
  fill(TX, 1);
  fill(TY, -1);

  print_mat(TA);
  print_mat(TX);
  print_mat(TY);

  gemv_mat(Tran, 10, TA, TX, 3, TY);

  print_mat(TA);
  print_mat(TX);
  print_mat(TY);

  free_mat(TA);
  free_mat(TX);
  free_mat(TY);
#endif
/**** COMPLEX ****/
#if _C

  CA = czeros(4, 3);
  CX = alloc_cmat(4);
  CY = alloc_cmat(3);

  cfill(CA, 3, 0);
  cfill(CX, -2, 0);
  cfill(CY, 4, -2);

  print_cmat(CA);
  print_cmat(CX);
  print_cmat(CY);

  gemv_cmat(NoTran, calpha, CA, CX, cbeta, CY);

  print_cmat(CA);
  print_cmat(CX);
  print_cmat(CY);

  free_cmat(CA);
  free_cmat(CX);
  free_cmat(CY);
#endif
/**** TRANSPOSE  ****/
#if _D
  TCA = czeros(3, 4);
  TCX = alloc_cmat(4);
  TCY = alloc_cmat(3);

  cfill(TCA, 0, 1);
  cfill(TCX, 1, 0);
  cfill(TCY, 4, -2);

  print_cmat(TCA);
  print_cmat(TCX);
  print_cmat(TCY);

  gemv_cmat(Tran, calpha, TCA, TCX, cbeta, TCY);

  print_cmat(TCA);
  print_cmat(TCX);
  print_cmat(TCY);

  free_cmat(TCA);
  free_cmat(TCX);
  free_cmat(TCY);
#endif
  /**** Conjugate transpose ****/
  TCA = czeros(3, 4);
  TCX = alloc_cmat(4);
  TCY = alloc_cmat(3);

  cfill(TCA, 0, 1.0);
  cfill(TCX, 1, 0);
  cfill(TCY, 4, -2);

  print_cmat(TCA);
  print_cmat(TCX);
  print_cmat(TCY);

  gemv_cmat(CTran, calpha, TCA, TCX, cbeta, TCY);

  print_cmat(TCA);
  print_cmat(TCX);
  print_cmat(TCY);

  free_cmat(TCA);
  free_cmat(TCX);
  free_cmat(TCY);
#if USE_CUDA
  finit();
#endif
  return 0;
}
