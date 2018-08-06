#include "mother.h"

int main() {
  MAT *A, *B;
  CMAT *CA, *CB;
  CTYPE ctemp = {5, -7};
  int i, j;
  init(0);

  printf("===== init =====\n");

  A = zeros(2, 3);
  CA = czeros(2, 3);

  B = alloc_mat(2, 3);
  CB = alloc_cmat(2, 3);

  printf("===== A,B ,CA,CB  =====\n");

  //	fill(B,10);
  //	cfill(CB,10,10);
  for (i = 0; i < A->d0 * A->d1 * A->d2; i++) A->data[i] = i;
  for (i = 0; i < B->d0 * B->d1 * B->d2; i++) B->data[i] = -i;
  for (i = 0; i < CA->d0 * CA->d1 * CA->d2; i++) {
    CA->data[i].re = i;
    CA->data[i].im = -i;
  }
  for (i = 0; i < CB->d0 * CB->d1 * CB->d2; i++) {
    CB->data[i].re = -i;
    CB->data[i].im = i;
  }

  print_mat(A);
  print_mat(B);
  print_cmat(CA);
  print_cmat(CB);
  printf("===== axpy(5,B,A) =====\n");

  axpy(5, B, A);
  caxpy(ctemp, CB, CA);

  printf("===== print A,CA =====\n");
  print_mat(A);
  print_cmat(CA);

  printf("===== copy(A,B) =====\n");

  copy(A, B);
  ccopy(CA, CB);

  printf("===== print B,CB =====\n");
  print_mat(B);
  print_cmat(CB);

  for (i = 0; i < A->d0 * A->d1 * A->d2; i++) A->data[i] = i;
  for (i = 0; i < B->d0 * B->d1 * B->d2; i++) B->data[i] = -i;
  for (i = 0; i < CA->d0 * CA->d1 * CA->d2; i++) {
    CA->data[i].re = i;
    CA->data[i].im = -i;
  }
  for (i = 0; i < CB->d0 * CB->d1 * CB->d2; i++) {
    CB->data[i].re = -i;
    CB->data[i].im = i;
  }

  printf("===== print A,B ,CA,CB =====\n");
  print_mat(A);
  print_mat(B);
  print_cmat(CA);
  print_cmat(CB);

  printf("===== asum(B,1) =====\n");
  printf("asum(B,1) = %.2f\n", asum(B, 1));
  printf("casum(CB,1) = %.2f\n", casum(CB, 1));

  /*
          printf("===== fill(AB, 2) =====\n");
          fill(A, 2);
          fill(B, 4);
          cfill(CA, 2, 1);
          cfill(CB, 4, 2);
          print_mat(A);
          print_mat(B);
  */
  printf("===== dot(A) =====\n");
  printf("dot(A, 1, B, 1) = %.2f\n\n", dot(A, 1, B, 1));

  printf("===== cdot(CA, A) =====\n");
  CTYPE temp = cdot(CA, 1, A, 1);
  printf("cdot(CA, 1, A, 1).re = %.2f\n", temp.re);
  printf("cdot(CA, 1, A, 1).im = %.2f\n\n", temp.im);

  printf("===== udot(CA, CB) =====\n");
  temp = udot(CA, 1, CB, 1);
  printf("udot(CA, 1, CB, 1).re = %.2f\n", temp.re);
  printf("udot(CA, 1, CB, 1).im = %.2f\n\n", temp.im);

  printf("===== swap(A, B) =====\n");
  swap(A, B);
  printf("#A\n");
  print_mat(A);
  printf("#B\n");
  print_mat(B);

  printf("===== cswap(CA, CB) =====\n");
  cswap(CA, CB);
  printf("#CA\n");
  print_cmat(CA);
  printf("#CB\n");
  print_cmat(CB);

  printf("===== amin(A) =====\n");
  print_mat(A);
  printf("A[%u] : %lf\n", amin(A), A->data[amin(A)]);

  printf("===== camin(CA) =====\n");
  print_mat(A);
  printf("CA[%u] : %lf %lf\n", camin(CA), CA->data[camin(CA)].re,
         CA->data[camin(CA)].im);

  printf("===== amax(A) =====\n");
  print_mat(A);
  printf("A[%u] : %lf\n", amax(A), A->data[amax(A)]);

  printf("===== camax(CA) =====\n");
  print_cmat(CA);
  printf("CA[%u] : %lf %lf\n", camax(CA), CA->data[camax(CA)].re,
         CA->data[camax(CA)].im);

  printf("===== cabs1(%lf,%lfi) =====\n", ctemp.re, ctemp.im);
  printf("cabs1(ctemp) : %lf\n", cabs1(ctemp));

  printf("===== nrm2(A) =====\n");
  printf("%lf\n", nrm2(A));

  printf("===== cnrm2(CA) =====\n");
  printf("%lf\n", cnrm2(CA));

  printf("==== print(A) scal(100.1,A) print(A) ====\n");
  print_mat(A);
  scal(100.1, A);
  print_mat(A);

  printf("==== print(CA) cscal(100.1,CA) print(CA) ====\n");
  print_cmat(CA);
  cscal(100.1, CA);
  print_cmat(CA);

  printf("==== print(CA) uscal(-10 + 10i  ,CA) print(CA) ====\n");
  print_cmat(CA);
  ctemp.re = -10;
  ctemp.im = 10;
  uscal(ctemp, CA);
  print_cmat(CA);

  free_mat(A);
  free_mat(B);

  free_cmat(CA);
  free_cmat(CB);

  finit();

  return 0;
}
