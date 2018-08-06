#include "mother.h"

#define _print 0

#define m 400
#define n 500
#define k 600

int main() {
  MAT *A, *B, *C;
  ITER i, j;
  init(16776960);

  for (j = 0; j < 10; j++) {
    A = mpalloc(sizeof(MAT));
    B = mpalloc(sizeof(MAT));
    C = mpalloc(sizeof(MAT));

    A->data = mpalloc(sizeof(DTYPE) * m * n);
    B->data = mpalloc(sizeof(DTYPE) * n * k);
    C->data = mpalloc(sizeof(DTYPE) * m * k);

    A->ndim = 1;
    B->ndim = 1;
    C->ndim = 1;

    A->d2 = 1;
    B->d2 = 1;
    C->d2 = 1;

    A->d1 = n;
    A->d0 = m;
    B->d1 = k;
    B->d0 = n;
    C->d1 = k;
    C->d0 = m;

    fill(A, 1);
    fill(B, 2);
    fill(C, 0);
#if _print
    print_mat(A);
    print_mat(B);
    print_mat(C);
#endif
    printf("[%d X %d] = [%d X %d] * [%d X %d]\n", m, k, m, n, n, k);
    matmul(A, B, C);
#if _print
    print_mat(A);
    print_mat(B);
    print_mat(C);
#endif
    mpfree(A->data);
    mpfree(A);
    mpfree(C->data);
    mpfree(B->data);
    mpfree(C);
    mpfree(B);

    A = mpalloc(sizeof(MAT));
    A->data = mpalloc(sizeof(DTYPE) * 20 * 100);
    B = mpalloc(sizeof(MAT));
    C = mpalloc(sizeof(MAT));
    C->data = mpalloc(sizeof(DTYPE) * 24 * 1);
    B->data = mpalloc(sizeof(DTYPE) * 10 * 133);

    mpfree(B->data);
    mpfree(A->data);
    mpfree(A);
    mpfree(B);
    mpfree(C->data);
    mpfree(C);
    //	A = mpalloc_mat(10);

    //	free_mpalloc_mat(10);
  }
  finit();
  return 0;
}
