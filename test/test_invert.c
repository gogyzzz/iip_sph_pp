#include "mother.h"

#define _print 0
int main() {
  CMAT *CA, *CB, *CC;
  ITER i;

  CA = alloc_cmat(8, 8);
  CB = alloc_cmat(8, 8);
  CC = alloc_cmat(8, 8);
  read_cmat("c8by8.bin", CA);
  read_cmat("c8by8_inv.bin", CC);
  cinvert(CA, CB);

  if (_print) {
    print_cmat(CA);
    print_cmat(CB);
    print_cmat(CC);
  }

  if (compare_cmat(CB, CC) == 1)
    printf("OK\n");
  else
    printf("SHIT\n");

  free_cmat(CA);
  free_cmat(CB);
  free_cmat(CC);
  return 0;
}
