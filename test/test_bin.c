#include "mother.h"

int main() {
  CMAT *CA;
  CMAT *CB;

  CA = czeros(4, 3);
  CB = czeros(4, 3);
  read_cmat("complex.bin", CA);
  print_cmat(CA);

  write_cmat("sec.bin", CA);
  read_cmat("sec.bin", CB);
  print_cmat(CB);

  free_cmat(CB);

  free_cmat(CA);

  return 0;
}
