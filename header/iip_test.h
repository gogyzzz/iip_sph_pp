#ifndef IIP_TEST_H
#define IIP_TEST_H

#include "iip_matrix.h"
#include "iip_type.h"

UINT _eqdd(DTYPE A, DTYPE B);
UINT compare_mat(MAT* A, MAT* B);
UINT compare_cmat(CMAT* A, CMAT* B);

void perform_test();
void do_test(char* filename);
void append_post(char* filename, const char* post, char* out);

#endif
