/*
 * ===========================================================
 *           Copyright (c) 2018, __IIPLAB__
 *                All rights reserved.
 *
 * This Source Code Form is subject to the terms of
 * the Mozilla Public License, v. 2.0.
 * If a copy of the MPL was not distributed with this file,
 *  You can obtain one at http://mozilla.org/MPL/2.0/.
 * ===========================================================
 */
#ifndef IIP_TEST_H
#define IIP_TEST_H

#include "iip_matrix.h"
#include "iip_type.h"

typedef struct FUNCS {
  char name[32];
  int param_cnt;
  void (*fp)();
} FUNCS;

FUNCS func_list[64];

/*check equal of two DTYPEs*/
UINT _eqdd(DTYPE A, DTYPE B);
/*compare 2 MAT struct*/
UINT compare_mat(MAT *A, MAT *B);
/*compare 2 CMAT struct*/
UINT compare_cmat(CMAT *A, CMAT *B);

void perform_test();
void init_list();
void preheat();

void test_verification(int heat, int print_flag, int compare);
// void test_performance(int heat, int print_flag);      //legacy
void test_viewlist();

/*char* operation : out = filename_post*/
void append_post(char *filename, const char *post, char *out);
void print_compare(char *func_name, int is_ctype, void *mat1, void *mat2);
void read_ans(int is_ctype, char *func_name, char *ans_name, void *data_d,
              void *data_c, void **result_mat);
void read_broad_ans(int is_ctype, int param_type, char *func_name,
                    char *ans_name, void *data_d, void *data_c,
                    void **result_mat);

#define PRGB_SIZE 30
void progressbar(int size);
void progress_update(DTYPE rate, int size);

#endif
