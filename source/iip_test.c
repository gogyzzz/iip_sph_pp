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
#include "iip_io.h"
#include "iip_test.h"
#include "iip_math.h"
#include "iip_matrix.h"
#include "iip_invert.h"
#include "iip_blas_lv1.h"
#include "iip_blas_lv3.h"
#include "iip_time.h"

static int is_init = 0;
const int func_list_size = 51;

UINT _eq(DTYPE A, DTYPE B) {
#if NTYPE == 0
  return fabsf((A + 1.0) - B) < (FZERO + 1.0) ? 1 : 0;
#elif NTYPE == 1
  return fabs((A + 1.0) - B) < (FZERO + 1.0) ? 1 : 0;
#endif
  return 0;
}

UINT compare_mat(MAT *A, MAT *B) {
  ITER i;
  for (i = 0; i < A->d0 * A->d1 * A->d2; i++)
    if (!_eq(A->data[i], B->data[i])) {
      printf("%.16lf != %.16lf \n",A->data[i],B->data[i]);
      return 0;}
  return 1;
}

UINT compare_cmat(CMAT *A, CMAT *B) {
  ITER i;
  for (i = 0; i < A->d0 * A->d1 * A->d2; i++) {
    if (!_eq(A->data[i].re, B->data[i].re)){ 
      printf("re %.16lf != %.16lf \n",A->data[i].re,B->data[i].re);
      return 0;}
    if (!_eq(A->data[i].im, B->data[i].im)){
      printf("im %.16lf != %.16lf \n",A->data[i].im,B->data[i].im);
      return 0;}
  }
  return 1;
}

void append_post(char *filename, const char *post, char *out) {
  char t[MAX_CHAR] = "../test_ans/";

  strcpy(out, filename);
  strtok(out, ".");
  strcat(t, out);
  strcat(t, post);
  strcat(t, ".bin");
  strcpy(out, t);
}


void test_verification(int heat, int print_flag, int compare){
  if(!is_init)
    init_list();
  printf("Verifying...\n");

/////////////////////////////////////////////

  int is_ctype = 0;     // 0 for DTYPE, 1 for CTYPE
  int i = 0, j = 0, calced = 0;
  long long total;
  char ans_path[256];

  /*** Input Data ***/
  MAT *A, *B, *C;
  MAT *A_, *B_, *C_;
  MAT *A_diagonal, *A_trace, *A_mul;
  MAT *B_mul;
  MAT *C_mul;
  MAT *A_broad[8] = {NULL, };
  MAT *C_broad[4] = {NULL, };

  CMAT *cA, *cB, *cC;
  CMAT *cA_, *cB_, *cC_;
  CMAT *cA_diagonal, *cA_trace, *cA_mul;
  CMAT *cB_mul;
  CMAT *cC_mul;
  CMAT *cA_broad[8] = {NULL, };
  CMAT *cC_broad[4] = {NULL, };
  
  void *pA_, *pB_, *pC_, *pA_diagonal, *pA_trace, *pA_mul, *pB_mul, *pC_mul, *result;
  void *pA_broad[8] = {NULL, };
  void *pC_broad[4] = {NULL, };

  /*** Answer Data ***/
  MAT *ans_A = NULL, *ans_mul = NULL, *ans_diagonal = NULL, *ans_trace = NULL;
  CMAT *ans_cA = NULL, *ans_cmul = NULL, *ans_cdiagonal = NULL, *ans_ctrace = NULL;

  MAT *ans_broad[4] = {NULL, };
  CMAT *ans_cbroad[4] = {NULL, };
  void *pAns = NULL;
  void *pAns_broad[4] = {NULL, };

  int mat_size = 4, mat_batch = 2;
  int broad_alpha = 1;
  int broad_beta = 2;
  int broad_gamma = 4;

  CTYPE temp_param = {0, 0}, temp_param2 = {0, 0};

  if(!is_init)
    init_list();
  if(print_flag)
    printf(" Measuring...\n");

  A = zeros(mat_size, mat_size, mat_batch);
  A_ = zeros(mat_size, mat_size, mat_batch);
  B = zeros(mat_size, mat_size, mat_batch);
  B_ = zeros(mat_size, mat_size, mat_batch);
  C = zeros(mat_size, mat_size, mat_batch);
  C_ = zeros(mat_size, mat_size, mat_batch);
  ans_A = zeros(mat_size, mat_size, mat_batch);

  cA = czeros(mat_size, mat_size, mat_batch);
  cA_ = czeros(mat_size, mat_size, mat_batch);
  cB = czeros(mat_size, mat_size, mat_batch);
  cB_ = czeros(mat_size, mat_size, mat_batch);
  cC = czeros(mat_size, mat_size, mat_batch);
  cC_ = czeros(mat_size, mat_size, mat_batch);
  ans_cA = czeros(mat_size, mat_size, mat_batch);

  A_broad[0] = zeros(broad_alpha, broad_alpha, broad_alpha);
  A_broad[1] = zeros(broad_alpha, broad_alpha, broad_beta);
  A_broad[2] = zeros(broad_alpha, broad_gamma, broad_alpha);
  A_broad[3] = zeros(broad_alpha, broad_gamma, broad_beta);
  A_broad[4] = zeros(broad_gamma, broad_alpha, broad_alpha);
  A_broad[5] = zeros(broad_gamma, broad_alpha, broad_beta);
  A_broad[6] = zeros(broad_gamma, broad_gamma, broad_alpha);
  A_broad[7] = zeros(broad_gamma, broad_gamma, broad_beta);

  for(i=0; i<8; i++){
    sprintf(ans_path, "../test_data/d_%d_%d_%d.bin", A_broad[i]->d0, A_broad[i]->d1, A_broad[i]->d2);
    read_mat(ans_path, A_broad[i]);
  }

  cA_broad[0] = czeros(broad_alpha, broad_alpha, broad_alpha);
  cA_broad[1] = czeros(broad_alpha, broad_alpha, broad_beta);
  cA_broad[2] = czeros(broad_alpha, broad_gamma, broad_alpha);
  cA_broad[3] = czeros(broad_alpha, broad_gamma, broad_beta);
  cA_broad[4] = czeros(broad_gamma, broad_alpha, broad_alpha);
  cA_broad[5] = czeros(broad_gamma, broad_alpha, broad_beta);
  cA_broad[6] = czeros(broad_gamma, broad_gamma, broad_alpha);
  cA_broad[7] = czeros(broad_gamma, broad_gamma, broad_beta);

  for(i=0; i<8; i++){
    sprintf(ans_path, "../test_data/c_%d_%d_%d.bin", cA_broad[i]->d0, cA_broad[i]->d1, cA_broad[i]->d2);
    read_cmat(ans_path, cA_broad[i]);
  }

  
  C_broad[0] = zeros(broad_alpha, broad_gamma, broad_beta);
  C_broad[1] = zeros(broad_gamma, broad_alpha, broad_beta);
  C_broad[2] = zeros(broad_gamma, broad_gamma, broad_beta);
  C_broad[3] = zeros(broad_gamma, broad_gamma, broad_alpha);

  ans_broad[0] = zeros(broad_alpha, broad_gamma, broad_beta);
  ans_broad[1] = zeros(broad_gamma, broad_alpha, broad_beta);
  ans_broad[2] = zeros(broad_gamma, broad_gamma, broad_beta);
  ans_broad[3] = zeros(broad_gamma, broad_gamma, broad_alpha);

  cC_broad[0] = czeros(broad_alpha, broad_gamma, broad_beta);
  cC_broad[1] = czeros(broad_gamma, broad_alpha, broad_beta);
  cC_broad[2] = czeros(broad_gamma, broad_gamma, broad_beta);
  cC_broad[3] = czeros(broad_gamma, broad_gamma, broad_alpha);

  ans_cbroad[0] = czeros(broad_alpha, broad_gamma, broad_beta);
  ans_cbroad[1] = czeros(broad_gamma, broad_alpha, broad_beta);
  ans_cbroad[2] = czeros(broad_gamma, broad_gamma, broad_beta);
  ans_cbroad[3] = czeros(broad_gamma, broad_gamma, broad_alpha);

  //set A
  read_mat("../test_data/d_4_4_2.bin", A);
  copy(A, A_);
  read_cmat("../test_data/c_4_4_2.bin", cA);
  ccopy(cA, cA_);
  //set B
  read_mat("../test_data/d_4_4_2.bin", B);
  copy(B, B_);
  read_cmat("../test_data/c_4_4_2.bin", cB);
  ccopy(cB, cB_);
  //set C
  copy(A, C);
  copy(C, C_);
  ccopy(cA, cC);
  ccopy(cC, cC_);
  //set A_*
  A_diagonal = zeros(mat_size, 1, mat_batch);
  A_trace = zeros(1, 1, mat_batch);
  ans_diagonal = zeros(mat_size, 1, mat_batch);
  ans_trace = zeros(1, 1, mat_batch);
  read_mat("../test_ans/d_4_4_2_Diagonal.bin", ans_diagonal);
  read_mat("../test_ans/d_4_4_2_Trace.bin", ans_trace);

  cA_diagonal = czeros(mat_size, 1, mat_batch);
  cA_trace = czeros(1, 1, mat_batch);
  ans_cdiagonal = czeros(mat_size, 1, mat_batch);
  ans_ctrace = czeros(1, 1, mat_batch);
  read_cmat("../test_ans/c_4_4_2_cDiagonal.bin", ans_cdiagonal);
  read_cmat("../test_ans/c_4_4_2_cTrace.bin", ans_ctrace);

  A_mul = zeros(1, 4, 2);
  B_mul = zeros(4, 4, 2);
  C_mul = zeros(1, 4, 2);
  read_mat("../test_data/d_1_4_2.bin", A_mul);
  read_mat("../test_data/d_4_4_2.bin", B_mul);
  ans_mul = zeros(1, 4, 2);
  read_mat("../test_ans/d_1_4_2_Gemm.bin", ans_mul);
  
  cA_mul = czeros(1, 4, 2);
  cB_mul = czeros(4, 4, 2);
  cC_mul = czeros(1, 4, 2);
  read_cmat("../test_data/c_1_4_2.bin", cA_mul);
  read_cmat("../test_data/c_4_4_2.bin", cB_mul);
  ans_cmul = czeros(1, 4, 2);
  read_cmat("../test_ans/c_1_4_2_cGemm.bin", ans_cmul);


  is_ctype = 0;
  total = 0;

  if(heat == 1){
    preheat();
    test_performance(0, 0);
  }

  if(print_flag){
    printf(" *** Input Data ***\n");
    print_mat(A);
  }

  for(i = 0; i<func_list_size; i++){
    is_ctype = i%2;

    // set data
    copy(A, A_);
    copy(B, B_);
    copy(C, C_);
    ccopy(cA, cA_);
    ccopy(cB, cB_);
    ccopy(cC, cC_);

    if (is_ctype == 0) {
      pA_ = A_;
      pB_ = B_;
      pC_ = C_;
      pA_diagonal = A_diagonal;
      pA_trace = A_trace;
      pA_mul = A_mul;
      pB_mul = B_mul;
      pC_mul = C_mul;
      for(j=0; j<8; j++){
        pA_broad[j] = A_broad[j];
      }
      for(j=0; j<4; j++){
        pC_broad[j] = C_broad[j];
        pAns_broad[j] = ans_broad[j];
      }
    } else {
      pA_ = cA_;
      pB_ = cB_;
      pC_ = cC_;
      pA_diagonal = cA_diagonal;
      pA_trace = cA_trace;
      pA_mul = cA_mul;
      pB_mul = cB_mul;
      pC_mul = cC_mul;
      for(j=0; j<8; j++){
        pA_broad[j] = cA_broad[j];
      }
      for(j=0; j<4; j++){
        pC_broad[j] = cC_broad[j];
        pAns_broad[j] = ans_cbroad[j];
      }
    }
    calced = 0;
    pAns = NULL;

    /* call appropriate func */

    // Require ONE matrix func
    if (func_list[i].param_cnt == 1) {
      calced = 1;
      stopwatch(0);
      func_list[i].fp(pA_);
      total += stopwatch(1);
      result = pA_;

      //read_ans(is_ctype, func_list[i].name, ans_path, ans_A, ans_cA, pAns);
    } else if (i == 24 || i == 25) {  // permute
      calced = 1;
      stopwatch(0);
      func_list[i].fp(pA_, 213);
      total += stopwatch(1);
      result = pA_;

      //read_ans(is_ctype, func_list[i].name, ans_path, ans_A, ans_cA, pAns);
    } else if (i == 34 || i == 35) {  // sum
      calced = 1;
      stopwatch(0);
      func_list[i].fp(pA_, 1);
      total += stopwatch(1);
      result = pA_;
	  printf("\n # skip asum\n");
	  continue;
    } else if (i == func_list_size - 9) {  // pow_mat
      calced = 1;
      pA_ = A_;
      stopwatch(0);
      func_list[i].fp(pA_, 2.0);
      total += stopwatch(1);
      result = pA_;
      is_ctype = 0;
    } else if (i == func_list_size - 8) {  // pow_cmat
      calced = 1;
      pA_ = cA_;
      stopwatch(0);
      func_list[i].fp(pA_, 2.0);
      total += stopwatch(1);
      result = pA_;
      is_ctype = 1;
    } else if (i == func_list_size - 7) {  // cpow_cmat
      calced = 1;
      pA_ = cA_;
      temp_param.re = 2.0;
      temp_param.im = 2.0;
      stopwatch(0);
      func_list[i].fp(pA_, temp_param);
      total += stopwatch(1);
      result = pA_;
      is_ctype = 1;
    } else if (i == func_list_size - 6) {  // add
      calced = 1;
      pA_ = A_;
      stopwatch(0);
      func_list[i].fp(10.0, pA_);
      total += stopwatch(1);
      result = pA_;
      is_ctype = 0;
    } else if (i == func_list_size - 5) {  // cadd
      calced = 1;
      pA_ = cA_;
      stopwatch(0);
      func_list[i].fp(10.0, pA_);
      total += stopwatch(1);
      result = pA_;
      is_ctype = 1;
    } else if (i == func_list_size - 4) {  // uadd
      calced = 1;
      pA_ = cA_;
      temp_param.re = 10.0;
      temp_param.im = 10.0;
      stopwatch(0);
      func_list[i].fp(temp_param, pA_);
      total += stopwatch(1);
      result = pA_;
      is_ctype = 1;
    } else if (i == func_list_size - 3) {  // scal
      calced = 1;
      pA_ = A_;
      stopwatch(0);
      func_list[i].fp(2.0, pA_);
      total += stopwatch(1);
      result = pA_;
      is_ctype = 0;
    } else if (i == func_list_size - 2) {  // cscal
      calced = 1;
      pA_ = cA_;
      stopwatch(0);
      func_list[i].fp(2.0, pA_);
      total += stopwatch(1);
      result = pA_;
      is_ctype = 1;
    } else if (i == func_list_size - 1) {  // uscal
      calced = 1;
      pA_ = cA_;
      temp_param.re = 2.0;
      temp_param.im = 2.0;
      stopwatch(0);
      func_list[i].fp(temp_param, pA_);
      total += stopwatch(1);
      result = pA_;
      is_ctype = 1;
    }
    //Require TWO matrix func (one input, one write buffer)
    else if(func_list[i].param_cnt == 2 && !(i == 22 || i == 23 || i == 26 || i == 27)){
      calced = 1;
      if(i == 30 || i == 31){         // diagonal
        stopwatch(0);
        func_list[i].fp(pA_, pA_diagonal);
        total += stopwatch(1);
        result = pA_diagonal;
        read_ans(is_ctype, func_list[i].name, ans_path, ans_diagonal, ans_cdiagonal, &pAns);
      }
      else if(i == 32 || i == 33){    // trace
        stopwatch(0);
        func_list[i].fp(pA_, pA_trace);
        total += stopwatch(1);
        result = pA_trace;
        read_ans(is_ctype, func_list[i].name, ans_path, ans_trace, ans_ctrace, &pAns);
      }
      else{                           // invert
        stopwatch(0);
        func_list[i].fp(pA_, pC_);
        total += stopwatch(1);
        result = pC_;
      }
    }
    //Require more than THREE matrix func (two input, one write buffer) (add, mul, div)_elements
    else if(func_list[i].param_cnt == 3){
      calced = 1;
      stopwatch(0);
      func_list[i].fp(pA_, pB_, pC_);
      total += stopwatch(1);
      result = pC_;

      switch(i){
        case 36:
        case 37:
          // Add_elements
          stopwatch(0);
          func_list[i].fp(pA_broad[1], pA_broad[3], pC_broad[0]);
          func_list[i].fp(pA_broad[1], pA_broad[5], pC_broad[1]);
          total += stopwatch(1);

          read_broad_ans(is_ctype, 0, func_list[i].name, ans_path, ans_broad[0], ans_cbroad[0], &pAns_broad[0]);
          read_broad_ans(is_ctype, 1, func_list[i].name, ans_path, ans_broad[1], ans_cbroad[1], &pAns_broad[1]);

          if(compare != 0 && print_flag != 0){
            print_compare(func_list[i].name, is_ctype, pC_broad[0], pAns_broad[0]);
            print_compare(func_list[i].name, is_ctype, pC_broad[1], pAns_broad[1]);
		  }
		  printf("=============================\n");
		  print_mat(pA_broad[1]);
		  printf("------------\n");
		  print_mat(pA_broad[5]);
		  printf("------------\n");
		  print_mat(pC_broad[1]);
		  printf("------------\n");
		  print_mat(pAns_broad[1]);
		  printf("------------\n");
		  printf("=============================\n");

          stopwatch(0);
          func_list[i].fp(pA_broad[3], pA_broad[1], pC_broad[0]);
          func_list[i].fp(pA_broad[5], pA_broad[0], pC_broad[1]);
          total += stopwatch(1);

          read_broad_ans(is_ctype, 2, func_list[i].name, ans_path, ans_broad[0], ans_cbroad[0], &pAns_broad[0]);
          read_broad_ans(is_ctype, 3, func_list[i].name, ans_path, ans_broad[1], ans_cbroad[1], &pAns_broad[1]);

          if(compare != 0 && print_flag != 0){
            print_compare(func_list[i].name, is_ctype, pC_broad[0], pAns_broad[0]);
            print_compare(func_list[i].name, is_ctype, pC_broad[1], pAns_broad[1]);
          }
          break;
        case 38:
        case 39:
          // Mul_elements
          //printf(" PASS MUL_ELEMENTS\n");
          stopwatch(0);
          func_list[i].fp(pA_broad[1], pA_broad[6], pC_broad[2]);
          total += stopwatch(1);
          read_broad_ans(is_ctype, 4, func_list[i].name, ans_path, ans_broad[2], ans_cbroad[2], &pAns_broad[2]);

          if(compare != 0 && print_flag != 0)
            print_compare(func_list[i].name, is_ctype, pC_broad[2], pAns_broad[2]);

          stopwatch(0);
          func_list[i].fp(pA_broad[7], pA_broad[0], pC_broad[2]);
          total += stopwatch(1);
          read_broad_ans(is_ctype, 5, func_list[i].name, ans_path, ans_broad[2], ans_cbroad[2], &pAns_broad[2]);

          if(compare != 0 && print_flag != 0)
            print_compare(func_list[i].name, is_ctype, pC_broad[2], pAns_broad[2]);

          stopwatch(0);
          func_list[i].fp(pA_broad[2], pA_broad[5], pC_broad[2]);
          total += stopwatch(1);
          read_broad_ans(is_ctype, 6, func_list[i].name, ans_path, ans_broad[2], ans_cbroad[2], &pAns_broad[2]);

          if(compare != 0 && print_flag != 0)
            print_compare(func_list[i].name, is_ctype, pC_broad[2], pAns_broad[2]);
          
          stopwatch(0);
          func_list[i].fp(pA_broad[4], pA_broad[3], pC_broad[2]);
          total += stopwatch(1);
          read_broad_ans(is_ctype, 7, func_list[i].name, ans_path, ans_broad[2], ans_cbroad[2], &pAns_broad[2]);

          if(compare != 0 && print_flag != 0)
            print_compare(func_list[i].name, is_ctype, pC_broad[2], pAns_broad[2]);

          stopwatch(0);
          func_list[i].fp(pA_broad[2], pA_broad[7], pC_broad[2]);
          total += stopwatch(1);
          read_broad_ans(is_ctype, 8, func_list[i].name, ans_path, ans_broad[2], ans_cbroad[2], &pAns_broad[2]);

          if(compare != 0 && print_flag != 0)
            print_compare(func_list[i].name, is_ctype, pC_broad[2], pAns_broad[2]);
          break;
        case 40:
        case 41:
          // Div_elements
          //printf(" PASS DIV_ELEMENTS\n");
          stopwatch(0);
          func_list[i].fp(pA_broad[4], pA_broad[6], pC_broad[3]);
          total += stopwatch(1);
          read_broad_ans(is_ctype, 9, func_list[i].name, ans_path, ans_broad[3], ans_cbroad[3], &pAns_broad[3]);

          if(compare != 0 && print_flag != 0)
            print_compare(func_list[i].name, is_ctype, pC_broad[3], pAns_broad[3]);

          stopwatch(0);
          func_list[i].fp(pA_broad[6], pA_broad[2], pC_broad[3]);
          total += stopwatch(1);
          read_broad_ans(is_ctype, 10, func_list[i].name, ans_path, ans_broad[3], ans_cbroad[3], &pAns_broad[3]);

          if(compare != 0 && print_flag != 0)
            print_compare(func_list[i].name, is_ctype, pC_broad[3], pAns_broad[3]);

          stopwatch(0);
          func_list[i].fp(pA_broad[6], pA_broad[4], pC_broad[3]);
          total += stopwatch(1);
          read_broad_ans(is_ctype, 11, func_list[i].name, ans_path, ans_broad[3], ans_cbroad[3], &pAns_broad[3]);

          if(compare != 0 && print_flag != 0)
            print_compare(func_list[i].name, is_ctype, pC_broad[3], pAns_broad[3]);
          
          stopwatch(0);
          func_list[i].fp(pA_broad[6], pA_broad[6], pC_broad[3]);
          total += stopwatch(1);
          read_broad_ans(is_ctype, 12, func_list[i].name, ans_path, ans_broad[3], ans_cbroad[3], &pAns_broad[3]);

          if(compare != 0 && print_flag != 0)
            print_compare(func_list[i].name, is_ctype, pC_broad[3], pAns_broad[3]);
          break;
        default:
          // nothing
          break;
      }

      // not implemented yet
      continue;
    }
    // only gemm (matmul)
    else if(func_list[i].param_cnt > 3){
      calced = 1;
      if(is_ctype == 0){
        stopwatch(0);
        func_list[i].fp(NoTran, NoTran, 1.0, pA_mul, pB_mul, 0.0, pC_mul);
        total += stopwatch(1);
        pAns = ans_mul;
      }
      if(is_ctype == 1){
        temp_param.re = 1.0;
        temp_param.im = 0.0;
        temp_param2.re = 0.0;
        temp_param2.im = 0.0;
        stopwatch(0);
        func_list[i].fp(NoTran, NoTran, temp_param, pA_mul, pB_mul, temp_param2, pC_mul);
        total += stopwatch(1);
        pAns = ans_cmul;
      }
      result = pC_mul;
    }
    // not using - repmat, reshape...
    else{
      printf(" passing function [%s]\n", func_list[i].name);
      // do nothing
      continue;
    }

    // print result
    if (calced != 0 && print_flag != 0) {
      if(compare != 0)
      {
        if(pAns == NULL)
          read_ans(is_ctype, func_list[i].name, ans_path, ans_A, ans_cA, &pAns);

        print_compare(func_list[i].name, is_ctype, result, pAns);
        //print_mat(pAns);
        //print_mat(result);
      }
      else{
        printf(" Done.\n");
      }
    }
  }

  if(print_flag)
    printf("\n *** Broad Cast Function ***\n");

  if(print_flag)
    printf("Done.\n *** Total Run time : %.3lfms\n", (double)total/1000.0);

/////////////////////////////////////////////

  /*** Release Memory ***/
  free_mat(A);  free_mat(B);  free_mat(C);
  free_mat(A_); free_mat(B_); free_mat(C_);
  free_mat(A_diagonal); free_mat(A_trace);  free_mat(A_mul);
  free_mat(B_mul);  free_mat(C_mul);
  
  free_cmat(cA);  free_cmat(cB);  free_cmat(cC);
  free_cmat(cA_); free_cmat(cB_); free_cmat(cC_);
  free_cmat(cA_diagonal); free_cmat(cA_trace);  free_cmat(cA_mul);
  free_cmat(cB_mul);  free_cmat(cC_mul);
}

void test_performance(int heat, int print_flag){
  int is_ctype = 0;     // 0 for DTYPE, 1 for CTYPE
  int i = 0, calced = 0;
  long long total;

  if(!is_init)
    init_list();
  if(print_flag)
    printf(" Measuring...\n");

  MAT *A, *B, *C;
  MAT *A_, *B_, *C_;
  MAT *A_diagonal, *A_trace, *A_mul;
  MAT *B_mul;
  MAT *C_mul;

  CMAT *cA, *cB, *cC;
  CMAT *cA_, *cB_, *cC_;
  CMAT *cA_diagonal, *cA_trace, *cA_mul;
  CMAT *cB_mul;
  CMAT *cC_mul;
  
  void *pA_, *pB_, *pC_, *pA_diagonal, *pA_trace, *pA_mul, *pB_mul, *pC_mul, *result;

  int mat_size = 4, mat_batch = 2;

  CTYPE temp_param = {0, 0}, temp_param2 = {0, 0};

  A = zeros(mat_size, mat_size, mat_batch);
  A_ = zeros(mat_size, mat_size, mat_batch);
  B = zeros(mat_size, mat_size, mat_batch);
  B_ = zeros(mat_size, mat_size, mat_batch);
  C = zeros(mat_size, mat_size, mat_batch);
  C_ = zeros(mat_size, mat_size, mat_batch);

  cA = czeros(mat_size, mat_size, mat_batch);
  cA_ = czeros(mat_size, mat_size, mat_batch);
  cB = czeros(mat_size, mat_size, mat_batch);
  cB_ = czeros(mat_size, mat_size, mat_batch);
  cC = czeros(mat_size, mat_size, mat_batch);
  cC_ = czeros(mat_size, mat_size, mat_batch);

  //set A
  read_mat("../test_data/d_4_4_2.bin", A);
  copy(A, A_);
  read_cmat("../test_data/c_4_4_2.bin", cA);
  ccopy(cA, cA_);
  //set B
  read_mat("../test_data/d_4_4_2.bin", B);
  copy(B, B_);
  read_cmat("../test_data/c_4_4_2.bin", cB);
  ccopy(cB, cB_);
  //set C
  copy(A, C);
  copy(C, C_);
  ccopy(cA, cC);
  ccopy(cC, cC_);
  //set A_*
  A_diagonal = zeros(mat_size, 1, mat_batch);
  A_trace = zeros(1, 1, mat_batch);

  cA_diagonal = czeros(mat_size, 1, mat_batch);
  cA_trace = czeros(1, 1, mat_batch);

  A_mul = zeros(2, 4, 1);
  B_mul = zeros(4, 6, 1);
  C_mul = zeros(2, 6, 1);
  read_mat("../test_data/d_2_4_1.bin", A_mul);
  read_mat("../test_data/d_4_6_1.bin", B_mul);
  
  cA_mul = czeros(2, 4, 1);
  cB_mul = czeros(4, 6, 1);
  cC_mul = czeros(2, 6, 1);
  read_cmat("../test_data/c_2_4_1.bin", cA_mul);
  read_cmat("../test_data/c_4_6_1.bin", cB_mul);


  is_ctype = 0;
  total = 0;

  if(heat == 1){
    preheat();
    test_performance(0, 0);
  }

  if(print_flag){
    printf(" *** Input Data ***\n");
    print_mat(A);
  }

  for(i = 0; i<func_list_size; i++){
    is_ctype = i%2;

    // set data
    copy(A, A_);
    copy(B, B_);
    copy(C, C_);
    ccopy(cA, cA_);
    ccopy(cB, cB_);
    ccopy(cC, cC_);

    if (is_ctype == 0) {
      pA_ = A_;
      pB_ = B_;
      pC_ = C_;
      pA_diagonal = A_diagonal;
      pA_trace = A_trace;
      pA_mul = A_mul;
      pB_mul = B_mul;
      pC_mul = C_mul;
    } else {
      pA_ = cA_;
      pB_ = cB_;
      pC_ = cC_;
      pA_diagonal = cA_diagonal;
      pA_trace = cA_trace;
      pA_mul = cA_mul;
      pB_mul = cB_mul;
      pC_mul = cC_mul;
    }
    calced = 0;

    /* call appropriate func */

    // Require ONE matrix func
    if (func_list[i].param_cnt == 1) {
      calced = 1;
      stopwatch(0);
      func_list[i].fp(pA_);
      total += stopwatch(1);
      result = pA_;
    } else if (i == 24 || i == 25) {  // permute
      calced = 1;
      stopwatch(0);
      func_list[i].fp(pA_, 213);
      total += stopwatch(1);
      result = pA_;
    } else if (i == 34 || i == 35) {  // sum
      calced = 1;
      stopwatch(0);
      func_list[i].fp(pA_, 1);
      total += stopwatch(1);
      result = pA_;
    } else if (i == func_list_size - 9) {  // pow_mat
      calced = 1;
      pA_ = A_;
      stopwatch(0);
      func_list[i].fp(pA_, 2.0);
      total += stopwatch(1);
      result = pA_;
    } else if (i == func_list_size - 8) {  // pow_cmat
      calced = 1;
      pA_ = cA_;
      stopwatch(0);
      func_list[i].fp(pA_, 2.0);
      total += stopwatch(1);
      result = pA_;
    } else if (i == func_list_size - 7) {  // cpow_cmat
      calced = 1;
      pA_ = cA_;
      temp_param.re = 2.0;
      temp_param.im = 2.0;
      stopwatch(0);
      func_list[i].fp(pA_, temp_param);
      total += stopwatch(1);
      result = pA_;
    } else if (i == func_list_size - 6) {  // add
      calced = 1;
      pA_ = A_;
      stopwatch(0);
      func_list[i].fp(10.0, pA_);
      total += stopwatch(1);
      result = pA_;
    } else if (i == func_list_size - 5) {  // cadd
      calced = 1;
      pA_ = cA_;
      stopwatch(0);
      func_list[i].fp(10.0, pA_);
      total += stopwatch(1);
      result = pA_;
    } else if (i == func_list_size - 4) {  // uadd
      calced = 1;
      pA_ = cA_;
      temp_param.re = 10.0;
      temp_param.im = 10.0;
      stopwatch(0);
      func_list[i].fp(temp_param, pA_);
      total += stopwatch(1);
      result = pA_;
    } else if (i == func_list_size - 3) {  // scal
      calced = 1;
      pA_ = A_;
      stopwatch(0);
      func_list[i].fp(5.0, pA_);
      total += stopwatch(1);
      result = pA_;
    } else if (i == func_list_size - 2) {  // cscal
      calced = 1;
      pA_ = cA_;
      stopwatch(0);
      func_list[i].fp(5.0, pA_);
      total += stopwatch(1);
      result = pA_;
    } else if (i == func_list_size - 1) {  // uscal
      calced = 1;
      pA_ = cA_;
      temp_param.re = 5.0;
      temp_param.im = 0.0;
      stopwatch(0);
      func_list[i].fp(temp_param, pA_);
      total += stopwatch(1);
      result = pA_;
    }
    //Require TWO matrix func (one input, one write buffer)
    else if(func_list[i].param_cnt == 2 && !(i == 22 || i == 23 || i == 26 || i == 27)){
      calced = 1;
      if(i == 30 || i == 31){         // diagonal
        stopwatch(0);
        func_list[i].fp(pA_, pA_diagonal);
        total += stopwatch(1);
        result = pA_diagonal;
      }
      else if(i == 32 || i == 33){    // trace
        stopwatch(0);
        func_list[i].fp(pA_, pA_trace);
        total += stopwatch(1);
        result = pA_trace;
      }
      else{                           // invert
        stopwatch(0);
        func_list[i].fp(pA_, pC_);
        total += stopwatch(1);
        result = pC_;
      }
    }
    //Require more than THREE matrix func (two input, one write buffer) (add, mul, div)_elements
    else if(func_list[i].param_cnt == 3){
      calced = 1;
      stopwatch(0);
      func_list[i].fp(pA_, pB_, pC_);
      total += stopwatch(1);
      result = pC_;
    }
    // only gemm (matmul)
    else if(func_list[i].param_cnt > 3){
      calced = 1;
      if(is_ctype == 0){
        stopwatch(0);
        func_list[i].fp(NoTran, NoTran, 1.0, pA_mul, pB_mul, 0.0, pC_mul);
        total += stopwatch(1);
      }
      if(is_ctype == 1){
        temp_param.re = 1.0;
        temp_param.im = 0.0;
        temp_param2.re = 0.0;
        temp_param2.im = 0.0;
        stopwatch(0);
        func_list[i].fp(NoTran, NoTran, temp_param, pA_mul, pB_mul, temp_param2, pC_mul);
        total += stopwatch(1);
      }
      result = pC_mul;
    }
    // not using - repmat, reshape...
    else{
      printf(" passing function [%s]\n", func_list[i].name);
      // do nothing
      continue;
    }

    // print result
    if (calced != 0 && print_flag != 0) {
      printf("\n # %s\n", func_list[i].name);
      print_mat(result);
    }
  }

  if(print_flag)
    printf("Done.\n *** Total Run time : %.3lfms\n", (double)total/1000.0);
  
  /*** Release Memory ***/
  free_mat(A);  free_mat(B);  free_mat(C);
  free_mat(A_); free_mat(B_); free_mat(C_);
  free_mat(A_diagonal); free_mat(A_trace);  free_mat(A_mul);
  free_mat(B_mul);  free_mat(C_mul);
  
  free_cmat(cA);  free_cmat(cB);  free_cmat(cC);
  free_cmat(cA_); free_cmat(cB_); free_cmat(cC_);
  free_cmat(cA_diagonal); free_cmat(cA_trace);  free_cmat(cA_mul);
  free_cmat(cB_mul);  free_cmat(cC_mul);
}

void test_viewlist(){
  if(!is_init)
    init_list();
  int i=0;
  printf("\n\n---- Func list ----\n");
  for(i=0; i<func_list_size-9; i++){
    printf("戍成 %s\n", func_list[i++].name);
    printf("弛戌式 %s\n", func_list[i].name);
    printf("弛\n");
  }
  i=func_list_size-9;
  printf("戍成 %s\n", func_list[i++].name);
  printf("弛戍式 %s\n", func_list[i++].name);
  printf("弛戌式 %s\n", func_list[i++].name);
  printf("弛\n");
  printf("戍成 %s\n", func_list[i++].name);
  printf("弛戍式 %s\n", func_list[i++].name);
  printf("弛戌式 %s\n", func_list[i++].name);
  printf("弛\n");
  printf("戌成 %s\n", func_list[i++].name);
  printf(" 戍式 %s\n", func_list[i++].name);
  printf(" 戌式 %s\n", func_list[i].name);

  printf("-------------------\n\n\n");
}

void preheat(){
  int i=0;
  double temp=0;
  long long times=0;

  printf(" CPU is warming up...\n");

  stopwatch(0);
  #pragma omp parallel for shared(temp) private(i)
  for(i=0; i<250000000; i++){
    temp += (double)i / 2.0;
    temp /= (double)i / 3.0;
    temp *= (double)i / 4.0;
  }
  times = stopwatch(1);

  printf(" Took %.3lfms to warm up.\n\n", (double)times / 1000.0);
}

void init_list(){
  is_init = 1;
  func_list[0].fp = invert;
  func_list[0].param_cnt = 2;
  strcpy(func_list[0].name, "Invert");
  func_list[1].fp = cinvert;
  func_list[1].param_cnt = 2;
  strcpy(func_list[1].name, "cInvert");
  
  func_list[2].fp = transpose;
  func_list[2].param_cnt = 1;
  strcpy(func_list[2].name, "Transpose");
  func_list[3].fp = ctranspose;
  func_list[3].param_cnt = 1;
  strcpy(func_list[3].name, "cTranspose");
  
  func_list[4].fp = sqrt_mat;
  func_list[4].param_cnt = 1;
  strcpy(func_list[4].name, "Sqrt");
  func_list[5].fp = sqrt_cmat;
  func_list[5].param_cnt = 1;
  strcpy(func_list[5].name, "cSqrt");
  
  func_list[6].fp = round_mat;
  func_list[6].param_cnt = 1;
  strcpy(func_list[6].name, "Round");
  func_list[7].fp = round_cmat;
  func_list[7].param_cnt = 1;
  strcpy(func_list[7].name, "cRound");
  
  func_list[8].fp = floor_mat;
  func_list[8].param_cnt = 1;
  strcpy(func_list[8].name, "Floor");
  func_list[9].fp = floor_cmat;
  func_list[9].param_cnt = 1;
  strcpy(func_list[9].name, "cFloor");
  
  func_list[10].fp = ceil_mat;
  func_list[10].param_cnt = 1;
  strcpy(func_list[10].name, "Ceil");
  func_list[11].fp = ceil_cmat;
  func_list[11].param_cnt = 1;
  strcpy(func_list[11].name, "cCeil");
  
  func_list[12].fp = log_mat;
  func_list[12].param_cnt = 1;
  strcpy(func_list[12].name, "Log");
  func_list[13].fp = log_cmat;
  func_list[13].param_cnt = 1;
  strcpy(func_list[13].name, "cLog");
  
  func_list[14].fp = log2_mat;
  func_list[14].param_cnt = 1;
  strcpy(func_list[14].name, "Log2");
  func_list[15].fp = log2_cmat;
  func_list[15].param_cnt = 1;
  strcpy(func_list[15].name, "cLog2");
  
  func_list[16].fp = log10_mat;
  func_list[16].param_cnt = 1;
  strcpy(func_list[16].name, "Log10");
  func_list[17].fp = log10_cmat;
  func_list[17].param_cnt = 1;
  strcpy(func_list[17].name, "cLog10");
  
  func_list[18].fp = exp_mat;
  func_list[18].param_cnt = 1;
  strcpy(func_list[18].name, "Exp");
  func_list[19].fp = exp_cmat;
  func_list[19].param_cnt = 1;
  strcpy(func_list[19].name, "cExp");
  
  func_list[20].fp = abs_mat;
  func_list[20].param_cnt = 1;
  strcpy(func_list[20].name, "Abs");
  func_list[21].fp = abs_cmat;
  func_list[21].param_cnt = 1;
  strcpy(func_list[21].name, "cAbs");
  
  func_list[22].fp = repmat;
  func_list[22].param_cnt = 2;
  strcpy(func_list[22].name, "Repmat");
  func_list[23].fp = crepmat;
  func_list[23].param_cnt = 2;
  strcpy(func_list[23].name, "cRepmat");
  
  func_list[24].fp = permute;
  func_list[24].param_cnt = 2;
  strcpy(func_list[24].name, "Permute213");
  func_list[25].fp = cpermute;
  func_list[25].param_cnt = 2;
  strcpy(func_list[25].name, "cPermute213");
  
  func_list[26].fp = reshape;
  func_list[26].param_cnt = 2;
  strcpy(func_list[26].name, "Reshape");
  func_list[27].fp = creshape;
  func_list[27].param_cnt = 2;
  strcpy(func_list[27].name, "cReshape");
  
  func_list[28].fp = gemm;
  func_list[28].param_cnt = 7;
  strcpy(func_list[28].name, "Gemm");
  func_list[29].fp = cgemm;
  func_list[29].param_cnt = 7;
  strcpy(func_list[29].name, "cGemm");
  
  func_list[30].fp = diagonal;
  func_list[30].param_cnt = 2;
  strcpy(func_list[30].name, "Diagonal");
  func_list[31].fp = cdiagonal;
  func_list[31].param_cnt = 2;
  strcpy(func_list[31].name, "cDiagonal");
  
  func_list[32].fp = trace;
  func_list[32].param_cnt = 2;
  strcpy(func_list[32].name, "Trace");
  func_list[33].fp = ctrace;
  func_list[33].param_cnt = 2;
  strcpy(func_list[33].name, "cTrace");
  
  func_list[34].fp = asum;
  func_list[34].param_cnt = 2;
  strcpy(func_list[34].name, "Asum");
  func_list[35].fp = casum;
  func_list[35].param_cnt = 2;
  strcpy(func_list[35].name, "cAsum");
  
  func_list[36].fp = add_elements;
  func_list[36].param_cnt = 3;
  strcpy(func_list[36].name, "Add Elements");
  func_list[37].fp = cadd_elements;
  func_list[37].param_cnt = 3;
  strcpy(func_list[37].name, "cAdd Elements");
  
  func_list[38].fp = mul_elements;
  func_list[38].param_cnt = 3;
  strcpy(func_list[38].name, "Mul Elements");
  func_list[39].fp = cmul_elements;
  func_list[39].param_cnt = 3;
  strcpy(func_list[39].name, "cMul Elements");
  
  func_list[40].fp = div_elements;
  func_list[40].param_cnt = 3;
  strcpy(func_list[40].name, "Div Elements");
  func_list[41].fp = cdiv_elements;
  func_list[41].param_cnt = 3;
  strcpy(func_list[41].name, "cDiv Elements");
  
  func_list[func_list_size-9].fp = pow_mat;
  func_list[func_list_size-9].param_cnt = 2;
  strcpy(func_list[func_list_size-9].name, "Pow Mat");
  func_list[func_list_size-8].fp = pow_cmat;
  func_list[func_list_size-8].param_cnt = 2;
  strcpy(func_list[func_list_size-8].name, "Pow cMat");
  func_list[func_list_size-7].fp = cpow_cmat;
  func_list[func_list_size-7].param_cnt = 2;
  strcpy(func_list[func_list_size-7].name, "cPow cMat");
  
  func_list[func_list_size-6].fp = add;
  func_list[func_list_size-6].param_cnt = 2;
  strcpy(func_list[func_list_size-6].name, "Add");
  func_list[func_list_size-5].fp = cadd;
  func_list[func_list_size-5].param_cnt = 2;
  strcpy(func_list[func_list_size-5].name, "cAdd");
  func_list[func_list_size-4].fp = uadd;
  func_list[func_list_size-4].param_cnt = 2;
  strcpy(func_list[func_list_size-4].name, "uAdd");
  
  func_list[func_list_size-3].fp = scal;
  func_list[func_list_size-3].param_cnt = 2;
  strcpy(func_list[func_list_size-3].name, "Scale");
  func_list[func_list_size-2].fp = cscal;
  func_list[func_list_size-2].param_cnt = 2;
  strcpy(func_list[func_list_size-2].name, "cScale");
  func_list[func_list_size-1].fp = uscal;
  func_list[func_list_size-1].param_cnt = 2;
  strcpy(func_list[func_list_size-1].name, "uScale");
}

void read_ans(int is_ctype, char *func_name, char *ans_name, void *data_d, void *data_c, void **result_mat){
  char ans_dtype[16] = "d_4_4_2_.bin";
  char ans_ctype[16] = "c_4_4_2_.bin";

  if(is_ctype == 0){
    append_post(ans_dtype, func_name, ans_name);
    read_mat(ans_name, data_d);
    *result_mat = data_d;
  }
  else if(is_ctype == 1){
    append_post(ans_ctype, func_name, ans_name);
    read_cmat(ans_name, data_c);
    *result_mat = data_c;
  }
}


void read_broad_ans(int is_ctype, int param_type, char *func_name, char *ans_name, void *data_d, void *data_c, void **result_mat){
  char ans_broad[32] = "Broadcast_.bin";

  sprintf(ans_broad, "Broadcast_%d_.bin", param_type);
  append_post(ans_broad, func_name, ans_name);

  if(is_ctype == 0){
    read_mat(ans_name, data_d);
    *result_mat = data_d;
  }
  else if(is_ctype == 1){
    read_cmat(ans_name, data_c);
    *result_mat = data_c;
  }
}

void print_compare(char *func_name, int is_ctype, void *mat1, void *mat2){
  // print result
  printf("\n # %s\n", func_name);
  if(is_ctype == 0){
    if(compare_mat(mat1, mat2) == 0)
      printf("  !! FAIL !!\n");
    else
      printf("  OK\n");
  }
  else{
    if(compare_cmat(mat1, mat2) == 0)
      printf("  !! FAIL !!\n");
    else
      printf("  OK\n");
  }
}
