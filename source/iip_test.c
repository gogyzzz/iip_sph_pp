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
const int func_list_size = 44;

UINT _eq(DTYPE A, DTYPE B) {
#if NTYPE == 0
  return fabsf(A - B) < FZERO ? 1 : 0;
#elif NTYPE == 1
  return fabs(A - B) < FZERO ? 1 : 0;
#endif
  return 0;
}

UINT compare_mat(MAT *A, MAT *B) {
  ITER i;
  for (i = 0; i < A->d0 * A->d1 * A->d2; i++)
    if (!_eq(A->data[i], B->data[i])) return 0;
  return 1;
}

UINT compare_cmat(CMAT *A, CMAT *B) {
  ITER i;
  for (i = 0; i < A->d0 * A->d1 * A->d2; i++) {
    if (!_eq(A->data[i].re, B->data[i].re)) return 0;
    if (!_eq(A->data[i].im, B->data[i].im)) return 0;
  }
  return 1;
}

void perform_test() {
  FILE *f;
  char cur[MAX_CHAR];
  f = NULL;
  f = fopen("../test/test.txt", "r");
  ASSERT_FILE(f, "../test/test.txt")
  while (fscanf(f, "%s\n", cur) != EOF)
  // while(fgets(cur,MAX_CHAR,f) != EOF)
  {
    printf("== %s == \n", cur);
    //do_test(cur);
  }

  fclose(f);
}

void append_post(char *filename, const char *post, char *out) {
  char t[MAX_CHAR] = "../test/";

  strcpy(out, filename);
  strtok(out, ".");
  strcat(t, out);
  strcat(t, post);
  strcat(t, ".bin");
  strcpy(out, t);
}


void test_verification(int heat){
  if(!is_init)
    init_list();
  printf("Verifying...\n");

  printf("Done.\n");
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
  A = zeros(6, 6, 1);
  B = zeros(6, 6, 1);
  C = NULL;

  read_mat("../test/d_6_6_1.bin", A);
  is_ctype = 0;
  copy(A, B);
  total = 0;

  if(heat == 1){
    preheat();
    test_performance(0, 0);
  }

  if(print_flag){
    printf(" *** Input Data ***\n");
    print_mat(A);
  }

  for(i = 0; i<func_list_size-6; i++){
    // continue for real / complex func
    if(is_ctype == 0 && i%2 != 0)
      continue;
    else if(is_ctype != 0 && i%2 == 0)
      continue;

    // set data
    copy(A, B);
    calced = 0;

    // call appropriate func
    if (func_list[i].param_cnt == 1) {
      calced = 1;
      stopwatch(0);
      func_list[i].fp(B);
      total += stopwatch(1);
    } else if (i == 4 || i == 5) {  // pow
      calced = 1;
      stopwatch(0);
      func_list[i].fp(B, 2.0);
      total += stopwatch(1);
    } else if (i == 25 || i == 26) {  // permute
      calced = 1;
      stopwatch(0);
      func_list[i].fp(B, 213);
      total += stopwatch(1);
    } else if (i == 35 || i == 36) {  // sum
      calced = 1;
      stopwatch(0);
      func_list[i].fp(B, 10.0);
      total += stopwatch(1);
    } else if (i >= 37) {  // scal
      calced = 1;
      stopwatch(0);
      func_list[i].fp(5.0, B);
      total += stopwatch(1);
    }

    // print result
    if (calced != 0 && print_flag != 0) {
      printf("\n # %s\n", func_list[i].name);
      print_mat(B);
    }
  }

  if(print_flag)
    printf("Done.\n *** Total Run time : %.3lfms\n", (double)total/1000.0);
}

void test_viewlist(){
  if(!is_init)
    init_list();
  int i=0;
  printf("\n\n---- Func list ----\n");
  for(i=0; i<func_list_size-6; i++){
    printf("戍成 %s\n", func_list[i++].name);
    printf("弛戌式 %s\n", func_list[i].name);
    printf("弛\n");
  }
  i=func_list_size-6;
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
  for(i=0; i<1000000000; i++){
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
  
  func_list[2].fp = trans;
  func_list[2].param_cnt = 1;
  strcpy(func_list[2].name, "Transpose");
  func_list[3].fp = ctrans;
  func_list[3].param_cnt = 1;
  strcpy(func_list[3].name, "cTranspose");
  
  func_list[4].fp = pow_mat;
  func_list[4].param_cnt = 2;
  strcpy(func_list[4].name, "Pow");
  func_list[5].fp = pow_cmat;
  func_list[5].param_cnt = 2;
  strcpy(func_list[5].name, "cPow");
  
  func_list[6].fp = sqrt_mat;
  func_list[6].param_cnt = 1;
  strcpy(func_list[6].name, "Sqrt");
  func_list[7].fp = sqrt_cmat;
  func_list[7].param_cnt = 1;
  strcpy(func_list[7].name, "cSqrt");
  
  func_list[8].fp = round_mat;
  func_list[8].param_cnt = 1;
  strcpy(func_list[8].name, "Round");
  func_list[9].fp = round_cmat;
  func_list[9].param_cnt = 1;
  strcpy(func_list[9].name, "cRound");
  
  func_list[10].fp = floor_mat;
  func_list[10].param_cnt = 1;
  strcpy(func_list[10].name, "Floor");
  func_list[11].fp = floor_cmat;
  func_list[11].param_cnt = 1;
  strcpy(func_list[11].name, "cFloor");
  
  func_list[12].fp = ceil_mat;
  func_list[12].param_cnt = 1;
  strcpy(func_list[12].name, "Ceil");
  func_list[13].fp = ceil_cmat;
  func_list[13].param_cnt = 1;
  strcpy(func_list[13].name, "cCeil");
  
  func_list[14].fp = log_mat;
  func_list[14].param_cnt = 1;
  strcpy(func_list[14].name, "Log");
  func_list[15].fp = log_cmat;
  func_list[15].param_cnt = 1;
  strcpy(func_list[15].name, "cLog");
  
  func_list[16].fp = log2_mat;
  func_list[16].param_cnt = 1;
  strcpy(func_list[16].name, "Log2");
  func_list[17].fp = log2_cmat;
  func_list[17].param_cnt = 1;
  strcpy(func_list[17].name, "cLog2");
  
  func_list[18].fp = log10_mat;
  func_list[18].param_cnt = 1;
  strcpy(func_list[18].name, "Log10");
  func_list[19].fp = log10_cmat;
  func_list[19].param_cnt = 1;
  strcpy(func_list[19].name, "cLog10");
  
  func_list[20].fp = exp_mat;
  func_list[20].param_cnt = 1;
  strcpy(func_list[20].name, "Exp");
  func_list[21].fp = exp_cmat;
  func_list[21].param_cnt = 1;
  strcpy(func_list[21].name, "cExp");
  
  func_list[22].fp = abs_mat;
  func_list[22].param_cnt = 1;
  strcpy(func_list[22].name, "Abs");
  func_list[23].fp = abs_cmat;
  func_list[23].param_cnt = 1;
  strcpy(func_list[23].name, "cAbs");
  
  func_list[24].fp = repmat;
  func_list[24].param_cnt = 2;
  strcpy(func_list[24].name, "Repeat");
  func_list[25].fp = crepmat;
  func_list[25].param_cnt = 2;
  strcpy(func_list[25].name, "cRepeat");
  
  func_list[26].fp = permute;
  func_list[26].param_cnt = 2;
  strcpy(func_list[26].name, "Permute");
  func_list[27].fp = cpermute;
  func_list[27].param_cnt = 2;
  strcpy(func_list[27].name, "cPermute");
  
  func_list[28].fp = reshape;
  func_list[28].param_cnt = 2;
  strcpy(func_list[28].name, "Reshape");
  func_list[29].fp = creshape;
  func_list[29].param_cnt = 2;
  strcpy(func_list[29].name, "cReshape");
  
  func_list[30].fp = gemm;
  func_list[30].param_cnt = 6;
  strcpy(func_list[30].name, "Gemm");
  func_list[31].fp = cgemm;
  func_list[31].param_cnt = 6;
  strcpy(func_list[31].name, "cGemm");
  
  func_list[32].fp = diagonal;
  func_list[32].param_cnt = 2;
  strcpy(func_list[32].name, "Diagonal");
  func_list[33].fp = cdiagonal;
  func_list[33].param_cnt = 2;
  strcpy(func_list[33].name, "cDiagonal");
  
  func_list[34].fp = trace;
  func_list[34].param_cnt = 2;
  strcpy(func_list[34].name, "Trace");
  func_list[35].fp = ctrace;
  func_list[35].param_cnt = 2;
  strcpy(func_list[35].name, "cTrace");
  
  func_list[36].fp = asum;
  func_list[36].param_cnt = 2;
  strcpy(func_list[36].name, "Asum");
  func_list[37].fp = casum;
  func_list[37].param_cnt = 2;
  strcpy(func_list[37].name, "cAsum");
  
  func_list[38].fp = scal;
  func_list[38].param_cnt = 2;
  strcpy(func_list[38].name, "Scale");
  func_list[39].fp = cscal;
  func_list[39].param_cnt = 2;
  strcpy(func_list[39].name, "cScale");
  func_list[40].fp = uscal;
  func_list[40].param_cnt = 2;
  strcpy(func_list[40].name, "uScale");
  
  func_list[41].fp = add;
  func_list[41].param_cnt = 3;
  strcpy(func_list[41].name, "Add");
  func_list[42].fp = cadd;
  func_list[42].param_cnt = 3;
  strcpy(func_list[42].name, "cAdd");
  func_list[43].fp = uadd;
  func_list[43].param_cnt = 3;
  strcpy(func_list[43].name, "uAdd");
}