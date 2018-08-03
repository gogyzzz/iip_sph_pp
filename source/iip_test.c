#include "iip_test.h"

UINT _eq(DTYPE A,DTYPE B){
#if NTYPE == 0
  return fabsf(A-B)<FZERO?1:0;
#elif NTYPE == 1
  return fabs(A-B)<FZERO?1:0;
#endif
  return 0;
}

UINT compare_mat(MAT*A,MAT*B){
ITER i;
  for(i=0;i<A->d0*A->d1*A->d2;i++)
    if(!_eq(A->data[i],B->data[i])) return 0;
  return 1;
}

UINT compare_cmat(CMAT*A,CMAT*B){
ITER i;
  for(i=0;i<A->d0*A->d1*A->d2;i++){
    if(!_eq(A->data[i].re,B->data[i].re)) return 0;
    if(!_eq(A->data[i].im,B->data[i].im)) return 0;
  }
  return 1;
}

void perform_test(){
  FILE * f;
  char cur[MAX_CHAR];
  f = NULL;
  f = fopen("../test/test.txt","r");
  if(f == NULL) ASSERT(NO_FILE)
  while(fscanf(f,"%s\n",cur) != EOF)
  //while(fgets(cur,MAX_CHAR,f) != EOF)
  {
    printf("== %s == \n",cur);
    do_test(cur);
  }
  
  fclose(f);

}

void do_test(char*name){
  UINT d0,d1,d2;
  char type;
  MAT*A,*B,*C;
  CMAT*CA,*CB,*CC;
  char filepath[MAX_CHAR] = "../test/";
  char testpath[MAX_CHAR] = "";
  char *t;
  FILE *f;
  
  UINT f_inv,f_tran;
  UINT sq=0;
  f_inv =1,f_tran=1;

  strcat(filepath,name);
  sscanf(name,"%c_%d_%d_%d\n",&type,&d0,&d1,&d2);
  
  if(type == 'd'){
    
    A = alloc_MAT(d0,d1,d2);
    B = alloc_MAT(d0,d1,d2);
    C = alloc_MAT(d0,d1,d2);
    
    read_mat(filepath,A); 
    //CHECK SQAURE
    if(d0 == d1){
      sq = 1;
    //READ original matrix to A
    //printf("file path : %s\n",filepath);
    if(!compare_mat(B,C))f_tran=0;
    post_append(name,"_inv",testpath); 
    //READ matlab inversed matrix to C
    read_mat(testpath,C);
    //invert A to B
    invert(A,B);
/*
    printf("==== B ====\n");
    print_MAT(B);
    printf("==== C ====\n");
    print_MAT(C);
*/
    //Compare B and C
    if(!compare_mat(B,C))f_inv=0;
    }

    //TRANSPOSE
    free_MAT(C);
    copy(A,B);
    C = alloc_MAT(d1,d0,d2);

    post_append(name,"_tran",testpath);
    read_mat(testpath,C);
    trans(B);
    if(!compare_mat(B,C))f_tran=0;
    free_MAT(A);
    free_MAT(B);
    free_MAT(C);
  }
  else if(type == 'c'){
  
  }
  else ASSERT(ARG_INVAL)

  if(sq){
    if(f_inv) printf("INVERT : PASS\n");
    else      printf("INVERT : FAIL\n");
  }

  if(f_tran)  printf("TRANS  : PASS\n");
  else        printf("TRANS  : FAIL\n");
}

void post_append(char*filename, const char*post,char* out){
 char t[MAX_CHAR]="../test/";
 
 strcpy(out,filename);
 strtok(out,"."); 
 strcat(t,out);
 strcat(t,post);
 strcat(t,".bin");
 strcpy(out,t);

}
