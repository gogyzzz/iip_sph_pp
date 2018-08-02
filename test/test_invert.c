#include "mother.h"

int main() {
  MAT *A2, *A3, *A4, *A5, *A6,*A7;
  MAT *B2, *B3, *B4, *B5, *B6,*B7;
  MAT *C2,*C3,*C4,*C5,*C6,*C7;
  ITER i;
  int *idx;
  int work;
  int seven = 7;
  UINT info;
  A2 = zeros(2, 2);
  B2 = zeros(2, 2);
  C2 = zeros(2,2);


  A3 = zeros(3, 3);
  B3 = zeros(3, 3);
  C3 = zeros(3,3);

  A4 = zeros(4, 4);
  B4 = zeros(4, 4);

  A5 = zeros(5, 5);
  B5 = zeros(5, 5);

  A6 = zeros(6, 6);
  B6 = zeros(6, 6);

  A7 = zeros(7,7);
  B7 = zeros(7,7);
  C7 = zeros(7,7);

  for (i = 0; i < 4; i++) A2->data[i] = i + 1;
  for (i = 0; i < 9; i++) A3->data[i] = i + 1 + (i + 1) * (i + 1);

  A4->data[0] = 5;
  A4->data[1] = 0;
  A4->data[2] = -1;
  A4->data[3] = 1;
  A4->data[4] = 4;
  A4->data[5] = 1;
  A4->data[6] = -1;
  A4->data[7] = 1;
  A4->data[8] = 2;
  A4->data[9] = -1;
  A4->data[10] = 3;
  A4->data[11] = -1;
  A4->data[12] = 1;
  A4->data[13] = -1;
  A4->data[14] = 0;
  A4->data[15] = 2;

  A5->data[0] = 1;
  A5->data[2] = 2;
  A5->data[4] = 5;
  A5->data[5] = -3;
  A5->data[8] = 4;
  A5->data[11] = -2;
  A5->data[14] = -5;
  A5->data[15] = -1;
  A5->data[18] = -4;
  A5->data[21] = 3;
  A5->data[24] = 6;

  A6->data[0] = 1;
  A6->data[2] = 2;
  A6->data[4] = 5;
  A6->data[6] = -3;
  A6->data[9] = 4;
  A6->data[11] = 1;
  A6->data[13] = 4;
  A6->data[15] = 6;
  A6->data[16] = 1;
  A6->data[17] = 2;
  A6->data[18] = 3;
  A6->data[19] = -2;
  A6->data[20] = 5;
  A6->data[22] = -5;
  A6->data[23] = 3;
  A6->data[24] = -1;
  A6->data[27] = -4;
  A6->data[29] = 4;
  A6->data[31] = 3;
  A6->data[34] = 6;
  A6->data[35] = 5;

  A7->data[ 0] = -8 ;
  A7->data[ 1] =  4;
  A7->data[ 2] = -50 ;
  A7->data[ 3] =  -82;
  A7->data[ 4] = 2 ;
  A7->data[ 5] = -56 ;
  A7->data[ 6] = 48 ;
  A7->data[ 7] =  -13;
  A7->data[ 8] = 7 ;
  A7->data[ 9] = -60 ;

  A7->data[ 10] = -98 ;
  A7->data[ 11] = 1 ;
  A7->data[ 12] = -66 ;
  A7->data[ 13] = 59 ;
  A7->data[ 14] =  15;
  A7->data[ 15] =  -6;
  A7->data[ 16] = 74 ;
  A7->data[ 17] = 118 ;
  A7->data[ 18] =  -2;
  A7->data[ 19] =  80;

  A7->data[ 20] = -70 ;
  A7->data[ 21] =  -16;
  A7->data[ 22] =  7;
  A7->data[ 23] =  -77;
  A7->data[ 24] =  -125;
  A7->data[ 25] = 2 ;
  A7->data[ 26] = -86 ;
  A7->data[ 27] =  75;
  A7->data[ 28] =  6;
  A7->data[ 9] = -2 ;

  A7->data[ 30] = 29 ;
  A7->data[ 31] =  47;
  A7->data[ 32] =  1;
  A7->data[ 33] =  32;
  A7->data[ 34] =  -28;
  A7->data[ 35] =  9;
  A7->data[ 36] =  -4;
  A7->data[ 37] =  43;
  A7->data[ 38] =  71;
  A7->data[ 39] =  -1;

  A7->data[ 40] =  50;
  A7->data[ 41] =  -42;
  A7->data[ 42] =  -3;
  A7->data[ 43] =  1;
  A7->data[ 44] =  -13;
  A7->data[ 45] =  -21;
  A7->data[ 46] =  0;
  A7->data[ 47] =  -14;
  A7->data[ 48] =  15;
  
  print_MAT(A2);
  wiki(A2, B2);
  print_MAT(B2);
  for (i = 0; i < 4; i++) A2->data[i] = i + 1;
  matmul(A2,B2,C2);
  print_MAT(C2);

  print_MAT(A3);
  wiki(A3, B3);
  print_MAT(B3);
  for (i = 0; i < 9; i++) A3->data[i] = i + 1 + (i + 1) * (i + 1);
  matmul(A3,B3,C3);
  print_MAT(C3);

  C4 = zeros(4,4);
  print_MAT(A4);
  printf("==== LU ====\n");
  fill(B4,0);
  wiki(A4, B4);

  print_MAT(B4);
  fill(A4,0); 
  A4->data[0] = 5;
  A4->data[1] = 0;
  A4->data[2] = -1;
  A4->data[3] = 1;
  A4->data[4] = 4;
  A4->data[5] = 1;
  A4->data[6] = -1;
  A4->data[7] = 1;
  A4->data[8] = 2;
  A4->data[9] = -1;
  A4->data[10] = 3;
  A4->data[11] = -1;
  A4->data[12] = 1;
  A4->data[13] = -1;
  A4->data[14] = 0;
  A4->data[15] = 2;


   matmul(A4,B4,C4);
   print_MAT(C4);
 
  printf("LAPACK : \n");  
  idx = (int*)malloc(sizeof(int)*4);
  LAPACK_dgetrf( &(A4->d0),&(A4->d0),(A4->data),&(A4->d0),idx,&info);
  printf("getrf info : %d\n",info);
  print_MAT(A4);
  for(i=0;i<4;i++)
    printf("%d ",idx[i]);
  printf("\n");
  free(idx);
   // invert_5by5(A5, B5);

    print_MAT(A5);
   // lu_invert(A5,B5);
    wiki(A5,B5);
    print_MAT(B5);
    fill(A5,0);
  A5->data[0] = 1;
  A5->data[2] = 2;
  A5->data[4] = 5;
  A5->data[5] = -3;
  A5->data[8] = 4;
  A5->data[11] = -2;
  A5->data[14] = -5;
  A5->data[15] = -1;
  A5->data[18] = -4;
  A5->data[21] = 3;
  A5->data[24] = 6;
  C5 = zeros(5,5);
  matmul(A5,B5,C5);
   print_MAT(C5);

    C6 = zeros(6,6);   
    print_MAT(A6);
    wiki(A6, B6);
    print_MAT(B6);
 fill(A6,0);
    A6->data[0] = 1;
  A6->data[2] = 2;
  A6->data[4] = 5;
  A6->data[6] = -3;
  A6->data[9] = 4;
  A6->data[11] = 1;
  A6->data[13] = 4;
  A6->data[15] = 6;
  A6->data[16] = 1;
  A6->data[17] = 2;
  A6->data[18] = 3;
  A6->data[19] = -2;
  A6->data[20] = 5;
  A6->data[22] = -5;
  A6->data[23] = 3;
  A6->data[24] = -1;
  A6->data[27] = -4;
  A6->data[29] = 4;
  A6->data[31] = 3;
  A6->data[34] = 6;
  A6->data[35] = 5;


    matmul(A6,B6,C6);
    print_MAT(C6);

  printf("MY IMPLEMENT\n");
  print_MAT(A7);
  wiki(A7,B7);
  print_MAT(B7);
A7->data[ 0] = -8 ;
  A7->data[ 1] =  4;
  A7->data[ 2] = -50 ;
  A7->data[ 3] =  -82;
  A7->data[ 4] = 2 ;
  A7->data[ 5] = -56 ;
  A7->data[ 6] = 48 ;
  A7->data[ 7] =  -13;
  A7->data[ 8] = 7 ;
  A7->data[ 9] = -60 ;

  A7->data[ 10] = -98 ;
  A7->data[ 11] = 1 ;
  A7->data[ 12] = -66 ;
  A7->data[ 13] = 59 ;
  A7->data[ 14] =  15;
  A7->data[ 15] =  -6;
  A7->data[ 16] = 74 ;
  A7->data[ 17] = 118 ;
  A7->data[ 18] =  -2;
  A7->data[ 19] =  80;

  A7->data[ 20] = -70 ;
  A7->data[ 21] =  -16;
  A7->data[ 22] =  7;
  A7->data[ 23] =  -77;
  A7->data[ 24] =  -125;
  A7->data[ 25] = 2 ;
  A7->data[ 26] = -86 ;
  A7->data[ 27] =  75;
  A7->data[ 28] =  6;
  A7->data[ 9] = -2 ;

  A7->data[ 30] = 29 ;
  A7->data[ 31] =  47;
  A7->data[ 32] =  1;
  A7->data[ 33] =  32;
  A7->data[ 34] =  -28;
  A7->data[ 35] =  9;
  A7->data[ 36] =  -4;
  A7->data[ 37] =  43;
  A7->data[ 38] =  71;
  A7->data[ 39] =  -1;

  A7->data[ 40] =  50;
  A7->data[ 41] =  -42;
  A7->data[ 42] =  -3;
  A7->data[ 43] =  1;
  A7->data[ 44] =  -13;
  A7->data[ 45] =  -21;
  A7->data[ 46] =  0;
  A7->data[ 47] =  -14;
  A7->data[ 48] =  15;
  


  matmul(A7,B7,C7);
  print_MAT(C7);

  printf("LAPACK : \n");  
  idx = (int*)malloc(sizeof(int)*8);
  LAPACK_dgetrf( &(seven),&(seven),(A7->data),&(seven),idx,&info);
  printf("getrf info : %d\n",info);
  print_MAT(A7);
  for(i=0;i<7;i++)
    printf("%d ",idx[i]);
  printf("\n");
  work = 7*7;
  fill(B7,0);
  LAPACK_dgetri(&(seven),A7->data,&(seven),idx,B7->data,&work,&info);
  printf("getri work : %d info :%d\n",work,info);
  print_MAT(A7);
  print_MAT(B7);

  copy(A7,B7);
  A7->data[ 0] = -8 ;
  A7->data[ 1] =  4;
  A7->data[ 2] = -50 ;
  A7->data[ 3] =  -82;
  A7->data[ 4] = 2 ;
  A7->data[ 5] = -56 ;
  A7->data[ 6] = 48 ;
  A7->data[ 7] =  -13;
  A7->data[ 8] = 7 ;
  A7->data[ 9] = -60 ;

  A7->data[ 10] = -98 ;
  A7->data[ 11] = 1 ;
  A7->data[ 12] = -66 ;
  A7->data[ 13] = 59 ;
  A7->data[ 14] =  15;
  A7->data[ 15] =  -6;
  A7->data[ 16] = 74 ;
  A7->data[ 17] = 118 ;
  A7->data[ 18] =  -2;
  A7->data[ 19] =  80;

  A7->data[ 20] = -70 ;
  A7->data[ 21] =  -16;
  A7->data[ 22] =  7;
  A7->data[ 23] =  -77;
  A7->data[ 24] =  -125;
  A7->data[ 25] = 2 ;
  A7->data[ 26] = -86 ;
  A7->data[ 27] =  75;
  A7->data[ 28] =  6;
  A7->data[ 9] = -2 ;

  A7->data[ 30] = 29 ;
  A7->data[ 31] =  47;
  A7->data[ 32] =  1;
  A7->data[ 33] =  32;
  A7->data[ 34] =  -28;
  A7->data[ 35] =  9;
  A7->data[ 36] =  -4;
  A7->data[ 37] =  43;
  A7->data[ 38] =  71;
  A7->data[ 39] =  -1;

  A7->data[ 40] =  50;
  A7->data[ 41] =  -42;
  A7->data[ 42] =  -3;
  A7->data[ 43] =  1;
  A7->data[ 44] =  -13;
  A7->data[ 45] =  -21;
  A7->data[ 46] =  0;
  A7->data[ 47] =  -14;
  A7->data[ 48] =  15;
  fill(C7,0);    
  matmul(A7,B7,C7);
  print_MAT(C7);



  free_MAT(A2);
  free_MAT(B2);
  free_MAT(C2);
  free_MAT(A3);
  free_MAT(B3);
  free_MAT(C3);
  free_MAT(A4);
  free_MAT(B4);
  free_MAT(C4);
  free_MAT(A5);
  free_MAT(B5);
  free_MAT(C5);
  free_MAT(A6);
  free_MAT(B6);

  free_MAT(A7);
  free_MAT(B7);
  free_MAT(C7);
  return 0;
}
