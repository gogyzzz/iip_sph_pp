#include "mother.h"

int main() {
  MAT *A2, *A3, *A4,*A5,*A6;
  MAT *B2, *B3, *B4,*B5,*B6;
  ITER i;

  A2 = zeros(2, 2);
  B2 = zeros(2, 2);

  A3 = zeros(3, 3);
  B3 = zeros(3, 3);

  A4 = zeros(4, 4);
  B4 = zeros(4, 4);
  
  A5 = zeros(5, 5);
  B5 = zeros(5, 5);
  
  A6 = zeros(6, 6);
  B6 = zeros(6, 6);

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

  A5->data[0] = 1  ;
  A5->data[2 ] = 2  ;
  A5->data[4 ] =   5;
  A5->data[5 ] =  -3 ;
  A5->data[8 ] =   4;
  A5->data[11 ] = -2  ;
  A5->data[ 14] =-5   ;
  A5->data[15 ] =  -1 ;
  A5->data[18 ] =  -4 ;
  A5->data[21 ] =  3 ;
  A5->data[ 24] = 6  ;

  A6->data[0 ] = 1  ;
  A6->data[2 ] = 2  ;
  A6->data[4 ] = 5  ;
  A6->data[6 ] = -3  ;
  A6->data[9 ] =  4 ;
  A6->data[11 ] =  1 ;
  A6->data[ 13] =  4 ;
  A6->data[ 15] =  6 ;
  A6->data[16 ] =  1 ;
  A6->data[17 ] =  2 ;
  A6->data[18 ] =  3 ;
  A6->data[19 ] =  -2 ;
  A6->data[20 ] =   5;
  A6->data[22 ] =   -5;
  A6->data[23 ] = 3  ;
  A6->data[24 ] = -1  ;
  A6->data[27 ] =  -4 ;
  A6->data[29 ] =  4 ;
  A6->data[31 ] =  3 ;
  A6->data[34 ] =  6 ;
  A6->data[35 ] =  5 ;


  print_MAT(A2);
  invert_2by2(A2, B2);
  print_MAT(B2);
 
  print_MAT(A3);
  invert_3by3(A3, B3);
  print_MAT(B3);
 
  print_MAT(A4);
  invert_4by4(A4, B4);
  print_MAT(B4);
 
  print_MAT(A5);
  invert_5by5(A5, B5);
  print_MAT(B5);
 
  print_MAT(A6);
  invert_6by6(A6, B6);
  print_MAT(B6);

  free_MAT(A2);
  free_MAT(B2);
  free_MAT(A3);
  free_MAT(B3);
  free_MAT(A4);
  free_MAT(B4);
  free_MAT(A5);
  free_MAT(B5);
  free_MAT(A6);
  free_MAT(B6);
  return 0;
}
