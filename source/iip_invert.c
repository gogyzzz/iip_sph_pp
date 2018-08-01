#include "iip_invert.h"

/**** get inverse of matrix ****/
void invert_2by2(MAT* mat, MAT* inv) {
  DTYPE det;

  det = mat->data[0] * mat->data[3] - mat->data[1] * mat->data[2];

  if (det > -FZERO && det < FZERO) ASSERT(DET_FAIL)
  det = 1 / det;
  inv->data[0] = mat->data[3] * det;
  inv->data[3] = mat->data[0] * det;
  inv->data[1] = -mat->data[2] * det;
  inv->data[2] = -mat->data[1] * det;
}

void invert_3by3(MAT* mat, MAT* inv) {
  DTYPE det;

  inv->data[0] = mat->data[4] * mat->data[8] - mat->data[7] * mat->data[5];
  inv->data[1] = mat->data[7] * mat->data[2] - mat->data[1] * mat->data[8];
  inv->data[2] = mat->data[1] * mat->data[5] - mat->data[4] * mat->data[2];

  det = mat->data[0] * inv->data[0] + mat->data[3] * inv->data[1] +
        mat->data[6] * inv->data[2];

  if (det > -FZERO && det < FZERO) ASSERT(DET_FAIL)

  inv->data[3] = mat->data[6] * mat->data[5] - mat->data[3] * mat->data[8];
  inv->data[4] = mat->data[0] * mat->data[6] - mat->data[6] * mat->data[2];
  inv->data[5] = mat->data[3] * mat->data[2] - mat->data[0] * mat->data[5];
  inv->data[6] = mat->data[3] * mat->data[7] - mat->data[6] * mat->data[4];
  inv->data[7] = mat->data[6] * mat->data[1] - mat->data[0] * mat->data[7];
  inv->data[8] = mat->data[0] * mat->data[4] - mat->data[3] * mat->data[1];

  det = 1 / det;
  inv->data[0] = inv->data[0] * det;
  inv->data[1] = inv->data[1] * det;
  inv->data[2] = inv->data[2] * det;
  inv->data[3] = inv->data[3] * det;
  inv->data[4] = inv->data[4] * det;
  inv->data[5] = inv->data[5] * det;
  inv->data[6] = inv->data[6] * det;
  inv->data[7] = inv->data[7] * det;
  inv->data[8] = inv->data[8] * det;
}

void invert_4by4(MAT* mat, MAT* inv) {
  DTYPE det;
  DTYPE t1, t2, t3, t4, t5;

  t1 = mat->data[10] * mat->data[15] - mat->data[14] * mat->data[11];
  t2 = mat->data[6] * mat->data[15] - mat->data[14] * mat->data[7];
  t3 = mat->data[6] * mat->data[11] - mat->data[10] * mat->data[7];

  inv->data[0] = mat->data[5] * t1 - mat->data[9] * t2 + mat->data[13] * t3;
  inv->data[4] = mat->data[8] * t2 - mat->data[4] * t1 - mat->data[12] * t3;

  t4 = mat->data[2] * mat->data[15] - mat->data[14] * mat->data[3];
  t5 = mat->data[2] * mat->data[11] - mat->data[10] * mat->data[3];

  inv->data[1] = mat->data[9] * t4 - mat->data[1] * t1 - mat->data[13] * t5;
  inv->data[5] = mat->data[0] * t1 - mat->data[8] * t4 + mat->data[12] * t5;

  t1 = mat->data[2] * mat->data[7] - mat->data[6] * mat->data[3];

  inv->data[2] = mat->data[1] * t2 - mat->data[5] * t4 + mat->data[13] * t1;
  inv->data[6] = mat->data[4] * t4 - mat->data[0] * t2 - mat->data[12] * t1;
  inv->data[3] = mat->data[5] * t5 - mat->data[1] * t3 - mat->data[9] * t1;

  det = mat->data[0] * inv->data[0] + mat->data[4] * inv->data[1] +
        mat->data[8] * inv->data[2] + mat->data[12] * inv->data[3];
#if DEBUG
  printf("det : %lf\n", det);
#endif
  if (det > -FZERO && det < FZERO) ASSERT(DET_FAIL)

  inv->data[7] = mat->data[0] * t3 - mat->data[4] * t5 + mat->data[8] * t1;

  t1 = mat->data[8] * mat->data[13] - mat->data[12] * mat->data[9];
  t2 = mat->data[4] * mat->data[13] - mat->data[12] * mat->data[5];
  t3 = mat->data[4] * mat->data[9] - mat->data[8] * mat->data[5];

  inv->data[8] = mat->data[7] * t1 - mat->data[11] * t2 + mat->data[15] * t3;
  inv->data[12] = mat->data[10] * t2 - mat->data[6] * t1 - mat->data[14] * t3;

  t4 = mat->data[0] * mat->data[13] - mat->data[12] * mat->data[1];
  t5 = mat->data[0] * mat->data[9] - mat->data[8] * mat->data[1];

  inv->data[9] = mat->data[11] * t4 - mat->data[3] * t1 - mat->data[15] * t5;
  inv->data[13] = mat->data[2] * t1 - mat->data[10] * t4 + mat->data[14] * t5;

  t1 = mat->data[0] * mat->data[5] - mat->data[4] * mat->data[1];

  inv->data[10] = mat->data[3] * t2 - mat->data[7] * t4 + mat->data[15] * t1;
  inv->data[14] = mat->data[6] * t4 - mat->data[2] * t2 - mat->data[14] * t1;
  inv->data[11] = mat->data[7] * t5 - mat->data[3] * t3 - mat->data[11] * t1;
  inv->data[15] = mat->data[2] * t3 - mat->data[6] * t5 + mat->data[10] * t1;

  det = 1. / det;

  inv->data[0] = inv->data[0] * det;
  inv->data[1] = inv->data[1] * det;
  inv->data[2] = inv->data[2] * det;
  inv->data[3] = inv->data[3] * det;
  inv->data[4] = inv->data[4] * det;
  inv->data[5] = inv->data[5] * det;
  inv->data[6] = inv->data[6] * det;
  inv->data[7] = inv->data[7] * det;
  inv->data[8] = inv->data[8] * det;
  inv->data[9] = inv->data[9] * det;
  inv->data[10] = inv->data[10] * det;
  inv->data[11] = inv->data[11] * det;
  inv->data[12] = inv->data[12] * det;
  inv->data[13] = inv->data[13] * det;
  inv->data[14] = inv->data[14] * det;
  inv->data[15] = inv->data[15] * det;
}

void invert_5by5(MAT* mat, MAT* inv) {
  DTYPE det;
  DTYPE t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16,
      t17, t18, t19, t20;

  t1 = mat->data[18] * mat->data[24] - mat->data[23] * mat->data[19];
  t2 = mat->data[13] * mat->data[24] - mat->data[23] * mat->data[14];
  t3 = mat->data[13] * mat->data[19] - mat->data[18] * mat->data[14];
  t4 = mat->data[8] * mat->data[24] - mat->data[23] * mat->data[9];
  t5 = mat->data[8] * mat->data[19] - mat->data[18] * mat->data[9];
  t6 = mat->data[8] * mat->data[14] - mat->data[13] * mat->data[9];
  t7 = mat->data[3] * mat->data[24] - mat->data[23] * mat->data[4];
  t8 = mat->data[3] * mat->data[19] - mat->data[18] * mat->data[4];
  t9 = mat->data[3] * mat->data[14] - mat->data[13] * mat->data[4];
  t10 = mat->data[3] * mat->data[9] - mat->data[8] * mat->data[4];

  t11 = mat->data[12] * t1 - mat->data[17] * t2 + mat->data[22] * t3;
  t12 = mat->data[7] * t1 - mat->data[17] * t4 + mat->data[22] * t5;
  t13 = mat->data[7] * t2 - mat->data[12] * t4 + mat->data[22] * t6;
  t14 = mat->data[7] * t3 - mat->data[12] * t5 + mat->data[17] * t6;
  t15 = mat->data[2] * t1 - mat->data[17] * t7 + mat->data[22] * t8;
  t16 = mat->data[2] * t2 - mat->data[12] * t7 + mat->data[22] * t9;
  t17 = mat->data[2] * t3 - mat->data[12] * t8 + mat->data[17] * t9;

  inv->data[0] = mat->data[6] * t11 - mat->data[11] * t12 +
                 mat->data[16] * t13 - mat->data[21] * t14;
  inv->data[5] = -mat->data[5] * t11 + mat->data[10] * t12 -
                 mat->data[15] * t13 + mat->data[20] * t14;
  inv->data[1] = -mat->data[1] * t11 + mat->data[11] * t15 -
                 mat->data[16] * t16 + mat->data[21] * t17;
  inv->data[6] = mat->data[0] * t11 - mat->data[10] * t15 +
                 mat->data[15] * t16 - mat->data[20] * t17;

  t18 = mat->data[2] * t4 - mat->data[7] * t7 + mat->data[22] * t10;
  t19 = mat->data[2] * t5 - mat->data[7] * t8 + mat->data[17] * t10;
  t20 = mat->data[2] * t6 - mat->data[7] * t9 + mat->data[12] * t10;

  inv->data[2] = mat->data[1] * t12 - mat->data[6] * t15 + mat->data[16] * t18 -
                 mat->data[21] * t19;
  inv->data[7] = -mat->data[0] * t12 + mat->data[5] * t15 -
                 mat->data[15] * t18 + mat->data[20] * t19;
  inv->data[3] = -mat->data[1] * t13 + mat->data[6] * t16 -
                 mat->data[11] * t18 + mat->data[21] * t20;
  inv->data[8] = mat->data[0] * t13 - mat->data[5] * t16 + mat->data[10] * t18 -
                 mat->data[20] * t20;
  inv->data[4] = mat->data[1] * t14 - mat->data[6] * t17 + mat->data[11] * t19 -
                 mat->data[16] * t20;
  inv->data[9] = -mat->data[0] * t14 + mat->data[5] * t17 -
                 mat->data[10] * t19 + mat->data[15] * t20;

  det = mat->data[0] * inv->data[0] + mat->data[5] * inv->data[1] +
        mat->data[10] * inv->data[2] + mat->data[15] * inv->data[3] +
        mat->data[20] * inv->data[4];

#if DEBUG
  printf("det : %lf\n", det);
#endif
  if (det > -FZERO && det < FZERO) ASSERT(DET_FAIL)
  t11 = mat->data[11] * t1 - mat->data[16] * t2 + mat->data[21] * t3;
  t12 = mat->data[6] * t1 - mat->data[16] * t4 + mat->data[21] * t5;
  t13 = mat->data[6] * t2 - mat->data[11] * t4 + mat->data[21] * t6;
  t14 = mat->data[6] * t3 - mat->data[11] * t5 + mat->data[16] * t6;
  t15 = mat->data[1] * t1 - mat->data[16] * t7 + mat->data[21] * t8;
  t16 = mat->data[1] * t2 - mat->data[11] * t7 + mat->data[21] * t9;
  t17 = mat->data[1] * t3 - mat->data[11] * t8 + mat->data[16] * t9;
  t18 = mat->data[1] * t4 - mat->data[6] * t7 + mat->data[21] * t10;
  t19 = mat->data[1] * t5 - mat->data[6] * t8 + mat->data[16] * t10;

  inv->data[10] = mat->data[5] * t11 - mat->data[10] * t12 +
                  mat->data[15] * t13 - mat->data[20] * t14;
  inv->data[11] = -mat->data[0] * t11 + mat->data[10] * t15 -
                  mat->data[15] * t16 + mat->data[20] * t17;
  inv->data[12] = mat->data[0] * t12 - mat->data[5] * t15 +
                  mat->data[15] * t18 - mat->data[20] * t19;

  t1 = mat->data[10] * mat->data[16] - mat->data[15] * mat->data[11];
  t2 = mat->data[5] * mat->data[16] - mat->data[15] * mat->data[6];
  t3 = mat->data[5] * mat->data[11] - mat->data[10] * mat->data[6];
  t4 = mat->data[0] * mat->data[16] - mat->data[15] * mat->data[1];
  t5 = mat->data[0] * mat->data[11] - mat->data[10] * mat->data[1];
  t6 = mat->data[0] * mat->data[6] - mat->data[5] * mat->data[1];
  t7 = mat->data[10] * mat->data[21] - mat->data[20] * mat->data[11];
  t8 = mat->data[5] * mat->data[21] - mat->data[20] * mat->data[6];
  t9 = mat->data[0] * mat->data[21] - mat->data[20] * mat->data[1];
  t10 = mat->data[15] * mat->data[21] - mat->data[20] * mat->data[16];

  t11 = mat->data[12] * t10 - mat->data[17] * t7 + mat->data[22] * t1;
  t12 = mat->data[7] * t10 - mat->data[17] * t8 + mat->data[22] * t2;
  t13 = mat->data[7] * t7 - mat->data[12] * t8 + mat->data[22] * t3;
  t14 = mat->data[7] * t1 - mat->data[12] * t2 + mat->data[17] * t3;
  t15 = mat->data[2] * t10 - mat->data[17] * t9 + mat->data[22] * t4;
  t16 = mat->data[2] * t7 - mat->data[12] * t9 + mat->data[22] * t5;
  t17 = mat->data[2] * t1 - mat->data[12] * t4 + mat->data[17] * t5;

  inv->data[15] = mat->data[9] * t11 - mat->data[14] * t12 +
                  mat->data[19] * t13 - mat->data[24] * t14;
  inv->data[20] = -mat->data[8] * t11 + mat->data[13] * t12 -
                  mat->data[18] * t13 + mat->data[23] * t14;
  inv->data[16] = -mat->data[4] * t11 + mat->data[14] * t15 -
                  mat->data[19] * t16 + mat->data[24] * t17;
  inv->data[21] = mat->data[3] * t11 - mat->data[13] * t15 +
                  mat->data[18] * t16 - mat->data[23] * t17;

  t18 = mat->data[2] * t8 - mat->data[7] * t9 + mat->data[22] * t6;
  t19 = mat->data[2] * t2 - mat->data[7] * t4 + mat->data[17] * t6;
  t20 = mat->data[2] * t3 - mat->data[7] * t5 + mat->data[12] * t6;

  inv->data[17] = mat->data[4] * t12 - mat->data[9] * t15 +
                  mat->data[19] * t18 - mat->data[24] * t19;
  inv->data[22] = -mat->data[3] * t12 + mat->data[8] * t15 -
                  mat->data[18] * t18 + mat->data[23] * t19;
  inv->data[18] = -mat->data[4] * t13 + mat->data[9] * t16 -
                  mat->data[14] * t18 + mat->data[24] * t20;
  inv->data[23] = mat->data[3] * t13 - mat->data[8] * t16 +
                  mat->data[13] * t18 - mat->data[23] * t20;
  inv->data[19] = mat->data[4] * t14 - mat->data[9] * t17 +
                  mat->data[14] * t19 - mat->data[19] * t20;
  inv->data[24] = -mat->data[3] * t14 + mat->data[8] * t17 -
                  mat->data[13] * t19 + mat->data[18] * t20;

  t11 = mat->data[8] * t7 - mat->data[13] * t8 + mat->data[23] * t3;
  t12 = mat->data[3] * t7 - mat->data[13] * t9 + mat->data[23] * t5;
  t13 = mat->data[3] * t8 - mat->data[8] * t9 + mat->data[23] * t6;
  t14 = mat->data[3] * t3 - mat->data[8] * t5 + mat->data[13] * t6;

  t15 = mat->data[8] * t1 - mat->data[13] * t2 + mat->data[18] * t3;
  t16 = mat->data[3] * t1 - mat->data[13] * t4 + mat->data[18] * t5;
  t17 = mat->data[3] * t2 - mat->data[8] * t4 + mat->data[18] * t6;

  inv->data[13] = mat->data[4] * t11 - mat->data[9] * t12 +
                  mat->data[14] * t13 - mat->data[24] * t14;
  inv->data[14] = -mat->data[4] * t15 + mat->data[9] * t16 -
                  mat->data[14] * t17 + mat->data[19] * t14;
  det = 1 / det;
  inv->data[0] = inv->data[0] * det;
  inv->data[1] = inv->data[1] * det;
  inv->data[2] = inv->data[2] * det;
  inv->data[3] = inv->data[3] * det;
  inv->data[4] = inv->data[4] * det;
  inv->data[5] = inv->data[5] * det;
  inv->data[6] = inv->data[6] * det;
  inv->data[7] = inv->data[7] * det;
  inv->data[8] = inv->data[8] * det;
  inv->data[9] = inv->data[9] * det;
  inv->data[10] = inv->data[10] * det;
  inv->data[11] = inv->data[11] * det;
  inv->data[12] = inv->data[12] * det;
  inv->data[13] = inv->data[13] * det;
  inv->data[14] = inv->data[14] * det;
  inv->data[15] = inv->data[15] * det;
  inv->data[16] = inv->data[16] * det;
  inv->data[17] = inv->data[17] * det;
  inv->data[18] = inv->data[18] * det;
  inv->data[19] = inv->data[19] * det;
  inv->data[20] = inv->data[20] * det;
  inv->data[21] = inv->data[21] * det;
  inv->data[22] = inv->data[22] * det;
  inv->data[23] = inv->data[23] * det;
  inv->data[24] = inv->data[24] * det;
}

void invert_6by6(MAT* mat, MAT* inv) {
  DTYPE det;
  DTYPE t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16,
      t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t29, t30, t31,
      t32, t33, t34, t35, t36, t37, t38, t39, t40, t41, t42, t43, t44, t45, t46,
      t47, t48, t49, t50;

  t1 = mat->data[28] * mat->data[35] - mat->data[34] * mat->data[29];
  t2 = mat->data[22] * mat->data[35] - mat->data[34] * mat->data[23];
  t3 = mat->data[22] * mat->data[29] - mat->data[28] * mat->data[23];
  t4 = mat->data[16] * mat->data[35] - mat->data[34] * mat->data[17];
  t5 = mat->data[16] * mat->data[29] - mat->data[28] * mat->data[17];
  t6 = mat->data[16] * mat->data[23] - mat->data[22] * mat->data[17];
  t7 = mat->data[10] * mat->data[35] - mat->data[34] * mat->data[11];
  t8 = mat->data[10] * mat->data[29] - mat->data[28] * mat->data[11];
  t9 = mat->data[10] * mat->data[23] - mat->data[22] * mat->data[11];
  t10 = mat->data[10] * mat->data[17] - mat->data[16] * mat->data[11];
  t11 = mat->data[4] * mat->data[35] - mat->data[34] * mat->data[5];
  t12 = mat->data[4] * mat->data[29] - mat->data[28] * mat->data[5];
  t13 = mat->data[4] * mat->data[23] - mat->data[22] * mat->data[5];
  t14 = mat->data[4] * mat->data[17] - mat->data[16] * mat->data[5];
  t15 = mat->data[4] * mat->data[11] - mat->data[10] * mat->data[5];

  t16 = mat->data[21] * t1 - mat->data[27] * t2 + mat->data[33] * t3;
  t17 = mat->data[15] * t1 - mat->data[27] * t4 + mat->data[33] * t5;
  t18 = mat->data[15] * t2 - mat->data[21] * t4 + mat->data[33] * t6;
  t19 = mat->data[15] * t3 - mat->data[21] * t5 + mat->data[27] * t6;
  t20 = mat->data[9] * t1 - mat->data[27] * t7 + mat->data[33] * t8;
  t21 = mat->data[9] * t2 - mat->data[21] * t7 + mat->data[33] * t9;
  t22 = mat->data[9] * t3 - mat->data[21] * t8 + mat->data[27] * t9;
  t23 = mat->data[9] * t4 - mat->data[15] * t7 + mat->data[33] * t10;
  t24 = mat->data[9] * t5 - mat->data[15] * t8 + mat->data[27] * t10;
  t25 = mat->data[9] * t6 - mat->data[15] * t9 + mat->data[21] * t10;
  t26 = mat->data[3] * t1 - mat->data[27] * t11 + mat->data[33] * t12;
  t27 = mat->data[3] * t2 - mat->data[21] * t11 + mat->data[33] * t13;
  t28 = mat->data[3] * t3 - mat->data[21] * t12 + mat->data[27] * t13;
  t29 = mat->data[3] * t4 - mat->data[15] * t11 + mat->data[33] * t14;
  t30 = mat->data[3] * t5 - mat->data[15] * t12 + mat->data[27] * t14;
  t31 = mat->data[3] * t6 - mat->data[15] * t13 + mat->data[21] * t14;
  t32 = mat->data[3] * t7 - mat->data[9] * t11 + mat->data[33] * t15;
  t33 = mat->data[3] * t8 - mat->data[9] * t12 + mat->data[27] * t15;
  t34 = mat->data[3] * t9 - mat->data[9] * t13 + mat->data[21] * t15;
  t35 = mat->data[3] * t10 - mat->data[9] * t14 + mat->data[15] * t15;

  t36 = mat->data[14] * t16 - mat->data[20] * t17 + mat->data[26] * t18 -
        mat->data[32] * t19;
  t37 = mat->data[8] * t16 - mat->data[20] * t20 + mat->data[26] * t21 -
        mat->data[32] * t22;
  t38 = mat->data[8] * t17 - mat->data[14] * t20 + mat->data[26] * t23 -
        mat->data[32] * t24;
  t39 = mat->data[8] * t18 - mat->data[14] * t21 + mat->data[20] * t23 -
        mat->data[32] * t25;
  t40 = mat->data[8] * t19 - mat->data[14] * t22 + mat->data[20] * t24 -
        mat->data[26] * t25;
  t41 = mat->data[2] * t16 - mat->data[20] * t26 + mat->data[26] * t27 -
        mat->data[32] * t28;
  t42 = mat->data[2] * t17 - mat->data[14] * t26 + mat->data[26] * t29 -
        mat->data[32] * t30;
  t43 = mat->data[2] * t18 - mat->data[14] * t27 + mat->data[20] * t29 -
        mat->data[32] * t31;
  t44 = mat->data[2] * t19 - mat->data[14] * t28 + mat->data[20] * t30 -
        mat->data[26] * t31;

  inv->data[0] = mat->data[7] * t36 - mat->data[13] * t37 +
                 mat->data[19] * t38 - mat->data[25] * t39 +
                 mat->data[31] * t40;
  inv->data[6] = -mat->data[6] * t36 + mat->data[12] * t37 -
                 mat->data[18] * t38 + mat->data[24] * t39 -
                 mat->data[30] * t40;
  inv->data[1] = -mat->data[1] * t36 + mat->data[13] * t41 -
                 mat->data[19] * t42 + mat->data[25] * t43 -
                 mat->data[31] * t44;
  inv->data[7] = mat->data[0] * t36 - mat->data[12] * t41 +
                 mat->data[18] * t42 - mat->data[24] * t43 +
                 mat->data[30] * t44;

  t45 = mat->data[2] * t20 - mat->data[8] * t26 + mat->data[26] * t32 -
        mat->data[32] * t33;
  t46 = mat->data[2] * t21 - mat->data[8] * t27 + mat->data[20] * t32 -
        mat->data[32] * t34;
  t47 = mat->data[2] * t22 - mat->data[8] * t28 + mat->data[20] * t33 -
        mat->data[26] * t34;
  t48 = mat->data[2] * t23 - mat->data[8] * t29 + mat->data[14] * t32 -
        mat->data[32] * t35;
  t49 = mat->data[2] * t24 - mat->data[8] * t30 + mat->data[14] * t33 -
        mat->data[26] * t35;

  inv->data[2] = mat->data[1] * t37 - mat->data[7] * t41 + mat->data[19] * t45 -
                 mat->data[25] * t46 + mat->data[31] * t47;
  inv->data[8] = -mat->data[0] * t37 + mat->data[6] * t41 -
                 mat->data[18] * t45 + mat->data[24] * t46 -
                 mat->data[30] * t47;
  inv->data[3] = -mat->data[1] * t38 + mat->data[7] * t42 -
                 mat->data[13] * t45 + mat->data[25] * t48 -
                 mat->data[31] * t49;
  inv->data[9] = mat->data[0] * t38 - mat->data[6] * t42 + mat->data[12] * t45 -
                 mat->data[24] * t48 + mat->data[30] * t49;

  t50 = mat->data[2] * t25 - mat->data[8] * t31 + mat->data[14] * t34 -
        mat->data[20] * t35;

  inv->data[4] = mat->data[1] * t39 - mat->data[7] * t43 + mat->data[13] * t46 -
                 mat->data[19] * t48 + mat->data[31] * t50;
  inv->data[10] = -mat->data[0] * t39 + mat->data[6] * t43 -
                  mat->data[12] * t46 + mat->data[18] * t48 -
                  mat->data[30] * t50;
  inv->data[5] = -mat->data[1] * t40 + mat->data[7] * t44 -
                 mat->data[13] * t47 + mat->data[19] * t49 -
                 mat->data[25] * t50;
  inv->data[11] = mat->data[0] * t40 - mat->data[6] * t44 +
                  mat->data[12] * t47 - mat->data[18] * t49 +
                  mat->data[24] * t50;

  det = mat->data[0] * inv->data[0] + mat->data[6] * inv->data[1] +
        mat->data[12] * inv->data[2] + mat->data[18] * inv->data[3] +
        mat->data[24] * inv->data[4] + mat->data[30] * inv->data[5];

#if DEBUG
  printf("det : %lf\n", det);
#endif
  if (det > -FZERO && det < FZERO) ASSERT(DET_FAIL)
  t36 = mat->data[13] * t16 - mat->data[19] * t17 + mat->data[25] * t18 -
        mat->data[31] * t19;
  t37 = mat->data[7] * t16 - mat->data[19] * t20 + mat->data[25] * t21 -
        mat->data[31] * t22;
  t38 = mat->data[7] * t17 - mat->data[13] * t20 + mat->data[25] * t23 -
        mat->data[31] * t24;
  t39 = mat->data[7] * t18 - mat->data[13] * t21 + mat->data[19] * t23 -
        mat->data[31] * t25;
  t40 = mat->data[7] * t19 - mat->data[13] * t22 + mat->data[19] * t24 -
        mat->data[25] * t25;
  t41 = mat->data[1] * t16 - mat->data[19] * t26 + mat->data[25] * t27 -
        mat->data[31] * t28;
  t42 = mat->data[1] * t17 - mat->data[13] * t26 + mat->data[25] * t29 -
        mat->data[31] * t30;
  t43 = mat->data[1] * t18 - mat->data[13] * t27 + mat->data[19] * t29 -
        mat->data[31] * t31;
  t44 = mat->data[1] * t19 - mat->data[13] * t28 + mat->data[19] * t30 -
        mat->data[25] * t31;
  t45 = mat->data[1] * t20 - mat->data[7] * t26 + mat->data[25] * t32 -
        mat->data[31] * t33;
  t46 = mat->data[1] * t21 - mat->data[7] * t27 + mat->data[19] * t32 -
        mat->data[31] * t34;
  t47 = mat->data[1] * t22 - mat->data[7] * t28 + mat->data[19] * t33 -
        mat->data[25] * t34;
  t48 = mat->data[1] * t23 - mat->data[7] * t29 + mat->data[13] * t32 -
        mat->data[31] * t35;
  t49 = mat->data[1] * t24 - mat->data[7] * t30 + mat->data[13] * t33 -
        mat->data[25] * t35;
  t50 = mat->data[1] * t25 - mat->data[7] * t31 + mat->data[13] * t34 -
        mat->data[19] * t35;

  inv->data[12] = mat->data[6] * t36 - mat->data[12] * t37 +
                  mat->data[18] * t38 - mat->data[24] * t39 +
                  mat->data[30] * t40;
  inv->data[13] = -mat->data[0] * t36 + mat->data[12] * t41 -
                  mat->data[18] * t42 + mat->data[24] * t43 -
                  mat->data[30] * t44;
  inv->data[14] = mat->data[0] * t37 - mat->data[6] * t41 +
                  mat->data[18] * t45 - mat->data[24] * t46 +
                  mat->data[30] * t47;
  inv->data[15] = -mat->data[0] * t38 + mat->data[6] * t42 -
                  mat->data[12] * t45 + mat->data[24] * t48 -
                  mat->data[30] * t49;
  inv->data[16] = mat->data[0] * t39 - mat->data[6] * t43 +
                  mat->data[12] * t46 - mat->data[18] * t48 +
                  mat->data[30] * t50;
  inv->data[17] = -mat->data[0] * t40 + mat->data[6] * t44 -
                  mat->data[12] * t47 + mat->data[18] * t49 -
                  mat->data[24] * t50;

  t1 = mat->data[18] * mat->data[25] - mat->data[24] * mat->data[19];
  t2 = mat->data[12] * mat->data[25] - mat->data[24] * mat->data[13];
  t3 = mat->data[12] * mat->data[19] - mat->data[18] * mat->data[13];
  t4 = mat->data[6] * mat->data[25] - mat->data[24] * mat->data[7];
  t5 = mat->data[6] * mat->data[19] - mat->data[18] * mat->data[7];
  t6 = mat->data[6] * mat->data[13] - mat->data[12] * mat->data[7];
  t7 = mat->data[0] * mat->data[25] - mat->data[24] * mat->data[1];
  t8 = mat->data[0] * mat->data[19] - mat->data[18] * mat->data[1];
  t9 = mat->data[0] * mat->data[13] - mat->data[12] * mat->data[1];
  t10 = mat->data[0] * mat->data[7] - mat->data[6] * mat->data[1];
  t11 = mat->data[18] * mat->data[31] - mat->data[30] * mat->data[19];
  t12 = mat->data[12] * mat->data[31] - mat->data[30] * mat->data[13];
  t13 = mat->data[6] * mat->data[31] - mat->data[30] * mat->data[7];
  t14 = mat->data[0] * mat->data[31] - mat->data[30] * mat->data[1];
  t15 = mat->data[24] * mat->data[31] - mat->data[30] * mat->data[25];

  t16 = mat->data[20] * t15 - mat->data[26] * t11 + mat->data[32] * t1;
  t17 = mat->data[14] * t15 - mat->data[26] * t12 + mat->data[32] * t2;
  t18 = mat->data[14] * t11 - mat->data[20] * t12 + mat->data[32] * t3;
  t19 = mat->data[14] * t1 - mat->data[20] * t2 + mat->data[26] * t3;
  t20 = mat->data[8] * t15 - mat->data[26] * t13 + mat->data[32] * t4;
  t21 = mat->data[8] * t11 - mat->data[20] * t13 + mat->data[32] * t5;
  t22 = mat->data[8] * t1 - mat->data[20] * t4 + mat->data[26] * t5;
  t23 = mat->data[8] * t12 - mat->data[14] * t13 + mat->data[32] * t6;
  t24 = mat->data[8] * t2 - mat->data[14] * t4 + mat->data[26] * t6;
  t25 = mat->data[8] * t3 - mat->data[14] * t5 + mat->data[20] * t6;
  t26 = mat->data[2] * t15 - mat->data[26] * t14 + mat->data[32] * t7;
  t27 = mat->data[2] * t11 - mat->data[20] * t14 + mat->data[32] * t8;
  t28 = mat->data[2] * t1 - mat->data[20] * t7 + mat->data[26] * t8;
  t29 = mat->data[2] * t12 - mat->data[14] * t14 + mat->data[32] * t9;
  t30 = mat->data[2] * t2 - mat->data[14] * t7 + mat->data[26] * t9;
  t31 = mat->data[2] * t3 - mat->data[14] * t8 + mat->data[20] * t9;
  t32 = mat->data[2] * t13 - mat->data[8] * t14 + mat->data[32] * t10;
  t33 = mat->data[2] * t4 - mat->data[8] * t7 + mat->data[26] * t10;
  t34 = mat->data[2] * t5 - mat->data[8] * t8 + mat->data[20] * t10;
  t35 = mat->data[2] * t6 - mat->data[8] * t9 + mat->data[14] * t10;

  t36 = mat->data[15] * t16 - mat->data[21] * t17 + mat->data[27] * t18 -
        mat->data[33] * t19;
  t37 = mat->data[9] * t16 - mat->data[21] * t20 + mat->data[27] * t21 -
        mat->data[33] * t22;
  t38 = mat->data[9] * t17 - mat->data[15] * t20 + mat->data[27] * t23 -
        mat->data[33] * t24;
  t39 = mat->data[9] * t18 - mat->data[15] * t21 + mat->data[21] * t23 -
        mat->data[33] * t25;
  t40 = mat->data[9] * t19 - mat->data[15] * t22 + mat->data[21] * t24 -
        mat->data[27] * t25;
  t41 = mat->data[3] * t16 - mat->data[21] * t26 + mat->data[27] * t27 -
        mat->data[33] * t28;
  t42 = mat->data[3] * t17 - mat->data[15] * t26 + mat->data[27] * t29 -
        mat->data[33] * t30;
  t43 = mat->data[3] * t18 - mat->data[15] * t27 + mat->data[21] * t29 -
        mat->data[33] * t31;
  t44 = mat->data[3] * t19 - mat->data[15] * t28 + mat->data[21] * t30 -
        mat->data[27] * t31;

  inv->data[24] = -mat->data[11] * t36 + mat->data[17] * t37 -
                  mat->data[23] * t38 + mat->data[29] * t39 -
                  mat->data[35] * t40;
  inv->data[30] = mat->data[10] * t36 - mat->data[16] * t37 +
                  mat->data[22] * t38 - mat->data[28] * t39 +
                  mat->data[34] * t40;
  inv->data[25] = mat->data[5] * t36 - mat->data[17] * t41 +
                  mat->data[23] * t42 - mat->data[29] * t43 +
                  mat->data[35] * t44;
  inv->data[31] = -mat->data[4] * t36 + mat->data[16] * t41 -
                  mat->data[22] * t42 + mat->data[28] * t43 -
                  mat->data[34] * t44;

  t45 = mat->data[3] * t20 - mat->data[9] * t26 + mat->data[27] * t32 -
        mat->data[33] * t33;
  t46 = mat->data[3] * t21 - mat->data[9] * t27 + mat->data[21] * t32 -
        mat->data[33] * t34;
  t47 = mat->data[3] * t22 - mat->data[9] * t28 + mat->data[21] * t33 -
        mat->data[27] * t34;
  t48 = mat->data[3] * t23 - mat->data[9] * t29 + mat->data[15] * t32 -
        mat->data[33] * t35;
  t49 = mat->data[3] * t24 - mat->data[9] * t30 + mat->data[15] * t33 -
        mat->data[27] * t35;

  inv->data[26] = -mat->data[5] * t37 + mat->data[11] * t41 -
                  mat->data[23] * t45 + mat->data[29] * t46 -
                  mat->data[35] * t47;
  inv->data[32] = mat->data[4] * t37 - mat->data[10] * t41 +
                  mat->data[22] * t45 - mat->data[28] * t46 +
                  mat->data[34] * t47;
  inv->data[27] = mat->data[5] * t38 - mat->data[11] * t42 +
                  mat->data[17] * t45 - mat->data[29] * t48 +
                  mat->data[35] * t49;
  inv->data[33] = -mat->data[4] * t38 + mat->data[10] * t42 -
                  mat->data[16] * t45 + mat->data[28] * t48 -
                  mat->data[34] * t49;

  t50 = mat->data[3] * t25 - mat->data[9] * t31 + mat->data[15] * t34 -
        mat->data[21] * t35;

  inv->data[28] = -mat->data[5] * t39 + mat->data[11] * t43 -
                  mat->data[17] * t46 + mat->data[23] * t48 -
                  mat->data[35] * t50;
  inv->data[34] = mat->data[4] * t39 - mat->data[10] * t43 +
                  mat->data[16] * t46 - mat->data[22] * t48 +
                  mat->data[34] * t50;
  inv->data[29] = mat->data[5] * t40 - mat->data[11] * t44 +
                  mat->data[17] * t47 - mat->data[23] * t49 +
                  mat->data[29] * t50;
  inv->data[35] = -mat->data[4] * t40 + mat->data[10] * t44 -
                  mat->data[16] * t47 + mat->data[22] * t49 -
                  mat->data[28] * t50;

  t36 = mat->data[16] * t16 - mat->data[22] * t17 + mat->data[28] * t18 -
        mat->data[34] * t19;
  t37 = mat->data[10] * t16 - mat->data[22] * t20 + mat->data[28] * t21 -
        mat->data[34] * t22;
  t38 = mat->data[10] * t17 - mat->data[16] * t20 + mat->data[28] * t23 -
        mat->data[34] * t24;
  t39 = mat->data[10] * t18 - mat->data[16] * t21 + mat->data[22] * t23 -
        mat->data[34] * t25;
  t40 = mat->data[10] * t19 - mat->data[16] * t22 + mat->data[22] * t24 -
        mat->data[28] * t25;
  t41 = mat->data[4] * t16 - mat->data[22] * t26 + mat->data[28] * t27 -
        mat->data[34] * t28;
  t42 = mat->data[4] * t17 - mat->data[16] * t26 + mat->data[28] * t29 -
        mat->data[34] * t30;
  t43 = mat->data[4] * t18 - mat->data[16] * t27 + mat->data[22] * t29 -
        mat->data[34] * t31;
  t44 = mat->data[4] * t19 - mat->data[16] * t28 + mat->data[22] * t30 -
        mat->data[28] * t31;
  t45 = mat->data[4] * t20 - mat->data[10] * t26 + mat->data[28] * t32 -
        mat->data[34] * t33;
  t46 = mat->data[4] * t21 - mat->data[10] * t27 + mat->data[22] * t32 -
        mat->data[34] * t34;
  t47 = mat->data[4] * t22 - mat->data[10] * t28 + mat->data[22] * t33 -
        mat->data[28] * t34;
  t48 = mat->data[4] * t23 - mat->data[10] * t29 + mat->data[16] * t32 -
        mat->data[34] * t35;
  t49 = mat->data[4] * t24 - mat->data[10] * t30 + mat->data[16] * t33 -
        mat->data[28] * t35;
  t50 = mat->data[4] * t25 - mat->data[10] * t31 + mat->data[16] * t34 -
        mat->data[22] * t35;

  inv->data[18] = mat->data[11] * t36 - mat->data[17] * t37 +
                  mat->data[23] * t38 - mat->data[29] * t39 +
                  mat->data[35] * t40;
  inv->data[19] = -mat->data[5] * t36 + mat->data[17] * t41 -
                  mat->data[23] * t42 + mat->data[29] * t43 -
                  mat->data[35] * t44;
  inv->data[20] = mat->data[5] * t37 - mat->data[11] * t41 +
                  mat->data[23] * t45 - mat->data[29] * t46 +
                  mat->data[35] * t47;
  inv->data[21] = -mat->data[5] * t38 + mat->data[11] * t42 -
                  mat->data[17] * t45 + mat->data[29] * t48 -
                  mat->data[35] * t49;
  inv->data[22] = mat->data[5] * t39 - mat->data[11] * t43 +
                  mat->data[17] * t46 - mat->data[23] * t48 +
                  mat->data[35] * t50;
  inv->data[23] = -mat->data[5] * t40 + mat->data[11] * t44 -
                  mat->data[17] * t47 + mat->data[23] * t49 -
                  mat->data[29] * t50;

  det = 1 / det;

  inv->data[0] = inv->data[0] * det;
  inv->data[1] = inv->data[1] * det;
  inv->data[2] = inv->data[2] * det;
  inv->data[3] = inv->data[3] * det;
  inv->data[4] = inv->data[4] * det;
  inv->data[5] = inv->data[5] * det;
  inv->data[6] = inv->data[6] * det;
  inv->data[7] = inv->data[7] * det;
  inv->data[8] = inv->data[8] * det;
  inv->data[9] = inv->data[9] * det;
  inv->data[10] = inv->data[10] * det;
  inv->data[11] = inv->data[11] * det;
  inv->data[12] = inv->data[12] * det;
  inv->data[13] = inv->data[13] * det;
  inv->data[14] = inv->data[14] * det;
  inv->data[15] = inv->data[15] * det;
  inv->data[16] = inv->data[16] * det;
  inv->data[17] = inv->data[17] * det;
  inv->data[18] = inv->data[18] * det;
  inv->data[19] = inv->data[19] * det;
  inv->data[20] = inv->data[20] * det;
  inv->data[21] = inv->data[21] * det;
  inv->data[22] = inv->data[22] * det;
  inv->data[23] = inv->data[23] * det;
  inv->data[24] = inv->data[24] * det;
  inv->data[25] = inv->data[25] * det;
  inv->data[26] = inv->data[26] * det;
  inv->data[27] = inv->data[27] * det;
  inv->data[28] = inv->data[28] * det;
  inv->data[29] = inv->data[29] * det;
  inv->data[30] = inv->data[30] * det;
  inv->data[31] = inv->data[31] * det;
  inv->data[32] = inv->data[32] * det;
  inv->data[33] = inv->data[33] * det;
  inv->data[34] = inv->data[34] * det;
  inv->data[35] = inv->data[35] * det;
}

/**** get determinant of matrix ****/
DTYPE det_2by2(MAT* mat) {
  DTYPE det;
#if DEBUF
  printf("%s\n", __func__);
#endif
  det = mat->data[0] * mat->data[3] - mat->data[1] * mat->data[2];
  return det;
}

DTYPE det_3by3(MAT* mat) {
  DTYPE det;
  DTYPE t0, t1, t2;

  t0 = mat->data[4] * mat->data[8] - mat->data[7] * mat->data[5];
  t1 = mat->data[7] * mat->data[2] - mat->data[1] * mat->data[8];
  t2 = mat->data[1] * mat->data[5] - mat->data[4] * mat->data[2];

  det = mat->data[0] * t0 + mat->data[3] * t1 + mat->data[6] * t2;

  return det;
}

DTYPE det_4by4(MAT* mat) {
  DTYPE det;
  DTYPE t1, t2, t3, t4, t5;
  DTYPE m0, m1, m2, m3;
#if DEBUF
  printf("%s\n", __func__);
#endif
  t1 = mat->data[10] * mat->data[15] - mat->data[14] * mat->data[11];
  t2 = mat->data[6] * mat->data[15] - mat->data[14] * mat->data[7];
  t3 = mat->data[6] * mat->data[11] - mat->data[10] * mat->data[7];

  m0 = mat->data[5] * t1 - mat->data[9] * t2 + mat->data[13] * t3;

  t4 = mat->data[2] * mat->data[15] - mat->data[14] * mat->data[3];
  t5 = mat->data[2] * mat->data[11] - mat->data[10] * mat->data[3];

  m1 = mat->data[9] * t4 - mat->data[1] * t1 - mat->data[13] * t5;

  t1 = mat->data[2] * mat->data[7] - mat->data[6] * mat->data[3];

  m2 = mat->data[1] * t2 - mat->data[5] * t4 + mat->data[13] * t1;
  m3 = mat->data[5] * t5 - mat->data[1] * t3 - mat->data[9] * t1;

  det = mat->data[0] * m0 + mat->data[4] * m1 + mat->data[8] * m2 +
        mat->data[12] * m3;
  return det;
}
DTYPE det_5by5(MAT* mat) {
  DTYPE det;
  DTYPE t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16,
      t17, t18, t19, t20;
  DTYPE m0, m1, m2, m3, m4;
#if DEBUF
  printf("%s\n", __func__);
#endif
  t1 = mat->data[18] * mat->data[24] - mat->data[23] * mat->data[19];
  t2 = mat->data[13] * mat->data[24] - mat->data[23] * mat->data[14];
  t3 = mat->data[13] * mat->data[19] - mat->data[18] * mat->data[14];
  t4 = mat->data[8] * mat->data[24] - mat->data[23] * mat->data[9];
  t5 = mat->data[8] * mat->data[19] - mat->data[18] * mat->data[9];
  t6 = mat->data[8] * mat->data[14] - mat->data[13] * mat->data[9];
  t7 = mat->data[3] * mat->data[24] - mat->data[23] * mat->data[4];
  t8 = mat->data[3] * mat->data[19] - mat->data[18] * mat->data[4];
  t9 = mat->data[3] * mat->data[14] - mat->data[13] * mat->data[4];
  t10 = mat->data[3] * mat->data[9] - mat->data[8] * mat->data[4];

  t11 = mat->data[12] * t1 - mat->data[17] * t2 + mat->data[22] * t3;
  t12 = mat->data[7] * t1 - mat->data[17] * t4 + mat->data[22] * t5;
  t13 = mat->data[7] * t2 - mat->data[12] * t4 + mat->data[22] * t6;
  t14 = mat->data[7] * t3 - mat->data[12] * t5 + mat->data[17] * t6;
  t15 = mat->data[2] * t1 - mat->data[17] * t7 + mat->data[22] * t8;
  t16 = mat->data[2] * t2 - mat->data[12] * t7 + mat->data[22] * t9;
  t17 = mat->data[2] * t3 - mat->data[12] * t8 + mat->data[17] * t9;

  m0 = mat->data[6] * t11 - mat->data[11] * t12 + mat->data[16] * t13 -
       mat->data[21] * t14;
  mat->data[15] * t13 + mat->data[20] * t14;
  m1 = -mat->data[1] * t11 + mat->data[11] * t15 - mat->data[16] * t16 +
       mat->data[21] * t17;
  mat->data[15] * t16 - mat->data[20] * t17;

  t18 = mat->data[2] * t4 - mat->data[7] * t7 + mat->data[22] * t10;
  t19 = mat->data[2] * t5 - mat->data[7] * t8 + mat->data[17] * t10;
  t20 = mat->data[2] * t6 - mat->data[7] * t9 + mat->data[12] * t10;

  m2 = mat->data[1] * t12 - mat->data[6] * t15 + mat->data[16] * t18 -
       mat->data[21] * t19;
  mat->data[15] * t18 + mat->data[20] * t19;
  m3 = -mat->data[1] * t13 + mat->data[6] * t16 - mat->data[11] * t18 +
       mat->data[21] * t20;
  mat->data[20] * t20;
  m4 = mat->data[1] * t14 - mat->data[6] * t17 + mat->data[11] * t19 -
       mat->data[16] * t20;
  mat->data[10] * t19 + mat->data[15] * t20;

  det = mat->data[0] * m0 + mat->data[5] * m1 + mat->data[10] * m2 +
        mat->data[15] * m3 + mat->data[20] * m4;

  return det;
}
DTYPE det_6by6(MAT* mat) {
  DTYPE det;
  DTYPE t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16,
      t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t29, t30, t31,
      t32, t33, t34, t35, t36, t37, t38, t39, t40, t41, t42, t43, t44, t45, t46,
      t47, t48, t49, t50;
  DTYPE m0, m1, m2, m3, m4, m5;
#if DEBUG
  printf("%s\n", __func__);
#endif

  t1 = mat->data[28] * mat->data[35] - mat->data[34] * mat->data[29];
  t2 = mat->data[22] * mat->data[35] - mat->data[34] * mat->data[23];
  t3 = mat->data[22] * mat->data[29] - mat->data[28] * mat->data[23];
  t4 = mat->data[16] * mat->data[35] - mat->data[34] * mat->data[17];
  t5 = mat->data[16] * mat->data[29] - mat->data[28] * mat->data[17];
  t6 = mat->data[16] * mat->data[23] - mat->data[22] * mat->data[17];
  t7 = mat->data[10] * mat->data[35] - mat->data[34] * mat->data[11];
  t8 = mat->data[10] * mat->data[29] - mat->data[28] * mat->data[11];
  t9 = mat->data[10] * mat->data[23] - mat->data[22] * mat->data[11];
  t10 = mat->data[10] * mat->data[17] - mat->data[16] * mat->data[11];
  t11 = mat->data[4] * mat->data[35] - mat->data[34] * mat->data[5];
  t12 = mat->data[4] * mat->data[29] - mat->data[28] * mat->data[5];
  t13 = mat->data[4] * mat->data[23] - mat->data[22] * mat->data[5];
  t14 = mat->data[4] * mat->data[17] - mat->data[16] * mat->data[5];
  t15 = mat->data[4] * mat->data[11] - mat->data[10] * mat->data[5];

  t16 = mat->data[21] * t1 - mat->data[27] * t2 + mat->data[33] * t3;
  t17 = mat->data[15] * t1 - mat->data[27] * t4 + mat->data[33] * t5;
  t18 = mat->data[15] * t2 - mat->data[21] * t4 + mat->data[33] * t6;
  t19 = mat->data[15] * t3 - mat->data[21] * t5 + mat->data[27] * t6;
  t20 = mat->data[9] * t1 - mat->data[27] * t7 + mat->data[33] * t8;
  t21 = mat->data[9] * t2 - mat->data[21] * t7 + mat->data[33] * t9;
  t22 = mat->data[9] * t3 - mat->data[21] * t8 + mat->data[27] * t9;
  t23 = mat->data[9] * t4 - mat->data[15] * t7 + mat->data[33] * t10;
  t24 = mat->data[9] * t5 - mat->data[15] * t8 + mat->data[27] * t10;
  t25 = mat->data[9] * t6 - mat->data[15] * t9 + mat->data[21] * t10;
  t26 = mat->data[3] * t1 - mat->data[27] * t11 + mat->data[33] * t12;
  t27 = mat->data[3] * t2 - mat->data[21] * t11 + mat->data[33] * t13;
  t28 = mat->data[3] * t3 - mat->data[21] * t12 + mat->data[27] * t13;
  t29 = mat->data[3] * t4 - mat->data[15] * t11 + mat->data[33] * t14;
  t30 = mat->data[3] * t5 - mat->data[15] * t12 + mat->data[27] * t14;
  t31 = mat->data[3] * t6 - mat->data[15] * t13 + mat->data[21] * t14;
  t32 = mat->data[3] * t7 - mat->data[9] * t11 + mat->data[33] * t15;
  t33 = mat->data[3] * t8 - mat->data[9] * t12 + mat->data[27] * t15;
  t34 = mat->data[3] * t9 - mat->data[9] * t13 + mat->data[21] * t15;
  t35 = mat->data[3] * t10 - mat->data[9] * t14 + mat->data[15] * t15;

  t36 = mat->data[14] * t16 - mat->data[20] * t17 + mat->data[26] * t18 -
        mat->data[32] * t19;
  t37 = mat->data[8] * t16 - mat->data[20] * t20 + mat->data[26] * t21 -
        mat->data[32] * t22;
  t38 = mat->data[8] * t17 - mat->data[14] * t20 + mat->data[26] * t23 -
        mat->data[32] * t24;
  t39 = mat->data[8] * t18 - mat->data[14] * t21 + mat->data[20] * t23 -
        mat->data[32] * t25;
  t40 = mat->data[8] * t19 - mat->data[14] * t22 + mat->data[20] * t24 -
        mat->data[26] * t25;
  t41 = mat->data[2] * t16 - mat->data[20] * t26 + mat->data[26] * t27 -
        mat->data[32] * t28;
  t42 = mat->data[2] * t17 - mat->data[14] * t26 + mat->data[26] * t29 -
        mat->data[32] * t30;
  t43 = mat->data[2] * t18 - mat->data[14] * t27 + mat->data[20] * t29 -
        mat->data[32] * t31;
  t44 = mat->data[2] * t19 - mat->data[14] * t28 + mat->data[20] * t30 -
        mat->data[26] * t31;

  m0 = mat->data[7] * t36 - mat->data[13] * t37 + mat->data[19] * t38 -
       mat->data[25] * t39 + mat->data[31] * t40;
  m1 = -mat->data[1] * t36 + mat->data[13] * t41 - mat->data[19] * t42 +
       mat->data[25] * t43 - mat->data[31] * t44;

  t45 = mat->data[2] * t20 - mat->data[8] * t26 + mat->data[26] * t32 -
        mat->data[32] * t33;
  t46 = mat->data[2] * t21 - mat->data[8] * t27 + mat->data[20] * t32 -
        mat->data[32] * t34;
  t47 = mat->data[2] * t22 - mat->data[8] * t28 + mat->data[20] * t33 -
        mat->data[26] * t34;
  t48 = mat->data[2] * t23 - mat->data[8] * t29 + mat->data[14] * t32 -
        mat->data[32] * t35;
  t49 = mat->data[2] * t24 - mat->data[8] * t30 + mat->data[14] * t33 -
        mat->data[26] * t35;

  m2 = mat->data[1] * t37 - mat->data[7] * t41 + mat->data[19] * t45 -
       mat->data[25] * t46 + mat->data[31] * t47;
  m3 = -mat->data[1] * t38 + mat->data[7] * t42 - mat->data[13] * t45 +
       mat->data[25] * t48 - mat->data[31] * t49;

  t50 = mat->data[2] * t25 - mat->data[8] * t31 + mat->data[14] * t34 -
        mat->data[20] * t35;

  m4 = mat->data[1] * t39 - mat->data[7] * t43 + mat->data[13] * t46 -
       mat->data[19] * t48 + mat->data[31] * t50;
  m5 = -mat->data[1] * t40 + mat->data[7] * t44 - mat->data[13] * t47 +
       mat->data[19] * t49 - mat->data[25] * t50;

  det = mat->data[0] * m0 + mat->data[6] * m1 + mat->data[12] * m2 +
        mat->data[18] * m3 + mat->data[24] * m4 + mat->data[30] * m5;

  return det;
}

/**** LAPACK ****/
void invert_nbyn(MAT*mat,MAT*inv){
  UINT n;
  UINT lwork,info;
  UINT *idx; 
  DTYPE* work;

#if DEBUG
printf("%s\n",__func__);
#endif
  n = mat->d0;
  lwork = n*n;
  
#if USE_LAPACK
  work = iip_malloc(sizeof(DTYPE)*n);
  idx = iip_malloc(sizeof(UINT)*n);
  #if NTYPE == 0 
  LAPACK_sgetrf( &n,&n,mat->data,&n,idx,&info);
  LAPACK_sgetri(&n,mat->data,&n,idx,work,&lwork,&info);
  #elif NTYPE == 1
  LAPACK_dgetrf( &n,&n,mat->data,&n,idx,&info);
  LAPACK_dgetri(&n,mat->data,&n,idx,work,&lwork,&info);
  #endif
  iip_free(idx);
  iip_free(work);

#elif USE_MKL
  idx = iip_malloc(sizeof(UINT)*n);
  #if NTYPE == 0 
  LAPACK_sgetrf(LAPACK_COL_MAJOR, &n,&n,mat->data,&n,idx;
  LAPACK_sgetri(LAPACK_COL_MAJOR,&n,mat->data,&n,idx);
  #elif NTYPE == 1
  LAPACK_dgetrf(LAPACK_COL_MAJOR, &n,&n,mat->data,&n,idx);
  LAPACK_dgetri(LAPACK_COL_MAJOR,&n,mat->data,&n,idx);
  #endif
  iip_free(idx);

#else
  printf("ERROR : You need to use 'LAPACK' or 'Intel MKL'\n");
  return;
#endif
}

DTYPE det_nbyn(MAT*mat){
  ITER i;
  DTYPE det;
  UINT info;
  UINT n;
  UINT * idx;
#if DEBUG
printf("%s\n",__func__);
#endif
  n = mat->d0;

#if USE_LAPACK
  idx = iip_malloc(sizeof(UINT)*n);

  #if NTYPE == 0
  LAPACK_sgetrf( &n,&n,mat->data,&n,idx,&info);
  #elif NTYPE == 1
  LAPACK_dgetrf( &n,&n,mat->data,&n,idx,&info);
  #endif

#elif USE_MKL
  idx = iip_malloc(sizeof(UINT)*n);

  #if NTYPE == 0
  LAPACK_sgetrf( LAPACK_COL_MAJOR,&n,&n,mat->data,&n,idx);
  #elif NTYPE == 1
  LAPACK_dgetrf( LAPACK_COL_MAJOR,&n,&n,mat->data,&n,idx);
  #endif

#else

  printf("ERROR : You need to use 'LAPACK' or 'Intel MKL'\n");
  return -1;
#endif
  det = 1; 
  for(i=0;i<mat->d0;i++){
    if(i+1 != idx[i])
      det *= - mat->data[i*n + i];
    else
      det *= mat->data[i*n + i];
  }
 iip_free(idx);
 return det;
}

/**** LU ****/
void lu_decomposition(DTYPE* X, UINT n, UINT* idx) {
  ITER i, j, k, imax, t;
  DTYPE big, dum, sum, temp;
  DTYPE* v;
  v = iip_malloc(sizeof(DTYPE) * (n));
  
  for(i=0;i<n;i++)
    idx[i] = i;

  for (j = 0; j < n; j++) {
    big = 0.0;
    for (i = 0; i < n; i++)
      if ((temp = fabs(X[i + n * j])) > big) big = temp;
    v[j] = 1. / big;
  }

  for (j = 0; j < n; j++) {
    for (i = 0; i < j; i++) {
      sum = X[i + n * j];
      for (k = 0; k < i; k++) sum -= X[i + n * k] * X[k + n * j];
      X[i + n * j] = sum;
    }
    big = 0.0;
    for (i = j; i < n; i++) {
      sum = X[i + n * j];
      for (k = 0; k < j; k++) sum -= X[i + n * k] * X[k + n * j];
      X[i + n * j] = sum;
      if ((dum = v[i] * fabs(sum)) >= big) {
        big = dum;
        imax = i;
      }
    }

    if (j != imax) {
      for (k = 0; k < n; k++) {
        dum = X[imax + n * k];
        X[imax + n * k] = X[j + n * k];
        X[j + n * k] = dum;
      }
      v[imax] = v[j];
    }
    t = idx[j];
    idx[j] = imax;
    idx[imax] = t;
    printf("idx[%d] = %d\nidx[%d] = %d\n",j,imax,imax,t);
    if (X[j + n * j] == 0.0) X[j + n * j] = FZERO;
    if ((j + 1) != n) {
      dum = 1.0 / X[j + n * j];
      for (i = j+1; i < n; i++) X[i + n * j] *= dum;
    }
  }
  

 iip_free(v);
}

void lu_backwardsubstitution(DTYPE* lu, UINT n, UINT* idx, DTYPE* B) {
  ITER i, ii, ip, j;
  DTYPE sum;

  ii = -1;
  for (i = 0; i < n; i++) {
    ip = idx[i];
    sum = B[ip];
    B[ip] = B[i];
    if (ii >=0)
      for (j = ii; j <= i - 1; j++) sum -= lu[i + n*j] * B[j];
    else if (sum)
      ii = i;
    B[i] = sum;
  }
  for (i = n - 1; i >= 0; i--) {
    sum = B[i];
    for (j = i+1; j < n; j++) sum -= lu[i + n * j] * B[j];
    B[i] = sum / lu[i + n * i];
  }
}

void lu_invert(MAT* X, MAT* Y) {
  UINT* idx;
  ITER i, j, k;
  DTYPE* col;
  UINT n, t = 0;

  n = X->d0;

  idx = iip_malloc(sizeof(UINT) * n);
  for(i=0;i<n;i++)idx[i]=i;
  col = iip_malloc(sizeof(DTYPE) * n);

  lu_decomposition(X->data, n, idx);
  printf("==== DECOMPOSITON ====\n");
  for (i = 0; i < n; i++) printf("%u ", idx[i]);
  printf("\nLU MATRIX : \n");
  print_MAT(X);
/*
  for (j = 0; j < n; j++) {
    for (i = 0; i < n; i++) col[i] = 0.0;
    col[j] = 1.0;

    for (k = 0; k < n; k++) printf("%lf ", col[k]);printf("\n");
    
    lu_backwardsubstitution(X->data, n, idx, col);
    for (i = 0; i < n; i++) Y->data[i + n * j] = col[i];
  }
*/
  wiki_invert(X->data,n,idx,Y->data);
  iip_free(idx);
  iip_free(col);
}

void cig_lu(DTYPE* A, int n) {
  DTYPE s;
  ITER i, j, k;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      if (j >= i) {  // U
        s = 0;
        for (k = 0; k < i; k++) s += A[k + j * n] * A[i + k * n];
        A[i + j * n] -= s;
      } else {  // L
        s = 0;
        for (k = 0; k < j; k++) s += A[k + j * n] * A[i + k * n];
        A[i + j * n] = (A[i + j * n] - s) / A[j + j * n];
      }
    }
  }
}

void cig_lusolve(DTYPE* LU, DTYPE* b, int n) {
  // Solve Ly = b
  ITER i, j;
  for (j = 0; j < n; j++) {
    for (i = j + 1; i < n; i++) b[i] -= LU[i + j * n] * b[j];
  }
  // Solve Ux = y
  for (i = n - 1; j >= 0; j--) {
    b[j] /= LU[j + j * n];
    for (i = 0; i < j; i++) b[i] -= LU[i + j * n] * b[j];
  }
}

void lu_dcmp(DTYPE* X, UINT n) {
  ITER i, j, k;
  for (i = 1; i < n; i++) {
    for (j = 0; j < n; j++) {
      if (i > j) {
        for (k = 0; k < j; k++) X[j * n + i] -= X[k * n + i] * X[j * n + k];
        X[j * n + i] /= X[j * n + j];
      } else {
        for (k = 0; k < i; k++) {
          X[j * n + i] -= X[k * n + i] * X[j * n + k];
        }
      }
    }
  }
}



void wiki_invert(DTYPE*X, UINT n,UINT* idx,DTYPE*Y){
  ITER i,j,k;
  /*
  X[0] = 5;
  X[1] = 0;
  X[2] = -0.2;
  X[3] = 0.2;
  X[4] = 4;
  X[5] = 1;
  X[6] = -.2;
  X[7] = .2;
  X[8] = 2;
  X[9] = -1;
  X[10] = 3.2;
  X[11] = -0.375;
  X[12] = 1;
  X[13] = -1;
  X[14] = 0;
  X[15] = 2;
  idx[0]=0;
  idx[1]=1;
  idx[2]=2;
  idx[3]=3;
  */
  for (j = 0; j < n; j++) {
        for (i = 0; i < n; i++) {
            if (idx[i] == j) 
                Y[i+ n*j] = 1.0;
            else
                Y[i+ n*j] = 0.0;

            for (k = 0; k < i; k++)
                Y[i+n*j] -= X[idx[i]+n*k] * Y[k+n*j];
        }

        for (i = n - 1; i >= 0; i--) {
            for (k = i + 1; k < n; k++)
                Y[i+n*j] -= X[idx[i]+n*k] * Y[k+n*j];

            Y[i+n*j] = Y[i+n*j] / X[idx[i]+n*i];
           // printf("%lf\n",X[i+n*i]);
        }
    }
}


void wiki_dcmp(DTYPE*X,UINT n,UINT*idx)
{
  ITER i, j, k, imax; 
  DTYPE maxA, *ptr, absA;
  DTYPE t;
    for (i = 0; i <n; i++)
       idx[i] = i; //Unit permutation matrix, P[N] initialized with N

    for(i=0;i<n;i++);

    for (i = 0; i < n; i++) {
        maxA = 0.0;
        imax = i;
        for (k = 0; k < n; k++)
            if ((absA = fabs(X[idx[k]+n*i])) > maxA) { 
                maxA = absA;
                t = X[idx[k]+ n*i];
                imax = idx[k];
            }
        if (maxA < FZERO){ printf("FAIL\n");return ;} //failure, matrix is degenerate

        if (imax != i) {
            //pivoting P
            j = idx[i];
            idx[i] = idx[imax];
            idx[imax] = j;
            printf("%lf : %d <->  %d\n",t,idx[imax],idx[i]);
        }
    }

 /*
        for (j = i + 1; j < n; j++) {
            X[idx[j] + n*i] /= X[idx[i] + n*i];

            for (k = i + 1; k < n; k++)
                X[idx[j] +n*k] -= X[idx[j] + n*i] * X[idx[i] + n*k];
        }
  */

    for (i = 0; i <n; i++){
     for (j = 0; j < n; j++) {
       if (i > j) {
         for (k = 0; k < j; k++) X[j*n + idx[i]] -= X[k * n + idx[i]] * X[j * n + idx[k]];
         X[j * n + idx[i]] /= X[j * n + idx[j]];
       } 
       else {
         for (k = 0; k < i; k++)
           X[j * n + idx[i]] -= X[k * n + idx[i]] * X[j * n + idx[k]];
       }
     }
  }

}


void wiki(MAT* X, MAT* Y) {
  UINT* idx;
  ITER i, j, k;
  UINT n, t = 0;

  n = X->d0;

  idx = iip_malloc(sizeof(UINT) * n);

  wiki_dcmp(X->data,n,idx); 
  printf("==== DECOMPOSITON ====\n");
  for (i = 0; i < n; i++) printf("%u ", idx[i]);
  printf("\nLU MATRIX : \n");
  print_MAT(X);
  
  wiki_invert(X->data,n,idx,Y->data);
  iip_free(idx);
}
