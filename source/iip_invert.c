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
#include "iip_invert.h"

/**** get inverse of matrix ****/

void invert(MAT* mat, MAT* inv) {
  ITER i;
  UINT size;
  ASSERT((mat->d0 == mat->d1), "This function requires SQUARE MATRIX.\n");
  ASSERT((inv->d0 == inv->d1), "This function requires SQUARE MATRIX.\n");
  if (mat->d0 != inv->d0) ASSERT_DIM_INVALID()

  size = mat->d0 * mat->d1;

  for (i = 0; i < mat->d2; i++) {
    if (mat->d0 == 2)
      invert_2by2(&(mat->data[i * size]), &(inv->data[i * size]));
    else if (mat->d0 == 3)
      invert_3by3(&(mat->data[i * size]), &(inv->data[i * size]));
    else if (mat->d0 == 4)
      invert_4by4(&(mat->data[i * size]), &(inv->data[i * size]));
    else if (mat->d0 == 5)
      invert_5by5(&(mat->data[i * size]), &(inv->data[i * size]));
    else if (mat->d0 == 6)
      invert_6by6(&(mat->data[i * size]), &(inv->data[i * size]));
    else if (mat->d0 > 6)
      invert_nbyn(&(mat->data[i * size]), &(inv->data[i * size]), mat->d0);
    else
      ASSERT_DIM_INVALID()
  }
}

void cinvert(CMAT* mat, CMAT* inv) {
  ITER i;
  UINT size;
  ASSERT((mat->d0 == mat->d1), "This function requires SQUARE MATRIX.\n");
  ASSERT((inv->d0 == inv->d1), "This function requires SQUARE MATRIX.\n");
  if (mat->d0 != inv->d0) ASSERT_DIM_INVALID()

  size = mat->d0 * mat->d1;
  for (i = 0; i < mat->d2; i++) {
    if (mat->d0 == 2)
      cinvert_2by2(&(mat->data[i * size]), &(inv->data[i * size]));
    else if (mat->d0 == 3)
      cinvert_3by3(&(mat->data[i * size]), &(inv->data[i * size]));
    else if (mat->d0 == 4)
      cinvert_4by4(&(mat->data[i * size]), &(inv->data[i * size]));
    else if (mat->d0 == 5)
      cinvert_5by5(&(mat->data[i * size]), &(inv->data[i * size]));
    else if (mat->d0 == 6)
      cinvert_6by6(&(mat->data[i * size]), &(inv->data[i * size]));
    else if (mat->d0 > 6)
      cinvert_nbyn(&(mat->data[i * size]), &(inv->data[i * size]), mat->d0);
    else
      ASSERT_DIM_INVALID()
  }
}

void invert_2by2(DTYPE* X, DTYPE* Y) {
  DTYPE det;
#if DEBUG
  printf("%s\n", __func__);
#endif

  det = X[0] * X[3] - X[2] * X[1];

  if (det > -FZERO && det < FZERO) ASSERT(0, "This matrix is singural.\n")
  det = 1 / det;
  Y[0] = X[3] * det;
  Y[3] = X[0] * det;
  Y[1] = -X[1] * det;
  Y[2] = -X[2] * det;
}

void cinvert_2by2(CTYPE* X, CTYPE* Y) {
  CTYPE det;
  DTYPE t;
#if DEBUG
  printf("%s\n", __func__);
#endif

  det.re = (X[0].re * X[3].re - X[0].im * X[3].im) -
           (X[2].re * X[1].re - X[2].im * X[1].im);
  det.im = (X[0].re * X[3].im + X[0].im * X[3].re) -
           (X[2].re * X[1].im + X[2].im * X[1].re);

#if NTYPE == 0
  if (cabsf(CXF(det)) < FZERO) ASSERT(0, "This matrix is singural.\n")
#elif NTYPE == 1
  if (cabs(CXD(det)) < FZERO) ASSERT(0, "This matrix is singural.\n")
#endif
  t = det.re;
  det.re = (det.re) / (det.re * det.re + det.im * det.im);
  det.im = -(det.im) / (t * t + det.im * det.im);
  Y[0].re = (X[3].re * det.re - X[3].im * det.im);
  Y[0].im = (X[3].re * det.im + X[3].im * det.re);

  Y[3].re = (X[0].re * det.re - X[0].im * det.im);
  Y[3].im = (X[0].re * det.im + X[0].im * det.re);

  Y[1].re = -(X[1].re * det.re - X[1].im * det.im);
  Y[1].im = -(X[1].re * det.im + X[1].im * det.re);

  Y[2].re = -(X[2].re * det.re - X[2].im * det.im);
  Y[2].im = -(X[2].re * det.im + X[2].im * det.re);
}

void invert_3by3(DTYPE* X, DTYPE* Y) {
  DTYPE det;
#if DEBUG
  printf("%s\n", __func__);
#endif

  Y[0] = X[4] * X[8] - X[7] * X[5];
  Y[1] = X[7] * X[2] - X[1] * X[8];
  Y[2] = X[1] * X[5] - X[4] * X[2];

  det = X[0] * Y[0] + X[3] * Y[1] + X[6] * Y[2];

  if (det > -FZERO && det < FZERO) ASSERT(0, "This matrix is singural.\n")

  Y[3] = X[6] * X[5] - X[3] * X[8];
  Y[4] = X[0] * X[8] - X[6] * X[2];
  Y[5] = X[3] * X[2] - X[0] * X[5];
  Y[6] = X[3] * X[7] - X[6] * X[4];
  Y[7] = X[6] * X[1] - X[0] * X[7];
  Y[8] = X[0] * X[4] - X[3] * X[1];

  det = 1 / det;
  Y[0] = Y[0] * det;
  Y[1] = Y[1] * det;
  Y[2] = Y[2] * det;
  Y[3] = Y[3] * det;
  Y[4] = Y[4] * det;
  Y[5] = Y[5] * det;
  Y[6] = Y[6] * det;
  Y[7] = Y[7] * det;
  Y[8] = Y[8] * det;
}

void cinvert_3by3(CTYPE* X, CTYPE* Y) {
  CTYPE det;
  DTYPE t;
#if DEBUG
  printf("%s\n", __func__);
#endif
  Y[0].re = (X[4].re * X[8].re - X[4].im * X[8].im) -
            (X[7].re * X[5].re - X[7].im * X[5].im);
  Y[0].im = (X[4].re * X[8].im + X[4].im * X[8].re) -
            (X[7].re * X[5].im + X[7].im * X[5].re);
  Y[1].re = (X[7].re * X[2].re - X[7].im * X[2].im) -
            (X[1].re * X[8].re - X[1].im * X[8].im);
  Y[1].im = (X[7].re * X[2].im + X[7].im * X[2].re) -
            (X[1].re * X[8].im + X[1].im * X[8].re);

  Y[2].re = (X[1].re * X[5].re - X[1].im * X[5].im) -
            (X[4].re * X[2].re - X[4].im * X[2].im);
  Y[2].im = (X[1].re * X[5].im + X[1].im * X[5].re) -
            (X[4].re * X[2].im + X[4].im * X[2].re);
  det.re = (X[0].re * Y[0].re - X[0].im * Y[0].im) +
           (X[3].re * Y[1].re - X[3].im * Y[1].im) +
           (X[6].re * Y[2].re - X[6].im * Y[2].im);
  det.im = (X[0].re * Y[0].im + X[0].im * Y[0].re) +
           (X[3].re * Y[1].im + X[3].im * Y[1].re) +
           (X[6].re * Y[2].im + X[6].im * Y[2].re);

#if NTYPE == 0
  if (cabsf(CXF(det)) < FZERO) ASSERT(0, "This matrix is singural.\n")
#elif NTYPE == 1
  if (cabs(CXD(det)) < FZERO) ASSERT(0, "This matrix is singural.\n")
#endif

  Y[3].re = (X[6].re * X[5].re - X[6].im * X[5].im) -
            (X[3].re * X[8].re - X[3].im * X[8].im);
  Y[3].im = (X[6].re * X[5].im + X[6].im * X[5].re) -
            (X[3].re * X[8].im + X[3].im * X[8].re);

  Y[4].re = (X[0].re * X[8].re - X[0].im * X[8].im) -
            (X[6].re * X[2].re - X[6].im * X[2].im);
  Y[4].im = (X[0].re * X[8].im + X[0].im * X[8].re) -
            (X[6].re * X[2].im + X[6].im * X[2].re);
  Y[5].re = (X[3].re * X[2].re - X[3].im * X[2].im) -
            (X[0].re * X[5].re - X[0].im * X[5].im);
  Y[5].im = (X[3].re * X[2].im + X[3].im * X[2].re) -
            (X[0].re * X[5].im + X[0].im * X[5].re);

  Y[6].re = (X[3].re * X[7].re - X[3].im * X[7].im) -
            (X[6].re * X[4].re - X[6].im * X[4].im);
  Y[6].im = (X[3].re * X[7].im + X[3].im * X[7].re) -
            (X[6].re * X[4].im + X[6].im * X[4].re);

  Y[7].re = (X[6].re * X[1].re - X[6].im * X[1].im) -
            (X[0].re * X[7].re - X[0].im * X[7].im);
  Y[7].im = (X[6].re * X[1].im + X[6].im * X[1].re) -
            (X[0].re * X[7].im + X[0].im * X[7].re);

  Y[8].re = (X[0].re * X[4].re - X[0].im * X[4].im) -
            (X[3].re * X[1].re - X[3].im * X[1].im);
  Y[8].im = (X[0].re * X[4].im + X[0].im * X[4].re) -
            (X[3].re * X[1].im + X[3].im * X[1].re);

  t = det.re;
  det.re = (det.re) / (det.re * det.re + det.im * det.im);
  det.im = -(det.im) / (t * t + det.im * det.im);

  t = Y[0].re;
  Y[0].re = (Y[0].re * det.re - Y[0].im * det.im);
  Y[0].im = (t * det.im + Y[0].im * det.re);

  t = Y[1].re;
  Y[1].re = (Y[1].re * det.re - Y[1].im * det.im);
  Y[1].im = (t * det.im + Y[1].im * det.re);

  t = Y[2].re;
  Y[2].re = (Y[2].re * det.re - Y[2].im * det.im);
  Y[2].im = (t * det.im + Y[2].im * det.re);

  t = Y[3].re;
  Y[3].re = (Y[3].re * det.re - Y[3].im * det.im);
  Y[3].im = (t * det.im + Y[3].im * det.re);

  t = Y[4].re;
  Y[4].re = (Y[4].re * det.re - Y[4].im * det.im);
  Y[4].im = (t * det.im + Y[4].im * det.re);

  t = Y[5].re;
  Y[5].re = (Y[5].re * det.re - Y[5].im * det.im);
  Y[5].im = (t * det.im + Y[5].im * det.re);

  t = Y[6].re;
  Y[6].re = (Y[6].re * det.re - Y[6].im * det.im);
  Y[6].im = (t * det.im + Y[6].im * det.re);

  t = Y[7].re;
  Y[7].re = (Y[7].re * det.re - Y[7].im * det.im);
  Y[7].im = (t * det.im + Y[7].im * det.re);

  t = Y[8].re;
  Y[8].re = (Y[8].re * det.re - Y[8].im * det.im);
  Y[8].im = (t * det.im + Y[8].im * det.re);
}

void invert_4by4(DTYPE* X, DTYPE* Y) {
  DTYPE det;
  DTYPE t1, t2, t3, t4, t5;
#if DEBUG
  printf("%s\n", __func__);
#endif

  t1 = X[10] * X[15] - X[14] * X[11];
  t2 = X[6] * X[15] - X[14] * X[7];
  t3 = X[6] * X[11] - X[10] * X[7];

  Y[0] = X[5] * t1 - X[9] * t2 + X[13] * t3;
  Y[4] = X[8] * t2 - X[4] * t1 - X[12] * t3;

  t4 = X[2] * X[15] - X[14] * X[3];
  t5 = X[2] * X[11] - X[10] * X[3];

  Y[1] = X[9] * t4 - X[1] * t1 - X[13] * t5;
  Y[5] = X[0] * t1 - X[8] * t4 + X[12] * t5;

  t1 = X[2] * X[7] - X[6] * X[3];

  Y[2] = X[1] * t2 - X[5] * t4 + X[13] * t1;
  Y[6] = X[4] * t4 - X[0] * t2 - X[12] * t1;
  Y[3] = X[5] * t5 - X[1] * t3 - X[9] * t1;

  det = X[0] * Y[0] + X[4] * Y[1] + X[8] * Y[2] + X[12] * Y[3];
  if (det > -FZERO && det < FZERO) ASSERT(0, "This matrix is singural.\n")

  Y[7] = X[0] * t3 - X[4] * t5 + X[8] * t1;

  t1 = X[8] * X[13] - X[12] * X[9];
  t2 = X[4] * X[13] - X[12] * X[5];
  t3 = X[4] * X[9] - X[8] * X[5];

  Y[8] = X[7] * t1 - X[11] * t2 + X[15] * t3;
  Y[12] = X[10] * t2 - X[6] * t1 - X[14] * t3;

  t4 = X[0] * X[13] - X[12] * X[1];
  t5 = X[0] * X[9] - X[8] * X[1];

  Y[9] = X[11] * t4 - X[3] * t1 - X[15] * t5;
  Y[13] = X[2] * t1 - X[10] * t4 + X[14] * t5;

  t1 = X[0] * X[5] - X[4] * X[1];

  Y[10] = X[3] * t2 - X[7] * t4 + X[15] * t1;
  Y[14] = X[6] * t4 - X[2] * t2 - X[14] * t1;
  Y[11] = X[7] * t5 - X[3] * t3 - X[11] * t1;
  Y[15] = X[2] * t3 - X[6] * t5 + X[10] * t1;

  det = 1. / det;

  Y[0] = Y[0] * det;
  Y[1] = Y[1] * det;
  Y[2] = Y[2] * det;
  Y[3] = Y[3] * det;
  Y[4] = Y[4] * det;
  Y[5] = Y[5] * det;
  Y[6] = Y[6] * det;
  Y[7] = Y[7] * det;
  Y[8] = Y[8] * det;
  Y[9] = Y[9] * det;
  Y[10] = Y[10] * det;
  Y[11] = Y[11] * det;
  Y[12] = Y[12] * det;
  Y[13] = Y[13] * det;
  Y[14] = Y[14] * det;
  Y[15] = Y[15] * det;
}

void cinvert_4by4(CTYPE* X, CTYPE* Y) {
  CTYPE det;
  CTYPE t1, t2, t3, t4, t5;
  DTYPE t;
  ITER i;
#if DEBUG
  printf("%s\n", __func__);
#endif
  t1.re = (X[10].re * X[15].re - X[10].im * X[15].im) -
          (X[14].re * X[11].re - X[14].im * X[11].im);
  t1.im = (X[10].re * X[15].im + X[10].im * X[15].re) -
          (X[14].re * X[11].im + X[14].im * X[11].re);

  t2.re = (X[6].re * X[15].re - X[6].im * X[15].im) -
          (X[14].re * X[7].re - X[14].im * X[7].im);
  t2.im = (X[6].re * X[15].im + X[6].im * X[15].re) -
          (X[14].re * X[7].im + X[14].im * X[7].re);

  t3.re = (X[6].re * X[11].re - X[6].im * X[11].im) -
          (X[10].re * X[7].re - X[10].im * X[7].im);
  t3.im = (X[6].re * X[11].im + X[6].im * X[11].re) -
          (X[10].re * X[7].im + X[10].im * X[7].re);

  Y[0].re = (X[5].re * t1.re - X[5].im * t1.im) -
            (X[9].re * t2.re - X[9].im * t2.im) +
            (X[13].re * t3.re - X[13].im * t3.im);
  Y[0].im = (X[5].re * t1.im + X[5].im * t1.re) -
            (X[9].re * t2.im + X[9].im * t2.re) +
            (X[13].re * t3.im + X[13].im * t3.re);

  Y[4].re = (X[8].re * t2.re - X[8].im * t2.im) -
            (X[4].re * t1.re - X[4].im * t1.im) -
            (X[12].re * t3.re - X[12].im * t3.im);
  Y[4].im = (X[8].re * t2.im + X[8].im * t2.re) -
            (X[4].re * t1.im + X[4].im * t1.re) -
            (X[12].re * t3.im + X[12].im * t3.re);

  t4.re = (X[2].re * X[15].re - X[2].im * X[15].im) -
          (X[14].re * X[3].re - X[14].im * X[3].im);
  t4.im = (X[2].re * X[15].im + X[2].im * X[15].re) -
          (X[14].re * X[3].im + X[14].im * X[3].re);

  t5.re = (X[2].re * X[11].re - X[2].im * X[11].im) -
          (X[10].re * X[3].re - X[10].im * X[3].im);
  t5.im = (X[2].re * X[11].im + X[2].im * X[11].re) -
          (X[10].re * X[3].im + X[10].im * X[3].re);

  Y[1].re = (X[9].re * t4.re - X[9].im * t4.im) -
            (X[1].re * t1.re - X[1].im * t1.im) -
            (X[13].re * t5.re - X[13].im * t5.im);
  Y[1].im = (X[9].re * t4.im + X[9].im * t4.re) -
            (X[1].re * t1.im + X[1].im * t1.re) -
            (X[13].re * t5.im + X[13].im * t5.re);

  Y[5].re = (X[0].re * t1.re - X[0].im * t1.im) -
            (X[8].re * t4.re - X[8].im * t4.im) +
            (X[12].re * t5.re - X[12].im * t5.im);
  Y[5].im = (X[0].re * t1.im + X[0].im * t1.re) -
            (X[8].re * t4.im + X[8].im * t4.re) +
            (X[12].re * t5.im + X[12].im * t5.re);

  t1.re = (X[2].re * X[7].re - X[2].im * X[7].im) -
          (X[6].re * X[3].re - X[6].im * X[3].im);
  t1.im = (X[2].re * X[7].im + X[2].im * X[7].re) -
          (X[6].re * X[3].im + X[6].im * X[3].re);

  Y[2].re = (X[1].re * t2.re - X[1].im * t2.im) -
            (X[5].re * t4.re - X[5].im * t4.im) +
            (X[13].re * t1.re - X[13].im * t1.im);
  Y[2].im = (X[1].re * t2.im + X[1].im * t2.re) -
            (X[5].re * t4.im + X[5].im * t4.re) +
            (X[13].re * t1.im + X[13].im * t1.re);

  Y[6].re = (X[4].re * t4.re - X[4].im * t4.im) -
            (X[0].re * t2.re - X[0].im * t2.im) -
            (X[12].re * t1.re - X[12].im * t1.im);
  Y[6].im = (X[4].re * t4.im + X[4].im * t4.re) -
            (X[0].re * t2.im + X[0].im * t2.re) -
            (X[12].re * t1.im + X[12].im * t1.re);

  Y[3].re = (X[5].re * t5.re - X[5].im * t5.im) -
            (X[1].re * t3.re - X[1].im * t3.im) -
            (X[9].re * t1.re - X[9].im * t1.im);
  Y[3].im = (X[5].re * t5.im + X[5].im * t5.re) -
            (X[1].re * t3.im + X[1].im * t3.re) -
            (X[9].re * t1.im + X[9].im * t1.re);

  det.re = (X[0].re * Y[0].re - X[0].im * Y[0].im) +
           (X[4].re * Y[1].re - X[4].im * Y[1].im) +
           (X[8].re * Y[2].re - X[8].im * Y[2].im) +
           (X[12].re * Y[3].re - X[12].im * Y[3].im);
  det.im = (X[0].re * Y[0].im + X[0].im * Y[0].re) +
           (X[4].re * Y[1].im + X[4].im * Y[1].re) +
           (X[8].re * Y[2].im + X[8].im * Y[2].re) +
           (X[12].re * Y[3].im + X[12].im * Y[3].re);

#if NTYPE == 0
  if (cabsf(CXF(det)) < FZERO) ASSERT(0, "This matrix is singural.\n")
#elif NTYPE == 1
  if (cabs(CXD(det)) < FZERO) ASSERT(0, "This matrix is singural.\n")
#endif

  Y[7].re = (X[0].re * t3.re - X[0].im * t3.im) -
            (X[4].re * t5.re - X[4].im * t5.im) +
            (X[8].re * t1.re - X[8].im * t1.im);
  Y[7].im = (X[0].re * t3.im + X[0].im * t3.re) -
            (X[4].re * t5.im + X[4].im * t5.re) +
            (X[8].re * t1.im + X[8].im * t1.re);

  t1.re = (X[8].re * X[13].re - X[8].im * X[13].im) -
          (X[12].re * X[9].re - X[12].im * X[9].im);
  t1.im = (X[8].re * X[13].im + X[8].im * X[13].re) -
          (X[12].re * X[9].im + X[12].im * X[9].re);

  t2.re = (X[4].re * X[13].re - X[4].im * X[13].im) -
          (X[12].re * X[5].re - X[12].im * X[5].im);
  t2.im = (X[4].re * X[13].im + X[4].im * X[13].re) -
          (X[12].re * X[5].im + X[12].im * X[5].re);

  t3.re = (X[4].re * X[9].re - X[4].im * X[9].im) -
          (X[8].re * X[5].re - X[8].im * X[5].im);
  t3.im = (X[4].re * X[9].im + X[4].im * X[9].re) -
          (X[8].re * X[5].im + X[8].im * X[5].re);

  Y[8].re = (X[7].re * t1.re - X[7].im * t1.im) -
            (X[11].re * t2.re - X[11].im * t2.im) +
            (X[15].re * t3.re - X[15].im * t3.im);
  Y[8].im = (X[7].re * t1.im + X[7].im * t1.re) -
            (X[11].re * t2.im + X[11].im * t2.re) +
            (X[15].re * t3.im + X[15].im * t3.re);

  Y[12].re = (X[10].re * t2.re - X[10].im * t2.im) -
             (X[6].re * t1.re - X[6].im * t1.im) -
             (X[14].re * t3.re - X[14].im * t3.im);
  Y[12].im = (X[10].re * t2.im + X[10].im * t2.re) -
             (X[6].re * t1.im + X[6].im * t1.re) -
             (X[14].re * t3.im + X[14].im * t3.re);

  t4.re = (X[0].re * X[13].re - X[0].im * X[13].im) -
          (X[12].re * X[1].re - X[12].im * X[1].im);
  t4.im = (X[0].re * X[13].im + X[0].im * X[13].re) -
          (X[12].re * X[1].im + X[12].im * X[1].re);

  t5.re = (X[0].re * X[9].re - X[0].im * X[9].im) -
          (X[8].re * X[1].re - X[8].im * X[1].im);
  t5.im = (X[0].re * X[9].im + X[0].im * X[9].re) -
          (X[8].re * X[1].im + X[8].im * X[1].re);

  Y[9].re = (X[11].re * t4.re - X[11].im * t4.im) -
            (X[3].re * t1.re - X[3].im * t1.im) -
            (X[15].re * t5.re - X[15].im * t5.im);
  Y[9].im = (X[11].re * t4.im + X[11].im * t4.re) -
            (X[3].re * t1.im + X[3].im * t1.re) -
            (X[15].re * t5.im + X[15].im * t5.re);

  Y[13].re = (X[2].re * t1.re - X[2].im * t1.im) -
             (X[10].re * t4.re - X[10].im * t4.im) +
             (X[14].re * t5.re - X[14].im * t5.im);
  Y[13].im = (X[2].re * t1.im + X[2].im * t1.re) -
             (X[10].re * t4.im + X[10].im * t4.re) +
             (X[14].re * t5.im + X[14].im * t5.re);

  t1.re = (X[0].re * X[5].re - X[0].im * X[5].im) -
          (X[4].re * X[1].re - X[4].im * X[1].im);
  t1.im = (X[0].re * X[5].im + X[0].im * X[5].re) -
          (X[4].re * X[1].im + X[4].im * X[1].re);

  Y[10].re = (X[3].re * t2.re - X[3].im * t2.im) -
             (X[7].re * t4.re - X[7].im * t4.im) +
             (X[15].re * t1.re - X[15].im * t1.im);
  Y[10].im = (X[3].re * t2.im + X[3].im * t2.re) -
             (X[7].re * t4.im + X[7].im * t4.re) +
             (X[15].re * t1.im + X[15].im * t1.re);

  Y[14].re = (X[6].re * t4.re - X[6].im * t4.im) -
             (X[2].re * t2.re - X[2].im * t2.im) -
             (X[14].re * t1.re - X[14].im * t1.im);
  Y[14].im = (X[6].re * t4.im + X[6].im * t4.re) -
             (X[2].re * t2.im + X[2].im * t2.re) -
             (X[14].re * t1.im + X[14].im * t1.re);

  Y[11].re = (X[7].re * t5.re - X[7].im * t5.im) -
             (X[3].re * t3.re - X[3].im * t3.im) -
             (X[11].re * t1.re - X[11].im * t1.im);
  Y[11].im = (X[7].re * t5.im + X[7].im * t5.re) -
             (X[3].re * t3.im + X[3].im * t3.re) -
             (X[11].re * t1.im + X[11].im * t1.re);

  Y[15].re = (X[2].re * t3.re - X[2].im * t3.im) -
             (X[6].re * t5.re - X[6].im * t5.im) +
             (X[10].re * t1.re - X[10].im * t1.im);
  Y[15].im = (X[2].re * t3.im + X[2].im * t3.re) -
             (X[6].re * t5.im + X[6].im * t5.re) +
             (X[10].re * t1.im + X[10].im * t1.re);
  t = det.re;
  det.re = (det.re) / (det.re * det.re + det.im * det.im);
  det.im = -(det.im) / (t * t + det.im * det.im);

  t = Y[0].re;
  Y[0].re = (Y[0].re * det.re - Y[0].im * det.im);
  Y[0].im = (t * det.im + Y[0].im * det.re);

  t = Y[1].re;
  Y[1].re = (Y[1].re * det.re - Y[1].im * det.im);
  Y[1].im = (t * det.im + Y[1].im * det.re);

  t = Y[2].re;
  Y[2].re = (Y[2].re * det.re - Y[2].im * det.im);
  Y[2].im = (t * det.im + Y[2].im * det.re);

  t = Y[3].re;
  Y[3].re = (Y[3].re * det.re - Y[3].im * det.im);
  Y[3].im = (t * det.im + Y[3].im * det.re);

  t = Y[4].re;
  Y[4].re = (Y[4].re * det.re - Y[4].im * det.im);
  Y[4].im = (t * det.im + Y[4].im * det.re);

  t = Y[5].re;
  Y[5].re = (Y[5].re * det.re - Y[5].im * det.im);
  Y[5].im = (t * det.im + Y[5].im * det.re);

  t = Y[6].re;
  Y[6].re = (Y[6].re * det.re - Y[6].im * det.im);
  Y[6].im = (t * det.im + Y[6].im * det.re);

  t = Y[7].re;
  Y[7].re = (Y[7].re * det.re - Y[7].im * det.im);
  Y[7].im = (t * det.im + Y[7].im * det.re);

  t = Y[8].re;
  Y[8].re = (Y[8].re * det.re - Y[8].im * det.im);
  Y[8].im = (t * det.im + Y[8].im * det.re);

  t = Y[9].re;
  Y[9].re = (Y[9].re * det.re - Y[9].im * det.im);
  Y[9].im = (t * det.im + Y[9].im * det.re);

  t = Y[10].re;
  Y[10].re = (Y[10].re * det.re - Y[10].im * det.im);
  Y[10].im = (t * det.im + Y[10].im * det.re);

  t = Y[11].re;
  Y[11].re = (Y[11].re * det.re - Y[11].im * det.im);
  Y[11].im = (t * det.im + Y[11].im * det.re);

  t = Y[12].re;
  Y[12].re = (Y[12].re * det.re - Y[12].im * det.im);
  Y[12].im = (t * det.im + Y[12].im * det.re);

  t = Y[13].re;
  Y[13].re = (Y[13].re * det.re - Y[13].im * det.im);
  Y[13].im = (t * det.im + Y[13].im * det.re);

  t = Y[14].re;
  Y[14].re = (Y[14].re * det.re - Y[14].im * det.im);
  Y[14].im = (t * det.im + Y[14].im * det.re);

  t = Y[15].re;
  Y[15].re = (Y[15].re * det.re - Y[15].im * det.im);
  Y[15].im = (t * det.im + Y[15].im * det.re);
}

void invert_5by5(DTYPE* X, DTYPE* Y) {
  DTYPE det;
  DTYPE t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16,
      t17, t18, t19, t20;
#if DEBUG
  printf("%s\n", __func__);
#endif

  t1 = X[18] * X[24] - X[23] * X[19];
  t2 = X[13] * X[24] - X[23] * X[14];
  t3 = X[13] * X[19] - X[18] * X[14];
  t4 = X[8] * X[24] - X[23] * X[9];
  t5 = X[8] * X[19] - X[18] * X[9];
  t6 = X[8] * X[14] - X[13] * X[9];
  t7 = X[3] * X[24] - X[23] * X[4];
  t8 = X[3] * X[19] - X[18] * X[4];
  t9 = X[3] * X[14] - X[13] * X[4];
  t10 = X[3] * X[9] - X[8] * X[4];

  t11 = X[12] * t1 - X[17] * t2 + X[22] * t3;
  t12 = X[7] * t1 - X[17] * t4 + X[22] * t5;
  t13 = X[7] * t2 - X[12] * t4 + X[22] * t6;
  t14 = X[7] * t3 - X[12] * t5 + X[17] * t6;
  t15 = X[2] * t1 - X[17] * t7 + X[22] * t8;
  t16 = X[2] * t2 - X[12] * t7 + X[22] * t9;
  t17 = X[2] * t3 - X[12] * t8 + X[17] * t9;

  Y[0] = X[6] * t11 - X[11] * t12 + X[16] * t13 - X[21] * t14;
  Y[5] = -X[5] * t11 + X[10] * t12 - X[15] * t13 + X[20] * t14;
  Y[1] = -X[1] * t11 + X[11] * t15 - X[16] * t16 + X[21] * t17;
  Y[6] = X[0] * t11 - X[10] * t15 + X[15] * t16 - X[20] * t17;

  t18 = X[2] * t4 - X[7] * t7 + X[22] * t10;
  t19 = X[2] * t5 - X[7] * t8 + X[17] * t10;
  t20 = X[2] * t6 - X[7] * t9 + X[12] * t10;

  Y[2] = X[1] * t12 - X[6] * t15 + X[16] * t18 - X[21] * t19;
  Y[7] = -X[0] * t12 + X[5] * t15 - X[15] * t18 + X[20] * t19;
  Y[3] = -X[1] * t13 + X[6] * t16 - X[11] * t18 + X[21] * t20;
  Y[8] = X[0] * t13 - X[5] * t16 + X[10] * t18 - X[20] * t20;
  Y[4] = X[1] * t14 - X[6] * t17 + X[11] * t19 - X[16] * t20;
  Y[9] = -X[0] * t14 + X[5] * t17 - X[10] * t19 + X[15] * t20;

  det = X[0] * Y[0] + X[5] * Y[1] + X[10] * Y[2] + X[15] * Y[3] + X[20] * Y[4];

  if (det > -FZERO && det < FZERO) ASSERT(0, "This matrix is singural.\n")
  t11 = X[11] * t1 - X[16] * t2 + X[21] * t3;
  t12 = X[6] * t1 - X[16] * t4 + X[21] * t5;
  t13 = X[6] * t2 - X[11] * t4 + X[21] * t6;
  t14 = X[6] * t3 - X[11] * t5 + X[16] * t6;
  t15 = X[1] * t1 - X[16] * t7 + X[21] * t8;
  t16 = X[1] * t2 - X[11] * t7 + X[21] * t9;
  t17 = X[1] * t3 - X[11] * t8 + X[16] * t9;
  t18 = X[1] * t4 - X[6] * t7 + X[21] * t10;
  t19 = X[1] * t5 - X[6] * t8 + X[16] * t10;

  Y[10] = X[5] * t11 - X[10] * t12 + X[15] * t13 - X[20] * t14;
  Y[11] = -X[0] * t11 + X[10] * t15 - X[15] * t16 + X[20] * t17;
  Y[12] = X[0] * t12 - X[5] * t15 + X[15] * t18 - X[20] * t19;

  t1 = X[10] * X[16] - X[15] * X[11];
  t2 = X[5] * X[16] - X[15] * X[6];
  t3 = X[5] * X[11] - X[10] * X[6];
  t4 = X[0] * X[16] - X[15] * X[1];
  t5 = X[0] * X[11] - X[10] * X[1];
  t6 = X[0] * X[6] - X[5] * X[1];
  t7 = X[10] * X[21] - X[20] * X[11];
  t8 = X[5] * X[21] - X[20] * X[6];
  t9 = X[0] * X[21] - X[20] * X[1];
  t10 = X[15] * X[21] - X[20] * X[16];

  t11 = X[12] * t10 - X[17] * t7 + X[22] * t1;
  t12 = X[7] * t10 - X[17] * t8 + X[22] * t2;
  t13 = X[7] * t7 - X[12] * t8 + X[22] * t3;
  t14 = X[7] * t1 - X[12] * t2 + X[17] * t3;
  t15 = X[2] * t10 - X[17] * t9 + X[22] * t4;
  t16 = X[2] * t7 - X[12] * t9 + X[22] * t5;
  t17 = X[2] * t1 - X[12] * t4 + X[17] * t5;

  Y[15] = X[9] * t11 - X[14] * t12 + X[19] * t13 - X[24] * t14;
  Y[20] = -X[8] * t11 + X[13] * t12 - X[18] * t13 + X[23] * t14;
  Y[16] = -X[4] * t11 + X[14] * t15 - X[19] * t16 + X[24] * t17;
  Y[21] = X[3] * t11 - X[13] * t15 + X[18] * t16 - X[23] * t17;

  t18 = X[2] * t8 - X[7] * t9 + X[22] * t6;
  t19 = X[2] * t2 - X[7] * t4 + X[17] * t6;
  t20 = X[2] * t3 - X[7] * t5 + X[12] * t6;

  Y[17] = X[4] * t12 - X[9] * t15 + X[19] * t18 - X[24] * t19;
  Y[22] = -X[3] * t12 + X[8] * t15 - X[18] * t18 + X[23] * t19;
  Y[18] = -X[4] * t13 + X[9] * t16 - X[14] * t18 + X[24] * t20;
  Y[23] = X[3] * t13 - X[8] * t16 + X[13] * t18 - X[23] * t20;
  Y[19] = X[4] * t14 - X[9] * t17 + X[14] * t19 - X[19] * t20;
  Y[24] = -X[3] * t14 + X[8] * t17 - X[13] * t19 + X[18] * t20;

  t11 = X[8] * t7 - X[13] * t8 + X[23] * t3;
  t12 = X[3] * t7 - X[13] * t9 + X[23] * t5;
  t13 = X[3] * t8 - X[8] * t9 + X[23] * t6;
  t14 = X[3] * t3 - X[8] * t5 + X[13] * t6;

  t15 = X[8] * t1 - X[13] * t2 + X[18] * t3;
  t16 = X[3] * t1 - X[13] * t4 + X[18] * t5;
  t17 = X[3] * t2 - X[8] * t4 + X[18] * t6;

  Y[13] = X[4] * t11 - X[9] * t12 + X[14] * t13 - X[24] * t14;
  Y[14] = -X[4] * t15 + X[9] * t16 - X[14] * t17 + X[19] * t14;
  det = 1 / det;
  Y[0] = Y[0] * det;
  Y[1] = Y[1] * det;
  Y[2] = Y[2] * det;
  Y[3] = Y[3] * det;
  Y[4] = Y[4] * det;
  Y[5] = Y[5] * det;
  Y[6] = Y[6] * det;
  Y[7] = Y[7] * det;
  Y[8] = Y[8] * det;
  Y[9] = Y[9] * det;
  Y[10] = Y[10] * det;
  Y[11] = Y[11] * det;
  Y[12] = Y[12] * det;
  Y[13] = Y[13] * det;
  Y[14] = Y[14] * det;
  Y[15] = Y[15] * det;
  Y[16] = Y[16] * det;
  Y[17] = Y[17] * det;
  Y[18] = Y[18] * det;
  Y[19] = Y[19] * det;
  Y[20] = Y[20] * det;
  Y[21] = Y[21] * det;
  Y[22] = Y[22] * det;
  Y[23] = Y[23] * det;
  Y[24] = Y[24] * det;
}

void cinvert_5by5(CTYPE* X, CTYPE* Y) {
  CTYPE det;
  DTYPE t;
  CTYPE t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16,
      t17, t18, t19, t20;
#if DEBUG
  printf("%s\n", __func__);
#endif

  t1.re = (X[18].re * X[24].re - X[18].im * X[24].im) -
          (X[23].re * X[19].re - X[23].im * X[19].im);
  t1.im = (X[18].re * X[24].im + X[18].im * X[24].re) -
          (X[23].re * X[19].im + X[23].im * X[19].re);

  t2.re = (X[13].re * X[24].re - X[13].im * X[24].im) -
          (X[23].re * X[14].re - X[23].im * X[14].im);
  t2.im = (X[13].re * X[24].im + X[13].im * X[24].re) -
          (X[23].re * X[14].im + X[23].im * X[14].re);

  t3.re = (X[13].re * X[19].re - X[13].im * X[19].im) -
          (X[18].re * X[14].re - X[18].im * X[14].im);
  t3.im = (X[13].re * X[19].im + X[13].im * X[19].re) -
          (X[18].re * X[14].im + X[18].im * X[14].re);

  t4.re = (X[8].re * X[24].re - X[8].im * X[24].im) -
          (X[23].re * X[9].re - X[23].im * X[9].im);
  t4.im = (X[8].re * X[24].im + X[8].im * X[24].re) -
          (X[23].re * X[9].im + X[23].im * X[9].re);

  t5.re = (X[8].re * X[19].re - X[8].im * X[19].im) -
          (X[18].re * X[9].re - X[18].im * X[9].im);
  t5.im = (X[8].re * X[19].im + X[8].im * X[19].re) -
          (X[18].re * X[9].im + X[18].im * X[9].re);

  t6.re = (X[8].re * X[14].re - X[8].im * X[14].im) -
          (X[13].re * X[9].re - X[13].im * X[9].im);
  t6.im = (X[8].re * X[14].im + X[8].im * X[14].re) -
          (X[13].re * X[9].im + X[13].im * X[9].re);

  t7.re = (X[3].re * X[24].re - X[3].im * X[24].im) -
          (X[23].re * X[4].re - X[23].im * X[4].im);
  t7.im = (X[3].re * X[24].im + X[3].im * X[24].re) -
          (X[23].re * X[4].im + X[23].im * X[4].re);

  t8.re = (X[3].re * X[19].re - X[3].im * X[19].im) -
          (X[18].re * X[4].re - X[18].im * X[4].im);
  t8.im = (X[3].re * X[19].im + X[3].im * X[19].re) -
          (X[18].re * X[4].im + X[18].im * X[4].re);

  t9.re = (X[3].re * X[14].re - X[3].im * X[14].im) -
          (X[13].re * X[4].re - X[13].im * X[4].im);
  t9.im = (X[3].re * X[14].im + X[3].im * X[14].re) -
          (X[13].re * X[4].im + X[13].im * X[4].re);

  t10.re = (X[3].re * X[9].re - X[3].im * X[9].im) -
           (X[8].re * X[4].re - X[8].im * X[4].im);
  t10.im = (X[3].re * X[9].im + X[3].im * X[9].re) -
           (X[8].re * X[4].im + X[8].im * X[4].re);

  t11.re = (X[12].re * t1.re - X[12].im * t1.im) -
           (X[17].re * t2.re - X[17].im * t2.im) +
           (X[22].re * t3.re - X[22].im * t3.im);
  t11.im = (X[12].re * t1.im + X[12].im * t1.re) -
           (X[17].re * t2.im + X[17].im * t2.re) +
           (X[22].re * t3.im + X[22].im * t3.re);

  t12.re = (X[7].re * t1.re - X[7].im * t1.im) -
           (X[17].re * t4.re - X[17].im * t4.im) +
           (X[22].re * t5.re - X[22].im * t5.im);
  t12.im = (X[7].re * t1.im + X[7].im * t1.re) -
           (X[17].re * t4.im + X[17].im * t4.re) +
           (X[22].re * t5.im + X[22].im * t5.re);

  t13.re = (X[7].re * t2.re - X[7].im * t2.im) -
           (X[12].re * t4.re - X[12].im * t4.im) +
           (X[22].re * t6.re - X[22].im * t6.im);
  t13.im = (X[7].re * t2.im + X[7].im * t2.re) -
           (X[12].re * t4.im + X[12].im * t4.re) +
           (X[22].re * t6.im + X[22].im * t6.re);

  t14.re = (X[7].re * t3.re - X[7].im * t3.im) -
           (X[12].re * t5.re - X[12].im * t5.im) +
           (X[17].re * t6.re - X[17].im * t6.im);
  t14.im = (X[7].re * t3.im + X[7].im * t3.re) -
           (X[12].re * t5.im + X[12].im * t5.re) +
           (X[17].re * t6.im + X[17].im * t6.re);

  t15.re = (X[2].re * t1.re - X[2].im * t1.im) -
           (X[17].re * t7.re - X[17].im * t7.im) +
           (X[22].re * t8.re - X[22].im * t8.im);
  t15.im = (X[2].re * t1.im + X[2].im * t1.re) -
           (X[17].re * t7.im + X[17].im * t7.re) +
           (X[22].re * t8.im + X[22].im * t8.re);

  t16.re = (X[2].re * t2.re - X[2].im * t2.im) -
           (X[12].re * t7.re - X[12].im * t7.im) +
           (X[22].re * t9.re - X[22].im * t9.im);
  t16.im = (X[2].re * t2.im + X[2].im * t2.re) -
           (X[12].re * t7.im + X[12].im * t7.re) +
           (X[22].re * t9.im + X[22].im * t9.re);

  t17.re = (X[2].re * t3.re - X[2].im * t3.im) -
           (X[12].re * t8.re - X[12].im * t8.im) +
           (X[17].re * t9.re - X[17].im * t9.im);
  t17.im = (X[2].re * t3.im + X[2].im * t3.re) -
           (X[12].re * t8.im + X[12].im * t8.re) +
           (X[17].re * t9.im + X[17].im * t9.re);

  Y[0].re = (X[6].re * t11.re - X[6].im * t11.im) -
            (X[11].re * t12.re - X[11].im * t12.im) +
            (X[16].re * t13.re - X[16].im * t13.im) -
            (X[21].re * t14.re - X[21].im * t14.im);
  Y[0].im = (X[6].re * t11.im + X[6].im * t11.re) -
            (X[11].re * t12.im + X[11].im * t12.re) +
            (X[16].re * t13.im + X[16].im * t13.re) -
            (X[21].re * t14.im + X[21].im * t14.re);

  Y[5].re = -(X[5].re * t11.re - X[5].im * t11.im) +
            (X[10].re * t12.re - X[10].im * t12.im) -
            (X[15].re * t13.re - X[15].im * t13.im) +
            (X[20].re * t14.re - X[20].im * t14.im);
  Y[5].im = -(X[5].re * t11.im + X[5].im * t11.re) +
            (X[10].re * t12.im + X[10].im * t12.re) -
            (X[15].re * t13.im + X[15].im * t13.re) +
            (X[20].re * t14.im + X[20].im * t14.re);

  Y[1].re = -(X[1].re * t11.re - X[1].im * t11.im) +
            (X[11].re * t15.re - X[11].im * t15.im) -
            (X[16].re * t16.re - X[16].im * t16.im) +
            (X[21].re * t17.re - X[21].im * t17.im);
  Y[1].im = -(X[1].re * t11.im + X[1].im * t11.re) +
            (X[11].re * t15.im + X[11].im * t15.re) -
            (X[16].re * t16.im + X[16].im * t16.re) +
            (X[21].re * t17.im + X[21].im * t17.re);

  Y[6].re = (X[0].re * t11.re - X[0].im * t11.im) -
            (X[10].re * t15.re - X[10].im * t15.im) +
            (X[15].re * t16.re - X[15].im * t16.im) -
            (X[20].re * t17.re - X[20].im * t17.im);
  Y[6].im = (X[0].re * t11.im + X[0].im * t11.re) -
            (X[10].re * t15.im + X[10].im * t15.re) +
            (X[15].re * t16.im + X[15].im * t16.re) -
            (X[20].re * t17.im + X[20].im * t17.re);

  t18.re = (X[2].re * t4.re - X[2].im * t4.im) -
           (X[7].re * t7.re - X[7].im * t7.im) +
           (X[22].re * t10.re - X[22].im * t10.im);
  t18.im = (X[2].re * t4.im + X[2].im * t4.re) -
           (X[7].re * t7.im + X[7].im * t7.re) +
           (X[22].re * t10.im + X[22].im * t10.re);

  t19.re = (X[2].re * t5.re - X[2].im * t5.im) -
           (X[7].re * t8.re - X[7].im * t8.im) +
           (X[17].re * t10.re - X[17].im * t10.im);
  t19.im = (X[2].re * t5.im + X[2].im * t5.re) -
           (X[7].re * t8.im + X[7].im * t8.re) +
           (X[17].re * t10.im + X[17].im * t10.re);

  t20.re = (X[2].re * t6.re - X[2].im * t6.im) -
           (X[7].re * t9.re - X[7].im * t9.im) +
           (X[12].re * t10.re - X[12].im * t10.im);
  t20.im = (X[2].re * t6.im + X[2].im * t6.re) -
           (X[7].re * t9.im + X[7].im * t9.re) +
           (X[12].re * t10.im + X[12].im * t10.re);

  Y[2].re = (X[1].re * t12.re - X[1].im * t12.im) -
            (X[6].re * t15.re - X[6].im * t15.im) +
            (X[16].re * t18.re - X[16].im * t18.im) -
            (X[21].re * t19.re - X[21].im * t19.im);
  Y[2].im = (X[1].re * t12.im + X[1].im * t12.re) -
            (X[6].re * t15.im + X[6].im * t15.re) +
            (X[16].re * t18.im + X[16].im * t18.re) -
            (X[21].re * t19.im + X[21].im * t19.re);

  Y[7].re = -(X[0].re * t12.re - X[0].im * t12.im) +
            (X[5].re * t15.re - X[5].im * t15.im) -
            (X[15].re * t18.re - X[15].im * t18.im) +
            (X[20].re * t19.re - X[20].im * t19.im);
  Y[7].im = -(X[0].re * t12.im + X[0].im * t12.re) +
            (X[5].re * t15.im + X[5].im * t15.re) -
            (X[15].re * t18.im + X[15].im * t18.re) +
            (X[20].re * t19.im + X[20].im * t19.re);

  Y[3].re = -(X[1].re * t13.re - X[1].im * t13.im) +
            (X[6].re * t16.re - X[6].im * t16.im) -
            (X[11].re * t18.re - X[11].im * t18.im) +
            (X[21].re * t20.re - X[21].im * t20.im);
  Y[3].im = -(X[1].re * t13.im + X[1].im * t13.re) +
            (X[6].re * t16.im + X[6].im * t16.re) -
            (X[11].re * t18.im + X[11].im * t18.re) +
            (X[21].re * t20.im + X[21].im * t20.re);

  Y[8].re = (X[0].re * t13.re - X[0].im * t13.im) -
            (X[5].re * t16.re - X[5].im * t16.im) +
            (X[10].re * t18.re - X[10].im * t18.im) -
            (X[20].re * t20.re - X[20].im * t20.im);
  Y[8].im = (X[0].re * t13.im + X[0].im * t13.re) -
            (X[5].re * t16.im + X[5].im * t16.re) +
            (X[10].re * t18.im + X[10].im * t18.re) -
            (X[20].re * t20.im + X[20].im * t20.re);

  Y[4].re = (X[1].re * t14.re - X[1].im * t14.im) -
            (X[6].re * t17.re - X[6].im * t17.im) +
            (X[11].re * t19.re - X[11].im * t19.im) -
            (X[16].re * t20.re - X[16].im * t20.im);
  Y[4].im = (X[1].re * t14.im + X[1].im * t14.re) -
            (X[6].re * t17.im + X[6].im * t17.re) +
            (X[11].re * t19.im + X[11].im * t19.re) -
            (X[16].re * t20.im + X[16].im * t20.re);

  Y[9].re = -(X[0].re * t14.re - X[0].im * t14.im) +
            (X[5].re * t17.re - X[5].im * t17.im) -
            (X[10].re * t19.re - X[10].im * t19.im) +
            (X[15].re * t20.re - X[15].im * t20.im);
  Y[9].im = -(X[0].re * t14.im + X[0].im * t14.re) +
            (X[5].re * t17.im + X[5].im * t17.re) -
            (X[10].re * t19.im + X[10].im * t19.re) +
            (X[15].re * t20.im + X[15].im * t20.re);

  det.re = (X[0].re * Y[0].re - X[0].im * Y[0].im) +
           (X[5].re * Y[1].re - X[5].im * Y[1].im) +
           (X[10].re * Y[2].re - X[10].im * Y[2].im) +
           (X[15].re * Y[3].re - X[15].im * Y[3].im) +
           (X[20].re * Y[4].re - X[20].im * Y[4].im);
  det.im = (X[0].re * Y[0].im + X[0].im * Y[0].re) +
           (X[5].re * Y[1].im + X[5].im * Y[1].re) +
           (X[10].re * Y[2].im + X[10].im * Y[2].re) +
           (X[15].re * Y[3].im + X[15].im * Y[3].re) +
           (X[20].re * Y[4].im + X[20].im * Y[4].re);

#if NTYPE == 0
  if (cabsf(CXF(det)) < FZERO) ASSERT(0, "This matrix is singural.\n")
#elif NTYPE == 1
  if (cabs(CXD(det)) < FZERO) ASSERT(0, "This matrix is singural.\n")
#endif

  t11.re = (X[11].re * t1.re - X[11].im * t1.im) -
           (X[16].re * t2.re - X[16].im * t2.im) +
           (X[21].re * t3.re - X[21].im * t3.im);
  t11.im = (X[11].re * t1.im + X[11].im * t1.re) -
           (X[16].re * t2.im + X[16].im * t2.re) +
           (X[21].re * t3.im + X[21].im * t3.re);

  t12.re = (X[6].re * t1.re - X[6].im * t1.im) -
           (X[16].re * t4.re - X[16].im * t4.im) +
           (X[21].re * t5.re - X[21].im * t5.im);
  t12.im = (X[6].re * t1.im + X[6].im * t1.re) -
           (X[16].re * t4.im + X[16].im * t4.re) +
           (X[21].re * t5.im + X[21].im * t5.re);

  t13.re = (X[6].re * t2.re - X[6].im * t2.im) -
           (X[11].re * t4.re - X[11].im * t4.im) +
           (X[21].re * t6.re - X[21].im * t6.im);
  t13.im = (X[6].re * t2.im + X[6].im * t2.re) -
           (X[11].re * t4.im + X[11].im * t4.re) +
           (X[21].re * t6.im + X[21].im * t6.re);

  t14.re = (X[6].re * t3.re - X[6].im * t3.im) -
           (X[11].re * t5.re - X[11].im * t5.im) +
           (X[16].re * t6.re - X[16].im * t6.im);
  t14.im = (X[6].re * t3.im + X[6].im * t3.re) -
           (X[11].re * t5.im + X[11].im * t5.re) +
           (X[16].re * t6.im + X[16].im * t6.re);

  t15.re = (X[1].re * t1.re - X[1].im * t1.im) -
           (X[16].re * t7.re - X[16].im * t7.im) +
           (X[21].re * t8.re - X[21].im * t8.im);
  t15.im = (X[1].re * t1.im + X[1].im * t1.re) -
           (X[16].re * t7.im + X[16].im * t7.re) +
           (X[21].re * t8.im + X[21].im * t8.re);

  t16.re = (X[1].re * t2.re - X[1].im * t2.im) -
           (X[11].re * t7.re - X[11].im * t7.im) +
           (X[21].re * t9.re - X[21].im * t9.im);
  t16.im = (X[1].re * t2.im + X[1].im * t2.re) -
           (X[11].re * t7.im + X[11].im * t7.re) +
           (X[21].re * t9.im + X[21].im * t9.re);

  t17.re = (X[1].re * t3.re - X[1].im * t3.im) -
           (X[11].re * t8.re - X[11].im * t8.im) +
           (X[16].re * t9.re - X[16].im * t9.im);
  t17.im = (X[1].re * t3.im + X[1].im * t3.re) -
           (X[11].re * t8.im + X[11].im * t8.re) +
           (X[16].re * t9.im + X[16].im * t9.re);

  t18.re = (X[1].re * t4.re - X[1].im * t4.im) -
           (X[6].re * t7.re - X[6].im * t7.im) +
           (X[21].re * t10.re - X[21].im * t10.im);
  t18.im = (X[1].re * t4.im + X[1].im * t4.re) -
           (X[6].re * t7.im + X[6].im * t7.re) +
           (X[21].re * t10.im + X[21].im * t10.re);

  t19.re = (X[1].re * t5.re - X[1].im * t5.im) -
           (X[6].re * t8.re - X[6].im * t8.im) +
           (X[16].re * t10.re - X[16].im * t10.im);
  t19.im = (X[1].re * t5.im + X[1].im * t5.re) -
           (X[6].re * t8.im + X[6].im * t8.re) +
           (X[16].re * t10.im + X[16].im * t10.re);

  Y[10].re = (X[5].re * t11.re - X[5].im * t11.im) -
             (X[10].re * t12.re - X[10].im * t12.im) +
             (X[15].re * t13.re - X[15].im * t13.im) -
             (X[20].re * t14.re - X[20].im * t14.im);
  Y[10].im = (X[5].re * t11.im + X[5].im * t11.re) -
             (X[10].re * t12.im + X[10].im * t12.re) +
             (X[15].re * t13.im + X[15].im * t13.re) -
             (X[20].re * t14.im + X[20].im * t14.re);

  Y[11].re = -(X[0].re * t11.re - X[0].im * t11.im) +
             (X[10].re * t15.re - X[10].im * t15.im) -
             (X[15].re * t16.re - X[15].im * t16.im) +
             (X[20].re * t17.re - X[20].im * t17.im);
  Y[11].im = -(X[0].re * t11.im + X[0].im * t11.re) +
             (X[10].re * t15.im + X[10].im * t15.re) -
             (X[15].re * t16.im + X[15].im * t16.re) +
             (X[20].re * t17.im + X[20].im * t17.re);

  Y[12].re = (X[0].re * t12.re - X[0].im * t12.im) -
             (X[5].re * t15.re - X[5].im * t15.im) +
             (X[15].re * t18.re - X[15].im * t18.im) -
             (X[20].re * t19.re - X[20].im * t19.im);
  Y[12].im = (X[0].re * t12.im + X[0].im * t12.re) -
             (X[5].re * t15.im + X[5].im * t15.re) +
             (X[15].re * t18.im + X[15].im * t18.re) -
             (X[20].re * t19.im + X[20].im * t19.re);

  t1.re = (X[10].re * X[16].re - X[10].im * X[16].im) -
          (X[15].re * X[11].re - X[15].im * X[11].im);
  t1.im = (X[10].re * X[16].im + X[10].im * X[16].re) -
          (X[15].re * X[11].im + X[15].im * X[11].re);

  t2.re = (X[5].re * X[16].re - X[5].im * X[16].im) -
          (X[15].re * X[6].re - X[15].im * X[6].im);
  t2.im = (X[5].re * X[16].im + X[5].im * X[16].re) -
          (X[15].re * X[6].im + X[15].im * X[6].re);

  t3.re = (X[5].re * X[11].re - X[5].im * X[11].im) -
          (X[10].re * X[6].re - X[10].im * X[6].im);
  t3.im = (X[5].re * X[11].im + X[5].im * X[11].re) -
          (X[10].re * X[6].im + X[10].im * X[6].re);

  t4.re = (X[0].re * X[16].re - X[0].im * X[16].im) -
          (X[15].re * X[1].re - X[15].im * X[1].im);
  t4.im = (X[0].re * X[16].im + X[0].im * X[16].re) -
          (X[15].re * X[1].im + X[15].im * X[1].re);

  t5.re = (X[0].re * X[11].re - X[0].im * X[11].im) -
          (X[10].re * X[1].re - X[10].im * X[1].im);
  t5.im = (X[0].re * X[11].im + X[0].im * X[11].re) -
          (X[10].re * X[1].im + X[10].im * X[1].re);

  t6.re = (X[0].re * X[6].re - X[0].im * X[6].im) -
          (X[5].re * X[1].re - X[5].im * X[1].im);
  t6.im = (X[0].re * X[6].im + X[0].im * X[6].re) -
          (X[5].re * X[1].im + X[5].im * X[1].re);

  t7.re = (X[10].re * X[21].re - X[10].im * X[21].im) -
          (X[20].re * X[11].re - X[20].im * X[11].im);
  t7.im = (X[10].re * X[21].im + X[10].im * X[21].re) -
          (X[20].re * X[11].im + X[20].im * X[11].re);

  t8.re = (X[5].re * X[21].re - X[5].im * X[21].im) -
          (X[20].re * X[6].re - X[20].im * X[6].im);
  t8.im = (X[5].re * X[21].im + X[5].im * X[21].re) -
          (X[20].re * X[6].im + X[20].im * X[6].re);

  t9.re = (X[0].re * X[21].re - X[0].im * X[21].im) -
          (X[20].re * X[1].re - X[20].im * X[1].im);
  t9.im = (X[0].re * X[21].im + X[0].im * X[21].re) -
          (X[20].re * X[1].im + X[20].im * X[1].re);

  t10.re = (X[15].re * X[21].re - X[15].im * X[21].im) -
           (X[20].re * X[16].re - X[20].im * X[16].im);
  t10.im = (X[15].re * X[21].im + X[15].im * X[21].re) -
           (X[20].re * X[16].im + X[20].im * X[16].re);

  t11.re = (X[12].re * t10.re - X[12].im * t10.im) -
           (X[17].re * t7.re - X[17].im * t7.im) +
           (X[22].re * t1.re - X[22].im * t1.im);
  t11.im = (X[12].re * t10.im + X[12].im * t10.re) -
           (X[17].re * t7.im + X[17].im * t7.re) +
           (X[22].re * t1.im + X[22].im * t1.re);

  t12.re = (X[7].re * t10.re - X[7].im * t10.im) -
           (X[17].re * t8.re - X[17].im * t8.im) +
           (X[22].re * t2.re - X[22].im * t2.im);
  t12.im = (X[7].re * t10.im + X[7].im * t10.re) -
           (X[17].re * t8.im + X[17].im * t8.re) +
           (X[22].re * t2.im + X[22].im * t2.re);

  t13.re = (X[7].re * t7.re - X[7].im * t7.im) -
           (X[12].re * t8.re - X[12].im * t8.im) +
           (X[22].re * t3.re - X[22].im * t3.im);
  t13.im = (X[7].re * t7.im + X[7].im * t7.re) -
           (X[12].re * t8.im + X[12].im * t8.re) +
           (X[22].re * t3.im + X[22].im * t3.re);

  t14.re = (X[7].re * t1.re - X[7].im * t1.im) -
           (X[12].re * t2.re - X[12].im * t2.im) +
           (X[17].re * t3.re - X[17].im * t3.im);
  t14.im = (X[7].re * t1.im + X[7].im * t1.re) -
           (X[12].re * t2.im + X[12].im * t2.re) +
           (X[17].re * t3.im + X[17].im * t3.re);

  t15.re = (X[2].re * t10.re - X[2].im * t10.im) -
           (X[17].re * t9.re - X[17].im * t9.im) +
           (X[22].re * t4.re - X[22].im * t4.im);
  t15.im = (X[2].re * t10.im + X[2].im * t10.re) -
           (X[17].re * t9.im + X[17].im * t9.re) +
           (X[22].re * t4.im + X[22].im * t4.re);

  t16.re = (X[2].re * t7.re - X[2].im * t7.im) -
           (X[12].re * t9.re - X[12].im * t9.im) +
           (X[22].re * t5.re - X[22].im * t5.im);
  t16.im = (X[2].re * t7.im + X[2].im * t7.re) -
           (X[12].re * t9.im + X[12].im * t9.re) +
           (X[22].re * t5.im + X[22].im * t5.re);

  t17.re = (X[2].re * t1.re - X[2].im * t1.im) -
           (X[12].re * t4.re - X[12].im * t4.im) +
           (X[17].re * t5.re - X[17].im * t5.im);
  t17.im = (X[2].re * t1.im + X[2].im * t1.re) -
           (X[12].re * t4.im + X[12].im * t4.re) +
           (X[17].re * t5.im + X[17].im * t5.re);

  Y[15].re = (X[9].re * t11.re - X[9].im * t11.im) -
             (X[14].re * t12.re - X[14].im * t12.im) +
             (X[19].re * t13.re - X[19].im * t13.im) -
             (X[24].re * t14.re - X[24].im * t14.im);
  Y[15].im = (X[9].re * t11.im + X[9].im * t11.re) -
             (X[14].re * t12.im + X[14].im * t12.re) +
             (X[19].re * t13.im + X[19].im * t13.re) -
             (X[24].re * t14.im + X[24].im * t14.re);

  Y[20].re = -(X[8].re * t11.re - X[8].im * t11.im) +
             (X[13].re * t12.re - X[13].im * t12.im) -
             (X[18].re * t13.re - X[18].im * t13.im) +
             (X[23].re * t14.re - X[23].im * t14.im);
  Y[20].im = -(X[8].re * t11.im + X[8].im * t11.re) +
             (X[13].re * t12.im + X[13].im * t12.re) -
             (X[18].re * t13.im + X[18].im * t13.re) +
             (X[23].re * t14.im + X[23].im * t14.re);

  Y[16].re = -(X[4].re * t11.re - X[4].im * t11.im) +
             (X[14].re * t15.re - X[14].im * t15.im) -
             (X[19].re * t16.re - X[19].im * t16.im) +
             (X[24].re * t17.re - X[24].im * t17.im);
  Y[16].im = -(X[4].re * t11.im + X[4].im * t11.re) +
             (X[14].re * t15.im + X[14].im * t15.re) -
             (X[19].re * t16.im + X[19].im * t16.re) +
             (X[24].re * t17.im + X[24].im * t17.re);

  Y[21].re = (X[3].re * t11.re - X[3].im * t11.im) -
             (X[13].re * t15.re - X[13].im * t15.im) +
             (X[18].re * t16.re - X[18].im * t16.im) -
             (X[23].re * t17.re - X[23].im * t17.im);
  Y[21].im = (X[3].re * t11.im + X[3].im * t11.re) -
             (X[13].re * t15.im + X[13].im * t15.re) +
             (X[18].re * t16.im + X[18].im * t16.re) -
             (X[23].re * t17.im + X[23].im * t17.re);

  t18.re = (X[2].re * t8.re - X[2].im * t8.im) -
           (X[7].re * t9.re - X[7].im * t9.im) +
           (X[22].re * t6.re - X[22].im * t6.im);
  t18.im = (X[2].re * t8.im + X[2].im * t8.re) -
           (X[7].re * t9.im + X[7].im * t9.re) +
           (X[22].re * t6.im + X[22].im * t6.re);

  t19.re = (X[2].re * t2.re - X[2].im * t2.im) -
           (X[7].re * t4.re - X[7].im * t4.im) +
           (X[17].re * t6.re - X[17].im * t6.im);
  t19.im = (X[2].re * t2.im + X[2].im * t2.re) -
           (X[7].re * t4.im + X[7].im * t4.re) +
           (X[17].re * t6.im + X[17].im * t6.re);

  t20.re = (X[2].re * t3.re - X[2].im * t3.im) -
           (X[7].re * t5.re - X[7].im * t5.im) +
           (X[12].re * t6.re - X[12].im * t6.im);
  t20.im = (X[2].re * t3.im + X[2].im * t3.re) -
           (X[7].re * t5.im + X[7].im * t5.re) +
           (X[12].re * t6.im + X[12].im * t6.re);

  Y[17].re = (X[4].re * t12.re - X[4].im * t12.im) -
             (X[9].re * t15.re - X[9].im * t15.im) +
             (X[19].re * t18.re - X[19].im * t18.im) -
             (X[24].re * t19.re - X[24].im * t19.im);
  Y[17].im = (X[4].re * t12.im + X[4].im * t12.re) -
             (X[9].re * t15.im + X[9].im * t15.re) +
             (X[19].re * t18.im + X[19].im * t18.re) -
             (X[24].re * t19.im + X[24].im * t19.re);

  Y[22].re = -(X[3].re * t12.re - X[3].im * t12.im) +
             (X[8].re * t15.re - X[8].im * t15.im) -
             (X[18].re * t18.re - X[18].im * t18.im) +
             (X[23].re * t19.re - X[23].im * t19.im);
  Y[22].im = -(X[3].re * t12.im + X[3].im * t12.re) +
             (X[8].re * t15.im + X[8].im * t15.re) -
             (X[18].re * t18.im + X[18].im * t18.re) +
             (X[23].re * t19.im + X[23].im * t19.re);

  Y[18].re = -(X[4].re * t13.re - X[4].im * t13.im) +
             (X[9].re * t16.re - X[9].im * t16.im) -
             (X[14].re * t18.re - X[14].im * t18.im) +
             (X[24].re * t20.re - X[24].im * t20.im);
  Y[18].im = -(X[4].re * t13.im + X[4].im * t13.re) +
             (X[9].re * t16.im + X[9].im * t16.re) -
             (X[14].re * t18.im + X[14].im * t18.re) +
             (X[24].re * t20.im + X[24].im * t20.re);

  Y[23].re = (X[3].re * t13.re - X[3].im * t13.im) -
             (X[8].re * t16.re - X[8].im * t16.im) +
             (X[13].re * t18.re - X[13].im * t18.im) -
             (X[23].re * t20.re - X[23].im * t20.im);
  Y[23].im = (X[3].re * t13.im + X[3].im * t13.re) -
             (X[8].re * t16.im + X[8].im * t16.re) +
             (X[13].re * t18.im + X[13].im * t18.re) -
             (X[23].re * t20.im + X[23].im * t20.re);

  Y[19].re = (X[4].re * t14.re - X[4].im * t14.im) -
             (X[9].re * t17.re - X[9].im * t17.im) +
             (X[14].re * t19.re - X[14].im * t19.im) -
             (X[19].re * t20.re - X[19].im * t20.im);
  Y[19].im = (X[4].re * t14.im + X[4].im * t14.re) -
             (X[9].re * t17.im + X[9].im * t17.re) +
             (X[14].re * t19.im + X[14].im * t19.re) -
             (X[19].re * t20.im + X[19].im * t20.re);

  Y[24].re = -(X[3].re * t14.re - X[3].im * t14.im) +
             (X[8].re * t17.re - X[8].im * t17.im) -
             (X[13].re * t19.re - X[13].im * t19.im) +
             (X[18].re * t20.re - X[18].im * t20.im);
  Y[24].im = -(X[3].re * t14.im + X[3].im * t14.re) +
             (X[8].re * t17.im + X[8].im * t17.re) -
             (X[13].re * t19.im + X[13].im * t19.re) +
             (X[18].re * t20.im + X[18].im * t20.re);

  t11.re = (X[8].re * t7.re - X[8].im * t7.im) -
           (X[13].re * t8.re - X[13].im * t8.im) +
           (X[23].re * t3.re - X[23].im * t3.im);
  t11.im = (X[8].re * t7.im + X[8].im * t7.re) -
           (X[13].re * t8.im + X[13].im * t8.re) +
           (X[23].re * t3.im + X[23].im * t3.re);

  t12.re = (X[3].re * t7.re - X[3].im * t7.im) -
           (X[13].re * t9.re - X[13].im * t9.im) +
           (X[23].re * t5.re - X[23].im * t5.im);
  t12.im = (X[3].re * t7.im + X[3].im * t7.re) -
           (X[13].re * t9.im + X[13].im * t9.re) +
           (X[23].re * t5.im + X[23].im * t5.re);

  t13.re = (X[3].re * t8.re - X[3].im * t8.im) -
           (X[8].re * t9.re - X[8].im * t9.im) +
           (X[23].re * t6.re - X[23].im * t6.im);
  t13.im = (X[3].re * t8.im + X[3].im * t8.re) -
           (X[8].re * t9.im + X[8].im * t9.re) +
           (X[23].re * t6.im + X[23].im * t6.re);

  t14.re = (X[3].re * t3.re - X[3].im * t3.im) -
           (X[8].re * t5.re - X[8].im * t5.im) +
           (X[13].re * t6.re - X[13].im * t6.im);
  t14.im = (X[3].re * t3.im + X[3].im * t3.re) -
           (X[8].re * t5.im + X[8].im * t5.re) +
           (X[13].re * t6.im + X[13].im * t6.re);

  t15.re = (X[8].re * t1.re - X[8].im * t1.im) -
           (X[13].re * t2.re - X[13].im * t2.im) +
           (X[18].re * t3.re - X[18].im * t3.im);
  t15.im = (X[8].re * t1.im + X[8].im * t1.re) -
           (X[13].re * t2.im + X[13].im * t2.re) +
           (X[18].re * t3.im + X[18].im * t3.re);

  t16.re = (X[3].re * t1.re - X[3].im * t1.im) -
           (X[13].re * t4.re - X[13].im * t4.im) +
           (X[18].re * t5.re - X[18].im * t5.im);
  t16.im = (X[3].re * t1.im + X[3].im * t1.re) -
           (X[13].re * t4.im + X[13].im * t4.re) +
           (X[18].re * t5.im + X[18].im * t5.re);

  t17.re = (X[3].re * t2.re - X[3].im * t2.im) -
           (X[8].re * t4.re - X[8].im * t4.im) +
           (X[18].re * t6.re - X[18].im * t6.im);
  t17.im = (X[3].re * t2.im + X[3].im * t2.re) -
           (X[8].re * t4.im + X[8].im * t4.re) +
           (X[18].re * t6.im + X[18].im * t6.re);

  Y[13].re = (X[4].re * t11.re - X[4].im * t11.im) -
             (X[9].re * t12.re - X[9].im * t12.im) +
             (X[14].re * t13.re - X[14].im * t13.im) -
             (X[24].re * t14.re - X[24].im * t14.im);
  Y[13].im = (X[4].re * t11.im + X[4].im * t11.re) -
             (X[9].re * t12.im + X[9].im * t12.re) +
             (X[14].re * t13.im + X[14].im * t13.re) -
             (X[24].re * t14.im + X[24].im * t14.re);

  Y[14].re = -(X[4].re * t15.re - X[4].im * t15.im) +
             (X[9].re * t16.re - X[9].im * t16.im) -
             (X[14].re * t17.re - X[14].im * t17.im) +
             (X[19].re * t14.re - X[19].im * t14.im);
  Y[14].im = -(X[4].re * t15.im + X[4].im * t15.re) +
             (X[9].re * t16.im + X[9].im * t16.re) -
             (X[14].re * t17.im + X[14].im * t17.re) +
             (X[19].re * t14.im + X[19].im * t14.re);

  det.re = (det.re) / (det.re * det.re + det.im * det.im);
  det.im = -(det.im) / (det.re * det.re + det.im * det.im);

  t = Y[0].re;
  Y[0].re = (Y[0].re * det.re - Y[0].im * det.im);
  Y[0].im = (t * det.im + Y[0].im * det.re);

  t = Y[1].re;
  Y[1].re = (Y[1].re * det.re - Y[1].im * det.im);
  Y[1].im = (t * det.im + Y[1].im * det.re);

  t = Y[2].re;
  Y[2].re = (Y[2].re * det.re - Y[2].im * det.im);
  Y[2].im = (t * det.im + Y[2].im * det.re);

  t = Y[3].re;
  Y[3].re = (Y[3].re * det.re - Y[3].im * det.im);
  Y[3].im = (t * det.im + Y[3].im * det.re);

  t = Y[4].re;
  Y[4].re = (Y[4].re * det.re - Y[4].im * det.im);
  Y[4].im = (t * det.im + Y[4].im * det.re);

  t = Y[5].re;
  Y[5].re = (Y[5].re * det.re - Y[5].im * det.im);
  Y[5].im = (t * det.im + Y[5].im * det.re);

  t = Y[6].re;
  Y[6].re = (Y[6].re * det.re - Y[6].im * det.im);
  Y[6].im = (t * det.im + Y[6].im * det.re);

  t = Y[7].re;
  Y[7].re = (Y[7].re * det.re - Y[7].im * det.im);
  Y[7].im = (t * det.im + Y[7].im * det.re);

  t = Y[8].re;
  Y[8].re = (Y[8].re * det.re - Y[8].im * det.im);
  Y[8].im = (t * det.im + Y[8].im * det.re);

  t = Y[9].re;
  Y[9].re = (Y[9].re * det.re - Y[9].im * det.im);
  Y[9].im = (t * det.im + Y[9].im * det.re);

  t = Y[10].re;
  Y[10].re = (Y[10].re * det.re - Y[10].im * det.im);
  Y[10].im = (t * det.im + Y[10].im * det.re);

  t = Y[11].re;
  Y[11].re = (Y[11].re * det.re - Y[11].im * det.im);
  Y[11].im = (t * det.im + Y[11].im * det.re);

  t = Y[12].re;
  Y[12].re = (Y[12].re * det.re - Y[12].im * det.im);
  Y[12].im = (t * det.im + Y[12].im * det.re);

  t = Y[13].re;
  Y[13].re = (Y[13].re * det.re - Y[13].im * det.im);
  Y[13].im = (t * det.im + Y[13].im * det.re);

  t = Y[14].re;
  Y[14].re = (Y[14].re * det.re - Y[14].im * det.im);
  Y[14].im = (t * det.im + Y[14].im * det.re);

  t = Y[15].re;
  Y[15].re = (Y[15].re * det.re - Y[15].im * det.im);
  Y[15].im = (t * det.im + Y[15].im * det.re);

  t = Y[16].re;
  Y[16].re = (Y[16].re * det.re - Y[16].im * det.im);
  Y[16].im = (t * det.im + Y[16].im * det.re);

  t = Y[17].re;
  Y[17].re = (Y[17].re * det.re - Y[17].im * det.im);
  Y[17].im = (t * det.im + Y[17].im * det.re);

  t = Y[18].re;
  Y[18].re = (Y[18].re * det.re - Y[18].im * det.im);
  Y[18].im = (t * det.im + Y[18].im * det.re);

  t = Y[19].re;
  Y[19].re = (Y[19].re * det.re - Y[19].im * det.im);
  Y[19].im = (t * det.im + Y[19].im * det.re);

  t = Y[20].re;
  Y[20].re = (Y[20].re * det.re - Y[20].im * det.im);
  Y[20].im = (t * det.im + Y[20].im * det.re);

  t = Y[21].re;
  Y[21].re = (Y[21].re * det.re - Y[21].im * det.im);
  Y[21].im = (t * det.im + Y[21].im * det.re);

  t = Y[22].re;
  Y[22].re = (Y[22].re * det.re - Y[22].im * det.im);
  Y[22].im = (t * det.im + Y[22].im * det.re);

  t = Y[23].re;
  Y[23].re = (Y[23].re * det.re - Y[23].im * det.im);
  Y[23].im = (t * det.im + Y[23].im * det.re);

  t = Y[24].re;
  Y[24].re = (Y[24].re * det.re - Y[24].im * det.im);
  Y[24].im = (t * det.im + Y[24].im * det.re);
}

void invert_6by6(DTYPE* X, DTYPE* Y) {
  DTYPE det;
  DTYPE t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16,
      t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t29, t30, t31,
      t32, t33, t34, t35, t36, t37, t38, t39, t40, t41, t42, t43, t44, t45, t46,
      t47, t48, t49, t50;
#if DEBUG
  printf("%s\n", __func__);
#endif

  t1 = X[28] * X[35] - X[34] * X[29];
  t2 = X[22] * X[35] - X[34] * X[23];
  t3 = X[22] * X[29] - X[28] * X[23];
  t4 = X[16] * X[35] - X[34] * X[17];
  t5 = X[16] * X[29] - X[28] * X[17];
  t6 = X[16] * X[23] - X[22] * X[17];
  t7 = X[10] * X[35] - X[34] * X[11];
  t8 = X[10] * X[29] - X[28] * X[11];
  t9 = X[10] * X[23] - X[22] * X[11];
  t10 = X[10] * X[17] - X[16] * X[11];
  t11 = X[4] * X[35] - X[34] * X[5];
  t12 = X[4] * X[29] - X[28] * X[5];
  t13 = X[4] * X[23] - X[22] * X[5];
  t14 = X[4] * X[17] - X[16] * X[5];
  t15 = X[4] * X[11] - X[10] * X[5];

  t16 = X[21] * t1 - X[27] * t2 + X[33] * t3;
  t17 = X[15] * t1 - X[27] * t4 + X[33] * t5;
  t18 = X[15] * t2 - X[21] * t4 + X[33] * t6;
  t19 = X[15] * t3 - X[21] * t5 + X[27] * t6;
  t20 = X[9] * t1 - X[27] * t7 + X[33] * t8;
  t21 = X[9] * t2 - X[21] * t7 + X[33] * t9;
  t22 = X[9] * t3 - X[21] * t8 + X[27] * t9;
  t23 = X[9] * t4 - X[15] * t7 + X[33] * t10;
  t24 = X[9] * t5 - X[15] * t8 + X[27] * t10;
  t25 = X[9] * t6 - X[15] * t9 + X[21] * t10;
  t26 = X[3] * t1 - X[27] * t11 + X[33] * t12;
  t27 = X[3] * t2 - X[21] * t11 + X[33] * t13;
  t28 = X[3] * t3 - X[21] * t12 + X[27] * t13;
  t29 = X[3] * t4 - X[15] * t11 + X[33] * t14;
  t30 = X[3] * t5 - X[15] * t12 + X[27] * t14;
  t31 = X[3] * t6 - X[15] * t13 + X[21] * t14;
  t32 = X[3] * t7 - X[9] * t11 + X[33] * t15;
  t33 = X[3] * t8 - X[9] * t12 + X[27] * t15;
  t34 = X[3] * t9 - X[9] * t13 + X[21] * t15;
  t35 = X[3] * t10 - X[9] * t14 + X[15] * t15;

  t36 = X[14] * t16 - X[20] * t17 + X[26] * t18 - X[32] * t19;
  t37 = X[8] * t16 - X[20] * t20 + X[26] * t21 - X[32] * t22;
  t38 = X[8] * t17 - X[14] * t20 + X[26] * t23 - X[32] * t24;
  t39 = X[8] * t18 - X[14] * t21 + X[20] * t23 - X[32] * t25;
  t40 = X[8] * t19 - X[14] * t22 + X[20] * t24 - X[26] * t25;
  t41 = X[2] * t16 - X[20] * t26 + X[26] * t27 - X[32] * t28;
  t42 = X[2] * t17 - X[14] * t26 + X[26] * t29 - X[32] * t30;
  t43 = X[2] * t18 - X[14] * t27 + X[20] * t29 - X[32] * t31;
  t44 = X[2] * t19 - X[14] * t28 + X[20] * t30 - X[26] * t31;

  Y[0] = X[7] * t36 - X[13] * t37 + X[19] * t38 - X[25] * t39 + X[31] * t40;
  Y[6] = -X[6] * t36 + X[12] * t37 - X[18] * t38 + X[24] * t39 - X[30] * t40;
  Y[1] = -X[1] * t36 + X[13] * t41 - X[19] * t42 + X[25] * t43 - X[31] * t44;
  Y[7] = X[0] * t36 - X[12] * t41 + X[18] * t42 - X[24] * t43 + X[30] * t44;

  t45 = X[2] * t20 - X[8] * t26 + X[26] * t32 - X[32] * t33;
  t46 = X[2] * t21 - X[8] * t27 + X[20] * t32 - X[32] * t34;
  t47 = X[2] * t22 - X[8] * t28 + X[20] * t33 - X[26] * t34;
  t48 = X[2] * t23 - X[8] * t29 + X[14] * t32 - X[32] * t35;
  t49 = X[2] * t24 - X[8] * t30 + X[14] * t33 - X[26] * t35;

  Y[2] = X[1] * t37 - X[7] * t41 + X[19] * t45 - X[25] * t46 + X[31] * t47;
  Y[8] = -X[0] * t37 + X[6] * t41 - X[18] * t45 + X[24] * t46 - X[30] * t47;
  Y[3] = -X[1] * t38 + X[7] * t42 - X[13] * t45 + X[25] * t48 - X[31] * t49;
  Y[9] = X[0] * t38 - X[6] * t42 + X[12] * t45 - X[24] * t48 + X[30] * t49;

  t50 = X[2] * t25 - X[8] * t31 + X[14] * t34 - X[20] * t35;

  Y[4] = X[1] * t39 - X[7] * t43 + X[13] * t46 - X[19] * t48 + X[31] * t50;
  Y[10] = -X[0] * t39 + X[6] * t43 - X[12] * t46 + X[18] * t48 - X[30] * t50;
  Y[5] = -X[1] * t40 + X[7] * t44 - X[13] * t47 + X[19] * t49 - X[25] * t50;
  Y[11] = X[0] * t40 - X[6] * t44 + X[12] * t47 - X[18] * t49 + X[24] * t50;

  det = X[0] * Y[0] + X[6] * Y[1] + X[12] * Y[2] + X[18] * Y[3] + X[24] * Y[4] +
        X[30] * Y[5];

  if (det > -FZERO && det < FZERO) ASSERT(0, "This matrix is singural.\n")

  t36 = X[13] * t16 - X[19] * t17 + X[25] * t18 - X[31] * t19;
  t37 = X[7] * t16 - X[19] * t20 + X[25] * t21 - X[31] * t22;
  t38 = X[7] * t17 - X[13] * t20 + X[25] * t23 - X[31] * t24;
  t39 = X[7] * t18 - X[13] * t21 + X[19] * t23 - X[31] * t25;
  t40 = X[7] * t19 - X[13] * t22 + X[19] * t24 - X[25] * t25;
  t41 = X[1] * t16 - X[19] * t26 + X[25] * t27 - X[31] * t28;
  t42 = X[1] * t17 - X[13] * t26 + X[25] * t29 - X[31] * t30;
  t43 = X[1] * t18 - X[13] * t27 + X[19] * t29 - X[31] * t31;
  t44 = X[1] * t19 - X[13] * t28 + X[19] * t30 - X[25] * t31;
  t45 = X[1] * t20 - X[7] * t26 + X[25] * t32 - X[31] * t33;
  t46 = X[1] * t21 - X[7] * t27 + X[19] * t32 - X[31] * t34;
  t47 = X[1] * t22 - X[7] * t28 + X[19] * t33 - X[25] * t34;
  t48 = X[1] * t23 - X[7] * t29 + X[13] * t32 - X[31] * t35;
  t49 = X[1] * t24 - X[7] * t30 + X[13] * t33 - X[25] * t35;
  t50 = X[1] * t25 - X[7] * t31 + X[13] * t34 - X[19] * t35;

  Y[12] = X[6] * t36 - X[12] * t37 + X[18] * t38 - X[24] * t39 + X[30] * t40;
  Y[13] = -X[0] * t36 + X[12] * t41 - X[18] * t42 + X[24] * t43 - X[30] * t44;
  Y[14] = X[0] * t37 - X[6] * t41 + X[18] * t45 - X[24] * t46 + X[30] * t47;
  Y[15] = -X[0] * t38 + X[6] * t42 - X[12] * t45 + X[24] * t48 - X[30] * t49;
  Y[16] = X[0] * t39 - X[6] * t43 + X[12] * t46 - X[18] * t48 + X[30] * t50;
  Y[17] = -X[0] * t40 + X[6] * t44 - X[12] * t47 + X[18] * t49 - X[24] * t50;

  t1 = X[18] * X[25] - X[24] * X[19];
  t2 = X[12] * X[25] - X[24] * X[13];
  t3 = X[12] * X[19] - X[18] * X[13];
  t4 = X[6] * X[25] - X[24] * X[7];
  t5 = X[6] * X[19] - X[18] * X[7];
  t6 = X[6] * X[13] - X[12] * X[7];
  t7 = X[0] * X[25] - X[24] * X[1];
  t8 = X[0] * X[19] - X[18] * X[1];
  t9 = X[0] * X[13] - X[12] * X[1];
  t10 = X[0] * X[7] - X[6] * X[1];
  t11 = X[18] * X[31] - X[30] * X[19];
  t12 = X[12] * X[31] - X[30] * X[13];
  t13 = X[6] * X[31] - X[30] * X[7];
  t14 = X[0] * X[31] - X[30] * X[1];
  t15 = X[24] * X[31] - X[30] * X[25];

  t16 = X[20] * t15 - X[26] * t11 + X[32] * t1;
  t17 = X[14] * t15 - X[26] * t12 + X[32] * t2;
  t18 = X[14] * t11 - X[20] * t12 + X[32] * t3;
  t19 = X[14] * t1 - X[20] * t2 + X[26] * t3;
  t20 = X[8] * t15 - X[26] * t13 + X[32] * t4;
  t21 = X[8] * t11 - X[20] * t13 + X[32] * t5;
  t22 = X[8] * t1 - X[20] * t4 + X[26] * t5;
  t23 = X[8] * t12 - X[14] * t13 + X[32] * t6;
  t24 = X[8] * t2 - X[14] * t4 + X[26] * t6;
  t25 = X[8] * t3 - X[14] * t5 + X[20] * t6;
  t26 = X[2] * t15 - X[26] * t14 + X[32] * t7;
  t27 = X[2] * t11 - X[20] * t14 + X[32] * t8;
  t28 = X[2] * t1 - X[20] * t7 + X[26] * t8;
  t29 = X[2] * t12 - X[14] * t14 + X[32] * t9;
  t30 = X[2] * t2 - X[14] * t7 + X[26] * t9;
  t31 = X[2] * t3 - X[14] * t8 + X[20] * t9;
  t32 = X[2] * t13 - X[8] * t14 + X[32] * t10;
  t33 = X[2] * t4 - X[8] * t7 + X[26] * t10;
  t34 = X[2] * t5 - X[8] * t8 + X[20] * t10;
  t35 = X[2] * t6 - X[8] * t9 + X[14] * t10;

  t36 = X[15] * t16 - X[21] * t17 + X[27] * t18 - X[33] * t19;
  t37 = X[9] * t16 - X[21] * t20 + X[27] * t21 - X[33] * t22;
  t38 = X[9] * t17 - X[15] * t20 + X[27] * t23 - X[33] * t24;
  t39 = X[9] * t18 - X[15] * t21 + X[21] * t23 - X[33] * t25;
  t40 = X[9] * t19 - X[15] * t22 + X[21] * t24 - X[27] * t25;
  t41 = X[3] * t16 - X[21] * t26 + X[27] * t27 - X[33] * t28;
  t42 = X[3] * t17 - X[15] * t26 + X[27] * t29 - X[33] * t30;
  t43 = X[3] * t18 - X[15] * t27 + X[21] * t29 - X[33] * t31;
  t44 = X[3] * t19 - X[15] * t28 + X[21] * t30 - X[27] * t31;

  Y[24] = -X[11] * t36 + X[17] * t37 - X[23] * t38 + X[29] * t39 - X[35] * t40;
  Y[30] = X[10] * t36 - X[16] * t37 + X[22] * t38 - X[28] * t39 + X[34] * t40;
  Y[25] = X[5] * t36 - X[17] * t41 + X[23] * t42 - X[29] * t43 + X[35] * t44;
  Y[31] = -X[4] * t36 + X[16] * t41 - X[22] * t42 + X[28] * t43 -

          X[34] * t44;

  t45 = X[3] * t20 - X[9] * t26 + X[27] * t32 - X[33] * t33;
  t46 = X[3] * t21 - X[9] * t27 + X[21] * t32 - X[33] * t34;
  t47 = X[3] * t22 - X[9] * t28 + X[21] * t33 - X[27] * t34;
  t48 = X[3] * t23 - X[9] * t29 + X[15] * t32 - X[33] * t35;
  t49 = X[3] * t24 - X[9] * t30 + X[15] * t33 - X[27] * t35;

  Y[26] = -X[5] * t37 + X[11] * t41 - X[23] * t45 + X[29] * t46 - X[35] * t47;
  Y[32] = X[4] * t37 - X[10] * t41 + X[22] * t45 - X[28] * t46 + X[34] * t47;
  Y[27] = X[5] * t38 - X[11] * t42 + X[17] * t45 - X[29] * t48 + X[35] * t49;
  Y[33] = -X[4] * t38 + X[10] * t42 - X[16] * t45 + X[28] * t48 - X[34] * t49;

  t50 = X[3] * t25 - X[9] * t31 + X[15] * t34 - X[21] * t35;

  Y[28] = -X[5] * t39 + X[11] * t43 - X[17] * t46 + X[23] * t48 - X[35] * t50;
  Y[34] = X[4] * t39 - X[10] * t43 + X[16] * t46 - X[22] * t48 + X[34] * t50;
  Y[29] = X[5] * t40 - X[11] * t44 + X[17] * t47 - X[23] * t49 + X[29] * t50;
  Y[35] = -X[4] * t40 + X[10] * t44 - X[16] * t47 + X[22] * t49 - X[28] * t50;

  t36 = X[16] * t16 - X[22] * t17 + X[28] * t18 - X[34] * t19;
  t37 = X[10] * t16 - X[22] * t20 + X[28] * t21 - X[34] * t22;
  t38 = X[10] * t17 - X[16] * t20 + X[28] * t23 - X[34] * t24;
  t39 = X[10] * t18 - X[16] * t21 + X[22] * t23 - X[34] * t25;
  t40 = X[10] * t19 - X[16] * t22 + X[22] * t24 - X[28] * t25;
  t41 = X[4] * t16 - X[22] * t26 + X[28] * t27 - X[34] * t28;
  t42 = X[4] * t17 - X[16] * t26 + X[28] * t29 - X[34] * t30;
  t43 = X[4] * t18 - X[16] * t27 + X[22] * t29 - X[34] * t31;
  t44 = X[4] * t19 - X[16] * t28 + X[22] * t30 - X[28] * t31;
  t45 = X[4] * t20 - X[10] * t26 + X[28] * t32 - X[34] * t33;
  t46 = X[4] * t21 - X[10] * t27 + X[22] * t32 - X[34] * t34;
  t47 = X[4] * t22 - X[10] * t28 + X[22] * t33 - X[28] * t34;
  t48 = X[4] * t23 - X[10] * t29 + X[16] * t32 - X[34] * t35;
  t49 = X[4] * t24 - X[10] * t30 + X[16] * t33 - X[28] * t35;
  t50 = X[4] * t25 - X[10] * t31 + X[16] * t34 - X[22] * t35;

  Y[18] = X[11] * t36 - X[17] * t37 + X[23] * t38 - X[29] * t39 + X[35] * t40;
  Y[19] = -X[5] * t36 + X[17] * t41 - X[23] * t42 + X[29] * t43 - X[35] * t44;
  Y[20] = X[5] * t37 - X[11] * t41 + X[23] * t45 - X[29] * t46 + X[35] * t47;
  Y[21] = -X[5] * t38 + X[11] * t42 - X[17] * t45 + X[29] * t48 - X[35] * t49;
  Y[22] = X[5] * t39 - X[11] * t43 + X[17] * t46 - X[23] * t48 + X[35] * t50;
  Y[23] = -X[5] * t40 + X[11] * t44 - X[17] * t47 + X[23] * t49 - X[29] * t50;

  det = 1 / det;

  Y[0] = Y[0] * det;
  Y[1] = Y[1] * det;
  Y[2] = Y[2] * det;
  Y[3] = Y[3] * det;
  Y[4] = Y[4] * det;
  Y[5] = Y[5] * det;
  Y[6] = Y[6] * det;
  Y[7] = Y[7] * det;
  Y[8] = Y[8] * det;
  Y[9] = Y[9] * det;
  Y[10] = Y[10] * det;
  Y[11] = Y[11] * det;
  Y[12] = Y[12] * det;
  Y[13] = Y[13] * det;
  Y[14] = Y[14] * det;
  Y[15] = Y[15] * det;
  Y[16] = Y[16] * det;
  Y[17] = Y[17] * det;
  Y[18] = Y[18] * det;
  Y[19] = Y[19] * det;
  Y[20] = Y[20] * det;
  Y[21] = Y[21] * det;
  Y[22] = Y[22] * det;
  Y[23] = Y[23] * det;
  Y[24] = Y[24] * det;
  Y[25] = Y[25] * det;
  Y[26] = Y[26] * det;
  Y[27] = Y[27] * det;
  Y[28] = Y[28] * det;
  Y[29] = Y[29] * det;
  Y[30] = Y[30] * det;
  Y[31] = Y[31] * det;
  Y[32] = Y[32] * det;
  Y[33] = Y[33] * det;
  Y[34] = Y[34] * det;
  Y[35] = Y[35] * det;
}

void cinvert_6by6(CTYPE* X, CTYPE* Y) {
  CTYPE det;
  DTYPE t;
  CTYPE t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16,
      t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t29, t30, t31,
      t32, t33, t34, t35, t36, t37, t38, t39, t40, t41, t42, t43, t44, t45, t46,
      t47, t48, t49, t50;
#if DEBUG
  printf("%s\n", __func__);
#endif

  t1.re = (X[28].re * X[35].re - X[28].im * X[35].im) -
          (X[34].re * X[29].re - X[34].im * X[29].im);
  t1.im = (X[28].re * X[35].im + X[28].im * X[35].re) -
          (X[34].re * X[29].im + X[34].im * X[29].re);

  t2.re = (X[22].re * X[35].re - X[22].im * X[35].im) -
          (X[34].re * X[23].re - X[34].im * X[23].im);
  t2.im = (X[22].re * X[35].im + X[22].im * X[35].re) -
          (X[34].re * X[23].im + X[34].im * X[23].re);

  t3.re = (X[22].re * X[29].re - X[22].im * X[29].im) -
          (X[28].re * X[23].re - X[28].im * X[23].im);
  t3.im = (X[22].re * X[29].im + X[22].im * X[29].re) -
          (X[28].re * X[23].im + X[28].im * X[23].re);

  t4.re = (X[16].re * X[35].re - X[16].im * X[35].im) -
          (X[34].re * X[17].re - X[34].im * X[17].im);
  t4.im = (X[16].re * X[35].im + X[16].im * X[35].re) -
          (X[34].re * X[17].im + X[34].im * X[17].re);

  t5.re = (X[16].re * X[29].re - X[16].im * X[29].im) -
          (X[28].re * X[17].re - X[28].im * X[17].im);
  t5.im = (X[16].re * X[29].im + X[16].im * X[29].re) -
          (X[28].re * X[17].im + X[28].im * X[17].re);

  t6.re = (X[16].re * X[23].re - X[16].im * X[23].im) -
          (X[22].re * X[17].re - X[22].im * X[17].im);
  t6.im = (X[16].re * X[23].im + X[16].im * X[23].re) -
          (X[22].re * X[17].im + X[22].im * X[17].re);

  t7.re = (X[10].re * X[35].re - X[10].im * X[35].im) -
          (X[34].re * X[11].re - X[34].im * X[11].im);
  t7.im = (X[10].re * X[35].im + X[10].im * X[35].re) -
          (X[34].re * X[11].im + X[34].im * X[11].re);

  t8.re = (X[10].re * X[29].re - X[10].im * X[29].im) -
          (X[28].re * X[11].re - X[28].im * X[11].im);
  t8.im = (X[10].re * X[29].im + X[10].im * X[29].re) -
          (X[28].re * X[11].im + X[28].im * X[11].re);

  t9.re = (X[10].re * X[23].re - X[10].im * X[23].im) -
          (X[22].re * X[11].re - X[22].im * X[11].im);
  t9.im = (X[10].re * X[23].im + X[10].im * X[23].re) -
          (X[22].re * X[11].im + X[22].im * X[11].re);

  t10.re = (X[10].re * X[17].re - X[10].im * X[17].im) -
           (X[16].re * X[11].re - X[16].im * X[11].im);
  t10.im = (X[10].re * X[17].im + X[10].im * X[17].re) -
           (X[16].re * X[11].im + X[16].im * X[11].re);

  t11.re = (X[4].re * X[35].re - X[4].im * X[35].im) -
           (X[34].re * X[5].re - X[34].im * X[5].im);
  t11.im = (X[4].re * X[35].im + X[4].im * X[35].re) -
           (X[34].re * X[5].im + X[34].im * X[5].re);

  t12.re = (X[4].re * X[29].re - X[4].im * X[29].im) -
           (X[28].re * X[5].re - X[28].im * X[5].im);
  t12.im = (X[4].re * X[29].im + X[4].im * X[29].re) -
           (X[28].re * X[5].im + X[28].im * X[5].re);

  t13.re = (X[4].re * X[23].re - X[4].im * X[23].im) -
           (X[22].re * X[5].re - X[22].im * X[5].im);
  t13.im = (X[4].re * X[23].im + X[4].im * X[23].re) -
           (X[22].re * X[5].im + X[22].im * X[5].re);

  t14.re = (X[4].re * X[17].re - X[4].im * X[17].im) -
           (X[16].re * X[5].re - X[16].im * X[5].im);
  t14.im = (X[4].re * X[17].im + X[4].im * X[17].re) -
           (X[16].re * X[5].im + X[16].im * X[5].re);

  t15.re = (X[4].re * X[11].re - X[4].im * X[11].im) -
           (X[10].re * X[5].re - X[10].im * X[5].im);
  t15.im = (X[4].re * X[11].im + X[4].im * X[11].re) -
           (X[10].re * X[5].im + X[10].im * X[5].re);

  t16.re = (X[21].re * t1.re - X[21].im * t1.im) -
           (X[27].re * t2.re - X[27].im * t2.im) +
           (X[33].re * t3.re - X[33].im * t3.im);
  t16.im = (X[21].re * t1.im + X[21].im * t1.re) -
           (X[27].re * t2.im + X[27].im * t2.re) +
           (X[33].re * t3.im + X[33].im * t3.re);

  t17.re = (X[15].re * t1.re - X[15].im * t1.im) -
           (X[27].re * t4.re - X[27].im * t4.im) +
           (X[33].re * t5.re - X[33].im * t5.im);
  t17.im = (X[15].re * t1.im + X[15].im * t1.re) -
           (X[27].re * t4.im + X[27].im * t4.re) +
           (X[33].re * t5.im + X[33].im * t5.re);

  t18.re = (X[15].re * t2.re - X[15].im * t2.im) -
           (X[21].re * t4.re - X[21].im * t4.im) +
           (X[33].re * t6.re - X[33].im * t6.im);
  t18.im = (X[15].re * t2.im + X[15].im * t2.re) -
           (X[21].re * t4.im + X[21].im * t4.re) +
           (X[33].re * t6.im + X[33].im * t6.re);

  t19.re = (X[15].re * t3.re - X[15].im * t3.im) -
           (X[21].re * t5.re - X[21].im * t5.im) +
           (X[27].re * t6.re - X[27].im * t6.im);
  t19.im = (X[15].re * t3.im + X[15].im * t3.re) -
           (X[21].re * t5.im + X[21].im * t5.re) +
           (X[27].re * t6.im + X[27].im * t6.re);

  t20.re = (X[9].re * t1.re - X[9].im * t1.im) -
           (X[27].re * t7.re - X[27].im * t7.im) +
           (X[33].re * t8.re - X[33].im * t8.im);
  t20.im = (X[9].re * t1.im + X[9].im * t1.re) -
           (X[27].re * t7.im + X[27].im * t7.re) +
           (X[33].re * t8.im + X[33].im * t8.re);

  t21.re = (X[9].re * t2.re - X[9].im * t2.im) -
           (X[21].re * t7.re - X[21].im * t7.im) +
           (X[33].re * t9.re - X[33].im * t9.im);
  t21.im = (X[9].re * t2.im + X[9].im * t2.re) -
           (X[21].re * t7.im + X[21].im * t7.re) +
           (X[33].re * t9.im + X[33].im * t9.re);

  t22.re = (X[9].re * t3.re - X[9].im * t3.im) -
           (X[21].re * t8.re - X[21].im * t8.im) +
           (X[27].re * t9.re - X[27].im * t9.im);
  t22.im = (X[9].re * t3.im + X[9].im * t3.re) -
           (X[21].re * t8.im + X[21].im * t8.re) +
           (X[27].re * t9.im + X[27].im * t9.re);

  t23.re = (X[9].re * t4.re - X[9].im * t4.im) -
           (X[15].re * t7.re - X[15].im * t7.im) +
           (X[33].re * t10.re - X[33].im * t10.im);
  t23.im = (X[9].re * t4.im + X[9].im * t4.re) -
           (X[15].re * t7.im + X[15].im * t7.re) +
           (X[33].re * t10.im + X[33].im * t10.re);

  t24.re = (X[9].re * t5.re - X[9].im * t5.im) -
           (X[15].re * t8.re - X[15].im * t8.im) +
           (X[27].re * t10.re - X[27].im * t10.im);
  t24.im = (X[9].re * t5.im + X[9].im * t5.re) -
           (X[15].re * t8.im + X[15].im * t8.re) +
           (X[27].re * t10.im + X[27].im * t10.re);

  t25.re = (X[9].re * t6.re - X[9].im * t6.im) -
           (X[15].re * t9.re - X[15].im * t9.im) +
           (X[21].re * t10.re - X[21].im * t10.im);
  t25.im = (X[9].re * t6.im + X[9].im * t6.re) -
           (X[15].re * t9.im + X[15].im * t9.re) +
           (X[21].re * t10.im + X[21].im * t10.re);

  t26.re = (X[3].re * t1.re - X[3].im * t1.im) -
           (X[27].re * t11.re - X[27].im * t11.im) +
           (X[33].re * t12.re - X[33].im * t12.im);
  t26.im = (X[3].re * t1.im + X[3].im * t1.re) -
           (X[27].re * t11.im + X[27].im * t11.re) +
           (X[33].re * t12.im + X[33].im * t12.re);

  t27.re = (X[3].re * t2.re - X[3].im * t2.im) -
           (X[21].re * t11.re - X[21].im * t11.im) +
           (X[33].re * t13.re - X[33].im * t13.im);
  t27.im = (X[3].re * t2.im + X[3].im * t2.re) -
           (X[21].re * t11.im + X[21].im * t11.re) +
           (X[33].re * t13.im + X[33].im * t13.re);

  t28.re = (X[3].re * t3.re - X[3].im * t3.im) -
           (X[21].re * t12.re - X[21].im * t12.im) +
           (X[27].re * t13.re - X[27].im * t13.im);
  t28.im = (X[3].re * t3.im + X[3].im * t3.re) -
           (X[21].re * t12.im + X[21].im * t12.re) +
           (X[27].re * t13.im + X[27].im * t13.re);

  t29.re = (X[3].re * t4.re - X[3].im * t4.im) -
           (X[15].re * t11.re - X[15].im * t11.im) +
           (X[33].re * t14.re - X[33].im * t14.im);
  t29.im = (X[3].re * t4.im + X[3].im * t4.re) -
           (X[15].re * t11.im + X[15].im * t11.re) +
           (X[33].re * t14.im + X[33].im * t14.re);

  t30.re = (X[3].re * t5.re - X[3].im * t5.im) -
           (X[15].re * t12.re - X[15].im * t12.im) +
           (X[27].re * t14.re - X[27].im * t14.im);
  t30.im = (X[3].re * t5.im + X[3].im * t5.re) -
           (X[15].re * t12.im + X[15].im * t12.re) +
           (X[27].re * t14.im + X[27].im * t14.re);

  t31.re = (X[3].re * t6.re - X[3].im * t6.im) -
           (X[15].re * t13.re - X[15].im * t13.im) +
           (X[21].re * t14.re - X[21].im * t14.im);
  t31.im = (X[3].re * t6.im + X[3].im * t6.re) -
           (X[15].re * t13.im + X[15].im * t13.re) +
           (X[21].re * t14.im + X[21].im * t14.re);

  t32.re = (X[3].re * t7.re - X[3].im * t7.im) -
           (X[9].re * t11.re - X[9].im * t11.im) +
           (X[33].re * t15.re - X[33].im * t15.im);
  t32.im = (X[3].re * t7.im + X[3].im * t7.re) -
           (X[9].re * t11.im + X[9].im * t11.re) +
           (X[33].re * t15.im + X[33].im * t15.re);

  t33.re = (X[3].re * t8.re - X[3].im * t8.im) -
           (X[9].re * t12.re - X[9].im * t12.im) +
           (X[27].re * t15.re - X[27].im * t15.im);
  t33.im = (X[3].re * t8.im + X[3].im * t8.re) -
           (X[9].re * t12.im + X[9].im * t12.re) +
           (X[27].re * t15.im + X[27].im * t15.re);

  t34.re = (X[3].re * t9.re - X[3].im * t9.im) -
           (X[9].re * t13.re - X[9].im * t13.im) +
           (X[21].re * t15.re - X[21].im * t15.im);
  t34.im = (X[3].re * t9.im + X[3].im * t9.re) -
           (X[9].re * t13.im + X[9].im * t13.re) +
           (X[21].re * t15.im + X[21].im * t15.re);

  t35.re = (X[3].re * t10.re - X[3].im * t10.im) -
           (X[9].re * t14.re - X[9].im * t14.im) +
           (X[15].re * t15.re - X[15].im * t15.im);
  t35.im = (X[3].re * t10.im + X[3].im * t10.re) -
           (X[9].re * t14.im + X[9].im * t14.re) +
           (X[15].re * t15.im + X[15].im * t15.re);

  t36.re = (X[14].re * t16.re - X[14].im * t16.im) -
           (X[20].re * t17.re - X[20].im * t17.im) +
           (X[26].re * t18.re - X[26].im * t18.im) -
           (X[32].re * t19.re - X[32].im * t19.im);
  t36.im = (X[14].re * t16.im + X[14].im * t16.re) -
           (X[20].re * t17.im + X[20].im * t17.re) +
           (X[26].re * t18.im + X[26].im * t18.re) -
           (X[32].re * t19.im + X[32].im * t19.re);

  t37.re = (X[8].re * t16.re - X[8].im * t16.im) -
           (X[20].re * t20.re - X[20].im * t20.im) +
           (X[26].re * t21.re - X[26].im * t21.im) -
           (X[32].re * t22.re - X[32].im * t22.im);
  t37.im = (X[8].re * t16.im + X[8].im * t16.re) -
           (X[20].re * t20.im + X[20].im * t20.re) +
           (X[26].re * t21.im + X[26].im * t21.re) -
           (X[32].re * t22.im + X[32].im * t22.re);

  t38.re = (X[8].re * t17.re - X[8].im * t17.im) -
           (X[14].re * t20.re - X[14].im * t20.im) +
           (X[26].re * t23.re - X[26].im * t23.im) -
           (X[32].re * t24.re - X[32].im * t24.im);
  t38.im = (X[8].re * t17.im + X[8].im * t17.re) -
           (X[14].re * t20.im + X[14].im * t20.re) +
           (X[26].re * t23.im + X[26].im * t23.re) -
           (X[32].re * t24.im + X[32].im * t24.re);

  t39.re = (X[8].re * t18.re - X[8].im * t18.im) -
           (X[14].re * t21.re - X[14].im * t21.im) +
           (X[20].re * t23.re - X[20].im * t23.im) -
           (X[32].re * t25.re - X[32].im * t25.im);
  t39.im = (X[8].re * t18.im + X[8].im * t18.re) -
           (X[14].re * t21.im + X[14].im * t21.re) +
           (X[20].re * t23.im + X[20].im * t23.re) -
           (X[32].re * t25.im + X[32].im * t25.re);

  t40.re = (X[8].re * t19.re - X[8].im * t19.im) -
           (X[14].re * t22.re - X[14].im * t22.im) +
           (X[20].re * t24.re - X[20].im * t24.im) -
           (X[26].re * t25.re - X[26].im * t25.im);
  t40.im = (X[8].re * t19.im + X[8].im * t19.re) -
           (X[14].re * t22.im + X[14].im * t22.re) +
           (X[20].re * t24.im + X[20].im * t24.re) -
           (X[26].re * t25.im + X[26].im * t25.re);

  t41.re = (X[2].re * t16.re - X[2].im * t16.im) -
           (X[20].re * t26.re - X[20].im * t26.im) +
           (X[26].re * t27.re - X[26].im * t27.im) -
           (X[32].re * t28.re - X[32].im * t28.im);
  t41.im = (X[2].re * t16.im + X[2].im * t16.re) -
           (X[20].re * t26.im + X[20].im * t26.re) +
           (X[26].re * t27.im + X[26].im * t27.re) -
           (X[32].re * t28.im + X[32].im * t28.re);

  t42.re = (X[2].re * t17.re - X[2].im * t17.im) -
           (X[14].re * t26.re - X[14].im * t26.im) +
           (X[26].re * t29.re - X[26].im * t29.im) -
           (X[32].re * t30.re - X[32].im * t30.im);
  t42.im = (X[2].re * t17.im + X[2].im * t17.re) -
           (X[14].re * t26.im + X[14].im * t26.re) +
           (X[26].re * t29.im + X[26].im * t29.re) -
           (X[32].re * t30.im + X[32].im * t30.re);

  t43.re = (X[2].re * t18.re - X[2].im * t18.im) -
           (X[14].re * t27.re - X[14].im * t27.im) +
           (X[20].re * t29.re - X[20].im * t29.im) -
           (X[32].re * t31.re - X[32].im * t31.im);
  t43.im = (X[2].re * t18.im + X[2].im * t18.re) -
           (X[14].re * t27.im + X[14].im * t27.re) +
           (X[20].re * t29.im + X[20].im * t29.re) -
           (X[32].re * t31.im + X[32].im * t31.re);

  t44.re = (X[2].re * t19.re - X[2].im * t19.im) -
           (X[14].re * t28.re - X[14].im * t28.im) +
           (X[20].re * t30.re - X[20].im * t30.im) -
           (X[26].re * t31.re - X[26].im * t31.im);
  t44.im = (X[2].re * t19.im + X[2].im * t19.re) -
           (X[14].re * t28.im + X[14].im * t28.re) +
           (X[20].re * t30.im + X[20].im * t30.re) -
           (X[26].re * t31.im + X[26].im * t31.re);

  Y[0].re = (X[7].re * t36.re - X[7].im * t36.im) -
            (X[13].re * t37.re - X[13].im * t37.im) +
            (X[19].re * t38.re - X[19].im * t38.im) -
            (X[25].re * t39.re - X[25].im * t39.im) +
            (X[31].re * t40.re - X[31].im * t40.im);
  Y[0].im = (X[7].re * t36.im + X[7].im * t36.re) -
            (X[13].re * t37.im + X[13].im * t37.re) +
            (X[19].re * t38.im + X[19].im * t38.re) -
            (X[25].re * t39.im + X[25].im * t39.re) +
            (X[31].re * t40.im + X[31].im * t40.re);

  Y[6].re = -(X[6].re * t36.re - X[6].im * t36.im) +
            (X[12].re * t37.re - X[12].im * t37.im) -
            (X[18].re * t38.re - X[18].im * t38.im) +
            (X[24].re * t39.re - X[24].im * t39.im) -
            (X[30].re * t40.re - X[30].im * t40.im);
  Y[6].im = -(X[6].re * t36.im + X[6].im * t36.re) +
            (X[12].re * t37.im + X[12].im * t37.re) -
            (X[18].re * t38.im + X[18].im * t38.re) +
            (X[24].re * t39.im + X[24].im * t39.re) -
            (X[30].re * t40.im + X[30].im * t40.re);

  Y[1].re = -(X[1].re * t36.re - X[1].im * t36.im) +
            (X[13].re * t41.re - X[13].im * t41.im) -
            (X[19].re * t42.re - X[19].im * t42.im) +
            (X[25].re * t43.re - X[25].im * t43.im) -
            (X[31].re * t44.re - X[31].im * t44.im);
  Y[1].im = -(X[1].re * t36.im + X[1].im * t36.re) +
            (X[13].re * t41.im + X[13].im * t41.re) -
            (X[19].re * t42.im + X[19].im * t42.re) +
            (X[25].re * t43.im + X[25].im * t43.re) -
            (X[31].re * t44.im + X[31].im * t44.re);

  Y[7].re = (X[0].re * t36.re - X[0].im * t36.im) -
            (X[12].re * t41.re - X[12].im * t41.im) +
            (X[18].re * t42.re - X[18].im * t42.im) -
            (X[24].re * t43.re - X[24].im * t43.im) +
            (X[30].re * t44.re - X[30].im * t44.im);
  Y[7].im = (X[0].re * t36.im + X[0].im * t36.re) -
            (X[12].re * t41.im + X[12].im * t41.re) +
            (X[18].re * t42.im + X[18].im * t42.re) -
            (X[24].re * t43.im + X[24].im * t43.re) +
            (X[30].re * t44.im + X[30].im * t44.re);

  t45.re = (X[2].re * t20.re - X[2].im * t20.im) -
           (X[8].re * t26.re - X[8].im * t26.im) +
           (X[26].re * t32.re - X[26].im * t32.im) -
           (X[32].re * t33.re - X[32].im * t33.im);
  t45.im = (X[2].re * t20.im + X[2].im * t20.re) -
           (X[8].re * t26.im + X[8].im * t26.re) +
           (X[26].re * t32.im + X[26].im * t32.re) -
           (X[32].re * t33.im + X[32].im * t33.re);

  t46.re = (X[2].re * t21.re - X[2].im * t21.im) -
           (X[8].re * t27.re - X[8].im * t27.im) +
           (X[20].re * t32.re - X[20].im * t32.im) -
           (X[32].re * t34.re - X[32].im * t34.im);
  t46.im = (X[2].re * t21.im + X[2].im * t21.re) -
           (X[8].re * t27.im + X[8].im * t27.re) +
           (X[20].re * t32.im + X[20].im * t32.re) -
           (X[32].re * t34.im + X[32].im * t34.re);

  t47.re = (X[2].re * t22.re - X[2].im * t22.im) -
           (X[8].re * t28.re - X[8].im * t28.im) +
           (X[20].re * t33.re - X[20].im * t33.im) -
           (X[26].re * t34.re - X[26].im * t34.im);
  t47.im = (X[2].re * t22.im + X[2].im * t22.re) -
           (X[8].re * t28.im + X[8].im * t28.re) +
           (X[20].re * t33.im + X[20].im * t33.re) -
           (X[26].re * t34.im + X[26].im * t34.re);

  t48.re = (X[2].re * t23.re - X[2].im * t23.im) -
           (X[8].re * t29.re - X[8].im * t29.im) +
           (X[14].re * t32.re - X[14].im * t32.im) -
           (X[32].re * t35.re - X[32].im * t35.im);
  t48.im = (X[2].re * t23.im + X[2].im * t23.re) -
           (X[8].re * t29.im + X[8].im * t29.re) +
           (X[14].re * t32.im + X[14].im * t32.re) -
           (X[32].re * t35.im + X[32].im * t35.re);

  t49.re = (X[2].re * t24.re - X[2].im * t24.im) -
           (X[8].re * t30.re - X[8].im * t30.im) +
           (X[14].re * t33.re - X[14].im * t33.im) -
           (X[26].re * t35.re - X[26].im * t35.im);
  t49.im = (X[2].re * t24.im + X[2].im * t24.re) -
           (X[8].re * t30.im + X[8].im * t30.re) +
           (X[14].re * t33.im + X[14].im * t33.re) -
           (X[26].re * t35.im + X[26].im * t35.re);

  Y[2].re = (X[1].re * t37.re - X[1].im * t37.im) -
            (X[7].re * t41.re - X[7].im * t41.im) +
            (X[19].re * t45.re - X[19].im * t45.im) -
            (X[25].re * t46.re - X[25].im * t46.im) +
            (X[31].re * t47.re - X[31].im * t47.im);
  Y[2].im = (X[1].re * t37.im + X[1].im * t37.re) -
            (X[7].re * t41.im + X[7].im * t41.re) +
            (X[19].re * t45.im + X[19].im * t45.re) -
            (X[25].re * t46.im + X[25].im * t46.re) +
            (X[31].re * t47.im + X[31].im * t47.re);

  Y[8].re = -(X[0].re * t37.re - X[0].im * t37.im) +
            (X[6].re * t41.re - X[6].im * t41.im) -
            (X[18].re * t45.re - X[18].im * t45.im) +
            (X[24].re * t46.re - X[24].im * t46.im) -
            (X[30].re * t47.re - X[30].im * t47.im);
  Y[8].im = -(X[0].re * t37.im + X[0].im * t37.re) +
            (X[6].re * t41.im + X[6].im * t41.re) -
            (X[18].re * t45.im + X[18].im * t45.re) +
            (X[24].re * t46.im + X[24].im * t46.re) -
            (X[30].re * t47.im + X[30].im * t47.re);

  Y[3].re = -(X[1].re * t38.re - X[1].im * t38.im) +
            (X[7].re * t42.re - X[7].im * t42.im) -
            (X[13].re * t45.re - X[13].im * t45.im) +
            (X[25].re * t48.re - X[25].im * t48.im) -
            (X[31].re * t49.re - X[31].im * t49.im);
  Y[3].im = -(X[1].re * t38.im + X[1].im * t38.re) +
            (X[7].re * t42.im + X[7].im * t42.re) -
            (X[13].re * t45.im + X[13].im * t45.re) +
            (X[25].re * t48.im + X[25].im * t48.re) -
            (X[31].re * t49.im + X[31].im * t49.re);

  Y[9].re = (X[0].re * t38.re - X[0].im * t38.im) -
            (X[6].re * t42.re - X[6].im * t42.im) +
            (X[12].re * t45.re - X[12].im * t45.im) -
            (X[24].re * t48.re - X[24].im * t48.im) +
            (X[30].re * t49.re - X[30].im * t49.im);
  Y[9].im = (X[0].re * t38.im + X[0].im * t38.re) -
            (X[6].re * t42.im + X[6].im * t42.re) +
            (X[12].re * t45.im + X[12].im * t45.re) -
            (X[24].re * t48.im + X[24].im * t48.re) +
            (X[30].re * t49.im + X[30].im * t49.re);

  t50.re = (X[2].re * t25.re - X[2].im * t25.im) -
           (X[8].re * t31.re - X[8].im * t31.im) +
           (X[14].re * t34.re - X[14].im * t34.im) -
           (X[20].re * t35.re - X[20].im * t35.im);
  t50.im = (X[2].re * t25.im + X[2].im * t25.re) -
           (X[8].re * t31.im + X[8].im * t31.re) +
           (X[14].re * t34.im + X[14].im * t34.re) -
           (X[20].re * t35.im + X[20].im * t35.re);

  Y[4].re = (X[1].re * t39.re - X[1].im * t39.im) -
            (X[7].re * t43.re - X[7].im * t43.im) +
            (X[13].re * t46.re - X[13].im * t46.im) -
            (X[19].re * t48.re - X[19].im * t48.im) +
            (X[31].re * t50.re - X[31].im * t50.im);
  Y[4].im = (X[1].re * t39.im + X[1].im * t39.re) -
            (X[7].re * t43.im + X[7].im * t43.re) +
            (X[13].re * t46.im + X[13].im * t46.re) -
            (X[19].re * t48.im + X[19].im * t48.re) +
            (X[31].re * t50.im + X[31].im * t50.re);

  Y[10].re = -(X[0].re * t39.re - X[0].im * t39.im) +
             (X[6].re * t43.re - X[6].im * t43.im) -
             (X[12].re * t46.re - X[12].im * t46.im) +
             (X[18].re * t48.re - X[18].im * t48.im) -
             (X[30].re * t50.re - X[30].im * t50.im);
  Y[10].im = -(X[0].re * t39.im + X[0].im * t39.re) +
             (X[6].re * t43.im + X[6].im * t43.re) -
             (X[12].re * t46.im + X[12].im * t46.re) +
             (X[18].re * t48.im + X[18].im * t48.re) -
             (X[30].re * t50.im + X[30].im * t50.re);

  Y[5].re = -(X[1].re * t40.re - X[1].im * t40.im) +
            (X[7].re * t44.re - X[7].im * t44.im) -
            (X[13].re * t47.re - X[13].im * t47.im) +
            (X[19].re * t49.re - X[19].im * t49.im) -
            (X[25].re * t50.re - X[25].im * t50.im);
  Y[5].im = -(X[1].re * t40.im + X[1].im * t40.re) +
            (X[7].re * t44.im + X[7].im * t44.re) -
            (X[13].re * t47.im + X[13].im * t47.re) +
            (X[19].re * t49.im + X[19].im * t49.re) -
            (X[25].re * t50.im + X[25].im * t50.re);

  Y[11].re = (X[0].re * t40.re - X[0].im * t40.im) -
             (X[6].re * t44.re - X[6].im * t44.im) +
             (X[12].re * t47.re - X[12].im * t47.im) -
             (X[18].re * t49.re - X[18].im * t49.im) +
             (X[24].re * t50.re - X[24].im * t50.im);
  Y[11].im = (X[0].re * t40.im + X[0].im * t40.re) -
             (X[6].re * t44.im + X[6].im * t44.re) +
             (X[12].re * t47.im + X[12].im * t47.re) -
             (X[18].re * t49.im + X[18].im * t49.re) +
             (X[24].re * t50.im + X[24].im * t50.re);

  det.re = (X[0].re * Y[0].re - X[0].im * Y[0].im) +
           (X[6].re * Y[1].re - X[6].im * Y[1].im) +
           (X[12].re * Y[2].re - X[12].im * Y[2].im) +
           (X[18].re * Y[3].re - X[18].im * Y[3].im) +
           (X[24].re * Y[4].re - X[24].im * Y[4].im) +
           (X[30].re * Y[5].re - X[30].im * Y[5].im);
  det.im = (X[0].re * Y[0].im + X[0].im * Y[0].re) +
           (X[6].re * Y[1].im + X[6].im * Y[1].re) +
           (X[12].re * Y[2].im + X[12].im * Y[2].re) +
           (X[18].re * Y[3].im + X[18].im * Y[3].re) +
           (X[24].re * Y[4].im + X[24].im * Y[4].re) +
           (X[30].re * Y[5].im + X[30].im * Y[5].re);

#if NTYPE == 0
  if (cabsf(CXF(det)) < FZERO) ASSERT(0, "This matrix is singural.\n")
#elif NTYPE == 1
  if (cabs(CXD(det)) < FZERO) ASSERT(0, "This matrix is singural.\n")
#endif

  t36.re = (X[13].re * t16.re - X[13].im * t16.im) -
           (X[19].re * t17.re - X[19].im * t17.im) +
           (X[25].re * t18.re - X[25].im * t18.im) -
           (X[31].re * t19.re - X[31].im * t19.im);
  t36.im = (X[13].re * t16.im + X[13].im * t16.re) -
           (X[19].re * t17.im + X[19].im * t17.re) +
           (X[25].re * t18.im + X[25].im * t18.re) -
           (X[31].re * t19.im + X[31].im * t19.re);

  t37.re = (X[7].re * t16.re - X[7].im * t16.im) -
           (X[19].re * t20.re - X[19].im * t20.im) +
           (X[25].re * t21.re - X[25].im * t21.im) -
           (X[31].re * t22.re - X[31].im * t22.im);
  t37.im = (X[7].re * t16.im + X[7].im * t16.re) -
           (X[19].re * t20.im + X[19].im * t20.re) +
           (X[25].re * t21.im + X[25].im * t21.re) -
           (X[31].re * t22.im + X[31].im * t22.re);

  t38.re = (X[7].re * t17.re - X[7].im * t17.im) -
           (X[13].re * t20.re - X[13].im * t20.im) +
           (X[25].re * t23.re - X[25].im * t23.im) -
           (X[31].re * t24.re - X[31].im * t24.im);
  t38.im = (X[7].re * t17.im + X[7].im * t17.re) -
           (X[13].re * t20.im + X[13].im * t20.re) +
           (X[25].re * t23.im + X[25].im * t23.re) -
           (X[31].re * t24.im + X[31].im * t24.re);

  t39.re = (X[7].re * t18.re - X[7].im * t18.im) -
           (X[13].re * t21.re - X[13].im * t21.im) +
           (X[19].re * t23.re - X[19].im * t23.im) -
           (X[31].re * t25.re - X[31].im * t25.im);
  t39.im = (X[7].re * t18.im + X[7].im * t18.re) -
           (X[13].re * t21.im + X[13].im * t21.re) +
           (X[19].re * t23.im + X[19].im * t23.re) -
           (X[31].re * t25.im + X[31].im * t25.re);

  t40.re = (X[7].re * t19.re - X[7].im * t19.im) -
           (X[13].re * t22.re - X[13].im * t22.im) +
           (X[19].re * t24.re - X[19].im * t24.im) -
           (X[25].re * t25.re - X[25].im * t25.im);
  t40.im = (X[7].re * t19.im + X[7].im * t19.re) -
           (X[13].re * t22.im + X[13].im * t22.re) +
           (X[19].re * t24.im + X[19].im * t24.re) -
           (X[25].re * t25.im + X[25].im * t25.re);

  t41.re = (X[1].re * t16.re - X[1].im * t16.im) -
           (X[19].re * t26.re - X[19].im * t26.im) +
           (X[25].re * t27.re - X[25].im * t27.im) -
           (X[31].re * t28.re - X[31].im * t28.im);
  t41.im = (X[1].re * t16.im + X[1].im * t16.re) -
           (X[19].re * t26.im + X[19].im * t26.re) +
           (X[25].re * t27.im + X[25].im * t27.re) -
           (X[31].re * t28.im + X[31].im * t28.re);

  t42.re = (X[1].re * t17.re - X[1].im * t17.im) -
           (X[13].re * t26.re - X[13].im * t26.im) +
           (X[25].re * t29.re - X[25].im * t29.im) -
           (X[31].re * t30.re - X[31].im * t30.im);
  t42.im = (X[1].re * t17.im + X[1].im * t17.re) -
           (X[13].re * t26.im + X[13].im * t26.re) +
           (X[25].re * t29.im + X[25].im * t29.re) -
           (X[31].re * t30.im + X[31].im * t30.re);

  t43.re = (X[1].re * t18.re - X[1].im * t18.im) -
           (X[13].re * t27.re - X[13].im * t27.im) +
           (X[19].re * t29.re - X[19].im * t29.im) -
           (X[31].re * t31.re - X[31].im * t31.im);
  t43.im = (X[1].re * t18.im + X[1].im * t18.re) -
           (X[13].re * t27.im + X[13].im * t27.re) +
           (X[19].re * t29.im + X[19].im * t29.re) -
           (X[31].re * t31.im + X[31].im * t31.re);

  t44.re = (X[1].re * t19.re - X[1].im * t19.im) -
           (X[13].re * t28.re - X[13].im * t28.im) +
           (X[19].re * t30.re - X[19].im * t30.im) -
           (X[25].re * t31.re - X[25].im * t31.im);
  t44.im = (X[1].re * t19.im + X[1].im * t19.re) -
           (X[13].re * t28.im + X[13].im * t28.re) +
           (X[19].re * t30.im + X[19].im * t30.re) -
           (X[25].re * t31.im + X[25].im * t31.re);

  t45.re = (X[1].re * t20.re - X[1].im * t20.im) -
           (X[7].re * t26.re - X[7].im * t26.im) +
           (X[25].re * t32.re - X[25].im * t32.im) -
           (X[31].re * t33.re - X[31].im * t33.im);
  t45.im = (X[1].re * t20.im + X[1].im * t20.re) -
           (X[7].re * t26.im + X[7].im * t26.re) +
           (X[25].re * t32.im + X[25].im * t32.re) -
           (X[31].re * t33.im + X[31].im * t33.re);

  t46.re = (X[1].re * t21.re - X[1].im * t21.im) -
           (X[7].re * t27.re - X[7].im * t27.im) +
           (X[19].re * t32.re - X[19].im * t32.im) -
           (X[31].re * t34.re - X[31].im * t34.im);
  t46.im = (X[1].re * t21.im + X[1].im * t21.re) -
           (X[7].re * t27.im + X[7].im * t27.re) +
           (X[19].re * t32.im + X[19].im * t32.re) -
           (X[31].re * t34.im + X[31].im * t34.re);

  t47.re = (X[1].re * t22.re - X[1].im * t22.im) -
           (X[7].re * t28.re - X[7].im * t28.im) +
           (X[19].re * t33.re - X[19].im * t33.im) -
           (X[25].re * t34.re - X[25].im * t34.im);
  t47.im = (X[1].re * t22.im + X[1].im * t22.re) -
           (X[7].re * t28.im + X[7].im * t28.re) +
           (X[19].re * t33.im + X[19].im * t33.re) -
           (X[25].re * t34.im + X[25].im * t34.re);

  t48.re = (X[1].re * t23.re - X[1].im * t23.im) -
           (X[7].re * t29.re - X[7].im * t29.im) +
           (X[13].re * t32.re - X[13].im * t32.im) -
           (X[31].re * t35.re - X[31].im * t35.im);
  t48.im = (X[1].re * t23.im + X[1].im * t23.re) -
           (X[7].re * t29.im + X[7].im * t29.re) +
           (X[13].re * t32.im + X[13].im * t32.re) -
           (X[31].re * t35.im + X[31].im * t35.re);

  t49.re = (X[1].re * t24.re - X[1].im * t24.im) -
           (X[7].re * t30.re - X[7].im * t30.im) +
           (X[13].re * t33.re - X[13].im * t33.im) -
           (X[25].re * t35.re - X[25].im * t35.im);
  t49.im = (X[1].re * t24.im + X[1].im * t24.re) -
           (X[7].re * t30.im + X[7].im * t30.re) +
           (X[13].re * t33.im + X[13].im * t33.re) -
           (X[25].re * t35.im + X[25].im * t35.re);

  t50.re = (X[1].re * t25.re - X[1].im * t25.im) -
           (X[7].re * t31.re - X[7].im * t31.im) +
           (X[13].re * t34.re - X[13].im * t34.im) -
           (X[19].re * t35.re - X[19].im * t35.im);
  t50.im = (X[1].re * t25.im + X[1].im * t25.re) -
           (X[7].re * t31.im + X[7].im * t31.re) +
           (X[13].re * t34.im + X[13].im * t34.re) -
           (X[19].re * t35.im + X[19].im * t35.re);

  Y[12].re = (X[6].re * t36.re - X[6].im * t36.im) -
             (X[12].re * t37.re - X[12].im * t37.im) +
             (X[18].re * t38.re - X[18].im * t38.im) -
             (X[24].re * t39.re - X[24].im * t39.im) +
             (X[30].re * t40.re - X[30].im * t40.im);
  Y[12].im = (X[6].re * t36.im + X[6].im * t36.re) -
             (X[12].re * t37.im + X[12].im * t37.re) +
             (X[18].re * t38.im + X[18].im * t38.re) -
             (X[24].re * t39.im + X[24].im * t39.re) +
             (X[30].re * t40.im + X[30].im * t40.re);

  Y[13].re = -(X[0].re * t36.re - X[0].im * t36.im) +
             (X[12].re * t41.re - X[12].im * t41.im) -
             (X[18].re * t42.re - X[18].im * t42.im) +
             (X[24].re * t43.re - X[24].im * t43.im) -
             (X[30].re * t44.re - X[30].im * t44.im);
  Y[13].im = -(X[0].re * t36.im + X[0].im * t36.re) +
             (X[12].re * t41.im + X[12].im * t41.re) -
             (X[18].re * t42.im + X[18].im * t42.re) +
             (X[24].re * t43.im + X[24].im * t43.re) -
             (X[30].re * t44.im + X[30].im * t44.re);

  Y[14].re = (X[0].re * t37.re - X[0].im * t37.im) -
             (X[6].re * t41.re - X[6].im * t41.im) +
             (X[18].re * t45.re - X[18].im * t45.im) -
             (X[24].re * t46.re - X[24].im * t46.im) +
             (X[30].re * t47.re - X[30].im * t47.im);
  Y[14].im = (X[0].re * t37.im + X[0].im * t37.re) -
             (X[6].re * t41.im + X[6].im * t41.re) +
             (X[18].re * t45.im + X[18].im * t45.re) -
             (X[24].re * t46.im + X[24].im * t46.re) +
             (X[30].re * t47.im + X[30].im * t47.re);

  Y[15].re = -(X[0].re * t38.re - X[0].im * t38.im) +
             (X[6].re * t42.re - X[6].im * t42.im) -
             (X[12].re * t45.re - X[12].im * t45.im) +
             (X[24].re * t48.re - X[24].im * t48.im) -
             (X[30].re * t49.re - X[30].im * t49.im);
  Y[15].im = -(X[0].re * t38.im + X[0].im * t38.re) +
             (X[6].re * t42.im + X[6].im * t42.re) -
             (X[12].re * t45.im + X[12].im * t45.re) +
             (X[24].re * t48.im + X[24].im * t48.re) -
             (X[30].re * t49.im + X[30].im * t49.re);

  Y[16].re = (X[0].re * t39.re - X[0].im * t39.im) -
             (X[6].re * t43.re - X[6].im * t43.im) +
             (X[12].re * t46.re - X[12].im * t46.im) -
             (X[18].re * t48.re - X[18].im * t48.im) +
             (X[30].re * t50.re - X[30].im * t50.im);
  Y[16].im = (X[0].re * t39.im + X[0].im * t39.re) -
             (X[6].re * t43.im + X[6].im * t43.re) +
             (X[12].re * t46.im + X[12].im * t46.re) -
             (X[18].re * t48.im + X[18].im * t48.re) +
             (X[30].re * t50.im + X[30].im * t50.re);

  Y[17].re = -(X[0].re * t40.re - X[0].im * t40.im) +
             (X[6].re * t44.re - X[6].im * t44.im) -
             (X[12].re * t47.re - X[12].im * t47.im) +
             (X[18].re * t49.re - X[18].im * t49.im) -
             (X[24].re * t50.re - X[24].im * t50.im);
  Y[17].im = -(X[0].re * t40.im + X[0].im * t40.re) +
             (X[6].re * t44.im + X[6].im * t44.re) -
             (X[12].re * t47.im + X[12].im * t47.re) +
             (X[18].re * t49.im + X[18].im * t49.re) -
             (X[24].re * t50.im + X[24].im * t50.re);

  t1.re = (X[18].re * X[25].re - X[18].im * X[25].im) -
          (X[24].re * X[19].re - X[24].im * X[19].im);
  t1.im = (X[18].re * X[25].im + X[18].im * X[25].re) -
          (X[24].re * X[19].im + X[24].im * X[19].re);

  t2.re = (X[12].re * X[25].re - X[12].im * X[25].im) -
          (X[24].re * X[13].re - X[24].im * X[13].im);
  t2.im = (X[12].re * X[25].im + X[12].im * X[25].re) -
          (X[24].re * X[13].im + X[24].im * X[13].re);

  t3.re = (X[12].re * X[19].re - X[12].im * X[19].im) -
          (X[18].re * X[13].re - X[18].im * X[13].im);
  t3.im = (X[12].re * X[19].im + X[12].im * X[19].re) -
          (X[18].re * X[13].im + X[18].im * X[13].re);

  t4.re = (X[6].re * X[25].re - X[6].im * X[25].im) -
          (X[24].re * X[7].re - X[24].im * X[7].im);
  t4.im = (X[6].re * X[25].im + X[6].im * X[25].re) -
          (X[24].re * X[7].im + X[24].im * X[7].re);

  t5.re = (X[6].re * X[19].re - X[6].im * X[19].im) -
          (X[18].re * X[7].re - X[18].im * X[7].im);
  t5.im = (X[6].re * X[19].im + X[6].im * X[19].re) -
          (X[18].re * X[7].im + X[18].im * X[7].re);

  t6.re = (X[6].re * X[13].re - X[6].im * X[13].im) -
          (X[12].re * X[7].re - X[12].im * X[7].im);
  t6.im = (X[6].re * X[13].im + X[6].im * X[13].re) -
          (X[12].re * X[7].im + X[12].im * X[7].re);

  t7.re = (X[0].re * X[25].re - X[0].im * X[25].im) -
          (X[24].re * X[1].re - X[24].im * X[1].im);
  t7.im = (X[0].re * X[25].im + X[0].im * X[25].re) -
          (X[24].re * X[1].im + X[24].im * X[1].re);

  t8.re = (X[0].re * X[19].re - X[0].im * X[19].im) -
          (X[18].re * X[1].re - X[18].im * X[1].im);
  t8.im = (X[0].re * X[19].im + X[0].im * X[19].re) -
          (X[18].re * X[1].im + X[18].im * X[1].re);

  t9.re = (X[0].re * X[13].re - X[0].im * X[13].im) -
          (X[12].re * X[1].re - X[12].im * X[1].im);
  t9.im = (X[0].re * X[13].im + X[0].im * X[13].re) -
          (X[12].re * X[1].im + X[12].im * X[1].re);

  t10.re = (X[0].re * X[7].re - X[0].im * X[7].im) -
           (X[6].re * X[1].re - X[6].im * X[1].im);
  t10.im = (X[0].re * X[7].im + X[0].im * X[7].re) -
           (X[6].re * X[1].im + X[6].im * X[1].re);

  t11.re = (X[18].re * X[31].re - X[18].im * X[31].im) -
           (X[30].re * X[19].re - X[30].im * X[19].im);
  t11.im = (X[18].re * X[31].im + X[18].im * X[31].re) -
           (X[30].re * X[19].im + X[30].im * X[19].re);

  t12.re = (X[12].re * X[31].re - X[12].im * X[31].im) -
           (X[30].re * X[13].re - X[30].im * X[13].im);
  t12.im = (X[12].re * X[31].im + X[12].im * X[31].re) -
           (X[30].re * X[13].im + X[30].im * X[13].re);

  t13.re = (X[6].re * X[31].re - X[6].im * X[31].im) -
           (X[30].re * X[7].re - X[30].im * X[7].im);
  t13.im = (X[6].re * X[31].im + X[6].im * X[31].re) -
           (X[30].re * X[7].im + X[30].im * X[7].re);

  t14.re = (X[0].re * X[31].re - X[0].im * X[31].im) -
           (X[30].re * X[1].re - X[30].im * X[1].im);
  t14.im = (X[0].re * X[31].im + X[0].im * X[31].re) -
           (X[30].re * X[1].im + X[30].im * X[1].re);

  t15.re = (X[24].re * X[31].re - X[24].im * X[31].im) -
           (X[30].re * X[25].re - X[30].im * X[25].im);
  t15.im = (X[24].re * X[31].im + X[24].im * X[31].re) -
           (X[30].re * X[25].im + X[30].im * X[25].re);

  t16.re = (X[20].re * t15.re - X[20].im * t15.im) -
           (X[26].re * t11.re - X[26].im * t11.im) +
           (X[32].re * t1.re - X[32].im * t1.im);
  t16.im = (X[20].re * t15.im + X[20].im * t15.re) -
           (X[26].re * t11.im + X[26].im * t11.re) +
           (X[32].re * t1.im + X[32].im * t1.re);

  t17.re = (X[14].re * t15.re - X[14].im * t15.im) -
           (X[26].re * t12.re - X[26].im * t12.im) +
           (X[32].re * t2.re - X[32].im * t2.im);
  t17.im = (X[14].re * t15.im + X[14].im * t15.re) -
           (X[26].re * t12.im + X[26].im * t12.re) +
           (X[32].re * t2.im + X[32].im * t2.re);

  t18.re = (X[14].re * t11.re - X[14].im * t11.im) -
           (X[20].re * t12.re - X[20].im * t12.im) +
           (X[32].re * t3.re - X[32].im * t3.im);
  t18.im = (X[14].re * t11.im + X[14].im * t11.re) -
           (X[20].re * t12.im + X[20].im * t12.re) +
           (X[32].re * t3.im + X[32].im * t3.re);

  t19.re = (X[14].re * t1.re - X[14].im * t1.im) -
           (X[20].re * t2.re - X[20].im * t2.im) +
           (X[26].re * t3.re - X[26].im * t3.im);
  t19.im = (X[14].re * t1.im + X[14].im * t1.re) -
           (X[20].re * t2.im + X[20].im * t2.re) +
           (X[26].re * t3.im + X[26].im * t3.re);

  t20.re = (X[8].re * t15.re - X[8].im * t15.im) -
           (X[26].re * t13.re - X[26].im * t13.im) +
           (X[32].re * t4.re - X[32].im * t4.im);
  t20.im = (X[8].re * t15.im + X[8].im * t15.re) -
           (X[26].re * t13.im + X[26].im * t13.re) +
           (X[32].re * t4.im + X[32].im * t4.re);

  t21.re = (X[8].re * t11.re - X[8].im * t11.im) -
           (X[20].re * t13.re - X[20].im * t13.im) +
           (X[32].re * t5.re - X[32].im * t5.im);
  t21.im = (X[8].re * t11.im + X[8].im * t11.re) -
           (X[20].re * t13.im + X[20].im * t13.re) +
           (X[32].re * t5.im + X[32].im * t5.re);

  t22.re = (X[8].re * t1.re - X[8].im * t1.im) -
           (X[20].re * t4.re - X[20].im * t4.im) +
           (X[26].re * t5.re - X[26].im * t5.im);
  t22.im = (X[8].re * t1.im + X[8].im * t1.re) -
           (X[20].re * t4.im + X[20].im * t4.re) +
           (X[26].re * t5.im + X[26].im * t5.re);

  t23.re = (X[8].re * t12.re - X[8].im * t12.im) -
           (X[14].re * t13.re - X[14].im * t13.im) +
           (X[32].re * t6.re - X[32].im * t6.im);
  t23.im = (X[8].re * t12.im + X[8].im * t12.re) -
           (X[14].re * t13.im + X[14].im * t13.re) +
           (X[32].re * t6.im + X[32].im * t6.re);

  t24.re = (X[8].re * t2.re - X[8].im * t2.im) -
           (X[14].re * t4.re - X[14].im * t4.im) +
           (X[26].re * t6.re - X[26].im * t6.im);
  t24.im = (X[8].re * t2.im + X[8].im * t2.re) -
           (X[14].re * t4.im + X[14].im * t4.re) +
           (X[26].re * t6.im + X[26].im * t6.re);

  t25.re = (X[8].re * t3.re - X[8].im * t3.im) -
           (X[14].re * t5.re - X[14].im * t5.im) +
           (X[20].re * t6.re - X[20].im * t6.im);
  t25.im = (X[8].re * t3.im + X[8].im * t3.re) -
           (X[14].re * t5.im + X[14].im * t5.re) +
           (X[20].re * t6.im + X[20].im * t6.re);

  t26.re = (X[2].re * t15.re - X[2].im * t15.im) -
           (X[26].re * t14.re - X[26].im * t14.im) +
           (X[32].re * t7.re - X[32].im * t7.im);
  t26.im = (X[2].re * t15.im + X[2].im * t15.re) -
           (X[26].re * t14.im + X[26].im * t14.re) +
           (X[32].re * t7.im + X[32].im * t7.re);

  t27.re = (X[2].re * t11.re - X[2].im * t11.im) -
           (X[20].re * t14.re - X[20].im * t14.im) +
           (X[32].re * t8.re - X[32].im * t8.im);
  t27.im = (X[2].re * t11.im + X[2].im * t11.re) -
           (X[20].re * t14.im + X[20].im * t14.re) +
           (X[32].re * t8.im + X[32].im * t8.re);

  t28.re = (X[2].re * t1.re - X[2].im * t1.im) -
           (X[20].re * t7.re - X[20].im * t7.im) +
           (X[26].re * t8.re - X[26].im * t8.im);
  t28.im = (X[2].re * t1.im + X[2].im * t1.re) -
           (X[20].re * t7.im + X[20].im * t7.re) +
           (X[26].re * t8.im + X[26].im * t8.re);

  t29.re = (X[2].re * t12.re - X[2].im * t12.im) -
           (X[14].re * t14.re - X[14].im * t14.im) +
           (X[32].re * t9.re - X[32].im * t9.im);
  t29.im = (X[2].re * t12.im + X[2].im * t12.re) -
           (X[14].re * t14.im + X[14].im * t14.re) +
           (X[32].re * t9.im + X[32].im * t9.re);

  t30.re = (X[2].re * t2.re - X[2].im * t2.im) -
           (X[14].re * t7.re - X[14].im * t7.im) +
           (X[26].re * t9.re - X[26].im * t9.im);
  t30.im = (X[2].re * t2.im + X[2].im * t2.re) -
           (X[14].re * t7.im + X[14].im * t7.re) +
           (X[26].re * t9.im + X[26].im * t9.re);

  t31.re = (X[2].re * t3.re - X[2].im * t3.im) -
           (X[14].re * t8.re - X[14].im * t8.im) +
           (X[20].re * t9.re - X[20].im * t9.im);
  t31.im = (X[2].re * t3.im + X[2].im * t3.re) -
           (X[14].re * t8.im + X[14].im * t8.re) +
           (X[20].re * t9.im + X[20].im * t9.re);

  t32.re = (X[2].re * t13.re - X[2].im * t13.im) -
           (X[8].re * t14.re - X[8].im * t14.im) +
           (X[32].re * t10.re - X[32].im * t10.im);
  t32.im = (X[2].re * t13.im + X[2].im * t13.re) -
           (X[8].re * t14.im + X[8].im * t14.re) +
           (X[32].re * t10.im + X[32].im * t10.re);

  t33.re = (X[2].re * t4.re - X[2].im * t4.im) -
           (X[8].re * t7.re - X[8].im * t7.im) +
           (X[26].re * t10.re - X[26].im * t10.im);
  t33.im = (X[2].re * t4.im + X[2].im * t4.re) -
           (X[8].re * t7.im + X[8].im * t7.re) +
           (X[26].re * t10.im + X[26].im * t10.re);

  t34.re = (X[2].re * t5.re - X[2].im * t5.im) -
           (X[8].re * t8.re - X[8].im * t8.im) +
           (X[20].re * t10.re - X[20].im * t10.im);
  t34.im = (X[2].re * t5.im + X[2].im * t5.re) -
           (X[8].re * t8.im + X[8].im * t8.re) +
           (X[20].re * t10.im + X[20].im * t10.re);

  t35.re = (X[2].re * t6.re - X[2].im * t6.im) -
           (X[8].re * t9.re - X[8].im * t9.im) +
           (X[14].re * t10.re - X[14].im * t10.im);
  t35.im = (X[2].re * t6.im + X[2].im * t6.re) -
           (X[8].re * t9.im + X[8].im * t9.re) +
           (X[14].re * t10.im + X[14].im * t10.re);

  t36.re = (X[15].re * t16.re - X[15].im * t16.im) -
           (X[21].re * t17.re - X[21].im * t17.im) +
           (X[27].re * t18.re - X[27].im * t18.im) -
           (X[33].re * t19.re - X[33].im * t19.im);
  t36.im = (X[15].re * t16.im + X[15].im * t16.re) -
           (X[21].re * t17.im + X[21].im * t17.re) +
           (X[27].re * t18.im + X[27].im * t18.re) -
           (X[33].re * t19.im + X[33].im * t19.re);

  t37.re = (X[9].re * t16.re - X[9].im * t16.im) -
           (X[21].re * t20.re - X[21].im * t20.im) +
           (X[27].re * t21.re - X[27].im * t21.im) -
           (X[33].re * t22.re - X[33].im * t22.im);
  t37.im = (X[9].re * t16.im + X[9].im * t16.re) -
           (X[21].re * t20.im + X[21].im * t20.re) +
           (X[27].re * t21.im + X[27].im * t21.re) -
           (X[33].re * t22.im + X[33].im * t22.re);

  t38.re = (X[9].re * t17.re - X[9].im * t17.im) -
           (X[15].re * t20.re - X[15].im * t20.im) +
           (X[27].re * t23.re - X[27].im * t23.im) -
           (X[33].re * t24.re - X[33].im * t24.im);
  t38.im = (X[9].re * t17.im + X[9].im * t17.re) -
           (X[15].re * t20.im + X[15].im * t20.re) +
           (X[27].re * t23.im + X[27].im * t23.re) -
           (X[33].re * t24.im + X[33].im * t24.re);

  t39.re = (X[9].re * t18.re - X[9].im * t18.im) -
           (X[15].re * t21.re - X[15].im * t21.im) +
           (X[21].re * t23.re - X[21].im * t23.im) -
           (X[33].re * t25.re - X[33].im * t25.im);
  t39.im = (X[9].re * t18.im + X[9].im * t18.re) -
           (X[15].re * t21.im + X[15].im * t21.re) +
           (X[21].re * t23.im + X[21].im * t23.re) -
           (X[33].re * t25.im + X[33].im * t25.re);

  t40.re = (X[9].re * t19.re - X[9].im * t19.im) -
           (X[15].re * t22.re - X[15].im * t22.im) +
           (X[21].re * t24.re - X[21].im * t24.im) -
           (X[27].re * t25.re - X[27].im * t25.im);
  t40.im = (X[9].re * t19.im + X[9].im * t19.re) -
           (X[15].re * t22.im + X[15].im * t22.re) +
           (X[21].re * t24.im + X[21].im * t24.re) -
           (X[27].re * t25.im + X[27].im * t25.re);

  t41.re = (X[3].re * t16.re - X[3].im * t16.im) -
           (X[21].re * t26.re - X[21].im * t26.im) +
           (X[27].re * t27.re - X[27].im * t27.im) -
           (X[33].re * t28.re - X[33].im * t28.im);
  t41.im = (X[3].re * t16.im + X[3].im * t16.re) -
           (X[21].re * t26.im + X[21].im * t26.re) +
           (X[27].re * t27.im + X[27].im * t27.re) -
           (X[33].re * t28.im + X[33].im * t28.re);

  t42.re = (X[3].re * t17.re - X[3].im * t17.im) -
           (X[15].re * t26.re - X[15].im * t26.im) +
           (X[27].re * t29.re - X[27].im * t29.im) -
           (X[33].re * t30.re - X[33].im * t30.im);
  t42.im = (X[3].re * t17.im + X[3].im * t17.re) -
           (X[15].re * t26.im + X[15].im * t26.re) +
           (X[27].re * t29.im + X[27].im * t29.re) -
           (X[33].re * t30.im + X[33].im * t30.re);

  t43.re = (X[3].re * t18.re - X[3].im * t18.im) -
           (X[15].re * t27.re - X[15].im * t27.im) +
           (X[21].re * t29.re - X[21].im * t29.im) -
           (X[33].re * t31.re - X[33].im * t31.im);
  t43.im = (X[3].re * t18.im + X[3].im * t18.re) -
           (X[15].re * t27.im + X[15].im * t27.re) +
           (X[21].re * t29.im + X[21].im * t29.re) -
           (X[33].re * t31.im + X[33].im * t31.re);

  t44.re = (X[3].re * t19.re - X[3].im * t19.im) -
           (X[15].re * t28.re - X[15].im * t28.im) +
           (X[21].re * t30.re - X[21].im * t30.im) -
           (X[27].re * t31.re - X[27].im * t31.im);
  t44.im = (X[3].re * t19.im + X[3].im * t19.re) -
           (X[15].re * t28.im + X[15].im * t28.re) +
           (X[21].re * t30.im + X[21].im * t30.re) -
           (X[27].re * t31.im + X[27].im * t31.re);

  Y[24].re = -(X[11].re * t36.re - X[11].im * t36.im) +
             (X[17].re * t37.re - X[17].im * t37.im) -
             (X[23].re * t38.re - X[23].im * t38.im) +
             (X[29].re * t39.re - X[29].im * t39.im) -
             (X[35].re * t40.re - X[35].im * t40.im);
  Y[24].im = -(X[11].re * t36.im + X[11].im * t36.re) +
             (X[17].re * t37.im + X[17].im * t37.re) -
             (X[23].re * t38.im + X[23].im * t38.re) +
             (X[29].re * t39.im + X[29].im * t39.re) -
             (X[35].re * t40.im + X[35].im * t40.re);

  Y[30].re = (X[10].re * t36.re - X[10].im * t36.im) -
             (X[16].re * t37.re - X[16].im * t37.im) +
             (X[22].re * t38.re - X[22].im * t38.im) -
             (X[28].re * t39.re - X[28].im * t39.im) +
             (X[34].re * t40.re - X[34].im * t40.im);
  Y[30].im = (X[10].re * t36.im + X[10].im * t36.re) -
             (X[16].re * t37.im + X[16].im * t37.re) +
             (X[22].re * t38.im + X[22].im * t38.re) -
             (X[28].re * t39.im + X[28].im * t39.re) +
             (X[34].re * t40.im + X[34].im * t40.re);

  Y[25].re = (X[5].re * t36.re - X[5].im * t36.im) -
             (X[17].re * t41.re - X[17].im * t41.im) +
             (X[23].re * t42.re - X[23].im * t42.im) -
             (X[29].re * t43.re - X[29].im * t43.im) +
             (X[35].re * t44.re - X[35].im * t44.im);
  Y[25].im = (X[5].re * t36.im + X[5].im * t36.re) -
             (X[17].re * t41.im + X[17].im * t41.re) +
             (X[23].re * t42.im + X[23].im * t42.re) -
             (X[29].re * t43.im + X[29].im * t43.re) +
             (X[35].re * t44.im + X[35].im * t44.re);

  Y[31].re = -(X[4].re * t36.re - X[4].im * t36.im) +
             (X[16].re * t41.re - X[16].im * t41.im) -
             (X[22].re * t42.re - X[22].im * t42.im) +
             (X[28].re * t43.re - X[28].im * t43.im) -
             (X[34].re * t44.re - X[34].im * t44.im);
  Y[31].im = -(X[4].re * t36.im + X[4].im * t36.re) +
             (X[16].re * t41.im + X[16].im * t41.re) -
             (X[22].re * t42.im + X[22].im * t42.re) +
             (X[28].re * t43.im + X[28].im * t43.re) -
             (X[34].re * t44.im + X[34].im * t44.re);

  t45.re = (X[3].re * t20.re - X[3].im * t20.im) -
           (X[9].re * t26.re - X[9].im * t26.im) +
           (X[27].re * t32.re - X[27].im * t32.im) -
           (X[33].re * t33.re - X[33].im * t33.im);
  t45.im = (X[3].re * t20.im + X[3].im * t20.re) -
           (X[9].re * t26.im + X[9].im * t26.re) +
           (X[27].re * t32.im + X[27].im * t32.re) -
           (X[33].re * t33.im + X[33].im * t33.re);

  t46.re = (X[3].re * t21.re - X[3].im * t21.im) -
           (X[9].re * t27.re - X[9].im * t27.im) +
           (X[21].re * t32.re - X[21].im * t32.im) -
           (X[33].re * t34.re - X[33].im * t34.im);
  t46.im = (X[3].re * t21.im + X[3].im * t21.re) -
           (X[9].re * t27.im + X[9].im * t27.re) +
           (X[21].re * t32.im + X[21].im * t32.re) -
           (X[33].re * t34.im + X[33].im * t34.re);

  t47.re = (X[3].re * t22.re - X[3].im * t22.im) -
           (X[9].re * t28.re - X[9].im * t28.im) +
           (X[21].re * t33.re - X[21].im * t33.im) -
           (X[27].re * t34.re - X[27].im * t34.im);
  t47.im = (X[3].re * t22.im + X[3].im * t22.re) -
           (X[9].re * t28.im + X[9].im * t28.re) +
           (X[21].re * t33.im + X[21].im * t33.re) -
           (X[27].re * t34.im + X[27].im * t34.re);

  t48.re = (X[3].re * t23.re - X[3].im * t23.im) -
           (X[9].re * t29.re - X[9].im * t29.im) +
           (X[15].re * t32.re - X[15].im * t32.im) -
           (X[33].re * t35.re - X[33].im * t35.im);
  t48.im = (X[3].re * t23.im + X[3].im * t23.re) -
           (X[9].re * t29.im + X[9].im * t29.re) +
           (X[15].re * t32.im + X[15].im * t32.re) -
           (X[33].re * t35.im + X[33].im * t35.re);

  t49.re = (X[3].re * t24.re - X[3].im * t24.im) -
           (X[9].re * t30.re - X[9].im * t30.im) +
           (X[15].re * t33.re - X[15].im * t33.im) -
           (X[27].re * t35.re - X[27].im * t35.im);
  t49.im = (X[3].re * t24.im + X[3].im * t24.re) -
           (X[9].re * t30.im + X[9].im * t30.re) +
           (X[15].re * t33.im + X[15].im * t33.re) -
           (X[27].re * t35.im + X[27].im * t35.re);

  Y[26].re = -(X[5].re * t37.re - X[5].im * t37.im) +
             (X[11].re * t41.re - X[11].im * t41.im) -
             (X[23].re * t45.re - X[23].im * t45.im) +
             (X[29].re * t46.re - X[29].im * t46.im) -
             (X[35].re * t47.re - X[35].im * t47.im);
  Y[26].im = -(X[5].re * t37.im + X[5].im * t37.re) +
             (X[11].re * t41.im + X[11].im * t41.re) -
             (X[23].re * t45.im + X[23].im * t45.re) +
             (X[29].re * t46.im + X[29].im * t46.re) -
             (X[35].re * t47.im + X[35].im * t47.re);

  Y[32].re = (X[4].re * t37.re - X[4].im * t37.im) -
             (X[10].re * t41.re - X[10].im * t41.im) +
             (X[22].re * t45.re - X[22].im * t45.im) -
             (X[28].re * t46.re - X[28].im * t46.im) +
             (X[34].re * t47.re - X[34].im * t47.im);
  Y[32].im = (X[4].re * t37.im + X[4].im * t37.re) -
             (X[10].re * t41.im + X[10].im * t41.re) +
             (X[22].re * t45.im + X[22].im * t45.re) -
             (X[28].re * t46.im + X[28].im * t46.re) +
             (X[34].re * t47.im + X[34].im * t47.re);

  Y[27].re = (X[5].re * t38.re - X[5].im * t38.im) -
             (X[11].re * t42.re - X[11].im * t42.im) +
             (X[17].re * t45.re - X[17].im * t45.im) -
             (X[29].re * t48.re - X[29].im * t48.im) +
             (X[35].re * t49.re - X[35].im * t49.im);
  Y[27].im = (X[5].re * t38.im + X[5].im * t38.re) -
             (X[11].re * t42.im + X[11].im * t42.re) +
             (X[17].re * t45.im + X[17].im * t45.re) -
             (X[29].re * t48.im + X[29].im * t48.re) +
             (X[35].re * t49.im + X[35].im * t49.re);

  Y[33].re = -(X[4].re * t38.re - X[4].im * t38.im) +
             (X[10].re * t42.re - X[10].im * t42.im) -
             (X[16].re * t45.re - X[16].im * t45.im) +
             (X[28].re * t48.re - X[28].im * t48.im) -
             (X[34].re * t49.re - X[34].im * t49.im);
  Y[33].im = -(X[4].re * t38.im + X[4].im * t38.re) +
             (X[10].re * t42.im + X[10].im * t42.re) -
             (X[16].re * t45.im + X[16].im * t45.re) +
             (X[28].re * t48.im + X[28].im * t48.re) -
             (X[34].re * t49.im + X[34].im * t49.re);

  t50.re = (X[3].re * t25.re - X[3].im * t25.im) -
           (X[9].re * t31.re - X[9].im * t31.im) +
           (X[15].re * t34.re - X[15].im * t34.im) -
           (X[21].re * t35.re - X[21].im * t35.im);
  t50.im = (X[3].re * t25.im + X[3].im * t25.re) -
           (X[9].re * t31.im + X[9].im * t31.re) +
           (X[15].re * t34.im + X[15].im * t34.re) -
           (X[21].re * t35.im + X[21].im * t35.re);

  Y[28].re = -(X[5].re * t39.re - X[5].im * t39.im) +
             (X[11].re * t43.re - X[11].im * t43.im) -
             (X[17].re * t46.re - X[17].im * t46.im) +
             (X[23].re * t48.re - X[23].im * t48.im) -
             (X[35].re * t50.re - X[35].im * t50.im);
  Y[28].im = -(X[5].re * t39.im + X[5].im * t39.re) +
             (X[11].re * t43.im + X[11].im * t43.re) -
             (X[17].re * t46.im + X[17].im * t46.re) +
             (X[23].re * t48.im + X[23].im * t48.re) -
             (X[35].re * t50.im + X[35].im * t50.re);

  Y[34].re = (X[4].re * t39.re - X[4].im * t39.im) -
             (X[10].re * t43.re - X[10].im * t43.im) +
             (X[16].re * t46.re - X[16].im * t46.im) -
             (X[22].re * t48.re - X[22].im * t48.im) +
             (X[34].re * t50.re - X[34].im * t50.im);
  Y[34].im = (X[4].re * t39.im + X[4].im * t39.re) -
             (X[10].re * t43.im + X[10].im * t43.re) +
             (X[16].re * t46.im + X[16].im * t46.re) -
             (X[22].re * t48.im + X[22].im * t48.re) +
             (X[34].re * t50.im + X[34].im * t50.re);

  Y[29].re = (X[5].re * t40.re - X[5].im * t40.im) -
             (X[11].re * t44.re - X[11].im * t44.im) +
             (X[17].re * t47.re - X[17].im * t47.im) -
             (X[23].re * t49.re - X[23].im * t49.im) +
             (X[29].re * t50.re - X[29].im * t50.im);
  Y[29].im = (X[5].re * t40.im + X[5].im * t40.re) -
             (X[11].re * t44.im + X[11].im * t44.re) +
             (X[17].re * t47.im + X[17].im * t47.re) -
             (X[23].re * t49.im + X[23].im * t49.re) +
             (X[29].re * t50.im + X[29].im * t50.re);

  Y[35].re = -(X[4].re * t40.re - X[4].im * t40.im) +
             (X[10].re * t44.re - X[10].im * t44.im) -
             (X[16].re * t47.re - X[16].im * t47.im) +
             (X[22].re * t49.re - X[22].im * t49.im) -
             (X[28].re * t50.re - X[28].im * t50.im);
  Y[35].im = -(X[4].re * t40.im + X[4].im * t40.re) +
             (X[10].re * t44.im + X[10].im * t44.re) -
             (X[16].re * t47.im + X[16].im * t47.re) +
             (X[22].re * t49.im + X[22].im * t49.re) -
             (X[28].re * t50.im + X[28].im * t50.re);

  t36.re = (X[16].re * t16.re - X[16].im * t16.im) -
           (X[22].re * t17.re - X[22].im * t17.im) +
           (X[28].re * t18.re - X[28].im * t18.im) -
           (X[34].re * t19.re - X[34].im * t19.im);
  t36.im = (X[16].re * t16.im + X[16].im * t16.re) -
           (X[22].re * t17.im + X[22].im * t17.re) +
           (X[28].re * t18.im + X[28].im * t18.re) -
           (X[34].re * t19.im + X[34].im * t19.re);

  t37.re = (X[10].re * t16.re - X[10].im * t16.im) -
           (X[22].re * t20.re - X[22].im * t20.im) +
           (X[28].re * t21.re - X[28].im * t21.im) -
           (X[34].re * t22.re - X[34].im * t22.im);
  t37.im = (X[10].re * t16.im + X[10].im * t16.re) -
           (X[22].re * t20.im + X[22].im * t20.re) +
           (X[28].re * t21.im + X[28].im * t21.re) -
           (X[34].re * t22.im + X[34].im * t22.re);

  t38.re = (X[10].re * t17.re - X[10].im * t17.im) -
           (X[16].re * t20.re - X[16].im * t20.im) +
           (X[28].re * t23.re - X[28].im * t23.im) -
           (X[34].re * t24.re - X[34].im * t24.im);
  t38.im = (X[10].re * t17.im + X[10].im * t17.re) -
           (X[16].re * t20.im + X[16].im * t20.re) +
           (X[28].re * t23.im + X[28].im * t23.re) -
           (X[34].re * t24.im + X[34].im * t24.re);

  t39.re = (X[10].re * t18.re - X[10].im * t18.im) -
           (X[16].re * t21.re - X[16].im * t21.im) +
           (X[22].re * t23.re - X[22].im * t23.im) -
           (X[34].re * t25.re - X[34].im * t25.im);
  t39.im = (X[10].re * t18.im + X[10].im * t18.re) -
           (X[16].re * t21.im + X[16].im * t21.re) +
           (X[22].re * t23.im + X[22].im * t23.re) -
           (X[34].re * t25.im + X[34].im * t25.re);

  t40.re = (X[10].re * t19.re - X[10].im * t19.im) -
           (X[16].re * t22.re - X[16].im * t22.im) +
           (X[22].re * t24.re - X[22].im * t24.im) -
           (X[28].re * t25.re - X[28].im * t25.im);
  t40.im = (X[10].re * t19.im + X[10].im * t19.re) -
           (X[16].re * t22.im + X[16].im * t22.re) +
           (X[22].re * t24.im + X[22].im * t24.re) -
           (X[28].re * t25.im + X[28].im * t25.re);

  t41.re = (X[4].re * t16.re - X[4].im * t16.im) -
           (X[22].re * t26.re - X[22].im * t26.im) +
           (X[28].re * t27.re - X[28].im * t27.im) -
           (X[34].re * t28.re - X[34].im * t28.im);
  t41.im = (X[4].re * t16.im + X[4].im * t16.re) -
           (X[22].re * t26.im + X[22].im * t26.re) +
           (X[28].re * t27.im + X[28].im * t27.re) -
           (X[34].re * t28.im + X[34].im * t28.re);

  t42.re = (X[4].re * t17.re - X[4].im * t17.im) -
           (X[16].re * t26.re - X[16].im * t26.im) +
           (X[28].re * t29.re - X[28].im * t29.im) -
           (X[34].re * t30.re - X[34].im * t30.im);
  t42.im = (X[4].re * t17.im + X[4].im * t17.re) -
           (X[16].re * t26.im + X[16].im * t26.re) +
           (X[28].re * t29.im + X[28].im * t29.re) -
           (X[34].re * t30.im + X[34].im * t30.re);

  t43.re = (X[4].re * t18.re - X[4].im * t18.im) -
           (X[16].re * t27.re - X[16].im * t27.im) +
           (X[22].re * t29.re - X[22].im * t29.im) -
           (X[34].re * t31.re - X[34].im * t31.im);
  t43.im = (X[4].re * t18.im + X[4].im * t18.re) -
           (X[16].re * t27.im + X[16].im * t27.re) +
           (X[22].re * t29.im + X[22].im * t29.re) -
           (X[34].re * t31.im + X[34].im * t31.re);

  t44.re = (X[4].re * t19.re - X[4].im * t19.im) -
           (X[16].re * t28.re - X[16].im * t28.im) +
           (X[22].re * t30.re - X[22].im * t30.im) -
           (X[28].re * t31.re - X[28].im * t31.im);
  t44.im = (X[4].re * t19.im + X[4].im * t19.re) -
           (X[16].re * t28.im + X[16].im * t28.re) +
           (X[22].re * t30.im + X[22].im * t30.re) -
           (X[28].re * t31.im + X[28].im * t31.re);

  t45.re = (X[4].re * t20.re - X[4].im * t20.im) -
           (X[10].re * t26.re - X[10].im * t26.im) +
           (X[28].re * t32.re - X[28].im * t32.im) -
           (X[34].re * t33.re - X[34].im * t33.im);
  t45.im = (X[4].re * t20.im + X[4].im * t20.re) -
           (X[10].re * t26.im + X[10].im * t26.re) +
           (X[28].re * t32.im + X[28].im * t32.re) -
           (X[34].re * t33.im + X[34].im * t33.re);

  t46.re = (X[4].re * t21.re - X[4].im * t21.im) -
           (X[10].re * t27.re - X[10].im * t27.im) +
           (X[22].re * t32.re - X[22].im * t32.im) -
           (X[34].re * t34.re - X[34].im * t34.im);
  t46.im = (X[4].re * t21.im + X[4].im * t21.re) -
           (X[10].re * t27.im + X[10].im * t27.re) +
           (X[22].re * t32.im + X[22].im * t32.re) -
           (X[34].re * t34.im + X[34].im * t34.re);

  t47.re = (X[4].re * t22.re - X[4].im * t22.im) -
           (X[10].re * t28.re - X[10].im * t28.im) +
           (X[22].re * t33.re - X[22].im * t33.im) -
           (X[28].re * t34.re - X[28].im * t34.im);
  t47.im = (X[4].re * t22.im + X[4].im * t22.re) -
           (X[10].re * t28.im + X[10].im * t28.re) +
           (X[22].re * t33.im + X[22].im * t33.re) -
           (X[28].re * t34.im + X[28].im * t34.re);

  t48.re = (X[4].re * t23.re - X[4].im * t23.im) -
           (X[10].re * t29.re - X[10].im * t29.im) +
           (X[16].re * t32.re - X[16].im * t32.im) -
           (X[34].re * t35.re - X[34].im * t35.im);
  t48.im = (X[4].re * t23.im + X[4].im * t23.re) -
           (X[10].re * t29.im + X[10].im * t29.re) +
           (X[16].re * t32.im + X[16].im * t32.re) -
           (X[34].re * t35.im + X[34].im * t35.re);

  t49.re = (X[4].re * t24.re - X[4].im * t24.im) -
           (X[10].re * t30.re - X[10].im * t30.im) +
           (X[16].re * t33.re - X[16].im * t33.im) -
           (X[28].re * t35.re - X[28].im * t35.im);
  t49.im = (X[4].re * t24.im + X[4].im * t24.re) -
           (X[10].re * t30.im + X[10].im * t30.re) +
           (X[16].re * t33.im + X[16].im * t33.re) -
           (X[28].re * t35.im + X[28].im * t35.re);

  t50.re = (X[4].re * t25.re - X[4].im * t25.im) -
           (X[10].re * t31.re - X[10].im * t31.im) +
           (X[16].re * t34.re - X[16].im * t34.im) -
           (X[22].re * t35.re - X[22].im * t35.im);
  t50.im = (X[4].re * t25.im + X[4].im * t25.re) -
           (X[10].re * t31.im + X[10].im * t31.re) +
           (X[16].re * t34.im + X[16].im * t34.re) -
           (X[22].re * t35.im + X[22].im * t35.re);

  Y[18].re = (X[11].re * t36.re - X[11].im * t36.im) -
             (X[17].re * t37.re - X[17].im * t37.im) +
             (X[23].re * t38.re - X[23].im * t38.im) -
             (X[29].re * t39.re - X[29].im * t39.im) +
             (X[35].re * t40.re - X[35].im * t40.im);
  Y[18].im = (X[11].re * t36.im + X[11].im * t36.re) -
             (X[17].re * t37.im + X[17].im * t37.re) +
             (X[23].re * t38.im + X[23].im * t38.re) -
             (X[29].re * t39.im + X[29].im * t39.re) +
             (X[35].re * t40.im + X[35].im * t40.re);

  Y[19].re = -(X[5].re * t36.re - X[5].im * t36.im) +
             (X[17].re * t41.re - X[17].im * t41.im) -
             (X[23].re * t42.re - X[23].im * t42.im) +
             (X[29].re * t43.re - X[29].im * t43.im) -
             (X[35].re * t44.re - X[35].im * t44.im);
  Y[19].im = -(X[5].re * t36.im + X[5].im * t36.re) +
             (X[17].re * t41.im + X[17].im * t41.re) -
             (X[23].re * t42.im + X[23].im * t42.re) +
             (X[29].re * t43.im + X[29].im * t43.re) -
             (X[35].re * t44.im + X[35].im * t44.re);

  Y[20].re = (X[5].re * t37.re - X[5].im * t37.im) -
             (X[11].re * t41.re - X[11].im * t41.im) +
             (X[23].re * t45.re - X[23].im * t45.im) -
             (X[29].re * t46.re - X[29].im * t46.im) +
             (X[35].re * t47.re - X[35].im * t47.im);
  Y[20].im = (X[5].re * t37.im + X[5].im * t37.re) -
             (X[11].re * t41.im + X[11].im * t41.re) +
             (X[23].re * t45.im + X[23].im * t45.re) -
             (X[29].re * t46.im + X[29].im * t46.re) +
             (X[35].re * t47.im + X[35].im * t47.re);

  Y[21].re = -(X[5].re * t38.re - X[5].im * t38.im) +
             (X[11].re * t42.re - X[11].im * t42.im) -
             (X[17].re * t45.re - X[17].im * t45.im) +
             (X[29].re * t48.re - X[29].im * t48.im) -
             (X[35].re * t49.re - X[35].im * t49.im);
  Y[21].im = -(X[5].re * t38.im + X[5].im * t38.re) +
             (X[11].re * t42.im + X[11].im * t42.re) -
             (X[17].re * t45.im + X[17].im * t45.re) +
             (X[29].re * t48.im + X[29].im * t48.re) -
             (X[35].re * t49.im + X[35].im * t49.re);

  Y[22].re = (X[5].re * t39.re - X[5].im * t39.im) -
             (X[11].re * t43.re - X[11].im * t43.im) +
             (X[17].re * t46.re - X[17].im * t46.im) -
             (X[23].re * t48.re - X[23].im * t48.im) +
             (X[35].re * t50.re - X[35].im * t50.im);
  Y[22].im = (X[5].re * t39.im + X[5].im * t39.re) -
             (X[11].re * t43.im + X[11].im * t43.re) +
             (X[17].re * t46.im + X[17].im * t46.re) -
             (X[23].re * t48.im + X[23].im * t48.re) +
             (X[35].re * t50.im + X[35].im * t50.re);

  Y[23].re = -(X[5].re * t40.re - X[5].im * t40.im) +
             (X[11].re * t44.re - X[11].im * t44.im) -
             (X[17].re * t47.re - X[17].im * t47.im) +
             (X[23].re * t49.re - X[23].im * t49.im) -
             (X[29].re * t50.re - X[29].im * t50.im);
  Y[23].im = -(X[5].re * t40.im + X[5].im * t40.re) +
             (X[11].re * t44.im + X[11].im * t44.re) -
             (X[17].re * t47.im + X[17].im * t47.re) +
             (X[23].re * t49.im + X[23].im * t49.re) -
             (X[29].re * t50.im + X[29].im * t50.re);

  det.re = (det.re) / (det.re * det.re + det.im * det.im);
  det.im = -(det.im) / (det.re * det.re + det.im * det.im);

  t = Y[0].re;
  Y[0].re = (Y[0].re * det.re - Y[0].im * det.im);
  Y[0].im = (t * det.im + Y[0].im * det.re);

  t = Y[1].re;
  Y[1].re = (Y[1].re * det.re - Y[1].im * det.im);
  Y[1].im = (t * det.im + Y[1].im * det.re);

  t = Y[2].re;
  Y[2].re = (Y[2].re * det.re - Y[2].im * det.im);
  Y[2].im = (t * det.im + Y[2].im * det.re);

  t = Y[3].re;
  Y[3].re = (Y[3].re * det.re - Y[3].im * det.im);
  Y[3].im = (t * det.im + Y[3].im * det.re);

  t = Y[4].re;
  Y[4].re = (Y[4].re * det.re - Y[4].im * det.im);
  Y[4].im = (t * det.im + Y[4].im * det.re);

  t = Y[5].re;
  Y[5].re = (Y[5].re * det.re - Y[5].im * det.im);
  Y[5].im = (t * det.im + Y[5].im * det.re);

  t = Y[6].re;
  Y[6].re = (Y[6].re * det.re - Y[6].im * det.im);
  Y[6].im = (t * det.im + Y[6].im * det.re);

  t = Y[7].re;
  Y[7].re = (Y[7].re * det.re - Y[7].im * det.im);
  Y[7].im = (t * det.im + Y[7].im * det.re);

  t = Y[8].re;
  Y[8].re = (Y[8].re * det.re - Y[8].im * det.im);
  Y[8].im = (t * det.im + Y[8].im * det.re);

  t = Y[9].re;
  Y[9].re = (Y[9].re * det.re - Y[9].im * det.im);
  Y[9].im = (t * det.im + Y[9].im * det.re);

  t = Y[10].re;
  Y[10].re = (Y[10].re * det.re - Y[10].im * det.im);
  Y[10].im = (t * det.im + Y[10].im * det.re);

  t = Y[11].re;
  Y[11].re = (Y[11].re * det.re - Y[11].im * det.im);
  Y[11].im = (t * det.im + Y[11].im * det.re);

  t = Y[12].re;
  Y[12].re = (Y[12].re * det.re - Y[12].im * det.im);
  Y[12].im = (t * det.im + Y[12].im * det.re);

  t = Y[13].re;
  Y[13].re = (Y[13].re * det.re - Y[13].im * det.im);
  Y[13].im = (t * det.im + Y[13].im * det.re);

  t = Y[14].re;
  Y[14].re = (Y[14].re * det.re - Y[14].im * det.im);
  Y[14].im = (t * det.im + Y[14].im * det.re);

  t = Y[15].re;
  Y[15].re = (Y[15].re * det.re - Y[15].im * det.im);
  Y[15].im = (t * det.im + Y[15].im * det.re);

  t = Y[16].re;
  Y[16].re = (Y[16].re * det.re - Y[16].im * det.im);
  Y[16].im = (t * det.im + Y[16].im * det.re);

  t = Y[17].re;
  Y[17].re = (Y[17].re * det.re - Y[17].im * det.im);
  Y[17].im = (t * det.im + Y[17].im * det.re);

  t = Y[18].re;
  Y[18].re = (Y[18].re * det.re - Y[18].im * det.im);
  Y[18].im = (t * det.im + Y[18].im * det.re);

  t = Y[19].re;
  Y[19].re = (Y[19].re * det.re - Y[19].im * det.im);
  Y[19].im = (t * det.im + Y[19].im * det.re);

  t = Y[20].re;
  Y[20].re = (Y[20].re * det.re - Y[20].im * det.im);
  Y[20].im = (t * det.im + Y[20].im * det.re);

  t = Y[21].re;
  Y[21].re = (Y[21].re * det.re - Y[21].im * det.im);
  Y[21].im = (t * det.im + Y[21].im * det.re);

  t = Y[22].re;
  Y[22].re = (Y[22].re * det.re - Y[22].im * det.im);
  Y[22].im = (t * det.im + Y[22].im * det.re);

  t = Y[23].re;
  Y[23].re = (Y[23].re * det.re - Y[23].im * det.im);
  Y[23].im = (t * det.im + Y[23].im * det.re);

  t = Y[24].re;
  Y[24].re = (Y[24].re * det.re - Y[24].im * det.im);
  Y[24].im = (t * det.im + Y[24].im * det.re);

  t = Y[25].re;
  Y[25].re = (Y[25].re * det.re - Y[25].im * det.im);
  Y[25].im = (t * det.im + Y[25].im * det.re);

  t = Y[26].re;
  Y[26].re = (Y[26].re * det.re - Y[26].im * det.im);
  Y[26].im = (t * det.im + Y[26].im * det.re);

  t = Y[27].re;
  Y[27].re = (Y[27].re * det.re - Y[27].im * det.im);
  Y[27].im = (t * det.im + Y[27].im * det.re);

  t = Y[28].re;
  Y[28].re = (Y[28].re * det.re - Y[28].im * det.im);
  Y[28].im = (t * det.im + Y[28].im * det.re);

  t = Y[29].re;
  Y[29].re = (Y[29].re * det.re - Y[29].im * det.im);
  Y[29].im = (t * det.im + Y[29].im * det.re);

  t = Y[30].re;
  Y[30].re = (Y[30].re * det.re - Y[30].im * det.im);
  Y[30].im = (t * det.im + Y[30].im * det.re);

  t = Y[31].re;
  Y[31].re = (Y[31].re * det.re - Y[31].im * det.im);
  Y[31].im = (t * det.im + Y[31].im * det.re);

  t = Y[32].re;
  Y[32].re = (Y[32].re * det.re - Y[32].im * det.im);
  Y[32].im = (t * det.im + Y[32].im * det.re);

  t = Y[33].re;
  Y[33].re = (Y[33].re * det.re - Y[33].im * det.im);
  Y[33].im = (t * det.im + Y[33].im * det.re);

  t = Y[34].re;
  Y[34].re = (Y[34].re * det.re - Y[34].im * det.im);
  Y[34].im = (t * det.im + Y[34].im * det.re);

  t = Y[35].re;
  Y[35].re = (Y[35].re * det.re - Y[35].im * det.im);
  Y[35].im = (t * det.im + Y[35].im * det.re);
}

/**** get determinant of matrix ****/
DTYPE det(MAT* mat) {
  if (mat->d0 == 2)
    return det_2by2(mat->data);
  else if (mat->d0 == 3)
    return det_3by3(mat->data);
  else if (mat->d0 == 4)
    return det_4by4(mat->data);
  else if (mat->d0 == 5)
    return det_5by5(mat->data);
  else if (mat->d0 == 6)
    return det_6by6(mat->data);
  else if (mat->d0 > 6)
    return det_nbyn(mat->data);

  ASSERT_DIM_INVALID()
}

CTYPE cdet(CMAT* mat) {
  if (mat->d0 == 2)
    return cdet_2by2(mat->data);
  else if (mat->d0 == 3)
    return cdet_3by3(mat->data);
  else if (mat->d0 == 4)
    return cdet_4by4(mat->data);
  else if (mat->d0 == 5)
    return cdet_5by5(mat->data);
  else if (mat->d0 == 6)
    return cdet_6by6(mat->data);
  else if (mat->d0 > 6)
    return cdet_nbyn(mat->data);

  ASSERT_DIM_INVALID()
}

DTYPE det_2by2(DTYPE* X) {
  DTYPE det;
#if DEBUG
  printf("%s\n", __func__);
#endif
  det = X[0] * X[3] - X[1] * X[2];
  return det;
}

CTYPE cdet_2by2(CTYPE* X) {
  CTYPE det;
#if DEBUG
  printf("%s\n", __func__);
#endif
  det.re = (X[0].re * X[3].re - X[0].im * X[3].im) -
           (X[1].re * X[2].re - X[1].im * X[2].im);
  det.im = (X[0].re * X[3].im + X[0].im * X[3].re) -
           (X[1].re * X[2].im + X[1].im * X[2].re);

  return det;
}

DTYPE det_3by3(DTYPE* X) {
  DTYPE det;
  DTYPE t0, t1, t2;
#if DEBUG
  printf("%s\n", __func__);
#endif

  t0 = X[4] * X[8] - X[7] * X[5];
  t1 = X[7] * X[2] - X[1] * X[8];
  t2 = X[1] * X[5] - X[4] * X[2];

  det = X[0] * t0 + X[3] * t1 + X[6] * t2;

  return det;
}

CTYPE cdet_3by3(CTYPE* X) {
  CTYPE det;
  CTYPE t0, t1, t2;
#if DEBUG
  printf("%s\n", __func__);
#endif

  t0.re = (X[4].re * X[8].re - X[4].im * X[8].im) -
          (X[7].re * X[5].re - X[7].im * X[5].im);
  t0.im = (X[4].re * X[8].im + X[4].im * X[8].re) -
          (X[7].re * X[5].im + X[7].im * X[5].re);

  t1.re = (X[7].re * X[2].re - X[7].im * X[2].im) -
          (X[1].re * X[8].re - X[1].im * X[8].im);
  t1.im = (X[7].re * X[2].im + X[7].im * X[2].re) -
          (X[1].re * X[8].im + X[1].im * X[8].re);

  t2.re = (X[1].re * X[5].re - X[1].im * X[5].im) -
          (X[4].re * X[2].re - X[4].im * X[2].im);
  t2.im = (X[1].re * X[5].im + X[1].im * X[5].re) -
          (X[4].re * X[2].im + X[4].im * X[2].re);

  det.re = (X[0].re * t0.re - X[0].im * t0.im) +
           (X[3].re * t1.re - X[3].im * t1.im) +
           (X[6].re * t2.re - X[6].im * t2.im);
  det.im = (X[0].re * t0.im + X[0].im * t0.re) +
           (X[3].re * t1.im + X[3].im * t1.re) +
           (X[6].re * t2.im + X[6].im * t2.re);

  return det;
}

DTYPE det_4by4(DTYPE* X) {
  DTYPE det;
  DTYPE t1, t2, t3, t4, t5;
  DTYPE m0, m1, m2, m3;
#if DEBUG
  printf("%s\n", __func__);
#endif
  t1 = X[10] * X[15] - X[14] * X[11];
  t2 = X[6] * X[15] - X[14] * X[7];
  t3 = X[6] * X[11] - X[10] * X[7];

  m0 = X[5] * t1 - X[9] * t2 + X[13] * t3;

  t4 = X[2] * X[15] - X[14] * X[3];
  t5 = X[2] * X[11] - X[10] * X[3];

  m1 = X[9] * t4 - X[1] * t1 - X[13] * t5;

  t1 = X[2] * X[7] - X[6] * X[3];

  m2 = X[1] * t2 - X[5] * t4 + X[13] * t1;
  m3 = X[5] * t5 - X[1] * t3 - X[9] * t1;

  det = X[0] * m0 + X[4] * m1 + X[8] * m2 + X[12] * m3;
  return det;
}

CTYPE cdet_4by4(CTYPE* X) {
  CTYPE det;
  CTYPE t1, t2, t3, t4, t5;
  CTYPE m0, m1, m2, m3;
#if DEBUG
  printf("%s\n", __func__);
#endif

  t1.re = (X[10].re * X[15].re - X[10].im * X[15].im) -
          (X[14].re * X[11].re - X[14].im * X[11].im);
  t1.im = (X[10].re * X[15].im + X[10].im * X[15].re) -
          (X[14].re * X[11].im + X[14].im * X[11].re);

  t2.re = (X[6].re * X[15].re - X[6].im * X[15].im) -
          (X[14].re * X[7].re - X[14].im * X[7].im);
  t2.im = (X[6].re * X[15].im + X[6].im * X[15].re) -
          (X[14].re * X[7].im + X[14].im * X[7].re);

  t3.re = (X[6].re * X[11].re - X[6].im * X[11].im) -
          (X[10].re * X[7].re - X[10].im * X[7].im);
  t3.im = (X[6].re * X[11].im + X[6].im * X[11].re) -
          (X[10].re * X[7].im + X[10].im * X[7].re);

  m0.re = (X[5].re * t1.re - X[5].im * t1.im) -
          (X[9].re * t2.re - X[9].im * t2.im) +
          (X[13].re * t3.re - X[13].im * t3.im);
  m0.im = (X[5].re * t1.im + X[5].im * t1.re) -
          (X[9].re * t2.im + X[9].im * t2.re) +
          (X[13].re * t3.im + X[13].im * t3.re);

  t4.re = (X[2].re * X[15].re - X[2].im * X[15].im) -
          (X[14].re * X[3].re - X[14].im * X[3].im);
  t4.im = (X[2].re * X[15].im + X[2].im * X[15].re) -
          (X[14].re * X[3].im + X[14].im * X[3].re);

  t5.re = (X[2].re * X[11].re - X[2].im * X[11].im) -
          (X[10].re * X[3].re - X[10].im * X[3].im);
  t5.im = (X[2].re * X[11].im + X[2].im * X[11].re) -
          (X[10].re * X[3].im + X[10].im * X[3].re);

  m1.re = (X[9].re * t4.re - X[9].im * t4.im) -
          (X[1].re * t1.re - X[1].im * t1.im) -
          (X[13].re * t5.re - X[13].im * t5.im);
  m1.im = (X[9].re * t4.im + X[9].im * t4.re) -
          (X[1].re * t1.im + X[1].im * t1.re) -
          (X[13].re * t5.im + X[13].im * t5.re);

  t1.re = (X[2].re * X[7].re - X[2].im * X[7].im) -
          (X[6].re * X[3].re - X[6].im * X[3].im);
  t1.im = (X[2].re * X[7].im + X[2].im * X[7].re) -
          (X[6].re * X[3].im + X[6].im * X[3].re);

  m2.re = (X[1].re * t2.re - X[1].im * t2.im) -
          (X[5].re * t4.re - X[5].im * t4.im) +
          (X[13].re * t1.re - X[13].im * t1.im);
  m2.im = (X[1].re * t2.im + X[1].im * t2.re) -
          (X[5].re * t4.im + X[5].im * t4.re) +
          (X[13].re * t1.im + X[13].im * t1.re);

  m3.re = (X[5].re * t5.re - X[5].im * t5.im) -
          (X[1].re * t3.re - X[1].im * t3.im) -
          (X[9].re * t1.re - X[9].im * t1.im);
  m3.im = (X[5].re * t5.im + X[5].im * t5.re) -
          (X[1].re * t3.im + X[1].im * t3.re) -
          (X[9].re * t1.im + X[9].im * t1.re);

  det.re = (X[0].re * m0.re - X[0].im * m0.im) +
           (X[4].re * m1.re - X[4].im * m1.im) +
           (X[8].re * m2.re - X[8].im * m2.im) +
           (X[12].re * m3.re - X[12].im * m3.im);
  det.im = (X[0].re * m0.im + X[0].im * m0.re) +
           (X[4].re * m1.im + X[4].im * m1.re) +
           (X[8].re * m2.im + X[8].im * m2.re) +
           (X[12].re * m3.im + X[12].im * m3.re);
  return det;
}
DTYPE det_5by5(DTYPE* X) {
  DTYPE det;
  DTYPE t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16,
      t17, t18, t19, t20;
  DTYPE m0, m1, m2, m3, m4;
#if DEBUG
  printf("%s\n", __func__);
#endif
  t1 = X[18] * X[24] - X[23] * X[19];
  t2 = X[13] * X[24] - X[23] * X[14];
  t3 = X[13] * X[19] - X[18] * X[14];
  t4 = X[8] * X[24] - X[23] * X[9];
  t5 = X[8] * X[19] - X[18] * X[9];
  t6 = X[8] * X[14] - X[13] * X[9];
  t7 = X[3] * X[24] - X[23] * X[4];
  t8 = X[3] * X[19] - X[18] * X[4];
  t9 = X[3] * X[14] - X[13] * X[4];
  t10 = X[3] * X[9] - X[8] * X[4];

  t11 = X[12] * t1 - X[17] * t2 + X[22] * t3;
  t12 = X[7] * t1 - X[17] * t4 + X[22] * t5;
  t13 = X[7] * t2 - X[12] * t4 + X[22] * t6;
  t14 = X[7] * t3 - X[12] * t5 + X[17] * t6;
  t15 = X[2] * t1 - X[17] * t7 + X[22] * t8;
  t16 = X[2] * t2 - X[12] * t7 + X[22] * t9;
  t17 = X[2] * t3 - X[12] * t8 + X[17] * t9;

  m0 = X[6] * t11 - X[11] * t12 + X[16] * t13 - X[21] * t14;
  X[15] * t13 + X[20] * t14;
  m1 = -X[1] * t11 + X[11] * t15 - X[16] * t16 + X[21] * t17;
  X[15] * t16 - X[20] * t17;

  t18 = X[2] * t4 - X[7] * t7 + X[22] * t10;
  t19 = X[2] * t5 - X[7] * t8 + X[17] * t10;
  t20 = X[2] * t6 - X[7] * t9 + X[12] * t10;

  m2 = X[1] * t12 - X[6] * t15 + X[16] * t18 - X[21] * t19;
  X[15] * t18 + X[20] * t19;
  m3 = -X[1] * t13 + X[6] * t16 - X[11] * t18 + X[21] * t20;
  X[20] * t20;
  m4 = X[1] * t14 - X[6] * t17 + X[11] * t19 - X[16] * t20;
  X[10] * t19 + X[15] * t20;

  det = X[0] * m0 + X[5] * m1 + X[10] * m2 + X[15] * m3 + X[20] * m4;

  return det;
}

CTYPE cdet_5by5(CTYPE* X) {
  CTYPE det;
  CTYPE t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16,
      t17, t18, t19, t20;
  CTYPE m0, m1, m2, m3, m4;
#if DEBUG
  printf("%s\n", __func__);
#endif

  t1.re = (X[18].re * X[24].re - X[18].im * X[24].im) -
          (X[23].re * X[19].re - X[23].im * X[19].im);
  t1.im = (X[18].re * X[24].im + X[18].im * X[24].re) -
          (X[23].re * X[19].im + X[23].im * X[19].re);

  t2.re = (X[13].re * X[24].re - X[13].im * X[24].im) -
          (X[23].re * X[14].re - X[23].im * X[14].im);
  t2.im = (X[13].re * X[24].im + X[13].im * X[24].re) -
          (X[23].re * X[14].im + X[23].im * X[14].re);

  t3.re = (X[13].re * X[19].re - X[13].im * X[19].im) -
          (X[18].re * X[14].re - X[18].im * X[14].im);
  t3.im = (X[13].re * X[19].im + X[13].im * X[19].re) -
          (X[18].re * X[14].im + X[18].im * X[14].re);

  t4.re = (X[8].re * X[24].re - X[8].im * X[24].im) -
          (X[23].re * X[9].re - X[23].im * X[9].im);
  t4.im = (X[8].re * X[24].im + X[8].im * X[24].re) -
          (X[23].re * X[9].im + X[23].im * X[9].re);

  t5.re = (X[8].re * X[19].re - X[8].im * X[19].im) -
          (X[18].re * X[9].re - X[18].im * X[9].im);
  t5.im = (X[8].re * X[19].im + X[8].im * X[19].re) -
          (X[18].re * X[9].im + X[18].im * X[9].re);

  t6.re = (X[8].re * X[14].re - X[8].im * X[14].im) -
          (X[13].re * X[9].re - X[13].im * X[9].im);
  t6.im = (X[8].re * X[14].im + X[8].im * X[14].re) -
          (X[13].re * X[9].im + X[13].im * X[9].re);

  t7.re = (X[3].re * X[24].re - X[3].im * X[24].im) -
          (X[23].re * X[4].re - X[23].im * X[4].im);
  t7.im = (X[3].re * X[24].im + X[3].im * X[24].re) -
          (X[23].re * X[4].im + X[23].im * X[4].re);

  t8.re = (X[3].re * X[19].re - X[3].im * X[19].im) -
          (X[18].re * X[4].re - X[18].im * X[4].im);
  t8.im = (X[3].re * X[19].im + X[3].im * X[19].re) -
          (X[18].re * X[4].im + X[18].im * X[4].re);

  t9.re = (X[3].re * X[14].re - X[3].im * X[14].im) -
          (X[13].re * X[4].re - X[13].im * X[4].im);
  t9.im = (X[3].re * X[14].im + X[3].im * X[14].re) -
          (X[13].re * X[4].im + X[13].im * X[4].re);

  t10.re = (X[3].re * X[9].re - X[3].im * X[9].im) -
           (X[8].re * X[4].re - X[8].im * X[4].im);
  t10.im = (X[3].re * X[9].im + X[3].im * X[9].re) -
           (X[8].re * X[4].im + X[8].im * X[4].re);

  t11.re = (X[12].re * t1.re - X[12].im * t1.im) -
           (X[17].re * t2.re - X[17].im * t2.im) +
           (X[22].re * t3.re - X[22].im * t3.im);
  t11.im = (X[12].re * t1.im + X[12].im * t1.re) -
           (X[17].re * t2.im + X[17].im * t2.re) +
           (X[22].re * t3.im + X[22].im * t3.re);

  t12.re = (X[7].re * t1.re - X[7].im * t1.im) -
           (X[17].re * t4.re - X[17].im * t4.im) +
           (X[22].re * t5.re - X[22].im * t5.im);
  t12.im = (X[7].re * t1.im + X[7].im * t1.re) -
           (X[17].re * t4.im + X[17].im * t4.re) +
           (X[22].re * t5.im + X[22].im * t5.re);

  t13.re = (X[7].re * t2.re - X[7].im * t2.im) -
           (X[12].re * t4.re - X[12].im * t4.im) +
           (X[22].re * t6.re - X[22].im * t6.im);
  t13.im = (X[7].re * t2.im + X[7].im * t2.re) -
           (X[12].re * t4.im + X[12].im * t4.re) +
           (X[22].re * t6.im + X[22].im * t6.re);

  t14.re = (X[7].re * t3.re - X[7].im * t3.im) -
           (X[12].re * t5.re - X[12].im * t5.im) +
           (X[17].re * t6.re - X[17].im * t6.im);
  t14.im = (X[7].re * t3.im + X[7].im * t3.re) -
           (X[12].re * t5.im + X[12].im * t5.re) +
           (X[17].re * t6.im + X[17].im * t6.re);

  t15.re = (X[2].re * t1.re - X[2].im * t1.im) -
           (X[17].re * t7.re - X[17].im * t7.im) +
           (X[22].re * t8.re - X[22].im * t8.im);
  t15.im = (X[2].re * t1.im + X[2].im * t1.re) -
           (X[17].re * t7.im + X[17].im * t7.re) +
           (X[22].re * t8.im + X[22].im * t8.re);

  t16.re = (X[2].re * t2.re - X[2].im * t2.im) -
           (X[12].re * t7.re - X[12].im * t7.im) +
           (X[22].re * t9.re - X[22].im * t9.im);
  t16.im = (X[2].re * t2.im + X[2].im * t2.re) -
           (X[12].re * t7.im + X[12].im * t7.re) +
           (X[22].re * t9.im + X[22].im * t9.re);

  t17.re = (X[2].re * t3.re - X[2].im * t3.im) -
           (X[12].re * t8.re - X[12].im * t8.im) +
           (X[17].re * t9.re - X[17].im * t9.im);
  t17.im = (X[2].re * t3.im + X[2].im * t3.re) -
           (X[12].re * t8.im + X[12].im * t8.re) +
           (X[17].re * t9.im + X[17].im * t9.re);

  m0.re = (X[6].re * t11.re - X[6].im * t11.im) -
          (X[11].re * t12.re - X[11].im * t12.im) +
          (X[16].re * t13.re - X[16].im * t13.im) -
          (X[21].re * t14.re - X[21].im * t14.im);
  m0.im = (X[6].re * t11.im + X[6].im * t11.re) -
          (X[11].re * t12.im + X[11].im * t12.re) +
          (X[16].re * t13.im + X[16].im * t13.re) -
          (X[21].re * t14.im + X[21].im * t14.re);

  m1.re = -(X[1].re * t11.re - X[1].im * t11.im) +
          (X[11].re * t15.re - X[11].im * t15.im) -
          (X[16].re * t16.re - X[16].im * t16.im) +
          (X[21].re * t17.re - X[21].im * t17.im);
  m1.im = -(X[1].re * t11.im + X[1].im * t11.re) +
          (X[11].re * t15.im + X[11].im * t15.re) -
          (X[16].re * t16.im + X[16].im * t16.re) +
          (X[21].re * t17.im + X[21].im * t17.re);

  t18.re = (X[2].re * t4.re - X[2].im * t4.im) -
           (X[7].re * t7.re - X[7].im * t7.im) +
           (X[22].re * t10.re - X[22].im * t10.im);
  t18.im = (X[2].re * t4.im + X[2].im * t4.re) -
           (X[7].re * t7.im + X[7].im * t7.re) +
           (X[22].re * t10.im + X[22].im * t10.re);

  t19.re = (X[2].re * t5.re - X[2].im * t5.im) -
           (X[7].re * t8.re - X[7].im * t8.im) +
           (X[17].re * t10.re - X[17].im * t10.im);
  t19.im = (X[2].re * t5.im + X[2].im * t5.re) -
           (X[7].re * t8.im + X[7].im * t8.re) +
           (X[17].re * t10.im + X[17].im * t10.re);

  t20.re = (X[2].re * t6.re - X[2].im * t6.im) -
           (X[7].re * t9.re - X[7].im * t9.im) +
           (X[12].re * t10.re - X[12].im * t10.im);
  t20.im = (X[2].re * t6.im + X[2].im * t6.re) -
           (X[7].re * t9.im + X[7].im * t9.re) +
           (X[12].re * t10.im + X[12].im * t10.re);

  m2.re = (X[1].re * t12.re - X[1].im * t12.im) -
          (X[6].re * t15.re - X[6].im * t15.im) +
          (X[16].re * t18.re - X[16].im * t18.im) -
          (X[21].re * t19.re - X[21].im * t19.im);
  m2.im = (X[1].re * t12.im + X[1].im * t12.re) -
          (X[6].re * t15.im + X[6].im * t15.re) +
          (X[16].re * t18.im + X[16].im * t18.re) -
          (X[21].re * t19.im + X[21].im * t19.re);

  m3.re = -(X[1].re * t13.re - X[1].im * t13.im) +
          (X[6].re * t16.re - X[6].im * t16.im) -
          (X[11].re * t18.re - X[11].im * t18.im) +
          (X[21].re * t20.re - X[21].im * t20.im);
  m3.im = -(X[1].re * t13.im + X[1].im * t13.re) +
          (X[6].re * t16.im + X[6].im * t16.re) -
          (X[11].re * t18.im + X[11].im * t18.re) +
          (X[21].re * t20.im + X[21].im * t20.re);

  m4.re = (X[1].re * t14.re - X[1].im * t14.im) -
          (X[6].re * t17.re - X[6].im * t17.im) +
          (X[11].re * t19.re - X[11].im * t19.im) -
          (X[16].re * t20.re - X[16].im * t20.im);
  m4.im = (X[1].re * t14.im + X[1].im * t14.re) -
          (X[6].re * t17.im + X[6].im * t17.re) +
          (X[11].re * t19.im + X[11].im * t19.re) -
          (X[16].re * t20.im + X[16].im * t20.re);

  det.re = (X[0].re * m0.re - X[0].im * m0.im) +
           (X[5].re * m1.re - X[5].im * m1.im) +
           (X[10].re * m2.re - X[10].im * m2.im) +
           (X[15].re * m3.re - X[15].im * m3.im) +
           (X[20].re * m4.re - X[20].im * m4.im);
  det.im = (X[0].re * m0.im + X[0].im * m0.re) +
           (X[5].re * m1.im + X[5].im * m1.re) +
           (X[10].re * m2.im + X[10].im * m2.re) +
           (X[15].re * m3.im + X[15].im * m3.re) +
           (X[20].re * m4.im + X[20].im * m4.re);

  return det;
}

DTYPE det_6by6(DTYPE* X) {
  DTYPE det;
  DTYPE t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16,
      t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t29, t30, t31,
      t32, t33, t34, t35, t36, t37, t38, t39, t40, t41, t42, t43, t44, t45, t46,
      t47, t48, t49, t50;
  DTYPE m0, m1, m2, m3, m4, m5;
#if DEBUG
  printf("%s\n", __func__);
#endif

  t1 = X[28] * X[35] - X[34] * X[29];
  t2 = X[22] * X[35] - X[34] * X[23];
  t3 = X[22] * X[29] - X[28] * X[23];
  t4 = X[16] * X[35] - X[34] * X[17];
  t5 = X[16] * X[29] - X[28] * X[17];
  t6 = X[16] * X[23] - X[22] * X[17];
  t7 = X[10] * X[35] - X[34] * X[11];
  t8 = X[10] * X[29] - X[28] * X[11];
  t9 = X[10] * X[23] - X[22] * X[11];
  t10 = X[10] * X[17] - X[16] * X[11];
  t11 = X[4] * X[35] - X[34] * X[5];
  t12 = X[4] * X[29] - X[28] * X[5];
  t13 = X[4] * X[23] - X[22] * X[5];
  t14 = X[4] * X[17] - X[16] * X[5];
  t15 = X[4] * X[11] - X[10] * X[5];

  t16 = X[21] * t1 - X[27] * t2 + X[33] * t3;
  t17 = X[15] * t1 - X[27] * t4 + X[33] * t5;
  t18 = X[15] * t2 - X[21] * t4 + X[33] * t6;
  t19 = X[15] * t3 - X[21] * t5 + X[27] * t6;
  t20 = X[9] * t1 - X[27] * t7 + X[33] * t8;
  t21 = X[9] * t2 - X[21] * t7 + X[33] * t9;
  t22 = X[9] * t3 - X[21] * t8 + X[27] * t9;
  t23 = X[9] * t4 - X[15] * t7 + X[33] * t10;
  t24 = X[9] * t5 - X[15] * t8 + X[27] * t10;
  t25 = X[9] * t6 - X[15] * t9 + X[21] * t10;
  t26 = X[3] * t1 - X[27] * t11 + X[33] * t12;
  t27 = X[3] * t2 - X[21] * t11 + X[33] * t13;
  t28 = X[3] * t3 - X[21] * t12 + X[27] * t13;
  t29 = X[3] * t4 - X[15] * t11 + X[33] * t14;
  t30 = X[3] * t5 - X[15] * t12 + X[27] * t14;
  t31 = X[3] * t6 - X[15] * t13 + X[21] * t14;
  t32 = X[3] * t7 - X[9] * t11 + X[33] * t15;
  t33 = X[3] * t8 - X[9] * t12 + X[27] * t15;
  t34 = X[3] * t9 - X[9] * t13 + X[21] * t15;
  t35 = X[3] * t10 - X[9] * t14 + X[15] * t15;

  t36 = X[14] * t16 - X[20] * t17 + X[26] * t18 - X[32] * t19;
  t37 = X[8] * t16 - X[20] * t20 + X[26] * t21 - X[32] * t22;
  t38 = X[8] * t17 - X[14] * t20 + X[26] * t23 - X[32] * t24;
  t39 = X[8] * t18 - X[14] * t21 + X[20] * t23 - X[32] * t25;
  t40 = X[8] * t19 - X[14] * t22 + X[20] * t24 - X[26] * t25;
  t41 = X[2] * t16 - X[20] * t26 + X[26] * t27 - X[32] * t28;
  t42 = X[2] * t17 - X[14] * t26 + X[26] * t29 - X[32] * t30;
  t43 = X[2] * t18 - X[14] * t27 + X[20] * t29 - X[32] * t31;
  t44 = X[2] * t19 - X[14] * t28 + X[20] * t30 - X[26] * t31;

  m0 = X[7] * t36 - X[13] * t37 + X[19] * t38 - X[25] * t39 + X[31] * t40;
  m1 = -X[1] * t36 + X[13] * t41 - X[19] * t42 + X[25] * t43 - X[31] * t44;

  t45 = X[2] * t20 - X[8] * t26 + X[26] * t32 - X[32] * t33;
  t46 = X[2] * t21 - X[8] * t27 + X[20] * t32 - X[32] * t34;
  t47 = X[2] * t22 - X[8] * t28 + X[20] * t33 - X[26] * t34;
  t48 = X[2] * t23 - X[8] * t29 + X[14] * t32 - X[32] * t35;
  t49 = X[2] * t24 - X[8] * t30 + X[14] * t33 - X[26] * t35;

  m2 = X[1] * t37 - X[7] * t41 + X[19] * t45 - X[25] * t46 + X[31] * t47;
  m3 = -X[1] * t38 + X[7] * t42 - X[13] * t45 + X[25] * t48 - X[31] * t49;

  t50 = X[2] * t25 - X[8] * t31 + X[14] * t34 - X[20] * t35;

  m4 = X[1] * t39 - X[7] * t43 + X[13] * t46 - X[19] * t48 + X[31] * t50;
  m5 = -X[1] * t40 + X[7] * t44 - X[13] * t47 + X[19] * t49 - X[25] * t50;

  det =
      X[0] * m0 + X[6] * m1 + X[12] * m2 + X[18] * m3 + X[24] * m4 + X[30] * m5;

  return det;
}

CTYPE cdet_6by6(CTYPE* X) {
  CTYPE det;
  CTYPE t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16,
      t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t29, t30, t31,
      t32, t33, t34, t35, t36, t37, t38, t39, t40, t41, t42, t43, t44, t45, t46,
      t47, t48, t49, t50;
  CTYPE m0, m1, m2, m3, m4, m5;
#if DEBUG
  printf("%s\n", __func__);
#endif

  t1.re = (X[28].re * X[35].re - X[28].im * X[35].im) -
          (X[34].re * X[29].re - X[34].im * X[29].im);
  t1.im = (X[28].re * X[35].im + X[28].im * X[35].re) -
          (X[34].re * X[29].im + X[34].im * X[29].re);

  t2.re = (X[22].re * X[35].re - X[22].im * X[35].im) -
          (X[34].re * X[23].re - X[34].im * X[23].im);
  t2.im = (X[22].re * X[35].im + X[22].im * X[35].re) -
          (X[34].re * X[23].im + X[34].im * X[23].re);

  t3.re = (X[22].re * X[29].re - X[22].im * X[29].im) -
          (X[28].re * X[23].re - X[28].im * X[23].im);
  t3.im = (X[22].re * X[29].im + X[22].im * X[29].re) -
          (X[28].re * X[23].im + X[28].im * X[23].re);

  t4.re = (X[16].re * X[35].re - X[16].im * X[35].im) -
          (X[34].re * X[17].re - X[34].im * X[17].im);
  t4.im = (X[16].re * X[35].im + X[16].im * X[35].re) -
          (X[34].re * X[17].im + X[34].im * X[17].re);

  t5.re = (X[16].re * X[29].re - X[16].im * X[29].im) -
          (X[28].re * X[17].re - X[28].im * X[17].im);
  t5.im = (X[16].re * X[29].im + X[16].im * X[29].re) -
          (X[28].re * X[17].im + X[28].im * X[17].re);

  t6.re = (X[16].re * X[23].re - X[16].im * X[23].im) -
          (X[22].re * X[17].re - X[22].im * X[17].im);
  t6.im = (X[16].re * X[23].im + X[16].im * X[23].re) -
          (X[22].re * X[17].im + X[22].im * X[17].re);

  t7.re = (X[10].re * X[35].re - X[10].im * X[35].im) -
          (X[34].re * X[11].re - X[34].im * X[11].im);
  t7.im = (X[10].re * X[35].im + X[10].im * X[35].re) -
          (X[34].re * X[11].im + X[34].im * X[11].re);

  t8.re = (X[10].re * X[29].re - X[10].im * X[29].im) -
          (X[28].re * X[11].re - X[28].im * X[11].im);
  t8.im = (X[10].re * X[29].im + X[10].im * X[29].re) -
          (X[28].re * X[11].im + X[28].im * X[11].re);

  t9.re = (X[10].re * X[23].re - X[10].im * X[23].im) -
          (X[22].re * X[11].re - X[22].im * X[11].im);
  t9.im = (X[10].re * X[23].im + X[10].im * X[23].re) -
          (X[22].re * X[11].im + X[22].im * X[11].re);

  t10.re = (X[10].re * X[17].re - X[10].im * X[17].im) -
           (X[16].re * X[11].re - X[16].im * X[11].im);
  t10.im = (X[10].re * X[17].im + X[10].im * X[17].re) -
           (X[16].re * X[11].im + X[16].im * X[11].re);

  t11.re = (X[4].re * X[35].re - X[4].im * X[35].im) -
           (X[34].re * X[5].re - X[34].im * X[5].im);
  t11.im = (X[4].re * X[35].im + X[4].im * X[35].re) -
           (X[34].re * X[5].im + X[34].im * X[5].re);

  t12.re = (X[4].re * X[29].re - X[4].im * X[29].im) -
           (X[28].re * X[5].re - X[28].im * X[5].im);
  t12.im = (X[4].re * X[29].im + X[4].im * X[29].re) -
           (X[28].re * X[5].im + X[28].im * X[5].re);

  t13.re = (X[4].re * X[23].re - X[4].im * X[23].im) -
           (X[22].re * X[5].re - X[22].im * X[5].im);
  t13.im = (X[4].re * X[23].im + X[4].im * X[23].re) -
           (X[22].re * X[5].im + X[22].im * X[5].re);

  t14.re = (X[4].re * X[17].re - X[4].im * X[17].im) -
           (X[16].re * X[5].re - X[16].im * X[5].im);
  t14.im = (X[4].re * X[17].im + X[4].im * X[17].re) -
           (X[16].re * X[5].im + X[16].im * X[5].re);

  t15.re = (X[4].re * X[11].re - X[4].im * X[11].im) -
           (X[10].re * X[5].re - X[10].im * X[5].im);
  t15.im = (X[4].re * X[11].im + X[4].im * X[11].re) -
           (X[10].re * X[5].im + X[10].im * X[5].re);

  t16.re = (X[21].re * t1.re - X[21].im * t1.im) -
           (X[27].re * t2.re - X[27].im * t2.im) +
           (X[33].re * t3.re - X[33].im * t3.im);
  t16.im = (X[21].re * t1.im + X[21].im * t1.re) -
           (X[27].re * t2.im + X[27].im * t2.re) +
           (X[33].re * t3.im + X[33].im * t3.re);

  t17.re = (X[15].re * t1.re - X[15].im * t1.im) -
           (X[27].re * t4.re - X[27].im * t4.im) +
           (X[33].re * t5.re - X[33].im * t5.im);
  t17.im = (X[15].re * t1.im + X[15].im * t1.re) -
           (X[27].re * t4.im + X[27].im * t4.re) +
           (X[33].re * t5.im + X[33].im * t5.re);

  t18.re = (X[15].re * t2.re - X[15].im * t2.im) -
           (X[21].re * t4.re - X[21].im * t4.im) +
           (X[33].re * t6.re - X[33].im * t6.im);
  t18.im = (X[15].re * t2.im + X[15].im * t2.re) -
           (X[21].re * t4.im + X[21].im * t4.re) +
           (X[33].re * t6.im + X[33].im * t6.re);

  t19.re = (X[15].re * t3.re - X[15].im * t3.im) -
           (X[21].re * t5.re - X[21].im * t5.im) +
           (X[27].re * t6.re - X[27].im * t6.im);
  t19.im = (X[15].re * t3.im + X[15].im * t3.re) -
           (X[21].re * t5.im + X[21].im * t5.re) +
           (X[27].re * t6.im + X[27].im * t6.re);

  t20.re = (X[9].re * t1.re - X[9].im * t1.im) -
           (X[27].re * t7.re - X[27].im * t7.im) +
           (X[33].re * t8.re - X[33].im * t8.im);
  t20.im = (X[9].re * t1.im + X[9].im * t1.re) -
           (X[27].re * t7.im + X[27].im * t7.re) +
           (X[33].re * t8.im + X[33].im * t8.re);

  t21.re = (X[9].re * t2.re - X[9].im * t2.im) -
           (X[21].re * t7.re - X[21].im * t7.im) +
           (X[33].re * t9.re - X[33].im * t9.im);
  t21.im = (X[9].re * t2.im + X[9].im * t2.re) -
           (X[21].re * t7.im + X[21].im * t7.re) +
           (X[33].re * t9.im + X[33].im * t9.re);

  t22.re = (X[9].re * t3.re - X[9].im * t3.im) -
           (X[21].re * t8.re - X[21].im * t8.im) +
           (X[27].re * t9.re - X[27].im * t9.im);
  t22.im = (X[9].re * t3.im + X[9].im * t3.re) -
           (X[21].re * t8.im + X[21].im * t8.re) +
           (X[27].re * t9.im + X[27].im * t9.re);

  t23.re = (X[9].re * t4.re - X[9].im * t4.im) -
           (X[15].re * t7.re - X[15].im * t7.im) +
           (X[33].re * t10.re - X[33].im * t10.im);
  t23.im = (X[9].re * t4.im + X[9].im * t4.re) -
           (X[15].re * t7.im + X[15].im * t7.re) +
           (X[33].re * t10.im + X[33].im * t10.re);

  t24.re = (X[9].re * t5.re - X[9].im * t5.im) -
           (X[15].re * t8.re - X[15].im * t8.im) +
           (X[27].re * t10.re - X[27].im * t10.im);
  t24.im = (X[9].re * t5.im + X[9].im * t5.re) -
           (X[15].re * t8.im + X[15].im * t8.re) +
           (X[27].re * t10.im + X[27].im * t10.re);

  t25.re = (X[9].re * t6.re - X[9].im * t6.im) -
           (X[15].re * t9.re - X[15].im * t9.im) +
           (X[21].re * t10.re - X[21].im * t10.im);
  t25.im = (X[9].re * t6.im + X[9].im * t6.re) -
           (X[15].re * t9.im + X[15].im * t9.re) +
           (X[21].re * t10.im + X[21].im * t10.re);

  t26.re = (X[3].re * t1.re - X[3].im * t1.im) -
           (X[27].re * t11.re - X[27].im * t11.im) +
           (X[33].re * t12.re - X[33].im * t12.im);
  t26.im = (X[3].re * t1.im + X[3].im * t1.re) -
           (X[27].re * t11.im + X[27].im * t11.re) +
           (X[33].re * t12.im + X[33].im * t12.re);

  t27.re = (X[3].re * t2.re - X[3].im * t2.im) -
           (X[21].re * t11.re - X[21].im * t11.im) +
           (X[33].re * t13.re - X[33].im * t13.im);
  t27.im = (X[3].re * t2.im + X[3].im * t2.re) -
           (X[21].re * t11.im + X[21].im * t11.re) +
           (X[33].re * t13.im + X[33].im * t13.re);

  t28.re = (X[3].re * t3.re - X[3].im * t3.im) -
           (X[21].re * t12.re - X[21].im * t12.im) +
           (X[27].re * t13.re - X[27].im * t13.im);
  t28.im = (X[3].re * t3.im + X[3].im * t3.re) -
           (X[21].re * t12.im + X[21].im * t12.re) +
           (X[27].re * t13.im + X[27].im * t13.re);

  t29.re = (X[3].re * t4.re - X[3].im * t4.im) -
           (X[15].re * t11.re - X[15].im * t11.im) +
           (X[33].re * t14.re - X[33].im * t14.im);
  t29.im = (X[3].re * t4.im + X[3].im * t4.re) -
           (X[15].re * t11.im + X[15].im * t11.re) +
           (X[33].re * t14.im + X[33].im * t14.re);

  t30.re = (X[3].re * t5.re - X[3].im * t5.im) -
           (X[15].re * t12.re - X[15].im * t12.im) +
           (X[27].re * t14.re - X[27].im * t14.im);
  t30.im = (X[3].re * t5.im + X[3].im * t5.re) -
           (X[15].re * t12.im + X[15].im * t12.re) +
           (X[27].re * t14.im + X[27].im * t14.re);

  t31.re = (X[3].re * t6.re - X[3].im * t6.im) -
           (X[15].re * t13.re - X[15].im * t13.im) +
           (X[21].re * t14.re - X[21].im * t14.im);
  t31.im = (X[3].re * t6.im + X[3].im * t6.re) -
           (X[15].re * t13.im + X[15].im * t13.re) +
           (X[21].re * t14.im + X[21].im * t14.re);

  t32.re = (X[3].re * t7.re - X[3].im * t7.im) -
           (X[9].re * t11.re - X[9].im * t11.im) +
           (X[33].re * t15.re - X[33].im * t15.im);
  t32.im = (X[3].re * t7.im + X[3].im * t7.re) -
           (X[9].re * t11.im + X[9].im * t11.re) +
           (X[33].re * t15.im + X[33].im * t15.re);

  t33.re = (X[3].re * t8.re - X[3].im * t8.im) -
           (X[9].re * t12.re - X[9].im * t12.im) +
           (X[27].re * t15.re - X[27].im * t15.im);
  t33.im = (X[3].re * t8.im + X[3].im * t8.re) -
           (X[9].re * t12.im + X[9].im * t12.re) +
           (X[27].re * t15.im + X[27].im * t15.re);

  t34.re = (X[3].re * t9.re - X[3].im * t9.im) -
           (X[9].re * t13.re - X[9].im * t13.im) +
           (X[21].re * t15.re - X[21].im * t15.im);
  t34.im = (X[3].re * t9.im + X[3].im * t9.re) -
           (X[9].re * t13.im + X[9].im * t13.re) +
           (X[21].re * t15.im + X[21].im * t15.re);

  t35.re = (X[3].re * t10.re - X[3].im * t10.im) -
           (X[9].re * t14.re - X[9].im * t14.im) +
           (X[15].re * t15.re - X[15].im * t15.im);
  t35.im = (X[3].re * t10.im + X[3].im * t10.re) -
           (X[9].re * t14.im + X[9].im * t14.re) +
           (X[15].re * t15.im + X[15].im * t15.re);

  t36.re = (X[14].re * t16.re - X[14].im * t16.im) -
           (X[20].re * t17.re - X[20].im * t17.im) +
           (X[26].re * t18.re - X[26].im * t18.im) -
           (X[32].re * t19.re - X[32].im * t19.im);
  t36.im = (X[14].re * t16.im + X[14].im * t16.re) -
           (X[20].re * t17.im + X[20].im * t17.re) +
           (X[26].re * t18.im + X[26].im * t18.re) -
           (X[32].re * t19.im + X[32].im * t19.re);

  t37.re = (X[8].re * t16.re - X[8].im * t16.im) -
           (X[20].re * t20.re - X[20].im * t20.im) +
           (X[26].re * t21.re - X[26].im * t21.im) -
           (X[32].re * t22.re - X[32].im * t22.im);
  t37.im = (X[8].re * t16.im + X[8].im * t16.re) -
           (X[20].re * t20.im + X[20].im * t20.re) +
           (X[26].re * t21.im + X[26].im * t21.re) -
           (X[32].re * t22.im + X[32].im * t22.re);

  t38.re = (X[8].re * t17.re - X[8].im * t17.im) -
           (X[14].re * t20.re - X[14].im * t20.im) +
           (X[26].re * t23.re - X[26].im * t23.im) -
           (X[32].re * t24.re - X[32].im * t24.im);
  t38.im = (X[8].re * t17.im + X[8].im * t17.re) -
           (X[14].re * t20.im + X[14].im * t20.re) +
           (X[26].re * t23.im + X[26].im * t23.re) -
           (X[32].re * t24.im + X[32].im * t24.re);

  t39.re = (X[8].re * t18.re - X[8].im * t18.im) -
           (X[14].re * t21.re - X[14].im * t21.im) +
           (X[20].re * t23.re - X[20].im * t23.im) -
           (X[32].re * t25.re - X[32].im * t25.im);
  t39.im = (X[8].re * t18.im + X[8].im * t18.re) -
           (X[14].re * t21.im + X[14].im * t21.re) +
           (X[20].re * t23.im + X[20].im * t23.re) -
           (X[32].re * t25.im + X[32].im * t25.re);

  t40.re = (X[8].re * t19.re - X[8].im * t19.im) -
           (X[14].re * t22.re - X[14].im * t22.im) +
           (X[20].re * t24.re - X[20].im * t24.im) -
           (X[26].re * t25.re - X[26].im * t25.im);
  t40.im = (X[8].re * t19.im + X[8].im * t19.re) -
           (X[14].re * t22.im + X[14].im * t22.re) +
           (X[20].re * t24.im + X[20].im * t24.re) -
           (X[26].re * t25.im + X[26].im * t25.re);

  t41.re = (X[2].re * t16.re - X[2].im * t16.im) -
           (X[20].re * t26.re - X[20].im * t26.im) +
           (X[26].re * t27.re - X[26].im * t27.im) -
           (X[32].re * t28.re - X[32].im * t28.im);
  t41.im = (X[2].re * t16.im + X[2].im * t16.re) -
           (X[20].re * t26.im + X[20].im * t26.re) +
           (X[26].re * t27.im + X[26].im * t27.re) -
           (X[32].re * t28.im + X[32].im * t28.re);

  t42.re = (X[2].re * t17.re - X[2].im * t17.im) -
           (X[14].re * t26.re - X[14].im * t26.im) +
           (X[26].re * t29.re - X[26].im * t29.im) -
           (X[32].re * t30.re - X[32].im * t30.im);
  t42.im = (X[2].re * t17.im + X[2].im * t17.re) -
           (X[14].re * t26.im + X[14].im * t26.re) +
           (X[26].re * t29.im + X[26].im * t29.re) -
           (X[32].re * t30.im + X[32].im * t30.re);

  t43.re = (X[2].re * t18.re - X[2].im * t18.im) -
           (X[14].re * t27.re - X[14].im * t27.im) +
           (X[20].re * t29.re - X[20].im * t29.im) -
           (X[32].re * t31.re - X[32].im * t31.im);
  t43.im = (X[2].re * t18.im + X[2].im * t18.re) -
           (X[14].re * t27.im + X[14].im * t27.re) +
           (X[20].re * t29.im + X[20].im * t29.re) -
           (X[32].re * t31.im + X[32].im * t31.re);

  t44.re = (X[2].re * t19.re - X[2].im * t19.im) -
           (X[14].re * t28.re - X[14].im * t28.im) +
           (X[20].re * t30.re - X[20].im * t30.im) -
           (X[26].re * t31.re - X[26].im * t31.im);
  t44.im = (X[2].re * t19.im + X[2].im * t19.re) -
           (X[14].re * t28.im + X[14].im * t28.re) +
           (X[20].re * t30.im + X[20].im * t30.re) -
           (X[26].re * t31.im + X[26].im * t31.re);

  m0.re = (X[7].re * t36.re - X[7].im * t36.im) -
          (X[13].re * t37.re - X[13].im * t37.im) +
          (X[19].re * t38.re - X[19].im * t38.im) -
          (X[25].re * t39.re - X[25].im * t39.im) +
          (X[31].re * t40.re - X[31].im * t40.im);
  m0.im = (X[7].re * t36.im + X[7].im * t36.re) -
          (X[13].re * t37.im + X[13].im * t37.re) +
          (X[19].re * t38.im + X[19].im * t38.re) -
          (X[25].re * t39.im + X[25].im * t39.re) +
          (X[31].re * t40.im + X[31].im * t40.re);

  m1.re = -(X[1].re * t36.re - X[1].im * t36.im) +
          (X[13].re * t41.re - X[13].im * t41.im) -
          (X[19].re * t42.re - X[19].im * t42.im) +
          (X[25].re * t43.re - X[25].im * t43.im) -
          (X[31].re * t44.re - X[31].im * t44.im);
  m1.im = -(X[1].re * t36.im + X[1].im * t36.re) +
          (X[13].re * t41.im + X[13].im * t41.re) -
          (X[19].re * t42.im + X[19].im * t42.re) +
          (X[25].re * t43.im + X[25].im * t43.re) -
          (X[31].re * t44.im + X[31].im * t44.re);

  t45.re = (X[2].re * t20.re - X[2].im * t20.im) -
           (X[8].re * t26.re - X[8].im * t26.im) +
           (X[26].re * t32.re - X[26].im * t32.im) -
           (X[32].re * t33.re - X[32].im * t33.im);
  t45.im = (X[2].re * t20.im + X[2].im * t20.re) -
           (X[8].re * t26.im + X[8].im * t26.re) +
           (X[26].re * t32.im + X[26].im * t32.re) -
           (X[32].re * t33.im + X[32].im * t33.re);

  t46.re = (X[2].re * t21.re - X[2].im * t21.im) -
           (X[8].re * t27.re - X[8].im * t27.im) +
           (X[20].re * t32.re - X[20].im * t32.im) -
           (X[32].re * t34.re - X[32].im * t34.im);
  t46.im = (X[2].re * t21.im + X[2].im * t21.re) -
           (X[8].re * t27.im + X[8].im * t27.re) +
           (X[20].re * t32.im + X[20].im * t32.re) -
           (X[32].re * t34.im + X[32].im * t34.re);

  t47.re = (X[2].re * t22.re - X[2].im * t22.im) -
           (X[8].re * t28.re - X[8].im * t28.im) +
           (X[20].re * t33.re - X[20].im * t33.im) -
           (X[26].re * t34.re - X[26].im * t34.im);
  t47.im = (X[2].re * t22.im + X[2].im * t22.re) -
           (X[8].re * t28.im + X[8].im * t28.re) +
           (X[20].re * t33.im + X[20].im * t33.re) -
           (X[26].re * t34.im + X[26].im * t34.re);

  t48.re = (X[2].re * t23.re - X[2].im * t23.im) -
           (X[8].re * t29.re - X[8].im * t29.im) +
           (X[14].re * t32.re - X[14].im * t32.im) -
           (X[32].re * t35.re - X[32].im * t35.im);
  t48.im = (X[2].re * t23.im + X[2].im * t23.re) -
           (X[8].re * t29.im + X[8].im * t29.re) +
           (X[14].re * t32.im + X[14].im * t32.re) -
           (X[32].re * t35.im + X[32].im * t35.re);

  t49.re = (X[2].re * t24.re - X[2].im * t24.im) -
           (X[8].re * t30.re - X[8].im * t30.im) +
           (X[14].re * t33.re - X[14].im * t33.im) -
           (X[26].re * t35.re - X[26].im * t35.im);
  t49.im = (X[2].re * t24.im + X[2].im * t24.re) -
           (X[8].re * t30.im + X[8].im * t30.re) +
           (X[14].re * t33.im + X[14].im * t33.re) -
           (X[26].re * t35.im + X[26].im * t35.re);

  m2.re = (X[1].re * t37.re - X[1].im * t37.im) -
          (X[7].re * t41.re - X[7].im * t41.im) +
          (X[19].re * t45.re - X[19].im * t45.im) -
          (X[25].re * t46.re - X[25].im * t46.im) +
          (X[31].re * t47.re - X[31].im * t47.im);
  m2.im = (X[1].re * t37.im + X[1].im * t37.re) -
          (X[7].re * t41.im + X[7].im * t41.re) +
          (X[19].re * t45.im + X[19].im * t45.re) -
          (X[25].re * t46.im + X[25].im * t46.re) +
          (X[31].re * t47.im + X[31].im * t47.re);

  m3.re = -(X[1].re * t38.re - X[1].im * t38.im) +
          (X[7].re * t42.re - X[7].im * t42.im) -
          (X[13].re * t45.re - X[13].im * t45.im) +
          (X[25].re * t48.re - X[25].im * t48.im) -
          (X[31].re * t49.re - X[31].im * t49.im);
  m3.im = -(X[1].re * t38.im + X[1].im * t38.re) +
          (X[7].re * t42.im + X[7].im * t42.re) -
          (X[13].re * t45.im + X[13].im * t45.re) +
          (X[25].re * t48.im + X[25].im * t48.re) -
          (X[31].re * t49.im + X[31].im * t49.re);

  t50.re = (X[2].re * t25.re - X[2].im * t25.im) -
           (X[8].re * t31.re - X[8].im * t31.im) +
           (X[14].re * t34.re - X[14].im * t34.im) -
           (X[20].re * t35.re - X[20].im * t35.im);
  t50.im = (X[2].re * t25.im + X[2].im * t25.re) -
           (X[8].re * t31.im + X[8].im * t31.re) +
           (X[14].re * t34.im + X[14].im * t34.re) -
           (X[20].re * t35.im + X[20].im * t35.re);

  m4.re = (X[1].re * t39.re - X[1].im * t39.im) -
          (X[7].re * t43.re - X[7].im * t43.im) +
          (X[13].re * t46.re - X[13].im * t46.im) -
          (X[19].re * t48.re - X[19].im * t48.im) +
          (X[31].re * t50.re - X[31].im * t50.im);
  m4.im = (X[1].re * t39.im + X[1].im * t39.re) -
          (X[7].re * t43.im + X[7].im * t43.re) +
          (X[13].re * t46.im + X[13].im * t46.re) -
          (X[19].re * t48.im + X[19].im * t48.re) +
          (X[31].re * t50.im + X[31].im * t50.re);

  m5.re = -(X[1].re * t40.re - X[1].im * t40.im) +
          (X[7].re * t44.re - X[7].im * t44.im) -
          (X[13].re * t47.re - X[13].im * t47.im) +
          (X[19].re * t49.re - X[19].im * t49.im) -
          (X[25].re * t50.re - X[25].im * t50.im);
  m5.im = -(X[1].re * t40.im + X[1].im * t40.re) +
          (X[7].re * t44.im + X[7].im * t44.re) -
          (X[13].re * t47.im + X[13].im * t47.re) +
          (X[19].re * t49.im + X[19].im * t49.re) -
          (X[25].re * t50.im + X[25].im * t50.re);

  det.re = (X[0].re * m0.re - X[0].im * m0.im) +
           (X[6].re * m1.re - X[6].im * m1.im) +
           (X[12].re * m2.re - X[12].im * m2.im) +
           (X[18].re * m3.re - X[18].im * m3.im) +
           (X[24].re * m4.re - X[24].im * m4.im) +
           (X[30].re * m5.re - X[30].im * m5.im);
  det.im = (X[0].re * m0.im + X[0].im * m0.re) +
           (X[6].re * m1.im + X[6].im * m1.re) +
           (X[12].re * m2.im + X[12].im * m2.re) +
           (X[18].re * m3.im + X[18].im * m3.re) +
           (X[24].re * m4.im + X[24].im * m4.re) +
           (X[30].re * m5.im + X[30].im * m5.re);
  return det;
}

/**** LAPACK *aq***/
void invert_nbyn(DTYPE* X, DTYPE* Y, UINT n) {
  UINT* idx;
  ITER i;
#if DEBUG
  printf("%s\n", __func__);
#endif
#if USE_CBLAS
  idx = mpalloc(sizeof(UINT) * n);

#pragma omp parallel for schedule(dynamic) shared(X, Y) private(i)
  for (i = 0; i < n * n; i++) Y[i] = X[i];

#if NTYPE == 0
  LAPACKE_sgetrf(LAPACK_COL_MAJOR, n, n, Y, n, idx);
  LAPACKE_sgetri(LAPACK_COL_MAJOR, n, Y, n, idx);
#elif NTYPE == 1
  LAPACKE_dgetrf(LAPACK_COL_MAJOR, n, n, Y, n, idx);
  LAPACKE_dgetri(LAPACK_COL_MAJOR, n, Y, n, idx);
#endif
  mpfree(idx);

#else
  printf("ERROR : 'OpenBLAS' or 'INTEL MKL' is required for this operation\n");
  return;
#endif
}

void cinvert_nbyn(CTYPE* X, CTYPE* Y, UINT n) {
  UINT* idx;
  ITER i;
#if DEBUG
  printf("%s\n", __func__);
#endif
#if USE_CBLAS
  idx = mpalloc(sizeof(UINT) * n);
#pragma omp parallel for schedule(dynamic) shared(X, Y) private(i)
  for (i = 0; i < n; i++) {
    Y[i].re = X[i].re;
    Y[i].im = X[i].im;
  }
#if NTYPE == 0
  LAPACKE_cgetrf(LAPACK_COL_MAJOR, n, n, Y, n, idx);
  LAPACKE_cgetri(LAPACK_COL_MAJOR, n, Y, n, idx);
#elif NTYPE == 1
  LAPACKE_zgetrf(LAPACK_COL_MAJOR, n, n, Y, n, idx);
  LAPACKE_zgetri(LAPACK_COL_MAJOR, n, Y, n, idx);
#endif
  mpfree(idx);

#else
  printf("ERROR : 'OpenBLAS' or 'INTEL MKL' is required for this operation\n");
  return;
#endif
}

DTYPE det_nbyn(MAT* mat) {
  ITER i;
  DTYPE det;
  UINT n;
  UINT* idx;
  MAT* tmat;
#if DEBUG
  printf("%s\n", __func__);
#endif
  n = mat->d0;
#if USE_CBLAS
  idx = mpalloc(sizeof(UINT) * n);
  tmat = mpalloc_mat(mat->d0, mat->d1, mat->d2);
  copy(mat, tmat);
#if NTYPE == 0
  LAPACKE_sgetrf(LAPACK_COL_MAJOR, n, n, mat->data, n, idx);
#elif NTYPE == 1
  LAPACKE_dgetrf(LAPACK_COL_MAJOR, n, n, mat->data, n, idx);
#endif

#else

  printf("ERROR : 'OpenBLAS' or 'INTEL MKL' is required for this operation\n");
  return -1;
#endif
  det = 1;
  for (i = 0; i < mat->d0; i++) {
    if (i + 1 != idx[i])
      det *= -mat->data[i * n + i];
    else
      det *= mat->data[i * n + i];
  }
  mpfree(idx);
  mpfree_mat(tmat);
  return det;
}

CTYPE cdet_nbyn(CMAT* mat) {
  ITER i;
  CTYPE det;
  DTYPE t;
  UINT n;
  UINT* idx;
  CMAT* tmat;
#if DEBUG
  printf("%s\n", __func__);
#endif
  n = mat->d0;
#if USE_CBLAS
  idx = mpalloc(sizeof(UINT) * n);
  tmat = mpalloc_cmat(mat->d0, mat->d1, mat->d2);
  ccopy(mat, tmat);
#if NTYPE == 0
  LAPACKE_cgetrf(LAPACK_COL_MAJOR, n, n, mat->data, n, idx);
#elif NTYPE == 1
  LAPACKE_zgetrf(LAPACK_COL_MAJOR, n, n, mat->data, n, idx);
#endif

#else

  printf("ERROR : 'OpenBLAS' or 'INTEL MKL' is required for this operation\n");
  det.re = -1;
  det.im = -1;
  return det;
#endif
  det.re = 1.0;
  for (i = 0; i < mat->d0; i++) {
    if (i + 1 != idx[i])
      CXMUL(det, -mat->data[i * n + i], t)
    else
      CXMUL(det, mat->data[i * n + i], t)
  }
  mpfree(idx);
  mpfree_cmat(tmat);
  return det;
}

/**** Work In Progress ****/

void wiki_invert(DTYPE* X, UINT n, UINT* idx, DTYPE* Y) {
  ITER i, j, k;
  for (j = 0; j < n; j++) {
    for (i = 0; i < n; i++) {
      if (idx[i] == j)
        Y[i + n * j] = 1.0;
      else
        Y[i + n * j] = 0.0;

      for (k = 0; k < i; k++) Y[i + n * j] -= X[idx[i] + n * k] * Y[k + n * j];
    }

    for (i = n - 1; i >= 0; i--) {
      for (k = i + 1; k < n; k++)
        Y[i + n * j] -= X[idx[i] + n * k] * Y[k + n * j];

      Y[i + n * j] = Y[i + n * j] / X[idx[i] + n * i];
      // printf("%lf\n",X[i+n*i]);
    }
  }
}

void my_dcmp(DTYPE* X, UINT n, UINT* idx) {
  ITER i, j, k, imax;
  DTYPE max, cur;
  DTYPE t;
  for (i = 0; i < n; i++)
    idx[i] = i;  // Unit permutation matrix, P[N] initialized with N

  for (i = 0; i < n; i++)
    ;

  for (i = 0; i < n; i++) {
    max = .0;
    imax = i;
    for (k = 0; k < n; k++)
      if ((cur = fabs(X[idx[k] + n * i])) > max) {
        max = cur;
        t = X[idx[k] + n * i];
        imax = idx[k];
      }
    if (max < FZERO) {
      printf("FAiL\n");
    }  // failure, matrix is degenerate

    if (imax != i) {
      // pivoting P
      j = idx[i];
      idx[i] = idx[imax];
      idx[imax] = j;
      printf("%lf : %d <->  %d\n", t, idx[imax], idx[i]);
    }
  }

  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      if (i > j) {
        for (k = 0; k < j; k++)
          X[j * n + idx[i]] -= X[k * n + idx[i]] * X[j * n + idx[k]];
        X[j * n + idx[i]] /= X[j * n + idx[j]];
      } else {
        for (k = 0; k < i; k++)
          X[j * n + idx[i]] -= X[k * n + idx[i]] * X[j * n + idx[k]];
      }
    }
  }
}

void my_invert(MAT* X, MAT* Y) {
  UINT* idx;
  ITER i, j, k;
  UINT n, t = 0;

  n = X->d0;

  idx = mpalloc(sizeof(UINT) * n);

  my_dcmp(X->data, n, idx);
  printf("==== DECOMPOSiTON ====\n");
  for (i = 0; i < n; i++) printf("%u ", idx[i]);
  printf("\nLU MATRiX : \n");
  print_mat(X);

  wiki_invert(X->data, n, idx, Y->data);
  mpfree(idx);
}
