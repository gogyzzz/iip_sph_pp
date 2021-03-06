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
#ifndef IIP_INVERT_H
#define IIP_INVERT_H

#include "iip_matrix.h"
#include "iip_type.h"

/**** get inverse matrix of mat  ****/
void invert(MAT* mat, MAT* inv);
void cinvert(CMAT* mat, CMAT* inv);

void invert_2by2(DTYPE* X, DTYPE* Y);
void cinvert_2by2(CTYPE* X, CTYPE* Y);

void invert_3by3(DTYPE* X, DTYPE* Y);
void cinvert_3by3(CTYPE* X, CTYPE* Y);

void invert_4by4(DTYPE* X, DTYPE* Y);
void cinvert_4by4(CTYPE* X, CTYPE* Y);

void invert_5by5(DTYPE* X, DTYPE* Y);
void cinvert_5by5(CTYPE* X, CTYPE* Y);

void invert_6by6(DTYPE* X, DTYPE* Y);
void cinvert_6by6(CTYPE* X, CTYPE* Y);

/**** get determinant of mat ****/
DTYPE det(MAT* mat);
CTYPE cdet(CMAT* mat);

DTYPE det_2by2(DTYPE* X);
CTYPE cdet_2by2(CTYPE* X);

DTYPE det_3by3(DTYPE* X);
CTYPE cdet_3by3(CTYPE* X);

DTYPE det_4by4(DTYPE* X);
CTYPE cdet_4by4(CTYPE* X);

DTYPE det_5by5(DTYPE* X);
CTYPE cdet_5by5(CTYPE* X);

DTYPE det_6by6(DTYPE* X);
CTYPE cdet_6by6(CTYPE* X);

/**** inversion by LAPACK ****/
void invert_nbyn(DTYPE* X, DTYPE* Y, UINT n);
void cinvert_nbyn(CTYPE* X, CTYPE* Y, UINT n);

DTYPE det_nbyn(MAT* mat);
CTYPE cdet_nbyn(CMAT* mat);

/**** Work In Progress ****/
void wiki_invert(DTYPE* X, UINT n, UINT* idx, DTYPE* Y);
void mp_dcmp(DTYPE* X, UINT n, UINT* idx);
void my_invert(MAT* X, MAT* Y);

#endif
