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

/*check equal of two DTYPEs*/
UINT _eqdd(DTYPE A, DTYPE B);
/*compare 2 MAT struct*/
UINT compare_mat(MAT* A, MAT* B);
/*compare 2 CMAT struct*/
UINT compare_cmat(CMAT* A, CMAT* B);

void perform_test();
void do_test(char* filename);

/*char* operation : out = filename_post*/
void append_post(char* filename, const char* post, char* out);

#endif
