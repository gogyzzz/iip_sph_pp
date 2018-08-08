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

/**** READ MATLAB .bin FILE ****/

void read_mat(const char* filename, MAT* mat) {
  FILE* f;
  ITER i = 0;
  f = NULL;
  f = fopen(filename, "rb");
  ASSERT_FILE(f, filename)
  while (fread(&(mat->data[i++]), sizeof(DTYPE), 1, f) > 0)
    ;
  fclose(f);
}

void read_cmat(const char* filename, CMAT* mat) {
  FILE* f;
  ITER i = 0;
  f = NULL;
  f = fopen(filename, "rb");
  ASSERT_FILE(f, filename)
  while (fread(&(mat->data[i++]), sizeof(CTYPE), 1, f) > 0)
    ;
  fclose(f);
}

/**** WRITE MATLAB .bin FILE ****/
void write_mat(const char* filename, MAT* mat) {
  FILE* f;
  ITER i = 0;
  f = NULL;
  f = fopen(filename, "wb");
  ASSERT_FILE(f, filename)
  for (i = 0; i < mat->d0 * mat->d1 * mat->d2; i++)
    fwrite(&(mat->data[i]), sizeof(DTYPE), 1, f);
  fclose(f);
}

void write_cmat(const char* filename, CMAT* mat) {
  FILE* f;
  ITER i = 0;
  f = NULL;
  f = fopen(filename, "wb");
  ASSERT_FILE(f, filename)
  for (i = 0; i < mat->d0 * mat->d1 * mat->d2; i++) {
    fwrite(&(mat->data[i].re), sizeof(CTYPE), 1, f);
  }
  fclose(f);
}
