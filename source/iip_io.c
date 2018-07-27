#include "iip_io.h"

/**** READ MATLAB .bin FILE ****/

void read_mat(const char* filename, MAT* mat) {
  FILE* f;
  ITER i = 0;
  f = NULL;
  f = fopen(filename, "rb");
  if (!f) ASSERT(NO_FILE);
  while (fread(&(mat->data[i++]), sizeof(DTYPE), 1, f) > 0)
    ;
  fclose(f);
}

void read_cmat(const char* filename, CMAT* mat) {
  FILE* f;
  ITER i = 0;
  f = NULL;
  f = fopen(filename, "rb");
  if (!f) ASSERT(NO_FILE);
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
  if (!f) ASSERT(NO_FILE);
  for (i = 0; i < mat->d0 * mat->d1 * mat->d2; i++)
    fwrite(&(mat->data[i]), sizeof(DTYPE), 1, f);
  fclose(f);
}

void write_cmat(const char* filename, CMAT* mat) {
  FILE* f;
  ITER i = 0;
  f = NULL;
  f = fopen(filename, "wb");
  if (!f) ASSERT(NO_FILE);
  for (i = 0; i < mat->d0 * mat->d1 * mat->d2; i++) {
    fwrite(&(mat->data[i].re), sizeof(CTYPE), 1, f);
  }
  fclose(f);
}
