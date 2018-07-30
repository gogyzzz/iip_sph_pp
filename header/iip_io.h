#ifndef IIP_IO_H
#define IIP_IO_H

#include "iip_type.h"

/*** READ MATLAB .bin FILE ****/
void read_mat(const char* filename, MAT* mat);
void read_cmat(const char* filename, CMAT* mat);

/*** WRTIE MATLAB .bin FILE ****/
void write_mat(const char* filename, MAT* mat);
void write_cmat(const char* filename, CMAT* mat);
#endif
