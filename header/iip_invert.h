#ifndef IIP_INVERT_H
#define IIP_INVERT_H

#include "iip_type.h"

/**** get inverse of matrix ****/
void invert_2by2(MAT* mat, MAT* inv);
void invert_3by3(MAT* mat, MAT* inv);
void invert_4by4(MAT* mat, MAT* inv);
void invert_5by5(MAT* mat, MAT* inv);
void invert_6by6(MAT* mat, MAT* inv);

/**** get determinant of matrix ****/
DTYPE det_2by2(MAT*mat);
DTYPE det_3by3(MAT*mat);
DTYPE det_4by4(MAT*mat);
DTYPE det_5by5(MAT*mat);
DTYPE det_6by6(MAT*mat);

#endif
