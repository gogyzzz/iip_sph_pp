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

/**** LAPACK ****/
void invert_nbyn(MAT*mat,MAT*inv);
DTYPE det_nbyn(MAT*mat);


/**** LU ****/
void lu_decomposition(DTYPE*X,UINT n,UINT* idx);

void lu_backwardsubstitution(DTYPE*lu,UINT n,UINT* idx, DTYPE*B);

void lu_invert(MAT*X,MAT*Y);

DTYPE lu_det(MAT*X);


void cig_lu(DTYPE* A, int n);

void cig_lusolve(DTYPE* LU, DTYPE* b, int n);

void lu_dcmp(DTYPE*X,UINT n);

void wiki_invert(DTYPE*X, UINT n,UINT* idx,DTYPE*Y);

void wiki_dcmp(DTYPE*X, UINT n,UINT*idx);

void wiki(MAT*X,MAT*Y);

#endif
