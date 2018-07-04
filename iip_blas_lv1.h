#include "mother.h"
/*********************
**** iip_blas_lv1 ****
**********************/

/* TODO
 * + common
 * asum
 * nrm2
 * rot
 * scal
 * swap
 * i?amax
 * i?amin
 * copy
 * -------axpy
 * + s,d
 * dot
 * rotm
 * rotmg
 * rotg
 * + c,z
 * dotc
 * dotu
 *
 * */
void mp_axpy(UINT,DTYPE,DTYPE*,UINT,DTYPE*,UINT);
void axpy(DTYPE,MAT*,MAT*);

void mp_caxpy(UINT,DTYPE,CTYPE*,UINT,CTYPE*,UINT);
void caxpy(DTYPE,CMAT*,CMAT*);


