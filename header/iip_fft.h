#ifndef IIP_FFT_H
#define IIP_FFT_H

#include "iip_type.h"

/*
 *
 *
 * */
void fft(MAT*in,CMAT*out);
/*
 *
 *
 * */
void hfft(MAT*in, CMAT*out);

/*
    Ooura's FFt - fft4g
    Copyright:
    Copyright(C) 1996-2001 Takuya OOURA
    email: ooura@mmm.t.u-tokyo.ac.jp
    download: http://momonga.t.u-tokyo.ac.jp/~ooura/fft.html
    You may use, copy, modify this code for any purpose and 
    without fee. You may distribute this ORIGINAL package.
*/
void cdft(int, int, double *, int *, double *);
void rdft(int, int, double *, int *, double *);
void ddct(int, int, double *, int *, double *);
void ddst(int, int, double *, int *, double *);
void dfct(int, double *, double *, int *, double *);
void dfst(int, double *, double *, int *, double *);

#endif
