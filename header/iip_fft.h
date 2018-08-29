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

#ifndef IIP_FFT_H
#define IIP_FFT_H

#include "iip_type.h"

/**** MKL HFFT ****/
/* N = in->d0
 * out->d0 = N/2 + 1
 */
/* example)
 * mkl_handle* my_handle;
 * MAT*A;
 * CMAT*C;
 * A = alloc_mat(1024);
 * C = alloc_cmat(1024/2 + 1);
 * read_mat("some_file",A);
 * my_handle = fft_handle(1024);
 * mkl_hfft(A,C); 
 * free_handle(my_handle);
 * */
/* Note)
 * If you want to parallelize mkl_fft.
 * you need to create multiple handles.
 * use each handle for each mkl_fft.
 * */
mkl_handle* fft_handle(UINT N);
void mkl_hfft(mkl_handle*handle,MAT*in,CMAT*out);
void mkl_hfft_col(mkl_handle*handle,DTYPE*in,CTYPE*out);

/**** MKL HIFFT ****/
/* N : out->d0
 * in ->d0 = N/2 + 1
 * */
mkl_handle* ifft_handle(UINT N);
void mkl_hifft(mkl_handle*handle,CMAT*in,MAT*out);
void mkl_hifft_col(mkl_handle*handle,CTYPE*in,DTYPE*out);

/**** MKL FFT & IFFT   ****/
/* use same handle with hfft.
 * */
/* example)
 * mkl_handle* my_handle;
 * MAT*A;
 * CMAT*C;
 * A = alloc_mat(1024);
 * C = alloc_cmat(1024);
 * read_mat("some_file",A);
 * my_handle = fft_handle(1024);
 * mkl_fft(A,C); 
 * free_handle(my_handle);
 * */

void mkl_fft(mkl_handle*handle,MAT*in,CMAT*out);
void mkl_ifft(mkl_handle*handle,CMAT*in,MAT*out);

void free_handle(mkl_handle*handle);

/*
    Ooura's FFt - fft4g.c
    Copyright(C) 1996-2001 Takuya OOURA
    email: ooura@mmm.t.u-tokyo.ac.jp
    download: http://momonga.t.u-tokyo.ac.jp/~ooura/fft.html
    You may use, copy, modify this code for any purpose and 
    without fee. You may distribute this ORIGINAL package.
*/
/* Fast Fourier Transform
 * perform fft on first dimension(d0).
 * d0 must be power of 2
 * Every fft function uses memory pool.
 * Be sure to initialize memory pool by init(<size>);
 * */
void fft(MAT*in,CMAT*out);
/*
 * increment in column for further application.
 * */
void ooura_fft_col(UINT N,DTYPE*in,CTYPE*out);

/* Inverse FFT */
void ifft(CMAT*in, MAT*out);
void ooura_ifft_col(UINT N,CTYPE*in,DTYPE*out);

/* Complex FFT */
void cfft(CMAT*in,MAT*out);
void ooura_cfft_col(UINT N,CTYPE*in,DTYPE*out);

/* Inverse Complex FFT */
void cifft(MAT*in,CMAT*out);
void ooura_cifft_col(UINT N,DTYPE*in,CTYPE*out);

/* Half FFT 
 *    | in    out
 * d0 | N     N/ 2 + 1
 *
 * */
void hfft(MAT*in, CMAT*out);
void ooura_hfft_col(UINT N,DTYPE*in, CTYPE*out);

/* Inverse Half FFT 
 *    | in    out
 * d0 | N     2*(N-1)
 * */
void hifft(CMAT*in, MAT*out);
void ooura_hifft_col(UINT N,CTYPE*in,DTYPE*out);

/*
-------- Complex DFT (Discrete Fourier Transform) --------
    [definition]
        <case1>
            X[k] = sum_j=0^n-1 x[j]*exp(2*pi*i*j*k/n), 0<=k<n
        <case2>
            X[k] = sum_j=0^n-1 x[j]*exp(-2*pi*i*j*k/n), 0<=k<n
        (notes: sum_j=0^n-1 is a summation from j=0 to n-1)
    [usage]
        <case1>
            ip[0] = 0; // first time only
            cdft(2*n, 1, a, ip, w);
        <case2>
            ip[0] = 0; // first time only
            cdft(2*n, -1, a, ip, w);
    [parameters]
        2*n            :data length (int)
                        n >= 1, n = power of 2
        a[0...2*n-1]   :input/output data (double *)
                        input data
                            a[2*j] = Re(x[j]), 
                            a[2*j+1] = Im(x[j]), 0<=j<n
                        output data
                            a[2*k] = Re(X[k]), 
                            a[2*k+1] = Im(X[k]), 0<=k<n
        ip[0...*]      :work area for bit reversal (int *)
                        length of ip >= 2+sqrt(n)
                        strictly, 
                        length of ip >= 
                            2+(1<<(int)(log(n+0.5)/log(2))/2).
                        ip[0],ip[1] are pointers of the cos/sin table.
        w[0...n/2-1]   :cos/sin table (double *)
                        w[],ip[] are initialized if ip[0] == 0.
    [remark]
        Inverse of 
            cdft(2*n, -1, a, ip, w);
        is 
            cdft(2*n, 1, a, ip, w);
            for (j = 0; j <= 2 * n - 1; j++) {
                a[j] *= 1.0 / n;
            }
*/
void cdft(int, int, double *, int *, double *);
/*
-------- Real DFT / Inverse of Real DFT --------
    [definition]
        <case1> RDFT
            R[k] = sum_j=0^n-1 a[j]*cos(2*pi*j*k/n), 0<=k<=n/2
            I[k] = sum_j=0^n-1 a[j]*sin(2*pi*j*k/n), 0<k<n/2
        <case2> IRDFT (excluding scale)
            a[k] = (R[0] + R[n/2]*cos(pi*k))/2 + 
                   sum_j=1^n/2-1 R[j]*cos(2*pi*j*k/n) + 
                   sum_j=1^n/2-1 I[j]*sin(2*pi*j*k/n), 0<=k<n
    [usage]
        <case1>
            ip[0] = 0; // first time only
            rdft(n, 1, a, ip, w);
        <case2>
            ip[0] = 0; // first time only
            rdft(n, -1, a, ip, w);
    [parameters]
        n              :data length (int)
                        n >= 2, n = power of 2
        a[0...n-1]     :input/output data (double *)
                        <case1>
                            output data
                                a[2*k] = R[k], 0<=k<n/2
                                a[2*k+1] = I[k], 0<k<n/2
                                a[1] = R[n/2]
                        <case2>
                            input data
                                a[2*j] = R[j], 0<=j<n/2
                                a[2*j+1] = I[j], 0<j<n/2
                                a[1] = R[n/2]
        ip[0...*]      :work area for bit reversal (int *)
                        length of ip >= 2+sqrt(n/2)
                        strictly, 
                        length of ip >= 
                            2+(1<<(int)(log(n/2+0.5)/log(2))/2).
                        ip[0],ip[1] are pointers of the cos/sin table.
        w[0...n/2-1]   :cos/sin table (double *)
                        w[],ip[] are initialized if ip[0] == 0.
    [remark]
        Inverse of 
            rdft(n, 1, a, ip, w);
        is 
            rdft(n, -1, a, ip, w);
            for (j = 0; j <= n - 1; j++) {
                a[j] *= 2.0 / n;
            }
        .
*/
void rdft(int, int, double *, int *, double *);
void ddct(int, int, double *, int *, double *);
void ddst(int, int, double *, int *, double *);
void dfct(int, double *, double *, int *, double *);
void dfst(int, double *, double *, int *, double *);

#endif
