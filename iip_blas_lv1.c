#include "iip_blas_lv1.h"

/*
 *
 * NTYPE == 0 -> float ->  <T> = real : s | complex : c
 * NTYPE == 1 -> double -> <T> = real : d | complex : z
 *  
 *
 */

/*********************************************
************ AXPY  ***************************
**********************************************/
/*
 *
 *  y = alpha * x + y
 *  
 *
 *  <?>axpy(integer N, DTYPE alpha, DTYPE *x, integer incx, DTYPE beta )
 * */
void axpy(DTYPE alpha, MAT *x, MAT *y)
{
	UINT size = x->d0 * x->d1 * x->d2;
#if DEBUG
	printf("%s\n", __func__);
#endif

#if USE_CBLAS

#if NTYPE == 0
	cblas_saxpy(size, alpha, x->data, 1, y->data, 1);

#elif NTYPE == 1
	cblas_daxpy(size, alpha, x->data, 1, y->data, 1);
#endif

#else
	mp_axpy(size, alpha, x->data, 1, y->data, 1);
#endif
}

void mp_axpy(UINT N, DTYPE alpha, DTYPE *X, UINT INCX, DTYPE *Y, UINT INCY)
{
	ITER i;

#if DEBUG
	printf("%s\n", __func__);
#endif

#pragma omp parallel for shared(X, Y) private(i)
	for (i = 0; i < N; i++)
	{
		Y[i * INCY] = X[i * INCX] * alpha + Y[i * INCY];
	}
}

void caxpy(CTYPE alpha, CMAT *x, CMAT *y)
{
	UINT size = x->d0 * x->d1 * x->d2;
#if DEBUG
	printf("%s\n", __func__);
#endif
#if USE_CBLAS
#if NTYPE == 0
	cblas_caxpy(size, &alpha, x->data, 1, y->data, 1);

#elif NTYPE == 1
	cblas_zaxpy(size, &alpha, x->data, 1, y->data, 1);
#endif

#else
	mp_caxpy(size, alpha, x->data, 1, y->data, 1);
#endif
}

void mp_caxpy(UINT N, CTYPE alpha, CTYPE *X, UINT INCX, CTYPE *Y, UINT INCY)
{
#if DEBUG
	printf("%s\n", __func__);
#endif
	ITER i;

#pragma omp parallel for shared(X, Y) private(i)
	for (i = 0; i < N; i++)
	{
		Y[i * INCY].re = X[i * INCX].re * alpha.re + Y[i * INCY].re;
		Y[i * INCY].im = X[i * INCX].im * alpha.im + Y[i * INCY].im;
	}
}

/*********************************************
************ COPY  ***************************
**********************************************/
/*  copies a vector, x, to a vector y,
 *
 *  y = x
 *
 *  uses unrolled loops for increments equal to 1.
 *
 *  <?>acopy(integer N, DTYPE* X, intefer INCX, DTYPE* Y, integer INCY)
 * */
void copy(MAT *src, MAT *des)
{
#if DEBUG
	printf("%s\n", __func__);
#endif
	UINT mat_size = src->d0 * src->d1 * src->d2;

	if (mat_size == 0)
	{
		printf("Wrong MAT size!\n");
		return 0;
	}

#if USE_CBLAS
#if NTYPE == 0
	cblas_scopy(mat_size, src->data, 1, des->data, 1);

#elif NTYPE == 1
	cblas_dcopy(mat_size, src->data, 1, des->data, 1);
#endif

//USE_BLAS = 0 -> just c implement
#else
	mp_copy(mat_size, src->data, 1, des->data, 1);
#endif
}

void mp_copy(UINT N, DTYPE *src, SINT src_inc, DTYPE *des, SINT des_inc)
{
	ITER iteration = 8;
	UINT repeat = N >> 3;
	UINT left = N & (UINT)(iteration - 1);
	UINT i = 0;
	UINT j = 0;

#pragma omp parallel for shared(des, src) private(j, i)
	for (j = 0; j < repeat; j++)
	{
		i = j * iteration;
		des[(i)*des_inc] = src[(i)*src_inc];
		des[(i + 1) * des_inc] = src[(i + 1) * src_inc];
		des[(i + 2) * des_inc] = src[(i + 2) * src_inc];
		des[(i + 3) * des_inc] = src[(i + 3) * src_inc];
		des[(i + 4) * des_inc] = src[(i + 4) * src_inc];
		des[(i + 5) * des_inc] = src[(i + 5) * src_inc];
		des[(i + 6) * des_inc] = src[(i + 6) * src_inc];
		des[(i + 7) * des_inc] = src[(i + 7) * src_inc];
	}

	for (j = 0; j < left; j++)
	{
		des[(i + j) * des_inc] = src[(i + j) * src_inc];
	}
}

void ccopy(CMAT *src, CMAT *des)
{
#if DEBUG
	printf("%s\n", __func__);
#endif
	UINT mat_size = src->d0 * src->d1 * src->d2;

	if (mat_size == 0)
	{
		printf("Wrong MAT size!\n");
		return 0;
	}

#if USE_CBLAS
/*
 *  <?>axpy(integer N, DTYPE alpha, DTYPE *x, integer incx, DTYPE beta )
 * */
//DTYPE = float
#if NTYPE == 0
	cblas_ccopy(mat_size, src->data, 1, des->data, 1);

//DTYPE = double
#elif NTYPE == 1
	cblas_zcopy(mat_size, src->data, 1, des->data, 1);
#endif

//USE_BLAS = 0 -> just c implement
#else
	mp_ccopy(mat_size, src->data, 1, des->data, 1);
#endif
}

void mp_ccopy(UINT N, CTYPE *src, SINT src_inc, CTYPE *des, SINT des_inc)
{
	ITER iteration = 8;
	UINT repeat = N >> 3;
	UINT left = N & (UINT)(iteration - 1);
	UINT i = 0, j = 0;

#pragma omp parallel for shared(des, src) private(j, i)
	for (j = 0; j < repeat; j++)
	{
		i = j * iteration;
		des[(i)*des_inc].re = src[(i)*src_inc].re;
		des[(i)*des_inc].im = src[(i)*src_inc].im;
		des[(i + 1) * des_inc].re = src[(i + 1) * src_inc].re;
		des[(i + 1) * des_inc].im = src[(i + 1) * src_inc].im;
		des[(i + 2) * des_inc].re = src[(i + 2) * src_inc].re;
		des[(i + 2) * des_inc].im = src[(i + 2) * src_inc].im;
		des[(i + 3) * des_inc].re = src[(i + 3) * src_inc].re;
		des[(i + 3) * des_inc].im = src[(i + 3) * src_inc].im;
		des[(i + 4) * des_inc].re = src[(i + 4) * src_inc].re;
		des[(i + 4) * des_inc].im = src[(i + 4) * src_inc].im;
		des[(i + 5) * des_inc].re = src[(i + 5) * src_inc].re;
		des[(i + 5) * des_inc].im = src[(i + 5) * src_inc].im;
		des[(i + 6) * des_inc].re = src[(i + 6) * src_inc].re;
		des[(i + 6) * des_inc].im = src[(i + 6) * src_inc].im;
		des[(i + 7) * des_inc].re = src[(i + 7) * src_inc].re;
		des[(i + 7) * des_inc].im = src[(i + 7) * src_inc].im;
	}

	for (j = 0; j < left; j++)
	{
		des[(i + j) * des_inc].re = src[(i + j) * src_inc].re;
		des[(i + j) * des_inc].im = src[(i + j) * src_inc].im;
	}
}

/*** Get sum of the magnitudes of elements of a vector ***/
DTYPE asum(MAT *mat, UINT inc)
{
#if DEBUG
	printf("%s\n", __func__);
#endif
	UINT mat_size = mat->d0 * mat->d1 * mat->d2;

	if (mat_size == 0)
	{
		printf("Wrong MAT size!\n");
		return 0;
	}

#if USE_CBLAS
/*
 *  <?>axpy(integer N, DTYPE alpha, DTYPE *x, integer incx, DTYPE beta )
 * */
//DTYPE = float
#if NTYPE == 0
	return cblas_sasum(mat_size, mat->data, inc);

//DTYPE = double
#elif NTYPE == 1
	return cblas_dasum(mat_size, mat->data, inc);
#endif

//USE_BLAS = 0 -> just c implement
#else
	return mp_asum(mat_size, mat->data, inc);
#endif
}
DTYPE mp_asum(UINT N, DTYPE *data, UINT inc)
{
	UINT i = 0;
	DTYPE sum = 0;

	for (i = 0; i < N; i++)
	{
		sum += (data[i * inc] < 0 ? (-data[i * inc]) : data[i * inc]);
	}

	return sum;
}

DTYPE casum(CMAT *mat, UINT inc)
{
#if DEBUG
	printf("%s\n", __func__);
#endif
	UINT mat_size = mat->d0 * mat->d1 * mat->d2;

	if (mat_size == 0)
	{
		printf("Wrong MAT size!\n");
		return 0;
	}

#if USE_CBLAS
/*
 *  <?>axpy(integer N, DTYPE alpha, DTYPE *x, integer incx, DTYPE beta )
 * */
//DTYPE = float
#if NTYPE == 0
	return cblas_scasum(mat_size, mat->data, inc);

//DTYPE = double
#elif NTYPE == 1
	return cblas_dzasum(mat_size, mat->data, inc);
#endif

//USE_BLAS = 0 -> just c implement
#else
	return mp_casum(mat_size, mat->data, inc);
#endif
}

DTYPE mp_casum(UINT N, CTYPE *data, UINT inc)
{
	UINT i = 0;
	DTYPE sum = 0;

	for (i = 0; i < N; i++)
	{
		sum += (data[i * inc].re < 0 ? (-data[i * inc].re) : data[i * inc].re);
		sum += (data[i * inc].im < 0 ? (-data[i * inc].im) : data[i * inc].im);
	}

	return sum;
}

/*** Get a vector-vector dot product ***/
DTYPE dot(MAT *src_x, UINT x_increment, MAT *src_y, UINT y_increment)
{
#if DEBUG
	printf("%s\n", __func__);
#endif
	UINT mat_size = src_x->d0 * src_x->d1 * src_x->d2;

	if (mat_size == 0)
	{
		printf("Wrong MAT size!\n");
		return 0;
	}

#if USE_CBLAS
/*
 *  <?>axpy(integer N, DTYPE alpha, DTYPE *x, integer incx, DTYPE beta )
 * */
//DTYPE = float
#if NTYPE == 0
	return cblas_sdot(N, src_x->data, x_increment, src_y->data, y_increment);

//DTYPE = double
#elif NTYPE == 1
	return cblas_ddot(N, src_x->data, x_increment, src_y->data, y_increment);
#endif

//USE_BLAS = 0 -> just c implement
#else
	return mp_dot(mat_size, src_x->data, x_increment, src_y->data, y_increment);
#endif
}

DTYPE mp_dot(UINT N, DTYPE *src_x, UINT x_inc, DTYPE *src_y, UINT y_inc)
{
	UINT i = 0;
	DTYPE dot = 0;

#pragma omp parallel for shared(src_x, src_y) private(i)
	for (i = 0; i < N; i++)
	{
		dot += src_x[i * x_inc] * src_y[i * y_inc];
	}

	return dot;
}

CTYPE cdot(CMAT *src_x, UINT x_increment, MAT *src_y, UINT y_increment){
#if DEBUG
	printf("%s\n", __func__);
#endif
	UINT mat_size = src_x->d0 * src_x->d1 * src_x->d2;
	CTYPE result = {0, 0};

	if (mat_size == 0)
	{
		printf("Wrong MAT size!\n");
		return result;
	}

#if USE_CBLAS
/*
 *  <?>axpy(integer N, DTYPE alpha, DTYPE *x, integer incx, DTYPE beta )
 * */
//DTYPE = float
#if NTYPE == 0
	cblas_cdotc_sub (N, src_x->data, x_increment, src_y->data, y_increment, &result);
	return result;

//DTYPE = double
#elif NTYPE == 1
	cblas_zdotc_sub (N, src_x->data, x_increment, src_y->data, y_increment, &result);
	return result;
#endif

//USE_BLAS = 0 -> just c implement
#else
	return mp_cdot(mat_size, src_x->data, x_increment, src_y->data, y_increment);
#endif
}

CTYPE mp_cdot(UINT N, CTYPE *src_x, UINT x_inc, DTYPE *src_y, UINT y_inc){
	UINT i = 0;
	CTYPE dot;

	dot.re = 0;
	dot.im = 0;

#pragma omp parallel for shared(src_x, src_y) private(i)
	for (i = 0; i < N; i++)
	{
		dot.re += src_x[i * x_inc].re * src_y[i * y_inc];
		dot.im += src_x[i * x_inc].im * src_y[i * y_inc];
	}

	return dot;
}


CTYPE udot(CMAT *src_x, UINT x_increment, CMAT *src_y, UINT y_increment){
#if DEBUG
	printf("%s\n", __func__);
#endif
	UINT mat_size = src_x->d0 * src_x->d1 * src_x->d2;
	CTYPE result = {0, 0};

	if (mat_size == 0)
	{
		printf("Wrong MAT size!\n");
		return result;
	}

#if USE_CBLAS
/*
 *  <?>axpy(integer N, DTYPE alpha, DTYPE *x, integer incx, DTYPE beta )
 * */
//DTYPE = float
#if NTYPE == 0
	cblas_cdotu_sub (N, src_x->data, x_increment, src_y->data, y_increment, &result);
	return result;

//DTYPE = double
#elif NTYPE == 1
	cblas_zdotu_sub (N, src_x->data, x_increment, src_y->data, y_increment, &result);
	return result;
#endif

//USE_BLAS = 0 -> just c implement
#else
	return mp_udot(mat_size, src_x->data, x_increment, src_y->data, y_increment);
#endif
}

CTYPE mp_udot(UINT N, CTYPE *src_x, UINT x_inc, CTYPE *src_y, UINT y_inc){
	UINT i = 0;
	CTYPE dot;

	dot.re = 0;
	dot.im = 0;

#pragma omp parallel for shared(src_x, src_y) private(i)
	for (i = 0; i < N; i++)
	{
		dot.re += src_x[i * x_inc].re * src_y[i * y_inc].re;
		dot.im += src_x[i * x_inc].re * src_y[i * y_inc].im;
		dot.re -= src_x[i * x_inc].im * src_y[i * y_inc].im;
		dot.im += src_x[i * x_inc].im * src_y[i * y_inc].re;
	}

	return dot;
}