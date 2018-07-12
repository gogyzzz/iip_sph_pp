#include "iip_blas_lv1.h"

/*
 *
 *  y = alpha * x + y
 *  
 *  */
void axpy(DTYPE alpha, MAT *x, MAT *y)
{
#if DEBUG
	printf("%s\n", __func__);
#endif
#if USE_CBLAS
/*
 *  <?>axpy(integer N, DTYPE alpha, DTYPE *x, integer incx, DTYPE beta )
 * */
//DTYPE = float
#if NTYPE == 0
	cblas_saxpy(x->d0, alpha, x->data, 1, y->data, 1);

//DTYPE = double
#elif NTYPE == 1
	cblas_daxpy(x->d0, alpha, x->data, 1, y->data, 1);
#endif

//USE_BLAS = 0 -> just c implement
#else
	mp_axpy(x->d0, alpha, x->data, 1, y->data, 1);
#endif
}

void mp_axpy(UINT N, DTYPE SA, DTYPE *SX, UINT INCX, DTYPE *SY, UINT INCY)
{
#if DEBUG
	printf("%s\n", __func__);
#endif

	ITER i;

#pragma omp parallel for shared(SX, SY) private(i)
	for (i = 0; i < N; i++)
	{
		SY[i * INCY] = SX[i * INCX] * SA + SY[i * INCY];
	}
}

void caxpy(DTYPE alpha, CMAT *x, CMAT *y)
{
#if DEBUG
	printf("%s\n", __func__);
#endif
#if USE_CBLAS
/*
 *  <?>axpy(integer N, DTYPE alpha, DTYPE *x, integer incx, DTYPE beta )
 * */
//DTYPE = float
#if NTYPE == 0
	cblas_zsaxpy(x->d0, &alpha, x->data, 1, y->data, 1);

//DTYPE = double
#elif NTYPE == 1
	cblas_caxpy(x->d0, &alpha, x->data, 1, y->data, 1);
#endif

//USE_BLAS = 0 -> just c implement
#else
	mp_caxpy(x->d0, alpha, x->data, 1, y->data, 1);
#endif
}

void mp_caxpy(UINT N, DTYPE SA, CTYPE *SX, UINT INCX, CTYPE *SY, UINT INCY)
{
#if DEBUG
	printf("%s\n", __func__);
#endif
	ITER i;

#pragma omp parallel for shared(SX, SY) private(i)
	for (i = 0; i < N; i++)
	{
		SY[i * INCY].re = SX[i * INCX].re * SA + SY[i * INCY].re;
		SY[i * INCY].im = SX[i * INCX].im * SA + SY[i * INCY].im;
	}
}

/*** Copy vector from src to des ***/
void copy(MAT *src, SINT src_increment, MAT *des, SINT des_increment)
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
	cblas_scopy(x->d0, &alpha, x->data, src_increment, y->data, des_increment);

//DTYPE = double
#elif NTYPE == 1
	cblas_dcopy(x->d0, &alpha, x->data, src_increment, y->data, des_increment);
#endif

//USE_BLAS = 0 -> just c implement
#else
	mp_copy(mat_size, src->data, src_increment, des->data, des_increment);
#endif
}

void mp_copy(UINT N, DTYPE *src, SINT src_inc, DTYPE *des, SINT des_inc)
{
	ITER iteration = 8;
	UINT repeat = N >> 3;
	UINT left = N & (UINT)(iteration - 1);
	UINT i = 0, j = 0;

	i = 0;
	for (j = 0; j < repeat; j++)
	{
		des[(i)*des_inc] = src[(i)*src_inc];
		des[(i + 1) * des_inc] = src[(i + 1) * src_inc];
		des[(i + 2) * des_inc] = src[(i + 2) * src_inc];
		des[(i + 3) * des_inc] = src[(i + 3) * src_inc];
		des[(i + 4) * des_inc] = src[(i + 4) * src_inc];
		des[(i + 5) * des_inc] = src[(i + 5) * src_inc];
		des[(i + 6) * des_inc] = src[(i + 6) * src_inc];
		des[(i + 7) * des_inc] = src[(i + 7) * src_inc];
		i += iteration;
	}

	for (j = 0; j < left; j++)
	{
		des[(i + j) * des_inc] = src[(i + j) * src_inc];
	}
}

void ccopy(CMAT *src, SINT src_increment, CMAT *des, SINT des_increment)
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
	cblas_ccopy(x->d0, &alpha, x->data, src_increment, y->data, des_increment);

//DTYPE = double
#elif NTYPE == 1
	cblas_zcopy(x->d0, &alpha, x->data, src_increment, y->data, des_increment);
#endif

//USE_BLAS = 0 -> just c implement
#else
	mp_ccopy(mat_size, src->data, src_increment, des->data, des_increment);
#endif
}

void mp_ccopy(UINT N, CTYPE *src, SINT src_inc, CTYPE *des, SINT des_inc)
{
	ITER iteration = 8;
	UINT repeat = N >> 3;
	UINT left = N & (UINT)(iteration - 1);
	UINT i = 0, j = 0;

	i = 0;
	for (j = 0; j < repeat; j++)
	{
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
		i += iteration;
	}

	for (j = 0; j < left; j++)
	{
		des[(i + j) * des_inc].re = src[(i + j) * src_inc].re;
		des[(i + j) * des_inc].im = src[(i + j) * src_inc].im;
	}
}

/*** Get sum of the magnitudes of elements of a vector ***/
DTYPE asum(UINT N, MAT *mat, UINT inc)
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
	return mp_asum(N, mat->data, inc);
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

DTYPE casum(UINT N, CMAT *mat, UINT inc)
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
	return mp_casum(N, mat->data, inc);
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
DTYPE dot(UINT N, MAT *src_x, UINT x_increment, MAT *src_y, UINT y_increment)
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
	return mp_dot(N, src_x->data, x_increment, src_y->data, y_increment);
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


DTYPE cdot(UINT N, CMAT *src_x, UINT x_increment, CMAT *src_y, UINT y_increment);
DTYPE mp_cdot(UINT N, CTYPE *src_x, UINT x_inc, CTYPE *src_y, UINT y_inc);