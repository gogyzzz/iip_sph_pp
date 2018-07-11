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
void axpy(DTYPE alpha, MAT* x,MAT* y)
{
UINT size = x->d0 * x->d1 * x->d2;
#if DEBUG
printf("%s\n",__func__);
#endif

#if USE_CBLAS

	#if NTYPE == 0
	cblas_saxpy(size,alpha,x->data,1,y->data,1);	

	#elif NTYPE == 1
	cblas_daxpy(size,alpha,x->data,1,y->data,1);
	#endif

#else
	mp_axpy(size, alpha, x->data,1,y->data,1);	
#endif
}

void mp_axpy(UINT N, DTYPE alpha, DTYPE* X, UINT INCX, DTYPE* Y, UINT INCY)
{
ITER i;

#if DEBUG
printf("%s\n",__func__);
#endif

#pragma omp parallel for shared(X,Y) private(i)
for(i=0;i<N;i++)
	{
		Y[i*INCY] = X[i*INCX] * alpha + Y[i*INCY];	
	}

}

void caxpy(CTYPE alpha, CMAT* x,CMAT* y)
{
UINT size = x->d0 * x->d1 * x-> d2;
#if DEBUG
printf("%s\n",__func__);
#endif
#if USE_CBLAS
	#if NTYPE == 0
	cblas_caxpy(size,&alpha,x->data,1,y->data,1);	

	#elif NTYPE == 1
	cblas_zaxpy(size,&alpha,x->data,1,y->data,1);
	#endif

#else
	mp_caxpy(size, alpha, x->data,1,y->data,1);	
#endif
}

void mp_caxpy(UINT N, CTYPE alpha, CTYPE* X, UINT INCX, CTYPE* Y, UINT INCY)
{
#if DEBUG
printf("%s\n",__func__);
#endif
ITER i;

#pragma omp parallel for shared(X,Y) private(i)
	for(i=0;i<N;i++)
	{
			Y[i*INCY].re = X[i*INCX].re * alpha.re + Y[i*INCY].re;
			Y[i*INCY].im = X[i*INCX].im * alpha.im + Y[i*INCY].im;
		
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
void copy(MAT* src, MAT* des){
#if DEBUG
	printf("%s\n",__func__);
#endif
	UINT mat_size = src->d0 * src->d1 * src->d2;

	if(mat_size == 0)
		printf("Wrong MAT size!\n");

#if USE_CBLAS
#if NTYPE == 0
	cblas_scopy(mat_size,src->data,1,des->data,1);	

#elif NTYPE == 1
	cblas_dcopy(mat_size,src->data,1,des->data,1);
#endif

//USE_BLAS = 0 -> just c implement
#else
	mp_copy(mat_size, src->data, 1, des->data, 1);	
#endif
}

void mp_copy(UINT N, DTYPE* src, SINT src_inc, DTYPE* des, SINT des_inc){
	ITER iteration = 8;
	UINT repeat = N>>3;
	UINT left = N&(UINT)(iteration-1);
	UINT i = 0;
	UINT j = 0;

#pragma omp parallel for shared(des,src) private(j)
	for(j = 0; j < repeat; j++){
		des[(i    ) * des_inc] = src[(i    ) * src_inc];
		des[(i + 1) * des_inc] = src[(i + 1) * src_inc];
		des[(i + 2) * des_inc] = src[(i + 2) * src_inc];
		des[(i + 3) * des_inc] = src[(i + 3) * src_inc];
		des[(i + 4) * des_inc] = src[(i + 4) * src_inc];
		des[(i + 5) * des_inc] = src[(i + 5) * src_inc];
		des[(i + 6) * des_inc] = src[(i + 6) * src_inc];
		des[(i + 7) * des_inc] = src[(i + 7) * src_inc];
		i+= iteration;
	}

	for(UINT j=0; j<left; j++){
		des[(i + j) * des_inc] = src[(i + j) * src_inc];
	}
}

void ccopy(CMAT* src, CMAT* des){
#if DEBUG
	printf("%s\n",__func__);
#endif
	UINT mat_size = src->d0 * src->d1 * src->d2;

	if(mat_size == 0)
		printf("Wrong MAT size!\n");

#if USE_CBLAS
/*
 *  <?>axpy(integer N, DTYPE alpha, DTYPE *x, integer incx, DTYPE beta )
 * */
//DTYPE = float
#if NTYPE == 0
	cblas_ccopy(mat_size,src->data,1,des->data,1);	

//DTYPE = double
#elif NTYPE == 1
	cblas_zcopy(mat_size,src->data,1,des->data,1);
#endif

//USE_BLAS = 0 -> just c implement
#else
	mp_ccopy(mat_size, src->data, 1, des->data, 1);	
#endif
}

void mp_ccopy(UINT N, CTYPE* src, SINT src_inc, CTYPE* des, SINT des_inc){
	ITER iteration = 8;
	UINT repeat = N>>3;
	UINT left = N&(UINT)(iteration-1);
	UINT i = 0, j=0;
	

#pragma omp parallel for shared(des,src) private(j)
	for(j = 0; j < repeat; j++){
		des[(i    ) * des_inc].re = src[(i    ) * src_inc].re;
		des[(i    ) * des_inc].im = src[(i    ) * src_inc].im;
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
		i+= iteration;
	}

	for(UINT j=0; j<left; j++){
		des[(i + j) * des_inc].re = src[(i + j) * src_inc].re;
		des[(i + j) * des_inc].im = src[(i + j) * src_inc].im;
	}
}
