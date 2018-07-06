#include "iip_blas_lv1.h"

/*
 *
 *  y = alpha * x + y
 *  
 *  */
void axpy(DTYPE alpha, MAT* x,MAT* y)
{
#if DEBUG
printf("%s\n",__func__);
#endif
#if USE_CBLAS
/*
 *  <?>axpy(integer N, DTYPE alpha, DTYPE *x, integer incx, DTYPE beta )
 * */
//DTYPE = float
#if NTYPE == 0
	cblas_saxpy(x->d0,alpha,x->data,1,y->data,1);	

//DTYPE = double
#elif NTYPE == 1
	cblas_daxpy(x->d0,alpha,x->data,1,y->data,1);
#endif

//USE_BLAS = 0 -> just c implement
#else
	mp_axpy(x->d0, alpha, x->data,1,y->data,1);	
#endif
}

void mp_axpy(UINT N, DTYPE SA, DTYPE* SX, UINT INCX, DTYPE* SY, UINT INCY)
{
#if DEBUG
printf("%s\n",__func__);
#endif

ITER i;

#pragma omp parallel for shared(SX,SY) private(i)
for(i=0;i<N;i++)
	{
		SY[i*INCY] = SX[i*INCX] * SA + SY[i*INCY];	
	}

}

void caxpy(DTYPE alpha, CMAT* x,CMAT* y)
{
#if DEBUG
printf("%s\n",__func__);
#endif
#if USE_CBLAS
/*
 *  <?>axpy(integer N, DTYPE alpha, DTYPE *x, integer incx, DTYPE beta )
 * */
//DTYPE = float
#if NTYPE == 0
	cblas_zsaxpy(x->d0,&alpha,x->data,1,y->data,1);	

//DTYPE = double
#elif NTYPE == 1
	cblas_caxpy(x->d0,&alpha,x->data,1,y->data,1);
#endif

//USE_BLAS = 0 -> just c implement
#else
	mp_caxpy(x->d0, alpha, x->data,1,y->data,1);	
#endif
}

void mp_caxpy(UINT N, DTYPE SA, CTYPE* SX, UINT INCX, CTYPE* SY, UINT INCY)
{
#if DEBUG
printf("%s\n",__func__);
#endif
ITER i;

#pragma omp parallel for shared(SX,SY) private(i)
	for(i=0;i<N;i++)
	{
			SY[i*INCY].re = SX[i*INCX].re * SA + SY[i*INCY].re;
			SY[i*INCY].im = SX[i*INCX].im * SA + SY[i*INCY].im;
		
	}
}


void copy(MAT* src, SINT src_increment, MAT* des, SINT des_increment){
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
	cblas_scopy(x->d0,&alpha,x->data,src_increment,y->data,des_increment);	

//DTYPE = double
#elif NTYPE == 1
	cblas_dcopy(x->d0,&alpha,x->data,src_increment,y->data,des_increment);
#endif

//USE_BLAS = 0 -> just c implement
#else
	mp_copy(mat_size, src->data, src_increment, des->data, des_increment);	
#endif
}

void mp_copy(UINT N, DTYPE* src, SINT src_inc, DTYPE* des, SINT des_inc){
	ITER iteration = 8;
	UINT repeat = N>>3;
	UINT left = N&(UINT)(iteration-1);
	UINT i = 0;

	for(UINT j = 0; j < repeat; j++){
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

void ccopy(CMAT* src, SINT src_increment, CMAT* des, SINT des_increment){
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
	cblas_ccopy(x->d0,&alpha,x->data,src_increment,y->data,des_increment);	

//DTYPE = double
#elif NTYPE == 1
	cblas_zcopy(x->d0,&alpha,x->data,src_increment,y->data,des_increment);
#endif

//USE_BLAS = 0 -> just c implement
#else
	mp_ccopy(mat_size, src->data, src_increment, des->data, des_increment);	
#endif
}

void mp_ccopy(UINT N, CTYPE* src, SINT src_inc, CTYPE* des, SINT des_inc){
	ITER iteration = 8;
	UINT repeat = N>>3;
	UINT left = N&(UINT)(iteration-1);
	UINT i = 0;

	for(UINT j = 0; j < repeat; j++){
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
