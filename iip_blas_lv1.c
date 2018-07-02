#include "mother.h"


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
 *
 *
 *
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
	for(ITER i=0;i<N;i++)
		SY[i*INCY] = SX[i*INCX] * SA + SY[i*INCY];	
		
}

