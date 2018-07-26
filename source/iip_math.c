#include "iip_math.h"


/**** SQUARE ROOT ****/
void sqrt_mat(MAT*mat)
{
	sqrt_mat_inc(mat->d0*mat->d1*mat->d2,mat->data,1);	
}

void sqrt_mat_inc(UINT size,DTYPE*X,ITER inc)
{
	ITER i;
#if DEBUG
printf("%s\n",__func__);
#endif

#pragma omp parallel for shared(X) private(i)
	for(i=0;i<size;i+=inc)
	{
#if NTYPE == 0
	X[i]=sqrtf(X[i]);				
#elif NTYPE == 1
	X[i]=sqrt(X[i]);
#endif
	}
}

void sqrt_cmat(CMAT*mat)
{
	sqrt_cmat_inc(mat->d0*mat->d1*mat->d2,mat->data,1);	
}

void sqrt_cmat_inc(UINT size,CTYPE*X,ITER incx)
{
	ITER i;
#if DEBUG
printf("%s\n",__func__);
#endif

#pragma omp parallel for shared(X) private(i)
	for(i=0; i < size; i+=incx)
	{
#if OS_WIN
	#if NTYPE == 0
	*(_Fcomplex*)(&X[i]) = csqrtf(*(_Fcomplex*)(&X[i])); 
	#elif NTYPE == 1
	*(_Dcomplex*)(&X[i]) = csqrt(*(_Dcomplex*)(&X[i])); 
	#endif
#elif OS_UNIX
	#if NTYPE == 0
	*(complex float*)(&X[i]) = csqrtf(*(complex float*)(&X[i])); 
	#elif NTYPE == 1
	*(complex double*)(&X[i]) = csqrt(*(complex double*)(&X[i])); 
	#endif
#endif
	}
}


/**** POWER ****/

void pow_mat(MAT*mat, DTYPE n)
{
	pow_mat_inc(mat->d0 * mat->d1 * mat->d2, mat->data, n, 1);
}

void pow_mat_inc(UINT size, DTYPE*X, DTYPE n, ITER incx)
{
	ITER i;
#if DEBUG
	printf("%s\n",__func__);
#endif

#pragma omp parallel for shared(X) private(i)
	for(i=0; i < size; i+=incx)
	{
#if OS_WIN
	#if NTYPE == 0
	*(_Fcomplex*)(&X[i]) = powf(*(_Fcomplex*)(&X[i]),n); 
	#elif NTYPE == 1
	*(_Dcomplex*)(&X[i]) = pow(*(_Dcomplex*)(&X[i]),n); 
	#endif
#elif OS_UNIX
	#if NTYPE == 0
	*(complex float*)(&X[i]) = powf(*(complex float*)(&X[i]),n); 
	#elif NTYPE == 1
	*(complex double*)(&X[i]) = pow(*(complex double*)(&X[i]),n); 
	#endif
#endif
	}
}
void pow_cmat(CMAT*mat, DTYPE n)
{
	pow_cmat_inc(mat->d0 * mat->d1 * mat->d2, mat->data, n, 1);
}

void pow_cmat_inc(UINT size, CTYPE*X, DTYPE n, ITER incx)
{
	ITER i;
	CTYPE temp;
	temp.re = n;
	temp.im = 0;
#if DEBUG
	printf("%s\n",__func__);
#endif

#pragma omp parallel for shared(X) private(i)
	for(i=0; i < size; i+=incx)
	{
#if OS_WIN
	#if NTYPE == 0
	*(_Fcomplex*)(&X[i]) = cpowf(*(_Fcomplex*)(&X[i]),*(_Fcomplex*)(&temp)); 
	#elif NTYPE == 1
	*(_Dcomplex*)(&X[i]) = cpow(*(_Dcomplex*)(&X[i]),*(_Dcomplex*)(&temp)); 
	#endif
#elif OS_UNIX
	#if NTYPE == 0
	*(complex float*)(&X[i]) = cpowf(*(complex float*)(&X[i]),*(complex float*)(&temp)); 
	#elif NTYPE == 1
	*(complex double*)(&X[i]) = cpow(*(complex double*)(&X[i]),*(complex double*)(&temp)); 
	#endif
#endif
	}
}

void cpow_cmat(CMAT*mat, CTYPE n)
{
	cpow_cmat_inc(mat->d0 * mat->d1 * mat->d2, mat->data, n, 1);
}

void cpow_cmat_inc(UINT size, CTYPE*X, CTYPE n, ITER incx)
{
	ITER i;
#if DEBUG
	printf("%s\n",__func__);
#endif

#pragma omp parallel for shared(X) private(i)
	for(i=0; i < size; i+=incx)
	{
#if OS_WIN
	#if NTYPE == 0
	*(_Fcomplex*)(&X[i]) = cpowf(*(_Fcomplex*)(&X[i]),*(_Fcomplex*)(&n)); 
	#elif NTYPE == 1
	*(_Dcomplex*)(&X[i]) = pow(*(_Dcomplex*)(&X[i]),*(_Dcomplex*)(&n)); 
	#endif
#elif OS_UNIX
	#if NTYPE == 0
	*(complex float*)(&X[i]) = cpowf(*(complex float*)(&X[i]),*(complex float*)(&n)); 
	#elif NTYPE == 1
	*(complex double*)(&X[i]) = pow(*(complex double*)(&X[i]),*(complex double*)(&n)); 

	#endif
#endif
	}
}

/**** Uniform distribution ****/

void randu(MAT*mat, DTYPE a, DTYPE b)
{	
	randu_inc(mat->d0 *mat->d1 * mat->d2,mat->data,a,b,1);
}

void randu_inc(UINT size, DTYPE*X,DTYPE a,DTYPE b, ITER incx)
{
	ITER i;
#if DEBUG
	printf("%s\n",__func__);
#endif
	srand(time(NULL));

#pragma omp parallel for shared(X) private(i)
	for(i=0;i<size;i+=incx)
	{
#if NTYPE == 0
	X[i]=((float)rand()/RAND_MAX )*(b-a)+a;
#elif NTYPE == 1
	X[i]=((double)rand()/RAND_MAX )*(b-a)+a;
#endif
	}
}

void crandu(CMAT*mat, DTYPE ra, DTYPE rb, DTYPE ia, DTYPE ib)
{	
	crandu_inc(mat->d0 *mat->d1 * mat->d2,mat->data,ra,rb,ia,ib,1);
}

void crandu_inc(UINT size, CTYPE*X,DTYPE ra,DTYPE rb,DTYPE ia,DTYPE ib, ITER incx)
{
	ITER i;
#if DEBUG
	printf("%s\n",__func__);
#endif
	srand(time(NULL));
#pragma omp parallel for shared(X) private(i)
	for(i=0;i<size;i+=incx)
	{
#if NTYPE == 0
	X[i].re=((float)rand()/RAND_MAX )*(rb-ra)+ra;
	X[i].im=((float)rand()/RAND_MAX )*(ib-ia)+ia;
#elif NTYPE == 1
	X[i].re=((double)rand()/RAND_MAX )*(rb-ra)+ra;
	X[i].im=((double)rand()/RAND_MAX )*(ib-ia)+ia;
#endif
	}
}

void randn(MAT*mat,DTYPE mean,DTYPE std)
{
	randn_inc(mat->d0*mat->d1*mat->d2,mat->data,mean,std,1);
}
/*
 * Using Box-Muller Transform, But I heard ziggurat is fastest
 * */
void randn_inc(UINT size,DTYPE*X,DTYPE mean,DTYPE std,ITER incx)
{
	ITER i;
	DTYPE u,v;
	DTYPE s;
#if DEBUG
	printf("%s\n",__func__);
#endif

	srand(time(NULL));
#pragma omp parallel for shared(X) private(i,u,v,s)
	for(i=0;i<size;i+=incx)
	{
		do
		{
			u = (2.0*rand()/RAND_MAX)-1;	
			v = (2.0*rand()/RAND_MAX)-1;
			s= u*u + v*v;	
		}
		//check if s is in unit circle or is zero
		while(s >=1.0 || s == 0);	
		
		//You can use u,v both
#if NTYPE == 0
		X[i]=std*(u*sqrtf(-2.0*logf(s)/s)) + mean;
#elif NTYPE == 1
		X[i]=std*(u*sqrt(-2.0*log(s)/s)) + mean;
#endif
	}	
}

void crandn(CMAT*mat,CTYPE mean,CTYPE std)
{
	crandn_inc(mat->d0*mat->d1*mat->d2,mat->data,mean,std,1);
}
void crandn_inc(UINT size,CTYPE*X,CTYPE mean,CTYPE std,ITER incx)
{
	ITER i;
	DTYPE u,v;
	DTYPE s;
#if DEBUG
	printf("%s\n",__func__);
#endif

	srand(time(NULL));
#pragma omp parallel for shared(X) private(i,u,v,s)
	for(i=0;i<size;i+=incx)
	{
		do
		{
			u = (2.0*rand()/RAND_MAX)-1;	
			v = (2.0*rand()/RAND_MAX)-1;
			s= u*u + v*v;	
		}
		//check if s is in unit circle or is zero
		while(s >=1.0 || s == 0);	
		
		//You can use u,v both
#if NTYPE == 0
		X[i].re=std.re*(u*sqrtf(-2.0*logf(s)/s)) + mean.re;
		X[i].im=std.im*(v*sqrtf(-2.0*logf(s)/s)) + mean.im;
#elif NTYPE == 1
		X[i].re=std.re*(u*sqrt(-2.0*log(s)/s)) + mean.re;
		X[i].im=std.im*(u*sqrt(-2.0*log(s)/s)) + mean.im;
#endif
	}	
}


/**** round ****/
void round_mat(MAT*mat)
{
	round_mat_inc(mat->d0*mat->d1*mat->d2,mat->data,1);
}
void round_mat_inc(UINT size, DTYPE*X,ITER incx)
{
	ITER i;
#pragma omp parallel for shared(X) private(i)
	for(i=0;i<size;i+=incx)
	{
#if NTYPE == 0
		X[i]=roundf(X[i]);
#elif NTYPE == 1
		X[i]=round(X[i]);
#endif
	}
}

void round_cmat(CMAT*mat)
{
	round_cmat_inc(mat->d0*mat->d1*mat->d2,mat->data,1);
}
void round_cmat_inc(UINT size, CTYPE*X,ITER incx)
{
	ITER i;
#pragma omp parallel for shared(X) private(i)
	for(i=0;i<size;i+=incx)
	{
#if NTYPE == 0
		X[i].re=roundf(X[i].re);
		X[i].im=roundf(X[i].im);
#elif NTYPE == 1
		X[i].re=round(X[i].re);
		X[i].im=round(X[i].im);
#endif
	}
}
/**** floor ****/
void floor_mat(MAT*mat)
{
	floor_mat_inc(mat->d0*mat->d1*mat->d2,mat->data,1);
}
void floor_mat_inc(UINT size, DTYPE*X,ITER incx)
{
	ITER i;
#pragma omp parallel for shared(X) private(i)
	for(i=0;i<size;i+=incx)
	{
#if NTYPE == 0
		X[i]=floorf(X[i]);
#elif NTYPE == 1
		X[i]=floor(X[i]);
#endif
	}
}

void floor_cmat(CMAT*mat)
{
	floor_cmat_inc(mat->d0*mat->d1*mat->d2,mat->data,1);
}
void floor_cmat_inc(UINT size, CTYPE*X,ITER incx)
{
	ITER i;
#pragma omp parallel for shared(X) private(i)
	for(i=0;i<size;i+=incx)
	{
#if NTYPE == 0
		X[i].re=floorf(X[i].re);
		X[i].im=floorf(X[i].im);
#elif NTYPE == 1
		X[i].re=floor(X[i].re);
		X[i].im=floor(X[i].im);
#endif
	}
}
/**** ceil ****/
void ceil_mat(MAT*mat)
{
	ceil_mat_inc(mat->d0*mat->d1*mat->d2,mat->data,1);
}
void ceil_mat_inc(UINT size, DTYPE*X,ITER incx)
{
	ITER i;
#pragma omp parallel for shared(X) private(i)
	for(i=0;i<size;i+=incx)
	{
#if NTYPE == 0
		X[i]=ceilf(X[i]);
#elif NTYPE == 1
		X[i]=ceil(X[i]);
#endif
	}
}

void ceil_cmat(CMAT*mat)
{
	ceil_cmat_inc(mat->d0*mat->d1*mat->d2,mat->data,1);
}
void ceil_cmat_inc(UINT size, CTYPE*X,ITER incx)
{
	ITER i;
#pragma omp parallel for shared(X) private(i)
	for(i=0;i<size;i+=incx)
	{
#if NTYPE == 0
		X[i].re=ceilf(X[i].re);
		X[i].im=ceilf(X[i].im);
#elif NTYPE == 1
		X[i].re=ceil(X[i].re);
		X[i].im=ceil(X[i].im);
#endif
	}
}
/**** log ****/
void log_mat(MAT*mat)
{
	log_mat_inc(mat->d0*mat->d1*mat->d2,mat->data,1);
}
void log_mat_inc(UINT size, DTYPE*X,ITER incx)
{
	ITER i;
#pragma omp parallel for shared(X) private(i)
	for(i=0;i<size;i+=incx)
	{
#if NTYPE == 0
		X[i]=logf(X[i]);
#elif NTYPE == 1
		X[i]=log(X[i]);
#endif
	}
}

void log_cmat(CMAT*mat)
{
	log_cmat_inc(mat->d0*mat->d1*mat->d2,mat->data,1);
}
void log_cmat_inc(UINT size, CTYPE*X,ITER incx)
{
	ITER i;
#pragma omp parallel for shared(X) private(i)
	for(i=0;i<size;i+=incx)
	{
#if NTYPE == 0
		CXF(X[i])=clogf( CXF(X[i]) );
#elif NTYPE == 1
		CXD(X[i])=clogf( CXD(X[i]) );
#endif
	}
}

void log2_mat(MAT*mat)
{
	log2_mat_inc(mat->d0*mat->d1*mat->d2,mat->data,1);
}
void log2_mat_inc(UINT size, DTYPE*X,ITER incx)
{
	ITER i;
#pragma omp parallel for shared(X) private(i)
	for(i=0;i<size;i+=incx)
	{
#if NTYPE == 0
		X[i]=log2f(X[i]);
#elif NTYPE == 1
		X[i]=log2(X[i]);
#endif
	}
}

void log2_cmat(CMAT*mat)
{
	log2_cmat_inc(mat->d0*mat->d1*mat->d2,mat->data,1);
}
void log2_cmat_inc(UINT size, CTYPE*X,ITER incx)
{
	ITER i;
#pragma omp parallel for shared(X) private(i)
	for(i=0;i<size;i+=incx)
	{
#if NTYPE == 0
		CXF(X[i])=clog2f( CXF(X[i]) );
#elif NTYPE == 1
		CXD(X[i])=clog2f( CXD(X[i]) );
#endif
	}
}

void log10_mat(MAT*mat)
{
	log10_mat_inc(mat->d0*mat->d1*mat->d2,mat->data,1);
}
void log10_mat_inc(UINT size, DTYPE*X,ITER incx)
{
	ITER i;
#pragma omp parallel for shared(X) private(i)
	for(i=0;i<size;i+=incx)
	{
#if NTYPE == 0
		X[i]=log10f(X[i]);
#elif NTYPE == 1
		X[i]=log10(X[i]);
#endif
	}
}

void log10_cmat(CMAT*mat)
{
	log10_cmat_inc(mat->d0*mat->d1*mat->d2,mat->data,1);
}
void log10_cmat_inc(UINT size, CTYPE*X,ITER incx)
{
	ITER i;
#pragma omp parallel for shared(X) private(i)
	for(i=0;i<size;i+=incx)
	{
#if NTYPE == 0
		CXF(X[i])=clog10f( CXF(X[i]) );
#elif NTYPE == 1
		CXD(X[i])=clog10f( CXD(X[i]) );
#endif
	}
}
/**** exp ****/
void exp_mat(MAT*mat)
{
	exp_mat_inc(mat->d0*mat->d1*mat->d2,mat->data,1);
}
void exp_mat_inc(UINT size, DTYPE*X,ITER incx)
{
	ITER i;
#pragma omp parallel for shared(X) private(i)
	for(i=0;i<size;i+=incx)
	{
#if NTYPE == 0
		X[i]=expf(X[i]);
#elif NTYPE == 1
		X[i]=exp(X[i]);
#endif
	}
}

void exp_cmat(CMAT*mat)
{
	exp_cmat_inc(mat->d0*mat->d1*mat->d2,mat->data,1);
}
void exp_cmat_inc(UINT size, CTYPE*X,ITER incx)
{
	ITER i;
#pragma omp parallel for shared(X) private(i)
	for(i=0;i<size;i+=incx)
	{
#if NTYPE == 0
		CXF(X[i])=cexpf( CXF(X[i]) );
#elif NTYPE == 1
		CXD(X[i])=cexpf( CXD(X[i]) );
#endif
	}
}
/**** abs ****/


