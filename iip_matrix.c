//#include "mother.h"
#include "iip_matrix.h"

/**** alloc_MAT ****/

MAT *alloc_MAT_1d(UINT d0)
{
	MAT *mat;
#if DEBUG
	printf("%s\n", __func__);
#endif

	mat = (MAT *)malloc(sizeof(MAT));
	mat->ndim = 0;
	mat->d0 = d0;
	mat->d1 = 1;
	mat->d2 = 1;

	mat->data = (DTYPE *)malloc(sizeof(DTYPE) * d0);

	return mat;
}
MAT *alloc_MAT_2d(UINT d0, UINT d1)
{
	MAT *mat;
#if DEBUG
	printf("%s\n", __func__);
#endif

	mat = (MAT *)malloc(sizeof(MAT));
	mat->ndim = 1;
	mat->d0 = d0;
	mat->d1 = d1;
	mat->d2 = 1;

	mat->data = (DTYPE *)malloc(sizeof(DTYPE) * d0 * d1);

	return mat;
}
MAT *alloc_MAT_3d(UINT d0, UINT d1, UINT d2)
{
	MAT *mat;
#if DEBUG
	printf("%s\n", __func__);
#endif

	mat = (MAT *)malloc(sizeof(MAT));
	mat->ndim = 2;
	mat->d0 = d0;
	mat->d1 = d1;
	mat->d2 = d2;

	mat->data = (DTYPE *)malloc(sizeof(DTYPE) * d0 * d1 * d2);

	return mat;
}

CMAT *alloc_CMAT_1d(UINT d0)
{
	CMAT *mat;
#if DEBUG
	printf("%s\n", __func__);
#endif

	mat = (CMAT *)malloc(sizeof(CMAT));
	mat->ndim = 0;
	mat->d0 = d0;
	mat->d1 = 1;
	mat->d2 = 1;

	mat->data = (CTYPE *)malloc(sizeof(CTYPE) * d0);

	return mat;
}
CMAT *alloc_CMAT_2d(UINT d0, UINT d1)
{
	CMAT *mat;
#if DEBUG
	printf("%s\n", __func__);
#endif

	mat = (CMAT *)malloc(sizeof(CMAT));
	mat->ndim = 1;
	mat->d0 = d0;
	mat->d1 = d1;
	mat->d2 = 1;

	mat->data = (CTYPE *)malloc(sizeof(CTYPE) * d0 * d1);

	return mat;
}
CMAT *alloc_CMAT_3d(UINT d0, UINT d1, UINT d2)
{
	CMAT *mat;
#if DEBUG
	printf("%s\n", __func__);
#endif

	mat = (CMAT *)malloc(sizeof(CMAT));
	mat->ndim = 2;
	mat->d0 = d0;
	mat->d1 = d1;
	mat->d2 = d2;

	mat->data = (CTYPE *)malloc(sizeof(CTYPE) * d0*d1*d2);

	return mat;
}

/**** zeros  ****/

MAT *zeros_1d(UINT d0)
{
#if DEBUG
	printf("%s\n", __func__);
#endif
	MAT *mat = (MAT *)malloc(sizeof(MAT));
	mat->ndim = 0;
	mat->d0 = d0;
	mat->d1 = 1;
	mat->d2 = 1;

	mat->data = (DTYPE *)malloc(sizeof(DTYPE) * d0);
	memset(mat->data, 0, sizeof(DTYPE) * d0);

	return mat;
}

MAT *zeros_2d(UINT d0, UINT d1)
{
#if DEBUG
	printf("%s\n", __func__);
#endif
	MAT *mat = (MAT *)malloc(sizeof(MAT));
	mat->ndim = 1;
	mat->d0 = d0;
	mat->d1 = d1;
	mat->d2 = 1;

	mat->data = (DTYPE *)malloc(sizeof(DTYPE) * d0 * d1);
	memset(mat->data, 0, sizeof(DTYPE) * d0 * d1);

	return mat;
}
MAT *zeros_3d(UINT d0, UINT d1, UINT d2)
{
#if DEBUG
	printf("%s\n", __func__);
#endif
	MAT *mat = (MAT *)malloc(sizeof(MAT));
	mat->ndim = 2;
	mat->d0 = d0;
	mat->d1 = d1;
	mat->d2 = d2;

	mat->data = (DTYPE *)malloc(sizeof(DTYPE) * d0 * d1 * d2);
	memset(mat->data, 0, sizeof(DTYPE) * d0 * d1 * d2);

	return mat;
}

CMAT *czeros_1d(UINT d0)
{
#if DEBUG
	printf("%s\n", __func__);
#endif
	CMAT *mat = (CMAT *)malloc(sizeof(CMAT));
	mat->ndim = 0;
	mat->d0 = d0;
	mat->d1 = 1;
	mat->d2 = 1;

	mat->data = (CTYPE *)malloc(sizeof(CTYPE) * d0);
	memset(mat->data, 0, sizeof(CTYPE) * d0);

	return mat;
}

CMAT *czeros_2d(UINT d0, UINT d1)
{
#if DEBUG
	printf("%s\n", __func__);
#endif
	CMAT *mat = (CMAT *)malloc(sizeof(CMAT));
	mat->ndim = 1;
	mat->d0 = d0;
	mat->d1 = d1;
	mat->d2 = 1;

	mat->data = (CTYPE *)malloc(sizeof(CTYPE) * d0 * d1);
	memset(mat->data, 0, sizeof(CTYPE) * d0 * d1);

	return mat;
}
CMAT *czeros_3d(UINT d0, UINT d1, UINT d2)
{
#if DEBUG
	printf("%s\n", __func__);
#endif
	CMAT *mat = (CMAT *)malloc(sizeof(CMAT));
	mat->ndim = 2;
	mat->d0 = d0;
	mat->d1 = d1;
	mat->d2 = d2;

	mat->data = (CTYPE *)malloc(sizeof(CTYPE) * d0 * d1 * d2);
	memset(mat->data, 0, sizeof(CTYPE) * d0 * d1 * d2);

	return mat;
}

/**** set ****/

void set_1d(MAT *mat, UINT idx0, DTYPE val)
{
#if DEBUG
	printf("%s\n", __func__);
#endif
	mat->data[idx0] = val;
}
void set_2d(MAT *mat, UINT idx0, UINT idx1, DTYPE val)
{
#if DEBUG
	printf("%s\n", __func__);
#endif
	mat->data[idx0 + (mat->d0) * idx1] = val;
	printf("i0: %d, i1: %d, val: %lf\n", idx0, idx1, val);
}
void set_3d(MAT *mat, UINT idx0, UINT idx1, UINT idx2, DTYPE val)
{
#if DEBUG
	printf("%s\n", __func__);
#endif
	mat->data[idx0 + (mat->d0) * idx1 + (mat->d0) * (mat->d1) * idx2] = val;
}

void cset_1d(CMAT *mat, UINT idx0, DTYPE re, DTYPE im)
{
#if DEBUG
	printf("%s\n", __func__);
#endif
	mat->data[idx0].re = re;
	mat->data[idx0].im = im;
}
void cset_2d(CMAT *mat, UINT idx0, UINT idx1, DTYPE re, DTYPE im)
{
#if DEBUG
	printf("%s\n", __func__);
#endif
	mat->data[idx0 + (mat->d0) * idx1].re = re;
	mat->data[idx0 + (mat->d0) * idx1].im = im;
}
void cset_3d(CMAT *mat, UINT idx0, UINT idx1, UINT idx2, DTYPE re, DTYPE im)
{
#if DEBUG
	printf("%s\n", __func__);
#endif
	mat->data[idx0 + (mat->d0) * idx1 + (mat->d0) * (mat->d1) * idx2].re = re;
	mat->data[idx0 + (mat->d0) * idx1 + (mat->d0) * (mat->d1) * idx2].im = im;
}

/**** fill ****/

void fill(MAT *mat, DTYPE val) // for real mat
{
	ITER i = 0;
	UINT len = mat->d0 * mat->d1 * mat->d2;
#if DEBUG
	printf("%s\n", __func__);
#endif
	for (i = 0; i < len; i++)
		mat->data[i] = val;
}

void cfill(CMAT *cmat, DTYPE re, DTYPE im) // for complex mat
{
	ITER i = 0;
	UINT len = cmat->d0 * cmat->d1 * cmat->d2;
#if DEBUG
	printf("%s\n", __func__);
#endif
	for (i = 0; i < len; i++)
	{
		cmat->data[i].re = re;
		cmat->data[i].im = im;
	}
}

/**** get ****/

DTYPE get_1d(MAT *mat, UINT idx0)
{
#if DEBUG
	printf("%s\n", __func__);
#endif
	return mat->data[idx0];
}

DTYPE get_2d(MAT *mat, UINT idx0, UINT idx1)
{
#if DEBUG
	printf("%s\n", __func__);
#endif
	return mat->data[idx0 + (mat->d0) * idx1];
}
DTYPE get_3d(MAT *mat, UINT idx0, UINT idx1, UINT idx2)
{
#if DEBUG
	printf("%s\n", __func__);
#endif
	return mat->data[idx0 + (mat->d0) * idx1 + (mat->d0) * (mat->d1) * idx0];
}

CTYPE cget_1d(CMAT *mat, UINT idx0)
{
#if DEBUG
	printf("%s\n", __func__);
#endif
	return mat->data[idx0];
}

CTYPE cget_2d(CMAT *mat, UINT idx0, UINT idx1)
{
#if DEBUG
	printf("%s\n", __func__);
#endif
	return mat->data[idx0 + (mat->d0) * idx1];
}
CTYPE cget_3d(CMAT *mat, UINT idx0, UINT idx1, UINT idx2)
{
#if DEBUG
	printf("%s\n", __func__);
#endif
	return mat->data[idx0 + (mat->d0) * idx1 + (mat->d0) * (mat->d1) * idx0];
}

/**** submat ****/

void submat_1d(MAT *mat, MAT *submat, ITER d0_st, ITER d0_ed)
{
#if DEBUG
	printf("%s\n", __func__);
#endif
	submat_3d(mat, submat, d0_st, d0_ed, -1, -1, -1, -1);
}

void submat_2d(MAT *mat, MAT *submat,
			   ITER d0_st, ITER d0_ed,
			   ITER d1_st, ITER d1_ed)
{
#if DEBUG
	printf("%s\n", __func__);
#endif
	submat_3d(mat, submat, d0_st, d0_ed, d1_st, d1_ed, -1, -1);
}

void submat_3d(MAT *mat, MAT *submat,
			   ITER d0_st, ITER d0_ed,
			   ITER d1_st, ITER d1_ed,
			   ITER d2_st, ITER d2_ed)
{
	ITER i, j, k;

#if DEBUG
	printf("%s\n", __func__);
#endif
	if (mat->d0 < submat->d0 || mat->d1 < submat->d1 || mat->d2 < submat->d2)
	{
		printf("error in \nmat->d0 < submat->d0 || mat->d1 < submat->d1 || mat->d2 < submat->d2\n");
		return;
	}

	if (d0_st == -1)
		d0_st = 0;
	if (d0_ed == -1)
		d0_ed = mat->d0;
	if (d1_st == -1)
		d1_st = 0;
	if (d1_ed == -1)
		d1_ed = mat->d1;
	if (d2_st == -1)
		d2_st = 0;
	if (d2_ed == -1)
		d2_ed = mat->d2;

	for (i = 0; d0_st + i < d0_ed; i++)
		for (j = 0; d1_st + j < d1_ed; j++)
			for (k = 0; d2_st + k < d2_ed; k++)
			{
				submat->data[i + j * (submat->d0) + k * (submat->d0 * submat->d1)] = mat->data[(d0_st + i) + (d1_st + j) * (mat->d0) + (d2_st + k) * (mat->d0 * mat->d1)];
				//printf("%d %d %d %d\n",i,j,k,i + j*(submat->d0) + k*(submat->d0*submat->d1));
				//printf("%d\n",(d0_st+i) + (d1_st+j)*(mat->d0) + (d2_st+k)*(mat->d0*mat->d1));
			}
	//printf("%d %d %d %d %d %d\n",d0_st, d0_ed, d1_st, d1_ed, d2_st, d2_ed);
}

void csubmat_1d(CMAT *mat, CMAT *submat, ITER d0_st, ITER d0_ed)
{
#if DEBUG
	printf("%s\n", __func__);
#endif
	csubmat_3d(mat, submat, d0_st, d0_ed, -1, -1, -1, -1);
}

void csubmat_2d(CMAT *mat, CMAT *submat,
				ITER d0_st, ITER d0_ed,
				ITER d1_st, ITER d1_ed)
{
#if DEBUG
	printf("%s\n", __func__);
#endif
	csubmat_3d(mat, submat, d0_st, d0_ed, d1_st, d1_ed, -1, -1);
}

void csubmat_3d(CMAT *mat, CMAT *submat,
				ITER d0_st, ITER d0_ed,
				ITER d1_st, ITER d1_ed,
				ITER d2_st, ITER d2_ed)
{
	ITER i, j, k;

#if DEBUG
	printf("%s\n", __func__);
#endif
	if (mat->d0 < submat->d0 || mat->d1 < submat->d1 || mat->d2 < submat->d2)
	{
		printf("error in \nmat->d0 < submat->d0 || mat->d1 < submat->d1 || mat->d2 < submat->d2\n");
		return;
	}
	if (d0_st == -1)
		d0_st = 0;
	if (d0_ed == -1)
		d0_ed = mat->d0;
	if (d1_st == -1)
		d1_st = 0;
	if (d1_ed == -1)
		d1_ed = mat->d1;
	if (d2_st == -1)
		d2_st = 0;
	if (d2_ed == -1)
		d2_ed = mat->d2;

	for (i = 0; d0_st + i < d0_ed; i++)
		for (j = 0; d1_st + j < d1_ed; j++)
			for (k = 0; d2_st + k < d2_ed; k++)
			{
				submat->data[i + j * (submat->d0) + k * (submat->d0 * submat->d1)].re = mat->data[(d0_st + i) + (d1_st + j) * (mat->d0) + (d2_st + k) * (mat->d0 * mat->d1)].re;
				submat->data[i + j * (submat->d0) + k * (submat->d0 * submat->d1)].im = mat->data[(d0_st + i) + (d1_st + j) * (mat->d0) + (d2_st + k) * (mat->d0 * mat->d1)].im;
			}
}


/**** transpose ****/
MAT* trans(MAT* mat)
{
	ITER i,j;
	UINT d0=mat->d0;
	UINT d1=mat->d1;
	UINT d2=mat->d2;
	
	MAT * t_mat;
#if DEBUG
	printf("%s\n", __func__);
#endif
	if(mat->ndim == 0)
	{
   t_mat = alloc_MAT(d1,d0);
#pragma omp parallel for shared(t_mat,mat) private(i)
	for(i=0;i<d0;i++)
		{
	 		t_mat->data[i] = mat->data[i];
		}	 
	}
	else if(mat->ndim == 1)
	{
		t_mat = alloc_MAT(d1,d0);
#pragma omp parallel for shared(t_mat,mat) private(i)
		for(i=0;i<d0 * d1  ; i++)
		{
			t_mat->data[i/d0 + i%d0*d1 ] = mat->data[i];
		}
	}
	else
	{
		t_mat = alloc_MAT(d1,d0,d2);
//일단은 3차원 배열을 2차원 배열의 batch로 여기고 구현함
#pragma omp parallel for shared(t_mat,mat) private(i,j)
		for(j=0;j< d2;j++)
		{
			for(i=0;i<d0 * d1; i++)
				{
					t_mat->data[j*d0*d1 + i/d0 + i%d0*d1 ] = mat->data[j*d0*d1 + i];
				}	
		}		
	}
	return t_mat;
}

CMAT* ctrans(CMAT* mat)
{
	ITER i,j;
	UINT d0=mat->d0;
	UINT d1=mat->d1;
	UINT d2=mat->d2;
	
	CMAT * t_mat;
#if DEBUG
	printf("%s\n", __func__);
#endif
	if(mat->ndim == 0)
	{
   t_mat = alloc_CMAT(d1,d0);
#pragma omp parallel for shared(t_mat,mat) private(i)
	for(i=0;i<d0;i++)
		{
	 		t_mat->data[i].re = mat->data[i].re;
	 		t_mat->data[i].im = mat->data[i].im;
		}	 
	}
	else if(mat->ndim == 1)
	{
		t_mat = alloc_CMAT(d1,d0);
#pragma omp parallel for shared(t_mat,mat) private(i)
		for(i=0;i<d0 * d1  ; i++)
		{
			t_mat->data[i/d0 + i%d0*d1 ].re = mat->data[i].re;
			t_mat->data[i/d0 + i%d0*d1 ].im = mat->data[i].im;
		}
	}
	else
	{
		t_mat = alloc_CMAT(d1,d0,d2);
//일단은 3차원 배열을 2차원 배열의 batch로 여기고 구현함
#pragma omp parallel for shared(t_mat,mat) private(i,j)
		for(j=0;j< d2;j++)
		{
			for(i=0;i<d0 * d1; i++)
				{
					t_mat->data[j*d0*d1 + i/d0 + i%d0*d1 ].re = mat->data[j*d0*d1 + i].re;
					t_mat->data[j*d0*d1 + i/d0 + i%d0*d1 ].im = mat->data[j*d0*d1 + i].im;
				}	
		}		
	}
	return t_mat;
}

CMAT* hermit(CMAT* mat)
{
	ITER i,j;
	UINT d0=mat->d0;
	UINT d1=mat->d1;
	UINT d2=mat->d2;
	
	CMAT * t_mat;
#if DEBUG
	printf("%s\n", __func__);
#endif
	if(mat->ndim == 0)
	{
   t_mat = alloc_CMAT(d1,d0);
#pragma omp parallel for shared(t_mat,mat) private(i)
	for(i=0;i<d0;i++)
		{
	 		t_mat->data[i].re = mat->data[i].re;
	 		t_mat->data[i].im = -(mat->data[i].im);
		}	 
	}
	else if(mat->ndim == 1)
	{
		t_mat = alloc_CMAT(d1,d0);
#pragma omp parallel for shared(t_mat,mat) private(i)
		for(i=0;i<d0 * d1  ; i++)
		{
			t_mat->data[i/d0 + i%d0*d1 ].re = mat->data[i].re;
			t_mat->data[i/d0 + i%d0*d1 ].im = -(mat->data[i].im);
		}
	}
	else
	{
		t_mat = alloc_CMAT(d1,d0,d2);
//일단은 3차원 배열을 2차원 배열의 batch로 여기고 구현함
#pragma omp parallel for shared(t_mat,mat) private(i,j)
		for(j=0;j< d2;j++)
		{
			for(i=0;i<d0 * d1; i++)
				{
					t_mat->data[j*d0*d1 + i/d0 + i%d0*d1 ].re = mat->data[j*d0*d1 + i].re;
					t_mat->data[j*d0*d1 + i/d0 + i%d0*d1 ].im = -(mat->data[j*d0*d1 + i].im);
				}	
		}		
	}
	return t_mat;
}

/**** miscellaneous  ****/
void free_MAT(MAT *mat)
{
#if DEBUG
	printf("%s\n", __func__);
#endif
	free(mat->data);
	free(mat);
}

void free_CMAT(CMAT *mat)
{
#if DEBUG
	printf("%s\n", __func__);
#endif
	free(mat->data);
	free(mat);
}

void print_MAT(MAT *mat)
{
	ITER k, j, i;
#if DEBUG
	printf("%s\n", __func__);
#endif
	for (k = 0; k < mat->d2; k++)
	{
		for (i = 0; i < mat->d0; i++)
		{
			for (j = 0; j < mat->d1; j++)
				printf("%.3lf ", mat->data[k * (mat->d1) * (mat->d0) + j * (mat->d0) + i]);
			printf("\n");
		}
		printf("\n");
	}
}

void print_CMAT(CMAT *mat)
{
	ITER k, j, i;
#if DEBUG
	printf("%s\n", __func__);
#endif
	for (k = 0; k < mat->d2; k++)
	{
		for (i = 0; i < mat->d0; i++)
		{
			for (j = 0; j < mat->d1; j++)
			{
				printf("%.3lf ", mat->data[k * (mat->d1) * (mat->d0) + j * (mat->d0) + i].re);
				printf("%.3lf| ", mat->data[k * (mat->d1) * (mat->d0) + j * (mat->d0) + i].im);
			}
			printf("\n");
		}
		printf("\n");
	}
}
