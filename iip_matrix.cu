#include "mother.h"


/**** zeros  ****/

MAT* zeros_1d(UINT d0)
{
#if DEBUG
printf("%s\n",__func__);
#endif
	MAT* mat = (MAT*)malloc(sizeof(MAT));
	mat->ndim = 0;
	mat->d0 = d0;
	mat->d1 = 1;
	mat->d2 = 1;
	cudaMalloc((void**)&(mat->data),sizeof(DTYPE) * d0);
    cudaMemset(mat->data,0,sizeof(DTYPE)*d0  );
	return mat;
}

MAT* zeros_2d(UINT d0,UINT d1)
{
#if DEBUG
printf("%s\n",__func__);
#endif
	MAT* mat = (MAT*)malloc(sizeof(MAT));
	mat->ndim = 1;
	mat->d0 = d0;
	mat->d1 = d1;
	mat->d2 = 1;
	
	cudaMalloc((void**)&(mat->data),sizeof(DTYPE) * d0 *d1);
    cudaMemset(mat->data,0,sizeof(DTYPE)*d0*d1);

	return mat;
}
MAT* zeros_3d(UINT d0,UINT d1,UINT d2)
{
#if DEBUG
printf("%s\n",__func__);
#endif
	MAT* mat = (MAT*)malloc(sizeof(MAT));
	mat->ndim = 2;
	mat->d0 = d0;
	mat->d1 = d1;
	mat->d2 = d2;

	cudaMalloc((void**)&(mat->data),sizeof(DTYPE) * d0 *d1 *d2);
    cudaMemset(mat->data,0,sizeof(DTYPE)*d0*d1*d2);
	return mat;
}


CMAT* czeros_1d(UINT d0)
{
#if DEBUG
printf("%s\n",__func__);
#endif
	CMAT* mat = (CMAT*)malloc(sizeof(CMAT));
	mat->ndim = 0;
	mat->d0 = d0;
	mat->d1 = 1;
	mat->d2 = 1;
	
	cudaMalloc((void**)&(mat->data),sizeof(CTYPE) * d0);
    cudaMemset(mat->data,0,sizeof(CTYPE)*d0);
	return mat;
}

CMAT* czeros_2d(UINT d0,UINT d1)
{
#if DEBUG
printf("%s\n",__func__);
#endif
	CMAT* mat = (CMAT*)malloc(sizeof(CMAT));
	mat->ndim = 1;
	mat->d0 = d0;
	mat->d1 = d1;
	mat->d2 = 1;
	
	cudaMalloc((void**)&(mat->data),sizeof(CTYPE) * d0 *d1);
    cudaMemset(mat->data,0,sizeof(CTYPE)*d0*d1);
	return mat;
}
CMAT* czeros_3d(UINT d0,UINT d1,UINT d2)
{
#if DEBUG
printf("%s\n",__func__);
#endif
	CMAT* mat = (CMAT*)malloc(sizeof(CMAT));
	mat->ndim = 2;
	mat->d0 = d0;
	mat->d1 = d1;
	mat->d2 = d2;
	
	cudaMalloc((void**)&(mat->data),sizeof(CTYPE) * d0 *d1*d2);
    cudaMemset(mat->data,0,sizeof(CTYPE)*d0*d1*d2);
	return mat;
}

/**** set ****/

void set_1d(MAT*mat, UINT idx0, DTYPE val)
{
#if DEBUG
printf("%s\n",__func__);
#endif
	mat->data[idx0] = val;
}
void set_2d(MAT*mat, UINT idx0, UINT idx1, DTYPE val)
{
#if DEBUG
printf("%s\n",__func__);
#endif
	mat->data[idx0 + (mat->d0)*idx1 ] = val;
  printf("i0: %d, i1: %d, val: %lf\n",idx0, idx1, val);

}
void set_3d(MAT*mat, UINT idx0, UINT idx1, UINT idx2, DTYPE val)
{
#if DEBUG
printf("%s\n",__func__);
#endif
	mat->data[idx0 + (mat->d0)*idx1 + (mat->d0) * (mat->d1) * idx2] = val;

}

void cset_1d(CMAT*mat, UINT idx0, DTYPE re, DTYPE im)
{
#if DEBUG
printf("%s\n",__func__);
#endif
	mat->data[idx0].re =re;
	mat->data[idx0].im =im;
}
void cset_2d(CMAT*mat, UINT idx0, UINT idx1, DTYPE re, DTYPE im)
{
#if DEBUG
printf("%s\n",__func__);
#endif
	mat->data[idx0 + (mat->d0)*idx1 ].re = re;
	mat->data[idx0 + (mat->d0)*idx1 ].im = im;

}
void cset_3d(CMAT*mat, UINT idx0, UINT idx1, UINT idx2, DTYPE re, DTYPE im)
{
#if DEBUG
printf("%s\n",__func__);
#endif
	mat->data[idx0 + (mat->d0)*idx1 + (mat->d0) * (mat->d1) * idx2].re = re;
	mat->data[idx0 + (mat->d0)*idx1 + (mat->d0) * (mat->d1) * idx2].im = im;

}

/**** fill ****/

void fill(MAT*mat, DTYPE val) // for real mat
{
  ITER i = 0;
  UINT len = mat->d0 * mat->d1 * mat->d2;
#if DEBUG
  printf("%s\n",__func__);
#endif
  for (i=0; i < len; i++)
    mat->data[i] = val;
}

void cfill(CMAT*cmat, DTYPE re, DTYPE im) // for complex mat
{
  ITER i = 0;
  UINT len = cmat->d0 * cmat->d1 * cmat->d2;
#if DEBUG
  printf("%s\n",__func__);
#endif
  for (i=0; i < len; i++)
  {
    cmat->data[i].re = re;
    cmat->data[i].im = im;
  }
}

/**** get ****/

DTYPE get_1d(MAT*mat, UINT idx0)
{
#if DEBUG
printf("%s\n",__func__);
#endif
	return mat->data[idx0];
}

DTYPE get_2d(MAT*mat, UINT idx0, UINT idx1)
{
#if DEBUG
printf("%s\n",__func__);
#endif
	return mat->data[idx0 +( mat->d0) * idx1];
}
DTYPE get_3d(MAT*mat, UINT idx0, UINT idx1, UINT idx2)
{
#if DEBUG
printf("%s\n",__func__);
#endif
	return mat->data[idx0 + (mat-> d0) * idx1 + (mat->d0) * (mat->d1) * idx0];
}

CTYPE cget_1d(CMAT*mat, UINT idx0)
{
#if DEBUG
printf("%s\n",__func__);
#endif
	return mat->data[idx0];
}

CTYPE cget_2d(CMAT*mat, UINT idx0, UINT idx1)
{
#if DEBUG
printf("%s\n",__func__);
#endif
	return mat->data[idx0 +( mat->d0) * idx1];
}
CTYPE cget_3d(CMAT*mat, UINT idx0, UINT idx1, UINT idx2)
{
#if DEBUG
printf("%s\n",__func__);
#endif
	return mat->data[idx0 + (mat-> d0) * idx1 + (mat->d0) * (mat->d1) * idx0];
}

/**** submat ****/

void submat_1d(MAT* mat, MAT* submat, ITER d0_st, ITER d0_ed)
{
#if DEBUG
  prITERf("%s\n",__func__);
#endif
  submat_3d(mat, submat, d0_st, d0_ed, -1, -1, -1, -1);
}


void submat_2d(MAT* mat, MAT* submat,
    ITER d0_st, ITER d0_ed,
    ITER d1_st, ITER d1_ed)
{
#if DEBUG
  prITERf("%s\n",__func__);
#endif
  submat_3d(mat, submat, d0_st, d0_ed, d1_st, d1_ed, -1, -1);
}


void submat_3d(MAT* mat, MAT* submat, 
    ITER d0_st, ITER d0_ed,
    ITER d1_st, ITER d1_ed,
    ITER d2_st, ITER d2_ed)
{
  ITER i, j, k;
  
#if DEBUG
  prITERf("%s\n",__func__);
#endif
  if (mat->d0 < submat->d0 || mat->d1 < submat->d1 || mat->d2 < submat->d2)
  {
    printf("error in \nmat->d0 < submat->d0 || mat->d1 < submat->d1 || mat->d2 < submat->d2\n");
    return;
  }

  if ( d0_st == -1 ) d0_st = 0;
  if ( d0_ed == -1 ) d0_ed = mat->d0;
  if ( d1_st == -1 ) d1_st = 0;
  if ( d1_ed == -1 ) d1_ed = mat->d1;
  if ( d2_st == -1 ) d2_st = 0;
  if ( d2_ed == -1 ) d2_ed = mat->d2;
  
  for (i = 0; d0_st+i < d0_ed; i++)
    for (j = 0; d1_st+j < d1_ed; j++)
      for (k = 0; d2_st+k < d2_ed; k++)
      {
        submat->data[i + j*(submat->d0) + k*(submat->d0*submat->d1)]
          = mat->data[(d0_st+i) + (d1_st+j)*(mat->d0) + (d2_st+k)*(mat->d0*mat->d1)];
        //printf("%d %d %d %d\n",i,j,k,i + j*(submat->d0) + k*(submat->d0*submat->d1));
        //printf("%d\n",(d0_st+i) + (d1_st+j)*(mat->d0) + (d2_st+k)*(mat->d0*mat->d1));
      }
  //printf("%d %d %d %d %d %d\n",d0_st, d0_ed, d1_st, d1_ed, d2_st, d2_ed);
}

void csubmat_1d(CMAT* mat, CMAT* submat, ITER d0_st, ITER d0_ed)
{
#if DEBUG
  prITERf("%s\n",__func__);
#endif
  csubmat_3d(mat, submat, d0_st, d0_ed, -1, -1, -1, -1);
}


void csubmat_2d(CMAT* mat, CMAT* submat,
    ITER d0_st, ITER d0_ed,
    ITER d1_st, ITER d1_ed)
{
#if DEBUG
  prITERf("%s\n",__func__);
#endif
  csubmat_3d(mat, submat, d0_st, d0_ed, d1_st, d1_ed, -1, -1);
}


void csubmat_3d(CMAT* mat, CMAT* submat, 
    ITER d0_st, ITER d0_ed,
    ITER d1_st, ITER d1_ed,
    ITER d2_st, ITER d2_ed)
{
  ITER i, j, k;
  
#if DEBUG
  prITERf("%s\n",__func__);
#endif
  if (mat->d0 < submat->d0 || mat->d1 < submat->d1 || mat->d2 < submat->d2)
  {
    printf("error in \nmat->d0 < submat->d0 || mat->d1 < submat->d1 || mat->d2 < submat->d2\n");
    return;
  }
  if ( d0_st == -1 ) d0_st = 0;
  if ( d0_ed == -1 ) d0_ed = mat->d0;
  if ( d1_st == -1 ) d1_st = 0;
  if ( d1_ed == -1 ) d1_ed = mat->d1;
  if ( d2_st == -1 ) d2_st = 0;
  if ( d2_ed == -1 ) d2_ed = mat->d2;
  
  for (i = 0; d0_st+i < d0_ed; i++)
    for (j = 0; d1_st+j < d1_ed; j++)
      for (k = 0; d2_st+k < d2_ed; k++)
      {
        submat->data[i + j*(submat->d0) + k*(submat->d0*submat->d1)].re = mat->data[(d0_st+i) + (d1_st+j)*(mat->d0) + (d2_st+k)*(mat->d0*mat->d1)].re;
        submat->data[i + j*(submat->d0) + k*(submat->d0*submat->d1)].im = mat->data[(d0_st+i) + (d1_st+j)*(mat->d0) + (d2_st+k)*(mat->d0*mat->d1)].im;
      }
}

/**** miscellaneous  ****/
void free_MAT(MAT *mat)
{
#if DEBUG
printf("%s\n",__func__);
#endif
	free(mat->data);
	free(mat);
}

void free_CMAT(CMAT *mat)
{
#if DEBUG
printf("%s\n",__func__);
#endif
	free(mat->data);
	free(mat);
}


void print_MAT(MAT* mat)
{
  ITER k, j, i;
#if DEBUG
printf("%s\n",__func__);
#endif
	for(k = 0; k < mat->d2; k++)
	{
		for(i=0; i < mat->d0 ; i++)	
		{
			for(j=0; j < mat->d1;j++)
				printf("%.3lf ",mat->data[k*(mat->d1)*(mat->d0) + j*(mat->d0) + i]);
			printf("\n");
		}
	printf("\n");		
	}
}

void print_CMAT(CMAT* mat)
{
  ITER k, j, i;
#if DEBUG
printf("%s\n",__func__);
#endif
	for(k = 0; k < mat->d2; k++)
	{
		for(i=0; i < mat->d0 ; i++)	
		{
			for(j=0; j < mat->d1;j++)
			{
				printf("%.3lf ",mat->data[k*(mat->d1)*(mat->d0) + j*(mat->d0) + i].re);
				printf("%.3lf| ",mat->data[k*(mat->d1)*(mat->d0) + j*(mat->d0) + i].im);
			
			}
			printf("\n");
		}
	printf("\n");		
	}
}
