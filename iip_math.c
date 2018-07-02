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
	
	mat->data = (DTYPE*)malloc(sizeof(DTYPE)*d0);
	memset(mat->data,0,sizeof(DTYPE)*d0);

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
	
	mat->data=(DTYPE*)malloc(sizeof(DTYPE)*d0*d1);
	memset(mat->data,0,sizeof(DTYPE)*d0*d1);

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
	
	mat->data=(DTYPE*)malloc(sizeof(DTYPE)*d0*d1*d2);
	memset(mat->data,0,sizeof(DTYPE)*d0*d1*d2);

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


void submat_1d(MAT* mat, MAT* submat, int d0_st, int d0_ed)
{
#if DEBUG
  printf("%s\n",__func__);
#endif
  submat_3d(mat, submat, d0_st, d0_ed, -1, -1, -1, -1);
}
void submat_2d(MAT* mat, MAT* submat,
    int d0_st, int d0_ed,
    int d1_st, int d1_ed)
{
#if DEBUG
  printf("%s\n",__func__);
#endif
  submat_3d(mat, submat, d0_st, d0_ed, d1_st, d1_ed, -1, -1);
}
void submat_3d(MAT* mat, MAT* submat, 
    int d0_st, int d0_ed,
    int d1_st, int d1_ed,
    int d2_st, int d2_ed)
{
  ITER i, j, k;
  
#if DEBUG
  printf("%s\n",__func__);
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

/**** miscellaneous  ****/
void free_MAT(MAT *mat)
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
				printf("%lf ",mat->data[k*(mat->d1)*(mat->d0) + j*(mat->d0) + i]);
			printf("\n");
		}
	printf("\n");		
	}
}

