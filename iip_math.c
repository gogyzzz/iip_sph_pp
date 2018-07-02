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

MAT* zeros_2d(UINT d1,UINT d0)
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
MAT* zeros_3d(UINT d2,UINT d1,UINT d0)
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

void set_1d(MAT*mat, UINT idx0, UINT val)
{
#if DEBUG
printf("%s\n",__func__);
#endif
	mat->data[idx0] = val;
}
void set_2d(MAT*mat, UINT idx1, UINT idx0, UINT val)
{
#if DEBUG
printf("%s\n",__func__);
#endif
	mat->data[idx0 + (mat->d0)*idx1 ] = val;

}
void set_3d(MAT*mat, UINT idx2, UINT idx1, UINT idx0,UINT val)
{
#if DEBUG
printf("%s\n",__func__);
#endif
	mat->data[idx0 + (mat->d0)*idx1 + (mat->d0) * (mat->d1) * idx2] = val;

}

/**** get ****/

DTYPE get_1d(MAT*mat, UINT idx0)
{
#if DEBUG
printf("%s\n",__func__);
#endif
	return mat->data[idx0];
}

DTYPE get_2d(MAT*mat, UINT idx1, UINT idx0)
{
#if DEBUG
printf("%s\n",__func__);
#endif
	return mat->data[idx0 +( mat->d0) * idx1];
}
DTYPE get_3d(MAT*mat, UINT idx2, UINT idx1, UINT idx0)
{
#if DEBUG
printf("%s\n",__func__);
#endif
	return mat->data[idx0 + (mat-> d0) * idx1 + (mat->d0) * (mat->d1) * idx0];
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
#if DEBUG
printf("%s\n",__func__);
#endif
	for(ITER k = 0; k < mat->d2; k++)
	{
		for(ITER j=0; j < mat->d1 ; j++)	
		{
			for(ITER i=0;i < mat->d0;i++)
				printf("%lf ",mat->data[k*(mat->d1)*(mat->d0) + j*(mat->d0)  + i]);
			printf("\n");
		}
	printf("\n");		
	}
}

