#include "mother.h"

MAT* zeros_1d(UINT d0)
{
	MAT* mat = (MAT*)malloc(sizeof(MAT));
	mat->ndim = 0;
	mat->d0 = d0;
	mat->d1 = -1;
	mat->d2 = -1;
	
	mat->data = (DTYPE*)malloc(sizeof(DTYPE)*d0);
	for(ITER i=0;i<d0;i++)
		mat->data[i] = 0;

	return mat;
}

MAT* zeros_2d(UINT d0,UINT d1)
{
	MAT* mat = (MAT*)malloc(sizeof(MAT));
	mat->ndim = 1;
	mat->d0 = d0;
	mat->d1 = d1;
	mat->d2 = -1;
	
	mat->data=(DTYPE*)malloc(sizeof(DTYPE)*d0*d1);
	for(ITER i=0;i<d0*d1;i++)
		mat->data[i] = 0;

	return mat;
}
MAT* zeros_3d(UINT d0,UINT d1,UINT d2)
{
	MAT* mat = (MAT*)malloc(sizeof(MAT));
	mat->ndim = 2;
	mat->d0 = d0;
	mat->d1 = d1;
	mat->d2 = d2;
	
	mat->data=(DTYPE*)malloc(sizeof(DTYPE)*d0*d1*d2);
	for(ITER i=0;i<d0*d1*d2;i++)
		mat->data[i] = 0;

	return mat;
}


void free_MAT(MAT *mat)
{
	free(mat->data);
	free(mat);
}

void print_MAT(MAT* mat)
{
	switch(mat->ndim)
	{
		case 0:
			for(ITER i=0;i<mat->d0;i++)
				printf("%lf ",mat->data[i]);
			printf("\n");
		break;
		case 1:
		for(ITER j=0; j < mat->d1 ; j++)	
		{	for(ITER i=0;i < mat->d0;i++)
				printf("%lf ",mat->data[j*(mat->d0) + i]);
			printf("\n");
		}		
		break;
		case 2:
		for(ITER k=0; k < mat->d2 ; k++)
		{
			for(ITER j=0; j < mat->d1 ; j++)	
			{	for(ITER i=0;i < mat->d0;i++)
					printf("%lf ",mat->data[k*(mat->d2)*(mat->d1) + j*(mat->d1)  + j*i]);
				printf("\n");
			}
			printf("\n");		
		}
		break;
	}
}


