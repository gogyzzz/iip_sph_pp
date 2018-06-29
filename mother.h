#include <stdlib.h>
#include <stdio.h>

#define DTYPE double
#define UINT int
#define ITER long

typedef struct MAT
{
	DTYPE* data;
	UINT ndim;
	UINT d0;
	UINT d1;
	UINT d2;
}MAT;



#define zeros_load(_1,_2,_3,_4,...) _4
#define zeros(...) zeros_load(__VA_ARGS__, zeros_3d,zeros_2d,zeros_1d)(__VA_ARGS__) 

MAT* zeros_1d(UINT);
MAT* zeros_2d(UINT,UINT);
MAT* zeros_3d(UINT,UINT,UINT);
void free_MAT(MAT*);

void print_MAT(MAT*);
