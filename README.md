# IIP speech preprocessing

C library for speech preprocessing.

## Requirement
cmake 3.10 or higher  
(OPTION)OpenBLAS  
(OPTION)intelMKL  
(OPTION)CUDA  
(OPTION)OpenMP  

## Installation

### CMAKE 3.10.3
+ ubuntu    
Use ubuntu\_cmake\_3\_10\_3\_installer.sh to install proper cmake  
+ Windows  
[32bit installer](https://cmake.org/files/v3.10/cmake-3.10.3-win64-x64.msi)   
[64bit installer](https://cmake.org/files/v3.10/cmake-3.10.3-win32-x86.msi)      

## Schedule

|date|content|
|---|---|
| 180629 | how to construct initial framework and first example |

## Coding style
snake\_case for most of code  
UPPERCASE for MACRO DATA TYPE, but snake\_case for macro overloaded function 

## Backend

- [openblas](https://github.com/gogyzzz/iip_sph_pp/wiki/OpenBLAS_%5Blinux%5D) (needed to be registered environment variable)
- [MKL](https://github.com/gogyzzz/iip_sph_pp/wiki/Intel_MKL_%EC%84%A4%EC%B9%98_%EB%B0%8F_%EC%82%AC%EC%9A%A9_%5BWINDOWS%5D) (needed to be registered environment variable)
- [cublas](https://developer.nvidia.com/cuda-downloads) (needed to be registered environment variable)
- standard c code with openMP 2.0 (cannot be higher than 2.0 for visual studio)

## FUNCTIONS
+ **Note** : Every function comes with complex, 1~3 Dimensions
+ math(matrix)
zeros - allocation
fill
set
get
submat
free
print
+ blas\_lv1
axpy : y = a\*x + y

+ blas\_lv2
genv : y = a\*A*x + b\*y

## Example 1 - Declaration matrix, Accessing and Operating

<details><summary>test_matrix.c</summary>

```c
#include "mother.h"

int main()
{
	
	MAT* a1, *sa1;
  MAT* a2, *sa2;
	MAT* a3, *sa3;

  a1 = zeros(4);
  a2 = zeros(3,4);
	a3 = zeros(2,3,4);

  sa1 = zeros(2);
  sa2 = zeros(2,2);
  sa3 = zeros(2,2,2);

	set(a1,2,1);
	print_MAT(a1);
	set(a2,1,2,1);
	print_MAT(a2);
	set(a3,0,1,2,1);
	print_MAT(a3);

	fill(a1,4);
	print_MAT(a1);
	fill(a2,5);
	print_MAT(a2);
	fill(a3,6);
	print_MAT(a3);
 
  DTYPE t1 = get(a1,0);
	printf("t1 : %lf\n",t1);	
 
  submat(a1,sa1,1,3);
  print_MAT(sa1);
  submat(a2,sa2,1,3,1,3);
  print_MAT(sa2);
  submat(a3,sa3,1,3,1,3,1,3);
  print_MAT(sa3);

  free_MAT(a1);
  free_MAT(a2);
	free_MAT(a3);
  free_MAT(sa1);
  free_MAT(sa2);
  free_MAT(sa3);
  
	CMAT* ca1, *sca1;
  CMAT* ca2, *sca2;
	CMAT* ca3, *sca3;

  ca1 = czeros(4);
	print_CMAT(ca1);
  ca2 = czeros(3,4);
	ca3 = czeros(2,3,4);

  sca1 = czeros(2);
  sca2 = czeros(2,2);
  sca3 = czeros(2,2,2);

	cset(ca1,2,1,1);
	cset(ca2,1,2,1,1);
	print_CMAT(ca2);
	cset(ca3,0,1,2,1,1);

	cfill(ca1,4,4);
	cfill(ca2,5,5);
	cfill(ca3,6,6);
	print_CMAT(ca3);
	
  csubmat(ca1,sca1,1,3);
  print_CMAT(sca1);
  csubmat(ca2,sca2,1,3,1,3);
  print_CMAT(sca2);
  csubmat(ca3,sca3,1,3,1,3,1,3);
  print_CMAT(sca3);

	free_CMAT(ca1);
  free_CMAT(ca2);
	free_CMAT(ca3);

	free_CMAT(sca1);
  free_CMAT(sca2);
	free_CMAT(sca3);
	return 0;
}



```
</details>

## Example 2 - BLAS LEVEL 1

<details><summary>test_blas_lv1.c</summary>

```C++

```
</details>

## EXample 3 - BLAS LEVEL 2

<details><summary>test_blas_lv2.c</summary>

```C
#include "mother.h"
//#include "header_for_test.h"

int main()
{

	MAT *A, *X, *Y;
	MAT *TA, *TX, *TY;

	CMAT *CA, *CX, *CY;
	CMAT *TCA, *TCX, *TCY;
	CTYPE calpha, cbeta;
#if USE_CUDA
  init();
#endif

  printf("%d\n",max_thread);
//	stopwatch(0);

	A = zeros(5, 4);
	X = alloc_MAT(5);
	Y = alloc_MAT(4);

	fill(A, 2);
	fill(X, 1);
	fill(Y, -1);

	print_MAT(A);
	print_MAT(X);
	print_MAT(Y);

	gemv(NoTran, 10, A, X, 3, Y);

	print_MAT(A);
	print_MAT(X);
	print_MAT(Y);

	free_MAT(A);
	free_MAT(X);
	free_MAT(Y);

	/**** TRANSPOSE  ****/

	TA = zeros(4, 5);
	TX = alloc_MAT(5);
	TY = alloc_MAT(4);

	fill(TA, 2);
	fill(TX, 1);
	fill(TY, -1);

	print_MAT(TA);
	print_MAT(TX);
	print_MAT(TY);

	gemv(Tran, 10, TA, TX, 3, TY);

	print_MAT(TA);
	print_MAT(TX);
	print_MAT(TY);

	free_MAT(TA);
	free_MAT(TX);
	free_MAT(TY);

	/**** COMPLEX ****/

	calpha.re = 5;
	calpha.im = 0;
	cbeta.re = 2;
	cbeta.im = 0;

	CA = czeros(4, 3);
	CX = alloc_CMAT(4);
	CY = alloc_CMAT(3);

	cfill(CA, 3, 0);
	cfill(CX, -2, 0);
	cfill(CY, 4, -2);

	print_CMAT(CA);
	print_CMAT(CX);
	print_CMAT(CY);

	cgemv(NoTran, calpha, CA, CX, cbeta, CY);

	print_CMAT(CA);
	print_CMAT(CX);
	print_CMAT(CY);

	free_CMAT(CA);
	free_CMAT(CX);
	free_CMAT(CY);

	/**** TRANSPOSE  ****/

	TCA = czeros(3, 4);
	TCX = alloc_CMAT(4);
	TCY = alloc_CMAT(3);

	cfill(TCA, 3, 0);
	cfill(TCX, -2, 0);
	cfill(TCY, 4, -2);

	print_CMAT(TCA);
	print_CMAT(TCX);
	print_CMAT(TCY);

	cgemv(Tran, calpha, TCA, TCX, cbeta, TCY);

	print_CMAT(TCA);
	print_CMAT(TCX);
	print_CMAT(TCY);

	free_CMAT(TCA);
	free_CMAT(TCX);
	free_CMAT(TCY);
#if USE_CUDA
  finit();
#endif
	return 0;
}
```
	
</details>

Result with DEBUG=1

<details><summary>blas_lv2.result</summary>

```
zeros_2d
alloc_MAT_1d
alloc_MAT_1d
fill
max thread : 1024
fill
max thread : 1024
fill
max thread : 1024
print_MAT
2.000 2.000 2.000 2.000 
2.000 2.000 2.000 2.000 
2.000 2.000 2.000 2.000 
2.000 2.000 2.000 2.000 
2.000 2.000 2.000 2.000 

print_MAT
1.000 
1.000 
1.000 
1.000 
1.000 

print_MAT
-1.000 
-1.000 
-1.000 
-1.000 

gemv
trans : 0 m : 4 n: 5 lda : 4
alpha : 10.000000 beta : 3.000000
print_MAT
2.000 2.000 2.000 2.000 
2.000 2.000 2.000 2.000 
2.000 2.000 2.000 2.000 
2.000 2.000 2.000 2.000 
2.000 2.000 2.000 2.000 

print_MAT
1.000 
1.000 
1.000 
1.000 
1.000 

print_MAT
97.000 
97.000 
97.000 
97.000 

free_MAT
free_MAT
free_MAT
zeros_2d
alloc_MAT_1d
alloc_MAT_1d
fill
max thread : 1024
fill
max thread : 1024
fill
max thread : 1024
print_MAT
2.000 2.000 2.000 2.000 2.000 
2.000 2.000 2.000 2.000 2.000 
2.000 2.000 2.000 2.000 2.000 
2.000 2.000 2.000 2.000 2.000 

print_MAT
1.000 
1.000 
1.000 
1.000 
1.000 

print_MAT
-1.000 
-1.000 
-1.000 
-1.000 

gemv
trans : 1 m : 5 n: 4 lda : 5
alpha : 10.000000 beta : 3.000000
print_MAT
2.000 2.000 2.000 2.000 2.000 
2.000 2.000 2.000 2.000 2.000 
2.000 2.000 2.000 2.000 2.000 
2.000 2.000 2.000 2.000 2.000 

print_MAT
1.000 
1.000 
1.000 
1.000 
1.000 

print_MAT
97.000 
97.000 
97.000 
97.000 

free_MAT
free_MAT
free_MAT
czeros_2d
alloc_CMAT_1d
alloc_CMAT_1d
cfill
cfill
cfill
print_CMAT
3.000 0.000|3.000 0.000|3.000 0.000|
3.000 0.000|3.000 0.000|3.000 0.000|
3.000 0.000|3.000 0.000|3.000 0.000|
3.000 0.000|3.000 0.000|3.000 0.000|

print_CMAT
-2.000 0.000|
-2.000 0.000|
-2.000 0.000|
-2.000 0.000|

print_CMAT
4.000 -2.000|
4.000 -2.000|
4.000 -2.000|

cgemv
trans : 0 m : 3 n: 4 lda : 3
alpha : 5.000000|0.000000 beta : 2.000000|0.000000
print_CMAT
3.000 0.000|3.000 0.000|3.000 0.000|
3.000 0.000|3.000 0.000|3.000 0.000|
3.000 0.000|3.000 0.000|3.000 0.000|
3.000 0.000|3.000 0.000|3.000 0.000|

print_CMAT
-2.000 0.000|
-2.000 0.000|
-2.000 0.000|
-2.000 0.000|

print_CMAT
-112.000 -4.000|
-112.000 -4.000|
-112.000 -4.000|

free_CMAT
free_CMAT
free_CMAT
czeros_2d
alloc_CMAT_1d
alloc_CMAT_1d
cfill
cfill
cfill
print_CMAT
3.000 0.000|3.000 0.000|3.000 0.000|3.000 0.000|
3.000 0.000|3.000 0.000|3.000 0.000|3.000 0.000|
3.000 0.000|3.000 0.000|3.000 0.000|3.000 0.000|

print_CMAT
-2.000 0.000|
-2.000 0.000|
-2.000 0.000|
-2.000 0.000|

print_CMAT
4.000 -2.000|
4.000 -2.000|
4.000 -2.000|

cgemv
trans : 1 m : 4 n: 3 lda : 4
alpha : 5.000000|0.000000 beta : 2.000000|0.000000
print_CMAT
3.000 0.000|3.000 0.000|3.000 0.000|3.000 0.000|
3.000 0.000|3.000 0.000|3.000 0.000|3.000 0.000|
3.000 0.000|3.000 0.000|3.000 0.000|3.000 0.000|

print_CMAT
-2.000 0.000|
-2.000 0.000|
-2.000 0.000|
-2.000 0.000|

print_CMAT
-112.000 -4.000|
-112.000 -4.000|
-112.000 -4.000|

free_CMAT
free_CMAT
free_CMAT

```

</details>


## References

[Discussion with SY](https://docs.google.com/document/d/1rCuWjxcCX7lz-VraY0BthAHz8QdSYxDVFVWy7HIMcDo/edit)

[Open audio library study document](https://github.com/kooBH/OpenAudioLibraryStudy)

## Log

## License



