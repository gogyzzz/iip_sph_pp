# IIP speech preprocessing

C library for speech preprocessing.

## TODO
+ blas level 1
+ more matrix operation
+ CUDA implement

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

- [openblas]() (needed to be registered environment variable)
- [cublas]() (needed to be registered environment variable)
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

## Example 1 - Declaration matrix, Accessing and Scalar addition

<details><summary>test\_matrix.c</summary>

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

	return 0;
}



```

</details>

## References

[Discussion with SY](https://docs.google.com/document/d/1rCuWjxcCX7lz-VraY0BthAHz8QdSYxDVFVWy7HIMcDo/edit)

[Open audio library study document](https://github.com/kooBH/OpenAudioLibraryStudy)

## Log

## License



