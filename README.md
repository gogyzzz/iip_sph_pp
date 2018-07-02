# IIP speech preprocessing

C library for speech preprocessing.

## TODO
+ OpenMP
+ Complex 

## Requirement
cmake 3.10 or higher

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

## Backend

- [openblas]() (needed to be registered environment variable)
- [cublas]() (needed to be registered environment variable)
- standard c code with openMP 2.0 (cannot be higher than 2.0 for visual studio)

## Example 1 - Declaration matrix, Accessing and Scalar addition

```c
#include <mother.h>

int main(void)
{

DTYPE alpha = 1.0;
DTYPE beta = 1.0; 

mat* a1 = zeros_1d(2); // 가칭, 매크로 함수로 오버로딩이 되면 이름을 같게 할 수 있다
mat* a2 = zeros_2d(2,2); // 가칭, 매크로 함수로 오버로딩이 되면 이름을 같게 할 수 있다

// print

set_1d(a1, 0, 0.8);
set_2d(a2,0,1,0.9); //매크로 함수로 오버로딩이 되면 이름을 같게 할 수 있다

// print

axpb_1d(alpha, a1, beta); // X = alpha*X + beta; // += 1.0
axpb_2d(alpha, a2, beta); // X = alpha*X + beta; // += 1.0

// print

free_mat(a1); // 가칭
free_mat(a2); // 가칭

return 0;
}
```

## References

[Discussion with SY](https://docs.google.com/document/d/1rCuWjxcCX7lz-VraY0BthAHz8QdSYxDVFVWy7HIMcDo/edit)

[Open audio library study document](https://github.com/kooBH/OpenAudioLibraryStudy)

## Log

## License
