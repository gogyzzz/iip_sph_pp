# IIP speech preprocessing

C library for speech preprocessing.

## TODO
+ OpenMP
+ blas level 1
+ more matrix operation

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

<details><summary>TEST PART 1</summary>
```c
//foo();
```

</details>

## References

[Discussion with SY](https://docs.google.com/document/d/1rCuWjxcCX7lz-VraY0BthAHz8QdSYxDVFVWy7HIMcDo/edit)

[Open audio library study document](https://github.com/kooBH/OpenAudioLibraryStudy)

## Log

## License
