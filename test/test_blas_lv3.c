#include "mother.h"

#define _print 1
	
#define _A 1  //N N
#define _B 0  //T N
#define _C 0  //T T
#define _D 0  //N T

int main()
{
	MAT*A,*B,*C;

#if _A
	A = alloc_MAT(3,4);
	B = alloc_MAT(4,5);
	C = alloc_MAT(3,5);

  fill(A,1.0); 
  fill(B,3.0); 
  fill(C,-0.1); 

	#if _print
	print_MAT(A);	
	print_MAT(B);	
	print_MAT(C);	
	#endif

	gemm(NoTran,NoTran,1,A,B,0,C);

	#if _print
	print_MAT(C);	
	#endif
  free_MAT(A);
  free_MAT(B);
  free_MAT(C);
#endif

#if _B
	A = alloc_MAT(4,3);
	B = alloc_MAT(4,5);
	C = alloc_MAT(3,5);

  fill(A,1.0); 
  fill(B,3.0); 
  fill(C,-0.1); 

	#if _print
	print_MAT(A);	
	print_MAT(B);	
	print_MAT(C);	
	#endif

	gemm(Tran,NoTran,1,A,B,0,C);

	#if _print
	print_MAT(C);	
	#endif
  free_MAT(A);
  free_MAT(B);
  free_MAT(C);
#endif

#if _C
	A = alloc_MAT(4,3);
	B = alloc_MAT(5,4);
	C = alloc_MAT(3,5);

  fill(A,1.0); 
  fill(B,3.0); 
  fill(C,-0.1); 

	#if _print
	print_MAT(A);	
	print_MAT(B);	
	print_MAT(C);	
	#endif

	gemm(Tran,Tran,1,A,B,0,C);

	#if _print
	print_MAT(C);	
	#endif
  free_MAT(A);
  free_MAT(B);
  free_MAT(C);
#endif

#if _D
	A = alloc_MAT(3,4);
	B = alloc_MAT(5,4);
	C = alloc_MAT(3,5);

  fill(A,1.0); 
  fill(B,3.0); 
  fill(C,-0.1); 

	#if _print
	print_MAT(A);	
	print_MAT(B);	
	print_MAT(C);	
	#endif

	gemm(NoTran,Tran,1,A,B,0,C);

	#if _print
	print_MAT(C);	
	#endif
  free_MAT(A);
  free_MAT(B);
  free_MAT(C);
#endif


	return 0;
}
