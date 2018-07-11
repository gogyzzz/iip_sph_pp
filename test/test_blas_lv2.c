#include "mother.h"
#include "header_for_test.h"

int main()
{

MAT*A,*X,*Y;
MAT*TA,*TX,*TY;

stopwatch(0);

A = zeros(5,4);
X = alloc_MAT(5);
Y = alloc_MAT(4);

fill(A,2);
fill(X,1);
fill(Y,-1);

print_MAT(A);
print_MAT(X);
print_MAT(Y);

gemv(NoTran,10,A,X,3,Y);

print_MAT(A);
print_MAT(X);
print_MAT(Y);

free_MAT(A);
free_MAT(X);
free_MAT(Y);

/**** TRANSPOSE  ****/

TA = zeros(4,5);
TX = alloc_MAT(5);
TY = alloc_MAT(4);

fill(TA,2);
fill(TX,1);
fill(TY,-1);

print_MAT(TA);
print_MAT(TX);
print_MAT(TY);

gemv(Tran,10,TA,TX,3,TY);

print_MAT(TA);
print_MAT(TX);
print_MAT(TY);

free_MAT(TA);
free_MAT(TX);
free_MAT(TY);


stopwatch(1);
return 0;
}
