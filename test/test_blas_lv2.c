#include "mother.h"
#include "header_for_test.h"

int main()
{

MAT*A,*X,*Y;
MAT*TA,*TX,*TY;

CMAT*CA,*CX,*CY;
CMAT*TCA,*TCX,*TCY;
CTYPE calpha,cbeta;

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


/**** COMPLEX ****/

calpha.re = 5;
calpha.im = 0;
cbeta.re = 2;
cbeta.im = 0;

CA = czeros(4,3);
CX = alloc_CMAT(4);
CY = alloc_CMAT(3);

cfill(CA,3,0);
cfill(CX,-2,0);
cfill(CY,4,-2);

print_CMAT(CA);
print_CMAT(CX);
print_CMAT(CY);

cgemv(NoTran,calpha,CA,CX,cbeta,CY);

print_CMAT(CA);
print_CMAT(CX);
print_CMAT(CY);

free_CMAT(CA);
free_CMAT(CX);
free_CMAT(CY);

/**** TRANSPOSE  ****/

TCA = czeros(3,4);
TCX = alloc_CMAT(4);
TCY = alloc_CMAT(3);

cfill(TCA,3,0);
cfill(TCX,-2,0);
cfill(TCY,4,-2);

print_CMAT(TCA);
print_CMAT(TCX);
print_CMAT(TCY);

cgemv(Tran,calpha,TCA,TCX,cbeta,TCY);

print_CMAT(TCA);
print_CMAT(TCX);
print_CMAT(TCY);

free_CMAT(TCA);
free_CMAT(TCX);
free_CMAT(TCY);


stopwatch(1);
return 0;
}
