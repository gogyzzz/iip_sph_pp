#include "mother.h"

int main()
{
	MAT*A;
	CMAT*CA;

	A= zeros(5,4);
	CA = czeros(5,4);


	
	printf("==== round ====\n");
	randu(A,-10.,10.);
	crandu(CA,-10.,10.,-10.,10.);
	print_MAT(A);
	print_CMAT(CA);
	round_mat(A);
	round_cmat(CA);
	print_MAT(A);
	print_CMAT(CA);

	printf("==== floor ====\n");
	randu(A,-10.,10.);
	crandu(CA,-10.,10.,-10.,10.);
	print_MAT(A);
	print_CMAT(CA);
	floor_mat(A);
	floor_cmat(CA);
	print_MAT(A);
	print_CMAT(CA);

	printf("==== ceil ====\n");
	randu(A,-10.,10.);
	crandu(CA,-10.,10.,-10.,10.);
	print_MAT(A);
	print_CMAT(CA);
	ceil_mat(A);
	ceil_cmat(CA);
	print_MAT(A);
	print_CMAT(CA);

	
	printf("==== log ====\n");
	randu(A,-10.,10.);
	crandu(CA,-10.,10.,-10.,10.);
	print_MAT(A);
	print_CMAT(CA);
	log_mat(A);
	log_cmat(CA);
	print_MAT(A);
	print_CMAT(CA);

	printf("==== log2 ====\n");
	randu(A,-10.,10.);
	crandu(CA,-10.,10.,-10.,10.);
	print_MAT(A);
	print_CMAT(CA);
	log2_mat(A);
	log2_cmat(CA);
	print_MAT(A);
	print_CMAT(CA);

	printf("==== log10 ====\n");
	randu(A,-10.,10.);
	crandu(CA,-10.,10.,-10.,10.);
	print_MAT(A);
	print_CMAT(CA);
	log10_mat(A);
	log10_cmat(CA);
	print_MAT(A);
	print_CMAT(CA);

	printf("==== exp ====\n");
	randu(A,-10.,10.);
	crandu(CA,-10.,10.,-10.,10.);
	print_MAT(A);
	print_CMAT(CA);
	exp_mat(A);
	exp_cmat(CA);
	print_MAT(A);
	print_CMAT(CA);

	return 0;
}
