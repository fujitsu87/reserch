#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include "../include/matrix.c"
#include "../include/random_number.c"
#include "../include/blas.c"
#include "../include/inverse.c"


int main(int argc, char *argv[])
{
    int n,k;
    gsl_matrix_complex *A=gsl_matrix_complex_calloc(2,2);
    gsl_matrix_complex *B=gsl_matrix_complex_calloc(2,2);
    gsl_matrix_complex *C=gsl_matrix_complex_calloc(2,2);
	
	gsl_complex z;z.dat[0]=1/sqrt(2);z.dat[1]=1/sqrt(2);
	gsl_matrix_complex_set(A,0,0,z);

	z.dat[0]=0;z.dat[1]=0;
	gsl_matrix_complex_set(A,0,1,z);

	z.dat[0]=0;z.dat[1]=0;
	gsl_matrix_complex_set(A,1,0,z);

	z.dat[0]=1/sqrt(2);z.dat[1]=1/sqrt(2);
	gsl_matrix_complex_set(A,1,1,z);
	
	PrintMatrix(stdout,2,2,A);
	
	InverseA(A,B);
	PrintMatrix(stdout,2,2,B);

	AB(A,B,C);
	PrintMatrix(stdout,2,2,C);

	GSLMatrixFree(A);
	GSLMatrixFree(B);
	GSLMatrixFree(C);
	
    return 0;
}