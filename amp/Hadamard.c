#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
// #include "../include/matrix.c"
// #include "../include/random_number.c"
// #include "../include/blas.c"
// #include "init_functions.c"
// #include "inspection.c"
// #include "../include/interleaver.c"

gsl_matrix* hadamard(int n)
{
	gsl_matrix* ans;
	int matrix_size = (int)pow(2,n);
	ans = gsl_matrix_calloc(matrix_size,matrix_size);
	if(n == 0)
	{
		gsl_matrix_set(ans,0,0,1);
	}
	else 
	{
		int div_size = matrix_size/2;

		gsl_matrix* pre_Hadamard;
		pre_Hadamard = gsl_matrix_calloc(div_size,div_size);
		pre_Hadamard = hadamard(n-1);
		int i,j;
		for(i = 0; i < matrix_size; ++i){
			for (j = 0; j < matrix_size; ++j)
			{
				if(i < div_size && j < div_size)gsl_matrix_set(ans,i,j,gsl_matrix_get(pre_Hadamard,i,j));
				else if(i < div_size && j >= div_size)gsl_matrix_set(ans,i,j,gsl_matrix_get(pre_Hadamard,i,j-div_size));
				else if(i >= div_size && j < div_size)gsl_matrix_set(ans,i,j,gsl_matrix_get(pre_Hadamard,i-div_size,j));
				else gsl_matrix_set(ans,i,j,-gsl_matrix_get(pre_Hadamard,i-div_size,j-div_size));
			}

		}
		free(pre_Hadamard);
	}
	return ans;
}

void GslMatrixRealInterleaver(int n,int m,gsl_matrix* x,gsl_matrix* y,int *P)
{
	int i,j;
	for(i = 0;i < n; i++)
	{
		for(j = 0; j < m; j++)
		{
 			gsl_matrix_set(y,i,j,gsl_matrix_get(x,P[i],j));
		}
	}
}

// int main(int argc, char *argv[])
// {
// 	int n = 4,i;
// 	int size = pow(2,n-1);
// 	gsl_matrix* original = gsl_matrix_calloc(size,size);
// 	gsl_matrix* interleaver = gsl_matrix_calloc(size,size);

// 	original = Hadamard(n);
//     PrintRealMatrix(stdout,size,size,original);
    
//     int *P = (int *)malloc(sizeof(int)*size);
//     RandomNumberInitialization (0);
//     Permutation(size,P);
//     GslMatrixRealInterleaver(size,size,original,interleaver,P);
    
//     PrintRealMatrix(stdout,size,size,interleaver);
//     for (i = 0; i < size; ++i)
//     {
//     	printf("%d\n",P[i]);
//     }
//     return 0;
// }

// int main()
// {
// 	// int Tpp = 64;
// 	int n = 0;
// 	int matrix_size = 1;
// 	while(1)
// 	{
// 		if(matrix_size >= Tp)
// 		{
// 			break;
// 		}
// 		matrix_size = matrix_size*2;
// 		++n;
// 	}
// 	printf("%d %d\n",n,(int)pow(2,n));
// }
